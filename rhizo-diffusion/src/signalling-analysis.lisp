(in-package :rhizo-diffusion)

;;; Analysis of signalling runs (PLAN.md, Phase 2 gate): analytic
;;; references and the measurements compared with them.
;;;
;;; Steady state around a single producer. A source of strength S (amount
;;; per unit depth per second) in the plane, with diffusivity D and
;;; first-order loss k, settles to
;;;
;;;     C(r) = S / (2 pi D) K0(r / lambda),   lambda = sqrt(D / k),
;;;
;;; the 2D fundamental solution of D lap C - k C = -S delta. On the slab
;;; grid of signalling.lisp (uniform in depth) a producer patch is a source
;;; of this kind. Sorption that only exchanges with the dissolved pool does
;;; not change it: at steady state the sorbed pool is at equilibrium.

(defun bessel-k0 (x)
  "Return the modified Bessel function of the second kind K0(X), X > 0,
from K0(x) = integral over t in [0, inf) of exp(-x cosh t), by the
trapezoidal rule with step 1/8. The integrand is analytic in the strip
|Im t| < pi/2 and decays doubly exponentially, so the rule's error is about
exp(-8 pi^2), far below double precision. Signals an error unless X > 0."
  (unless (plusp x)
    (error "BESSEL-K0: argument must be positive, got ~a" x))
  (let ((x (float x 1d0)) (h 0.125d0))
    (* h (+ (* 0.5d0 (exp (- x)))
            (loop for i from 1
                  for term = (exp (- (* x (cosh (* i h)))))
                  sum term into total of-type double-float
                  until (< term (* 1d-18 (+ total (exp (- x)))))
                  finally (return total))))))

(defun point-source-steady-state (source-rate diffusivity decay-rate r)
  "Return the steady concentration at distance R from a source of
SOURCE-RATE (amount per unit depth per unit time) in the plane, with
DIFFUSIVITY and first-order DECAY-RATE: S / (2 pi D) K0(r / lambda), lambda
= sqrt(D / k). See the file header."
  (* (/ source-rate (* 2 pi diffusivity))
     (bessel-k0 (/ r (sqrt (/ diffusivity decay-rate))))))

(defun species-column-change (before after s)
  "Return max_i |AFTER[i,S] - BEFORE[i]| / max_i |AFTER[i,S]|, the relative
change of one species from BEFORE (a column, as from STATE-COLUMN) to
column S of the state AFTER; 0 if that column is all zero."
  (let ((change 0d0) (scale 0d0))
    (dotimes (i (length before))
      (setf change (max change (abs (- (aref after i s) (aref before i))))
            scale (max scale (abs (aref after i s)))))
    (if (zerop scale) 0d0 (/ change scale))))

(defun step-to-steady-state! (world species-name &key (tolerance 1d-10) (max-steps 100000))
  "Step WORLD with STEP-MECHANISTIC! until the relative change of
SPECIES-NAME over one whole step (SPECIES-COLUMN-CHANGE) falls below
TOLERANCE, or MAX-STEPS steps have been taken. Returns (VALUES steps
last-change). Mutates WORLD."
  (let ((s (gethash species-name (mechanistic-world-species-index world)))
        (steps 0)
        (change 0d0))
    (call-with-kernel
     (lambda ()
       (loop
         ;; a copy: substeps may update the state array in place
         (let ((before (state-column (mechanistic-world-state world) s)))
           (step-mechanistic! world)
           (incf steps)
           (setf change (species-column-change before (mechanistic-world-state world) s)))
         (when (or (< change tolerance) (>= steps max-steps))
           (return)))))
    (values steps change)))

(defun single-producer-steady-state (tessellation &key centre parameters (density 5d14)
                                                    r-min r-max (dt-factor 0.5d0)
                                                    (tolerance 1d-10))
  "Measure the steady AHL profile around a single producer (PLAN.md, Phase
2 gate). DENSITY producers per litre (default 5e14: cells of 1 um^3 filling
half the volume, as in Hense et al. 2012) occupy the cell of TESSELLATION
nearest CENTRE; no degraders. The AHL world with PARAMETERS (scheme :IMEX)
is stepped to steady state (STEP-TO-STEADY-STATE! with TOLERANCE) with dt =
DT-FACTOR / (hydrolysis rate + SORPTION-RELAXATION-RATE) - below 1 keeps the
explicit reactions positive, and :IMEX's steady state does not depend on dt
\(measured: the same error to 4 digits at 0.1 and 0.5) - and compared with
POINT-SOURCE-STEADY-STATE at every cell whose position is R-MIN to R-MAX
from the producer's. The source rate is the producer cell's area times its
production rate at the final concentration, so induction is allowed;
hydrolysis must be the only loss, as K0 assumes. Returns a plist:
  :L2           relative L2 error against K0 on the annulus
  :N-COMPARED   cells on the annulus
  :DECAY-LENGTH sqrt(D / k)
  :SOURCE-RATE  S, amount per unit depth per second
  :RESIDUAL     relative change of AHL in the last step
  :MASS-BALANCE drift of SIGNAL-POOL-BALANCE over produced AHL
  :STEPS        steps taken"
  (flet ((p (name) (parameter-value parameters name)))
    (let* ((source (nearest-cell tessellation centre))
           (origin (cell-position tessellation source))
           (d (p :diffusivity))
           (k (p :hydrolysis-rate))
           (world (make-ahl-world tessellation
                                  :parameters parameters
                                  :dt (/ dt-factor (+ k (sorption-relaxation-rate parameters)))
                                  :initial-state (lambda (id)
                                                   (when (equal id source)
                                                     `((:producer . ,density))))))
           (balance-0 (signal-pool-balance world :ahl)))
      (multiple-value-bind (steps residual) (step-to-steady-state! world :ahl :tolerance tolerance)
        (let* ((c-source (aref (mechanistic-world-state world)
                               (cell-id->index (mechanistic-world-cell-index world) source)
                               (gethash :ahl (mechanistic-world-species-index world))))
               (source-rate (* (cell-area tessellation source) density
                               (+ (p :basal-production)
                                  (* (p :induced-production)
                                     (hill-activation c-source (p :induction-threshold)
                                                      (p :hill-coefficient))))))
               (on-annulus (lambda (point) (<= r-min (point-distance point origin) r-max))))
          (list :l2 (relative-l2-error world :ahl
                                       (lambda (point)
                                         (point-source-steady-state source-rate d k
                                                                    (point-distance point origin)))
                                       :where on-annulus)
                :n-compared (count-if on-annulus (cell-ids tessellation)
                                      :key (lambda (id) (cell-position tessellation id)))
                :decay-length (sqrt (/ d k))
                :source-rate source-rate
                :residual residual
                :mass-balance (/ (abs (- (signal-pool-balance world :ahl) balance-0))
                                 (species-mass world :ahl-produced))
                :steps steps))))))
