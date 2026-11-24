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

(defun bessel-k (order x)
  "Return the modified Bessel function of the second kind K_ORDER(X), X >
0, from K_v(x) = integral over t in [0, inf) of exp(-x cosh t) cosh(v t), by
the trapezoidal rule with step 1/8. The integrand is analytic in the strip
|Im t| < pi/2 and decays doubly exponentially, so the rule's error is about
exp(-8 pi^2), far below double precision. Signals an error unless X > 0."
  (unless (plusp x)
    (error "BESSEL-K: argument must be positive, got ~a" x))
  (let ((x (float x 1d0)) (h 0.125d0))
    (* h (+ (* 0.5d0 (exp (- x)))
            (loop for i from 1
                  for term = (* (exp (- (* x (cosh (* i h))))) (cosh (* order i h)))
                  sum term into total of-type double-float
                  until (< term (* 1d-18 (+ total (exp (- x)))))
                  finally (return total))))))

(defun bessel-k0 (x)
  "Return K0(X), X > 0 (BESSEL-K of order 0)."
  (bessel-k 0 x))

(defun bessel-k1 (x)
  "Return K1(X) = -K0'(X), X > 0 (BESSEL-K of order 1)."
  (bessel-k 1 x))

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

(defun step-to-steady-state! (world species-name &key (tolerance 1d-10) (max-steps 100000)
                                                     stop-when)
  "Step WORLD with STEP-MECHANISTIC! until the relative change of
SPECIES-NAME over one whole step (SPECIES-COLUMN-CHANGE) falls below
TOLERANCE, until the predicate STOP-WHEN (a function of WORLD, if given) is
true after a step, or until MAX-STEPS steps have been taken. Returns
\(VALUES steps last-change outcome), OUTCOME being :STEADY, :STOPPED or
:MAX-STEPS - the last is not a steady state. Mutates WORLD."
  (let ((s (gethash species-name (mechanistic-world-species-index world)))
        (steps 0)
        (change 0d0)
        (outcome nil))
    (call-with-kernel
     (lambda ()
       (loop
         ;; a copy: substeps may update the state array in place
         (let ((before (state-column (mechanistic-world-state world) s)))
           (step-mechanistic! world)
           (incf steps)
           (setf change (species-column-change before (mechanistic-world-state world) s)))
         (setf outcome (cond ((and stop-when (funcall stop-when world)) :stopped)
                             ((< change tolerance) :steady)
                             ((>= steps max-steps) :max-steps)))
         (when outcome (return)))))
    (values steps change outcome)))

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

;;; Quorum induction of a uniform population (PLAN.md, Phase 2 gate).
;;;
;;; Producers at a uniform density rho have no gradients, so the whole
;;; domain behaves as one well-mixed cell: dC/dt = f(C), with
;;;
;;;     f(C) = rho (alpha + beta H(C)) - k C - kcat rho_d C / (Km + C),
;;;
;;; H the Hill activation. Steady states are the roots of f. Sorption that
;;; only exchanges with the dissolved pool slows the approach (by up to
;;; 1 + K') but does not move them. Solving f = 0 for the density gives
;;; rho(C) = loss(C) / (alpha + beta H(C)), the density at which C is a
;;; steady state. Where the feedback is strong enough, rho(C) rises to a
;;; local maximum and falls to a local minimum: a fold (saddle-node). Above
;;; the maximum, rho_on, the low (uninduced) state no longer exists and the
;;; population switches on; below the minimum, rho_off, the induced state
;;; no longer exists. Between them both are stable (hysteresis). Without a
;;; fold the response is graded, and the switch density is taken where the
;;; steady AHL concentration equals the induction threshold.

(defun uniform-ahl-balance (parameters producer-density c &key (degrader-density 0d0))
  "Return f(C), the net rate (nM/s) of AHL production minus hydrolysis and
biodegradation at concentration C in a well-mixed population of
PRODUCER-DENSITY producers and DEGRADER-DENSITY degraders (cells/L), with
PARAMETERS. See the comment above."
  (flet ((p (name) (parameter-value parameters name)))
    (- (* producer-density (+ (p :basal-production)
                              (* (p :induced-production)
                                 (hill-activation c (p :induction-threshold) (p :hill-coefficient)))))
       (* (p :hydrolysis-rate) c)
       (/ (* (p :biodegradation-kcat) degrader-density c) (+ (p :biodegradation-km) c)))))

(defun bisect-root (f lo hi &key (iterations 100))
  "Return a root of F in [LO, HI], where F(LO) and F(HI) differ in sign, by
ITERATIONS bisections (to double precision for the default)."
  (let ((f-lo (funcall f lo)))
    (loop repeat iterations
          for mid = (/ (+ lo hi) 2d0)
          for f-mid = (funcall f mid)
          do (if (eq (minusp f-mid) (minusp f-lo))
                 (setf lo mid f-lo f-mid)
                 (setf hi mid))
          finally (return (/ (+ lo hi) 2d0)))))

(defun uniform-steady-states (parameters producer-density &key (degrader-density 0d0)
                                                            (grid-points 4000))
  "Return the steady AHL concentrations, in increasing order, of a
well-mixed population (UNIFORM-AHL-BALANCE = 0) with PARAMETERS: one when
the response is monostable, three in the bistable range (low, unstable,
induced). Roots are bracketed on GRID-POINTS log-spaced concentrations up to
the bound rho (alpha + beta) / k, then bisected; two roots closer than one
grid step (0.7% at the default) are missed. Requires a hydrolysis rate > 0."
  (let* ((k (parameter-value parameters :hydrolysis-rate))
         (c-max (/ (* producer-density (+ (parameter-value parameters :basal-production)
                                          (parameter-value parameters :induced-production)))
                   k))
         (f (lambda (c) (uniform-ahl-balance parameters producer-density c
                                             :degrader-density degrader-density)))
         (grid (cons 0d0 (loop for i from 0 to grid-points
                               collect (* c-max (expt 10d0 (* -12d0 (- 1d0 (/ i grid-points)))))))))
    (unless (plusp k)
      (error "UNIFORM-STEADY-STATES: needs a positive hydrolysis rate"))
    (loop for (lo hi) on grid
          while hi
          when (not (eq (minusp (funcall f lo)) (minusp (funcall f hi))))
            collect (bisect-root f lo hi))))

(defun golden-section-maximum (f lo hi &key (tolerance 1d-12))
  "Return the X in [LO, HI] maximising F, assumed unimodal there, by
golden-section search to a bracket of TOLERANCE times (|LO| + |HI|). Near
a smooth maximum F changes with the square of the distance, so X is good to
about sqrt(double epsilon), 1.5e-8 relative, but F(X) to round-off."
  (let* ((ratio (/ (- (sqrt 5d0) 1d0) 2d0))
         (a lo) (b hi)
         (x1 (- b (* ratio (- b a)))) (x2 (+ a (* ratio (- b a))))
         (f1 (funcall f x1)) (f2 (funcall f x2)))
    (loop while (> (- b a) (* tolerance (+ (abs lo) (abs hi))))
          do (if (< f1 f2)
                 (setf a x1 x1 x2 f1 f2 x2 (+ a (* ratio (- b a))) f2 (funcall f x2))
                 (setf b x2 x2 x1 f2 f1 x1 (- b (* ratio (- b a))) f1 (funcall f x1))))
    (/ (+ a b) 2d0)))

(defun quorum-switch-density (parameters &key (degrader-density 0d0) (grid-points 2000))
  "Return the producer density (cells/L) at which a well-mixed population
switches to the induced state, from the steady-state density curve rho(C)
\(see the comment above), scanned on GRID-POINTS log-spaced concentrations
from 1e-6 to 1e6 times the induction threshold and refined by
GOLDEN-SECTION-MAXIMUM. Returns a plist:
  (:KIND :BISTABLE :ON rho_on :OFF rho_off :C-ON C at the upper fold
   :C-OFF C at the lower fold), or
  (:KIND :GRADED :HALF rho at which the steady AHL is the threshold)."
  (flet ((p (name) (parameter-value parameters name)))
    (let* ((threshold (p :induction-threshold))
           (rho (lambda (log-c)
                  (let ((c (exp log-c)))
                    (/ (- (uniform-ahl-balance parameters 0d0 c :degrader-density degrader-density))
                       (+ (p :basal-production)
                          (* (p :induced-production)
                             (hill-activation c threshold (p :hill-coefficient))))))))
           (logs (loop for i from 0 to grid-points
                       collect (+ (log (* 1d-6 threshold))
                                  (* (/ i grid-points) (log 1d12)))))
           (values (mapcar rho logs))
           (peak (loop for (l0 l1 l2) on logs
                       for (v0 v1 v2) on values
                       while l2
                       when (and (> v1 v0) (>= v1 v2)) return (list l0 l2)))
           (trough (and peak
                        (loop for (l0 l1 l2) on logs
                              for (v0 v1 v2) on values
                              while l2
                              when (and (> l1 (second peak)) (< v1 v0) (<= v1 v2))
                                return (list l0 l2)))))
      (if (and peak trough)
          (let ((log-on (golden-section-maximum rho (first peak) (second peak)))
                (log-off (golden-section-maximum (lambda (l) (- (funcall rho l)))
                                                 (first trough) (second trough))))
            (list :kind :bistable
                  :on (funcall rho log-on) :off (funcall rho log-off)
                  :c-on (exp log-on) :c-off (exp log-off)))
          (list :kind :graded :half (funcall rho (log threshold)))))))

(defun well-mixed-ahl-world (parameters producer-density degrader-density dt-factor)
  "Return an AHL world with PARAMETERS on WELL-MIXED-TESSELLATION holding
PRODUCER-DENSITY producers and DEGRADER-DENSITY degraders per litre and no
AHL, using the Strang scheme (RK4 reactions) with dt = DT-FACTOR / (the sum
of the hydrolysis, first-order biodegradation and sorption relaxation
rates). One cell has no diffusion, so Strang's steady-state bias (a
splitting artefact) cannot arise, and RK4 is stable at about 10 times the
explicit Euler step."
  (flet ((p (name) (parameter-value parameters name)))
    (make-ahl-world
     (well-mixed-tessellation)
     :parameters parameters
     :diffusion-scheme :implicit
     :dt (/ dt-factor (+ (p :hydrolysis-rate)
                         (/ (* (p :biodegradation-kcat) degrader-density) (p :biodegradation-km))
                         (sorption-relaxation-rate parameters)))
     :initial-state `((:producer . ,producer-density) (:degrader . ,degrader-density)))))

(defun well-mixed-ahl (world)
  "Return the AHL concentration of a world over WELL-MIXED-TESSELLATION."
  (aref (mechanistic-world-state world) 0 (gethash :ahl (mechanistic-world-species-index world))))

(defun uniform-steady-ahl (parameters producer-density &key (degrader-density 0d0)
                                                         (dt-factor 2d0) (tolerance 1d-12)
                                                         (max-steps 100000))
  "Return the steady AHL concentration a well-mixed population of
PRODUCER-DENSITY producers and DEGRADER-DENSITY degraders per litre reaches
from no AHL, simulated (WELL-MIXED-AHL-WORLD with DT-FACTOR) until AHL
changes by less than TOLERANCE in a step. Second value: the steps taken.
Signals an error if MAX-STEPS pass first."
  (let ((world (well-mixed-ahl-world parameters producer-density degrader-density dt-factor)))
    (multiple-value-bind (steps change outcome)
        (step-to-steady-state! world :ahl :tolerance tolerance :max-steps max-steps)
      (unless (eq outcome :steady)
        (error "UNIFORM-STEADY-AHL: not steady after ~d steps (last change ~,1e)" steps change))
      (values (well-mixed-ahl world) steps))))

(defun uniform-induction-p (parameters producer-density &key (degrader-density 0d0)
                                                          (dt-factor 2d0) (tolerance 1d-9)
                                                          (max-steps 1000000))
  "Return true if a simulated well-mixed population (WELL-MIXED-AHL-WORLD)
starting from no AHL becomes induced: its AHL passes the induction
threshold. The run stops as soon as it does, or when AHL changes by less
than TOLERANCE in a step below it (false): starting from no AHL, it climbs
toward the low steady state from below, and the default leaves it within
about 1e-6 of it for densities 1% from the fold. Near the fold both
outcomes are slow - the low state is approached at a rate falling as the
square root of the distance to the fold, and just above it the run crawls
through the bottleneck the fold leaves - so signals an error if MAX-STEPS
pass without a decision."
  (let ((world (well-mixed-ahl-world parameters producer-density degrader-density dt-factor))
        (threshold (parameter-value parameters :induction-threshold)))
    (multiple-value-bind (steps change outcome)
        (step-to-steady-state! world :ahl :tolerance tolerance :max-steps max-steps
                                          :stop-when (lambda (w) (> (well-mixed-ahl w) threshold)))
      (ecase outcome
        (:stopped t)
        (:steady (> (well-mixed-ahl world) threshold))
        (:max-steps
         (error "UNIFORM-INDUCTION-P: undecided after ~d steps (last change ~,1e)" steps change))))))

(defun simulated-switch-density (parameters &key lo hi (tolerance 1d-3) (degrader-density 0d0))
  "Return the producer density at which simulated well-mixed populations
switch on (UNIFORM-INDUCTION-P), found by bisection between LO (must stay
uninduced) and HI (must switch) until the bracket is within TOLERANCE of
LO. Signals an error if LO and HI do not bracket the switch."
  (flet ((induced-p (rho) (uniform-induction-p parameters rho :degrader-density degrader-density)))
    (when (or (induced-p lo) (not (induced-p hi)))
      (error "SIMULATED-SWITCH-DENSITY: ~a and ~a do not bracket the switch" lo hi))
    (loop while (> (- hi lo) (* tolerance lo))
          for mid = (/ (+ lo hi) 2d0)
          do (if (induced-p mid) (setf hi mid) (setf lo mid)))
    (/ (+ lo hi) 2d0)))
