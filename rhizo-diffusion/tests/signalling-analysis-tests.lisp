(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Analytic references and measurements for the Phase 2 gate.

(test bessel-k0-matches-tabulated-values
  ;; Abramowitz & Stegun, Table 9.8 / standard tables
  (loop for (x k0) in '((0.1d0 2.4270690247020166d0)
                        (0.5d0 0.9244190712276659d0)
                        (1d0 0.42102443824070834d0)
                        (2d0 0.11389387274953344d0)
                        (5d0 0.0036910983340425942d0)
                        (10d0 1.778006231616917d-5))
        do (is (< (relative-difference (bessel-k0 x) k0) 1d-12)
               "K0(~a): expected ~a, got ~a" x k0 (bessel-k0 x)))
  ;; small-argument asymptote -ln(x/2) - Euler's gamma, error O(x^2 ln x)
  (is (< (abs (- (bessel-k0 1d-4) (- (- (log 5d-5)) 0.5772156649015329d0))) 1d-7))
  (signals error (bessel-k0 0d0)))

(test point-source-steady-state-is-the-k0-profile
  ;; S / (2 pi D) K0(r / lambda), lambda = sqrt(D / k): with S = 2 pi D and
  ;; r = lambda the value is K0(1)
  (let ((d 490d0) (k 0.5d0))
    (is (< (relative-difference (point-source-steady-state (* 2 pi d) d k (sqrt (/ d k)))
                                (bessel-k0 1d0))
           1d-14))
    (is (< (relative-difference (point-source-steady-state 3d0 d k 10d0)
                                (* (/ 3d0 (* 2 pi d)) (bessel-k0 (/ 10d0 (sqrt (/ d k))))))
           1d-14))))

(test signalling-worlds-default-to-imex
  ;; exact steady states (see STEP-MECHANISTIC!)
  (is (eq :imex (mechanistic-world-diffusion-scheme
                 (make-ahl-world (single-cell-tessellation 1d0) :dt 1d0))))
  (is (eq :imex (mechanistic-world-diffusion-scheme
                 (make-root-exudate-world (root-line-tessellation) #'root-cell-p :dt 1d0)))))

(defun k0-gate-parameters (decay-length &key (k-prime 4d0))
  "AHL literature parameters with the hydrolysis rate set so that the decay
length sqrt(D / k) is DECAY-LENGTH, sorption K' = K-PRIME exchanging at the
decay rate, and no induction (a constant source)."
  (let* ((d (parameter-value (ahl-literature-parameters) :diffusivity))
         (k (/ d (* decay-length decay-length))))
    (with-parameter-values (ahl-literature-parameters)
      `((:hydrolysis-rate . ,k) (:desorption-rate . ,k) (:induced-production . 0d0)
        (:kd . ,(/ (* k-prime 0.3d0) 1.3d0))))))

(test single-producer-steady-state-matches-the-k0-profile
  ;; PLAN.md Phase 2 gate. One producer cell at the centre of a 400 um
  ;; square (1,600 cells, spacing 10 um), decay length 30 um, compared on the
  ;; annulus 30-90 um (1 to 3 decay lengths; the nearest boundary image adds
  ;; under 0.1% there). Measured 2026-11-24, mesh seed 1: 4.5%. On 6,400
  ;; cells it is 1.7%; scripts/phase2-gate.lisp prints the refinement study.
  (let ((run (single-producer-steady-state
              (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1))
              :centre '(200d0 . 200d0) :parameters (k0-gate-parameters 30d0 :k-prime 1d0)
              :r-min 30d0 :r-max 90d0)))
    (is (< (getf run :l2) 0.05d0) "relative L2 error against K0: ~a" (getf run :l2))
    (is (< (getf run :residual) 1d-8) "not at steady state: last-step change ~a" (getf run :residual))
    (is (< (getf run :mass-balance) 1d-10) "pool balance drift ~a" (getf run :mass-balance))
    (is (< (abs (- (getf run :decay-length) 30d0)) 1d-9))
    (is (< 30 (getf run :n-compared)))))

(test species-column-change-is-relative-to-the-new-state
  (let ((after (make-array '(3 2) :element-type 'double-float
                                  :initial-contents '((1d0 0d0) (4d0 0d0) (2d0 0d0)))))
    (is (= 0.25d0 (species-column-change (dvec 1 3 2) after 0)))
    (is (= 0d0 (species-column-change (dvec 0 0 0) after 1)) "an all-zero column has not changed")))

(test step-to-steady-state-stops-at-the-fixed-point
  ;; Two-cell source world under :IMEX converges to c0 = 1.2, c1 = 0.8.
  (let ((world (two-cell-source-world :imex 1d0)))
    (multiple-value-bind (steps change outcome) (step-to-steady-state! world :c :tolerance 1d-12
                                                                                :max-steps 1000)
      (is (< steps 1000) "stopped after ~a steps, last change ~a" steps change)
      (is (< change 1d-12))
      (is (eq :steady outcome))
      (is (< (abs (- (cell-species-value world "0" :c) 1.2d0)) 1d-10)))))

(test step-to-steady-state-says-why-it-stopped
  (multiple-value-bind (steps change outcome)
      (step-to-steady-state! (two-cell-source-world :imex 1d0) :c :max-steps 3)
    (declare (ignore change))
    (is (= 3 steps))
    (is (eq :max-steps outcome) "running out of steps is not a steady state"))
  (let ((world (two-cell-source-world :imex 1d0)))
    (multiple-value-bind (steps change outcome)
        (step-to-steady-state! world :c :stop-when (lambda (w) (> (cell-species-value w "0" :c) 1d0)))
      (declare (ignore change))
      (is (eq :stopped outcome))
      (is (< steps 100))
      (is (> (cell-species-value world "0" :c) 1d0)))))
