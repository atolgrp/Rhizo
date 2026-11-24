(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Calling distance and its sensitivity (PLAN.md, Phase 2 gate).

(test bessel-k1-matches-tabulated-values
  (loop for (x k1) in '((0.1d0 9.853844780870606d0)
                        (1d0 0.6019072301972346d0)
                        (2d0 0.13986588181652243d0)
                        (5d0 0.004044613445452164d0))
        do (is (< (relative-difference (bessel-k1 x) k1) 1d-12)
               "K1(~a): expected ~a, got ~a" x k1 (bessel-k1 x)))
  (signals error (bessel-k1 0d0)))

(test steady-calling-distance-is-where-the-k0-profile-meets-the-threshold
  (let* ((s 3d4) (d 490d0) (k 1d-3) (threshold 70d0)
         (r (steady-calling-distance s d k threshold)))
    (is (< (relative-difference (point-source-steady-state s d k r) threshold) 1d-10))))

(test steady-calling-distance-elasticities-match-finite-differences
  (let* ((s 3d4) (d 490d0) (k 1d-3) (threshold 70d0)
         (e (steady-calling-distance-elasticities s d k threshold))
         (f 1.001d0))
    (flet ((fd (fn) (/ (log (/ (funcall fn f) (funcall fn (/ f)))) (* 2 (log f)))))
      (is (< (abs (- (getf e :source-rate)
                     (fd (lambda (x) (steady-calling-distance (* s x) d k threshold)))))
             1d-5))
      (is (< (abs (- (getf e :threshold)
                     (fd (lambda (x) (steady-calling-distance s d k (* threshold x))))))
             1d-5))
      (is (< (abs (- (getf e :diffusivity)
                     (fd (lambda (x) (steady-calling-distance s (* d x) k threshold)))))
             1d-5))
      (is (< (abs (- (getf e :decay-rate)
                     (fd (lambda (x) (steady-calling-distance s d (* k x) threshold)))))
             1d-5))
      (is (= -0.5d0 (getf e :decay-rate))))))

(test point-source-transient-reduces-to-known-limits
  (let ((s 2d0) (d 490d0) (r 30d0))
    ;; no decay: S / (4 pi D) E1(R r^2 / (4 D t)); E1(1) = 0.21938393439552027
    (dolist (retardation '(1d0 5d0))
      (let ((time (/ (* retardation r r) (* 4 d))))
        (is (< (relative-difference (point-source-transient s d 0d0 retardation r time)
                                    (* (/ s (* 4 pi d)) 0.21938393439552027d0))
               1d-9))))
    ;; long times: the K0 steady state, whatever the retardation
    (dolist (retardation '(1d0 5d0))
      (is (< (relative-difference (point-source-transient s d 0.5d0 retardation r 1d4)
                                  (point-source-steady-state s d 0.5d0 r))
             1d-9)))
    (is (= 0d0 (point-source-transient s d 0.5d0 1d0 r 1d-6)) "far ahead of the front")))

(test transient-calling-distance-meets-the-threshold-at-that-time
  (let* ((s 3d4) (d 490d0) (k 1d-3) (retardation 10d0) (time 3600d0) (threshold 70d0)
         (r (transient-calling-distance s d k retardation threshold time)))
    (is (< (relative-difference (point-source-transient s d k retardation r time) threshold) 1d-8))
    (is (< r (steady-calling-distance s d k threshold)) "the transient front lags the steady one")))

(test log-elasticity-of-a-power-law-is-its-exponent
  (is (< (abs (- 2.5d0 (log-elasticity (lambda (x) (expt x 2.5d0)) 3d0))) 1d-12))
  (is (null (log-elasticity #'identity 0d0)) "no elasticity about zero"))

(test radial-profile-crossing-interpolates-in-log-concentration
  ;; profile c = 100 exp(-r / 10): crossing 10 at r = 10 ln 10
  (let ((profile (loop for r from 1d0 to 60d0 by 2d0 collect (cons r (* 100 (exp (- (/ r 10))))))))
    (is (< (abs (- (profile-crossing profile 10d0) (* 10 (log 10d0)))) 1d-9))
    (is (null (profile-crossing profile 1d-9)) "never crosses inside the profile")))

(test measured-calling-distance-matches-the-analytic-one
  ;; A single producer cell on 1,600 cells (spacing 10 um), decay length
  ;; 30 um, no sorption; the threshold is the analytic concentration at
  ;; 60 um, so the calling distance should come out at 60 um.
  (let* ((parameters (k0-gate-parameters 30d0 :k-prime 0d0))
         (tess (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1)))
         (density 5d14)
         (source (nearest-cell tess '(200d0 . 200d0)))
         (s-rate (* (cell-area tess source) density (parameter-value parameters :basal-production)))
         (k (parameter-value parameters :hydrolysis-rate))
         (threshold (point-source-steady-state s-rate 490d0 k 60d0))
         (world (make-ahl-world tess :parameters parameters :dt (/ 0.5d0 k)
                                     :initial-state (lambda (id)
                                                      (when (equal id source)
                                                        `((:producer . ,density)))))))
    (step-to-steady-state! world :ahl)
    (let ((r (calling-distance world :ahl (cell-position tess source) threshold)))
      (is (< (abs (- r 60d0)) 3d0) "measured calling distance ~a um, analytic 60 um" r))))

(test default-calling-scenario-reaches-beyond-the-patch
  ;; a calling distance inside the producer patch would be meaningless
  (let ((patch-radius (sqrt (/ (scenario-value nil :patch-area) pi))))
    (dolist (scenario '(() (:degrader-density 1d9)))
      (let ((result (calling-distance-scenario (ahl-literature-parameters) scenario)))
        (is (> (getf result :steady) (* 2 patch-radius)))
        (is (> (getf result :transient) (* 2 patch-radius)))))))

(test calling-distance-sensitivity-signs
  ;; The analytic table for the literature parameters.
  (flet ((rows (scenario) (calling-distance-sensitivity (ahl-literature-parameters)
                                                        :scenario scenario)))
    (let ((abiotic (rows '()))
          (biotic (rows '(:degrader-density 1d9))))
      (flet ((steady (rows name) (third (assoc name rows)))
             (transient (rows name) (fourth (assoc name rows))))
        ;; sorption cannot move the steady calling distance ...
        (is (< (abs (steady abiotic :kd)) 1d-9))
        ;; ... but with slow (abiotic) loss it governs the day-one one, which
        ;; then scales as sqrt(D t / R) - as strongly as time itself
        (is (< (transient abiotic :kd) -0.3d0))
        (is (< (abs (+ (transient abiotic :kd) (transient abiotic :time))) 0.05d0))
        ;; with fast biodegradation day one is near steady state already
        (is (< (abs (transient biotic :kd)) 0.1d0))
        ;; r ~ k^(-1/2) at steady state: hydrolysis alone, or shared
        (is (< (abs (+ 0.5d0 (steady abiotic :hydrolysis-rate))) 1d-6))
        (is (< (abs (+ 0.5d0 (steady biotic :hydrolysis-rate) (steady biotic :degrader-density)))
               1d-3))
        (is (null (steady abiotic :degrader-density)) "no elasticity about zero degraders")
        (dolist (rows (list abiotic biotic))
          (is (minusp (steady rows :receiver-threshold)))
          (is (plusp (steady rows :producer-density)))
          (is (plusp (transient rows :producer-density))))))))
