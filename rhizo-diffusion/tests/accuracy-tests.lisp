(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun two-cell-world (values &key (area 1.0d0) (distance 4.0d0))
  "A diffusion-free world on TWO-CELL-TESSELLATION (cells at (0,0) and
\(DISTANCE,0)) with species :C set to the two VALUES."
  (make-mechanistic-world
   (two-cell-tessellation area 1.0d0 distance)
   :species (list (make-species :name :c :diffusivity 0d0)) :reactions '() :dt 1d0
   :initial-state (lambda (id) (list (cons :c (if (string= id "0") (first values) (second values)))))))

(defun mean-of (key plists)
  (/ (reduce #'+ plists :key (lambda (p) (getf p key))) (length plists)))

(test gaussian-field-is-the-normalised-heat-kernel
  (let ((g (gaussian-field '(1 . 2) 3.0d0)))
    (is (< (abs (- (funcall g '(1d0 . 2d0)) (/ 1 (* 2 pi 9)))) 1d-15))
    (is (< (abs (- (funcall g '(4d0 . 2d0)) (/ (exp -0.5d0) (* 2 pi 9)))) 1d-15))
    ;; unit integral: midpoint sum over +-8 widths at spacing 0.1
    (is (< (abs (- 1d0 (* 0.01d0 (loop for i from -240 below 240
                                       sum (loop for j from -240 below 240
                                                 sum (funcall g (cons (+ 1.05d0 (* 0.1d0 i))
                                                                      (+ 2.05d0 (* 0.1d0 j)))))))))
           1d-9))))

(test species-mass-and-moments-by-hand
  ;; areas 2, cells at (0,0) and (4,0), values 10 and 10
  (let ((world (two-cell-world '(10d0 10d0) :area 2d0)))
    (is (= 40d0 (species-mass world :c)))
    (is (= 8d0 (species-mean-squared-distance world :c '(0d0 . 0d0))))
    (is (= 4d0 (species-mean-squared-distance world :c '(2d0 . 0d0))))
    (signals error (species-mass world :no-such-species))))

(test relative-l2-error-by-hand
  ;; sqrt(((1-1)^2 + (3-1)^2) / (1 + 1)) = sqrt(2)
  (let ((world (two-cell-world '(1d0 3d0))))
    (is (= (sqrt 2d0) (relative-l2-error world :c (lambda (p) (declare (ignore p)) 1d0))))
    (is (= 0d0 (relative-l2-error world :c (lambda (p) (if (zerop (car p)) 1d0 3d0)))))))

(test observed-order-in-2d
  ;; 4x the cells halves the cell side
  (is (< (abs (- 2d0 (observed-order 0.4d0 0.1d0 4))) 1d-12))
  (is (< (abs (- 1d0 (observed-order 0.2d0 0.1d0 4))) 1d-12))
  (is (< (abs (- 2d0 (observed-order 0.9d0 0.1d0 9))) 1d-12)))

(test gaussian-spreading-error-reports-its-measurements
  (let ((result (gaussian-spreading-error
                 (voronoi-grid-tessellation (make-nearest-seed-grid 80 80 40 :seed 1))
                 :centre '(40d0 . 40d0) :sigma-start 4d0 :sigma-end 8d0)))
    (is (plusp (getf result :steps)))
    (is (< (getf result :mass-error) 1d-12))
    (is (< 0.5d0 (getf result :d-eff) 1.5d0))
    (is (< (getf result :l2) 1d0))))

(test fv-gaussian-accuracy-on-nearest-seed-grids
  ;; PLAN.md Phase 1, step 5 - formerly scripts/fv-accuracy.lisp. A Gaussian
  ;; spreads from width W/20 to W/10 with D = 1 on nearest-seed grids of
  ;; about 225 pixels per cell; 3 meshes per size, seeds 1-3, because single
  ;; meshes scatter (0.9%-3.5% error at 1,600 cells). Measured 2026-11-24:
  ;;   400 cells   9.1%  D_eff 0.998
  ;;   1600 cells  2.3%  D_eff 1.003
  (let ((coarse (grid-gaussian-errors #'make-nearest-seed-grid 300 400 '(1 2 3)))
        (fine (grid-gaussian-errors #'make-nearest-seed-grid 600 1600 '(1 2 3))))
    (dolist (run (append coarse fine))
      (is (< (getf run :mass-error) 1d-12)
          "mass must be conserved to round-off (seed ~a): ~a" (getf run :seed) (getf run :mass-error)))
    (is (< (mean-of :l2 fine) 0.03d0)
        "error at 1,600 cells should be under 3%: ~a" (mean-of :l2 fine))
    (let ((order (observed-order (mean-of :l2 coarse) (mean-of :l2 fine) 4)))
      (is (< 1.8d0 order 2.5d0) "convergence should be second order: observed ~a" order))
    (is (< (abs (- 1d0 (mean-of :d-eff fine))) 0.01d0)
        "D_eff at 1,600 cells should be within 1% of D = 1: ~a" (mean-of :d-eff fine))
    (is (< (abs (- 1d0 (mean-of :d-eff coarse))) 0.01d0)
        "D_eff at 400 cells should be within 1% of D = 1: ~a" (mean-of :d-eff coarse))))

;;; ------------------------------------------------ isotropy and retardation

(test pooled-species-mass-and-moment
  (let ((world (make-mechanistic-world
                (two-cell-tessellation 2.0d0 1.0d0 4.0d0)
                :species (list (make-species :name :a :diffusivity 0d0)
                               (make-species :name :b :diffusivity 0d0))
                :reactions '() :dt 1d0
                :initial-state (lambda (id) (if (string= id "0") '((:a . 3d0)) '((:b . 1d0)))))))
    ;; :a holds 6 at (0,0), :b holds 2 at (4,0)
    (is (= 8d0 (species-mass world '(:a :b))))
    (is (= 6d0 (species-mass world :a)))
    (is (= 4d0 (species-mean-squared-distance world '(:a :b) '(0d0 . 0d0)))) ; (6*0 + 2*16) / 8
    (is (= 0d0 (species-mean-squared-distance world :a '(0d0 . 0d0))))))

(test point-source-arrival-time-solves-for-the-threshold
  (let* ((r 10d0) (peak (/ 1 (* pi (exp 1d0) r r))))
    (is (< (abs (- peak (point-source-concentration 1d0 1d0 r 25d0))) 1d-15)
        "the peak at r = 10 is at t = r^2 / 4D = 25")
    (let ((arrival (point-source-arrival-time 1d0 1d0 r (* 0.5d0 peak))))
      (is (< arrival 25d0))
      (is (< (abs (- (point-source-concentration 1d0 1d0 r arrival) (* 0.5d0 peak))) (* 1d-9 peak))))
    (is (null (point-source-arrival-time 1d0 1d0 r (* 1.01d0 peak))) "above the peak: never")))

(test compass-points-are-equidistant
  (let ((points (compass-points '(5d0 . 7d0) 3d0)))
    (is (equal '(:e :n :w :s :ne :nw :sw :se) (mapcar #'car points)))
    (dolist (p points)
      (is (< (abs (- 3d0 (point-distance '(5d0 . 7d0) (cdr p)))) 1d-12)))
    (is (equal '(8d0 . 7d0) (cdr (assoc :e points))))))

(test anisotropy-index-by-hand
  (is (< (abs (- (/ 2d0 11d0)
                 (anisotropy-index '((:e . 10) (:n . 10) (:w . 10) (:s . 10)
                                     (:ne . 12) (:nw . 12) (:sw . 12) (:se . 12)))))
         1d-12))
  (is (zerop (anisotropy-index (mapcar (lambda (d) (cons (car d) 5)) *compass-directions*)))))

(test threshold-arrival-times-interpolates-the-crossing
  ;; Two cells, explicit: C1 = 5 (1 - (1 - 2 k dt)^n), k = 0.25, dt = 0.1.
  ;; Cell "0" starts above the threshold, so it arrives at time 0.
  (let* ((dt 0.1d0) (threshold 2d0)
         (world (make-mechanistic-world
                 (two-cell-tessellation 1.0d0 1.0d0 4.0d0)
                 :species (list (make-species :name :c :diffusivity 1d0)) :reactions '() :dt dt
                 :initial-state (lambda (id) (if (string= id "0") '((:c . 10d0)) '((:c . 0d0))))))
         (arrivals (threshold-arrival-times world :c '((:far . "1") (:source . "0")) threshold 1000))
         (c (lambda (n) (* 5 (- 1 (expt (- 1 (* 2 0.25d0 dt)) n)))))
         (n (loop for n from 1 when (>= (funcall c n) threshold) return n))
         (expected (* dt (+ (1- n) (/ (- threshold (funcall c (1- n)))
                                      (- (funcall c n) (funcall c (1- n))))))))
    (is (< (abs (- expected (cdr (assoc :far arrivals)))) 1d-12))
    (is (= 0d0 (cdr (assoc :source arrivals))))
    (is (null (cdr (assoc :never (threshold-arrival-times
                                   world :c '((:never . "1")) 100d0 5)))))))

(test diffusion-front-is-isotropic-on-nearest-seed-grids
  ;; Phase 1 gate. The geometry of the published isotropy study
  ;; (voronoi-ca-isotropy/README.md): 200x200, one seed per 10 px, probes at
  ;; r = 49, its seeds 20210728-31. That study's Voronoi CA front has an
  ;; anisotropy index of 1.62% at this radius. Here a point source diffuses
  ;; and each probe's threshold-crossing time is compared with the analytic
  ;; one at that probe's own distance. Measured 2026-11-24: index 0.81%
  ;; (2.95% raw, probe placement included), CoV 1.84%, arrival 0.986 of
  ;; analytic (the zero-flux boundary and a one-cell source: assumed causes).
  (let* ((runs (loop for seed from 20210728 to 20210731
                     collect (point-source-isotropy
                              (voronoi-grid-tessellation (make-nearest-seed-grid 200 200 4000 :seed seed))
                              :centre '(100d0 . 100d0) :radius 49d0 :scheme :explicit)))
         (normalised (mean-of :normalised runs)))
    (is (< normalised 0.0162d0) "anisotropy index ~,2f% should be under the published 1.62%"
        (* 100 normalised))
    (is (< (mean-of :raw runs) 0.05d0) "raw index ~,2f%" (* 100 (mean-of :raw runs)))
    (is (< (mean-of :cov runs) 0.03d0) "CoV ~,2f%" (* 100 (mean-of :cov runs)))
    (is (< (abs (- 1d0 (mean-of :mean-ratio runs))) 0.03d0)
        "mean arrival ~,4f of analytic" (mean-of :mean-ratio runs))))

(test retardation-matches-one-plus-partition-coefficient
  ;; Phase 1 gate: R = 1 + rho_b K_d / theta. With K' = rho_b K_d / theta = 4,
  ;; total mass must spread at D / 5. Measured 2026-11-24: R = 5.0010
  ;; explicit, 5.0017 implicit at 5x the explicit step; mass to 2e-15.
  (let ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 300 300 400 :seed 3))))
    (dolist (spec '((:explicit 1) (:implicit 5)))
      (let ((run (retardation-measurement tess 4d0 :centre '(150d0 . 150d0)
                                                   :sigma-start 15d0 :sigma-end 30d0
                                                   :scheme (first spec) :time-step-factor (second spec))))
        (is (= 5d0 (getf run :r-expected)))
        (is (< (abs (- 1 (/ (getf run :r-measured) 5d0))) 0.005d0)
            "~a: R = ~,4f, expected 5" (first spec) (getf run :r-measured))
        (is (< (getf run :mass-error) 1d-12) "~a: mass error ~a" (first spec) (getf run :mass-error))))))
