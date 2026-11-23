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
