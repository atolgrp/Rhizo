(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(test isoperimetric-quotient-of-a-square-is-pi-over-4
  ;; a large axis-aligned square: area=n^2, perimeter=4n -> Q = 4*pi*n^2/(16n^2) = pi/4
  (is (< (abs (- (isoperimetric-quotient 100.0d0 40.0d0) (/ pi 4))) 1d-9)))

(test isoperimetric-quotient-is-scale-invariant-for-similar-shapes
  ;; a 2x bigger square (same shape) must give the same Q
  (is (< (abs (- (isoperimetric-quotient 100.0d0 40.0d0)
                 (isoperimetric-quotient 400.0d0 80.0d0)))
         1d-9)))

(test cell-shape-metrics-areas-sum-to-grid-size
  (let* ((grid (make-voronoi-grid 30 30 25 :seed 1))
         (metrics (cell-shape-metrics grid)))
    (is (= (* 30 30) (reduce #'+ (mapcar #'cadr metrics))))))

(test mean-isoperimetric-quotient-is-between-0-and-1
  (let ((q (mean-isoperimetric-quotient (make-voronoi-grid 60 60 100 :seed 2))))
    (is (< 0.0d0 q 1.0d0))))

(test isotropy-coefficient-of-variation-basics
  (is (= 0.0d0 (isotropy-coefficient-of-variation '((:n . 5) (:e . 5) (:s . 5) (:w . 5)))))
  (is (> (isotropy-coefficient-of-variation '((:n . 1) (:e . 10) (:s . 5) (:w . 5))) 0.0d0))
  (is (null (isotropy-coefficient-of-variation '((:n . 5) (:e . nil))))))
