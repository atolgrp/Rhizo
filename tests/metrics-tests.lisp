(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(test compass-measurement-points-are-valid-cells
  (let* ((grid (make-voronoi-grid 100 100 1500 :seed 10))
         (tess (voronoi-grid-tessellation grid))
         (points (compass-measurement-points grid))
         (valid-ids (cell-ids tess)))
    (is (= 8 (length points)))
    (is (set-equal '(:n :s :e :w :ne :nw :se :sw) (mapcar #'car points)))
    (dolist (p points)
      (is (member (cdr p) valid-ids)))))

(test front-propagation-is-approximately-isotropic
  ;; The whole point of alternating 4-/8-connected growth is that the front
  ;; should reach the 4 cardinal and 4 diagonal measurement points at
  ;; comparable ticks, not have the diagonals lag far behind the way a
  ;; fixed-neighborhood CA grid would. This is an automated stand-in for
  ;; the manual CSV/gnuplot comparison done by hand in the old notebook.
  (let* ((grid (make-voronoi-grid 150 150 3000 :seed 42))
         (points (compass-measurement-points grid))
         (world (make-world (voronoi-grid-tessellation grid) :seed 43))
         (center-id (voronoi-grid-id-at grid 75 75)))
    (seed-molecules-at! world (list center-id))
    (let* ((arrivals (run-front-propagation-with-metrics world points :max-ticks 500))
           (ticks (mapcar #'cdr arrivals)))
      (is (every #'identity ticks) "every direction should be reached within max-ticks")
      (is (<= (/ (reduce #'max ticks) (max 1 (reduce #'min ticks))) 3.0)
          "arrival ticks spread too much across directions: ~a" arrivals))))
