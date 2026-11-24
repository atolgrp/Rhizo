(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Rectangular grids and outer-boundary geometry (PLAN.md, Phase 5: the
;;; Richards benchmark is a 1D column, and flow needs boundary faces).

(test rectangular-tessellation-has-the-finite-volume-geometry
  ;; 3 columns x 2 rows of 2 x 1 cells; row-major ids, y grows downward
  (let ((tess (make-rectangular-tessellation 3 2 :dx 2d0 :dy 1d0)))
    (is (= 6 (cell-count tess)))
    (is (every (lambda (id) (= 2d0 (cell-area tess id))) (cell-ids tess)))
    (is (equal '(1d0 . 0.5d0) (cell-position tess 0)))
    (is (equal '(5d0 . 1.5d0) (cell-position tess 5)))
    (is (null (set-exclusive-or '(1 3) (neighbors-of tess 0))))
    (is (null (set-exclusive-or '(0 2 4) (neighbors-of tess 1))))
    (is (= 1d0 (edge-length tess 0 1)) "side-by-side cells share a vertical face of length dy")
    (is (= 2d0 (edge-length tess 0 3)) "stacked cells share a horizontal face of length dx")
    (is (= 2d0 (face-distance tess 0 1)))
    (is (= 1d0 (face-distance tess 0 3)))))

(test rectangular-boundary-faces-cover-each-side
  (let ((tess (make-rectangular-tessellation 3 2 :dx 2d0 :dy 1d0)))
    (is (equal '((0 2d0 0.5d0) (1 2d0 0.5d0) (2 2d0 0.5d0))
               (sort (copy-list (boundary-faces tess :top)) #'< :key #'first)))
    (is (equal '(3 4 5) (sort (mapcar #'first (boundary-faces tess :bottom)) #'<)))
    (is (equal '((0 1d0 1d0) (3 1d0 1d0))
               (sort (copy-list (boundary-faces tess :left)) #'< :key #'first)))
    (is (equal '(2 5) (sort (mapcar #'first (boundary-faces tess :right)) #'<)))
    (signals error (boundary-faces tess :front))))

(test column-tessellation-is-a-one-cell-wide-rectangular-grid
  (let ((tess (make-rectangular-tessellation 1 100 :dx 1d0 :dy 1d0)))
    (is (= 100 (cell-count tess)))
    (is (equal '((0 1d0 0.5d0)) (boundary-faces tess :top)))
    (is (equal '((99 1d0 0.5d0)) (boundary-faces tess :bottom)))
    (is (= 100 (length (boundary-faces tess :left))))))

(test every-cell-perimeter-is-accounted-for
  ;; interior face lengths plus boundary lengths add up to the perimeters
  (flet ((total-perimeter (tess)
           (+ (loop for id in (cell-ids tess)
                    sum (loop for nb in (neighbors-of tess id) sum (edge-length tess id nb)))
              (loop for side in '(:top :bottom :left :right)
                    sum (reduce #'+ (boundary-faces tess side) :key #'second)))))
    (let ((tess (make-rectangular-tessellation 4 3 :dx 1.5d0 :dy 0.5d0)))
      (is (< (abs (- (total-perimeter tess) (* 12 2 (+ 1.5d0 0.5d0)))) 1d-12)))))

(test nearest-seed-boundary-faces-follow-the-raster-edges
  (let* ((grid (make-nearest-seed-grid 60 40 30 :seed 7))
         (tess (voronoi-grid-tessellation grid))
         (ids (voronoi-grid-ids grid)))
    (flet ((side-length (side) (reduce #'+ (boundary-faces tess side) :key #'second)))
      (is (= 60d0 (side-length :top)))
      (is (= 60d0 (side-length :bottom)))
      (is (= 40d0 (side-length :left)))
      (is (= 40d0 (side-length :right))))
    (dolist (face (boundary-faces tess :top))
      (destructuring-bind (id length distance) face
        (is (= length (loop for x below 60 count (= id (aref ids x 0)))))
        (is (= distance (+ 0.5d0 (cdr (cell-position tess id))))
            "perpendicular distance from the seed to the edge y = -1/2")))
    (dolist (face (boundary-faces tess :right))
      (destructuring-bind (id length distance) face
        (declare (ignore length))
        (is (= distance (- 59.5d0 (car (cell-position tess id)))))))))

(test boundary-faces-need-recorded-geometry
  (signals error (boundary-faces (well-mixed-tessellation) :top)))
