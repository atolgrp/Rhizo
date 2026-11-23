(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun brute-force-nearest-seed (x y seed-points)
  "Reference answer: index of the nearest seed by scanning all of them,
lowest index on ties."
  (let ((best 0))
    (loop for k from 1 below (length seed-points)
          for (sx . sy) = (aref seed-points k)
          for (bx . by) = (aref seed-points best)
          when (< (+ (expt (- x sx) 2) (expt (- y sy) 2))
                  (+ (expt (- x bx) 2) (expt (- y by) 2)))
            do (setf best k))
    best))

(test seed-bin-size-gives-about-one-seed-per-bin
  (is (= 10 (seed-bin-size 100 100 100)))
  (is (= 15 (seed-bin-size 600 600 1600)))
  (is (= 1 (seed-bin-size 3 3 50)) "never below one pixel")
  (is (= 2 (seed-bin-size 3 1 1)) "rounds up"))

(test seed-bins-hold-each-seed-once-in-index-order
  (let* ((seeds (vector '(0 . 0) '(9 . 9) '(1 . 2) '(5 . 0)))
         (bins (seed-bins seeds 5 10 10)))
    (is (equal '(2 2) (array-dimensions bins)))
    (is (equal '(0 2) (aref bins 0 0)))
    (is (equal '(3) (aref bins 1 0)))
    (is (equal '(1) (aref bins 1 1)))
    (is (null (aref bins 0 1)))))

(test ring-bins-is-the-square-ring-clipped-to-the-array
  (is (equal '((2 . 2)) (ring-bins 2 2 0 5 5)))
  (is (= 8 (length (ring-bins 2 2 1 5 5))))
  (is (= 16 (length (ring-bins 2 2 2 5 5))))
  (is (set-equal '((1 . 0) (0 . 1) (1 . 1)) (ring-bins 0 0 1 5 5) :test #'equal))
  (is (null (ring-bins 0 0 5 5 5)) "a ring wholly outside the array is empty"))

(test nearest-seed-index-matches-brute-force-everywhere
  ;; Every pixel of three random layouts, including a sparse one where the
  ;; ring search has to go several bins out.
  (dolist (spec '((37 23 12 1) (40 40 50 2) (31 17 3 3)))
    (destructuring-bind (width height n seed) spec
      (let* ((seeds (coerce (rhizo-diffusion::random-distinct-points
                             width height n (sb-ext:seed-random-state seed))
                            'simple-vector))
             (size (seed-bin-size width height n))
             (bins (seed-bins seeds size width height))
             (mismatches (loop for x below width
                               sum (loop for y below height
                                         count (/= (nearest-seed-index x y seeds bins size)
                                                   (brute-force-nearest-seed x y seeds))))))
        (is (zerop mismatches) "~a mismatched pixels on ~ax~a with ~a seeds"
            mismatches width height n)))))

(test nearest-seed-index-breaks-ties-to-the-lowest-index
  ;; Pixel (2, 0) is 2 from both seeds; whichever bin each lands in, seed 0 wins.
  (let ((seeds (vector '(4 . 0) '(0 . 0))))
    (dolist (size '(1 2 5))
      (let ((bins (seed-bins seeds size 5 1)))
        (is (= 0 (nearest-seed-index 2 0 seeds bins size)))))))

(test nearest-seed-grid-cells-contain-their-seeds-and-tile-the-grid
  (let* ((grid (make-nearest-seed-grid 50 40 20 :seed 5))
         (tess (voronoi-grid-tessellation grid))
         (seeds (voronoi-grid-seed-points grid)))
    (is (= 20 (cell-count tess)))
    (is (= (* 50 40) (reduce #'+ (mapcar (lambda (id) (cell-area tess id)) (cell-ids tess)))))
    (loop for id from 0 for (x . y) across seeds
          do (is (= id (voronoi-grid-id-at grid x y)))
             (is (equal (cons (float x 1d0) (float y 1d0)) (cell-position tess id))))))

(test nearest-seed-grid-uses-the-same-seeds-as-the-ca-grid
  (is (equalp (voronoi-grid-seed-points (make-voronoi-grid 50 40 20 :seed 5))
              (voronoi-grid-seed-points (make-nearest-seed-grid 50 40 20 :seed 5)))))

(test nearest-seed-grid-is-reproducible-from-its-seed
  (let ((a (make-nearest-seed-grid 30 30 9 :seed 42))
        (b (make-nearest-seed-grid 30 30 9 :seed 42)))
    (is (= 42 (voronoi-grid-seed a)))
    (is (equalp (rhizo-diffusion::voronoi-grid-ids a) (rhizo-diffusion::voronoi-grid-ids b)))
    (is (integerp (voronoi-grid-seed (make-nearest-seed-grid 10 10 2)))
        "a seed is generated and recorded when none is given")))

(defun ids-array (columns)
  "A fixnum raster from a list of columns: (nth x COLUMNS) is column x,
listed from y = 0 upward."
  (make-array (list (length columns) (length (first columns)))
              :element-type 'fixnum :initial-contents columns))

(test boundary-step-counts-separates-x-and-y-steps
  ;;   y=1  0 1 1        cell 0 meets cell 1 across x twice ((0,1)|(1,1)
  ;;   y=0  0 0 1        and (1,0)|(2,0)) and across y once ((1,0)/(1,1))
  (let ((counts (boundary-step-counts (ids-array '((0 0) (0 1) (1 1))) 3 2)))
    (is (= 1 (hash-table-count counts)))
    (is (equal '(2 . 1) (gethash '(0 . 1) counts)))))

(test boundary-step-counts-ignores-corner-only-contact
  ;;   y=1  2 1          cells 0 and 1 touch only at the centre corner
  ;;   y=0  0 2
  (let ((counts (boundary-step-counts (ids-array '((0 2) (2 1))) 2 2)))
    (is (null (gethash '(0 . 1) counts)))
    (is (equal '(1 . 1) (gethash '(0 . 2) counts)))
    (is (equal '(1 . 1) (gethash '(1 . 2) counts)))))

(test face-length-from-steps-is-euclidean
  (is (eql 5.0d0 (face-length-from-steps 3 4)))
  (is (eql 7.0d0 (face-length-from-steps 7 0)))
  (is (< (abs (- (* 9 (sqrt 2d0)) (face-length-from-steps 9 9))) 1d-12)))

(test nearest-seed-tessellation-drops-corner-only-neighbours
  (let ((tess (nearest-seed-tessellation (ids-array '((0 2) (2 1)))
                                         (vector '(0 . 0) '(1 . 1) '(1 . 0)) 2 2)))
    (is (set-equal '(2) (neighbors-of tess 0)))
    (is (set-equal '(2) (neighbors-of tess 1)))
    (is (set-equal '(0 1) (neighbors-of tess 2)))
    (is (= (sqrt 2d0) (edge-length tess 0 2)))))

(test nearest-seed-tessellation-45-degree-face
  ;; Seeds (2,2) and (7,7) on 10x10: the bisector is x + y = 9, and pixels
  ;; on it are tied, so they go to seed 0. The face is a 9 x 9 staircase:
  ;; 18 pixel edges, true length 9*sqrt(2). Cell 0 has the 55 pixels with
  ;; x + y <= 9.
  (let* ((seeds (vector '(2 . 2) '(7 . 7)))
         (tess (nearest-seed-tessellation (nearest-seed-ids 10 10 seeds) seeds 10 10)))
    (is (< (abs (- (* 9 (sqrt 2d0)) (edge-length tess 0 1))) 1d-12))
    (is (< (abs (- (* 5 (sqrt 2d0)) (face-distance tess 0 1))) 1d-12))
    (is (= (face-distance tess 0 1) (edge-weight tess 0 1)) "hop cost is seed distance here")
    (is (< (abs (- 1.8d0 (face-conductance tess 0 1))) 1d-12))
    (is (= 55 (cell-area tess 0)))
    (is (= 45 (cell-area tess 1)))))

(test nearest-seed-tessellation-axis-aligned-face-is-unchanged
  ;; A vertical face is the one case the staircase count already got right.
  (let* ((seeds (vector '(2 . 5) '(7 . 5)))
         (tess (nearest-seed-tessellation (nearest-seed-ids 10 10 seeds) seeds 10 10)))
    (is (= 10.0d0 (edge-length tess 0 1)))
    (is (= 5.0d0 (face-distance tess 0 1)))))

(test nearest-seed-grid-has-symmetric-positive-face-geometry
  (let ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 60 60 25 :seed 3))))
    (dolist (a (cell-ids tess))
      (dolist (b (neighbors-of tess a))
        (is (plusp (edge-length tess a b)) "4-connected neighbours share a real face")
        (is (= (edge-length tess a b) (edge-length tess b a)))
        (is (= (face-distance tess a b) (face-distance tess b a)))))))

(test nearest-seed-grid-pixels-are-closer-to-their-own-seed
  ;; The defining property, checked on the finished grid: no pixel is
  ;; strictly closer to another seed than to its own.
  (let* ((grid (make-nearest-seed-grid 45 35 15 :seed 9))
         (seeds (voronoi-grid-seed-points grid)))
    (is (zerop (loop for x below 45
                     sum (loop for y below 35
                               count (/= (voronoi-grid-id-at grid x y)
                                         (brute-force-nearest-seed x y seeds))))))))
