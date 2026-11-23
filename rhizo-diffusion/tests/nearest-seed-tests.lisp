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

(test nearest-seed-grid-pixels-are-closer-to-their-own-seed
  ;; The defining property, checked on the finished grid: no pixel is
  ;; strictly closer to another seed than to its own.
  (let* ((grid (make-nearest-seed-grid 45 35 15 :seed 9))
         (seeds (voronoi-grid-seed-points grid)))
    (is (zerop (loop for x below 45
                     sum (loop for y below 35
                               count (/= (voronoi-grid-id-at grid x y)
                                         (brute-force-nearest-seed x y seeds))))))))
