(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun tiny-square-adjacency ()
  "A hand-built 2x2 ring: 0-1-2-3-0, with distinct weights so EDGE-WEIGHT
can be checked precisely. Small enough to check everything by hand."
  (let ((adjacency (make-hash-table :test 'eql)))
    (add-edge! adjacency 0 1 1.0d0)
    (add-edge! adjacency 1 2 2.0d0)
    (add-edge! adjacency 2 3 3.0d0)
    (add-edge! adjacency 3 0 4.0d0)
    adjacency))

(test hash-graph-tessellation-basic-protocol
  (let ((tess (make-hash-graph-tessellation (tiny-square-adjacency))))
    (is (= 4 (cell-count tess)))
    (is (set-equal '(0 1 2 3) (cell-ids tess)))
    (is (set-equal '(1 3) (neighbors-of tess 0)))
    (is (set-equal '(0 2) (neighbors-of tess 1)))))

(test hash-graph-tessellation-edge-weight
  (let ((tess (make-hash-graph-tessellation (tiny-square-adjacency))))
    (is (= 1.0d0 (edge-weight tess 0 1)))
    (is (= 1.0d0 (edge-weight tess 1 0))) ; symmetric
    (is (= 3.0d0 (edge-weight tess 2 3)))
    (is (= 4.0d0 (edge-weight tess 3 0)))))

(test add-edge!-is-idempotent
  (let ((adjacency (make-hash-table :test 'eql)))
    (add-edge! adjacency 0 1 5.0d0)
    (add-edge! adjacency 0 1 999.0d0) ; must not overwrite or duplicate
    (add-edge! adjacency 1 0 999.0d0)
    (is (= 1 (length (gethash 0 adjacency))))
    (is (= 1 (length (gethash 1 adjacency))))
    (is (= 5.0d0 (cdr (assoc 1 (gethash 0 adjacency)))))))

(test set-edge-length!-accumulates-across-calls
  ;; Unlike ADD-EDGE!, repeat calls for an already-recorded pair ADD to the
  ;; length - needed for both real callers: the CA grid calls it once per
  ;; shared 4-connected pixel pair, and a pair of GIS regions touching along
  ;; more than one boundary stretch should get the true total.
  (let ((edge-lengths (make-hash-table :test 'equal)))
    (set-edge-length! edge-lengths 0 1 1.0d0)
    (set-edge-length! edge-lengths 0 1 1.0d0)
    (set-edge-length! edge-lengths 1 0 1.0d0)
    (is (= 3.0d0 (cdr (assoc 1 (gethash 0 edge-lengths)))))
    (is (= 3.0d0 (cdr (assoc 0 (gethash 1 edge-lengths))))))) ; symmetric

(test hash-graph-tessellation-without-geometry-signals-clear-errors
  ;; A tessellation built without CELL-AREAS/EDGE-LENGTHS - every topological
  ;; fixture in this file, and every caller before the mechanistic engine
  ;; needed real geometry - must fail loudly if CELL-AREA/EDGE-LENGTH are
  ;; called on it, not silently return some default that would hide a real
  ;; gap. Same "fail loud, not quiet" call already made for
  ;; LINEAR-REGRESSION's zero-variance case.
  (let ((tess (make-hash-graph-tessellation (tiny-square-adjacency))))
    (signals error (cell-area tess 0))
    (signals error (edge-length tess 0 1))
    (signals error (face-distance tess 0 1))))

(defun tiny-square-positions (id)
  "Corners of a 3 x 4 rectangle for the 0-1-2-3 ring: sides 3 and 4."
  (ecase id (0 '(0 . 0)) (1 '(3 . 0)) (2 '(3 . 4)) (3 '(0 . 4))))

(test point-distance-is-double-precision-for-integer-points
  ;; Pixel seeds are integer conses; SQRT of an integer would otherwise be
  ;; a single-float, losing about 8 digits before any flux is computed.
  (is (eql 5.0d0 (point-distance '(0 . 0) '(3 . 4))))
  (is (typep (point-distance '(0 . 0) '(1 . 1)) 'double-float))
  (is (= (sqrt 2.0d0) (point-distance '(0 . 0) '(1 . 1)))))

(test face-distances-from-positions-is-symmetric-euclidean
  (let ((distances (face-distances-from-positions (tiny-square-adjacency)
                                                  #'tiny-square-positions)))
    (flet ((d (a b) (cdr (assoc b (gethash a distances)))))
      (is (= 3.0d0 (d 0 1) (d 1 0)))
      (is (= 4.0d0 (d 1 2) (d 2 1)))
      (is (= 3.0d0 (d 2 3) (d 3 2)))
      (is (= 4.0d0 (d 3 0) (d 0 3)))
      (is (null (d 0 2)) "only adjacent pairs get a distance")
      (is (= 2 (length (gethash 0 distances))) "no duplicate entries"))))

(test face-distances-from-positions-rejects-coincident-positions
  (signals error (face-distances-from-positions (tiny-square-adjacency)
                                                (lambda (id) (declare (ignore id)) '(1 . 1)))))

(test face-distance-is-independent-of-edge-weight
  ;; The point of FACE-DISTANCE: the ring's EDGE-WEIGHTs are 1, 2, 3, 4,
  ;; but its face distances come from the rectangle's corners.
  (let ((tess (make-hash-graph-tessellation
               (tiny-square-adjacency)
               :face-distances (face-distances-from-positions (tiny-square-adjacency)
                                                              #'tiny-square-positions))))
    (is (= 1.0d0 (edge-weight tess 0 1)))
    (is (= 3.0d0 (face-distance tess 0 1)))
    (is (= 2.0d0 (edge-weight tess 1 2)))
    (is (= 4.0d0 (face-distance tess 1 2)))
    (is (= 4.0d0 (face-distance tess 2 1)))
    (signals error (face-distance tess 0 2) "non-adjacent pair has no face distance")))

(test hash-graph-tessellation-with-geometry
  (let ((tess (make-hash-graph-tessellation
               (tiny-square-adjacency)
               :cell-areas (alexandria:alist-hash-table '((0 . 5.0d0) (1 . 6.0d0)) :test 'equal)
               :edge-lengths (let ((h (make-hash-table :test 'equal)))
                               (set-edge-length! h 0 1 2.5d0)
                               h))))
    (is (= 5.0d0 (cell-area tess 0)))
    (is (= 6.0d0 (cell-area tess 1)))
    (is (= 2.5d0 (edge-length tess 0 1)))
    (is (= 2.5d0 (edge-length tess 1 0)))))

(test cell-index-basic
  (let* ((tess (make-hash-graph-tessellation (tiny-square-adjacency)))
         (index (make-cell-index tess)))
    (is (= 4 (cell-index-count index)))
    (is (set-equal '(0 1 2 3) (loop for i from 0 below 4 collect (index->cell-id index i))))
    (dotimes (i 4)
      (is (= i (cell-id->index index (index->cell-id index i)))))
    (signals error (cell-id->index index 999))))

(test ca-voronoi-grid-geometry-basic
  ;; width=3, height=1, n-points=3: every pixel is its own seed point, so
  ;; growth does nothing and the result is fully deterministic regardless of
  ;; random seed - three unit-area cells in a row, each pair of horizontal
  ;; neighbors sharing exactly one unit of boundary.
  (let* ((grid (make-voronoi-grid 3 1 3 :seed 1))
         (tess (voronoi-grid-tessellation grid))
         (ids (cell-ids tess)))
    (is (= 3 (length ids)))
    (is (= 3.0d0 (reduce #'+ (mapcar (lambda (id) (cell-area tess id)) ids))))
    (let ((pairs (remove-duplicates
                  (loop for a in ids append (loop for b in (neighbors-of tess a) collect (sort (list a b) #'<)))
                  :test #'equal)))
      (is (= 2 (length pairs)) "a 1x3 row of singleton cells should have exactly 2 adjacent pairs")
      (dolist (pair pairs)
        (is (= 1.0d0 (edge-length tess (first pair) (second pair))))
        (is (= 1.0d0 (edge-length tess (second pair) (first pair))))
        ;; every pixel is its own seed, so neighbouring seeds are 1 apart
        (is (= 1.0d0 (face-distance tess (first pair) (second pair))))))))

(test ca-voronoi-grid-keeps-uniform-edge-weight-and-gets-seed-distances
  ;; The CA-grown builder's EDGE-WEIGHT must stay 1.0 everywhere so the
  ;; isotropy study reproduces (CLAUDE.md, rule 1). FACE-DISTANCE is added
  ;; next to it: seed points are distinct integer pixels, so every face
  ;; distance is at least 1 and, on a 40x40 grid with 16 seeds, almost
  ;; all are well above 1.
  (let* ((tess (voronoi-grid-tessellation (make-voronoi-grid 40 40 16 :seed 11)))
         (pairs (loop for a in (cell-ids tess)
                      append (loop for b in (neighbors-of tess a) collect (cons a b))))
         (distances (mapcar (lambda (p) (face-distance tess (car p) (cdr p))) pairs)))
    (is (every (lambda (p) (= 1.0d0 (edge-weight tess (car p) (cdr p)))) pairs))
    (is (every (lambda (d) (>= d 1.0d0)) distances))
    (is (every (lambda (p) (= (face-distance tess (car p) (cdr p))
                              (face-distance tess (cdr p) (car p))))
               pairs))
    (is (> (/ (reduce #'+ distances) (length distances)) 3.0d0)
        "mean seed spacing on a 40x40 grid with 16 cells should be several pixels")))

(test ca-voronoi-growth-is-unchanged-since-the-isotropy-study
  ;; Fingerprint of one grown grid, measured on commit 8989ba1 (before any
  ;; Phase 1 change). CA growth and adjacency must stay bit-identical so the
  ;; published isotropy study reproduces (CLAUDE.md, rule 1).
  (let* ((grid (make-voronoi-grid 60 50 30 :seed 20261123))
         (tess (voronoi-grid-tessellation grid)))
    (is (= 58986553 (loop for x below 60
                          sum (loop for y below 50
                                    sum (* (voronoi-grid-id-at grid x y) (+ 1 x (* 60 y)))))))
    (is (= 71 (/ (loop for id in (cell-ids tess) sum (length (neighbors-of tess id))) 2)))))

(test ca-voronoi-grid-keeps-its-seed-points
  (let* ((grid (make-voronoi-grid 40 40 16 :seed 11))
         (tess (voronoi-grid-tessellation grid))
         (seeds (voronoi-grid-seed-points grid)))
    (is (typep seeds 'simple-vector))
    (is (= 16 (length seeds)))
    (loop for id from 0 for (x . y) across seeds
          do (is (= id (voronoi-grid-id-at grid x y)) "seed ~a must lie in its own cell" id)
             (is (equal (cons (float x 1d0) (float y 1d0)) (cell-position tess id))))
    ;; FACE-DISTANCE is exactly the distance between the two seeds
    (dolist (a (cell-ids tess))
      (dolist (b (neighbors-of tess a))
        (is (= (point-distance (aref seeds a) (aref seeds b)) (face-distance tess a b)))))))

(test positions-table-converts-to-double-floats
  (let ((table (positions-table '(0 1) (lambda (id) (if (= id 0) '(1 . 2) '(3/2 . 4))))))
    (is (equal '(1.0d0 . 2.0d0) (gethash 0 table)))
    (is (equal '(1.5d0 . 4.0d0) (gethash 1 table)))
    (is (= 2 (hash-table-count table)))))

(test cell-position-signals-clear-errors
  (let ((bare (make-hash-graph-tessellation (tiny-square-adjacency)))
        (partial (make-hash-graph-tessellation
                  (tiny-square-adjacency)
                  :cell-positions (positions-table '(0) #'tiny-square-positions))))
    (signals error (cell-position bare 0))
    (is (equal '(0.0d0 . 0.0d0) (cell-position partial 0)))
    (signals error (cell-position partial 1))))

(test nearest-cell-picks-the-closest-value-point
  (let ((tess (make-hash-graph-tessellation
               (tiny-square-adjacency)
               :cell-positions (positions-table '(0 1 2 3) #'tiny-square-positions))))
    (is (= 0 (nearest-cell tess '(0.1d0 . 0.2d0))))
    (is (= 1 (nearest-cell tess '(2.9d0 . -5.0d0))))
    (is (= 2 (nearest-cell tess '(10.0d0 . 10.0d0))))
    (is (= 3 (nearest-cell tess '(0.0d0 . 4.0d0))) "exactly on a value point")))

(test nearest-cell-on-a-ca-grid-finds-the-seed-cell
  ;; Standing on a seed, the nearest value point is that seed itself.
  (let* ((grid (make-voronoi-grid 40 40 16 :seed 11))
         (tess (voronoi-grid-tessellation grid)))
    (loop for id from 0 for seed across (voronoi-grid-seed-points grid)
          do (is (= id (nearest-cell tess seed))))))
