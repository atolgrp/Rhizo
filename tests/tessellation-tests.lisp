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
    (signals error (edge-length tess 0 1))))

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
        (is (= 1.0d0 (edge-length tess (second pair) (first pair))))))))
