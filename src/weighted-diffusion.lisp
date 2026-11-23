(in-package :rhizo-diffusion)

;;; Weighted shortest-path ("arrival time") diffusion via Dijkstra, for
;;; tessellations where edges have real, non-uniform weights (GIS polygons
;;; of varying size - see gis-tessellation.lisp). This is a different
;;; shape of computation from the CA-Voronoi tick-based front-propagation
;;; in diffusion.lisp (a one-shot shortest-path result over the whole
;;; tessellation, not incremental stepping), so it gets its own function
;;; rather than being bolted onto WORLD/STEP-FRONT-PROPAGATION!. Backed by
;;; cl-heap's priority queue rather than a hand-rolled heap.
;;;
;;; Over a uniform-weight tessellation (every CA-Voronoi grid, weight 1.0
;;; on every edge - see voronoi-ca.lisp), this is exactly equivalent to
;;; breadth-first search: the distance this returns for a cell equals the
;;; tick at which STEP-FRONT-PROPAGATION! would first reach it. See
;;; tests/weighted-diffusion-tests.lisp for that cross-check.

(defun compute-weighted-arrival-times (tessellation source-ids)
  "Dijkstra from SOURCE-IDS (a list of cell ids, all starting at distance
0) over TESSELLATION, using NEIGHBORS-OF + EDGE-WEIGHT. Returns a
hash-table of cell-id -> shortest cumulative distance from the nearest
source. Cells unreachable from every source are absent from the table."
  (let ((distances (make-hash-table :test 'equal))
        (finalized (make-hash-table :test 'equal))
        (queue (make-instance 'cl-heap:priority-queue)))
    (dolist (id source-ids)
      (setf (gethash id distances) 0.0d0)
      (cl-heap:enqueue queue id 0.0d0))
    (loop until (zerop (cl-heap:queue-size queue))
          do (let ((id (cl-heap:dequeue queue)))
               ;; the same id can be enqueued more than once (whenever a
               ;; shorter path to it is found before it's been finalized);
               ;; only the first, shortest-distance pop needs to relax its
               ;; edges - later pops of the same id are stale and skipped
               (unless (gethash id finalized)
                 (setf (gethash id finalized) t)
                 (let ((d (gethash id distances)))
                   (dolist (nb (neighbors-of tessellation id))
                     (let ((candidate (+ d (edge-weight tessellation id nb))))
                       (when (or (not (gethash nb distances))
                                 (< candidate (gethash nb distances)))
                         (setf (gethash nb distances) candidate)
                         (cl-heap:enqueue queue nb candidate))))))))
    distances))
