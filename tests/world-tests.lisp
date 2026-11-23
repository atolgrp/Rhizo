(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; These two guard exactly the bug classes found in the old code:
;;;  - a builder parameter silently ignored (old *voronoi-points* hardcoded
;;;    regardless of the `vp' argument)
;;;  - setup work redundantly duplicated (old make-voronoi-graph-circles
;;;    re-running initialize-voronoi-points on top of what the caller
;;;    already did, doubling *color-dict*/*neighbors-dict*)

(test voronoi-grid-honors-requested-point-count
  (dolist (n '(10 50 200))
    (let ((grid (make-voronoi-grid 40 40 n :seed 1)))
      (is (= n (cell-count (voronoi-grid-tessellation grid)))
          "requested ~a points but got ~a cells" n (cell-count (voronoi-grid-tessellation grid))))))

(test voronoi-grid-has-no-duplicate-or-missing-ids
  (let* ((n 80)
         (grid (make-voronoi-grid 40 40 n :seed 2))
         (ids (cell-ids (voronoi-grid-tessellation grid))))
    (is (= n (length ids)))
    (is (= n (length (remove-duplicates ids))))
    (is (every (lambda (id) (<= 0 id (1- n))) ids))))

(test voronoi-grid-is-fully-connected
  ;; Every cell should be reachable from every other cell via NEIGHBORS-OF -
  ;; a tessellation with an isolated or disconnected region would indicate a
  ;; broken growth/adjacency step.
  (let* ((grid (make-voronoi-grid 40 40 60 :seed 3))
         (tess (voronoi-grid-tessellation grid))
         (all-ids (cell-ids tess))
         (start (first all-ids))
         (visited (make-hash-table :test 'eql))
         (queue (list start)))
    (setf (gethash start visited) t)
    (loop while queue
          do (let ((id (pop queue)))
               (dolist (nb (neighbors-of tess id))
                 (unless (gethash nb visited)
                   (setf (gethash nb visited) t)
                   (push nb queue)))))
    (is (= (length all-ids) (hash-table-count visited)))))

(test voronoi-grid-deterministic-given-same-seed
  (let ((grid-a (make-voronoi-grid 30 30 40 :seed 99))
        (grid-b (make-voronoi-grid 30 30 40 :seed 99)))
    (dotimes (x 30)
      (dotimes (y 30)
        (is (= (voronoi-grid-id-at grid-a x y) (voronoi-grid-id-at grid-b x y)))))))

(test voronoi-grid-adjacency-is-symmetric
  ;; if B is a neighbor of A, A must be a neighbor of B
  (let* ((grid (make-voronoi-grid 30 30 40 :seed 4))
         (tess (voronoi-grid-tessellation grid)))
    (dolist (id (cell-ids tess))
      (dolist (nb (neighbors-of tess id))
        (is (member id (neighbors-of tess nb))
            "~a lists ~a as a neighbor, but not vice versa" id nb)))))

(test make-world-gives-each-world-its-own-random-state
  (let* ((grid (make-voronoi-grid 20 20 20 :seed 5))
         (tess (voronoi-grid-tessellation grid))
         (world-a (make-world tess :seed 111))
         (world-b (make-world tess :seed 111)))
    (is (= 111 (world-seed world-a) (world-seed world-b)))
    ;; same seed value, but must not be the SAME random-state object -
    ;; mutating one (by drawing from it) must never affect the other.
    (is (not (eq (world-random-state world-a) (world-random-state world-b))))))
