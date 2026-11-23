(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun hand-built-weighted-graph ()
  "0 --2-- 1 --3-- 2 --1-- 3, plus a direct-but-expensive shortcut 0--10--3.
Shortest path from 0 to 3 must go via 1 and 2 (cost 6), not the direct
edge (cost 10) - a real check that Dijkstra explores, not just follows
the cheapest single edge."
  (let ((adjacency (make-hash-table :test 'eql)))
    (add-edge! adjacency 0 1 2.0d0)
    (add-edge! adjacency 1 2 3.0d0)
    (add-edge! adjacency 2 3 1.0d0)
    (add-edge! adjacency 0 3 10.0d0)
    (make-hash-graph-tessellation adjacency)))

(test weighted-arrival-times-hand-checked
  (let* ((tess (hand-built-weighted-graph))
         (distances (compute-weighted-arrival-times tess (list 0))))
    (is (= 0.0d0 (gethash 0 distances)))
    (is (= 2.0d0 (gethash 1 distances)))
    (is (= 5.0d0 (gethash 2 distances)))
    (is (= 6.0d0 (gethash 3 distances)) "should take the 0-1-2-3 path (cost 6), not the direct 0-3 edge (cost 10)")))

(test weighted-arrival-times-multiple-sources
  ;; with two sources, each cell's distance should be to whichever source is nearer
  (let* ((tess (hand-built-weighted-graph))
         (distances (compute-weighted-arrival-times tess (list 0 3))))
    (is (= 0.0d0 (gethash 0 distances)))
    (is (= 0.0d0 (gethash 3 distances)))
    (is (= 1.0d0 (gethash 2 distances)) "2 is closer to source 3 (cost 1) than to source 0 (cost 5)")
    (is (= 2.0d0 (gethash 1 distances)) "1 is closer to source 0 (cost 2) than to source 3 (cost 4 via 2)")))

(test weighted-diffusion-agrees-with-tick-based-front-propagation
  ;; On a uniform-weight (CA-Voronoi) tessellation, Dijkstra distance to a
  ;; cell must equal the tick at which STEP-FRONT-PROPAGATION! first
  ;; reaches it - the two diffusion modes describe the same process when
  ;; every edge costs 1.
  (let* ((grid (make-voronoi-grid 25 25 30 :seed 9))
         (tess (voronoi-grid-tessellation grid))
         (start (first (cell-ids tess)))
         (world (make-world tess :seed 9))
         (arrival-tick (make-hash-table :test 'eql)))
    (seed-molecules-at! world (list start))
    (setf (gethash start arrival-tick) 0)
    (dotimes (i (length (cell-ids tess)))
      (step-front-propagation! world)
      (loop for id in (cell-ids tess)
            when (and (front-active-p world id) (not (gethash id arrival-tick)))
              do (setf (gethash id arrival-tick) (world-tick world))))
    (let ((distances (compute-weighted-arrival-times tess (list start))))
      (dolist (id (cell-ids tess))
        (is (= (gethash id arrival-tick) (round (gethash id distances)))
            "cell ~a: tick-based=~a weighted=~a" id (gethash id arrival-tick) (gethash id distances))))))

(test weighted-arrival-times-on-gis-fixture
  ;; A-B-E and A-D-E are two equal-cost paths (both 3.5) by construction -
  ;; see gis-tests.lisp for the ground-truth layout/centroids this relies on.
  (let* ((regions (parse-geojson-regions *fixture-path*))
         (tess (make-gis-tessellation regions))
         (a-id (region-id (find "A" regions :key #'region-id :test #'string=)))
         (e-id (region-id (find "E" regions :key #'region-id :test #'string=)))
         (distances (compute-weighted-arrival-times tess (list a-id))))
    (is (< (abs (- (gethash e-id distances) 3.5d0)) 1d-6))))

(test weighted-arrival-times-with-a-freshly-constructed-equal-string-source
  ;; Same bug class as gis-tests.lisp's
  ;; LOOKUPS-WORK-WITH-A-FRESHLY-CONSTRUCTED-EQUAL-STRING-ID, isolated for
  ;; this function specifically: a caller naturally passes their own
  ;; string literal as a source id, not the exact object cl-json allocated
  ;; while parsing the fixture. DISTANCES must be keyed by content, or this
  ;; silently returns a near-empty result (only the source itself, nothing
  ;; reachable) instead of an error.
  (let* ((regions (parse-geojson-regions *fixture-path*))
         (tess (make-gis-tessellation regions))
         (distances (compute-weighted-arrival-times tess (list (copy-seq "A")))))
    (is (= 6 (hash-table-count distances)) "all 6 regions should be reachable from A")
    (is (< (abs (- (gethash (copy-seq "E") distances) 3.5d0)) 1d-6))))
