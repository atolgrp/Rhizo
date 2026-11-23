(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Thresholds below are set with real margin under the empirically
;;; observed range (linearity R^2 0.95-0.995, distance-correlation r
;;; 0.93-0.98 across several seeds/grid sizes when this was calibrated) -
;;; loose enough not to be flaky on a different seed, tight enough to
;;; catch a genuine regression in either property.

(test linear-regression-basics
  (multiple-value-bind (slope intercept r2) (linear-regression '(1.0d0 2.0d0 3.0d0 4.0d0) '(2.0d0 4.0d0 6.0d0 8.0d0))
    (is (< (abs (- slope 2.0d0)) 1d-9))
    (is (< (abs intercept) 1d-9))
    (is (< (abs (- r2 1.0d0)) 1d-9))))

(test linear-regression-signals-a-clear-error-on-zero-x-variance
  ;; Found via a real bug: feeding all-identical x values used to crash
  ;; deep in a (/ 0.0d0 0.0d0) floating-point trap instead of failing
  ;; clearly. A caller doing this (e.g. every cell reached on the same
  ;; tick) should get an understandable error, not a cryptic FP trap.
  (signals error (linear-regression '(5.0d0 5.0d0 5.0d0) '(1.0d0 2.0d0 3.0d0))))

(defun build-test-world-and-track-arrivals (width height n-points seed)
  "Shared setup for the two physical-sanity checks below: build a grid,
seed the front at its center, run to completion, and return (VALUES CURVE
AVG-AREA TICK-OF-ARRIVAL CENTROIDS CENTER-ID)."
  (let* ((grid (make-voronoi-grid width height n-points :seed seed))
         (tess (voronoi-grid-tessellation grid))
         (avg-area (/ (float (* width height) 1.0d0) (cell-count tess)))
         (centroids (cell-centroids grid))
         (center-id (voronoi-grid-id-at grid (floor width 2) (floor height 2)))
         (world (make-world tess :seed seed))
         (tick-of-arrival (make-hash-table :test 'equal))
         (curve-points (list (cons 0 0))))
    (seed-molecules-at! world (list center-id))
    (setf (gethash center-id tick-of-arrival) 0)
    (setf (first curve-points) (cons 0 (hash-table-count (world-frontier world))))
    (loop until (zerop (hash-table-count (world-molecules world)))
          do (step-front-propagation! world)
             (push (cons (world-tick world) (hash-table-count (world-frontier world))) curve-points)
             (loop for id being the hash-keys of (world-molecules world)
                   unless (gethash id tick-of-arrival)
                     do (setf (gethash id tick-of-arrival) (world-tick world))))
    (values (nreverse curve-points) avg-area tick-of-arrival centroids center-id)))

(test front-propagation-is-a-constant-speed-wave
  ;; STEP-FRONT-PROPAGATION! is a graph BFS: the front advances exactly
  ;; one hop everywhere, every tick. That's a constant-velocity wave
  ;; (radius ~ tick), not Fickian diffusion (radius ~ sqrt(tick)) - if
  ;; this stopped being linear, the tessellation's connectivity would be
  ;; doing something structurally wrong (e.g. a strong directional bias
  ;; that speeds up/slows down growth as it goes).
  (multiple-value-bind (curve avg-area) (build-test-world-and-track-arrivals 100 100 1000 42)
    (multiple-value-bind (slope r2) (front-growth-linearity curve avg-area)
      (is (> r2 0.90) "wavefront should expand at a linear (constant-speed) rate, got R^2=~a" r2)
      (is (> slope 0) "wavefront speed should be positive"))))

(test graph-distance-correlates-with-euclidean-distance
  ;; The critical property for any spatial interpretation of arrival
  ;; times/hops: "N hops away" should mean roughly "N*(some constant)
  ;; pixels/meters away," not something arbitrary - otherwise results
  ;; wouldn't translate to real geography once GIS polygons are involved.
  (multiple-value-bind (curve avg-area tick-of-arrival centroids center-id)
      (build-test-world-and-track-arrivals 100 100 1000 42)
    (declare (ignore curve avg-area))
    (multiple-value-bind (r r2 n) (distance-correlation tick-of-arrival (lambda (id) (gethash id centroids)) center-id)
      (declare (ignore r2))
      (is (> n 500) "expected most of the ~1000 cells to be reachable and checked")
      (is (> r 0.85) "graph-hop distance should correlate strongly with true Euclidean distance, got r=~a" r))))
