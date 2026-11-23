(in-package :rhizo-diffusion)

;;; Physical-sanity validation: does the simulation actually behave like a
;;; sensible spatial process, independent of the isotropy question already
;;; tested elsewhere?
;;;
;;; Two properties checked here:
;;;
;;;  1. STEP-FRONT-PROPAGATION! is a graph BFS - every tick, the front
;;;     advances exactly one hop everywhere at once. That is a constant-
;;;     velocity wave, not Fickian diffusion (which spreads as sqrt(t), not
;;;     t) - so the area covered should grow as tick^2, i.e. effective
;;;     radius should grow LINEARLY with tick. FRONT-GROWTH-CURVE measures
;;;     it; LINEAR-REGRESSION fits a line and reports R^2.
;;;
;;;  2. Graph-hop (or weighted-Dijkstra) distance between two cells should
;;;     track their true spatial (Euclidean) distance reasonably well - if
;;;     it didn't, "N hops/units away" wouldn't mean anything geographically
;;;     once real GIS polygons are involved. DISTANCE-CORRELATION checks
;;;     this by correlating graph distance against Euclidean centroid
;;;     distance across every reachable cell from a source.

(defun linear-regression (xs ys)
  "Ordinary least-squares fit of YS = SLOPE*XS + INTERCEPT. Returns
(VALUES SLOPE INTERCEPT R-SQUARED)."
  (let* ((mx (alexandria:mean xs)) (my (alexandria:mean ys))
         (sxy (reduce #'+ (mapcar (lambda (x y) (* (- x mx) (- y my))) xs ys)))
         (sxx (reduce #'+ (mapcar (lambda (x) (expt (- x mx) 2)) xs)))
         (syy (reduce #'+ (mapcar (lambda (y) (expt (- y my) 2)) ys))))
    (when (zerop sxx)
      (error "linear-regression: every x value is identical (~a) - no line can be fit through zero variance in x" mx))
    (let* ((slope (/ sxy sxx))
           (intercept (- my (* slope mx)))
           (r2 (if (zerop syy) 1.0d0 (/ (* sxy sxy) (* sxx syy)))))
      (values slope intercept r2))))

(defun pearson-correlation (xs ys)
  (nth-value 2 (linear-regression xs ys)))

(defun cell-centroids (grid)
  "Return a hash-table of cell-id -> (x . y) pixel centroid for every cell
in GRID, computed directly from the pixel assignment (mean x, mean y of
every pixel with that id) - the CA-grown-grid analog of GEOJSON.LISP's
RING-CENTROID for real polygons."
  (let ((sum-x (make-hash-table :test 'equal))
        (sum-y (make-hash-table :test 'equal))
        (count (make-hash-table :test 'equal))
        (width (voronoi-grid-width grid))
        (height (voronoi-grid-height grid)))
    (dotimes (x width)
      (dotimes (y height)
        (let ((id (voronoi-grid-id-at grid x y)))
          (incf (gethash id sum-x 0) x)
          (incf (gethash id sum-y 0) y)
          (incf (gethash id count 0)))))
    (let ((centroids (make-hash-table :test 'equal)))
      (loop for id being the hash-keys of count
            do (setf (gethash id centroids)
                     (cons (/ (float (gethash id sum-x) 1.0d0) (gethash id count))
                           (/ (float (gethash id sum-y) 1.0d0) (gethash id count)))))
      centroids)))

(defun euclidean-distance (p1 p2)
  (sqrt (+ (expt (- (car p1) (car p2)) 2) (expt (- (cdr p1) (cdr p2)) 2))))

(defun distance-correlation (distances centroid-of source-id)
  "DISTANCES is a hash-table of cell-id -> graph distance (ticks, or
Dijkstra weight) from some source, e.g. from RUN-FRONT-PROPAGATION-WITH-
METRICS-style tracking or COMPUTE-WEIGHTED-ARRIVAL-TIMES. CENTROID-OF is a
function cell-id -> (x . y). Returns (VALUES R R-SQUARED N) - the Pearson
correlation between graph distance and true Euclidean distance from
SOURCE-ID's centroid, across every cell in DISTANCES except the source
itself, and how many points that correlation is over."
  (let* ((source-centroid (funcall centroid-of source-id))
         (pairs (loop for id being the hash-keys of distances using (hash-value d)
                      unless (equal id source-id)
                        collect (cons d (euclidean-distance source-centroid (funcall centroid-of id))))))
    (let ((graph-ds (mapcar #'car pairs)) (euclid-ds (mapcar #'cdr pairs)))
      (multiple-value-bind (slope intercept r2) (linear-regression graph-ds euclid-ds)
        (declare (ignore slope intercept))
        (values (sqrt r2) r2 (length pairs))))))

(defun front-growth-curve (world &key (max-ticks 1000))
  "Run WORLD's front-propagation to completion (or MAX-TICKS), recording
the cumulative number of reached cells after every tick. Returns a list of
(tick . cumulative-cell-count), tick 0 included."
  (let ((curve (list (cons 0 (hash-table-count (world-frontier world))))))
    (loop until (or (zerop (hash-table-count (world-molecules world)))
                     (>= (world-tick world) max-ticks))
          do (step-front-propagation! world)
             (push (cons (world-tick world) (hash-table-count (world-frontier world))) curve))
    (nreverse curve)))

(defun front-growth-linearity (curve avg-cell-area &key (skip-first 2) (skip-last-fraction 1/10))
  "Convert a FRONT-GROWTH-CURVE (tick . cumulative-cell-count) into
effective radius (sqrt(area-covered/pi)) at each tick, and fit a line
against tick. Excludes the first SKIP-FIRST ticks (discrete/small-sample
noise right after seeding) and the last SKIP-LAST-FRACTION of ticks
(growth saturates as the front runs out of grid to cover, which is a
finite-canvas edge effect, not a property of the propagation itself).
Returns (VALUES SLOPE R-SQUARED) - SLOPE is the wavefront's effective
speed in pixels/tick; R-SQUARED close to 1 means the front is expanding at
a genuinely constant speed (a wave), not accelerating/decelerating."
  (let* ((n (length curve))
         (last-index (- n 1 (round (* n skip-last-fraction))))
         (usable (subseq curve skip-first (max (1+ skip-first) last-index)))
         (ticks (mapcar (lambda (p) (float (car p) 1.0d0)) usable))
         (radii (mapcar (lambda (p) (sqrt (/ (* (cdr p) avg-cell-area) pi))) usable)))
    (multiple-value-bind (slope intercept r2) (linear-regression ticks radii)
      (declare (ignore intercept))
      (values slope r2))))
