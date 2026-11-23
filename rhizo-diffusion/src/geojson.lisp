(in-package :rhizo-diffusion)

;;; Minimal GeoJSON ingestion: a FeatureCollection of Polygon/MultiPolygon
;;; features becomes a list of plain REGION structs (id, exterior ring,
;;; centroid). Interior rings (holes) are ignored - a deliberate
;;; simplification at water-district/census-tract scale, called out here
;;; rather than silently mishandled. GeoJSON is plain JSON, so this needs
;;; nothing more than cl-json - no GIS-specific dependency.

(defstruct region
  id        ; from the feature's "id" property if present, else a generated index
  ring      ; list of (x . y) conses, the exterior boundary (no closing duplicate)
  centroid) ; (x . y) cons

(defun coordinates->ring (coordinates)
  "COORDINATES is GeoJSON's list of [x y] pairs for one ring. Returns a
list of (x . y) conses with the closing point (equal to the first) dropped."
  (let ((points (mapcar (lambda (pt) (cons (first pt) (second pt))) coordinates)))
    (if (and (rest points) (equal (first points) (car (last points))))
        (butlast points)
        points)))

(defun polygon-signed-area (ring)
  "The shoelace-formula signed area of RING (a list of (x . y) conses) -
positive for a counter-clockwise ring, negative for clockwise. Factored out
so POLYGON-AREA and RING-CENTROID share one implementation rather than two
that could silently diverge."
  (let ((next (append (rest ring) (list (first ring))))
        (area 0.0d0))
    (loop for (x0 . y0) in ring
          for (x1 . y1) in next
          do (incf area (- (* x0 y1) (* x1 y0))))
    (/ area 2.0d0)))

(defun polygon-area (ring)
  "The (always nonnegative) area enclosed by RING, regardless of winding
direction - the mechanistic reaction-diffusion engine's per-cell area (a 2D
proxy for volume; see tessellation.lisp CELL-AREA)."
  (abs (polygon-signed-area ring)))

(defun ring-centroid (ring)
  "Centroid of a simple polygon ring (list of (x . y) conses) via the
standard signed-area centroid formula - correct for non-convex shapes,
unlike a plain average of vertices."
  (let* ((n (length ring))
         (next (append (rest ring) (list (first ring))))
         (area (polygon-signed-area ring))
         (cx 0.0d0) (cy 0.0d0))
    (loop for (x0 . y0) in ring
          for (x1 . y1) in next
          for cross = (- (* x0 y1) (* x1 y0))
          do (incf cx (* (+ x0 x1) cross))
             (incf cy (* (+ y0 y1) cross)))
    (if (zerop area)
        ;; degenerate ring (e.g. duplicate/collinear points) - fall back
        ;; to a plain vertex average rather than dividing by zero
        (cons (/ (reduce #'+ (mapcar #'car ring)) n)
              (/ (reduce #'+ (mapcar #'cdr ring)) n))
        (cons (/ cx (* 6.0d0 area)) (/ cy (* 6.0d0 area))))))

(defun geometry->rings (geometry)
  "Return the exterior ring of each polygon part in GEOMETRY (a decoded
GeoJSON geometry alist) - one ring for a Polygon, one per part for a
MultiPolygon. Only the first (exterior) ring of each polygon is used."
  (let ((type (cdr (assoc :type geometry)))
        (coordinates (cdr (assoc :coordinates geometry))))
    (cond
      ((string= type "Polygon") (list (coordinates->ring (first coordinates))))
      ((string= type "MultiPolygon")
       (mapcar (lambda (polygon) (coordinates->ring (first polygon))) coordinates))
      (t (error "geometry->rings: unsupported geometry type ~s" type)))))

(defun parse-geojson-regions (path)
  "Parse a GeoJSON FeatureCollection at PATH into a list of REGION structs.
A feature's id comes from its \"id\" property if present, else its index
in the collection. (A MultiPolygon feature with more than one part is a
known limitation: every part gets the same id, so parts collide in the
adjacency table rather than being distinguished - fine for the
single-polygon-per-feature data this project targets, not yet handled for
genuine multi-part features.)"
  (let* ((json (with-open-file (stream path) (cl-json:decode-json stream)))
         (features (cdr (assoc :features json))))
    (loop for feature in features
          for index from 0
          for props = (cdr (assoc :properties feature))
          for feature-id = (and props (cdr (assoc :id props)))
          append (loop for ring in (geometry->rings (cdr (assoc :geometry feature)))
                       collect (make-region :id (or feature-id index)
                                            :ring ring
                                            :centroid (ring-centroid ring))))))
