(in-package :rhizo-diffusion)

;;; Build a HASH-GRAPH-TESSELLATION from real polygons (REGIONs, from
;;; geojson.lisp) via shared-edge contiguity: two regions are adjacent iff
;;; they share a full boundary edge. Detected in one pass over every
;;; region's boundary edges (canonicalized and hashed), not by comparing
;;; every pair of polygons - O(total vertices), not O(regions^2). This is
;;; the same "contiguity" approach GIS tools use (e.g. PySAL's rook
;;; contiguity); it deliberately does NOT count polygons that only touch
;;; at a single corner point as adjacent.
;;;
;;; Edge weight is centroid-to-centroid Euclidean distance - a simple,
;;; standard spatial proxy for "distance the front must cross." A real
;;; hydrology/road-network distance would be more accurate for something
;;; like PFAS transport via water systems; this project doesn't have that
;;; data yet, so this is flagged as a simplification, not hidden as exact.

(defparameter *coordinate-rounding-digits* 6
  "Boundary coordinates are rounded to this many decimal digits before
being used as shared-edge keys, so floating-point noise between two
polygons' nominally-identical shared vertices doesn't hide a real
adjacency. Increase if real data has coarser precision than that; decrease
if it's dense enough that distinct-but-nearby vertices start colliding.")

(defun round-coordinate (x)
  (let ((scale (expt 10 *coordinate-rounding-digits*)))
    (/ (round (* x scale)) (float scale 1.0d0))))

(defun canonical-edge-key (p0 p1)
  "An order-independent key for the boundary edge between points P0 and
P1, so the same physical edge - walked in opposite directions by two
neighboring polygons' rings - hashes identically."
  (let ((a (cons (round-coordinate (car p0)) (round-coordinate (cdr p0))))
        (b (cons (round-coordinate (car p1)) (round-coordinate (cdr p1)))))
    (if (or (< (car a) (car b)) (and (= (car a) (car b)) (< (cdr a) (cdr b))))
        (cons a b)
        (cons b a))))

(defun region-edges (region)
  "Return a list of (EDGE-KEY . LENGTH) for every boundary edge of
REGION's ring - EDGE-KEY via CANONICAL-EDGE-KEY (order-independent, so the
same physical edge - walked in opposite directions by two neighboring
regions - hashes identically), LENGTH the edge's real Euclidean length,
needed for the mechanistic engine's finite-volume flux (see
tessellation.lisp EDGE-LENGTH)."
  (let ((ring (region-ring region)))
    (loop for p0 in ring
          for p1 in (append (rest ring) (list (first ring)))
          collect (cons (canonical-edge-key p0 p1) (point-distance p0 p1)))))

(defun centroid-distance (region-a region-b)
  (point-distance (region-centroid region-a) (region-centroid region-b)))

(defun build-gis-adjacency (regions)
  "Build adjacency (EDGE-WEIGHT = centroid distance), per-region area, and
per-edge length in one pass over REGIONS' boundaries. Returns (VALUES
ADJACENCY CELL-AREAS EDGE-LENGTHS) - MAKE-GIS-TESSELLATION feeds all three
into MAKE-HASH-GRAPH-TESSELLATION. Extra return values are backward
compatible with any caller that only wanted ADJACENCY, as before."
  (let ((adjacency (make-hash-table :test 'equal))
        (cell-areas (make-hash-table :test 'equal))
        (edge-lengths (make-hash-table :test 'equal))
        (edge->region-ids (make-hash-table :test 'equal))
        (edge->length (make-hash-table :test 'equal))
        (region-by-id (make-hash-table :test 'equal)))
    (dolist (region regions)
      (setf (gethash (region-id region) region-by-id) region)
      (setf (gethash (region-id region) cell-areas) (polygon-area (region-ring region)))
      (ensure-cell! adjacency (region-id region))
      (loop for (key . length) in (region-edges region)
            do (pushnew (region-id region) (gethash key edge->region-ids) :test 'equal)
               (setf (gethash key edge->length) length)))
    (loop for key being the hash-keys of edge->region-ids using (hash-value region-ids)
          when (> (length region-ids) 1)
            do (loop for (a . rest) on region-ids
                     do (dolist (b rest)
                          (add-edge! adjacency a b
                                     (centroid-distance (gethash a region-by-id)
                                                         (gethash b region-by-id)))
                          (set-edge-length! edge-lengths a b (gethash key edge->length)))))
    (values adjacency cell-areas edge-lengths)))

(defun make-gis-tessellation (regions)
  "Return a HASH-GRAPH-TESSELLATION over REGIONS with full finite-volume
geometry. Each region's value is taken to sit at its centroid, so
FACE-DISTANCE is centroid-to-centroid distance - the same number EDGE-WEIGHT
holds here, but read through its own accessor (see tessellation.lisp).
Assumed, not measured: on irregular polygons the centroid-to-centroid
segment is not in general perpendicular to the shared face, so the
two-point flux is only approximately consistent there."
  (let* ((by-id (alexandria:alist-hash-table
                 (mapcar (lambda (r) (cons (region-id r) r)) regions) :test 'equal))
         (centroids (positions-table (mapcar #'region-id regions)
                                     (lambda (id) (region-centroid (gethash id by-id))))))
    (multiple-value-bind (adjacency cell-areas edge-lengths) (build-gis-adjacency regions)
      (make-fv-tessellation adjacency :cell-areas cell-areas :edge-lengths edge-lengths
                                      :cell-positions centroids))))

(defun tessellation-from-geojson (path)
  "Parse the GeoJSON FeatureCollection at PATH and build a tessellation
from it via shared-edge contiguity. Returns (VALUES TESSELLATION REGIONS)
so callers can still look up each region's centroid/ring afterward."
  (let ((regions (parse-geojson-regions path)))
    (values (make-gis-tessellation regions) regions)))
