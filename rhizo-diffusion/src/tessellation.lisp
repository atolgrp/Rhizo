(in-package :rhizo-diffusion)

;;; The tessellation protocol.
;;;
;;; Diffusion and metrics code never touches a grid, a random point set,
;;; a polygon, or anything specific to how a tessellation was built - they
;;; only ever call these four functions. That is deliberate: it means a
;;; tessellation built by growing a CA-Voronoi diagram (VORONOI-GRID, in
;;; voronoi-ca.lisp) and one built from real GIS polygons (GIS-TESSELLATION.
;;; LISP) can be handed to the exact same diffusion/metrics code unchanged.

(defgeneric cell-ids (tessellation)
  (:documentation "Return a fresh list of every cell id in TESSELLATION."))

(defgeneric neighbors-of (tessellation cell-id)
  (:documentation "Return the list of ids adjacent to CELL-ID in TESSELLATION."))

(defgeneric cell-count (tessellation)
  (:documentation "Return the number of cells in TESSELLATION."))

(defgeneric edge-weight (tessellation a b)
  (:documentation "Return the weight (distance/cost of crossing) of the
edge between adjacent cells A and B in TESSELLATION. Symmetric: (EDGE-WEIGHT
tess a b) and (EDGE-WEIGHT tess b a) are always equal."))

;;; Geometry: needed by the mechanistic reaction-diffusion engine
;;; (mechanistic-step.lisp) for real finite-volume flux, on top of what
;;; topological diffusion (front-propagation, Dijkstra) already uses.
;;; The two-point flux across a face is D * EDGE-LENGTH / FACE-DISTANCE * dC:
;;; EDGE-LENGTH is the length of the shared boundary, FACE-DISTANCE the
;;; distance between the two points where the cells hold their values.
;;;
;;; FACE-DISTANCE is deliberately NOT EDGE-WEIGHT. EDGE-WEIGHT is the cost
;;; of one hop for front propagation and Dijkstra; the CA-grown grid sets it
;;; to 1.0 on every edge, which is right for hop counting and the published
;;; isotropy study but wrong as a flux distance (real seed spacing on a
;;; 200x200 grid with 400 cells is 3-22 pixels - PLAN.md, Phase 1). Keeping
;;; the two apart lets each caller read the quantity it actually means.
;;;
;;; Deliberately NOT defaulted to 1.0 for tessellations without geometry
;;; data (e.g. hand-built test fixtures for the topological-only tests) -
;;; a silent default here would be exactly the shape of bug the eql/equal
;;; mixup was: something that looks like it works but silently computes on
;;; the wrong data. See the HASH-GRAPH-TESSELLATION methods below.

(defgeneric cell-area (tessellation cell-id)
  (:documentation "Return the area of CELL-ID in TESSELLATION (a 2D proxy
for volume, used to convert flux into a concentration rate)."))

(defgeneric edge-length (tessellation a b)
  (:documentation "Return the length of the shared boundary between
adjacent cells A and B in TESSELLATION. Symmetric, like EDGE-WEIGHT."))

(defgeneric face-distance (tessellation a b)
  (:documentation "Return the distance between the points at which
adjacent cells A and B hold their values - the d of the two-point
finite-volume flux D * EDGE-LENGTH / d * (C_b - C_a). Seed-to-seed distance
for Voronoi grids, centroid distance for GIS polygons. Always positive and
symmetric. Not the same as EDGE-WEIGHT - see the comment above."))

(defgeneric cell-position (tessellation cell-id)
  (:documentation "Return the point (x . y), as double-floats, at which
CELL-ID holds its value: its seed on a Voronoi grid, its centroid for a GIS
polygon. FACE-DISTANCE is measured between these points, so this is where
the two-point flux is consistent. Initial conditions, source placement and
comparisons with analytic solutions evaluate fields here."))

;;; The one concrete implementation used everywhere in this project: an
;;; explicit adjacency table, id -> alist of (neighbor-id . weight). Any
;;; builder that can produce such a table (CA growth today; GIS polygons;
;;; a Delaunay triangulation, etc. later) can be wrapped in this class.
;;;
;;; CA-grown cells are roughly equal size, so voronoi-ca.lisp gives every
;;; edge a uniform weight of 1.0 - which makes weighted diffusion exactly
;;; equivalent to the plain hop-counting it replaces there. Real GIS
;;; polygons vary hugely in size, so gis-tessellation.lisp gives each edge
;;; a real distance. Same representation either way.

(defclass hash-graph-tessellation ()
  ((adjacency :initarg :adjacency :reader tessellation-adjacency
              :documentation "EQL hash-table: cell-id -> alist of (neighbor-id . weight).")
   (cell-areas :initarg :cell-areas :initform nil :reader tessellation-cell-areas
               :documentation "EQUAL hash-table: cell-id -> area, or NIL if
this tessellation was built without geometry data - see CELL-AREA.")
   (edge-lengths :initarg :edge-lengths :initform nil :reader tessellation-edge-lengths
                 :documentation "EQUAL hash-table: cell-id -> alist of
(neighbor-id . length), or NIL if built without geometry data - see EDGE-LENGTH.")
   (face-distances :initarg :face-distances :initform nil :reader tessellation-face-distances
                   :documentation "EQUAL hash-table: cell-id -> alist of
(neighbor-id . distance), or NIL if built without geometry data - see
FACE-DISTANCE and FACE-DISTANCES-FROM-POSITIONS.")
   (cell-positions :initarg :cell-positions :initform nil :reader tessellation-cell-positions
                   :documentation "EQUAL hash-table: cell-id -> (x . y)
double-float point, or NIL if built without geometry data - see
CELL-POSITION.")))

(defun make-hash-graph-tessellation (adjacency &key cell-areas edge-lengths face-distances
                                                    cell-positions)
  "Return a HASH-GRAPH-TESSELLATION over ADJACENCY. CELL-AREAS,
EDGE-LENGTHS, FACE-DISTANCES and CELL-POSITIONS are optional - omit them (as
every topological-only caller does) to get a tessellation that still
supports every topological operation (NEIGHBORS-OF, EDGE-WEIGHT, front
propagation, Dijkstra), just not the finite-volume geometry. Builders
derive FACE-DISTANCES from CELL-POSITIONS (FACE-DISTANCES-FROM-POSITIONS),
so the two agree by construction."
  (make-instance 'hash-graph-tessellation
                 :adjacency adjacency :cell-areas cell-areas :edge-lengths edge-lengths
                 :face-distances face-distances :cell-positions cell-positions))

(defmethod cell-ids ((tess hash-graph-tessellation))
  (loop for id being the hash-keys of (tessellation-adjacency tess) collect id))

(defmethod neighbors-of ((tess hash-graph-tessellation) cell-id)
  (mapcar #'car (gethash cell-id (tessellation-adjacency tess))))

(defmethod cell-count ((tess hash-graph-tessellation))
  (hash-table-count (tessellation-adjacency tess)))

(defmethod edge-weight ((tess hash-graph-tessellation) a b)
  (cdr (assoc b (gethash a (tessellation-adjacency tess)) :test 'equal)))

(defmethod cell-area ((tess hash-graph-tessellation) cell-id)
  (let ((areas (tessellation-cell-areas tess)))
    (unless areas
      (error "CELL-AREA: this tessellation was built without geometry data - ~
it can only be used for topological diffusion (front-propagation, Dijkstra), ~
not the mechanistic reaction-diffusion engine."))
    (multiple-value-bind (area foundp) (gethash cell-id areas)
      (unless foundp
        (error "CELL-AREA: no area recorded for cell ~s" cell-id))
      area)))

(defmethod edge-length ((tess hash-graph-tessellation) a b)
  (let ((lengths (tessellation-edge-lengths tess)))
    (unless lengths
      (error "EDGE-LENGTH: this tessellation was built without geometry data - ~
it can only be used for topological diffusion (front-propagation, Dijkstra), ~
not the mechanistic reaction-diffusion engine."))
    ;; No entry for this specific pair, even though the table itself
    ;; exists, is a legitimate zero, not a gap: BUILD-ADJACENCY's
    ;; 8-connected adjacency detection (see voronoi-ca.lisp) counts two
    ;; only-diagonally-touching cells as neighbors (needed for the
    ;; isotropy-averaging growth rule), but they share a single corner
    ;; point, not a boundary - zero shared length, correctly zero flux.
    (let ((entry (assoc b (gethash a lengths) :test 'equal)))
      (if entry (cdr entry) 0.0d0))))

(defmethod face-distance ((tess hash-graph-tessellation) a b)
  (let ((distances (tessellation-face-distances tess)))
    (unless distances
      (error "FACE-DISTANCE: this tessellation was built without geometry data - ~
it can only be used for topological diffusion (front-propagation, Dijkstra), ~
not the mechanistic reaction-diffusion engine."))
    ;; Unlike EDGE-LENGTH, a missing pair is never a legitimate value: every
    ;; pair of neighbours has two value points some positive distance apart.
    (let ((entry (assoc b (gethash a distances) :test 'equal)))
      (unless entry
        (error "FACE-DISTANCE: no distance recorded between cells ~s and ~s" a b))
      (cdr entry))))

(defmethod cell-position ((tess hash-graph-tessellation) cell-id)
  (let ((positions (tessellation-cell-positions tess)))
    (unless positions
      (error "CELL-POSITION: this tessellation was built without geometry data."))
    (multiple-value-bind (point foundp) (gethash cell-id positions)
      (unless foundp
        (error "CELL-POSITION: no position recorded for cell ~s" cell-id))
      point)))

(defun nearest-cell (tessellation point)
  "Return the id of the cell of TESSELLATION whose CELL-POSITION is closest
to POINT (x . y). On a nearest-seed Voronoi grid this is the cell that
contains POINT; on other tessellations it is the cell whose value point is
nearest. Used to place sources. Ties go to the first cell in CELL-IDS order."
  (let ((ids (cell-ids tessellation)))
    (reduce (lambda (best id)
              (if (< (point-distance point (cell-position tessellation id))
                     (point-distance point (cell-position tessellation best)))
                  id best))
            (rest ids) :initial-value (first ids))))

(defun positions-table (ids position-of)
  "Return a fresh EQUAL hash-table cell-id -> (x . y) for every id in IDS,
with each point from (FUNCALL POSITION-OF id) converted to double-floats."
  (let ((table (make-hash-table :test 'equal :size (length ids))))
    (dolist (id ids table)
      (let ((p (funcall position-of id)))
        (setf (gethash id table) (cons (float (car p) 1.0d0) (float (cdr p) 1.0d0)))))))

(defun point-distance (p0 p1)
  "Return the Euclidean distance between points P0 and P1, each (x . y),
as a double-float (integer pixel coordinates are converted first, so the
square root is never taken in single precision)."
  (let ((dx (float (- (car p1) (car p0)) 1.0d0))
        (dy (float (- (cdr p1) (cdr p0)) 1.0d0)))
    (sqrt (+ (* dx dx) (* dy dy)))))

(defun face-distances-from-positions (adjacency position-of)
  "Return a fresh EQUAL hash-table cell-id -> alist of (neighbor-id .
distance), with one symmetric entry for every adjacent pair in ADJACENCY.
The distance is POINT-DISTANCE between (FUNCALL POSITION-OF a) and
\(FUNCALL POSITION-OF b) - the points where the two cells hold their
values. Signals an error if two adjacent cells share a position, because a
zero distance would make the two-point flux infinite."
  (let ((distances (make-hash-table :test 'equal)))
    (maphash (lambda (a neighbors)
               (ensure-cell! distances a)
               (dolist (b (mapcar #'car neighbors))
                 (let ((d (point-distance (funcall position-of a) (funcall position-of b))))
                   (unless (plusp d)
                     (error "FACE-DISTANCES-FROM-POSITIONS: adjacent cells ~s and ~s ~
share the position ~s" a b (funcall position-of a)))
                   (add-edge! distances a b d))))
             adjacency)
    distances))

;;; Shared by every tessellation builder, so an edge is never accidentally
;;; recorded in only one direction.

(defun ensure-cell! (adjacency id)
  "Make sure ID has a (possibly empty) entry in ADJACENCY."
  (unless (nth-value 1 (gethash id adjacency))
    (setf (gethash id adjacency) '())))

(defun add-edge! (adjacency a b weight)
  "Record a two-way edge between A and B with WEIGHT in ADJACENCY (a
hash-table id -> alist of (neighbor-id . weight)). Idempotent - calling it
again for a pair that's already connected does nothing.

Cell ids are integers for CA-grown tessellations but strings for GIS ones
(taken straight from GeoJSON feature properties), so every comparison here
is :TEST 'EQUAL, never the default EQL - EQL happens to work for integers
but silently fails to match two distinct string objects with the same
characters, which is exactly the shape of bug a caller doing
\(compute-weighted-arrival-times tess (list \"some-id\")\) with a freshly
constructed string would hit. ADJACENCY itself must also be built with
:TEST 'EQUAL for the same reason - see voronoi-ca.lisp/gis-tessellation.lisp."
  (ensure-cell! adjacency a)
  (ensure-cell! adjacency b)
  (unless (assoc b (gethash a adjacency) :test 'equal)
    (push (cons b weight) (gethash a adjacency)))
  (unless (assoc a (gethash b adjacency) :test 'equal)
    (push (cons a weight) (gethash b adjacency))))

(defun set-edge-length! (edge-lengths a b length)
  "Record a two-way edge LENGTH between A and B in EDGE-LENGTHS (same shape
ADD-EDGE! uses: id -> alist of (neighbor-id . length)). Unlike ADD-EDGE!,
repeat calls for an already-recorded pair ADD to the length rather than
being a no-op - both builders call this once per unit of shared boundary
they find (once per 4-connected pixel pair for the CA grid, once per shared
polygon edge for GIS regions), so a pair touching along more than one
stretch of boundary gets its true total contact length, not just the first
piece found."
  (ensure-cell! edge-lengths a)
  (ensure-cell! edge-lengths b)
  (let ((entry-ab (assoc b (gethash a edge-lengths) :test 'equal)))
    (if entry-ab
        (incf (cdr entry-ab) length)
        (push (cons b length) (gethash a edge-lengths))))
  (let ((entry-ba (assoc a (gethash b edge-lengths) :test 'equal)))
    (if entry-ba
        (incf (cdr entry-ba) length)
        (push (cons a length) (gethash b edge-lengths)))))

;;; A dense 0..N-1 index over a tessellation's cells, in whatever order
;;; CELL-IDS returns them. The mechanistic engine (mechanistic-world.lisp,
;;; mechanistic-step.lisp) stores per-tick state as plain
;;; (ARRAY DOUBLE-FLOAT (n-cells n-species)) arrays keyed by this index,
;;; rather than a hash-table per tick - needed for PMAP-CELL-RANGE! (see
;;; parallel.lisp) to split the cell range into contiguous chunks across
;;; worker threads.

(defstruct (cell-index (:constructor %make-cell-index))
  ids        ; (simple-vector n) - index -> cell-id
  id->index) ; EQUAL hash-table: cell-id -> index

(defun make-cell-index (tessellation)
  (let* ((ids (coerce (cell-ids tessellation) 'simple-vector))
         (id->index (make-hash-table :test 'equal :size (length ids))))
    (loop for i from 0 below (length ids)
          do (setf (gethash (aref ids i) id->index) i))
    (%make-cell-index :ids ids :id->index id->index)))

(defun cell-index-count (cell-index)
  (length (cell-index-ids cell-index)))

(defun index->cell-id (cell-index i)
  (aref (cell-index-ids cell-index) i))

(defun cell-id->index (cell-index id)
  (or (gethash id (cell-index-id->index cell-index))
      (error "CELL-ID->INDEX: ~s is not a cell in this CELL-INDEX" id)))
