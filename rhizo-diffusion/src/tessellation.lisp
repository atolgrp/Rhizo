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
;;; EDGE-WEIGHT stays centroid distance; EDGE-LENGTH is the length of the
;;; shared boundary itself - flux needs both: D * edge-length/distance * dC.
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
(neighbor-id . length), or NIL if built without geometry data - see EDGE-LENGTH.")))

(defun make-hash-graph-tessellation (adjacency &key cell-areas edge-lengths)
  "CELL-AREAS and EDGE-LENGTHS are optional - omit them (as every existing
caller before the mechanistic engine did) to get a tessellation that still
supports every topological operation (NEIGHBORS-OF, EDGE-WEIGHT, front
propagation, Dijkstra), just not CELL-AREA/EDGE-LENGTH."
  (make-instance 'hash-graph-tessellation
                 :adjacency adjacency :cell-areas cell-areas :edge-lengths edge-lengths))

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
