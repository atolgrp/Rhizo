(in-package :rhizo-diffusion)

;;; CA-grown Voronoi tessellation.
;;;
;;; Port of the "isotropy-11" algorithm from RhizoCA/diff30.org: scatter
;;; random seed points on a grid, then flood-fill outward from each one
;;; generation-by-generation, ALTERNATING between 4-connected (von Neumann)
;;; and 8-connected (Moore) neighborhoods on even/odd generations. That
;;; alternation is what removes the diagonal-bias artifact a fixed
;;; neighborhood would have - it was validated empirically in the original
;;; notebook (front-propagation arrival times at N/E/S/W vs NW/NE/SW/SE
;;; converge) and again here (see tests and the validation script).
;;;
;;; Differences from the original, both deliberate:
;;;  - No global specials. Growth state lives in local arrays/tables owned
;;;    by this function call, so two calls (even concurrent ones) never
;;;    interfere with each other.
;;;  - The final adjacency graph is built in ONE pass over the grid
;;;    (O(width*height)), not one pass PER SEED POINT
;;;    (O(points*width*height), the confirmed bottleneck in the old code).

(defparameter *neighbor-offsets-4*
  '((-1 . 0) (1 . 0) (0 . -1) (0 . 1)))

(defparameter *neighbor-offsets-8*
  (append *neighbor-offsets-4* '((-1 . -1) (-1 . 1) (1 . -1) (1 . 1))))

(defstruct (voronoi-grid (:constructor %make-voronoi-grid))
  tessellation ; a HASH-GRAPH-TESSELLATION - what diffusion/metrics actually use
  width
  height
  ids          ; (array fixnum (width height)) - kept around so metrics can
               ; map a grid position (e.g. "the cell 1/4 of the way east of
               ; center") to a cell id. Diffusion code never looks at this.
  seed         ; integer seed used to grow this grid, recorded for reproducibility
  seed-points) ; simple-vector of (x . y) integer pixels; element I is the seed
               ; of cell I, where that cell's value sits (CELL-POSITION)

(declaim (inline in-bounds-p))
(defun in-bounds-p (x y width height)
  (and (<= 0 x) (< x width) (<= 0 y) (< y height)))

(defun grid-neighbor-coords (x y width height offsets)
  (loop for (dx . dy) in offsets
        for nx = (+ x dx)
        for ny = (+ y dy)
        when (in-bounds-p nx ny width height)
          collect (cons nx ny)))

(defun random-distinct-points (width height n random-state)
  "Return a vector of N distinct (x . y) conses drawn uniformly from a
WIDTH x HEIGHT grid, using RANDOM-STATE (never the global *random-state*,
so this is safe to call concurrently for different worlds)."
  (let ((seen (make-hash-table :test 'equal))
        (points (make-array n :fill-pointer 0))
        (max-attempts (* 50 (max n 1))))
    (loop repeat max-attempts
          until (= (fill-pointer points) n)
          do (let ((p (cons (random width random-state) (random height random-state))))
               (unless (gethash p seen)
                 (setf (gethash p seen) t)
                 (vector-push p points))))
    (unless (= (fill-pointer points) n)
      (warn "random-distinct-points: only found ~a of ~a requested distinct points (grid too small/dense?)"
            (fill-pointer points) n))
    points))

(defun grow-voronoi-ids (width height seed-points random-state)
  "Flood-fill grid of size WIDTH x HEIGHT outward from SEED-POINTS (a
vector of (x . y) conses; point at index I gets id I). Returns the
completed (array fixnum (width height)) with every cell assigned.
Alternates 4-/8-connected neighborhoods by generation parity, and
processes each generation's assigned cells in random order so that
boundary ties between two equidistant regions are broken randomly rather
than systematically."
  (let ((current (make-array (list width height) :element-type 'fixnum :initial-element -1))
        (safety-cap (* 4 (+ width height))))
    (loop for i from 0
          for (x . y) across seed-points
          do (setf (aref current x y) i))
    (loop with generation = 0
          with unassigned = (- (* width height) (length seed-points))
          while (and (plusp unassigned) (< generation safety-cap))
          do (let ((next (alexandria:copy-array current))
                    (offsets (if (evenp generation) *neighbor-offsets-4* *neighbor-offsets-8*))
                    (assigned-cells '()))
                (dotimes (x width)
                  (dotimes (y height)
                    (when (/= -1 (aref current x y))
                      (push (cons x y) assigned-cells))))
                ;; alexandria:shuffle only reads the global *random-state*, so bind it
                ;; locally to this call's own state - keeps growth reproducible/parallel-safe.
                (dolist (cell (let ((*random-state* random-state))
                                (alexandria:shuffle assigned-cells)))
                  (let* ((x (car cell)) (y (cdr cell)) (id (aref current x y)))
                    (dolist (nb (grid-neighbor-coords x y width height offsets))
                      (when (= -1 (aref next (car nb) (cdr nb)))
                        (setf (aref next (car nb) (cdr nb)) id)
                        (decf unassigned)))))
                (setf current next)
                (incf generation))
          finally (when (plusp unassigned)
                    (warn "grow-voronoi-ids: ~a cells still unassigned after ~a generations"
                          unassigned generation)))
    current))

(defun build-adjacency (ids width height)
  "Single pass over the grid: for every cell, look at its 8-connected
neighbors, and record a two-way edge whenever a neighbor has a different
id. O(width*height), independent of how many seed points were used -
this is the fix for the old O(points*width*height) bottleneck.

Every edge gets a uniform weight of 1.0 - CA-grown cells are roughly equal
size, so this makes weighted diffusion (see weighted-diffusion.lisp)
exactly equivalent to the plain hop-counting front-propagation.lisp
already does, and leaves the previously-validated isotropy result
unchanged."
  (let ((adjacency (make-hash-table :test 'equal)))
    (dotimes (x width)
      (dotimes (y height)
        (let ((id (aref ids x y)))
          (ensure-cell! adjacency id)
          (dolist (nb (grid-neighbor-coords x y width height *neighbor-offsets-8*))
            (let ((nb-id (aref ids (car nb) (cdr nb))))
              (when (/= nb-id id)
                (add-edge! adjacency id nb-id 1.0d0)))))))
    adjacency))

(defparameter *neighbor-offsets-4-forward* '((1 . 0) (0 . 1))
  "Only the 'right' and 'up' halves of *NEIGHBOR-OFFSETS-4* - visiting every
adjacent pixel PAIR exactly once (not once from each side) is what
BUILD-GEOMETRY needs for edge-length counting, since SET-EDGE-LENGTH!
accumulates on every call rather than being idempotent like ADD-EDGE! -
calling it once per side would silently double-count every shared
boundary.")

(defun build-geometry (ids width height)
  "Single pass alongside BUILD-ADJACENCY: computes each cell's area (pixel
count) and the shared-boundary length between adjacent cells, for the
mechanistic reaction-diffusion engine's finite-volume flux (see
tessellation.lisp CELL-AREA/EDGE-LENGTH, mechanistic-step.lisp).

Edge length is counted only over 4-connected neighbor pairs, not the
8-connected set BUILD-ADJACENCY uses for detecting adjacency - two
diagonally-touching pixels share a single corner point, not an edge, so
counting them would overstate contact length. This is a standard discrete
perimeter-counting approximation from raster landscape analysis. Returns
(VALUES CELL-AREAS EDGE-LENGTHS)."
  (let ((cell-areas (make-hash-table :test 'equal))
        (edge-lengths (make-hash-table :test 'equal)))
    (dotimes (x width)
      (dotimes (y height)
        (let ((id (aref ids x y)))
          (incf (gethash id cell-areas 0.0d0) 1.0d0)
          (dolist (nb (grid-neighbor-coords x y width height *neighbor-offsets-4-forward*))
            (let ((nb-id (aref ids (car nb) (cdr nb))))
              (when (/= nb-id id)
                (set-edge-length! edge-lengths id nb-id 1.0d0)))))))
    (values cell-areas edge-lengths)))

(defun make-voronoi-grid (width height n-points &key seed)
  "Build a CA-grown Voronoi tessellation on a WIDTH x HEIGHT grid with
N-POINTS seed regions. SEED is an integer used to seed this call's own
random-state (generated and recorded if not supplied), so a run is fully
reproducible from the returned VORONOI-GRID-SEED, and independent calls
never share random state - safe to run concurrently.

Growth, adjacency and the uniform EDGE-WEIGHT of 1.0 are exactly as in
the published isotropy study. Added alongside them, for the mechanistic
engine only: the seed points are kept, each cell's CELL-POSITION is its
seed, and FACE-DISTANCE is seed-to-seed distance (PLAN.md, Phase 1, steps
1-2). Edge lengths are still staircase pixel counts here - see PLAN.md,
Phase 1, step 4."
  (let* ((seed (or seed (random most-positive-fixnum (make-random-state t))))
         (random-state (sb-ext:seed-random-state seed))
         (seed-points (coerce (random-distinct-points width height n-points random-state)
                              'simple-vector))
         (ids (grow-voronoi-ids width height seed-points random-state))
         (adjacency (build-adjacency ids width height))
         (positions (positions-table (alexandria:iota (length seed-points))
                                     (lambda (id) (aref seed-points id)))))
    (multiple-value-bind (cell-areas edge-lengths) (build-geometry ids width height)
      (%make-voronoi-grid
       :tessellation (make-hash-graph-tessellation
                      adjacency :cell-areas cell-areas :edge-lengths edge-lengths
                                :cell-positions positions
                                :face-distances (face-distances-from-positions
                                                 adjacency (lambda (id) (gethash id positions))))
       :width width :height height :ids ids :seed seed :seed-points seed-points))))

(defun voronoi-grid-id-at (grid x y)
  (aref (voronoi-grid-ids grid) x y))
