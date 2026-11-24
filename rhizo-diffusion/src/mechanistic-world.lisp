(in-package :rhizo-diffusion)

;;; MECHANISTIC-WORLD: the state container for the reaction-diffusion engine
;;; - a third "shape of computation" alongside WORLD (tick-based front
;;; propagation, world.lisp) and the one-shot Dijkstra of
;;; weighted-diffusion.lisp, following the same pattern already established
;;; in this codebase of giving a genuinely different kind of computation its
;;; own type rather than overloading an existing one.
;;;
;;; Per-tick state lives in one dense (ARRAY DOUBLE-FLOAT (n-cells
;;; n-species)) array, indexed via the tessellation's CELL-INDEX
;;; (tessellation.lisp) rather than per-cell hash-tables - needed so
;;; PMAP-CELL-RANGE! (parallel.lisp) can split the cell range into
;;; contiguous chunks across worker threads in mechanistic-step.lisp.

(defstruct (mechanistic-world (:constructor %make-mechanistic-world))
  tessellation    ; must have real CELL-AREA/EDGE-LENGTH geometry
  cell-index      ; dense id <-> array-index mapping (tessellation.lisp)
  stencil         ; FLUX-STENCIL in CELL-INDEX order, built once (flux-stencil.lisp)
  species-list
  species-index   ; name -> array column, shared with every REACTION's rate-fn
  reactions
  n-species
  state           ; (array double-float (n-cells n-species)) - current tick
  diffusivities   ; (array double-float (n-cells n-species)) - D of each species in
                  ; each cell; faces use the harmonic mean (flux-stencil.lisp)
  capacities      ; (array double-float (n-cells n-species)) - storage capacity kappa:
                  ; the state is an amount per cell volume m, and diffusion moves
                  ; it down gradients of m / kappa (e.g. kappa = water content,
                  ; m / kappa = pore-water concentration). 1 by default.
  dt
  (diffusion-scheme :explicit) ; :EXPLICIT, :IMPLICIT or :IMEX - see STEP-MECHANISTIC!
  (tick 0 :type fixnum)
  (face-fluxes nil)       ; NIL (no flow) or, per stencil face, the water flux (volume per
                          ; time) out of the cell through it - see SET-FLOW-FIELD!
  (boundary-outflows nil) ; NIL or, per cell, water leaving through outer-boundary faces
  exports)                ; per species: cumulative amount carried out with that water

(defun species-field-array (cell-index species field default valid-p what)
  "Return a fresh (array double-float (n-cells n-species)) in CELL-INDEX
order holding, for each species, (FUNCALL DEFAULT species) in every cell -
except for species named in FIELD, an alist (species-name . function of a
cell-id returning that cell's value). Signals an error, naming WHAT, on a
value failing VALID-P or an unknown species name."
  (let* ((n-cells (cell-index-count cell-index))
         (array (make-array (list n-cells (length species)) :element-type 'double-float)))
    (dolist (entry field)
      (unless (find (car entry) species :key #'species-name)
        (error "~a: unknown species ~s" what (car entry))))
    (loop for sp in species
          for s from 0
          for fn = (cdr (assoc (species-name sp) field))
          do (dotimes (i n-cells)
               (let ((value (float (if fn
                                       (funcall fn (index->cell-id cell-index i))
                                       (funcall default sp))
                                   1.0d0)))
                 (unless (funcall valid-p value)
                   (error "~a: invalid value ~a for ~s in cell ~s"
                          what value (species-name sp) (index->cell-id cell-index i)))
                 (setf (aref array i s) value))))
    array))

(defun diffusivity-array (cell-index species diffusivity-field)
  "Return a fresh (array double-float (n-cells n-species)) of diffusivities
in CELL-INDEX order: each species' SPECIES-DIFFUSIVITY in every cell, except
for species named in DIFFUSIVITY-FIELD, an alist (species-name . function of
a cell-id returning that cell's diffusivity). Signals an error on a negative
value or an unknown species name."
  (species-field-array cell-index species diffusivity-field #'species-diffusivity
                       (lambda (d) (not (minusp d))) "DIFFUSIVITY-ARRAY"))

(defun capacity-array (cell-index species capacity-field)
  "Return a fresh (array double-float (n-cells n-species)) of storage
capacities in CELL-INDEX order: 1 for every species in every cell, except
for species named in CAPACITY-FIELD, an alist (species-name . function of a
cell-id). Signals an error on a value that is not positive or an unknown
species name."
  (species-field-array cell-index species capacity-field (constantly 1d0) #'plusp
                       "CAPACITY-ARRAY"))

(defun make-mechanistic-world (tessellation &key species reactions dt initial-state
                                                  (diffusion-scheme :explicit)
                                                  diffusivity-field capacity-field)
  "Build a fresh MECHANISTIC-WORLD over TESSELLATION. SPECIES/REACTIONS
define the reaction network (see reaction-network.lisp; examples/ holds a
complete network built this way). DT is
the tick size - see STABLE-TIMESTEP for how to choose a safe one for a
given TESSELLATION/SPECIES pair. INITIAL-STATE is either an alist
\(species-name . value), applied to every cell, or a function of one
cell-id returning such an alist - so a source-region seeding pattern (like
SEED-MOLECULES-AT! in diffusion.lisp) is possible. DIFFUSION-SCHEME is
:EXPLICIT (DIFFUSION-SUBSTEP!, DT must respect STABLE-TIMESTEP) or :IMPLICIT
\(IMPLICIT-DIFFUSION-SUBSTEP!, backward Euler, stable at any DT), both
Strang-split with RK4 reactions, or :IMEX (IMEX Euler, not split: exact
steady states, explicit Euler reactions - see STEP-MECHANISTIC!).
DIFFUSIVITY-FIELD gives species a per-cell diffusivity (see
DIFFUSIVITY-ARRAY); other species use their SPECIES-DIFFUSIVITY everywhere.
SET-SPECIES-DIFFUSIVITY! changes it during a run. CAPACITY-FIELD likewise
gives species a per-cell storage capacity (CAPACITY-ARRAY, default 1;
SET-SPECIES-CAPACITY!)."
  (unless (member diffusion-scheme '(:explicit :implicit :imex))
    (error "MAKE-MECHANISTIC-WORLD: unknown diffusion scheme ~s" diffusion-scheme))
  (let* ((cell-index (make-cell-index tessellation))
         (n-cells (cell-index-count cell-index))
         (species-index (make-species-index species))
         (n-species (length species))
         (state (make-array (list n-cells n-species)
                             :element-type 'double-float :initial-element 0.0d0)))
    (dotimes (i n-cells)
      (let* ((id (index->cell-id cell-index i))
             (values (if (functionp initial-state) (funcall initial-state id) initial-state))
             (row (make-cell-state species-index values)))
        (dotimes (s n-species)
          (setf (aref state i s) (aref row s)))))
    (%make-mechanistic-world
     :tessellation tessellation :cell-index cell-index
     :stencil (make-flux-stencil tessellation cell-index)
     :species-list species :species-index species-index :reactions reactions
     :n-species n-species :state state :dt (float dt 1.0d0) :tick 0
     :diffusivities (diffusivity-array cell-index species diffusivity-field)
     :capacities (capacity-array cell-index species capacity-field)
     :diffusion-scheme diffusion-scheme
     :exports (make-array n-species :element-type 'double-float :initial-element 0d0))))

(defun set-species-column! (world array species-name field valid-p what)
  "Set column SPECIES-NAME of ARRAY, a per-cell per-species array of WORLD,
to (FUNCALL FIELD cell-id) in every cell. Signals an error, naming WHAT, on
an unknown species or a value failing VALID-P. Returns WORLD."
  (let ((index (mechanistic-world-cell-index world))
        (s (or (gethash species-name (mechanistic-world-species-index world))
               (error "~a: unknown species ~s" what species-name))))
    (dotimes (i (cell-index-count index) world)
      (let ((value (float (funcall field (index->cell-id index i)) 1.0d0)))
        (unless (funcall valid-p value)
          (error "~a: invalid value ~a in cell ~s" what value (index->cell-id index i)))
        (setf (aref array i s) value)))))

(defun set-species-diffusivity! (world species-name field)
  "Set the diffusivity of SPECIES-NAME in every cell of WORLD to (FUNCALL
FIELD cell-id), e.g. when water content changes (PLAN.md, Phase 3). Changes
only WORLD's diffusivity array. Signals an error on a negative value.
Returns WORLD."
  (set-species-column! world (mechanistic-world-diffusivities world) species-name field
                       (lambda (d) (not (minusp d))) "SET-SPECIES-DIFFUSIVITY!"))

(defun set-species-capacity! (world species-name field)
  "Set the storage capacity of SPECIES-NAME in every cell of WORLD to
\(FUNCALL FIELD cell-id), e.g. the water content (PLAN.md, Phase 3). The
amounts held are unchanged, so a smaller capacity raises m / kappa: solute
stays behind as water leaves. Signals an error on a value that is not
positive. Returns WORLD."
  (set-species-column! world (mechanistic-world-capacities world) species-name field
                       #'plusp "SET-SPECIES-CAPACITY!"))

(defun world-stable-timestep (world &key (safety-factor 0.9d0))
  "Return the explicit-diffusion stability bound for WORLD's own per-cell
diffusivities and capacities: the smallest STENCIL-STABLE-TIMESTEP over its
species. Infinity if nothing diffuses."
  (let ((stencil (mechanistic-world-stencil world))
        (diffusivities (mechanistic-world-diffusivities world))
        (capacities (mechanistic-world-capacities world)))
    (reduce #'min (loop for s below (mechanistic-world-n-species world)
                        collect (stencil-stable-timestep stencil (state-column diffusivities s)
                                                         :safety-factor safety-factor
                                                         :capacity (state-column capacities s)))
            :initial-value sb-ext:double-float-positive-infinity)))

(defun field-initial-state (tessellation fields)
  "Return a function of one cell-id, suitable as MAKE-MECHANISTIC-WORLD's
:INITIAL-STATE, giving the alist (species-name . value) obtained by
evaluating each field at the cell's CELL-POSITION. FIELDS is an alist
\(species-name . function of a point (x . y)). Sampling at the value point
rather than averaging over the cell is what a cell-centred scheme assumes;
it is second-order accurate for smooth fields."
  (lambda (id)
    (let ((point (cell-position tessellation id)))
      (mapcar (lambda (entry) (cons (car entry) (funcall (cdr entry) point))) fields))))

;;; Water flow through the world (PLAN.md, Phase 5). A flow field is a
;;; water flux per stencil face, positive out of the cell, antisymmetric
;;; (the flux out of i into j is minus that out of j into i), plus each
;;; cell's water outflow through the outer boundary. Water entering
;;; through the boundary (rain) carries no solute. RICHARDS-STEP produces
;;; exactly such fields; the transport step advects every dissolved species
;;; with them (ADVECTION-DIFFUSION-COLUMN).

(defun set-flow-field! (world face-fluxes boundary-outflows)
  "Set WORLD's water flow: FACE-FLUXES, a double-float vector over its
stencil faces (volume per time out of the cell through each face; see
above), and BOUNDARY-OUTFLOWS, one non-negative value per cell, or NIL for
none. FACE-FLUXES NIL removes the flow. Changes only WORLD. Returns
WORLD."
  (let ((n-faces (length (flux-stencil-neighbors (mechanistic-world-stencil world))))
        (n-cells (cell-index-count (mechanistic-world-cell-index world))))
    (when (and face-fluxes (/= (length face-fluxes) n-faces))
      (error "SET-FLOW-FIELD!: ~d face fluxes for ~d faces" (length face-fluxes) n-faces))
    (when boundary-outflows
      (unless (= (length boundary-outflows) n-cells)
        (error "SET-FLOW-FIELD!: ~d boundary outflows for ~d cells" (length boundary-outflows) n-cells))
      (when (some #'minusp boundary-outflows)
        (error "SET-FLOW-FIELD!: a boundary outflow is negative")))
    (setf (mechanistic-world-face-fluxes world)
          (and face-fluxes (coerce face-fluxes '(simple-array double-float (*))))
          (mechanistic-world-boundary-outflows world)
          (and face-fluxes boundary-outflows
               (coerce boundary-outflows '(simple-array double-float (*)))))
    world))

(defun species-exported (world species)
  "Return the cumulative amount of SPECIES (one name or a list summed)
carried out of WORLD with water leaving through its outer boundary."
  (loop for name in (alexandria:ensure-list species)
        sum (aref (mechanistic-world-exports world)
                  (or (gethash name (mechanistic-world-species-index world))
                      (error "SPECIES-EXPORTED: unknown species ~s" name)))))
