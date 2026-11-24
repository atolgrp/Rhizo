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
  dt
  (diffusion-scheme :explicit) ; :EXPLICIT or :IMPLICIT - see STEP-MECHANISTIC!
  (tick 0 :type fixnum))

(defun diffusivity-array (cell-index species diffusivity-field)
  "Return a fresh (array double-float (n-cells n-species)) of diffusivities
in CELL-INDEX order: each species' SPECIES-DIFFUSIVITY in every cell, except
for species named in DIFFUSIVITY-FIELD, an alist (species-name . function of
a cell-id returning that cell's diffusivity). Signals an error on a negative
value or an unknown species name."
  (let* ((n-cells (cell-index-count cell-index))
         (array (make-array (list n-cells (length species)) :element-type 'double-float)))
    (dolist (entry diffusivity-field)
      (unless (find (car entry) species :key #'species-name)
        (error "DIFFUSIVITY-ARRAY: unknown species ~s" (car entry))))
    (loop for sp in species
          for s from 0
          for field = (cdr (assoc (species-name sp) diffusivity-field))
          do (dotimes (i n-cells)
               (let ((d (float (if field
                                   (funcall field (index->cell-id cell-index i))
                                   (species-diffusivity sp))
                               1.0d0)))
                 (when (minusp d)
                   (error "DIFFUSIVITY-ARRAY: negative diffusivity ~a for ~s in cell ~s"
                          d (species-name sp) (index->cell-id cell-index i)))
                 (setf (aref array i s) d))))
    array))

(defun make-mechanistic-world (tessellation &key species reactions dt initial-state
                                                  (diffusion-scheme :explicit)
                                                  diffusivity-field)
  "Build a fresh MECHANISTIC-WORLD over TESSELLATION. SPECIES/REACTIONS
define the reaction network (see reaction-network.lisp; examples/ holds a
complete network built this way). DT is
the tick size - see STABLE-TIMESTEP for how to choose a safe one for a
given TESSELLATION/SPECIES pair. INITIAL-STATE is either an alist
\(species-name . value), applied to every cell, or a function of one
cell-id returning such an alist - so a source-region seeding pattern (like
SEED-MOLECULES-AT! in diffusion.lisp) is possible. DIFFUSION-SCHEME is
:EXPLICIT (DIFFUSION-SUBSTEP!, DT must respect STABLE-TIMESTEP) or :IMPLICIT
\(IMPLICIT-DIFFUSION-SUBSTEP!, backward Euler, stable at any DT).
DIFFUSIVITY-FIELD gives species a per-cell diffusivity (see
DIFFUSIVITY-ARRAY); other species use their SPECIES-DIFFUSIVITY everywhere.
SET-SPECIES-DIFFUSIVITY! changes it during a run."
  (unless (member diffusion-scheme '(:explicit :implicit))
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
     :diffusion-scheme diffusion-scheme)))

(defun set-species-diffusivity! (world species-name field)
  "Set the diffusivity of SPECIES-NAME in every cell of WORLD to (FUNCALL
FIELD cell-id), e.g. when water content changes (PLAN.md, Phase 3). Changes
only WORLD's diffusivity array. Signals an error on a negative value.
Returns WORLD."
  (let* ((index (mechanistic-world-cell-index world))
         (s (or (gethash species-name (mechanistic-world-species-index world))
                (error "SET-SPECIES-DIFFUSIVITY!: unknown species ~s" species-name)))
         (array (mechanistic-world-diffusivities world)))
    (dotimes (i (cell-index-count index) world)
      (let ((d (float (funcall field (index->cell-id index i)) 1.0d0)))
        (when (minusp d)
          (error "SET-SPECIES-DIFFUSIVITY!: negative diffusivity ~a in cell ~s"
                 d (index->cell-id index i)))
        (setf (aref array i s) d)))))

(defun world-stable-timestep (world &key (safety-factor 0.9d0))
  "Return the explicit-diffusion stability bound for WORLD's own per-cell
diffusivities: the smallest STENCIL-STABLE-TIMESTEP over its species.
Infinity if nothing diffuses."
  (let ((stencil (mechanistic-world-stencil world))
        (diffusivities (mechanistic-world-diffusivities world)))
    (reduce #'min (loop for s below (mechanistic-world-n-species world)
                        collect (stencil-stable-timestep stencil (state-column diffusivities s)
                                                         :safety-factor safety-factor))
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
