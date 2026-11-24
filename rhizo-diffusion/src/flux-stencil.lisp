(in-package :rhizo-diffusion)

;;; The finite-volume diffusion operator, precomputed once per world.
;;;
;;; DIFFUSION-SUBSTEP! used to walk the tessellation protocol on every tick:
;;; a fresh NEIGHBORS-OF list, two alist lookups for the face conductance,
;;; and a hash lookup per neighbour - none of which changes during a run.
;;; A FLUX-STENCIL holds the same numbers in compressed sparse row form: for
;;; cell I, the faces OFFSETS[I] .. OFFSETS[I+1]-1, each with the dense index
;;; of the neighbour across it and its conductance. It is built through the
;;; protocol, so the physics still never sees a specific grid type.
;;;
;;; Conductances are stored NOT divided by cell area, so the stencil is the
;;; symmetric matrix K of the two-point flux (K_ij = K_ji). Areas are kept
;;; separately. The backward-Euler system (A + dt D K) c' = A c (A the
;;; diagonal of areas, PLAN.md Phase 1 step 6) is then symmetric positive
;;; definite, which conjugate gradient needs.

(defun face-conductance (tessellation a b)
  "Return the geometric conductance EDGE-LENGTH / FACE-DISTANCE of the face
between adjacent cells A and B - the two-point flux across it, per unit
diffusivity and unit concentration difference. Zero for cells that touch
only at a corner (EDGE-LENGTH 0)."
  (/ (edge-length tessellation a b) (face-distance tessellation a b)))

(defstruct (flux-stencil (:constructor %make-flux-stencil))
  (offsets nil :type (simple-array fixnum (*)))           ; n + 1 row starts
  (neighbors nil :type (simple-array fixnum (*)))         ; dense index across each face
  (conductances nil :type (simple-array double-float (*))) ; FACE-CONDUCTANCE of each face
  (areas nil :type (simple-array double-float (*))))      ; CELL-AREA of each cell

(defun make-flux-stencil (tessellation cell-index)
  "Return a fresh FLUX-STENCIL for TESSELLATION, with rows in CELL-INDEX
order and each row's faces in NEIGHBORS-OF order. Faces of zero conductance
(corner-only contacts) are left out, since they carry no flux."
  (let* ((n (cell-index-count cell-index))
         (rows (loop for i below n
                     for id = (index->cell-id cell-index i)
                     collect (loop for nb in (neighbors-of tessellation id)
                                   for g = (face-conductance tessellation id nb)
                                   when (plusp g)
                                     collect (cons (cell-id->index cell-index nb)
                                                   (float g 1.0d0)))))
         (n-faces (reduce #'+ rows :key #'length))
         (offsets (make-array (1+ n) :element-type 'fixnum))
         (neighbors (make-array n-faces :element-type 'fixnum))
         (conductances (make-array n-faces :element-type 'double-float))
         (areas (make-array n :element-type 'double-float)))
    (loop with k = 0
          for i from 0
          for row in rows
          do (setf (aref offsets i) k
                   (aref areas i) (float (cell-area tessellation (index->cell-id cell-index i))
                                         1.0d0))
             (loop for (j . g) in row
                   do (setf (aref neighbors k) j
                            (aref conductances k) g)
                      (incf k))
          finally (setf (aref offsets n) k))
    (%make-flux-stencil :offsets offsets :neighbors neighbors
                        :conductances conductances :areas areas)))

(declaim (inline harmonic-mean)) ; called once per face per step in the loops below
(defun harmonic-mean (a b)
  "Return the harmonic mean 2ab / (a + b) of the non-negative numbers A and
B, as a double-float; 0 when either is 0. This is the face diffusivity
between two cells of diffusivity A and B on a Voronoi grid: the face bisects
the seed-to-seed segment, so the flux crosses two equal half-lengths in
series, and series resistances add. A cell of zero diffusivity is a barrier."
  (if (or (zerop a) (zerop b))
      0.0d0
      (/ (* 2.0d0 a b) (+ a b))))

(defun stencil-total-conductance (stencil diffusivity i)
  "Return the sum over cell I's faces in STENCIL of face diffusivity times
conductance, the face diffusivity being the HARMONIC-MEAN of DIFFUSIVITY (a
double-float vector, one value per cell) on either side. This is the
diagonal entry of the diffusion matrix for row I."
  (let ((offsets (flux-stencil-offsets stencil))
        (neighbors (flux-stencil-neighbors stencil))
        (conductances (flux-stencil-conductances stencil)))
    (declare (type (simple-array double-float (*)) diffusivity))
    (loop for k from (aref offsets i) below (aref offsets (1+ i))
          sum (* (aref conductances k)
                 (harmonic-mean (aref diffusivity i) (aref diffusivity (aref neighbors k))))
            of-type double-float)))

(defun stencil-net-flux (stencil x diffusivity i)
  "Return the net diffusive inflow into cell I of the field X (a
\(simple-array double-float (n-cells)), one species): the sum over faces of
D_face * conductance * (x_j - x_i), with D_face the HARMONIC-MEAN of
DIFFUSIVITY (one value per cell) on either side; i.e. -(K x)_i. Divide by
the cell's area to get dx_i/dt."
  (let ((offsets (flux-stencil-offsets stencil))
        (neighbors (flux-stencil-neighbors stencil))
        (conductances (flux-stencil-conductances stencil))
        (xi (aref x i))
        (di (aref diffusivity i)))
    (declare (type (simple-array double-float (*)) x diffusivity)
             (type double-float xi di))
    (loop for k from (aref offsets i) below (aref offsets (1+ i))
          for j = (aref neighbors k)
          sum (* (aref conductances k) (harmonic-mean di (aref diffusivity j)) (- (aref x j) xi))
            of-type double-float)))

(defun stencil-stable-timestep (stencil diffusivity &key (safety-factor 0.9d0) capacity)
  "Return SAFETY-FACTOR times the smallest area_i kappa_i / (2 *
STENCIL-TOTAL-CONDUCTANCE) over the cells of STENCIL with DIFFUSIVITY (one
value per cell) and CAPACITY kappa (one value per cell, default 1): the
explicit-diffusion stability bound. Cells with no conducting face do not
constrain it; infinity if none conducts."
  (let ((areas (flux-stencil-areas stencil)))
    (* safety-factor
       (reduce #'min (loop for i below (length areas)
                           for g = (stencil-total-conductance stencil diffusivity i)
                           when (plusp g)
                             collect (/ (* (aref areas i) (if capacity (aref capacity i) 1d0))
                                        (* 2.0d0 g)))
               :initial-value sb-ext:double-float-positive-infinity))))

(defun intensity-column (column capacity)
  "Return COLUMN (amounts per cell volume) divided element-wise by CAPACITY
- the quantity diffusion equalises, e.g. a pore-water concentration - as a
fresh vector; COLUMN itself when CAPACITY is NIL."
  (if capacity
      (map '(simple-array double-float (*)) #'/ column capacity)
      column))

(defun uniform-column (n value)
  "Return a fresh (simple-array double-float (N)) filled with VALUE."
  (make-array n :element-type 'double-float :initial-element (float value 1.0d0)))

(defun state-column (state s)
  "Return a fresh (simple-array double-float (n-cells)) copy of column S
\(one species) of STATE, an (array double-float (n-cells n-species))."
  (let* ((n (array-dimension state 0))
         (column (make-array n :element-type 'double-float)))
    (dotimes (i n column)
      (setf (aref column i) (aref state i s)))))
