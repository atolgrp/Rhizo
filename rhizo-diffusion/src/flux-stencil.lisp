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

(defun stencil-total-conductance (stencil i)
  "Return the sum of the conductances of every face of cell I in STENCIL."
  (let ((offsets (flux-stencil-offsets stencil))
        (conductances (flux-stencil-conductances stencil)))
    (loop for k from (aref offsets i) below (aref offsets (1+ i))
          sum (aref conductances k) of-type double-float)))

(defun stencil-net-flux (stencil state i s)
  "Return the net diffusive inflow into cell I of column S of STATE (an
\(array double-float (n-cells n-species))), per unit diffusivity:
sum over faces of conductance * (C_j - C_i). Divide by the cell's area and
multiply by the diffusivity to get dC_i/dt."
  (let ((offsets (flux-stencil-offsets stencil))
        (neighbors (flux-stencil-neighbors stencil))
        (conductances (flux-stencil-conductances stencil))
        (ci (aref state i s)))
    (declare (type (array double-float (* *)) state) (type double-float ci))
    (loop for k from (aref offsets i) below (aref offsets (1+ i))
          sum (* (aref conductances k) (- (aref state (aref neighbors k) s) ci))
            of-type double-float)))
