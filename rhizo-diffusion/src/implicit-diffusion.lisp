(in-package :rhizo-diffusion)

;;; Implicit (backward Euler) diffusion - PLAN.md, Phase 1, step 6.
;;;
;;; The explicit step is stable only below STABLE-TIMESTEP, which the
;;; smallest cell sets (a 1,600-cell nearest-seed grid has cells from 8 to
;;; over 400 pixels). A run that must span slow diffusion and a fast event
;;; cannot afford that. Backward Euler is stable at any step:
;;;
;;;     (A + dt K) c' = A c
;;;
;;; with A the diagonal of cell areas and K the symmetric two-point-flux
;;; matrix held by the FLUX-STENCIL, (K x)_i = sum_j D_ij g_ij (x_i - x_j),
;;; D_ij the harmonic mean of the two cells' diffusivities. The
;;; matrix is symmetric positive definite, so it is solved by conjugate
;;; gradient with a Jacobi (diagonal) preconditioner. It is also an
;;; M-matrix, so the step obeys a discrete maximum principle: no new extrema
;;; and no negative concentrations, at any dt.
;;;
;;; Measured, not assumed: mass. Every row of K sums to zero, so
;;; sum_i A_i c'_i = sum_i b_i - sum_i r_i, where r is the solver residual.
;;; Mass is conserved only as well as the linear system is solved, hence the
;;; tight default tolerance, and an error rather than a silent result if CG
;;; does not converge. Backward Euler is first order in time; Strang
;;; splitting with it is first order overall.

(defparameter *implicit-tolerance* 1.0d-15
  "Default relative residual ||b - M x|| / ||b|| at which CONJUGATE-GRADIENT
stops. Mass is conserved to about this relative accuracy per step (see file
header). Measured 2026-11-24 on a 1,600-cell nearest-seed grid: 1e-13 left
4e-15 mass error per step (1.4e-11 over 4,466 steps); 1e-15 leaves 8e-16 to
1e-15 per step, i.e. round-off, for 1-5 extra CG iterations. CG still
converged at 3e-16 there.")

(defun vector-dot (u v)
  "Return the dot product of double-float vectors U and V."
  (declare (type (simple-array double-float (*)) u v))
  (loop for i below (length u) sum (* (aref u i) (aref v i)) of-type double-float))

(defun conjugate-gradient (apply-matrix b x0 &key inverse-diagonal
                                                  (tolerance *implicit-tolerance*)
                                                  (max-iterations 10000))
  "Solve M x = B for a symmetric positive definite M given as APPLY-MATRIX,
a function from a double-float vector to a fresh vector M v. Starts from X0
and preconditions with INVERSE-DIAGONAL (a vector of 1/M_ii; Jacobi) when
given. Stops when ||B - M x|| <= TOLERANCE * ||B||. Returns (VALUES X
ITERATIONS RELATIVE-RESIDUAL), with X fresh; B and X0 are not modified.
Signals an error if MAX-ITERATIONS pass without convergence."
  (declare (type (simple-array double-float (*)) b x0))
  (let* ((n (length b))
         (b-norm (sqrt (vector-dot b b)))
         (x (copy-seq x0))
         (r (let ((mx (funcall apply-matrix x0)))
              (map '(simple-array double-float (*)) #'- b mx)))
         (precondition (lambda (v) (if inverse-diagonal
                                       (map '(simple-array double-float (*)) #'* inverse-diagonal v)
                                       (copy-seq v))))
         (z (funcall precondition r))
         (p (copy-seq z))
         (rz (vector-dot r z)))
    (declare (type (simple-array double-float (*)) x r z p) (type double-float rz))
    (when (zerop b-norm)
      (return-from conjugate-gradient
        (values (make-array n :element-type 'double-float :initial-element 0d0) 0 0d0)))
    (loop for iteration from 0
          for residual = (/ (sqrt (vector-dot r r)) b-norm)
          when (<= residual tolerance)
            do (return-from conjugate-gradient (values x iteration residual))
          when (>= iteration max-iterations)
            do (error "CONJUGATE-GRADIENT: no convergence after ~d iterations ~
(relative residual ~,2e, tolerance ~,2e)" iteration residual tolerance)
          do (let* ((mp (funcall apply-matrix p))
                    (alpha (/ rz (vector-dot p mp))))
               (declare (type (simple-array double-float (*)) mp) (type double-float alpha))
               (dotimes (i n)
                 (incf (aref x i) (* alpha (aref p i)))
                 (decf (aref r i) (* alpha (aref mp i))))
               (setf z (funcall precondition r))
               (let ((rz-new (vector-dot r z)))
                 (dotimes (i n)
                   (setf (aref p i) (+ (aref z i) (* (/ rz-new rz) (aref p i)))))
                 (setf rz rz-new))))))

(defun implicit-matrix-product (stencil dt diffusivity x &optional capacity)
  "Return a fresh vector (A kappa + DT K) X, with A the cell areas, kappa the
storage CAPACITY (one value per cell, default 1) and K the two-point-flux
matrix of STENCIL weighted by face diffusivity (harmonic mean of
DIFFUSIVITY, one value per cell) - the backward Euler matrix."
  (let* ((areas (flux-stencil-areas stencil))
         (n (length areas))
         (y (make-array n :element-type 'double-float)))
    (dotimes (i n y)
      (setf (aref y i) (- (* (aref areas i) (if capacity (aref capacity i) 1d0) (aref x i))
                          (* dt (stencil-net-flux stencil x diffusivity i)))))))

(defun implicit-inverse-diagonal (stencil dt diffusivity &optional capacity)
  "Return a fresh vector of 1 / (A_i kappa_i + DT * STENCIL-TOTAL-CONDUCTANCE):
the Jacobi preconditioner for IMPLICIT-MATRIX-PRODUCT with the same
arguments."
  (let* ((areas (flux-stencil-areas stencil))
         (inverse (make-array (length areas) :element-type 'double-float)))
    (dotimes (i (length areas) inverse)
      (setf (aref inverse i)
            (/ 1.0d0 (+ (* (aref areas i) (if capacity (aref capacity i) 1d0))
                        (* dt (stencil-total-conductance stencil diffusivity i))))))))

(defun implicit-diffusion-column (stencil dt diffusivity column
                                  &key (tolerance *implicit-tolerance*) capacity)
  "Return the amounts after one backward Euler diffusion step of DT from
COLUMN (one species' amounts per cell volume), with DIFFUSIVITY that
species' value in each cell and CAPACITY its storage capacity kappa (one
value per cell; NIL for 1): kappa x', where x' solves (A kappa + DT K) x' =
A COLUMN, x = m / kappa being what diffusion equalises. Second value: the
number of CG iterations."
  (let* ((areas (flux-stencil-areas stencil))
         (b (map '(simple-array double-float (*)) #'* areas column)))
    (multiple-value-bind (solution iterations)
        (conjugate-gradient (lambda (v) (implicit-matrix-product stencil dt diffusivity v capacity))
                            b (intensity-column column capacity)
                            :inverse-diagonal (implicit-inverse-diagonal stencil dt diffusivity
                                                                         capacity)
                            :tolerance tolerance)
      (values (if capacity (map '(simple-array double-float (*)) #'* capacity solution) solution)
              iterations))))

(defun implicit-diffusion-substep! (world dt &key (tolerance *implicit-tolerance*))
  "Advance every mobile species of WORLD by DT of backward Euler
finite-volume diffusion (see the file header), stable at any DT, with the
world's per-cell diffusivities. Species with zero diffusivity in every cell
pass through unchanged. Returns WORLD, and as
second value the largest number of CG iterations any species needed."
  (let* ((stencil (mechanistic-world-stencil world))
         (old-state (mechanistic-world-state world))
         (new-state (alexandria:copy-array old-state))
         (diffusivities (mechanistic-world-diffusivities world))
         (capacities (mechanistic-world-capacities world))
         (max-iterations 0))
    (dotimes (s (mechanistic-world-n-species world))
      (let ((diffusivity (state-column diffusivities s)))
        (when (some #'plusp diffusivity)
          (multiple-value-bind (column iterations)
              (implicit-diffusion-column stencil dt diffusivity (state-column old-state s)
                                         :tolerance tolerance
                                         :capacity (state-column capacities s))
            (setf max-iterations (max max-iterations iterations))
            (dotimes (i (length column))
              (setf (aref new-state i s) (aref column i)))))))
    (setf (mechanistic-world-state world) new-state)
    (values world max-iterations)))
