(in-package :rhizo-diffusion)

;;; Advancing a MECHANISTIC-WORLD by one tick, via Strang (operator)
;;; splitting: a half-step of local reaction kinetics, a full step of
;;; diffusive transport, then another reaction half-step. Reaction (a system
;;; of local nonlinear ODEs, no neighbor reads) and diffusion (a spatial
;;; linear flux, no reaction terms) are different numerical processes;
;;; splitting them this way is the standard way to combine the two without
;;; either one degrading the other's accuracy, and it keeps each half's own
;;; parallel structure clean (see REACTION-SUBSTEP!/DIFFUSION-SUBSTEP!).

(defun stable-timestep (tessellation species &key (safety-factor 0.9d0))
  "A conservative CFL-style bound on the explicit diffusion timestep: for
every cell and every mobile species, the classic explicit-diffusion
stability condition

    area / (2 * diffusivity * sum-of-neighbor-conductances)

(conductance of an edge = edge-length / distance) is computed, and the
minimum across every cell/species is returned, scaled by SAFETY-FACTOR.
Immobile species (diffusivity 0.0d0) don't constrain anything and are
skipped. This is a heuristic bound - validated empirically by the CFL
boundary test in tests/mechanistic-tests.lisp (state stays bounded at/under
it, diverges meaningfully above it), not a proof of stability for the
coupled reaction-diffusion system as a whole, since reaction kinetics can
independently impose their own stiffness constraint that this doesn't see."
  (let ((bound sb-ext:double-float-positive-infinity))
    (dolist (sp species)
      (when (plusp (species-diffusivity sp))
        (dolist (id (cell-ids tessellation))
          (let ((area (cell-area tessellation id))
                (conductance 0.0d0))
            (dolist (nb (neighbors-of tessellation id))
              (incf conductance (/ (edge-length tessellation id nb)
                                    (edge-weight tessellation id nb))))
            (when (plusp conductance)
              (setf bound (min bound
                                (/ area (* 2.0d0 (species-diffusivity sp) conductance)))))))))
    (* safety-factor bound)))

(defparameter *reaction-underflow-epsilon* 1.0d-9
  "How far below zero RK4-CELL-REACTION! will tolerate before snapping a
result back to 0.0d0 - concentrations shouldn't go negative, and RK4 can
overshoot by a tiny floating-point amount near zero even when the true
dynamics wouldn't (e.g. a fast Monod/Michaelis-Menten term nearly exhausting
a pool within one step). Deliberately NOT a blanket (MAX 0.0d0 ...): a large
negative excursion is the actual, useful symptom of a genuinely unstable
timestep (see STABLE-TIMESTEP and its CFL boundary test), and silently
floor-clamping every negative value would hide that instability instead of
surfacing it - the same 'fail loud, not quiet' call made elsewhere in this
codebase.")

(defun rk4-cell-reaction! (state i n-species species-index reactions dt)
  "Advance row I of STATE by DT of local reaction kinetics via a fixed
4-stage Runge-Kutta step - deterministic and simple rather than an adaptive
stiff-ODE solver, which is more than v1 needs (see the plan's rationale;
revisit if the mass-conservation/analytic tests ever reveal real kinetics
too stiff for this)."
  (let ((y0 (make-array n-species :element-type 'double-float))
        (k1 (make-array n-species :element-type 'double-float))
        (k2 (make-array n-species :element-type 'double-float))
        (k3 (make-array n-species :element-type 'double-float))
        (k4 (make-array n-species :element-type 'double-float))
        (tmp (make-array n-species :element-type 'double-float)))
    (dotimes (s n-species) (setf (aref y0 s) (aref state i s)))
    (net-reaction-derivative! k1 y0 species-index reactions)
    (dotimes (s n-species) (setf (aref tmp s) (+ (aref y0 s) (* 0.5d0 dt (aref k1 s)))))
    (net-reaction-derivative! k2 tmp species-index reactions)
    (dotimes (s n-species) (setf (aref tmp s) (+ (aref y0 s) (* 0.5d0 dt (aref k2 s)))))
    (net-reaction-derivative! k3 tmp species-index reactions)
    (dotimes (s n-species) (setf (aref tmp s) (+ (aref y0 s) (* dt (aref k3 s)))))
    (net-reaction-derivative! k4 tmp species-index reactions)
    (dotimes (s n-species)
      (let ((y (+ (aref y0 s)
                  (* (/ dt 6.0d0)
                     (+ (aref k1 s) (* 2.0d0 (aref k2 s)) (* 2.0d0 (aref k3 s)) (aref k4 s))))))
        (setf (aref state i s)
              (if (and (minusp y) (> y (- *reaction-underflow-epsilon*))) 0.0d0 y)))))
  state)

(defun reaction-substep! (world dt)
  "Advance every cell's local reaction state by DT (see RK4-CELL-REACTION!)
- purely local per cell, no neighbor reads, so safe to run in parallel via
PMAP-CELL-RANGE!."
  (let ((state (mechanistic-world-state world))
        (species-index (mechanistic-world-species-index world))
        (reactions (mechanistic-world-reactions world))
        (n-species (mechanistic-world-n-species world))
        (n-cells (cell-index-count (mechanistic-world-cell-index world))))
    (pmap-cell-range!
     n-cells
     (lambda (i) (rk4-cell-reaction! state i n-species species-index reactions dt))))
  world)

(defun diffusion-substep! (world dt)
  "Advance every mobile species by DT of explicit finite-volume diffusion:
flux from neighbor J into cell I is D * edge-length(I,J)/distance(I,J) *
(C_J - C_I), summed over neighbors and divided by cell I's own area.
Immobile species (diffusivity 0.0d0) pass through unchanged. Reads only the
tick's starting state (double-buffered via a fresh array, then swapped in)
so it's safe to parallelize with PMAP-CELL-RANGE! even though every cell's
new value depends on its neighbors' old ones - the standard safe way to do
a synchronous CA-style update."
  (let* ((tess (mechanistic-world-tessellation world))
         (cell-index (mechanistic-world-cell-index world))
         (n-cells (cell-index-count cell-index))
         (n-species (mechanistic-world-n-species world))
         (species-list (mechanistic-world-species-list world))
         (old-state (mechanistic-world-state world))
         (new-state (make-array (list n-cells n-species) :element-type 'double-float)))
    (dotimes (i n-cells)
      (dotimes (s n-species) (setf (aref new-state i s) (aref old-state i s))))
    (pmap-cell-range!
     n-cells
     (lambda (i)
       (let ((id (index->cell-id cell-index i)))
         (dolist (nb-id (neighbors-of tess id))
           (let ((j (cell-id->index cell-index nb-id))
                 (conductance (/ (edge-length tess id nb-id) (edge-weight tess id nb-id)))
                 (area (cell-area tess id)))
             (loop for sp in species-list
                   for s from 0
                   for d = (species-diffusivity sp)
                   when (plusp d)
                     do (incf (aref new-state i s)
                              (* dt (/ (* d conductance) area)
                                 (- (aref old-state j s) (aref old-state i s))))))))))
    (setf (mechanistic-world-state world) new-state))
  world)

(defun step-mechanistic! (world)
  "Advance WORLD by one tick of its own DT via Strang splitting - reaction
half-step, full diffusion step, reaction half-step. See the file header."
  (let ((dt (mechanistic-world-dt world)))
    (reaction-substep! world (/ dt 2.0d0))
    (diffusion-substep! world dt)
    (reaction-substep! world (/ dt 2.0d0)))
  (incf (mechanistic-world-tick world))
  world)
