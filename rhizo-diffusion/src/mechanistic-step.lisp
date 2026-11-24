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

(see FACE-CONDUCTANCE, MAKE-FLUX-STENCIL) is computed, and the minimum
across every cell/species is returned, scaled by SAFETY-FACTOR. The largest
diffusivity always gives the smallest bound, so only it is used. This is
for species' uniform diffusivities; a world with a per-cell diffusivity
field uses WORLD-STABLE-TIMESTEP.
Immobile species (diffusivity 0.0d0) don't constrain anything and are
skipped. This is a heuristic bound - validated empirically by the CFL
boundary test in tests/mechanistic-tests.lisp (state stays bounded at/under
it, diverges meaningfully above it), not a proof of stability for the
coupled reaction-diffusion system as a whole, since reaction kinetics can
independently impose their own stiffness constraint that this doesn't see."
  (let* ((stencil (make-flux-stencil tessellation (make-cell-index tessellation)))
         (max-d (reduce #'max species :key #'species-diffusivity :initial-value 0.0d0)))
    (stencil-stable-timestep stencil (uniform-column (length (flux-stencil-areas stencil)) max-d)
                             :safety-factor safety-factor)))

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

(declaim (inline snap-underflow))
(defun snap-underflow (y)
  "Return Y, or 0.0d0 if Y is negative by less than
*REACTION-UNDERFLOW-EPSILON* (round-off near an exhausted pool; see there)."
  (if (and (minusp y) (> y (- *reaction-underflow-epsilon*))) 0.0d0 y))

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
        (setf (aref state i s) (snap-underflow y)))))
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

(defun euler-reaction-substep! (world dt)
  "Advance every cell's local reaction state by one explicit Euler step of
DT: c + DT f(c), f the NET-REACTION-DERIVATIVE! at the current state, with
SNAP-UNDERFLOW. Purely local per cell, so run in parallel. The reaction
half of the :IMEX scheme (STEP-MECHANISTIC!); positive and stable only for
DT below 1 / (the fastest reaction rate)."
  (let ((state (mechanistic-world-state world))
        (species-index (mechanistic-world-species-index world))
        (reactions (mechanistic-world-reactions world))
        (n-species (mechanistic-world-n-species world)))
    (pmap-cell-range!
     (cell-index-count (mechanistic-world-cell-index world))
     (lambda (i)
       (let ((y (make-array n-species :element-type 'double-float))
             (f (make-array n-species :element-type 'double-float)))
         (dotimes (s n-species) (setf (aref y s) (aref state i s)))
         (net-reaction-derivative! f y species-index reactions)
         (dotimes (s n-species)
           (setf (aref state i s) (snap-underflow (+ (aref y s) (* dt (aref f s))))))))))
  world)

(defun diffusion-substep! (world dt)
  "Advance every mobile species by DT of explicit finite-volume diffusion:
flux from neighbor J into cell I is D_face * FACE-CONDUCTANCE(I,J) *
\(x_J - x_I), x = m / kappa the amount over the storage capacity (the
world's CAPACITIES) and D_face the harmonic mean of the two cells'
diffusivities, summed over neighbors (STENCIL-NET-FLUX on the world's
precomputed FLUX-STENCIL) and divided by cell I's own area. Species with zero
diffusivity in every cell pass through unchanged. DT must respect
WORLD-STABLE-TIMESTEP. Reads only the
tick's starting state (double-buffered via a fresh array, then swapped in)
so it's safe to parallelize with PMAP-CELL-RANGE! even though every cell's
new value depends on its neighbors' old ones - the standard safe way to do
a synchronous CA-style update."
  (when (mechanistic-world-face-fluxes world)
    (error "DIFFUSION-SUBSTEP!: advection needs the :IMPLICIT or :IMEX scheme (implicit upwinding)"))
  (let* ((stencil (mechanistic-world-stencil world))
         (areas (flux-stencil-areas stencil))
         (n-cells (length areas))
         (n-species (mechanistic-world-n-species world))
         (diffusivities (mechanistic-world-diffusivities world))
         (capacities (mechanistic-world-capacities world))
         (old-state (mechanistic-world-state world))
         (new-state (make-array (list n-cells n-species) :element-type 'double-float)))
    (dotimes (i n-cells)
      (dotimes (s n-species) (setf (aref new-state i s) (aref old-state i s))))
    (dotimes (s n-species)
      (let ((diffusivity (state-column diffusivities s)))
        (when (some #'plusp diffusivity)
          (let ((column (intensity-column (state-column old-state s) (state-column capacities s))))
            (pmap-cell-range!
             n-cells
             (lambda (i)
               (incf (aref new-state i s)
                     (* (/ dt (aref areas i)) (stencil-net-flux stencil column diffusivity i)))))))))
    (setf (mechanistic-world-state world) new-state))
  world)

(defun scheme-diffusion-substep! (world dt)
  "Advance WORLD by DT of diffusion with the world's own DIFFUSION-SCHEME:
DIFFUSION-SUBSTEP! for :EXPLICIT, IMPLICIT-DIFFUSION-SUBSTEP! for :IMPLICIT
and :IMEX (whose diffusion is backward Euler). Returns WORLD."
  (ecase (mechanistic-world-diffusion-scheme world)
    (:explicit (diffusion-substep! world dt))
    ((:implicit :imex) (implicit-diffusion-substep! world dt)))
  world)

;;; The :IMEX scheme is not split. One explicit Euler reaction step followed
;;; by one backward Euler diffusion step is exactly IMEX Euler,
;;;
;;;     (A + dt K) c' = A (c + dt f(c)),
;;;
;;; whose fixed point is K c = A f(c): the steady state of the discrete
;;; equations, at any dt. Strang splitting does not have this property: its
;;; fixed point is off by O(dt), and by much more for modes that diffusion
;;; damps within one step (measured 2026-11-24: 11% error in the steady
;;; profile around a single producer at k dt = 0.1). Both are first order
;;; in time for transients (backward Euler); Strang's RK4 tolerates stiffer
;;; reactions.

(defun step-mechanistic! (world)
  "Advance WORLD by one tick of its own DT. For DIFFUSION-SCHEME :EXPLICIT
or :IMPLICIT, by Strang splitting - RK4 reaction half-step, full diffusion
step, reaction half-step (see the file header). For :IMEX, by IMEX Euler -
EULER-REACTION-SUBSTEP! then IMPLICIT-DIFFUSION-SUBSTEP!, which keeps
steady states exact (see above)."
  (let ((dt (mechanistic-world-dt world)))
    (if (eq (mechanistic-world-diffusion-scheme world) :imex)
        (progn (euler-reaction-substep! world dt)
               (implicit-diffusion-substep! world dt))
        (progn (reaction-substep! world (/ dt 2.0d0))
               (scheme-diffusion-substep! world dt)
               (reaction-substep! world (/ dt 2.0d0)))))
  (incf (mechanistic-world-tick world))
  world)

(defun run-mechanistic! (world n-ticks)
  "Advance WORLD by N-TICKS calls to STEP-MECHANISTIC! and return WORLD.
One lparallel kernel is held for the whole run (CALL-WITH-KERNEL), instead
of one being started and stopped by every PMAP-CELL-RANGE! call - three per
tick."
  (call-with-kernel (lambda () (dotimes (k n-ticks) (step-mechanistic! world))))
  world)
