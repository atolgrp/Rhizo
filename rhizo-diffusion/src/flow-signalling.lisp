(in-package :rhizo-diffusion)

;;; Flow-driven water content and directional signalling (PLAN.md, Phase 5).
;;;
;;; Coupling. Each world step first advances Richards flow over the step
;;; (RUN-RICHARDS, adaptive substeps), then sets the world's water content
;;; to the new theta (SET-WATER-CONTENT!: capacities and diffusivities) and
;;; its flow field to the step's mean water fluxes (SET-FLOW-FIELD!), then
;;; takes the ordinary transport step. The transport therefore moves solute
;;; with exactly the water that moved: per cell, A (theta' - theta) = -dt x
;;; (net outflow), the same equation as the solute balance at uniform
;;; concentration. Assumed: rain carries no signal, and a cell's boundary
;;; faces are either all inflow or all outflow (true unless one cell
;;; touches both the top and the bottom of the domain).
;;;
;;; Directional calling distance. With a uniform Darcy flux q (bulk form,
;;; signalling.lisp: amounts per soil volume, flux D_s grad c + q c on the
;;; pore-water concentration c), the steady point source becomes
;;;
;;;   c = S / (2 pi D_s) exp(q x / (2 D_s)) K0(r sqrt(q^2 / (4 D_s^2) + k / D_s)),
;;;
;;; x measured along the flow. Writing k_eff = k + q^2 / (4 D_s), this is
;;; exp(q x / (2 D_s)) times the still-water solution with decay k_eff, and
;;; the same factorisation holds for the transient point source with
;;; capacity (retardation) C: the Galilean shift (x - q s / C) expands to
;;; exactly that factor and extra decay. So both directional solutions
;;; reuse POINT-SOURCE-STEADY-STATE and POINT-SOURCE-TRANSIENT. A signalling
;;; link i -> j exists when i's signal at j reaches j's threshold; with
;;; flow, a link can exist one way only.

(defun richards-mean-fluxes (start end)
  "Return (values FACE-FLUXES BOUNDARY-OUTFLOWS): the mean water flux
through each stencil face between Richards states START and END, and each
cell's mean outflow through boundary faces (net boundary inflow counted as
no outflow), as fresh vectors - a flow field for SET-FLOW-FIELD!."
  (let ((dt (- (richards-state-time end) (richards-state-time start))))
    (values (map '(simple-array double-float (*)) (lambda (v1 v0) (/ (- v1 v0) dt))
                 (richards-state-face-volumes end) (richards-state-face-volumes start))
            (map '(simple-array double-float (*)) (lambda (v1 v0) (max 0d0 (/ (- v0 v1) dt)))
                 (richards-state-boundary-volumes end) (richards-state-boundary-volumes start)))))

(defun step-flow-transport! (world problem state relation &key (blocked-p (constantly nil))
                                                              richards-options)
  "Advance WORLD by one step of its DT with Richards flow: PROBLEM from
STATE (which must be at WORLD's time and share its cell index) to the end
of the step, then water content, flow field and transport (see the file
header), with relative diffusivity RELATION and BLOCKED-P as for
SET-WATER-CONTENT!. RICHARDS-OPTIONS is a plist passed to RUN-RICHARDS.
Changes WORLD; returns the new Richards state."
  (let* ((dt (mechanistic-world-dt world))
         (index (mechanistic-world-cell-index world))
         (time (* (mechanistic-world-tick world) dt)))
    (unless (eq (richards-problem-cell-index problem) index)
      (error "STEP-FLOW-TRANSPORT!: the flow problem must use the world's cell index"))
    (unless (<= (abs (- (richards-state-time state) time)) (* 1d-9 (max 1d0 time)))
      (error "STEP-FLOW-TRANSPORT!: flow state at time ~a, world at ~a"
             (richards-state-time state) time))
    (let* ((next (apply #'run-richards problem state (+ time dt)
                        (append richards-options
                                (list :dt (or (richards-state-next-dt state) dt)))))
           (theta (richards-state-water-contents next)))
      (multiple-value-bind (faces outflows) (richards-mean-fluxes state next)
        (set-water-content! world (lambda (id) (aref theta (cell-id->index index id))) relation
                            :blocked-p blocked-p)
        (set-flow-field! world faces outflows))
      (step-mechanistic! world)
      next)))

(defun run-with-flow! (world problem state n-steps relation &key (blocked-p (constantly nil))
                                                                 richards-options observer)
  "Advance WORLD by N-STEPS calls to STEP-FLOW-TRANSPORT! from Richards
STATE; OBSERVER, if given, is called with WORLD and the flow state after
each step. Returns the final Richards state."
  (call-with-kernel
   (lambda ()
     (dotimes (k n-steps state)
       (setf state (step-flow-transport! world problem state relation :blocked-p blocked-p
                                                                      :richards-options richards-options))
       (when observer (funcall observer world state))))))

;;; Directional point sources (see the file header).

(defun advective-decay-rate (decay-rate flux diffusivity)
  "Return k + q^2 / (4 D): the decay rate that, times exp(q x / 2D), turns
the still-water point source into the one in a uniform Darcy FLUX q."
  (+ decay-rate (/ (* flux flux) (* 4 diffusivity))))

(defun advective-point-source-steady-state (source-rate diffusivity decay-rate flux along across)
  "Return the steady concentration at offset (ALONG, ACROSS) from a point
source of SOURCE-RATE, ALONG measured in the direction of the uniform
Darcy FLUX, with DIFFUSIVITY and DECAY-RATE (bulk form; file header)."
  (* (exp (/ (* flux along) (* 2 diffusivity)))
     (point-source-steady-state source-rate diffusivity
                                (advective-decay-rate decay-rate flux diffusivity)
                                (sqrt (+ (* along along) (* across across))))))

(defun advective-point-source-transient (source-rate diffusivity decay-rate capacity flux along across
                                         time)
  "Return the concentration at offset (ALONG, ACROSS) and TIME after a
point source of SOURCE-RATE switches on in a uniform Darcy FLUX, with
DIFFUSIVITY, DECAY-RATE and storage CAPACITY (local-equilibrium
retardation; POINT-SOURCE-TRANSIENT)."
  (* (exp (/ (* flux along) (* 2 diffusivity)))
     (point-source-transient source-rate diffusivity
                             (advective-decay-rate decay-rate flux diffusivity) capacity
                             (sqrt (+ (* along along) (* across across))) time)))

(defun directional-calling-distance (concentration threshold angle scale)
  "Return the distance, along the ray at ANGLE from the flow direction, at
which CONCENTRATION - a function of (along across) - falls to THRESHOLD;
SCALE is a typical length for the root search (CROSSING-DISTANCE)."
  (crossing-distance (lambda (r) (funcall concentration (* r (cos angle)) (* r (sin angle))))
                     threshold scale))

(defun flow-signal-field (parameters flux scenario measure)
  "Return (values CONCENTRATION THRESHOLD SCALE) for a producer patch of
CALLING-DISTANCE-SCENARIO with PARAMETERS and SCENARIO in a uniform Darcy
FLUX: CONCENTRATION a function of (along across) - the steady field for
MEASURE :STEADY, at the scenario's :TIME for :TRANSIENT - or NIL where the
signal cannot move (no connected water)."
  (let* ((inputs (calling-distance-scenario parameters scenario))
         (s (getf inputs :source-rate)) (d (getf inputs :diffusivity))
         (k (getf inputs :decay-rate)) (capacity (getf inputs :capacity))
         (time (scenario-value scenario :time)))
    (values (unless (zerop d)
              (ecase measure
                (:steady (lambda (along across)
                           (advective-point-source-steady-state s d k flux along across)))
                (:transient (lambda (along across)
                              (advective-point-source-transient s d k capacity flux along across
                                                                time)))))
            (getf inputs :threshold)
            (if (zerop d) 1d0 (sqrt (/ d (advective-decay-rate k flux d)))))))

(defun flow-signalling-edges (positions parameters flux scenario measure
                              &key (direction '(0d0 . 1d0)))
  "Return the directed signalling links (i . j) among equal producer
patches at POSITIONS under a uniform Darcy FLUX along the unit vector
DIRECTION (default +y, downward): i -> j when i's signal (FLOW-SIGNAL-FIELD)
at j reaches the receiver threshold. Pairs farther apart than the
downstream range cannot link and are not evaluated."
  (multiple-value-bind (concentration threshold scale)
      (flow-signal-field parameters flux scenario measure)
    (if (null concentration)
        '()
        (let ((r-max (directional-calling-distance concentration threshold 0d0 scale))
              (points (coerce positions 'simple-vector))
              (ux (car direction)) (uy (cdr direction)))
          (loop for i below (length points)
                nconc (loop for j below (length points)
                            for dx = (- (car (aref points j)) (car (aref points i)))
                            for dy = (- (cdr (aref points j)) (cdr (aref points i)))
                            when (and (/= i j)
                                      (<= (sqrt (+ (* dx dx) (* dy dy))) r-max)
                                      (>= (funcall concentration (+ (* dx ux) (* dy uy))
                                                   (- (* dy ux) (* dx uy)))
                                          threshold))
                              collect (cons i j)))))))

(defun one-way-fraction (edges)
  "Return the fraction of directed EDGES (i . j) whose reverse (j . i) is
not among them; 0 when there are none."
  (if (null edges)
      0
      (let ((table (make-hash-table :test 'equal)))
        (dolist (e edges) (setf (gethash e table) t))
        (/ (count-if-not (lambda (e) (gethash (cons (cdr e) (car e)) table)) edges)
           (length edges)))))

(defun directed-spanning-p (positions edges from-p to-p)
  "Return true if a directed path along EDGES leads from a node whose
position satisfies FROM-P to one whose position satisfies TO-P."
  (let* ((n (length positions))
         (points (coerce positions 'simple-vector))
         (out (make-array n :initial-element '()))
         (seen (make-array n :initial-element nil))
         (queue (loop for i below n when (funcall from-p (aref points i)) collect i)))
    (dolist (e edges) (push (cdr e) (aref out (car e))))
    (dolist (i queue) (setf (aref seen i) t))
    (loop while queue
          do (let ((i (pop queue)))
               (when (funcall to-p (aref points i)) (return-from directed-spanning-p t))
               (dolist (j (aref out i))
                 (unless (aref seen j)
                   (setf (aref seen j) t)
                   (push j queue)))))
    nil))

(defun flow-graph-connectivity (positions edges width height margin)
  "Return the connectivity of the directed signalling graph with nodes at
POSITIONS and EDGES on a WIDTH x HEIGHT domain with flow along +y, as a
plist: :LARGEST-WEAK-FRACTION, :LARGEST-STRONG-FRACTION, :MEAN-DEGREE and
:SPANNING-ACROSS (a weak component within MARGIN of the left and right
edges) from GRAPH-CONNECTIVITY; :SPANNING-ALONG (the same, top and bottom);
:SPANNING-ALONG-DIRECTED (a directed path from within MARGIN of the top to
within MARGIN of the bottom: relay downstream) and :ONE-WAY-FRACTION."
  (let* ((base (graph-connectivity positions edges width margin))
         (flipped (mapcar (lambda (p) (cons (cdr p) (car p))) positions))
         (weak (weak-components (length positions) edges)))
    (list :largest-weak-fraction (getf base :largest-weak-fraction)
          :largest-strong-fraction (getf base :largest-strong-fraction)
          :mean-degree (getf base :mean-degree)
          :spanning-across (getf base :spanning)
          :spanning-along (if (some (lambda (c) (spanning-p c flipped 0d0 height margin)) weak) t nil)
          :spanning-along-directed (directed-spanning-p positions edges
                                                        (lambda (p) (<= (cdr p) margin))
                                                        (lambda (p) (>= (cdr p) (- height margin))))
          :one-way-fraction (one-way-fraction edges))))

;;; Steady rain-fed water content. Under steady infiltration q into a deep
;;; uniform soil with free drainage the head gradient is zero, so the flux
;;; is carried by gravity alone: K(theta) = q. Below Ks this fixes theta
;;; from the rain rate - the simplest flow-driven water content.

(defun steady-drainage-head (soil flux)
  "Return the pressure head at which SOIL's conductivity equals FLUX (unit
hydraulic gradient), by bisection in log |h|; 0 when FLUX is at or above
Ks (saturated)."
  (if (>= flux (van-genuchten-ks soil))
      0d0
      (- (exp (bisect-root (lambda (u) (- (vg-conductivity soil (- (exp u))) flux))
                           (log (/ 1d-9 (van-genuchten-alpha soil)))
                           (log (/ 1d9 (van-genuchten-alpha soil))))))))

(defun steady-drainage-water-content (soil flux)
  "Return the water content of SOIL carrying FLUX under unit gradient
\(STEADY-DRAINAGE-HEAD)."
  (vg-water-content soil (steady-drainage-head soil flux)))

(defun flow-soil-parameters (parameters soil theta)
  "Return PARAMETERS (a signalling set, e.g. AHL-LITERATURE-PARAMETERS)
with :WATER-CONTENT THETA and the soil made consistent with the hydraulic
SOIL: :POROSITY = theta_s (Millington-Quirk) and :BULK-DENSITY = 2.65 (1 -
theta_s) g/cm^3 (sorption), replacing the flagged defaults."
  (let ((porosity (van-genuchten-theta-s soil)))
    (mapcar (lambda (p)
              (case (model-parameter-name p)
                (:water-content (make-model-parameter :water-content (float theta 1d0) "-"
                                                      :source "scenario: flow-driven water content"))
                (:porosity (make-model-parameter :porosity porosity "-"
                                                 :source "theta_s of the hydraulic soil (van Genuchten)"))
                (:bulk-density (make-model-parameter
                                :bulk-density (* 2.65d0 (- 1d0 porosity)) "g/cm^3"
                                :source "derived: (1 - theta_s) x 2.65 g/cm^3, the particle density of quartz-dominated mineral soil"
                                :note "porosity taken as theta_s"))
                (t p)))
            parameters)))

;;; Prescribed uniform flow and directional measurement.

(defun uniform-flow-field (world flux &key (direction '(0d0 . 1d0)))
  "Return (values FACE-FLUXES BOUNDARY-OUTFLOWS) of a uniform Darcy FLUX
along the unit vector DIRECTION (default +y, downward) on WORLD's cells: q
x EDGE-LENGTH x the component along DIRECTION of the unit vector between
value points, per stencil face, and q L max(n . DIRECTION, 0) over each
cell's boundary faces of outward normal n. Exactly divergence-free where
cell faces close exactly (rectangular grids); on pixelated Voronoi cells
only approximately - use a Richards solution there."
  (let* ((tess (mechanistic-world-tessellation world))
         (index (mechanistic-world-cell-index world))
         (stencil (mechanistic-world-stencil world))
         (offsets (flux-stencil-offsets stencil))
         (neighbors (flux-stencil-neighbors stencil))
         (conductances (flux-stencil-conductances stencil))
         (n (cell-index-count index))
         (positions (map 'vector (lambda (i) (cell-position tess (index->cell-id index i)))
                         (alexandria:iota n)))
         (faces (make-array (length neighbors) :element-type 'double-float))
         (outflows (make-array n :element-type 'double-float :initial-element 0d0))
         (ux (float (car direction) 1d0)) (uy (float (cdr direction) 1d0)))
    (dotimes (i n)
      (loop for k from (aref offsets i) below (aref offsets (1+ i))
            for pi* = (aref positions i)
            for pj = (aref positions (aref neighbors k))
            do (setf (aref faces k)
                     (* flux (aref conductances k)
                        (+ (* ux (- (car pj) (car pi*))) (* uy (- (cdr pj) (cdr pi*))))))))
    (loop for (side nx . ny) in '((:top 0d0 . -1d0) (:bottom 0d0 . 1d0)
                                  (:left -1d0 . 0d0) (:right 1d0 . 0d0))
          for out = (max 0d0 (+ (* nx ux) (* ny uy)))
          when (plusp out)
            do (loop for (id length) in (boundary-faces tess side)
                     do (incf (aref outflows (cell-id->index index id)) (* flux length out))))
    (values faces outflows)))

(defun sector-predicate (origin angle half-width &key (direction '(0d0 . 1d0)))
  "Return a predicate of a point: true when the point lies within
HALF-WIDTH (radians) of the ray from ORIGIN at ANGLE from the flow
DIRECTION (0 downstream, pi upstream, pi/2 across). For directional
profiles (RADIAL-PROFILE's KEEP)."
  (let ((ux (float (car direction) 1d0)) (uy (float (cdr direction) 1d0)))
    (lambda (point)
      (let* ((dx (- (car point) (car origin))) (dy (- (cdr point) (cdr origin)))
             (along (+ (* dx ux) (* dy uy))) (across (- (* dy ux) (* dx uy))))
        (and (or (/= dx 0) (/= dy 0))
             (let ((difference (- (atan across along) angle)))
               (<= (abs (atan (sin difference) (cos difference))) half-width)))))))
