(in-package :rhizo-diffusion)

;;; Richards flow on the TPFA cells (PLAN.md, Phase 5).
;;;
;;; Mixed form, backward Euler in time, on the same two-point flux stencil
;;; as diffusion (flux-stencil.lisp). For cell i of area A_i,
;;;
;;;   R_i = A_i (theta_i(h) - theta_i^n) / dt
;;;         - sum_faces g_ij K_ij (H_j - H_i) - sum_boundary F_b = 0,
;;;
;;; with H = h + z the hydraulic head (z the elevation of the value point),
;;; g_ij = EDGE-LENGTH / FACE-DISTANCE and K_ij the face conductivity: the
;;; arithmetic mean of the two cells' K (the choice of Celia, Bouloutas &
;;; Zarba 1990) or the upstream one. Interior fluxes cancel in pairs, so
;;; sum_i R_i dt = (storage change) - (boundary inflow): the step conserves
;;; water exactly as well as R is driven to zero. That is Celia et al.'s
;;; argument for the mixed form; the :HEAD form, A_i C(h) (h - h^n) / dt,
;;; is kept only to reproduce their comparison - it does not conserve mass.
;;;
;;; Each step is solved by Newton's method on the full Jacobian (BICGSTAB:
;;; it is not symmetric), until max_i |R_i| dt / A_i - the water-content
;;; error per cell - is below a tolerance. RUN-RICHARDS adapts the step:
;;; it grows after easy steps and retries smaller after a failed one.
;;;
;;; Boundary conditions act on the faces BOUNDARY-FACES reports:
;;;   (:head h)            fixed pressure head at the face (two-point flux)
;;;   (:flux q)            imposed inflow per unit face length (q < 0: out)
;;;   (:free-drainage)     unit hydraulic gradient: outflow K(h_i) (bottom only)
;;;   (:atmospheric q h-min)
;;;                        flux q, unless the soil cannot deliver it with the
;;;                        surface at h-min (evaporation) or take it in with
;;;                        the surface at 0 (infiltration; excess runs off)
;;; h and q are numbers or functions of time. Elevation is a function of a
;;; point; by default -y, since y grows downward (depth).
;;;
;;; Assumed: a rigid soil, no hysteresis in theta(h), no ponded storage,
;;; an air phase at atmospheric pressure (Richards' assumption).

(defun depth-elevation (point)
  "Return the elevation -y of POINT (x . y): y is depth on every grid here."
  (- (cdr point)))

(defun time-value (value time)
  "Return VALUE at TIME: (FUNCALL VALUE TIME) for a function, VALUE
otherwise, as a double-float."
  (float (if (functionp value) (funcall value time) value) 1d0))

(defstruct (richards-boundary (:constructor %make-richards-boundary))
  (cell 0 :type fixnum)           ; dense index of the cell
  (length 0d0 :type double-float) ; face length
  (distance 0d0 :type double-float) ; value point to face
  (elevation 0d0 :type double-float) ; of the face point
  (kind :head :type keyword)
  (value 0d0)                     ; head or flux: a number or a function of time
  (limit 0d0))                    ; :ATMOSPHERIC: the lowest surface head

(defun boundary-point (point side distance)
  "Return the point DISTANCE from POINT toward the boundary SIDE (y down)."
  (destructuring-bind (x . y) point
    (ecase side
      (:top (cons x (- y distance)))
      (:bottom (cons x (+ y distance)))
      (:left (cons (- x distance) y))
      (:right (cons (+ x distance) y)))))

(defun make-richards-boundaries (tessellation cell-index spec elevation)
  "Return a fresh list of RICHARDS-BOUNDARY, one per boundary face of
TESSELLATION on the side named by SPEC, a list (side kind . arguments) -
see the file header. ELEVATION is a function of a point."
  (destructuring-bind (side kind &rest arguments) spec
    (let ((expected (case kind (:head 1) (:flux 1) (:free-drainage 0) (:atmospheric 2)
                      (t (error "RICHARDS boundary: unknown condition ~s" kind)))))
      (unless (= expected (length arguments))
        (error "RICHARDS boundary ~s: ~s takes ~d argument~:p" spec kind expected))
      (when (and (eq kind :free-drainage) (not (eq side :bottom)))
        (error "RICHARDS boundary: free drainage is a bottom condition, not ~s" side))
      (mapcar (lambda (face)
                (destructuring-bind (id length distance) face
                  (%make-richards-boundary
                   :cell (cell-id->index cell-index id)
                   :length (float length 1d0) :distance (float distance 1d0)
                   :elevation (float (funcall elevation
                                              (boundary-point (cell-position tessellation id)
                                                              side distance))
                                     1d0)
                   :kind kind
                   :value (if arguments (first arguments) 0d0)
                   :limit (if (eq kind :atmospheric) (float (second arguments) 1d0) 0d0))))
              (boundary-faces tessellation side)))))

(defstruct (richards-problem (:constructor %make-richards-problem))
  tessellation
  cell-index
  stencil      ; FLUX-STENCIL: conductances g_ij and areas
  elevations   ; (simple-array double-float (n)): z of each value point
  soils        ; simple-vector of VAN-GENUCHTEN, one per cell
  boundaries   ; list of RICHARDS-BOUNDARY
  face-mean    ; :ARITHMETIC or :UPSTREAM
  form)        ; :MIXED or :HEAD

(defun make-richards-problem (tessellation &key soil soil-field boundaries cell-index
                                                (elevation #'depth-elevation)
                                                (face-mean :arithmetic) (form :mixed))
  "Return a RICHARDS-PROBLEM on TESSELLATION. SOIL is one VAN-GENUCHTEN
soil for every cell, or SOIL-FIELD a function of a cell id returning one.
BOUNDARIES is a list of (side kind . arguments) (see the file header);
faces not named are no-flux. CELL-INDEX defaults to a fresh one (pass a
world's to couple flow to its transport). ELEVATION is a function of a
point. FACE-MEAN is :ARITHMETIC or :UPSTREAM; FORM :MIXED or :HEAD."
  (unless (member face-mean '(:arithmetic :upstream))
    (error "MAKE-RICHARDS-PROBLEM: unknown face mean ~s" face-mean))
  (unless (member form '(:mixed :head))
    (error "MAKE-RICHARDS-PROBLEM: unknown form ~s" form))
  (unless (alexandria:xor soil soil-field)
    (error "MAKE-RICHARDS-PROBLEM: give exactly one of SOIL and SOIL-FIELD"))
  (let* ((cell-index (or cell-index (make-cell-index tessellation)))
         (n (cell-index-count cell-index))
         (elevations (make-array n :element-type 'double-float))
         (soils (make-array n)))
    (dotimes (i n)
      (let ((id (index->cell-id cell-index i)))
        (setf (aref elevations i) (float (funcall elevation (cell-position tessellation id)) 1d0)
              (svref soils i) (if soil soil (funcall soil-field id)))))
    (%make-richards-problem
     :tessellation tessellation :cell-index cell-index
     :stencil (make-flux-stencil tessellation cell-index)
     :elevations elevations :soils soils
     :boundaries (mapcan (lambda (spec)
                           (make-richards-boundaries tessellation cell-index spec elevation))
                         boundaries)
     :face-mean face-mean :form form)))

(defstruct (richards-state (:constructor %make-richards-state))
  (time 0d0 :type double-float)
  heads            ; (simple-array double-float (n)) pressure heads
  water-contents   ; theta(h) in each cell
  face-volumes     ; per stencil face: cumulative volume out of the cell through it
  boundary-volumes ; per cell: cumulative net inflow through boundary faces
  next-dt)         ; the step RUN-RICHARDS would take next, or NIL

(defun richards-water-contents (problem heads)
  "Return a fresh vector of each cell's theta at HEADS under PROBLEM's soils."
  (let ((soils (richards-problem-soils problem))
        (theta (make-array (length heads) :element-type 'double-float)))
    (dotimes (i (length heads) theta)
      (setf (aref theta i) (vg-water-content (svref soils i) (aref heads i))))))

(defun make-richards-state (problem heads)
  "Return a RICHARDS-STATE at time 0 for PROBLEM with pressure HEADS: a
number for every cell, a function of a cell id, or a vector in the
problem's cell-index order. No water has crossed any face yet."
  (let* ((index (richards-problem-cell-index problem))
         (n (cell-index-count index))
         (h (make-array n :element-type 'double-float)))
    (dotimes (i n)
      (setf (aref h i) (float (cond ((functionp heads) (funcall heads (index->cell-id index i)))
                                    ((vectorp heads) (aref heads i))
                                    (t heads))
                              1d0)))
    (%make-richards-state
     :time 0d0 :heads h :water-contents (richards-water-contents problem h)
     :face-volumes (make-array (length (flux-stencil-neighbors (richards-problem-stencil problem)))
                               :element-type 'double-float :initial-element 0d0)
     :boundary-volumes (make-array n :element-type 'double-float :initial-element 0d0))))

(defun face-conductivity (face-mean k-i k-j dk-i dk-j dh)
  "Return (values K_face dK_face/dh_i dK_face/dh_j) for a face between
cells of conductivity K-I and K-J (derivatives DK-I, DK-J), with hydraulic
head difference DH = H_j - H_i: their mean for :ARITHMETIC, the upstream
cell's for :UPSTREAM."
  (ecase face-mean
    (:arithmetic (values (* 0.5d0 (+ k-i k-j)) (* 0.5d0 dk-i) (* 0.5d0 dk-j)))
    (:upstream (if (plusp dh) (values k-j 0d0 dk-j) (values k-i dk-i 0d0)))))

(defun boundary-head-inflow (face-mean boundary h-b soil h-i k-i dk-i z-i)
  "Return (values F dF/dh_i): the two-point inflow through BOUNDARY with
pressure head H-B at the face into a cell of SOIL at head H-I (K-I, DK-I,
elevation Z-I)."
  (let* ((g (/ (richards-boundary-length boundary) (richards-boundary-distance boundary)))
         (dh (- (+ h-b (richards-boundary-elevation boundary)) (+ h-i z-i))))
    (multiple-value-bind (k-face dk-face)
        (face-conductivity face-mean k-i (vg-conductivity soil h-b) dk-i 0d0 dh)
      (values (* g k-face dh) (* g (- (* dk-face dh) k-face))))))

(defun boundary-inflow (face-mean boundary time soil h-i k-i dk-i z-i)
  "Return (values F dF/dh_i): the inflow through BOUNDARY at TIME into its
cell, of SOIL at head H-I (K-I, DK-I, elevation Z-I) - see the file header."
  (let ((length (richards-boundary-length boundary))
        (value (richards-boundary-value boundary)))
    (ecase (richards-boundary-kind boundary)
      (:head (boundary-head-inflow face-mean boundary (time-value value time) soil h-i k-i dk-i z-i))
      (:flux (values (* length (time-value value time)) 0d0))
      (:free-drainage (values (- (* length k-i)) (- (* length dk-i))))
      (:atmospheric
       (let ((potential (* length (time-value value time))))
         (multiple-value-bind (low d-low)
             (boundary-head-inflow face-mean boundary (richards-boundary-limit boundary)
                                   soil h-i k-i dk-i z-i)
           (multiple-value-bind (high d-high)
               (boundary-head-inflow face-mean boundary 0d0 soil h-i k-i dk-i z-i)
             (cond ((< potential low) (values low d-low))
                   ((> potential high) (values high d-high))
                   (t (values potential 0d0))))))))))

(defun richards-system (problem h h-old theta-old dt time)
  "Return (values RESIDUAL DIAGONAL OFF-DIAGONAL FACE-FLUXES
BOUNDARY-FLUXES) of PROBLEM's backward Euler step of DT ending at TIME,
at heads H, from heads H-OLD with water contents THETA-OLD: the residual R
of the file header, the Jacobian dR/dh as a CSR matrix on the problem's
stencil (APPLY-CSR-MATRIX), the volumetric flux out of each cell through
each stencil face, and each cell's net inflow through boundary faces.
All fresh vectors."
  (let* ((stencil (richards-problem-stencil problem))
         (areas (flux-stencil-areas stencil))
         (offsets (flux-stencil-offsets stencil))
         (neighbors (flux-stencil-neighbors stencil))
         (conductances (flux-stencil-conductances stencil))
         (z (richards-problem-elevations problem))
         (soils (richards-problem-soils problem))
         (face-mean (richards-problem-face-mean problem))
         (n (length areas))
         (n-faces (length neighbors))
         (theta (make-array n :element-type 'double-float))
         (capacity (make-array n :element-type 'double-float))
         (k (make-array n :element-type 'double-float))
         (dk (make-array n :element-type 'double-float))
         (residual (make-array n :element-type 'double-float))
         (diagonal (make-array n :element-type 'double-float))
         (off (make-array n-faces :element-type 'double-float))
         (face-fluxes (make-array n-faces :element-type 'double-float))
         (boundary-fluxes (make-array n :element-type 'double-float :initial-element 0d0)))
    (declare (type (simple-array double-float (*)) h h-old theta-old z areas conductances))
    (dotimes (i n)
      (let ((soil (svref soils i)) (hi (aref h i)))
        (setf (aref theta i) (vg-water-content soil hi)
              (aref capacity i) (vg-capacity soil hi)
              (aref k i) (vg-conductivity soil hi)
              (aref dk i) (vg-conductivity-derivative soil hi))))
    (dotimes (i n)
      (let ((a/dt (/ (aref areas i) dt)))
        (ecase (richards-problem-form problem)
          (:mixed (setf (aref residual i) (* a/dt (- (aref theta i) (aref theta-old i)))))
          (:head (setf (aref residual i) (* a/dt (aref capacity i) (- (aref h i) (aref h-old i))))))
        (setf (aref diagonal i) (* a/dt (aref capacity i)))))
    (dotimes (i n)
      (let ((head-i (+ (aref h i) (aref z i))))
        (loop for f from (aref offsets i) below (aref offsets (1+ i))
              for j = (aref neighbors f)
              for g = (aref conductances f)
              for dh = (- (+ (aref h j) (aref z j)) head-i)
              do (multiple-value-bind (k-face dk-i dk-j)
                     (face-conductivity face-mean (aref k i) (aref k j) (aref dk i) (aref dk j) dh)
                   (let ((inflow (* g k-face dh)))
                     (decf (aref residual i) inflow)
                     (setf (aref face-fluxes f) (- inflow))
                     (incf (aref diagonal i) (* g (- k-face (* dk-i dh))))
                     (setf (aref off f) (- (* g (+ k-face (* dk-j dh))))))))))
    (dolist (boundary (richards-problem-boundaries problem))
      (let ((i (richards-boundary-cell boundary)))
        (multiple-value-bind (inflow d-inflow)
            (boundary-inflow face-mean boundary time (svref soils i)
                             (aref h i) (aref k i) (aref dk i) (aref z i))
          (decf (aref residual i) inflow)
          (decf (aref diagonal i) d-inflow)
          (incf (aref boundary-fluxes i) inflow))))
    (values residual diagonal off face-fluxes boundary-fluxes)))

(define-condition richards-nonconvergence (error)
  ((time :initarg :time :reader richards-nonconvergence-time)
   (iterations :initarg :iterations :reader richards-nonconvergence-iterations)
   (residual :initarg :residual :reader richards-nonconvergence-residual))
  (:report (lambda (c stream)
             (format stream "RICHARDS-STEP: no convergence by time ~a after ~d Newton iterations (water-content residual ~,2e)"
                     (richards-nonconvergence-time c) (richards-nonconvergence-iterations c)
                     (richards-nonconvergence-residual c))))
  (:documentation "Signalled by RICHARDS-STEP when Newton's method does not
reach its tolerance; RUN-RICHARDS retries with a smaller step."))

(defun safe-inverse-diagonal (diagonal)
  "Return a fresh vector of 1 / DIAGONAL_i, with 1 where DIAGONAL_i is 0:
a Jacobi preconditioner that cannot divide by zero."
  (map '(simple-array double-float (*)) (lambda (d) (if (zerop d) 1d0 (/ 1d0 d))) diagonal))

(defun richards-step (problem state dt &key (tolerance 1d-10) (max-iterations 15)
                                            (linear-tolerance 1d-10))
  "Return (values NEW-STATE ITERATIONS): STATE advanced by one backward
Euler step of DT under PROBLEM, solved by Newton's method until max_i
|R_i| DT / A_i <= TOLERANCE (a water-content error per cell). NEW-STATE is
fresh; its face and boundary volumes include this step's fluxes. Signals
RICHARDS-NONCONVERGENCE after MAX-ITERATIONS Newton iterations;
LINEAR-SOLVER-FAILURE and arithmetic errors pass through."
  (let* ((stencil (richards-problem-stencil problem))
         (areas (flux-stencil-areas stencil))
         (n (length areas))
         (h-old (richards-state-heads state))
         (theta-old (richards-state-water-contents state))
         (time (+ (richards-state-time state) dt))
         (h (copy-seq h-old)))
    (loop for iteration from 0
          do (multiple-value-bind (residual diagonal off face-fluxes boundary-fluxes)
                 (richards-system problem h h-old theta-old dt time)
               (let ((err (loop for i below n
                                maximize (abs (/ (* dt (aref residual i)) (aref areas i))))))
                 (when (<= err tolerance)
                   (return-from richards-step
                     (values (%make-richards-state
                              :time time :heads h
                              :water-contents (richards-water-contents problem h)
                              :face-volumes (vector-axpy dt face-fluxes
                                                         (richards-state-face-volumes state))
                              :boundary-volumes (vector-axpy dt boundary-fluxes
                                                             (richards-state-boundary-volumes state))
                              :next-dt (richards-state-next-dt state))
                             iteration)))
                 (when (>= iteration max-iterations)
                   (error 'richards-nonconvergence :time time :iterations iteration :residual err))
                 (let ((delta (bicgstab (lambda (v) (apply-csr-matrix stencil diagonal off v))
                                        (map '(simple-array double-float (*)) #'- residual)
                                        (make-array n :element-type 'double-float
                                                      :initial-element 0d0)
                                        :inverse-diagonal (safe-inverse-diagonal diagonal)
                                        :tolerance linear-tolerance :max-iterations 2000)))
                   (dotimes (i n) (incf (aref h i) (aref delta i)))))))))

(defun run-richards (problem state t-end &key dt (dt-min 1d-8) (dt-max most-positive-double-float)
                                              (growth 1.25d0) (shrink 0.5d0)
                                              (easy-iterations 3) (hard-iterations 8)
                                              (max-iterations 15) (tolerance 1d-10) observer)
  "Return (values STATE STATS): STATE advanced to time T-END under PROBLEM
by RICHARDS-STEPs with an adaptive step. The first step is DT, or the
step the state carries (RICHARDS-STATE-NEXT-DT). After a step that took
at most EASY-ITERATIONS Newton iterations the step grows by GROWTH (up to
DT-MAX); after one that took HARD-ITERATIONS or more it shrinks by
SHRINK. A failed step (no Newton or linear convergence, or an arithmetic
error) is retried SHRINK times smaller; an error is signalled below
DT-MIN. The last step is shortened to land exactly on T-END. OBSERVER, if
given, is called with each accepted state. STATS is a plist (:STEPS
:REJECTED :NEWTON-ITERATIONS :FINAL-DT). STATE itself is not modified."
  (let ((dt (float (or dt (richards-state-next-dt state)
                       (error "RUN-RICHARDS: no initial step DT")) 1d0))
        (t-end (float t-end 1d0))
        (steps 0) (rejected 0) (newton 0))
    (loop while (< (richards-state-time state) t-end)
          do (let* ((remaining (- t-end (richards-state-time state)))
                    (last-step (<= remaining dt))
                    (step (if last-step remaining dt)))
               (handler-case
                   (multiple-value-bind (next iterations)
                       (richards-step problem state step :tolerance tolerance
                                                         :max-iterations max-iterations)
                     (when last-step (setf (richards-state-time next) t-end))
                     (setf state next)
                     (incf steps)
                     (incf newton iterations)
                     (unless last-step
                       (cond ((<= iterations easy-iterations) (setf dt (min dt-max (* dt growth))))
                             ((>= iterations hard-iterations) (setf dt (max dt-min (* dt shrink))))))
                     (when observer (funcall observer state)))
                 ((or richards-nonconvergence linear-solver-failure arithmetic-error) ()
                   (incf rejected)
                   (setf dt (* shrink step))
                   (when (< dt dt-min)
                     (error "RUN-RICHARDS: the step fell below ~a at time ~a"
                            dt-min (richards-state-time state)))))))
    (let ((final (copy-richards-state state)))
      (setf (richards-state-next-dt final) dt)
      (values final (list :steps steps :rejected rejected :newton-iterations newton
                          :final-dt dt)))))

;;; Water balance.

(defun richards-storage (problem state)
  "Return the water held in PROBLEM's cells at STATE: sum of area x theta."
  (let ((areas (flux-stencil-areas (richards-problem-stencil problem)))
        (theta (richards-state-water-contents state)))
    (loop for i below (length areas) sum (* (aref areas i) (aref theta i)))))

(defun richards-cumulative-inflow (state)
  "Return the net volume that has entered through boundary faces up to STATE."
  (reduce #'+ (richards-state-boundary-volumes state)))

(defun richards-mass-balance-error (problem start end &key absolute)
  "Return the water-balance error between states START and END of PROBLEM:
(storage change - net boundary inflow) / net inflow, i.e. 1 minus the mass
balance ratio of Celia et al. (1990); the plain difference when ABSOLUTE.
Measured, not assumed: both terms come from the states."
  (let ((stored (- (richards-storage problem end) (richards-storage problem start)))
        (inflow (- (richards-cumulative-inflow end) (richards-cumulative-inflow start))))
    (if absolute
        (- stored inflow)
        (/ (- stored inflow) inflow))))

(defun richards-mass-balance-ratio (problem start end)
  "Return Celia et al.'s (1990) mass balance ratio between states START and
END of PROBLEM: storage change over net boundary inflow (1 when exact)."
  (/ (- (richards-storage problem end) (richards-storage problem start))
     (- (richards-cumulative-inflow end) (richards-cumulative-inflow start))))
