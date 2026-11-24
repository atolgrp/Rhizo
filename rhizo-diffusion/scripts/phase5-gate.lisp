;;; Phase 5 gate report (PLAN.md): Richards flow.
;;;
;;; 1. Celia, Bouloutas & Zarba (1990) infiltration: 100 cm of sand,
;;;    h = -1000 cm initially, -75 cm at the top, -1000 cm at the bottom,
;;;    profiles at 1 day. Gate: mass-balance error under 0.1%.
;;;    a. reference: mixed form, dz = 0.1 cm, adaptive steps
;;;    b. mixed vs head form at dz = 1 cm and fixed steps (Celia et al.'s
;;;       comparison): mass balance ratio and wetting-front depth
;;;    c. the same column on an irregular nearest-seed Voronoi slab
;;;       (horizontal mean against the 1D reference)
;;; 2. Repeat of Phase 4 with flow-driven water content and directional
;;;    (advective) edges - see the second half of this file.
;;; Every random grid and layout is seeded; seeds are printed.
;;;
;;; Run with:
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase5-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

;;; RHIZO_GATE_QUICK=1 shrinks every section (smoke test of the script only;
;;; its numbers are not the gate's).
(defparameter *quick* (uiop:getenv "RHIZO_GATE_QUICK"))
(defmacro quick-or (small full) `(if *quick* ,small ,full))

(defun celia-problem (tess &key (form :mixed))
  (make-richards-problem tess :soil (van-genuchten-soil (celia-sand-parameters))
                              :boundaries '((:top :head -75d0) (:bottom :head -1000d0))
                              :form form))

(defun column-profile (problem state)
  "List of (depth theta) of a one-cell-wide column, top to bottom."
  (let ((index (richards-problem-cell-index problem))
        (tess (richards-problem-tessellation problem)))
    (sort (loop for i below (cell-index-count index)
                collect (list (cdr (cell-position tess (index->cell-id index i)))
                              (aref (richards-state-water-contents state) i)))
          #'< :key #'first)))

(defun front-depth (profile)
  "Depth at which theta crosses halfway between its top and bottom values,
by linear interpolation in PROFILE ((depth theta) top to bottom)."
  (let* ((top (second (first profile))) (bottom (second (car (last profile))))
         (mid (* 0.5d0 (+ top bottom))))
    (loop for ((y0 t0) (y1 t1)) on profile
          while y1
          when (and (>= t0 mid) (< t1 mid))
            do (return (+ y0 (* (- y1 y0) (/ (- t0 mid) (- t0 t1))))))))

(defun profile-at (profile depth)
  "Theta of PROFILE at DEPTH by linear interpolation."
  (loop for ((y0 t0) (y1 t1)) on profile
        while y1
        when (<= y0 depth y1)
          do (return (+ t0 (* (- t1 t0) (/ (- depth y0) (- y1 y0)))))
        finally (return (second (car (last profile))))))

(defun timed (thunk)
  (let ((t0 (get-internal-real-time)))
    (multiple-value-call #'values
      (funcall thunk)
      (/ (- (get-internal-real-time) t0) internal-time-units-per-second 1.0))))

(format t "~&== 1. Celia, Bouloutas & Zarba (1990) infiltration, profiles at 1 day ==~%")
(format t "  sand: theta_r 0.102, theta_s 0.368, alpha 0.0335 /cm, n 2, Ks 0.00922 cm/s; 100 cm column~%")
(defparameter *day* 86400d0)

(defparameter *reference*
  (let* ((problem (celia-problem (quick-or (make-rectangular-tessellation 1 200 :dy 0.5d0)
                                            (make-rectangular-tessellation 1 1000 :dy 0.1d0))))
         (start (make-richards-state problem -1000d0)))
    (multiple-value-bind (end stats seconds)
        (timed (lambda () (run-richards problem start *day* :dt 0.1d0 :dt-max 60d0)))
      (let ((profile (column-profile problem end)))
        (format t "  a. reference, dz 0.1 cm, adaptive (dt <= 60 s): ~d steps, ~d rejected, ~,1f s~%"
                (getf stats :steps) (getf stats :rejected) seconds)
        (format t "     mass balance ratio ~,8f (error ~,2e); front at ~,2f cm; inflow ~,4f cm^3/cm^2~%"
                (richards-mass-balance-ratio problem start end)
                (richards-mass-balance-error problem start end)
                (front-depth profile) (richards-cumulative-inflow end))
        profile))))
(finish-output)

(defun profile-error (profile)
  "Largest |theta - theta_reference| over PROFILE's depths."
  (loop for (y th) in profile maximize (abs (- th (profile-at *reference* y)))))

(format t "  b. dz 1 cm, fixed steps (halved after a failed Newton solve, then grown back):~%")
(format t "     form   dt(s)  steps  rejected   MB ratio      MB error   front(cm)  max|dtheta| vs a.~%")
(dolist (form '(:mixed :head))
  (dolist (dt '(10d0 60d0 600d0 1800d0))
    (let* ((problem (celia-problem (make-rectangular-tessellation 1 100 :dy 1d0) :form form))
           (start (make-richards-state problem -1000d0)))
      (multiple-value-bind (end stats)
          (run-richards problem start *day* :dt dt :dt-max dt :growth 2d0
                                            :easy-iterations 1000 :hard-iterations 1000
                                            :max-iterations 50)
        (let ((profile (column-profile problem end)))
          (format t "     ~6a ~6,0f ~6d ~6d   ~11,8f  ~10,2e  ~8,2f   ~8,4f~%"
                  form dt (getf stats :steps) (getf stats :rejected)
                  (richards-mass-balance-ratio problem start end)
                  (richards-mass-balance-error problem start end)
                  (front-depth profile) (profile-error profile))))
      (finish-output))))

(let* ((problem (celia-problem (make-rectangular-tessellation 1 100 :dy 1d0)))
       (start (make-richards-state problem -1000d0)))
  (multiple-value-bind (end stats seconds)
      (timed (lambda () (run-richards problem start *day* :dt 1d0)))
    (let ((profile (column-profile problem end)))
      (format t "     mixed, adaptive (no dt cap): ~d steps (~d rejected, ~d Newton iterations, ~,2f s), final dt ~,0f s~%"
              (getf stats :steps) (getf stats :rejected) (getf stats :newton-iterations) seconds
              (getf stats :final-dt))
      (format t "       MB error ~,2e, front ~,2f cm, max|dtheta| ~,4f~%"
              (richards-mass-balance-error problem start end) (front-depth profile)
              (profile-error profile)))))
(finish-output)

(format t "  c. the same column on a nearest-seed Voronoi slab (20 cm x 100 cm, 0.5 cm pixels):~%")
(dolist (n (quick-or '(400) '(400 1600 6400)))
  (let* ((grid (make-nearest-seed-grid 40 200 n :seed 11))
         (tess (voronoi-grid-tessellation grid))
         ;; lengths in pixels of 0.5 cm: 2 length units per cm (heads too)
         (problem (make-richards-problem tess
                                         :soil (van-genuchten-soil (celia-sand-parameters)
                                                                   :length-scale 2d0)
                                         :boundaries '((:top :head -150d0) (:bottom :head -2000d0))))
         (start (make-richards-state problem -2000d0)))
    (multiple-value-bind (end stats seconds)
        (timed (lambda () (run-richards problem start *day* :dt 1d0 :dt-max 600d0)))
      (let* ((index (richards-problem-cell-index problem))
             (bins (make-array 50 :initial-element nil)))
        ;; horizontal mean of theta in 2 cm depth bins
        (dotimes (i (cell-index-count index))
          (let ((depth-cm (* 0.5d0 (+ 0.5d0 (cdr (cell-position tess (index->cell-id index i)))))))
            (push (aref (richards-state-water-contents end) i)
                  (aref bins (min 49 (floor depth-cm 2))))))
        (let ((profile (loop for b below 50
                             when (aref bins b)
                               collect (list (+ 1d0 (* 2 b))
                                             (/ (reduce #'+ (aref bins b)) (length (aref bins b)))))))
          (format t "     ~4d cells (seed 11): ~d steps, ~,1f s; MB error ~,2e; front ~,2f cm (1D reference ~,2f); max|dtheta| of 2 cm means ~,4f~%"
                  n (getf stats :steps) seconds
                  (richards-mass-balance-error problem start end)
                  (front-depth profile) (front-depth *reference*) (profile-error profile)))))
    (finish-output)))

;;; ------------------------------------------------------------------
;;; Data for the report's graphs: written as JSON when RHIZO_REPORT_DIR is
;;; set (objects are alists, arrays vectors).

(defparameter *report-dir* (uiop:getenv "RHIZO_REPORT_DIR"))

(defun dump-json (name object)
  (when *report-dir*
    (with-open-file (s (merge-pathnames (format nil "~a.json" name)
                                        (uiop:ensure-directory-pathname *report-dir*))
                       :direction :output :if-exists :supersede :if-does-not-exist :create)
      (cl-json:encode-json object s))))

(defun vec (list) (coerce list 'vector))

(dump-json "celia-reference" `((:depth . ,(vec (mapcar #'first *reference*)))
                               (:theta . ,(vec (mapcar #'second *reference*)))))
(dump-json "celia-forms"
           (vec (loop for form in '(:mixed :head)
                      nconc (loop for dt in '(60d0 600d0 1800d0)
                                  collect (let* ((problem (celia-problem (make-rectangular-tessellation 1 100 :dy 1d0)
                                                                         :form form))
                                                 (start (make-richards-state problem -1000d0))
                                                 (end (run-richards problem start *day* :dt dt :dt-max dt :growth 2d0
                                                                    :easy-iterations 1000 :hard-iterations 1000
                                                                    :max-iterations 50))
                                                 (profile (column-profile problem end)))
                                            `((:form . ,(string-downcase form)) (:dt . ,dt)
                                              (:mass-balance-error . ,(richards-mass-balance-error problem start end))
                                              (:depth . ,(vec (mapcar #'first profile)))
                                              (:theta . ,(vec (mapcar #'second profile)))))))))

;;; ------------------------------------------------------------------
(format t "~&== 2. Phase 4 repeated with flow-driven water content and directional links ==~%")
(format t "  soil: sandy clay loam (Carsel & Parrish 1988), um units; signalling porosity = theta_s 0.39,~%")
(format t "  bulk density 2.65 (1 - theta_s); Millington-Quirk; AHL literature parameters; flow along +y (down)~%")
(defparameter *loam* (van-genuchten-soil (sandy-clay-loam-parameters) :length-scale 1d4))
(defparameter *fluxes* '(1d-5 3d-5 1d-4 3d-4 1d-3 2d-3 3d-3 5d-3 1d-2 3d-2 0.1d0 0.3d0 1d0 2d0 3d0))
(defparameter *biotic* '(:degrader-density 1d9))
(defparameter *abiotic* '())

(defun loam-parameters (theta)
  (flow-soil-parameters (ahl-literature-parameters) *loam* theta))

(format t "~&  2a. steady rain q: water content (K(theta) = q), D_s, Peclet over the still range, directional ranges (um)~%")
(format t "      q(um/s) mm/day  theta   D_s  | biotic steady: Pe  down  across  up | abiotic day one: Pe  down  across  up~%")
(defparameter *range-table*
  (loop for q in *fluxes*
        collect (let* ((theta (steady-drainage-water-content *loam* q))
                       (p (loam-parameters theta)))
                  (flet ((ranges (scenario measure)
                           (multiple-value-bind (c thr scale) (flow-signal-field p q scenario measure)
                             (let ((still (getf (calling-distance-scenario p scenario) measure))
                                   (d (getf (calling-distance-scenario p scenario) :diffusivity)))
                               (list* (edge-peclet q d still)
                                      (mapcar (lambda (a) (directional-calling-distance c thr a scale))
                                              (list 0d0 (/ pi 2) pi)))))))
                    (let ((b (ranges *biotic* :steady)) (a (ranges *abiotic* :transient))
                          (d (getf (calling-distance-scenario p *biotic*) :diffusivity)))
                      (format t "      ~8,5f ~6,2f  ~5,3f ~6,1f | ~6,3f ~6,1f ~6,1f ~6,1f | ~6,3f ~7,1f ~6,1f ~6,1f~%"
                              q (* q 86.4d0) theta d (first b) (second b) (third b) (fourth b)
                              (first a) (second a) (third a) (fourth a))
                      (finish-output)
                      `((:flux . ,q) (:theta . ,theta) (:diffusivity . ,d)
                        (:biotic . ,(vec b)) (:abiotic . ,(vec a))))))))
(dump-json "flow-ranges" (vec *range-table*))

(format t "~&  2b. engine check: one induced producer under Richards-driven steady rain, q = 1 um/s, 1e9 degraders/L~%")
(format t "      1600 um square (nearest-seed, seed 1); top flux, free drainage; steady state; ranges by 30-degree sectors~%")
(let* ((q 1d0)
       (theta (steady-drainage-water-content *loam* q))
       (base (loam-parameters theta))
       (p (with-parameter-values base
            `((:basal-production . ,(+ (parameter-value base :basal-production)
                                       (parameter-value base :induced-production)))
              (:induced-production . 0d0))))
       (analytic (multiple-value-bind (c thr scale) (flow-signal-field base q *biotic* :steady)
                   (mapcar (lambda (a) (directional-calling-distance c thr a scale))
                           (list 0d0 (/ pi 2) pi)))))
  (dolist (n (quick-or '(400) '(1600 6400)))
    (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 1600 1600 n :seed 1)))
           (centre '(800d0 . 500d0))
           (source (nearest-cell tess centre))
           (density (/ (* 1000d0 5d14) (cell-area tess source)))
           (dt (/ 0.5d0 (+ (parameter-value p :hydrolysis-rate) (sorption-relaxation-rate p))))
           (world (make-ahl-world tess :parameters p :dt dt
                                       :initial-state (lambda (id)
                                                        `((:degrader . 1d9)
                                                          ,@(when (equal id source)
                                                              `((:producer . ,density)))))))
           (index (mechanistic-world-cell-index world))
           (problem (make-richards-problem tess :soil *loam* :cell-index index
                                                :boundaries `((:top :flux ,q) (:bottom :free-drainage))))
           (h (steady-drainage-head *loam* q))
           (start (make-richards-state problem h))
           (steady (run-richards problem start 4000d0 :dt 10d0))
           (t0 (get-internal-real-time)))
      (multiple-value-bind (faces outflows)
          (let ((later (run-richards problem steady 5000d0)))
            (richards-mean-fluxes steady later))
        (set-water-content! world (let ((th (richards-state-water-contents steady)))
                                    (lambda (id) (aref th (cell-id->index index id))))
                            (relative-diffusivity-model :millington-quirk p))
        (set-flow-field! world faces outflows))
      (let ((balance0 (signal-pool-balance world :ahl)))
        (multiple-value-bind (steps change outcome) (step-to-steady-state! world :ahl :max-steps 200000)
          (let ((simulated (mapcar (lambda (a)
                                     (calling-distance world (pore-water :ahl) (cell-position tess source) 70d0
                                                       :keep (sector-predicate (cell-position tess source) a
                                                                               (/ pi 12))))
                                   (list 0d0 (/ pi 2) pi))))
            (format t "      ~4d cells: ~d steps (~a, change ~,1e), ~,0f s; theta ~,3f~%"
                    n steps outcome change
                    (/ (- (get-internal-real-time) t0) internal-time-units-per-second) theta)
            (loop for label in '("downstream" "across" "upstream")
                  for s in simulated for a in analytic
                  do (format t "        ~10a simulated ~6,1f um  analytic ~6,1f um  (~5,1f%)~%"
                             label s a (* 100 (- (/ s a) 1))))
            (format t "        mass balance (signal pools + exported, over produced): ~,1e; exported ~,3e of ~,3e produced~%"
                    (/ (abs (- (signal-pool-balance world :ahl) balance0)) (species-mass world :ahl-produced))
                    (species-exported world :ahl) (species-mass world :ahl-produced))
            (finish-output)))))))

(format t "~&  2c. percolation under steady rain: 600 patches, 4 mm square, seeds 1-20; margin 400 um~%")
(format t "      fraction of layouts spanning across the flow (x), along it (y), directed downstream (y, relay);~%")
(format t "      mean largest weak / strong component; one-way links. 'still' = same theta, no advection.~%")
(defparameter *seeds* (quick-or '(1 2) (alexandria:iota 20 :start 1)))
(defun mean (xs) (/ (reduce #'+ xs) (max 1 (length xs))))
(defun fraction-true (xs) (/ (count-if #'identity xs) (max 1 (length xs))))
(defparameter *percolation*
  (loop for (label scenario measure) in `(("1e9 degraders/L, steady" ,*biotic* :steady)
                                          ("abiotic, day one" ,*abiotic* :transient))
        collect
        (progn
          (format t "      ~a~%" label)
          (format t "      q(um/s)  theta   Pe   | flow: across along directed  weak strong one-way | still: across along weak~%")
          (cons label
                (vec
                 (loop for q in (quick-or '(1d-3 1d0) *fluxes*)
                       collect
                       (let* ((samples (mapcar (lambda (seed)
                                                 (flow-connectivity-sample
                                                  (ahl-literature-parameters) *loam* q
                                                  (random-patch-positions 600 4000d0 4000d0 :seed seed)
                                                  4000d0 4000d0 400d0 scenario measure))
                                               *seeds*))
                              (flow (mapcar (lambda (s) (getf s :flow)) samples))
                              (still (mapcar (lambda (s) (getf s :still)) samples)))
                         (flet ((f (rows key) (float (fraction-true (mapcar (lambda (r) (getf r key)) rows))))
                                (m (rows key) (float (mean (mapcar (lambda (r) (getf r key)) rows)))))
                           (format t "      ~8,5f ~5,3f ~6,3f | ~5,2f ~5,2f ~5,2f  ~5,2f ~5,2f ~5,2f | ~5,2f ~5,2f ~5,2f~%"
                                   q (getf (first samples) :theta) (getf (first samples) :peclet)
                                   (f flow :spanning-across) (f flow :spanning-along) (f flow :spanning-along-directed)
                                   (m flow :largest-weak-fraction) (m flow :largest-strong-fraction)
                                   (m flow :one-way-fraction)
                                   (f still :spanning-across) (f still :spanning-along)
                                   (m still :largest-weak-fraction))
                           (finish-output)
                           `((:flux . ,q) (:theta . ,(getf (first samples) :theta))
                             (:peclet . ,(getf (first samples) :peclet))
                             (:flow-across . ,(f flow :spanning-across)) (:flow-along . ,(f flow :spanning-along))
                             (:flow-directed . ,(f flow :spanning-along-directed))
                             (:flow-weak . ,(m flow :largest-weak-fraction))
                             (:flow-strong . ,(m flow :largest-strong-fraction))
                             (:one-way . ,(m flow :one-way-fraction))
                             (:still-across . ,(f still :spanning-across))
                             (:still-along . ,(f still :spanning-along))
                             (:still-weak . ,(m still :largest-weak-fraction)))))))))))
(dump-json "flow-percolation" *percolation*)

(format t "~&  2d. hysteresis through weather (tagged simulation, Richards-driven theta):~%")
(format t "      1.2 mm square, 1600 cells (nearest-seed, seed 1), 12 patches (layout seed 3), 1e9 degraders/L,~%")
(format t "      constant induced production; top: rain 2 mm/h for 2 h, then evaporation 4 mm/day (limit -15000 cm);~%")
(format t "      bottom held at -50 cm (moist subsoil); dt 30 s; 36 h; samples every 30 min~%")
(let* ((side 1200d0)
       (tess (voronoi-grid-tessellation (make-nearest-seed-grid 1200 1200 1600 :seed 1)))
       (positions (random-patch-positions 12 side side :seed 3))
       (nodes (mapcar (lambda (p) (nearest-cell tess p)) positions))
       (cell-positions (mapcar (lambda (id) (cell-position tess id)) nodes))
       (weather (rain-evaporation-schedule :rain-rate (/ 2000d0 3600d0) :rain-time 7200d0
                                           :evaporation-rate (/ 4000d0 86400d0) :period 86400d0))
       (base (loam-parameters 0.3d0))
       (parameters (with-parameter-values base
                     `((:basal-production . ,(+ (parameter-value base :basal-production)
                                                (parameter-value base :induced-production)))
                       (:induced-production . 0d0))))
       (runs '()))
  (dolist (advect '(t nil))
    (let* ((world (make-tagged-ahl-world tess nodes :parameters parameters :dt 30d0 :feedback nil
                                                    :degrader-density 1d9 :producer-density 5d14))
           (problem (make-richards-problem tess :soil *loam* :cell-index (mechanistic-world-cell-index world)
                                                :boundaries `((:top :atmospheric ,weather -1.5d8)
                                                              (:bottom :head -5d5))))
           (state (make-richards-state problem (lambda (id) (+ -5d5 -1200d0 (cdr (cell-position tess id))))))
           (balance0 (tagged-pool-balance world))
           (t0 (get-internal-real-time)))
      (multiple-value-bind (samples end)
          (connectivity-over-flow world problem state nodes cell-positions :millington-quirk
                                  (quick-or 240 (round (* 36 3600) 30)) 60 70d0 side side :advect advect)
        (declare (ignore end))
        (push (cons advect samples) runs)
        (format t "~&      ~:[without advection (theta from the flow only)~;with advection~]: ~,0f s, mass balance ~,1e~%"
                advect (/ (- (get-internal-real-time) t0) internal-time-units-per-second)
                (/ (abs (- (tagged-pool-balance world) balance0))
                   (max 1d-300 (species-mass world :tags-produced))))
        (let ((second-cycle (or (remove-if (lambda (s) (< (getf s :time) 86400d0)) samples)
                                samples))) ; quick mode has no second cycle
          (format t "        hysteresis gap (second cycle): mean degree ~,3f, largest weak ~,3f, one-way ~,3f~%"
                  (hysteresis-gap second-cycle :key :mean-degree)
                  (hysteresis-gap second-cycle :key :largest-weak-fraction)
                  (hysteresis-gap second-cycle :key :one-way-fraction)))
        (format t "        time(h) theta  inflow(um^2/s)  weak  strong  degree  one-way  across along directed~%")
        (dolist (s samples)
          (when (zerop (mod (round (getf s :time)) 7200))
            (format t "        ~5,1f  ~5,3f  ~10,3f   ~5,3f  ~5,3f  ~6,2f  ~6,3f    ~3a   ~3a   ~3a~%"
                    (/ (getf s :time) 3600) (getf s :theta) (getf s :inflow)
                    (float (getf s :largest-weak-fraction)) (float (getf s :largest-strong-fraction))
                    (float (getf s :mean-degree)) (float (getf s :one-way-fraction))
                    (if (getf s :spanning-across) "yes" "no") (if (getf s :spanning-along) "yes" "no")
                    (if (getf s :spanning-along-directed) "yes" "no"))))
        (finish-output))))
  (dump-json "flow-hysteresis"
             (vec (loop for (advect . samples) in runs
                        collect `((:advect . ,(if advect t :false))
                                  (:samples . ,(vec (mapcar (lambda (s)
                                                              `((:time . ,(getf s :time)) (:theta . ,(getf s :theta))
                                                                (:inflow . ,(getf s :inflow))
                                                                (:weak . ,(float (getf s :largest-weak-fraction)))
                                                                (:strong . ,(float (getf s :largest-strong-fraction)))
                                                                (:degree . ,(float (getf s :mean-degree)))
                                                                (:one-way . ,(float (getf s :one-way-fraction)))))
                                                            samples))))))))

(format t "~&  2e. Peclet and Damkohler numbers of the links (biotic steady, layout seed 1, 600 patches):~%")
(dolist (q '(1d-3 0.1d0 1d0))
  (let* ((theta (steady-drainage-water-content *loam* q))
         (p (loam-parameters theta))
         (positions (random-patch-positions 600 4000d0 4000d0 :seed 1))
         (edges (flow-signalling-edges positions p q *biotic* :steady))
         (d (getf (calling-distance-scenario p *biotic*) :diffusivity))
         (lengths (sort (mapcar (lambda (e) (point-distance (nth (car e) positions) (nth (cdr e) positions)))
                                edges)
                        #'<))
         (median (if lengths (nth (floor (length lengths) 2) lengths) 0d0)))
    (format t "      q ~6,3f um/s (theta ~5,3f): ~d links, median length ~,0f um; Pe median ~,3f (max ~,3f); Da median ~,2f~%"
            q theta (length edges) median (edge-peclet q d median)
            (if lengths (edge-peclet q d (car (last lengths))) 0d0)
            (edge-damkohler p theta median :degrader-density 1d9))))
(finish-output)
