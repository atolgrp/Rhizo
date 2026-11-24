;;; Phase 4 gate report (PLAN.md): the signalling graph and its percolation.
;;;
;;; 1. Critical water content theta_c, with its spread over random patch
;;;    layouts (seeds) and over parameter ensembles (relative-diffusivity
;;;    relation, scenario, K_d of two soils), against the infinite-domain
;;;    prediction from continuum percolation.
;;; 2. Hysteresis: a tagged simulation (each patch's AHL its own species)
;;;    through two daily wet/dry cycles, with and without quorum feedback.
;;; 3. Damkohler numbers of the edges at theta_c (no flow: Peclet = 0).
;;; 4. Qualitative comparison with Ziadne et al. (2026, Pest Manag Sci,
;;;    doi:10.1002/ps.71085).
;;; Every layout and grid is seeded; seeds are printed.
;;;
;;; Run with (about 25 minutes):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase4-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defparameter *seeds* (alexandria:iota 20 :start 1))
(defparameter *n* 600)
(defparameter *side* 4000d0)

(defun mean-sd (xs)
  (let* ((n (length xs)) (mean (/ (reduce #'+ xs) n)))
    (values mean (sqrt (/ (reduce #'+ xs :key (lambda (x) (expt (- x mean) 2))) (max 1 (1- n)))))))

(format t "~&== 1. Critical water content theta_c: ~d patches on ~,1f mm square, seeds ~a-~a ==~%"
        *n* (/ *side* 1000) (first *seeds*) (car (last *seeds*)))
(format t "  critical connection radius (continuum percolation, mean degree 4.512): ~,1f um~%"
        (critical-connection-radius *n* (* *side* *side*)))
(format t "  ~16a ~9a ~18a ~6a  ~8a ~15a ~a~%" "scenario" "measure" "relation" "K_d" "predicted" "layouts mean+-sd" "no span")
(dolist (scenario '(("1e9 degraders/L" (:degrader-density 1d9) :steady)
                    ("abiotic" () :transient)))
  (dolist (relation '(:millington-quirk :papendick-campbell :mullins-sommer))
    (dolist (kd '(33.61d0 84.31d0))
      (let* ((parameters (with-parameter-values (ahl-literature-parameters) `((:kd . ,kd))))
             (scen (list* :relative-diffusivity relation (second scenario)))
             (measure (third scenario))
             (predicted (handler-case (predicted-percolation-water-content
                                       parameters *n* (* *side* *side*) :scenario scen :measure measure)
                          (error () nil)))
             (failures 0)
             (thetas (loop for seed in *seeds*
                           for theta = (handler-case
                                           (percolation-water-content
                                            parameters :n *n* :width *side* :height *side* :seed seed
                                                       :scenario scen :measure measure :tolerance 1d-4)
                                         (error () (incf failures) nil))
                           when theta collect theta)))
        (multiple-value-bind (mean sd) (if thetas (mean-sd thetas) (values 0 0))
          (format t "  ~16a ~9a ~18a ~6,2f  ~:[   -   ~;~:*~8,3f~] ~6,3f +- ~5,3f   ~d~%"
                  (first scenario) measure relation kd predicted mean sd failures))
        (finish-output)))))

(format t "~&== 2. Hysteresis through wet/dry cycles (tagged simulation) ==~%")
(let* ((side 600d0)
       (tess (voronoi-grid-tessellation (make-nearest-seed-grid 600 600 900 :seed 1)))
       (positions (random-patch-positions 12 side side :seed 3))
       (nodes (mapcar (lambda (p) (nearest-cell tess p)) positions))
       (cell-positions (mapcar (lambda (id) (cell-position tess id)) nodes))
       (base (ahl-literature-parameters))
       (theta (wet-dry-cycle :theta-wet 0.4d0 :theta-dry 0.15d0 :wet-time 21600d0
                             :drying-time 21600d0 :dry-time 21600d0 :wetting-time 21600d0))
       (dt 20d0)
       (steps-per-day (round 86400 dt)))
  (format t "  12 patches (layout seed 3) on a 600 um square, 900 cells (mesh 1); 1e9 degraders/L;~%")
  (format t "  theta 0.40 (6 h) -> 0.15 (6 h ramp, 6 h dry) -> 0.40 (6 h ramp); dt ~,0f s; samples every 30 min~%" dt)
  (dolist (feedback '(nil t))
    (let* ((parameters (if feedback
                           base
                           (with-parameter-values base
                             `((:basal-production . ,(+ (parameter-value base :basal-production)
                                                        (parameter-value base :induced-production)))
                               (:induced-production . 0d0)))))
           (world (make-tagged-ahl-world tess nodes :parameters parameters :dt dt :feedback feedback
                                                    :degrader-density 1d9
                                                    :producer-density 5d14))
           (t0 (get-internal-real-time))
           (balance-0 (tagged-pool-balance world))
           (samples (connectivity-over-cycle world nodes cell-positions
                                             (lambda (time) (funcall theta time))
                                             :millington-quirk (* 2 steps-per-day) 90
                                             (parameter-value base :receiver-threshold) side)))
      (let ((second-cycle (remove-if (lambda (s) (< (getf s :time) 86400d0)) samples)))
        (format t "~&  ~:[constant induced production~;quorum feedback~]: ~,0f s, mass balance ~,1e~%"
                feedback (/ (- (get-internal-real-time) t0) internal-time-units-per-second)
                (/ (abs (- (tagged-pool-balance world) balance-0)) (species-mass world :tags-produced)))
        (format t "    second cycle: hysteresis gap, largest weak component ~,3f; mean degree ~,3f~%"
                (hysteresis-gap second-cycle)
                (hysteresis-gap second-cycle :key :mean-degree))
        (format t "    time(h) theta  weak  strong  degree  spans~%")
        (dolist (s second-cycle)
          (when (zerop (mod (round (getf s :time)) 7200))
            (format t "    ~6,1f  ~4,2f  ~5,3f  ~5,3f  ~6,2f  ~a~%"
                    (/ (getf s :time) 3600) (getf s :theta)
                    (float (getf s :largest-weak-fraction)) (float (getf s :largest-strong-fraction))
                    (float (getf s :mean-degree)) (getf s :spanning))))
        (finish-output)))))

(format t "~&== 3. Damkohler numbers k L^2 / D_s of the edges at theta_c (1e9 degraders/L, MQ, seed 1) ==~%")
(let* ((parameters (ahl-literature-parameters))
       (scen '(:relative-diffusivity :millington-quirk :degrader-density 1d9))
       (theta-c (percolation-water-content parameters :n *n* :width *side* :height *side* :seed 1
                                                      :scenario scen :measure :steady))
       (r-c (range-at-water parameters theta-c scen :steady))
       (positions (random-patch-positions *n* *side* *side* :seed 1))
       (edges (signalling-edges positions (lambda (i j r) (declare (ignore i j)) (< r r-c))))
       (das (sort (mapcar (lambda (e) (edge-damkohler parameters theta-c
                                                      (point-distance (nth (car e) positions)
                                                                      (nth (cdr e) positions))
                                                      :degrader-density 1d9))
                          edges)
                  #'<)))
  (format t "  theta_c ~,3f, calling distance ~,1f um, ~d edges; Da median ~,3f, 90th percentile ~,3f, max ~,3f~%"
          theta-c r-c (length edges) (nth (floor (length das) 2) das)
          (nth (floor (* 9 (length das)) 10) das) (car (last das)))
  (format t "  Peclet number: 0 for every edge (no flow until Phase 5)~%"))
(finish-output)

(format t "~&== 4. Qualitative comparison with Ziadne et al. 2026 (GR24, 0.05-14.2 mm) ==~%")
(let* ((p (ahl-literature-parameters))
       (range (lambda (params scenario theta)
                (range-at-water params theta (list* :relative-diffusivity :millington-quirk scenario)
                                :transient))))
  (format t "  wetter soil lengthens range (theirs: yes, more in coarse soils):~%")
  (format t "    day-one range, abiotic, MQ: theta 0.15 ~,0f um, 0.30 ~,0f um, 0.45 ~,0f um~%"
          (funcall range p '() 0.15d0) (funcall range p '() 0.3d0) (funcall range p '() 0.45d0))
  (format t "  biological degradation shortens range (theirs: sterilisation lengthens it):~%")
  (format t "    theta 0.30: abiotic ~,0f um, 1e9 degraders/L ~,0f um~%"
          (funcall range p '() 0.3d0) (funcall range p '(:degrader-density 1d9) 0.3d0))
  (format t "  stronger sorption (finer texture) shortens range (theirs: clay < loamy sand < sand):~%")
  (format t "    theta 0.30, abiotic: K_d 33.6 L/kg ~,0f um, 84.3 L/kg ~,0f um~%"
          (funcall range p '() 0.3d0)
          (funcall range (with-parameter-values p '((:kd . 84.31d0))) '() 0.3d0))
  (format t "  The model's ranges (0.1-0.9 mm, AHL, 2D slab, local equilibrium) fall inside their 0.05-14.2 mm (GR24, columns).~%"))
(finish-output)
