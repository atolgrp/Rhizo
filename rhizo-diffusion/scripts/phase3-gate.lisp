;;; Phase 3 gate report (PLAN.md): signalling range against water content,
;;; its sensitivity to K_d, and the decision point.
;;;
;;; Analytic tables use the bulk form of CALLING-DISTANCE-SCENARIO with the
;;; literature AHL parameters (P. putida IsoF, 3-oxo-C10-HSL) and the
;;; default scenario (a 1000 um^2 induced microcolony; abiotic loss, or
;;; 1e9 degraders per litre of soil). Two simulations check the analytic
;;; ranges against the engine. Every grid is seeded; seeds are printed.
;;;
;;; Run with (about 20 minutes):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase3-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defparameter *thetas* '(0.10d0 0.15d0 0.20d0 0.25d0 0.26d0 0.30d0 0.35d0 0.40d0 0.45d0))
(defparameter *relations* '(:millington-quirk :papendick-campbell :mullins-sommer))
(defparameter *scenarios* '(("abiotic" ()) ("1e9 degraders/L" (:degrader-density 1d9))))

(format t "~&== 0. Relative diffusivity xi(theta) = D_s / D_0 (soil porosity ~,3f, flagged) ==~%"
        (parameter-value (soil-parameters) :porosity))
(format t "  theta  ~{~20a~}~%" *relations*)
(dolist (theta *thetas*)
  (format t "  ~5,2f  ~{~20,5f~}~%" theta
          (mapcar (lambda (r) (funcall (relative-diffusivity-model r (soil-parameters)) theta))
                  *relations*)))
(format t "  Mullins-Sommer threshold (sandy clay loam, Chou et al. 2012): ~,4f~%"
        (mullins-sommer-threshold (parameter-value (soil-parameters) :mullins-sommer-delta)
                                  (parameter-value (soil-parameters) :mullins-sommer-mu)))
(finish-output)

(format t "~&== 1. Calling distance (um) against water content ==~%")
(dolist (scenario *scenarios*)
  (dolist (relation *relations*)
    (let ((rows (range-versus-water (ahl-literature-parameters) *thetas*
                                    :scenario (list* :relative-diffusivity relation
                                                     (second scenario)))))
      (format t "  ~a, ~a~%    theta   ~{~8,2f~}~%    steady  ~{~8,1f~}~%    day one ~{~8,1f~}~%"
              (first scenario) relation *thetas*
              (mapcar (lambda (r) (getf r :steady)) rows)
              (mapcar (lambda (r) (getf r :transient)) rows)))))
(finish-output)

(format t "~&== 2. Sensitivity to K_d: day-one range (um) for K_d x 0.1, 1, 10 (Millington-Quirk) ==~%")
(let ((kd (parameter-value (ahl-literature-parameters) :kd)))
  (dolist (scenario *scenarios*)
    (format t "  ~a~%    theta      ~{~8,2f~}~%" (first scenario) *thetas*)
    (dolist (factor '(0.1d0 1d0 10d0))
      (let ((rows (range-versus-water (with-parameter-values (ahl-literature-parameters)
                                        (list (cons :kd (* factor kd))))
                                      *thetas*
                                      :scenario (list* :relative-diffusivity :millington-quirk
                                                       (second scenario)))))
        (format t "    K_d x~4,1f ~{~8,1f~}~%" factor (mapcar (lambda (r) (getf r :transient)) rows))))))
(finish-output)

(format t "~&== 3. Decision point: elasticities at the transition theta* (range = half its theta=0.45 value) ==~%")
(format t "  ~18a ~18a ~9a  ~6a  ~8a ~8a ~7a  ~a~%"
        "scenario" "relation" "measure" "theta*" "e_theta" "e_Kd" "ratio" "verdict")
(dolist (scenario *scenarios*)
  (dolist (relation *relations*)
    (dolist (measure '(:steady :transient))
      (let* ((scen (list* :relative-diffusivity relation (second scenario)))
             (theta* (water-transition (ahl-literature-parameters) :scenario scen :measure measure))
             (s (water-sensitivity (ahl-literature-parameters) theta* :scenario scen
                                                                      :measure measure)))
        (format t "  ~18a ~18a ~9a  ~6,3f  ~8,3f ~8,3f ~7,3f  ~a~%"
                (first scenario) relation measure theta*
                (getf s :theta-elasticity) (getf s :kd-elasticity) (getf s :ratio)
                (if (getf s :co-gated) "sorption-and-hydration co-gated" "hydration-gated"))))))
(finish-output)

(format t "~&== 4a. Simulated vs analytic steady range (scaled: lambda 30 um at theta 0.3, MQ, no sorption) ==~%")
(let* ((base (ahl-literature-parameters))
       (xi (relative-diffusivity-model :millington-quirk base))
       (k (/ (* 490d0 (funcall xi 0.3d0)) (* 0.3d0 900d0)))
       (parameters (with-parameter-values base `((:hydrolysis-rate . ,k) (:kd . 0d0)
                                                 (:induced-production . 0d0)))))
  (dolist (n '(1600 6400))
    (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 n :seed 1)))
           (source (nearest-cell tess '(200d0 . 200d0)))
           (s-rate (* (cell-area tess source) 5d14 (parameter-value parameters :basal-production)))
           (threshold (point-source-steady-state s-rate (* 490d0 (funcall xi 0.3d0)) (* 0.3d0 k) 60d0)))
      (dolist (theta '(0.2d0 0.3d0 0.45d0))
        (multiple-value-bind (simulated mass-balance)
            (single-producer-calling-distance
             tess (with-parameter-values parameters `((:water-content . ,theta)))
             :centre '(200d0 . 200d0) :threshold threshold
             :relative-diffusivity :millington-quirk)
          (let ((analytic (steady-calling-distance s-rate (* 490d0 (funcall xi theta)) (* theta k)
                                                   threshold)))
            (format t "  ~4d cells (mesh 1)  theta ~4,2f  simulated ~6,2f um  analytic ~6,2f um  (~5,2f%)  mass ~,1e~%"
                    n theta simulated analytic (* 100 (- (/ simulated analytic) 1)) mass-balance))
          (finish-output))))))

(format t "~&== 4b. Simulated (kinetic sorption) vs analytic (local equilibrium) range at 6 h, MQ, abiotic ==~%")
(format t "  2500 cells on 1.6 mm (spacing 32 um, mesh 1); the 1000 um^2 induced patch in one cell~%")
(defparameter *range-grid* (voronoi-grid-tessellation (make-nearest-seed-grid 1600 1600 2500 :seed 1)))
(defparameter *range-centre* '(800d0 . 800d0))

(defun kinetic-range (parameters &key (dt-factor 0.5d0) (time (* 6 3600d0)) (desorption-factor 1d0))
  "Simulated range of the constant induced patch after TIME with PARAMETERS."
  (let* ((source (nearest-cell *range-grid* *range-centre*))
         (p (with-parameter-values parameters
              `((:basal-production . ,(+ (parameter-value parameters :basal-production)
                                         (parameter-value parameters :induced-production)))
                (:induced-production . 0d0)
                (:desorption-rate . ,(* desorption-factor (parameter-value parameters :desorption-rate)))))))
    (single-producer-calling-distance *range-grid* p :centre *range-centre* :threshold 70d0 :time time
                                      :density (/ (* 1000d0 5d14) (cell-area *range-grid* source))
                                      :dt-factor dt-factor :relative-diffusivity :millington-quirk)))

(defun lea-range (parameters &key (time (* 6 3600d0)))
  (getf (calling-distance-scenario parameters `(:relative-diffusivity :millington-quirk :time ,time))
        :transient))

(let ((base (ahl-literature-parameters)))
  (dolist (case '((1d0 0.5d0) (1d0 0.1d0) (10d0 0.5d0)))
    (destructuring-bind (desorption-factor dt-factor) case
      (multiple-value-bind (sim mass-balance)
          (kinetic-range base :dt-factor dt-factor :desorption-factor desorption-factor)
        (format t "  sorption rate x~4,1f, dt x rate ~4,2f: simulated ~6,1f um, local-equilibrium analytic ~6,1f um (~5,1f%), mass ~,1e~%"
                desorption-factor dt-factor sim (lea-range base) (* 100 (- (/ sim (lea-range base)) 1))
                mass-balance)
        (finish-output))))
  (format t "  (24 h on 10000 cells, measured 2026-11-24: simulated 738.0 um vs analytic 656.3 um, +12.4%)~%"))
(finish-output)

(format t "~&== 4c. Elasticities at theta 0.3, 6 h: kinetic sorption (simulated) vs local equilibrium (analytic) ==~%")
(let* ((base (ahl-literature-parameters))
       (f 1.5d0)
       (kd (parameter-value base :kd))
       (theta (parameter-value base :water-content)))
  (flet ((sim-e (name value)
           (/ (log (/ (kinetic-range (with-parameter-values base (list (cons name (* value f)))))
                      (kinetic-range (with-parameter-values base (list (cons name (/ value f)))))))
              (* 2 (log f))))
         (lea-e (name value)
           (/ (log (/ (lea-range (with-parameter-values base (list (cons name (* value f)))))
                      (lea-range (with-parameter-values base (list (cons name (/ value f)))))))
              (* 2 (log f)))))
    (let ((e-kd-sim (sim-e :kd kd)) (e-theta-sim (sim-e :water-content theta)))
      (format t "  K_d:   simulated ~6,3f  analytic ~6,3f~%" e-kd-sim (lea-e :kd kd))
      (format t "  theta: simulated ~6,3f  analytic ~6,3f~%" e-theta-sim (lea-e :water-content theta))
      (format t "  ratio |e_Kd| / |e_theta|: simulated ~,3f  analytic ~,3f~%"
              (/ (abs e-kd-sim) (abs e-theta-sim))
              (/ (abs (lea-e :kd kd)) (abs (lea-e :water-content theta)))))))
(finish-output)
