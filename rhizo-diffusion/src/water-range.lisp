(in-package :rhizo-diffusion)

;;; Signalling range against water content (PLAN.md, Phase 3 gate).
;;;
;;; Water content enters the calling distance three ways (bulk form,
;;; CALLING-DISTANCE-SCENARIO): through the effective diffusivity
;;; D_s = D_0 xi(theta), through the dissolved amount exposed to hydrolysis
;;; (k theta), and through the capacity theta + rho_b K_d (retardation
;;; 1 + rho_b K_d / theta). Sorption cannot move the steady range, only how
;;; soon it is reached, so its influence is on the range attained within a
;;; finite time - one wet spell.
;;;
;;; PLAN.md's decision point compares the range's sensitivity to K_d with
;;; that to theta "near the transition". The transition is taken here as the
;;; water content theta* at which the range has fallen to half its value in
;;; wet soil (WATER-TRANSITION); the sensitivities are elasticities there
;;; (WATER-SENSITIVITY). If |e_Kd| > 2 |e_theta| the narrative becomes
;;; "sorption-and-hydration co-gated".

(defun range-versus-water (parameters thetas &key scenario)
  "Return one plist (:THETA theta :STEADY r :TRANSIENT r) per water content
in THETAS: the analytic calling distances of CALLING-DISTANCE-SCENARIO with
PARAMETERS' :WATER-CONTENT set to theta, in SCENARIO."
  (mapcar (lambda (theta)
            (let ((result (calling-distance-scenario
                           (with-parameter-values parameters (list (cons :water-content theta)))
                           scenario)))
              (list :theta theta :steady (getf result :steady) :transient (getf result :transient))))
          thetas))

(defun water-transition (parameters &key scenario (measure :transient) (theta-wet 0.45d0)
                                      (theta-dry 0.02d0) (fraction 0.5d0))
  "Return the water content theta* at which the calling distance MEASURE
\\(:STEADY or :TRANSIENT, RANGE-VERSUS-WATER) has fallen to FRACTION of its
value at THETA-WET, by bisection between THETA-DRY and THETA-WET. Signals an
error if the range at THETA-DRY is not below that target."
  (flet ((range (theta) (getf (first (range-versus-water parameters (list theta)
                                                         :scenario scenario))
                              measure)))
    (let ((target (* fraction (range theta-wet))))
      (unless (< (range theta-dry) target)
        (error "WATER-TRANSITION: the range at theta = ~a is not below ~a" theta-dry target))
      (bisect-root (lambda (theta) (- (range theta) target)) theta-dry theta-wet))))

(defun water-sensitivity (parameters theta &key scenario (measure :transient) (factor 1.01d0))
  "Return the sensitivity of the calling distance MEASURE to water content
and to K_d at water content THETA, as a plist :THETA, :RANGE, the
elasticities :THETA-ELASTICITY and :KD-ELASTICITY (LOG-ELASTICITY with
FACTOR, small so as not to cross a relation's threshold), :RATIO |e_Kd| /
|e_theta| and :CO-GATED (ratio above 2, PLAN.md's decision threshold)."
  (let ((at-theta (with-parameter-values parameters (list (cons :water-content theta)))))
    (flet ((range-with (name)
             (lambda (x)
               (getf (calling-distance-scenario
                      (with-parameter-values at-theta (list (cons name x))) scenario)
                     measure))))
      (let* ((e-theta (log-elasticity (range-with :water-content) theta :factor factor))
             (e-kd (log-elasticity (range-with :kd) (parameter-value parameters :kd) :factor factor))
             (ratio (/ (abs e-kd) (abs e-theta))))
        (list :theta theta
              :range (getf (calling-distance-scenario at-theta scenario) measure)
              :theta-elasticity e-theta :kd-elasticity e-kd
              :ratio ratio :co-gated (> ratio 2d0))))))

(defun single-producer-calling-distance (tessellation parameters &key centre (density 5d14)
                                                                   threshold time
                                                                   (relative-diffusivity
                                                                    :millington-quirk)
                                                                   (dt-factor 0.5d0))
  "Return the calling distance measured in a simulation: DENSITY producers
per litre of soil in the cell of TESSELLATION nearest CENTRE, an AHL world
with PARAMETERS (water content :WATER-CONTENT, relation
RELATIVE-DIFFUSIVITY, scheme :IMEX), run to steady state - or for TIME
seconds when TIME is given - and measured where the pore-water AHL falls to
THRESHOLD (CALLING-DISTANCE). The step is DT-FACTOR / (hydrolysis rate +
SORPTION-RELAXATION-RATE), shortened to divide TIME evenly. Second value:
the drift of SIGNAL-POOL-BALANCE over the AHL produced - mass balance,
pass/fail for every transport run."
  (let* ((source (nearest-cell tessellation centre))
         (dt (/ dt-factor (+ (parameter-value parameters :hydrolysis-rate)
                             (sorption-relaxation-rate parameters))))
         (steps (and time (ceiling time dt)))
         (world (make-ahl-world tessellation :parameters parameters
                                             :relative-diffusivity relative-diffusivity
                                             :dt (if time (/ time steps) dt)
                                             :initial-state (lambda (id)
                                                              (when (equal id source)
                                                                `((:producer . ,density)))))))
    (let ((balance-0 (signal-pool-balance world :ahl)))
      (if time
          (run-mechanistic! world steps)
          (step-to-steady-state! world :ahl))
      (values (calling-distance world (pore-water :ahl) (cell-position tessellation source) threshold)
              (/ (abs (- (signal-pool-balance world :ahl) balance-0))
                 (species-mass world :ahl-produced))))))
