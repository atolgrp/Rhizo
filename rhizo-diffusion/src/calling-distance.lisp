(in-package :rhizo-diffusion)

;;; Calling distance and its sensitivity (PLAN.md, Phase 2 gate).
;;;
;;; The calling distance of a producer is the distance out to which
;;; receivers are induced: where the AHL concentration falls to the
;;; receivers' detection threshold. Measured here three ways:
;;;
;;; - at steady state, analytically: the K0 profile of signalling-analysis
;;;   meets the threshold (STEADY-CALLING-DISTANCE). Its elasticities have
;;;   closed forms (STEADY-CALLING-DISTANCE-ELASTICITIES).
;;; - at a finite time after production starts, analytically
;;;   (TRANSIENT-CALLING-DISTANCE), from the point-source solution with
;;;   decay and retardation:
;;;
;;;     C(r, t) = integral over s in [0, t] of
;;;               S / (4 pi D s) exp(-R r^2 / (4 D s) - k s / R) ds,
;;;
;;;   for R dC/dt = D lap C - k C + S delta: sorption at local equilibrium
;;;   (retardation R = 1 + K') and decay of the dissolved signal only. As t
;;;   grows it tends to the K0 steady state, which does not depend on R.
;;; - from a simulated world (CALLING-DISTANCE), on its radial profile.
;;;
;;; With literature parameters steady state can be far away: hydrolysis
;;; alone (half-life 125 h) and K' = 146 give R / k of about three years,
;;; so the transient calling distance is the one that can be compared with
;;; observations (Gantner et al. 2006: 4-78 um on roots).

(defun crossing-distance (concentration threshold scale)
  "Return the distance r at which CONCENTRATION, a function of distance
that decreases from infinity at r -> 0, equals THRESHOLD: bisection in log r
between 1e-12 SCALE and the first SCALE * 2^n at which it is below
THRESHOLD."
  (let ((hi (loop for r = (float scale 1d0) then (* 2 r)
                  until (< (funcall concentration r) threshold)
                  finally (return r))))
    (exp (bisect-root (lambda (u) (- (funcall concentration (exp u)) threshold))
                      (log (* 1d-12 scale)) (log hi)))))

(defun steady-calling-distance (source-rate diffusivity decay-rate threshold)
  "Return the distance at which the steady K0 profile of a source of
SOURCE-RATE (POINT-SOURCE-STEADY-STATE) equals THRESHOLD."
  (crossing-distance (lambda (r) (point-source-steady-state source-rate diffusivity decay-rate r))
                     threshold (sqrt (/ diffusivity decay-rate))))

(defun steady-calling-distance-elasticities (source-rate diffusivity decay-rate threshold)
  "Return the elasticities d ln r_c / d ln p of the STEADY-CALLING-DISTANCE
r_c, as a plist over :SOURCE-RATE, :THRESHOLD, :DIFFUSIVITY and
:DECAY-RATE. With u = r_c / lambda and g = u K1(u) / K0(u): 1/g, -1/g,
1/2 - 1/g and -1/2. Where r_c is much shorter than lambda, g is small and
the calling distance is very sensitive to source strength and threshold
\(the profile is nearly logarithmic there)."
  (let* ((u (/ (steady-calling-distance source-rate diffusivity decay-rate threshold)
               (sqrt (/ diffusivity decay-rate))))
         (g (/ (* u (bessel-k1 u)) (bessel-k0 u))))
    (list :source-rate (/ g) :threshold (- (/ g))
          :diffusivity (- 0.5d0 (/ g)) :decay-rate -0.5d0)))

(defun point-source-transient (source-rate diffusivity decay-rate retardation r time)
  "Return the concentration at distance R > 0 and TIME after a source of
SOURCE-RATE (amount per unit depth per unit time) switches on in the plane,
with DIFFUSIVITY, first-order DECAY-RATE of the dissolved signal and
RETARDATION (see the file header). The integral is taken over ln s by
composite Simpson's rule with 4000 intervals, from where the integrand is
below exp(-50) of its scale; zero if that is after TIME."
  (let* ((a (/ (* retardation r r) (* 4 diffusivity)))
         (b (/ decay-rate retardation))
         (u-lo (log (/ a 50d0)))
         (u-hi (log time))
         (n 4000))
    (if (>= u-lo u-hi)
        0d0
        (let ((h (/ (- u-hi u-lo) n)))
          (flet ((g (u) (exp (- (- (* a (exp (- u)))) (* b (exp u))))))
            (* (/ source-rate (* 4 pi diffusivity))
               (/ h 3d0)
               (+ (g u-lo) (g u-hi)
                  (loop for i from 1 below n
                        sum (* (if (oddp i) 4 2) (g (+ u-lo (* i h))))))))))))

(defun transient-calling-distance (source-rate diffusivity decay-rate retardation threshold time)
  "Return the distance at which POINT-SOURCE-TRANSIENT equals THRESHOLD at
TIME."
  (crossing-distance (lambda (r) (point-source-transient source-rate diffusivity decay-rate
                                                         retardation r time))
                     threshold (sqrt (/ (* 4 diffusivity time) retardation))))

(defun log-elasticity (fn x &key (factor 1.1d0))
  "Return the central-difference elasticity d ln FN / d ln x at X,
ln(FN(X FACTOR) / FN(X / FACTOR)) / (2 ln FACTOR); NIL when X is zero."
  (unless (zerop x)
    (/ (log (/ (funcall fn (* x factor)) (funcall fn (/ x factor))))
       (* 2 (log factor)))))

(defun radial-profile (world species origin bin-width)
  "Return the radial profile of SPECIES in WORLD about the point ORIGIN: a
list of (r . c), one per non-empty annulus [i w, (i+1) w) of width
BIN-WIDTH w, with r its centre and c the area-weighted mean concentration of
the cells whose CELL-POSITION lies in it, in increasing r."
  (let ((bins (make-hash-table)))
    (loop for (area position value) in (species-cells world species)
          for i = (floor (point-distance position origin) bin-width)
          do (let ((bin (or (gethash i bins) (setf (gethash i bins) (list 0d0 0d0)))))
               (incf (first bin) area)
               (incf (second bin) (* area value))))
    (sort (loop for i being the hash-keys of bins using (hash-value (area mass))
                collect (cons (* (+ i 0.5d0) bin-width) (/ mass area)))
          #'< :key #'car)))

(defun profile-crossing (profile threshold)
  "Return the r at which PROFILE, a list of (r . c) in increasing r, first
falls below THRESHOLD, interpolating ln c linearly between the two points
either side (c linearly if the outer one is not positive); NIL if it never
does."
  (loop for ((r0 . c0) (r1 . c1)) on profile
        while r1
        when (and (>= c0 threshold) (< c1 threshold))
          return (if (plusp c1)
                     (+ r0 (* (- r1 r0) (/ (log (/ c0 threshold)) (log (/ c0 c1)))))
                     (+ r0 (* (- r1 r0) (/ (- c0 threshold) (- c0 c1)))))))

(defun calling-distance (world species origin threshold &key bin-width)
  "Return the calling distance measured in WORLD: where the RADIAL-PROFILE
of SPECIES about ORIGIN crosses THRESHOLD (PROFILE-CROSSING). BIN-WIDTH
defaults to the mean cell spacing, sqrt(total area / cells). NIL if the
profile never falls below THRESHOLD."
  (let ((tess (mechanistic-world-tessellation world)))
    (profile-crossing
     (radial-profile world species origin
                     (or bin-width
                         (sqrt (/ (reduce #'+ (cell-ids tess) :key (lambda (id) (cell-area tess id)))
                                  (cell-count tess)))))
     threshold)))

(defun default-calling-scenario ()
  "Return the default scenario plist of CALLING-DISTANCE-SCENARIO - choices
for the simulated experiment, not literature values:
  :PATCH-AREA 1000 um^2 - a producer microcolony; at :PRODUCER-DENSITY it
      holds 500 cells per um of slab depth. Much smaller and, with the
      literature rates, AHL never reaches the 70 nM threshold beyond the
      patch's own edge (100 um^2 gives 0.7 um).
  :PRODUCER-DENSITY 5e14 cells/L - cells of 1 um^3 filling half the
      volume, as in Hense et al. 2012.
  :DEGRADER-DENSITY 0 - abiotic loss only.
  :TIME 86400 s - one day after production starts.
  :INDUCED T - the patch produces at the induced rate."
  (list :patch-area 1000d0 :producer-density 5d14 :degrader-density 0d0 :time 86400d0
        :induced t))

(defun scenario-value (scenario name)
  "Return NAME's value in the plist SCENARIO, or its DEFAULT-CALLING-SCENARIO
value."
  (getf scenario name (getf (default-calling-scenario) name)))

(defun calling-distance-scenario (parameters &optional scenario)
  "Return the analytic calling distances of a producer patch with
PARAMETERS in SCENARIO (a plist over the keys of DEFAULT-CALLING-SCENARIO,
plus :RELATIVE-DIFFUSIVITY, default :MILLINGTON-QUIRK; missing keys take
their defaults), as a plist :STEADY, :TRANSIENT (at :TIME), and the inputs
:SOURCE-RATE, :DIFFUSIVITY, :DECAY-RATE, :CAPACITY, :RETARDATION,
:DECAY-LENGTH, :THRESHOLD. All in the bulk form of signalling.lisp, with
theta = :WATER-CONTENT:
  S = patch area x producer density (per soil volume) x (alpha + beta when
      :INDUCED, else alpha);
  D_s = D_0 xi(theta);  k = theta x hydrolysis + kcat rho_d / Km;
  capacity theta + rho_b K_d (retardation 1 + K').
Receivers sense at :RECEIVER-THRESHOLD (pore water). Assumed:
biodegradation by :DEGRADER-DENSITY degraders is linear (AHL far below Km =
1.5 uM); sorption is at local equilibrium; the patch's production is fixed
\(no feedback from its own AHL); receivers do not take up AHL; the slab is
2D (see signalling.lisp), so ranges exceed those around a point source in
3D. Where D_s is zero (no connected water) both ranges are zero."
  (flet ((p (name) (parameter-value parameters name))
         (v (name) (scenario-value scenario name)))
    (let* ((theta (p :water-content))
           (s (* (v :patch-area) (v :producer-density)
                 (+ (p :basal-production) (if (v :induced) (p :induced-production) 0d0))))
           (d (* (p :diffusivity)
                 (funcall (relative-diffusivity-model (getf scenario :relative-diffusivity
                                                            :millington-quirk)
                                                      parameters)
                          theta)))
           (k (+ (* theta (p :hydrolysis-rate))
                 (/ (* (p :biodegradation-kcat) (v :degrader-density)) (p :biodegradation-km))))
           (capacity (+ theta (sorption-coefficient parameters)))
           (threshold (p :receiver-threshold)))
      (if (zerop d)
          ;; no connected water (e.g. below the Mullins-Sommer threshold)
          (list :steady 0d0 :transient 0d0
                :source-rate s :diffusivity d :decay-rate k :capacity capacity
                :retardation (/ capacity theta) :decay-length 0d0 :threshold threshold)
          (list :steady (steady-calling-distance s d k threshold)
                :transient (transient-calling-distance s d k capacity threshold (v :time))
                :source-rate s :diffusivity d :decay-rate k :capacity capacity
                :retardation (/ capacity theta)
                :decay-length (sqrt (/ d k)) :threshold threshold)))))

(defun calling-distance-sensitivity (parameters
                                     &key scenario (factor 1.1d0)
                                       (parameter-names '(:basal-production :induced-production
                                                          :diffusivity :hydrolysis-rate
                                                          :biodegradation-kcat :biodegradation-km
                                                          :receiver-threshold :kd :bulk-density
                                                          :water-content))
                                       (scenario-names '(:producer-density :patch-area
                                                         :degrader-density :time)))
  "Return the sensitivity table of CALLING-DISTANCE-SCENARIO (with PARAMETERS
and SCENARIO): one row (name base-value steady-elasticity
transient-elasticity) per name in PARAMETER-NAMES (model parameters) and
SCENARIO-NAMES (scenario variables), each elasticity a LOG-ELASTICITY with
FACTOR. A scenario variable at zero (e.g. no degraders) has NIL
elasticities."
  (flet ((row (name base distances-at)
           (list name base
                 (log-elasticity (lambda (x) (getf (funcall distances-at x) :steady)) base
                                 :factor factor)
                 (log-elasticity (lambda (x) (getf (funcall distances-at x) :transient)) base
                                 :factor factor))))
    (append
     (loop for name in parameter-names
           collect (let ((name name))
                     (row name (parameter-value parameters name)
                          (lambda (x)
                            (calling-distance-scenario
                             (with-parameter-values parameters (list (cons name x))) scenario)))))
     (loop for name in scenario-names
           collect (let ((name name))
                     (row name (scenario-value scenario name)
                          (lambda (x) (calling-distance-scenario parameters
                                                                 (list* name x scenario)))))))))
