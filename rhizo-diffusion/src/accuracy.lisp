(in-package :rhizo-diffusion)

;;; Accuracy of transport against analytic solutions.
;;;
;;; Mass conservation alone does not show that diffusion is right: the
;;; CA-grown geometry conserved mass to round-off while diffusing about 10x
;;; too fast (PLAN.md, Phase 1). These functions measure a run against a
;;; known solution instead. Every field is sampled at CELL-POSITION, the
;;; point where a cell-centred scheme holds its value.
;;;
;;; The reference problem is a Gaussian spreading under pure diffusion in
;;; 2D. An initial Gaussian of width s0 stays Gaussian, with width
;;; s(t)^2 = s0^2 + 2 D t, and its mean squared distance from the centre
;;; grows as <r^2> = 2 s^2, i.e. by 4 D t. Measured: the relative L2 error
;;; against that solution, and an effective diffusivity from the growth of
;;; <r^2>. Assumed: the domain is large enough that the zero-flux boundary
;;; does not matter (the default widths keep the Gaussian 5 widths from the
;;; nearest edge).

(defun gaussian-field (centre sigma)
  "Return the function of a point (x . y) giving the 2D Gaussian of width
SIGMA centred at CENTRE, normalised to unit integral over the plane: the
fundamental solution of the diffusion equation with SIGMA^2 = 2 D t."
  (let ((cx (float (car centre) 1.0d0))
        (cy (float (cdr centre) 1.0d0))
        (two-s2 (* 2.0d0 sigma sigma)))
    (lambda (p)
      (let ((dx (- (car p) cx)) (dy (- (cdr p) cy)))
        (/ (exp (- (/ (+ (* dx dx) (* dy dy)) two-s2))) (* pi two-s2))))))

(defun species-cells (world species-name)
  "Return a list of (area position value) for every cell of WORLD, with
VALUE the concentration of SPECIES-NAME and POSITION the CELL-POSITION."
  (let* ((tess (mechanistic-world-tessellation world))
         (index (mechanistic-world-cell-index world))
         (state (mechanistic-world-state world))
         (s (or (gethash species-name (mechanistic-world-species-index world))
                (error "SPECIES-CELLS: unknown species ~s" species-name))))
    (loop for i below (cell-index-count index)
          for id = (index->cell-id index i)
          collect (list (cell-area tess id) (cell-position tess id) (aref state i s)))))

(defun species-mass (world species)
  "Return the total amount of SPECIES in WORLD - one species name, or a list
of names summed together (e.g. a mobile and a sorbed pool): sum of
concentration times cell area."
  (loop for name in (alexandria:ensure-list species)
        sum (loop for (area nil value) in (species-cells world name)
                  sum (* area value))))

(defun species-mean-squared-distance (world species centre)
  "Return the mass-weighted mean squared distance from the point CENTRE of
SPECIES in WORLD - one species name, or a list of names pooled together -
each cell's mass taken to sit at its CELL-POSITION."
  (let ((cells (mapcan (lambda (name) (species-cells world name))
                       (alexandria:ensure-list species))))
    (/ (loop for (area position value) in cells
             sum (* area value (expt (point-distance position centre) 2)))
       (loop for (area nil value) in cells sum (* area value)))))

(defun relative-l2-error (world species-name field &key (where (constantly t)))
  "Return the area-weighted relative L2 error of SPECIES-NAME in WORLD
against FIELD, a function of a point evaluated at each CELL-POSITION:
sqrt(sum A (c - f)^2 / sum A f^2), over the cells whose CELL-POSITION
satisfies the predicate WHERE (default: every cell)."
  (let ((cells (remove-if-not (lambda (cell) (funcall where (second cell)))
                              (species-cells world species-name))))
    (sqrt (/ (loop for (area position value) in cells
                   sum (* area (expt (- value (funcall field position)) 2)))
             (loop for (area position nil) in cells
                   sum (* area (expt (funcall field position) 2)))))))

(defun observed-order (coarse-error fine-error cell-ratio)
  "Return the observed order of convergence between two runs in 2D, where
the fine mesh has CELL-RATIO times as many cells, so a cell side
SQRT(CELL-RATIO) times smaller: log(coarse/fine) / log(sqrt(cell-ratio)).
Second order gives 2."
  ;; FLOAT first: (SQRT 4) of an integer is a single-float
  (/ (log (/ coarse-error fine-error)) (log (sqrt (float cell-ratio 1.0d0)))))

(defun gaussian-spreading-error (tessellation &key centre sigma-start sigma-end
                                                    (diffusivity 1.0d0)
                                                    (scheme :explicit)
                                                    (time-step-factor 1))
  "Diffuse a Gaussian on TESSELLATION from width SIGMA-START to SIGMA-END
with DIFFUSIVITY, in equal steps of at most TIME-STEP-FACTOR times
STABLE-TIMESTEP. SCHEME is :EXPLICIT (DIFFUSION-SUBSTEP!; a factor above 1
is unstable) or :IMPLICIT (IMPLICIT-DIFFUSION-SUBSTEP!). Returns a plist:
  :L2         relative L2 error against the analytic Gaussian at the end
  :D-EFF      effective diffusivity, from the growth of the mean squared
              distance (4 D t in 2D)
  :MASS-ERROR relative change in total mass (pass/fail: round-off only)
  :STEPS      number of explicit steps taken"
  (let* ((duration (/ (- (* sigma-end sigma-end) (* sigma-start sigma-start))
                      (* 2.0d0 diffusivity)))
         (species (list (make-species :name :c :diffusivity (float diffusivity 1.0d0))))
         (steps (ceiling duration (* time-step-factor (stable-timestep tessellation species))))
         (dt (/ duration steps))
         (world (make-mechanistic-world
                 tessellation :species species :reactions '() :dt dt :diffusion-scheme scheme
                              :initial-state (field-initial-state
                                              tessellation
                                              (list (cons :c (gaussian-field centre sigma-start))))))
         (mass-0 (species-mass world :c))
         (r2-0 (species-mean-squared-distance world :c centre)))
    (call-with-kernel
     (lambda ()
       (dotimes (k steps) (scheme-diffusion-substep! world dt))))
    (list :l2 (relative-l2-error world :c (gaussian-field centre sigma-end))
          :d-eff (/ (- (species-mean-squared-distance world :c centre) r2-0) (* 4.0d0 duration))
          :mass-error (/ (abs (- (species-mass world :c) mass-0)) (abs mass-0))
          :steps steps)))

;;; Isotropy of the diffusion front (PLAN.md, Phase 1 gate). A second
;;; moment cannot detect lattice anisotropy - any 4-fold symmetric spread
;;; has an isotropic second-moment tensor - so the front is timed along
;;; compass directions instead, as in the published isotropy study
;;; (voronoi-ca-isotropy/README.md): a point source, probes at equal
;;; distance along the 4 axial and 4 diagonal directions, and the time at
;;; which each probe's concentration first crosses a threshold.

(defun point-source-concentration (mass diffusivity r time)
  "Return the concentration at distance R and TIME after releasing MASS at
a point in the infinite plane, with DIFFUSIVITY:
MASS / (4 pi D t) * exp(-r^2 / (4 D t))."
  (/ (* mass (exp (- (/ (* r r) (* 4.0d0 diffusivity time)))))
     (* 4.0d0 pi diffusivity time)))

(defun point-source-arrival-time (mass diffusivity r threshold)
  "Return the first time at which POINT-SOURCE-CONCENTRATION at distance R
reaches THRESHOLD, found by bisection to a relative 1e-12 on (0, r^2/4D],
over which the concentration rises to its peak. NIL if THRESHOLD exceeds
the peak value MASS / (pi e r^2)."
  (let ((t-peak (/ (* r r) (* 4.0d0 diffusivity))))
    (when (<= threshold (point-source-concentration mass diffusivity r t-peak))
      (loop with lo = 0.0d0 and hi = t-peak
            repeat 200
            for mid = (/ (+ lo hi) 2.0d0)
            do (if (< (point-source-concentration mass diffusivity r mid) threshold)
                   (setf lo mid)
                   (setf hi mid))
            until (< (- hi lo) (* 1d-12 t-peak))
            finally (return (/ (+ lo hi) 2.0d0))))))

(defparameter *compass-directions*
  '((:e 1 0) (:n 0 1) (:w -1 0) (:s 0 -1)
    (:ne 1 1) (:nw -1 1) (:sw -1 -1) (:se 1 -1))
  "The 4 axial then 4 diagonal compass directions, as (key dx dy).")

(defun compass-points (centre radius)
  "Return an alist (direction . point) of the 8 points at exactly RADIUS
from CENTRE along *COMPASS-DIRECTIONS* (diagonals offset by radius/sqrt 2
on each axis)."
  (loop for (key dx dy) in *compass-directions*
        for scale = (/ radius (sqrt (float (+ (* dx dx) (* dy dy)) 1.0d0)))
        collect (cons key (cons (+ (car centre) (* scale dx)) (+ (cdr centre) (* scale dy))))))

(defun threshold-arrival-times (world species-name probes threshold max-steps)
  "Step WORLD by its own DT with SCHEME-DIFFUSION-SUBSTEP! until every
probe's concentration of SPECIES-NAME has reached THRESHOLD, or MAX-STEPS
pass. PROBES is an alist (key . cell-id). Returns an alist (key . time),
each time linearly interpolated within the step where the threshold was
crossed, or NIL for a probe never reached. Mutates WORLD."
  (let* ((index (mechanistic-world-cell-index world))
         (s (gethash species-name (mechanistic-world-species-index world)))
         (dt (mechanistic-world-dt world))
         (rows (mapcar (lambda (p) (cons (car p) (cell-id->index index (cdr p)))) probes))
         (arrivals (mapcar (lambda (p) (cons (car p) nil)) probes)))
    (flet ((value (row) (aref (mechanistic-world-state world) (cdr row) s)))
      ;; a probe already at the threshold before any step arrives at time 0
      (loop for row in rows
            when (>= (value row) threshold)
              do (setf (cdr (assoc (car row) arrivals)) 0.0d0))
      (call-with-kernel
       (lambda ()
         (loop for step from 1 to max-steps
               for before = (mapcar #'value rows)
               do (scheme-diffusion-substep! world dt)
                  (loop for row in rows
                        for c0 in before
                        for c1 = (value row)
                        for entry = (assoc (car row) arrivals)
                        when (and (null (cdr entry)) (>= c1 threshold))
                          do (setf (cdr entry)
                                   (* dt (+ (1- step) (/ (- threshold c0) (- c1 c0))))))
               until (every #'cdr arrivals)))))
    arrivals))

(defun anisotropy-index (arrivals)
  "Return |mean(axial) - mean(diagonal)| / mean(all) for ARRIVALS, an alist
\(direction . time) over *COMPASS-DIRECTIONS* - the index of the published
isotropy study, as a fraction (multiply by 100 for its percentages)."
  (flet ((mean (keys) (/ (reduce #'+ keys :key (lambda (k) (cdr (assoc k arrivals))))
                         (length keys))))
    (let ((axial (mean '(:e :n :w :s)))
          (diagonal (mean '(:ne :nw :sw :se))))
      (/ (abs (- axial diagonal)) (mean '(:e :n :w :s :ne :nw :sw :se))))))

(defun point-source-isotropy (tessellation &key centre radius (diffusivity 1.0d0)
                                                (scheme :implicit) dt)
  "Release unit mass in the cell of TESSELLATION nearest CENTRE and time the
diffusion front at 8 compass probes at RADIUS from that cell's position
\(each probe is the cell nearest the compass point). The threshold is half
the analytic peak concentration at RADIUS. Returns a plist:
  :RAW         ANISOTROPY-INDEX of the measured arrival times
  :NORMALISED  ANISOTROPY-INDEX of measured / analytic arrival time, the
               analytic time taken at each probe cell's own distance, so
               probe placement noise drops out
  :COV         ISOTROPY-COEFFICIENT-OF-VARIATION of the normalised times
  :MEAN-RATIO  mean of measured / analytic arrival time (accuracy)
  :ARRIVALS    the measured alist (direction . time)
DT defaults to STABLE-TIMESTEP; SCHEME is :EXPLICIT or :IMPLICIT."
  (let* ((source (nearest-cell tessellation centre))
         (origin (cell-position tessellation source))
         (probes (mapcar (lambda (p) (cons (car p) (nearest-cell tessellation (cdr p))))
                         (compass-points origin radius)))
         (species (list (make-species :name :c :diffusivity (float diffusivity 1.0d0))))
         (dt (or dt (stable-timestep tessellation species)))
         (threshold (* 0.5d0 (/ 1.0d0 (* pi (exp 1.0d0) radius radius))))
         (world (make-mechanistic-world
                 tessellation :species species :reactions '() :dt dt :diffusion-scheme scheme
                              :initial-state (lambda (id)
                                               (if (equal id source)
                                                   (list (cons :c (/ 1.0d0 (cell-area tessellation id))))
                                                   '()))))
         (t-peak (/ (* radius radius) (* 4.0d0 diffusivity)))
         (arrivals (threshold-arrival-times world :c probes threshold (ceiling (* 2 t-peak) dt)))
         (ratios (mapcar (lambda (probe arrival)
                           (cons (car probe)
                                 (/ (cdr arrival)
                                    (point-source-arrival-time
                                     1.0d0 diffusivity
                                     (point-distance origin (cell-position tessellation (cdr probe)))
                                     threshold))))
                         probes arrivals)))
    (list :raw (anisotropy-index arrivals)
          :normalised (anisotropy-index ratios)
          :cov (isotropy-coefficient-of-variation ratios)
          :mean-ratio (/ (reduce #'+ ratios :key #'cdr) (length ratios))
          :arrivals arrivals)))

(defun retardation-measurement (tessellation partition-coefficient
                                &key centre sigma-start sigma-end (diffusivity 1.0d0)
                                  (scheme :explicit) (time-step-factor 1))
  "Measure the retardation factor of linear equilibrium sorption on
TESSELLATION. A Gaussian of total mass is split at equilibrium between a
mobile pool (DIFFUSIVITY) and an immobile sorbed pool, sorbed =
PARTITION-COEFFICIENT * mobile, coupled by LINEAR-REVERSIBLE-RATE, and run
with STEP-MECHANISTIC! until the total spreads from SIGMA-START to SIGMA-END
at the predicted rate D / R. PARTITION-COEFFICIENT is the dimensionless
K' = rho_b K_d / theta, so R = 1 + K' is the retardation factor
R = 1 + rho_b K_d / theta (PLAN.md, Phase 1 gate).

Assumed: local equilibrium. The sorption rate constant is set so the
relaxation rate k (1 + K') is 0.5 / dt - fast against diffusion across a
cell, and inside RK4's stability limit. Returns a plist :R-EXPECTED,
:R-MEASURED (D / D_eff, D_eff from the growth of the total's mean squared
distance), :MASS-ERROR and :STEPS."
  (let* ((r-expected (+ 1.0d0 partition-coefficient))
         (species (list (make-species :name :mobile :diffusivity (float diffusivity 1.0d0))
                        (make-species :name :sorbed :diffusivity 0.0d0)))
         (duration (/ (* r-expected (- (* sigma-end sigma-end) (* sigma-start sigma-start)))
                      (* 2.0d0 diffusivity)))
         (steps (ceiling duration (* time-step-factor (stable-timestep tessellation species))))
         (dt (/ duration steps))
         (rate (/ 0.5d0 (* r-expected dt)))
         (initial (gaussian-field centre sigma-start))
         (world (make-mechanistic-world
                 tessellation :species species :dt dt :diffusion-scheme scheme
                 :reactions (list (make-reaction
                                   :name :sorption
                                   :stoichiometry '((:mobile . -1.0d0) (:sorbed . 1.0d0))
                                   :rate-fn (linear-reversible-rate rate :mobile :sorbed
                                                                    partition-coefficient)))
                 :initial-state (lambda (id)
                                  (let ((total (funcall initial (cell-position tessellation id))))
                                    (list (cons :mobile (/ total r-expected))
                                          (cons :sorbed (/ (* partition-coefficient total)
                                                           r-expected)))))))
         (pools '(:mobile :sorbed))
         (mass-0 (species-mass world pools))
         (r2-0 (species-mean-squared-distance world pools centre)))
    (run-mechanistic! world steps)
    (let ((d-eff (/ (- (species-mean-squared-distance world pools centre) r2-0)
                    (* 4.0d0 duration))))
      (list :r-expected r-expected
            :r-measured (/ diffusivity d-eff)
            :mass-error (/ (abs (- (species-mass world pools) mass-0)) mass-0)
            :steps steps))))

(defun grid-gaussian-errors (grid-builder width n-cells mesh-seeds)
  "Return one GAUSSIAN-SPREADING-ERROR plist per seed in MESH-SEEDS, each on
a fresh grid of WIDTH x WIDTH pixels and N-CELLS cells made by
GRID-BUILDER (MAKE-NEAREST-SEED-GRID or MAKE-VORONOI-GRID, called as
\(GRID-BUILDER width width n-cells :seed seed)). The Gaussian is centred on
the square and spreads from width WIDTH/20 to WIDTH/10 with D = 1. The mesh
seed is added to each plist as :SEED."
  (mapcar (lambda (seed)
            (list* :seed seed
                   (gaussian-spreading-error
                    (voronoi-grid-tessellation (funcall grid-builder width width n-cells :seed seed))
                    :centre (cons (/ width 2.0d0) (/ width 2.0d0))
                    :sigma-start (/ width 20.0d0) :sigma-end (/ width 10.0d0))))
          mesh-seeds))
