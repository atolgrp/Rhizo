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

(defun species-mass (world species-name)
  "Return the total amount of SPECIES-NAME in WORLD: sum of concentration
times cell area."
  (loop for (area nil value) in (species-cells world species-name)
        sum (* area value)))

(defun species-mean-squared-distance (world species-name centre)
  "Return the mass-weighted mean squared distance of SPECIES-NAME in WORLD
from the point CENTRE, each cell's mass taken to sit at its CELL-POSITION."
  (let ((cells (species-cells world species-name)))
    (/ (loop for (area position value) in cells
             sum (* area value (expt (point-distance position centre) 2)))
       (loop for (area nil value) in cells sum (* area value)))))

(defun relative-l2-error (world species-name field)
  "Return the area-weighted relative L2 error of SPECIES-NAME in WORLD
against FIELD, a function of a point evaluated at each CELL-POSITION:
sqrt(sum A (c - f)^2 / sum A f^2)."
  (let ((cells (species-cells world species-name)))
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
