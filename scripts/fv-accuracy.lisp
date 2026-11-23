;;; Finite-volume diffusion accuracy on raster Voronoi grids.
;;;
;;; Measures the edge-distance bug found on 2026-11-23 and the fix adopted
;;; for Phase 1 (see PLAN.md). A Gaussian of width W/20 diffuses with D = 1
;;; until its width doubles; the result is compared with the analytic
;;; solution. Three geometries are compared on the same seeds:
;;;
;;;   current        - what MAKE-VORONOI-GRID gives today: every edge 1.0
;;;                    apart, face length = staircase pixel count
;;;   ca-grown+fix   - same CA-grown cells, seed-to-seed distance,
;;;                    face length sqrt(nx^2 + ny^2)
;;;   nearest+fix    - cells assigned to the nearest seed (true pixelated
;;;                    Voronoi), same two fixes
;;;
;;; Expected (measured 2026-11-23, 3 meshes per row, ~225 px per cell):
;;;
;;;   cells  current          ca-grown+fix     nearest+fix
;;;   100    92%, D = 10.2    46%, D = 1.02    46%, D = 1.00
;;;   400    92%, D = 10.3    8.5%, D = 1.02   9.1%, D = 1.00
;;;   1600   92%, D = 10.3    2.3%, D = 1.02   2.3%, D = 1.00
;;;
;;; The fixed geometries converge at second order (error /4 per 4x cells);
;;; the current one does not converge. Mass is conserved in every case.
;;;
;;; Run with (about 3 minutes):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/fv-accuracy.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

;;; ---------------------------------------------------------------- geometry

(defun squared-distance (p q)
  (+ (expt (- (car p) (car q)) 2) (expt (- (cdr p) (cdr q)) 2)))

(defun nearest-seed (seeds x y)
  "Index of the seed in SEEDS (a vector of (x . y)) closest to pixel (X, Y)."
  (let ((pixel (cons x y)))
    (reduce (lambda (best i)
              (if (< (squared-distance pixel (aref seeds i))
                     (squared-distance pixel (aref seeds best)))
                  i best))
            (alexandria:iota (length seeds)))))

(defun nearest-seed-ids (width height seeds)
  "A fresh WIDTH x HEIGHT id array: each pixel belongs to its nearest seed."
  (let ((ids (make-array (list width height) :element-type 'fixnum)))
    (dotimes (x width ids)
      (dotimes (y height)
        (setf (aref ids x y) (nearest-seed seeds x y))))))

(defun pixel-pairs (width height)
  "Every 4-connected pixel pair, once: ((x . y) (nx . ny) axis), axis :x or :y."
  (loop for x below width
        nconc (loop for y below height
                    nconc (append (when (< (1+ x) width) (list (list (cons x y) (cons (1+ x) y) :x)))
                                  (when (< (1+ y) height) (list (list (cons x y) (cons x (1+ y)) :y)))))))

(defun boundary-step-counts (ids width height)
  "EQUAL hash-table (a . b), a < b -> (nx . ny): how many boundary steps
between cells a and b cross in x and in y. A straight face drawn in pixels
is a staircase of exactly these steps."
  (reduce (lambda (counts pair)
            (destructuring-bind ((x . y) (nx . ny) axis) pair
              (let ((a (aref ids x y)) (b (aref ids nx ny)))
                (if (= a b)
                    counts
                    (let* ((key (cons (min a b) (max a b)))
                           (old (gethash key counts '(0 . 0))))
                      (setf (gethash key counts)
                            (if (eq axis :x)
                                (cons (1+ (car old)) (cdr old))
                                (cons (car old) (1+ (cdr old)))))
                      counts)))))
          (pixel-pairs width height)
          :initial-value (make-hash-table :test 'equal)))

(defun cell-pixel-areas (ids)
  "EQUAL hash-table cell-id -> pixel count."
  (reduce (lambda (areas id) (incf (gethash id areas 0d0)) areas)
          (make-array (array-total-size ids) :displaced-to ids :element-type 'fixnum)
          :initial-value (make-hash-table :test 'equal)))

(defun face-geometry (mode seeds a b counts)
  "(values distance length) for the face between A and B under MODE."
  (destructuring-bind (nx . ny) counts
    (ecase mode
      (:current (values 1d0 (float (+ nx ny) 1d0)))
      (:fixed (values (sqrt (float (squared-distance (aref seeds a) (aref seeds b)) 1d0))
                      (sqrt (float (+ (* nx nx) (* ny ny)) 1d0)))))))

(defun raster-tessellation (ids width height seeds mode)
  "A HASH-GRAPH-TESSELLATION over the cells of IDS, with face geometry per MODE."
  (let ((adjacency (make-hash-table :test 'equal))
        (lengths (make-hash-table :test 'equal))
        (areas (cell-pixel-areas ids)))
    (maphash (lambda (id area) (declare (ignore area)) (ensure-cell! adjacency id)) areas)
    (maphash (lambda (key counts)
               (multiple-value-bind (distance len) (face-geometry mode seeds (car key) (cdr key) counts)
                 (add-edge! adjacency (car key) (cdr key) distance)
                 (set-edge-length! lengths (car key) (cdr key) len)))
             (boundary-step-counts ids width height))
    (make-hash-graph-tessellation adjacency :cell-areas areas :edge-lengths lengths)))

;;; ---------------------------------------------------------------- experiment

(defun gaussian (centre sigma)
  "Normalised 2D Gaussian of width SIGMA at CENTRE, as a function of a point."
  (lambda (p) (/ (exp (- (/ (squared-distance p centre) (* 2 sigma sigma))))
                 (* 2 pi sigma sigma))))

(defun cell-values (world)
  "List of (cell-id . value) for WORLD's single species."
  (let ((index (mechanistic-world-cell-index world))
        (state (mechanistic-world-state world)))
    (loop for i below (cell-index-count index)
          collect (cons (index->cell-id index i) (aref state i 0)))))

(defun second-moment (tess positions centre values)
  "(mass . mean squared distance from CENTRE) of VALUES over TESS."
  (let ((mass (loop for (id . c) in values sum (* c (cell-area tess id))))
        (r2 (loop for (id . c) in values
                  sum (* c (cell-area tess id) (squared-distance (gethash id positions) centre)))))
    (cons mass (/ r2 mass))))

(defun diffuse-gaussian (tess positions width)
  "Diffuse a Gaussian from sigma W/20 to W/10 with D = 1. Returns a plist:
:l2 relative error against the analytic solution, :d-eff from second-moment
growth, :mass-error."
  (let* ((centre (cons (/ width 2d0) (/ width 2d0)))
         (s0 (/ width 20d0)) (s1 (/ width 10d0))
         (t-end (/ (- (* s1 s1) (* s0 s0)) 2d0))
         (species (list (make-species :name :c :diffusivity 1d0)))
         (n-steps (ceiling t-end (stable-timestep tess species)))
         (initial (gaussian centre s0))
         (world (make-mechanistic-world
                 tess :species species :reactions nil :dt (/ t-end n-steps)
                      :initial-state (lambda (id) (list (cons :c (funcall initial (gethash id positions)))))))
         (before (second-moment tess positions centre (cell-values world))))
    (dotimes (k n-steps) (step-mechanistic! world))
    (let* ((values (cell-values world))
           (after (second-moment tess positions centre values))
           (exact (gaussian centre s1))
           (err2 (loop for (id . c) in values
                       sum (* (cell-area tess id) (expt (- c (funcall exact (gethash id positions))) 2))))
           (norm2 (loop for (id . c) in values
                        sum (* (cell-area tess id) (expt (funcall exact (gethash id positions)) 2)))))
      (list :l2 (sqrt (/ err2 norm2))
            :d-eff (/ (- (cdr after) (cdr before)) (* 4 t-end)) ; <r^2> grows as 4Dt in 2D
            :mass-error (/ (abs (- (car after) (car before))) (car before))))))

(defun seed-positions (seeds)
  "EQUAL hash-table cell-id -> seed point as double floats. Cell values are
taken to sit at their seed, where the two-point flux is consistent."
  (let ((table (make-hash-table :test 'equal)))
    (loop for i from 0 for (x . y) across seeds
          do (setf (gethash i table) (cons (float x 1d0) (float y 1d0))))
    table))

(defun mean-of (key results)
  (/ (reduce #'+ results :key (lambda (r) (getf r key))) (length results)))

(defun report (label results)
  (format t "  ~14a error ~5,1f%   D_eff ~5,2f   mass error ~,1e~%"
          label (* 100 (mean-of :l2 results)) (mean-of :d-eff results) (mean-of :mass-error results)))

(defun compare-geometries (width n-cells mesh-seeds)
  (let ((runs (mapcar
               (lambda (s)
                 (let* ((rs (sb-ext:seed-random-state s))
                        (seeds (random-distinct-points width width n-cells rs))
                        (grown (grow-voronoi-ids width width seeds rs))
                        (nearest (nearest-seed-ids width width seeds))
                        (positions (seed-positions seeds)))
                   (flet ((run (ids mode)
                            (diffuse-gaussian (raster-tessellation ids width width seeds mode) positions width)))
                     (list :current (run grown :current)
                           :grown-fixed (run grown :fixed)
                           :nearest-fixed (run nearest :fixed)))))
               mesh-seeds)))
    (format t "~&~a cells on ~ax~a (~a px/cell), ~a meshes, true D = 1:~%"
            n-cells width width (round (* width width) n-cells) (length mesh-seeds))
    (report "current" (mapcar (lambda (r) (getf r :current)) runs))
    (report "ca-grown+fix" (mapcar (lambda (r) (getf r :grown-fixed)) runs))
    (report "nearest+fix" (mapcar (lambda (r) (getf r :nearest-fixed)) runs))
    (finish-output)))

(mapc (lambda (spec) (compare-geometries (first spec) (second spec) '(1 2 3)))
      '((150 100) (300 400) (600 1600)))
