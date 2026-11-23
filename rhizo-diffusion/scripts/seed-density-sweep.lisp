;;; For a FIXED canvas size, sweep the number of Voronoi seed points and
;;; measure, at each density:
;;;   1. cell-shape regularity (mean isoperimetric quotient) - purely
;;;      geometric, no diffusion involved. Answers "how irregular/blocky
;;;      are the polygons themselves."
;;;   2. front-propagation isotropy (coefficient of variation of arrival
;;;      ticks across the 8 compass directions) - the thing actually
;;;      validated earlier. Answers "does diffusion still behave
;;;      isotropically at this resolution."
;;;
;;; Several replicate grids (different seeds) are run per density value to
;;; separate a real trend from run-to-run noise, using the existing
;;; parallel replicate runner.
;;;
;;; Run with:
;;;   sbcl --non-interactive --load scripts/seed-density-sweep.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

(defparameter *sweep-width* 300)
(defparameter *sweep-height* 300)
(defparameter *n-points-values* '(25 50 100 200 400 800 1600 3200 6400 12800 25600 45000
                                   60000 68000 72000 75000 78000 80000 82000 84000 86000 88000))
(defparameter *replicates-per-density* 6)
(defparameter *out-dir* (merge-pathnames "../renders/" *load-pathname*))
(ensure-directories-exist *out-dir*)

(defun sweep-task (n-points seed)
  (let* ((grid (make-voronoi-grid *sweep-width* *sweep-height* n-points :seed seed))
         (tess (voronoi-grid-tessellation grid))
         (cell-count (cell-count tess))
         (shape-q (mean-isoperimetric-quotient grid))
         (points (compass-measurement-points grid))
         (center-id (voronoi-grid-id-at grid (floor *sweep-width* 2) (floor *sweep-height* 2)))
         (world (make-world tess :seed seed)))
    (seed-molecules-at! world (list center-id))
    (let* ((arrivals (run-front-propagation-with-metrics world points :max-ticks 5000))
           (cv (isotropy-coefficient-of-variation arrivals)))
      (list :n-points n-points :seed seed :cell-count cell-count
            :avg-area (/ (float (* *sweep-width* *sweep-height*) 1.0d0) cell-count)
            :shape-q shape-q :isotropy-cv (or cv 0.0d0)))))

(defun mean (xs) (/ (reduce #'+ xs) (float (length xs) 1.0d0)))
(defun sample-stddev (xs)
  (let ((m (mean xs)))
    (sqrt (/ (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) xs)) (float (max 1 (1- (length xs))) 1.0d0)))))

(format t "~%Sweeping ~a density values x ~a replicates on a ~ax~a grid...~%"
        (length *n-points-values*) *replicates-per-density* *sweep-width* *sweep-height*)

(defparameter *raw-rows* '())
(defparameter *summary-rows* '())

(dolist (n *n-points-values*)
  (multiple-value-bind (results) (run-replicates *replicates-per-density*
                                                  (lambda (seed) (sweep-task n seed))
                                                  :workers 4 :base-seed (+ 100000 n))
    (dolist (r results)
      (push (list (getf r :n-points) (getf r :seed) (getf r :cell-count)
                  (getf r :avg-area) (getf r :shape-q) (getf r :isotropy-cv))
            *raw-rows*))
    (let* ((areas (mapcar (lambda (r) (getf r :avg-area)) results))
           (qs (mapcar (lambda (r) (getf r :shape-q)) results))
           (cvs (mapcar (lambda (r) (getf r :isotropy-cv)) results))
           (avg-diameter (sqrt (mean areas))))
      (push (list n avg-diameter (mean qs) (sample-stddev qs) (mean cvs) (sample-stddev cvs))
            *summary-rows*)
      (format t "n=~6a  avg-diameter=~5,1fpx  shape-Q=~,3f (+-~,3f)  isotropy-CV=~,3f (+-~,3f)~%"
              n avg-diameter (mean qs) (sample-stddev qs) (mean cvs) (sample-stddev cvs)))))

(setf *summary-rows* (sort *summary-rows* #'< :key #'first))
(setf *raw-rows* (sort *raw-rows* #'< :key #'first))

(write-csv (merge-pathnames "seed-density-sweep-raw.csv" *out-dir*)
           '("n_points" "seed" "cell_count" "avg_cell_area" "shape_q" "isotropy_cv")
           *raw-rows*)
(write-csv (merge-pathnames "seed-density-sweep-summary.csv" *out-dir*)
           '("n_points" "avg_cell_diameter_px" "mean_shape_q" "std_shape_q" "mean_isotropy_cv" "std_isotropy_cv")
           *summary-rows*)

(format t "~%wrote ~a and ~a~%"
        (merge-pathnames "seed-density-sweep-raw.csv" *out-dir*)
        (merge-pathnames "seed-density-sweep-summary.csv" *out-dir*))
