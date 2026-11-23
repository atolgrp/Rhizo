;;; Does the seed-density "breakdown" threshold found on a 300x300 canvas
;;; (renders/seed-density-sweep-summary.csv) depend on canvas size? Sweeps
;;; the SAME set of target average cell diameters across several canvas
;;; sizes, so results are comparable by diameter (a size-invariant
;;; quantity) rather than by raw seed count (which is not).
;;;
;;; Kept deliberately lighter than the single-size sweep (fewer
;;; replicates, fewer parallel workers, coarser floor on diameter) - this
;;; machine has ~3.4GB free RAM, and the finest/densest grids are the most
;;; memory-hungry part of the whole exercise.
;;;
;;; Run with:
;;;   sbcl --dynamic-space-size 3072 --non-interactive --load scripts/seed-density-sweep-multisize.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

(defparameter *canvas-sizes* '(150 300 600))
(defparameter *target-diameters* '(40.0d0 20.0d0 10.0d0 5.0d0 2.5d0 1.8d0 1.4d0 1.2d0 1.1d0))
(defparameter *replicates* 4)
(defparameter *workers* 2)
(defparameter *out-dir* (merge-pathnames "../renders/" *load-pathname*))
(ensure-directories-exist *out-dir*)

(defun sweep-task (width height n-points seed)
  (let* ((grid (make-voronoi-grid width height n-points :seed seed))
         (tess (voronoi-grid-tessellation grid))
         (cell-count (cell-count tess))
         (shape-q (mean-isoperimetric-quotient grid))
         (points (compass-measurement-points grid))
         (center-id (voronoi-grid-id-at grid (floor width 2) (floor height 2)))
         (world (make-world tess :seed seed)))
    (seed-molecules-at! world (list center-id))
    (let* ((arrivals (run-front-propagation-with-metrics world points :max-ticks 5000))
           (cv (isotropy-coefficient-of-variation arrivals)))
      (list :canvas width :n-points n-points :seed seed :cell-count cell-count
            :avg-area (/ (float (* width height) 1.0d0) cell-count)
            :shape-q shape-q :isotropy-cv (or cv 0.0d0)))))

(defun mean (xs) (/ (reduce #'+ xs) (float (length xs) 1.0d0)))
(defun sample-stddev (xs)
  (let ((m (mean xs)))
    (sqrt (/ (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) xs)) (float (max 1 (1- (length xs))) 1.0d0)))))

(defparameter *summary-rows* '())

(dolist (size *canvas-sizes*)
  (format t "~%~%=== canvas ~ax~a ===~%" size size)
  (dolist (d *target-diameters*)
    (let* ((area (* size size))
           (n-points (max 4 (round (/ area (* d d))))))
      (multiple-value-bind (results) (run-replicates *replicates*
                                                      (lambda (seed) (sweep-task size size n-points seed))
                                                      :workers *workers* :base-seed (+ (* size 1000000) (round (* d 100))))
        (let* ((areas (mapcar (lambda (r) (getf r :avg-area)) results))
               (qs (mapcar (lambda (r) (getf r :shape-q)) results))
               (cvs (mapcar (lambda (r) (getf r :isotropy-cv)) results))
               (avg-diameter (sqrt (mean areas))))
          (push (list size d n-points avg-diameter (mean qs) (sample-stddev qs) (mean cvs) (sample-stddev cvs))
                *summary-rows*)
          (format t "target-d=~5,2f  n=~7a  actual-d=~5,2fpx  shape-Q=~,3f (+-~,3f)  isotropy-CV=~,3f (+-~,3f)~%"
                  d n-points avg-diameter (mean qs) (sample-stddev qs) (mean cvs) (sample-stddev cvs)))))))

(setf *summary-rows* (sort *summary-rows* #'< :key (lambda (r) (+ (* (first r) 1000000) (- 100 (fourth r))))))

(write-csv (merge-pathnames "seed-density-sweep-multisize.csv" *out-dir*)
           '("canvas_size" "target_diameter_px" "n_points" "actual_avg_diameter_px"
             "mean_shape_q" "std_shape_q" "mean_isotropy_cv" "std_isotropy_cv")
           *summary-rows*)

(format t "~%~%wrote ~a~%" (merge-pathnames "seed-density-sweep-multisize.csv" *out-dir*))
