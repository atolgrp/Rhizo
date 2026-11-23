;;; Export raw data for the two physical-sanity checks (see
;;; tests/validation-tests.lisp and src/validation.lisp) at a
;;; representative resolution, for charting.
;;;
;;; Run with:
;;;   sbcl --non-interactive --load scripts/validate-physical-behavior.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

(defparameter *width* 300)
(defparameter *height* 300)
(defparameter *n-points* 4000)
(defparameter *seed* 77)
(defparameter *out-dir* (merge-pathnames "../renders/" *load-pathname*))
(ensure-directories-exist *out-dir*)

(let* ((grid (make-voronoi-grid *width* *height* *n-points* :seed *seed*))
       (tess (voronoi-grid-tessellation grid))
       (avg-area (/ (float (* *width* *height*) 1.0d0) (cell-count tess)))
       (centroids (cell-centroids grid))
       (center-id (voronoi-grid-id-at grid (floor *width* 2) (floor *height* 2)))
       (world (make-world tess :seed *seed*))
       (tick-of-arrival (make-hash-table :test 'equal))
       (curve-points (list (cons 0 0))))
  (seed-molecules-at! world (list center-id))
  (setf (gethash center-id tick-of-arrival) 0)
  (setf (first curve-points) (cons 0 (hash-table-count (world-frontier world))))
  (loop until (zerop (hash-table-count (world-molecules world)))
        do (step-front-propagation! world)
           (push (cons (world-tick world) (hash-table-count (world-frontier world))) curve-points)
           (loop for id being the hash-keys of (world-molecules world)
                 unless (gethash id tick-of-arrival)
                   do (setf (gethash id tick-of-arrival) (world-tick world))))
  (let ((curve (nreverse curve-points)))
    (multiple-value-bind (slope r2) (front-growth-linearity curve avg-area)
      (format t "wavefront speed=~,3f px/tick  R^2=~,5f~%" slope r2)
      (write-csv (merge-pathnames "wavefront-growth-curve.csv" *out-dir*)
                 '("tick" "cumulative_cells" "effective_radius_px")
                 (mapcar (lambda (p) (list (car p) (cdr p) (sqrt (/ (* (cdr p) avg-area) pi)))) curve)))
    (multiple-value-bind (r r2 n) (distance-correlation tick-of-arrival (lambda (id) (gethash id centroids)) center-id)
      (format t "distance correlation r=~,4f R^2=~,4f n=~a~%" r r2 n)
      (let* ((source-centroid (gethash center-id centroids))
             (rows (loop for id being the hash-keys of tick-of-arrival using (hash-value tick)
                         unless (equal id center-id)
                           collect (list tick (euclidean-distance source-centroid (gethash id centroids))))))
        (write-csv (merge-pathnames "distance-correlation.csv" *out-dir*)
                   '("graph_hops" "euclidean_distance_px")
                   rows)))))

(format t "~%wrote wavefront-growth-curve.csv and distance-correlation.csv to ~a~%" *out-dir*)
