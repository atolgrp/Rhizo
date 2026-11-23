;;; Render a CA-grown Voronoi mosaic and a front-propagation sequence to
;;; PNG, for visual comparison against the old show-voronoi/
;;; show-front-propagation SDL windows from RhizoCA/diff30.org.
;;;
;;; Run with:
;;;   sbcl --non-interactive --load scripts/render-comparison.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

(defparameter *width* 800)
(defparameter *height* 800)
(defparameter *n-points* 6000)
(defparameter *seed* 1000)
(defparameter *out-dir* (merge-pathnames "../renders/" *load-pathname*))
(defparameter *n-frames* 6)

(ensure-directories-exist *out-dir*)

(format t "~%Building a ~ax~a CA-grown Voronoi grid with ~a points (seed ~a)...~%"
        *width* *height* *n-points* *seed*)
(defparameter *grid* (make-voronoi-grid *width* *height* *n-points* :seed *seed*))

(defparameter *mosaic-path* (merge-pathnames "voronoi-mosaic.png" *out-dir*))
(render-voronoi-grid-png *grid* *mosaic-path*)
(format t "wrote ~a~%" *mosaic-path*)

(defparameter *tess* (voronoi-grid-tessellation *grid*))
(defparameter *center-id* (voronoi-grid-id-at *grid* (floor *width* 2) (floor *height* 2)))
(defparameter *world* (make-world *tess* :seed *seed*))
(seed-molecules-at! *world* (list *center-id*))

;; Run to completion first to find out how many ticks it actually takes,
;; then re-run from scratch saving a frame at evenly-spaced ticks - simpler
;; than guessing a tick count up front.
(defparameter *total-ticks* 0)
(loop until (zerop (hash-table-count (world-molecules *world*)))
      do (step-front-propagation! *world*)
         (incf *total-ticks*))
(format t "front-propagation fully covered the grid in ~a ticks~%" *total-ticks*)

(defparameter *world* (make-world *tess* :seed *seed*))
(seed-molecules-at! *world* (list *center-id*))
(defparameter *frame-ticks*
  (remove-duplicates (loop for i from 0 to *n-frames*
                            collect (round (* i (/ *total-ticks* *n-frames*))))))

(render-front-propagation-png *grid* *world*
                               (merge-pathnames (format nil "front-tick-~3,'0d.png" 0) *out-dir*)
                               :source-ids (list *center-id*))
(format t "wrote front-tick-000.png (tick 0)~%")

(dotimes (i *total-ticks*)
  (step-front-propagation! *world*)
  (when (member (1+ i) *frame-ticks*)
    (let ((path (merge-pathnames (format nil "front-tick-~3,'0d.png" (1+ i)) *out-dir*)))
      (render-front-propagation-png *grid* *world* path :source-ids (list *center-id*))
      (format t "wrote ~a (tick ~a)~%" (file-namestring path) (1+ i)))))

(format t "~%Done. Files are in ~a~%" *out-dir*)
