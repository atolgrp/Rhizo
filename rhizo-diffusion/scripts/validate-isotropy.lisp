;;; End-to-end validation: reproduce the isotropy result the whole project
;;; exists to demonstrate (front reaches N/E/S/W and NW/NE/SW/SE at
;;; comparable ticks), across many replicate runs, and confirm the
;;; replicate sweep is actually running in parallel and not just claiming to.
;;;
;;; Run with:
;;;   sbcl --non-interactive --load scripts/validate-isotropy.lisp

(ql:quickload (list "alexandria" "bordeaux-threads" "lparallel" "csv-parser") :silent t)
(load (merge-pathnames "../rhizo-diffusion.asd" *load-pathname*))
(asdf:load-system "rhizo-diffusion" :verbose nil)

(in-package :rhizo-diffusion)

(defparameter *grid-width* 200)
(defparameter *grid-height* 200)
(defparameter *n-points* 4000)
(defparameter *n-replicates* 10)
(defparameter *directions* '(:n :e :s :w :ne :nw :se :sw))

(defun run-one-isotropy-replicate (seed)
  "Build a fresh grid+world from SEED, measure front-propagation arrival
ticks in all 8 compass directions from center. Returns (seed . arrivals)."
  (let* ((grid (make-voronoi-grid *grid-width* *grid-height* *n-points* :seed seed))
         (points (compass-measurement-points grid))
         (world (make-world (voronoi-grid-tessellation grid) :seed (1+ seed)))
         (center-id (voronoi-grid-id-at grid (floor *grid-width* 2) (floor *grid-height* 2))))
    (seed-molecules-at! world (list center-id))
    (cons seed (run-front-propagation-with-metrics world points :max-ticks 1000))))

(format t "~%=== Running ~a replicates SEQUENTIALLY (grid ~ax~a, ~a points) ===~%"
        *n-replicates* *grid-width* *grid-height* *n-points*)
(let ((start (get-internal-real-time)))
  (loop for i below *n-replicates* collect (run-one-isotropy-replicate i))
  (format t "sequential wall time: ~,2f s~%"
          (/ (- (get-internal-real-time) start) internal-time-units-per-second 1.0)))

(format t "~%=== Running ~a replicates in PARALLEL (4 workers) ===~%" *n-replicates*)
(let (results base-seed seeds (start (get-internal-real-time)))
  (multiple-value-setq (results base-seed seeds)
    (run-replicates *n-replicates* #'run-one-isotropy-replicate :workers 4 :base-seed 2026))
  (format t "parallel wall time: ~,2f s  (base-seed ~a, used to derive ~a per-replicate seeds: ~a)~%"
          (/ (- (get-internal-real-time) start) internal-time-units-per-second 1.0)
          base-seed (length seeds) seeds)

  (format t "~%=== Aggregating arrival ticks across replicates ===~%")
  (let ((rows (loop for (seed . arrivals) in results
                     collect (cons seed (mapcar (lambda (d) (cdr (assoc d arrivals))) *directions*)))))
    (write-csv "csv/isotropy-validation.csv"
               (cons "seed" (mapcar #'string-downcase (mapcar #'symbol-name *directions*)))
               rows)
    (dolist (row rows)
      (format t "seed ~4a  ~{~a: ~3a~^  ~}~%" (first row)
              (loop for d in *directions* for v in (rest row) append (list d v))))
    (let* ((per-direction (loop for i from 1 to (length *directions*)
                                 collect (mapcar (lambda (r) (nth i r)) rows)))
           (means (mapcar (lambda (col) (/ (reduce #'+ col) (float (length col)))) per-direction)))
      (format t "~%mean arrival tick per direction: ~{~a=~,1f~^  ~}~%"
              (loop for d in *directions* for m in means append (list d m)))
      (format t "max/min ratio of means (isotropy check, want close to 1.0): ~,2f~%"
              (/ (reduce #'max means) (reduce #'min means))))))
