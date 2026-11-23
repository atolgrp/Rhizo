(in-package :rhizo-diffusion)

;;; Isotropy measurement: does the front reach the diagonal directions
;;; (NW/NE/SW/SE) at a different tick than the cardinal ones (N/E/S/W)?
;;; That comparison is the entire reason this project exists - it's the
;;; empirical check on whether the alternating-neighborhood CA growth
;;; actually removed the anisotropy a fixed neighborhood would have.
;;;
;;; This is deliberately grid-specific (it only makes sense to ask "the
;;; cell northeast of center" when there is a notion of 2D position), so
;;; it depends on VORONOI-GRID rather than the generic tessellation
;;; protocol - unlike diffusion.lisp, which knows nothing about grids.
;;;
;;; One deliberate simplification vs. the original notebook: the old code
;;; tracked "peak number of active neighbors" at each measurement point,
;;; which conflated front arrival with local density and was fiddly to get
;;; right. Here we track first-arrival tick directly (the tick at which a
;;; measurement point first joins WORLD-FRONTIER) - same scientific
;;; question, simpler to compute and to test.

(defun compass-measurement-points (grid &key (radius-fraction 1/4))
  "Return an alist ((:n . id) (:e . id) (:s . id) (:w . id)
(:ne . id) (:nw . id) (:se . id) (:sw . id)) of the cell ids at
RADIUS-FRACTION of GRID's half-width from its center, in each of the 8
compass directions."
  (let* ((width (voronoi-grid-width grid))
         (height (voronoi-grid-height grid))
         (cx (floor width 2))
         (cy (floor height 2))
         (radius (max 1 (round (* radius-fraction (min width height) 1/2))))
         (diag (max 1 (round (/ radius (sqrt 2))))))
    (flet ((id-at (x y)
             (voronoi-grid-id-at grid
                                  (alexandria:clamp x 0 (1- width))
                                  (alexandria:clamp y 0 (1- height)))))
      (list (cons :n  (id-at cx (- cy radius)))
            (cons :s  (id-at cx (+ cy radius)))
            (cons :e  (id-at (+ cx radius) cy))
            (cons :w  (id-at (- cx radius) cy))
            (cons :ne (id-at (+ cx diag) (- cy diag)))
            (cons :nw (id-at (- cx diag) (- cy diag)))
            (cons :se (id-at (+ cx diag) (+ cy diag)))
            (cons :sw (id-at (- cx diag) (+ cy diag)))))))

(defun run-front-propagation-with-metrics (world measurement-points &key (max-ticks 1000))
  "Step WORLD's front forward until every id in MEASUREMENT-POINTS (an
alist as returned by COMPASS-MEASUREMENT-POINTS) has been reached, or
MAX-TICKS is hit. Returns an alist ((:n . tick) (:e . tick) ...) of the
first tick at which the front reached each measurement point (NIL for any
not reached within MAX-TICKS)."
  (let ((arrival (make-hash-table :test 'eq)))
    (flet ((record-new-arrivals ()
             (loop for (direction . id) in measurement-points
                   when (and (front-active-p world id) (not (gethash direction arrival)))
                     do (setf (gethash direction arrival) (world-tick world)))))
      (record-new-arrivals)
      (loop until (or (= (hash-table-count arrival) (length measurement-points))
                       (>= (world-tick world) max-ticks))
            do (step-front-propagation! world)
               (record-new-arrivals)))
    (loop for (direction . nil) in measurement-points
          collect (cons direction (gethash direction arrival)))))

(defun isotropy-coefficient-of-variation (arrivals)
  "Summarize an arrival-time alist (as returned by
RUN-FRONT-PROPAGATION-WITH-METRICS) as one scalar: the coefficient of
variation (stddev/mean) of the arrival ticks across directions. 0.0 means
every direction arrived at exactly the same tick (perfect isotropy);
larger values mean more directional spread. NIL if fewer than two
directions were reached (not enough data to say anything)."
  (let ((ticks (remove nil (mapcar #'cdr arrivals))))
    (when (>= (length ticks) 2)
      (let* ((n (float (length ticks) 1.0d0))
             (mean (/ (reduce #'+ ticks) n))
             (variance (/ (reduce #'+ (mapcar (lambda (x) (expt (- x mean) 2)) ticks)) n)))
        (if (zerop mean) 0.0d0 (/ (sqrt variance) mean))))))
