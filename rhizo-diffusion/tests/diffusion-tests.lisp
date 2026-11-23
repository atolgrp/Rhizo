(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun small-test-world (&key (seed 1))
  (let* ((grid (make-voronoi-grid 25 25 30 :seed seed)))
    (values (make-world (voronoi-grid-tessellation grid) :seed seed) grid)))

(test seed-molecules-at!-marks-active-and-visited
  (multiple-value-bind (world grid) (small-test-world)
    (declare (ignore grid))
    (let ((start (first (cell-ids (world-tessellation world)))))
      (seed-molecules-at! world (list start))
      (is (gethash start (world-molecules world)))
      (is (front-active-p world start)))))

(test front-propagation-is-monotonic
  ;; the set of ever-visited cells must never shrink from one tick to the next
  (multiple-value-bind (world grid) (small-test-world)
    (declare (ignore grid))
    (let ((start (first (cell-ids (world-tessellation world)))))
      (seed-molecules-at! world (list start))
      (let ((previous-count (hash-table-count (world-frontier world))))
        (dotimes (i 40)
          (step-front-propagation! world)
          (let ((count (hash-table-count (world-frontier world))))
            (is (>= count previous-count))
            (setf previous-count count)))))))

(test front-propagation-reaches-every-connected-cell
  (multiple-value-bind (world grid) (small-test-world)
    (declare (ignore grid))
    (let* ((tess (world-tessellation world))
           (all-ids (cell-ids tess))
           (start (first all-ids)))
      (seed-molecules-at! world (list start))
      (dotimes (i (* 2 (length all-ids)))
        (step-front-propagation! world))
      (is (= (length all-ids) (hash-table-count (world-frontier world)))))))

(test front-propagation-self-terminates
  ;; once every cell has been visited, further steps must be safe no-ops
  ;; (an empty MOLECULES front, not an error or infinite growth)
  (multiple-value-bind (world grid) (small-test-world)
    (declare (ignore grid))
    (let* ((all-ids (cell-ids (world-tessellation world)))
           (start (first all-ids)))
      (seed-molecules-at! world (list start))
      (dotimes (i (* 3 (length all-ids)))
        (step-front-propagation! world))
      (is (= 0 (hash-table-count (world-molecules world))))
      ;; calling it again on an already-exhausted front must not error
      (finishes (step-front-propagation! world)))))
