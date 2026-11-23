(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun tiny-replicate (seed)
  "A cheap, deterministic-given-SEED replicate: build a small grid+world
from SEED, run a few diffusion steps, return how many cells the front has
reached - enough to tell two replicates apart without doing real work."
  (let* ((grid (make-voronoi-grid 15 15 20 :seed seed))
         (world (make-world (voronoi-grid-tessellation grid) :seed seed))
         (start (first (cell-ids (voronoi-grid-tessellation grid)))))
    (seed-molecules-at! world (list start))
    (dotimes (i 3) (step-front-propagation! world))
    (list :seed seed :visited (hash-table-count (world-frontier world)))))

(test run-replicates-returns-n-results
  (let ((results (run-replicates 6 #'tiny-replicate :workers 4 :base-seed 1)))
    (is (= 6 (length results)))))

(test run-replicates-reproducible-given-same-base-seed
  (multiple-value-bind (results-a base-a seeds-a) (run-replicates 5 #'tiny-replicate :workers 4 :base-seed 123)
    (multiple-value-bind (results-b base-b seeds-b) (run-replicates 5 #'tiny-replicate :workers 4 :base-seed 123)
      (is (= base-a base-b 123))
      (is (equal seeds-a seeds-b))
      (is (equal results-a results-b)))))

(test run-replicates-uses-distinct-seeds-per-replicate
  (multiple-value-bind (results base seeds) (run-replicates 8 #'tiny-replicate :workers 4 :base-seed 7)
    (declare (ignore base))
    (is (= 8 (length (remove-duplicates seeds))))
    (is (= 8 (length (remove-duplicates results :test #'equal))))))
