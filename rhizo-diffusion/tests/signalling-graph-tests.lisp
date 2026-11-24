(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; The signalling graph (PLAN.md, Phase 4).

(test random-patch-positions-are-seeded-and-inside
  (let ((a (random-patch-positions 50 300d0 200d0 :seed 7))
        (b (random-patch-positions 50 300d0 200d0 :seed 7))
        (c (random-patch-positions 50 300d0 200d0 :seed 8)))
    (is (= 50 (length a)))
    (is (equal a b) "same seed, same patches")
    (is (not (equal a c)))
    (is (every (lambda (p) (and (<= 0d0 (car p) 300d0) (<= 0d0 (cdr p) 200d0))) a))))

(test signalling-edges-follow-the-reach-predicate
  ;; reach: within 1.5 of the source, and only from source 0 at distance 2
  (let* ((positions '((0d0 . 0d0) (1d0 . 0d0) (3d0 . 0d0) (2d0 . 0d0)))
         (edges (signalling-edges positions (lambda (i j r)
                                              (or (< r 1.5d0) (and (= i 0) (= j 3)))))))
    (is (equal '((0 . 1) (0 . 3) (1 . 0) (1 . 3) (2 . 3) (3 . 1) (3 . 2))
               (sort (copy-list edges) (lambda (a b) (or (< (car a) (car b))
                                                         (and (= (car a) (car b))
                                                              (< (cdr a) (cdr b)))))))
        "directed: 0 reaches 3 but 3 does not reach 0")))

(test weak-and-strong-components
  ;; 0 <-> 1 -> 2, 3 alone, 4 <-> 5
  (let ((edges '((0 . 1) (1 . 0) (1 . 2) (4 . 5) (5 . 4))))
    (is (equal '((0 1 2) (3) (4 5)) (weak-components 6 edges)))
    (is (equal '((0 1) (2) (3) (4 5)) (strong-components 6 edges)))
    (is (= 3 (largest-component-size (weak-components 6 edges))))
    (is (= 2 (largest-component-size (strong-components 6 edges))))))

(test spanning-component-touches-both-sides
  ;; a chain across x in [0, 10] spans; a chain ending at x = 6 does not
  (let ((positions '((0.5d0 . 1d0) (3d0 . 1d0) (6d0 . 1d0) (9.5d0 . 1d0))))
    (is (spanning-p '(0 1 2 3) positions 0d0 10d0 1d0))
    (is (not (spanning-p '(0 1 2) positions 0d0 10d0 1d0)))))

(test critical-connection-radius-is-the-continuum-percolation-point
  ;; mean degree n pi r^2 = 4 eta_c = 4.5123 (Mertens & Moore 2012)
  (let ((r (critical-connection-radius 100 1d4)))
    (is (< (abs (- (* (/ 100 1d4) pi r r) (* 4 1.12808737d0))) 1d-9))))

;;; Percolation against water content (PLAN.md, Phase 4 gate)

(test percolation-at-range-limits
  (let ((none (percolation-at-range 50 100d0 100d0 0d0 :seed 3))
        (all (percolation-at-range 50 100d0 100d0 500d0 :seed 3)))
    (is (= (/ 1 50) (getf none :largest-weak-fraction)))
    (is (not (getf none :spanning)))
    (is (= 0 (getf none :mean-degree)))
    (is (= 1 (getf all :largest-weak-fraction)))
    (is (= 1 (getf all :largest-strong-fraction)))
    (is (getf all :spanning))
    (is (= 49 (getf all :mean-degree)))))

(defparameter *percolation-scenario*
  '(:relative-diffusivity :millington-quirk :degrader-density 1d9)
  "Equal induced patches, 1e9 degraders per litre: range rises with theta.")

(test predicted-percolation-water-content-meets-the-critical-radius
  (let* ((parameters (ahl-literature-parameters))
         (theta (predicted-percolation-water-content parameters 400 (* 4000d0 4000d0)
                                                     :scenario *percolation-scenario*
                                                     :measure :steady)))
    (is (< 0.1d0 theta 0.45d0))
    (is (< (abs (- (getf (first (range-versus-water parameters (list theta)
                                                    :scenario *percolation-scenario*))
                         :steady)
                   (critical-connection-radius 400 (* 4000d0 4000d0))))
           1d-6))))

(test percolation-water-content-is-where-the-layout-first-spans
  (let* ((parameters (ahl-literature-parameters))
         (args (list :n 300 :width 2800d0 :height 2800d0 :seed 11
                     :scenario *percolation-scenario* :measure :steady))
         (theta-c (apply #'percolation-water-content parameters args)))
    (flet ((spans (theta)
             (getf (percolation-at-range
                    300 2800d0 2800d0
                    (getf (first (range-versus-water parameters (list theta)
                                                     :scenario *percolation-scenario*))
                          :steady)
                    :seed 11)
                   :spanning)))
      (is (< 0.1d0 theta-c 0.45d0))
      (is (spans (+ theta-c 1d-4)))
      (is (not (spans (- theta-c 1d-4)))))))
