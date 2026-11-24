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

;;; Tagged signal: each patch's AHL as its own species (PLAN.md, Phase 4)

(defun tagged-test-grid ()
  (voronoi-grid-tessellation (make-nearest-seed-grid 120 120 144 :seed 5)))

(defun fast-decay-parameters ()
  "Literature AHL parameters with a decay length of about 30 um under
Millington-Quirk at theta 0.3, no sorption, no induction."
  (let* ((base (ahl-literature-parameters))
         (xi (funcall (relative-diffusivity-model :millington-quirk base) 0.3d0)))
    (with-parameter-values base `((:hydrolysis-rate . ,(/ (* 490d0 xi) (* 0.3d0 900d0)))
                                  (:kd . 0d0) (:induced-production . 0d0)))))

(test a-single-tag-is-the-ahl-field
  ;; one patch, linear production: the tag is exactly the AHL of an AHL world
  (let* ((tess (tagged-test-grid))
         (node (nearest-cell tess '(60d0 . 60d0)))
         (parameters (fast-decay-parameters))
         (tagged (make-tagged-ahl-world tess (list node) :parameters parameters :dt 0.5d0
                                                         :producer-density 1d14))
         (plain (make-ahl-world tess :parameters parameters :dt 0.5d0
                                     :initial-state (lambda (id) (when (equal id node)
                                                                   '((:producer . 1d14)))))))
    (run-mechanistic! tagged 40)
    (run-mechanistic! plain 40)
    (let ((a (species-cells tagged :tag-0)) (b (species-cells plain :ahl)))
      (is (< (reduce #'max (mapcar (lambda (x y) (abs (- (third x) (third y)))) a b))
             (* 1d-12 (reduce #'max (mapcar #'third b))))))))

(test tags-add-up-to-the-joint-field
  ;; two patches: each tag equals that patch's field alone (linearity)
  (let* ((tess (tagged-test-grid))
         (nodes (list (nearest-cell tess '(40d0 . 60d0)) (nearest-cell tess '(80d0 . 60d0))))
         (parameters (fast-decay-parameters))
         (both (make-tagged-ahl-world tess nodes :parameters parameters :dt 0.5d0))
         (second-alone (make-tagged-ahl-world tess (list (second nodes)) :parameters parameters
                                                                        :dt 0.5d0)))
    (run-mechanistic! both 30)
    (run-mechanistic! second-alone 30)
    (is (< (reduce #'max (mapcar (lambda (x y) (abs (- (third x) (third y))))
                                 (species-cells both :tag-1) (species-cells second-alone :tag-0)))
           1d-12))))

(test tagged-edges-read-each-tag-at-each-patch
  (let* ((tess (tagged-test-grid))
         (nodes (list (nearest-cell tess '(40d0 . 60d0)) (nearest-cell tess '(70d0 . 60d0))
                      (nearest-cell tess '(110d0 . 110d0))))
         (world (make-tagged-ahl-world tess nodes :parameters (fast-decay-parameters) :dt 0.5d0)))
    (step-to-steady-state! world :tag-0)
    (let* ((c (lambda (tag node)
                (/ (aref (mechanistic-world-state world)
                         (cell-id->index (mechanistic-world-cell-index world) node)
                         (gethash (tag-name tag) (mechanistic-world-species-index world)))
                   0.3d0)))
           ;; threshold between what patch 0 delivers to patch 1 and to patch 2
           (threshold (sqrt (* (funcall c 0 (second nodes)) (funcall c 0 (third nodes)))))
           (edges (tagged-edges world nodes threshold)))
      (is (member '(0 . 1) edges :test #'equal))
      (is (not (member '(0 . 2) edges :test #'equal))))))

(test tagged-world-keeps-its-pool-balance
  (let* ((tess (tagged-test-grid))
         (nodes (list (nearest-cell tess '(40d0 . 60d0)) (nearest-cell tess '(80d0 . 60d0))))
         (world (make-tagged-ahl-world tess nodes :dt 50d0 :feedback t))
         (before (tagged-pool-balance world)))
    (run-mechanistic! world 40)
    (is (plusp (species-mass world :tags-produced)))
    (is (< (abs (- (tagged-pool-balance world) before))
           (* 1d-12 (species-mass world :tags-produced))))))

;;; Connectivity over wet/dry cycles (PLAN.md, Phase 4)

(test graph-connectivity-summarises-a-graph
  ;; 0 <-> 1 -> 2 across [0, 10]; 3 alone
  (let ((c (graph-connectivity '((0.5d0 . 0d0) (5d0 . 0d0) (9.5d0 . 0d0) (5d0 . 9d0))
                               '((0 . 1) (1 . 0) (1 . 2)) 10d0 1d0)))
    (is (= 3/4 (getf c :largest-weak-fraction)))
    (is (= 1/2 (getf c :largest-strong-fraction)))
    (is (getf c :spanning))
    (is (= 3/4 (getf c :mean-degree)))))

(test hysteresis-gap-compares-the-two-branches
  ;; drying: fraction = theta; wetting: fraction = theta - 0.1
  (let* ((drying (loop for theta from 0.4d0 downto 0.1d0 by 0.05d0
                       collect (list :theta theta :largest-weak-fraction theta)))
         (wetting (loop for theta from 0.1d0 to 0.4d0 by 0.05d0
                        collect (list :theta theta :largest-weak-fraction (- theta 0.1d0)))))
    (is (< (abs (- 0.1d0 (hysteresis-gap (append drying wetting)))) 1d-9))
    (is (< (hysteresis-gap (append drying (reverse drying))) 1d-12) "same path both ways: none")))

(test edge-damkohler-is-decay-over-diffusion
  ;; Da = k L^2 / D_s with k = theta k_h (+ biodegradation), D_s = D_0 xi(theta)
  (let* ((parameters (ahl-literature-parameters))
         (xi (funcall (relative-diffusivity-model :millington-quirk parameters) 0.3d0)))
    (is (< (relative-difference (edge-damkohler parameters 0.3d0 100d0)
                                (/ (* 0.3d0 (parameter-value parameters :hydrolysis-rate) 1d4)
                                   (* 490d0 xi)))
           1d-12))))

(test connectivity-over-cycle-samples-the-schedule
  (let* ((tess (tagged-test-grid))
         (nodes (list (nearest-cell tess '(40d0 . 60d0)) (nearest-cell tess '(80d0 . 60d0))))
         (positions (mapcar (lambda (id) (cell-position tess id)) nodes))
         (world (make-tagged-ahl-world tess nodes :parameters (fast-decay-parameters) :dt 0.5d0))
         (theta (wet-dry-cycle :theta-wet 0.4d0 :theta-dry 0.2d0
                               :wet-time 2d0 :drying-time 4d0 :dry-time 2d0 :wetting-time 4d0))
         (samples (connectivity-over-cycle world nodes positions
                                           (lambda (time) (funcall theta time))
                                           :millington-quirk 24 4 1d-6 120d0)))
    (is (= 6 (length samples)))
    (is (equal '(2d0 4d0 6d0 8d0 10d0 12d0) (mapcar (lambda (s) (getf s :time)) samples)))
    (is (every (lambda (s) (= (getf s :theta) (funcall theta (getf s :time)))) samples))
    (is (every (lambda (s) (<= 0 (getf s :largest-weak-fraction) 1)) samples))))
