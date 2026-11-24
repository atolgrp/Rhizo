(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; IMEX Euler: implicit diffusion with the reaction rate as an explicit
;;; source, not split - its fixed point is exactly the steady state.

(defun two-cell-source-world (scheme dt)
  "Two cells (area 1, face conductance 1) with species :C (D = 1, decay
0.5) and an immobile source :S = 1 in cell \"0\" producing :C at rate [:S].
Steady state: c0 = 1.2, c1 = 0.8 (production 1 = 0.5 (c0 + c1))."
  (make-mechanistic-world
   (two-cell-tessellation 1d0 1d0 1d0)
   :species (list (make-species :name :c :diffusivity 1d0) (make-species :name :s))
   :reactions (list (make-reaction :name :source :stoichiometry '((:c . 1d0))
                                   :rate-fn (mass-action-rate 1d0 :s))
                    (make-reaction :name :decay :stoichiometry '((:c . -1d0))
                                   :rate-fn (first-order-decay-rate 0.5d0 :c)))
   :dt dt :diffusion-scheme scheme
   :initial-state (lambda (id) (if (string= id "0") '((:s . 1d0)) '()))))

(test imex-steady-state-is-exact-at-large-steps
  (let ((world (two-cell-source-world :imex 1d0))) ; k dt = 0.5
    (dotimes (i 200) (step-mechanistic! world))
    (is (< (abs (- (cell-species-value world "0" :c) 1.2d0)) 1d-12))
    (is (< (abs (- (cell-species-value world "1" :c) 0.8d0)) 1d-12))))

(test strang-splitting-steady-state-is-biased-at-large-steps
  ;; Why :IMEX exists: the split schemes' fixed point is off by O(dt).
  (let ((world (two-cell-source-world :implicit 1d0)))
    (dotimes (i 200) (step-mechanistic! world))
    (is (> (abs (- (cell-species-value world "1" :c) 0.8d0)) 0.01d0))))

(test imex-transient-is-forward-euler-in-the-reactions
  ;; One well-mixed cell: c_n = (q / k) (1 - (1 - k dt)^n) exactly.
  (let* ((k 0.5d0) (dt 0.1d0) (n 30)
         (world (make-mechanistic-world
                 (single-cell-tessellation 1d0)
                 :species (list (make-species :name :c :diffusivity 1d0) (make-species :name :s))
                 :reactions (list (make-reaction :name :source :stoichiometry '((:c . 1d0))
                                                 :rate-fn (mass-action-rate 1d0 :s))
                                  (make-reaction :name :decay :stoichiometry '((:c . -1d0))
                                                 :rate-fn (first-order-decay-rate k :c)))
                 :dt dt :diffusion-scheme :imex :initial-state '((:s . 2d0)))))
    (dotimes (i n) (step-mechanistic! world))
    (is (< (abs (- (cell-species-value world "only" :c)
                   (* (/ 2d0 k) (- 1d0 (expt (- 1d0 (* k dt)) n)))))
           1d-12))))

(test imex-without-reactions-is-backward-euler-diffusion
  (let ((imex (nearest-seed-gaussian-world :imex 5d0))
        (implicit (nearest-seed-gaussian-world :implicit 5d0)))
    (dotimes (i 3)
      (step-mechanistic! imex)
      (step-mechanistic! implicit))
    (is (< (reduce #'max (map 'list (lambda (a b) (abs (- a b)))
                              (sb-ext:array-storage-vector (mechanistic-world-state imex))
                              (sb-ext:array-storage-vector (mechanistic-world-state implicit))))
           1d-14))))

(test imex-conserves-the-ahl-pool-balance
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 60 60 30 :seed 4)))
         (producer-cells (subseq (cell-ids tess) 0 5))
         (world (make-ahl-world
                 tess :dt 50d0 :diffusion-scheme :imex ; dt x sorption relaxation (0.01/s) < 1
                      :initial-state (lambda (id)
                                       `((:ahl . 10d0)
                                         (:producer . ,(if (member id producer-cells) 1d12 0d0))
                                         (:degrader . ,(if (member id producer-cells) 0d0 1d8))))))
         (before (signal-pool-balance world :ahl)))
    (run-mechanistic! world 50)
    (let ((produced (species-mass world :ahl-produced)))
      (is (plusp produced))
      (is (< (abs (- (signal-pool-balance world :ahl) before)) (* 1d-12 (+ produced (abs before))))))))

(test snap-underflow-clears-round-off-but-not-real-negatives
  (is (= 0d0 (snap-underflow -1d-12)))
  (is (= -1d-3 (snap-underflow -1d-3)) "a large negative is a symptom to keep visible")
  (is (= 2.5d0 (snap-underflow 2.5d0))))
