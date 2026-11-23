(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(test net-reaction-derivative-stoichiometry-bookkeeping
  (let* ((species (list (make-species :name :a :diffusivity 0.0d0)
                         (make-species :name :b :diffusivity 0.0d0)))
         (index (make-species-index species))
         (reactions (list (make-reaction :name :a->b
                                          :stoichiometry '((:a . -1.0d0) (:b . 1.0d0))
                                          :rate-fn (first-order-decay-rate 0.5d0 :a))))
         (state (make-cell-state index '((:a . 2.0d0) (:b . 0.0d0)))))
    (let ((deriv (net-reaction-derivative state index reactions)))
      (is (= -1.0d0 (state-value deriv index :a)))
      (is (= 1.0d0 (state-value deriv index :b))))))

(test net-reaction-derivative-sums-multiple-reactions-per-species
  (let* ((species (list (make-species :name :a :diffusivity 0.0d0)))
         (index (make-species-index species))
         (reactions (list (make-reaction :name :decay-1
                                          :stoichiometry '((:a . -1.0d0))
                                          :rate-fn (first-order-decay-rate 0.5d0 :a))
                           (make-reaction :name :decay-2
                                          :stoichiometry '((:a . -1.0d0))
                                          :rate-fn (first-order-decay-rate 0.25d0 :a))))
         (state (make-cell-state index '((:a . 4.0d0)))))
    ;; two independent decay pathways acting on the same pool: rates add
    (is (= -3.0d0 (state-value (net-reaction-derivative state index reactions) index :a)))))

(test first-order-decay-matches-closed-form-exponential
  ;; RK4-CELL-REACTION! integrating pure first-order decay (no other
  ;; reactions, no diffusion) should track the exact analytic solution
  ;; a(t) = a0 * exp(-k*t) closely - the strongest single piece of evidence
  ;; the reaction integrator itself is numerically correct, independent of
  ;; anything spatial.
  (let* ((k 0.3d0) (a0 10.0d0) (dt 0.01d0) (n-steps 500)
         (species (list (make-species :name :a :diffusivity 0.0d0)))
         (index (make-species-index species))
         (reactions (list (make-reaction :name :decay
                                          :stoichiometry '((:a . -1.0d0))
                                          :rate-fn (first-order-decay-rate k :a))))
         (state (make-array (list 1 1) :element-type 'double-float :initial-element a0)))
    (dotimes (i n-steps)
      (rk4-cell-reaction! state 0 1 index reactions dt))
    (let* ((total-time (* dt n-steps))
           (expected (* a0 (exp (- (* k total-time))))))
      (is (< (abs (- (aref state 0 0) expected)) 1d-6)
          "expected ~a after t=~a, got ~a" expected total-time (aref state 0 0)))))

(test linear-reversible-rate-drives-toward-equilibrium
  (let* ((species (list (make-species :name :mobile :diffusivity 0.0d0)
                         (make-species :name :sorbed :diffusivity 0.0d0)))
         (index (make-species-index species))
         (rate-fn (linear-reversible-rate 1.0d0 :mobile :sorbed 2.0d0)))
    ;; equilibrium is sorbed = 2*mobile; mobile=1,sorbed=0 is below it -> net sorption (positive)
    (is (> (funcall rate-fn (make-cell-state index '((:mobile . 1.0d0) (:sorbed . 0.0d0))) index)
           0.0d0))
    ;; mobile=1,sorbed=5 is above equilibrium -> net desorption (negative)
    (is (< (funcall rate-fn (make-cell-state index '((:mobile . 1.0d0) (:sorbed . 5.0d0))) index)
           0.0d0))
    ;; exactly at equilibrium -> zero net rate
    (is (< (abs (funcall rate-fn (make-cell-state index '((:mobile . 1.0d0) (:sorbed . 2.0d0))) index))
           1d-12))))

(test monod-rate-saturates-with-substrate
  (let* ((species (list (make-species :name :biomass :diffusivity 0.0d0)
                         (make-species :name :substrate :diffusivity 0.0d0)))
         (index (make-species-index species))
         (rate-fn (monod-rate 1.0d0 :biomass :substrate 1.0d0)))
    (let ((low (funcall rate-fn (make-cell-state index '((:biomass . 1.0d0) (:substrate . 0.1d0))) index))
          (high (funcall rate-fn (make-cell-state index '((:biomass . 1.0d0) (:substrate . 100.0d0))) index)))
      (is (< low high) "rate should increase with more substrate")
      (is (< high 1.05d0) "rate should approach mu-max=1.0 as substrate -> infinity, not exceed it"))))
