(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; AHL quorum-sensing species set (PLAN.md, Phase 2).

(defun relative-difference (a b)
  (/ (abs (- a b)) (max (abs a) (abs b))))

(defun ahl-test-parameters (&rest overrides)
  "AHL-LITERATURE-PARAMETERS with OVERRIDES, a plist name value."
  (with-parameter-values (ahl-literature-parameters)
    (loop for (name value) on overrides by #'cddr collect (cons name value))))

(test ahl-literature-parameters-carry-sources-and-units
  (let ((parameters (ahl-literature-parameters)))
    (is (every #'model-parameter-unit parameters))
    (is (equal '(:bulk-density) (flagged-parameters parameters))
        "only the soil bulk density lacks a literature source")
    ;; unit conversions of the published values
    (is (= 490d0 (parameter-value parameters :diffusivity)))           ; 4.9e-10 m^2/s
    (is (< (relative-difference (parameter-value parameters :basal-production)
                                (/ (* 2.3d-19 1d9) 3600d0))
           1d-12))                                                    ; mol/cell/h -> nmol/cell/s
    (is (< (relative-difference (/ (log 2d0) (parameter-value parameters :hydrolysis-rate) 3600d0)
                                125.0d0)
           1d-3))                                                     ; half-life 125 h
    (is (= 70d0 (parameter-value parameters :induction-threshold)))))

(test ahl-partition-coefficient-is-bulk-density-times-kd-over-water-content
  (is (< (relative-difference
          6d0 (ahl-partition-coefficient (ahl-test-parameters :bulk-density 1.5d0 :kd 2d0
                                                              :water-content 0.5d0)))
         1d-15)))

(test ahl-reaction-rates-assemble-every-process
  ;; d[AHL]/dt in one cell = production - hydrolysis - biodegradation - net sorption
  (let* ((parameters (ahl-test-parameters :bulk-density 1d0 :kd 2d0 :water-content 0.5d0
                                          :desorption-rate 1d-3))
         (species (ahl-species parameters))
         (index (make-species-index species))
         (c 35d0) (s 50d0) (producers 1d10) (degraders 1d8)
         (state (make-cell-state index `((:ahl . ,c) (:ahl-sorbed . ,s)
                                         (:producer . ,producers) (:degrader . ,degraders))))
         (rates (net-reaction-derivative state index (ahl-reactions parameters)))
         (p (lambda (name) (parameter-value parameters name)))
         (production (* producers (+ (funcall p :basal-production)
                                     (* (funcall p :induced-production)
                                        (hill-activation c 70d0 2.5d0)))))
         (hydrolysis (* (funcall p :hydrolysis-rate) c))
         (biodegradation (/ (* (funcall p :biodegradation-kcat) degraders c)
                            (+ (funcall p :biodegradation-km) c)))
         (sorption (* 1d-3 (- (* 4d0 c) s)))) ; K' = 1 * 2 / 0.5 = 4
    (is (< (relative-difference (state-value rates index :ahl)
                                (- production hydrolysis biodegradation sorption))
           1d-12))
    (is (< (relative-difference (state-value rates index :ahl-sorbed) sorption) 1d-12))
    (is (< (relative-difference (state-value rates index :ahl-produced) production) 1d-12))
    (is (< (relative-difference (state-value rates index :ahl-lost) (+ hydrolysis biodegradation))
           1d-12))
    (is (zerop (state-value rates index :producer)) "biomass is fixed in Phase 2")))

(test ahl-world-basal-production-matches-closed-form
  ;; One well-mixed cell, no feedback, no sorption, no degraders:
  ;; dC/dt = rho alpha - k C, so C(t) = (rho alpha / k) (1 - exp(-k t)).
  (let* ((k 1d-4) (rho 1d11)
         (parameters (ahl-test-parameters :induced-production 0d0 :kd 0d0 :hydrolysis-rate k))
         (dt 500d0) (steps 40)
         (world (make-ahl-world (single-cell-tessellation 1d0)
                                :parameters parameters :dt dt
                                :initial-state `((:producer . ,rho))))
         (alpha (parameter-value parameters :basal-production)))
    (dotimes (i steps) (step-mechanistic! world))
    (is (< (relative-difference (cell-species-value world "only" :ahl)
                                (* (/ (* rho alpha) k) (- 1d0 (exp (- (* k dt steps))))))
           1d-7))))

(test ahl-sorption-relaxes-to-the-partition-equilibrium
  (let* ((parameters (ahl-test-parameters :hydrolysis-rate 0d0 :basal-production 0d0
                                          :induced-production 0d0))
         (k-prime (ahl-partition-coefficient parameters))
         (world (make-ahl-world (single-cell-tessellation 1d0) :parameters parameters :dt 10d0
                                :initial-state '((:ahl . 100d0)))))
    (dotimes (i 200) (step-mechanistic! world))
    (let ((c (cell-species-value world "only" :ahl))
          (s (cell-species-value world "only" :ahl-sorbed)))
      (is (< (relative-difference (/ s c) k-prime) 1d-6))
      (is (< (relative-difference (+ c s) 100d0) 1d-12) "sorption moves AHL, never makes it"))))

(test ahl-pool-balance-holds-to-round-off
  ;; Production, hydrolysis, biodegradation, sorption and implicit
  ;; diffusion together: AHL + sorbed + lost - produced never changes.
  (let* ((grid (make-nearest-seed-grid 60 60 30 :seed 4))
         (tess (voronoi-grid-tessellation grid))
         (ids (cell-ids tess))
         (producer-cells (subseq ids 0 5))
         (world (make-ahl-world
                 tess :dt 200d0 ; RK4 half-step 100 s against sorption relaxing at 0.01/s
                      :initial-state (lambda (id)
                                       `((:ahl . 10d0)
                                         (:producer . ,(if (member id producer-cells) 1d12 0d0))
                                         (:degrader . ,(if (member id producer-cells) 0d0 1d8))))))
         (before (ahl-pool-balance world)))
    (run-mechanistic! world 50)
    (let ((produced (species-mass world :ahl-produced)))
      (is (plusp produced))
      (is (plusp (species-mass world :ahl-lost)))
      (is (< (abs (- (ahl-pool-balance world) before)) (* 1d-12 (+ produced (abs before))))
          "pool balance moved from ~a to ~a" before (ahl-pool-balance world)))))

(test sorption-relaxation-rate-is-desorption-rate-times-one-plus-k-prime
  ;; K' = 1 * 2 / 0.5 = 4, so k_des (1 + K') = 1e-3 * 5
  (is (< (relative-difference
          5d-3 (sorption-relaxation-rate (ahl-test-parameters :bulk-density 1d0 :kd 2d0
                                                              :water-content 0.5d0
                                                              :desorption-rate 1d-3)))
         1d-15)))
