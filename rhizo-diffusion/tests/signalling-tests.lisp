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
    (is (equal '(:bulk-density :porosity) (flagged-parameters parameters))
        "only the soil bulk density (and the porosity derived from it) lack a source")
    ;; unit conversions of the published values
    (is (= 490d0 (parameter-value parameters :diffusivity)))           ; 4.9e-10 m^2/s
    (is (< (relative-difference (parameter-value parameters :basal-production)
                                (/ (* 2.3d-19 1d9) 3600d0))
           1d-12))                                                    ; mol/cell/h -> nmol/cell/s
    (is (< (relative-difference (/ (log 2d0) (parameter-value parameters :hydrolysis-rate) 3600d0)
                                125.0d0)
           1d-3))                                                     ; half-life 125 h
    (is (= 70d0 (parameter-value parameters :induction-threshold)))))

(test partition-coefficient-is-bulk-density-times-kd-over-water-content
  (is (< (relative-difference
          6d0 (partition-coefficient (ahl-test-parameters :bulk-density 1.5d0 :kd 2d0
                                                              :water-content 0.5d0)))
         1d-15)))

(test ahl-reaction-rates-assemble-every-process
  ;; Pools are amounts per soil volume; chemistry reads the pore-water
  ;; concentration c = m / theta. d[AHL]/dt per soil volume = production -
  ;; hydrolysis (of the dissolved amount) - biodegradation - net sorption.
  (let* ((parameters (ahl-test-parameters :bulk-density 1d0 :kd 2d0 :water-content 0.5d0
                                          :desorption-rate 1d-3))
         (species (ahl-species parameters))
         (index (make-species-index species))
         (theta 0.5d0) (m 17.5d0) (c (/ m theta)) (s 50d0) (producers 1d10) (degraders 1d8)
         (state (make-cell-state index `((:ahl . ,m) (:ahl-sorbed . ,s) (:water . ,theta)
                                         (:producer . ,producers) (:degrader . ,degraders))))
         (rates (net-reaction-derivative state index (ahl-reactions parameters)))
         (p (lambda (name) (parameter-value parameters name)))
         (production (* producers (+ (funcall p :basal-production)
                                     (* (funcall p :induced-production)
                                        (hill-activation c 70d0 2.5d0)))))
         (hydrolysis (* (funcall p :hydrolysis-rate) m))
         (biodegradation (/ (* (funcall p :biodegradation-kcat) degraders c)
                            (+ (funcall p :biodegradation-km) c)))
         (sorption (* 1d-3 (- (* 2d0 c) s)))) ; rho_b K_d = 1 * 2, against the pore water
    (is (< (relative-difference (state-value rates index :ahl)
                                (- production hydrolysis biodegradation sorption))
           1d-12))
    (is (< (relative-difference (state-value rates index :ahl-sorbed) sorption) 1d-12))
    (is (< (relative-difference (state-value rates index :ahl-produced) production) 1d-12))
    (is (< (relative-difference (state-value rates index :ahl-lost) (+ hydrolysis biodegradation))
           1d-12))
    (is (zerop (state-value rates index :producer)) "biomass is fixed in Phase 2")
    (is (zerop (state-value rates index :water)) "water is prescribed, not reacted")))

(test ahl-world-tracks-water
  ;; theta in every cell; dissolved AHL stored in it with D_0 xi(theta)
  (let* ((parameters (ahl-literature-parameters))
         (tess (two-cell-tessellation 1d0 1d0 1d0))
         (theta (lambda (id) (if (string= id "0") 0.2d0 0.4d0))))
    (flet ((world (&rest args) (apply #'make-ahl-world tess :parameters parameters :dt 1d0 args))
           (at (world array id name)
             (aref array (cell-id->index (mechanistic-world-cell-index world) id)
                   (gethash name (mechanistic-world-species-index world)))))
      (let ((w (world)))
        (is (= 0.3d0 (cell-species-value w "0" :water)) "default: the parameters' water content")
        (is (< (relative-difference (at w (mechanistic-world-diffusivities w) "0" :ahl)
                                    (* 490d0 (funcall (relative-diffusivity-model
                                                       :millington-quirk parameters)
                                                      0.3d0)))
               1d-12)
            "default relation: Millington-Quirk"))
      (let ((w (world :water-content theta :relative-diffusivity :no-tortuosity)))
        (is (= 0.4d0 (cell-species-value w "1" :water)))
        (is (= 0.4d0 (at w (mechanistic-world-capacities w) "1" :ahl)))
        (is (< (relative-difference (at w (mechanistic-world-diffusivities w) "0" :ahl) (* 0.2d0 490d0))
               1d-12))))))

(test ahl-world-basal-production-matches-closed-form
  ;; One well-mixed cell, no feedback, no sorption, no degraders:
  ;; dC/dt = rho alpha - k C, so C(t) = (rho alpha / k) (1 - exp(-k t)).
  (let* ((k 1d-4) (rho 1d11)
         (parameters (ahl-test-parameters :induced-production 0d0 :kd 0d0 :hydrolysis-rate k))
         (dt 500d0) (steps 40)
         (world (make-ahl-world (single-cell-tessellation 1d0)
                                :parameters parameters :dt dt
                                :diffusion-scheme :implicit ; Strang: RK4 reactions
                                :initial-state `((:producer . ,rho))))
         (alpha (parameter-value parameters :basal-production)))
    (dotimes (i steps) (step-mechanistic! world))
    (is (< (relative-difference (cell-species-value world "only" :ahl)
                                (* (/ (* rho alpha) k) (- 1d0 (exp (- (* k dt steps))))))
           1d-7))))

(test ahl-sorption-relaxes-to-the-partition-equilibrium
  (let* ((parameters (ahl-test-parameters :hydrolysis-rate 0d0 :basal-production 0d0
                                          :induced-production 0d0))
         (k-prime (partition-coefficient parameters))
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
                      :diffusion-scheme :implicit
                      :initial-state (lambda (id)
                                       `((:ahl . 10d0)
                                         (:producer . ,(if (member id producer-cells) 1d12 0d0))
                                         (:degrader . ,(if (member id producer-cells) 0d0 1d8))))))
         (before (signal-pool-balance world :ahl)))
    (run-mechanistic! world 50)
    (let ((produced (species-mass world :ahl-produced)))
      (is (plusp produced))
      (is (plusp (species-mass world :ahl-lost)))
      (is (< (abs (- (signal-pool-balance world :ahl) before)) (* 1d-12 (+ produced (abs before))))
          "pool balance moved from ~a to ~a" before (signal-pool-balance world :ahl)))))

(test sorption-relaxation-rate-is-desorption-rate-times-one-plus-k-prime
  ;; K' = 1 * 2 / 0.5 = 4, so k_des (1 + K') = 1e-3 * 5
  (is (< (relative-difference
          5d-3 (sorption-relaxation-rate (ahl-test-parameters :bulk-density 1d0 :kd 2d0
                                                              :water-content 0.5d0
                                                              :desorption-rate 1d-3)))
         1d-15)))

;;; Root exudate as a boundary source (PLAN.md, Phase 2)

(defun root-line-tessellation ()
  "Three cells in a row at x = 0, 2, 4, each of area 4: root cell \"r\",
soil cell \"a\" sharing a face of length 2 with it, and soil cell \"b\"
beyond \"a\" (face length 3), not touching the root."
  (let ((adjacency (make-hash-table :test 'equal))
        (areas (make-hash-table :test 'equal))
        (lengths (make-hash-table :test 'equal)))
    (add-edge! adjacency "r" "a" 1d0)
    (add-edge! adjacency "a" "b" 1d0)
    (dolist (id '("r" "a" "b")) (setf (gethash id areas) 4d0))
    (set-edge-length! lengths "r" "a" 2d0)
    (set-edge-length! lengths "a" "b" 3d0)
    (make-fv-tessellation adjacency :cell-areas areas :edge-lengths lengths
                                    :cell-positions (positions-table '("r" "a" "b")
                                                                     (lambda (id)
                                                                       (cons (* 2 (position id '("r" "a" "b")
                                                                                            :test #'string=))
                                                                             0))))))

(defun root-cell-p (id) (string= id "r"))

(test signal-pool-names-follow-the-signal
  (is (eq :exudate-sorbed (signal-pool :exudate :sorbed)))
  (is (eq :ahl-lost (signal-pool :ahl :lost))))

(test root-perimeter-sums-faces-between-root-and-soil
  (is (= 2d0 (root-perimeter (root-line-tessellation) #'root-cell-p))))

(test root-surface-source-spreads-flux-over-adjacent-soil-cells
  ;; flux 5 over the root-soil face of length 2, into soil cell "a" of area 4
  (let ((source (root-surface-source (root-line-tessellation) #'root-cell-p 5d0)))
    (is (= 2.5d0 (funcall source "a")))
    (is (= 0d0 (funcall source "b")) "a cell not touching the root gets nothing")
    (is (= 0d0 (funcall source "r")) "the root itself is not soil")))

(test exudate-literature-parameters-flag-what-has-no-source
  (let ((parameters (exudate-literature-parameters)))
    (is (every #'model-parameter-unit parameters))
    (is (equal '(:exudation-flux :kd :desorption-rate :bulk-density :porosity)
               (flagged-parameters parameters)))
    ;; Stokes-Einstein scaling of the AHL diffusivity by molar mass, 269.34 -> 330.38 g/mol
    (is (< (relative-difference (parameter-value parameters :diffusivity)
                                (* 490d0 (expt (/ 269.34d0 330.38d0) (/ 1d0 3d0))))
           1d-3))
    ;; half-life 1.5 d
    (is (< (relative-difference (/ (log 2d0) (parameter-value parameters :decay-rate)) (* 1.5d0 86400))
           1d-12))))

(test root-exudate-world-releases-flux-times-perimeter
  ;; No decay, no sorption: total exudate produced after time t is
  ;; flux * root perimeter * t, all of it in the soil, none in the root.
  (let* ((parameters (with-parameter-values (exudate-literature-parameters)
                       '((:decay-rate . 0d0) (:kd . 0d0) (:exudation-flux . 5d0))))
         (tess (root-line-tessellation))
         (world (make-root-exudate-world tess #'root-cell-p :parameters parameters :dt 0.01d0)))
    (dotimes (i 100) (step-mechanistic! world))
    (is (< (relative-difference (species-mass world :exudate-produced) (* 5d0 2d0 1d0)) 1d-12))
    (is (< (relative-difference (species-mass world :exudate) (* 5d0 2d0 1d0)) 1d-9))
    (is (= 0d0 (cell-species-value world "r" :exudate)))
    (is (plusp (cell-species-value world "b" :exudate)) "exudate diffuses on into the soil")
    (is (< (abs (signal-pool-balance world :exudate)) 1d-12))))

(test bulk-unit-helpers
  (let ((parameters (ahl-test-parameters :bulk-density 1.5d0 :kd 2d0 :water-content 0.4d0)))
    (is (= 3d0 (sorption-coefficient parameters)))
    (is (equal '(:per :ahl :water) (pore-water :ahl)))
    (is (= 0.4d0 (funcall (water-content-field nil parameters) "any")))
    (is (= 0.25d0 (funcall (water-content-field 0.25d0 parameters) "any")))
    (is (= 0.1d0 (funcall (water-content-field (lambda (id) (declare (ignore id)) 0.1d0) parameters)
                          "any")))))
