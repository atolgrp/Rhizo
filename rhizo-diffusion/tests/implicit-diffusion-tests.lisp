(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun dvec (&rest xs)
  (make-array (length xs) :element-type 'double-float
                          :initial-contents (mapcar (lambda (x) (float x 1d0)) xs)))

(defun spd-3x3-product (v)
  "M v for M = [[4 1 0] [1 3 1] [0 1 2]], symmetric positive definite."
  (dvec (+ (* 4 (aref v 0)) (aref v 1))
        (+ (aref v 0) (* 3 (aref v 1)) (aref v 2))
        (+ (aref v 1) (* 2 (aref v 2)))))

(defun nearest-seed-gaussian-world (scheme dt &key (width 300) (n-cells 400) (seed 3))
  "World on a nearest-seed grid holding a Gaussian of width W/20 at the centre."
  (let ((tess (voronoi-grid-tessellation (make-nearest-seed-grid width width n-cells :seed seed))))
    (make-mechanistic-world
     tess :species (list (make-species :name :c :diffusivity 1d0)) :reactions '()
          :dt dt :diffusion-scheme scheme
          :initial-state (field-initial-state
                          tess (list (cons :c (gaussian-field (cons (/ width 2d0) (/ width 2d0))
                                                              (/ width 20d0))))))))

(test vector-dot-by-hand
  (is (= 32d0 (vector-dot (dvec 1 2 3) (dvec 4 5 6))))
  (is (= 0d0 (vector-dot (dvec) (dvec)))))

(test conjugate-gradient-solves-a-small-spd-system
  ;; b = M (1 2 3) = (6 10 8); CG is exact in at most n = 3 steps
  (dolist (inverse-diagonal (list nil (dvec 1/4 1/3 1/2)))
    (multiple-value-bind (x iterations residual)
        (conjugate-gradient #'spd-3x3-product (dvec 6 10 8) (dvec 0 0 0)
                            :inverse-diagonal inverse-diagonal)
      (is (< (reduce #'max (map 'list (lambda (a b) (abs (- a b))) x (dvec 1 2 3))) 1d-13))
      (is (<= iterations 4))
      (is (<= residual 1d-15)))))

(test conjugate-gradient-leaves-its-inputs-alone
  (let ((b (dvec 6 10 8)) (x0 (dvec 1 1 1)))
    (conjugate-gradient #'spd-3x3-product b x0)
    (is (equalp (dvec 6 10 8) b))
    (is (equalp (dvec 1 1 1) x0))))

(test conjugate-gradient-zero-right-hand-side-gives-zero
  (multiple-value-bind (x iterations) (conjugate-gradient #'spd-3x3-product (dvec 0 0 0) (dvec 5 5 5))
    (is (equalp (dvec 0 0 0) x))
    (is (= 0 iterations))))

(test conjugate-gradient-signals-non-convergence
  (signals error (conjugate-gradient #'spd-3x3-product (dvec 6 10 8) (dvec 0 0 0)
                                     :max-iterations 1)))

(test implicit-matrix-product-and-diagonal-by-hand
  ;; areas 2, conductance 3/4 = 0.75, coefficient 0.5, x = (10, 0):
  ;; (A + 0.5 K) x = (20 + 0.5 * 0.75 * 10, 0 - 0.5 * 0.75 * 10)
  (let* ((tess (two-cell-tessellation 2.0d0 3.0d0 4.0d0))
         (index (make-cell-index tess))
         (stencil (make-flux-stencil tess index))
         (x (make-array 2 :element-type 'double-float :initial-element 0d0)))
    (setf (aref x (cell-id->index index "0")) 10d0)
    (let ((y (implicit-matrix-product stencil 0.5d0 x)))
      (is (= 23.75d0 (aref y (cell-id->index index "0"))))
      (is (= -3.75d0 (aref y (cell-id->index index "1")))))
    (is (equalp (dvec (/ 1 2.375d0) (/ 1 2.375d0)) (implicit-inverse-diagonal stencil 0.5d0)))))

(test implicit-two-cell-matches-backward-euler-exactly
  ;; Backward Euler on the two-compartment problem: the difference shrinks
  ;; by 1 / (1 + 2 k dt) per step, k = D L / (d A) = 1 * 1 / (4 * 1) = 0.25.
  ;; Run through STEP-MECHANISTIC!, so the world's scheme is honoured.
  (let* ((dt 0.5d0) (n 20)
         (world (make-mechanistic-world
                 (two-cell-tessellation 1.0d0 1.0d0 4.0d0)
                 :species (list (make-species :name :c :diffusivity 1d0)) :reactions '()
                 :dt dt :diffusion-scheme :implicit
                 :initial-state (lambda (id) (if (string= id "0") '((:c . 10d0)) '((:c . 0d0)))))))
    (run-mechanistic! world n)
    (let ((c0 (cell-species-value world "0" :c)) (c1 (cell-species-value world "1" :c)))
      (is (< (abs (- (- c0 c1) (/ 10d0 (expt (+ 1 (* 2 0.25d0 dt)) n)))) 1d-12))
      (is (< (abs (- 10d0 (+ c0 c1))) 1d-13)))))

(test implicit-diffusion-conserves-mass-far-beyond-the-explicit-limit
  (let* ((probe (nearest-seed-gaussian-world :implicit 1d0))
         (dt (* 100 (stable-timestep (mechanistic-world-tessellation probe)
                                     (mechanistic-world-species-list probe))))
         (world (nearest-seed-gaussian-world :implicit dt))
         (mass-0 (species-mass world :c)))
    (dotimes (k 20) (implicit-diffusion-substep! world dt))
    (is (< (/ (abs (- (species-mass world :c) mass-0)) mass-0) 1d-13))))

(test implicit-diffusion-obeys-the-maximum-principle-at-any-step
  ;; At 1,000x the explicit limit: no negative values, no new maximum.
  (let* ((probe (nearest-seed-gaussian-world :implicit 1d0))
         (dt (* 1000 (stable-timestep (mechanistic-world-tessellation probe)
                                      (mechanistic-world-species-list probe))))
         (world (nearest-seed-gaussian-world :implicit dt))
         (max-0 (reduce #'max (state-column (mechanistic-world-state world) 0))))
    (dotimes (k 5)
      (implicit-diffusion-substep! world dt)
      (let ((column (state-column (mechanistic-world-state world) 0)))
        (is (>= (reduce #'min column) (* -1d-14 max-0)))
        (is (<= (reduce #'max column) (* (+ 1 1d-14) max-0)))))))

(test implicit-diffusion-leaves-immobile-species-alone
  (let* ((tess (two-cell-tessellation 1.0d0 1.0d0 1.0d0))
         (world (make-mechanistic-world
                 tess :species (list (make-species :name :c :diffusivity 1d0)
                                     (make-species :name :fixed :diffusivity 0d0))
                      :reactions '() :dt 10d0 :diffusion-scheme :implicit
                      :initial-state (lambda (id) (if (string= id "0")
                                                      '((:c . 4d0) (:fixed . 7d0))
                                                      '((:c . 0d0) (:fixed . 1d0)))))))
    (multiple-value-bind (result iterations) (implicit-diffusion-substep! world 10d0)
      (is (eq world result))
      (is (plusp iterations)))
    (is (= 7d0 (cell-species-value world "0" :fixed)))
    (is (= 1d0 (cell-species-value world "1" :fixed)))))

(test make-mechanistic-world-rejects-an-unknown-scheme
  (signals error (make-mechanistic-world (two-cell-tessellation 1d0 1d0 1d0)
                                         :species (list (make-species :name :c :diffusivity 1d0))
                                         :reactions '() :dt 1d0 :diffusion-scheme :crank-nicolson)))

(test implicit-and-explicit-agree-where-both-are-stable
  ;; Phase 1 gate. Same grid, same Gaussian, same dt (the explicit limit);
  ;; both schemes are first order in time with opposite-signed error, so
  ;; their difference is small but not zero. Measured 2026-11-24, 400 cells,
  ;; seed 3: L2 7.70% vs 8.04% against the analytic solution.
  (let* ((probe (nearest-seed-gaussian-world :explicit 1d0))
         (dt (stable-timestep (mechanistic-world-tessellation probe)
                              (mechanistic-world-species-list probe)))
         (explicit (nearest-seed-gaussian-world :explicit dt))
         (implicit (nearest-seed-gaussian-world :implicit dt)))
    (run-mechanistic! explicit 100)
    (run-mechanistic! implicit 100)
    (let* ((areas (flux-stencil-areas (mechanistic-world-stencil explicit)))
           (ce (state-column (mechanistic-world-state explicit) 0))
           (ci (state-column (mechanistic-world-state implicit) 0))
           (difference (sqrt (/ (loop for i below (length areas)
                                      sum (* (aref areas i) (expt (- (aref ce i) (aref ci i)) 2)))
                                (loop for i below (length areas)
                                      sum (* (aref areas i) (expt (aref ce i) 2)))))))
      (is (< difference 0.01d0) "explicit and implicit fields differ by ~,3f%" (* 100 difference)))))

(test implicit-gaussian-accuracy-at-ten-times-the-explicit-step
  ;; The implicit step keeps the spatial accuracy criteria of step 5 at 10x
  ;; the explicit limit. Measured 2026-11-24: 2.43% mean error at 1,600
  ;; cells (explicit at 1x: 2.34%), D_eff 1.003.
  (let ((fine (mapcar (lambda (seed)
                        (gaussian-spreading-error
                         (voronoi-grid-tessellation (make-nearest-seed-grid 600 600 1600 :seed seed))
                         :centre '(300d0 . 300d0) :sigma-start 30d0 :sigma-end 60d0
                         :scheme :implicit :time-step-factor 10))
                      '(1 2 3))))
    (is (< (mean-of :l2 fine) 0.03d0) "error at 1,600 cells: ~a" (mean-of :l2 fine))
    (is (< (abs (- 1d0 (mean-of :d-eff fine))) 0.01d0) "D_eff: ~a" (mean-of :d-eff fine))
    (dolist (run fine)
      (is (< (getf run :mass-error) 1d-12) "mass error ~a" (getf run :mass-error)))))
