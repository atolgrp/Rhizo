(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Storage capacity (PLAN.md, Phase 3): a species' state is an amount per
;;; cell volume m, and diffusion moves it down gradients of m / kappa - the
;;; pore-water concentration, when kappa is the water content.

(defun two-cell-capacity-world (scheme dt &key (capacities '(1d0 3d0)) (values '(4d0 0d0)))
  "Two cells (area 1, conductance 1), species :C with D = 1 and capacity
CAPACITIES in cells \"0\" and \"1\", amounts VALUES."
  (make-mechanistic-world
   (two-cell-tessellation 1d0 1d0 1d0)
   :species (list (make-species :name :c :diffusivity 1d0)) :reactions '()
   :dt dt :diffusion-scheme scheme
   :capacity-field (list (cons :c (lambda (id) (if (string= id "0")
                                                   (first capacities)
                                                   (second capacities)))))
   :initial-state (lambda (id) `((:c . ,(if (string= id "0") (first values) (second values)))))))

(test capacity-array-defaults-to-one
  (let ((world (two-cell-world '(1d0 2d0))))
    (is (every (lambda (x) (= x 1d0))
               (sb-ext:array-storage-vector (mechanistic-world-capacities world))))))

(test diffusion-equalises-amount-over-capacity
  ;; total 4 shared so that m0 / 1 = m1 / 3: m0 = 1, m1 = 3
  (dolist (scheme '(:explicit :implicit :imex))
    (let ((world (two-cell-capacity-world scheme (if (eq scheme :explicit) 0.05d0 1d0))))
      (dotimes (i 2000) (step-mechanistic! world))
      (is (< (abs (- (cell-species-value world "0" :c) 1d0)) 1d-9) "~a: m0 ~a" scheme
          (cell-species-value world "0" :c))
      (is (< (abs (- (cell-species-value world "1" :c) 3d0)) 1d-9))
      (is (< (abs (- (species-mass world :c) 4d0)) 1d-12)))))

(test implicit-step-with-capacity-is-backward-euler-exactly
  ;; x = m / kappa; (A kappa + dt K) x' = A m. With kappa = (1, 3), dt = 0.5,
  ;; K = [[1 -1] [-1 1]]: [[1.5 -0.5] [-0.5 3.5]] x' = (4, 0), so
  ;; x' = (4 * 3.5, 4 * 0.5) / (1.5 * 3.5 - 0.25) = (2.8, 0.4); m' = (2.8, 1.2).
  (let ((world (two-cell-capacity-world :implicit 0.5d0)))
    (step-mechanistic! world)
    (is (< (abs (- (cell-species-value world "0" :c) 2.8d0)) 1d-13))
    (is (< (abs (- (cell-species-value world "1" :c) 1.2d0)) 1d-13))))

(test explicit-step-with-capacity-moves-flux-of-amount-over-capacity
  ;; flux into "1" = D g (m0/k0 - m1/k1) = 4; dt = 0.05: m = (3.8, 0.2)
  (let ((world (two-cell-capacity-world :explicit 0.05d0)))
    (diffusion-substep! world 0.05d0)
    (is (< (abs (- (cell-species-value world "0" :c) 3.8d0)) 1d-14))
    (is (< (abs (- (cell-species-value world "1" :c) 0.2d0)) 1d-14))))

(test stable-timestep-scales-with-capacity
  ;; area * kappa / (2 D g): capacities 1 and 3 -> the smaller, 0.5
  (is (< (abs (- (world-stable-timestep (two-cell-capacity-world :explicit 0.01d0)
                                        :safety-factor 1d0)
                 0.5d0))
         1d-14))
  (is (< (abs (- (world-stable-timestep (two-cell-capacity-world :explicit 0.01d0
                                                                 :capacities '(2d0 3d0))
                                        :safety-factor 1d0)
                 1d0))
         1d-14)))

(test set-species-capacity-rejects-non-positive-values
  (let ((world (two-cell-capacity-world :implicit 1d0)))
    (set-species-capacity! world :c (constantly 2d0))
    (is (every (lambda (x) (= x 2d0))
               (sb-ext:array-storage-vector (mechanistic-world-capacities world))))
    (signals error (set-species-capacity! world :c (constantly 0d0)))
    (signals error (set-species-capacity! world :no-such-species (constantly 1d0)))))

(test intensity-column-divides-by-capacity
  (is (equalp (dvec 2 0.5) (intensity-column (dvec 4 1.5) (dvec 2 3))))
  (let ((column (dvec 1 2)))
    (is (eq column (intensity-column column nil)) "no capacity: the column itself")))

(test capacity-array-validates-its-field
  (let* ((tess (two-cell-tessellation 1d0 1d0 1d0))
         (index (make-cell-index tess))
         (species (list (make-species :name :c :diffusivity 1d0) (make-species :name :s))))
    (is (equalp #2A((0.3d0 1d0) (0.3d0 1d0))
                (capacity-array index species (list (cons :c (constantly 0.3d0))))))
    (signals error (capacity-array index species (list (cons :c (constantly -1d0)))))
    (signals error (capacity-array index species (list (cons :nope (constantly 1d0)))))))
