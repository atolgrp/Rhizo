(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun single-cell-tessellation (area)
  "A tessellation with one cell \"only\" and no edges - isolates local
reaction dynamics from diffusion entirely, for tests that only care about
the reaction-network side of the mechanistic engine."
  (well-mixed-tessellation :area area :id "only"))

(defun two-cell-tessellation (area edge-len distance &key (hop-weight distance))
  "Two cells of AREA joined by one face of length EDGE-LEN, with value
points DISTANCE apart. HOP-WEIGHT is the topological EDGE-WEIGHT, which the
finite-volume engine must ignore; it defaults to DISTANCE."
  (let ((adjacency (make-hash-table :test 'equal))
        (cell-areas (make-hash-table :test 'equal))
        (edge-lengths (make-hash-table :test 'equal))
        (positions (positions-table '("0" "1")
                                    (lambda (id) (if (string= id "0") '(0 . 0) (cons distance 0))))))
    (add-edge! adjacency "0" "1" hop-weight)
    (setf (gethash "0" cell-areas) area)
    (setf (gethash "1" cell-areas) area)
    (set-edge-length! edge-lengths "0" "1" edge-len)
    (make-fv-tessellation adjacency :cell-areas cell-areas :edge-lengths edge-lengths
                                    :cell-positions positions)))

(test field-initial-state-samples-fields-at-cell-positions
  ;; Cell "0" sits at (0, 0) and cell "1" at (4, 0).
  (let* ((tess (two-cell-tessellation 1.0d0 1.0d0 4.0d0))
         (init (field-initial-state tess (list (cons :a (lambda (p) (+ 1 (car p))))
                                               (cons :b (lambda (p) (declare (ignore p)) 7))))))
    (is (equal '((:a . 1.0d0) (:b . 7)) (funcall init "0")))
    (is (equal '((:a . 5.0d0) (:b . 7)) (funcall init "1")))
    (let ((world (make-mechanistic-world
                  tess :species (list (make-species :name :a :diffusivity 0d0)
                                      (make-species :name :b :diffusivity 0d0))
                       :reactions '() :dt 1d0 :initial-state init)))
      (is (= 5.0d0 (cell-species-value world "1" :a)))
      (is (= 7.0d0 (cell-species-value world "0" :b))))))

(test face-conductance-is-edge-length-over-face-distance
  (let ((tess (two-cell-tessellation 1.0d0 3.0d0 2.0d0 :hop-weight 1.0d0)))
    (is (= 1.5d0 (face-conductance tess "0" "1")))
    (is (= 1.5d0 (face-conductance tess "1" "0")))))

(test diffusion-reads-face-distance-not-edge-weight
  ;; The defect of 2026-11-23 in miniature: hop weight 1.0, true distance 4.0.
  ;; Exchange rate k = D * L / (d * area) must use d = 4, so the two cells'
  ;; difference decays as exp(-2kt) with k = 1 * 1 / (4 * 1) = 0.25, not 1.0.
  ;; STABLE-TIMESTEP must likewise scale with d: 1 / (2 * D * L/d) = 2.0.
  (let* ((tess (two-cell-tessellation 1.0d0 1.0d0 4.0d0 :hop-weight 1.0d0))
         (species (list (make-species :name :c :diffusivity 1.0d0)))
         (dt 0.01d0) (n-steps 200)
         (world (make-mechanistic-world
                 tess :species species :reactions '() :dt dt
                      :initial-state (lambda (id) (if (string= id "0") '((:c . 10.0d0)) '((:c . 0.0d0)))))))
    (is (< (abs (- (stable-timestep tess species :safety-factor 1.0d0) 2.0d0)) 1d-12))
    (dotimes (i n-steps) (step-mechanistic! world))
    (let ((difference (- (cell-species-value world "0" :c) (cell-species-value world "1" :c)))
          (expected (* 10.0d0 (exp (* -2.0d0 0.25d0 dt n-steps)))))
      ;; explicit Euler's own error here is 0.009; k = 1.0 would give 0.18 vs 3.68
      (is (< (abs (- difference expected)) 0.02d0)
          "C0 - C1: expected ~a (k from face distance) got ~a" expected difference))))

(defun cell-species-value (world cell-id species-name)
  (aref (mechanistic-world-state world)
        (cell-id->index (mechanistic-world-cell-index world) cell-id)
        (gethash species-name (mechanistic-world-species-index world))))

(defun total-mechanistic-mass (world)
  "Sum of concentration*area over every cell and species - the quantity
STEP-MECHANISTIC! is designed to conserve exactly (see
mechanistic-mass-is-conserved-over-many-ticks below for why)."
  (let* ((tess (mechanistic-world-tessellation world))
         (cell-index (mechanistic-world-cell-index world))
         (state (mechanistic-world-state world))
         (n-cells (cell-index-count cell-index))
         (n-species (mechanistic-world-n-species world))
         (total 0.0d0))
    (dotimes (i n-cells)
      (let ((area (cell-area tess (index->cell-id cell-index i))))
        (dotimes (s n-species) (incf total (* area (aref state i s))))))
    total))

(test two-cell-diffusion-matches-analytic-equilibration
  ;; With no reactions, two equal-area cells connected by one edge reduce to
  ;; the textbook two-compartment linear ODE: d(C1-C2)/dt = -2k(C1-C2) where
  ;; k = D*edge-length/(distance*area) - solvable in closed form, so this is
  ;; a real numerical-accuracy check on DIFFUSION-SUBSTEP!, not just a
  ;; sanity smell test.
  (let* ((tess (two-cell-tessellation 1.0d0 1.0d0 1.0d0))
         (d 0.5d0) (k 0.5d0) ; k = D*edge-len/(distance*area) = 0.5*1/(1*1)
         (dt 0.01d0) (n-steps 500)
         (species (list (make-species :name :c :diffusivity d)))
         (world (make-mechanistic-world
                 tess :species species :reactions '() :dt dt
                      :initial-state (lambda (id)
                                       (if (string= id "0") '((:c . 10.0d0)) '((:c . 0.0d0)))))))
    (dotimes (i n-steps) (step-mechanistic! world))
    (let* ((c1 (cell-species-value world "0" :c))
           (c2 (cell-species-value world "1" :c))
           (total-time (* dt n-steps))
           (analytic-diff (* 10.0d0 (exp (* -2.0d0 k total-time))))
           (expected-c1 (/ (+ 10.0d0 analytic-diff) 2.0d0))
           (expected-c2 (/ (- 10.0d0 analytic-diff) 2.0d0)))
      (is (< (abs (- c1 expected-c1)) 0.05d0) "c1: expected ~a got ~a" expected-c1 c1)
      (is (< (abs (- c2 expected-c2)) 0.05d0) "c2: expected ~a got ~a" expected-c2 c2))))

(test stable-timestep-diverges-above-bound-stays-bounded-below
  ;; STABLE-TIMESTEP is a heuristic (see its docstring) rather than an exact
  ;; threshold, so this checks well inside/outside it rather than right at
  ;; the boundary: a timestep at half the computed bound should stay
  ;; comfortably bounded, and one at 5x it should visibly diverge - proof
  ;; the stability guard is measuring something real, not just present.
  (let* ((tess (two-cell-tessellation 1.0d0 1.0d0 1.0d0))
         (species (list (make-species :name :c :diffusivity 0.5d0)))
         (bound (stable-timestep tess species :safety-factor 1.0d0))
         (initial-state (lambda (id) (if (string= id "0") '((:c . 10.0d0)) '((:c . 0.0d0))))))
    (let ((world (make-mechanistic-world tess :species species :reactions '()
                                               :dt (* 0.5d0 bound) :initial-state initial-state)))
      (dotimes (i 50) (step-mechanistic! world))
      (is (< (abs (cell-species-value world "0" :c)) 100.0d0)
          "well under the stability bound, state should stay close to its starting scale"))
    (let ((world (make-mechanistic-world tess :species species :reactions '()
                                               :dt (* 5.0d0 bound) :initial-state initial-state)))
      (dotimes (i 50) (step-mechanistic! world))
      (is (> (abs (cell-species-value world "0" :c)) 1000.0d0)
          "well over the stability bound, explicit diffusion should visibly blow up"))))
