(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun chain-tessellation (n)
  "N unit-area cells in a row at (0,0), (1,0), ...; neighbours share a face
of length 1, so every conductance is 1."
  (let ((adjacency (make-hash-table :test 'equal))
        (areas (make-hash-table :test 'equal))
        (lengths (make-hash-table :test 'equal)))
    (dotimes (i n)
      (ensure-cell! adjacency i)
      (setf (gethash i areas) 1d0)
      (when (plusp i)
        (add-edge! adjacency (1- i) i 1d0)
        (set-edge-length! lengths (1- i) i 1d0)))
    (make-fv-tessellation adjacency :cell-areas areas :edge-lengths lengths
                                    :cell-positions (positions-table (alexandria:iota n)
                                                                     (lambda (i) (cons i 0))))))

(defun chain-world (diffusivities initial scheme dt)
  "World on CHAIN-TESSELLATION with species :C, per-cell DIFFUSIVITIES and
INITIAL concentrations (lists, one entry per cell)."
  (make-mechanistic-world
   (chain-tessellation (length diffusivities))
   :species (list (make-species :name :c :diffusivity 0d0)) :reactions '()
   :dt dt :diffusion-scheme scheme
   :diffusivity-field (list (cons :c (lambda (i) (nth i diffusivities))))
   :initial-state (lambda (i) (list (cons :c (nth i initial))))))

(defun chain-values (world)
  (loop for i below (cell-count (mechanistic-world-tessellation world))
        collect (cell-species-value world i :c)))

(defun random-diffusivity-world (scheme dt &key (seed 7))
  "Nearest-seed grid (400 cells) with a Gaussian of :C and a diffusivity
drawn per cell, uniformly from [0.1, 2.1), with a seeded random state."
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 300 300 400 :seed 3)))
         (random-state (sb-ext:seed-random-state seed))
         (field (make-hash-table :test 'equal)))
    (dolist (id (sort (cell-ids tess) #'<))
      (setf (gethash id field) (+ 0.1d0 (random 2d0 random-state))))
    (make-mechanistic-world
     tess :species (list (make-species :name :c :diffusivity 1d0)) :reactions '()
          :dt dt :diffusion-scheme scheme
          :diffusivity-field (list (cons :c (lambda (id) (gethash id field))))
          :initial-state (field-initial-state
                          tess (list (cons :c (gaussian-field '(150d0 . 150d0) 15d0)))))))

(test harmonic-mean-by-hand
  (is (= 2d0 (harmonic-mean 2d0 2d0)))
  (is (= 1.5d0 (harmonic-mean 1d0 3d0)))
  (is (= 1.5d0 (harmonic-mean 3d0 1d0)))
  (is (= 0d0 (harmonic-mean 0d0 5d0)) "a zero-diffusivity cell is a barrier")
  (is (= 0d0 (harmonic-mean 0d0 0d0)))
  (is (typep (harmonic-mean 1 3) 'double-float)))

(test diffusivity-array-defaults-and-overrides
  (let* ((tess (chain-tessellation 3))
         (index (make-cell-index tess))
         (species (list (make-species :name :a :diffusivity 0.5d0)
                        (make-species :name :b :diffusivity 0d0)))
         (array (diffusivity-array index species (list (cons :b (lambda (i) (* 2 i)))))))
    (dotimes (k 3)
      (let ((i (index->cell-id index k)))
        (is (= 0.5d0 (aref array k 0)) "unlisted species keep their own diffusivity")
        (is (= (* 2d0 i) (aref array k 1)))))
    (signals error (diffusivity-array index species (list (cons :b (lambda (i) (declare (ignore i)) -1)))))
    (signals error (diffusivity-array index species (list (cons :nope (lambda (i) i)))))))

(test set-species-diffusivity!-changes-one-species-only
  (let* ((world (make-mechanistic-world
                 (chain-tessellation 3)
                 :species (list (make-species :name :a :diffusivity 0.5d0)
                                (make-species :name :b :diffusivity 0.25d0))
                 :reactions '() :dt 0.1d0)))
    (is (eq world (set-species-diffusivity! world :a (lambda (i) (+ 1 i)))))
    (let ((index (mechanistic-world-cell-index world))
          (array (mechanistic-world-diffusivities world)))
      (dotimes (k 3)
        (is (= (+ 1d0 (index->cell-id index k)) (aref array k 0)))
        (is (= 0.25d0 (aref array k 1)))))
    (signals error (set-species-diffusivity! world :a (lambda (i) (- i 5))))
    (signals error (set-species-diffusivity! world :nope (lambda (i) i)))))

(test two-cells-exchange-at-the-harmonic-mean-diffusivity
  ;; D = 1 and 3 give a face diffusivity of 1.5. Unit areas, face length 1,
  ;; distance 1: k = 1.5. The difference C0 - C1 shrinks by exactly
  ;; (1 - 2 k dt) per explicit step and 1 / (1 + 2 k dt) per implicit step.
  (let ((dt 0.1d0) (n 10) (k 1.5d0))
    (let ((explicit (chain-world '(1d0 3d0) '(10d0 0d0) :explicit dt))
          (implicit (chain-world '(1d0 3d0) '(10d0 0d0) :implicit dt)))
      (run-mechanistic! explicit n)
      (run-mechanistic! implicit n)
      (destructuring-bind (c0 c1) (chain-values explicit)
        (is (< (abs (- (- c0 c1) (* 10 (expt (- 1 (* 2 k dt)) n)))) 1d-12)))
      (destructuring-bind (c0 c1) (chain-values implicit)
        (is (< (abs (- (- c0 c1) (/ 10 (expt (+ 1 (* 2 k dt)) n)))) 1d-12))))))

(test a-zero-diffusivity-cell-is-a-barrier
  ;; Chain 0-1-2-3-4 with D = 0 in cell 2: nothing ever reaches cells 3
  ;; and 4, and cell 2 keeps its own value, under either scheme.
  (dolist (scheme '(:explicit :implicit))
    (let ((world (chain-world '(1d0 1d0 0d0 1d0 1d0) '(4d0 0d0 7d0 0d0 0d0) scheme 0.2d0)))
      (run-mechanistic! world 50)
      (destructuring-bind (c0 c1 c2 c3 c4) (chain-values world)
        (is (< (abs (- 4d0 (+ c0 c1))) 1d-13) "~a: mass stays left of the barrier" scheme)
        (is (< (abs (- c0 c1)) 1d-3) "~a: the left pair has equilibrated" scheme)
        (is (= 7d0 c2) "~a: the barrier cell is untouched" scheme)
        (is (= 0d0 c3 c4) "~a: nothing crosses" scheme)))))

(test a-uniform-field-reproduces-the-species-diffusivity
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 100 100 40 :seed 2)))
         (init (field-initial-state tess (list (cons :c (gaussian-field '(50d0 . 50d0) 6d0))))))
    (flet ((world (species-d field)
             (make-mechanistic-world tess :species (list (make-species :name :c :diffusivity species-d))
                                          :reactions '() :dt 0.5d0 :initial-state init
                                          :diffusivity-field field)))
      (let ((by-species (world 0.7d0 nil))
            (by-field (world 0d0 (list (cons :c (lambda (id) (declare (ignore id)) 0.7d0))))))
        (run-mechanistic! by-species 30)
        (run-mechanistic! by-field 30)
        (is (equalp (mechanistic-world-state by-species) (mechanistic-world-state by-field)))))))

(test world-stable-timestep-follows-the-field
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 100 100 40 :seed 2)))
         (species (list (make-species :name :c :diffusivity 1d0)))
         (world (make-mechanistic-world tess :species species :reactions '() :dt 1d0)))
    (is (< (abs (- (world-stable-timestep world) (stable-timestep tess species))) 1d-12))
    (set-species-diffusivity! world :c (lambda (id) (declare (ignore id)) 2d0))
    (is (< (abs (- (world-stable-timestep world) (/ (stable-timestep tess species) 2))) 1d-12))
    (set-species-diffusivity! world :c (lambda (id) (declare (ignore id)) 0d0))
    (is (= sb-ext:double-float-positive-infinity (world-stable-timestep world)))))

(test heterogeneous-field-conserves-mass-and-stays-bounded
  ;; Both schemes, on a random per-cell field. Explicit at its own stability
  ;; bound, implicit at 50x that bound.
  (let* ((probe (random-diffusivity-world :explicit 1d0))
         (dt (world-stable-timestep probe)))
    (dolist (spec (list (list :explicit dt 200) (list :implicit (* 50 dt) 20)))
      (destructuring-bind (scheme step n) spec
        (let* ((world (random-diffusivity-world scheme step))
               (mass-0 (species-mass world :c))
               (max-0 (reduce #'max (state-column (mechanistic-world-state world) 0))))
          (run-mechanistic! world n)
          (let ((column (state-column (mechanistic-world-state world) 0)))
            (is (< (/ (abs (- (species-mass world :c) mass-0)) mass-0) 1d-13)
                "~a: mass error ~a" scheme (/ (abs (- (species-mass world :c) mass-0)) mass-0))
            (is (>= (reduce #'min column) (* -1d-14 max-0)) "~a: no negative values" scheme)
            (is (<= (reduce #'max column) (* (+ 1 1d-14) max-0)) "~a: no new maximum" scheme)))))))

(test heterogeneous-field-explicit-and-implicit-agree
  (let* ((probe (random-diffusivity-world :explicit 1d0))
         (dt (world-stable-timestep probe))
         (explicit (random-diffusivity-world :explicit dt))
         (implicit (random-diffusivity-world :implicit dt)))
    (run-mechanistic! explicit 100)
    (run-mechanistic! implicit 100)
    (let* ((areas (flux-stencil-areas (mechanistic-world-stencil explicit)))
           (ce (state-column (mechanistic-world-state explicit) 0))
           (ci (state-column (mechanistic-world-state implicit) 0))
           (difference (sqrt (/ (loop for i below (length areas)
                                      sum (* (aref areas i) (expt (- (aref ce i) (aref ci i)) 2)))
                                (loop for i below (length areas)
                                      sum (* (aref areas i) (expt (aref ce i) 2)))))))
      (is (< difference 0.01d0) "fields differ by ~,3f%" (* 100 difference)))))
