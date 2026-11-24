(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Advective solute transport with implicit upwinding (PLAN.md, Phase 5).
;;; A species is dissolved - carried by water - when its free-water
;;; diffusivity is above zero (the convention of SET-WATER-CONTENT!), so the
;;; test solute diffuses, if only slightly.

(defun face-flux-vector (world flux-of)
  "Test helper: the vector, in WORLD's stencil face order, of (FUNCALL
FLUX-OF i j) - the water flux out of cell index I into J."
  (let* ((stencil (mechanistic-world-stencil world))
         (offsets (flux-stencil-offsets stencil))
         (neighbors (flux-stencil-neighbors stencil))
         (fluxes (make-array (length neighbors) :element-type 'double-float)))
    (dotimes (i (length (flux-stencil-areas stencil)) fluxes)
      (loop for k from (aref offsets i) below (aref offsets (1+ i))
            do (setf (aref fluxes k) (float (funcall flux-of i (aref neighbors k)) 1d0))))))

(defun one-species-world (tess &key (diffusivity 1d-9) (dt 1d0) initial capacity
                                    (scheme :implicit))
  (make-mechanistic-world tess :species (list (make-species :name :a :diffusivity diffusivity)
                                              (make-species :name :fixed :diffusivity 0d0))
                               :reactions '() :dt dt :diffusion-scheme scheme
                               :initial-state initial
                               :capacity-field (when capacity `((:a . ,(constantly capacity))))))

(defun index-of (world id)
  (cell-id->index (mechanistic-world-cell-index world) id))

(test circulating-flow-keeps-a-uniform-solute-uniform
  ;; a closed loop 0 -> 1 -> 3 -> 2 -> 0 on a 2 x 2 grid: divergence-free
  (let* ((world (one-species-world (make-rectangular-tessellation 2 2)
                                   :initial '((:a . 3d0) (:fixed . 1d0))))
         (loop-ids '((0 . 1) (1 . 3) (3 . 2) (2 . 0)))
         (ids (mechanistic-world-cell-index world)))
    (flet ((flux (i j)
             (let ((a (index->cell-id ids i)) (b (index->cell-id ids j)))
               (cond ((member (cons a b) loop-ids :test #'equal) 0.4d0)
                     ((member (cons b a) loop-ids :test #'equal) -0.4d0)
                     (t 0d0)))))
      (set-flow-field! world (face-flux-vector world #'flux) nil)
      (run-mechanistic! world 5)
      (is (every (lambda (id) (< (abs (- 3d0 (aref (mechanistic-world-state world) (index-of world id) 0)))
                                 1d-11))
                 '(0 1 2 3)))
      (is (= 0d0 (species-exported world :a))))))

(test upwinding-moves-solute-downstream-only
  (let* ((world (one-species-world (make-rectangular-tessellation 2 2)
                                   :initial (lambda (id) `((:a . ,(if (= id 0) 1d0 0d0)) (:fixed . 0d0)))))
         (loop-ids '((0 . 1) (1 . 3) (3 . 2) (2 . 0)))
         (ids (mechanistic-world-cell-index world)))
    (flet ((flux (i j)
             (let ((a (index->cell-id ids i)) (b (index->cell-id ids j)))
               (cond ((member (cons a b) loop-ids :test #'equal) 0.4d0)
                     ((member (cons b a) loop-ids :test #'equal) -0.4d0)
                     (t 0d0)))))
      (set-flow-field! world (face-flux-vector world #'flux) nil)
      (step-mechanistic! world)
      (let ((state (mechanistic-world-state world)))
        ;; implicit upwinding carries some mass several hops in one step,
        ;; falling off with the number of hops downstream: 1, then 3, then 2
        (is (> (aref state (index-of world 1) 0) (aref state (index-of world 3) 0)
               (aref state (index-of world 2) 0) 0d0))
        (is (> (aref state (index-of world 1) 0) 0.1d0) "one hop downstream gains most")
        (is (< (abs (- 1d0 (species-mass world :a))) 1d-12) "closed loop: mass conserved")
        (is (= 0d0 (aref state (index-of world 0) 1)) "an immobile species is not carried")))))

(test outflow-through-a-boundary-exports-solute
  ;; two stacked unit cells, water 0 -> 1 -> out of the bottom at rate Q,
  ;; diffusivity D across their unit face: the exact backward Euler upwind
  ;; solution of the 2 x 2 system
  ;;   (1 + c + dt Q) x0 - c x1 = 1,   -(c + dt Q) x0 + (1 + c + dt Q) x1 = 0.5
  ;; with c = dt D.
  (let* ((q 0.3d0) (dt 2d0) (d 0.2d0) (c (* dt d))
         (world (one-species-world (make-rectangular-tessellation 1 2) :dt dt :diffusivity d
                                   :initial (lambda (id) `((:a . ,(if (= id 0) 1d0 0.5d0))
                                                            (:fixed . 0d0)))))
         (i0 (index-of world 0)) (i1 (index-of world 1))
         (outflow (make-array 2 :element-type 'double-float :initial-element 0d0)))
    (setf (aref outflow i1) q)
    (set-flow-field! world (face-flux-vector world (lambda (i j) (if (= i i0) q (if (= j i0) (- q) 0d0))))
                     outflow)
    (step-mechanistic! world)
    (let* ((a11 (+ 1d0 c (* dt q))) (a12 (- c)) (a21 (- (+ c (* dt q)))) (a22 (+ 1d0 c (* dt q)))
           (det (- (* a11 a22) (* a12 a21)))
           (x0 (/ (- (* 1d0 a22) (* a12 0.5d0)) det))
           (x1 (/ (- (* a11 0.5d0) (* a21 1d0)) det))
           (state (mechanistic-world-state world)))
      (is (< (abs (- (aref state i0 0) x0)) 1d-12))
      (is (< (abs (- (aref state i1 0) x1)) 1d-12))
      (is (< (abs (- (species-exported world :a) (* dt q x1))) 1d-12))
      (is (< (abs (- (pool-balance world '(:a) '()) 1.5d0)) 1d-12)
          "pool balances count exported solute: mass + exports = 1.5"))))

(test advected-pulse-moves-at-the-pore-velocity
  ;; Gaussian pulse in a row of unit cells, capacity 0.5, Darcy flux 0.5:
  ;; pore velocity 1. Implicit upwind: the centroid moves at exactly v; the
  ;; variance grows at 2 (D/kappa + v dx/2 + v^2 dt/2) - measured
  ;; numerical dispersion included, not hidden.
  (let* ((n 200) (dt 0.5d0) (steps 100) (kappa 0.5d0) (d 0.5d0) (v 1d0)
         (world (one-species-world (make-rectangular-tessellation n 1) :diffusivity d :dt dt
                                   :capacity kappa
                                   :initial (lambda (id) `((:a . ,(exp (- (/ (expt (- (+ id 0.5d0) 50d0) 2)
                                                                             50d0))))
                                                            (:fixed . 0d0)))))
         (ids (mechanistic-world-cell-index world))
         (outflow (make-array n :element-type 'double-float :initial-element 0d0)))
    (setf (aref outflow (index-of world (1- n))) (* v kappa))
    (set-flow-field! world (face-flux-vector world (lambda (i j)
                                                     (let ((a (index->cell-id ids i))
                                                           (b (index->cell-id ids j)))
                                                       (if (> b a) (* v kappa) (- (* v kappa))))))
                     outflow)
    (flet ((moments ()
             (let* ((state (mechanistic-world-state world))
                    (cells (loop for i below n collect (cons (+ 0.5d0 (index->cell-id ids i))
                                                             (aref state i 0))))
                    (mass (reduce #'+ cells :key #'cdr))
                    (mean (/ (reduce #'+ cells :key (lambda (c) (* (car c) (cdr c)))) mass)))
               (values mean (/ (reduce #'+ cells :key (lambda (c) (* (expt (- (car c) mean) 2) (cdr c))))
                               mass)))))
      (multiple-value-bind (mean0 var0) (moments)
        (run-mechanistic! world steps)
        (multiple-value-bind (mean1 var1) (moments)
          (let ((time (* dt steps)))
            (is (< (abs (- (- mean1 mean0) (* v time))) 1d-6) "centroid moved ~a" (- mean1 mean0))
            (let ((expected (* 2 time (+ (/ d kappa) (* 0.5d0 v 1d0) (* 0.5d0 v v dt)))))
              (is (< (abs (- (- var1 var0) expected)) (* 0.03d0 expected))
                  "variance grew ~a, expected ~a" (- var1 var0) expected))))))))

(test flow-needs-an-implicit-transport-scheme
  (let ((world (one-species-world (make-rectangular-tessellation 2 1) :scheme :explicit
                                  :diffusivity 0.1d0 :dt 0.1d0
                                  :initial '((:a . 1d0) (:fixed . 0d0)))))
    (set-flow-field! world (face-flux-vector world (lambda (i j) (if (< i j) 0.1d0 -0.1d0))) nil)
    (signals error (step-mechanistic! world))))

(test flow-field-is-checked-and-can-be-cleared
  (let ((world (one-species-world (make-rectangular-tessellation 2 1)
                                  :initial '((:a . 1d0) (:fixed . 0d0)))))
    (signals error (set-flow-field! world (uniform-column 7 0d0) nil) "wrong length")
    (signals error (set-flow-field! world (face-flux-vector world (constantly 0d0))
                                    (uniform-column 2 -1d0))
      "negative outflow")
    (set-flow-field! world nil nil)
    (is (null (mechanistic-world-face-fluxes world)))))
