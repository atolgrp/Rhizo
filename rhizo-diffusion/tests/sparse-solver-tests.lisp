(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Non-symmetric sparse solves (PLAN.md, Phase 5: the Newton Jacobian of
;;; Richards' equation and implicit upwind advection are not symmetric).

(defun dense-solve (matrix b)
  "Test helper: solve MATRIX x = B by Gaussian elimination with partial
pivoting; MATRIX is a list of rows, B a list. Returns a list."
  (let* ((n (length b))
         (a (make-array (list n (1+ n)))))
    (loop for row in matrix for i from 0
          do (loop for v in row for j from 0 do (setf (aref a i j) (float v 1d0)))
             (setf (aref a i n) (float (nth i b) 1d0)))
    (dotimes (k n)
      (let ((p (loop with best = k for i from k below n
                     when (> (abs (aref a i k)) (abs (aref a best k))) do (setf best i)
                     finally (return best))))
        (dotimes (j (1+ n)) (rotatef (aref a k j) (aref a p j)))
        (loop for i from (1+ k) below n
              do (let ((f (/ (aref a i k) (aref a k k))))
                   (loop for j from k to n do (decf (aref a i j) (* f (aref a k j))))))))
    (let ((x (make-list n :initial-element 0d0)))
      (loop for i from (1- n) downto 0
            do (setf (nth i x) (/ (- (aref a i n) (loop for j from (1+ i) below n
                                                          sum (* (aref a i j) (nth j x))))
                                  (aref a i i))))
      x)))

(defun csr-test-system (stencil diagonal off-diagonal)
  "Test helper: the dense rows of the CSR matrix (STENCIL, DIAGONAL,
OFF-DIAGONAL) of APPLY-CSR-MATRIX."
  (let* ((n (length diagonal))
         (rows (loop for i below n collect (make-list n :initial-element 0d0)))
         (offsets (flux-stencil-offsets stencil))
         (neighbors (flux-stencil-neighbors stencil)))
    (dotimes (i n rows)
      (setf (nth i (nth i rows)) (aref diagonal i))
      (loop for k from (aref offsets i) below (aref offsets (1+ i))
            do (incf (nth (aref neighbors k) (nth i rows)) (aref off-diagonal k))))))

(defun upwind-test-matrix (stencil)
  "Test helper: (DIAGONAL OFF-DIAGONAL) of a diagonally dominant,
non-symmetric matrix on STENCIL - diffusion plus a one-sided coupling."
  (let* ((n (length (flux-stencil-areas stencil)))
         (offsets (flux-stencil-offsets stencil))
         (neighbors (flux-stencil-neighbors stencil))
         (diagonal (make-array n :element-type 'double-float :initial-element 1d0))
         (off (make-array (length neighbors) :element-type 'double-float)))
    (dotimes (i n (list diagonal off))
      (loop for k from (aref offsets i) below (aref offsets (1+ i))
            for j = (aref neighbors k)
            for w = (+ 0.5d0 (if (< j i) 0.7d0 0d0)) ; stronger pull from lower ids
            do (setf (aref off k) (- w))
               (incf (aref diagonal i) w)))))

(test apply-csr-matrix-equals-the-dense-product
  (let* ((tess (make-rectangular-tessellation 3 2))
         (stencil (make-flux-stencil tess (make-cell-index tess))))
    (destructuring-bind (diagonal off) (upwind-test-matrix stencil)
      (let* ((x (coerce '(1d0 -2d0 0.5d0 3d0 0d0 -1d0) '(simple-array double-float (*))))
             (rows (csr-test-system stencil diagonal off))
             (y (apply-csr-matrix stencil diagonal off x)))
        (loop for row in rows for i from 0
              do (is (< (abs (- (aref y i) (loop for a in row for xj across x sum (* a xj))))
                        1d-14)))))))

(test bicgstab-solves-a-non-symmetric-system
  (let* ((tess (make-rectangular-tessellation 5 4))
         (stencil (make-flux-stencil tess (make-cell-index tess)))
         (n 20))
    (destructuring-bind (diagonal off) (upwind-test-matrix stencil)
      (let* ((b (make-array n :element-type 'double-float))
             (inverse (map '(simple-array double-float (*)) (lambda (d) (/ 1d0 d)) diagonal)))
        (dotimes (i n) (setf (aref b i) (float (- (mod (* 7 i) 5) 2) 1d0)))
        (multiple-value-bind (x iterations residual)
            (bicgstab (lambda (v) (apply-csr-matrix stencil diagonal off v))
                      b (make-array n :element-type 'double-float :initial-element 0d0)
                      :inverse-diagonal inverse :tolerance 1d-13)
          (is (<= residual 1d-13))
          (is (< 0 iterations n))
          (loop for xi across x for expected in (dense-solve (csr-test-system stencil diagonal off)
                                                             (coerce b 'list))
                do (is (< (abs (- xi expected)) 1d-11))))))))

(test bicgstab-agrees-with-conjugate-gradient-on-a-symmetric-system
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 40 40 50 :seed 3)))
         (stencil (make-flux-stencil tess (make-cell-index tess)))
         (n (length (flux-stencil-areas stencil)))
         (d (uniform-column n 1d0))
         (coefficients (face-coefficients stencil 5d0 d))
         (diagonal (capacity-diagonal stencil nil))
         (apply (lambda (v) (apply-implicit-matrix stencil coefficients diagonal v)))
         (b (make-array n :element-type 'double-float)))
    (dotimes (i n) (setf (aref b i) (float (mod i 3) 1d0)))
    (let ((x-cg (conjugate-gradient apply b (copy-seq b) :tolerance 1d-14))
          (x-bi (bicgstab apply b (copy-seq b) :tolerance 1d-14)))
      (is (every (lambda (u v) (< (abs (- u v)) 1d-10)) x-cg x-bi)))))

(test bicgstab-returns-zero-for-a-zero-right-hand-side
  (let ((x (bicgstab #'copy-seq (uniform-column 4 0d0) (uniform-column 4 3d0))))
    (is (every #'zerop x))))

(test linear-solvers-signal-their-failure
  ;; b must not be constant: the test matrix has unit row sums, so x = 1
  ;; would solve b = 1 at once
  (let* ((tess (make-rectangular-tessellation 6 6))
         (stencil (make-flux-stencil tess (make-cell-index tess)))
         (b (make-array 36 :element-type 'double-float)))
    (dotimes (i 36) (setf (aref b i) (float (- (mod (* 7 i) 5) 2) 1d0)))
    (destructuring-bind (diagonal off) (upwind-test-matrix stencil)
      (signals linear-solver-failure
        (bicgstab (lambda (v) (apply-csr-matrix stencil diagonal off v))
                  b (uniform-column 36 0d0) :tolerance 1d-15 :max-iterations 2))
      (signals linear-solver-failure
        (conjugate-gradient (lambda (v) (apply-csr-matrix stencil diagonal off v))
                            b (uniform-column 36 0d0) :tolerance 1d-15 :max-iterations 1)))))

(test jacobi-preconditioner-scales-by-the-inverse-diagonal
  (let ((precondition (jacobi-preconditioner (coerce '(2d0 0.5d0) '(simple-array double-float (*)))))
        (v (coerce '(3d0 4d0) '(simple-array double-float (*)))))
    (is (equalp #(6d0 2d0) (funcall precondition v)))
    (is (equalp #(3d0 4d0) (funcall (jacobi-preconditioner nil) v)) "NIL means no preconditioning")
    (is (not (eq v (funcall (jacobi-preconditioner nil) v))) "always a fresh vector")))
