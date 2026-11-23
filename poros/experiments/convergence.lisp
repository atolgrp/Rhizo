;;;; convergence.lisp --- Phase 0, step 2: does Femlisp converge at the expected rate?
;;;;
;;;; Method of manufactured solutions. On the unit square take
;;;;
;;;;     u(x,y) = sin(pi x) sin(pi y)
;;;;
;;;; which vanishes on the boundary, so homogeneous Dirichlet conditions are
;;;; exact. Then -laplacian u = 2 pi^2 u, so the source term is known, and
;;;;
;;;;     integral of u over (0,1)^2  =  (2/pi)^2  =  4/pi^2  =  0.405284734...
;;;;
;;;; exactly. Comparing the computed integral against that number tests the
;;;; whole chain at once -- problem definition, assembly, solve, and
;;;; FE-INTEGRATE -- against a value that owes nothing to the code under test.
;;;;
;;;; A correct Lagrange discretisation of order p should drive this functional
;;;; error down at a rate of roughly h^(p+1).

(require :asdf)
(load (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname)))
(let ((*standard-output* (make-broadcast-stream))
      (*error-output* (make-broadcast-stream)))
  (handler-bind ((warning #'muffle-warning))
    (funcall (intern "QUICKLOAD" "QL") "femlisp" :silent t)))

(in-package :fl.application)

(defparameter *exact-integral* (/ 4d0 (* pi pi))
  "Integral of sin(pi x) sin(pi y) over the unit square.")

(defun manufactured-problem ()
  "-laplacian u = 2 pi^2 sin(pi x) sin(pi y), u = 0 on the boundary."
  (cdr-model-problem
   2
   :source (lambda (x)
             (* 2d0 pi pi
                (sin (* pi (aref x 0)))
                (sin (* pi (aref x 1)))))
   :dirichlet 0d0))

(defun solve-at-level (order levels)
  "Solve on a uniformly refined mesh and return the solution vector."
  (let ((bb (blackboard :problem (manufactured-problem)
                        :fe-class (lagrange-fe order)
                        :success-if `(>= :nr-levels ,levels)
                        :output 0 :plot nil)))
    (let ((*standard-output* (make-broadcast-stream)))
      (solve bb))
    (getbb bb :solution)))

(defun scalar-value (x)
  "Unwrap whatever FE-INTEGRATE hands back into a double float.

For a scalar problem it returns a SIMPLE-VECTOR with one entry per solution
component, each entry a 1x1 matlisp matrix -- i.e. #(#M((0.4052...))). Neither
layer is a number, so both have to be peeled."
  (cond ((numberp x) (float x 1d0))
        ((and (vectorp x) (not (stringp x)) (plusp (length x)))
         (scalar-value (aref x 0)))
        (t (float (fl.matlisp:mref x 0 0) 1d0))))

(defun integral-of (asv)
  "Integral of the FE function over the domain, as a double float."
  (scalar-value (fe-integrate asv)))

(defun run (&key (order 1) (max-levels 5))
  (format t "~&~%manufactured solution: u = sin(pi x) sin(pi y) on the unit square~%")
  (format t "exact integral of u    = ~,12f~%" *exact-integral*)
  (format t "Lagrange order         = ~d~%~%" order)
  (format t "~6a ~18a ~14a ~8a~%" "levels" "computed integral" "abs error" "rate")
  (format t "~a~%" (make-string 52 :initial-element #\-))
  (let ((prev nil))
    (loop for levels from 2 to max-levels do
      (handler-case
          (let* ((asv (solve-at-level order levels))
                 (val (integral-of asv))
                 (err (abs (- val *exact-integral*))))
            (format t "~6d ~18,12f ~14,3e ~8a~%" levels val err
                    (if (and prev (> err 0))
                        (format nil "~,2f" (/ (log (/ prev err)) (log 2d0)))
                        "--"))
            (setf prev err))
        (error (e)
          (format t "~6d  FAILED: ~a~%" levels
                  (substitute #\Space #\Newline (princ-to-string e)))
          (return))))))

(run :order 1 :max-levels 5)
(run :order 2 :max-levels 4)
(format t "~&~%done.~%")
