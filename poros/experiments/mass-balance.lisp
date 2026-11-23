;;;; mass-balance.lisp --- Phase 0, step 3: the measurement that decides the architecture.
;;;;
;;;; Transport above a flow solver is only meaningful if the scheme conserves
;;;; mass. This measures whether it does.
;;;;
;;;; Setup: pure diffusion on the unit square, no source, and NO-FLUX boundaries
;;;; (passing :dirichlet NIL leaves the natural boundary condition, which for
;;;; this operator is zero normal flux). Nothing enters or leaves, so the total
;;;; mass must not change.
;;;;
;;;; Initial condition:
;;;;
;;;;     u0(x,y) = 1 + 1/2 cos(pi x) cos(pi y)
;;;;
;;;; chosen because the cosine terms integrate to zero over (0,1), so
;;;;
;;;;     integral of u0 = 1   exactly.
;;;;
;;;; The L2 projection onto a Lagrange space preserves that integral exactly --
;;;; constants lie in the test space, so the projection error is orthogonal to 1
;;;; -- which means any drift away from 1 is a conservation failure of the
;;;; time-stepping scheme, not an artefact of representing the initial state.
;;;;
;;;; Diffusion drives u towards the uniform state 1, so the solution changes a
;;;; great deal while its integral must not change at all. That is what makes
;;;; this a sharp test rather than a tautology.
;;;;
;;;; WHAT THIS DOES AND DOES NOT SETTLE
;;;;
;;;; Continuous Galerkin is GLOBALLY conservative by construction: summing all
;;;; equations telescopes the interior fluxes, so this test is expected to pass.
;;;; It does not establish LOCAL conservation -- that fluxes balance on each
;;;; cell -- which is the property porous-media transport actually needs, and
;;;; which standard Lagrange elements do not have. A pass here is necessary,
;;;; not sufficient. See local-conservation.lisp for the discriminating test.

(require :asdf)
(load (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname)))
(let ((*standard-output* (make-broadcast-stream))
      (*error-output* (make-broadcast-stream)))
  (handler-bind ((warning #'muffle-warning))
    (funcall (intern "QUICKLOAD" "QL") "femlisp" :silent t)))

(in-package :fl.application)

(defparameter *exact-mass* 1d0
  "Integral of u0 over the unit square, exactly 1 by construction.")

(defun scalar-value (x)
  "Unwrap FE-INTEGRATE's #(#M((v))) into a double float."
  (cond ((numberp x) (float x 1d0))
        ((and (vectorp x) (not (stringp x)) (plusp (length x)))
         (scalar-value (aref x 0)))
        (t (float (fl.matlisp:mref x 0 0) 1d0))))

(defun u0 (x)
  (+ 1d0 (* 0.5d0 (cos (* pi (aref x 0))) (cos (* pi (aref x 1))))))

(defun no-flux-diffusion-problem ()
  "Pure diffusion, no source, no Dirichlet data -> natural (no-flux) boundary."
  (cdr-model-problem 2
                     :initial #'u0
                     :source nil
                     :dirichlet nil))

(defun run-steps (&key (order 2) (levels 3) (time-step 0.01d0) (steps 20))
  "Integrate forward and return the total mass at the end."
  (let ((bb (blackboard :problem (no-flux-diffusion-problem)
                        :fe-class (lagrange-fe order)
                        :number-of-levels levels
                        :time-step time-step
                        :number-of-time-steps steps
                        :plot-mesh nil :plot nil :output 0)))
    (let ((*standard-output* (make-broadcast-stream)))
      (solve bb))
    (scalar-value (fe-integrate (getbb bb :solution)))))

(defun main ()
  (format t "~&~%no-flux diffusion, u0 = 1 + 1/2 cos(pi x) cos(pi y)~%")
  (format t "exact total mass (conserved) = ~,15f~%~%" *exact-mass*)
  (format t "~6a ~7a ~8a ~20a ~14a ~12a~%"
          "order" "levels" "steps" "total mass" "abs error" "rel error")
  (format t "~a~%" (make-string 72 :initial-element #\-))
  (dolist (cfg '((2 3 0.01d0 1) (2 3 0.01d0 10) (2 3 0.01d0 50)
                 (1 3 0.01d0 50) (2 4 0.01d0 50) (2 3 0.001d0 200)))
    (destructuring-bind (order levels dt steps) cfg
      (handler-case
          (let* ((m (run-steps :order order :levels levels
                               :time-step dt :steps steps))
                 (err (abs (- m *exact-mass*))))
            (format t "~6d ~7d ~8d ~20,15f ~14,3e ~12,3e~%"
                    order levels steps m err (/ err *exact-mass*)))
        (error (e)
          (format t "~6d ~7d ~8d  FAILED: ~a~%" order levels steps
                  (let ((s (substitute #\Space #\Newline (princ-to-string e))))
                    (subseq s 0 (min 60 (length s))))))))
    (finish-output))
  (format t "~%Global conservation is necessary but NOT sufficient -- see the~%")
  (format t "header of this file and local-conservation.lisp.~%"))

(main)
