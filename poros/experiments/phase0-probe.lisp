;;;; phase0-probe.lisp --- Phase 0, step 1.
;;;;
;;;; Establishes the facts the plan rests on, so they can be re-checked on your
;;;; machine rather than taken on trust:
;;;;
;;;;   * Femlisp loads from Quicklisp on this SBCL
;;;;   * it solves a 2D convection-diffusion-reaction model problem
;;;;   * which FL.* packages and entry points are actually present
;;;;
;;;; Run:  sbcl --dynamic-space-size 4096 --script experiments/phase0-probe.lisp
;;;;
;;;; Verified 2026-10-01 on SBCL 2.2.9 / Quicklisp dist 2026-01-01:
;;;; femlisp loaded in ~17 s and the solve completed.
;;;;
;;;; NEXT, and the actual point of Phase 0: this only shows that Femlisp runs.
;;;; It does NOT show that a mass-conservative scheme is expressible, which is
;;;; the question that decides the architecture. See PLAN.md Phase 0, step 3 --
;;;; the no-flux mass-balance test is the measurement that matters.

(require :asdf)

(let ((ql (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname))))
  (if (probe-file ql)
      (load ql)
      (error "Quicklisp not found at ~a" ql)))

(defun quiet-load (system)
  (let ((start (get-internal-real-time)))
    (let ((*standard-output* (make-broadcast-stream))
          (*error-output* (make-broadcast-stream)))
      (handler-bind ((warning #'muffle-warning))
        (funcall (intern "QUICKLOAD" "QL") system :silent t)))
    (/ (- (get-internal-real-time) start) internal-time-units-per-second)))

(format t "~&== loading femlisp ==~%")
(finish-output)
(format t "   loaded in ~,1f s~%" (quiet-load "femlisp"))

(format t "~&~%== FL.* packages present ==~%")
(let ((names (sort (remove-if-not
                    (lambda (n) (and (> (length n) 3) (string= "FL." (subseq n 0 3))))
                    (mapcar #'package-name (list-all-packages)))
                   #'string<)))
  (format t "   ~d packages: ~{~a~^ ~}~%" (length names) names))

(format t "~&~%== entry points ==~%")
(dolist (spec '(("FL.CDR" "CDR-MODEL-PROBLEM")
                ("FL.STRATEGY" "SOLVE")
                ("FL.PROBLEM" "BLACKBOARD")
                ("FL.DOMAINS" "N-CUBE-DOMAIN")
                ("FL.MESH" "UNIFORMLY-REFINED-HIERARCHICAL-MESH")))
  (let* ((pkg (find-package (first spec)))
         (sym (and pkg (find-symbol (second spec) pkg))))
    (format t "   ~16a ~36a ~a~%" (first spec) (second spec)
            (cond ((null pkg) "PACKAGE MISSING")
                  ((null sym) "symbol missing")
                  (t "present")))))

(format t "~&~%== solve a 2D CDR model problem ==~%")
(finish-output)
(handler-case
    (let* ((make  (find-symbol "CDR-MODEL-PROBLEM" "FL.CDR"))
           (solve (find-symbol "SOLVE" "FL.STRATEGY"))
           (board (find-symbol "BLACKBOARD" "FL.PROBLEM")))
      (if (and make solve board)
          (let ((bb (funcall board :problem (funcall make 2)
                                   :success-if '(> :nr-levels 3))))
            (let ((*standard-output* (make-broadcast-stream)))
              (funcall solve bb))
            (format t "   SOLVE COMPLETED~%"))
          (format t "   entry points missing; cannot attempt~%")))
  (error (e)
    (format t "   SOLVE FAILED: ~a~%"
            (substitute #\Space #\Newline (princ-to-string e)))))

(format t "~&~%Phase 0 step 1 complete. Step 3 (mass balance) is the real gate.~%")
