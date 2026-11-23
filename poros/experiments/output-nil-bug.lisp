;;;; output-nil-bug.lisp --- Verifying the :OUTPUT NIL crash and its fix.
;;;;
;;;; Three parts:
;;;;   1. reproduce the crash and confirm the root cause in isolation
;;;;   2. check the proposed fix is behaviour-preserving for every valid value
;;;;   3. patch the live image and confirm :OUTPUT NIL then solves
;;;;
;;;; Run: sbcl --dynamic-space-size 4096 --script output-nil-bug.lisp

(load (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname)))
(let ((*standard-output* (make-broadcast-stream))
      (*error-output* (make-broadcast-stream)))
  (handler-bind ((warning #'muffle-warning))
    (funcall (intern "QUICKLOAD" "QL") "femlisp" :silent t)))

(in-package :fl.application)

(defun say (fmt &rest args) (apply #'format t fmt args) (finish-output))

;;; ------------------------------------------------------------------
;;; Part 1 -- root cause, in isolation
;;; ------------------------------------------------------------------
;;; OUTPUT-P's fallback branch is
;;;
;;;   (or (member *output-depth* '(t :all :infinity))
;;;       (<= *iteration-depth* *output-depth*))
;;;
;;; so a NIL *output-depth* reaches <= as its second argument.

(say "~&=== 1. root cause ===~%")
(say "   reproducing the arithmetic directly:~%")
(dolist (depth '(0 1 t :all :infinity nil))
  (let ((fl.iteration::*iteration-depth* 0)
        (fl.iteration::*output-depth* depth))
    (handler-case
        (let ((r (or (member fl.iteration::*output-depth* '(t :all :infinity))
                     (<= fl.iteration::*iteration-depth*
                         fl.iteration::*output-depth*))))
          (say "     *output-depth* = ~8s -> ~a~%" depth (and r t)))
      (error (e)
        (say "     *output-depth* = ~8s -> ~a: ~a~%" depth (type-of e)
             (substitute #\Space #\Newline (princ-to-string e)))))))

;;; ------------------------------------------------------------------
;;; Part 2 -- is the fix behaviour-preserving?
;;; ------------------------------------------------------------------
;;; Proposed: guard the comparison with REALP, so anything that is neither a
;;; real nor an unlimited-sentinel means "no output" instead of erroring.

(defun fallback-original (iter-depth out-depth)
  (or (member out-depth '(t :all :infinity))
      (<= iter-depth out-depth)))

(defun fallback-patched (iter-depth out-depth)
  (or (member out-depth '(t :all :infinity))
      (and (realp out-depth)
           (<= iter-depth out-depth))))

(say "~&=== 2. is the fix behaviour-preserving? ===~%")
(say "   ~14a ~14a ~14a ~14a ~a~%"
     "*output-depth*" "*iter-depth*" "original" "patched" "verdict")
(say "   ~a~%" (make-string 66 :initial-element #\-))
(let ((differences 0) (fixed 0))
  (dolist (out-depth '(0 1 2 t :all :infinity nil))
    (dolist (iter-depth '(0 1 2))
      (let* ((orig (handler-case (list :ok (and (fallback-original iter-depth out-depth) t))
                     (error (e) (list :error (type-of e)))))
             (new  (handler-case (list :ok (and (fallback-patched iter-depth out-depth) t))
                     (error (e) (list :error (type-of e)))))
             (verdict (cond ((equal orig new) "same")
                            ((and (eq (first orig) :error) (eq (first new) :ok))
                             (incf fixed) "FIXED")
                            (t (incf differences) "DIFFERS"))))
        (say "   ~14s ~14d ~14a ~14a ~a~%" out-depth iter-depth
             (if (eq (first orig) :error) (second orig) (second orig))
             (if (eq (first new) :error) (second new) (second new))
             verdict))))
  (say "~%   behaviour changed only where it previously crashed: ~a~%"
       (if (zerop differences) "YES" "NO -- REVIEW"))
  (say "   cases fixed: ~d   unintended differences: ~d~%" fixed differences))

;;; ------------------------------------------------------------------
;;; Part 3 -- patch the live image and re-test end to end
;;; ------------------------------------------------------------------

(say "~&=== 3. end-to-end, before and after the patch ===~%")

(defun try-solve (label)
  (handler-case
      (let ((bb (blackboard :problem (cdr-model-problem 2)
                            :success-if '(> :nr-levels 2)
                            :output nil :plot nil)))
        (let ((*standard-output* (make-broadcast-stream))) (solve bb))
        (say "   ~22a :output nil -> OK (~a)~%" label
             (type-of (getbb bb :solution))))
    (error (e) (say "   ~22a :output nil -> ~a: ~a~%" label (type-of e)
                    (substitute #\Space #\Newline (princ-to-string e))))))

(try-solve "before patch")

;; The fix, applied to the live generic function.
(defmethod fl.iteration::output-p ((iter fl.iteration::<iteration>) blackboard)
  (declare (ignore blackboard))
  (if (slot-boundp iter 'fl.iteration::output)
      (slot-value iter 'fl.iteration::output)
      (or (member fl.iteration::*output-depth* '(t :all :infinity))
          (and (realp fl.iteration::*output-depth*)
               (<= fl.iteration::*iteration-depth*
                   fl.iteration::*output-depth*)))))

(try-solve "after patch")

;; And confirm the previously-working values still work.
(say "~%   regression -- values that already worked:~%")
(dolist (val '(0 1 t))
  (handler-case
      (let ((bb (blackboard :problem (cdr-model-problem 2)
                            :success-if '(> :nr-levels 2)
                            :output val :plot nil)))
        (let ((*standard-output* (make-broadcast-stream))) (solve bb))
        (say "     :output ~10s -> OK~%" val))
    (error (e) (say "     :output ~10s -> ~a~%" val (type-of e)))))
(handler-case
    (let ((bb (blackboard :problem (cdr-model-problem 2)
                          :success-if '(> :nr-levels 2) :plot nil)))
      (let ((*standard-output* (make-broadcast-stream))) (solve bb))
      (say "     :output absent     -> OK~%"))
  (error (e) (say "     :output absent     -> ~a~%" (type-of e))))

(say "~&~%done.~%")
