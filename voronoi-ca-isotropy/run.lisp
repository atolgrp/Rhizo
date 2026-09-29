;;;; run.lisp --- Entry point. No Quicklisp, no external dependencies.
;;;;
;;;;   sbcl --script run.lisp             # 200x200, 5 replicates per condition
;;;;   sbcl --script run.lisp 300 10      # 300x300, 10 replicates per condition

(setf *load-verbose* nil *compile-verbose* nil *compile-print* nil)

(let ((here (directory-namestring *load-truename*)))
  (handler-bind ((warning #'muffle-warning))
    (dolist (file '("package" "utils" "grid" "voronoi"
                    "propagation" "ppm" "measurement" "experiment"))
      (load (merge-pathnames (format nil "src/~a.lisp" file) here)))))

(voronoi-ca:main)
