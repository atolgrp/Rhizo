;;;; utils.lisp --- Small list / numeric helpers.
;;;;
;;;; These are faithful to the original notebook implementations except where
;;;; noted. Where an implementation was changed for performance, the change is
;;;; order- and value-preserving, so results are bit-identical.

(in-package #:voronoi-ca)

(defun unique (list)
  "Remove duplicates from LIST, keeping the FIRST occurrence of each element.

Behaviourally identical to the original O(n^2) `member'-based version, but uses a
hash table so it is O(n). This matters: the tessellation seeds a list of up to
40,000 coordinate pairs through this function."
  (let ((seen (make-hash-table :test #'equal))
        (out '()))
    (dolist (x list (nreverse out))
      (unless (gethash x seen)
        (setf (gethash x seen) t)
        (push x out)))))

(defun flatten (tree)
  "Flatten TREE into a list of atoms."
  (cond ((null tree) nil)
        ((atom tree) (list tree))
        (t (loop for branch in tree appending (flatten branch)))))

(defun coord-pair-p (list)
  "True when LIST looks like a single (X Y) coordinate pair."
  (cond ((null list) nil)
        ((and (atom (car list)) (atom (cadr list))) t)
        (t nil)))

(defun flatten-to-pairs (tree)
  "Flatten TREE down to the level of (X Y) coordinate pairs."
  (cond ((null tree) nil)
        ((coord-pair-p tree) (list tree))
        (t (append (flatten-to-pairs (car tree))
                   (flatten-to-pairs (cdr tree))))))

(defun pick-random-element (sequence)
  "Return a uniformly random element of SEQUENCE."
  (nth (random (length sequence)) sequence))

(defun mean (list)
  "Arithmetic mean of LIST, or NIL when empty."
  (let ((n (length list)))
    (unless (zerop n)
      (float (/ (reduce #'+ list) n)))))

(defun copy-2d-array (array)
  "Return a fresh copy of the 2-D ARRAY, preserving element type."
  (let ((new (make-array (array-dimensions array)
                         :element-type (array-element-type array))))
    (dotimes (i (array-total-size array) new)
      (setf (row-major-aref new i) (row-major-aref array i)))))

(defun write-csv-row (stream values)
  "Write VALUES to STREAM as one comma-separated line."
  (loop for v in values
        for first = t then nil
        do (unless first (write-char #\, stream))
           (princ v stream))
  (terpri stream))

(defun seed-rng (seed)
  "Make the global random state deterministic for SEED.

Determinism is what makes the regression test against the original notebook
possible, and what makes published results reproducible."
  #+sbcl (setf *random-state* (sb-ext:seed-random-state seed))
  #-sbcl (progn
           (warn "Deterministic seeding is only implemented for SBCL; ~
                  results will not be reproducible on this implementation.")
           (setf *random-state* (make-random-state t)))
  seed)
