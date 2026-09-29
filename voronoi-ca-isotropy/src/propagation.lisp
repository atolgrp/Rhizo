;;;; propagation.lisp --- Front propagation across the region graph.
;;;;
;;;; This is the cellular automaton proper. State lives on the Voronoi regions,
;;;; not on the lattice cells: a region is either part of the advancing front,
;;;; already visited, or untouched. One time step advances the front to every
;;;; adjacent region that has not been reached before.
;;;;
;;;; Because the automaton runs on the region adjacency graph rather than on the
;;;; lattice, the directional bias of the lattice is not inherited directly --
;;;; that is the hypothesis the experiment tests.

(in-package #:voronoi-ca)

(defvar *front* '()
  "Region IDs currently on the advancing front, most recently added first.")

(defvar *front-set* (make-hash-table :test #'eql)
  "Membership mirror of *FRONT*. Kept exactly in sync so membership tests are
O(1) instead of a linear scan; the list is retained because downstream results
depend on its order.")

(defvar *visited* '()
  "Region IDs the front has already passed through.")

(defvar *visited-set* (make-hash-table :test #'eql)
  "Membership mirror of *VISITED*.")

(defun reset-front ()
  "Clear all front state."
  (setf *front* '()
        *visited* '()
        *front-set* (make-hash-table :test #'eql)
        *visited-set* (make-hash-table :test #'eql)))

(defun seed-front (radius)
  "Seed the front at the centre of the lattice, expanded RADIUS graph hops out.
Returns the list of seeded region IDs."
  (let* ((centre (list (floor (1- (array-dimension *world* 0)) 2)
                       (floor (1- (array-dimension *world* 1)) 2)))
         (centre-id (aref *world* (first centre) (second centre)))
         (out (list centre-id)))
    (dotimes (i radius)
      (dolist (id (flatten (unique out)))
        (push (region-neighbours id) out)))
    (let ((ids (flatten (unique out))))
      (setf *front* ids)
      (setf *front-set* (make-hash-table :test #'eql))
      (dolist (id ids) (setf (gethash id *front-set*) t))
      ids)))

(defun step-front ()
  "Advance the front by one time step.

Every region adjacent to a front region, and not already on the front or
visited, becomes the new front. Regions processed in this step drop out.

This is a restructuring of the original loop, which repeatedly rebuilt the front
list with REMOVE while iterating it -- quadratic in the front size, and ~60% of
run time after the neighbour-table fix. Because every element of the old front
is removed by the end of the sweep, the surviving list is exactly the newly
pushed regions in reverse push order, so accumulating them directly is
equivalent. Membership is checked against the live set, as before."
  (let ((new '()))
    (dolist (cell *front*)
      (dolist (nb (region-neighbours cell))
        (unless (or (gethash nb *front-set*)
                    (gethash nb *visited-set*))
          (push nb new)
          (setf (gethash nb *front-set*) t)
          (push nb *visited*)
          (setf (gethash nb *visited-set*) t)))
      (remhash cell *front-set*))
    (setf *front* new)))

(declaim (inline on-front-p))
(defun on-front-p (id)
  "True when region ID is currently on the front."
  (gethash id *front-set*))

(defun active-neighbour-count (id)
  "How many of region ID's neighbours are currently on the front.

This is the measurement probe: a measurement point registers the front's passage
by counting how many of its adjacent regions are active."
  (let ((n 0))
    (dolist (nb (region-neighbours id) n)
      (when (on-front-p nb) (incf n)))))
