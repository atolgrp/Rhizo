;;;; grid.lisp --- The square lattice the tessellation is drawn on.
;;;;
;;;; The world is a 2-D fixnum array. Each cell holds the integer ID of the
;;;; Voronoi region it belongs to, or 0 if it is not yet assigned.

(in-package #:voronoi-ca)

(defparameter *width* 200
  "Width of the lattice in cells.")

(defparameter *height* 200
  "Height of the lattice in cells.")

(defvar *world* nil
  "The 2-D fixnum array holding a region ID per lattice cell.")

(defvar *counter* 0
  "Global cell-update counter. Its parity alternates the neighbourhood used
during tessellation growth, which is what keeps the growing regions from
becoming diamond-shaped.")

(defun make-world (&optional (width *width*) (height *height*))
  "Allocate a fresh zeroed world array."
  (make-array (list width height) :element-type 'fixnum :initial-element 0))

(defun clear-world (world)
  "Set every cell of WORLD to 0."
  (dotimes (i (array-total-size world) world)
    (setf (row-major-aref world i) 0)))

(defun in-bounds (world coords)
  "Keep only those (X Y) pairs in COORDS that lie inside WORLD."
  (remove nil
          (loop for (x y) in coords
                collect (when (array-in-bounds-p world x y) (list x y)))))

(defun von-neumann-neighbours (i j world distance)
  "The 4 orthogonal neighbours of (I J) at DISTANCE, clipped to WORLD."
  (let ((d distance))
    (when (array-in-bounds-p world i j)
      (in-bounds world
                 (list (list (- i d) j)
                       (list i (- j d))
                       (list i (+ j d))
                       (list (+ i d) j))))))

(defun moore-neighbours (i j world distance)
  "The 8 surrounding neighbours of (I J) at DISTANCE, clipped to WORLD."
  (let ((d distance))
    (when (array-in-bounds-p world i j)
      (in-bounds world
                 (list (list (- i d) (- j d))
                       (list (- i d) j)
                       (list (- i d) (+ j d))
                       (list i (- j d))
                       (list i (+ j d))
                       (list (+ i d) (- j d))
                       (list (+ i d) j)
                       (list (+ i d) (+ j d)))))))

(defun cell-id (world coords)
  "Region ID stored at the (X Y) pair COORDS."
  (aref world (first coords) (second coords)))

(defun all-cells (&optional (world *world*))
  "A list of every (X Y) coordinate in WORLD, in row-major order."
  (let ((out '()))
    (dotimes (x (array-dimension world 0))
      (dotimes (y (array-dimension world 1))
        (push (list x y) out)))
    (nreverse out)))
