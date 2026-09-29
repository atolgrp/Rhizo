;;;; ppm.lisp --- Dependency-free image output.
;;;;
;;;; The original rendered the lattice through lispbuilder-sdl, which binds
;;;; SDL 1.2 and needs a display. Writing binary PPM instead costs nothing, has
;;;; no dependencies, runs headless, and any image tool will convert it.

(in-package #:voronoi-ca)

(defun colour-lookup-table ()
  "Map region ID -> RGB triple, built from the region registry."
  (let ((table (make-hash-table :test #'eql)))
    (dolist (entry *region-ids* table)
      (let ((id (first entry))
            (colour (second entry)))
        (unless (gethash id table)
          (setf (gethash id table) colour))))))

(defun wash-out (rgb amount)
  "Blend RGB towards white by AMOUNT (0..1)."
  (mapcar (lambda (c) (round (+ (* c (- 1 amount)) (* 255 amount))))
          (or rgb '(128 128 128))))

(defun write-tessellation-ppm (path &key (scale 1)
                                         (front-colour '(214 40 40))
                                         (swept-colour '(70 90 120)))
  "Write the current propagation state to PATH as a binary PPM (P6).

Three layers, which together show both the medium and the result:

  * unvisited cells keep their region's random colour, washed out, so the
    granularity of the tessellation stays visible but recessive;
  * cells the front has already swept are filled in SWEPT-COLOUR -- the shape
    of this region is the experiment's result, square on a lattice and round
    on a Voronoi tessellation;
  * cells on the current front are drawn in FRONT-COLOUR.

SCALE magnifies each lattice cell to a SCALE x SCALE block."
  (let* ((w (array-dimension *world* 0))
         (h (array-dimension *world* 1))
         (colours (colour-lookup-table))
         (out-w (* w scale))
         (out-h (* h scale)))
    (with-open-file (stream path :direction :output
                                 :element-type '(unsigned-byte 8)
                                 :if-exists :supersede
                                 :if-does-not-exist :create)
      (let ((header (format nil "P6~%~d ~d~%255~%" out-w out-h)))
        (loop for ch across header do (write-byte (char-code ch) stream)))
      (dotimes (row out-h)
        (dotimes (col out-w)
          (let* ((i (floor col scale))
                 (j (floor row scale))
                 (id (aref *world* i j))
                 (rgb (cond ((on-front-p id) front-colour)
                            ((gethash id *visited-set*) swept-colour)
                            (t (wash-out (gethash id colours) 0.62)))))
            (write-byte (mod (or (first rgb) 0) 256) stream)
            (write-byte (mod (or (second rgb) 0) 256) stream)
            (write-byte (mod (or (third rgb) 0) 256) stream)))))
    path))
