;;;; measurement.lisp --- Measurement points and the single-trial loop.
;;;;
;;;; Sixteen probes are placed on two concentric circles around the source, at
;;;; the four axial and four diagonal compass directions. All eight probes on a
;;;; circle are at the same Euclidean distance from the source -- the diagonal
;;;; ones are offset by radius/sqrt(2) in each axis -- so any systematic
;;;; difference between axial and diagonal arrival times is lattice anisotropy
;;;; and not a difference in distance.

(in-package #:voronoi-ca)

(defparameter +probe-names+
  '("N" "E" "S" "W" "NW" "NE" "SW" "SE"
    "N2" "E2" "S2" "W2" "NW2" "NE2" "SW2" "SE2")
  "Probe labels, in CSV column order.")

(defparameter +axial-probes+ '(:n :e :s :w)
  "The four axial directions.")

(defun probe-coordinates ()
  "Return the 16 probe coordinates as a list of (X Y), in CSV column order.

Ring 1 sits at radius = width/4, ring 2 at radius = 4*(width/10)."
  (let* ((end-x (1- (array-dimension *world* 0)))
         (end-y (1- (array-dimension *world* 1)))
         (half-x (floor end-x 2))
         (half-y (floor end-y 2))
         (r1 (floor end-x 4))
         (r2 (* (floor end-x 10) 4)))
    (flet ((ring (r)
             ;; diagonal offset places the diagonal probes on the same circle
             (let ((d (round (/ r (sqrt 2)))))
               (list :n  (list half-x (- half-y r))
                     :e  (list (+ half-x r) half-y)
                     :s  (list half-x (+ half-y r))
                     :w  (list (- half-x r) half-y)
                     :nw (list (- half-x d) (- half-y d))
                     :ne (list (+ half-x d) (- half-y d))
                     :sw (list (- half-x d) (+ half-y d))
                     :se (list (+ half-x d) (+ half-y d))))))
      (let ((a (ring r1))
            (b (ring r2)))
        (list (getf a :n) (getf a :e) (getf a :s) (getf a :w)
              (getf a :nw) (getf a :ne) (getf a :sw) (getf a :se)
              (getf b :n) (getf b :e) (getf b :s) (getf b :w)
              (getf b :nw) (getf b :ne) (getf b :sw) (getf b :se))))))

(defun measurement-ids ()
  "Region IDs under each of the 16 probes, in CSV column order."
  (mapcar (lambda (xy) (aref *world* (first xy) (second xy)))
          (probe-coordinates)))

(defun run-trial (&key (width *width*) (height *height*) (n-seeds *n-seeds*)
                       snapshot (snapshot-scale 1))
  "Run one complete trial and return the 16 per-probe front-passage times.

Returns two values: the list of times in CSV column order, and the number of
distinct regions actually realised by the tessellation.

The recorded time is the LAST step at which the probe had any active neighbour,
i.e. the moment the front finishes passing the probe. The original notebook
named these variables `peak' but reset the running maximum on every iteration,
so it never recorded a peak. The behaviour is preserved here because it is
applied identically to all sixteen probes, which leaves the axial/diagonal
comparison valid; the name has been corrected to match what is computed."
  (let ((*width* width)
        (*height* height)
        (*n-seeds* n-seeds))
    (setf *world* (make-world width height)
          *counter* 0
          *region-ids* '()
          *region-colours* '()
          *seed-points* '())
    (setf *cell-list* (all-cells *world*))
    ;; scatter seeds and colours (colour draws are part of the RNG sequence)
    (push (make-seeds *n-seeds* *world*) *seed-points*)
    (push (random-colour-list *n-seeds*) *region-colours*)
    (build-tessellation)
    (build-neighbour-table)
    (reset-front)
    (seed-front 1)
    (let* ((probes (measurement-ids))
           (n-probes (length probes))
           ;; -1 means "the front never reached this probe"
           (last-active (make-array n-probes :initial-element -1))
           ;; Run until the front dies out rather than for a fixed number of
           ;; steps. The original stopped at width/3, which is fine for a coarse
           ;; tessellation but too few steps for the front to cross a true
           ;; lattice, where it advances only one cell per step -- the outer
           ;; ring would silently never be reached.
           (max-steps (* 2 width))
           (step 0)
           (snapped nil))
      (loop while (and *front* (< step max-steps)) do
        (step-front)
        (loop for p in probes
              for k from 0
              do (when (> (active-neighbour-count p) 0)
                   (setf (aref last-active k) step)))
        ;; Snapshot the instant the front first touches ring 1. Taking it at a
        ;; fixed step instead would compare the two conditions at different
        ;; stages of propagation, since the Voronoi front crosses the lattice
        ;; in far fewer steps.
        (when (and snapshot (not snapped)
                   (loop for p in (subseq probes 0 8)
                         thereis (> (active-neighbour-count p) 0)))
          (write-tessellation-ppm snapshot :scale snapshot-scale)
          (setf snapped t))
        (incf step))
      (values (coerce last-active 'list)
              (hash-table-count *neighbour-table*)))))
