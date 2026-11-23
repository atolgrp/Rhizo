(in-package :rhizo-diffusion)

;;; PNG rendering for CA-grown Voronoi grids, via zpng (pure Lisp, no
;;; external C library) - the old code rendered through lispbuilder-sdl +
;;; SDL_gfx, and SDL_gfx's underlying C library has since been dropped
;;; from this OS entirely (see the session notes on getting the original
;;; notebook running). Rendering straight to PNG files sidesteps that dead
;;; dependency altogether and needs no live display, which also makes it
;;; usable from a batch/headless run.

(defun id->color (id)
  "A deterministic pseudo-random RGB color for ID - same id always gets
the same color, unlike the old code's mutable *color-dict* (a separate
random color chosen once and remembered per run). Works for any id type
(integers from CA growth, strings from GIS) since it hashes via SXHASH."
  (let ((h (logand (* (logior 1 (sxhash id)) 2654435761) #xFFFFFF)))
    (list (ldb (byte 8 0) h) (ldb (byte 8 8) h) (ldb (byte 8 16) h))))

(defun render-voronoi-grid-png (grid path)
  "Write a PATH PNG coloring every pixel of GRID by its Voronoi cell id -
directly comparable to the old show-voronoi window. Returns PATH."
  (let* ((width (voronoi-grid-width grid))
         (height (voronoi-grid-height grid))
         (png (make-instance 'zpng:png :width width :height height :color-type :truecolor))
         (data (zpng:data-array png)))
    (dotimes (x width)
      (dotimes (y height)
        (destructuring-bind (r g b) (id->color (voronoi-grid-id-at grid x y))
          (setf (aref data y x 0) r (aref data y x 1) g (aref data y x 2) b))))
    (zpng:write-png png path)
    path))

(defparameter *background-color* '(0 0 0))
(defparameter *frontier-color* '(55 55 55)
  "Dim gray for cells the front has already passed through - not in the
old notebook's rendering (it only ever showed the current instant's front,
not a trail), but a cheap, useful addition: it makes the overall shape of
the covered region visible at a glance, which is exactly what the isotropy
question is about.")
(defparameter *active-color* '(255 255 255))
(defparameter *source-color* '(50 120 255))

(defun render-front-propagation-png (grid world path &key source-ids)
  "Write a PATH PNG of WORLD's current front-propagation state over
GRID's cell layout: SOURCE-IDS in blue, the current tick's active front in
white, previously-reached cells in dim gray, everything else black.
Directly comparable to the old show-front-propagation/show-voronoi-
diffusion windows. Returns PATH."
  (let* ((width (voronoi-grid-width grid))
         (height (voronoi-grid-height grid))
         (png (make-instance 'zpng:png :width width :height height :color-type :truecolor))
         (data (zpng:data-array png)))
    (dotimes (x width)
      (dotimes (y height)
        (let* ((id (voronoi-grid-id-at grid x y))
               (color (cond ((member id source-ids :test 'equal) *source-color*)
                            ((gethash id (world-molecules world)) *active-color*)
                            ((gethash id (world-frontier world)) *frontier-color*)
                            (t *background-color*))))
          (destructuring-bind (r g b) color
            (setf (aref data y x 0) r (aref data y x 1) g (aref data y x 2) b)))))
    (zpng:write-png png path)
    path))
