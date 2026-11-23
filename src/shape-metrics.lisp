(in-package :rhizo-diffusion)

;;; Cell-shape regularity, measured directly from the pixel grid - no
;;; diffusion involved. This is a separate question from front-propagation
;;; isotropy: front-propagation isotropy asks "does a diffusing front reach
;;; all directions at the same rate," which is a macroscopic question about
;;; many cells together and can stay fine even if individual cells are a
;;; bit blocky. Shape regularity asks "how close is each individual cell to
;;; a nice, roughly-round polygon," which is what actually breaks down
;;; first as cells shrink toward the pixel grid's own resolution - there
;;; just aren't enough pixels per cell left for the alternating-neighborhood
;;; growth to average out the underlying grid's 8-direction bias.
;;;
;;; The isoperimetric quotient Q = 4*pi*Area/Perimeter^2 is 1 for a perfect
;;; circle and smaller for anything else (an axis-aligned block, however
;;; large, tops out around pi/4 ~ 0.785 - so this measures irregularity
;;; beyond that, not just "is it round" in an absolute sense that unfairly
;;; penalizes blocky-but-consistent shapes).

(defun cell-shape-metrics (grid)
  "Return an alist of cell-id -> (area . perimeter) for every cell in
GRID, computed directly from its pixel assignment. AREA is the pixel
count; PERIMETER counts 4-connected pixel-edge boundaries (an edge between
two orthogonally adjacent pixels of different ids, or a pixel on the
grid's outer edge)."
  (let ((area (make-hash-table :test 'equal))
        (perimeter (make-hash-table :test 'equal))
        (width (voronoi-grid-width grid))
        (height (voronoi-grid-height grid)))
    (dotimes (x width)
      (dotimes (y height)
        (let ((id (voronoi-grid-id-at grid x y)))
          (incf (gethash id area 0))
          (dolist (d '((1 . 0) (-1 . 0) (0 . 1) (0 . -1)))
            (let ((nx (+ x (car d))) (ny (+ y (cdr d))))
              (if (in-bounds-p nx ny width height)
                  (when (/= id (voronoi-grid-id-at grid nx ny))
                    (incf (gethash id perimeter 0)))
                  (incf (gethash id perimeter 0))))))))
    (loop for id being the hash-keys of area
          collect (cons id (cons (gethash id area) (gethash id perimeter))))))

(defun isoperimetric-quotient (area perimeter)
  (if (zerop perimeter) 0.0d0 (/ (* 4.0d0 pi area) (* perimeter perimeter))))

(defun mean-isoperimetric-quotient (grid)
  "Average isoperimetric quotient across every cell in GRID - a single
scalar summarizing how round/regular the tessellation's cells are overall.
Higher is more regular; a perfect circle would be 1.0, a large clean
square tessellation tops out around 0.785, and small pixel-quantized
cells score progressively lower and more erratically."
  (let ((metrics (cell-shape-metrics grid)))
    (/ (reduce #'+ (mapcar (lambda (m) (isoperimetric-quotient (float (cadr m) 1.0d0)
                                                                (float (cddr m) 1.0d0)))
                           metrics))
       (float (length metrics) 1.0d0))))
