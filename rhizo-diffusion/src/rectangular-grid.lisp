(in-package :rhizo-diffusion)

;;; Rectangular finite-volume grids (PLAN.md, Phase 5).
;;;
;;; The Richards benchmark of Celia, Bouloutas & Zarba (1990) is a 1D
;;; vertical column; a grid one cell wide is that column, and the same
;;; builder gives structured 2D grids. Cells are numbered row-major (id = j
;;; NX + i), hold their values at their centres, and y grows downward, as on
;;; the pixel grids, so row 0 is the :TOP. On a rectangle the face between
;;; two neighbours is perpendicular to the segment between their centres,
;;; so the two-point flux is consistent, as on a Voronoi grid.

(defun make-rectangular-tessellation (nx ny &key (dx 1d0) (dy 1d0))
  "Return the finite-volume HASH-GRAPH-TESSELLATION of NX x NY rectangular
cells of size DX x DY: cell j NX + i at ((i + 1/2) DX, (j + 1/2) DY), area
DX DY, faces of length DY between horizontal and DX between vertical
neighbours, and BOUNDARY-FACES on all four sides (distance DY/2 or DX/2).
EDGE-WEIGHT is the centre-to-centre distance."
  (let* ((dx (float dx 1d0)) (dy (float dy 1d0))
         (ids (alexandria:iota (* nx ny)))
         (position-of (lambda (id) (multiple-value-bind (j i) (floor id nx)
                                     (cons (* (+ i 0.5d0) dx) (* (+ j 0.5d0) dy)))))
         (positions (positions-table ids position-of))
         (adjacency (make-hash-table :test 'equal))
         (areas (make-hash-table :test 'equal))
         (edge-lengths (make-hash-table :test 'equal)))
    (dolist (id ids)
      (ensure-cell! adjacency id)
      (setf (gethash id areas) (* dx dy))
      (multiple-value-bind (j i) (floor id nx)
        (when (< (1+ i) nx)
          (add-edge! adjacency id (1+ id) dx)
          (set-edge-length! edge-lengths id (1+ id) dy))
        (when (< (1+ j) ny)
          (add-edge! adjacency id (+ id nx) dy)
          (set-edge-length! edge-lengths id (+ id nx) dx))))
    (make-fv-tessellation
     adjacency :cell-areas areas :edge-lengths edge-lengths :cell-positions positions
               :boundary-faces
               (list (cons :top (loop for i below nx collect (list i dx (* 0.5d0 dy))))
                     (cons :bottom (loop for i below nx
                                         collect (list (+ i (* nx (1- ny))) dx (* 0.5d0 dy))))
                     (cons :left (loop for j below ny collect (list (* j nx) dy (* 0.5d0 dx))))
                     (cons :right (loop for j below ny
                                        collect (list (+ (* j nx) (1- nx)) dy (* 0.5d0 dx))))))))
