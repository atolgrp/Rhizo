;;; Finite-volume diffusion accuracy on raster Voronoi grids - report.
;;;
;;; The pass/fail check is now a FiveAM test,
;;; FV-GAUSSIAN-ACCURACY-ON-NEAREST-SEED-GRIDS (tests/accuracy-tests.lisp),
;;; built on the functions in src/accuracy.lisp. This script prints the same
;;; measurement as a table for both grid builders, at three sizes. A
;;; Gaussian of width W/20 diffuses with D = 1 until its width doubles and is
;;; compared with the analytic solution; 3 meshes per row, ~225 px per cell.
;;;
;;; History. On 2026-11-23 this script found that the finite-volume flux
;;; used EDGE-WEIGHT (1.0 on every CA edge) as the cell-to-cell distance:
;;; 92% error and D_eff 10.2 at every resolution, mass still conserved.
;;; Phase 1 steps 1-4 fixed it (FACE-DISTANCE, seed positions, the
;;; nearest-seed builder, Euclidean face lengths). Measured 2026-11-24:
;;;
;;;   cells  ca-grown (staircase faces)   nearest-seed
;;;   100    35.0%, D_eff 1.33            46.1%, D_eff 1.005
;;;   400    11.4%, D_eff 1.34             9.1%, D_eff 0.998
;;;   1600   13.7%, D_eff 1.33             2.3%, D_eff 1.003
;;;
;;; Nearest-seed converges at second order; CA-grown does not, because its
;;; staircase face lengths overstate diagonal faces by 4/pi on average.
;;; The CA-grown grid stays for the isotropy study, not reaction-diffusion.
;;;
;;; Run with (about 40 seconds; the CA-grown 600x600 growth is the slow part):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/fv-accuracy.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defun report-row (label grid-builder width n-cells)
  "Print the mean error, mean D_eff and worst mass error over meshes 1-3."
  (let ((runs (grid-gaussian-errors grid-builder width n-cells '(1 2 3))))
    (flet ((mean (key) (/ (reduce #'+ runs :key (lambda (r) (getf r key))) (length runs))))
      (format t "  ~13a ~5d cells  error ~5,1f%   D_eff ~5,3f   worst mass error ~,1e~%"
              label n-cells (* 100 (mean :l2)) (mean :d-eff)
              (reduce #'max runs :key (lambda (r) (getf r :mass-error))))
      (finish-output))))

(format t "~&Gaussian, D = 1, width W/20 -> W/10, meshes 1-3:~%")
(dolist (size '((150 100) (300 400) (600 1600)))
  (destructuring-bind (width n-cells) size
    (report-row "ca-grown" #'make-voronoi-grid width n-cells)
    (report-row "nearest-seed" #'make-nearest-seed-grid width n-cells)))
