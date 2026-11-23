;;; Phase 1 gate report (PLAN.md): every gate measurement in one run.
;;;
;;; The pass/fail versions are FiveAM tests (tests/accuracy-tests.lisp,
;;; tests/implicit-diffusion-tests.lisp); this script runs the fuller
;;; versions for the record - 10 replicates and both radii for isotropy, as
;;; in the published study (voronoi-ca-isotropy/README.md: 200x200, one seed
;;; per 10 px, probes at r = 49 and 76, seeds 20210728-37; its Voronoi CA
;;; front has an anisotropy index of 1.62% and 2.55%).
;;;
;;; Every run is seeded; seeds are printed with the results.
;;;
;;; Run with (about 4 minutes):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase1-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defparameter *isotropy-seeds* (alexandria:iota 10 :start 20210728))

(defun mean-of (key plists)
  (/ (reduce #'+ plists :key (lambda (p) (getf p key))) (length plists)))

(defun ca-front-anisotropy (seed radius)
  "Anisotropy index of the CA hop front on a CA-grown grid, as in the
published study: probes are the cells containing the 8 compass points at
RADIUS from the centre; arrival is the tick the front first reaches them."
  (let* ((grid (make-voronoi-grid 200 200 4000 :seed seed))
         (probes (mapcar (lambda (p) (cons (car p) (voronoi-grid-id-at grid (round (cadr p))
                                                                         (round (cddr p)))))
                         (compass-points '(100 . 100) radius)))
         (world (make-world (voronoi-grid-tessellation grid) :seed seed)))
    (seed-molecules-at! world (list (voronoi-grid-id-at grid 100 100)))
    (anisotropy-index (run-front-propagation-with-metrics world probes :max-ticks 1000))))

(format t "~&== 1. Gaussian accuracy on nearest-seed grids (D = 1, meshes 1-3) ==~%")
(let ((rows (mapcar (lambda (size)
                      (cons (second size)
                            (grid-gaussian-errors #'make-nearest-seed-grid (first size) (second size) '(1 2 3))))
                    '((150 100) (300 400) (600 1600)))))
  (dolist (row rows)
    (format t "  ~5d cells  error ~5,2f%  D_eff ~,4f  worst mass error ~,1e~%"
            (car row) (* 100 (mean-of :l2 (cdr row))) (mean-of :d-eff (cdr row))
            (reduce #'max (cdr row) :key (lambda (r) (getf r :mass-error)))))
  (format t "  observed order 400 -> 1600 cells: ~,2f~%"
          (observed-order (mean-of :l2 (cdr (second rows))) (mean-of :l2 (cdr (third rows))) 4)))
(finish-output)

(format t "~&== 2. Isotropy, 200x200, 4000 cells, seeds ~a-~a ==~%"
        (first *isotropy-seeds*) (car (last *isotropy-seeds*)))
(dolist (radius '(49d0 76d0))
  (let ((ca (mapcar (lambda (s) (ca-front-anisotropy s radius)) *isotropy-seeds*))
        (fv (mapcar (lambda (s)
                      (point-source-isotropy
                       (voronoi-grid-tessellation (make-nearest-seed-grid 200 200 4000 :seed s))
                       :centre '(100d0 . 100d0) :radius radius :scheme :explicit))
                    *isotropy-seeds*)))
    (format t "  r = ~2d  CA hop front (CA-grown):        anisotropy ~5,2f%~%"
            (round radius) (* 100 (/ (reduce #'+ ca) (length ca))))
    (format t "          FV diffusion (nearest-seed): anisotropy ~5,2f% raw, ~5,2f% normalised, CoV ~5,2f%, arrival ~,4f of analytic~%"
            (* 100 (mean-of :raw fv)) (* 100 (mean-of :normalised fv)) (* 100 (mean-of :cov fv))
            (mean-of :mean-ratio fv))
    (finish-output)))

(format t "~&== 3. Retardation R = 1 + rho_b K_d / theta (400 cells, mesh seed 3) ==~%")
(let ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 300 300 400 :seed 3))))
  (dolist (k '(0d0 1d0 4d0 9d0))
    (dolist (spec '((:explicit 1) (:implicit 5)))
      (let ((run (retardation-measurement tess k :centre '(150d0 . 150d0) :sigma-start 15d0 :sigma-end 30d0
                                                 :scheme (first spec) :time-step-factor (second spec))))
        (format t "  K' = ~3,1f  ~9a  R expected ~5,2f  measured ~,4f  (~,3f%)  mass error ~,1e~%"
                k (first spec) (getf run :r-expected) (getf run :r-measured)
                (* 100 (- (/ (getf run :r-measured) (getf run :r-expected)) 1)) (getf run :mass-error))
        (finish-output)))))

(format t "~&== 4. Implicit vs explicit (Gaussian, meshes 1-3, 1600 cells) ==~%")
(dolist (spec '((:explicit 1) (:implicit 1) (:implicit 10) (:implicit 100)))
  (let ((runs (mapcar (lambda (seed)
                        (gaussian-spreading-error
                         (voronoi-grid-tessellation (make-nearest-seed-grid 600 600 1600 :seed seed))
                         :centre '(300d0 . 300d0) :sigma-start 30d0 :sigma-end 60d0
                         :scheme (first spec) :time-step-factor (second spec)))
                      '(1 2 3))))
    (format t "  ~9a x~3d explicit step: error ~5,2f%  D_eff ~,4f  worst mass error ~,1e  ~5d steps~%"
            (first spec) (second spec) (* 100 (mean-of :l2 runs)) (mean-of :d-eff runs)
            (reduce #'max runs :key (lambda (r) (getf r :mass-error))) (getf (first runs) :steps))
    (finish-output)))
