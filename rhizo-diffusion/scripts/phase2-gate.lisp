;;; Phase 2 gate report (PLAN.md): every gate measurement in one run.
;;;
;;; The pass/fail versions are FiveAM tests (tests/signalling-analysis-tests,
;;; tests/quorum-tests, tests/calling-distance-tests); this script runs the
;;; fuller versions for the record. Every grid is seeded; seeds are printed.
;;;
;;; Run with (about 5 minutes):
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase2-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defun gate-parameters (decay-length k-prime &key (induced 0d0))
  "AHL literature parameters with hydrolysis set for DECAY-LENGTH (um),
sorption K' = K-PRIME exchanging at the decay rate, induced rate INDUCED."
  (let* ((d (parameter-value (ahl-literature-parameters) :diffusivity))
         (k (/ d (* decay-length decay-length))))
    (with-parameter-values (ahl-literature-parameters)
      `((:hydrolysis-rate . ,k) (:desorption-rate . ,k) (:induced-production . ,induced)
        (:kd . ,(/ (* k-prime (parameter-value (ahl-literature-parameters) :water-content))
                   (parameter-value (ahl-literature-parameters) :bulk-density)))))))

(defun print-parameters (title parameters)
  (format t "~&~a~%" title)
  (dolist (row (parameter-table-rows parameters))
    (destructuring-bind (name value unit source note) row
      (declare (ignore note))
      (format t "  ~22a ~12,4e ~14a ~a~%" name value unit
              (subseq source 0 (min 70 (length source)))))))

(print-parameters "== Parameters: AHL (3-oxo-C10-HSL, P. putida IsoF) ==" (ahl-literature-parameters))
(format t "  flagged (no source): ~a~%" (flagged-parameters (ahl-literature-parameters)))
(print-parameters "== Parameters: root exudate (5-deoxystrigol placeholder) ==" (exudate-literature-parameters))
(format t "  flagged (no source): ~a~%" (flagged-parameters (exudate-literature-parameters)))
(finish-output)

(format t "~&== 1. Steady state around a single producer vs S/(2 pi D) K0(r/lambda) ==~%")
(format t "  400 um square, lambda 30 um, annulus 30-90 um (lambda to 3 lambda), :IMEX~%")
(let ((rows (loop for (n seeds) in '((400 (1 2 3)) (1600 (1 2 3)) (6400 (1 2 3)))
                  collect (list n (mapcar (lambda (seed)
                                            (single-producer-steady-state
                                             (voronoi-grid-tessellation
                                              (make-nearest-seed-grid 400 400 n :seed seed))
                                             :centre '(200d0 . 200d0)
                                             :parameters (gate-parameters 30d0 1d0)
                                             :relative-diffusivity :no-tortuosity
                                             :r-min 30d0 :r-max 90d0))
                                          seeds)
                                seeds))))
  (loop for (n runs seeds) in rows
        for l2 = (/ (reduce #'+ runs :key (lambda (r) (getf r :l2))) (length runs))
        do (format t "  ~5d cells (meshes ~a)  L2 ~5,2f%  worst mass ~,1e  worst residual ~,1e  steps ~a~%"
                   n seeds (* 100 l2)
                   (reduce #'max runs :key (lambda (r) (getf r :mass-balance)))
                   (reduce #'max runs :key (lambda (r) (getf r :residual)))
                   (mapcar (lambda (r) (getf r :steps)) runs)))
  (flet ((l2-of (n) (let ((runs (second (assoc n rows))))
                      (/ (reduce #'+ runs :key (lambda (r) (getf r :l2))) (length runs)))))
    (format t "  observed order 400 -> 1600 cells (meshes 1-3): ~,2f~%" (observed-order (l2-of 400) (l2-of 1600) 4))
    (format t "  observed order 1600 -> 6400 cells (meshes 1-3): ~,2f~%" (observed-order (l2-of 1600) (l2-of 6400) 4))))
(dolist (k-prime '(0d0 4d0))
  (let ((run (single-producer-steady-state
              (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1))
              :centre '(200d0 . 200d0) :parameters (gate-parameters 30d0 k-prime)
              :relative-diffusivity :no-tortuosity :r-min 30d0 :r-max 90d0)))
    (format t "  K' = ~3,1f (1600 cells, mesh 1): L2 ~5,3f%  steps ~d  (sorption must not move it)~%"
            k-prime (* 100 (getf run :l2)) (getf run :steps))))
(let ((run (single-producer-steady-state
            (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1))
            :centre '(200d0 . 200d0) :parameters (gate-parameters 30d0 1d0 :induced 6.389d-13)
            :relative-diffusivity :no-tortuosity :density 5d16 :r-min 30d0 :r-max 90d0)))
  (format t "  with induction (dense producer, source self-induced): L2 ~5,3f%  source ~,3e nM um^2/s~%"
          (* 100 (getf run :l2)) (getf run :source-rate)))
(finish-output)

(format t "~&== 2. Quorum induction of a uniform population ==~%")
(dolist (degraders '(0d0 1d8 1d9))
  (let ((switch (quorum-switch-density (ahl-literature-parameters) :degrader-density degraders)))
    (format t "  degraders ~8,1e /L: ~a  on ~,4e /L (AHL ~,2f nM at the fold)  off ~,4e /L~%"
            degraders (getf switch :kind) (getf switch :on) (getf switch :c-on) (getf switch :off))))
(let* ((p (with-parameter-values (ahl-literature-parameters) '((:kd . 0d0))))
       (predicted (getf (quorum-switch-density p) :on))
       (simulated (simulated-switch-density p :lo (* 0.95d0 predicted) :hi (* 1.05d0 predicted)
                                              :tolerance 1d-3)))
  (format t "  simulated switch (no sorption, bisection to 0.1%): ~,4e /L, predicted ~,4e (~,3f%)~%"
          simulated predicted (* 100 (- (/ simulated predicted) 1))))
(let* ((p (gate-parameters 30d0 4d0))
       (p (with-parameter-values (ahl-literature-parameters)
            `((:kd . ,(parameter-value p :kd)))))
       (on (getf (quorum-switch-density p) :on)))
  (format t "  with sorption K' = 4: 1% below the fold induced: ~a; 1% above: ~a~%"
          (uniform-induction-p p (* 0.99d0 on)) (uniform-induction-p p (* 1.01d0 on))))
(let ((p (ahl-literature-parameters)))
  (format t "  literature K' = ~,1f: steady state takes R / k = ~,1f years (hydrolysis only)~%"
          (partition-coefficient p)
          (/ (+ 1 (partition-coefficient p)) (parameter-value p :hydrolysis-rate) 86400 365.25)))
(finish-output)

(format t "~&== 3. Calling distance: analytic tables (literature parameters) ==~%")
(dolist (scenario '((:relative-diffusivity :no-tortuosity)
                    (:relative-diffusivity :no-tortuosity :degrader-density 1d9)))
  (let ((base (calling-distance-scenario (ahl-literature-parameters) scenario)))
    (format t "~&  scenario ~a (defaults: ~a)~%" scenario (default-calling-scenario))
    (format t "    steady ~,1f um, day-one ~,1f um; lambda ~,1f um, R ~,1f, relaxation time ~,2f h~%"
            (getf base :steady) (getf base :transient) (getf base :decay-length)
            (getf base :retardation) (/ (getf base :capacity) (getf base :decay-rate) 3600))
    (format t "    ~22a ~12a  ~8a  ~8a~%" "parameter" "base" "steady" "day one")
    (dolist (row (calling-distance-sensitivity (ahl-literature-parameters) :scenario scenario))
      (destructuring-bind (name value steady transient) row
        (format t "    ~22a ~12,4e  ~:[    -   ~;~:*~8,3f~]  ~:[    -   ~;~:*~8,3f~]~%"
                name value steady transient)))))
(finish-output)

(format t "~&== 3b. Calling distance: simulated vs analytic (scaled: lambda 30 um, 1600 cells) ==~%")
(defun simulated-calling-distance (tess parameters density threshold source)
  (let* ((k (parameter-value parameters :hydrolysis-rate))  ; of the dissolved amount
         (world (make-ahl-world tess :parameters parameters :dt (/ 0.5d0 k)
                                     :relative-diffusivity :no-tortuosity
                                     :initial-state (lambda (id)
                                                      (when (equal id source)
                                                        `((:producer . ,density)))))))
    (step-to-steady-state! world :ahl)
    (calling-distance world '(:per :ahl :water) (cell-position tess source) threshold)))
(let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1)))
       (source (nearest-cell tess '(200d0 . 200d0)))
       (p (gate-parameters 30d0 0d0))
       (theta (parameter-value p :water-content))
       (d (* theta (parameter-value p :diffusivity)))  ; bulk form, free water
       (k (* theta (parameter-value p :hydrolysis-rate)))
       (density 5d14)
       (s (* (cell-area tess source) density (parameter-value p :basal-production)))
       (threshold (point-source-steady-state s d k 60d0))
       (analytic (steady-calling-distance-elasticities s d k threshold))
       (f 1.2d0))
  (format t "  base: simulated ~,2f um, analytic 60.00 um (mesh 1)~%"
          (simulated-calling-distance tess p density threshold source))
  (let ((e-source (/ (log (/ (simulated-calling-distance tess p (* density f) threshold source)
                             (simulated-calling-distance tess p (/ density f) threshold source)))
                     (* 2 (log f))))
        (e-decay (/ (log (/ (simulated-calling-distance
                             tess (with-parameter-values
                                      p `((:hydrolysis-rate . ,(* (/ k theta) f))))
                             density threshold source)
                            (simulated-calling-distance
                             tess (with-parameter-values
                                      p `((:hydrolysis-rate . ,(/ (/ k theta) f))))
                             density threshold source)))
                    (* 2 (log f)))))
    (format t "  elasticity to source strength: simulated ~,3f, analytic ~,3f~%" e-source
            (getf analytic :source-rate))
    (format t "  elasticity to decay rate:      simulated ~,3f, analytic ~,3f~%" e-decay
            (getf analytic :decay-rate))))
(finish-output)
