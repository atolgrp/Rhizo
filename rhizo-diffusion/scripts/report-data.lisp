;;; Data for the report's graphs of Phases 2-4, written as JSON into
;;; RHIZO_REPORT_DIR (Phase 5's come from scripts/phase5-gate.lisp). Every
;;; grid and layout is seeded. About 10 minutes.
;;;
;;;   RHIZO_REPORT_DIR=dir/ sbcl --dynamic-space-size 4096 --non-interactive --load scripts/report-data.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defparameter *report-dir* (or (uiop:getenv "RHIZO_REPORT_DIR") (error "set RHIZO_REPORT_DIR")))

(defun vec (list) (coerce list 'vector))

(defun dump-json (name object)
  (with-open-file (s (merge-pathnames (format nil "~a.json" name)
                                      (uiop:ensure-directory-pathname *report-dir*))
                     :direction :output :if-exists :supersede :if-does-not-exist :create)
    (cl-json:encode-json object s))
  (format t "wrote ~a~%" name)
  (finish-output))

;;; Phase 2: steady profile around one producer against K0 (the Phase 3
;;; engine-check setup: decay length 30 um at theta 0.3, no sorption).
(let* ((base (ahl-literature-parameters))
       (xi (relative-diffusivity-model :millington-quirk base))
       (theta 0.3d0)
       (d (* 490d0 (funcall xi theta)))
       (k (/ d (* theta 900d0)))
       (parameters (with-parameter-values base `((:hydrolysis-rate . ,k) (:kd . 0d0)
                                                 (:induced-production . 0d0) (:water-content . ,theta))))
       (tess (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1)))
       (source (nearest-cell tess '(200d0 . 200d0)))
       (density 5d14)
       (s-rate (* (cell-area tess source) density (parameter-value parameters :basal-production)))
       (world (make-ahl-world tess :parameters parameters
                                   :dt (/ 0.5d0 (+ k (sorption-relaxation-rate parameters)))
                                   :initial-state (lambda (id) (when (equal id source)
                                                                 `((:producer . ,density)))))))
  (step-to-steady-state! world :ahl)
  (let ((profile (radial-profile world (pore-water :ahl) (cell-position tess source) 10d0)))
    (dump-json "phase2-k0"
               `((:r . ,(vec (mapcar #'car profile)))
                 (:simulated . ,(vec (mapcar #'cdr profile)))
                 (:analytic . ,(vec (mapcar (lambda (p) (point-source-steady-state s-rate d (* theta k) (car p)))
                                            profile)))))))

;;; Phase 2: well-mixed quorum response - steady AHL against producer density.
(let* ((parameters (ahl-literature-parameters))
       (densities (loop for e from 7d0 to 9.5d0 by 0.025d0 collect (expt 10d0 e))))
  (dump-json "phase2-quorum"
             `((:switch-density . ,(quorum-switch-density parameters))
               (:threshold . ,(parameter-value parameters :induction-threshold))
               (:points . ,(vec (mapcan (lambda (rho)
                                          (mapcar (lambda (c) (vec (list rho c)))
                                                  (uniform-steady-states parameters rho)))
                                        densities))))))

;;; Phase 2: calling-distance elasticities.
(let ((rows (calling-distance-sensitivity (ahl-literature-parameters))))
  (dump-json "phase2-sensitivity"
             (vec (mapcar (lambda (row)
                            (destructuring-bind (name base steady transient) row
                              `((:name . ,(string-downcase name)) (:base . ,base)
                                (:steady . ,steady) (:transient . ,transient))))
                          rows))))

;;; Phase 3: calling distance against water content, three relations, two scenarios.
(let ((thetas (loop for th from 0.05d0 to 0.4501d0 by 0.01d0 collect th)))
  (dump-json "phase3-ranges"
             (vec (loop for (label scenario) in '(("abiotic" ()) ("1e9 degraders/L" (:degrader-density 1d9)))
                        nconc (loop for relation in '(:millington-quirk :papendick-campbell :mullins-sommer)
                                    collect (let ((rows (range-versus-water (ahl-literature-parameters) thetas
                                                                            :scenario (list* :relative-diffusivity
                                                                                             relation scenario))))
                                              `((:scenario . ,label) (:relation . ,(string-downcase relation))
                                                (:theta . ,(vec thetas))
                                                (:steady . ,(vec (mapcar (lambda (r) (getf r :steady)) rows)))
                                                (:day-one . ,(vec (mapcar (lambda (r) (getf r :transient)) rows)))))))))
  ;; K_d sensitivity of the day-one range (Millington-Quirk)
  (let ((kd (parameter-value (ahl-literature-parameters) :kd)))
    (dump-json "phase3-kd"
               (vec (loop for factor in '(0.1d0 1d0 10d0)
                          collect (let ((rows (range-versus-water
                                               (with-parameter-values (ahl-literature-parameters)
                                                 (list (cons :kd (* factor kd))))
                                               thetas :scenario '(:relative-diffusivity :millington-quirk))))
                                    `((:factor . ,factor) (:theta . ,(vec thetas))
                                      (:day-one . ,(vec (mapcar (lambda (r) (getf r :transient)) rows))))))))))

;;; Phase 4: percolation against water content, 600 patches on 4 mm, 20 layouts
;;; (Millington-Quirk, 1e9 degraders/L, steady; the Phase 4 soil defaults).
(let ((thetas (loop for th from 0.18d0 to 0.3201d0 by 0.005d0 collect th))
      (seeds (alexandria:iota 20 :start 1)))
  (dump-json "phase4-percolation"
             (vec (mapcar (lambda (theta)
                            (let* ((p (with-parameter-values (ahl-literature-parameters)
                                        `((:water-content . ,theta))))
                                   (rows (mapcar (lambda (seed)
                                                   (let ((positions (random-patch-positions 600 4000d0 4000d0
                                                                                            :seed seed)))
                                                     (graph-connectivity positions
                                                                         (flow-signalling-edges
                                                                          positions p 0d0 '(:degrader-density 1d9)
                                                                          :steady)
                                                                         4000d0 400d0)))
                                                 seeds)))
                              `((:theta . ,theta)
                                (:range . ,(getf (calling-distance-scenario p '(:degrader-density 1d9)) :steady))
                                (:spanning . ,(float (/ (count-if (lambda (r) (getf r :spanning)) rows) 20)))
                                (:weak . ,(float (/ (reduce #'+ rows :key (lambda (r) (getf r :largest-weak-fraction)))
                                                    20)))
                                (:degree . ,(float (/ (reduce #'+ rows :key (lambda (r) (getf r :mean-degree))) 20))))))
                          thetas))))
