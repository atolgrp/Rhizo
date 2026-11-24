(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Signalling range against water content (PLAN.md, Phase 3 gate).

(test range-versus-water-falls-as-the-soil-dries
  ;; The range reached in a day falls as the soil dries, as does the steady
  ;; range under biodegradation. With slow abiotic loss alone the steady
  ;; range (years away) instead peaks at intermediate moisture: in its
  ;; logarithmic regime wetter soil dilutes the signal faster.
  (let ((parameters (ahl-literature-parameters))
        (thetas '(0.15d0 0.2d0 0.3d0 0.35d0 0.4d0 0.45d0)))
    (flet ((ranges (scenario measure)
             (mapcar (lambda (row) (getf row measure))
                     (range-versus-water parameters thetas :scenario scenario))))
      (let ((abiotic '(:relative-diffusivity :millington-quirk))
            (biotic '(:relative-diffusivity :millington-quirk :degrader-density 1d9)))
        (dolist (case (list (list abiotic :transient) (list biotic :transient) (list biotic :steady)))
          (let ((r (apply #'ranges case)))
            (is (every #'< r (rest r)) "~a: ~a" case r)))
        (let ((r (ranges abiotic :steady)))
          (is (< (first r) (fourth r)) "rises from dry ...")
          (is (> (fourth r) (sixth r)) "... and falls again toward saturation"))))))

(test mullins-sommer-range-vanishes-below-its-threshold
  ;; sandy clay loam threshold 0.2505: no connected water, no signalling
  (let* ((parameters (ahl-literature-parameters))
         (rows (range-versus-water parameters '(0.24d0 0.3d0)
                                   :scenario '(:relative-diffusivity :mullins-sommer))))
    (is (= 0d0 (getf (first rows) :steady)))
    (is (= 0d0 (getf (first rows) :transient)))
    (is (plusp (getf (second rows) :steady)))
    (is (plusp (getf (second rows) :transient)))))

(test water-transition-halves-the-wet-range
  (let ((parameters (ahl-literature-parameters)))
    (dolist (relation '(:millington-quirk :mullins-sommer))
      (dolist (measure '(:steady :transient))
        (let* ((scenario (list :relative-diffusivity relation :degrader-density 1d9))
               (theta* (water-transition parameters :scenario scenario :measure measure
                                                    :theta-wet 0.45d0)))
          (flet ((range (theta) (getf (first (range-versus-water parameters (list theta)
                                                                 :scenario scenario))
                                      measure)))
            (is (< 0.1d0 theta* 0.45d0))
            (is (< (abs (- (range theta*) (* 0.5d0 (range 0.45d0)))) (* 1d-6 (range 0.45d0)))
                "~a ~a: range at theta* ~a" relation measure theta*)))))))

(test water-sensitivity-reports-theta-and-kd-elasticities
  (let* ((parameters (ahl-literature-parameters))
         (scenario '(:relative-diffusivity :millington-quirk))
         (steady (water-sensitivity parameters 0.3d0 :scenario scenario :measure :steady))
         (transient (water-sensitivity parameters 0.3d0 :scenario scenario :measure :transient)))
    (is (< (abs (getf steady :kd-elasticity)) 1d-9) "sorption cannot move the steady range")
    (is (zerop (getf steady :ratio)))
    (is (< (getf transient :kd-elasticity) -0.3d0) "but it slows the day-one range")
    (is (plusp (getf transient :theta-elasticity)))
    (is (< (abs (- (getf transient :ratio)
                   (/ (abs (getf transient :kd-elasticity)) (abs (getf transient :theta-elasticity)))))
           1d-12))
    (is (eq (getf transient :co-gated) (> (getf transient :ratio) 2d0)))))

(test simulated-range-follows-the-water-content
  ;; Engine check: one producer cell, no sorption, Millington-Quirk; the
  ;; hydrolysis rate gives a decay length of 30 um at theta = 0.3 and the
  ;; threshold a 60 um range there. At theta = 0.45 the simulated range
  ;; must follow the analytic one (bulk form).
  (let* ((relation :millington-quirk)
         (base (ahl-literature-parameters))
         (xi (relative-diffusivity-model relation base))
         (k (/ (* 490d0 (funcall xi 0.3d0)) (* 0.3d0 900d0)))
         (parameters (with-parameter-values base `((:hydrolysis-rate . ,k) (:kd . 0d0)
                                                   (:induced-production . 0d0))))
         (tess (voronoi-grid-tessellation (make-nearest-seed-grid 400 400 1600 :seed 1)))
         (source (nearest-cell tess '(200d0 . 200d0)))
         (density 5d14)
         (s-rate (* (cell-area tess source) density (parameter-value parameters :basal-production)))
         (threshold (point-source-steady-state s-rate (* 490d0 (funcall xi 0.3d0)) (* 0.3d0 k) 60d0)))
    (dolist (theta '(0.3d0 0.45d0))
      (let ((analytic (steady-calling-distance s-rate (* 490d0 (funcall xi theta)) (* theta k)
                                               threshold))
            (simulated (single-producer-calling-distance
                        tess (with-parameter-values parameters `((:water-content . ,theta)))
                        :centre '(200d0 . 200d0) :density density :threshold threshold
                        :relative-diffusivity relation)))
        (is (< (abs (- simulated analytic)) (* 0.05d0 analytic))
            "theta ~a: simulated ~a um, analytic ~a um" theta simulated analytic)))))
