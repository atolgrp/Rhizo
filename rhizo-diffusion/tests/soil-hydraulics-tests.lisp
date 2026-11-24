(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Van Genuchten-Mualem soil hydraulic relations (PLAN.md, Phase 5).

(defun celia-sand ()
  (van-genuchten-soil (celia-sand-parameters)))

(test van-genuchten-water-content-follows-the-retention-curve
  ;; theta = theta_r + (theta_s - theta_r) [1 + (alpha |h|)^n]^-m, m = 1 - 1/n
  (let ((soil (celia-sand)))
    (is (= 0.368d0 (vg-water-content soil 0d0)) "saturated at h = 0")
    (is (= 0.368d0 (vg-water-content soil 5d0)) "and under positive pressure")
    (let ((expected (+ 0.102d0 (* 0.266d0 (expt (+ 1d0 (expt (* 0.0335d0 1000d0) 2)) -0.5d0)))))
      (is (< (abs (- (vg-water-content soil -1000d0) expected)) 1d-15)))
    (is (< (vg-water-content soil -1d6) 0.1021d0) "tends to theta_r when very dry")))

(test van-genuchten-conductivity-matches-celia-et-al-form
  ;; Celia et al. (1990) and Kavetski et al. (2001, eq. 37) write the Mualem
  ;; conductivity in terms of h: K = Ks {1 - (a|h|)^(n-1) [1 + (a|h|)^n]^-m}^2
  ;; / [1 + (a|h|)^n]^(m/2). The implementation works through Se instead.
  (let ((soil (celia-sand)))
    (dolist (h '(-1d0 -10d0 -75d0 -300d0 -1000d0))
      (let* ((x (* 0.0335d0 (- h)))
             (u (+ 1d0 (expt x 2)))
             (celia (* 0.00922d0 (/ (expt (- 1d0 (* x (expt u -0.5d0))) 2) (expt u 0.25d0)))))
        (is (< (abs (- (vg-conductivity soil h) celia)) (* 1d-12 celia))
            "h ~a: ~a vs ~a" h (vg-conductivity soil h) celia)))
    (is (= 0.00922d0 (vg-conductivity soil 0d0)))
    (is (= 0.00922d0 (vg-conductivity soil 3d0)))))

(test van-genuchten-derivatives-match-finite-differences
  (dolist (soil (list (celia-sand) (van-genuchten-soil (sandy-clay-loam-parameters))))
    (dolist (h '(-0.5d0 -5d0 -75d0 -500d0 -5000d0))
      (let ((dh (* 1d-6 (abs h))))
        (flet ((fd (fn) (/ (- (funcall fn soil (+ h dh)) (funcall fn soil (- h dh))) (* 2 dh))))
          (let ((c (vg-capacity soil h)) (c-fd (fd #'vg-water-content))
                (dk (vg-conductivity-derivative soil h)) (dk-fd (fd #'vg-conductivity)))
            (is (< (abs (- c c-fd)) (* 1d-6 (abs c-fd))) "C at h ~a: ~a vs ~a" h c c-fd)
            (is (< (abs (- dk dk-fd)) (* 1d-5 (abs dk-fd))) "dK/dh at h ~a: ~a vs ~a" h dk dk-fd)
            (is (plusp c))
            (is (plusp dk))))))
    (is (= 0d0 (vg-capacity soil 0d0)) "no storage change when saturated")
    (is (= 0d0 (vg-conductivity-derivative soil 1d0)))))

(test van-genuchten-head-inverts-the-retention-curve
  (let ((soil (van-genuchten-soil (sandy-clay-loam-parameters))))
    (dolist (h '(-1d0 -50d0 -330d0 -15000d0))
      (is (< (abs (- (vg-head soil (vg-water-content soil h)) h)) (* 1d-9 (abs h)))))
    (is (= 0d0 (vg-head soil 0.39d0)) "saturated")
    (signals error (vg-head soil 0.1d0) "no finite head at residual water content")))

(test van-genuchten-soil-converts-length-units
  ;; parameter sets are in cm and s; a run in um scales alpha and Ks
  (let ((cm (celia-sand))
        (um (van-genuchten-soil (celia-sand-parameters) :length-scale 1d4)))
    (is (< (abs (- (van-genuchten-alpha um) 0.0335d-4)) 1d-18))
    (is (< (abs (- (van-genuchten-ks um) 92.2d0)) 1d-12))
    (is (< (abs (- (vg-water-content um -1d6) (vg-water-content cm -100d0))) 1d-15)
        "the same physical head gives the same water content")))

(test soil-hydraulic-parameter-sets-are-sourced
  (dolist (parameters (list (celia-sand-parameters) (sandy-clay-loam-parameters)))
    (is (null (flagged-parameters parameters)))
    (is (every #'model-parameter-unit parameters)))
  (is (= 0.00922d0 (parameter-value (celia-sand-parameters) :ks)))
  (is (= 1.48d0 (parameter-value (sandy-clay-loam-parameters) :n)))
  (signals error (make-van-genuchten :theta-r 0.3d0 :theta-s 0.2d0 :alpha 0.01d0 :n 2d0 :ks 1d0))
  (signals error (make-van-genuchten :theta-r 0.1d0 :theta-s 0.4d0 :alpha 0.01d0 :n 1d0 :ks 1d0)))
