(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Soil water and solute diffusion (PLAN.md, Phase 3).

(test millington-quirk-relative-diffusivity
  (let ((xi (millington-quirk 0.5d0)))
    (is (< (abs (- (funcall xi 0.3d0) (/ (expt 0.3d0 10/3) 0.25d0))) 1d-15))
    (is (< (abs (- (funcall xi 0.5d0) (expt 0.5d0 4/3))) 1d-15) "saturated: phi^(4/3)")
    (is (= 0d0 (funcall xi 0d0)))))

(test papendick-campbell-relative-diffusivity
  (is (< (abs (- (funcall (papendick-campbell 2.8d0 3d0) 0.3d0) (* 2.8d0 0.027d0))) 1d-15)))

(test mullins-sommer-relative-diffusivity-has-a-threshold
  ;; delta theta^2 + mu theta, zero at and below theta_th = -mu / delta
  (let ((xi (mullins-sommer 2.1679d0 -0.5431d0)))
    (is (< (abs (- (funcall xi 0.4d0) (+ (* 2.1679d0 0.16d0) (* -0.5431d0 0.4d0)))) 1d-15))
    (is (= 0d0 (funcall xi 0.2d0)) "below the threshold water films do not connect")
    (is (= 0d0 (funcall xi 0.25d0)))
    (is (< (abs (- (mullins-sommer-threshold 2.1679d0 -0.5431d0) 0.2505d0)) 1d-4)
        "the threshold Chou et al. (2012) quote for their sandy clay loam")))

(test soil-parameters-carry-the-diffusion-model-constants
  (let ((parameters (soil-parameters)))
    (is (every #'model-parameter-unit parameters))
    (is (equal '(:bulk-density :porosity) (flagged-parameters parameters)))
    (is (= 2.1679d0 (parameter-value parameters :mullins-sommer-delta)))
    (is (= -0.5431d0 (parameter-value parameters :mullins-sommer-mu)))
    (is (= 2.8d0 (parameter-value parameters :papendick-campbell-coefficient)))))

(test relative-diffusivity-model-builds-each-relation-from-parameters
  (let ((parameters (soil-parameters)))
    (is (< (abs (- (funcall (relative-diffusivity-model :millington-quirk parameters) 0.3d0)
                   (funcall (millington-quirk (parameter-value parameters :porosity)) 0.3d0)))
           1d-15))
    (is (< (abs (- (funcall (relative-diffusivity-model :papendick-campbell parameters) 0.3d0)
                   (* 2.8d0 (expt 0.3d0 3))))
           1d-15))
    (is (= 0d0 (funcall (relative-diffusivity-model :mullins-sommer parameters) 0.2d0)))
    (signals error (relative-diffusivity-model :no-such-model parameters))))
