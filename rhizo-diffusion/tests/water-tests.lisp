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

(test no-tortuosity-relative-diffusivity-is-the-water-content
  ;; D_s = theta D_0: pore-water diffusivity equals that in free water
  (is (= 0.3d0 (funcall (relative-diffusivity-model :no-tortuosity (soil-parameters)) 0.3d0))))

(defun water-test-world ()
  "Two cells; :C dissolved (free-water D 2), :S immobile, :WATER the water
content."
  (make-mechanistic-world
   (two-cell-tessellation 1d0 1d0 1d0)
   :species (list (make-species :name :c :diffusivity 2d0) (make-species :name :s)
                  (water-species))
   :reactions '() :dt 1d0 :diffusion-scheme :imex
   :initial-state '((:c . 3d0) (:s . 1d0))))

(test set-water-content-sets-water-capacity-and-diffusivity
  (let ((world (water-test-world))
        (theta (lambda (id) (if (string= id "0") 0.2d0 0.4d0)))
        (xi (lambda (th) (* th th))))
    (set-water-content! world theta xi)
    (let ((index (mechanistic-world-cell-index world))
          (species-index (mechanistic-world-species-index world)))
      (flet ((at (array id name) (aref array (cell-id->index index id) (gethash name species-index))))
        (is (= 0.4d0 (at (mechanistic-world-state world) "1" :water)))
        (is (= 0.2d0 (at (mechanistic-world-capacities world) "0" :c)) "dissolved: capacity theta")
        (is (= 1d0 (at (mechanistic-world-capacities world) "0" :s)) "immobile: unchanged")
        (is (< (abs (- (at (mechanistic-world-diffusivities world) "1" :c) (* 2d0 0.16d0))) 1d-15)
            "D_0 xi(theta)")
        (is (= 3d0 (at (mechanistic-world-state world) "0" :c)) "amounts are not touched")))))

(test set-water-content-keeps-blocked-cells-walled
  (let ((world (water-test-world)))
    (set-water-content! world (constantly 0.3d0) #'identity :blocked-p (lambda (id) (string= id "0")))
    (is (= 0d0 (aref (mechanistic-world-diffusivities world)
                     (cell-id->index (mechanistic-world-cell-index world) "0")
                     (gethash :c (mechanistic-world-species-index world)))))
    (signals error (set-water-content! world (constantly 0d0) #'identity))))

(test drying-concentrates-solute-and-conserves-it
  ;; Halving theta leaves every amount in place, so the pore-water
  ;; concentration m / theta doubles and the total is unchanged.
  (let ((world (water-test-world)))
    (set-water-content! world (constantly 0.4d0) #'identity)
    (let ((mass (species-mass world :c)))
      (set-water-content! world (constantly 0.2d0) #'identity)
      (is (= mass (species-mass world :c)))
      (is (< (abs (- 15d0 (concentration-of (let ((row (make-array 3 :element-type 'double-float)))
                                              (dotimes (s 3 row)
                                                (setf (aref row s) (aref (mechanistic-world-state world) 0 s))))
                                            (mechanistic-world-species-index world)
                                            '(:per :c :water))))
             1d-12)))))
