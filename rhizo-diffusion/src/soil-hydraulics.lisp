(in-package :rhizo-diffusion)

;;; Soil hydraulic relations for Richards flow (PLAN.md, Phase 5).
;;;
;;; Van Genuchten (1980) retention with the Mualem (1976) conductivity model
;;; and pore-connectivity exponent l = 1/2, the form Celia, Bouloutas &
;;; Zarba (1990) use for their infiltration benchmark. With x = alpha |h|,
;;; u = 1 + x^n and m = 1 - 1/n, for pressure head h < 0:
;;;
;;;     Se    = u^-m                         effective saturation
;;;     theta = theta_r + (theta_s - theta_r) Se
;;;     K     = Ks u^(-m/2) (1 - x^(n-1) u^-m)^2
;;;
;;; and theta = theta_s, K = Ks for h >= 0. C = dtheta/dh and dK/dh are
;;; written in x as well, derived by hand and checked against finite
;;; differences in the tests. For n < 2, dK/dh grows without bound as h -> 0
;;; (a known property of the Mualem form), although K itself stays finite.
;;;
;;; Assumed, not measured: one retention curve (no hysteresis between
;;; wetting and drying), and a rigid soil.
;;;
;;; Parameter sets are in cm and s, as published. VAN-GENUCHTEN-SOIL converts
;;; them to the length unit of a run (the signalling worlds use um).

(defstruct (van-genuchten (:constructor %make-van-genuchten (theta-r theta-s alpha n ks)))
  (theta-r 0d0 :type double-float)  ; residual water content
  (theta-s 0d0 :type double-float)  ; saturated water content
  (alpha 0d0 :type double-float)    ; 1 / length
  (n 0d0 :type double-float)        ; > 1
  (ks 0d0 :type double-float))      ; saturated conductivity, length / time

(defun make-van-genuchten (&key theta-r theta-s alpha n ks)
  "Return a VAN-GENUCHTEN soil with these parameters (as double-floats).
Signals an error unless 0 <= THETA-R < THETA-S <= 1, ALPHA > 0, N > 1 and
KS > 0."
  (unless (and (<= 0 theta-r) (< theta-r theta-s) (<= theta-s 1)
               (plusp alpha) (> n 1) (plusp ks))
    (error "MAKE-VAN-GENUCHTEN: invalid parameters theta_r ~a theta_s ~a alpha ~a n ~a Ks ~a"
           theta-r theta-s alpha n ks))
  (%make-van-genuchten (float theta-r 1d0) (float theta-s 1d0) (float alpha 1d0)
                       (float n 1d0) (float ks 1d0)))

(declaim (inline vg-m))
(defun vg-m (soil)
  "Return the van Genuchten exponent m = 1 - 1/n of SOIL."
  (- 1d0 (/ 1d0 (van-genuchten-n soil))))

(defun vg-effective-saturation (soil h)
  "Return the effective saturation Se = [1 + (alpha |H|)^n]^-m of SOIL at
pressure head H; 1 for H >= 0."
  (if (>= h 0)
      1d0
      (expt (+ 1d0 (expt (* (van-genuchten-alpha soil) (- h)) (van-genuchten-n soil)))
            (- (vg-m soil)))))

(defun vg-water-content (soil h)
  "Return the volumetric water content of SOIL at pressure head H."
  (+ (van-genuchten-theta-r soil)
     (* (- (van-genuchten-theta-s soil) (van-genuchten-theta-r soil))
        (vg-effective-saturation soil h))))

(defun vg-capacity (soil h)
  "Return the specific moisture capacity C = dtheta/dh of SOIL at pressure
head H: (theta_s - theta_r) alpha m n x^(n-1) u^(-m-1); 0 for H >= 0."
  (if (>= h 0)
      0d0
      (let* ((alpha (van-genuchten-alpha soil)) (n (van-genuchten-n soil)) (m (vg-m soil))
             (x (* alpha (- h)))
             (u (+ 1d0 (expt x n))))
        (* (- (van-genuchten-theta-s soil) (van-genuchten-theta-r soil))
           alpha m n (expt x (- n 1d0)) (expt u (- -1d0 m))))))

(defun vg-conductivity (soil h)
  "Return the Mualem hydraulic conductivity of SOIL at pressure head H:
Ks u^(-m/2) (1 - x^(n-1) u^-m)^2; Ks for H >= 0."
  (if (>= h 0)
      (van-genuchten-ks soil)
      (let* ((n (van-genuchten-n soil)) (m (vg-m soil))
             (x (* (van-genuchten-alpha soil) (- h)))
             (u (+ 1d0 (expt x n)))
             (f (- 1d0 (* (expt x (- n 1d0)) (expt u (- m))))))
        (* (van-genuchten-ks soil) (expt u (* -0.5d0 m)) f f))))

(defun vg-conductivity-derivative (soil h)
  "Return dK/dh of SOIL at pressure head H (see the file header); 0 for
H >= 0."
  (if (>= h 0)
      0d0
      (let* ((alpha (van-genuchten-alpha soil)) (n (van-genuchten-n soil)) (m (vg-m soil))
             (x (* alpha (- h)))
             (u (+ 1d0 (expt x n)))
             (f (- 1d0 (* (expt x (- n 1d0)) (expt u (- m)))))
             (df/dx (- (* (- n 1d0) (expt x (- n 2d0)) (expt u (- -1d0 m)))))
             (dk/dx (* (van-genuchten-ks soil) (expt u (* -0.5d0 m))
                       (+ (/ (* -0.5d0 m n (expt x (- n 1d0)) f f) u)
                          (* 2d0 f df/dx)))))
        ;; x = alpha |h| = -alpha h
        (* (- alpha) dk/dx))))

(defun vg-head (soil theta)
  "Return the pressure head at which SOIL holds water content THETA: 0 at
or above theta_s. Signals an error at or below theta_r, where no finite
head gives THETA."
  (let ((theta-r (van-genuchten-theta-r soil)) (theta-s (van-genuchten-theta-s soil)))
    (cond ((>= theta theta-s) 0d0)
          ((<= theta theta-r)
           (error "VG-HEAD: water content ~a is at or below theta_r ~a" theta theta-r))
          (t (let ((se (/ (- theta theta-r) (- theta-s theta-r))))
               (- (/ (expt (- (expt se (/ -1d0 (vg-m soil))) 1d0) (/ 1d0 (van-genuchten-n soil)))
                     (van-genuchten-alpha soil))))))))

;;; Parameter sets (cm, s).

(defun celia-sand-parameters ()
  "Return the van Genuchten-Mualem parameters of the sand in the Celia,
Bouloutas & Zarba (1990) infiltration benchmark, in cm and s."
  (let ((source "Celia, Bouloutas & Zarba 1990, Water Resour Res 26:1483 (benchmark soil); also Kavetski et al. 2001, Adv Water Resour 24:595"))
    (list (make-model-parameter :theta-r 0.102d0 "-" :source source)
          (make-model-parameter :theta-s 0.368d0 "-" :source source)
          (make-model-parameter :alpha 0.0335d0 "1/cm" :source source)
          (make-model-parameter :n 2d0 "-" :source source)
          (make-model-parameter :ks 0.00922d0 "cm/s" :source source))))

(defun sandy-clay-loam-parameters ()
  "Return mean van Genuchten-Mualem parameters for a sandy clay loam, in cm
and s: the texture whose diffusion constants Phase 3 uses (Chou et al.
2012)."
  (let ((source "Carsel & Parrish 1988, Water Resour Res 24:755, Table 3 class means (sandy clay loam)"))
    (list (make-model-parameter :theta-r 0.100d0 "-" :source source)
          (make-model-parameter :theta-s 0.39d0 "-" :source source)
          (make-model-parameter :alpha 0.059d0 "1/cm" :source source)
          (make-model-parameter :n 1.48d0 "-" :source source)
          (make-model-parameter :ks (/ 31.44d0 86400d0) "cm/s"
                                :source source
                                :note "31.44 cm/day as reproduced in the HYDRUS-1D soil catalogue; theta_r, theta_s, alpha and n were checked against a second secondary source, Ks was not"))))

(defun van-genuchten-soil (parameters &key (length-scale 1d0))
  "Return the VAN-GENUCHTEN soil of PARAMETERS (a set in cm and s, such as
CELIA-SAND-PARAMETERS) in a length unit LENGTH-SCALE times smaller than a
cm (1d4 for um): alpha / LENGTH-SCALE, Ks x LENGTH-SCALE."
  (flet ((p (name) (parameter-value parameters name)))
    (make-van-genuchten :theta-r (p :theta-r) :theta-s (p :theta-s)
                        :alpha (/ (p :alpha) length-scale) :n (p :n)
                        :ks (* (p :ks) length-scale))))
