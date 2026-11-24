(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Mixed-form Richards equation on the TPFA cells (PLAN.md, Phase 5).

(defun column-heads (problem state)
  "Test helper: the heads of a one-cell-wide column, top to bottom."
  (let ((index (richards-problem-cell-index problem)))
    (loop for id from 0 below (cell-index-count index)
          collect (aref (richards-state-heads state) (cell-id->index index id)))))

(defun steady-infiltration-profile (soil h-top h-bottom depth points)
  "Test helper, independent of the finite-volume code: the steady 1D head
profile between heads H-TOP at y = 0 and H-BOTTOM at y = DEPTH (y down),
from Darcy's law with constant downward flux q, dh/dy = 1 - q / K(h),
integrated by RK4 with q found by bisection. Returns (values q heads), the
heads at the y values in POINTS."
  ;; trial fluxes that are too large drive h toward -infinity, where K
  ;; underflows: clamp h at -1e5 cm, below any profile of interest
  (labels ((slope (q h) (- 1d0 (/ q (vg-conductivity soil (max h -1d5)))))
           (integrate (q y-end)
             (let* ((steps 4000) (dy (/ y-end steps)) (h h-top))
               (dotimes (k steps h)
                 (let* ((k1 (slope q h)) (k2 (slope q (+ h (* 0.5d0 dy k1))))
                        (k3 (slope q (+ h (* 0.5d0 dy k2)))) (k4 (slope q (+ h (* dy k3)))))
                   (setf h (+ h (* (/ dy 6d0) (+ k1 (* 2 k2) (* 2 k3) k4))))
                   (when (< h -1d5) (return -1d5)))))))
    (let ((lo 0d0) (hi (van-genuchten-ks soil)))
      (loop repeat 60 ; 2^-60: below double precision of q
            do (let ((mid (* 0.5d0 (+ lo hi))))
                 (if (> (integrate mid depth) h-bottom) (setf lo mid) (setf hi mid))))
      (let ((q (* 0.5d0 (+ lo hi))))
        (values q (mapcar (lambda (y) (integrate q y)) points))))))

(defun celia-column (n-cells &key (form :mixed) (face-mean :arithmetic))
  "Test helper: the Celia et al. (1990) problem on N-CELLS cells over 100 cm."
  (make-richards-problem (make-rectangular-tessellation 1 n-cells :dx 1d0 :dy (/ 100d0 n-cells))
                         :soil (van-genuchten-soil (celia-sand-parameters))
                         :boundaries '((:top :head -75d0) (:bottom :head -1000d0))
                         :form form :face-mean face-mean))

(test hydrostatic-column-stays-at-rest
  ;; H = h + z constant: no flux anywhere, whatever the step
  (let* ((tess (make-rectangular-tessellation 1 20 :dy 5d0))
         (problem (make-richards-problem tess :soil (van-genuchten-soil (celia-sand-parameters))
                                              :boundaries '((:bottom :head -20d0))))
         (state (make-richards-state problem (lambda (id) (- -20d0 (- 100d0 (* 5d0 (+ id 0.5d0))))))))
    (let ((next (run-richards problem state 3600d0 :dt 3600d0)))
      (is (every (lambda (a b) (< (abs (- a b)) 1d-9))
                 (richards-state-heads state) (richards-state-heads next)))
      (is (< (abs (richards-mass-balance-error problem state next :absolute t)) 1d-12)))))

(test unit-gradient-drainage-is-steady
  ;; uniform h with inflow K(h) at the top and free drainage at the bottom
  (let* ((soil (van-genuchten-soil (sandy-clay-loam-parameters)))
         (h -80d0)
         (problem (make-richards-problem (make-rectangular-tessellation 1 10 :dy 2d0)
                                         :soil soil
                                         :boundaries `((:top :flux ,(vg-conductivity soil h))
                                                       (:bottom :free-drainage))))
         (state (run-richards problem (make-richards-state problem h) 86400d0 :dt 100d0)))
    (is (every (lambda (x) (< (abs (- x h)) 1d-6)) (richards-state-heads state)))))

(test steady-infiltration-matches-darcy-quadrature
  ;; Independent reference: Darcy's law integrated as an ODE, not the FV
  ;; equations. A smooth profile (wetter at depth, 50 cm of Celia sand):
  ;; the arithmetic face mean converges at second order, the upstream one
  ;; at first. (A drier bottom with a large flux makes a boundary layer
  ;; thinner than the cells - measured 2026-11-24: -50/-200 cm gives 8.7,
  ;; 4.4, 2.2, 1.0% at 25-200 cells, all in the last cell, where the
  ;; gradient is about -325 cm/cm. That case is pre-asymptotic, not wrong.)
  (let ((soil (van-genuchten-soil (celia-sand-parameters))))
    (flet ((max-error (n face-mean)
             (let* ((dy (/ 50d0 n))
                    (problem (make-richards-problem (make-rectangular-tessellation 1 n :dy dy)
                                                    :soil soil :face-mean face-mean
                                                    :boundaries '((:top :head -80d0)
                                                                  (:bottom :head -60d0))))
                    (state (run-richards problem (make-richards-state problem -60d0) 3d6
                                         :dt 10d0 :dt-max 1d5))
                    (reference (nth-value 1 (steady-infiltration-profile
                                             soil -80d0 -60d0 50d0
                                             (loop for j below n collect (* dy (+ j 0.5d0)))))))
               (loop for h in (column-heads problem state) for r in reference
                     maximize (abs (/ (- h r) r))))))
      (let ((a25 (max-error 25 :arithmetic)) (a50 (max-error 50 :arithmetic))
            (u25 (max-error 25 :upstream)) (u50 (max-error 50 :upstream)))
        (is (< a50 1d-4) "arithmetic, 50 cells (1 cm): max relative head error ~a" a50)
        (is (< 3.5d0 (/ a25 a50) 4.5d0) "arithmetic: second order, ~a -> ~a" a25 a50)
        (is (< 1.7d0 (/ u25 u50) 2.3d0) "upstream: first order, ~a -> ~a" u25 u50)))))

(test richards-jacobian-matches-finite-differences
  (let* ((tess (make-rectangular-tessellation 3 3 :dx 2d0 :dy 1.5d0))
         (soil (van-genuchten-soil (sandy-clay-loam-parameters)))
         (h-old (make-array 9 :element-type 'double-float))
         (h (make-array 9 :element-type 'double-float)))
    (dotimes (i 9)
      (setf (aref h-old i) (- -40d0 (* 13d0 i))
            (aref h i) (- -35d0 (* 11d0 (mod (* 4 i) 9)))))
    (dolist (face-mean '(:arithmetic :upstream))
      (dolist (form '(:mixed :head))
        (let* ((problem (make-richards-problem
                         tess :soil soil :face-mean face-mean :form form
                              :boundaries '((:top :atmospheric -1d-4 -1000d0)
                                            (:bottom :free-drainage)
                                            (:left :head -30d0)
                                            (:right :flux 2d-5))))
               (theta-old (map '(simple-array double-float (*))
                               (lambda (x) (vg-water-content soil x)) h-old)))
          (multiple-value-bind (r0 diagonal off) (richards-system problem h h-old theta-old 50d0 0d0)
            (declare (ignore r0))
            (dotimes (j 9)
              (let* ((e (make-array 9 :element-type 'double-float :initial-element 0d0))
                     (eps (* 1d-6 (abs (aref h j)))))
                (setf (aref e j) 1d0)
                (let ((column (apply-csr-matrix (richards-problem-stencil problem) diagonal off e))
                      (hp (copy-seq h)) (hm (copy-seq h)))
                  (incf (aref hp j) eps) (decf (aref hm j) eps)
                  (let ((rp (richards-system problem hp h-old theta-old 50d0 0d0))
                        (rm (richards-system problem hm h-old theta-old 50d0 0d0)))
                    (dotimes (i 9)
                      (let ((fd (/ (- (aref rp i) (aref rm i)) (* 2 eps)))
                            (an (aref column i)))
                        (unless (eq form :head) ; the head form drops dC/dh (quasi-Newton)
                          (is (< (abs (- fd an)) (+ 1d-9 (* 1d-5 (abs fd))))
                              "~a ~a dR~a/dh~a: analytic ~a, finite difference ~a"
                              face-mean form i j an fd))))))))))))))

(test closed-column-takes-in-exactly-the-imposed-flux
  (let* ((problem (make-richards-problem (make-rectangular-tessellation 1 25 :dy 2d0)
                                         :soil (van-genuchten-soil (celia-sand-parameters))
                                         :boundaries '((:top :flux 1d-3))))
         (start (make-richards-state problem -500d0))
         (end (run-richards problem start 7200d0 :dt 10d0)))
    (is (< (abs (- (richards-storage problem end) (richards-storage problem start) 7.2d0)) 1d-8)
        "1e-3 cm/s over 1 cm for 2 h adds 7.2 cm^2")
    (is (< (abs (- (richards-cumulative-inflow end) 7.2d0)) 1d-10))))

(test atmospheric-boundary-limits-evaporation-from-dry-soil
  ;; demand 1e-4 cm/s; a dry sand cannot supply it with the surface at -1e4 cm
  (let* ((soil (van-genuchten-soil (celia-sand-parameters)))
         (problem (make-richards-problem (make-rectangular-tessellation 1 20 :dy 1d0)
                                         :soil soil
                                         :boundaries '((:top :atmospheric -1d-4 -1d4))))
         (wet (run-richards problem (make-richards-state problem -20d0) 600d0 :dt 10d0))
         (dry (run-richards problem (make-richards-state problem -3000d0) 600d0 :dt 10d0)))
    (is (< (abs (- (richards-cumulative-inflow wet) -0.06d0)) 1d-9) "wet soil: at the demand")
    (is (> (richards-cumulative-inflow dry) -0.06d0) "dry soil: below the demand")
    (is (minusp (richards-cumulative-inflow dry)))))

(test celia-infiltration-conserves-mass-in-mixed-form-only
  ;; Celia et al. (1990): the mixed form conserves mass at any step; the
  ;; head form does not. Coarse version of the gate (scripts/phase5-gate.lisp).
  (let* ((mixed (celia-column 50))
         (head (celia-column 50 :form :head))
         (start-m (make-richards-state mixed -1000d0))
         (start-h (make-richards-state head -1000d0))
         (end-m (run-richards mixed start-m 3600d0 :dt 10d0 :dt-max 600d0))
         ;; fixed 10 min steps: no growth, and no shrinking for slow
         ;; (quasi-Newton) convergence - only after an outright failure
         (end-h (run-richards head start-h 3600d0 :dt 600d0 :growth 1d0
                                                  :hard-iterations 1000 :max-iterations 50)))
    (is (< (abs (richards-mass-balance-error mixed start-m end-m)) 1d-6)
        "mixed form: ~a" (richards-mass-balance-error mixed start-m end-m))
    (is (> (abs (richards-mass-balance-error head start-h end-h)) 1d-3)
        "head form at 10 min steps: ~a" (richards-mass-balance-error head start-h end-h))
    (is (plusp (richards-cumulative-inflow end-m)) "water enters at the top")))

(test run-richards-adapts-the-step-and-lands-on-the-end-time
  (let* ((problem (celia-column 25))
         (start (make-richards-state problem -1000d0)))
    (multiple-value-bind (end stats) (run-richards problem start 1800d0 :dt 1d0 :dt-max 300d0)
      (is (= 1800d0 (richards-state-time end)))
      (is (> (richards-state-next-dt end) 1d0) "the step grew")
      (is (< (getf stats :steps) 1800))
      (is (plusp (getf stats :newton-iterations))))
    ;; a first step far too large is rejected and retried smaller
    (multiple-value-bind (end stats) (run-richards problem start 600d0 :dt 1d6 :max-iterations 6)
      (is (= 600d0 (richards-state-time end)))
      (is (plusp (getf stats :rejected))))))

(test richards-boundaries-take-their-elevation-from-the-face
  (let* ((problem (make-richards-problem (make-rectangular-tessellation 1 4 :dy 10d0)
                                         :soil (van-genuchten-soil (celia-sand-parameters))
                                         :boundaries '((:top :head -1d0) (:bottom :head -2d0))))
         (faces (richards-problem-boundaries problem)))
    ;; y down: the top face is at y = 0 (elevation 0), the bottom at y = 40
    (is (every #'= '(0d0 -40d0) (sort (mapcar #'richards-boundary-elevation faces) #'>)))
    (is (every (lambda (f) (= 5d0 (richards-boundary-distance f))) faces))))
