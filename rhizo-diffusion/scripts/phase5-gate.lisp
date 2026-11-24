;;; Phase 5 gate report (PLAN.md): Richards flow.
;;;
;;; 1. Celia, Bouloutas & Zarba (1990) infiltration: 100 cm of sand,
;;;    h = -1000 cm initially, -75 cm at the top, -1000 cm at the bottom,
;;;    profiles at 1 day. Gate: mass-balance error under 0.1%.
;;;    a. reference: mixed form, dz = 0.1 cm, adaptive steps
;;;    b. mixed vs head form at dz = 1 cm and fixed steps (Celia et al.'s
;;;       comparison): mass balance ratio and wetting-front depth
;;;    c. the same column on an irregular nearest-seed Voronoi slab
;;;       (horizontal mean against the 1D reference)
;;; 2. Repeat of Phase 4 with flow-driven water content and directional
;;;    (advective) edges - see the second half of this file.
;;; Every random grid and layout is seeded; seeds are printed.
;;;
;;; Run with:
;;;   sbcl --dynamic-space-size 4096 --non-interactive --load scripts/phase5-gate.lisp

(ql:quickload "rhizo-diffusion" :silent t)

(in-package :rhizo-diffusion)

(defun celia-problem (tess &key (form :mixed))
  (make-richards-problem tess :soil (van-genuchten-soil (celia-sand-parameters))
                              :boundaries '((:top :head -75d0) (:bottom :head -1000d0))
                              :form form))

(defun column-profile (problem state)
  "List of (depth theta) of a one-cell-wide column, top to bottom."
  (let ((index (richards-problem-cell-index problem))
        (tess (richards-problem-tessellation problem)))
    (sort (loop for i below (cell-index-count index)
                collect (list (cdr (cell-position tess (index->cell-id index i)))
                              (aref (richards-state-water-contents state) i)))
          #'< :key #'first)))

(defun front-depth (profile)
  "Depth at which theta crosses halfway between its top and bottom values,
by linear interpolation in PROFILE ((depth theta) top to bottom)."
  (let* ((top (second (first profile))) (bottom (second (car (last profile))))
         (mid (* 0.5d0 (+ top bottom))))
    (loop for ((y0 t0) (y1 t1)) on profile
          while y1
          when (and (>= t0 mid) (< t1 mid))
            do (return (+ y0 (* (- y1 y0) (/ (- t0 mid) (- t0 t1))))))))

(defun profile-at (profile depth)
  "Theta of PROFILE at DEPTH by linear interpolation."
  (loop for ((y0 t0) (y1 t1)) on profile
        while y1
        when (<= y0 depth y1)
          do (return (+ t0 (* (- t1 t0) (/ (- depth y0) (- y1 y0)))))
        finally (return (second (car (last profile))))))

(defun timed (thunk)
  (let ((t0 (get-internal-real-time)))
    (multiple-value-call #'values
      (funcall thunk)
      (/ (- (get-internal-real-time) t0) internal-time-units-per-second 1.0))))

(format t "~&== 1. Celia, Bouloutas & Zarba (1990) infiltration, profiles at 1 day ==~%")
(format t "  sand: theta_r 0.102, theta_s 0.368, alpha 0.0335 /cm, n 2, Ks 0.00922 cm/s; 100 cm column~%")
(defparameter *day* 86400d0)

(defparameter *reference*
  (let* ((problem (celia-problem (make-rectangular-tessellation 1 1000 :dy 0.1d0)))
         (start (make-richards-state problem -1000d0)))
    (multiple-value-bind (end stats seconds)
        (timed (lambda () (run-richards problem start *day* :dt 0.1d0 :dt-max 60d0)))
      (let ((profile (column-profile problem end)))
        (format t "  a. reference, dz 0.1 cm, adaptive (dt <= 60 s): ~d steps, ~d rejected, ~,1f s~%"
                (getf stats :steps) (getf stats :rejected) seconds)
        (format t "     mass balance ratio ~,8f (error ~,2e); front at ~,2f cm; inflow ~,4f cm^3/cm^2~%"
                (richards-mass-balance-ratio problem start end)
                (richards-mass-balance-error problem start end)
                (front-depth profile) (richards-cumulative-inflow end))
        profile))))
(finish-output)

(defun profile-error (profile)
  "Largest |theta - theta_reference| over PROFILE's depths."
  (loop for (y th) in profile maximize (abs (- th (profile-at *reference* y)))))

(format t "  b. dz 1 cm, fixed steps (halved after a failed Newton solve, then grown back):~%")
(format t "     form   dt(s)  steps  rejected   MB ratio      MB error   front(cm)  max|dtheta| vs a.~%")
(dolist (form '(:mixed :head))
  (dolist (dt '(10d0 60d0 600d0 1800d0))
    (let* ((problem (celia-problem (make-rectangular-tessellation 1 100 :dy 1d0) :form form))
           (start (make-richards-state problem -1000d0)))
      (multiple-value-bind (end stats)
          (run-richards problem start *day* :dt dt :dt-max dt :growth 2d0
                                            :easy-iterations 1000 :hard-iterations 1000
                                            :max-iterations 50)
        (let ((profile (column-profile problem end)))
          (format t "     ~6a ~6,0f ~6d ~6d   ~11,8f  ~10,2e  ~8,2f   ~8,4f~%"
                  form dt (getf stats :steps) (getf stats :rejected)
                  (richards-mass-balance-ratio problem start end)
                  (richards-mass-balance-error problem start end)
                  (front-depth profile) (profile-error profile))))
      (finish-output))))

(let* ((problem (celia-problem (make-rectangular-tessellation 1 100 :dy 1d0)))
       (start (make-richards-state problem -1000d0)))
  (multiple-value-bind (end stats seconds)
      (timed (lambda () (run-richards problem start *day* :dt 1d0)))
    (let ((profile (column-profile problem end)))
      (format t "     mixed, adaptive (no dt cap): ~d steps (~d rejected, ~d Newton iterations, ~,2f s), final dt ~,0f s~%"
              (getf stats :steps) (getf stats :rejected) (getf stats :newton-iterations) seconds
              (getf stats :final-dt))
      (format t "       MB error ~,2e, front ~,2f cm, max|dtheta| ~,4f~%"
              (richards-mass-balance-error problem start end) (front-depth profile)
              (profile-error profile)))))
(finish-output)

(format t "  c. the same column on a nearest-seed Voronoi slab (20 cm x 100 cm, 0.5 cm pixels):~%")
(dolist (n '(400 1600 6400))
  (let* ((grid (make-nearest-seed-grid 40 200 n :seed 11))
         (tess (voronoi-grid-tessellation grid))
         ;; lengths in pixels of 0.5 cm: 2 length units per cm (heads too)
         (problem (make-richards-problem tess
                                         :soil (van-genuchten-soil (celia-sand-parameters)
                                                                   :length-scale 2d0)
                                         :boundaries '((:top :head -150d0) (:bottom :head -2000d0))))
         (start (make-richards-state problem -2000d0)))
    (multiple-value-bind (end stats seconds)
        (timed (lambda () (run-richards problem start *day* :dt 1d0 :dt-max 600d0)))
      (let* ((index (richards-problem-cell-index problem))
             (bins (make-array 50 :initial-element nil)))
        ;; horizontal mean of theta in 2 cm depth bins
        (dotimes (i (cell-index-count index))
          (let ((depth-cm (* 0.5d0 (+ 0.5d0 (cdr (cell-position tess (index->cell-id index i)))))))
            (push (aref (richards-state-water-contents end) i)
                  (aref bins (min 49 (floor depth-cm 2))))))
        (let ((profile (loop for b below 50
                             when (aref bins b)
                               collect (list (+ 1d0 (* 2 b))
                                             (/ (reduce #'+ (aref bins b)) (length (aref bins b)))))))
          (format t "     ~4d cells (seed 11): ~d steps, ~,1f s; MB error ~,2e; front ~,2f cm (1D reference ~,2f); max|dtheta| of 2 cm means ~,4f~%"
                  n (getf stats :steps) seconds
                  (richards-mass-balance-error problem start end)
                  (front-depth profile) (front-depth *reference*) (profile-error profile)))))
    (finish-output)))
