(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Quorum induction of a uniform population (PLAN.md, Phase 2 gate)

(test uniform-ahl-balance-is-production-minus-losses
  (let* ((p (ahl-literature-parameters))
         (rho 1d9) (rho-d 1d8) (c 35d0)
         (v (lambda (name) (parameter-value p name))))
    (is (< (relative-difference
            (uniform-ahl-balance p rho c :degrader-density rho-d)
            (- (* rho (+ (funcall v :basal-production)
                         (* (funcall v :induced-production) (hill-activation c 70d0 2.5d0))))
               (* (funcall v :hydrolysis-rate) c)
               (/ (* (funcall v :biodegradation-kcat) rho-d c) (+ (funcall v :biodegradation-km) c))))
           1d-12))))

(test uniform-steady-states-are-roots-of-the-balance
  (let ((p (ahl-literature-parameters)))
    (dolist (rho '(1d7 1d8 2.5d8 1d9 1d11))
      (let ((roots (uniform-steady-states p rho))
            (scale (* rho (parameter-value p :basal-production))))
        (is (member (length roots) '(1 3)) "~a steady states at density ~a" (length roots) rho)
        (is (equal roots (sort (copy-list roots) #'<)))
        (dolist (c roots)
          (is (< (abs (uniform-ahl-balance p rho c)) (* 1d-9 scale))))))))

(test quorum-switch-density-brackets-the-bistable-range
  (let* ((p (ahl-literature-parameters))
         (switch (quorum-switch-density p)))
    (is (eq :bistable (getf switch :kind)))
    (is (< (getf switch :off) (getf switch :on)))
    (is (= 3 (length (uniform-steady-states p (* 0.99d0 (getf switch :on))))))
    (is (= 1 (length (uniform-steady-states p (* 1.01d0 (getf switch :on))))))
    (is (= 1 (length (uniform-steady-states p (* 0.99d0 (getf switch :off))))))
    ;; With n = 1 the feedback cannot make a fold: the switch is graded, and
    ;; its density is where the steady AHL concentration is the threshold.
    (let* ((p1 (ahl-test-parameters :hill-coefficient 1d0))
           (graded (quorum-switch-density p1)))
      (is (eq :graded (getf graded :kind)))
      (is (< (relative-difference (first (uniform-steady-states p1 (getf graded :half))) 70d0)
             1d-9)))))

(test well-mixed-tessellation-is-one-cell-without-faces
  (let ((tess (well-mixed-tessellation :area 2.5d0 :id "pot")))
    (is (equal '("pot") (cell-ids tess)))
    (is (null (neighbors-of tess "pot")))
    (is (= 2.5d0 (cell-area tess "pot")))))

(test golden-section-maximum-finds-the-peak
  (is (< (abs (- 2d0 (golden-section-maximum (lambda (x) (- (expt (- x 2d0) 2))) 0d0 5d0))) 1d-8))
  ;; near a maximum f differs from its peak by the square of the distance,
  ;; so the location resolves only to about sqrt(double epsilon) = 1.5e-8
  (is (< (abs (- (/ pi 2) (golden-section-maximum #'sin 0d0 3d0))) 1d-7)))

(test quorum-induction-switches-at-the-predicted-density
  ;; PLAN.md Phase 2 gate: a uniform population has no gradients, so one
  ;; well-mixed cell is the whole domain. Literature parameters, no sorption.
  (let* ((p (ahl-test-parameters :kd 0d0))
         (on (getf (quorum-switch-density p) :on))
         (below (uniform-steady-ahl p (* 0.97d0 on)))
         (above (uniform-steady-ahl p (* 1.03d0 on))))
    (is (< below 70d0) "3% below the predicted density AHL reached ~a nM" below)
    (is (> above 70d0) "3% above the predicted density AHL reached only ~a nM" above)
    (is (< (relative-difference below (first (uniform-steady-states p (* 0.97d0 on)))) 1d-6))
    (is (< (relative-difference above (first (last (uniform-steady-states p (* 1.03d0 on))))) 1d-6))))

(test quorum-switch-is-unmoved-by-sorption
  ;; Sorption (K' = 4) slows the approach fivefold but cannot move the fold:
  ;; 1% either side of the predicted density.
  (let* ((p (ahl-test-parameters :kd (/ (* 4d0 0.3d0) 1.3d0)))
         (on (getf (quorum-switch-density p) :on)))
    (is (null (uniform-induction-p p (* 0.99d0 on))))
    (is (uniform-induction-p p (* 1.01d0 on)))))

(test simulated-switch-density-finds-the-fold
  (let* ((p (ahl-test-parameters :kd 0d0))
         (predicted (getf (quorum-switch-density p) :on))
         (simulated (simulated-switch-density p :lo (* 0.9d0 predicted) :hi (* 1.1d0 predicted)
                                                :tolerance 5d-3)))
    (is (< (relative-difference simulated predicted) 1d-2)
        "simulated switch at ~a cells/L, predicted ~a" simulated predicted)))

(test uniform-steady-ahl-refuses-an-unconverged-run
  (signals error (uniform-steady-ahl (ahl-literature-parameters) 1d9 :max-steps 2)))

(test uniform-induction-p-classifies-either-side-of-the-fold
  (let* ((p (ahl-test-parameters :kd 0d0))
         (on (getf (quorum-switch-density p) :on)))
    (is (null (uniform-induction-p p (* 0.97d0 on))))
    (is (uniform-induction-p p (* 1.03d0 on)))))
