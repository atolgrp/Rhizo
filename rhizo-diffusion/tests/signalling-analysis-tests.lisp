(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Analytic references and measurements for the Phase 2 gate.

(test bessel-k0-matches-tabulated-values
  ;; Abramowitz & Stegun, Table 9.8 / standard tables
  (loop for (x k0) in '((0.1d0 2.4270690247020166d0)
                        (0.5d0 0.9244190712276659d0)
                        (1d0 0.42102443824070834d0)
                        (2d0 0.11389387274953344d0)
                        (5d0 0.0036910983340425942d0)
                        (10d0 1.778006231616917d-5))
        do (is (< (relative-difference (bessel-k0 x) k0) 1d-12)
               "K0(~a): expected ~a, got ~a" x k0 (bessel-k0 x)))
  ;; small-argument asymptote -ln(x/2) - Euler's gamma, error O(x^2 ln x)
  (is (< (abs (- (bessel-k0 1d-4) (- (- (log 5d-5)) 0.5772156649015329d0))) 1d-7))
  (signals error (bessel-k0 0d0)))

(test point-source-steady-state-is-the-k0-profile
  ;; S / (2 pi D) K0(r / lambda), lambda = sqrt(D / k): with S = 2 pi D and
  ;; r = lambda the value is K0(1)
  (let ((d 490d0) (k 0.5d0))
    (is (< (relative-difference (point-source-steady-state (* 2 pi d) d k (sqrt (/ d k)))
                                (bessel-k0 1d0))
           1d-14))
    (is (< (relative-difference (point-source-steady-state 3d0 d k 10d0)
                                (* (/ 3d0 (* 2 pi d)) (bessel-k0 (/ 10d0 (sqrt (/ d k))))))
           1d-14))))
