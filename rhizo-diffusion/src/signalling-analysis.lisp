(in-package :rhizo-diffusion)

;;; Analysis of signalling runs (PLAN.md, Phase 2 gate): analytic
;;; references and the measurements compared with them.
;;;
;;; Steady state around a single producer. A source of strength S (amount
;;; per unit depth per second) in the plane, with diffusivity D and
;;; first-order loss k, settles to
;;;
;;;     C(r) = S / (2 pi D) K0(r / lambda),   lambda = sqrt(D / k),
;;;
;;; the 2D fundamental solution of D lap C - k C = -S delta. On the slab
;;; grid of signalling.lisp (uniform in depth) a producer patch is a source
;;; of this kind. Sorption that only exchanges with the dissolved pool does
;;; not change it: at steady state the sorbed pool is at equilibrium.

(defun bessel-k0 (x)
  "Return the modified Bessel function of the second kind K0(X), X > 0,
from K0(x) = integral over t in [0, inf) of exp(-x cosh t), by the
trapezoidal rule with step 1/8. The integrand is analytic in the strip
|Im t| < pi/2 and decays doubly exponentially, so the rule's error is about
exp(-8 pi^2), far below double precision. Signals an error unless X > 0."
  (unless (plusp x)
    (error "BESSEL-K0: argument must be positive, got ~a" x))
  (let ((x (float x 1d0)) (h 0.125d0))
    (* h (+ (* 0.5d0 (exp (- x)))
            (loop for i from 1
                  for term = (exp (- (* x (cosh (* i h)))))
                  sum term into total of-type double-float
                  until (< term (* 1d-18 (+ total (exp (- x)))))
                  finally (return total))))))

(defun point-source-steady-state (source-rate diffusivity decay-rate r)
  "Return the steady concentration at distance R from a source of
SOURCE-RATE (amount per unit depth per unit time) in the plane, with
DIFFUSIVITY and first-order DECAY-RATE: S / (2 pi D) K0(r / lambda), lambda
= sqrt(D / k). See the file header."
  (* (/ source-rate (* 2 pi diffusivity))
     (bessel-k0 (/ r (sqrt (/ diffusivity decay-rate))))))
