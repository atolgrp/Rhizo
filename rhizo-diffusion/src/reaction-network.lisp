(in-package :rhizo-diffusion)

;;; A small, reusable chemical-reaction-network (CRN) engine: SPECIES hold a
;;; name and a diffusivity, REACTIONs hold a stoichiometry and a rate law,
;;; and NET-REACTION-DERIVATIVE! sums rate*stoichiometry across reactions to
;;; get each species' local (non-spatial) rate of change - the standard CRN
;;; ODE right-hand side, nothing bespoke.
;;;
;;; Deliberately knows nothing about PFAS, GIS, or rhizospheres - a cell's
;;; state is just a vector of numbers keyed by species name via a shared
;;; SPECIES-INDEX. Concrete species and reactions are defined elsewhere
;;; (examples/pfas-rhizobiome.lisp); a different project reuses this file
;;; unchanged by defining a different species/reaction list.

(defstruct species
  name         ; a keyword, e.g. :substrate - the key used everywhere below
  (diffusivity 0.0d0 :type double-float)) ; 0.0d0 for immobile/sorbed species

(defstruct reaction
  name
  stoichiometry ; alist of (species-name . signed-coefficient) - negative
                ; entries are consumed, positive ones produced
  rate-fn)      ; (function (state species-index) -> rate). Nonnegative for
                ; an irreversible reaction (Monod growth, decay); signed for
                ; a reversible one (linear equilibrium partitioning) - the
                ; sign combined with STOICHIOMETRY determines each species'
                ; contribution, exactly like a standard rate law.

(defun make-species-index (species-list)
  "EQ hash-table: species-name -> dense array index, in the order
SPECIES-LIST lists them. Shared by every cell-state vector and every
reaction's rate-fn, so a species name always means the same array slot."
  (let ((index (make-hash-table :test 'eq :size (length species-list))))
    (loop for sp in species-list
          for i from 0
          do (setf (gethash (species-name sp) index) i))
    index))

(declaim (inline state-value))
(defun state-value (state species-index name)
  (aref state (or (gethash name species-index)
                   (error "STATE-VALUE: unknown species ~s" name))))

(defun concentration-of (state species-index reference)
  "Return the value a rate law reads for REFERENCE in STATE: for a species
name, its STATE-VALUE; for (:PER name water-name), name's value divided by
water-name's - an amount per soil volume over the water content, i.e. a
pore-water concentration (PLAN.md, Phase 3). Signals an error for anything
else."
  (cond ((symbolp reference) (state-value state species-index reference))
        ((and (consp reference) (eq (first reference) :per) (= (length reference) 3))
         (/ (state-value state species-index (second reference))
            (state-value state species-index (third reference))))
        (t (error "CONCENTRATION-OF: not a species reference: ~s" reference))))

(defun make-cell-state (species-index &optional initial-values)
  "A fresh (SIMPLE-ARRAY DOUBLE-FLOAT (n)) cell-state vector, one slot per
species in SPECIES-INDEX, all zero except any overridden by INITIAL-VALUES
(an alist of species-name -> value)."
  (let* ((n (hash-table-count species-index))
         (state (make-array n :element-type 'double-float :initial-element 0.0d0)))
    (dolist (entry initial-values)
      (setf (aref state (or (gethash (car entry) species-index)
                             (error "MAKE-CELL-STATE: unknown species ~s" (car entry))))
            (float (cdr entry) 1.0d0)))
    state))

(defun net-reaction-derivative! (output state species-index reactions)
  "Write into OUTPUT (a DOUBLE-FLOAT vector the same length as STATE) the
local (non-spatial) rate of change of every species under REACTIONS alone -
sum over reactions of RATE * stoichiometric-coefficient, per species.
Zeroes OUTPUT first. Returns OUTPUT. Destructive/output-array-taking so the
RK4 substeps in mechanistic-step.lisp don't allocate a fresh vector per
cell per substep across a whole tessellation, every tick."
  (fill output 0.0d0)
  (dolist (rx reactions)
    (let ((rate (funcall (reaction-rate-fn rx) state species-index)))
      (dolist (entry (reaction-stoichiometry rx))
        (let ((i (or (gethash (car entry) species-index)
                      (error "NET-REACTION-DERIVATIVE!: unknown species ~s in reaction ~s"
                             (car entry) (reaction-name rx)))))
          (incf (aref output i) (* rate (cdr entry)))))))
  output)

(defun net-reaction-derivative (state species-index reactions)
  "Non-destructive convenience wrapper over NET-REACTION-DERIVATIVE!, for
tests and small examples where per-call allocation doesn't matter."
  (net-reaction-derivative! (make-array (length state) :element-type 'double-float
                                                        :initial-element 0.0d0)
                             state species-index reactions))

;;; Rate-law constructors. MONOD-RATE and MICHAELIS-MENTEN-RATE are the same
;;; saturating (Hill-type) functional form under two different domain names -
;;; one implementation, shared, so they can't silently diverge.

;;; Every rate law below reads species through CONCENTRATION-OF, so a
;;; species name or a (:PER name water-name) pore-water reference can be
;;; given wherever a name is.

(defun %saturating-rate (vmax catalyst-name limiting-name half-saturation)
  (lambda (state species-index)
    (let ((catalyst (concentration-of state species-index catalyst-name))
          (limiting (max 0.0d0 (concentration-of state species-index limiting-name))))
      (/ (* vmax catalyst limiting) (+ half-saturation limiting)))))

(defun monod-rate (mu-max biomass-name substrate-name half-saturation)
  "Monod growth rate: MU-MAX * [BIOMASS-NAME] * [SUBSTRATE-NAME] /
(HALF-SATURATION + [SUBSTRATE-NAME]). HALF-SATURATION must be > 0."
  (%saturating-rate mu-max biomass-name substrate-name half-saturation))

(defun michaelis-menten-rate (k-cat enzyme-name substrate-name km)
  "Michaelis-Menten rate: K-CAT * [ENZYME-NAME] * [SUBSTRATE-NAME] /
(KM + [SUBSTRATE-NAME]). Used here for microbially-mediated contaminant
degradation, where ENZYME-NAME is a biomass pool acting as the catalyst.
KM must be > 0."
  (%saturating-rate k-cat enzyme-name substrate-name km))

(defun linear-reversible-rate (rate-constant mobile-name sorbed-name partition-coefficient)
  "Net rate transferring mass from MOBILE-NAME to SORBED-NAME, driving the
pair toward the linear equilibrium [SORBED-NAME] = PARTITION-COEFFICIENT *
[MOBILE-NAME] (e.g. a sorption equilibrium: mobile contaminant <-> sorbed
contaminant). Positive when the mobile pool exceeds equilibrium (net
sorption); negative when the sorbed pool exceeds it (net desorption) - a
reversible rate law, unlike the one-directional MONOD-RATE/
MICHAELIS-MENTEN-RATE above."
  (lambda (state species-index)
    (let ((mobile (concentration-of state species-index mobile-name))
          (sorbed (concentration-of state species-index sorbed-name)))
      (* rate-constant (- (* partition-coefficient mobile) sorbed)))))

(defun mass-action-rate (rate-constant &rest species-names)
  "Mass-action rate law: RATE-CONSTANT times the product of the
concentrations of SPECIES-NAMES (zero-order, just RATE-CONSTANT, when none
are given). A species named here but not consumed by the reaction's
stoichiometry acts as a catalyst - e.g. a fixed root source driving
exudation."
  (lambda (state species-index)
    (reduce #'* species-names
            :key (lambda (name) (concentration-of state species-index name))
            :initial-value rate-constant)))

(defun first-order-decay-rate (rate-constant species-name)
  "Rate = RATE-CONSTANT * [SPECIES-NAME] - standard first-order decay, used
for microbial biomass death and abiotic hydrolysis. The one-species case of
MASS-ACTION-RATE."
  (mass-action-rate rate-constant species-name))

(defun hill-activation (x threshold hill-coefficient)
  "Return the Hill activation x^n / (THRESHOLD^n + x^n), n =
HILL-COEFFICIENT: 0 at no signal, 1/2 at X = THRESHOLD, rising to 1. X at or
below zero (round-off) gives 0. THRESHOLD must be positive."
  (if (plusp x)
      (let ((xn (expt (float x 1.0d0) hill-coefficient)))
        (/ xn (+ (expt (float threshold 1.0d0) hill-coefficient) xn)))
      0.0d0))

(defun induced-production-rate (basal induced biomass-name signal-name threshold hill-coefficient)
  "Production rate of a signal by quorum-sensing producers:
[BIOMASS-NAME] * (BASAL + INDUCED * HILL-ACTIVATION([SIGNAL-NAME], THRESHOLD,
HILL-COEFFICIENT)). BASAL and INDUCED are per unit biomass; the induced
term is the positive feedback of autoinduction (luxI expression rising once
the signal passes the LuxR threshold), the form fitted by Fekete et al.
\(2010) - see src/signalling.lisp."
  (lambda (state species-index)
    (* (concentration-of state species-index biomass-name)
       (+ basal (* induced (hill-activation (concentration-of state species-index signal-name)
                                            threshold hill-coefficient))))))
