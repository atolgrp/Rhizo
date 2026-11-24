(in-package :rhizo-diffusion)

;;; Signalling chemistry (PLAN.md, Phase 2): AHL quorum sensing.
;;;
;;; The species set for the thesis, replacing the PFAS demonstration set
;;; (now examples/pfas-rhizobiome.lisp). It is a reaction network like any
;;; other: species, reactions built from the rate laws in
;;; reaction-network.lisp, and a parameter list whose every entry carries
;;; its source (parameters.lisp).
;;;
;;; Model. Producers and receivers are patch-level biomass: an immobile
;;; cell density in each grid cell, fixed during a run. Producers make AHL
;;; at a basal rate plus an induced rate switched on by AHL itself through
;;; a Hill function (autoinduction), as fitted for Pseudomonas putida IsoF
;;; by Fekete et al. (2010). Receivers only listen: a receiver is induced
;;; where AHL exceeds its detection threshold. AHL is lost by abiotic
;;; hydrolysis (lactonolysis, first order) and by biotic degradation
;;; (lactonase activity of degrader biomass, Michaelis-Menten), and
;;; partitions at a finite rate onto the soil (linear sorption, K_d).
;;;
;;; Units. Length um, time s, concentration nM (nmol per litre of pore
;;; water), cell density cells per litre of pore water. Concentrations and
;;; densities are per volume of water, so receivers compare the state
;;; directly with their threshold; sorbed AHL is likewise expressed per
;;; volume of water, with equilibrium sorbed / dissolved = K' = rho_b K_d /
;;; theta (dimensionless; rho_b bulk density, theta volumetric water
;;; content).
;;;
;;; Assumed, not measured:
;;; - The 2D grid is a slab uniform in depth, so a producer patch is a
;;;   column of producers and the steady profile around it is K0, not 1/r.
;;; - Water content is uniform and constant (Phase 3 makes it a field), and
;;;   the diffusivity is that in free water (soil tortuosity is Phase 3).
;;; - Only dissolved AHL is hydrolysed or biodegraded; sorbed AHL is
;;;   protected. With this, sorption does not change the steady state, only
;;;   how fast it is reached.
;;; - Biomass does not grow or die during a run.
;;;
;;; Mass ledger. Two immobile bookkeeping pools, AHL-PRODUCED and AHL-LOST,
;;; record cumulative production and loss, so AHL + AHL-SORBED + AHL-LOST
;;; - AHL-PRODUCED is invariant (AHL-POOL-BALANCE). It changes only if
;;; transport or the reaction bookkeeping fails to conserve mass, which
;;; makes mass balance a pass/fail check for runs with sources and sinks.

(defun ahl-literature-parameters ()
  "Return a fresh list of the MODEL-PARAMETERs of the AHL set, in the units
of the file header, each with its source. The signal is
3-oxo-decanoyl-homoserine lactone (3-oxo-C10-HSL) of Pseudomonas putida
IsoF, a tomato-rhizosphere isolate, unless a note names an analogue."
  (list
   (make-model-parameter
    :diffusivity 490d0 "um^2/s"
    :source "Hense et al. 2012, Sensors 12:4156, doi:10.3390/s120404156, Parameter Choice section, citing their ref. 33: 4.9e-10 m^2/s"
    :note "AHL in free water; the same paper's molar-mass relation gives 9e-10 m^2/s for 3-oxo-C6-HSL. Soil tortuosity is applied in Phase 3.")
   (make-model-parameter
    :basal-production (/ (* 2.3d-19 1d9) 3600d0) "nmol/(cell s)"
    :source "Fekete et al. 2010, FEMS Microbiol Ecol 72:22, doi:10.1111/j.1574-6941.2009.00828.x, parameter table: alpha = 2.3e-19 mol per cell per h")
   (make-model-parameter
    :induced-production (/ (* 2.3d-18 1d9) 3600d0) "nmol/(cell s)"
    :source "Fekete et al. 2010, parameter table: beta = 2.3e-18 mol per cell per h"
    :note "Added to the basal rate at full induction, so a fully induced cell makes 11 alpha.")
   (make-model-parameter
    :induction-threshold 70d0 "nM"
    :source "Fekete et al. 2010, parameter table: A_thresh = 70 nmol/L"
    :note "Their Hill function acts on the LuxR-AHL complex; here on AHL, taking the complex at quasi-steady state in proportion to AHL.")
   (make-model-parameter
    :hill-coefficient 2.5d0 "1"
    :source "Fekete et al. 2010, parameter table: n = 2.5")
   (make-model-parameter
    :receiver-threshold 70d0 "nM"
    :source "Assumed equal to the P. putida IsoF autoinduction threshold (Fekete et al. 2010)"
    :note "Strain-specific: Hense et al. 2012 use about 10 nM for Vibrio fischeri LuxR.")
   (make-model-parameter
    :hydrolysis-rate (/ 0.005545d0 3600d0) "1/s"
    :source "Fekete et al. 2010, parameter table: abiotic degradation gamma = 0.005545 per h"
    :note "Half-life 125 h in culture medium near pH 7. Lactonolysis speeds up with pH and temperature (Yates et al. 2002, Infect Immun 70:5635).")
   (make-model-parameter
    :biodegradation-kcat (/ (/ 13.4d0 4.6d5) 3600d0) "nmol/(cell s)"
    :source "Derived from Wang & Leadbetter 2005, Appl Environ Microbiol 71:1291, doi:10.1128/AEM.71.3.1291-1299.2005 (abstract): Vmax 13.4 nmol/h per g fresh soil divided by a most probable number of 4.6e5 degrading cells per g"
    :note "Hexanoyl-HSL in a turf soil. An MPN undercounts cells, so this per-cell rate is an upper bound.")
   (make-model-parameter
    :biodegradation-km 1500d0 "nM"
    :source "Wang & Leadbetter 2005 (abstract): apparent Km 1.5 uM")
   (make-model-parameter
    :kd 33.61d0 "L/kg"
    :source "Sheng et al. 2018, RSC Adv, doi:10.1039/c7ra10421a, isotherm table: Freundlich K_F of C10-HSL on an oxisol"
    :note "Analogue C10-HSL (no 3-oxo group). Freundlich exponent 0.85, taken as linear. On an alfisol K_F = 84.31.")
   (make-model-parameter
    :desorption-rate (/ (/ 0.33d0 (+ 1d0 (* 0.01d0 33.61d0))) 3600d0) "1/s"
    :source "Derived from Sheng et al. 2018, kinetics table: pseudo-first-order k = 0.33 per h for C10-HSL on the oxisol"
    :note "Batch of 50 mg soil in 5 mL (0.01 kg/L) relaxes at k_des (1 + 0.01 K_d), so k_des = 0.33 / (1 + 0.336) per h.")
   (make-model-parameter
    :bulk-density 1.3d0 "kg/L"
    :note "Typical mineral-soil value, not from a measurement; set to the bulk density of the soil modelled.")
   (make-model-parameter
    :water-content 0.3d0 "1"
    :source "Scenario choice: uniform volumetric water content (Phase 3 makes it a field)")))

(defun ahl-partition-coefficient (parameters)
  "Return the dimensionless equilibrium ratio K' = rho_b K_d / theta of
sorbed to dissolved AHL, both per volume of water, from the :BULK-DENSITY,
:KD and :WATER-CONTENT of PARAMETERS. The retardation factor is 1 + K'."
  (/ (* (parameter-value parameters :bulk-density) (parameter-value parameters :kd))
     (parameter-value parameters :water-content)))

(defun sorption-relaxation-rate (parameters)
  "Return the rate k_des (1 + K') at which dissolved and sorbed AHL relax
to equilibrium, from PARAMETERS. RK4 is stable for a reaction half-step
dt/2 below about 2.8 divided by this rate."
  (* (parameter-value parameters :desorption-rate) (+ 1d0 (ahl-partition-coefficient parameters))))

(defun ahl-species (parameters)
  "Return the species of the AHL set: dissolved AHL (diffusing with
PARAMETERS' :DIFFUSIVITY), sorbed AHL, producer, receiver and degrader
biomass, and the two immobile ledger pools AHL-PRODUCED and AHL-LOST."
  (cons (make-species :name :ahl :diffusivity (float (parameter-value parameters :diffusivity) 1d0))
        (mapcar (lambda (name) (make-species :name name :diffusivity 0d0))
                '(:ahl-sorbed :producer :receiver :degrader :ahl-produced :ahl-lost))))

(defun ahl-reactions (parameters)
  "Return the reactions of the AHL set with PARAMETERS: induced production
by producers, hydrolysis, biodegradation by degraders, and linear sorption
at rate :DESORPTION-RATE toward K' (AHL-PARTITION-COEFFICIENT). Production
also feeds AHL-PRODUCED and every loss feeds AHL-LOST (see file header)."
  (flet ((p (name) (parameter-value parameters name)))
    (list
     (make-reaction
      :name :ahl-production
      :stoichiometry '((:ahl . 1d0) (:ahl-produced . 1d0))
      :rate-fn (induced-production-rate (p :basal-production) (p :induced-production)
                                        :producer :ahl
                                        (p :induction-threshold) (p :hill-coefficient)))
     (make-reaction
      :name :ahl-hydrolysis
      :stoichiometry '((:ahl . -1d0) (:ahl-lost . 1d0))
      :rate-fn (first-order-decay-rate (p :hydrolysis-rate) :ahl))
     (make-reaction
      :name :ahl-biodegradation
      :stoichiometry '((:ahl . -1d0) (:ahl-lost . 1d0))
      :rate-fn (michaelis-menten-rate (p :biodegradation-kcat) :degrader :ahl
                                      (p :biodegradation-km)))
     (make-reaction
      :name :ahl-sorption
      :stoichiometry '((:ahl . -1d0) (:ahl-sorbed . 1d0))
      :rate-fn (linear-reversible-rate (p :desorption-rate) :ahl :ahl-sorbed
                                       (ahl-partition-coefficient parameters))))))

(defun make-ahl-world (tessellation &key (parameters (ahl-literature-parameters)) dt initial-state
                                         (diffusion-scheme :implicit) diffusivity-field)
  "Return a fresh MECHANISTIC-WORLD of the AHL set over TESSELLATION with
PARAMETERS (default AHL-LITERATURE-PARAMETERS). DT, INITIAL-STATE,
DIFFUSION-SCHEME and DIFFUSIVITY-FIELD are as for MAKE-MECHANISTIC-WORLD;
the scheme defaults to :IMPLICIT because the processes span seconds
\(sorption) to days (hydrolysis). Biomass is placed through INITIAL-STATE
\(:PRODUCER, :RECEIVER, :DEGRADER in cells per litre). DT must keep RK4
stable against SORPTION-RELAXATION-RATE."
  (make-mechanistic-world tessellation :species (ahl-species parameters)
                                       :reactions (ahl-reactions parameters)
                                       :dt dt :initial-state initial-state
                                       :diffusion-scheme diffusion-scheme
                                       :diffusivity-field diffusivity-field))

(defun pool-balance (world plus minus)
  "Return the total amount (SPECIES-MASS) of the species named in PLUS
minus that of those in MINUS, in WORLD - an accounting identity such as
stored + lost - produced, which a conserving run keeps constant."
  (- (species-mass world plus) (species-mass world minus)))

(defun ahl-pool-balance (world)
  "Return AHL + AHL-SORBED + AHL-LOST - AHL-PRODUCED in WORLD (POOL-BALANCE).
Invariant in a conserving run of the AHL set; see the file header."
  (pool-balance world '(:ahl :ahl-sorbed :ahl-lost) '(:ahl-produced)))
