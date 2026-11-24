(in-package :rhizo-diffusion)

;;; Signalling chemistry (PLAN.md, Phase 2): AHL quorum sensing, and root
;;; exudate released from the root surface.
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
;;; Units. Length um, time s. Every pool is an amount per litre of SOIL
;;; (nmol/L for signal, cells/L for biomass), so drying changes no pool
;;; (PLAN.md, Phase 3). Chemistry responds to the pore-water concentration
;;; c = m / theta in nM (nmol per litre of water), read in the rate laws as
;;; (:PER :AHL :WATER); theta is the immobile species :WATER, set by
;;; SET-WATER-CONTENT! (water.lisp), which also makes theta the storage
;;; capacity of dissolved signal and D_0 xi(theta) its diffusivity. Sorbed
;;; signal per soil volume equilibrates with the pore water at rho_b K_d c,
;;; i.e. sorbed / dissolved = K' = rho_b K_d / theta (PARTITION-COEFFICIENT;
;;; rho_b bulk density).
;;;
;;; Assumed, not measured:
;;; - The 2D grid is a slab uniform in depth, so a producer patch is a
;;;   column of producers and the steady profile around it is K0, not 1/r.
;;; - Water content is prescribed (no flow); the relative diffusivity
;;;   xi(theta) is a chosen relation (water.lisp), Millington-Quirk unless
;;;   another is given; Phase 2's free-water assumption is :NO-TORTUOSITY.
;;; - Only dissolved AHL is hydrolysed or biodegraded; sorbed AHL is
;;;   protected. With this, sorption does not change the steady state, only
;;;   how fast it is reached.
;;; - Biomass does not grow or die during a run.
;;;
;;; Every signal (AHL, exudate) has the same four pools, named after it:
;;; the dissolved signal itself (diffusing), SIGNAL-SORBED, and two immobile
;;; bookkeeping pools SIGNAL-PRODUCED and SIGNAL-LOST recording cumulative
;;; production and loss. Then signal + sorbed + lost - produced is invariant
;;; (SIGNAL-POOL-BALANCE). It changes only if transport or the reaction
;;; bookkeeping fails to conserve mass, which makes mass balance a pass/fail
;;; check for runs with sources and sinks. Decay and sorption are built the
;;; same way for every signal (SIGNAL-FATE-REACTIONS).

(defun ahl-literature-parameters ()
  "Return a fresh list of the MODEL-PARAMETERs of the AHL set, in the units
of the file header, each with its source. The signal is
3-oxo-decanoyl-homoserine lactone (3-oxo-C10-HSL) of Pseudomonas putida
IsoF, a tomato-rhizosphere isolate, unless a note names an analogue.
The SOIL-PARAMETERS come last."
  (list*
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
   (soil-parameters)))

(defun partition-coefficient (parameters)
  "Return the dimensionless equilibrium ratio K' = rho_b K_d / theta of
sorbed to dissolved signal, both per volume of water, from the
:BULK-DENSITY, :KD and :WATER-CONTENT of PARAMETERS. The retardation factor
is 1 + K'."
  (/ (* (parameter-value parameters :bulk-density) (parameter-value parameters :kd))
     (parameter-value parameters :water-content)))

(defun sorption-coefficient (parameters)
  "Return rho_b K_d (dimensionless: kg/L times L/kg) from PARAMETERS: the
sorbed amount per soil volume in equilibrium with a unit pore-water
concentration."
  (* (parameter-value parameters :bulk-density) (parameter-value parameters :kd)))

(defun sorption-relaxation-rate (parameters)
  "Return the rate k_des (1 + K') at which dissolved and sorbed signal relax
to equilibrium, from PARAMETERS. RK4 is stable for a reaction half-step
dt/2 below about 2.8 divided by this rate."
  (* (parameter-value parameters :desorption-rate) (+ 1d0 (partition-coefficient parameters))))

(defun signal-pool (signal pool)
  "Return the keyword naming POOL (:SORBED, :PRODUCED or :LOST) of SIGNAL,
e.g. :AHL-SORBED for (SIGNAL-POOL :AHL :SORBED)."
  (intern (format nil "~a-~a" (symbol-name signal) (symbol-name pool)) :keyword))

(defun signal-species (signal diffusivity)
  "Return the four species of SIGNAL (see file header): SIGNAL itself with
DIFFUSIVITY, and its immobile SORBED, PRODUCED and LOST pools."
  (cons (make-species :name signal :diffusivity (float diffusivity 1d0))
        (mapcar (lambda (pool) (make-species :name (signal-pool signal pool) :diffusivity 0d0))
                '(:sorbed :produced :lost))))

(defun production-reaction (signal rate-fn)
  "Return a reaction producing SIGNAL at RATE-FN, recorded in its PRODUCED
pool."
  (make-reaction :name (signal-pool signal :production)
                 :stoichiometry `((,signal . 1d0) (,(signal-pool signal :produced) . 1d0))
                 :rate-fn rate-fn))

(defun loss-reaction (signal name rate-fn)
  "Return a reaction NAME removing dissolved SIGNAL at RATE-FN, recorded in
its LOST pool."
  (make-reaction :name name
                 :stoichiometry `((,signal . -1d0) (,(signal-pool signal :lost) . 1d0))
                 :rate-fn rate-fn))

(defun pore-water (signal)
  "Return the rate-law reference to SIGNAL's pore-water concentration,
\(:PER signal :WATER)."
  (list :per signal :water))

(defun signal-fate-reactions (signal decay-rate parameters)
  "Return the decay and sorption reactions of SIGNAL: first-order decay of
the dissolved amount at DECAY-RATE (into its LOST pool), and linear
sorption at PARAMETERS' :DESORPTION-RATE toward SORPTION-COEFFICIENT times
the pore-water concentration."
  (list (loss-reaction signal (signal-pool signal :decay) (first-order-decay-rate decay-rate signal))
        (make-reaction
         :name (signal-pool signal :sorption)
         :stoichiometry `((,signal . -1d0) (,(signal-pool signal :sorbed) . 1d0))
         :rate-fn (linear-reversible-rate (parameter-value parameters :desorption-rate)
                                          (pore-water signal) (signal-pool signal :sorbed)
                                          (sorption-coefficient parameters)))))

(defun ahl-species (parameters)
  "Return the species of the AHL set: the SIGNAL-SPECIES of :AHL (free-water
diffusivity PARAMETERS' :DIFFUSIVITY), immobile producer, receiver and
degrader biomass, and the WATER-SPECIES."
  (append (signal-species :ahl (parameter-value parameters :diffusivity))
          (mapcar (lambda (name) (make-species :name name :diffusivity 0d0))
                  '(:producer :receiver :degrader))
          (list (water-species))))

(defun ahl-reactions (parameters)
  "Return the reactions of the AHL set with PARAMETERS: induced production
by producers, biodegradation by degraders, and the SIGNAL-FATE-REACTIONS
\(hydrolysis at :HYDROLYSIS-RATE, sorption)."
  (flet ((p (name) (parameter-value parameters name)))
    (list* (production-reaction :ahl (induced-production-rate
                                      (p :basal-production) (p :induced-production)
                                      :producer (pore-water :ahl)
                                      (p :induction-threshold) (p :hill-coefficient)))
           (loss-reaction :ahl :ahl-biodegradation
                          (michaelis-menten-rate (p :biodegradation-kcat) :degrader (pore-water :ahl)
                                                 (p :biodegradation-km)))
           (signal-fate-reactions :ahl (p :hydrolysis-rate) parameters))))

(defun water-content-field (water-content parameters)
  "Return a function of a cell id giving its water content: WATER-CONTENT
itself if it is a function, a constant function if it is a number, and the
constant PARAMETERS' :WATER-CONTENT if it is NIL."
  (cond ((functionp water-content) water-content)
        (water-content (constantly (float water-content 1d0)))
        (t (constantly (parameter-value parameters :water-content)))))

(defun make-ahl-world (tessellation &key (parameters (ahl-literature-parameters)) dt initial-state
                                         (diffusion-scheme :imex) water-content
                                         (relative-diffusivity :millington-quirk))
  "Return a fresh MECHANISTIC-WORLD of the AHL set over TESSELLATION with
PARAMETERS (default AHL-LITERATURE-PARAMETERS), its water content set by
SET-WATER-CONTENT! from WATER-CONTENT (a number, a function of a cell id,
or NIL for the parameters' :WATER-CONTENT; see WATER-CONTENT-FIELD) and the
relation named RELATIVE-DIFFUSIVITY (RELATIVE-DIFFUSIVITY-MODEL; default
Millington-Quirk, the Phase 3 baseline). DT, INITIAL-STATE and
DIFFUSION-SCHEME are as for MAKE-MECHANISTIC-WORLD; the scheme defaults to
:IMEX, which keeps steady states exact and is stable for diffusion at any
DT. Biomass is placed through INITIAL-STATE (:PRODUCER, :RECEIVER,
:DEGRADER in cells per litre of soil). DT times the fastest reaction rate -
often SORPTION-RELAXATION-RATE - must stay below 1 for :IMEX's explicit
reactions (about 5.6 for the Strang schemes' RK4)."
  (set-water-content!
   (make-mechanistic-world tessellation :species (ahl-species parameters)
                                        :reactions (ahl-reactions parameters)
                                        :dt dt :initial-state initial-state
                                        :diffusion-scheme diffusion-scheme)
   (water-content-field water-content parameters)
   (relative-diffusivity-model relative-diffusivity parameters)))

(defun pool-balance (world plus minus)
  "Return the total amount (SPECIES-MASS plus SPECIES-EXPORTED, i.e.
including what water carried out of the domain) of the species named in
PLUS minus that of those in MINUS, in WORLD - an accounting identity such
as stored + lost - produced, which a conserving run keeps constant."
  (- (+ (species-mass world plus) (species-exported world plus))
     (+ (species-mass world minus) (species-exported world minus))))

(defun signal-pool-balance (world signal)
  "Return SIGNAL + its SORBED and LOST pools - its PRODUCED pool, in WORLD
\(POOL-BALANCE). Invariant in a conserving run; see the file header."
  (pool-balance world (list signal (signal-pool signal :sorbed) (signal-pool signal :lost))
                (list (signal-pool signal :produced))))

;;; Root exudate (PLAN.md, Phase 2): a root releases a signal from its
;;; surface - a placeholder for a strigolactone or flavonoid. Root cells
;;; are not soil: the exudate cannot diffuse into them (their diffusivity
;;; is zero, so the harmonic face mean makes every root face a no-flux
;;; wall), and each soil cell touching the root receives the exudation flux
;;; over its faces with the root. That is a prescribed-flux (Neumann)
;;; condition on the root surface, written as a source in the first layer
;;; of soil cells.
;;;
;;; The root surface is the set of cell faces between root and soil. On a
;;; pixelated Voronoi grid a straight root's surface is longer than the
;;; straight line (faces zigzag), so the flux is per unit of this perimeter;
;;; divide a total rate by ROOT-PERIMETER to release exactly that total.
;;;
;;; The flux is in nmol/L um/s (amount per unit root surface per second,
;;; the amount per litre of soil), so flux * face length / cell area is a
;;; rate per soil volume in nmol/L/s.

(defun root-perimeter (tessellation root-p)
  "Return the total length of the faces between cells satisfying the
predicate ROOT-P (root) and cells that do not (soil) in TESSELLATION."
  (loop for id in (cell-ids tessellation)
        when (funcall root-p id)
          sum (loop for nb in (neighbors-of tessellation id)
                    unless (funcall root-p nb)
                      sum (edge-length tessellation id nb))))

(defun root-surface-source (tessellation root-p flux)
  "Return a function of a cell id giving the exudation rate per soil volume
\(nmol/L/s) in that cell: FLUX times the length of its faces with root cells
\(ROOT-P), divided by its CELL-AREA, for soil cells; zero for root cells."
  (lambda (id)
    (if (funcall root-p id)
        0d0
        (/ (* flux (loop for nb in (neighbors-of tessellation id)
                         when (funcall root-p nb)
                           sum (edge-length tessellation id nb)))
           (cell-area tessellation id)))))

(defun exudate-literature-parameters ()
  "Return a fresh list of the MODEL-PARAMETERs of the root exudate, a
strigolactone placeholder modelled on 5-deoxystrigol (C19H22O5, 330.38
g/mol), in the units of the file header. Flagged where no source was found
\(2026-11-24): the exudation flux, K_d and sorption rate."
  (let ((ahl-d (parameter-value (ahl-literature-parameters) :diffusivity)))
    (list*
     (make-model-parameter
      :diffusivity (* ahl-d (expt (/ 269.34d0 330.38d0) (/ 1d0 3d0))) "um^2/s"
      :source "Derived: the AHL diffusivity (Hense et al. 2012) scaled by (M_AHL / M_5DS)^(1/3), Stokes-Einstein for molecules of equal density; 3-oxo-C10-HSL C14H23NO4 269.34 g/mol, 5-deoxystrigol 330.38 g/mol")
     (make-model-parameter
      :decay-rate (/ (log 2d0) (* 1.5d0 86400d0)) "1/s"
      :source "Akiyama et al. 2010 (Plant Cell Physiol): 5-deoxystrigol half-life 1.5 d in water, as cited by Klein et al. 2024, Front Microbiol, doi:10.3389/fmicb.2024.1470469"
      :note "Strigolactones hydrolyse faster above pH 7.5; no soil half-life found.")
     (make-model-parameter
      :detection-threshold 0.1d0 "nM"
      :source "Rial et al. 2018, Phytochem Anal, doi:10.1002/pca.2796 (abstract): strigolactones are active below 1e-10 M"
      :note "Germination of parasitic-plant seeds; arbuscular mycorrhizal hyphae may respond lower.")
     (make-model-parameter
      :exudation-flux 0.05d0 "nmol/L um/s"
      :note "Placeholder: gives a steady dissolved concentration J / sqrt(D k) of about 1 nM at a flat root surface, ten times the detection threshold. Replace with a measured exudation rate.")
     (make-model-parameter
      :kd 33.61d0 "L/kg"
      :note "Placeholder: the AHL analogue value (C10-HSL on an oxisol, Sheng et al. 2018), both being lactones of similar hydrophobicity. No measured strigolactone K_d found.")
     (make-model-parameter
      :desorption-rate (parameter-value (ahl-literature-parameters) :desorption-rate) "1/s"
      :note "Placeholder: the AHL analogue value (Sheng et al. 2018). No strigolactone sorption kinetics found.")
     (soil-parameters))))

(defun exudate-species (parameters)
  "Return the species of the root exudate set: the SIGNAL-SPECIES of
:EXUDATE (free-water diffusivity PARAMETERS' :DIFFUSIVITY), the immobile
:ROOT-SOURCE (each cell's exudation rate per soil volume), and the
WATER-SPECIES."
  (append (signal-species :exudate (parameter-value parameters :diffusivity))
          (list (make-species :name :root-source :diffusivity 0d0) (water-species))))

(defun exudate-reactions (parameters)
  "Return the reactions of the root exudate set: release at the rate held
in :ROOT-SOURCE, and the SIGNAL-FATE-REACTIONS (decay at :DECAY-RATE,
sorption) with PARAMETERS."
  (cons (production-reaction :exudate (mass-action-rate 1d0 :root-source))
        (signal-fate-reactions :exudate (parameter-value parameters :decay-rate) parameters)))

(defun make-root-exudate-world (tessellation root-p &key (parameters (exudate-literature-parameters))
                                                         dt initial-state
                                                         (diffusion-scheme :imex) water-content
                                                         (relative-diffusivity :millington-quirk))
  "Return a fresh MECHANISTIC-WORLD of the root exudate set over
TESSELLATION, the cells satisfying ROOT-P being root: their :ROOT-SOURCE
comes from ROOT-SURFACE-SOURCE with PARAMETERS' :EXUDATION-FLUX, and they
are no-flux walls for the exudate (SET-WATER-CONTENT!'s BLOCKED-P).
WATER-CONTENT and RELATIVE-DIFFUSIVITY are as for MAKE-AHL-WORLD.
INITIAL-STATE is an alist or a function of a cell id, as for
MAKE-MECHANISTIC-WORLD, and must not set :ROOT-SOURCE. The scheme defaults
to :IMEX, as for MAKE-AHL-WORLD."
  (let ((source (root-surface-source tessellation root-p (parameter-value parameters :exudation-flux))))
    (set-water-content!
     (make-mechanistic-world
      tessellation :species (exudate-species parameters) :reactions (exudate-reactions parameters)
                   :dt dt :diffusion-scheme diffusion-scheme
                   :initial-state (lambda (id)
                                    (acons :root-source (funcall source id)
                                           (if (functionp initial-state)
                                               (funcall initial-state id)
                                               initial-state))))
     (water-content-field water-content parameters)
     (relative-diffusivity-model relative-diffusivity parameters)
     :blocked-p root-p)))
