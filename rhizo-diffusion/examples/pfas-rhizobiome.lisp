(in-package :rhizo-diffusion.examples)

;;; Example reaction network: PFAS degradation in a two-guild rhizobiome.
;;;
;;; This was the engine's first ("v1") species set. PLAN.md, Phase 2,
;;; replaced it with signalling chemistry (src/signalling.lisp) and moved it
;;; here, outside the engine, to show that the engine is general: nothing in
;;; src/ knows about these species, and this file uses only the public
;;; reaction-network and mechanistic-world API.
;;;
;;; Its parameters are illustrative, not taken from the literature, so this
;;; set is a demonstration, not a model of a real system.

;;; The v1 concrete species/reaction set: two microbial guilds competing for
;;; one shared substrate pool, one of which (the "degrader") also
;;; co-metabolically degrades a contaminant that partitions between a
;;; mobile (diffusing) and sorbed (immobile) pool at linear equilibrium. A
;;; MINERALIZED sink pool tracks degraded contaminant so total system mass
;;; (all six pools summed) stays a meaningful, checkable quantity - see
;;; tests/pfas-example-tests.lisp. All stoichiometric coefficients are +/-1 -
;;; i.e. every reaction here is an exact 1:1 mass conversion between pools
;;; (yield = 1.0), the simplest choice that keeps the whole-system mass
;;; conservation check a plain sum rather than a yield-weighted one. A real
;;; PFAS/rhizobiome model would give growth reactions yield < 1 (not all
;;; consumed substrate becomes biomass) - that only changes the parameters
;;; passed to MAKE-REACTION's stoichiometry, not this engine.

(defun make-pfas-rhizobiome-species (&key (diffusivity-contaminant 0.1d0)
                                           (diffusivity-substrate 0.1d0))
  (list (make-species :name :contaminant-mobile
                       :diffusivity (float diffusivity-contaminant 1.0d0))
        (make-species :name :contaminant-sorbed :diffusivity 0.0d0)
        (make-species :name :substrate :diffusivity (float diffusivity-substrate 1.0d0))
        (make-species :name :biomass-degrader :diffusivity 0.0d0)
        (make-species :name :biomass-competitor :diffusivity 0.0d0)
        (make-species :name :mineralized :diffusivity 0.0d0)))

(defun make-pfas-rhizobiome-reactions (&key (mu-max-degrader 0.5d0) (mu-max-competitor 0.5d0)
                                             (ks-degrader 1.0d0) (ks-competitor 1.0d0)
                                             (k-deg 0.2d0) (km-deg 1.0d0)
                                             (k-sorb 0.1d0) (kd 2.0d0)
                                             (decay-rate-degrader 0.05d0)
                                             (decay-rate-competitor 0.05d0))
  (list
   (make-reaction
    :name :sorption
    :stoichiometry '((:contaminant-mobile . -1.0d0) (:contaminant-sorbed . 1.0d0))
    :rate-fn (linear-reversible-rate k-sorb :contaminant-mobile :contaminant-sorbed kd))
   (make-reaction
    :name :degrader-growth
    :stoichiometry '((:substrate . -1.0d0) (:biomass-degrader . 1.0d0))
    :rate-fn (monod-rate mu-max-degrader :biomass-degrader :substrate ks-degrader))
   (make-reaction
    :name :competitor-growth
    :stoichiometry '((:substrate . -1.0d0) (:biomass-competitor . 1.0d0))
    :rate-fn (monod-rate mu-max-competitor :biomass-competitor :substrate ks-competitor))
   (make-reaction
    :name :contaminant-degradation
    :stoichiometry '((:contaminant-mobile . -1.0d0) (:mineralized . 1.0d0))
    :rate-fn (michaelis-menten-rate k-deg :biomass-degrader :contaminant-mobile km-deg))
   (make-reaction
    :name :degrader-decay
    :stoichiometry '((:biomass-degrader . -1.0d0) (:substrate . 1.0d0))
    :rate-fn (first-order-decay-rate decay-rate-degrader :biomass-degrader))
   (make-reaction
    :name :competitor-decay
    :stoichiometry '((:biomass-competitor . -1.0d0) (:substrate . 1.0d0))
    :rate-fn (first-order-decay-rate decay-rate-competitor :biomass-competitor))))

(defun make-pfas-rhizobiome-world (tessellation &key (dt 0.05d0) initial-state
                                                      species-args reaction-args)
  "Convenience constructor for the v1 species/reaction set above.
SPECIES-ARGS/REACTION-ARGS are plists forwarded to
MAKE-PFAS-RHIZOBIOME-SPECIES/-REACTIONS, for parameter sweeps. DT's default
is a conservative starting point for the default parameters above, not a
computed bound - call STABLE-TIMESTEP (mechanistic-step.lisp) directly for
a real one once diffusivities, geometry, or reaction rates change."
  (let ((species (apply #'make-pfas-rhizobiome-species species-args))
        (reactions (apply #'make-pfas-rhizobiome-reactions reaction-args)))
    (make-mechanistic-world tessellation :species species :reactions reactions
                                          :dt dt :initial-state initial-state)))
