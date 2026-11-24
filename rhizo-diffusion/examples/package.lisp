(defpackage :rhizo-diffusion.examples
  (:use :cl :rhizo-diffusion)
  (:documentation "Example reaction networks built only on the public
rhizo-diffusion API, kept outside the engine (PLAN.md, Phase 2).")
  (:export
   ;; PFAS degradation in a two-guild rhizobiome (the former v1 set)
   #:make-pfas-rhizobiome-species
   #:make-pfas-rhizobiome-reactions
   #:make-pfas-rhizobiome-world))
