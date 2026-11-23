;;;; poros.asd --- System definition.
;;;;
;;;; Deliberately minimal. The module split below is the INTENDED architecture
;;;; from PLAN.md, not a description of working code -- only the probe in
;;;; experiments/ runs today. Phase 0 may still change this structure, since
;;;; whether Poros supplies its own finite-volume layer is unresolved.

(asdf:defsystem #:poros
  :description "Porous-media and reactive-transport layer for Common Lisp."
  :author "Apostolos Almpanis"
  :license "MIT"
  :version "0.0.1"
  :depends-on (#:femlisp)
  :serial t
  :components ((:module "src"
                :serial t
                :components ((:file "package")))))

;; Planned, not yet present -- see PLAN.md:
;;
;;   poros/flow       Phase 1  Richards equation, van Genuchten-Mualem,
;;                             mass-conservative mixed form
;;   poros/transport  Phase 2  advection-dispersion, sorption, degradation,
;;                             operator splitting
;;   poros/network    Phase 4  time-varying directed signalling graph,
;;                             percolation detection
