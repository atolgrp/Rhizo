(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; The PFAS degrader/competitor set lives in examples/ (PLAN.md, Phase 2):
;;; it shows the engine is general, but it is not part of the engine.

(test engine-exports-no-pfas-symbols
  (let ((pfas (loop for symbol being the external-symbols of :rhizo-diffusion
                    when (search "PFAS" (symbol-name symbol)) collect symbol)))
    (is (null pfas) "the engine package should not export the PFAS example: ~s" pfas)))

(test mechanistic-mass-is-conserved-over-many-ticks
  ;; Both halves of STEP-MECHANISTIC! conserve total system mass by
  ;; construction: diffusion moves mass between cells symmetrically (never
  ;; creating or destroying it, per species - the flux cell I gains is
  ;; exactly the flux cell J loses), and every PFAS-example reaction (see
  ;; examples/pfas-rhizobiome.lisp) is a 1:1 stoichiometric exchange between two of
  ;; the six pools, so a cell's own total across species is unchanged by
  ;; reactions too. If this test ever drifts, one of those two properties
  ;; broke.
  (let* ((grid (make-voronoi-grid 12 12 15 :seed 7))
         (tess (voronoi-grid-tessellation grid))
         (world (make-pfas-rhizobiome-world
                 tess :dt 0.02d0
                      :initial-state (lambda (id)
                                       (declare (ignore id))
                                       '((:contaminant-mobile . 5.0d0) (:substrate . 3.0d0)
                                         (:biomass-degrader . 0.5d0) (:biomass-competitor . 0.5d0))))))
    (let ((initial-mass (total-mechanistic-mass world)))
      (dotimes (i 200) (step-mechanistic! world))
      (let ((final-mass (total-mechanistic-mass world)))
        (is (< (abs (/ (- final-mass initial-mass) initial-mass)) 1d-6)
            "total system mass should be conserved: initial=~a final=~a" initial-mass final-mass)))))

(test degrader-guild-suppresses-competitor-growth
  ;; The actual multi-guild ecological interaction the PFAS example
  ;; (examples/pfas-rhizobiome.lisp) was built to demonstrate: both guilds draw
  ;; down the same shared SUBSTRATE pool via independent Monod terms, so an
  ;; active degrader guild should measurably slow the competitor's growth
  ;; relative to a control where the degrader never activates - starting it
  ;; at zero biomass is a stable fixed point of its own Monod growth term
  ;; (rate is proportional to biomass, so zero stays zero), making a clean
  ;; apples-to-apples control using the exact same reaction network.
  (let ((tess (single-cell-tessellation 1.0d0)))
    (flet ((final-competitor-biomass (initial-degrader-biomass)
             (let ((world (make-pfas-rhizobiome-world
                           tess :dt 0.05d0
                                :initial-state
                                (lambda (id)
                                  (declare (ignore id))
                                  (list (cons :substrate 5.0d0)
                                        (cons :biomass-competitor 0.1d0)
                                        (cons :biomass-degrader initial-degrader-biomass))))))
               (dotimes (i 300) (step-mechanistic! world))
               (cell-species-value world "only" :biomass-competitor))))
      (let ((control (final-competitor-biomass 0.0d0))
            (treatment (final-competitor-biomass 0.3d0)))
        (is (> control treatment)
            "competitor biomass should be lower competing with an active degrader guild: control=~a treatment=~a"
            control treatment)))))
