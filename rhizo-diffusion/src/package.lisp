(defpackage :rhizo-diffusion
  (:use :cl)
  (:export
   ;; tessellation protocol - implement these three on any new representation
   ;; (CA-grown grid, GIS polygons, a Delaunay triangulation, ...) to make it
   ;; usable by everything below.
   #:cell-ids
   #:neighbors-of
   #:cell-count
   #:edge-weight
   #:cell-area
   #:edge-length
   #:face-distance
   #:cell-position
   #:nearest-cell
   #:well-mixed-tessellation

   ;; concrete protocol implementation backed by an explicit adjacency table
   #:hash-graph-tessellation
   #:make-hash-graph-tessellation
   #:ensure-cell!
   #:add-edge!
   #:set-edge-length!
   #:face-distances-from-positions
   #:positions-table
   #:make-fv-tessellation

   ;; dense cell-id <-> array-index mapping, for the mechanistic engine
   #:cell-index
   #:make-cell-index
   #:cell-index-count
   #:index->cell-id
   #:cell-id->index

   ;; CA-grown Voronoi builder
   #:voronoi-grid
   #:voronoi-grid-tessellation
   #:voronoi-grid-width
   #:voronoi-grid-height
   #:voronoi-grid-id-at
   #:voronoi-grid-seed
   #:voronoi-grid-seed-points
   #:make-voronoi-grid

   ;; nearest-seed Voronoi builder - the grid for reaction-diffusion
   #:seed-bin-size
   #:seed-bins
   #:ring-bins
   #:nearest-seed-index
   #:nearest-seed-ids
   #:boundary-step-counts
   #:face-length-from-steps
   #:nearest-seed-tessellation
   #:make-nearest-seed-grid

   ;; world: the explicit, non-global state of one simulation instance
   #:world
   #:make-world
   #:world-tessellation
   #:world-molecules
   #:world-frontier
   #:world-random-state
   #:world-seed
   #:world-tick

   ;; diffusion
   #:seed-molecules-at!
   #:step-front-propagation!
   #:front-active-p

   ;; metrics
   #:compass-measurement-points
   #:run-front-propagation-with-metrics
   #:isotropy-coefficient-of-variation

   ;; io
   #:write-csv

   ;; parallel replicate runner
   #:run-replicates

   ;; GeoJSON ingestion
   #:region
   #:make-region
   #:region-id
   #:region-ring
   #:region-centroid
   #:ring-centroid
   #:parse-geojson-regions

   ;; GIS tessellation builder
   #:canonical-edge-key
   #:centroid-distance
   #:point-distance
   #:build-gis-adjacency
   #:make-gis-tessellation
   #:tessellation-from-geojson
   #:polygon-area

   ;; weighted (Dijkstra) diffusion
   #:compute-weighted-arrival-times

   ;; PNG rendering
   #:id->color
   #:render-voronoi-grid-png
   #:render-front-propagation-png

   ;; cell-shape regularity metrics
   #:cell-shape-metrics
   #:isoperimetric-quotient
   #:mean-isoperimetric-quotient

   ;; physical-sanity validation
   #:linear-regression
   #:pearson-correlation
   #:cell-centroids
   #:euclidean-distance
   #:distance-correlation
   #:front-growth-curve
   #:front-growth-linearity

   ;; reaction network (species, reactions, rate laws)
   #:species
   #:make-species
   #:species-name
   #:species-diffusivity
   #:reaction
   #:make-reaction
   #:reaction-name
   #:reaction-stoichiometry
   #:reaction-rate-fn
   #:make-species-index
   #:state-value
   #:concentration-of
   #:make-cell-state
   #:net-reaction-derivative!
   #:net-reaction-derivative
   #:monod-rate
   #:michaelis-menten-rate
   #:linear-reversible-rate
   #:first-order-decay-rate
   #:mass-action-rate
   #:hill-activation
   #:induced-production-rate

   ;; model parameters with their literature sources
   #:model-parameter
   #:make-model-parameter
   #:model-parameter-name
   #:model-parameter-value
   #:model-parameter-unit
   #:model-parameter-source
   #:model-parameter-note
   #:find-parameter
   #:parameter-value
   #:flagged-parameters
   #:with-parameter-values
   #:parameter-table-rows

   ;; mechanistic world: state container for the reaction-diffusion engine
   #:mechanistic-world
   #:make-mechanistic-world
   #:mechanistic-world-tessellation
   #:mechanistic-world-cell-index
   #:mechanistic-world-stencil
   #:mechanistic-world-species-list
   #:mechanistic-world-species-index
   #:mechanistic-world-reactions
   #:mechanistic-world-n-species
   #:mechanistic-world-state
   #:mechanistic-world-dt
   #:mechanistic-world-tick
   #:mechanistic-world-diffusion-scheme
   #:mechanistic-world-diffusivities
   #:mechanistic-world-capacities
   #:species-field-array
   #:diffusivity-array
   #:capacity-array
   #:set-species-column!
   #:set-species-diffusivity!
   #:set-species-capacity!
   #:world-stable-timestep
   #:field-initial-state

   ;; precomputed finite-volume diffusion operator
   #:face-conductance
   #:flux-stencil
   #:make-flux-stencil
   #:flux-stencil-offsets
   #:flux-stencil-neighbors
   #:flux-stencil-conductances
   #:flux-stencil-areas
   #:stencil-total-conductance
   #:stencil-net-flux
   #:harmonic-mean
   #:stencil-stable-timestep
   #:intensity-column
   #:uniform-column
   #:state-column

   ;; implicit (backward Euler) diffusion
   #:*implicit-tolerance*
   #:vector-dot
   #:conjugate-gradient
   #:implicit-matrix-product
   #:implicit-inverse-diagonal
   #:face-coefficients
   #:capacity-diagonal
   #:apply-implicit-matrix
   #:implicit-diffusion-column
   #:implicit-diffusion-substep!

   ;; mechanistic stepping
   #:stable-timestep
   #:rk4-cell-reaction!
   #:reaction-substep!
   #:snap-underflow
   #:euler-reaction-substep!
   #:diffusion-substep!
   #:scheme-diffusion-substep!
   #:step-mechanistic!
   #:run-mechanistic!

   ;; accuracy against analytic solutions
   #:gaussian-field
   #:species-cells
   #:species-mass
   #:species-mean-squared-distance
   #:relative-l2-error
   #:observed-order
   #:gaussian-spreading-error
   #:grid-gaussian-errors
   #:point-source-concentration
   #:point-source-arrival-time
   #:*compass-directions*
   #:compass-points
   #:threshold-arrival-times
   #:anisotropy-index
   #:point-source-isotropy
   #:retardation-measurement

   ;; soil water and solute diffusion
   #:millington-quirk
   #:papendick-campbell
   #:mullins-sommer-threshold
   #:mullins-sommer
   #:soil-parameters
   #:relative-diffusivity-model
   #:no-tortuosity
   #:water-species
   #:set-water-content!
   #:wet-dry-period
   #:wet-dry-cycle
   #:run-with-water!

   ;; signalling chemistry: AHL quorum sensing and root exudate
   #:partition-coefficient
   #:sorption-coefficient
   #:sorption-relaxation-rate
   #:pore-water
   #:water-content-field
   #:signal-pool
   #:signal-species
   #:production-reaction
   #:loss-reaction
   #:signal-fate-reactions
   #:pool-balance
   #:signal-pool-balance
   #:ahl-literature-parameters
   #:ahl-species
   #:ahl-reactions
   #:make-ahl-world
   #:root-perimeter
   #:root-surface-source
   #:exudate-literature-parameters
   #:exudate-species
   #:exudate-reactions
   #:make-root-exudate-world

   ;; analysis of signalling runs (Phase 2 gate)
   #:bessel-k
   #:bessel-k0
   #:bessel-k1
   #:point-source-steady-state
   #:species-column-change
   #:step-to-steady-state!
   #:single-producer-steady-state
   #:uniform-ahl-balance
   #:bisect-root
   #:uniform-steady-states
   #:golden-section-maximum
   #:quorum-switch-density
   #:well-mixed-ahl-world
   #:well-mixed-ahl
   #:uniform-steady-ahl
   #:uniform-induction-p
   #:simulated-switch-density

   ;; calling distance and its sensitivity (Phase 2 gate)
   #:crossing-distance
   #:steady-calling-distance
   #:steady-calling-distance-elasticities
   #:point-source-transient
   #:transient-calling-distance
   #:log-elasticity
   #:radial-profile
   #:profile-crossing
   #:calling-distance
   #:default-calling-scenario
   #:scenario-value
   #:calling-distance-scenario
   #:calling-distance-sensitivity

   ;; signalling range against water content (Phase 3 gate)
   #:range-versus-water
   #:water-transition
   #:water-sensitivity
   #:single-producer-calling-distance

   ;; the signalling graph and its percolation (Phase 4)
   #:*disc-percolation-filling-factor*
   #:random-patch-positions
   #:signalling-edges
   #:group-components
   #:weak-components
   #:strong-components
   #:largest-component-size
   #:spanning-p
   #:critical-connection-radius
   #:graph-connectivity
   #:percolation-at-range
   #:connectivity-over-cycle
   #:branch-value
   #:hysteresis-gap
   #:edge-damkohler
   #:range-at-water
   #:predicted-percolation-water-content
   #:percolation-water-content

   ;; tagged signal: each patch's AHL as its own species (Phase 4)
   #:tag-name
   #:node-name
   #:tagged-species
   #:total-tag-concentration
   #:tagged-reactions
   #:make-tagged-ahl-world
   #:tagged-edges
   #:tagged-pool-balance

   ;; soil hydraulics: van Genuchten-Mualem relations (Phase 5)
   #:van-genuchten
   #:make-van-genuchten
   #:van-genuchten-p
   #:van-genuchten-theta-r
   #:van-genuchten-theta-s
   #:van-genuchten-alpha
   #:van-genuchten-n
   #:van-genuchten-ks
   #:vg-m
   #:vg-effective-saturation
   #:vg-water-content
   #:vg-capacity
   #:vg-conductivity
   #:vg-conductivity-derivative
   #:vg-head
   #:celia-sand-parameters
   #:sandy-clay-loam-parameters
   #:van-genuchten-soil

   ;; outer-boundary geometry and rectangular grids (Phase 5)
   #:boundary-faces
   #:make-rectangular-tessellation
   #:voronoi-grid-ids
   #:raster-boundary-faces

   ;; non-symmetric sparse solves (Phase 5)
   #:linear-solver-failure
   #:linear-solver-failure-solver
   #:linear-solver-failure-iterations
   #:linear-solver-failure-residual
   #:jacobi-preconditioner
   #:bicgstab
   #:apply-csr-matrix
   #:vector-axpy
   #:vector-norm
   #:residual-vector

   ;; Richards flow: mixed form, Newton, adaptive steps (Phase 5)
   #:depth-elevation
   #:time-value
   #:richards-boundary
   #:richards-boundary-cell
   #:richards-boundary-length
   #:richards-boundary-distance
   #:richards-boundary-elevation
   #:richards-boundary-kind
   #:make-richards-boundaries
   #:richards-problem
   #:make-richards-problem
   #:richards-problem-tessellation
   #:richards-problem-cell-index
   #:richards-problem-stencil
   #:richards-problem-elevations
   #:richards-problem-soils
   #:richards-problem-boundaries
   #:richards-problem-face-mean
   #:richards-problem-form
   #:richards-state
   #:make-richards-state
   #:richards-state-time
   #:richards-state-heads
   #:richards-state-water-contents
   #:richards-state-face-volumes
   #:richards-state-boundary-volumes
   #:richards-state-next-dt
   #:richards-water-contents
   #:face-conductivity
   #:boundary-inflow
   #:richards-system
   #:richards-nonconvergence
   #:richards-step
   #:run-richards
   #:richards-storage
   #:richards-cumulative-inflow
   #:richards-mass-balance-error
   #:richards-mass-balance-ratio
   #:cell-id-
   #:index-

   ;; advective transport: flow fields and exports (Phase 5)
   #:set-flow-field!
   #:species-exported
   #:mechanistic-world-face-fluxes
   #:mechanistic-world-boundary-outflows
   #:mechanistic-world-exports
   #:advection-diffusion-column
   #:safe-inverse-diagonal

   ;; flow-driven signalling: coupling and directional links (Phase 5)
   #:richards-mean-fluxes
   #:step-flow-transport!
   #:run-with-flow!
   #:advective-decay-rate
   #:advective-point-source-steady-state
   #:advective-point-source-transient
   #:directional-calling-distance
   #:flow-signal-field
   #:flow-signalling-edges
   #:one-way-fraction
   #:directed-spanning-p
   #:flow-graph-connectivity

   ;; steady rain-fed water content (Phase 5)
   #:steady-drainage-head
   #:steady-drainage-water-content
   #:flow-soil-parameters

   ;; prescribed uniform flow and directional profiles (Phase 5)
   #:uniform-flow-field
   #:sector-predicate

   ;; parallel per-cell update (single-node, shared-memory)
   #:call-with-kernel
   #:pmap-cell-range!))
