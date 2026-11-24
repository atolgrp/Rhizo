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
   #:diffusivity-array
   #:set-species-diffusivity!
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
   #:uniform-column
   #:state-column

   ;; implicit (backward Euler) diffusion
   #:*implicit-tolerance*
   #:vector-dot
   #:conjugate-gradient
   #:implicit-matrix-product
   #:implicit-inverse-diagonal
   #:implicit-diffusion-column
   #:implicit-diffusion-substep!

   ;; mechanistic stepping
   #:stable-timestep
   #:rk4-cell-reaction!
   #:reaction-substep!
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

   ;; parallel per-cell update (single-node, shared-memory)
   #:call-with-kernel
   #:pmap-cell-range!))
