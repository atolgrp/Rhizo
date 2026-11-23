# Poros — project brief

*A standing brief. Read this first. It is also the prompt to paste into a fresh
chat to bring an assistant up to speed on this project in one message.*

## What this is

**Poros** (πόρος — passage, pore) is a porous-media and reactive-transport layer
for Common Lisp. It does not exist yet; building it is the project.

It is being built as a **general framework, publishable in its own right**, and
only then used for the scientific question that motivated it (below). That
ordering is deliberate and was chosen with the downside understood: a framework
is slower to a first result than a one-off model, but it is reusable, citable,
and the thing the ecosystem is actually missing.

## Why Common Lisp

Not nostalgia, and not a constraint to be worked around. Python has FiPy,
PorePy, OpenPNM; C++ has DuMuX and OpenFOAM; Julia has a growing stack. Common
Lisp has **Femlisp** — a general finite-element framework — and nothing above it
for porous media. That gap is real and narrow enough for one person to fill.

The honest counter-argument, which should be restated rather than forgotten:
the audience for a Common Lisp scientific library is small. This is the main
strategic risk and it is accepted knowingly, not overlooked.

## The scientific target (the reason the framework exists)

Simulate **signal exchange between plant roots and the rhizosphere microbial
community**, where the effective range of a chemical signal is gated by soil
water.

The specific, falsifiable claim to be tested once the framework exists:

> The rhizosphere signalling network undergoes a **percolation transition at a
> critical water content**, with hysteresis across wetting–drying cycles. Below
> θ_c the network fragments into isolated root–microbe islands; above it a
> spanning component forms. This mechanistically explains the measured
> soil-to-soil variation in signalling distance.

Signalling is modelled as a **time-varying directed graph**, not only a
concentration field: nodes are root segments and microbial patches; a directed
edge exists at time *t* if a molecule released at the source reaches the
receiver above its detection threshold before degrading. Advection makes edges
directional — during infiltration, downgradient nodes become reachable while
upgradient ones are cut off.

Three dimensionless groups govern every edge: **Péclet** (advection vs
diffusion, swinging ~4 orders of magnitude between quiescent and flooding
regimes), **Damköhler** (degradation vs transport), and **hydration
connectivity** θ/θ_c.

### Anchors and constraints established by literature review

- Effective signalling distance for strigolactones measured at **0.05–14.2 mm
  across soils** (Ziadne, Gafni & Paporisch 2026, *Pest Manag. Sci.*,
  10.1002/ps.71085). Unexplained mechanistically — this is the dataset to
  reproduce. Caveat: synthetic racemic GR24, not natural strigolactones.
- GR24 dissipates below detection within 24 h in live soil (Wallach et al.,
  *Weed Research* 65:e12619, 10.1111/wre.12619). No DT50 in days exists.
- **Sorption is the main scientific threat to the hypothesis.** If effective
  range is set by the sorption isotherm rather than water-phase connectivity,
  the percolation claim weakens. Sorption must be a first-class term from day
  one, not a later correction. Demonstrated for flavonoids
  (10.1126/sciadv.aax8254); materially controls AHLs (Sheng et al. 2018,
  *RSC Adv.* 8:9364, 10.1039/c7ra10421a).
- **AHLs are the best-parameterised first target** (measured Kd and isotherms,
  abiotic hydrolysis half-lives, biosensor detection thresholds).
  Strigolactones are the best-motivated but worst-parameterised — no measured
  aqueous diffusion coefficient or Kd exists in the literature.
- Closest prior art: Bickel & Or 2023 (*Commun. Biol.* 6:322,
  10.1038/s42003-023-04703-7) on moisture-gated microbial interaction ranges;
  Jung & Meile 2020 (*PeerJ* 8:e9942) on quorum-sensing calling distance under
  flow. Neither builds a percolating signalling graph. The novelty is the
  synthesis and the falsifiable θ_c prediction, not the ingredients.
- Likely convergent competitor: Jülich (Schnepf/Vereecken), CPlantBox + DuMuX.

## What has already been established empirically

Tested on SBCL 2.2.9 against Quicklisp dist 2026-01-01:

- `femlisp` **loads in ~17 s and solves a 2D convection–diffusion–reaction model
  problem out of the box.** Last upstream update 2025-06-22. Modified BSD.
- It provides 37 `FL.*` packages including `FL.MESH` (hierarchical meshes),
  `FL.DISCRETIZATION`, `FL.MULTIGRID`, `FL.CDR` (convection–diffusion–reaction),
  `FL.ELLSYS`, `FL.NAVIER-STOKES`, `FL.PARALLEL`, `FL.STRATEGY` (adaptive),
  `FL.ALIEN` (FFI), `FL.LAPACK`.
- `magicl` and `lla` also load (dense linear algebra over BLAS/LAPACK).
- Dependencies pulled: alexandria, bordeaux-threads, cffi, cl-ppcre, closer-mop,
  lparallel, cl-mpi, trees, uiop. `libgd` is optional at load time.

**Consequence: the FEM and solver machinery does not need to be written.** The
missing layer is porous-media physics — variably saturated flow, constitutive
relations, multiphase, sorption, unsaturated effective diffusivity, and the
reactive-transport coupling. That is what Poros is.

## The crux question — answered (2026-10-05)

Femlisp is a **finite element** framework and ships **continuous Lagrange
elements only**. A search of the whole distribution finds no Raviart–Thomas,
BDM, Nédélec, discontinuous Galerkin or finite-volume discretisation, and no
conservative flux reconstruction. It does define the abstractions `<fe>`,
`<scalar-fe>`, `<vector-fe>` and `<vector-fe-discretization>`, so the extension
point for a new element family exists and is already used by its Navier–Stokes
and elasticity modules.

Global mass conservation was measured and passes to round-off (1e-16 to 1e-14),
but that was expected and is not sufficient: continuous Galerkin is globally
conservative by construction and locally conservative never. Since Darcy
velocity is −K∇h, a continuous head field gives a velocity that does not balance
cell to cell, and transport driven by it gains or loses solute at every face.

**Decision: implement lowest-order Raviart–Thomas in mixed form on Femlisp's
existing `<vector-fe>`, and contribute it upstream to Femlisp.** Full reasoning
and the measurements are in PHASE-0-DECISION.md.

## Working agreements

- **Upstream first.** If something general is missing from Femlisp, contribute it
  to Femlisp rather than building a private copy inside Poros. The test is
  whether the piece is useful to someone who has never heard of the rhizosphere:
  mixed finite elements are, van Genuchten is not. General numerics go upstream;
  porous-media physics stays here. Develop in Poros so the work is usable at
  once, upstream it once validated, and never block on review.
- **Validate against published benchmarks, not intuition.** Every solver phase
  has a named benchmark and a numeric gate before the next phase starts.
- **Every run leaves a log.** Capture the environment, the sha256 of the script
  and the output, as in `experiments/*.log`. A result that exists only in a
  terminal is not reproducible evidence.
- **Mass balance is reported for every transport run.** A run that does not
  conserve mass to stated tolerance is a failed run, not a result. Note that
  *global* conservation is cheap and *local* conservation is the one that
  matters — see PHASE-0-DECISION.md.
- **No duplicate code.** The previous incarnation of this project was a 39k-line
  org notebook with 192 function names defined 1,762 times. Not again.
- **Deterministic by construction.** Seed every stochastic run and record the
  seed in the output.
- **Say what is measured vs what is assumed.** Parameters without a literature
  source are flagged in the code, not quietly defaulted.
- Prior work worth reading before arguing about language choice: N. Neuss,
  *On Using Common Lisp for Scientific Computing*.

## What this project is not

- Not a cellular automaton. The predecessor was; the lattice-anisotropy study
  that came out of it is at github.com/atolgrp/Rhizo and is closed.
- Not a DuMuX competitor. It is a Common Lisp porous-media layer, and the
  comparison that matters is correctness against benchmarks, not feature count.
- Not expected to be finished quickly. See the timeline realism note in PLAN.md.
