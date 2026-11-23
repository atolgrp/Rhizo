# Build plan

Read [PROMPT.md](PROMPT.md) first. This is the staged plan, with a validation
gate at the end of every phase. The gates are the point: each one is a place
where the project can be stopped, redirected, or published, rather than drifting
into permanent infrastructure work.

## Timeline realism, stated up front

To the end of Phase 3 — a documented, tested, citable framework — is realistically
**6 to 12 months of part-time work**. It will not be ready for a PhD application
submitted this cycle. The artifact for that is the completed lattice-anisotropy
study already published at github.com/atolgrp/Rhizo. This project is what the
application can point at as *current direction*, with Phase 0 and ideally Phase 1
visibly underway.

Anyone who tells you a general porous-media framework is a few weekends has not
written one. DuMuX has roughly twenty years and dozens of contributors behind it.
Poros is viable because it is a **thin, well-validated layer on Femlisp**, not
because the problem is small.

---

## Phase 0 — Settle the discretisation question — **DONE (2026-10-05)**

**Outcome: option (b), amended — build on Femlisp, and contribute the missing
discretisation upstream to Femlisp rather than keeping it private.**

Measured: convergence correct (order 1 → h², order 2 → h⁴ on a manufactured
solution); global mass conservation to round-off (1e-16 … 1e-14); and Femlisp
ships **continuous Lagrange elements only** — no RT, BDM, Nédélec, DG or FV
anywhere in the distribution, and no conservative flux reconstruction. Since
Darcy velocity is −K∇h, a continuous head field yields a velocity that does not
balance cell to cell, so transport above it gains or loses solute at every face.

Full note with the numbers: **PHASE-0-DECISION.md**. Logs: `experiments/*.log`.

Remaining to close Phase 0 honestly: **quantify the flux imbalance** for a Darcy
problem discretised with Lagrange. That turns a textbook statement into a number
for this code and becomes the regression test the RT implementation must beat.

### Upstream policy

General numerics go to Femlisp; porous-media physics stays in Poros. The test is
whether a piece helps someone who has never heard of the rhizosphere.

| Upstream to Femlisp | Stays in Poros |
|---|---|
| Raviart–Thomas / mixed elements | Richards equation |
| Conservative flux reconstruction | van Genuchten–Mualem |
| `:output nil` crash (first patch, ready now) | Sorption, retardation, D_eff(θ) |
| | Signalling network / percolation |

Develop here so the work is usable immediately; upstream once validated; never
block on review.

---

## Phase 1 — Variably saturated flow (`poros/flow`)

**6–10 weeks.** The genuinely missing physics: Richards' equation.

- Constitutive relations: van Genuchten θ(h) and Mualem K(θ), with Brooks–Corey
  as an alternative. Hysteresis deferred to Phase 4 but the interface must not
  preclude it.
- Mass-conservative time discretisation. Use the mixed form — the head-based
  form is well known to lose mass badly on infiltration fronts.
- Newton iteration with line search and adaptive time stepping.
- Boundary conditions: Dirichlet, Neumann, seepage face, atmospheric.

**Gate — reproduce two published benchmarks:**
1. **Celia, Bouloutas & Zarba (1990)** infiltration into dry soil. This is *the*
   standard Richards test and exists precisely because naive schemes fail it.
2. A drainage/evaporation benchmark with a reported solution.

Required: mass balance error < 0.1% over the run, and the wetting-front position
matching published values within plotting accuracy. **If Celia cannot be
reproduced, do not proceed to Phase 2** — every transport result above a broken
flow solver is worthless.

---

## Phase 2 — Reactive transport (`poros/transport`)

**6–10 weeks.** Solute movement on the flow field from Phase 1.

- Advection–diffusion–dispersion with velocity from the Richards solve.
- Effective diffusivity as a function of water content. Millington–Quirk as
  baseline, **with the expectation of replacing it** — it performed worst of the
  tested tortuosity models in Chou, Hu & Dane (2012, 10.1029/2011WR011653).
  Make the relation pluggable from the start.
- Sorption: linear (Kd), Freundlich, Langmuir; retardation factor. **First-class,
  not an afterthought** — see the sorption risk in PROMPT.md.
- First-order and Monod degradation.
- Operator splitting (Strang) with a quasi-static flow field between events,
  which is what makes the slow-diffusion/fast-flooding timescale separation
  tractable.
- Stiff reaction sub-steps via a proper integrator rather than shrinking the
  global step.

**Gate:**
- Ogata–Banks analytic solution for 1D advection–dispersion reproduced to
  stated tolerance.
- Retardation verified against the analytic R = 1 + ρKd/θ.
- Mass balance reported and within tolerance on every test.
- A Péclet sweep showing the scheme does not oscillate or smear unacceptably
  across the 4-order range the science demands.

---

## Phase 3 — Release the framework

**~4 weeks. This is the publication.**

- Documentation: tutorial, worked examples, API reference.
- Test suite covering all gates above, runnable with one command.
- CI on SBCL, and at least one other implementation if feasible.
- `poros.asd` split into `poros/core`, `poros/flow`, `poros/transport`.
- Submit to Quicklisp.
- **Software paper — JOSS is the natural venue.** JOSS wants working, documented,
  tested, useful software rather than novelty, which fits exactly.

**Gate:** someone who is not you can install it and run the tutorial unaided.
Test this on an actual person before submitting.

---

## Phase 4 — The science

**The original question, now on your own framework.**

- Signalling-network layer: nodes, directed edges, time-varying graph, percolation
  and spanning-cluster detection.
- Parameterise for **AHLs first** — the only signal class with a complete enough
  parameter set.
- Sorption sensitivity sweep. **Decision point:** if effective range is more than
  ~2× as sensitive to Kd as to θ near θ_c, the narrative changes from
  "hydration-gated" to "sorption-and-hydration co-gated". That is a finding, not
  a failure, but it must be reported honestly either way.
- Wetting–drying hysteresis, which is where the novel prediction lives.
- Reproduce the Ziadne et al. (2026) cross-soil signalling distances, including
  the ordering (clay < loamy sand < sand) and the sterilisation effect.

**Gate for a 3D run:** the 1D/2D model reproduces the Ziadne ordering within
experimental error. Only then does the full 3D domain earn its compute.

---

## Phase 5 — 3D and root architecture

Only after Phase 4 produces a result worth scaling.

- 3D domains, 1D-radial-rhizosphere embedded in a coarse 3D grid (the Mai et al.
  2019 pattern — mass-conservative and the established approach).
- Root architecture, either imported from CPlantBox output or generated natively.
- Microbial agents in the near-root shell only. At 10⁸–10⁹ cells per gram,
  individual agents are tractable nowhere else, and that limitation should be a
  deliberate modelling choice rather than a discovered surprise.

---

## Standing risks

| Risk | Mitigation |
|---|---|
| Infrastructure forever, no science | Hard gates; Phase 3 publishes before Phase 4 starts |
| Femlisp bus factor is 1 (Neuss) | BSD-licensed and vendorable; Phase 0 establishes how deep the coupling goes |
| Small Common Lisp audience | Accepted knowingly. JOSS + clean benchmarks is the best available answer |
| Sorption dominates, hypothesis weakens | Built in from Phase 2; explicit decision threshold in Phase 4 |
| Scooped by Jülich/Schnepf | Their work is 3D continuum; the percolating-graph framing is not theirs. Speed matters more than secrecy — publish Phase 3 |
| No measured D or Kd for strigolactones | Lead with AHLs; strigolactones via explicit sensitivity analysis over plausible ranges |
