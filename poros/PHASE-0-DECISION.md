# Phase 0 decision note

**Date:** 2026-10-05
**Question:** can Poros build on Femlisp's discretisation, or does it need its own?
**Answer:** build on Femlisp, and contribute the missing discretisation *to* Femlisp
rather than keeping it private.

Evidence is in `experiments/*.log`, each log carrying its environment and the
sha256 of the script that produced it. Everything below was measured on
SBCL 2.2.9, Quicklisp dist 2026-01-01, `femlisp-20250622-git`.

---

## What was measured

### 1. Femlisp runs, and it is being driven correctly

Loads from Quicklisp in ~17 s cold. Solves a 2D convection–diffusion–reaction
model problem. 37 `FL.*` packages, including `FL.MESH`, `FL.DISCRETIZATION`,
`FL.MULTIGRID`, `FL.CDR`, `FL.NAVIER-STOKES`, `FL.PARALLEL`, `FL.STRATEGY`.
(`phase0-probe.log`)

### 2. Convergence is correct — method of manufactured solutions

u = sin(πx)sin(πy) on the unit square, so ∫u = 4/π² = 0.405284734… exactly.
Comparing the computed integral to that number tests problem definition,
assembly, solve and `FE-INTEGRATE` together against a value that owes nothing
to the code under test. (`convergence.log`)

| order | levels | computed ∫u | abs error | rate |
|---|---|---|---|---|
| 1 | 3 | 0.383569401357 | 2.172e-2 | 2.18 |
| 1 | 4 | 0.400022450959 | 5.262e-3 | 2.04 |
| 1 | 5 | 0.403979298262 | 1.305e-3 | **2.01** |
| 2 | 3 | 0.405492534250 | 2.078e-4 | 3.82 |
| 2 | 4 | 0.405298026086 | 1.329e-5 | **3.97** |

Order 1 → h², order 2 → h⁴. Both are what theory predicts for this functional.
The discretisation is sound and the harness is trustworthy.

### 3. Global mass conservation — passes to round-off

Pure diffusion, no source, no-flux boundaries (`:dirichlet nil` leaves the
natural condition, which is zero normal flux). Initial condition
u₀ = 1 + ½cos(πx)cos(πy), whose integral over the unit square is exactly 1
because the cosine terms vanish. The L2 projection onto a Lagrange space
preserves that integral exactly — constants are in the test space — so any
drift is a failure of the time-stepping, not of representing u₀.
(`mass-balance.log`)

| order | levels | steps | total mass | rel error |
|---|---|---|---|---|
| 2 | 3 | 1 | 0.999999999999999 | 6.7e-16 |
| 2 | 3 | 50 | 0.999999999999995 | 4.7e-15 |
| 1 | 3 | 50 | 1.000000000000007 | 7.3e-15 |
| 2 | 4 | 50 | 0.999999999999935 | 6.5e-14 |
| 2 | 3 | 200 | 1.000000000000031 | 3.1e-14 |

Round-off across every configuration. **This was the expected result and it is
not sufficient.** Continuous Galerkin is globally conservative by construction:
summing all equations telescopes the interior fluxes. The test confirms the
time-stepping does not leak mass; it says nothing about local conservation.

### 4. The element inventory — this is what decides it

Femlisp defines the abstractions `<fe>`, `<scalar-fe>`, `<vector-fe>`,
`<fe-discretization>`, `<vector-fe-discretization>` — but **`lagrange.lisp` is
the only concrete element family in the tree**. A search of the entire
distribution, source and documentation, finds no Raviart–Thomas, BDM, Nédélec,
discontinuous Galerkin or finite-volume discretisation, and no conservative
flux reconstruction.

So Femlisp offers continuous Lagrange elements and nothing locally conservative.

---

## Why local conservation is non-negotiable here

Darcy velocity is **v = −K∇h**. With a continuous Lagrange head field, ∇h is
discontinuous across element faces, so the velocity reconstructed from it does
not balance cell by cell: what leaves one cell is not what enters its neighbour.
Driving transport with such a field injects or destroys solute at every face.
The error does not vanish under refinement — it is a property of the scheme, not
a truncation term — and it is precisely why mixed FEM and finite volume dominate
porous-media codes.

A globally conservative transport run above a non-conservative velocity field
can still report total mass to round-off while being locally wrong everywhere.
Measurement 3 would not catch that. This is the gap.

---

## Decision

**Implement a locally conservative discretisation, and offer it upstream to
Femlisp rather than keeping it inside Poros.**

Concretely: lowest-order Raviart–Thomas in mixed form, implemented as a new
concrete subclass of Femlisp's existing `<vector-fe>`. The extension point
already exists and is already used by the Navier–Stokes and elasticity modules,
so this fits the architecture rather than fighting it.

This is option (b) from PLAN.md, with the user's amendment: the general numerics
go upstream, not into a private layer.

### Division of labour

| Goes upstream to Femlisp | Stays in Poros |
|---|---|
| Raviart–Thomas / mixed elements | Richards equation |
| Conservative flux reconstruction | van Genuchten–Mualem relations |
| The `:output nil` bug (below) | Sorption isotherms, retardation |
| | Effective diffusivity θ-relations |
| | The signalling-network / percolation layer |

The test is whether a piece is useful to someone who has never heard of the
rhizosphere. Mixed elements are; van Genuchten is not.

### Why upstream rather than private

- It is the contribution that actually fills the ecosystem gap. A private RT
  implementation inside Poros helps one project; in Femlisp it helps everyone
  using Common Lisp for PDEs.
- It is more citable. "Contributed mixed finite elements to Femlisp" is a
  stronger claim than "wrote some elements for my own code".
- It reduces Poros to what it should be — a porous-media physics layer — rather
  than a half-finished FEM framework competing with the one it depends on.

The cost is honest: upstream work moves at the maintainer's pace, and a patch
may be rejected or substantially revised. Mitigation is to develop in Poros
first so the work is usable immediately, and upstream it once validated, rather
than blocking on review.

---

## First upstream contribution, available now

`:output nil` crashes `solve` with `The value NIL is not of type REAL`, thrown
deep in the solver with no indication of the culprit. `:output 0` is the correct
way to silence output. `nil` is the obvious thing to write and the error gives
no clue.

Small, self-contained, clearly a bug, and reproducible in four lines. A good
first patch — it opens a channel with the maintainer before anything large is
proposed.

---

## Next

1. Send the `:output nil` fix upstream.
2. **Quantify the gap before building the fix.** Solve a Darcy problem with
   Lagrange, reconstruct the velocity, and measure the flux imbalance across
   element faces. This turns "CG is not locally conservative" from a textbook
   statement into a number for this code, and becomes the regression test the
   RT implementation must beat. It is the honest completion of Phase 0.
3. Implement lowest-order RT on `<vector-fe>`; verify the same imbalance is at
   round-off.
4. Then Phase 1 — Richards — on a discretisation that conserves locally.
