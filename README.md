# Rhizo

Cellular automata modelling of the rhizobiome — signal exchange between plant
roots and the microbial community in the surrounding soil.

A cellular automaton on a square lattice does not spread isotropically. A
diffusion front advances faster along the lattice diagonals than along its axes,
so a process that ought to be radially symmetric comes out square-shaped. For a
model of chemical signalling in soil that is not cosmetic: the artifact
determines which microbes a signal is predicted to reach.

**How much of that bias survives if the automaton runs on a Voronoi tessellation
instead of on the lattice? It does not. Directional bias falls from 34% to under
3%.**

![Anisotropy result](voronoi-ca-isotropy/results/anisotropy.png)

The left panel is the mechanism: front speed against compass direction. On the
lattice it traces a square — the front reaches the diagonals about 41% sooner
than the axes. On the Voronoi tessellation it traces a circle.

![Swept region](voronoi-ca-isotropy/results/tessellations.png)

## What is here

### [`voronoi-ca-isotropy/`](voronoi-ca-isotropy/) — a self-contained experiment

A controlled measurement of lattice anisotropy and of how much a Voronoi
tessellation removes. Runs in about a minute, needs only SBCL, has no external
Lisp dependencies and no display:

```sh
cd voronoi-ca-isotropy
sbcl --script run.lisp 200 10
python3 analysis/analyse.py
```

Every trial is seeded and the seed is recorded in the output, so any row of the
results can be reproduced exactly. Full method, results and limitations are in
[its own README](voronoi-ca-isotropy/README.md).

The headline numbers, from a 200×200 lattice with 10 replicates per condition:

| condition | ring 1 | ring 2 |
|---|---|---|
| square lattice | 34.15% | 34.38% |
| Voronoi tessellation | 1.62% | 2.55% |

The lattice bias does not decay with distance, which is the reason it matters —
it is a systematic geometric error rather than noise, so a longer path will not
average it away. The lattice result also has an analytic prediction: under a
Moore neighbourhood a front travels in the Chebyshev metric, so the axial to
diagonal ratio should be √2. Measured, it is 1.4118 and 1.4151 at the two radii,
within 0.4% at both. That makes the measurement a test of the apparatus and not
only a report of its output.

### `rhizo.org` — exploratory notebook

An org-mode lab notebook on the multi-receiver signalling model: signal sources,
receivers, and the distance over which a signal is received. Literate and
exploratory, in the sense that it records a line of thinking rather than
presenting a finished program.

## Status and direction

The isotropy experiment above is complete and stands on its own. The wider
modelling question — signal exchange between roots and the rhizobiome under
realistic soil water conditions — is a larger problem in three dimensions, and
this work is a component of it rather than a foundation for it:

- The result here is 2D. It does not carry to 3D automatically: no
  three-dimensional Bravais lattice with a single-speed neighbour set has an
  isotropic fourth-rank tensor, which is why hexagonal geometry rescues the 2D
  case and nothing so cheap rescues the 3D one.
- The tessellation is not free. Regions average about 15.5 face-neighbours in 3D
  against 6 on a cubic grid, memory access becomes indirect, and the smallest
  cell sets the explicit stability limit. For pure diffusion an isotropic
  27-point Laplacian is far cheaper; for flow, D3Q19 lattice Boltzmann is
  rank-4 isotropic on a plain cubic grid.

So this repository quantifies what a Voronoi layer buys, which is the question
worth settling before deciding whether to pay for it.

## Licence

MIT — see [LICENSE](voronoi-ca-isotropy/LICENSE).
