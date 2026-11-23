# Poros

A porous-media and reactive-transport layer for Common Lisp, built on
[Femlisp](http://www.femlisp.org/).

**Status: shelved (2026-11-23).** The rhizobiome simulator is being built on
`rhizo-diffusion` instead: its Voronoi-cell diffusion is already a two-point-flux
finite-volume scheme, so it conserves mass locally, which is the property
Phase 0 found missing in Femlisp. See `rhizo-diffusion/PLAN.md`. Poros remains a
possible long-term Common Lisp framework; the Phase 0 findings below still stand.

**Phase 0 status.** Nothing works yet. The architecture question that decides
everything downstream — whether a locally mass-conservative scheme can be
expressed in Femlisp, or whether Poros needs its own finite-volume layer — is
still open and is the first thing being settled.

- **[PROMPT.md](PROMPT.md)** — what this is, why, and what has been established.
  Read this first. It doubles as a brief to paste into a fresh chat.
- **[PLAN.md](PLAN.md)** — staged plan with a validation gate at each phase.

## Why

Python has FiPy, PorePy and OpenPNM; C++ has DuMuX and OpenFOAM. Common Lisp has
Femlisp, a capable general finite-element framework, and nothing above it for
porous media. Poros is that layer: variably saturated flow, constitutive
relations, sorption, unsaturated effective diffusivity, and reactive transport.

It exists to support a model of root–microbe chemical signalling in soil, but is
being built to stand on its own.

## Verifying the foundation

```sh
sbcl --dynamic-space-size 4096 --script experiments/phase0-probe.lisp
```

Confirms Femlisp loads from Quicklisp and solves a 2D convection–diffusion–
reaction problem. Verified 2026-10-01 on SBCL 2.2.9 with Quicklisp dist
2026-01-01: loaded in ~17 s, solve completed.

## Licence

MIT.
