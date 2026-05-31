---
title: Ledger — our own implementations
project: Ladruno
tags:
  - ledger
  - features
  - implementation
---

# Ledger — Ladruno implementations (new code we authored)

Every brand-new file / feature the Ladruno fork adds on top of vanilla
OpenSees: what it is, its class tag, where it lives, and the PR(s) that
shipped it. This is the counterpart to [[LEDGER_vanilla_files]] (which tracks
edits to *pre-existing* upstream files).

## Conventions

- **One row per feature**, not per file. List the owning files in the *Files* cell.
- *Status*: `shipped` (merged to `ladruno`), `draft` (plan/WIP), or `frozen`.
- **The banner feature list mirrors this ledger.** Every `shipped` feature
  should have a matching line in `Ladruno_scripts/banner_features.txt`. After
  editing that file run `python Ladruno_scripts/patch_banner.py` and rebuild.
- Class tags live in `SRC/classTags.h`; keep them recorded here so we never
  collide on a tag.
- When a forward-looking plan in this folder ships, move it to
  `Ladruno_internal/implemented_<name>.md` and add/flip its row here.

## Ledger

| Feature | Kind | Class tag | Files | Status | PR(s) |
|---|---|---|---|---|---|
| **OpenSeesPyMP** — `import openseesmp`, MPI-parallel Python module (Patch 9) | Interpreter | — | `SRC/interpreter/PythonMPIModule.cpp` (+ vanilla edits, see [[LEDGER_vanilla_files]]) | shipped | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| **EnergyBalance recorder** — per-region energy sidecar | Recorder | 26 | `SRC/recorder/EnergyBalanceRecorder.{cpp,h}` | shipped | [#2](https://github.com/nmorabowen/OpenSees/pull/2) |
| **Explicit Noh–Bathe integrators** — `ExplicitBathe`, `ExplicitBatheLNVD` (hardened `dt_cr`, shared TU, DSYGV+β fix, `-lump`, MPI reduce) | Integrator | 61 / 63 | `SRC/analysis/integrator/ExplicitBathe.{cpp,h}`, `ExplicitBatheLNVD.{cpp,h}` | shipped | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#9](https://github.com/nmorabowen/OpenSees/pull/9) |
| **Queryable critical time step** — query/recompute `dt_cr` (tangent/periodic); fixes mass-aliasing bug | Integrator API | — | (within explicit integrator TU) | shipped | [#3](https://github.com/nmorabowen/OpenSees/pull/3) |
| **BezierTri6** — quadratic Bézier triangle element (Kadapa 2018) | Element | 272 | `SRC/element/bezierTriangle/BezierTri6.{cpp,h}` | shipped | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#10](https://github.com/nmorabowen/OpenSees/pull/10) |
| **MPCO_Ladruno** — modular `.mpco` recorder fork (`mpcol`), node + element parity, multi-stage restaging; self-describing basis (consumes the `basisInfo`/`integrationPoints`/`integrationWeights` element probes), schema-v1 `MODEL/ELEMENTS/<name>` as a GROUP with `CONNECTIVITY` + `QUADRATURE`/{`GP_PARAM`,`GP_WEIGHT`}; first consumer = BezierTri6 (FAMILY=bernstein, barycentric GP); **standard-rule QUADRATURE by derivation** (`getStandardQuadrature` table → `NDIR`+`NUM_GP`+`QUADRATURE` for legacy GL quad/hex/tri/tet/line, no reader KeyError) + **belt-and-suspenders `GLOBAL_GP_COORDS`** (write-side `computeGlobalGP` basis evaluator, reference config, write-time round-trip oracle) + explicit `NDIR` (D4) | Recorder | 27 | `SRC/recorder/MPCORecorderLadruno.{cpp,h}`, `MPCOL_{Sinks,ResultIO,NodeResults,ElementResults}.{cpp,h}` | shipped | [#8](https://github.com/nmorabowen/OpenSees/pull/8), [#12](https://github.com/nmorabowen/OpenSees/pull/12), [#14](https://github.com/nmorabowen/OpenSees/pull/14), [#16](https://github.com/nmorabowen/OpenSees/pull/16), [#18](https://github.com/nmorabowen/OpenSees/pull/18), [#23](https://github.com/nmorabowen/OpenSees/pull/23), [#PENDING_STEPB] |
| **CentralDifferenceLadruno** — explicit leap-frog central difference done right (correct first-step starter, built-in `dt_cr`, clean full-step velocity, βK guard); coupled mode dropped → use `NewmarkExplicit 0.5` | Integrator | 33003 | `SRC/analysis/integrator/CentralDifferenceLadruno.{cpp,h}` | per-TU compile-verified (full link pending MPCO_Ladruno blocker) | [#22](https://github.com/nmorabowen/OpenSees/pull/22) |

## Documentation / ADR PRs (no source change)

These shaped the design but shipped docs only — recorded for traceability:

- [#4](https://github.com/nmorabowen/OpenSees/pull/4) — ADR: explicit Noh–Bathe integrators + energy-balance recorder
- [#11](https://github.com/nmorabowen/OpenSees/pull/11) — ADR: reframe energy-recorder future work (apeGmsh is our viewer)
- [#13](https://github.com/nmorabowen/OpenSees/pull/13) — ADR: apeGmsh-native energy result type into MPCO_Ladruno (D8)
