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
| **EnergyBalance recorder** — per-region energy sidecar; energy math now lifted into the shared `EnergyBalanceKernel.h` (consumed by both this recorder and MPCO_Ladruno's `EnergyBalanceSource`, ADR D8) | Recorder | 33000 | `SRC/recorder/EnergyBalanceRecorder.{cpp,h}`, `SRC/recorder/EnergyBalanceKernel.h` | shipped | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#24](https://github.com/nmorabowen/OpenSees/pull/24) |
| **Explicit Noh–Bathe integrators** — `ExplicitBathe`, `ExplicitBatheLNVD` (hardened `dt_cr`, shared TU, DSYGV+β fix, `-lump`, MPI reduce) | Integrator | 33000 / 33002 | `SRC/analysis/integrator/ExplicitBathe.{cpp,h}`, `ExplicitBatheLNVD.{cpp,h}` | shipped | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#9](https://github.com/nmorabowen/OpenSees/pull/9) |
| **Queryable critical time step** — query/recompute `dt_cr` (tangent/periodic); fixes mass-aliasing bug | Integrator API | — | (within explicit integrator TU) | shipped | [#3](https://github.com/nmorabowen/OpenSees/pull/3) |
| **BezierTri6** — quadratic Bézier triangle element (Kadapa 2018) | Element | 33000 | `SRC/element/bezierTriangle/BezierTri6.{cpp,h}` | shipped | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#10](https://github.com/nmorabowen/OpenSees/pull/10) |
| **MPCO_Ladruno** — modular `.mpco` recorder fork (`mpcol`), node + element parity, multi-stage restaging; self-describing basis (consumes the `basisInfo`/`integrationPoints`/`integrationWeights` element probes), schema-v1 `MODEL/ELEMENTS/<name>` as a GROUP with `CONNECTIVITY` + `QUADRATURE`/{`GP_PARAM`,`GP_WEIGHT`}; first consumer = BezierTri6 (FAMILY=bernstein, barycentric GP); **standard-rule QUADRATURE by derivation** (`getStandardQuadrature` table → `NDIR`+`NUM_GP`+`QUADRATURE` for legacy GL quad/hex/tri/tet/line, no reader KeyError; higher-order quad9/tet10/hex20) + **belt-and-suspenders `GLOBAL_GP_COORDS`** (write-side `computeGlobalGP` basis evaluator incl. quad9/hex8/tet10 + hex20 serendipity, reference config, write-time round-trip oracle) + explicit `NDIR` (D4); **chunked extensible `[T×nIds×nComp]` time-series (ADR D3)** — `DATA` per result is one chunked+shuffle+deflate dataset growing one slab/step + `STEP[T]`/`TIME[T]` axes, replacing per-step `DATA/STEP_<k>` (kills 10⁵–10⁶-dataset blowup, compresses the time axis, O(1) slice reads); reader handles chunked-or-legacy transparently; **energy result type (ADR D8): whole-model `RESULTS/ON_DOMAIN/energyBalance` + per-region `ON_REGIONS` (new `ResultFamily::OnRegions`) via `mpcol::EnergyBalanceSource` calling the shared `EnergyBalanceKernel.h`; `-G energy <regionTag...>` verb; `MODEL/SETS/SET_<tag>` writer** | Recorder | 33001 | `SRC/recorder/MPCORecorderLadruno.{cpp,h}`, `MPCOL_{Sinks,ResultIO,NodeResults,ElementResults,DomainResults}.{cpp,h}`, `EnergyBalanceKernel.h` | shipped | [#8](https://github.com/nmorabowen/OpenSees/pull/8), [#12](https://github.com/nmorabowen/OpenSees/pull/12), [#14](https://github.com/nmorabowen/OpenSees/pull/14), [#16](https://github.com/nmorabowen/OpenSees/pull/16), [#18](https://github.com/nmorabowen/OpenSees/pull/18), [#23](https://github.com/nmorabowen/OpenSees/pull/23), [#24](https://github.com/nmorabowen/OpenSees/pull/24), [#29](https://github.com/nmorabowen/OpenSees/pull/29), [#32](https://github.com/nmorabowen/OpenSees/pull/32), [#33](https://github.com/nmorabowen/OpenSees/pull/33), [#36](https://github.com/nmorabowen/OpenSees/pull/36) |
| **CentralDifferenceLadruno** — explicit leap-frog central difference done right (correct first-step starter, built-in `dt_cr`, clean full-step velocity, βK guard); coupled mode dropped → use `NewmarkExplicit 0.5` | Integrator | 33003 | `SRC/analysis/integrator/CentralDifferenceLadruno.{cpp,h}` | per-TU compile-verified (full link pending MPCO_Ladruno blocker) | [#22](https://github.com/nmorabowen/OpenSees/pull/22) |
| **Ladruno brick element(s)** — our own higher-order hexahedral element(s) (planned), sibling to BezierTri6 on the solid side. Will implement the [[mpco_ladruno_element_contract]] (`basisInfo`/`integrationPoints`/`integrationWeights` self-declaration) so MPCO_Ladruno records geometry + GP coords with zero recorder edits; non-negative lumped mass for explicit dynamics, B-bar/assumed-strain against volumetric locking. Scope/order TBD (serendipity 20N vs Lagrange 27N vs enriched/Bézier). | Element | TBD (≥33000 band) | (plan: `Ladruno_implementation/`, TBD `SRC/element/`) | draft | — |

## Documentation / ADR PRs (no source change)

These shaped the design but shipped docs only — recorded for traceability:

- [#4](https://github.com/nmorabowen/OpenSees/pull/4) — ADR: explicit Noh–Bathe integrators + energy-balance recorder
- [#11](https://github.com/nmorabowen/OpenSees/pull/11) — ADR: reframe energy-recorder future work (apeGmsh is our viewer)
- [#13](https://github.com/nmorabowen/OpenSees/pull/13) — ADR: apeGmsh-native energy result type into MPCO_Ladruno (D8)
