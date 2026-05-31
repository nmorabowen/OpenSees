---
title: MPCO_Ladruno — session handoff (cold-resume map)
project: Ladruno
status: in-progress
owner: nmora
updated: 2026-05-30
tags:
  - handoff
  - recorder
---

# MPCO_Ladruno — session handoff

Cold-resume map for the `mpcoLadruno` recorder. Deep detail in the `.claude`
memory `project_mpco_ladruno.md`; design in [[03_mpco_ladruno]] (orig ADR) +
[[07_adr_post_review_storage]] (post-review decisions); on-disk format in
[[mpco_ladruno_schema_v1]]; in-flight work in [[06_quadrature_global_gp_plan]];
gotchas in [[LEDGER_quirks]].

## What it is
A sibling recorder (`recorder mpcoLadruno` → `.ladruno` HDF5), apeGmsh-native,
forked from the frozen `MPCORecorder` (untouched) and wrapped in `namespace mpcol`
(ODR-clean). Splits `ResultSource` (what to compute) from `ResultSink` (how to
persist). Value channels reproduce the frozen recorder to 1e-12.

## Merged on `ladruno` (done + verified)
- Recorder Phases 1–3 (#8); multi-stage restaging fix C1/C2 (#12); element-result
  parity gate (#14); HDF5-DIAG noise fix (#16); BezierTri6 wiring +
  **bucket-as-GROUP + `basisInfo` self-declaration probe + QUADRATURE for custom
  rules** (#18).
- **Verified (fresh build):** nodal 80/80, multistage 108/108 (2 stages),
  element quad 96/96 + beam 144/144, bezier 72/72, `SECTION_ASSIGNMENTS`
  byte-identical to frozen, zero HDF5-DIAG. Harness pytest green.

## Architecture review (done) → decisions
6-lens adversarial review = **proceed-with-changes, do NOT pivot the format**
(see [[07_adr_post_review_storage]]). Architect decisions locked:
**D2** belt-and-suspenders `GLOBAL_GP_COORDS` (+ write-time oracle); **D3** chunked
extensible `[T×nIds×nComp]` time-series replacing per-step `DATA/STEP_<k>` (do
early); **D4** explicit `NDIR`; **D5** shared validator + hold `FORMAT_VERSION=1`.

## DONE — Steps A + B (merged/PR'd)
**Step A (PR #23, merged):** `getStandardQuadrature(...)` table in
`MPCOL_ElementResults.h`. **Step B (PR [#29](https://github.com/nmorabowen/OpenSees/pull/29),
branch `feature/mpco-step-b-global-gp` off current ladruno):**
- `computeGlobalGP(...)` write-side basis evaluator (line2/quad4/quad9/tri3/tet4/hex8,
  node orders verified vs element sources) in `MPCOL_ElementResults.h`.
- `writeModelElements` rewired: custom→`getStandardQuadrature`→none; writes
  `NDIR`+`NUM_GP`+`QUADRATURE` for ALL resolved rules + `GLOBAL_GP_COORDS`
  `[nElem×nGP*ndim]` (2-D).
- `ladruno_format.py`: QUADRATURE-tolerant (warn), NDIR-authoritative (D4),
  GLOBAL_GP_COORDS-aware.
- New gate `standard_quad_{model,check}.py` with **write-time round-trip oracle**
  (ALL PASS ≤1e-12). No regression: 80/80·96/96·144/144·72/72·108/108·pytest 10/10.

## Resume (next session) — Step D/E
`"continue Step D of MPCO_Ladruno: higher-order GLOBAL_GP_COORDS"` —
1. **Step D — higher-order.** Add 9N quad / 20N serendipity / 27N Lagrange hex,
   10N tet, 6N Lagrange tri shape functions to BOTH `computeGlobalGP` (C++) and an
   importable `ladruno_basis.py` oracle (bary/tri/tet + order-2/serendipity). Node
   orderings MUST be read from each element source (Twenty_Node_Brick uses external
   `shp3d.h`/`shp3dv.h`). These are ungateable until their rules are tabulated (Step E),
   so land D+E together.
2. **Step E — remaining rules + fixtures.** Tabulate Hex_GL_3 (27pt — which 27-node
   brick? none currently maps), Tri_GL_2/2B/2C (shells, e.g. ShellDKGT), Tet_GL_2
   (read TenNodeTetrahedron.cpp) in `getStandardQuadrature`; add `NDIR`+`GLOBAL_GP_COORDS`
   +a simplex group to `make_synthetic.py`; extend the pytest harness.
3. **Step F (D3, separate later step) — chunked time-series** `[T×nIds×nComp]`
   replacing per-step `DATA/STEP_<k>` (ADR D3). Then KIND/LOCAL_AXES, shared validator,
   parallel.

## Build / run recipe
- Worktree build tree is warm. Incremental: `cmd /c "Ladruno_scripts\setup_env.bat
  && cmake --build build\build\Release --target OpenSeesPy && copy /y
  build\build\Release\OpenSeesPy.dll dist\bin\opensees.pyd"` via the **PowerShell
  tool** (~2 min; NOT `build.bat` — its reconfigure forces a 1807-obj full rebuild).
- Run models with BUILD python `…\pythoncore-3.12-64\python.exe` (script does
  `add_dll_directory` + `sys.path`→`<wt>\dist\bin`); validate/oracle with venv python
  (h5py). `LADRUNO_OPENSEES_QUIET=1`; `grep -a` (banner has binary bytes).
- `gh pr create` needs `--repo nmorabowen/OpenSees` (else targets upstream).
