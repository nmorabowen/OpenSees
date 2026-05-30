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

## In progress — branch `feature/mpco-quadrature-global-gp`
**Step 2+3:** standard-rule `QUADRATURE` (close reader KeyError on quad/hex/tri/tet
GL elements) + belt-and-suspenders `GLOBAL_GP_COORDS` + `NDIR`. Full higher-order
scope. Plan + verified GP-ordering spec in [[06_quadrature_global_gp_plan]] and the
memory. **Step A DONE:** `getStandardQuadrature(...)` in `MPCOL_ElementResults.h`
(compiles+links clean). **Resume at Step B.**

## Resume (next session)
`"continue Step B of MPCO_Ladruno standard-rule QUADRATURE"` —
1. **Step B:** add `computeGlobalGP(...)` (linear quad/hex tensor in each element's
   own node order + tri/tet bary + line) to `MPCOL_ElementResults.h`; rewire
   `MPCORecorderLadruno.cpp::writeModelElements` (~606–639) to resolve
   `(num_gp,ndir,gp_param,gp_weight)` via custom→`getStandardQuadrature`→none, write
   `NDIR`/`NUM_GP`/`QUADRATURE` for **all** resolved rules, and compute+write
   `GLOBAL_GP_COORDS[nElem×nGP*ndim]` (2-D, reshape doc'd; reuse the CONNECTIVITY
   loop's `node->getCrds()`); incremental build; h5py-inspect a quad+brick+tri model.
2. **Step C:** Python — `ladruno_basis.py` (bary/tri/tet), `ladruno_format.py`
   (read `NDIR`; `GLOBAL_GP_COORDS` optional; **QUADRATURE tolerant — warn not
   KeyError**); new `standard_quad_{model,check}.py` with the round-trip oracle.
3. **Step D:** higher-order (9N/27N Lagrange + 20N serendipity, C++ + Python).
4. **Step E:** fixture + harness self-tests + remaining rules (27N hex, tri 3/4-pt,
   tet GL2 — read sources or rely on the tolerant reader).
5. **Step F:** full regression gate (80/80·96/96·144/144·72/72) → merge PR → `ladruno`.

## Build / run recipe
- Worktree build tree is warm. Incremental: `cmd /c "Ladruno_scripts\setup_env.bat
  && cmake --build build\build\Release --target OpenSeesPy && copy /y
  build\build\Release\OpenSeesPy.dll dist\bin\opensees.pyd"` via the **PowerShell
  tool** (~2 min; NOT `build.bat` — its reconfigure forces a 1807-obj full rebuild).
- Run models with BUILD python `…\pythoncore-3.12-64\python.exe` (script does
  `add_dll_directory` + `sys.path`→`<wt>\dist\bin`); validate/oracle with venv python
  (h5py). `LADRUNO_OPENSEES_QUIET=1`; `grep -a` (banner has binary bytes).
- `gh pr create` needs `--repo nmorabowen/OpenSees` (else targets upstream).
