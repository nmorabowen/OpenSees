---
title: Ladruno — session handoff (cold-resume map)
project: Ladruno
status: in-progress
owner: nmora
updated: 2026-05-30
tags:
  - handoff
  - recorder
---

# Ladruno — session handoff

Cold-resume map for the `ladruno` recorder. Deep detail in the `.claude`
memory `project_mpco_ladruno.md`; design in [[03_ladruno_recorder]] (orig ADR) +
[[07_adr_post_review_storage]] (post-review decisions); on-disk format in
[[ladruno_schema_v1]]; in-flight work in [[06_quadrature_global_gp_plan]];
gotchas in [[LEDGER_quirks]].

## What it is
A sibling recorder (`recorder ladruno` → `.ladruno` HDF5), apeGmsh-native,
forked from the frozen `MPCORecorder` (untouched) and wrapped in `namespace ladruno`
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
`Ladruno_ElementResults.h`. **Step B (PR [#29](https://github.com/nmorabowen/OpenSees/pull/29),
branch `feature/mpco-step-b-global-gp` off current ladruno):**
- `computeGlobalGP(...)` write-side basis evaluator (line2/quad4/quad9/tri3/tet4/hex8,
  node orders verified vs element sources) in `Ladruno_ElementResults.h`.
- `writeModelElements` rewired: custom→`getStandardQuadrature`→none; writes
  `NDIR`+`NUM_GP`+`QUADRATURE` for ALL resolved rules + `GLOBAL_GP_COORDS`
  `[nElem×nGP*ndim]` (2-D).
- `ladruno_format.py`: QUADRATURE-tolerant (warn), NDIR-authoritative (D4),
  GLOBAL_GP_COORDS-aware.
- New gate `standard_quad_{model,check}.py` with **write-time round-trip oracle**
  (ALL PASS ≤1e-12). No regression: 80/80·96/96·144/144·72/72·108/108·pytest 10/10.

## Step D DONE — PR #32 (quad9+tet10) + PR #33 (hex20)
Higher-order GLOBAL_GP_COORDS for all three source-verifiable elements:
- **quad9 (NineNodeQuad, Quad_GL_3)** — gated (rule+shape fn already shipped Steps A/B). [#32]
- **tet10 (TenNodeTetrahedron, Tet_GL_2 4-pt α/β)** — NEW `Tet_GL_2` rule + tet10 shape
  fn (node-8/9 swap) in `computeGlobalGP`. Round-trip 1.1e-16. [#32]
- **hex20 (Twenty_Node_Brick, Hex_GL_3 27-pt)** — traced `shp3dv` `brcshl`: GP order
  `b·(2·RA,2·SA,2·TA)` over the serendipity node pattern = element `materialPointers[L]`
  order; NEW `Hexahedron_GaussLegendre_3` rule + 20N serendipity basis in `computeGlobalGP`.
  Round-trip 2.2e-16. [#33]
- Gate `standard_quad_{model,check}.py` covers quad4/tri3/hex8/quad9/tet10/hex20; all
  CONFORMANT; no regression (80/80·96/96·144/144·72/72·108/108·pytest 10/10).

## D3 chunked time-series DONE — PR #36
`StreamingSink` now writes one chunked+shuffle+deflate `DATA[T×nIds×nComp]` dataset per
result + `STEP[T]`/`TIME[T]` axes (was per-step `DATA/STEP_<k>`). New
`Ladruno_Hdf5.h` `createTimeSeries3d`/`appendSlab3d`/`appendDouble1d`/`appendInt1d`. Reader
`ladruno_format.iter_step_slices` reads chunked-or-legacy transparently → chunked
`.ladruno` still diffs 1e-12 vs per-step `.mpco`. `make_synthetic.py` emits chunked.
Full regression green.

## KIND + LOCAL_AXES DONE — PR #38
- **KIND**: `-kind transient|static|eigen` option (default static), auto-`eigen` when a
  modal result is requested. Written at `MODEL_STAGE/KIND`.
- **LOCAL_AXES**: `writeModelLocalAxes()` writes per-class
  `MODEL/LOCAL_AXES/<classTag>-<name>/{ID,FRAME[E×4 quaternion]}` from each element's
  `"localAxes"` response (`quatFromMat`); never a silent identity. **ElasticBeam3d** wired
  (new `"localAxes"` response id 30 → `theCoordTransf->getLocalAxes`). Gate
  `localaxes_{model,check}.py` (diagonal beam → non-identity frame; quaternion→axes vs
  element axis). NOTE: frozen `quatFromMat` stores the **transpose** convention → axes are
  the ROWS of the reconstructed matrix (checker documents this).

## Resume (next session) — finish the recorder
1. **Remaining-beam localAxes** — replicate the ElasticBeam3d `"localAxes"` response (id →
   `crdTransf->getLocalAxes`, 9 packed cosines) on ElasticBeam2d, DispBeamColumn2d/3d,
   ForceBeamColumn2d/3d (identical pattern; 2D fills z=(0,0,1)). Each is a vanilla edit →
   LEDGER_vanilla_files row.
2. **Shared validator + CI round-trip oracle (D5)** → then **freeze `FORMAT_VERSION=1`**.
3. **Step E cleanup** — importable `ladruno_basis.py` (bary/tri/tet + quad9/tet10/hex20);
   `make_synthetic.py` NDIR/GLOBAL_GP_COORDS/simplex; `COLUMN_MAP/COMP_NAMES` authoritative.
4. **Tier 2 — parallel** (`sendSelf`/`recvSelf` + `part-N` + `Allreduce`); **Tier 3 —
   envelopes** (wire `EnvelopeSink`).
- **Minor:** geom-derived `ORDER` for 9N/10N/20N reports linear (NDIR authoritative).

## Build / run recipe
- Worktree build tree is warm. Incremental: `cmd /c "Ladruno_scripts\setup_env.bat
  && cmake --build build\build\Release --target OpenSeesPy && copy /y
  build\build\Release\OpenSeesPy.dll dist\bin\opensees.pyd"` via the **PowerShell
  tool** (~2 min; NOT `build.bat` — its reconfigure forces a 1807-obj full rebuild).
- Run models with BUILD python `…\pythoncore-3.12-64\python.exe` (script does
  `add_dll_directory` + `sys.path`→`<wt>\dist\bin`); validate/oracle with venv python
  (h5py). `LADRUNO_OPENSEES_QUIET=1`; `grep -a` (banner has binary bytes).
- `gh pr create` needs `--repo nmorabowen/OpenSees` (else targets upstream).
