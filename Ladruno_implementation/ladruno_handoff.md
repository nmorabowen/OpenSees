---
title: Ladruno — session handoff (cold-resume map)
project: Ladruno
status: in-progress
owner: nmora
updated: 2026-06-01
tags:
  - handoff
  - recorder
  - element
  - finite-strain
---

# Ladruno — session handoff

Cold-resume map. Three tracks below: the **`ladruno` recorder** (original), the
**solid element / finite-strain stack** (LadrunoBrick + SolidTransformation +
LogStrain), and **Track 3 — LadrunoDispBeamColumn** (regularized disp-based frame +
embedded softening hinge; COMPLETE through the coupled biaxial hinge as of 2026-06-18,
see its §"CURRENT STATE"). Recorder deep detail in the `.claude` memory
`project_mpco_ladruno.md`; design in [[03_ladruno_recorder]] (orig ADR) +
[[07_adr_post_review_storage]] (post-review decisions); on-disk format in
[[ladruno_schema_v1]]; in-flight work in [[06_quadrature_global_gp_plan]]; gotchas in
[[LEDGER_quirks]].

---

# Track 2 — Solid element / finite strain (DONE, merged + CI-green 2026-06-01)

Deep detail in `.claude` memory `project_solid_transformation_wrapper.md` +
`project_finite_strain_wrapper.md`; plans [[09_ladruno_brick]],
[[solid_transformation_wrapper]], [[09_finite_strain_material_wrapper]].

## Shipped (all merged to `ladruno`, Zone-A + manifest gates green)
- **LadrunoBrick** (`ELE_TAG 33002`) — unified 8-node hex, `-formulation {std|bbar|uri|eas}`
  + `-hourglass {stiffness|physical|viscous}` (#69/#72/#75). std/bbar reproduce upstream
  `Brick`/`bbarBrick`; eas = SSPbrick bbar+condensed-EAS.
- **SolidTransformation** geometry-method layer — `-geom {linear|finite}` axis, orthogonal
  to `-formulation`. v1 linear (identity, #71) + **v3 finite (updated-Lagrangian, #76)**.
- **LogStrainNDMaterial** (`ND_TAG 33010`) — Hencky finite-strain adaptor (#70), broker +
  contract-lock (#73).

## The load-bearing fact (don't re-derive — it cost a 2nd build + an adversarial sweep)
The spatial constitutive modulus `c = (1/2J)[D:L:B]` is **NOT minor-symmetric in (k,l)**
(factor B eq.14.102 isn't k↔l symmetric), so it **cannot pass through a 6×6 Voigt
matrix** — `NDMaterial::getTangent()` (6×6) is a LOSSY projection. The finite element
gets the **full 4th-order tensor** via `FiniteStrainNDMaterial::getSpatialTangentTensor(double
c[3][3][3][3])` (+ kernel `logstrain_kernel::spatial_tangent_full`), then forms
`a_ijkl = c_ijkl − σ_il δ_jk` (the −σδ geometric term is element-owned) and assembles
`K_{(a,i)(b,k)} = ∫ (∂Nₐ/∂xⱼ) a_ijkl (∂N_b/∂xₗ) dv` with the FULL nodal gradients. This is
the ONLY split that yields a symmetric, FD-consistent tangent. Arbiter:
`tests/test_ladrunoBrick_finite.py` (FD-tangent + symmetry + homogeneous patch + objectivity
+ reduce-to-linear, 6/6).

## Resume (next session) — finite follow-ups (all fresh work, no loose ends)
- **corot (v2):** `SolidTransformationCorot` (Higham polar R + EICR K_geo + rigid-rotation
  patch test). Will need dedicated core buffers (the v1 in-place identity `globalize` aliasing
  won't survive Tᵀ K T).
- **plastic inner materials** through `-geom finite` (seam-3 state protocol; v3 LogStrain is
  elastic-exact only).
- **bbar+finite = F-bar**, and **`-geom finite` + damping** (currently rejected at the factory).
- **Banner:** `Ladruno_scripts/banner_features.txt` still lacks a `-geom finite` / LogStrain
  line (CLAUDE.md banner-sync rule) — a small follow-up PR (edit + `patch_banner.py` + rebuild).

## Process lessons from the v3 saga (avoid repeating)
- After resolving a merge, **PUSH it** and verify `git merge-base --is-ancestor origin/ladruno
  origin/<branch>` — a locally-resolved-but-unpushed merge left #76 showing CONFLICTING for two
  sessions.
- GitHub paused **auto-triggered** Actions runs (manual `workflow_dispatch` still worked) —
  likely an account Actions usage cap; check billing if PR checks stop populating.
- `ladruno` has **no branch protection** → a red/missing CI badge does NOT block merge; a real
  merge conflict does.

---

# Track 1 — `ladruno` recorder

Cold-resume map for the `ladruno` recorder (the original subject of this doc).

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

---

# Track 3 — LadrunoDispBeamColumn (regularized disp-based frame) (DONE through ADR 34 PR-4a, merged + 88/88 green 2026-06-18)

Design in [[32_ladruno_dispbeamcolumn_regularization_adr]] (2D), [[33_ladruno_dispbeamcolumn3d_hinge_adr]] (3D hinge), [[34_ladruno_cohesive_hinge_biaxial_adr]] (coupled material); gotchas in [[LEDGER_quirks]].

## CURRENT STATE (2026-06-18): the embedded-hinge feature is COMPLETE end-to-end.
All of Stage 0–1 + Stage 2 (2D hinge, 3D strong-axis, 3D biaxial block-diagonal, 3D coupled
interaction-surface) is shipped to `ladruno`, **88/88** (2D+3D element + scalar/biaxial/coupled hinge
+ cohesive material), CI green (Zone-A Ubuntu + classTag/manifest gates). 0 open PRs.
Session PRs: #274 (PR-3b biaxial), #275 (PR-4a coupled material), #276 (adversarial-review
dispositions), #280 (manifest dispatch.tcl fix). Then a 32-agent adversarial review of #274+#275 —
its dispositions are in [[34_ladruno_cohesive_hinge_biaxial_adr]] §"Adversarial review dispositions"
(net: 2 doc nits applied, the 2 "must-fix"es were a false positive + a convention conflict; see quirks).

## The element's four hinge modes (one ndm-dispatched `LadrunoDispBeamColumn` command)
- 2D `-hinge $matTag` (scalar α); 3D `-hinge $matTag` (strong-axis α_z, rank-1 condensation);
- 3D `-hinge $mZ -hingeY $mY` (BIAXIAL block-diagonal: two scalar laws, true coupled 2×2 K_αα via the
  bulk ks(MZ,MY), eigenvalue-floored inverse, rank-2 condensation);
- 3D `-hingeBiaxial $ndMatTag` (COUPLED: one `LadrunoCohesiveHingeBiaxial` order-2 NDMaterial carrying
  the Mz–My interaction ellipse; the full 2×2 cohesive tangent enters K_αα). Exclusivity:
  `-hingeY` requires `-hinge`; `-hingeBiaxial` excludes `-hinge`/`-hingeY`; any hinge excludes `-nl`.

## Shipped (all merged to `ladruno`; PRs #251 #254 #255 #258 #260 + review-fix; Stage-2 cohesive material #264; 2D hinge #267/#269; 3D PR-3a #271; PR-3b #274; PR-4a #275)
- **LadrunoDispBeamColumn2d** (`ELE_TAG 33013`) + **LadrunoDispBeamColumn3d** (`ELE_TAG 33014`) —
  displacement-based fiber frame, clones of `DispBeamColumn2d/3d`. ONE ndm-dispatched command
  `LadrunoDispBeamColumn` (ndm2/ndf3 → 2D, ndm3/ndf6 → 3D). Reduces to stock bit-identically.
- **Tier-1 `lch` channel** — `current_section_lch = wt[i]*crdTransf->getInitialLength()` (REFERENCE
  length) set INSIDE the `update()` IP loop right before `setTrialSectionDeformation`;
  `getCharacteristicLength()` override; `-lch {ip|element|<value>}` (`isfinite`-guarded). Makes
  crack-band/auto-reg materials (ASDConcrete1D `-autoRegularization`, ASDSteel1D, LadrunoUniaxialJ2+
  Lemaitre) regularize over the localizing IP, not the whole element — fixes the `LEDGER_quirks §59`
  pathology stock `DispBeamColumn` has. Mesh-objectivity + correct-band validated.
- **Large displacement** — via the existing Corotational `CrdTransf` (validated == stock); no
  element-side geometric code.
- **`-nl`** — ½θ² (2D) / ½(θz²+θy²) (3D) bowing strain (`DispBeamColumnNL2d/NL3d`). Default linear.

## Load-bearing facts (don't re-derive)
- **3-SITE registration** for a Ladruno element: `classTags.h` + `FEM_ObjectBrokerAllClasses.cpp` +
  `interpreter/OpenSeesElementCommands.cpp` `functionMap` (OpenSeesPy) **AND**
  `element/TclElementCommands.cpp` `ladrunoElementTable` (the standalone Tcl `OpenSees.exe`). Missing
  the last → builds/links clean but `element ... not known` ONLY in the exe. Enumerate sites:
  `grep -rln "OPS_LadrunoIMKBeam" SRC/ | grep -v ladrunoIMKBeam/`.
- The `lch` assignment MUST stay inside the IP loop (once-only material latch). REFERENCE length.
- 3D `getTangentStiff` INLINES kb+q (unlike 2D which calls `getBasicStiff`); the `-nl` path overwrites
  kb via `getBasicStiff` + adds bowing to q. The damping stiffness-multiplier was MOVED from 3D
  `getBasicStiff` → `getInitialStiff` (mirroring 2D) so the `-nl`/linear tangents stay consistent
  under `-damp`.
- Tier-1 only helps `lch`-CONSUMING materials; non-regularizing ones (Concrete02) stay mesh-dependent →
  that is what Tier-2 addresses.

## Resume (next session) — recommended NEXT, in priority order
1. **Consistent off-radial onset tangent for the coupled law (HIGHEST VALUE).** The documented v1
   limitation: `LadrunoCohesiveHingeBiaxial`'s 2×2 tangent is **frozen-mix** (exact only on radial
   jump paths) and sign-discontinuous at the elliptical onset (r=1). On deeply weak-axis-dominant
   paths driven by INDIRECT control (a non-radial jump path), the inner Newton transiently misses at
   onset → the element falls back to best-effort + warn (see quirk). Fix: a consistent tangent that
   carries ∂(mode-mix)/∂α (or a line-searched inner solve that re-evaluates the residual). This is
   ALSO what would let the strict ADR-33 "cut the global step on inner non-convergence" behavior hold
   robustly (currently best-effort, deliberately — see quirk). Without it, prescribe rotations
   (radial control) for deep coupled softening.
2. **B-K / power-law mode-mix** for `Gf_mix(w)` (currently linear `w_z Gf_z + w_y Gf_y` in
   `effectiveLaw`); reduces to the per-axis Gf on the pure axes regardless, so this is a refinement.
3. **Torsional jump `α_t` (`-torsion`)** — a 3rd condensed channel; needs a torsional cohesive-law
   concept (fiber rupture in torsion is unstandardized — research-y, lower confidence).
4. **`-nl` + hinge cross-terms** — lift the v1 parser ban (2D+3D); the ½θ² bowing couples the jump
   into the axial channel, so the condensation gets cross-terms (deferred algebra).
5. **`ladruno_drive` snap-back collapse gate (2D)** — blocked on the RESERVED dissipation arc-length
   [[22_ladruno_dissipation_arclength_adr]] (or trace with `LadrunoIndirectControl` 33006).
- Low-sev nice-to-have: 3D `getTangentStiff` recomputes the linear kb then discards it when `-nl`
  (perf only; guard `if(!nlGeom)`).

## Load-bearing facts (don't re-derive / don't "fix" back)
- **Registration sites.** A Ladruno ELEMENT needs 4 sites (`classTags.h`,
  `FEM_ObjectBrokerAllClasses.cpp`, `interpreter/OpenSeesElementCommands.cpp` functionMap,
  `element/TclElementCommands.cpp` ladrunoElementTable). A Ladruno **nDMaterial** (the coupled hinge)
  needs 4: `classTags.h`, broker `getNewNDMaterial` case+include, `OpenSeesNDMaterialCommands.cpp`
  map+extern, `material/nD/CMakeLists.txt`. **PLUS a `Ladruno_implementation/testbed/manifest.yaml`
  row** for every new classTag — the `classTag + manifest` CI gate (`ci/check_manifest.py`) goes RED
  without it (it failed after #275; cf. quirks "manifest row"). Mirror an existing ND row.
- **Best-effort inner Newton (quirk — DON'T return -1).** `solveHingeJump`/`solveHingeJumpBiaxial`
  return 0 + warn on maxIter, NOT a failure code. Returning -1 (to "cut the step" per the ADR's literal
  text) makes softening fragile and conflicts with OpenSees convention (`ForceBeamColumn`) — it broke
  4 coupled tests at the onset. The global `NormDispIncr` test is the arbiter. See [[LEDGER_quirks]].
- **Eigenvalue floor is DIAGONAL-only (quirk).** The 2×2 condensation floor scales by
  `|bulk_zz|+|bulk_yy|+|Czz|+|Cyy|` — folding in the off-diagonal `Czy` over-floors near onset and
  destabilizes the radial inner Newton (broke 4 tests when tried). See [[LEDGER_quirks]].
- **Deep biaxial/coupled softening test recipe: PRESCRIBE the rotations.** A free DOF under a dead
  moment SNAPS at hinge activation (a single-element control artifact, not a tangent bug — the strong
  axis drives to full break cleanly either way). The robust driver is `sp`-prescribed nodal rotations
  ramped via LoadControl (`_prescribe` helper in the test files); for finite-rotation gates prescribe
  skew transverse displacements (member rotates + hinge opens, no dead-load snap). Corotational coupled
  must stay strong-axis-DOMINANT (weak-axis-dominant indirect paths hit the unstable onset, item #1).
- **FE_Datastore ID collision.** Two same-size IDs sent with the same dbTag+commitTag COLLIDE
  (key = dbTag/commitTag/size); the second overwrites the first. The Z/Y hinge metadata is sent as ONE
  combined ID (size 2 strong-only / size 4 biaxial); the coupled nD material's classTag/dbTag ride in
  the element `data` Vector (23→26) to dodge it. (Caught by the biaxial DB-roundtrip gate.)
- The Stage-0/1 facts still hold: `lch` assignment stays INSIDE the IP loop (REFERENCE length); 3D
  `getTangentStiff` INLINES kb+q; the `-damp` stiffness-multiplier lives in `getInitialStiff`.

## Build / run recipe
- `Ladruno_scripts\build.bat OpenSeesPy OpenSees`. Editing `classTags.h` forces a wide recompile; a
  new `add_subdirectory`/source needs the build re-run once. Tests (PYTHONPATH/PATH = `dist\bin`):
  `python -m pytest tests/test_ladrunoDispBeamColumn2d_element.py tests/test_ladrunoDispBeamColumn3d_element.py tests/test_ladrunoDispBeamColumn2d_hinge.py tests/test_ladrunoDispBeamColumn3d_hinge.py tests/test_ladrunoDispBeamColumn3d_hinge_biaxial.py tests/test_ladrunoDispBeamColumn3d_hinge_coupled.py tests/test_ladrunoCohesiveHinge_material.py`
  (88/88). CI gates locally: `python ci/check_manifest.py && python ci/check_classtags.py`.
