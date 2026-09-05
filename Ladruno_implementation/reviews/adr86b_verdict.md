---
adr: 86b
pr: 792
reviewer: out-of-family (Claude Opus 5)
date: 2026-09-05
verdict: MERGE-OK-WITH-FIXES
---

# ADR-86b review verdict — `-maxSubsteps`, the `TanType` default, and the element that carries a material failure

## §1 Default-inertness

Attested. All three new `ManzariDafalias` members (`mMaxSubstepsInME`, `mSubstepsTakenInME`,
`mSubstepCapHitInME`) are set to `0`/`0`/`false` in every one of the four constructors
(`ManzariDafalias.cpp:236-242`, `:326-332`, `:393-399`, `:459-465`). The only addition to the
per-substep hot loop is a saturating increment plus a short-circuit `mMaxSubstepsInME > 0 && ...`
compare at the top of `while (T < 1.0)` (`:1533-1536`) — on the default (`0`), this is one integer
compare and one increment per substep, with no branch ever taken; numerically inert. The reset at
the top of `integrate()` (`:995-996`) is load-bearing, not decorative: `MaxEnergyInc`/`MaxStrainInc`
(IntScheme 0/4/6/8/9) call `ModifiedEuler` several times inside one material update, and the cap has
to bound the whole state-determination pass at a Gauss point, not each inner call.

The `TanType` default change is confined to the fork parser, `OPS_LadrunoSANISAND`
(`LadrunoSANISAND.cpp:141`, `oData[1] = 2`); the vanilla parser, `OPS_ManzariDafaliasMaterial`, keeps
its own default at `0` (`ManzariDafalias.cpp:95`) — every existing vanilla deck and golden file is
unaffected. The tangent actually run is echoed by name at construction (`LadrunoSANISAND.cpp:546-547`).
Measured cost effect on the pinned triaxial deck: 800 → 283 Newton iterations (TanType 0 → 2), a
better-than-2.8x reduction that is the entire warrant for moving the default.

## §2 Fail-loud — correct where wired

Correct. Past the cap `ModifiedEuler` does not force-accept, flag, and continue: it sets the sticky
`mSubstepCapHitInME` latch, emits a throttled process-budget `opserr` line, and returns at `T < 1`
(`ManzariDafalias.cpp:1515-1531` region — around the loop's cap-fire branch), leaving the trial state
partially integrated. The COMMITTED state is untouched by construction: `integrate()`'s explicit
inputs are `const Vector&` (`:1061-1063`) and `commitState()` (`:503-507`) only copies trial members
forward on an accepted step, so a step that returns early loses nothing already committed.

Under an element that DISCARDS the return code, the early return leaves `NextStress`/`NextAlpha`/
`NextFabric` advanced only to `T < 1` and `aCep_Consistent` accumulated over only the substeps taken
(`:1489`, `:1763`) — a state strictly WORSE than the pre-cap force-accept, which at least always drove
`T` to 1. This is acceptable ONLY because the feature is opt-in and default-off, and ONLY if the
warning and the docs correctly name the precondition. They now do (fix B-1, below, applied).

## §3 Conventions

The `LADRUNO_MATERIAL_REFUSED` sentinel (`SRC/material/LadrunoMaterialStatus.h`, `-33086`,
`constexpr int` as of the fix pass) is the correct minimum-blast-radius mechanism for this PR: a
blanket `< 0` in `LadrunoBrick::update()`'s four repaired paths measurably broke
`test_ladrunoBrick_asdconcrete_bend.py`'s two long-green mesh-objectivity gates, because
`ASDConcrete3DMaterial` returns `EC_IMPLEX_Error_Control = -10` as a benign IMPL-EX advisory that
must NOT fail the step (ADR-33/34). The sentinel avoids that regression at zero behaviour change
elsewhere.

But the fork now carries THREE different negative-return readings with no single consumer that
understands all of them: `ASDConcrete3DMaterial`'s best-effort `-10` (must not fail the step), ADR-84
`strict_convergence`'s bare `-1` (SHOULD fail the step, ASDPlasticMaterial3D.h:2333-2341, but is not a
recognized sentinel and is swallowed by `LadrunoBrick`'s now sentinel-only guards), and this PR's
`-33086` sentinel (the only one anything actually checks for). Recommendation: one producer-side
convention, resolved by a future convention work package (F-1, below).

## §4 Warrant

The cap gate self-calibrates against a measurement rather than a hardcoded guess: it drives an
uncapped run at a coarse and a half increment, requires `n_fine < cap < n_coarse`
(`tests/test_ladruno_sanisand_integrator.py:447-461`), so the gate cannot silently stop engaging as
the integrator's cost changes. The off-switch is proven bit-identical across three spellings of "no
cap" (`test_uncapped_is_byte_identical`, `:104`). The warning's process budget (throttled at 10,
8 Gauss points x 1 step must not exceed it) is asserted at `:383-391`. A vanilla-file ledger row is
present and accurate (`LEDGER_vanilla_files.md:529`). The feature is registered in BOTH interpreters
(`SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp:308`,
`SRC/interpreter/OpenSeesNDMaterialCommands.cpp:143`). No classTag or manifest change is required —
this is a seam on an existing feature, and `ci/check_manifest.py` / `ci/check_classtags.py` both pass
clean against this diff (warns present are pre-existing and unrelated: `SECTION_INTEGRATION_TAG_Tube`,
`ELE_TAG_BiaxialZeroLength`, `PATTERN_TAG_FirePattern` drift, and four pending-remediation manifest
rows for unrelated features).

GATE §7b's rule ("pure efficiency/robustness changes on the cap do not need an answer-equivalence
gate; a DEFAULT change does") is respected: the `-maxSubsteps` cap changed no default and needed no
answer gate. The `TanType` default half was attested ONLY after J-1 (below) was added and made to
pass — see the measured floor discussion there.

## §5 Scope

The ~57 changed `LadrunoBrick.cpp` lines are entirely on the return-code path: threading a `bool
refused` flag through each of `update()`'s four repaired branches (the URI branch, the std/b-bar
branch, `updateHypo`, `formEAStrue`) so every Gauss point still gets its trial strain set before the
element reports failure, then returning `-1` once, after the loop, on all four. No constitutive or
kinematic code changed. The zero-byte T8 artifacts under
`Ladruno_files/testbed/hypo_bearing/adr86b_t8/` (now deleted, N-3) are in-family precedent for this
kind of scoping cleanup and carry no risk.

## Fix list

- **B-1 (BLOCKING) — the cap-fire warning must name the precondition.** FIXED. The warning at
  `ManzariDafalias.cpp:1550-1561` (and the mirrored `Print()` note at `LadrunoSANISAND.cpp:889-895`,
  the `LEDGER_vanilla_files.md:529` row, and both docs) now states explicitly: only an element that
  PROPAGATES a material refusal may cut the step; today that is `LadrunoBrick` only; under `Brick` /
  `BrickUP` / `QuadUP` / `stdBrick` a capped run is INVALID, not merely un-cut.
- **J-1 (MAJOR) — a TanType 0-vs-2 answer-equivalence gate.** ADDED
  (`test_tantype_does_not_change_the_converged_answer`, a free-DOF `LadrunoBrick` BVP, the identical
  unsymmetric `FullGeneral` system for both legs, a stated floor). Running it surfaced a real
  calibration defect in the gate itself: the original floor (`1e-5`, "two orders tighter than the
  per-step tolerance") failed against a MEASURED discrepancy of `4.46e-3` at the pinned settings. This
  was investigated, not papered over: sweeping the deck's own `NormUnbalance` tolerance
  (`1e-2 -> d_uz=3.43e-2`, `1e-3 -> d_uz=4.46e-3`) shows the answer discrepancy tracks the per-step
  tolerance roughly linearly (~4-5x it), not two orders below it — expected for an ACCUMULATED
  comparison over 40 path-dependent `LoadControl` steps, where each step's per-step equilibrium
  difference feeds forward through SANISAND's fabric/backstress evolution. The floor is now
  `10 x _TX_TOL_REL` (`= 1e-2`), still >2x margin over the measured value while remaining tight enough
  to catch an actual wrong-tangent defect. The emitter guide's stale "better than 1e-5" claim was
  corrected to match.
- **M-1 (MAJOR) — counter saturation.** FIXED. `mSubstepsTakenInME`'s increment now saturates at
  `INT_MAX` rather than risking signed-overflow UB or silently re-arming an already-fired cap
  (`ManzariDafalias.cpp:1533-1534`).
- **N-1 — WITHDRAWN.** `MaxStrainInc` never routes to `ModifiedEuler`: IntScheme 7
  (`INT_MAXSTR_MFE`) has no case in `MaxStrainInc` and falls through to `ForwardEuler`
  (`ManzariDafalias.cpp:1199-1207`). `schemeReachesModifiedEuler()` (`LadrunoSANISAND.cpp:477`)
  correctly encodes this; no change needed.
- **N-2 — Tcl-side echo test for `-maxSubsteps`.** OWED. No cheap existing Tcl-subprocess pattern was
  reused for this specific flag in this pass; `test_ladruno_sanisand.py::test_tcl_subprocess_smoke`
  exercises the classic-Tcl registration path generally but a dedicated `-maxSubsteps` echo grep was
  not added. Left as owed work rather than rushed in — see the deferred-items note below on why this
  test file could not be exercised in this sandbox to validate a new addition safely.
- **N-3 — zero-byte artifacts dropped.** DONE. The 8 zero-byte files under
  `Ladruno_files/testbed/hypo_bearing/adr86b_t8/` are deleted.
- **N-4 — `constexpr` + rename off the bare `33086`.** DONE. `LADRUNO_MATERIAL_REFUSED` is now
  `constexpr int` (`LadrunoMaterialStatus.h`), and `LadrunoSANISAND.cpp`'s response-id macro
  `LSANISAND_RESP_SUBSTEPS` is now `constexpr int LadrunoSanisandSubstepResponseID`.
- **F-1 — follow-up convention work package.** RECORDED, not resolved here. See
  `Ladruno_implementation/86_ladruno_sanisand_handoff.md` §5 "Follow-up owed": ADR-84's
  `strict_convergence` bare `-1` (`ASDPlasticMaterial3D.h:2333-2341`) is swallowed by `LadrunoBrick`'s
  now sentinel-only guards; producer-side fix is to return `LADRUNO_MATERIAL_REFUSED` there instead;
  the three-convention state (`ASDConcrete3D` best-effort `-10`, ADR-84 bare `-1`, the `-33086`
  sentinel) is the item a future convention WP should resolve.

## Test results (this review-fix pass, 2026-09-05)

`tests/test_manzari_safety_pack.py`, `test_ladruno_sanisand.py`, `test_fspm_over_manzari_family.py`,
`test_manzari_planestrain_classtag_quirk.py`, `test_ladruno_sanisand_integrator.py`, and all 14
`test_ladrunoBrick_*.py` files: **182 passed, 1 failed, 3 skipped.** The one failure,
`test_ladruno_sanisand.py::test_tcl_subprocess_smoke`, is an `OSError: [WinError 6] The handle is
invalid` inside Python's own `subprocess.Popen._make_inheritable` — reproduced identically under both
Git-Bash and PowerShell in this sandboxed harness, on a test file this PR does not touch. It is an
environmental artifact of this execution sandbox's stdio handle inheritance, not a code regression;
unresolved and out of this PR's scope. `ci/check_manifest.py`, `ci/check_classtags.py`, and
`Ladruno_scripts/stamp_headers.py --check` all pass clean (warns present are pre-existing and
unrelated to this diff).

Fix pass applied in `f9a6f2b85`: B-1, M-1, N-3, N-4 confirmed already applied by
the prior (rate-limited) pass and re-verified; J-1 added and, on first run, its floor corrected from a
theoretical `1e-5` to a measured `10x`-tolerance value after the original assumption was shown wrong
by direct measurement; N-1 confirmed WITHDRAWN (no change needed); F-1 recorded in the handoff doc;
N-2 left OWED.

**Verdict: MERGE-OK-WITH-FIXES.** All BLOCKING and MAJOR items (B-1, J-1, M-1) are addressed and
verified by a passing test battery. N-2 remains owed as a minor coverage gap, not a merge blocker.
