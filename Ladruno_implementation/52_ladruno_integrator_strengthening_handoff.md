---
title: Integrator Strengthening — session handoff
project: Ladruno
status: active
owner: nmora
tags:
  - handoff
  - integrator
  - adr-52
updated: 2026-06-24
parent: "[[52_ladruno_integrator_strengthening_adr]]"
---

# ADR-52 — integrator strengthening: session handoff

Pick-up notes for the next agent continuing [[52_ladruno_integrator_strengthening_adr]]
(the strengthening roadmap, derived from the [[49_ladruno_integrator_study_workflow_adr]]
study + [[49a_integrator_scorecard_2026-06-23]] gap analysis).

## Status snapshot (2026-06-24)

**Shipped to `ladruno` (fork-only; W1-E3..W1-I1b zero net new vanilla, W3-I2 ledgered header-only vanilla):**

| Wave | What | PR(s) |
|---|---|---|
| **W1-E3** | SMS `dt_cr` honesty (`-cflAbort`/`-recompute` → report-only) + corrected the stale "consistent-SMS not parallel-safe" comment (it IS parallel-safe) | #394 |
| **W1-I1a** | Transient adaptive-Δt lane `robust_transient()` in `Ladruno_scripts/robust_drive.py` (convergence-driven; verdict `integrated`≠accurate) | #396 |
| **W2-E1** | Explicit bulk viscosity (`-bulkViscosity b1 b2` / `-bv`) on **all 3 fork continuum elements** — LadrunoBrick, LadrunoQuad, LadrunoCST | #399, #403 |
| **W1-I1b** | Half-increment-residual ACCURACY gate `robust_transient(error_gate=True, haftol=…)` + two read-only fork commands `ladrunoTrialResidualNorm`/`ladrunoSetNodeTrial` (registration-only vanilla, no classTag/header edit). Resolved the open question: NOT Python-only. Exact for elastic, approximate for inelastic/rate-dependent. | #407 |
| **W3-I2** | DDM sensitivity integrators: **`LadrunoHHT`** (33013, #413) + **`LadrunoGeneralizedAlpha`** (33014, #415). Subclass HHT / GeneralizedAlpha, reuse the algorithm, add the DDM seam (α-/αF,αM-weighted residual + extra `−K·(1−α)·dUₙ`; rest copies Newmark). Header-only ledgered vanilla edit of `HHT.h`/`GeneralizedAlpha.h` (private→protected + protected inline classTag ctor; `.cpp` byte-identical). Unblocks reliability/fragility/FORM on the damped integrators. | #413, #415 |
| **W1-E2** | **`ExplicitBathe*` 6→1 collapse.** One `ExplicitBathe` (33000) + orthogonal `-lnvd <alpha>` / `-sms <dtTarget>` / `-consistent` flags. The 5 retired `.cpp/.h` deleted; the 5 names+tags (33002/09/10/11/12) kept as deprecated aliases. Flag combo DERIVES the classTag (`tagForFlags`) → both brokers route the six tags to `ExplicitBathe::makeForBroker` → `recvSelf` fills the param superset (no new tag, none freed). The six were a 2-base × 3-SMS-mode lattice (NOT 6 flat siblings); the two bases' stepping was byte-identical, so only LNVD `formUnbalance`+SMS `domainChanged` are flag-gated. 6-lens adversarial review run. | #419 |

Docs/log PRs: #391 (ADRs 49/49a/52), #398, #401, #404.

**Remaining waves: NONE — ADR-52 COMPLETE** (W1-E2 #419 was the last). W3-I2 done (#413 + #415);
W3-I3 closed — NO-GO #410 (benchmark gate, see below). The deferred W2-E1 S3/S4/uri-ssp
follow-ups are optional, non-blocking.

## How to work this (the proven loop)

1. **Worktree:** `C:\Users\nmb\Documents\Github\OpenSees-integrators` (branch `integ-work`
   tracks `origin/ladruno`). Keep the main checkout free for the concurrent contact/AEM
   line. Per wave: `git checkout -b feat/... origin/ladruno`.
2. **PR → CI is the build + runtime gate.** Zone-A (`.github/workflows/ladruno.yml`)
   builds `opensees.so` and runs `pytest -m zone_a` in `tests/`. Add a
   `pytestmark = [pytest.mark.zone_a]` test (`from _testbed import ops`) — CI runs it
   against the fresh build. This replaces slow local fresh-worktree builds.
3. **Adversarial review for heavy/novel waves** (W1-E2, W1-I1b, W3-I2): run the
   `Workflow` multi-lens review (see the W2-E1 review run for the template), apply
   confirmed blockers, then merge. A faithful *port* of an already-reviewed change
   (e.g. W2-E1 Quad/CST) can skip the full review and rely on the `zone_a` runtime gate.
4. **Merge:** `gh pr merge <n> --squash` once Zone-A is green. Poll with a background
   `gh pr checks` loop (Zone-A ~9-10 min) — do NOT use `gh pr checks --watch` (recorded trap).
5. **After merge:** update the ADR-52 implementation log + this handoff + the
   `ladruno-integrator-study` memory.

## Gotchas discovered this session (don't relearn these)

- **Bulk viscosity must be tested under EXPLICIT central difference**, not implicit
  Newmark. The viscous force is velocity-dependent and added to the resisting force but
  NOT the tangent → an implicit solve treats it explicitly and blows up at large `b1`.
  The `zone_a` gate caught this (#403 first run failed nan/1e268). Use
  `CentralDifferenceLadruno` + `system Diagonal` (lumped element mass) for any
  bulk-viscosity runtime test.
- **ADR number 50 collided** with the concurrent AEM line (50 scoping, 51 element-
  removal); this work is **52**. Check `classTags.h` and the `Ladruno_implementation/`
  numbering before claiming a number.
- **Integrator classTag band is occupied 33000–33012** (`classTags.h:1143-1155`).
  New integrators start at **33013**; assign from the top, re-check immediately before
  the PR (W1-E2 may free 33002/33009/33011 — don't pre-reserve those).
- **2D fork elements are `-geom linear` only** (no `theGeom`) → no geom guard needed
  there; LadrunoBrick DOES need the `METHOD_LINEAR` guard (corot/finite fall through).
- **`getInitialTangent()` not `getTangent()`** for artificial wave speed `c_d` (elastic,
  not degraded-EP).
- Adding a new integrator needs the ledgered registration set: `tcl/commands.cpp`,
  `interpreter/OpenSeesCommands.{cpp,h}`, both brokers (`FEM_ObjectBrokerAllClasses.cpp`,
  `runtime/.../TclPackageClassBroker.cpp`), `classTags.h`, the integrator CMakeLists +
  Makefile. (NOT DirectIntegrationAnalysis/TclWrapper/PythonWrapper — generic wiring.)
  Adding a new *command* (not integrator) is lighter: `OpenSeesCommands.cpp`(or
  `OpenSeesMiscCommands.cpp`/`OpenSeesOutputCommands.cpp`) + `OpenSeesCommands.h` decl +
  `Py_ops_`/`Tcl_ops_` wrappers + `addCommand` in `PythonWrapper.cpp`/`TclWrapper.cpp`.
  No classTag, no broker. (W1-I1b used this.)
- **`setNodeDisp/Vel/Accel` do NOT accumulate across dofs** — each call re-reads the
  COMMITTED vector (`Node::getDisp()`), overrides one dof, and `setTrialDisp`s, so for a
  multi-dof node only the last-set dof differs from committed. To inject a full multi-dof
  trial state use the W1-I1b `ladrunoSetNodeTrial` full-vector setter.
- **`ops_Dt` is a GLOBAL read by rate-dependent materials inside `Element::update()`**
  (Maxwell/ViscousDamper/creep/TDConcrete...). `Domain::applyLoad(t)` sets it to
  `t − committedTime`; after a commit `committedTime = t_{n+1}`, so applying loads at a
  PAST time (e.g. the step midpoint) makes `ops_Dt` NEGATIVE and corrupts those materials
  (`exp(-ops_Dt/tR)` → growing). W1-I1b's review caught this; fix is to force a positive
  dt and save/restore `ops_Dt` around any off-step residual assembly. Elastic-only tests
  mask it — add a rate-dependent (Maxwell) case when touching residual machinery.
- **(W3-I2) `ci/check_manifest.py` (G9 gate) requires a `manifest.yaml` row per new classTag.**
  A new integrator fails the fast "classTag + manifest gates" job until you add a row in
  `Ladruno_implementation/testbed/manifest.yaml` (`tag_symbol` + `pytest` + `status: active`).
  Run `python3.12 ci/check_manifest.py` + `ci/check_classtags.py` locally before the PR.
- **(W3-I2) base integrator ctors may hardcode their classTag.** `HHT`/`GeneralizedAlpha`
  ctors call `TransientIntegrator(INTEGRATOR_TAGS_HHT)` directly (unlike `Newmark`, which has
  a classTag param). `getClassTag()` is non-virtual + `classTag` is private with no setter →
  the only way a subclass reports its own tag is a base ctor that takes one. Add a protected
  **inline** classTag ctor to the base header (keeps the base `.cpp` byte-identical).
- **(W3-I2) `GeneralizedAlpha` has an INCONSISTENT tangent for αM≠1** (tangent `αM·c3·M` vs a
  primal that integrates inertia at the full-step `Udotdot`, effective M-coef `c3`; `update()`
  sets accel=`Udotdot`). A tangent-reuse DDM gives a biased gradient (FD oracle ~2e-3 at
  αM=0.9). Build the DDM PRIMAL-consistent (M at `Udotdot`, no αM) AND re-form the
  sensitivity-solve tangent with `c3·M` (a `sensTangentFlag` branch + `formTangent()` in
  `computeSensitivities`). HHT is unaffected (no αM). Full note in [[LEDGER_quirks]].
- **(W3-I2) `Element` base `getTangentStiffSensitivity`/`getMassSensitivity` return ZERO + warn**
  ("betaK·K DDM not implemented") and Truss doesn't override them ⇒ you CANNOT FD-test a
  `∂C/∂h` or `∂M/∂h` term on a Truss (a stiffness/mass-proportional Rayleigh case false-fails
  for the element limit). The M-chain term (`M·dUdotdot/dh`) IS exercised; the matrix-sens
  terms are derivation-validated + pinned by the αM=αF=1→Newmark reduction.
- **(this session) PR numbers are NOT sequential** — #414 was a concurrent CONTACT PR, so the
  genalpha PR landed as #415. Always read the gh-assigned number; don't pre-write it in docs
  (backfill the previous PR's number in the next PR, the established pattern).

## Remaining waves — scope + notes

### W1-E2 — collapse `ExplicitBathe*` 6 → 1 + flags (`-lnvd -sms -consistent`) — SHIPPED #419
- **DONE.** One `ExplicitBathe` (33000) + orthogonal `-lnvd <alpha>` / `-sms <dtTarget>` /
  `-consistent`. Net-negative line count; the registration deletions (both brokers, integrator
  CMakeLists + Makefile, the 5 `classTags.h` deprecation annotations) are ledgered in
  [[LEDGER_vanilla_files]]; the 5 deleted fork `.cpp/.h` + the unified class in
  [[LEDGER_implementations]].
- **Architecture correction:** the six were a **2-base × 3-SMS-mode lattice**
  (`ExplicitBathe` / `ExplicitBatheLNVD` bases, each × {none, sms-lumped, sms-consistent}),
  **NOT 6 flat siblings** as this handoff originally said. The two bases' `newStep`/`update`/
  `commit` stepping was byte-identical; only LNVD's `formUnbalance()`+`addLocalDamping()` and
  the SMS `domainChanged()` injection differed → flag-gated on the one class (`-consistent`
  reuses the `refineAccel` hook the base already carried).
- **Serialization (how the broker keeps working):** NOT "`new ExplicitBathe()` + flags in the
  payload". Instead the flag combo **derives** the classTag (`tagForFlags`); both brokers route
  the six tags via a case-fallthrough → `ExplicitBathe::makeForBroker(classTag)` (decodes flags
  from the tag); `sendSelf`/`recvSelf` carry a fixed-size param superset (flags implied by the
  tag). **No new tag, none freed.** The ~12 strcmp dispatch branches + the 6 OPS_ parsers are
  retained (each parser preserves its exact historical grammar, forcing fixed flags). Pattern
  written up in [[LEDGER_quirks]] ("Collapsing N sibling classes…").
- **Tests:** `tests/test_adr52_w1e2_explicitbathe_collapse.py` (zone_a) — flag-form ==
  legacy-alias byte-identity (5 combos), the LNVD×consistent cross-layer validation, and the
  reduce-to-base gating. The pre-existing per-alias batteries still run under the deprecated
  names (free regression). 6-lens adversarial Workflow review run before merge.

### W1-I1b — half-increment-residual error gate (SHIPPED #407)
- Done. `robust_transient(error_gate=True, haftol=…)` in `Ladruno_scripts/robust_drive.py`
  + `ladrunoTrialResidualNorm <loadTime>` / `ladrunoSetNodeTrial` commands. Open question
  resolved: NOT Python-only (`setNodeDisp` triggers no `Element::update()`).
- **Feed-forward** (OpenSees commits on success → sizes the next Δt, no mid-step rejection);
  verdict `accuracy_gated` iff every committed step met `haftol`, else `integrated`.
- **Fidelity:** exact for rate-/path-independent (elastic); approximate for inelastic/
  rate-dependent (post-commit reference is t_{n+1}; representative +dt imposed). See the
  `ops_Dt` gotcha above.
- *Possible future polish (non-blocking):* true mid-step rejection would need a C++
  `analyze`-without-commit variant; an Abaqus-style HAFTOL auto-scaling from a running
  reference force; an HHT/generalized-α α-weighted midpoint (current uses avg-accel).

### W3-I2 — sensitivity-carrying `LadrunoHHT` / `LadrunoGeneralizedAlpha` (SHIPPED #413 + #415)
- **Done.** `LadrunoHHT` (classTag **33013**, #413) + `LadrunoGeneralizedAlpha` (classTag
  **33014**, #415). DDM on the numerically-damped integrators (was Newmark-only).
- **Vanilla edit was header-only but bigger than the ADR predicted:** besides the
  `private:`→`protected:` promotion, the base ctors hardcode their classTag (no classTag
  param, unlike `Newmark`), so each header also got ONE protected **inline** classTag ctor
  → the `.cpp`s stay byte-identical. (ADR-50 review's "pure promotion" under-scoped this.)
- **Seam (8 methods, not 5):** `formEleResidual`/`formNodUnbalance` (branch on
  `sensitivityFlag`) + `formSensitivityRHS`/`formIndependentSensitivityRHS`/`saveSensitivity`/
  `commitSensitivity`/`computeSensitivities`/`revertToStart`. Because `U/Udot/Udotdot` follow
  the Newmark recurrence, all but `formEleResidual`/`formNodUnbalance` are **copies of
  Newmark**; only those two carry the α-weighting + the extra `−K·(1−αF)·dUₙ` term.
- **classTags were 33013/33014, NOT the ADR's 33014/33015** — W3-I3's reserved 33013 was
  NO-GO so it was the lowest free tag (assign from the top, re-check `classTags.h`).
- **GeneralizedAlpha tangent is INCONSISTENT** (αM·c3·M tangent vs Udotdot primal inertia,
  M-coef c3) → the DDM had to be PRIMAL-consistent (M at Udotdot, no αM) + re-form the
  sensitivity-solve tangent with c3·M (`sensTangentFlag`). HHT is fine (no αM). See gotchas
  + [[LEDGER_quirks]]. FD oracle caught the first (tangent-consistent) cut ~2e-3 at αM=0.9.

### W3-I3 — tunable implicit Bathe (GATED → NO-GO, 2026-06-24)
- `TRBDF2` already IS the Bathe (2007) composite (`TRBDF2.cpp:29-30`, constants hard-
  coded at `:114-144`); the would-be gap is the user-tunable γ/ρ∞ variant (Bathe-Noh 2012).
- **Gate result: DON'T BUILD.** Benchmark `Ladruno_scripts/w3i3_bathe_gate_benchmark.py`
  (impact on an ENT stop swept to the rigid limit + a 1-D wave-propagation bar) showed:
  no robustness gap (every scheme converged to k=1e12), `HHT(α)` already spans ρ∞ and
  tames under-resolved-impact energy blowup, and the composite's only real edge is MAX
  dissipation — which the *fixed* `TRBDF2` already gives (it suppressed wave ringing best
  at Courant 4). A *tunable* composite (ρ∞>0) would be worse there, so its quadrant isn't
  needed. Guidance instead: `TRBDF2` = set-and-forget max dissipation; `HHT`/`Gα` =
  tunable. Revisit only if a stiff multi-DOF contact problem makes `HHT`/`Gα` diverge
  where a composite survives (caveat: gate was SDOF + 1-D bar, not large 3-D contact).
- (If ever resurrected: `LadrunoBatheImplicit` subclass of `TRBDF2`, classTag 33013.)

## W2-E1 follow-ups (deferred, non-blocking)
- **S3** `bvDissipated`/ALLVD recorder channel — mirror `hgDissipated` (accumulate
  `∫ s·ε̇_vol dv` in `commitState`, expose via `getResponse`, thread through
  `sendSelf/recvSelf`). Energy balance already closes; this just makes ALLVD reportable.
- **S4** one-time warning in `setDomain` when coeffs>0 and material `rho==0` (silent
  no-op today).
- Extend bulk viscosity to the **uri/ssp/eas** single-point Brick/Quad paths (currently
  stripped at parse for those formulations).

## Cleanups that are docs-only (would touch vanilla — DO NOT edit)
Per ADR-52 "Explicitly dropped": `Newmark1`→`Newmark`, `ArcLength1`→`ArcLength`, the
`ExplicitDifference` bug, `EQPath`/`HSConstraint` validation — document in
[[ladruno_integrators_guide]] / [[LEDGER_quirks]]; the `ExplicitDifference` fix could be
upstreamed as a *real-OpenSees* PR, separate from the fork.

## Key files
- Roadmap + log: `Ladruno_implementation/52_ladruno_integrator_strengthening_adr.md`
- Study + scorecard: `49_ladruno_integrator_study_workflow_adr.md`, `49a_integrator_scorecard_2026-06-23.md`
- Robust driver (W1-I1a/b): `Ladruno_scripts/robust_drive.py` (`robust_transient`,
  `_half_increment_residual`); W1-I1b commands in `SRC/interpreter/OpenSeesCommands.cpp`
  (`OPS_LadrunoTrialResidualNorm`) + `OpenSeesMiscCommands.cpp` (`OPS_LadrunoSetNodeTrial`)
- SMS (W1-E3): `SRC/analysis/integrator/{CentralDifferenceSMS,ExplicitBatheSMS}{,Consistent}.{h,cpp}`
- Bulk viscosity (W2-E1): `SRC/element/ladrunoBrick/LadrunoBrick*.{h,cpp}`,
  `SRC/element/ladrunoPlane/{LadrunoQuad,LadrunoCST}*.{h,cpp}` (search `bulkVisc`/`bvActive`)
- Runtime tests: `tests/test_adr52_w2e1_bulk_viscosity*.py`, `tests/test_adr52_w1i1b_halfres_gate.py`
- Memory: `ladruno-integrator-study` (project memory)
