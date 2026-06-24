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
| **W3-I2** | DDM sensitivity integrators: **`LadrunoHHT`** (33013, #413) + **`LadrunoGeneralizedAlpha`** (33014, #414). Subclass HHT / GeneralizedAlpha, reuse the algorithm, add the DDM seam (α-/αF,αM-weighted residual + extra `−K·(1−α)·dUₙ`; rest copies Newmark). Header-only ledgered vanilla edit of `HHT.h`/`GeneralizedAlpha.h` (private→protected + protected inline classTag ctor; `.cpp` byte-identical). Unblocks reliability/fragility/FORM on the damped integrators. | #413, #414 |

Docs/log PRs: #391 (ADRs 49/49a/52), #398, #401, #404.

**Remaining waves:** W1-E2 only (details below). **W3-I2 done** (#413 + #414). **W3-I3 closed — NO-GO** (benchmark
gate, see below).

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

## Remaining waves — scope + notes

### W1-E2 — collapse `ExplicitBathe*` 6 → 1 + flags (`-lnvd -sms -consistent`)
- **Most invasive.** Net-negative line count but NON-zero vanilla (registration
  deletions across the standard files — rebase-fragile, ledger them).
- The 6 are **sibling** classes (all `: TransientIntegrator`), not a chain. `-lnvd`
  routes through a unified `formUnbalance()` override; `-sms`/`-consistent` select the
  mass path in `domainChanged`.
- **Keep the broker cases** (route retired tags → `new ExplicitBathe()` and extend
  `sendSelf/recvSelf` to carry the `{lnvd,sms,consistent}` flag set decoded from the
  incoming tag) or saved-DB/parallel `recvSelf` breaks. No alias facility exists —
  during a deprecation window the ~12 strcmp dispatch branches are retained.
- **Do the adversarial review.** Add `zone_a` tests proving each flag combo == its old
  dedicated class (byte-identity).

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

### W3-I2 — sensitivity-carrying `LadrunoHHT` / `LadrunoGeneralizedAlpha`
- **Requires 2 ledgered vanilla base-header edits**: promote the needed members in
  `HHT.h` / `GeneralizedAlpha.h` from `private:` → `protected:` (no `.cpp`/algorithm
  change), marked `// Ladruno: protected for sensitivity subclass`, recorded in
  [[LEDGER_vanilla_files]]. Pure subclassing is impossible (confirmed by the ADR-50
  review) — Newmark only works because it's already `protected`.
- The seam is **FIVE virtual overrides**, not three: `formSensitivityRHS`,
  `saveSensitivity`, `commitSensitivity`, **plus** `formEleResidual`/`formNodUnbalance`
  branching on `sensitivityFlag` (re-derive the α-weighted M/C terms; the Newmark
  pattern at `Newmark.cpp:577-747`).
- classTags 33014 (`LadrunoHHT`), 33015 (`LadrunoGeneralizedAlpha`) — re-confirm.

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
