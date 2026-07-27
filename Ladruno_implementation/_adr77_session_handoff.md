---
title: ADR-77 session handoff — implicit transient study SHIPPED; next = post-ship reviews
project: Ladruno
status: handoff — ADR-77 complete and merged (#650); REVIEW WAVE DONE 2026-07-27 (see §1a); two parked lanes still queued
owner: nmora
tags:
  - handoff
  - adr
  - transient
  - review
---

# ADR-77 session handoff (written 2026-07-27, post-merge)

> **State: [#650](https://github.com/nmorabowen/OpenSees/pull/650) is MERGED to
> `ladruno`** (2026-07-27, Zone-A green, landed-proof checked). ADR-77 is
> **complete — no open study items**. The full record:
> [[77_ladruno_implicit_transient_adr]] (decisions + gates) and
> [[77a_c0_t0_results_2026-07-26]] (every measurement, §1–§6i). The PR's comment
> thread carries the stage-by-stage writeups. The owner's stated next step:
> **continue the reviews in a new session.**

## 1. What #650 shipped (the review surface)

| unit | what changed | review angle |
|---|---|---|
| **C0-6 + DDM follow-through** | `GeneralizedAlpha::update()` now passes `Ualphadotdot` (was discarding it — every αM≠1 result before the fix is invalid); `LadrunoGeneralizedAlpha`'s DDM re-derived to match (αM-weighted multiplicator, `dM/dh` at `Ualphadotdot`, sensitivity Jacobian `αM·c3·M`) | the DDM algebra was re-derived by me under CI pressure — an independent re-derivation of `formEleResidual`'s a2..a8/αM terms against the generalized-α sensitivity literature is the highest-value single review |
| **C0 wave (15 vanilla files)** | `deltaT` guards ×4; `HALL_TANGENT` branches ×3 (was silent zero/stiffness-free tangents); dead `Newmark` static removed; P3 scopes for 8 EquiSolnAlgo | the HALL branches transposed Newmark's `c1*cFactor/c1*iFactor` pattern — verify that transposition is right for *Collocation's* c1 semantics specifically (its c1 differs from Newmark's) |
| **LadrunoBrick mass cache** | guard-checked per-instance `Mi`, `-noMassCache`; T2-measured −12.6% w/ Rayleigh | reviewed heavily in-session; low residual risk |
| **Mass-cache extension ×5** | shared `SRC/element/LadrunoMassCache.h` → BezierTet10/Tri6/Quad/LST/SolidShell | **only A/B + battery tested, no independent code review yet**; two bugs were already caught and fixed in-session (see §3) — a 5-lens adversarial review of the helper + the 5 integrations is warranted, fork-standard practice |
| **T0 instrumentation** | `dof.tangent`/`dof.residual` scopes in `TransientIntegrator.cpp`/`IncrementalIntegrator.cpp` | trivial |

Review-wave shape that fits the fork's precedent (ADR-41/57 pattern): one
adversarial multi-lens review per row above, patch-in-place fixes, ledger rows
amended. All five units are on `ladruno` now — review against HEAD.

## 1a. Review-wave OUTCOME (2026-07-27)

- **DDM algebra (priority 1): CONFIRMED by independent re-derivation.**
  Differentiating the Newmark update relations against the post-C0-6 residual
  `R = F(t+αF·dt) − P(Uα) − C·U̇α − M·Üα` reproduces every implemented
  constant (`a2=−c3`, `a3=−c2/γ`, `a4=1−1/(2β)`, `a6=−c2`, `a7=1−γ/β`,
  `a8=dt(1−γ/(2β))`), the αM/αF-weighted multiplicators, the `−K(1−αF)dUn`
  term (`addK_Force` uses the consistent tangent at the trial state = Uα
  pre-commit — correct), the `dM/dh@Üα` / `dC/dh@U̇α` states, and the LHS
  `αF·K+αF·c2·C+αM·c3·M` = base tangent. `saveSensitivity` is Newmark-exact
  as claimed. ONE fix: the `computeSensitivities` comment still described the
  pre-C0-6 design (claimed M-coef `c3` vs the base's `αM·c3` — the code says
  `αM·c3`); rewritten. The tangent re-form itself is KEPT — its live value is
  refreshing a stale/initial factorization left by ModifiedNewton/`-initial`
  primal solves, not a coefficient difference.
- **Collocation HALL branch (priority 3): CONFIRMED, no change.** The
  transposition is structurally byte-for-byte vanilla `WilsonTheta` — the
  exact precedent for a t+θΔt-substep integrator (`c1=1.0`, C/M outside the
  chain); `c1` plays the same effective-tangent stiffness-coefficient role in
  both. BackwardEuler/Newmark1 additions spot-checked against their vanilla
  precedents too — clean.
- **Mass-cache extension (priority 2): 5-lens × 5 integrations + helper — NO
  serial-run-corrupting bug; lifecycle hardening applied.** The in-session
  Quad/LST 4-site side-effect fix is verified complete. Patched: recvSelf
  cache invalidation ×6 (sig-exempt `massType`/`quadz`/`formulation`/`cMass`
  can flip on restore-into-live-object — brick donor had the gap too),
  helper null-node tolerance (`setDomain(0)` removal flows / pre-setDomain
  getMass no longer segfault where uncached code survived), contract-scope +
  dangling-hit-pointer + serial-only `-noMassCache` notes in the header, and
  do-not-cache breadcrumbs at the four bare-`getMass()` sites in the
  cache-less CST/CSTPair. Three new [[LEDGER_quirks]] entries. Known
  accepted non-issues: setDomain-snapshot geometry staleness (controlPts /
  SolidShell `X` / SSP J-tables) is pre-existing element behavior the cache
  faithfully preserves; the Bezier/plane miss path still returns the class
  static (per-instance only on hits) — an ADR-75b nuance, not a bug.

## 2. Parked lanes (decided, not started)

1. **LadrunoUP mass review** — the u-p element was *deliberately excluded* from
   the mass-cache rollout (its `buildSolidMass` inputs include the fluid side;
   a pattern-match cache would be reckless). Needs its own read before any cache.
2. **ADR-75b de-statication brief** — now the *only* route to both threaded
   assembly AND cross-element batching/vectorization (V0 closed `/arch:AVX2`:
   no win, no drift — do not add the flag). Blockers list = ADR-75b §5.4.
3. **ADR-49a scorecard caveat** — ✅ DONE 2026-07-27: the "GeneralizedAlpha =
   gold standard" row now carries the post-C0-6-only note (pre-fix αM≠1
   results invalid, #650).

## 3. Traps banked THIS session (do not rediscover — full text in WORKFLOW_GOTCHAS + memory)

- **A/B on/off equality cannot detect a broken default** — both arms identically
  wrong passes. The Zone-A batteries are the correctness oracle; run them for
  any element-side change. (Caught the brace-less-`if` zero-mass bug.)
- **Bare `this->getMass();` side-effect callers**: some elements call it only to
  refill class-static `K` then read `K` directly. Any cache/refactor that skips
  the formation breaks them silently. Grep `^\s*this->getMass\(\);\s*$` before
  touching any element's mass path. (Quad/LST had 4 sites; fixed.)
- **When a primal scheme changes, grep for DDM/sensitivity subclasses that
  mirror it** — `LadrunoGeneralizedAlpha` had faithfully modelled the C0-6 bug
  and broke the moment the primal was fixed. `LadrunoHHT` mirrors HHT the same
  way; HHT's primal was NOT changed, so it is fine — but the rule generalizes.
- **The Bash heredoc layer eats one backslash level even quoted** — inline
  python patch scripts with `\n` in anchors silently never match. Write scripts
  to a file; build backslashes with `chr(92)`.
- **Regression thresholds must sit above the oracle's reproducibility floor**
  (threaded PARDISO ~1e-12 over a history norm; 1–2 ulp on a single uz;
  wall-clock ±12%). Control-run the same build against itself before reading
  any small drift as a regression.
- CI: a `synchronize` push occasionally fails to trigger the workflow — check
  `gh run list --branch <b>` matches HEAD sha; empty-commit to retrigger.

## 4. Environment facts for the next session

- **This worktree** (`pardisio-profiling-0a03b1`) has `dist/bin/opensees.pyd`
  built from the merged state (mass-cache ext + DDM fix, no `/arch` flag) — the
  worktree branch is now `claude/adr77-session-handoff`; branch fresh from
  `origin/ladruno` for any new work (stranded-commit rule).
- Bench stack: `system Pardiso`, `MKL_NUM_THREADS=4` operating point (T0b),
  1 thread for bit-exact oracles (P1f). `numberer LadrunoParallel*` is
  **parallel-only** (serial silently falls back to RCM — T0 gotcha).
- Harnesses live in `Ladruno_files/testbed/perf/phase77/` — T0/T0b/T1/T2/T3
  sweeps, `c06_genalpha_verify.py`, `g2_mass_cache_verify.py`,
  `g2_mass_cache_ext_verify.py`, `c0_wave_regression.py` (all runnable against
  any future build; baselines committed).
- Recommended default for implicit transient work on this stack:
  **`algorithm KrylovNewton`** (T1); `rayleigh betaKinit` over `betaKcurr`
  where acceptable (1.86x cheaper, T2).
