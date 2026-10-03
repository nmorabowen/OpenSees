---
wp: LEGACY
title: "ManzariDafalias's substepped ModifiedEuler has NO iteration cap — only a dT_min floor — so one analyze(1) on a softening BVP can take tens of minutes"
legacy_seq: 360
---
## `ManzariDafalias`'s substepped `ModifiedEuler` has NO iteration cap — only a `dT_min` floor — so one `analyze(1)` on a softening BVP can take tens of minutes

**Found 2026-09-05, ADR-90 WP-A2. It is what stopped that study from reaching an answer.**

`ManzariDafalias::ModifiedEuler` integrates the stress update over pseudo-time `T` in
`while (T < 1.0)` with an adaptive substep `dT` bounded below by
`dT_min = 1e-6` (`ManzariDafalias.cpp:1380`, `:1543`, `:1663`). There is **no substep COUNT
cap**. At `dT == dT_min` the loop force-accepts the substep and advances `T`
(`:1649-1663`), so it does terminate — but the bound is **10^6 return maps per Gauss point per
state-determination pass**, and a Newton ladder that tries three algorithms will repeat that up
to 125 times in one load step.

Measured on a strip footing on softening `LadrunoSANISAND` (`LadrunoBrick -formulation bbar`,
200-782 elements, three meshes, two densities): as the plastic zone develops, single `analyze(1)`
calls start costing **11 minutes**, then **20-28 minutes**, on every mesh and both densities at
once. The deepest leg reached `s/B = 0.0228` of a `0.25` target in 40 minutes of push.

- **It does not look like a hang from outside.** The process burns 100 % CPU, the engine log stops
  growing (nothing is failing, so nothing is logged), and the analysis eventually returns a
  CONVERGED step. Only the curve file's mtime tells you.
- **It is NOT the stepping controller, and the controller's own diagnostics prove it.** Every leg
  in that campaign used **0 of its 80** pinned subdivisions and ended with its step **800-12800x
  above** the `DS_MIN` floor. A wall-clock stop under those conditions is a `WALL` seizure whose
  cause is the constitutive integrator, not the guard -- report it that way.
- **A wall-clock budget cannot interrupt it.** A budget checked at the top of a stepping loop only
  binds once the current `analyze(1)` RETURNS, so a leg can overshoot its budget by the length of
  one seized step. Size the budget with that in mind, or watch the curve file.
- This is the boundary-value-problem form of what ADR 90's A0 note measures on a 1-D bar: the
  rate-independent problem needs 4000 -> 64000 load steps across four meshes to converge at all,
  while every regularized run finishes at 250 steps on every mesh. The cost explosion IS the
  ill-posedness.

> **ADDRESSED in WP-86b (ADR-86b, PR pending) — opt-in, default unchanged.** `nDMaterial
> LadrunoSANISAND` gains **`-maxSubsteps N`**, wired to a new vanilla flag seam
> `ManzariDafalias::mMaxSubstepsInME` read at the top of `ModifiedEuler`'s `while (T < 1.0)`.
> **Default `0` = UNCAPPED = exactly the behaviour described above**, so nothing changes for a deck
> that does not ask, and vanilla `ManzariDafalias` stays bit-identical.
> Past the cap the integrator does **NOT** force-accept — force-accepting is what hid the cost in
> the first place. It flags, prints one throttled `opserr` line naming tag/`T`/`dT` (PROCESS budget
> of 10), and returns; the committed state is untouched (`integrate()` writes only trial members),
> `setTrialStrain` returns `-1`, and the step FAILS so a subdivision controller finally has
> something to react to. Precedent: ADR-84's `strict_convergence`.
> **Size the cap from a measurement, not a guess:** `eleResponse <ele> material <gp> substeps`
> returns `[substeps_taken, cap_hit]` for the last update at that point.
> **THE PRECONDITION, and it is the sharp edge of this feature.** A cap is only safe under an
> element that **propagates** a material refusal. The capped update
> returns at `T < 1`, so it leaves a **partially-integrated** trial stress/`alpha`/`fabric` and a
> partial `aCep_Consistent`. An element that discards the return code assembles that partial state
> and reports convergence, which is strictly **WORSE** than the un-capped force-accept it replaces
> (that at least always drove `T` to 1). Nothing checks this at run time — a
> material cannot see its element — so the default `0` (which cannot reach the branch) is the only
> thing standing between a user and that state.
> **CORRECTED 2026-09-14 (WP-99 / F7).** This paragraph used to end "today **`LadrunoBrick` only**"
> and to list "`Brick`, `BrickUP` / `QuadUP`, `stdBrick`" as the discarders. Both halves were
> wrong: `QuadUP` (`FourNodeQuadUP.cpp:419`) does `ret += ...->setTrialStrain(...)` and therefore
> **propagates**, and `stdBrick` **is** `Brick` under its Tcl name (`TclBrickCommand.cpp:210`). The
> audited lists (`9c2f964ea`) are in the dedicated entry
> "`Domain::commit()` discards element commit returns" at the end of this ledger; the short form is
> *propagate any nonzero*: `LadrunoBrick20`, `LadrunoQuad`/`CST`/`LST`, `BezierTet10`/`Tri6`,
> `FourNodeQuad`, `FourNodeQuadUP`; *propagate only the sentinel*: `LadrunoBrick` (ADR-33/34);
> *discard*: `Brick` (= `stdBrick`), `BbarBrick`, `BrickUP`, `SSPbrick`, `SSPquad`,
> `LadrunoSolidShell`. And **at commit time no element propagates anything**, which is the hole
> WP-99's latch closes.
> **One more thing that will bite.** The cap bounds one whole `integrate()`, not one `ModifiedEuler`
> call — `MaxEnergyInc`/`MaxStrainInc` (IntScheme 0/4/6/8/9) call `ModifiedEuler` several times
> inside one material update, which is why the counter is reset in `integrate()` and not at the top
> of `ModifiedEuler`. (IntScheme 7 is *called* `INT_MAXSTR_MFE` and does NOT reach `ModifiedEuler`:
> `MaxStrainInc` has no case for it and falls through to `ForwardEuler` — read the switch, not the
> name.)
> **A `-maxSubsteps` cap of 1 is VACUOUS.** `ModifiedEuler`'s error-controlled stepper always tries
> `dT = 1` (the whole increment) FIRST; that first attempt almost always exceeds `TolE` on anything
> but a trivial elastic step, so it is rejected and retried at a smaller `dT` — the loop body runs at
> least twice (the failed `dT=1` attempt, then the first real substep) before a single "substep" has
> been accepted. So `mSubstepsTakenInME` is `>= 2` for essentially every plastic update regardless of
> how cheap it is, and a cap of `1` fires on every single one of them, not just the pathological ones
> a cap is meant to catch. **Workaround:** never hand-pick a small cap — measure the deck's own cost
> first (`eleResponse <ele> material <gp> substeps` after an uncapped run, exactly as
> `tests/test_ladruno_sanisand_integrator.py`'s gates do) and set the cap below THAT, not below some
> assumed-cheap constant like `1` or `2`. Learned 2026-09-05, ADR-86b review-fix pass.
