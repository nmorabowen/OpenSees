---
title: LadrunoSANISAND — session handoff (after PR-1, PR-2 and PR-3)
project: Ladruno
status: active
owner: nmora
tags:
  - handoff
  - adr
  - material
  - sanisand
  - manzari-dafalias
---

# LadrunoSANISAND handoff

State after PR-1 ([#767](https://github.com/nmorabowen/OpenSees/pull/767), merged), PR-2
([#768](https://github.com/nmorabowen/OpenSees/pull/768), merged) and PR-3. Spec is
[[86_ladruno_sanisand_adr]]; read it first, **including its correction boxes** — several of its
original numbers were refuted by measurement and corrected in place.

---

## 1. Read this before you touch anything

Four traps cost real time in the first two PRs. Each is cheap to avoid and expensive to rediscover.

**A plain `pytest` tests the WRONG BINARY.** `opensees_venv` carries `ladruno_opensees.pth`, which at
interpreter startup puts `C:\Program Files\Ladruno\OpenSees\bin` on `sys.path[0]` and eagerly
imports `opensees`. Measured: `test_manzari_safety_pack.py` is *2 failed, 1 passed* against the
installed engine and *3 passed* against the worktree one. Use a rebinding plugin and **confirm the
bound path on stderr every run**. Never modify anything under `C:\Program Files\Ladruno\`.

**Your binary can be stale, and green tests will not tell you.** A rebuild that finishes in under a
second did nothing. A header change to `ManzariDafalias.h` triggers a **~1500-TU, >10-minute**
rebuild. At the end of PR-2 the working binary was ahead of `SRC/` and nobody noticed until the
final rebuild recompiled 1492 TUs. **Rebuild from committed HEAD before recording any number.**

**CI does not run what you think it runs.** No workflow sets `--runslow` or `LADRUNO_RUN_SLOW`, so
anything marked `@pytest.mark.slow` **never runs in CI**. Before PR-1's fix, `_opensees_exe()` looked
only in `dist/bin`, which CI never creates — so the classic-Tcl gate skipped on every push and had
*never once executed*. Marking a test slow is equivalent to deleting it from CI.

**A skip and a pass look identical in a summary line.** `WORKFLOW_GOTCHAS` §7 records the G1 Tcl
suite reporting `OK (0 checks passed)` for months. When you arm a gate, expect it to go red — that
is evidence it exists.

---

## 2. The architecture, in one page

`LadrunoSANISAND : public ManzariDafalias`. Not a copy, not a wrapper — **inheritance**. When a deck
uses it, `integrate()` / `ModifiedEuler()` are `ManzariDafalias`'s functions, the same compiled code
vanilla runs. Our ~600 lines contain **no physics**.

What it changes: two `protected` members (`ManzariDafalias.h:162-163`) that ~30 sites in the physics
read at run time. Vanilla's `initialize()` hardcodes them; we let it run and then overwrite —
***win the last write***. That is also why G1 is bit-identical: given vanilla's constants, it is
vanilla's code on vanilla's data.

The six overrides exist only to stop the base taking those values back (re-`initialize()`,
`revertToStart` mid-analysis, a wire format with no slot for `p_r`, `getCopy` hardcoding vanilla
wrapper types at every Gauss point, and `Print` so the record can state what it ran).

**The ceiling.** Inheritance reaches `protected` data and `virtual` methods. `ModifiedEuler`,
`integrate` and `elastic_integrator` are **non-virtual** and dispatched through a pointer-to-member
explicitly typed to `ManzariDafalias`, so a subclass cannot substitute them — and you **cannot** fix
that by adding `virtual`, because `ManzariDafaliasRO` (`ManzariDafaliasRO.h:93-97`) *shadows*
`initialize()` and two `GetElasticModuli` overloads with identical signatures. Promoting them would
silently convert those shadows into overrides and change every existing Ramberg-Osgood deck.

### The three lanes (D9)

| lane | for | vanilla decks |
|---|---|---|
| **1. subclass** | anything reachable via protected data | untouched |
| **2. vanilla bugfix** | genuine **errors** | move — disclosed |
| **3. flag seam in vanilla** | **opinions**, and calibration-moving errors | bit-identical |

**The test is: error, or opinion?** A flag seam is a protected member on the base, default = vanilla,
read at the defect site, set only from a derived constructor. Two exist: `mHonorTolRInME` and
`mUseCurrentVoidRatioInG`. **Neither is wired to the subclass or the parser yet.**

**The fork tripwire.** The subclass stops being viable the day we change the *formulation* rather
than fix an error — the `alpha` re-set, `Elastic2Plastic`'s `M_c` repair, D5a's sigmoid shape, §7.3's
`m_e_init`. All are non-virtual and entangled, and changing them in vanilla imposes a modelling
opinion on every existing deck. **D5a is queued for PR-3, so this decision is one PR away.** Open it
consciously, not by drift.

---

## 3. What shipped

**PR-1** — the class, both wrappers, all six overrides, five registration sites, ND tags
33019/33020/33021, the battery, ledgers, manifest, banner, header stamps. Vanilla untouched.

**PR-2** — five vanilla commits (`mTolR` seam, clamp made observable, `D_factor` non-dimensionalised
at four sites, void-ratio interpolant fixed at three, clamp `p` self-consistency), one declined as
ambiguous, one flag seam, and a confine-first test deck.

### Constants (re-measured on an engine rebuilt from PR-2 HEAD)

| | ramp deck | confined deck |
|---|---|---|
| G1 equivalence | **0.0 bit-identical** | 0.0 bit-identical |
| non-vacuity (1.01 vs 0) | 5.579874e-02 | **9.325164e-02** |
| isolate-the-constant (1.01 vs 0.5) | 2.692933e-02 | **4.241113e-02** |
| PlaneStrain sensitivity | 3.177732e-02 | **1.091772e-01** |
| PlaneStrain mirror margin | 0.1866 | **19.157** |
| `\|eps_pl\|` at `p_r = 0` | **0.0 exactly** | 5.675e-03 |
| Outside Bounding / leg | **8** | **0** |

**The interpolant fix (PR-2 commit 4) moves `IntScheme 45` decks by ~0.6 %.** RK4 does not move (its loop reads a separate `CurVoidRatio`); **RK45 does** (11 in-loop reads of `NextVoidRatio`, 0 of `CurVoidRatio`). Measured by A/B rebuild on a NON-isochoric triaxial: reldiff **6.33e-3**. PR-2 initially reported zero because the confined test deck is exactly isochoric and the two anchors differ by exactly trace(d_eps) — a measurement that was true on its deck and falsely generalised. Caught by adversarial review before merge. **Lesson: an inertness measurement is only as general as the deck's strain path.**

Low-`p` gate: vanilla **+18.1258 %**, `p_r = 0` **+0.0754 %**, invariant `err·p_end ~= p_r` on both.
`-Presidual` plastic threshold: between **0.0230 and 0.0240 kPa**.

### ADR-86b (WP-86b, 2026-09-05) — the integrator follow-up ADR-90 GATE U asked for

PR-1/2/3 made the material *honest*. **86b makes it usable on a boundary-value problem**, which is
where every one of these defects had been invisible: the fork's whole SANISAND battery is
zero-free-DOF material-point decks, so a bad tangent costs nothing and a non-terminating integrator
never gets the chance to seize. ADR-90's GATE U put the material into a strip-footing BVP and all
three showed up at once (`_adr90_tau0_qu_band.md`).

| # | what | kind |
|---|---|---|
| **T1** | `-maxSubsteps N` — a substep-COUNT cap on `ModifiedEuler`, default **0 = uncapped**, wired through a new base seam `mMaxSubstepsInME` exactly as `-honorTolR` is wired through `mHonorTolRInME` | lane 3 (flag seam in vanilla) + subclass |
| **T2** | `LadrunoSANISAND`'s parser default for **`TanType` 0 → 2**. Vanilla's parser keeps 0 | subclass only |
| **T3** | convergence-test guidance (`NormUnbalance`, not `NormDispIncr`) | docs — emitter guide §6 |
| **T4** | the `-Presidual` decision a deep push owes | docs — emitter guide §7 and §5 below |

**T1's design point, and it is the whole point.** Vanilla, at `dT == dT_min`, **force-accepts** the
substep and advances `T` (`ManzariDafalias.cpp:1649-1663`). So the integrator never *fails* — it
just takes an unbounded amount of time, and GATE U measured what that looks like from outside: 0 of
80 pinned subdivisions used on every leg, terminal steps 6400–25000x above the `DS_MIN` floor, and
one `analyze(1)` costing **2056 s** (59 % of that leg's entire budget in a single step). The cap
therefore does **not** force-accept: past `N` substeps it flags, prints one throttled line naming
tag/`T`/`dT`, and returns, `setTrialStrain` returns `-1`, and the committed state is untouched. The
precedent is ADR-84's `strict_convergence`; the alternative (force-accept + warn) reproduces exactly
the silence being fixed.

**T1 needed an element change too, and the obvious version of it was wrong.** `LadrunoBrick::update()`
**discarded** `setTrialStrain`'s return code on four of its five paths (std/b-bar, SSP, and both URI
branches) — only `updateHypo`/`updateFinite`/EAS propagated — so a material failure could not reach
the analysis at all. The first cut made all four honour any `< 0`, and that **broke two long-green
gates**: `test_ladrunoBrick_asdconcrete_bend.py`'s two mesh-objectivity tests died at
`Domain::update - domain failed in update`, load factor 605.

`ASDConcrete3DMaterial::setTrialStrain` does return negative, but it means *"an inner iteration
missed, here is my best state"* — the case ADR-33/34 already settled the other way (return 0 with
the last iterate and a loud warning; a failure code there makes softening analyses fragile).
OpenSees has no convention separating that from *"I did not integrate this increment"*, so ADR-86b
introduced one: **`LADRUNO_MATERIAL_REFUSED`** (`SRC/material/LadrunoMaterialStatus.h`, `-33086`).
`LadrunoBrick` propagates **only that exact value**, every other non-zero code keeps the treatment
it had, and adopting it cost zero behaviour change anywhere else. `updateHypo`/`updateFinite` keep
their pre-existing `< 0` tests, untouched.

`stdBrick` still swallows return codes (upstream) — a capped material inside a `stdBrick` still
reports success, and that is why the ADR-86b gates use `LadrunoBrick`.

**T2 is a two-parser divergence, on purpose.** Vanilla `OPS_ManzariDafaliasMaterial` also defaults
`oData[1] = 0` (`ManzariDafalias.cpp:93`) — verified at source, **not** changed, because every
existing vanilla deck and golden file depends on it. The Ladruno parser now defaults 2 and **echoes
the tangent it will run** at construction, so the divergence cannot be silent. `mCep_Consistent` is
unsymmetric under non-associated flow: pair it with an unsymmetric solver.

**Wire format widened** `Vector(4)` → `Vector(5)` (`data(4) = mMaxSubsteps`). Fork-internal: both
ends of any supported channel are the same build. It does break reading a datastore written by a
pre-86b engine — loudly, as a size mismatch, not silently.

**New diagnostic response:** `eleResponse <ele> material <gp> substeps` → `[substeps_taken,
cap_hit]` for the last update. Defined on `LadrunoSANISAND`, so both wrappers inherit it.

---

## 4. Two mechanisms you must know

**The radial freeze.** Every *elastic* step re-sets `NextAlpha = GetDevPart(NextStress)/p` with `p`
including `m_Presidual` (`ManzariDafalias.cpp:974-977`). On a perfectly proportional path with
`p_r = 0`, `s/p` is constant along the ray, so `alpha` is frozen at exactly the current stress ratio
and `f = ||s − p·alpha|| − √(2/3)·m·p` collapses to `−√(2/3)·m·p < 0` **forever**. The material never
yields, on any integration scheme. Not a defect of the `p_r = 0` default — `NTUASand02` ships
`p_r = 0` with identical code and is used for staged foundation work where the stress direction
changes constantly. It is a property of **proportional decks**, and a user can hit it with a
single-element calibration deck. Pinned by `test_radial_ramp_with_pr0_never_yields`.

**`M_c` inflation.** On the 3D ramp deck every leg — including vanilla's — starts the plastic stage
outside the bounding surface, and `Elastic2Plastic` silently raises `M_c` from 1.3309 to 1.99878, 8×
per leg. Any η/M^b measured there is against a surface the material stopped using. Zero on the
confined deck, and now asserted.

---

## 5. Owed work

**DONE in PR-3** (kept here so the next reader can see what moved, and what the doing of it taught)
- ~~Wire `-honorTolR 1` through to `mHonorTolRInME`~~ — wired from `applyLadrunoConstants()`.
  **Learned:** the seam is read at exactly ONE site, inside `ModifiedEuler`, so the flag is inert
  on IntScheme 2/3/4/5/6/7/8/9/45 (measured 0.0 on 3 and 45). Read the dispatch, not the names:
  IntScheme 7 is called `INT_MAXSTR_MFE` and reaches `ForwardEuler`.
- ~~A behavioural `-Pmin` gate~~ — shipped. **The predicted mechanism was wrong.** This section
  said the gate was writable "via `Elastic2Plastic`'s `3*m_Pmin` reset" and via PR-2's clamp
  diagnostic. Neither fires: on the gate's deck `Elastic2Plastic` does not trigger (`p` at the
  stage flip stays above `m_Pmin`) and **zero clamp warnings are emitted**. What fires is a
  silent whole-tensor reset to `m_Pmin*mI1` at a site PR-2 did not instrument.
- ~~`LadrunoSANISAND3D::getCopy(void)` is uncovered~~ — covered, and proven by mutation.
  **But only the 3D twin.** `LadrunoSANISANDPlaneStrain::getCopy(void)` is still uncovered and
  the same route does not reach it: `InitStrainNDMaterial::getCopy(const char*)` delegates to
  the no-argument form only for `"ThreeDimensional"` (`:292-293`); the `"PlaneStrain"` branch
  calls `theMaterial->getCopy(type)`, the string form (`:310`). No deck-reachable caller of the
  no-argument form was found on that lane. Recorded in the manifest row, not claimed closed.

**Still owed. The first item is the one to fix first.**

- **CI verifies none of this ADR's headline claim.** `test_presidual_is_the_low_p_defect` is the
  ONLY test that measures the +18.13 % vs +0.075 % departure from the model's own bounding-surface
  identity -- the silent-wrong-answer defect the whole ADR exists to fix -- and it is
  `@pytest.mark.slow`, so it runs on **no CI push** (no workflow passes `--runslow` or sets
  `LADRUNO_RUN_SLOW`). Measured by adversarial review: of the file's 17 tests, **15 execute in CI**;
  this one and the MP smoke do not. The marker is honestly disclosed and has a real reason (ADR
  risk 6: its `p_r = 0` leg fails at 400 steps and needs 1200), but disclosure is not coverage.
  Options, cheapest first: run the slow tier on a nightly/self-hosted job; or split out a
  cheaper-but-still-diagnostic leg at a less extreme `p0` that can carry the identity assertion in
  the default tier. **Do not simply unmark it** -- that reintroduces the convergence fragility the
  marker exists for.
- Add a **subprocess** gate for `-honorTolR 1` on IntScheme 3. Its inertness is measured once and
  argued from the dispatch, but nothing pins it, and it cannot be pinned in-process:
  `ManzariDafalias`'s scheme-3/5 warning is behind a once-per-process static latch that
  `test_manzari_safety_pack.py::test_scheme3_warns_no_error_control` must observe first.
- `-Presidual`/`-Pmin` **parser error paths** have no test (only `-honorTolR`'s do). Same defect
  family as ADR sec.7.1.
- `revertToStart`'s **`ops_InitialStateAnalysis` branch** is untested; only the plain mid-life
  `reset()` path is exercised, though the guard is the reason the override exists (sec.4.4).
- The **analytical Jacobian** (`JacoType`) branch is never exercised on this class.

**Also owed, and PR-3 added the first two**
- **Instrument the other `m_Pmin` resets.** PR-2's throttled clamp diagnostic covers
  `ModifiedEuler:1378` and `RungeKutta45:1902`. **Three** more sites rebuild the whole stress
  tensor from `m_Pmin` and say nothing: `explicit_integrator:1074/1078`,
  `Stress_Correction:2555/2557/2624` (a **live** `p = m_Pmin + m_Presidual` store — the twin PR-2
  found dead in `ModifiedEuler` is live here), and `BackwardEuler_CPPM:2213/2229`. Same throttle
  shape as PR-2's: a PROCESS budget, not per instance.
  **Note which way the asymmetry runs:** the two INSTRUMENTED sites preserve the deviator
  (`GetDevPart(NextStress) + m_Pmin*mI1`); all three UNINSTRUMENTED ones write a purely isotropic
  tensor and zero `alpha`. The warning covers the gentler rebuild.
  *(A draft of this list also named `Stress_Correction:2608`. It is dead code — inside the
  `if (false)` opening at `:2559` — and was one logical site counted twice. Caught by adversarial
  review; do not re-add it.)*
- **The fifth `D_factor` sigmoid** (`NewtonSol2:3556-3563`) is dimensionally worse than the four
  PR-2 repaired and its steepness is proportional to `m_Pmin`. It is **dead code** today
  (`NewtonIter3` has no caller anywhere in `SRC/`), so it is recorded rather than fixed — but
  PR-2's "all four sites" is incomplete about source, and this is why.
- A **silent-elasticity diagnostic**: override `commitState()`, and if the plastic stage is active
  while `mDGamma` has been identically zero for N consecutive commits, warn (throttled, process
  budget — *not* per instance; `getCopy` makes every Gauss point an instance). Catches the radial
  freeze for any user on any deck. Entirely subclass-side. **Considered for PR-3 and deliberately
  not taken** (owner's scope call), not blocked.
- Extend the echo/`Print`: `p_r = 0` costs more than cohesion — it also removes a **numerical
  regulariser** (the low-`p` leg fails at 400 steps where vanilla survives; needs 1200) and permits
  the radial degeneracy. PR-3 got part of the way: the echo and `Print` now name the honoured
  `TolE` and the true D5a/D5b state, but neither yet warns about the regulariser.

**FORK TRIPWIRE — investigated in PR-3, NOT taken. Read the memo before touching either.**
[[86_ladruno_sanisand_pr3_tripwire_memo]] measures both and takes no decision.
- **D5a** — should the `D_factor` sigmoid exist at all once `p_r = 0`? **This question is no longer
  neutral for us.** `GetStateDependent` hands the sigmoid a `p` that INCLUDES `m_Presidual`, the
  sigmoid's half-suppression point is 1.0500 kPa and vanilla's `p_r` is 1.01 kPa — so vanilla's
  residual pressure was **bounding `D_factor` from below at 0.4278**, and our `p_r = 0` default
  drops that floor to `4.830e-4`, a factor of **886**. Measured on one deck: min `D_factor` 0.7227
  vs 0.0016821. Setting `p_r = 0` does not only remove a cohesion; it un-masks a dilatancy
  suppressor. The memo lists the four experiments that would settle it.
- **§7.3** `m_e_init` — seam is open, unwired. Wiring it moves a calibrated quantity: `G0 = 264.32`
  was fitted against the frozen form. Quantified: `dG/G = 2.489 · d tr(eps)`, i.e. **~2.5 % per 1 %
  volumetric strain**, and **no single `G0` rescale reproduces the old model** because the
  correction is state-dependent. Note the current battery **cannot see this change at all** — the
  confine-first deck is exactly isochoric, so `e` never leaves `e_init` — which means a
  non-isochoric gate is a prerequisite for the seam, not a follow-up to it.

**Separate PR — `SAniSandMS`**
- §7.1: `numData = numArgs - 19` at `:134` and `numData -= 5` at `:143` silently drop `TolF`/`TolR`
  for 4–5 optionals. The ADR's line numbers are **exact**; a proposed "correction" to `:132`/`:141`
  was itself wrong.
- Its `D_factor` trigger is `0.001*P_atm`, **not** `0.05` — over that window the factor runs
  4.8e-4…1.0e-3, a ~1000× kill of dilatancy. A *shape* defect on top of the units one.
- Its `RungeKutta4:1550` carries the same interpolant bug.

**Do NOT re-propose (D8):** upstreaming. Considered and declined at fork level. Prof. Gorini is to be
informed; the upstream-PR call is **his**.

---

## 5b. `-Presidual` on a deep push — the decision, written up and NOT taken (ADR-86b T4)

ADR-90 GATE U's owed-work list item 2 says: *"Settle the low-p floor before any deeper push."*
ADR-86b **does not change the default**, and this section is the write-up of what the owner and the
TIMs consumer have to decide, so that the choice is made rather than inherited.

**What the default is, and what it costs.** `-Presidual 0.0` (a cohesionless sand has no cohesion)
leaves the low-`p` floor at `-Pmin` alone, `1e-3 * P_atm = 0.101` kPa. Both `ModifiedEuler` and
`RungeKutta45` **rebuild the stress** at that floor — deviator preserved, pressure pasted in — and
say so in `opserr`. A clamped integration point is no longer the model's answer.

**Measured (GATE U, strip footing on `LadrunoBrick -formulation bbar`):** the clamp first fires at
`s/B ≈ 0.0153`, on the **coarsest** mesh (`h0 = 1.0`) and the **dense** density (`e_init = 0.60`,
`psi = -0.22`) **only**. The loose (`e_init = 0.6944`) legs and both finer meshes fired zero clamps
over the same and deeper settlements. It is past every number that WP reported — and it will not be
past the numbers of the next campaign, because ADR-86b's whole purpose is to let a leg get further.

**Why this is a decision and not a bugfix, and why it is not one parameter.** `p_residual` is read
at ~30 mean-stress sites, including the one that feeds `GetPSI()`, so it displaces `psi` and with it
`M^b`, `M^d` and the dilatancy. And per the PR-3 tripwire memo it **also bounds the `D_factor`
dilatancy sigmoid from below**: vanilla's 1.01 kPa held `D_factor ≥ 0.4278`, while `p_r = 0` drops
that floor to 4.83e-4 — a factor of **886**. Setting `p_r = 0` does not only remove a cohesion; it
un-masks a dilatancy suppressor (ADR-86 **D5a**, still open). Changing it back is therefore a
modelling move on a **calibrated** soil, not a numerical tidy-up.

**The three answers, and what each buys:**

| option | buys | costs |
|---|---|---|
| **A.** declared non-zero `-Presidual` | a numerical regulariser at low `p` (measured: the `p_r = 0` low-`p` leg fails at 400 steps where vanilla survives) and no clamp | an apparent cohesion `c = p_r·tan(phi)`, a shifted `M^b`, and it re-engages the D5a sigmoid floor — all on Gorini's calibrated set |
| **B.** a small free-surface surcharge | keeps the shallow points confined; changes the BVP, not the material; auditable in one line of the deck | it is a different boundary-value problem from the one TIMs may have specified |
| **C.** accept the clamp, disclose the count | no modelling change at all | every reported result carries "N integration points were set by the clamp" |

**Recommendation (fork-side, non-binding on the consumer): B for a footing campaign, C never
silently.** A surcharge is the cheapest honest fix and does not touch a calibrated constant. A is
defensible but must not be taken without re-opening D5a, because it moves the dilatancy floor by
~886x as a side effect.

**Whichever is chosen, the count must be measured, not hoped for.** The message is `opserr`-only:
capture it (`ops.logFile(path, '-noEcho')` → run → redirect away → count `"CLAMPING"`) and report
the number alongside the result. Emitter guide §7 carries the same decision in consumer language.

---

## 6. Method that worked

**Adversarial review by a different model found something real every single time** — a false
"self-arming" coverage claim, a CI gate that had never executed, three undefended overrides, a false
"8 events" claim, and two wrong premises of the orchestrator's. ADR §9's staffing note is not
ceremony.

**Mutation testing is what separates a gate from a decoration.** Three overrides had no defence at
all and read as covered. And the naive PlaneStrain assertion — "the stress should be compressive" —
*passes the mutant*, because the sign-flipped material still returns compression (−0.73 kPa). Only
the **deviator sign** discriminates. Reading the code would not have told you that.

**Name the deck, the leg and the constants when claiming a diagnostic fires.** Two agents reported
"8 events" about two unrelated diagnostics in one PR. Number collisions are how false claims survive
review.

**G1 is a pair, not a test.** Under the `getCopy(const char*)` mutation — every Gauss point silently
vanilla — G1 *passed*, comparing vanilla with vanilla. The non-vacuity companion caught it. The ADR
asked only for G1; the companion was added on suspicion and is the half that cannot be fooled.

---

## Log

- 2026-08-27 — Written at the end of the PR-2 session.
- 2026-09-05 — **ADR-86b** (WP-86b, branch `wp/86b-sanisand-integrator`): §3's new subsection
  records T1 (`-maxSubsteps`) and T2 (`TanType` default 0 → 2); §5b writes up the `-Presidual`
  decision without taking it. Driven by ADR-90 GATE U (`_adr90_tau0_qu_band.md` §5, §6.3), whose
  §8.1 close-out named this work as the actionable remainder.
