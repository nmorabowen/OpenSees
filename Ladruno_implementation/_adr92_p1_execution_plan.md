---
title: "ADR 92 / P1 — execution plan: `-implex` on LadrunoSANISAND3D"
project: Ladruno
type: execution-plan
status: "OPEN — W1-W9 C++ WRITTEN and syntax-clean, NOT BUILT; every gate in §2 pending the owner's build; D0 discharged at CP1"
priority: high
owner: nmora
orchestrator: "Opus 5 session on `wp/92c-implex-p1`"
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[_adr92_cp1_surcharge_results]]"
tags: [adr, sanisand, implex, p1, plan]
updated: 2026-09-05
---

# ADR 92 / P1 — execution plan

**Why this is open.** D0 gated the C++ on #792's T8. T8 reported (the cap works, no
plateau), and CP1 supplied what T8 lacked: a measurement of *where the wall clock goes*.
Only ~30–35 steps per leg converge on plain Newton; 61–83 % need a fallback rung; **89–93 %
of all Newton iterations are spent in rungs that fail and are discarded**
(`_adr92_cp1_surcharge_results` §6). A linear step has no ladder. That is IMPL-EX's
structural claim and the one thing P0 could not test.

**What P1 is NOT.** It is not a licence to assume the cost win. §5's BVP gate states the
prediction so it can fail, and the table gets reported either way.

---

## 1. Work items

| | item | notes |
|---|---|---|
| **W1** | Extrapolated stress, **incremental**: `σ~ = σ_n + Ce(p_n):(Δε − f·Δε_p)`, tangent `Ce(p_n)` | `G`, `K` **frozen at the committed `p`** — the base already reads them there (`:1008`, `:2223`). The total-strain form is 78 % wrong at `p0 = 5` at a *zero* increment (P0 §8); it is the negative control, not an alternative. |
| **W2** | **Floor clamp on `σ~`** | `σ~ = dev(σ~) + p_min·I1` when `tr(σ~)/3 < p_min`. Without it P0 measured `−1.37 kPa` on a `+0.0101` state. See the open question in §4 — clamp vs bounding `f`. |
| **W3** | History: one new committed vector `Δε_p` = `eps_p(n) − eps_p(n−1)`, from `mEpsilon − mEpsilonE` | Integrator-agnostic (D1). **No vanilla hook** — that is the zero-footprint claim to protect. |
| **W4** | Companion return at `commitState`, **scheme 1 + `-maxSubsteps` required** | D3 as reversed: scheme 2's low-`p` Newton is a literal `errFlag = 0` (`:2264`) and 58–74 % of corner calls are `ModifiedEuler` anyway. Refuse schemes 3/5/7/8/9. |
| **W5** | Stage-0 inertness | `mElastFlag == 0` ⇒ `-implex` is inert; history initialises at the stage flip, not before. Gravity and the `LoadControl 0.0` hold must be **bit-identical**. |
| **W6** | `-implexDt {pseudo\|strain\|user}` + guards | `dt = 0` holds; refuse solved-load-factor integrators. See the June `LEDGER_quirks` entry — and note ADR-92 D2 *departs* from its clamp-and-degrade fix on purpose (the ADR says why). |
| **W7** | **`-implexControl`** through `LADRUNO_MATERIAL_REFUSED` (`-33086`, #792) | Moved up from P2: P0 measured A unusable from `Δε = 5e-4` at `p0 = 5 kPa`, and the corner is where the campaign lives. **Reuse the sentinel, propagate by exact value only** — a blanket `< 0` broke two green ASDConcrete gates once already. |
| **W8** | `implexError` / `avgImplexError` responses | `ASDConcrete3DMaterial.cpp:2073-2077` is the template. |
| **W9** | The six ADR-86 overrides | `getCopy` ×2, `sendSelf`/`recvSelf`, `revertToLastCommit`, `Print`. The wire grows; the FSPM `getCopy` lesson applies. |

**Not in P1:** plane strain, the cyclic/reversal test (both P2); vanilla `ManzariDafalias`;
any change to `LadrunoUP`.

## 2. Gates, ordered by what they would refute

1. **BVP gate (decisive).** CP1's deck (`h1.0_e0.6944`, `Q = 10`, cap 1000) with
   `-implex -implexControl`, same ladder decomposition as CP1 §6.
   **Prediction: past-rung-1 61 % → near zero; failed-rung iterations 89 % → single digits.**
   If the ladder still fires at CP1's rate, the cost case is refuted and P1 closes as a
   correctness-only feature. Report the table either way.
2. **Oracle parity.** C++ `implex_A` reproduces `adr92_p0_oracle/sanisand_implex_oracle.py`'s
   `implex_A` to `1e-8` on the recorded binary paths — G0's discipline one level up.
3. **Tangent identity.** Returned tangent == numerical `dσ~/dε` to machine precision (P0
   measured `3.5e-11`).
4. **Floor clamp.** `tr(σ~)/3 ≥ p_min` at every iteration of every step on the P0 G3 path;
   the unclamped oracle (`−1.37 kPa`) is the negative control.
5. **Byte-identity** with `-implex` unset, on every existing SANISAND deck; **stage-0
   inertness** bit-identical.
6. **Symmetry.** Assembled tangent symmetric to round-off; `system Pardiso -matrixType sym`
   reproduces the general solver (ADR-75 P1d) — drained `LadrunoBrick` legs only,
   `LadrunoUP` is unsymmetric regardless (ADR 71).
7. **`dt` guards** and the refusal path: over-tolerance returns the sentinel, the element
   propagates only that value, subdivision engages, the committed state is unchanged.

## 3. Staffing and the build dependency

Per ADR 92 §9. The orchestrator holds the context; agents get a self-contained brief.

| work | agent | model / effort |
|---|---|---|
| W1–W9 C++ | `general-purpose` in this worktree | Opus, high |
| tests (`tests/test_ladruno_sanisand_implex.py`) | `general-purpose` | Sonnet, medium — the author does not write the tests |
| adversarial pass on W2/W5/W7 and the §2.1 reading | `general-purpose` | Fable, one pass |
| `/code-review high` | slash command | — |
| **build** (`Ladruno_scripts\build.bat OpenSees OpenSeesPy`), Esmeralda, merge | **the human** | — |

**The build is the hard dependency in the chain.** No gate above can run until the owner
builds; the C++ can be written and read before that, and gate 2 (oracle parity) is the first
thing to run once a binary exists.

## 4. Open questions the C++ must not silently decide

- **Clamp vs bound.** §1 W2 clamps `σ~`. The `/code-review` of P0 argued the clamp repairs
  `tr(σ~)` only and leaves an equally over-extrapolated deviator, and that bounding `f` so
  the *whole* `σ~` stays admissible is the deeper fix. **P1 implements the clamp and
  measures the deviator error against the oracle**; if it is the same order as the pressure
  error was, the bound replaces the clamp and gate 4 grows a deviatoric clause.
- **Where `-implexControl` computes its implicit.** ASD runs the implicit on *every*
  `setTrialStrain` when control is on (`:1665-1684`) — the cost IMPL-EX exists to remove.
  The alternative is an a-posteriori check at commit that shrinks the *next* step (one-step
  lag, cannot refuse the step it measured). **Decide on measurement, not preference**, and
  record which in the results memo; the BVP gate is the instrument.

### 4.1 What P1's C++ decided, and on what

*(Written by the P1 C++ author, 2026-09-05, `wp/92c-implex-p1`, on a branch that has never
been compiled. Neither decision is closed by this note — both name the measurement that
would overturn them.)*

**Q1 — clamp vs bound: the clamp is IMPLEMENTED, and the instrument for the decision now
lives in the binary rather than only in the oracle.** `ladrunoImplexMeasureError()` splits
the error at every commit into a deviatoric and a volumetric leg over the same denominator,
on the exact identity `‖Δσ‖² = ‖Δdev‖² + 3(Δp)²`, and both are readable from a deck through
the new `implexDetail` response together with the clamp's per-point fire flag and count. So
gate 4 need not be argued from the oracle: **the C++ reports the two legs itself, at every
Gauss point, on the real path.** The rule is unchanged — if the deviatoric leg comes back
the same order the pressure leg was before the clamp, the clamp is the wrong fix and a bound
on `f` replaces it.

**One argument for the bound the plan did not carry, found by writing the code.** The clamp
is a *projection*, so wherever it fires `d(σ~)/d(ε)` is the deviatoric projection of `Ce`,
not `Ce`. The material keeps handing out `Ce` (symmetric, positive definite — the point of
the exercise), so the tangent is knowingly inexact at exactly the Gauss points sitting on
the floor, i.e. the free-surface ring. **Bounding `f` keeps `σ~` admissible with no
projection and keeps the tangent identity exact everywhere.** Consequence for gate 3 either
way: *the tangent-identity test must be run on a path where the clamp is idle*, and a test
that reports machine precision on a clamped path is measuring nothing.

**Q2 — where `-implexControl` computes its implicit: IN-STEP (the ASD site), and not on
preference.** Two structural reasons, neither of which the a-posteriori site can answer at
any cost:

1. **Only the in-step site can refuse the step it measured.** P0 measured `implexError` at
   `1.26` on the `p0 = 5 kPa`, `Δε = 5e-4` row, where the answer is not inaccurate but wrong
   (`q/p` 0.092 against 2.07). A one-step-lag check commits that step and shrinks the next.
2. **Only the in-step site can see a COMPANION FAILURE before the step commits.** If the
   `-maxSubsteps` cap fires inside the commit-time return there is no good move left:
   committing that state is garbage, committing `σ~` instead manufactures history, and
   refusing at commit leaves the analysis believing a step advanced that did not. With the
   control on, the same failure is discovered in `setTrialStrain` and the step is refused
   through `LADRUNO_MATERIAL_REFUSED` like any other. (The control-off path still handles
   it — loudly, by refusing the commit and naming this as the reason to turn the control
   on.) **This reason is new; it is not in the ADR.**

The cost objection is the real one and it is **not** settled here: with the control on the
constitutive cost is **two** full returns per step — the in-step probe plus the commit-time
companion — because P1 deliberately does **not** cache the probe's state. A cache would be
correct only under an exact strain comparison plus a seven-field state save, and a bug there
is a silent wrong answer; the 2x is visible, bounded and reversible. **The falsifier is
gate 1.** IMPL-EX's structural claim is that the global step is linear, so 2 returns/step
sits against CP1's up-to-125 state-determination passes. If the BVP gate shows the ladder
still firing at CP1's rate, the claim is refuted and this is the first thing to move —
either to the cache, or to the a-posteriori site.

### 4.2 Three things the plan got wrong when the code met it

1. **`-implexDt strain` silently destroys the tangent identity unless `f` is frozen.** W6
   lists three `dt` sources as if they were interchangeable. They are not: with a
   strain-based source `f = (‖Δε‖/‖Δε_n‖)·α` is a *function of the trial strain*, so
   `d(σ~)/d(ε) ≠ Ce` by a term nobody would look for — and the run still *looks* like
   IMPL-EX (one global iteration). P1 therefore computes `f` **once per step**, at the first
   trial call after a commit or a `revertToLastCommit`, and freezes it. Under `pseudo` this
   is a no-op (`ops_Dt` is constant within a step), which is exactly why the trap stays
   invisible until someone picks another source. **Gate 3 must be run on a non-`pseudo`
   source or it does not test this.** Recorded in `LEDGER_quirks`.
2. **W6's "refuse solved-load-factor integrators" is not implementable as written.** A
   material has no handle on the integrator — there is no `ops_TheActiveIntegrator` beside
   `ops_Dt` / `ops_TheActiveElement`. P1 enforces D2 **by symptom**: a *negative* pseudo-time
   increment is a load factor that has turned round past a limit point, and that step is
   refused with a sentence. **The gap is real and is not papered over: a `DisplacementControl`
   run BEFORE its limit point has `dλ > 0` and is NOT caught.** Such a deck must pass
   `-implexDt user` (or `strain`). If the campaign ever puts `DisplacementControl` under
   `-implex`, this needs a real seam — and that seam would be the first thing in ADR 92 to
   cost a `LEDGER_vanilla_files` row.
3. **Gate 2's parity reference does not cover the clamp.** The P0 oracle's
   `Implex.extrapolate` has **no** floor clamp — the clamp is the item G3 *discovered* — so
   1e-8 parity is meaningful only where the C++ clamp stays idle. On the G3 corner path the
   two *must* differ, and a parity run there that passes is evidence the clamp is not firing,
   not evidence of parity. Either the test agent adds the clamp to a copy of the oracle, or
   gate 2 declares its paths clamp-free and gate 4 covers the rest.

### 4.3 One scope deviation, deliberate

§1 says plane strain is not in P1. `LadrunoSANISANDPlaneStrain::setTrialStrain` is
nevertheless routed through the same `ladrunoTrialUpdate()` seam, because the operator is
dimension-agnostic (`mEpsilon`, `mSigma` and the frozen `Ce` are 6-vectors in both wrappers;
the plane-strain view just reads three of the six back out) and because the alternative —
leaving that line calling `integrate()` — makes `-implex` **silently inert** on a
plane-strain deck, which is the defect class this whole file family exists to prevent (see
its own `-honorTolR` / `-maxSubsteps` inertness warnings). It is **untested at P1**; P2 owns
its gate.

## Log

- 2026-09-05 (later) — W1-W9 written on `wp/92c-implex-p1`. `g++ -fsyntax-only` clean on all
  three touched translation units; **nothing has been built or run** — every gate in §2 waits
  on the owner's `Ladruno_scripts\build.bat OpenSees OpenSeesPy`, and gate 2 (oracle parity)
  is the first to run once a binary exists. Both §4 questions answered in §4.1, three plan
  errors recorded in §4.2, one scope deviation in §4.3. Vanilla footprint: **zero** (only
  `LadrunoSANISAND{,3D,PlaneStrain}` touched); no new classTag.
- 2026-09-05 — Opened. D0 discharged on CP1's ladder decomposition; `-implexControl` moved
  up from P2; branch `wp/92c-implex-p1` cut from `ladruno`.
