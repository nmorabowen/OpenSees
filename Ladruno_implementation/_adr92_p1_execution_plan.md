---
title: "ADR 92 / P1 — execution plan: `-implex` on LadrunoSANISAND3D"
project: Ladruno
type: execution-plan
status: "OPEN — C++ not started; D0 discharged at CP1"
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

## Log

- 2026-09-05 — Opened. D0 discharged on CP1's ladder decomposition; `-implexControl` moved
  up from P2; branch `wp/92c-implex-p1` cut from `ladruno`.
