# BLUE TEAM 2 — verification of the red measurement attack on the ADR-92 P1 BVP gate

Target: PR #798 / `wp/92c-implex-p1`, results commit `70f6bf0e8`.
Worktree read: `C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr92-p1-implex`.
Read-only. Every number below was recomputed from the committed artifacts.

Scoreboard: **CONFIRMED 6** (F1, F3, F5, F7*, F8, F10, F11) — F7 partial on one bullet —
**PARTIAL 5** (F2, F4, F6, F9, F12), **REFUTED 0 outright**. Two red bullets are
materially wrong (F2's cost inference, F7's ADR-79 comparison) and one blue finding
is stronger than anything red filed (see §13).

---

## F1 — the gate returns PARTIAL, not PREDICTION MET — **CONFIRMED**

Four invocations of the committed `adr92_bvp_gate.py`, all against
`--baseline adr92_bvp/implex_OFF`:

| `--implex` argument | past rung 1 | failed-rung iters | VERDICT |
|---|---|---|---|
| `adr92_bvp` (whole tree — `load()` does **not** filter on the `implex` flag, so `implex_OFF` is counted as a treatment leg too; `implex_ON` has no `.json`, only `FAILED.txt`, so it is skipped) | 48.4 → **16.1 %** | 81.6 → **55.5 %** | **PARTIAL** |
| both IMPL-EX arms only (red's staging) | 48.4 → **0.0 %** | 81.6 → **42.5 %** | **PARTIAL** |
| `implex_ctl` alone (**the registered arm**) | 48.4 → **0.0 %** | 81.6 → **85.0 %** | **PARTIAL** |
| `implex_noctl` alone (**unregistered arm**) | 48.4 → **0.0 %** | 81.6 → **0.0 %** | **PREDICTION MET** |

The PASS branch is `if ip <= PASS_PAST_RUNG1 and ifs <= PASS_FAILED_SHARE` (`:136`,
both 10.0). "PREDICTION MET" is reachable **only** on arm B alone, i.e. only by
excluding the arm the pre-registration named. It is not legitimate as a gate verdict.
Two qualifiers red did not give: (a) the memo's frontmatter is scoped —
`status: "PREDICTION MET **on the ladder claim**"` — and clause 1 (past rung 1 = 0.0 %)
does hold on every invocation including the registered arm and the mean; (b) red's
"exactly as committed, on exactly the committed payloads" is loose — no committed
directory contains only the two IMPL-EX arms, so their run required staging (I
reproduced it by copying `implex_ctl` + `implex_noctl` into a scratch dir; numbers match
to the digit). The instrument fact stands: the committed script never prints
PREDICTION MET on the registered configuration, and it prints "do not round it toward
either."

## F2 — "not one step left rung 1" — **PARTIAL** (literal sentence false; cost inference refuted)

Warning counts per algorithm, `grep -c "WARNING <algo>::solveCurrentStep"`:

| leg | NewtonRaphson | NewtonLineSearch | AcceleratedNewton | `CTestNormUnbalance … failed to converge` |
|---|---|---|---|---|
| `implex_OFF` | 2 | 10 | 0 | **35** |
| `implex_ctl` | **25** | **25** | **25** | **0** |
| `implex_noctl` | 0 | 0 | 0 | 0 |
| `implex_ON` | 7 | 7 | 7 | 0 |

25/25/25 confirmed; `nfail = 75 = 3 × nsub = 3 × 25`, `nrelax = 0`. `steps = len(rows)`
and `rows.append` sits **after** the failure `continue` (driver `:753-756`), so the
denominator is survivorship — confirmed. Attempts table reproduces from the JSONs:
OFF 64 attempts / 31 past rung 1 (48.4 %) / 47 failed rungs = 0.734 per attempt;
ctl 110 + 25 = **135** / 25 (18.5 %) / 75 = **0.556**; noctl **142** / 0 / 0.
So the memo's unqualified "**Not one step in either IMPL-EX leg left rung 1**" and
"0 of 252 steps across both arms needed a fallback rung" are **false as written** —
25 steps in arm A ran rungs 2 and 3 and failed both.

But red's *inference* does not follow, and the evidence is in the same logs red counted.
Every one of arm A's 75 rung failures is `solveCurrentStep() -the Integrator failed in
update()` — the element propagating `LADRUNO_MATERIAL_REFUSED` — and arm A has **zero**
`CTestNormUnbalance::test() - failed to converge` lines, against 35 in the control.
The rungs did not grind; they aborted on the first element update, at roughly one
iteration each, not at their 25/40/60 caps. "Exhausts all three rungs every time it
does" is true of the *attempt* and false of the *cost*, and "IMPL-EX cuts failed rungs
per attempt by 24 %" charges refusals as ladder work — the exact error the memo's §2
disclaims and the gate's new `refused` split was written to fix. Red's own "WHAT WOULD
REFUTE ME" clause ("evidence that the warnings are emitted without the rung actually
being attempted") is answered halfway: the rung was entered, no Newton iteration failed.
Red's log table is also misleading for the control: those warnings are message-throttled
(`ManzariDafalias::ModifiedEuler` and `LadrunoBrick::update` both cap at exactly 10),
so OFF's 2/10/0 cannot be read as OFF using fewer rungs than ctl — OFF's real count is
31 rung-1 failures + 16 rung-2 failures from the JSON. Red did not draw that inference,
but the table invites it.

**Verdict:** the sentence must be withdrawn or qualified. The "18.5 % of attempts still
fires the ladder" reframing must not replace it.

## F3 — three pre-registration deviations — **CONFIRMED**, one bullet sharper than red filed

`_adr92_p1_execution_plan.md:50-53`, verbatim:

> 1. **BVP gate (decisive).** CP1's deck (`h1.0_e0.6944`, `Q = 10`, cap 1000) with
>    `-implex -implexControl`, same ladder decomposition as CP1 §6.
>    **Prediction: past-rung-1 61 % → near zero; failed-rung iterations 89 % → single digits.**

1. **Arm.** Registered `-implex -implexControl`; every headline sentence comes from
   `implex_noctl` (`implex_control = None`). Confirmed.
2. **Cap.** `max_substeps: 20000` in all three payloads vs registered `cap 1000` — a 20x
   change. **PARTIAL on "undisclosed"**: the *value* is stated in the memo's own table
   header and commit message; what is undisclosed is that it *deviates*. The cap moves
   the answer more than red said: the `adr86b_t8` `h1.0_e0.6944` legs give
   cap1000 → `s/B` 0.05875, `q_u` **1433.4 kPa**, mode BUDGET-free WALL; cap10000 →
   `s/B` 0.04780, `q_u` **1035.3 kPa** — a 38 % swing in `q_u`, not just depth.
   Caveat against red: those T8 legs carry `surcharge_kpa: None`, so they are the same
   *mesh and void ratio*, **not** "this exact deck" (the gate deck has `Q = 10`).
   The plan's own Log (2026-09-06) does say "at an adequate substep cap" and "the cap
   was the variable" — the cap's *sensitivity* is disclosed on the branch; its
   *re-registration at 20000* is not.
3. **Baseline.** Registered 61 % / 89 % (CP1); used 48.4 % / 81.6 % (a new `implex_OFF`
   leg). Sharper than red: the memo's §1 does not merely report against the new
   baseline, it **restates the pre-registration itself with the substituted number** —
   memo §1: *"past rung 1 falls from ~48 % to near zero"* against the plan's and the
   gate docstring's *"falls from 61 %"*. That is an altered quotation of a pre-registered
   prediction, in the sentence that claims it was "written before the run".

Pre-registration timing verified: `adr92_bvp_gate.py` first committed `5e5ec4db7`
2026-09-06 00:44:59; first leg started 01:29:30/32. The thresholds and `past1` were
fixed before the run.
On the registered arm the second clause fails outright: 85.0 %, against a registered
89 % and a measured 81.6 % — i.e. arm A's failed-rung share is **above** its own
measured control.

## F4 — "no ladder" is a theorem about the operator — **PARTIAL**

The affinity argument is **correct for this deck**. Driver: `LadrunoBrick` with
`-geom linear` and `-formulation bbar` (`:597-599`), one `nDMaterial`, constant body
force `-b 0 0 -GAMMA`, prescribed settlement via `ops.sp(t, 3, 1.0)` under a `Linear`
series with `LoadControl` (`:685-687`), `constraints Transformation`. No contact, no
corotational/finite geometry, no u-p coupling, no follower load. With `dσ~/dε = Ce(p_n)`
constant in the step the residual is affine in the free DOFs. `implex_noctl` records
`nfail = nsub = nrelax = 0` and an 8-line engine log. ADR §2 item 2 reads
"**A global step that is linear.** Newton converges in one iteration on a frozen
operator; the ladder never fires…" and carries no measurement marker (unlike item 1's
"**Unmeasured.**" and item 3's "measured at P3") — so "a derivation" is a fair
characterisation of the *mathematics*.

Two corrections. (a) The affinity is not unconditional even on this deck: the P1 floor
clamp `σ~ = dev(σ~) + p_min·I1` is a kink that would make the step piecewise-linear if
it fired. It did not — `n_clamping = 0` in all three payloads — which is a *measurement*,
not a theorem. (b) "A leg that cannot fail" is false of the implementation: the first
attempt, `implex_ON`, ran `-implex` on this deck and produced **0 converged steps**
(`FAILED.txt`: "the push never started"), 34 272 refusal lines, and 7/7/7 rung failures.
The gate demonstrably could and did return a negative on an `-implex` leg. So the
implementation-level claim was open even though the operator-level claim was not.
Red's "measured on the boundary-value problem for the first time" objection stands as
rhetoric; "could not have come out otherwise" does not.

## F5 — the memo mis-cites which ADR claim it upgrades — **CONFIRMED**

`grep -n "warranted by CP1\|unproven at BVP\|BVP" 92_ladruno_sanisand_implex_adr.md`
→ **no matches** (exit 1). `git show --stat 70f6bf0e8` lists 38 files; the ADR is not
among them. So the memo's §3 quotation *"warranted by CP1 but unproven at BVP level"*
is nowhere on the branch. ADR §2 item **1** is *"**Bounded work per step.** The
expensive return runs **once per committed step** instead of once per state-determination
pass — up to 125 of them at the seizure. That is a ~100x cut in constitutive work …
**Unmeasured.**"* — its unit is constitutive returns per step. Item **2** is the linear
step / no ladder. The gate measured item 2. No count of returns, substeps or
state-determination passes exists in any shipped artifact (`max_substeps` is an input;
`n_material_refused` is absent by the memo's own admission; `implexError` is absent per
F8). The memo (and the commit message: "ADR sec.2 item 1 moves from 'warranted, unproven
at BVP level' to proven for the ladder term") moves an unmeasured cost claim on evidence
for a different item. The trailing hedge "for the ladder term" is doing all the work and
the item number is wrong.

## F6 — the checkable overlap was never checked — **PARTIAL** (numbers reproduce; two comparisons are range-mismatched)

Overlay reproduces to the digit (linear interpolation of each leg's
`a2_h1.0_e0.6944_curve.csv` on `s_over_B`):

| s/B | q_OFF | q_noctl | q_ctl | noctl/OFF | ctl/OFF | ds_noctl |
|---|---|---|---|---|---|---|
| 0.00001 | 1.34 | 2.74 | 0.95 | 2.050 | 0.715 | 0.020 |
| 0.00100 | 31.88 | 35.10 | 33.05 | 1.101 | 1.037 | 0.320 |
| 0.00500 | 131.84 | 161.58 | 139.76 | 1.226 | 1.060 | 1.280 |
| 0.01000 | 248.85 | 322.05 | — | 1.294 | — | 2.560 |
| 0.02000 | 477.98 | 657.65 | — | **1.376** | — | 5.000 |
| 0.04000 | 926.47 | 933.00 | — | 1.007 | — | 5.000 |
| 0.05530 | 1261.86 | 1273.05 | — | 1.009 | — | 5.000 |

Per-step over the overlap: **noctl vs OFF mean |dev| 11.36 %, max +104.98 % (step 1),
max excluding step 1 +38.06 % at s/B = 0.0203**; **ctl vs OFF mean |dev| 4.80 %,
max positive +8.67 % at s/B = 0.00032**. Non-monotone confirmed (+37.6 % at 0.02 →
+0.25 % at 0.0378). Field L2 (volume-weighted relative L2 of `eps_q_p`):
**6.58 / 7.20 / 8.22 / 9.35 / 7.35 %** at s/B 0.002/0.005/0.01/0.02/0.04 — matches.
Element 120 (`xc = 1.5, zc = -0.5`) at the 0.02 checkpoint: `eps_q_p` 0.011343 (OFF)
vs 0.010074 (noctl) = **−11.19 %**, at `q_foot` 484.69 → 669.14 kPa = **+38.05 %**.
Static check after step 5: OFF max|q_foot+q_base| = **0.3893 kPa** (0.0548 % of q),
noctl **22.3398 kPa** (0.5063 % of q).

**Comparability of `q_foot`: yes.** Same driver process, same extraction
(`ops.reactions()`, `−(Σ nodeReaction(foot,3) − r0_foot)/area`, `:751-753`), and the
field checkpoints are at *identical* actual settlements (0.002180 / 0.005060 / 0.010180 /
0.020300 / 0.040300 in both legs). Red's "no defensible reason they are not comparable"
holds.

**Two range-mismatches in red's presentation.** (i) The static-check pairing compares
maxima taken over different ranges: noctl's 22.34 kPa occurs at s/B = 0.25, far outside
the overlap. Restricted to the overlap (step ≥ 5, s/B ≤ 0.0553) the numbers invert on
the relative axis — OFF max 0.3893 kPa / **3.99 %** of q, noctl 13.96 kPa / **2.89 %**.
"9.2x worse relative" is an artefact. (ii) The 4.80 % vs 11.36 % headline compares
ctl over s/B ≤ 0.005 against noctl over s/B ≤ 0.0553. On the matched window
(s/B ≤ 0.00504): noctl mean |dev| **11.16 %**, ctl **4.80 %** — the ordering survives,
so the conclusion holds — but red reports noctl's step-1 outlier (+105 %) and suppresses
ctl's (−28.5 %); ctl's max |deviation| is **28.53 %**, not 8.67 %.

## F7 — the deep 78 % — **CONFIRMED on 4 of 5 bullets; bullet 4 REFUTED**

Confirmed: **78 of 78** steps past s/B = 0.0553 carry `ds_mm = 5.000` (= `ds_max`);
`implexError` never recorded (F8); the checkable error is non-monotone (F6); and the
final step is **q-decreasing by 65.65 kPa** (4478.303 → 4412.654) at s/B 0.2478 → 0.2500
— only the second q-decreasing step in the leg (the other is step 2, F10c), and the
"TARGET" endpoint sits on it. (Pedantic: the `ds_mm` column records the *post-growth*
`ds`, so the final clipped step actually advanced 4.4 mm; immaterial, 77/78 are exact.)

**Bullet 4 is REFUTED as filed.** ADR-79 §9's 1108–1152 kPa is a **different problem**:
a *square* footing on a Drucker-Prager cone (PDMY-derived), non-associated, 14.5 B
clearance / 10 B depth, no surcharge. ADR-92's deck is a *plane-strain strip* footing on
`LadrunoSANISAND` at `e = 0.6944` with a 10 kPa surcharge on the R3 60×20 m domain.
Different footing shape, different material, different confinement — the two capacities
are not commensurable and "~4x any measured collapse load" does not follow. The
internal comparison points the other way: the control leg is at 1261.9 kPa at s/B
0.0553 with `t_last = 10 815 kPa/m` still climbing; a tangent extension of the *control*
to s/B 0.2478 gives ≈ 5 430 kPa, so IMPL-EX's 4 478 kPa is **below**, not above, what
the implicit leg's own terminal slope projects. The fork's own T8 legs on this mesh
reach `q_u` 1035–2177 kPa by s/B 0.03–0.06 without plateau. The deep range is
unverified — that is bullet 1–3 and it stands — but it is not "~4x any measured
collapse load".

## F8 — `implexError` in no artifact — **CONFIRMED**

`grep -c implexError` on the driver (both the branch copy and the copy actually used):
**0 / 0**. `grep -ril implexError adr92_bvp/`: three engine logs, one hit each — all the
banner line "READING HAZARD (ADR 92 section 8) … and implexError must be printed beside
the verdict." Zero numeric values anywhere. ADR §8, verbatim:

> **A reading hazard, and it is the serious one.** An IMPL-EX curve satisfies equilibrium
> with the *extrapolated* stress. A limit point called on an IMPL-EX leg must be confirmed
> by the implicit material up to the last settlement the implicit solver reaches, and
> `implexError` must be printed beside every WP1 verdict.

The material computes and exposes it (`ladrunoImplexMeasureError`, W8 in the plan); the
driver was never taught to read it. The gate does not satisfy the ADR's own reporting
mandate. Red's "re-run condition, not a disclaimer" is a judgement call, not a
measurement; the measurement is unambiguous.

## F9 — "the mechanism works, the tolerance as set does not" — **PARTIAL**

Confirmed: over the matched window ctl is the more accurate arm (mean |dev| 4.80 % vs
11.16 %), and ADR §8 predicts the seizure verbatim —

> The corner Gauss point sees 10–100x the nominal, so `-implexControl` at `0.05` will
> halve the step exactly where the wall was. Correct behaviour, bounded by the reduction
> limit — **IMPL-EX trades the wall for a cost there rather than removing it.** … so at
> the corner `-implexControl` is a requirement, not an option.

The memo does not cite that prediction. Confirmed: the two arms are **not** "one flag
apart" in trajectory — ctl's first successful step is `ds = 0.00125 mm` (curve row 1:
`s_m = 1.25e-06`), four halvings below `ds_base = 0.02 mm`, at which OFF and noctl both
start.

Not confirmed: (a) ctl's error is not uniformly small — max |dev| **28.53 %** on its
first step, and it collapses at s/B 0.005, so "reproduces the implicit curve to <9 %" is
red's own suppression of the outlier they charged noctl for; (b) "nothing in the data
supports a usable middle" is an argument, not a measurement — no tolerance other than
0.05 was run, so the data is equally silent for and against a middle, which is exactly
what the memo says ("it has not been found") and what its owed-item 1 proposes to
measure; (c) "loosening the tolerance admits more extrapolation error by construction"
is true and does not show that the admitted error is unacceptable at 0.2, because
`implexError` was never recorded (F8). Red is right that the memo reads a fulfilled
prediction as a mis-set knob; red overreaches in declaring the sweep pointless before it
is run.

## F10 — `tail = 95.9 %` rests on a four-point fit — **CONFIRMED**

Driver `:844-845`: `n0 = max(4, len(s) // 50)`; all three legs have `n0 = 4`.

| leg | first 4 `q` (kPa) | fit window (m) | `t_init` |
|---|---|---|---|
| OFF | 1.336, 1.836, 2.295, 2.768 | 2e-5 … 8e-5 | **23 773.3** |
| noctl | **2.739, 2.231**, 2.595, 2.978 | 2e-5 … 8e-5 | **5 404.2** |
| ctl | **0.462, 0.349**, 0.378, 0.425 | 1.25e-6 … 5e-6 | **−6 709.3** |

Matched-window refit, `s ≤ 1 mm` (19 / 19 / 42 points): **16 745.0 / 17 581.6 /
18 848.6 kPa/m** — the three legs agree to 12.6 %. Recomputed tails:
`t_last/t_init(OFF)` = **21.80 %**; `t_last/`matched-noctl = **29.48 %**; against the
reported 95.92 % — inflated **3.3–4.4x**. The estimator returns a negative initial
stiffness on ctl (`tail_pct = −240.03`) and fits the three legs over *different*
settlement windows (ctl's window is 16x shorter), so it is not even a like-for-like
statistic. Step-2 `q`-drop confirmed: 2.739 → 2.231 (**−18.6 %**) under increasing
prescribed settlement, while the control is monotone over the same four steps.
The qualitative claim survives at any denominator (21.8 % ≫ `PLATEAU_FRAC = 2 %`); the
number 95.9 %, quoted three times, does not.

## F11 — provenance — **CONFIRMED**

1. `build: 'any'` in both IMPL-EX payloads, curve headers and field headers.
   `EXPECTED_BUILD` (`:234`) is a single variable used both for the assert (`:290-291`)
   and as the recorded `build` field (`:842`, `:861`, `:925`, `:944`) — confirmed
   mechanism. The gate prints it: `builds: implex ['any'], baseline ['80e65a4de']`.
2. Banner, first line of both IMPL-EX `run.log`s, verbatim:
   `*** LADRUNO_A2_EXPECT_BUILD=any: build hash NOT pinned (running 8fe4f5630...); these
   numbers are not comparable to the reported campaign's ***`. The memo does not quote it.
3. **4896 confirmed 7x low.** `implex_ON/a2_h1.0_e0.6944_engine.log` contains **34 272**
   lines "`-implex REFUSES this step -- the pseudo-time increment is negative`", made of
   **exactly 4896 occurrences of each of 7 distinct `dt` values** (−1e-05, −2e-05,
   −5e-06, −2.5e-06, −1.25e-06, −6.25e-07, −3.125e-07) — one per subdivision level.
   34272 = 7 × 4896. The memo's "fired 4896 times before" is the per-level count.
   `grep -rn 4896` over `Ladruno_implementation/` and `Ladruno_files/` finds it nowhere
   else. (Fair to the memo: 4896 is a real count in the log, just not the run's.)
4. Cross-build: baseline `80e65a4de`, treatment `8fe4f5630` + uncommitted `3c788778f`.
   `implex_ON` also ran on `80e65a4de` — so the "before/after" pair is itself a
   two-binary comparison. The memo does not say the baseline and treatment legs are
   different binaries.

## F12 — the residuals — **PARTIAL** (four of six confirmed, two overstated)

- **Disclaimer coverage** — the three disclaimers are silent on the PARTIAL verdict,
  arm A's 75 rung failures, the overlap deviation, the three deviations, the missing
  `implexError` and the `t_init` defect. Confirmed as a matter of text.
- **Cost model edited in the results commit** — confirmed. `git diff 5e5ec4db7 70f6bf0e8
  -- adr92_bvp_gate.py` adds `refused` and cites "measured 2026-09-06, arm A came back
  85.0 %" in the new comment. `past1` and both thresholds are untouched. The edit is
  inert on this data (`n_material_refused` absent ⇒ `refused = 0`). All confirmed.
  Blue note: the edit's *rationale* is now independently supported — arm A has zero
  convergence failures (F2) — so "it changes only the narrative" understates it; it
  fixes a real mis-costing that this data happens not to exercise.
- **Wall-clock at matched settlement** — reproduces exactly:
  **11.0x** at s/B 0.002, **24.9x** at 0.005, **37.8x** at 0.010, **51.1x** at 0.020,
  **75.7x** at 0.040, **88.2x** at 0.0553; 0.606 s/step (noctl) vs 37.53 s/step (OFF).
  The memo forbids the number and the number is one interpolation away. Confirmed.
- **Concurrency** — `run.log` dates: OFF 01:29:32, ON 01:29:30, ctl 02:12:02,
  noctl 02:12:05; ctl's JSON stamps 02:12:23 and noctl's 02:13:32. So ctl and noctl
  overlapped for ~18 of noctl's 86 s — real contamination of the ratios above.
  **Overstated for the baseline**: `implex_ON` exited within the same minute it started
  (`run.log`: `exit=0 01:29`, 0 converged steps), so OFF ran essentially alone for its
  2400 s and its WALL termination is not meaningfully contaminated by it.
- **Cited baseline not in the PR** — confirmed. `_adr92_cp1_surcharge_results.md` is not
  in `git ls-files Ladruno_implementation` on this branch (it lands on `6eca6d774`,
  `wp/92b-cp1-surcharge`), yet the memo lists it under `related:` and quotes its
  61 %/89 % throughout.
- **Field dumps cannot localise the error** — columns are exactly
  `ele,xc,zc,hx,hz,vol,eps_q_p,eta_grav`; `eta_grav` is bit-identical between the
  s/B 0.002 and s/B 0.15 checkpoints, i.e. it is the frozen gravity state. No `p`, `q`,
  `η/M`, `e` or `α` at the reported settlement. Confirmed.

## 13. What red missed — the reproducibility hole is larger than F11 says

The driver that produced every leg is **not on the branch**. `run.log` and the JSONs
name it as
`…\worktrees\ladruno-sanisand-implex-adr-31cfcf\Ladruno_files\testbed\hypo_bearing\sanisand_tau0_band.py`;
that worktree is on branch **`wp/92b-cp1-surcharge`** and the file there is
**uncommitted** (`git status` → ` M `). The copy committed on `wp/92c-implex-p1` has no
`--implex`, no `--implexControl`, no `--surcharge`, no `--xlim/--zbot` — `run_leg`'s
signature ends at `verbose=True` and the `nDMaterial` call ends at `-maxSubsteps`
(diff verified). PR #798 therefore ships data that **cannot be regenerated from PR #798**,
by anyone, including its author after that worktree is pruned. Red cited driver line
numbers (`:754`, `:845`, `:234`) from the uncommitted copy without noticing they do not
resolve on the branch under review.

---

## What the data licenses

**On the ladder claim:** with the error control OFF, every one of 142 converged steps
solved on rung 1 with zero subdivisions, which is what an affine residual on a frozen
`Ce(p_n)` must do and is now shown to hold in the implementation (the clamp never fired,
`n_clamping = 0`; and the same implementation failed outright once, so the check was
not vacuous). With the control ON at 0.05, every converged step also solved on rung 1,
but 25 further step attempts entered and abandoned all three rungs on a material refusal
— no Newton iteration failed to converge in that arm. So the ladder-as-*iteration-cost*
is gone on both arms; the ladder-as-*control-flow* still runs, futilely, on the
registered one, and the memo's unqualified "not one step left rung 1" is false while its
§2 "not a ladder cost" is right.

**On the deep curve:** nothing past s/B = 0.055 is licensed. Over the 0–0.055 overlap
the uncontrolled leg differs from the implicit leg by 0.25–38 % (mean 11.4 %,
non-monotone, peaking at s/B 0.02), carries 6.6–9.4 % less plastic strain in volume-
weighted L2 while carrying up to 38 % more load, and all 78 deeper steps ran at `ds_max`
with `implexError` never recorded — so the error there is unmeasured, not small, and the
one number the ADR mandates beside the verdict does not exist. `q_u = 4478 kPa` is
neither absurd (it is below what the control's own terminal tangent projects) nor
supported; it is a solver trace, and `tail = 95.9 %` must be withdrawn in favour of
21.8–29.5 % on a defensible denominator.
