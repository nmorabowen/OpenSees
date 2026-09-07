# RED TEAM 2 — measurement attack on the ADR-92 P1 "BVP gate"

Target: PR #798 / `wp/92c-implex-p1`, results commit `70f6bf0e8`.
Worktree read: `C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr92-p1-implex`.
All numbers below were recomputed from the committed artifacts; every command is cited.
Read-only: nothing was edited, built, or re-run in OpenSees.

---

## Overlay table (the memo contains no such table)

`q_foot` (kPa) at matched settlement, linear interpolation of each leg's
`a2_h1.0_e0.6944_curve.csv` on `s_over_B`. `ds_noctl` is the step size the
uncontrolled arm was running at that point.

| s/B | q_OFF | q_noctl | q_ctl | noctl/OFF | ctl/OFF | ds_noctl (mm) |
|---|---|---|---|---|---|---|
| 0.00001 | 1.34 | 2.74 | 0.95 | **2.050** | 0.715 | 0.020 |
| 0.00060 | 20.64 | 22.90 | 21.45 | 1.109 | 1.039 | 0.160 |
| 0.00100 | 31.88 | 35.10 | 33.05 | 1.101 | 1.037 | 0.320 |
| 0.00200 | 58.16 | 63.75 | 59.91 | 1.096 | 1.030 | 0.640 |
| 0.00300 | 83.27 | 91.34 | 86.36 | 1.097 | 1.037 | 0.640 |
| 0.00500 | 131.84 | 161.58 | 139.76 | 1.226 | 1.060 | 1.280 |
| 0.01000 | 248.85 | 322.05 | — | 1.294 | — | 2.560 |
| 0.02000 | 477.98 | 657.65 | — | **1.376** | — | 5.000 |
| 0.03000 | 703.55 | 743.36 | — | 1.057 | — | 5.000 |
| 0.04000 | 926.47 | 933.00 | — | 1.007 | — | 5.000 |
| 0.05530 | 1261.86 | 1273.05 | — | 1.009 | — | 5.000 |

Summary statistics over the overlap `s/B <= 0.0553` (step-by-step, not interpolated
onto a coarse grid): **noctl vs OFF: mean |deviation| 11.36 %, max +104.98 %** (step 1),
**max +38.1 % at s/B = 0.0203**. **ctl vs OFF: mean |deviation| 4.80 %, max +8.67 %.**

---

## 1. The gate's own instrument returns **PARTIAL**, not "PREDICTION MET" — BLOCKER

**CLAIM:** memo frontmatter, `status: "PREDICTION MET on the ladder claim"`; title
"the ladder is gone"; callout "**The ladder claim is CONFIRMED.**"

**EVIDENCE.** I ran the gate script exactly as committed, on exactly the committed
payloads (`--implex` = both IMPL-EX arms, `--baseline` = `implex_OFF`):

```
  past rung 1          48.4 %  ->     0.0 %
  failed-rung iters    81.6 %  ->    42.5 %

  VERDICT: PARTIAL -- between the prediction and a refutation. Report the
  table as measured; do not round it toward either.
```

The PASS branch in `adr92_bvp_gate.py:main` is `if ip <= PASS_PAST_RUNG1 and ifs <=
PASS_FAILED_SHARE` — i.e. **both** clauses, 10 % each. `ifs` = mean(85.0, 0.0) = 42.5 %.
The script therefore refuses to print "PREDICTION MET" on this data, and it prints an
explicit instruction against exactly what the memo did: *"do not round it toward either."*

The memo escapes this in §2 by declaring the failing half of the criterion "not
meaningful" — but that half was pre-registered (`_adr92_p1_execution_plan.md:52`,
"failed-rung iterations 89 % → single digits") and its input (`n_material_refused`)
still does not exist, so the disqualification is an assertion, not a measurement.

**WHAT THE DATA LICENSES:** "The BVP gate returned PARTIAL. One of its two
pre-registered clauses was met on the counted subset; the other was not, and the
authors believe (without the instrument to show it) that the failing clause is an
artefact." Not "PREDICTION MET", not "the ladder is gone".

**WHAT WOULD REFUTE ME:** a run of the committed gate on the committed data that
prints "PREDICTION MET", or a pre-run (pre-`70f6bf0e8`) written statement that clause 2
was advisory.

---

## 2. "Not one step in either IMPL-EX leg left rung 1" is **false** — the engine log the PR ships records 25 rung-2 and 25 rung-3 failures — BLOCKER

**CLAIM:** "**Not one step in either IMPL-EX leg left rung 1.**" and §3 MAY-say:
"`-implex` removes the solver ladder entirely: 0 of 252 steps across both arms needed
a fallback rung".

**EVIDENCE.** Counting the algorithm-failure lines in the committed engine logs:

```
grep -c "WARNING <algo>::solveCurrentStep" adr92_bvp/*/a2_h1.0_e0.6944_engine.log
                    OFF   ctl   ON
NewtonRaphson         2    25     7
NewtonLineSearch     10    25     7
AcceleratedNewton     0    25     7      (= KrylovNewton, rung 3)
```

The `implex_ctl` arm entered **rung 2 twenty-five times and rung 3 twenty-five times**.
That matches its own payload exactly: `nfail = 75 = 3 x nsub = 3 x 25`. The literal log
sequence for the first event is `NewtonRaphson` fail → `NewtonLineSearch` fail →
`AcceleratedNewton` fail → subdivide. The driver's own verdict string agrees:
`"step collapsed to the DS_MIN floor ... (every ladder rung failed at ds = 0.000156 mm)"`.

The "0 %" is a **survivorship denominator**. `steps = len(rows)` counts only *successful*
steps (`sanisand_tau0_band.py:754`), so the gate's `past1 = (steps - rung1)/steps`
structurally cannot see a step that failed all three rungs. The control leg loses
nothing to this (`nsub = 0`); the ctl arm loses 25 attempts. Recomputed on **step
attempts**, which is the apples-to-apples denominator:

| leg | attempts | attempts using a fallback rung | failed rungs / attempt |
|---|---|---|---|
| OFF (control) | 64 | 31 (48.4 %) | 47/64 = **0.734** |
| implex + control 0.05 | 135 | 25 (**18.5 %**) | 75/135 = **0.556** |
| implex, control OFF | 142 | 0 (0.0 %) | **0.000** |

With the error control on — the configuration the pre-registered gate specified, and
the one the material's own banner calls *"a requirement, not an option"* at a
low-confinement corner — IMPL-EX cuts failed rungs per attempt by 24 %, not to zero.

**WHAT THE DATA LICENSES:** "With `-implexControl` OFF the ladder never fires. With it
ON at 0.05 the ladder still fires on 18.5 % of attempts and exhausts all three rungs
every time it does." The unqualified sentence must be withdrawn.

**WHAT WOULD REFUTE ME:** evidence that `AcceleratedNewton::solveCurrentStep` /
`NewtonLineSearch::solveCurrentStep` warnings are emitted without the rung actually
being attempted.

---

## 3. The headline arm is **not** the pre-registered configuration; three undisclosed protocol deviations — BLOCKER

**CLAIM:** memo §1, "Written before the run, as a prediction with a refutation clause".

**EVIDENCE.** `_adr92_p1_execution_plan.md:50-53` registers the gate as:

> "**BVP gate (decisive).** CP1's deck (`h1.0_e0.6944`, `Q = 10`, **cap 1000**) with
> `-implex -implexControl`, same ladder decomposition as CP1 §6. **Prediction:
> past-rung-1 61 % → near zero; failed-rung iterations 89 % → single digits.**"

Deviations, none disclosed in the memo:

1. **Configuration.** The registered arm is `-implex -implexControl`. That arm
   (`implex_ctl`) *seized at s/B 0.005*. The headline sentences ("the ladder is gone",
   "reached TARGET", "142 steps, zero subdivisions") come from `implex_noctl`, an
   **unregistered** arm with the control off.
2. **Substep cap.** Registered `cap 1000`; run at `-maxSubsteps 20000` (`max_substeps:
   20000` in all three payloads) — a 20x change to a variable the sibling commit
   `5e5ec4db7` itself titles *"the substep cap is a CONTROLLED variable"*. The
   `adr86b_t8` archive in this same tree holds `cap1000` and `cap10000` legs of this
   exact deck reaching different depths (0.0588 vs 0.0478), so the cap demonstrably
   moves the result.
3. **Baseline.** Registered baseline: CP1's legs at 61 % / 89 %. Used baseline: a new
   `implex_OFF` leg at 48.4 % / 81.6 %. The memo reports the delta against the new
   baseline while §1 and the gate docstring quote the registered one.

Also: on the **registered** arm, the second registered clause fails outright —
89 % → 85.0 % is a 4.5 % relative reduction, not "single digits".

**WHAT THE DATA LICENSES:** the registered gate was run in a modified configuration and
its decisive result comes from an arm that was never registered. The pre-registration
protects the `past-rung-1` metric and the 10/40 % thresholds (verified: `git log` shows
`adr92_bvp_gate.py` first committed `5e5ec4db7` 00:44:59, first run started 01:29:32) —
it does not protect the deck, the arm, or the baseline.

**WHAT WOULD REFUTE ME:** a pre-`70f6bf0e8` commit or note authorising the control-OFF
arm as a registered arm and `-maxSubsteps 20000` as the registered deck.

---

## 4. "No ladder" on the control-OFF arm is a **theorem about the operator**, not a measurement of IMPL-EX's merit — MAJOR

**CLAIM:** "A frozen SPD operator makes the step linear and a linear step has no ladder
to fail — IMPL-EX's structural claim, **measured on the boundary-value problem for the
first time**."

**EVIDENCE.** `LadrunoSANISAND::ladrunoImplexFreezeTangent` writes `Ce(p_n)` into
`mCe`, `mCep` and `mCep_Consistent`, and ADR §2 states the response as
`sigma~ = sigma_n + Ce(p_n) : ((eps_{n+1} - eps_n) - d_eps_p~)` with
`d sigma~/d eps = Ce(p_n)` *constant in the step*. The deck has no contact, no
geometric nonlinearity, and one material. **The residual is therefore affine in the
free displacements and Newton converges exactly, in one iteration, for every step, for
any deck, at any settlement, no matter how wrong `sigma~` is.** The engine banner even
says so: *"TanType 2 is INERT under -implex."*

Consistently, `implex_noctl` records `nfail = 0, nsub = 0, nrelax = 0` and its engine
log is 8 lines with zero warnings. A leg that cannot fail is not evidence that it
succeeded. ADR §2 already asserts this a priori as item **2** ("A global step that is
linear. Newton converges in one iteration on a frozen operator; the ladder never
fires") — a derivation, not an open prediction.

The non-trivial version of the claim — *does the linear step still converge when the
material is allowed to police itself?* — **was** measured, in the `ctl` arm, and the
answer is no (finding 2: 25 three-rung exhaustions, seizure at 10x shallower).

**WHAT THE DATA LICENSES:** "the implementation behaves as the operator predicts."
Not "IMPL-EX's structural claim measured for the first time" as if it could have come
out otherwise on the arm that produced the headline.

**WHAT WOULD REFUTE ME:** a deck-level nonlinearity (contact, corotational, u-p
coupling, load-dependent BC) in the `sanisand_tau0_band.py` R3 deck that could make the
frozen-`Ce` step nonlinear. I found none.

---

## 5. The memo mis-cites which ADR claim it upgrades; the claim it names (constitutive work, ~100x) was **not measured at all** — MAJOR

**CLAIM:** §3 MAY-say: "The ADR §2 item-1 cost claim, *'warranted by CP1 but unproven at
BVP level'*, is now **proven at BVP level for the ladder term**."

**EVIDENCE.**
- The quoted phrase **does not exist in the ADR**: `grep -n "warranted by CP1\|unproven
  at BVP\|BVP" 92_ladruno_sanisand_implex_adr.md` returns **nothing**. And
  `git show --stat 70f6bf0e8` shows the ADR was **not modified** by the results commit,
  so the memo is quoting a sentence that is nowhere on the branch.
- ADR §2 item **1** is *"**Bounded work per step.** The expensive return runs once per
  committed step instead of once per state-determination pass — up to 125 of them at the
  seizure. That is a ~100x cut in constitutive work ... **Unmeasured.**"* Its unit is
  *constitutive returns per step*. The gate measured *ladder rungs*, which is ADR §2
  item **2**.
- **No count of constitutive returns, substeps, or state-determination passes exists in
  any artifact of this gate.** The payload has `max_substeps` (an input) and nothing
  else; `n_material_refused` is admittedly missing; `implexError` is absent (finding 8).

**WHAT THE DATA LICENSES:** ADR §2 item 2 (already a derivation) is confirmed on the
uncontrolled arm. ADR §2 item 1 remains exactly as the ADR marks it: **Unmeasured**. The
memo's sentence moves an unmeasured cost claim using evidence for a different item.

**WHAT WOULD REFUTE ME:** an ADR revision on this branch containing the quoted phrase,
or a returns-per-step count in the shipped data.

---

## 6. The verification that *was* possible was not done: over the checkable range the uncontrolled arm is 11.4 % off on average and +38 % at worst — MAJOR

**CLAIM:** §3 MUST-NOT: "That the deeper reach is physically right. It is an equilibrium
of the **extrapolated** stress and must be confirmed on the implicit material — **which
currently cannot get past `s/B = 0.055` to do the confirming.**"

**EVIDENCE.** The implicit material *did* reach s/B = 0.0553, which overlaps 22 % of the
IMPL-EX curve — and that overlap was never compared. Overlay table above; recomputed
per-step over the overlap:

- `implex_noctl` vs `implex_OFF`: **mean |deviation| 11.36 %**, **max +38.1 % at
  s/B = 0.0203**, +105 % on step 1, and the deviation is **non-monotone** (rises to
  +37.6 % at s/B 0.02, falls back to +0.7 % at s/B 0.04).
- Plastic-strain fields at the matched checkpoints (`*_field_sB*.csv`, volume-weighted
  relative L2 of `eps_q_p`): 6.6 % / 7.2 % / 8.2 % / **9.4 %** / 7.4 % at
  s/B = 0.002 / 0.005 / 0.01 / 0.02 / 0.04. At the footing-edge element (`ele 120`,
  x = 1.5, z = -0.5) the IMPL-EX arm carries **11.2 % less** plastic strain at
  s/B = 0.02 while carrying **37.6 % more** load — under-plasticised and over-stressed,
  the signature of an unpoliced elastic extrapolation.
- Global static check `max|q_foot + q_base|` after step 5: OFF **0.389 kPa**
  (0.055 % of q), noctl **22.34 kPa** (0.506 % of q) — 9.2x worse relative,
  57x worse absolute.

The framing "cannot get past 0.055 to do the confirming" describes work that was
available and skipped, not work that was impossible.

**WHAT THE DATA LICENSES:** "Over the 0–0.055 s/B range where both legs exist, the
control-OFF IMPL-EX curve differs from the implicit curve by 0.7–38 %, mean 11.4 %,
with no monotone trend. Nothing in that range bounds the error beyond it."

**WHAT WOULD REFUTE ME:** a defensible reason the two legs' `q_foot` are not comparable
at matched settlement (they use the same nodes, the same reaction sum, the same
`r0_foot` subtraction, and the same `ds` schedule from step 1 to step 43).

---

## 7. The deep 78 % of the curve has **no evidentiary value**, and the memo's §3 does not say so — MAJOR

**CLAIM:** the table's "**termination TARGET @ 0.25000**" and "the first leg in this
entire campaign ... to finish its push", plus §4.3's framing that "what licenses the
other 78 % of the curve is an open question".

**EVIDENCE.** It is not an open question on this data; it is settled negatively.
1. All **78 of 78** steps past s/B = 0.0553 ran at `ds = 5.000 mm` (`ds_max`) — the
   largest step in the schedule, and the regime in which the checkable part of the curve
   shows the +28–38 % excursion (overlay table, s/B 0.010–0.023).
2. `implexError` was never recorded anywhere (finding 8), so there is no per-step bound
   on the error in that range.
3. The error in the checkable range is non-monotone, so "it came back to +0.7 % by
   s/B 0.04" licenses no extrapolation.
4. `q_u = 4478 kPa` at s/B = 0.2478 on a `Q = 10 kPa` surcharge with all 200/200
   elements yielded and the curve still climbing at ~5183 kPa/m. For context, the fork's
   own ADR-79 §9 bearing campaign measured 1108–1152 kPa against a 1525 kPa Davis
   anchor. A number ~4x any previously measured collapse load, with the whole domain
   plastic and no plateau, is not a datum awaiting a licence — it is unsupported.
5. The last recorded step of the leg is **q-decreasing by 65.6 kPa** (4478.3 → 4412.7),
   i.e. the "TARGET" endpoint is on a descending step.

**WHAT THE DATA LICENSES:** "The IMPL-EX leg produced 78 unverifiable points at the
largest step size in the schedule, with no error measurement. They should be reported as
a solver trace, not as a load–settlement curve, and no `q` past s/B = 0.055 should
appear in the ADR."

**WHAT WOULD REFUTE ME:** an `implexError` trace over the deep steps showing it bounded
and small, or a coarser-step convergence study on the deep range.

---

## 8. ADR §7/§8's own mandate was violated: `implexError` is not in a single artifact — MAJOR

**CLAIM:** memo disclaimer 2 — "**That leg ran with `-implexControl` OFF**, so nothing
policed the extrapolation error"; and §3 MUST-NOT closing "That asymmetry is itself a
finding and is the ADR §8 reading hazard, live."

**EVIDENCE.** ADR §8 states: *"a limit point called on an IMPL-EX leg must be confirmed
by the implicit material ... and `implexError` **must be printed beside every WP1
verdict**."* The material computes it at every commit
(`ladrunoImplexMeasureError`) and exposes it via `implexDetail`. Yet:

```
grep -c implexError sanisand_tau0_band.py            -> 0
grep -ril implexError adr92_bvp/                     -> only the banner text
grep -c implexError implex_noctl/..._engine.log      -> 1  (the READING HAZARD sentence)
```

Zero numeric `implexError` values were recorded in this gate. The memo's disclaimer
understates the problem: the error was not merely *unpoliced*, it was **not measured**,
on a leg whose entire claim to interest is how far it went. The driver was never taught
to read the response the C++ was built to provide.

**WHAT THE DATA LICENSES:** the gate does not satisfy the ADR's own reporting mandate,
and the single number that would have quantified the hazard was available and skipped.
This is a re-run condition, not a disclaimer.

**WHAT WOULD REFUTE ME:** an `implexError`/`implexDetail` record anywhere in the PR.

---

## 9. "The mechanism works, the tolerance as set does not" inverts what the two arms show — MAJOR

**CLAIM:** "the two arms bracket the truth: **the mechanism works, the tolerance as set
does not.** The operating point is somewhere between 'refuses everything' and 'polices
nothing', and it has not been found."

**EVIDENCE.**
- The **controlled** arm is the accurate one: mean |deviation| from the implicit control
  **4.80 %**, max **8.67 %** (overlay table: 1.030–1.060 through s/B 0.005). The
  **uncontrolled** arm is 11.36 % mean / 38.1 % max. The control is not
  malfunctioning — it is refusing exactly the steps whose error is large, and the arm it
  produced is 2.4x closer to the implicit answer.
- The seizure is **ADR §8's own prediction, fulfilled**: *"The corner Gauss point sees
  10–100x the nominal, so `-implexControl` at `0.05` will halve the step exactly where
  the wall was. Correct behaviour, bounded by the reduction limit — **IMPL-EX trades the
  wall for a cost there rather than removing it.**"* The memo reads the confirmation of
  a predicted failure mode as a mis-set knob and does not cite the prediction.
- Nothing in the data supports a usable middle. Loosening the tolerance admits *more*
  extrapolation error by construction; the memo proposes 0.2 / 0.5 / 1.0, i.e. 4x–20x
  the tolerance whose refusals were already producing a 4.8 %-accurate curve. There is
  no measurement that a tolerance exists which both reaches depth and bounds the error.
- The "bracket" is also not a bracket: the two arms differ in **two** flags and, from
  step 1, in step-size trajectory (ctl's first successful step is `ds = 0.00125 mm`,
  4 halvings below `ds_base`; OFF and noctl both start at `0.02 mm`). "One flag apart"
  in the memo's table header is wrong for the ctl column.

**WHAT THE DATA LICENSES:** "With the error control on, `-implex` reproduces the
implicit curve to <9 % but cannot advance past s/B 0.005 at the free-surface corner —
which is ADR §8's predicted behaviour. That is evidence against IMPL-EX's applicability
to this corner, not evidence of a mis-set tolerance."

**WHAT WOULD REFUTE ME:** a tolerance sweep showing an operating point that both reaches
s/B > 0.05 and keeps `implexError` bounded. It has not been run (memo §4 item 1).

---

## 10. `tail = 95.9 %` is quoted three times and rests on a four-point fit that is noise — MAJOR

**CLAIM:** "**It still did not plateau.** `tail = 95.9 %` of the initial slope at
`s/B = 0.25`." (repeated in the table and in the commit message).

**EVIDENCE.** `sanisand_tau0_band.py:845`: `n0 = max(4, len(s)//50)`; every leg here has
`n0 = 4`. `t_init` is therefore a straight-line fit through the **first four converged
steps**, where `q` is 0.3–3 kPa:

| leg | first 4 q (kPa) | window (m) | `t_init` |
|---|---|---|---|
| OFF | 1.336, 1.836, 2.295, 2.768 | 2e-5 … 8e-5 | 23773 |
| noctl | **2.739, 2.231**, 2.595, 2.978 | 2e-5 … 8e-5 | 5404 |
| ctl | **0.462, 0.349**, 0.378, 0.425 | 1.25e-6 … 5e-6 | **−6709** |

Two things follow.

(a) **The estimator is unfit.** It produced a *negative* initial stiffness for the ctl
arm and a `tail_pct = −240 %`. The memo does not retire the statistic; it blanks the ctl
cell with "—" and keeps quoting the same statistic's output for noctl.

(b) **The 4.4x `t_init` gap is an artefact, not a different problem.** Refitting all
three legs on a *matched* window, `s <= 1 mm` (19 / 19 / 42 points): **16745 / 17582 /
18849 kPa/m** — the three arms agree to 12 %. So the "IMPL-EX is solving a different
problem" reading of `t_init` is not supported; what *is* supported is that the memo's
denominator is noise. Recomputing the headline with a defensible denominator:
`t_last / t_init(OFF) = 5183/23773 = **21.8 %**`; with the matched-window slope,
`5183/17582 = **29.5 %**`. The memo's 95.9 % is inflated **3.3–4.4x**.

(c) The noise has a mechanism worth naming: the first IMPL-EX step has
`d_eps_p(n) = 0`, so `sigma~` is a **pure elastic predictor** — `q` = 2.05x the control's
(2.739 vs 1.336) — and the second step then **drops** `q` by 18.6 % under *increasing*
settlement (2.739 → 2.231). The control leg is monotone over the same 4 steps. A
`q`-decreasing step under monotone prescribed settlement is non-physical and it sits at
the exact points the plateau statistic is built from.

**WHAT THE DATA LICENSES:** "the leg did not plateau" (true at any denominator: 21.8 %
still exceeds `PLATEAU_FRAC = 2 %`). The **number** 95.9 % must be withdrawn, and
`t_init` / `tail_pct` retired or redefined on a matched window.

**WHAT WOULD REFUTE ME:** an argument that a 4-point fit spanning 60 µm of settlement on
a 2 m footing is a meaningful "initial tangent".

---

## 11. Provenance: the driver printed its own non-comparability banner, and the memo's verification number contradicts the shipped log — MAJOR

**CLAIM:** "**Provenance, stated exactly.** Engine built from `wp/92c-implex-p1` after
the D2 sign-guard fix; the binary **contains** that fix but **reports** hash
`8fe4f5630` ... The fix was verified **behaviourally, not by hash**: the negative-`dt`
refusal fired **4896** times before and **0** times after."

**EVIDENCE.**
1. **The payloads carry no provenance at all.** `build: 'any'` in both IMPL-EX legs.
   `EXPECTED_BUILD` (`sanisand_tau0_band.py:234`) is an env var used *both* for the
   assert *and* as the recorded `build` field, so `LADRUNO_A2_EXPECT_BUILD=any` disables
   the pin **and** stamps the literal string `any` into every JSON, curve header and
   field header. The gate script prints the consequence itself:
   `builds: implex ['any'], baseline ['80e65a4de']`.
2. **The driver emitted an explicit self-invalidation the memo does not quote.** First
   line of both IMPL-EX `run.log`s:
   `*** LADRUNO_A2_EXPECT_BUILD=any: build hash NOT pinned (running 8fe4f5630...);
   these numbers are not comparable to the reported campaign's ***`
3. **The 4896 does not match the shipped artifact.** The only "before" leg in the PR is
   `implex_ON`, whose engine log contains **34272** occurrences of
   `-implex REFUSES this step -- the pseudo-time increment is negative`
   (`grep -o ... | sort | uniq -c`). 34272 = 7 x 4896, i.e. 4896 is the count for **one**
   of the 7 subdivision attempts, not the run. `grep -rn 4896` over
   `Ladruno_implementation` and `Ladruno_files` finds the figure nowhere except the fix
   commit message. A behavioural-provenance argument whose one number is 7x off the
   shipped evidence is not "stated exactly".
4. **Cross-build control.** Baseline `implex_OFF` = `80e65a4de`; treatment arms =
   `8fe4f5630` + an uncommitted `3c788778f`. Three commits apart. `3c788778f` touches
   `LadrunoSANISAND.cpp` only inside `ladrunoImplexTrial`, so the risk is low — but the
   comparison is nonetheless between two different binaries and the memo does not say so.

**WHAT THE DATA LICENSES:** "the IMPL-EX legs were run on an unpinned binary that the
driver itself flagged as non-comparable, and the behavioural check reported in the memo
does not reproduce from the committed log."

**WHAT WOULD REFUTE ME:** a re-run with `LADRUNO_A2_EXPECT_BUILD` pinned to a real
hash, on the same binary as the baseline.

---

## 12. The disclaimers are calibrated on the axes the authors chose; the cost model was edited post-run; and n = 1 — MINOR (three residuals)

**CLAIM:** "**Three things this does NOT say, each of which must be quoted with it**";
and §2 "Rather than manufacture a number, arm A's failed-rung column is reported as not
meaningful."

**EVIDENCE.**
- **The disclaimers are silent on every axis found above**: the ctl arm's 75 failed
  rungs, the 11.4 %/38 % curve deviation in the checkable overlap, the three
  pre-registration deviations, the PARTIAL verdict, the missing `implexError`, and the
  `t_init` defect. Each of these is a *stronger* negative than the three that were
  disclosed. They are careful, but they are careful about what the authors already knew.
- **The cost model was changed in the results commit itself.** `git diff 5e5ec4db7
  70f6bf0e8 -- adr92_bvp_gate.py` adds the `refused` split, in a direction that lowers
  the IMPL-EX arm's failure share, with the measured 85.0 % quoted in the new comment as
  the motivation. The `past1` metric and the 10/40 % thresholds are untouched (good), but
  the clause the verdict fails on was edited after seeing it fail. And the edit is inert:
  `n_material_refused` still does not exist, so `refused = 0` and the number is unchanged
  — it changes only the narrative around it.
- **The memo's forbidden number is computable and large.** §3 MUST-NOT: "Anything about
  wall-clock speedup: the arms had different terminations and different work, and **no
  timing comparison at matched settlement was made**." The curve CSVs carry `wall_s` per
  row; the comparison takes one interpolation: **11.0x at s/B 0.002, 24.9x at 0.005,
  37.8x at 0.010, 51.1x at 0.020, 75.7x at 0.040, 88.2x at 0.0553** (0.606 s/step vs
  37.5 s/step). The honest pairing is "75x faster at s/B = 0.04, for an answer that
  deviates 0.7–38 % over that range" — a sentence the memo has neither half of.
- **n = 1, and the arms were not isolated.** One deck, one void ratio, one surcharge,
  one tolerance, one mesh, one `ds_max`. `ctl` and `noctl` started at **02:12:02** and
  **02:12:05** (their `run.log` headers) and both finished by 02:13:32 — two threaded
  PARDISO processes on one machine, so neither IMPL-EX arm's `wall_s` is a clean
  measurement. `implex_ON` (which wrote 12.4 MB of warnings) started at 01:29:30, two
  seconds before the baseline `implex_OFF` at 01:29:32, whose 2400 s WALL termination is
  the entire reason it stopped at s/B 0.055.
- **The cited baseline is not in the PR.** `_adr92_cp1_surcharge_results.md` is listed in
  the memo's `related:` and its 61 %/89 %/125-passes numbers are quoted throughout, but
  the file does not exist on this branch (it is on `6eca6d774`). A reviewer of #798
  cannot check the baseline the prediction was written against.
- **The field dumps cannot localise the error.** `*_field*.csv` columns are
  `ele,xc,zc,hx,hz,vol,eps_q_p,eta_grav` — `eta_grav` is the *gravity* state, constant.
  No `p`, `q`, `eta/M`, `e` or `alpha` at the reported settlement, so the requested
  corner-state divergence cannot be computed from the shipped artifacts at all.

**WHAT THE DATA LICENSES:** a re-run, not a rewrite of the prose.

**WHAT WOULD REFUTE ME:** a second deck / void ratio / surcharge reproducing the ladder
result on the *registered* arm, plus process-isolated timings.

---

## Bottom line

The one thing this gate demonstrates is that **the implementation behaves the way a
frozen-`Ce` operator must**: with the error control off, the global step is linear and
therefore cannot fail. That is ADR §2 item 2, which the ADR already asserts as a
derivation.

Everything the memo adds to that is unsupported by its own data:
the gate returns **PARTIAL**; the ladder is **not** gone on the registered
(`-implexControl`) arm, where it fires on 18.5 % of attempts and exhausts all three
rungs each time; the headline arm was never registered and the deck's substep cap was
changed 20x; the cited ADR claim (§2 item 1, constitutive work) was **not measured**;
`implexError` — mandated by ADR §8 — was **never recorded**; `tail = 95.9 %` is inflated
3.3–4.4x by a four-point fit that returns a negative slope on the sibling arm; and the
overlap where verification *was* possible shows the headline curve **11.4 % off on
average and +38 % at worst**, which the memo never computed.
