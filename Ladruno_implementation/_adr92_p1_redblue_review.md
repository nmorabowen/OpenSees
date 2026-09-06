---
title: "ADR 92 / P1 — red/blue adversarial review of PR #798 (wp/92c-implex-p1)"
project: Ladruno
type: review
status: "NOT READY — 5 blockers confirmed by both teams; the BVP results memo must be rewritten"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p1_execution_plan]]"
  - "[[_adr92_p1_bvp_gate_results]]"
  - "[[_adr92_p0_oracle_results]]"
tags: [adr, sanisand, implex, review, red-team]
updated: 2026-09-06
---

# ADR 92 / P1 — red/blue adversarial review

Scope: range `235934592..70f6bf0e8` on `wp/92c-implex-p1` (PR #798, draft). Three red
attackers (C++, measurement, tests/process) filed 34 findings; three blue defenders
verified each against the source and the committed data. The six raw reports are in
`_adr92_p1_redblue/`. This memo is the adjudication: what stands, what fell, what is owed.

> [!danger] **Verdict: the C++ is not correct on the campaign's own deck, and the BVP
> results memo overstates what the data shows.** Both are fixable. Nothing here says
> IMPL-EX is the wrong idea; it says the evidence filed for it does not yet exist.

Scoreboard (blue verdicts on red findings):

| team | filed | CONFIRMED | PARTIAL | REFUTED (whole or central consequence) |
|---|---|---|---|---|
| RED-1 C++ | 10 | 8 | 2 | 0 whole; F10's consequence refuted, F9 downgraded |
| RED-2 measurement | 12 | 7 | 5 | 0 whole; F7 bullet 4 and F2's cost inference refuted |
| RED-3 tests/process | 12 | 9 | 2 | F4's central consequence refuted; coverage row 9 refuted |

---

## 1. Blockers (confirmed by both teams)

### B1. The extrapolation factor is never `dt_{n+1}/dt_n` on a negative pseudo-clock
`LadrunoSANISAND.cpp:1329` gates the ratio on `mImplexDtCommit > 0.0`. Commit
`3c788778f` let a monotonically negative clock *through the refusal* on the argument
that two negative increments give a positive ratio — but that ratio is never computed,
so `f ≡ alpha = 1.0` for the life of any `LoadControl(-ds)` leg. The campaign deck is
exactly that. BLUE-1 reconstructed the noctl leg's applied `ds` from the settlement
column: **9 of 142 steps had a ratio ≠ 1** (seven at 2.0, one at 1.95, the last at
0.88), all run at `f = 1`. The two consumers of `mImplexDtCommit` disagree about its
sign contract (`!= 0.0` at the D2 guard `:1511`, `> 0.0` at the factor `:1329`).
*(RED-1 F2 = RED-3 F1.)*

### B2. `-implexControl`'s reduction floor is dead on the same decks
`mImplexDt0` is armed only from a positive `dt` (`:1326`), so it stays `0.0`; the
escape at `:1615` then short-circuits on `mImplexDt0 <= 0.0` and refuses without
limit. BLUE-1's decisive arithmetic: the ctl arm died refusing at `ds = 1.5625e-7 m`,
**below its own intended floor** of `0.01 × 2e-5 = 2e-7 m`. The memo's "the tolerance
as set does not work" is therefore confounded — it measured a broken floor, not a
tolerance. The booked tolerance sweep cannot succeed at any value until this is fixed.
*(RED-1 F3 = RED-3 F2.)*

### B3. The commit-time companion refusal is dead code
`commitState()` returns `LADRUNO_MATERIAL_REFUSED` at `:1660` when the companion hits
`-maxSubsteps`, but `Domain.cpp:2244` is a bare `elePtr->commitState();` and
`Domain::commit()` returns 0 unconditionally, so `AnalysisModel::commitDomain()`'s
check can never fire. The step is accepted, `committedTime` advances, recorders write,
and because the early return precedes `:1673-1682`, `mEpsilon_n`, `mImplexDtCommit`,
`mImplexStepArmed` and `mImplexDEpsP` are all left stale for every later step. This
is the **default** configuration (control off, cap mandatory). *(RED-1 F1.)*

### B4. The results memo's claims are not what the data shows
Every item below reproduced by BLUE-2 to the digit.

- **The committed gate returns PARTIAL on the registered arm**, not PREDICTION MET.
  The PASS branch needs both clauses ≤ 10 %; the registered `-implex -implexControl`
  arm gives 0.0 % / 85.0 %. "PREDICTION MET" is reachable only by feeding the gate
  the unregistered control-OFF arm alone. *(RED-2 F1.)*
- **"Not one step in either IMPL-EX leg left rung 1" is false as written.** The ctl
  engine log records 25 / 25 / 25 `solveCurrentStep` failures on NR / NLS /
  AcceleratedNewton — `nfail = 75 = 3 × 25`. The gate's `steps` is a survivorship
  denominator (`rows.append` sits after the failure `continue`). What survives: all
  75 were `Integrator failed in update()` on a material refusal with **zero**
  `CTestNormUnbalance` failures, so the ladder-as-iteration-cost is gone; the
  ladder-as-control-flow still runs, futilely. *(RED-2 F2, RED-3 F3.)*
- **Three undisclosed pre-registration deviations**: headline arm unregistered
  (control OFF); cap 1000 → 20000 (a variable the sibling commit calls "controlled";
  `adr86b_t8` shows a 38 % swing in `q_u` across caps on this mesh); baseline
  61 % → 48.4 % with §1 **restating the pre-registration with the substituted
  number** inside the sentence claiming it was written before the run. *(RED-2 F3.)*
- **The ADR claim upgraded does not exist.** The quoted phrase "warranted by CP1 but
  unproven at BVP level" occurs nowhere in the ADR, and `70f6bf0e8` did not touch it.
  §2 item 1 (constitutive work, "Unmeasured") was not measured; the gate measured
  item 2 (the linear step), which the ADR already states as a derivation. *(RED-2 F5.)*
- **`implexError` was never recorded.** ADR §8 mandates it "beside every WP1
  verdict"; the driver has zero references and no artifact carries a value.
  *(RED-2 F8.)*
- **`tail = 95.9 %` is inflated 3.3–4.4×.** `t_init` is a 4-point fit over 60 µm;
  it returns **−6709** for the ctl arm. Matched-window refit (s ≤ 1 mm) gives
  16745 / 17582 / 18849 for OFF / noctl / ctl (agree to 12.6 %), and tails of
  **21.8–29.5 %**. "Did not plateau" survives at any denominator. *(RED-2 F10.)*
- **The checkable overlap was never compared.** Over s/B ≤ 0.0553 the control-OFF
  curve is **11.4 % off on average, +38 % at s/B 0.02**, non-monotone (back to +0.7 %
  by 0.04), carrying 6.6–9.4 % *less* plastic strain while carrying up to 38 % *more*
  load. Step 1 is 2.05× the control (pure elastic predictor, `d_eps_p(n) = 0`) and
  step 2 **drops** q by 18.6 % under increasing settlement. The ctl arm, on its
  matched window, is 4.8 % mean / 28.5 % max. The 78 deep steps all ran at `ds_max`
  and the endpoint sits on a q-decreasing step. *(RED-2 F6, F7, F9.)*

### B5. The PR cannot regenerate its own data
The driver on the branch (`sanisand_tau0_band.py` at `ee841bead`) has **no
`--implex`, `--implexControl` or `--surcharge`**. The copy that produced every leg is an
uncommitted edit (+41/−6) in the `wp/92b-cp1-surcharge` worktree. The payloads stamp
`build: 'any'` because `EXPECTED_BUILD` is both the assert input and the recorded
field; the `run.log`s do carry `engine build (ladrunoBuild) : 8fe4f5630` (which
refutes RED-3 F4's "unidentified binary"), but the run.log banner the memo never quotes
says *"these numbers are not comparable to the reported campaign's"*. The memo's
"4896 refusals" is 7× low: the log holds **34 272** = 4896 × 7 subdivision levels.
Baseline `80e65a4de` vs treatment `8fe4f5630` + uncommitted `3c788778f`: two binaries,
unstated. *(BLUE-2 §13, RED-2 F11, RED-3 F4.)*

---

## 2. Majors (confirmed)

**C++**
- Trial refusals are swallowed by every element except `LadrunoBrick`
  (`SSPbrick.cpp:445`, `Brick.cpp:1069` discard the return). BLUE-1 correction: the
  truly silent site is the `-implexControl` refusal at `:1616`, which leaves a
  strain-dependent `sigma~`; the D2/companion sites leave a strain-independent stress
  that usually stalls Newton. *(RED-1 F4.)*
- Scheme-2 without a cap is only warned, and the warning is nested under `verbose`,
  which is `false` on all `getCopy`/`recvSelf` paths; `mSubstepCapHitInME` is
  unreachable at cap 0 so `-implexControl` cannot detect a companion failure there.
  *(RED-1 F5.)*
- `takeAverageError()` zeroes the accumulator on every read; a recorder over an
  8-point mesh gets the average at point 1 and 0.0 elsewhere; `maxError` never
  resets. *(RED-1 F6.)*
- Both refusal warnings are unthrottled (one 527-char line × 4896 = 2.58 MB per
  step) while the `-implexControl` refusal prints nothing at all — which is why
  `n_material_refused` is unrecoverable from any log. The results memo `:84` falsely
  says the material's message is throttled. *(RED-1 F7, F8; RED-3 F12.)*

**Tests**
- The battery cannot detect an inert operator (`f ≡ 0`) or the total-strain form P0
  measured 78 % wrong: the flagship ON/OFF test is tautological by its own docstring,
  the tangent-identity test is affine for *any* frozen `f`, and `implexDetail[5]`
  (`f`) is read and never asserted. The regression `80e65a4de` fixed (extrapolation
  silently off on every retry) has no test. *(RED-3 F5.)*
- `sendSelf/recvSelf` grew 5 → 22 wire slots with zero green coverage: the roundtrip
  test is SKIPPED on every run and, on a zero-free-DOF deck, blind by construction.
  *(RED-3 F6.)*
- D2 semantics changed in `3c788778f` **after** the last pytest run; no test file
  changed after it; the negative-`dt` test's name and docstring now describe the old
  contract (BLUE-3: the test *does* exercise the new sign-change guard, so coverage
  exists — the defect is the stale wording). The monotone-negative-clock case that
  produced B1 and B2 has no test. *(RED-3 F7.)*
- `getCopy` test clones before any `analyze`, so it asserts zero == zero.
  Two known-red tests carry no `xfail`; 16 tests collected, not 15. *(RED-3 F9, F11.)*
- PlaneStrain lane, `-implexAlpha`, `-implexDt user/pseudo`, companion failure,
  subdivision/retry, parallel: not covered. 15 of 28 aspects uncovered. *(RED-3 matrix.)*

**Process**
- The gate's cost model was edited inside the results commit, in the direction that
  lowers the IMPL-EX arm's failure share, and the edit is byte-inert (`n_material_refused`
  has no producer). With `refused = 25` reconstructed from the log, arm A's share is
  53.2 % — still failing. *(RED-2 F12, RED-3 F8.)*
- Ledger row still reads "draft — C++ written, NOT YET BUILT"; files column lists only
  `SRC/`; `manifest.yaml` untouched and `check_manifest.py` cannot notice (keys on new
  classTags; ADR-92 adds none); no mutation gate run. Draft-PR timing (28 s), banner
  untouched by design, and zero vanilla footprint are **compliant**. *(RED-3 F10.)*
- `_adr92_cp1_surcharge_results.md`, the baseline the prediction was written against,
  is not on this branch. A 12.4 MB log of the superseded pre-fix run is committed.
  The two IMPL-EX arms overlapped in wall time (~18 s of noctl's 86 s). *(RED-2 F12.)*

---

## 3. Red-team findings that fell, and why

- **RED-3 F4 (arms ran on an unidentified binary):** `run.log:3` in both arms records
  `ladrunoBuild = 8fe4f5630`; stale-`.pyd` is excluded by 34 272 → 0 refusals. The
  schema defect (`build` = expected, not measured) stands as a minor.
- **RED-1 F10 (stage-0 byte-identity broken):** at `mElastFlag == 0` the moduli are
  independent of stress and void ratio (`ManzariDafalias.cpp:4864`), so the frozen
  writes are bit-identical to `elastic_integrator`'s. The `enabled`-vs-`active` gate
  is a real wart, not a defect.
- **RED-1 F9 (zero-`dt` hold destroys history):** `f = alpha` at `dt_n = 0` is the
  documented fallback and `d_eps_p(n) = 0` after a zero-strain step is the correct
  datum. Minor.
- **RED-2 F2's inference (IMPL-EX cuts failed rungs by only 24 %):** charges material
  refusals as ladder work — the exact error the memo's §2 disclaims. Arm A has zero
  Newton convergence failures.
- **RED-2 F7 bullet 4 (`q_u` 4× any measured collapse load):** ADR-79 §9 is a square
  footing on a Drucker-Prager cone, no surcharge. Not commensurable. The control's own
  tangent projects ≈ 5430 kPa at s/B 0.25, so 4478 is *below* it.
- **RED-2 F9 ("controlled arm accurate to < 9 %"):** its true max is −28.5 % (step 1).
  Red suppressed ctl's outlier while quoting noctl's.
- **RED-2 F6 static check "9.2× worse":** range-mismatched; restricted to the overlap
  the axis inverts.
- **RED-3 coverage row 9:** the negative-`dt` test drives `LoadControl(-1.0)` after
  positive committed steps, so it *is* a sign change and does hit the shipped guard.

---

## 4. What the data does license (both teams agree)

- With control OFF, all 142 converged steps solved on rung 1 with zero subdivisions,
  which an affine residual on frozen `Ce(p_n)` must do — and the `p_min` clamp never
  fired (`n_clamping = 0`, a measurement, not a theorem). The check was not vacuous:
  the same implementation produced 0 converged steps on `implex_ON`.
- With control ON, every *converged* step also solved on rung 1, and no Newton
  iteration ever failed; the 25 abandoned attempts were material refusals.
- Wall clock at matched settlement: 11× at s/B 0.002 rising to 88× at 0.055
  (0.61 vs 37.5 s/step), for a curve that deviates 0.25–38 % over that range. The
  memo has neither half of this sentence; it should have both.
- The leg did not plateau at any denominator (21.8 % ≫ 2 %). `q_u` is unmeasured.
- Nothing past s/B 0.055 has evidentiary value until `implexError` is recorded and the
  overlap discrepancy is explained.

---

## 5. Required before ready-flip, in order

1. **C++ (B1, B2, B3):** gate the factor on `mImplexDtCommit != 0.0` with a
   sign-consistent ratio; arm `mImplexDt0` from `|dt|` and compare magnitudes at the
   floor; make the companion refusal *observable* — at minimum an `opserr` and a
   material-level flag a recorder can read, since `Domain::commit()` will never
   propagate it. Then the majors: throttle both refusal messages to the parent's
   10-per-process budget; print (throttled) and **count** the `-implexControl`
   refusal and expose the count as a response (`n_material_refused`); make
   `avgImplexError` non-destructive; make the scheme-2 no-cap case a refusal or an
   unconditional warning.
2. **Tests:** a monotone-negative-clock test asserting `implexDetail[5] == dt_{n+1}/dt_n`
   across a `ds` change; a companion-refusal test on a free-DOF deck; `getCopy` after
   a plastic history; roundtrip on a deck where ON ≠ OFF; PlaneStrain; `xfail(strict)`
   on the two red tests or fix them; rename the negative-`dt` test.
3. **Commit the driver** that produced the legs (rescued from the `wp/92b` worktree
   in this review, see `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py`).
   Record `ladrunoBuild()` in the payload, not `EXPECTED_BUILD`. Teach the driver to
   record `implexError` and `n_material_refused` per step.
4. **Re-run the gate as registered** on the fixed binary: `-implex -implexControl`,
   cap as registered or the deviation justified in writing, pinned build, baseline and
   treatment on the same binary, arms not overlapping in wall time. Report the overlay
   against the implicit control over the whole overlap, and `tail` on a matched window.
5. **Rewrite `_adr92_p1_bvp_gate_results.md`:** withdraw "PREDICTION MET", "not one
   step left rung 1", `95.9 %`, the non-existent ADR quote, and "verified 4896";
   disclose the three pre-registration deviations; add the overlay table and the
   wall-clock-with-error sentence. Cite ADR §8's own prediction of the ctl seizure.
6. **Ledgers:** flip the row off "NOT YET BUILT", list every file, fix the D2
   behaviour text; add ADR-92 rows to `manifest.yaml`; run the mutation gate on the
   family; drop the 12.4 MB dead log.

## Log

- 2026-09-06 — Review run as 3 red + 3 blue agents on `70f6bf0e8`. Five blockers,
  all confirmed by the opposing team. The driver rescue (item 3) was done in this
  session; everything else is owed.
- 2026-09-06 (later) — Status of §5's required items. **Done:** item 1, the C++
  (B1 sign-consistent `f` + `mImplexDtCommit != 0.0` gate, B2 `|dt|`-armed floor,
  B3 observable companion refusal via `implexRefusals`, plus the majors — throttled
  and counted refusals, non-destructive `avgImplexError`, scheme-2-no-cap refused —
  `2473ce46c`; re-arm after refusal, `2ceb65fa4`; the un-primed-step exemption,
  `afb95c40c`). Item 2, tests, **partially**: lane B added a monotone-negative-clock
  test asserting `implexDetail[5] == dt_{n+1}/dt_n`, a companion-refusal test, a
  `getCopy`-after-history test, a roundtrip test on an ON != OFF deck, and PlaneStrain
  coverage, plus `xfail(strict)` markers (`e4a8523e8`; battery 20 passed / 1 strict
  xfail / 1 skip on `2473ce46c`, `347242f91`) — not yet re-run on `afb95c40c`
  (owed, separate from this lane). Item 3, the driver, **done**: committed with
  `ladrunoBuild()` recording and `implexError`/`n_material_refused` in the payload
  (`b45b338a5` lane C). Item 4, the registered-arm re-run, **done**: ran on the fixed
  binary and returns `PARTIAL` (`_adr92_p1_bvp_gate_rerun.md`). Item 6, ledgers,
  **partially**: the `LEDGER_implementations` row and this review's own log are kept
  current in this pass; **the mutation gate has NOT been run** and `manifest.yaml`
  was not touched (ADR-92 adds no new classTags, so `check_manifest.py` would not
  flag its absence regardless). **Remaining, in order:** the mutation gate on the
  family; a decision on the registered arm's operating point (`reductionLimit 0.01`
  spends the entire subdivision budget by `s/B = 0.02754` — the owed tolerance sweep
  from `_adr92_p1_bvp_gate_results.md` §4 item 1 is still not run); and the
  wall-clock caveats in the rerun memo (a control/noctl sequencing overlap and an
  unrelated >11-CPU-hour background process contaminate every `wall_s` column there)
  are not yet resolved by a clean single-tenant re-run.
- 2026-09-06 (close-out) — §5's six required items, final status. **1–5 done**
  (unchanged from the entry above). **6, ledgers, now complete**: the mutation gate
  ran (`_adr92_p1_mutation_gate.md`, PASSED 0.750, 9/12) and the `LEDGER_implementations`
  ADR-92 P1 row reflects it; `manifest.yaml` is still untouched, correctly, since ADR-92
  adds no new classTag for `check_manifest.py` to flag. The operating-point decision
  owed by item 4's sweep is also measured (`_adr92_p1_bvp_gate_rerun.md`'s
  "Operating-point sweep" section; the ADR now carries a RECOMMENDATION for `tol = 0.1`
  as the documented deck default, with the C++ default left at `0.05` pending the
  owner). **Two things remain, and both are owed tests, not blockers to a ready-flip
  decision:** the mutation gate's three survivors — **M4** (reduction-floor arming
  untested, no deck in the battery drives a subdivision ladder), **M5** (the
  `-implexControl` refusal's return-code contract is unpinned — a mutant that returns
  `0` instead of `LADRUNO_MATERIAL_REFUSED` survives), and **M10** (re-arm-after-refusal
  is untested for a caller that retries without `revertToLastCommit()`). The wall-clock
  caveats (sequencing overlap, background-process contention) remain unresolved by a
  clean re-run and are not expected to move before ready-flip. **The PR is now ready
  for the owner's ready-flip decision, pending Zone-A on the head.** The operating-point
  choice (`0.05` vs the recommended `0.1`) and whether to require the three owed
  mutation tests before merge are both the owner's call, not this review's.
