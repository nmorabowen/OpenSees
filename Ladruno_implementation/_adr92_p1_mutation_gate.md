# ADR-92 P1 — IMPL-EX mutation gate (ADR-87 D2)

**Verdict: GATE PASSED. Score 0.750 (9 of 12 mutants killed), floor 0.60.**

Suite under test: `tests/test_ladruno_sanisand_implex.py` — 22 collected,
20 passed / 1 skipped / 1 xfailed at baseline.
Feature row: `ND_TAG_LadrunoSANISAND` (`Ladruno_implementation/testbed/manifest.yaml`).
Run 2026-09-06 on `wp/92c-implex-p1`.

This closes the last outstanding item of `_adr92_p1_redblue_review.md` §5
("mutation gate") and answers RED-3 F5 directly: **the post-fix battery DOES
detect an inert operator and DOES detect the total-strain form** — the two
things RED-3 measured it could not (M1 → 3 kills, M3 → 4 kills).

---

## 1. Why this is a hand-mutant lane and not a `mutation_gate.py` family run

`SRC/Ladruno_mutation.h` expresses exactly three modes — `ZERO`, `SCALE`,
`IDENT` — applied at the element/section force and tangent accessors. Every
ADR-92 P1 mechanism is **semantic**: a frozen extrapolation factor, a sign gate
on the pseudo-time clock, an incremental-vs-total stress form, a refusal
sentinel, a counter. None of them is expressible as "zero the resultant", and
`SANISAND` is deliberately **not** in `GATED_FAMILIES` (registering paths there
documents the intended evidence suite; running it before the call sites exist
would score 0.0 and mean nothing).

So this gate uses the **other** mutation lane this fork already runs on this
same feature row — the one ADR-86 PR-3 used for `-honorTolR`, `-Pmin`,
`getCopy(void)` and the PlaneStrain sign flip, whose verdicts are recorded in
the manifest as prose (`PROVEN by mutation: ... fails THIS TEST AND ONLY THIS
TEST`). One minimal source edit at a time in
`SRC/material/nD/LadrunoSANISAND.cpp`, rebuilt, battery re-run, restored.

**Scoring is `mutation_gate.py score`'s own semantics, unchanged**: only
baseline-**passing** tests are detectors (20 of 22 — the strict `xfail` and the
skipped roundtrip carry no information and are excluded from the denominator);
a detector that goes **red** under a mutant is a kill. The score reported here
is the fraction of **mutants killed**, which is the right denominator for a set
of hand-picked semantic mutants; the per-mutant detector lists are in §3.

### The identity check, which is not optional

`mutation_gate.py --expect` exists because a silently failed mutant build makes
the harness re-run the previous binary and report the exact inversion of the
truth. A hand-edited mutant reports `ladrunoMutation() == "none"` by
construction, so `--expect` cannot serve here. The same guarantee was obtained
three ways on **every one of the 13 builds**:

1. **The TU is confirmed recompiled by name** in the build log
   (`Select-String LadrunoSANISAND`) — never inferred from `FASTBUILD: OK`,
   which reports the *link*. (LEDGER_quirks: the `Copy-Item` stale-mtime trap.)
2. **`dist\bin\opensees.pyd` mtime strictly advanced and its MD5 changed** on
   every build (all 13 hashes distinct — §3).
3. **Restore writes the original bytes back** and stamps a fresh mtime — never
   `git checkout --` (which restores from the index and would silently revert
   uncommitted work) and never `Copy-Item` (which preserves the source mtime and
   leaves ninja compiling the mutant).

Anchor hygiene: every mutation is applied by a script that **refuses** unless
its anchor matches the expected number of sites (1, except M10's 3). All 12
anchors were dry-run applied and reverted before the first build, so a broken
anchor could never be mistaken for a vacuous test.

**One environment note.** `HEAD` advanced from `0d93d3c07` to `6f52a30bb`
mid-gate — lane C committed the `-implexControl` operating-point sweep to the
same branch. `git diff 0d93d3c07 6f52a30bb -- SRC/` is **empty**, so the
baseline stayed valid; the only visible effect is that `ladrunoBuild()` reports
`6f52a30bb` from M3 onward instead of `0d93d3c07`.

---

## 2. Kill table

| # | mutant (one minimal edit) | build | verdict | detectors that went RED |
|---|---|---|---|---|
| M1 | extrapolation inert — `mImplexFactor = 0` after it is computed | ok | **KILLED** | `test_negative_monotone_clock_runs_the_spec_factor`, `test_implexcontrol_refuses_past_tolerance_and_leaves_committed_state_unchanged`, `test_companion_refusal_at_commit_is_observable` |
| M2 | factor sign-blind — restore the pre-B1 `mImplexDtCommit > 0.0` gate | ok | **KILLED** | `test_negative_monotone_clock_runs_the_spec_factor` |
| M3 | total-strain form — drop `- mEpsilon_n` from the `sigma~` assembly | ok | **KILLED** | `test_implexcontrol_refuses_past_tolerance_...`, `test_implexcontrol_refusal_is_counted_and_reported`, `test_implex_sign_change_in_pseudo_dt_is_refused_...`, `test_db_roundtrip_on_a_deck_where_on_differs_from_off` |
| M4 | reduction floor dead — `mImplexDt0` never armed | ok | *survived* | — |
| M5 | `-implexControl` refusal never returns the sentinel (`return 0`) | ok | *survived* | — |
| M6 | commit-time companion refusal not counted (`noteRefusalCompanion` removed) | ok | **KILLED** | `test_companion_refusal_at_commit_is_observable` |
| M7 | un-primed-step exemption removed (`implexPrimed = true`) | ok | **KILLED** | `test_implexcontrol_refuses_past_tolerance_and_leaves_committed_state_unchanged` |
| M8 | `p_min` clamp on `sigma~` removed | ok | **KILLED** | `test_floor_clamp_fires_and_implexDetail_reports_it` |
| M9 | `avgImplexError` destructive read restored (RED-1 F6) | ok | **KILLED** | `test_implexcontrol_refusal_is_counted_and_reported` |
| M10 | re-arm after refusal removed (all 3 refusal sites) | ok | *survived* | — |
| M11 | scheme-2 no-cap refusal removed | ok | **KILLED** | `test_scheme2_without_maxsubsteps_is_refused_under_implex` |
| M12 | D2 sign-change guard removed | ok | **KILLED** | `test_implex_sign_change_in_pseudo_dt_is_refused_and_leaves_committed_state_unchanged` |

**9 / 12 killed = 0.750.** Floor **0.60** — the ADR-87 D2 `ZERO` floor
`CONTINUUM` and `SHELLMOD` are held to. (`mutation_gate.py`'s registered
`SANISAND` `ZERO` floor is the looser 0.50; the score clears both, so no floor
is being chosen after the fact.) All 12 builds succeeded; none was cut.

---

## 3. Provenance of each build

The baseline was measured on the pre-existing `afb95c40c` pyd (20 passed /
1 skipped / 1 xfailed; `SRC/` is identical between `afb95c40c` and `HEAD`, both
intervening commits being docs and measurement data). Per-mutant pyd md5 below,
all distinct, all with `LadrunoSANISAND` confirmed recompiled in the build log:

```
M1  89A8320D6A0C20EEBD81442C73562F9F      M7  909C8CABC6AC28B36D9D8449657925BA
M2  8D367028B211612C55A69B5DF6100134      M8  1A4A3DA18859516BD7706207CB479171
M3  AEEFDFEE864C1D7E0F074543F47083D6      M9  CE9338A0BBF8BF42452B474446EB88BC
M4  DDBB4E9D7B25D6716F7258921F84FC71      M10 9076BE1F017FD2CB3F9C213544692093
M5  2969DB3465C842E8FE2C9FA8DA27C094      M11 01B6A4344899D8534C70FCF5CB239177
M6  8D7CFFEA3DBD9EA17EF5CBA2AB6DD15C      M12 38B3463F7AA68B4A328BD22A6ADD044B

RESTORED (== HEAD 6f52a30bb)  366E8FCA33AD6617E326C7B37C1EEEC8
```

Post-gate: `git status` shows **no change under `SRC/`**, the pyd was rebuilt
from the restored tree, `ops.ladrunoBuild()` returns
`6f52a30bb4cae4d057f1007455b19f012dff7df9` (== `HEAD`),
`ops.ladrunoMutation()` returns `none`, and the battery is back at
**20 passed / 1 skipped / 1 xfailed**.

---

## 4. The three survivors, audited one at a time

A survivor is the actionable output of this tool, more than the score is. None
of the three is a "test that could not notice"; all three are **real coverage
gaps in behaviour the C++ deliberately implements**.

### M4 — the `-implexControl` reduction floor is never armed. **REAL GAP.**

`mImplexDt0` records `|dt|` from the first non-zero step and is the *only* thing
that stops the control refusing without limit once a driver has cut its
increment as far as it can. Killing the arming line changes nothing in this
battery because **no deck here drives a subdivision ladder**: the refusal tests
refuse on their first attempt and stop, so `mImplexDt0` is never consulted at a
second, smaller `|dt|`. The floor is exactly the B2 repair whose absence killed
the BVP `ctl` arm (`_adr92_p1_bvp_gate_results.md`), and the evidence for it
today is a BVP campaign log, not a test.

*Owed:* a test that refuses, halves `ds`, refuses again, and asserts the refusal
**stops** at `reductionLimit * |dt_0|` — i.e. that `analyze()` finally returns 0
rather than refusing forever. It also needs to assert the **magnitude**
semantics (`LoadControl(-ds)`), since a signed `mImplexDt0` was half of the
original B2 contradiction.

### M5 — the `-implexControl` refusal returns `0` instead of the sentinel. **REAL GAP, and the sharpest one.**

`test_implexcontrol_refuses_past_tolerance_...` and
`test_implexcontrol_refusal_is_counted_and_reported` both stay green with the
sentinel deleted, because M5 leaves the other two thirds of the refusal contract
intact: `mSigma = mSigma_n` (so Newton still stalls and `analyze()` still
returns non-zero) and `noteRefusalControl()` (so `implexRefusals` still counts).
The battery therefore pins the refusal's **symptoms** and not its **return-code
contract** — and the return code is precisely what `LadrunoBrick` propagates and
what a driver would use to cut the step deliberately rather than by stalling.

This is the same defect shape RED-1 F4 named at this exact site (the file's one
genuinely silent wrong answer) wearing the opposite hat: the *fix* is now
partly unpinned.

*Owed:* assert the sentinel is observable — a `LadrunoBrick` deck whose
`analyze()` failure is attributable to the propagated
`LADRUNO_MATERIAL_REFUSED` rather than to a stalled residual (e.g. one Newton
iteration allowed, or the refusal count compared against the element's own
failure path).

### M10 — a refused trial no longer re-arms the step. **REAL GAP.**

The re-arm (`mImplexStepArmed = true` at all three refusal sites) exists so a
retry measures **its own** increment rather than the first attempt's frozen `dt`
and `f`, whatever the caller does. It survives here because
`StaticAnalysis.cpp:185` calls `Domain::revertToLastCommit()` on a failed
`solveCurrentStep`, and `revertToLastCommit()` re-arms anyway — so on every deck
in this battery the removed line is genuinely redundant. It is *not* redundant
for the driver the fix was written for: one that halves its increment and
re-analyses **without** a revert.

*Owed:* a test that refuses and then re-drives the same material without a
`revertToLastCommit()`, asserting `implexDetail[5]` (`f`) tracks the **retry's**
`dt` ratio and not the first attempt's. This is the one survivor whose gap is
about a caller the battery does not have.

---

## 5. Detectors that killed no mutant

Twelve of the twenty detectors killed nothing. Audited individually, they fall
in three groups and only the last is a gap:

- **Parse-time refusals — cannot notice deleted physics, and that is correct.**
  `test_implex_refuses_unsupported_schemes[3|5|7|8|9]`,
  `test_implex_scheme1_requires_maxsubsteps`,
  `test_implex_maxsubsteps_zero_is_also_refused`. Exactly the SHELLMOD
  precedent (three parse-time refusals among its four survivors). Their scheme-2
  twin, `test_scheme2_without_maxsubsteps_is_refused_under_implex`, **does** kill
  M11 — because M11 is itself a parser mutation.
- **Paths no `-implex` mutation can reach, by construction.**
  `test_implex_unset_reproduces_recorded_sensitivity` (the `-implex` **off**
  ADR-86b regression check — its whole job is that the OFF path is untouched),
  and `test_stage0_inertness_gravity_and_hold_is_bit_identical`
  (`mElastFlag == 0`, where the operator is inert by design).
- **Green for reasons the review already named.**
  `test_implex_on_matches_off_on_a_zero_free_dof_deck` is tautological by its own
  docstring (a leak detector, not a physics gate — RED-3 F5) and survived all 12
  as predicted; `test_getcopy_after_plastic_history_shares_options_not_history`
  and `test_plane_strain_implex_smoke` are structural/smoke assertions that no
  mutant here perturbs. The PlaneStrain lane remains the coverage gap RED-3's
  matrix recorded.

Not detectors, and correctly outside the denominator:
`test_tangent_identity_frozen_ce_on_strain_dt_source` (strict `xfail`) and
`test_implex_db_roundtrip_carries_flags_and_history` (skipped on every run —
its replacement `test_db_roundtrip_on_a_deck_where_on_differs_from_off` is a
real detector and kills M3).

---

## 6. What the gate settles, and what it does not

**Settles.** RED-3 F5's charge is answered on the record: the battery kills the
inert operator (M1) and the total-strain form (M3), it pins the frozen-factor
sign rule (M2), the D2 sign-change guard (M12), the un-primed-step exemption
(M7), the `p_min` clamp (M8), the non-destructive average (M9), the commit-time
companion counter (M6) and the scheme-2 cap refusal (M11).

**Does not settle.** Three refusal-contract mechanisms are unpinned (§4) and
should be recorded as owed tests, not papered over: the reduction floor's
arming, the `-implexControl` sentinel itself, and the re-arm on a retry without a
revert. `SCALE` mode was not run: this lane's mutants are semantic and have no
scale analogue.
