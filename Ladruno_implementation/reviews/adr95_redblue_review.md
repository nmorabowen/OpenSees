# ADR-95 red/blue review — note 95 (Prandtl-Reissner quadratic root cause)

**VERDICT: MERGE-OK-WITH-EDITS.** The causal claim (corner-branch tangent defect in vanilla
`DruckerPrager.cpp`, no element touched) is solid and independently reproduced on H20 and both
tet bases (P1/P3/P4, cross-checked against source and against the CSV/log data on disk). Three
factual errors and two disclosure gaps in the note need fixing before the PR flips to ready; none
of them changes the verdict.

## Findings

| # | claim | red attack | blue answer | status |
|---|---|---|---|---|
| 1 | wall caused by corner branch, not effect of it | is the corner state a committed equilibrium or a failed-trial artifact? | `qpd_h20uri_h1.0_p1_branch.csv`: s/B 0.011125 has `n_branch3=0`, detAmin −0.064 (healthy); s/B 0.0111375 (= `q_max`/`end s/B` of the terminated leg, i.e. the last **converged** step) has `n_branch3=4`, detAmin −1.4e7…−6.6e7. `p1_h20uri.log` shows the corner read at "AT THE WALL" comes from the domain's last committed state, and the step that fails is the *next* one. Cause precedes effect. | HOLDS |
| 2 | "material only, nothing in any element changed" | check `git log --stat` for `SRC/element` on the branch | `git diff ladruno...wp/95 -- SRC/element` and `git log --stat` on the three ADR-95 commits (`56310fd39`,`038888971`,`31322a47a`) touch only `DruckerPrager.{h,cpp}` (+new test files, +python harness). The `ASDPlasticMaterial3D` diff visible in the full branch diff is from `feb358fda` (`Merge origin/ladruno`, wp/94a work), not ADR-95 authorship. | HOLDS |
| 3 | fix is correct; f1-only branch change is only the trial-norm denominator | derive corner residual/tangent; check FD test exercises γ1>0; check f1-only delta is denominator-only | `git diff 31322a47a~1 31322a47a -w` (real diff ≈141 lines, see note on diff size below) shows exactly: index-driven `if(Jact(0)==1)/if(Jact(1)==1)` replacing the dead `else if (Jact(i)==2)`, and `norm_eta_trial` latched before the return, guarded `>1e-13`, used only in the `devSoft` term. For `Jact=(1,0)` the old loop's `i=1` pass was already a no-op, so the ONLY delta on the cone is the norm swap — confirmed by inspection, matches the note's claim. `test_hardened_corner_is_non_degenerate` asserts `gamma0>1e-9 AND gamma1>1e-9 AND ||eta||>1e-2` before `test_corner_tangent_matches_finite_difference[k=0..5]` FDs that same non-degenerate corner state — the FD genuinely exercises γ1>0, not the degenerate apex. | HOLDS |
| 4 | fix doesn't silently move gate #722 beyond its bands | check the numbers in §4 against the actual gate bands | Bands are `_MEASURED={1.0:1.0842,...}` × (1±0.03) in `test_r3_prandtl_collapse_gate.py` → h0=1.0 band is exactly 1.051674–1.116726 as quoted; fastest-leg pre/post (1.0849417→1.0849561) and linear control (1.7e-6 rel) both sit inside it. **But** only the h0=1.0 non-associated leg was re-run post-fix — `_adr95_p0_results.md` §6 and `_adr95_p4_results.md` §7/§11 both say the other two resolutions and the associated-flow control ("owed before the PR flips to ready") were never run against the fixed material, and note 95 does not repeat that caveat anywhere. | EDIT |
| 5 | tet transfer (P3) and fixed-material tet plateaus are the same event / real capacities | check termination modes, tail slopes, corner-GP identity | `_adr95_p3_results.md`: pre-fix `tet10`/`beziertet10` FLOOR at the SAME element/GP (ele 534, gp 3) as H20's corner; `p4_tet10_long.log`/`p4_bezstd_long.log`/`p4_bezbbar_long.log` all show `MODE=TARGET`, `nfail=0`, `nsub=0/200`, tails ≤0.02% — genuine flat plateaus, not budget artifacts. Corner-GP counts (181/170/0 at the final station) and `detAmin_min` back to the −0.06…−0.08 ambient floor confirm the same defect signature, now healed. | HOLDS |
| 6 | P2 (SY knob) is a legitimate confirmatory/falsifying experiment for H1 | is "reach moves monotonically with T" independent evidence, or the same mechanism restated? | `_adr95_p4_results.md` §9, verbatim: "P2's SY ladder is **superseded** as a diagnostic: raising SY moved the wall because it moved T out of reach of the heave zone, which is the same mechanism, not an independent knob." Note 95 §3 reports the P2 result as a clean confirmation without carrying this demotion forward. | EDIT |
| 7 | H4 (element rank) correctly demoted to secondary | does the rank data actually support "not the cause"? | `_adr95_h4_rank_results.md`: H4 REFUTED for LadrunoBrick20-std (+4 over predicted) and LadrunoBrick-bbar (+1); the wall's own signature (flat σ_min, abrupt collapse, removed by material fix alone) is inconsistent with a rank story. Note §5 states this accurately. | HOLDS |
| 8 | no walled/allowance number is quoted as a capacity, and quoted numbers are accurate | scan §0 table against the source JSON/logs | Two defects found. (a) §0 quotes tail = 0.11% for BOTH the `tet10` and `beziertet10 std` long legs and 0.02% for `beziertet10bbar`; the actual `tail_pct` in `tpd_tet10_h1.0_p4long.json` / `tpd_beziertet10_h1.0_p4long.json` / `tpd_beziertet10bbar_h1.0_p4long.json` is 0.0180% / 0.0180% / 0.0035% — off by ~6x, ~6x, ~5.7x. Still flat enough that CAPACITY=yes is not in question, but the printed numbers are wrong. (b) The pre-fix `BezierTet10 -bbar` cell reads "TARGET s/B 0.02, 0.713 at matched s/B 0.008" — this conflates two different ratios into one number. `tpd_beziertet10bbar_p3_run.log`/P3 §2a: ratio AT s/B=0.02 is **0.9726** (still hardening); the **0.7128** figure is the separate matched-s/B=0.008 reading from the P3 comparison table. As written the cell misstates the s/B=0.02 capacity-adjacent number by 36%. No FLOOR/WALL number anywhere is misreported as a capacity — that rule holds; these are two accuracy bugs inside otherwise-correctly-labeled cells. | EDIT |
| 9 | note 95 doesn't contradict notes 82/83 without saying so | cross-check the abrupt-tangent numbers and DR framing | Note 95 §2's σ_min/cond numbers (2.16e-4→2.29e-7, s/B 0.01115→0.01118) match note 82 §7.3 exactly; the "DR walked through because no tangent formed" line matches note 83 §5.2/§7 item 7. No contradiction found. | HOLDS |
| 10 | ledger/banner obligations met | check `LEDGER_vanilla_files.md`, `LEDGER_quirks.md`, banner | Both ledgers carry full ADR-95 rows (vanilla file: 5 rows for P0+P4+response; quirks: a ~2000-word entry dated 2026-09-07 covering the two-surface model, the dead arm, the tangent bug, and the fix). No `LEDGER_implementations.md` row and no banner line — correct, since P4 shipped a bug fix to a vanilla file, not a new Ladruno feature/option (plan §5 only required a ledger row "if P4 ships an option"). | HOLDS |

## Required edits to the note

1. **§0 table, tail percentages** — replace `tail 0.11 %` (tet10 row) with `tail 0.018 %`;
   replace `tail 0.11 %` (BezierTet10 std row) with `tail 0.018 %`; replace `tail 0.02 %`
   (BezierTet10 -bbar row) with `tail 0.0035 %`. Source: `tpd_{tet10,beziertet10,beziertet10bbar}_h1.0_p4long.json` → `tail_pct`.
2. **§0 table, BezierTet10 -bbar pre-fix cell** — replace
   `TARGET s/B 0.02, 0.713 at matched s/B 0.008 (no corner GP in range)` with
   `TARGET s/B 0.02, 0.9726 of exact (still hardening, tail 6.1 %, not a plateau); at matched
   s/B 0.008, ratio 0.7128 (no corner GP in range)`.
3. **§0 table, corner-GP-count phrasing** — the quoted "40–181" / "120–203" ranges are the
   spread over the last ~6 sampled stations near the plateau, not the full-run min/max (true
   min is 4 at both legs' first station). Reword to `"corner GP count fluctuates 40–181 (resp.
   120–203) over the plateau's sampled stations"` so a reader doesn't read it as a monotone count.
4. **§3 (the knob)**, append: "P4 §9 supersedes this as an independent diagnostic: raising SY
   moves T, the same quantity the corner defect triggers on — the monotone reach is consistent
   with H1 but is not a second, independent confirmation of it. H1's confirmation rests on
   P1/P3/P4's branch/state evidence, not on P2."
5. **§4 or §7**, add a line: "Only the fastest gate leg (h0 = 1.0, non-associated) has been
   re-run against the fixed material. The other two resolutions and the associated-flow control
   (`test_r3_prandtl_collapse_gate.py`, ~61 min total) are still owed before the PR flips to
   ready (P0 §6, P4 §7/§11)."

## Minor, not required
`git show 31322a47a --stat` reports ~2380 changed lines in `DruckerPrager.cpp`; this is a CRLF
line-ending flip on a file the parent commit stored as LF (`git diff -w` shows the real change is
130 insertions/11 deletions ≈ 141 lines, matching "~200 lines"). Worth one line in the PR
description so a reviewer isn't alarmed by the diff size; not a note defect.
