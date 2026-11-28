# WP-260 Review Protocol — Gemini Authors, Claude Reviews

Every WP-260 PR gets a Claude review before ready-flip. This file is what the author prepares and what the reviewer checks, so neither side has to guess. It does not replace ADR-87 D4 (owner-run out-of-family review with a committed verdict, ADR §7 / D7).

## 0. Evidence rules (both sides)

The reviews of #1158 and #1160 found claims that went further than their evidence: test names that did not exist, "0-ULP oracle" tests that compared the new code with itself, "all green" on a run whose conclusion was `failure`, a gate reported as passed that the script did not implement, and an "owner-approved" pick with no owner comment. The code was mostly right; the reports were not. These rules make a report checkable without trusting its author.

1. **A claim needs an artifact.** Every claim in a review request, handoff or board row points to something that exists at the stated SHA: a test name in the branch, a committed file (`results/<sha9>/…`, an MD5 list written by a script in the repo), or a CI run URL. A claim with no artifact counts as not made. A claim whose artifact is missing or says otherwise is an `important` finding, even when the code is correct.
2. **Name what an oracle compares against.** A `G-BYTE` proof names both builds (base SHA and head SHA) or both paths (switch on and off). A `G-EQUIV` proof names an independent reference: an analytic solution, a numpy oracle, or a dense path behind a test-only toggle. The new code compared with itself is a regression test, not an oracle. Label it as one.
3. **Copy status from the tool, never summarize it.** CI state is what `gh pr checks <pr>` / `gh run view <id> --json jobs` print, job by job. When the run conclusion is `failure` because of a known non-required job, say exactly that. "Passed" means the check exists and ran. A gate the script does not implement is listed under **Known gaps**.
4. **Approvals and decisions carry a link.** "Owner-approved", "owner decision" or "amendment accepted" counts only with a link to the owner's comment. Without one, write "proposed, awaiting owner".
5. **Measure on a quiet machine.** No timing run while a build, a test suite or another benchmark runs on the same host. Record a load probe per run. A record whose deep-profiler overhead is below 0.95×, or whose scope call counts differ between modes, is void.
6. **One package per PR, conservative defaults.** A PR implements one menu package. A change that alters results (`G-EQUIV`, `G-OPT`) ships default-off unless a committed independent oracle (rule 2) backs it and the owner accepts the default on the `program` issue. Files outside the brief's **Owned Files**, and any change to a brief target, need an amendment the owner accepts on the `program` issue **before** the code is written.
7. **Ideas go to issues, not code.** A new lever found while working (e.g. the `enforceSPs` double update) becomes a GitHub issue: the hypothesis, the phase and expected `p` / `s`, the cheapest measurement that could prove it wrong, and the evidence so far. No code until a measurement supports it (plan §0.5) and the owner picks it. Ideas never live in untracked folders.

## 1. Author: the review request

Post on the draft PR only when the warrant package is green (fork-checklist §5 items 1–6):

```markdown
### Review request — 260<x> — <head sha>
- **What changed:** <one line per change, with file::symbol>
- **Gate per change:** <G-BYTE | G-EQUIV | G-OPT> + the proof, in the evidence table below
- **Contracts touched:** <#1a/#1b/#2/#3/#4/#5/#6 or none>
- **Measured:** the before/after table (§3 D below), D1 decks, thread count, SHA of both builds, load probe
- **Selection:** the plan §0.5 numbers this WP was picked on (p, s, projected G) and the link to the owner's pick
- **Scope:** files outside **Owned Files** and brief deviations, each with the link to the accepted amendment
- **CI:** `gh pr checks` output for <head sha>, job by job (rule 3)
- **Known gaps:** <everything not done, not proved, or proved only by a self-comparison>

| # | Claim | Gate | Artifact (test name / committed file / run URL) | Compared against | Command to reproduce |
|---|---|---|---|---|---|
```

Every claim in the request, the PR body and the handoff has a row. Rows are what the reviewer checks; prose without a row is not a claim (rule 1).

## 2. Reviewer: what Claude does

1. Reviews the PR **at the requested SHA** (a later push needs a new request).
2. Checks every row of the evidence table: the test exists and runs what its row says, the file is committed, the run belongs to the SHA. A row that fails is an `important` finding (§0 rule 1).
3. **Reproduces one headline claim independently**: rebuilds the head SHA in its own worktree and re-runs one MD5 pair or one interleaved timing. If it does not reproduce, the review stops there with `changes needed`.
4. Checks every item of §3, reading the code paths, not only the diff (callers of a changed function, subclasses of a changed base).
5. Posts **one** PR review comment: numbered findings (severity `important` / `nit`, file:line, what goes wrong, fix) and a verdict:
   - `ready` — no open `important` findings;
   - `changes needed` — at least one `important` finding;
   - `stop — ADR amendment needed` — the WP conflicts with a contract, a prior ADR closure, or the selection rule.
6. Re-checks after the author answers and posts an updated verdict.
7. Updates the evidence scorecard: one running comment on the `program` issue (#1156) with, per review, the claims made, claims verified, claims reproduced, and claims found false (with the finding number). The owner reads it to decide how much of a WP to delegate.

## 3. Review checklist

### A. Scope
- Only the brief's **owned files** are edited; none of the files open PRs #1147/#1153/#1154 edit (ADR §2.2) while they are open.
- The WP was **picked** under plan §0.5 (or is a declared ride-along), and its **Blocked by** is satisfied.
- Nothing re-opens a lane a prior ADR closed (ADR §2.3) without an owner decision.

### B. Correctness gates (ADR-68)
- `G-BYTE`: before/after recorder outputs bit-identical (MD5 + signbit) on the D1 decks **and** the relevant Zone-A tests, same thread count; for threaded code also at ≥ 2 thread counts.
- `G-EQUIV`: an **independent** oracle (not the new code path compared with itself), a fail-closed/refusal test, and Linux TSan for threaded code.
- `G-OPT`: off by default, a refusal when combined with an option it can't honour (notably `-deterministic`), and a documented accuracy envelope.
- Anything inside the WP-132 `-deterministic` block re-proves CNR bit identity run-to-run.

### C. Contracts (ADR §4)
- Signatures match the ADR exactly; defaults leave un-upgraded classes unchanged.
- #1a: a test that a tangent requested after a stress-only trial equals the full-trial tangent bit for bit, on **every** channel the element uses (`getTangent`, `getSpatialTangentTensor`).
- #1b: no second extrapolation on a material with its own `-implex`.
- #2: `TRBDF2`/`TRBDF3` keys differ between legs; reuse keyed on `(c1, c2, c3)`.
- #3: a test per row of the commit-side-effect table (porous overlay advances once per step, contact never in augment mode, recorders silent).
- #4: no subclass bypasses the NVI `setSize`; a plan from an old SOE never matches a new one.

### D. Performance evidence
- Same machine, same build flags, same thread count, before and after; both `ladrunoBuild` SHAs stated.
- At least 3 interleaved runs per side; median and spread reported.
- The PR reports the **phase** speedup `s` and the **whole-run** gain on each D1 deck, and checks the whole-run gain against `1 / ((1 − p) + p / s)`.
- Phase reading follows the profiler quirks: `elem.tangent` includes the `addA` scatter; deep profiles are 1-thread only.
- No claim exceeds what was measured; a measured `G` below half the projection gets a "why" line.

### E. Fork rules
- `// Ladruno ...` marks on every upstream edit and a `vanilla/` fragment row per file; an `implementations/` fragment with `class_tags` for new classes; a `quirks/` fragment for anything learned.
- Registrations complete (fork-checklist §4): classTags band, manifest, command table, broker, `stamp_headers` GLOBS, banner line when shipped.
- MKL edits: the `MKL 2024.2 compat (WP-148)` job is green; no symbol newer than oneMKL 2024.2.2.
- All 5 targets built from one SHA (`opensees.pyd` not stale); Zone-A green.

### F. Threading
- Each `ladrunoThreadSafe*()` opt-in is justified for the whole call graph (element + materials + sections); no function-scope statics or shared scratch left on the threaded path.
- No `opserr` or Python callbacks from worker threads.
- The thread count follows the fork's existing thread-count policy; no new OpenMP runtime.

## 4. Program-level review

After `260-w0` and after each merge, Claude reviews the phase-share table and the re-ranked menu (plan §0.5). Findings go on the `program` issue (#1156); script nits go on the `260-w0` PR. The checks follow the `260-w0` brief:
- **Decks (§1):** each deck is the D1 configuration (or an owner amendment recorded before the run), at production solver settings (factor reuse, `Linear -factorOnce` on a linear deck), and the implicit decks are run at two sizes, with one ≥ 10⁵ DOF.
- **Runs (§2):** one SHA across all records; `1t_deep`, `1t_coarse` and `nt_coarse` all present; no non-default `OMP_*`/`KMP_*` override; the "threaded?" column agrees with the `ladrunoThreads` announcement; step 1 dropped; spread ≤ 5 % or 5 reps.
- **Scopes (§3):** `elem_by_type` and `soe.addA` disjoint; residual within ±0.5 %; no unexplained remainder above 5 %.
- **Ranking (§4):** `p` limited to the phases each package changes; `s` per rule 2 with a source; Wave-3, unranked and not-phase-rankable packages kept out of the pick.
- **Pick (§5):** does it follow from the numbers?
