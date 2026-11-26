# WP-260 Review Protocol — Gemini Authors, Claude Reviews

Every WP-260 PR gets a Claude review before ready-flip. This file is what the author prepares and what the reviewer checks, so neither side has to guess. It does not replace ADR-87 D4 (owner-run out-of-family review with a committed verdict, ADR §7 / D7).

## 1. Author: the review request

Post on the draft PR only when the warrant package is green (fork-checklist §5 items 1–6):

```markdown
### Review request — 260<x> — <head sha>
- **What changed:** <one line per change, with file::symbol>
- **Gate per change:** <G-BYTE | G-EQUIV | G-OPT> + the proof (test name, MD5 list, oracle)
- **Contracts touched:** <#1a/#1b/#2/#3/#4/#5/#6 or none>
- **Measured:** the before/after table (§3 D below), D1 decks, thread count, SHA of both builds
- **Selection:** the plan §0.5 numbers this WP was picked on (p, s, projected G)
- **Known gaps:** <what is not done or not proved>
```

## 2. Reviewer: what Claude does

1. Reviews the PR **at the requested SHA** (a later push needs a new request).
2. Checks every item of §3, reading the code paths, not only the diff (callers of a changed function, subclasses of a changed base).
3. Posts **one** PR review comment: numbered findings (severity `important` / `nit`, file:line, what goes wrong, fix) and a verdict:
   - `ready` — no open `important` findings;
   - `changes needed` — at least one `important` finding;
   - `stop — ADR amendment needed` — the WP conflicts with a contract, a prior ADR closure, or the selection rule.
4. Re-checks after the author answers and posts an updated verdict.

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

After `260-w0` and after each merge, Claude reviews the phase-share table and the re-ranked menu (plan §0.5): are the scopes read correctly, is `s` conservative, does the proposed pick follow from the numbers? Findings go on the `program` issue (#1156), not on a PR.
