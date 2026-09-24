# WP-115 — Agents surface pilot: task guides + quirk-pattern lints

Revision 1. Not yet adversarially reviewed.

Status: **built; draft PR #850. Merge after #841** — L2 stays red until #841 lands, because
SANISAND's process-wide IMPL-EX globals are the live instance it exists to catch (2026-09-23).

Scoped 2026-09-23. Branch `wp/115-agents-surface`, cut from `ladruno` @ `79e062367`. Source: a
read of how `basecamp/omarchy` organizes its repo for agents (`AGENTS.md` + short task-triggered
guides in `agents/skills/` + a metadata lint, `omarchy commands --check`, that fails CI).

## Problem

The fork rarely loses a lesson. It fails to put the lesson in front of the agent at the moment of
action. `LEDGER_quirks.md` is 7,619 lines with 470 headings, and it records its own repeats (32
recurrence/rediscovery notes as of 2026-09-23):

- **Rayleigh P-clobber (#562).** The quirk was already in the ledger, and the ADR-70 P4a author
  reasoned about it in a code comment. It still recurred across the whole plane family
  (Quad/CST/LST finite + CSTPair). The ledger's own lesson: *"When a quirk names a pattern, grep
  for the pattern, don't reason about the instance."*
- **State that survives `wipe()`**, at least five times: the global `getCommitTag()` counter
  (quirks, "`getCommitTag()` is a GLOBAL monotonic counter"), `EQ_Constraints` leaking across
  models ("`wipe` / `Domain::clearAll()` did NOT clear EQ_Constraints"), the Profiler singleton
  ("The Profiler is a process-global singleton"), the ADR-69 energy registry mixing state between
  models, and open PR #841 (SANISAND's process-wide IMPL-EX diagnostic ledger).

A 7,600-line ledger is a good archive and a poor warning system.

## Shape

1. **`AGENTS.md`** holds the working rules. `CLAUDE.md` becomes one line, `@AGENTS.md`, so Codex
   and other agents read the same rules Claude does. The content moves verbatim; the only
   addition is a table routing to the task guides.
   *Accept:* every rule heading of the old `CLAUDE.md` is still present.
2. **Two task guides** in `.claude/skills/<name>/SKILL.md`: `ladruno-new-element` and
   `ladruno-new-material`. Claude Code loads them by description; the `AGENTS.md` table lists
   their paths for every other agent. One copy, no sync script. Each is a checklist of at most
   100 lines whose items link into `LEDGER_quirks.md` entries rather than copying them.
   `.gitignore` ignored all of `.claude/`; it now re-includes `.claude/skills/` only.
3. **`ci/check_quirk_patterns.py`**, a static gate in `ladruno.yml` (runs on drafts). It scans
   only fork-authored files (those carrying the `LADRUNO-HEADER-START` stamp).
   - **L1, Rayleigh (fails).** Flags adding `getRayleighDampingForces()` into anything that is
     not a function-local vector. A same-line or preceding-line comment
     `// ladruno-lint: rayleigh-ok <reason>` waives one site; the reason is mandatory.
   - **L2, wipe-state (fails).** Every process-wide singleton/registry in a fork-authored file
     (an `instance()` accessor) must either be reset from `Domain::clearAll()` or carry
     `// ladruno-lint: wipe-ok <reason>` at its declaration.
   - *Acceptance, as a mutation gate:* each lint must catch the incident that motivated it. L1
     run on `4a975edee` (parent of the #562 fix) must flag the four plane elements. L2 run on
     `ladruno` before #841 must flag SANISAND. A lint that misses its own incident does not ship.
4. **Ledgers + measurement.** `LEDGER_implementations` row (kind: tooling); pointer from
   `WORKFLOW_GOTCHAS`; the baseline above.
   *Review:* after about 10 work packages, count adversarial/post-merge findings that match an
   entry already in the ledger. If that count does not drop, stop investing in guides.

## Rejected approaches

- **Lint every mutable `static`.** 123 candidates in 29 fork files, mostly one-shot warning
  latches and scratch arrays. The noise would get ignored. L2 targets singletons only.
- **Timestamped ADR filenames (Omarchy's `migrations/<unix-ts>.sh`).** ADR numbers are
  load-bearing in branch names, commits and memory, and only two prefixes (`78_`, `19_`) hold
  more than one ADR. The other shared prefixes are deliberate companion-doc families.
- **Split or rewrite `LEDGER_quirks.md`.** Out of scope; the guides point into it.
- **Run the lints on vanilla code.** `SRC` has ~240 upstream Rayleigh call sites. The
  vanilla-footprint rule says upstream bugs are raised with the owner first.
- **A command index/router for `Ladruno_scripts/`.** Deferred; worth less than the quirk lints.
- **Pilot in apeGmsh.** No recent activity in its checkout and no recurrence baseline.
- **Guides in Omarchy's `agents/skills/`.** Claude Code would not load them automatically;
  `.claude/skills/` plus the `AGENTS.md` table reaches both Claude and other agents.

## Results (2026-09-23)

- **Step 1.** `AGENTS.md` holds the rules; every heading of the old `CLAUDE.md` is present.
  `upstream_pr_campaign.md`'s never-ships list gains `AGENTS.md` and the lint.
- **Step 2.** Guides: `ladruno-new-element` (78 lines), `ladruno-new-material` (79 lines); all 46
  ledger pointers resolve (enforced by L3). The #588 degeneracy-guard lesson was only in code
  comments and a commit message, so it gained a `LEDGER_quirks` entry for the guide to point at.
- **Step 3.** `ci/check_quirk_patterns.py` + a 10-case self-test, last in the `static-gates` job
  (the required job name is unchanged). On the branch it reports one finding: LadrunoSANISAND.
- **Mutation acceptance, all passed** (lint run on `git archive` trees):

  | Run | Tree | Expected | Got |
  |---|---|---|---|
  | A1 | `4a975edee` (parent of #562) | L1 flags Quad, CST, LST, CSTPair | 8 plane sites flagged (2 per element) |
  | A1b | `f89687274` (the #562 fix) | plane family clean | clean |
  | A2 | `ladruno` @ `79e062367` | L2 flags SANISAND | flagged |
  | A2b | #841 head `e900f49e0` | SANISAND passes | passes |

- **Waivers (6, comment-only, each checked by reading the code).** `BezierTet10`, `BezierTri6`,
  `LadrunoIMKBeam`, `LadrunoIMKBeam2d`: the shared buffer is written only by
  `getResistingForce`/`getResistingForceIncInertia`, which nothing on the Rayleigh path
  (`getMass`/`getTangentStiff`/`getInitialStiff`) calls. Safe today, but it depends on that
  staying true; converting them to the snapshot idiom would remove the dependency (C++ change, not
  in this WP). `MassScalingEnergyRegistry`: owner-scoped, cleared in each publisher's destructor,
  and `wipe` deletes the integrator. `Profiler`: survives `wipe` by current design; whether it
  should reset is left open.
- **Step 4.** `LEDGER_implementations` row; the `LEDGER_quirks` conventions now say to enforce
  greppable quirks in the gate and point to them from the guides. That replaces the planned
  `WORKFLOW_GOTCHAS` pointer: the convention belongs where quirks are written.

## Step 5 — snapshot conversion of the four waived sites (2026-09-23)

The owner asked for the four waived Rayleigh sites to be converted, after an adversarial review.

**Review (two independent Opus reviewers, read-only).**
- *Safety/conversion:* the "safe today" claim HOLDS at all four sites (every path from `Element::getRayleighDampingForces()` traced). The conversion is bit-identical if operation order is kept: Bezier is ((f−Q)+M·a)+R, IMK is ((f−Q)+R)+m·a. Copying the CST template onto IMK would have changed the last bit.
- *Concurrency:* no race today or after. Only `Element::update()` runs in the WP-107 OpenMP loop (`Domain.cpp` ~2693); assembly is serial and none of the four classes is allowlisted.
- *Coverage:* none of the four had a transient Rayleigh test, and the existing plane gate (tiny βK) cannot see a dropped Rayleigh force.
- *Lint:* one CRITICAL hole (a snapshot taken after `const Vector &v = getRayleighDampingForces()` passed) and several MAJOR ones (plain/copy binds, pointer and `this->` targets, one-line `if` bodies, unqualified L2 hook matching, any `X::clearAll` counted as a hook).

**Done.**
- Conversion in `BezierTet10`, `BezierTri6`, `LadrunoIMKBeam`, `LadrunoIMKBeam2d`: function-local `static Vector res` seeded before the Rayleigh call, `return res`. Waivers removed.
- `tests/test_rayleigh_inertia_bezier_imk.py` (zone_a, 34 cases). IMK is differential against `elasticBeamColumn` carrying the same lumped masses as nodal masses (agreement 1.2e-14); Bezier uses the self-validating overshoot rig (1.93–1.97) plus a 5% damping leg (ratios 0.914–0.931, band 0.80–0.97).
- `ci/check_quirk_patterns.py` rewritten as a small C++ scanner; every review hole has a self-test case (33 cases). Historical acceptance re-run and unchanged.

**Evidence** (pyd rebuilt per row; sources restored and a full 5-target rebuild after):

| Build | Expected | Result |
|---|---|---|
| original (pre-conversion) code | bit-identical to converted | 4,382 / 4,382 recorded displacements identical over 36 runs |
| original + re-entry hazard in `getTangentStiff` | βK legs fail | 12 failed — exactly the βK legs (the only path through `getTangentStiff`) |
| converted + same hazard | all pass | 34 passed |
| converted, Rayleigh add dropped | every leg with element damping fails | 26 failed; the 8 passes have nothing to detect (static, tiny-βK by design, IMK αM with nodal mass) |
| converted, inertia add dropped | every leg with element mass fails | 20 failed; the 14 passes have no element mass |
| final committed code, full rebuild | all pass, bit-identical | 148 passed (new + Bezier + IMK + plane-dynamics + response tokens); bit-identical to converted |

**Found along the way — owner decision needed (not fixed; each changes results):**
1. **CRITICAL — Bezier ground-motion inertia has the wrong sign.** `BezierTet10::addInertiaLoadToUnbalance` / `BezierTri6` build `Q += +M·a_g`; the vanilla and LadrunoBrick convention is `−M·a_g`. Proven by running it: under a constant +2.0 ground acceleration, a rigid-body probe gives relative acceleration −2.000 for vanilla `quad`/`stdBrick` and +2.000 for both Bezier elements. Present since BezierTet10 was added (2026-05-30); no test ran a Bezier element under `UniformExcitation`.
2. **MAJOR — `betaKc` damping frozen on the IMK beams.** `LadrunoIMKBeam(2d)::commitState` never calls `Element::commitState()`, so `Kc` keeps the tangent from when `rayleigh` ran — zero if that was before the first step.
3. **Upstream — vanilla `ElasticBeam2d` subtracts the ground-motion Q twice** with element `-mass` (response exactly 2× the same beam with nodal masses). `ElasticBeam3d` is correct. Recorded in `LEDGER_quirks`; not fixed (vanilla-footprint rule).

## Open questions

- Whether L2's site list stays small enough to hand-classify as fork singletons grow.
- Should `wipe` reset the Profiler? (Waived as current design.)
- ~~Convert the four waived Rayleigh sites?~~ Done (Step 5).
- Fix the three findings in Step 5 (Bezier ground-motion sign, IMK `betaKc`, upstream `ElasticBeam2d`)?
