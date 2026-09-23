# WP-115 — Agents surface pilot: task guides + quirk-pattern lints

Revision 1. Not yet adversarially reviewed.

Status: **in progress on the branch; draft PR open** (2026-09-23).

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

## Open questions

- Whether L2's site list stays small enough to hand-classify as fork singletons grow.
