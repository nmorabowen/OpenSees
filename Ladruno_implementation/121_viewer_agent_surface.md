# WP-121 — Viewer tools join the agents surface: a task guide + a ledger gate

Revision 1. Not yet adversarially reviewed.

Status: **built; draft PR #856.** No merge-order constraint: V1 reads a PR's own diff, so the
stale ledger rows listed under "Live incidents" do not turn it red.

Scoped 2026-09-25. Branch `wp/121-viewer-agent-surface`, cut from `ladruno` @ `bc5c33453`. Extends
WP-115 ([[115_agents_surface_pilot]]) from elements and materials to the viewer tools only:
`Ladruno_tools/profiler_viewer/` (`ProfilerResults`, FastAPI `profiler_api`, React/TypeScript/Vite
`frontend/`, `launch.py`, `profiler_monitor`, the `*_smoke_*` pairs) and `Ladruno_tools/monitor_viewer/`
(`MonitorReader`, `monitor_server` + `monitor_page.html`, `monitor_view`). The rest of the fork is out
of scope. Method: the user-level `agent-surface` playbook, phases 0, 2 and 3 (phase 1, `AGENTS.md`,
already exists).

## Problem

The viewers have 14 PRs (#35 to #487, 2026-05-30 to 2026-07-05) and no activity since. Their lessons
live mostly outside the repo: two `LEDGER_quirks` entries, commit messages (#484, #487), plan-doc
logs (`06_profiler.md`, `08_analysis_monitor.md`), and session memory
(`project_profiler.md`, `project_analysis_monitor.md`). Three lessons were written down and still
bit again:

| # | Lesson (where it was written) | Recurred | Strength |
|---|---|---|---|
| R1 | "prior PRs ... missed LEDGER_implementations' row; check that row explicitly" — session memory, 2026-05-31, after PR #56 had to backfill #53/#55. The standing `CLAUDE.md` ledger rule said the same. | 2026-07-05: #485 added `profiler_monitor.py`, #487 added all of `monitor_viewer/`; neither touched any ledger. | strong |
| R2 | Verification on `make_sample.py` data only: P8 (#51, 2026-05-31) was checked panel by panel on the sample; the Series tab was empty on real runs until #52. Recorded as a symptom, not as a rule. | 2026-07-04, #484: flame graph empty and every `share` ~1e8 on 100% of engine-written files, all tests green. | soft |
| R3 | "preview_screenshot kept TIMING OUT ... assert via DOM queries, not screenshots" — session memory, 2026-05-31 (#51). | 2026-07-04 (#484): "preview_screenshot still flaky". | soft |

Baseline for the WP-115 review: **3 viewer recurrences** (1 strong, 2 soft) as of 2026-09-25.

Single occurrences, not recurrences: #487's three self-review bugs (500 instead of 503 on a
missing sink, a Follow button that did nothing, `--watch` giving up before the writer started) and
the Vite template's TypeScript strictness traps (`tsc` reports them).

## Shape

1. **One task guide**, `.claude/skills/ladruno-viewer-tools/SKILL.md` (89 lines), built on the
   playbook's "a UI needs its own visual-verification guide": the contracts other code reads (the C++
   writer, `api.ts`, the monitor sink, apeGmsh's `profiler.py`), the data to check against, the
   checks to run (none of them run in CI), reference-vs-candidate verification driven through the
   DOM, stopping what you launched by PID, and the ledger items.
   *Accept:* 89 lines ≤ 100; every `Quirks: "..."` pointer resolves (L3 of
   `ci/check_quirk_patterns.py` scans all `.claude/skills/*/SKILL.md`, so it now covers this guide).
2. **One row in the `AGENTS.md` task-guide table.** Nothing else in `AGENTS.md` changes.
   *Accept:* `git diff AGENTS.md` is one added line.
3. **Lessons into the repo.** R2 and R3 lived only in memory and a commit message. They become two
   `LEDGER_quirks` entries (appended at the end, per the `merge=union` convention) for the guide to
   point at. R1 is enforced instead (step 4); its evidence is this document.
4. **V1, `ci/check_viewer_ledger.py`:** a change that ADDS a file under `Ladruno_tools/` must also
   modify `LEDGER_implementations.md`. It reads the PR diff three-dot (`base...head`, the change since
   the merge base), so a multi-commit PR whose ledger edit landed in a later commit passes.
   Modification-only changes are not checked; renames are not additions, copies are. Two steps
   last in the `static-gates` job: the self-test (every event) and the gate (`pull_request` only,
   with the event's base and head SHAs). The required job name is unchanged.
   *Accept (mutation gate):* V1 flags the incident PRs and passes their fix and every compliant PR
   (Results).
5. **Ledgers:** a `LEDGER_implementations` row (test infrastructure), a `ci/README.md` row, and
   `ci/check_viewer_ledger.py` on the never-ships list in `upstream_pr_campaign.md`. No vanilla file
   is touched; `.gitignore` already re-includes `.claude/skills/`.

## Rejected approaches

- **A code-pattern rule for #484** (`x or 1` denominators, the icicle's `wall_ms <= 0` early return).
  The fix `acb932935` keeps both idioms (`profiler_results.py` `eff_root_ns or 1`, `wall_ns or 1`;
  `Icicle.tsx` unchanged); it fixed the data, not the idiom. A rule that flags the fix does not ship.
- **Rules for #487's three bugs.** They were fixed before the first commit, so there is no pre-fix
  commit to accept against.
- **A TypeScript lint for `verbatimModuleSyntax` / `erasableSyntaxOnly`.** `tsc -b` in
  `npm run build` already fails on them.
- **A tree-state ledger check** ("every viewer file is named in the ledger"). The profiler row lists
  `Ladruno_tools/profiler_viewer/*`, so it passes today with the row stale. It cannot see the incident.
- **Checking each commit instead of the PR.** It flags #58's launcher commit (`2da2d1d08`) although
  the same PR carried the ledger note (`6b87fa6b9`). The PR is the unit AGENTS.md sets.
- **Scoping V1 to the two existing directories.** #487 created `monitor_viewer/`, a new sibling. A
  rule scoped to `profiler_viewer/` before it would have missed its own incident. `Ladruno_tools/`
  holds only the viewers today, so the wider scope adds no noise.
- **Running V1 on push, schedule or dispatch.** There is no PR base; what lands on `ladruno` already
  passed as a PR.
- **Carrying #487's "h5py raises a bare `OSError`, not `FileNotFoundError`" into the ledger.** It
  does not reproduce: h5py 3.16.0 raises `FileNotFoundError` (an `OSError` subclass) for a missing
  file, with and without `swmr=True`. The comment at `monitor_reader.py:41` is wrong on this
  version; the `os.path.exists` guard it motivated is harmless. Older h5py not tested.
- **A second guide for the monitor viewer.** Both viewers share the stack (h5py, FastAPI, a browser
  page) and the verification loop.
- **Pixel-diff screenshot baselines.** No browser in CI, and screenshots were the flaky part (R3).
- **Backfilling the stale profiler/monitor rows here.** Out of scope for this WP; listed below.

## Results (2026-09-25)

**V1 survey — every PR that touched `Ladruno_tools/`** (squash PRs: `c^...c`; merge PRs:
`M^1...M^2`):

| PR | Files added under `Ladruno_tools/` | Ledger edited | V1 |
|---|---|---|---|
| #35 | 6 (`profiler_viewer/` created) | no (row first added by #46) | **flagged** |
| #46, #49, #50, #51, #52 | 2, 2, 3, 24, 2 | yes | pass |
| #47, #48, #55, #484 | 0 | — | pass (nothing added) |
| #53 | 2 (`header_smoke_*`) | no (backfilled by #56) | **flagged** |
| #58 | 3 (launcher) | yes, in a later commit of the same PR | pass |
| #485 | 2 (`profiler_monitor.py`, `monitor_smoke.py`) | no | **flagged** |
| #487 | 6 (`monitor_viewer/` created) | no | **flagged** |

Four flags, each a real miss of the `AGENTS.md` rule; no false positive.

**Mutation acceptance** (the checker run on real git ranges):

| Run | Range | Expected | Got |
|---|---|---|---|
| A1 | #487 `0ec96549a^...0ec96549a` | flag | 6 findings, exit 1 |
| A2 | #485 `e7e403763^...e7e403763` | flag | 2 findings, exit 1 |
| A3 | #53 `d549987e7...3e3b017e8` | flag | 2 findings, exit 1 |
| A3b | #53 with its fix #56, `d549987e7...6796db41d` | pass | 0, exit 0 |
| A3c | the fix #56 alone | pass | 0, exit 0 |
| A4 | #58 as a PR, `f23d64c99...d9b09882a` | pass | 0, exit 0 |
| A4b | #58's first commit alone, `2da2d1d08^...2da2d1d08` | flag (why the PR is the unit) | 3 findings, exit 1 |
| A5 | this branch vs `origin/ladruno` | pass | 0, exit 0 (0.7 s) |

**Self-test** `ci/test_check_viewer_ledger.py`: 15 cases — the #487 and #485/#53 shapes, the #50
compliant shape, ledger added, modification-only (#484/#55) silent, out-of-scope files,
directory-bounded scope, a new sibling tool directory, rename vs copy, deletions, `-z` parsing,
and three real-git cases (#58's multi-commit PR, a base that moved on after the fork point, an
unresolvable base is exit 2 rather than a pass). Each of six one-line mutations of the checker is
killed by at least one case:

| Mutation | Result |
|---|---|
| ledger check removed | 3 failed |
| three-dot → two-dot diff | 1 failed |
| renames counted as additions | 1 failed |
| copies not counted | 1 failed |
| scope not directory-bounded (`Ladruno_tools` prefix) | 1 failed |
| unresolvable base returns 0 | 1 failed |

**Other gates on this branch:** `python ci/check_quirk_patterns.py` 0 findings (L1–L6, L3 now
including the new guide); its self-test 52 passed; `ruff check` clean on both new files (the fork
has no ruff config; run for hygiene). No C++ build: docs and CI tooling only.

**Viewer checks run while writing the guide** (h5py 3.16 + numpy, no FastAPI on this machine):
`test_contract.py` OK, `monitor_smoke.py` OK, `test_monitor_view.py` reader cases OK. Not run:
`test_api.py` and the `test_monitor_view.py` API cases (need FastAPI), `npm run build` (needs
`npm install`), and any live browser session.

## Live incidents — merge order

V1 checks new PRs only, so these do not block this PR. They are the live instances of R1 — the
ledger rows are still stale:

- `LEDGER_implementations.md`, **Stack profiler** row: cites PRs up to #58; nothing on #484 (the
  untimed-root backfill) or #485 (`profiler_monitor.py`, `monitor_smoke.py`).
- `LEDGER_implementations.md`, **Analysis monitor** row: cites only #64 and still lists the viewer
  ("React live mode") as a follow-up, although #487 shipped `Ladruno_tools/monitor_viewer/`.
- `06_profiler.md` and `08_analysis_monitor.md` logs mention none of #484, #485, #487.

Suggested fix: a small docs-only WP that brings the two rows and two logs current.

## Open questions

- **`make_sample.py` still writes a timed top node**, so the frontend README's dev fixture cannot
  show the #484 class of bug. Add an engine-shaped run to it (tooling change, its own WP)?
- **No viewer test runs in CI.** `test_contract.py` and `monitor_smoke.py` need only h5py and numpy,
  which `static-gates` already installs; `test_api.py` / `test_monitor_view.py` need FastAPI and
  httpx. Worth a step?
- Widen V1 beyond `Ladruno_tools/` (the `AGENTS.md` rule is fork-wide)? Not measured outside the
  viewers.
- R3's cause (screenshot timeouts in the preview tool) was not re-tested here.
- Measurement: after ~10 viewer WPs, count review or post-merge findings that match an existing
  ledger entry against the baseline of 3.
