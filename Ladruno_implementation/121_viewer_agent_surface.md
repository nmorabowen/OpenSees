# WP-121 — Viewer tools join the agents surface: a task guide + a ledger gate

Revision 2 — adversarial review (Opus), 2026-09-25: verdict "ship-with-changes". All six findings
were reproduced before acting and all were applied. Dispositions are in the table under
"Revision 2". Revision 1 was not reviewed.

Status: **built; draft PR #856.** No merge-order constraint. V1 checks a change's own diff, so the
stale ledger rows listed under "Live incidents" do not turn it red.

Scoped 2026-09-25. Branch `wp/121-viewer-agent-surface`, cut from `ladruno` @ `bc5c33453`. Extends
WP-115 ([[115_agents_surface_pilot]]) from elements and materials to the viewer tools only:
- `Ladruno_tools/profiler_viewer/`: `ProfilerResults`, FastAPI `profiler_api`, the
  React/TypeScript/Vite `frontend/`, `launch.py`, `profiler_monitor`, and the `*_smoke_*` pairs.
- `Ladruno_tools/monitor_viewer/`: `MonitorReader`, `monitor_server` + `monitor_page.html`, and
  `monitor_view`.

The rest of the fork is out of scope. Method: the user-level `agent-surface` playbook, phases 0, 2
and 3. Phase 1 (`AGENTS.md`) already exists.

## Problem

The viewers have 14 PRs (#35 to #487, 2026-05-30 to 2026-07-05) and no activity since. Their lessons
live mostly outside the repo:
- two `LEDGER_quirks` entries;
- commit messages (#484, #487);
- plan-doc logs (`06_profiler.md`, `08_analysis_monitor.md`);
- session memory (`project_profiler.md`, `project_analysis_monitor.md`).

Three lessons were written down and still bit again:

| # | Lesson (where it was written) | Recurred | Strength |
|---|---|---|---|
| R1 | "prior PRs ... missed LEDGER_implementations' row; check that row explicitly" — session memory, 2026-05-31, after PR #56 had to backfill #53/#55. The standing `CLAUDE.md` ledger rule said the same. | 2026-07-05: #485 added `profiler_monitor.py`, #487 added all of `monitor_viewer/`; neither touched any ledger. | strong |
| R2 | Verification on `make_sample.py` data only: P8 (#51, 2026-05-31) was checked panel by panel on the sample; the Series tab was empty on real runs until #52. Recorded as a symptom, not as a rule. | 2026-07-04, #484: flame graph empty and every `share` ~1e8 on 100% of engine-written files, all tests green. | soft |
| R3 | "preview_screenshot kept TIMING OUT ... assert via DOM queries, not screenshots" — session memory, 2026-05-31 (#51). | 2026-07-04 (#484): "preview_screenshot still flaky". | soft |

Baseline for the WP-115 review: **3 viewer recurrences** (1 strong, 2 soft) as of 2026-09-25.

Single occurrences, not recurrences:
- #487's three self-review bugs: 500 instead of 503 on a missing sink, a Follow button that did
  nothing, and `--watch` giving up before the writer started.
- The Vite template's TypeScript strictness traps (`tsc` reports them).

## Shape

1. **One task guide**, `.claude/skills/ladruno-viewer-tools/SKILL.md` (≤ 100 lines). It follows the
   playbook's rule that a UI needs its own visual-verification guide, and covers:
   - the contracts other code reads: the C++ writer, `api.ts`, the monitor sink, apeGmsh's
     `profiler.py`;
   - the data to check against;
   - the checks to run (none of them run in CI);
   - reference-vs-candidate verification driven through the DOM;
   - stopping what you launched, by PID;
   - the ledger items.

   *Accept:* ≤ 100 lines, and every `Quirks: "..."` pointer resolves. L3 of
   `ci/check_quirk_patterns.py` scans all `.claude/skills/*/SKILL.md`, so it now covers this guide.
2. **One row in the `AGENTS.md` task-guide table.** Nothing else in `AGENTS.md` changes.
   *Accept:* `git diff AGENTS.md` is one added line.
3. **Lessons into the repo.** R2 and R3 lived only in memory and a commit message. They become two
   `LEDGER_quirks` entries for the guide to point at, appended at the end per the `merge=union`
   convention. R1 is enforced instead (step 4); its evidence is this document.
4. **V1, `ci/check_viewer_ledger.py` (Revision 2 rule).** A change that brings a new **source**
   file into a tool directory `Ladruno_tools/<tool>/` must **add a line** to
   `LEDGER_implementations.md` that **names `Ladruno_tools/<tool>`**. In other words, it edits that
   tool's row.
   - *New:* status A or C, or a rename whose destination is in a different tool directory than its
     source. That covers a move into `Ladruno_tools/` and a move between tools. A rename inside one
     tool is not new.
   - *Source:* `.py .pyw .ts .tsx .js .jsx .mjs .cjs .html .css .bat .cmd .sh .ps1`, compared
     case-insensitively.
   - *Row edit:* an added line that names the tool, bounded so that `profiler_viewer2` and
     `profiler_viewer.old` do not match. A line that only re-adds or re-spaces a removed line does
     not count.
   - *Range:* three-dot `base...head`, the change since the merge base, so a multi-commit PR whose
     row edit landed in a later commit passes.
   - *Silent by design:* modification-only changes, and whether the edited row actually describes
     the new file.
   - *Exit codes:* 0 clean, 1 finding, 2 the range cannot be evaluated, 3 git cannot run.
   - *CI:* two steps last in `static-gates`, **on every event** (Revision 2). On a PR the range is
     the PR's base and head SHAs; on any other event it is `origin/ladruno...HEAD`. The required job
     name is unchanged.

   *Accept (mutation gate):* V1 flags the incident PRs and passes every compliant PR (Results).
5. **Ledgers:**
   - a `LEDGER_implementations` row (test infrastructure);
   - a `ci/README.md` row;
   - `ci/check_viewer_ledger.py` added to the never-ships list in `upstream_pr_campaign.md`.

   No vanilla file is touched, and `.gitignore` already re-includes `.claude/skills/`.

## Revision 2 — adversarial review dispositions

Every finding below was reproduced by running it before being acted on:

| # | Finding (severity) | Reproduced | Disposition |
|---|---|---|---|
| 1 | **Bypass (MAJOR):** the gate ran on `pull_request` only, but `ladruno` needs no PR and a `workflow_dispatch` run also produces the required check. | 6189f5291 has two suites, `pull_request` run 35946791470 and `workflow_dispatch` run 35946815553. Protection on `ladruno`: required checks `Zone-A (Ubuntu)` and `classTag + manifest gates`, no review requirement, `enforce_admins` on. The dispatch run's checkout fetches `+refs/heads/*:refs/remotes/origin/*`, so `origin/ladruno` resolves. | `if:` dropped; the base/head fall back to `origin/ladruno` / `HEAD` on every non-PR event. Simulated locally: HEAD == `origin/ladruno` (a push to `ladruno` or the schedule) gives an empty diff and exit 0, also when `ladruno` has already moved on; an unmerged branch vs `origin/ladruno` (dispatch) gives the same verdicts as the PR-style range (constructed table). The Revision 1 claim "what lands on `ladruno` already passed as a PR" was false and is withdrawn. |
| 2 | **Any ledger change passed (MAJOR),** including whitespace-only edits or deleting the ledger; every WP-era PR adds its own row. | 14/14 of the latest WP PRs (#840–#854) edit `LEDGER_implementations.md` (31 of the last 39). Revision 1 passed all four bypass shapes (constructed table). | The rule now needs an added line naming `Ladruno_tools/<tool>` that is not a re-add or whitespace variant of a removed line. The survey is unchanged: the same 4 flagged, 10 pass, including all six compliant PRs. |
| 3 | **Moves into scope were silent (MINOR):** `git mv scripts/x.py Ladruno_tools/...` gave 0 findings. `--no-renames` and ignoring the failed diff both survived all 15 tests. | All three reproduced (0 findings; 15/15 passed under both mutations). | A rename counts as new when its destination tool directory differs from its source's; a rename within one tool does not. New real-git tests for both. The two diff calls now share one helper that raises on failure, and the no-merge-base test asserts the message. Both mutations are killed. |
| 4 | **No waiver (MINOR):** a `.gitignore` or fixture under `Ladruno_tools/` was flagged. | #35's `.gitignore` and `README.md` were flagged in Revision 1. | Scoped to source files, with no waiver. In L1–L6 a waiver covers code that matches the pattern but is correct. Here the fix, editing the tool's row, is always available, so a new source file with no row edit is never a legitimate exception. Non-source files (`.gitignore`, docs, JSON/config, fixtures, images) are out of scope. |
| 5 | **"Fix passes" runs were constructed (NIT).** | Yes: A3b spans #53+#55+#56, A3c adds nothing, and #485/#487 were never fixed. | The acceptance tables below say which runs are real PRs and which are constructed. |
| 6 | **Missing git gave a traceback with exit 1 (NIT),** indistinguishable from a finding. | `PATH=/nonexistent` produced a traceback and exit 1. | Now exit 3 with a one-line message; tested. The loud exit 2 on an unresolvable range (a shallow clone, a missing base) deliberately departs from the playbook's "stay silent". That rule covers one unreadable case in a readable tree; here the whole gate cannot run, and a silent pass would switch it off unnoticed. |

## Rejected approaches

- **A code-pattern rule for #484** (`x or 1` denominators, the icicle's `wall_ms <= 0` early return).
  The fix `acb932935` keeps both idioms: `profiler_results.py` still has `eff_root_ns or 1` and
  `wall_ns or 1`, and `Icicle.tsx` is unchanged. It fixed the data, not the idiom. A rule that flags
  the fix does not ship.
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
- **Running V1 on pull requests only (Revision 1).** Refuted by review finding 1: a dispatch run
  produces the required check without V1.
- **Accepting any `LEDGER_implementations.md` change (Revision 1).** Refuted by review finding 2:
  every WP PR edits the ledger for its own row.
- **Scoping V1 to new tool directories only** (the other option for finding 4). It would miss #53
  and #485, which added files to an existing tool: half the incidents.
- **A waiver mechanism** (e.g. a commit-message trailer). Compliance is always possible (see finding
  4), so a waiver would only be a way around the rule.
- **Carrying #487's "h5py raises a bare `OSError`, not `FileNotFoundError`" into the ledger.** It
  does not reproduce: h5py 3.16.0 raises `FileNotFoundError` (an `OSError` subclass) for a missing
  file, with and without `swmr=True`. The comment at `monitor_reader.py:41` is wrong on this version,
  and the `os.path.exists` guard it motivated is harmless. Older h5py was not tested.
- **A second guide for the monitor viewer.** Both viewers share the stack (h5py, FastAPI, a browser
  page) and the verification loop.
- **Pixel-diff screenshot baselines.** No browser in CI, and screenshots were the flaky part (R3).
- **Backfilling the stale profiler/monitor rows here.** Out of scope for this WP; listed below.

## Results (Revision 2, 2026-09-25)

### V1 survey

Every PR that touched `Ladruno_tools/`, run with the Revision 2 rule. Squash PRs use `c^...c`; merge
PRs use `M^1...M^2`:

| PR | New source files | Added ledger line names the tool | V1 |
|---|---|---|---|
| #35 | 4 (`profiler_viewer/` created) | no; the row was first added by #46 | **flagged** |
| #46, #49, #50, #51, #52 | 2, 2, 2, 16, 2 | yes | pass |
| #47, #48, #55, #484 | 0 | — | pass (nothing new) |
| #53 | 2 (`header_smoke_*`) | no; backfilled by #56 | **flagged** |
| #58 | 3 (launcher) | yes, in a later commit of the same PR | pass |
| #485 | 2 (`profiler_monitor.py`, `monitor_smoke.py`) | no | **flagged** |
| #487 | 5 (`monitor_viewer/` created; the README is not source) | no | **flagged** |

Four flags, each a real miss of the `AGENTS.md` rule, and no false positive. The results match
Revision 1, so the stricter row rule costs nothing on this history.

### Acceptance on real ranges

| Run | Range | Kind | Expected | Got |
|---|---|---|---|---|
| A1 | #487 `0ec96549a^...0ec96549a` | real incident | flag | 5 findings, exit 1 |
| A2 | #485 `e7e403763^...e7e403763` | real incident | flag | 2 findings, exit 1 |
| A3 | #53 `d549987e7...3e3b017e8` | real incident | flag | 2 findings, exit 1 |
| A4 | #58 as a PR, `f23d64c99...d9b09882a` | real compliant PR (multi-commit) | pass | 0, exit 0 |
| A4b | #58's first commit alone, `2da2d1d08^...2da2d1d08` | a real commit, not a PR | flag (shows why the PR is the unit) | 3 findings, exit 1 |
| A5 | #46, #49, #50, #51, #52 | real compliant PRs | pass | 0 each |
| A6 | #53..#56, `d549987e7...6796db41d` | **constructed**: spans the #53, #55 and #56 merges, not one PR | pass | 0, exit 0 |
| A7 | #56 alone | **trivial**: adds nothing, so it proves nothing about V1 | pass | 0, exit 0 |
| A8 | this branch vs `origin/ladruno` | real | pass | 0, exit 0 |

#485 and #487 were **never fixed**; their rows are still stale ("Live incidents"). So no real
fix-side run exists for them, and the constructed replays below stand in.

### Constructed replays

Each replay is #485's or #487's real tree with a synthetic `LEDGER_implementations.md`, committed
with plumbing in a scratch clone. The last column is the "dispatch" range, `origin/ladruno...head`.

| Tree | Ledger variant | Revision 1 | Revision 2 | Revision 2, dispatch range |
|---|---|---|---|---|
| #487 | a WP-style row that names no viewer | pass (bypass) | **flag** | flag |
| #487 | a row naming `Ladruno_tools/monitor_viewer/*` | pass | pass | pass |
| #487 | whitespace-only edit of the profiler row | pass (bypass) | **flag** | flag |
| #487 | ledger deleted | pass (bypass) | **flag** | flag |
| #487 | real edit of the *profiler* row (wrong tool) | pass (bypass) | **flag** | flag |
| #485 | a WP-style row that names no viewer | pass (bypass) | **flag** | flag |
| #485 | a row naming the *monitor* viewer (wrong tool) | pass (bypass) | **flag** | flag |
| #485 | whitespace-only edit of the profiler row | pass (bypass) | **flag** | flag |
| #485 | ledger deleted | pass (bypass) | **flag** | flag |
| #485 | real edit of the profiler row | pass | pass | pass |

Non-PR events, simulated on a clone after `git fetch`:
- `origin/ladruno...origin/ladruno` (a push to `ladruno`, or the schedule): 0 findings, exit 0.
- `origin/ladruno...origin/ladruno~3` (`ladruno` moved on before the job ran): 0 findings, exit 0.

### Self-test

`ci/test_check_viewer_ledger.py` has **27 cases**. The Revision 1 shapes are kept (rewritten for
the new rule), and one case was added for every review finding:
- which files count: source suffixes (case-insensitive), non-source files out of scope, renames
  within and across tools, a move into scope, copies, deletions;
- what counts as a row edit: another row, the other tool's row, whitespace-only, a moved row, the
  ledger deleted, bounded tool names, a top-level file;
- the `-z` and unified-diff parsers;
- real git: a multi-commit PR, a rename within a tool, a move into scope, a base that moved on, an
  unresolvable base (exit 2), no merge base (exit 2), and git missing (exit 3).

**Mutation gate: 18/18 one-line mutations of the checker killed**, each by a real assertion
failure. The runner retries a run that ends in a pytest *error* rather than a failure, and does not
count it as a kill; one batch had transient errors that were not failures.

| Mutation | Failed |
|---|---|
| `-M` → `--no-renames` (review survivor) | 1 |
| ignore a failed diff / no merge base (review survivor) | 1 |
| diff return code unchecked | 1 |
| row check removed | 4 |
| any added ledger line counts | 4 |
| any ledger change counts (Revision 1 rule) | 7 |
| whitespace-only edit counts | 2 |
| tool name matched unbounded | 1 |
| every file is source | 2 |
| suffix match case-sensitive | 1 |
| tool = whole path, not directory | 5 |
| three-dot → two-dot diff | 2 |
| every rename is new | 2 |
| no rename is new | 2 |
| copies not counted | 1 |
| scope not directory-bounded | 17 |
| unresolvable range returns 0 | 2 |
| git missing returns 1 | 1 |

### Other gates

- `python ci/check_quirk_patterns.py`: 0 findings (L1–L6; L3 now includes the new guide).
- The quirk lint's self-test: 52 passed.
- `ruff check` is clean on both new files. The fork has no ruff config; this was run for hygiene.
- No C++ build: docs and CI tooling only.
- CI (`static-gates` on the PR): Revision 1 run 36194815317 passed, with the self-test at 0.13 s and
  the gate under 1 s. Revision 2's run is linked from the PR.

### Viewer checks run while writing the guide

On h5py 3.16 + numpy, with no FastAPI on this machine:
- **Passed:** `test_contract.py`, `monitor_smoke.py`, and the `test_monitor_view.py` reader cases.
- **Not run:** `test_api.py` and the `test_monitor_view.py` API cases (they need FastAPI),
  `npm run build` (it needs `npm install`), and any live browser session.

## Live incidents — merge order

V1 checks new changes only, so these do not block this PR. They are the live instances of R1 — the
ledger rows are still stale:

- `LEDGER_implementations.md`, **Stack profiler** row: cites PRs up to #58, with nothing on #484 (the
  untimed-root backfill) or #485 (`profiler_monitor.py`, `monitor_smoke.py`).
- `LEDGER_implementations.md`, **Analysis monitor** row: cites only #64 and still lists the viewer
  ("React live mode") as a follow-up, although #487 shipped `Ladruno_tools/monitor_viewer/`. No line
  names `Ladruno_tools/monitor_viewer`, so under Revision 2 the next PR that adds a monitor-viewer
  source file must add one.
- The `06_profiler.md` and `08_analysis_monitor.md` logs mention none of #484, #485, #487.

Suggested fix: a small docs-only WP that brings the two rows and two logs current.

## Open questions

- **`make_sample.py` still writes a timed top node**, so the frontend README's dev fixture cannot
  show the #484 class of bug. Add an engine-shaped run to it (a tooling change, its own WP)?
- **No viewer test runs in CI.** `test_contract.py` and `monitor_smoke.py` need only h5py and numpy,
  which `static-gates` already installs; `test_api.py` / `test_monitor_view.py` also need FastAPI and
  httpx. Worth a step?
- **V1 cannot judge the content of a row edit.** Any non-whitespace edit to a line naming the tool
  passes, for example appending one word to the row.
- **Widen V1 beyond `Ladruno_tools/`?** The `AGENTS.md` rule is fork-wide, but V1 was not measured
  outside the viewers.
- R3's cause (screenshot timeouts in the preview tool) was not re-tested here.
- **Measurement:** after ~10 viewer WPs, count review or post-merge findings that match an existing
  ledger entry, against the baseline of 3.
