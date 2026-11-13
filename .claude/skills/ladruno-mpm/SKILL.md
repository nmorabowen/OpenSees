---
name: ladruno-mpm
description: >
  Run and implement the Ladruno MPM program: WP-232, a general Material Point Method driver for
  the OpenSees fork, with apeGmsh as its front and back end. Use when the user runs
  "/ladruno-mpm <wp>" (232a, 232b, 232c, 232d, R, 232e, 232f…) or asks to start, continue,
  orchestrate, review or unblock an MPM work package; and before writing or reviewing MPM code:
  the explicit engine and its LoadPattern host, the `mpm` commands, particle elements, P2G / G2P,
  GIMP, MUSL, PIC / FLIP, the incremental hypo-kernel entry, the particle-set and `.ladruno`
  particle contracts, lattice seeding, or the apeGmshViewer points layer. Holds the work-package
  map, the model and effort for each role, the dispatch briefs, the frozen design and the walls.
---

# `/ladruno-mpm <wp>` — run one MPM work package

**Authority:** the ADR `Ladruno_implementation/232_ladruno_mpm_general_driver_adr.md` and the plan
`Ladruno_implementation/_wp232_orchestration_plan.md`. This skill is the operating layer on top of
them. When the skill and the ADR disagree, the ADR wins; fix the skill in the same PR.

**Two modes.**
- **Orchestrator** — the session that runs `/ladruno-mpm <wp>`. It plans, dispatches, verifies and
  records. It never edits `SRC/`, `tests/` or apeGmsh code itself.
- **Worker** — an agent dispatched with a brief from [references/wp-briefs.md](references/wp-briefs.md).
  It implements one WP in its own worktree and follows the worker rules below,
  [references/architecture.md](references/architecture.md) and
  [references/fork-checklist.md](references/fork-checklist.md).

## The map

```
WP-232 (ADR, plan, this skill) ── the owner merges it first
  │
  ├─ Lane C: 232a contracts + apeGmsh ADR ──► owner accepts the apeGmsh ADR ──► 232d apeGmsh: seeding, readers, points layer
  │            └─ Contract-2 schema ─────────────────┐
  └─ Lane N: 232b numpy oracle ──► interface frozen ──► 232c fork slice 1 ──► R: release + FORK_PIN ──► 232e apeGmsh bridge
then 232f+ = the ADR's §9 phases, each started with /ladruno-mpm 232f, 232g, …
```

- **232a and 232b start together.** Neither waits for the other.
- **232c starts when 232b's interface is frozen**, not when 232b finishes. The interface is frozen
  when 232b has pushed `tests/mpm_interface.md` and a first `tests/mpm_cases.txt`, and the
  orchestrator has reviewed them and said so on the program issue. 232c's `.ladruno` writer also
  waits for 232a's Contract-2 schema.
- **232d starts when the owner accepts the apeGmsh ADR.** It builds against 232a's fixture files,
  not against a fork build.
- **R needs 232c merged; 232e needs R.**
- **Merge order on `ladruno`:** WP-232 → 232b → 232a (fork half) → 232c. The owner merges every PR.

## Roster — set `model` and `effort` on every Agent call

| Role | `model` · `effort` | WPs |
|---|---|---|
| Orchestrator (this session's own setting) | fable · high; max when adjudicating a review, declaring M0b, or deciding the D-list | all |
| Contracts and the apeGmsh ADR draft | opus · high | 232a |
| Numpy oracle | opus · high | 232b |
| Fork C++ (never the oracle's agent) | opus · high | 232c, 232f+ |
| apeGmsh readers (Python, TypeScript) and the three.js layer | opus · high | 232d |
| apeGmsh seeding plumbing; bridge primitives | sonnet · medium | 232d, 232e |
| Critic at ready-flip, from the other family | fable · high for Opus or Sonnet authors; opus · high for Fable authors | every WP |
| Measurement harnesses | sonnet · medium | M0b |
| Read-only sweeps | the `Explore` agent, "very thorough" | any |

This skill is the instruction that authorizes `effort` (the Agent tool sets it only when
instructions ask). At boot, confirm the Agent tool still exposes `effort`; if it does not, write
"effort: harness default" in the handoff.

## Orchestrator loop

### 1. Boot (at most 8 tool calls)

1. **Duplicates.** `ListAgents` for a live session on this WP;
   `gh pr list --repo nmorabowen/OpenSees --search "WP-<n> in:title"` (and `--repo
   nmorabowen/apeGmsh` for 232d / 232e); `gh issue list --label in-flight`. A live duplicate →
   stop and report it.
2. **Program issue.** `gh issue list --repo nmorabowen/OpenSees --label program --search "MPM in:title"`,
   then `gh issue view <n> --comments`. No program issue → stop: opening it is the owner's
   go/no-go item 2.
3. **The WP's draft PR**, if it exists: `gh pr view <n> --comments`. The last `### Handoff` is
   where you resume.
4. **Gates** (the map above). A red gate → post a `blocked: <gate>` handoff and stop.
5. **Fresh refs.** `git fetch origin`; read fork code from `origin/ladruno` and apeGmsh code from
   its `origin/main`.

### 2. Dispatch

- Copy the WP's brief from [references/wp-briefs.md](references/wp-briefs.md), fill in the PR
  number and the last handoff's **Next**, and dispatch it with `subagent_type: general-purpose`,
  the roster's `model` and `effort`, and `isolation: "worktree"`.
- Send independent workers in one message (parallel Agent calls). Subagents run in the background
  and notify you when done; never poll them.
- Add `in-flight` to the WP's PR or issue.

### 3. Verify — never accept a worker's word

- Read the files and logs the report names. Re-run the cheapest gate yourself: the oracle cases,
  and the `ci/check_*.py` gates in [references/fork-checklist.md](references/fork-checklist.md).
- **You launch every build**, in the worker's worktree, with the WMI recipe in the fork checklist
  and a bounded wait on the log. A subagent that starts its own background build is never woken
  by it.
- Open every probe with `ops.ladrunoBuild()` (Python) or `ladrunoBuild` (Tcl) and compare the hash
  with the worker's head, on both `OpenSees.exe` and `dist\bin\opensees.pyd`.
- Verify every review finding in-tree before accepting it, and log a declined fix with its reason
  (precedent: ADR §14.3, M6).

### 4. Review — once per WP, at ready-flip

- One critic from the other family (roster), briefed with the diff and the ADR sections it
  implements — not with your conclusions.
- `/code-review <level>`: `low` for docs, `medium` for the oracle, `high` for C++ and contracts.
- 232c also needs the out-of-family verdict (ADR-87 D4), which the owner runs: a `human-gate`.

### 5. Land

- Flip to ready only when the WP's done-when holds and its warrant package is in the PR (fork
  checklist, "Warrant package").
- A ready-flip starts no Zone-A run: `gh workflow run ladruno.yml --ref <branch>`, then read CI
  with the `ccd_pr` tools (`get_status`; `bind_pr` if the PR is not listed). Never poll or
  schedule CI checks yourself.
- **Who merges.** The owner — unless the program issue records that the owner has delegated fork
  merges to the orchestrator (done for this program on 2026-11-13). Then merge a fork WP PR with
  `gh pr merge <n> --merge` (a merge commit; `--squash` only for a small fix), and only when every
  required check is green and the WP's done-when holds. Never `--auto`, never `--admin`, never
  around a red check. Without a delegation: add `human-gate` and say what is ready.
- **apeGmsh PRs are never merged from here** — see "apeGmsh — only through its T coordinator".

### 6. Record and chain

1. Post the handoff (`Ladruno_internal/ISSUES.md` §3 template) on the WP's PR, one line on the
   program issue, and refresh the WP's row on board #937. apeGmsh's board #1203 belongs to its T
   coordinator: report to it instead.
2. Commit and push everything before queuing anything: a started chip can run in this worktree and
   commit what you left uncommitted.
3. Queue the next fork WP with `spawn_task`, prompt `/ladruno-mpm <next>`. apeGmsh work (232a's
   ADR, 232d, 232e, R's pin bump) is queued with the T coordinator, never by chip.
4. A failure class seen for the second time → open a `lesson` issue.

### Stop when

the WP is done; a human gate is reached; a gate is red; or your context is getting large — then post
the handoff and queue a continuation chip for the same WP.

## apeGmsh — only through its T coordinator

apeGmsh runs its own program with one live coordinator session (`/apegmsh-program T`; coordinator
issue #1202, board #1203). Every apeGmsh PR, merge, board row and label for this program is
arranged with it (the owner's rule, 2026-11-13).

- Find it with `ListAgents` (a live session titled "T coordinator") and message it with
  `SendMessage` (load the tool with ToolSearch `select:SendMessage`). Introduce the work, hand it
  the brief or the draft, and report to it at every landing, blocker and handoff, as its protocol
  asks of every program session.
- If no T session is live, post on #1202 and wait. Never land apeGmsh work around it.
- A worker isolated in a fork worktree cannot write in the apeGmsh checkout. apeGmsh changes are
  made by a session in apeGmsh, as T arranges; fork workers only draft apeGmsh text in the fork
  branch.

## Worker rules (every brief tells the worker to read this section)

1. **Your worktree only.** Branch `wp/<n>-<slug>` from `origin/ladruno` (apeGmsh: its own branch
   rule, from `origin/main`). Open the draft PR on your first push and commit continuously.
   Report your worktree path.
2. **The design is frozen.** Code against [references/architecture.md](references/architecture.md);
   the ADR is the authority. A design objection goes in your report as a proposed ADR amendment —
   never as a silent deviation.
3. **Quote, don't recall.** Before claiming "X works with Y", "X reuses Y" or "Y is missing", open
   X's code and its own scope or limits section and quote it. Open every ADR number on the branch
   you cite it from.
4. **The walls** (anchors in architecture.md). W1: integrators ignore node writes between steps.
   W4: recorders cache pointers. W9: `ops_TheActiveElement` is a live global — set it to 0 around
   particle updates. W10: `analyze()` sub-steps — call `analyzeStep`. W11: G2P before recorders.
   W12: stage switches walk Domain elements only. W15: never a zero-DOF element in a Domain.
   W16: `StaticAnalysis` has no `analyzeStep`.
5. **Licences.** Never open GPL, LGPL or unlicensed source (Anura3D, Karamelo, AMPLE, GeoTaichi,
   Matter, NairnMPM). Kratos (BSD-4) may enter the fork with its acknowledgement, never apeGmsh.
6. **Independence.** The oracle (232b) reads nothing under `SRC/`. The C++ worker (232c) reads
   `tests/mpm_interface.md` and `tests/mpm_cases.txt`, never `tests/mpm_reference.py`.
7. **Builds.** Batch your C++ edits, then ask the orchestrator for one build of all five targets.
   Do not start `build.bat` in the background yourself.
8. **Ledgers.** Per-WP fragments under `Ladruno_implementation/ledger/`; never edit `LEDGER_*.md`;
   run `python3.12 ci/check_ledger_fragments.py` before every push.
9. **Report in at most 300 words:** worktree path and head sha; what is done; what is verified
   (commands and results); what is blocked; the next action.

## Speed rules

- Two lanes from day one (232a and 232b); 232c on the frozen interface; 232d on fixture files.
- One review per WP, at ready-flip — not one per draft iteration.
- One build per round, all five targets.
- Owner decisions batched: D1–D7 in one sitting. The apeGmsh ADR acceptance is the one
  mid-stream owner gate besides merges.
- A worker on its third fix round for the same finding stops; the orchestrator takes it at max
  effort.
