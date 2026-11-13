---
title: "WP-232 orchestration plan — work packages, models, effort, and protocol for the general MPM driver"
project: Ladruno
status: draft — awaiting the owner's approval (§5). Third review folded (ADR §14.3)
orchestrator: "the authoring session (Opus 5.5), the owner's choice on 2026-11-13; for a fresh session the recommendation stays Fable @ high, max at the decision points"
companion: 232_ladruno_mpm_general_driver_adr
updated: 2026-11-13
---

# WP-232 orchestration plan

ADR-87 requires that agent assignment and effort per work package be decided at kickoff and
recorded; `Ladruno_internal/ISSUES.md` requires that live state live on GitHub. This plan does the
first and names the GitHub homes for the second. It supersedes the ADR-98 orchestration draft.

## 0. Two mechanics to keep straight

- **Model and effort are both settable per spawned agent** in the current harness: the `Agent` tool
  takes a `model` and an `effort` (low … max), and effort is to be set when the owner or project
  instructions ask for it — once approved, this plan is that instruction. This corrects the first
  draft, which said the tool took a model only (ADR §14.3, m11). It is harness-dependent: at each
  kickoff, check that the tool still exposes `effort` before relying on it. Agent definitions under
  `.claude/agents/` are not needed for this; the repo has none and ignores them (`.gitignore`:
  `.claude/*` except `skills/`). For a PR, effort is also the `/code-review` level.
- **GitHub holds the state.** The MPM `program` issue holds the plan's live status; each WP's draft
  PR holds its state and handoffs (header `### Handoff — WP-<nnn> — <date> — <model> @ <effort>`);
  the fork board (#937) holds the fork's row, and apeGmsh's program board (#1203) holds the apeGmsh
  slices. Private memory holds only pointers.

## 1. Roles

| Role | Model · effort | Why |
|---|---|---|
| Orchestrator — sequencing, synthesis, decisions, launching builds and long runs | **Fable · high**; **max** at the decision points (adjudicating a review, declaring the M0b result, the D-list) | house precedent (ADR-90, ADR-94); a fresh session is independent of the author's context. For this program the owner chose the authoring Opus session (2026-11-13) |
| Critics — the next review of the folded ADR, and each family PR | **Fable · high** lead; Opus lenses added for family PRs | Fable found every blocker in this lineage; critics must not share the author's context |
| Numpy oracles — the incremental kernel entry, the explicit engine reference, two-phase checks | **Opus · high** | correctness-critical numerics; independence protocol (the oracle reads nothing under `SRC/`) |
| Fork C++ — the explicit engine and its `LoadPattern` host, the `mpm` commands, particle elements, the kernel entry, the `.ladruno` particle group | **Opus · high**, one agent per WP continued across review rounds | integrator interaction, commit ordering, three-engine registration, a format bump |
| Contract design — the two contracts and the apeGmsh ADR draft | **Opus · high**, with the fork's and apeGmsh's checklists in the prompt | a contract error costs both repos; the ADR is accepted only by apeGmsh's maintainer |
| apeGmsh — the `.ladruno` particle readers (Python and TypeScript), the three.js points layer | **Opus · high** | TypeScript feature work in apeGmshViewer under a new ADR, not pattern-following (ADR §14.3, B2) |
| apeGmsh — lattice seeding and emit plumbing, under the accepted ADR | **Sonnet · medium** | follows the ADR and apeGmsh's checklists |
| Measurement harnesses (M0b) | **Sonnet · medium**; the orchestrator launches the runs | written once, run many times |
| Read-only sweeps | **Explore (Sonnet)**, "very thorough" | keeps the orchestrator's context lean |
| Out-of-family verdict (ADR-87 D4) | **SuperGrok**, run by the owner, at the first family PR | the citable review |
| Haiku | not used | nothing here is mechanical enough |

apeGmsh slices run under apeGmsh's own program protocol (the `apegmsh-program` skill, board #1203,
its pinned workers and cross-family review). The models above are this plan's recommendation there,
not an override of that protocol.

## 2. Work packages

Numbers after 232 are claimed through the program issue, never by guessing the next free one.

| WP | Repo | Content | Working agent | `/code-review` | Gate |
|---|---|---|---|---|---|
| **232** | fork | the ADR, this plan, the quirk fragments, the `ladruno-mpm` skill; three reviews folded | orchestrator + Fable critic | `low` (docs) | owner sign-off on §5; merged first, so every later branch carries the skill |
| **232a** | both | **the apeGmsh ADR** (particle-set stream, particle-cloud results) and **the two contracts** as field-by-field schemas (ADR §7.2–7.3) | Opus · high | `high` | the apeGmsh ADR **accepted by its maintainer**; both repos' checklists ticked; a fixture file each way |
| **232b** | fork | numpy oracle: the explicit MUSL engine and the incremental kernel entry (objectivity included) | Opus · high | `medium` | vibrating bar, self-weight column (with local damping), rigid rotation, two-body impact reproduced |
| **232c** | fork | **slice 1:** the explicit engine hosted as a `LoadPattern`, `mpm` as one `LadrunoCommandTable.h` row, the admitted materials, the `materialState` stage path, critical Δt, local damping, the `.ladruno` particle group with its `FORMAT_VERSION` bump and format document; manifest + mutation gate + guide | Opus · high (a different agent from 232b's) | `high` + out-of-family verdict | the S1 rows of ADR §8 green against 232b; mutants red; three-engine parity |
| **232d** | apeGmsh | lattice seeding writing the particle set at emit; the `.ladruno` particle readers in Python and TypeScript; the three.js points layer; the visual check | Opus · high (readers, viewer); Sonnet · medium (seeding plumbing) | per apeGmsh rules, through its T coordinator | the accepted apeGmsh ADR; apeGmsh's checklists; the visual check |
| **R** | both | a fork release containing `mpm` (the release recipe), then apeGmsh's `FORK_PIN` bumped — all four keys together | the owner, with the orchestrator | — | the release published; `live-fork` green on the new pin |
| **232e** | apeGmsh | bridge primitives for the `mpm` commands: `ladruno_fork` marker, capability-map row, probes that return `None` when inconclusive, staged runs deck-only | Sonnet · medium | per apeGmsh rules, through its T coordinator | **only after 232c has merged into `ladruno` and R is done** |
| 232f+ | fork / both | the phase order in the ADR's §9 — implicit via A, M0b measurements, B-spline, anti-locking, soils, two-phase, contact, the seam, the handoff, MPM → FE remeshing, parallel | re-decided at each kickoff | — | per WP |

## 3. Orchestrator protocol

The operating layer is the `ladruno-mpm` skill (`.claude/skills/ladruno-mpm/`): a session runs
`/ladruno-mpm <wp>`, and the skill holds the boot sequence, the dispatch briefs, the worker rules,
the fork checklist and the frozen design. The rules below are its summary.

1. Boot from the board and the program issue; claim with `in-flight`; post a handoff before stopping.
2. Every WP opens as a draft PR on day one; merges follow rule 8; never `--auto`.
3. Launch builds and long runs yourself (WMI-launched build, a bounded foreground loop on the
   artifact's mtime); never let a subagent background its own build.
4. Read result files, not transcripts; verify every review finding in-tree before accepting it, and
   say so when a proposed fix is declined (ADR §14.3, M6).
5. Before writing or accepting any claim of the form "works with X", "reuses X" or "X is missing",
   open X's code and its own scope or limits section, and quote it beside the claim.
6. Every ADR number cited is opened on the branch it is cited from; apeGmsh claims cite apeGmsh
   `origin/main` at a named commit.
7. Keep the oracle author and the C++ author different agents, meeting only at a frozen interface.
8. **Merges (owner, 2026-11-13: "handle the merges, I trust you"):** the orchestrator merges fork
   WP PRs — merge commits, only on green required checks, never `--auto` or `--admin`. Recorded on
   the program issue.
9. **apeGmsh (owner, 2026-11-13):** every apeGmsh PR, merge and board change goes through
   apeGmsh's T coordinator session (`/apegmsh-program T`, #1202, #1203). The owner approved its
   routing the same day: apeGmsh chain M, nmorabowen/apeGmsh#2972 — M1 the ADR, M2 = 232d,
   M3 = the `FORK_PIN` bump + 232e.

## 4. Standing rules for every agent (verbatim in the prompt)

1. Work only in your own worktree. Read fork code from `origin/ladruno` when the checkout may be stale.
2. Licence rule: never open GPL / LGPL / unlicensed source (Anura3D, Karamelo, AMPLE, GeoTaichi,
   Matter, NairnMPM). Kratos code may enter the fork, never apeGmsh.
3. Integrators ignore node writes between steps (W1), and `analyzeStep` re-seeds them only when the
   domain stamp moved; mode A raises `Domain::domainChange()` every step.
4. No stock recorder may watch an entity that can be removed (W4). `LadrunoRecorder` re-maps only
   when its topology fingerprint changes (WP-165), so it does not follow particle re-pointing either;
   particle output is the `.ladruno` particle group.
5. Never put a zero-node element in a Domain: the first analysis set up there exits the process (W15).
6. A non-element host sets `ops_TheActiveElement` to 0 on every worker thread around material
   updates (W9).
7. Ledger entries are per-WP fragments under `Ladruno_implementation/ledger/`; never edit `LEDGER_*.md`.
8. Every reported number names the build hash, the deck, the output path and the date.

## 5. Go / no-go items for the owner

Owner's "go" on 2026-11-13 ("go, orchestrate from here"): items 1–4 are settled as below.

1. ~~Approve §1 and §2, or reassign.~~ Approved.
2. ~~Open the MPM `program` issue claiming WP-232 onward; add the row on #937 and a link on
   #1203.~~ Done by the orchestrator; the #1203 link goes through apeGmsh's T coordinator.
3. ~~Choose the orchestrator session.~~ The owner chose the authoring session (Opus), not a fresh
   Fable one.
4. ~~Approve contract-first co-design across the two repos.~~ Approved.
5. Decide ADR D7: amend ADR-86 D2, moving the FE → MPM handoff to phase 9. **Open.**
6. As apeGmsh's maintainer, accept (or reject) the apeGmsh ADR when 232a drafts it — an agent
   cannot. **Pending 232a.**
