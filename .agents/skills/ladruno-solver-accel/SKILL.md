---
name: ladruno-solver-accel
description: >
  Orchestrate and implement the Ladruno Solver & Wave Acceleration program (WP-260) across
  implicit and explicit nonlinear mechanics, material return mapping, global SOE/assembly,
  explicit wave engines, and modern integrators. Use whenever starting, continuing, building,
  benchmarking, or reviewing any WP-260 work package.
---

# `ladruno-solver-accel` — WP-260 Operating Guide

**Authority:**
- [ADR `260_ladruno_solver_acceleration_adr.md`](../../../Ladruno_implementation/260_ladruno_solver_acceleration_adr.md)
- [Plan `_wp260_orchestration_plan.md`](../../../Ladruno_implementation/_wp260_orchestration_plan.md)
- Detailed technical blueprints live in the `program` issue and local scoping docs. When this skill and the ADR disagree, the ADR wins.

**Dual-Agent Protocol:**
- Live state lives **only** on GitHub: the `WP-260` `program` issue (#1156), pinned board `#937`, and each WP's Draft PR (`### Handoff` comments).
- **Gemini / Antigravity** authors WPs inside `.agents/worktrees/wp260*` (never touching `.claude/` or `CLAUDE.md`), and **Claude** runs the out-of-family adversarial review (ADR-87 D4) on each PR.

---

## Orchestrator & Worker Loop

### 1. Boot
1. Read pinned board `#937` (`gh issue view 937 --repo nmorabowen/OpenSees`) and the `WP-260` `program` issue (`gh issue list --repo nmorabowen/OpenSees --label program --search "WP-260 in:title"`).
2. Check if the target WP already has an open PR (`gh pr list --repo nmorabowen/OpenSees --search "WP-260 in:title"`) and read its last `### Handoff` comment.
3. Read [references/contracts.md](references/contracts.md) (locked C++ contracts, mixed-precision policy, Tier-3 traps) and [references/fork-checklist.md](references/fork-checklist.md).
4. Fetch fresh `origin/ladruno` (`git fetch origin ladruno`) and cut an isolated in-workspace worktree:
   `git worktree add .agents/worktrees/wp260<x> -b wp/260<x>-<slug> origin/ladruno`.

### 2. Implement (Disjoint WPs)
- Copy the WP's brief from [references/wp-briefs.md](references/wp-briefs.md).
- Respect live de-duplication: never re-implement `WP-253` (#1147), `WP-258` (#1153), or `WP-259` (#1154).
- Mark every edited upstream file with `// Ladruno ...` and write per-WP ledger fragments under `Ladruno_implementation/ledger/{implementations,vanilla,quirks}/WP-260<x>-<slug>.md` (never edit `LEDGER_*.md`).
- If adding new fork classes/commands, complete all registrations in [references/fork-checklist.md](references/fork-checklist.md).

### 3. Build & Verify (Serialized Local Build)
1. **Static gates first:** Run all `ci/check_*.py` + `stamp_headers.py --check` from [references/fork-checklist.md](references/fork-checklist.md).
2. **Local 5-target build (one worktree at a time on Windows):**
   ```cmd
   cmd.exe /c "cd /d <worktree> && call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat"
   ```
   Confirm `ops.ladrunoBuild()` matches the worktree commit SHA on both `dist\bin\OpenSees.exe` and `dist\bin\opensees.pyd`.
3. **ADR-68 Gates + Profiler Numbers:**
   - Run targeted `G-BYTE` / `G-EQUIV` / `G-OPT` pytest & Tcl tests (`pytest -m zone_a`).
   - Run `-perStep` + `-deep` before/after comparison against `origin/ladruno` and paste exact timings into the PR description (never claim unmeasured speedups).

### 4. Draft PR, Zone-A & Claude Review
1. Push branch and open a **Draft PR** targeting `ladruno` on day one (`in-flight` label).
2. Trigger Zone-A Ubuntu CI (`gh workflow run ladruno.yml --ref <branch>`) and verify status via `gh pr checks <pr>` / `gh run view <id>`.
3. Submit the green Draft PR for **Claude adversarial review**, verify every finding in-tree, fix or log declined items with rationale, re-verify gates, and post a `Ladruno_internal/ISSUES.md` §3 `### Handoff` comment.
4. Flip Ready for Review only when the entire warrant package is green; the owner merges (`--merge`, never `--auto` or `--admin`).

### Stop Rules
Stop and post a `### Handoff` comment when:
- A gate or benchmark reveals a contract conflict requiring an ADR update,
- An owner decision (`D1–D4` or `human-gate`) is reached, or
- A WP is ready for Claude review / owner merge.
