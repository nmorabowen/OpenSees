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
- When this skill and the ADR disagree, the ADR wins. When the ADR and a prior measured ADR (67, 68, 75, 77) disagree, the prior ADR wins until the owner decides otherwise (ADR §2.3).

**Dual-Agent Protocol:**
- Live state lives **only** on GitHub: the `WP-260` `program` issue (#1156), pinned board `#937`, and each WP's Draft PR (`### Handoff` comments).
- **Gemini / Antigravity** authors WPs inside `.agents/worktrees/wp260*` (never touching `.claude/` or `CLAUDE.md`). **Claude** runs an adversarial bug-finding review on each Draft PR. That is **not** the ADR-87 D4 review: D4 is an owner-run out-of-family review of a handoff pack, with a committed verdict file (ADR §7).

Quirks: "`HHTGeneralizedExplicit`) add `C` to the LHS tangent", "PARDISO `iparm[27]=1` requires `float*` buffers", "`TRBDF2` uses different effective-stiffness weights on its two legs", "An MKL symbol newer than esmeralda's oneMKL 2024.2 compiles on Windows and breaks the Linux PARDISO build", "Threaded MKL PARDISO is NOT byte-reproducible run-to-run", "A held-load augmentation loop that re-enters `Domain::commit()` SILENTLY corrupts the recorder stream".

---

## Orchestrator & Worker Loop

### 1. Boot
1. Read pinned board `#937` (`gh issue view 937 --repo nmorabowen/OpenSees`) and the `WP-260` `program` issue (`gh issue view 1156 --repo nmorabowen/OpenSees --comments`).
2. Check if the target WP already has an open PR (`gh pr list --repo nmorabowen/OpenSees --search "WP-260 in:title"`) and read its last `### Handoff` comment.
3. Read [references/contracts.md](references/contracts.md) (contracts, mixed-precision policy, Tier-3 traps) and [references/fork-checklist.md](references/fork-checklist.md).
4. **Is this WP picked?** Only `260-w0` and the packages the owner picked from its ranking may start (plan §0.5: phase share ≥ 10 %, projected whole-run gain ≥ 1.10×, top 2–3 at a time). Check the latest pick on the `program` issue; if the WP isn't on it, stop.
5. Check the WP's **Blocked by** line in [references/wp-briefs.md](references/wp-briefs.md): if it names an open PR (#1147/#1153/#1154) or an owner decision (D4/D6), stop and post a handoff instead of starting.
6. Fetch fresh `origin/ladruno` (`git fetch origin ladruno`) and cut an isolated in-workspace worktree:
   `git worktree add .agents/worktrees/wp260<x> -b wp/260<x>-<slug> origin/ladruno`.

### 2. Implement (Disjoint WPs)
- Copy the WP's brief from [references/wp-briefs.md](references/wp-briefs.md).
- Respect live de-duplication: never re-implement `WP-253` (#1147), `WP-258` (#1153), or `WP-259` (#1154), and never edit a file they edit while they are open (ADR §2.2 table).
- Mark every edited upstream file with `// Ladruno ...` and write per-WP ledger fragments under `Ladruno_implementation/ledger/{implementations,vanilla,quirks}/WP-260<x>-<slug>.md` (never edit `LEDGER_*.md`).
- If adding new fork classes/commands, complete all registrations in [references/fork-checklist.md](references/fork-checklist.md).

### 3. Build & Verify (Serialized Local Build)
1. **Static gates first:** run the exact list in [references/fork-checklist.md §1](references/fork-checklist.md) (the same commands `.github/workflows/ladruno.yml` runs on drafts). Do **not** glob `ci/check_*.py`: some need arguments (`check_mkl_compat.py --build-dir`, `check_viewer_ledger.py --base/--head`, `check_tcl_results.py`).
2. **MKL edits** (anything under `SRC/system_of_eqn/linearSOE/pardiso/`, `_PARDISO` / `_LADRUNO_MKL_FEAST` blocks, `mkl_rci.h`): every routine, macro and `iparm` value must exist in **oneMKL 2024.2.2** (esmeralda). Prove it with `ci/check_mkl_compat.py` (fork-checklist §1b); your local 2025.1 is not evidence.
3. **Local 5-target build (one worktree at a time on Windows):**
   ```cmd
   cmd.exe /c "cd /d <worktree> && call Ladruno_scripts\setup_env.bat && Ladruno_scripts\build.bat"
   ```
   Confirm `ladrunoBuild` matches the worktree commit SHA on both `dist\bin\OpenSees.exe` and `dist\bin\opensees.pyd`.
4. **ADR-68 Gates + Profiler Numbers:**
   - Run targeted `G-BYTE` / `G-EQUIV` / `G-OPT` pytest & Tcl tests (`pytest -m zone_a`). Anything inside the WP-132 `-deterministic` block of `PARDISOGenLinSolver.cpp` re-proves CNR bit identity.
   - Profile before/after against `origin/ladruno` on the D1 decks the WP was picked on (method as in the `260-w0` brief: 1-thread `-deep` for phase shares, production threads for wall time, 3 interleaved runs, median). Report the phase speedup `s` and the whole-run gain, and paste the table into the PR (never claim unmeasured speedups).

### 4. Draft PR, Zone-A & Review
1. Push branch and open a **Draft PR** targeting `ladruno` on day one (`in-flight` label).
2. Trigger Zone-A Ubuntu CI (`gh workflow run ladruno.yml --ref <branch>`) and verify status via `gh pr checks <pr>` / `gh run view <id>`.
3. Post a `### Review request — 260<x> — <sha>` comment (template and checklist: [references/review-protocol.md](references/review-protocol.md)). Claude reviews at that SHA and posts numbered findings + a verdict. Answer every finding by number (`fixed in <sha>` / `declined: <rationale>`), re-verify gates, and post a `Ladruno_internal/ISSUES.md` §3 `### Handoff` comment.
4. For a family-level PR (per owner decision D7): write `Ladruno_implementation/reviews/handoff_260<x>.md` (ADR section, key diffs, specific questions) for the owner's ADR-87 D4 out-of-family run; the verdict lands as `Ladruno_implementation/reviews/adr260<x>_verdict.md`.
5. After merge, post the measured whole-run gain on the D1 decks to the `program` issue so the menu is re-ranked before the next pick (plan §0.5 rule 5).
6. Flip Ready for Review only when the entire warrant package is green and Claude's verdict is `ready`; the owner merges (`--merge`, never `--auto` or `--admin`).

### Stop Rules
Stop and post a `### Handoff` comment when:
- A gate or benchmark reveals a contract conflict requiring an ADR update,
- The WP turns out to re-open something a prior ADR closed (ADR §2.3),
- An owner decision (`D1–D7` or `human-gate`) is reached, or
- A WP is ready for Claude review / owner merge.
