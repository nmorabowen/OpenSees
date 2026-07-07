---
title: "Modal-analysis family (ADR 45) — session handoff: P-A + P-B(serial) COMPLETE"
project: Ladruno
type: handoff
related:
  - "[[45_ladruno_modal_family_roadmap_adr]]"
  - "[[46_ladruno_complex_modal_adr]]"
  - "[[43_ladruno_feast_eigensolver_adr]]"
  - "[[LadrunoComplexEigen_guide]]"
  - "[[LEDGER_quirks]]"
updated: 2026-07-07
---

# Modal family handoff — 2026-07-07

## State: umbrella phases P-A and P-B(serial) SHIPPED; gates G-A + G-B CLOSED

| PR | What | Gate |
|---|---|---|
| [#505](https://github.com/nmorabowen/OpenSees/pull/505) | Umbrella/docs review fixes (added the missing G-E gate, D1 reword, ADR-47 collision) | — |
| [#506](https://github.com/nmorabowen/OpenSees/pull/506) | ADR 46 P0 — reduced-pencil QZ kernel (`complexEigen -qz`), 33019 reserved | numpy companion-form oracle @1e-10 |
| [#507](https://github.com/nmorabowen/OpenSees/pull/507) | ADR 46 P1 — Route-A closed-form Rayleigh + Domain getters, 33019 ACTIVE | classical ζ/ω_d @1e-8 |
| [#510](https://github.com/nmorabowen/OpenSees/pull/510) | ADR 46 P2 — Route-B assembled projection (`LadrunoDampingAssembler`), DEFAULT | 2-DOF non-classical @1e-8, B==A @1e-10 |
| [#514](https://github.com/nmorabowen/OpenSees/pull/514) | ADR 46 P3 — phased mode shapes + recorder + banner = v1 ship | **G-A**: stick model vs log-decrement @5% |
| [#515](https://github.com/nmorabowen/OpenSees/pull/515) | ADR 43 P1 — `eigen -feast fmin fmax` (33022/33023) | ARPACK parity + band-targeting + ΦᵀMΦ==I |
| [#517](https://github.com/nmorabowen/OpenSees/pull/517) | ADR 43 P2 — `-certify` Sturm/inertia certificate + banner | **G-B**: certified==m, edge-refusal |

Batteries: `tests/test_ladrunoComplexEigen.py` (21) + `tests/test_feastEigen.py` (10,
`skipif` non-win32 — FEAST needs MKL, Windows/oneAPI build only). Five Opus
adversarial gates ran; every finding is fixed-with-regression or in
[[LEDGER_quirks]]. User guide: [[LadrunoComplexEigen_guide]].

## What the next session picks up (in umbrella order)

### 1. The D2 spike — RUN 2026-07-07, decision made (see ADR 43 §5.2/§9 R1, umbrella §5 D2,
[[feast_d2_spike/README]]). No MPI toy build was needed — reading
`MumpsParallelSolver`/`MumpsParallelSOE` plus the local MUMPS build state settled it:
1. **`MumpsParallelSolver` hardcodes `MPI_COMM_WORLD`** — the `mpi_comm` constructor
   arg is accepted and silently discarded (`MumpsParallelSolver.cpp:54-64,97,104-105`).
   Bounded plumbing bug, now [[LEDGER_quirks|ledgered]].
2. **ADR 43 §5.2's "real solver per conjugate pair" claim was wrong** — refuted against
   the FEAST v3/v4 User Guides and MKL's `?feast_srci` reference (Opus
   cross-check): every contour solve is genuinely complex. The fork's local MUMPS is
   real-only (`arith=d`, no `zmumps.lib`) anyway.
3. Numerically verified (numpy, 1e-15 rel. error) that a symmetric real $2n\times2n$
   block-augmented system reproduces the complex solve exactly and can be factored
   LDLᵀ (`SYM=2`) with the *existing* real `dmumps` — zero new dependency.

**Decision: stay on MKL FEAST, no PFEAST vendoring.** P3 (next session, per ADR 43's
phased roadmap) splits into P3a (sub-communicator plumbing fix) → P3b (new symmetric
2n×2n block-real inner SOE) → P3c (`MPI_Comm_split` orchestration + `Allreduce`-sum
projector), then P4 SP/MP build-flag unification unchanged. P3a/b/c and P4 are all
core/parallel-build edits → full adversarial gate (per
[[feedback_adversarial_gate_when]]).

### 2. Demand-driven lanes (build when a project asks)
- **ADR 42 buckling** (33021): rides serial eigen; `-shift` exposure is
  parser-only (`ArpackSOE` plumbing already complete but unreachable —
  `OpenSeesCommands.cpp` hard-zeros it). Gate G-D.
- **ADR 44 frequency domain** (33024): consumes ζ_k from `complexEigen` and
  reuses `LadrunoDampingAssembler` verbatim (umbrella D1). NOTE from the P0
  audit: upstream `ResponseSpectrumAnalysis` performs NO combination — CQC/
  SRSS is an *addition*, favoring the fork-sibling route (D4). Gate G-F
  includes the byte-identical `-combine`-absent regression.
- **ADR 46 deferred slice**: `modalProperties -complex` columns + complex
  participation factors (scope-trimmed at P3, documented in the ADR).

### 3. P-E (complex contours at scale) — only after P-A (done) + P-C; gate G-E
(serial-projection parity + conjugate completeness) is in the umbrella §6.

## Build/test recipes that actually work (this machine)

- Build: `cmd /c Ladruno_scripts\build.bat OpenSeesPy` **from PowerShell**
  (from Bash the cmd wrapper silently no-ops — check the pyd mtime, not the
  exit code). Full log capture: `*>&1 | Out-File` (a `Select-Object` pipe
  discards the diagnostics you'll need).
- Tests in a worktree: `pythoncore-3.12-64\python.exe -S <runner>` with
  `PMI_RANK=0` set BEFORE `site.main()` (stale-pyd `.pth` hijack), paths
  inserted after, and `assert opensees.__file__` contains the worktree.
  Runner template: the scratchpad `run_p1_tests.py` pattern (in-repo tests
  import `_testbed`).
- PR flow: repo auto-merge is DISABLED → background `until gh pr checks`
  watcher, then `gh pr merge --squash`. `ladruno` moves mid-CI on busy days:
  on a merge conflict, `git merge origin/ladruno` locally (ledgers union-merge
  clean; GitHub ignores `merge=union` server-side) and re-watch.

## Load-bearing gotchas (full entries in [[LEDGER_quirks]])

`wipe()` = `clearAll` on the SAME Domain (new domain-level state must
self-reset there; `theEigenvalues` leaked across wipe for years);
`getEigenvalues()`/`getEigenvectors()` exit() when unset (use the additive
`getNum*` probes); `region -rayleigh` bypasses the domain factor copy;
**`Truss` defaults `-doRayleigh` 0** (whole family, not just zeroLength);
Element/Node `getDamp()`/`getMass()` share static scratch (deep-copy first);
recorder files need `-precision 15` for 1e-8 asserts; PARDISO perturbed-pivot
guard needs a COUPLED model to test (isolated diagonal zeros are scaled away).
