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

### 1. The D2 spike (BEFORE any P-C work — decision gate, not a PR)
Question: can MKL FEAST **RCI** (`dfeast_srci`) run its contour solves through a
`MumpsParallelSOE` on an `MPI_Comm_split` sub-communicator? Launch package:
- Toy: ~1k-DOF frame, 2–4 ranks, split into 2 sub-comms; per sub-comm factor
  `(z_j M - K)` with MUMPS (**complex** shift — MUMPS must run its complex
  build, or exploit the conjugate-pair trick in ADR 43 §5.2 if it survives
  scrutiny; this is exactly what P1 dodged by using `dfeast_scsrgv`).
- Outcome A (works) → P3 per-contour MPI + P4 SP/MP build-flag unification
  (`_PARALLEL_PROCESSING` vs `_PARALLEL_INTERPRETERS` surgery — full
  adversarial gate, touches the parallel mains).
- Outcome B (MKL resists) → vendor PFEAST 4.0 into `OTHER/FEAST/` behind a
  CMake gate (ARPACK precedent); pin MPI-int/threading ABI in
  [[Ladruno_internal]]'s compilation journal.

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
