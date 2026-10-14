---
wp: WP-190
title: "MUMPS stats: wall-clock analysis / factorization / solve timings"
date: 2026-10-14
status: "draft"
pr: "#978"
section: "table"
---
| **WP-190 -- MUMPS phase timings in the `-stats` block** -- `system Mumps ... -stats` now also prints `t_analysis s (max)` (only on the print that follows an analysis), `t_factor s (max)` (JOB=2, per factorization) and `t_solve s (max, n)` (JOB=3 wall time and call count since the previous print), each a max over ranks via one `MPI_Allreduce` on the solver communicator. Existing stats lines are byte-identical (`compare_validation.py` parses them). With stats on, the factorization is issued as JOB=2 then JOB=3 (MUMPS defines JOB=5 as exactly these two) so the phases can be timed; stats off keeps the single JOB=5, no `MPI_Wtime`, no reduction. A factorization skipped by WP-188 reuse prints no block (stale-time safe); the solves of the last factor's lifetime are flushed at teardown as `MUMPS solve tail:` (rank 0 local, non-collective). `printStats` now rides `sendSelf`/`recvSelf` so SP subordinates join the collective. | Solver diagnostic | -- (no new tag) | `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}`; `tests/test_wp190_mumps_phase_timings.py` | draft -- PR open, not merged | [#978](https://github.com/nmorabowen/OpenSees/pull/978) |
