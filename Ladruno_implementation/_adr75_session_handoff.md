---
title: ADR-75 session handoff — sparse-direct solver lane (PARDISO desktop / MUMPS cluster)
project: Ladruno
status: handoff — P0/P1/P1c/P2/P2b SHIPPED; P1b + Lane 3 + cluster validation OPEN
owner: nmora
relates: 75_ladruno_sparse_direct_strategy_adr
---

# ADR-75 handoff (session of 2026-07-24/25)

Everything below is **merged to `ladruno`**. Pick up from "What to do next".

## What shipped

| PR | Merge | What |
|---|---|---|
| [#617](https://github.com/nmorabowen/OpenSees/pull/617) | `e8fda0e` | ADR-75 scoping (3 lanes + Abaqus/Kratos/LS-DYNA precedent) |
| [#618](https://github.com/nmorabowen/OpenSees/pull/618) | `57c5b52` | Adversarial-review corrections to the ADR |
| [#619](https://github.com/nmorabowen/OpenSees/pull/619) | `3d36da5` | **P0** trade study (portfolio confirmed) + Lane-B bench harness |
| [#620](https://github.com/nmorabowen/OpenSees/pull/620) | `434b841` | Lane-B baseline measured (UmfPack 22.711 s = P1 gate) |
| [#621](https://github.com/nmorabowen/OpenSees/pull/621) | `c9494db` | Clarify the baseline's MUMPS row is serial-only |
| [#622](https://github.com/nmorabowen/OpenSees/pull/622) | `f5c43db` | **P1a** PARDISO compiles + factorization reuse |
| [#623](https://github.com/nmorabowen/OpenSees/pull/623) | `4144a65` | **P1b** wire `system Pardiso`, threaded MKL — gate PASSED |
| [#624](https://github.com/nmorabowen/OpenSees/pull/624) | `6921f93` | **P1c** scaling study + capability wall |
| [#625](https://github.com/nmorabowen/OpenSees/pull/625) | `2bb348c` | **P2** MUMPS `-BLR` + oneAPI-2026 toolchain fixes |
| [#626](https://github.com/nmorabowen/OpenSees/pull/626) | `434b1a6` | **P2b** `-stats` (INFOG/RINFOG) + BLR memory verdict |

## Headline results (all measured, not estimated)

**Desktop PARDISO (`system Pardiso`) — the big win.**
- Lane B, same binary: **1.71× UmfPack @4 threads**, 1.76× @8, **1.19× even single-threaded**.
- **The win COMPOUNDS with size: 1.61× (11.5k DOF) → 2.15× (26k) → 3.40× (51k)** — UmfPack scales
  ~O(N²), PARDISO ~O(N^1.45).
- **Capability, not just speed: UmfPack OOM'd at 86,490 DOF; PARDISO solved it in 30.4 s and
  136,080 DOF in 68.6 s.** Largest single-machine model raised ≥1.6× in DOF; true ceiling untested.
- **Bit-identical to UmfPack at every thread count** ⇒ no FP-determinism problem in the solver lane.
- **4 threads is the practical default** (1→8 only buys 1.58×; memory-bandwidth-bound).

**MUMPS BLR (`-BLR <eps>`) — shipped, but did NOT pay at any size reachable here.**
- At ~32k DOF/np2: `eps=1e-4` → factor entries −21.8%, BLR flops −45%, **but peak MB/proc only −4.6%**,
  and **3.17× slower**. `eps=1e-9` is strictly worse (+8.4% MB/proc, 1.70× slower).
- **Peak memory is dominated by the active frontal/working space, not stored factors** — so "BLR saves
  memory" is largely false of the number that decides whether a model fits (`INFOG(21)`).
- Nonlinear trap: an inexact solve degrades Newton ⇒ more iterations ⇒ looser eps can be *slower*.
- This bounds the **small end only**. BLR's payoff is governed by **front size** (≈ nested-dissection
  top separator, ~N^(2/3) DOF; front memory ~N^(4/3)) — i.e. by **model size, not rank count**.
  ~20 MB fronts here vs ~1.7 GB at 1M DOF, ~36 GB at 10M DOF.

## Build state (IMPORTANT)

- This worktree's `dist/` has **both** rebuilt modules: serial `dist/bin/opensees.pyd` **with PARDISO**
  (477 pardiso symbols, was 0) and `dist/openseesmp/openseesmp.pyd` **with `-BLR`/`-stats`**.
- **The SHARED checkout `C:\Users\nmb\Documents\Github\OpenSees\dist\` still has the OLD binaries** —
  rebuild there before anyone expects `system Pardiso` outside this worktree.
- ⚠ **oneAPI 2026.1 ships NO Fortran** (no `ifx.exe`, no `ifconsol.lib`). An update re-pointed
  `compiler\latest` 2025.1 → 2026.1 mid-session and broke MP builds three ways. **All three fixed** in
  `Ladruno_scripts/setup_env.bat` + `build.bat`; written up as **`BUILD_GOTCHAS.md §8`**. Serial builds
  are unaffected, which is why it hides until an MP build.

## What to do next (ranked)

1. **Cluster validation of BLR — the open question, ~10 min, NO code needed.** Run a production deck
   twice and compare: `system Mumps -ICNTL14 200 -stats` vs `... -BLR 1e-8 -stats`. Compare
   **`INFOG(21)` (MB/proc)**, total wall (must include the whole Newton loop), and the answer vs
   full-rank. Start near `1e-8`. If `INFOG(21)` doesn't move, BLR is not the memory lever at your scale
   either and the next candidate is **MUMPS out-of-core `ICNTL(22)`** (trades memory for I/O, not
   accuracy — not yet exposed by the fork; small addition next to `-BLR`).
2. **P1b — symmetric PARDISO (`mtype ±2`). ← CHOSEN AS THE NEXT SESSION'S TASK.**
   **Self-contained brief: [[_adr75_p1b_brief]] — start there, it needs no re-derivation.** Needs **upper-triangle SOE storage** in
   `PARDISOGenLinSOE` (a real change; `MumpsSOE`'s `matType != 0` branch is the template). Targets
   memory *and* time — attractive now that memory is the measured constraint. **Must be measured, not
   assumed**: `SparseSYM` (a symmetric solver) is **2.10× SLOWER** than unsymmetric UmfPack, so
   symmetric ≠ automatically better.
3. **Lane 3 — threaded element assembly.** The only lever for the ~34% PARDISO cannot touch (Amdahl:
   measured 1.76× vs a predicted ~2.2× ceiling). Deserves **its own sub-ADR**. Order: (a) scope
   `elem.update` (ADR-40b's #1 instrumentation gap); (b) de-`static` element/material scratch; (c)
   explicit diagonal-SOE private-buffer reduction **with an ordered variant from day one**;
   (d) implicit — **freeze-sparsity + atomic-scatter (Kratos), NOT graph coloring** (see ADR-75 §4);
   (e) thread the `update` loop. Gate each stage on a measured >40% element fraction.
4. **Rebuild the shared checkout** so `system Pardiso` is available outside this worktree.
5. **Version-control the perf skill.** `~/.claude/skills/opensees-performance/SKILL.md` was updated with
   all of the above but **is not a git repo**. Memory says its home should be
   `.claude/skills/opensees-performance/` in the fork — copy + commit it.

## Gotchas banked this session (don't rediscover)

- `theSOE` is a **member of `OpenSeesCommands`** (`OpenSeesCommands.h:178`), not a file-scope global —
  a free `OPS_*` factory cannot assign it. The serial `#else` branch of `OPS_MumpsSolver()` still does,
  which is a **third proof serial MUMPS has never compiled** (with missing `libseq` + CMake gating).
- `~LinearSOE()` **deletes the solver**, so by the time a solver destructor runs the SOE has already
  freed `A/B/X/rowStartA/colA`. Never dereference `theSOE` in a solver destructor (cache what you need).
- The `system Mumps` option loop was `> 1` on the premise **every option takes a value**; a bare flag
  (`-stats`) at the end was silently dropped. Now `> 0`.
- `setup_env.bat` must **not** `setlocal` (it exports env to its caller) ⇒ **no delayed expansion
  `!VAR!`**; write those blocks flat with labels.
- CMake caches the **resolved absolute** compiler path; `-DCMAKE_Fortran_COMPILER=ifx` does **not**
  override a stale cached path. `build.bat` now deletes such a cache.
- A BLR run that gives a **bit-identical answer at a loose eps means BLR never engaged** — don't read it
  as a passing test.
- `sweep_p1.sh` originally `mv`'d the locked baseline JSON and destroyed it (twice). Now `cp`.

## Reproduce the benches
```bash
cd Ladruno_files/testbed/perf/phase1
./sweep_p1.sh                                    # UmfPack vs PARDISO, 1/2/4/8 threads
OPS_SCALE_SIZES=15,20,25,30 MKL_NUM_THREADS=4 python3.12 -S laneB_scaling.py
BLR_NX=20 BLR_NY=20 BLR_NZ=24 ../../../../dist/openseesmp/mpiexec.exe -n 2 \
  python3.12 -S mp_blr_smoke.py ../../../../dist/openseesmp
```
Results: `RESULTS_laneB_baseline.md`, `RESULTS_p1_pardiso.md`, `RESULTS_p1c_scaling.md`,
`RESULTS_p2_blr.md`.

## Standing context
The production regime is **huge solid nonlinear models** (not fiber frames). ADR-75 §1 carries the
correction: Lane B is the *primary* lane, the 66% solve fraction is a **floor** that grows with N, and
no gate measured at ≤136k DOF should be treated as production-representative.
