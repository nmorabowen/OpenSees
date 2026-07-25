---
title: ADR-75 P1b work brief — symmetric PARDISO (mtype ±2) + a memory instrument
project: Ladruno
status: ready to start (self-contained; no re-derivation needed)
owner: nmora
relates: 75_ladruno_sparse_direct_strategy_adr, _adr75_session_handoff
---

# P1b — symmetric PARDISO: work brief

**Goal.** Let `system Pardiso` factor symmetric tangents with MKL PARDISO `mtype = ±2` instead of the
current hardcoded `11` (real unsymmetric). Targets **factor memory ~halved** (store one triangle, not
two) and possibly time.

**Why now.** P1c measured **memory** as the binding constraint on huge solid nonlinear models —
UmfPack hit a hard OOM wall at 86,490 DOF. Symmetric storage attacks exactly that, and unlike BLR it
is **exact** (no accuracy loss, no Newton penalty). It is the most promising remaining desktop lever.

## What has to change

1. **`PARDISOGenLinSOE` — upper-triangle CSR storage (the real work).** It currently builds a FULL
   unsymmetric CSR (`setSize` fills `colA`/`rowStartA` over the whole adjacency; `addA` scatters every
   `(i,j)`). PARDISO `mtype ±2` requires **only the upper triangle**, with the diagonal present.
   **Template: `MumpsSOE`'s `matType != 0` branch** already does exactly this pattern —
   `newNNZ = (nnz - size)/2 + size` in `setSize`, and `addA` skipping `row < col`. Mirror it.
   Note MumpsSOE is *row/col*-flipped relative to PARDISO's CSR; keep the triangle consistent with
   PARDISO's expectation (upper, 1-based, diagonal included, column indices ascending per row).
2. **Solver: make `mtype` selectable** (it is already a member, hardcoded to 11 in the ctor).
   - `2` = real **symmetric positive definite**
   - `-2` = real **symmetric indefinite** ← **the safe default for structural tangents**; softening /
     large-deformation tangents are routinely indefinite and `mtype=2` would fail on a negative pivot.
   - `11` = unsymmetric (keep as the default; contact / non-associated flow genuinely need it).
3. **Interpreter option**, e.g. `system Pardiso -sym` (→ −2), `-spd` (→ 2), default unsymmetric.
   ⚠ The `system Pardiso` parser currently takes **no options at all** — adding one means adding a
   parse loop. Copy the `system Mumps` loop shape, but note **its bound is `> 0`, not `> 1`** (see the
   handoff gotcha: the `> 1` premise breaks for bare flags).
4. **Also add a `-stats` for PARDISO** (small, high value — the analogue of the MUMPS work already
   shipped). MKL reports memory in `iparm` after the symbolic/numeric phases:
   - `iparm[14]` peak memory (KB) for the symbolic phase
   - `iparm[15]` permanent memory (KB) from the symbolic phase
   - `iparm[16]` peak memory (KB) for the numeric factorization + solve
   - `iparm[17]` number of nonzeros in the factors (it is currently set to **0** by our P1a change —
     set it back to **−1** *only when stats are requested*, since −1 makes MKL compute it)
   Without this there is **no way to verify the memory claim** — exactly the gap that made the MUMPS
   BLR result unmeasurable until P2b.

## Gates (do NOT skip — symmetric ≠ automatically better)

- **Correctness first:** symmetric result must match the unsymmetric run to ~1e-12 on Lane B.
- **Memory:** `iparm[16]` (peak numeric) must drop materially vs `mtype=11`. **This is the headline
  claim; if it does not drop, P1b has failed regardless of timing.**
- **Time:** re-run `sweep_p1.sh` and `laneB_scaling.py`. **Do not assume a 2× win** — the fork's own
  `SparseSYM` (a symmetric solver) is **2.10× SLOWER** than unsymmetric UmfPack, so implementation
  quality dominates the storage-format advantage.
- **Capability:** re-run the scaling sweep past n=35 (136k DOF) and see whether the solvable ceiling
  moves. That is the result that matters most for production.
- **Refusal path:** a symmetric solver on a genuinely unsymmetric tangent gives silently wrong
  answers. `LadrunoConcrete3D` (non-associated) and the contact tangents are unsymmetric — `-sym`
  must be opt-in, documented as "your tangent must be symmetric", and ideally warn when a
  known-unsymmetric material/element is present.

## Files
- `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSOE.{h,cpp}` (triangle storage) — VANILLA
- `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.{h,cpp}` (mtype, stats) — VANILLA
- `SRC/interpreter/OpenSeesCommands.cpp` (`OPS_PARDISOGenLinSolver` option parsing) — VANILLA
- Template to copy: `SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp` (`matType != 0` branches)
- Ledger row required in `Ladruno_implementation/LEDGER_vanilla_files.md` (same PR).

## Build + measure
```bash
Ladruno_scripts\build.bat OpenSeesPy          # serial target only; ~15 min with a warm cache
cd Ladruno_files/testbed/perf/phase1 && ./sweep_p1.sh
OPS_SCALE_SIZES=25,30,35,40 MKL_NUM_THREADS=4 python3.12 -S laneB_scaling.py
```
Baselines to beat are in `RESULTS_p1_pardiso.md` (17.775 s UmfPack / 10.396 s PARDISO @4 threads on
Lane B) and `RESULTS_p1c_scaling.md` (the 86k-DOF UmfPack OOM wall; PARDISO 30.4 s there, 68.6 s at
136k).
