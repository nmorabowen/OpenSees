# T-MPI — parallel shared-node lumped-SMS validation

**Result: PASS, bit-identical.** Lumped selective mass scaling (ADR 36) is
**parallel-correct**. The earlier handoff worry that a partition-boundary node
would get only rank-local ΔM ("shared-node masses desync across ranks") is
**wrong for the lumped path**.

## Why it is correct (architecture)

- A parallel build auto-swaps the explicit `system Diagonal` for a *distributed*
  diagonal SOE: `DistributedDiagonalSOE` (OpenSeesSP, `_PARALLEL_PROCESSING`) or
  `system MPIDiagonal` → `MPIDiagonalSOE` (OpenSeesMP, `_PARALLEL_INTERPRETERS`).
  See `SRC/tcl/commands.cpp:3213`.
- Both solvers **sum the shared-DOF diagonal across ranks at solve time**:
  - `DistributedDiagonalSolver::solve()` — gather→P0→`*vectShared += otherShared`→broadcast.
  - `MPIDiagonalSolver::solve()` — `intersectionsAB()` accumulates `A` (and `B`) on the
    first solve; subsequent solves reuse the factored `A` and reduce only `B`.
- `CentralDifferenceLadruno` obtains M through `formTangent → SOE → solve`, i.e. it reads
  the **reduced** diagonal, NOT raw `Node::getMass()`.
- Each element lives wholly on one rank ⇒ per-element `dt_e`/scale `s` are rank-locally
  correct. SMS injects `(s−1)·m` into each rank's local copy of the node via `setMass`; the
  solve-time reduction then sums all ranks' contributions to the exact physical total.
- `dtTarget` is a user input (identical on all ranks) ⇒ no dt desync.

**Consistent (Olovsson) is the opposite:** `consistentPCG`/`consistentMatVec`
(`SRC/analysis/integrator/LadrunoMassScaling.h`) run a matrix-free Jacobi-PCG with
**rank-local** inner products (`res ^ z`, `p ^ Ap`) and no shared-DOF `M̄` matvec
exchange — NOT parallel-safe. A real distributed CG is required.

## The test

`sms_bar_mp.py` — a 1-D fixed-free truss bar, 20 elements, with a fine 4-element
zone (elements 9..12, `h=0.1` vs bulk `h=1.0`) straddling the central node 11.
Under a 2-rank manual partition the bar splits at node 11 (SHARED). Elements 10
(rank 0) and 11 (rank 1) are both below `dtTarget=0.004`, so `CentralDifferenceSMS`
injects fictitious mass into node 11 **from both ranks**. Correct behaviour needs
the MPIDiagonal solver to sum those across the partition.

- `np=1`: whole bar on one rank → trusted full-model reference (all injection local).
- `np=2`: split at node 11 → cross-rank ΔM reduction exercised.

If the tip-displacement histories match, the shared-node ΔM was reduced correctly.

## Run

`CentralDifferenceSMS` is NOT in the legacy Tcl `integrator` parser (it lives in the
interpreter/openseespy layer), so the validation uses **OpenSeesPyMP**, not the Tcl
`OpenSeesMP.exe`. The `openseesmp.pyd` is built to `dist/openseesmp/`.

```powershell
# build once (incremental; core objects are already compiled):
cmd /c "Ladruno_scripts\build.bat OpenSeesPyMP"
# run np=1 vs np=2 and compare:
powershell -ExecutionPolicy Bypass -File 'Ladruno_implementation\mass_scaling_mpi\run.ps1'
```

## Measured (this bundle)

```
steps compared : 150
peak |tipDisp| : 3.221400e-03
final ref/par  : 2.3511100000e-03  2.3511100000e-03
max |abs diff| : 0.000e+00
max  rel diff  : 0.000e+00
PASS: np=2 matches np=1 -> shared-node DeltaM reduced correctly
```

The integrators also print a (false-alarm) CFL warning because the per-element dt_cr
eigensolve cannot see the injected nodal mass — the run is in fact stable (bounded,
bit-identical). This is the documented dt_cr-vs-augmentation blind spot (LEDGER_quirks).

## Files

- `sms_bar_mp.py` — the openseespy MP model (handles np=1 and np=2).
- `sms_bar_mp.tcl` — Tcl twin; **does not run** (`CentralDifferenceSMS` absent from the
  legacy Tcl `integrator` parser) — kept to document that gap.
- `compare.py` — reads `tip_np1.out` / `tip_np2.out`, prints PASS/FAIL.
- `run.ps1` — runs np=1 then np=2 via Intel `mpiexec` and compares.
