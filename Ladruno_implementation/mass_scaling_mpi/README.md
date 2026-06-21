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

## Consistent (Olovsson) parallel — V5

The lumped result above needs **no** MPI work (the distributed/MPI diagonal solver already
reduces shared-node mass). The **consistent (Olovsson)** path is different: it solves
`M_tilde a = r` by a matrix-free Jacobi-PCG, and the serial `consistentPCG` is rank-local
(local inner products, no shared-DOF `M_bar` exchange) — wrong at partition boundaries. V5
adds a **distributed PCG** (`consistentParPCG` in `LadrunoMassScaling.h`, driven by
`LadrunoConsistentRefine.h`) keyed on one weight `w_i = 1/multiplicity_i`:

- matvec: GLOBAL (replicated) lumped diagonal applied WEIGHTED + off-diagonal `M_bar_e` in
  FULL, then `assembleSharedSum` across ranks ⇒ diagonal collapses to full-once, off-diagonal
  accumulates over all element-owning ranks;
- inner products: same `w` + `globalReduceSum` (all-reduce) ⇒ each shared DOF counted once;
- all CG control scalars are global ⇒ identical iteration count on every rank ⇒ the per-iter
  collectives stay in lockstep (no deadlock). At `np=1` it collapses to the serial PCG.

Because the integrators live in the shared `OpenSeesLIB` (compiled once, linked into serial
AND MP) they cannot reference the MP-only `MPIDiagonalSOE`, so dispatch goes through new
`LinearSOE` base virtuals that `MPIDiagonalSOE` overrides.

```powershell
cmd /c "Ladruno_scripts\build.bat OpenSeesPyMP"   # (also OpenSeesPy for the serial reference)
powershell -ExecutionPolicy Bypass -File 'Ladruno_implementation\mass_scaling_mpi\run_consistent.ps1'
```

Measured: `np=1` vs `np=2` match to recorder precision (`max |abs diff| 0`, final 1.746e-3),
and the **gold cross-check** — MPI `np=2` vs the serial `DiagonalSOE`+`consistentPCG`
reference (`sms_bar_serial_consistent.py`) — also matches to recorder precision. The
consistent final (1.746e-3) differs from the lumped final (2.35e-3), confirming it is
genuinely Olovsson, not a lumped fallback. Serial Zone-A 36/36 stays green.

> Non-vacuity: a stale binary once made BOTH ranks diverge to the same overflow, which a
> naive A/B compare "passed". `compare*.py` now reject non-finite / `|disp|>1.0` output, and
> the serial cross-check rules out a bug shared by `np=1` and `np=2`.

## Files

- `sms_bar_mp.py` — the openseespy MP model, LUMPED `CentralDifferenceSMS` (np=1 and np=2).
- `sms_bar_mp_consistent.py` — MP model, CONSISTENT `CentralDifferenceSMSConsistent`.
- `sms_bar_serial_consistent.py` — serial reference (`opensees.pyd`, `system Diagonal`).
- `sms_bar_mp.tcl` — Tcl twin; **does not run** (`CentralDifferenceSMS` absent from the
  legacy Tcl `integrator` parser) — kept to document that gap.
- `compare.py` / `compare_consistent.py` — read the `tip*_np*.out` files, print PASS/FAIL
  (both reject diverged/non-finite output).
- `run.ps1` / `run_consistent.ps1` — run np=1 then np=2 via Intel `mpiexec` and compare.
