---
wp: WP-240
title: "MPIDiagonalSOE keeps the A of its first solve; zeroA and addA are then no-ops"
pr: ""
date: 2026-11-19
---
### `system MPIDiagonal` ignores `zeroA()` / `addA()` after its first solve, so a rebuilt mass never reaches the solve

- **Bites:** an integrator that reassembles the diagonal during a run under OpenSeesMP
  (`LadrunoDynamicRelaxation` rebuilding `M*` at kinetic-energy peaks or with `-recompute`).
  `MPIDiagonalSolver::solve` inverts `A` in place at the first solve and sets `isAfactored`;
  from then on `MPIDiagonalSOE::zeroA` does nothing and `MPIDiagonalSOE::addA` returns 0 at
  once. Nothing is printed. The integrator's own copy of the mass changes, the solve keeps the
  old one, and a kinetic energy or residual computed from the copy no longer matches the march.
  This holds on one process as well as on several.
- **Also:** after a solve, `B` of that system holds the neighbours' contributions at the
  interface DOFs (the solver adds them before dividing), and the `A` of a process is the sum of
  all processes' shares there. A process's own row sums are only its share. So a residual
  computed as `M* a` from the local mass is wrong at an interface; read `getB()`.
- **And:** `OpenSees.exe` (serial) answers `getNP`, so a Tcl deck cannot tell
  that it runs under OpenSeesMP from `[info commands getNP]`. Pass the engine in from the
  driver (`tests/tcl/wp240_dr_mp.tcl`, `WP240_PARALLEL`).
- **Do instead:** ask `theSOE->isDistributedDiagonal()` (true only for `MPIDiagonalSOE`) and
  do not rebuild the mass when it is true. `LadrunoDynamicRelaxation::domainChanged` switches
  off the refresh there and says so once (WP-240). The serial run that matches a parallel one
  is the one with `-noAutoRefresh`.
