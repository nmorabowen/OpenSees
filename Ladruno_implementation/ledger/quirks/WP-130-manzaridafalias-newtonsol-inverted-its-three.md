---
wp: WP-130
title: "ManzariDafalias::NewtonSol inverted its three 6x6 blocks on Matrix::matrixWork, the process-wide scratch; NewtonIter kept function-scope STATIC work arrays (WP…"
legacy_seq: 483
---
### `ManzariDafalias::NewtonSol` inverted its three 6x6 blocks on `Matrix::matrixWork`, the process-wide scratch; `NewtonIter` kept function-scope STATIC work arrays (WP-130)
- **Bites:** a threaded state-determination (WP-107 / F19) on IntScheme 2 shares one DGETRF/DGETRI
  work and pivot array across every Gauss point (a silently WRONG inverse), and anything inverting
  a matrix larger than 400 doubles REALLOCATES it under the others (use-after-free). The
  `static Vector sol/R/R2/dX/norms/aux` and `static Matrix jaco/jInv` in `NewtonIter` are also
  shared and sized by whichever call ran first -- but `NewtonIter` has NO caller (IntScheme 2 runs
  `NewtonIter2 -> NewtonSol`), so the statics were never the live hazard; WP-131's inventory (D6/D9)
  says so. The trap is fixing the dead statics and calling IntScheme 2 re-entrant.
- **Fixed (WP-130, #868):** `NewtonSol` calls a file-static `ladrunoInvertLocal` (same LAPACK pair,
  pivot/work on the stack, lwork 400) -- bit-identical to `Matrix::Invert` for n = 6 (unblocked
  DGETRI), pinned on seven IntScheme-2 decks including a free-DOF Newton deck whose iteration counts
  are compared too. `NewtonIter`'s statics are locals. `ladrunoThreadSafeUpdate()` still answers
  false: re-entrancy is MEASURED in WP-131, not argued here. Other `Matrix::Invert/Solve` callers
  in the file (the dead `NewtonSol2`/`NewtonSol_negP`, `NewtonIter`) still use the shared scratch.
