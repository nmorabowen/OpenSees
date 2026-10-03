---
wp: ADR-20
title: "LadrunoStabilizedUnbalance"
pr: "#184"
status: "shipped — Zone-A ST battery (graceful fallback + true-vs-polluted residual)"
section: "table"
legacy_seq: 61
---
| **LadrunoStabilizedUnbalance** — true-equilibrium `NormUnbalance` convergence test (classTag 33000, CONVERGENCE_TEST registry), **ADR-20 §8 follow-up #4**. In the `LadrunoArcLength -stabilize` mode the SOE residual a stock `CTestNormUnbalance` norms is the **f_v-polluted** `λp − f_int − f_v`, so Newton is satisfied while the TRUE static unbalance is still nonzero by `‖f_v‖`. This test norms `‖λp − f_int‖` instead (a **stricter** criterion driving the artificial force below tol ⇒ a genuine equilibrium, not a regularized one). **Zero analysis-core surgery:** `setEquiSolnAlgo` captures the active integrator via `EquiSolnAlgo::getIncrementalIntegratorPtr`; if it is a stabilizing `LadrunoArcLength` it reads the new `getStabilizedTrueResidual()` accessor, else **degrades gracefully** to the stock SOE `‖B‖` (drop-in for `NormUnbalance`). Python `test LadrunoStabilizedUnbalance $tol $maxIter <$printFlag $normType>`. | Convergence test | 33000 | `SRC/convergenceTest/LadrunoStabilizedUnbalance.{cpp,h}`, `SRC/analysis/integrator/LadrunoArcLength.{cpp,h}` (added `getStabilizedTrueResidual` accessor), `tests/test_ladrunoStabilizedUnbalance_test.py` | shipped — Zone-A ST battery (graceful fallback + true-vs-polluted residual) | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
