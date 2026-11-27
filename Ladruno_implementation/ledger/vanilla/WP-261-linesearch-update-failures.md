---
wp: WP-261
title: "Line searches no longer swallow a failed Integrator update"
pr: "#1168"
date: 2026-11-27
files: ["SRC/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.cpp", "SRC/analysis/algorithm/equiSolnAlgo/BisectionLineSearch.cpp", "SRC/analysis/algorithm/equiSolnAlgo/RegulaFalsiLineSearch.cpp"]
---
| `SRC/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.cpp` | WP-261: a negative `theLineSearch->search()` return (an Integrator `update`/`formUnbalance` failure inside the search, e.g. an element whose local solve did not converge) now fails the step with -4. Upstream ignored it and ran the convergence test on the failed trial state. | #1168 |
| `SRC/analysis/algorithm/equiSolnAlgo/BisectionLineSearch.cpp` | WP-261: the no-bracket reset now checks `update`/`formUnbalance` and returns -1, like every other update in the search. Upstream ignored both. | #1168 |
| `SRC/analysis/algorithm/equiSolnAlgo/RegulaFalsiLineSearch.cpp` | WP-261: as BisectionLineSearch. | #1168 |
