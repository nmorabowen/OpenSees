---
wp: WP-240
title: "LadrunoDynamicRelaxation: MPI reductions and ladrunoDR relax wiring"
date: 2026-11-19
pr: "#1130"
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-240` (two lines): the `integrator LadrunoDynamicRelaxation` arm calls `Ladruno::drInstallMPIReductions()` (a no-op outside `_PARALLEL_INTERPRETERS`), so `OpenSeesMP.exe` reduces the kinetic energy and the residual over the processes; the classic `ladrunoDR` bridge passes `theTransientAnalysis` to `Ladruno::drCommand` so that `ladrunoDR relax` can march. Ported from jaabell/OpenSees `3c2f585b5` and `9a429237e` (J. Abell). | [#1130](https://github.com/nmorabowen/OpenSees/pull/1130) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-240` (two lines): the same call in the `LadrunoDynamicRelaxation` arm of `OPS_Integrator` (this is the engine of `openseesmp.pyd`, which never runs `commands.cpp`; upstream has the hook in `commands.cpp` only), and `OPS_LadrunoDRCmd` passes `cmds->getTransientAnalysis()`. | [#1130](https://github.com/nmorabowen/OpenSees/pull/1130) |
