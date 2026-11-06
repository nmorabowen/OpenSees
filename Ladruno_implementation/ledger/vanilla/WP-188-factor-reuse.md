---
wp: WP-188
title: "KrylovNewton -reuseFactor hooks in AcceleratedNewton and the two algorithm parsers"
date: 2026-10-11
pr: "#975"
files: ["`SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.h`", "`SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/analysis/algorithm/equiSolnAlgo/CMakeLists.txt`", "`CMakeLists.txt`"]
---
| `SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.h` | `// Ladruno WP-188`: `setFactorReuse()`, a `domainChanged()` override, a `LadrunoFactorReuse *theReuse` member (null = off) and `friend class LadrunoFactorReuse`. Upstream status: not filed (candidate after the 4D validation). | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
| `SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp` | `// Ladruno WP-188`: `solveCurrentStep` returns `theReuse->solveCurrentStep(*this)` when the option is on, else runs the vanilla body unchanged; `domainChanged()` invalidates a kept factor (off: the inherited no-op); destructor deletes `theReuse` (prints the final statistics); `Print` appends them. Not filed upstream. | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-188`: `algorithm KrylovNewton` parses `-reuseFactor/-reuseIter/-reuseRate/-reuseReport` (fail closed) and calls `setFactorReuse`. This is the parser `OpenSees.exe`, `OpenSeesSP.exe` **and `OpenSeesMP.exe`** use. Not filed upstream. | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-188`: the same options in `OPS_KrylovNewton` (openseespy, openseesmp). Not filed upstream. | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
| `SRC/analysis/algorithm/equiSolnAlgo/CMakeLists.txt` | `# Ladruno WP-188`: `LadrunoFactorReuse.{cpp,h}` in OPS_Analysis. | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
| `CMakeLists.txt` | `# Ladruno WP-188`: `LadrunoFactorReuseMPI.cpp` added to `OPS_MPI_PER_TARGET_SOURCES` (the reuse decisions' MPI_Allreduce must see each target's `_PARALLEL_*` define). | [#975](https://github.com/nmorabowen/OpenSees/pull/975) |
