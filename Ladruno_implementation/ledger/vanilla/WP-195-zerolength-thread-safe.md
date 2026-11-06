---
wp: WP-195
title: "zeroLength + Parallel + Viscous audited and on the WP-192 thread-safe opt-in"
date: 2026-10-14
pr: "#983"
files: ["`SRC/element/zeroLength/ZeroLength.h`", "`SRC/element/zeroLength/ZeroLength.cpp`", "`SRC/material/uniaxial/ParallelMaterial.h`", "`SRC/material/uniaxial/ParallelMaterial.cpp`", "`SRC/material/uniaxial/ViscousMaterial.h`"]
---
| `SRC/element/zeroLength/ZeroLength.h` | `// Ladruno WP-195`: the class-wide `ZeroLengthM*/V*` are `static thread_local`; `ladrunoM()/ladrunoV()` accessors; opt-in `ladrunoThreadSafeUpdate/Form`. Not filed upstream. | [#983](https://github.com/nmorabowen/OpenSees/pull/983) |
| `SRC/element/zeroLength/ZeroLength.cpp` | `// Ladruno WP-195`: the pool definitions `thread_local`; every state-path read of `theMatrix/theVector` goes through `ladrunoM()/ladrunoV()` (this thread's copy, selected by `numDOF` exactly as `setDomain` pairs them); `displaySelf` statics `thread_local`; opt-in = no Damping object and every uniaxial material (incl. the `-doRayleigh 2` damping copies) opted in. Not filed upstream. | [#983](https://github.com/nmorabowen/OpenSees/pull/983) |
| `SRC/material/uniaxial/ParallelMaterial.h` | `// Ladruno WP-195`: opt-in declarations. Not filed upstream. | [#983](https://github.com/nmorabowen/OpenSees/pull/983) |
| `SRC/material/uniaxial/ParallelMaterial.cpp` | `// Ladruno WP-195`: opt-in = every child material's answer (the state path touches only the instance and its children). Not filed upstream. | [#983](https://github.com/nmorabowen/OpenSees/pull/983) |
| `SRC/material/uniaxial/ViscousMaterial.h` | `// Ladruno WP-195`: opt-in true (no shared scratch on the state path). Not filed upstream. | [#983](https://github.com/nmorabowen/OpenSees/pull/983) |
