---
wp: WP-239
title: "Material thread races: ASDSteel1D Newton history and IMPL-EX accumulator, shell-section getType latch, thread-static waivers"
date: 2026-11-19
files: ["`SRC/material/uniaxial/ASDSteel1DMaterial.cpp`", "`SRC/material/section/LayeredShellFiberSection.cpp`", "`SRC/material/section/ElasticMembranePlateSection.cpp`", "`SRC/material/uniaxial/ASDConcrete1DMaterial.cpp`", "`SRC/material/nD/ASDConcrete3DMaterial.cpp`", "`SRC/material/section/FiberSection3d.cpp`"]
---
| `SRC/material/uniaxial/ASDSteel1DMaterial.cpp` | `// Ladruno WP-239`: the `-buckling -implex` implicit-correction Newton history (`U_iterative`, `U_iterative_el`, `Rnorm_iterative`) goes from plain function-scope `static` to `static thread_local` (written and read within one call, so serial is bit-identical; shared, two points committing on different WP-192 threads picked each other's iterate). `commitState`'s IMPL-EX error goes through `ladrunoRecordImplexError`: per-thread partials inside a threaded region, folded by `ladrunoFlushImplexError` (the ASDConcrete1D WP-192 pattern this class missed). Not filed upstream. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
| `SRC/material/section/LayeredShellFiberSection.cpp` | `// Ladruno WP-239`: `getType()` drops its process-wide `static bool initialized` latch and fills the `thread_local` `array` on every call (the latch filled only the first calling thread's copy). Not filed upstream. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
| `SRC/material/section/ElasticMembranePlateSection.cpp` | `// Ladruno WP-239`: same `getType()` latch fix. Not filed upstream. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
| `SRC/material/uniaxial/ASDConcrete1DMaterial.cpp` | `// ladruno-lint: thread-ok` on the IMPL-EX partials registry (only touched under its mutex). Comment only. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
| `SRC/material/nD/ASDConcrete3DMaterial.cpp` | `// ladruno-lint: thread-ok` on the IMPL-EX partials registry (only touched under its mutex). Comment only. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
| `SRC/material/section/FiberSection3d.cpp` | `// ladruno-lint: thread-ok` on `getEnergy`'s `fiberArea[10000]` (recorder-only; the WP-196 energy sweep never calls it). Comment only. | [#1129](https://github.com/nmorabowen/OpenSees/pull/1129) |
