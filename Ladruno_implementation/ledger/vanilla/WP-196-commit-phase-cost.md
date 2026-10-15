---
wp: WP-196
title: "Commit-phase profiler scopes; EnergyBalance opt-in predicate on Element and ASDShellQ4"
date: 2026-10-14
pr: "#984"
files: ["`SRC/domain/domain/Domain.cpp`", "`SRC/element/Element.h`", "`SRC/element/shell/ASDShellQ4.h`", "`SRC/element/shell/ASDShellQ4.cpp`"]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno WP-196`: coarse profiler scopes `commit.nodes` / `commit.elements` / `commit.record` and one `rec.<kind>` scope per recorder in `commit()` (one branch each when the profiler is off). Not filed upstream. | [#984](https://github.com/nmorabowen/OpenSees/pull/984) |
| `SRC/element/Element.h` | `// Ladruno WP-196`: new opt-in virtual `ladrunoEnergyRayleighOnly()` (default false): the element's damping is exactly the base Rayleigh combination and its post-commit evaluations are idempotent, so the energy kernel may skip the WP-184 check. Not filed upstream. | [#984](https://github.com/nmorabowen/OpenSees/pull/984) |
| `SRC/element/shell/ASDShellQ4.h` | `// Ladruno WP-196`: declares `ladrunoEnergyRayleighOnly()`. Not filed upstream. | [#984](https://github.com/nmorabowen/OpenSees/pull/984) |
| `SRC/element/shell/ASDShellQ4.cpp` | `// Ladruno WP-196`: `ladrunoEnergyRayleighOnly()` = no Damping object and no AGQI (`-noeas`): the AGQI condensation state keeps the kind of the last evaluation. Not filed upstream. | [#984](https://github.com/nmorabowen/OpenSees/pull/984) |
