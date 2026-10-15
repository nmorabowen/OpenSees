---
wp: WP-195
title: "stdBrick element tangent cache keyed on the raw bits of the material tangents and the geometry"
date: 2026-10-14
pr: "TBD"
files: ["`SRC/element/brick/Brick.h`", "`SRC/element/brick/Brick.cpp`"]
---
| `SRC/element/brick/Brick.h` | `// Ladruno WP-195`: per-element tangent cache members (`ladrunoKcache` buffer, state, capacity, consecutive-miss count, on/off, hit/build counters), `ladrunoInvalidateTangentCache()`, `formResidAndTangentCached()`. Not filed upstream. | TBD |
| `SRC/element/brick/Brick.cpp` | `// Ladruno WP-195`: `getTangentStiff()` goes through `formResidAndTangentCached()` (no Damping object, cache on, < 16 consecutive misses): the residual by `formResidAndTangent(0)` itself, the eight material tangents' raw bits and `xl`'s bits compared with the ones that built the cached 24x24 stiffness; equal => the cached matrix, else the stiffness rebuilt by the same loop nest as `formResidAndTangent(1)` and re-cached (bit-identical by construction). Invalidated by `revertToStart`, `setDomain`, `setDamping`, `setParameter`, `updateParameter`, `recvSelf`. Switch `LADRUNO_BRICK_TANGENT_CACHE=0` read at construction; element response `ladrunoTangentCache` = [enabled, state, hits, builds, misses]. Not filed upstream. | TBD |
