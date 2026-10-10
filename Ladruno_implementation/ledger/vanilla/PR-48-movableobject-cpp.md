---
wp: PR-48
title: "48 -- 1 vanilla row(s)"
pr: "#48"
files: ["`SRC/actor/actor/MovableObject.cpp`"]
table: "main"
legacy_seq: [117]
---
| `SRC/actor/actor/MovableObject.cpp` | `// Ladruno P4`: TaggedObject live-component census — `OPS_PROFILE_CENSUS_BORN(classTag)` in both ctor bodies, `OPS_PROFILE_CENSUS_DIED(classTag)` in the dtor (+`#include <profiler/ProfilerMacros.h>`). Seam is `MovableObject` (not `TaggedObject`): `classTag` is a plain member valid through the dtor, whereas `TaggedObject::getClassTag()` is virtual and unsafe in a base ctor/dtor. Runtime-gated on `mem()`; off-path cost is one branch on this hot path. Counts every MovableObject by raw classTag (shared integer space → superset; viewer filters by classTag band). Additive. | [#48](https://github.com/nmorabowen/OpenSees/pull/48) |
