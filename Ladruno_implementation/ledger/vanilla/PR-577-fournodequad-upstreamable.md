---
wp: PR-577
title: "577 -- upstreamable-table row(s)"
pr: "#577"
files: ["`SRC/element/fourNodeQuad/{FourNodeQuad,FourNodeQuad3d,NineNodeQuad,EightNodeQuad,SixNodeTri}.cpp` + `SRC/element/triangle/Tri31.cpp`"]
table: "upstreamable"
legacy_seq: [322]
---
| `SRC/element/fourNodeQuad/{FourNodeQuad,FourNodeQuad3d,NineNodeQuad,EightNodeQuad,SixNodeTri}.cpp` + `SRC/element/triangle/Tri31.cpp` | `// Ladruno` **upstream bug fix**: these continuum elements gate the mass matrix on `if (rho == 0) use material density; else use element rho`, but `sendSelf`/`recvSelf` never serialized the element-level `rho`, and the no-arg (broker) ctor left it UNINITIALIZED (only `FourNodeQuad3d`'s ctor already zeroed it). On `database`/`restore` (and OpenSeesMP send) the element `rho` was garbage heap → garbage≠0 hijacked the mass → NON-deterministic, indefinite M (negative generalized eigenvalues), diverging transient restarts. Fix (A+B, strictly additive): append `rho` to the send/recv data Vector (+1 slot) and add `rho(0.0)` to each blank ctor. Committed nodal disp round-trips fine, so it was invisible to a plain nodeDisp `database_roundtrip`; caught by `tests/test_quad_tri_rho_db_restart.py` (continues the transient across the restart). Vanilla bug — affects all mainline users doing DB/parallel dynamic restarts; upstream candidate. | [#577](https://github.com/nmorabowen/OpenSees/pull/577) |
