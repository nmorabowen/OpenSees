---
wp: PR-317
title: "317 -- 1 vanilla row(s)"
pr: "#317"
files: ["`SRC/domain/node/Node.{h,cpp}`"]
table: "main"
legacy_seq: [149]
---
| `SRC/domain/node/Node.{h,cpp}` | `// Ladruno` (ADR-30 P4): add a lazily-allocated per-node constraint tie-force buffer `projTieForce` (in-class `= 0` init → safe across all 6 ctors without touching their init lists; freed in `~Node`) + `getProjectionTieForce()` / `setProjectionTieForce()`, mirroring the `reaction` slot. CentralDifferenceLadruno scatters `M(a_raw−a_proj)` here each commit so the node-based LadrunoRecorder can record a `constraintTieForce` field. Additive; Node.h size change ⇒ recompile-all, no behavior touched. | [#317](https://github.com/nmorabowen/OpenSees/pull/317) |
