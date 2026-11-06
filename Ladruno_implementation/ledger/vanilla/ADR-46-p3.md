---
wp: ADR-46
title: "ADR46 P3 -- 1 vanilla row(s)"
files: ["`SRC/recorder/NodeRecorder.cpp`"]
table: "main"
legacy_seq: [258]
---
| `SRC/recorder/NodeRecorder.cpp` | `// Ladruno` ADR46 P3: `complexEigenRe<N>` / `complexEigenIm<N>` response channels — parse branches checked BEFORE the `eigen` strncmp; record branches in dataFlag bands 200000/210000+mode inserted BEFORE the open-ended `>=3000` sensitivity catch-all; null-safe (unset complex storage records 0.0). Magnitude/phase = post-processing of the pair. | ADR46 P3 |
