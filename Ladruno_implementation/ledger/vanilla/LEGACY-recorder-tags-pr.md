---
wp: LEGACY
title: "recorder-tags PR -- 1 vanilla row(s)"
files: ["`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`"]
table: "main"
legacy_seq: [276]
---
| `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` | `// Ladruno` material-response labeling: `setResponse` now emits `NdMaterialOutput` + `ResponseType` XML for every recognized token (`stress`/`strain`/`pstrain`/`eqpstrain`/`PStress`/`J2Stress`/`VolStrain`/`J2Strain` + internal-variable fallback), so `-E material.<token>` records with real column names (`eqpstrain`, `epsP11..epsP13`, …) instead of the recorder's generic `C1..CN`. Also guards the unknown-token path: `pos < 0 || iv_size <= 0` now returns 0 (a mistyped token like `plasticStrain` — the real spelling is `pstrain` — records nothing instead of a malformed `Vector(-1)` bucket) + `argc < 1` guard. `+#include <cstdio>`. Upstreamable (ASDP is stock, by jaabell/Petracca). NOTE: the reported "empty group" was really a token error (`-E eqpstrain` bare never reaches the material; use `-E material.eqpstrain`); the tri6n/tet10 elements already forward `material`/`integrPoint` requests, no element change needed. | recorder-tags PR |
