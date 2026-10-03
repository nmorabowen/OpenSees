---
wp: LEGACY
title: "wp/94d -- upstreamable-table row(s)"
files: ["`SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/HoekBrown_YF.h`"]
table: "upstreamable"
legacy_seq: [521]
---
| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/HoekBrown_YF.h` | `// Ladruno (ADR-94 wp/94d)`: ported jaabell/ASDP `60d9b9b23`'s composite `max(f_shear, f_tension)` Hoek-Brown yield function (arg clamped to 0 before `pow`, continuous at the tensile apex) into `YIELD_FUNCTION`, replacing the discontinuous `if (arg>0) ... else ...` split this tree had live (kept as a historical comment, matching jaabell's own convention). ADR-94 H10a/M4 measured the old split locking a uniaxial-tension path onto `sigci*s` (~2.4x = `mb` times the textbook tensile strength `sigma_t=-s*sigci/mb`) and then stalling (`analyze()==-3`). `CHECK_APEX_REGION`/`APEX_STRESS`/`yf_has_apex` kept LIVE from our pre-port tree (jaabell's carries them as dead comments too; the `ASDPlasticMaterial3D.h` apex call site is dead code on both trees today, so this is cosmetic-only). Compression paths are bit-identical before/after (never cross `arg<=0`). See `Ladruno_implementation/_adr94_hb_drift.md`. | wp/94d |
