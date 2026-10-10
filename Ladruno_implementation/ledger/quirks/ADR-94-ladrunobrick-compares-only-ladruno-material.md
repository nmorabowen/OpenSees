---
wp: ADR-94
title: "LadrunoBrick compares ONLY == LADRUNO_MATERIAL_REFUSED — a bare -1 from a material is treated as success (ADR-94 B2)"
legacy_seq: 389
---
### `LadrunoBrick` compares ONLY `== LADRUNO_MATERIAL_REFUSED` — a bare `-1` from a material is treated as success (ADR-94 B2)
- **Bites:** every material failure path that returns a plain `-1` (ASDP: 13 of its 15 failure sites, including the NaN guard and the singular-tangent guard; only two `Backward_Euler` sites return the sentinel). `stdBrick` drops every code. Only `TenNodeTetrahedron` (`success += ...`) propagates both. Measured: Drucker-Prager commits NaN stress with `analyze() == 0` on `LadrunoBrick` while the material prints `NaN!` and returns `-1`.
- **Why:** the ADR-86b review-fix deliberately narrowed the host check to the sentinel so that only "commit guaranteed unchanged" refusals abort; nobody widened the material side to match.
- **Rule:** a material that wants to be heard by `LadrunoBrick` must return the sentinel, not `-1`. Fix direction is widening the material's sites, not loosening the host to `< 0`.
- **Status (wp/94a):** FIXED for ASDPlasticMaterial3D — all 13 bare-`-1` sites now return `LADRUNO_MATERIAL_REFUSED` (the host side is unchanged, still sentinel-only, by design). The rule above still stands for **every other** fork material: `stdBrick` drops all codes, and `LadrunoBrick` hears only the sentinel. A NaN or singular-tangent return from ASDP now aborts the step on `LadrunoBrick` where it used to be committed as success — the intended behaviour change, and the first thing to suspect if an ASDP model that "worked" starts failing steps.
