---
wp: LEGACY
title: "quad8 mortar dual masses are NEGATIVE at corners (−A/12) by construction — sign-aware guards, do not \"fix\""
legacy_seq: 262
---
### quad8 mortar dual masses are NEGATIVE at corners (−A/12) by construction — sign-aware guards, do not "fix"
- **Bites:** the serendipity quad8 corner shape functions integrate to **−A/12** over the facet (mids +A/3; the total is still A), so the ADR-62 P2.1 dual condensation `Aᵉ=diag(∫N)(Dᵉ)⁻¹` legitimately produces NEGATIVE diagonal dual masses at every quad8 corner node. The sign cancels exactly in `P = Mdual/Ddual` (partition of unity `P·1=1` holds algebraically), but any guard, recorder, or future reader that assumes `Ddual > 0` — the shipped `Ddual[I] <= 1e-300` "uncovered node" refusal did exactly this — false-refuses every valid quad8 tie. Same trap wherever a rowsum `∫N_I ≥ 0` assumption hides: the pre-ADR-78 coverage ratio (`cover/fullCov ≥ 1−1e-3` FLIPS for negative rowsums) and the `|gap|/cover` normalization.
- **Rule:** for serendipity bases, node-wise rowsum measures are SIGNED; guards must be sign-free (areas, L1 integrals) or sign-aware (`|Ddual|`). ADR-78 D3 unified the mortar-tie guards on per-facet AREA coverage + per-facet `∫|g_N|/area`.
- **Workaround/status:** ✅ shipped that way (ADR-78). *2026-08-04.*
