---
wp: LEGACY
title: "tri6 SLAVE facets are structurally incompatible with the dual mortar — corner ∫N = 0, refused by name"
legacy_seq: 263
---
### tri6 SLAVE facets are structurally incompatible with the dual mortar — corner ∫N = 0, refused by name
- **Bites:** on the reference triangle the tri6 CORNER shape functions integrate to exactly ZERO (the three midsides carry the whole area). The dual scaling divides by the per-facet corner rowsum ⇒ division by zero, and no tolerance rescues it — it is structural, not conditioning. The rowsum-based coverage machinery is equally meaningless there. quad8 corners (−A/12 ≠ 0) are fine.
- **Rule:** `LadrunoTie -mortar` refuses `npsS == 6` in BOTH bases (ADR-78 D2, mirroring apeGmsh ADR 0086 v1); tri6 MASTER facets are fully supported. Remedy in the message: swap master/slave, or put quad8/hex20 faces on the slave side. Revisit only if a tet10-interface user materializes.
- **Workaround/status:** ✅ named refusal shipped (ADR-78). *2026-08-04.*
