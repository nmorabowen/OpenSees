---
wp: WP-155
title: "Mortar contact: one tie around a CLOSED surface silently binds the far side; a per-commit Uzawa makes a LINEAR tie step-count dependent (WP-155, 2026-09-29)"
date: 2026-09-29
legacy_seq: 1
---
### Mortar contact: one tie around a CLOSED surface silently binds the far side; a per-commit Uzawa makes a LINEAR tie step-count dependent (WP-155, 2026-09-29)
- **Bites:** (1) the mortar broad phase is brute force and the kernel clip accepts anti-parallel facets (`|cos|`), so on a closed skin (a pile) each facet is also tied to the ANTIPODAL facets — R1 measured one cylinder tie 3.65× too stiff, with no warning. (2) the shipped default augments λ once per `Domain::commit()`, so even a linear penalty tie gives a different answer for 1 vs 2 vs 5 load steps (exactly the 1-D Uzawa recursion: 2.2e-3 / 2.1e-3 / 2.04e-3 on the split column).
- **Also:** a `-adjust`ed node sits at p = 0 EXACTLY, which the shipped `pr < 0` mask treats as open — the first iterate then has no interface stiffness (a floating pile). WP-155 takes the closed tangent branch there. And a frictionless contact on a near-circular interface has NO torsional stiffness: a test slice must pin the rigid rotation or Newton chatters.
- **Also:** the `-augment` gate must cover EVERY Uzawa in `commit()`: the ADR-57 E6 edge-edge λ_N was missed at first (review #897 finding 1).
- **Status:** `-maxGap d` (pairing guard), `-augment request|never`, `-adjust`/`-gapOffset` — [[155_pile_contact_r05]]. Defaults unchanged.
