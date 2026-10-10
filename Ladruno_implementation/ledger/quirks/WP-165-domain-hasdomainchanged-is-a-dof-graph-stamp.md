---
wp: WP-165
title: "Domain::hasDomainChanged() is a DOF-graph stamp, not an analysis-stage marker — key anything \"per stage\" on the node/element set (WP-165)"
date: 2026-10-03
---
### `Domain::hasDomainChanged()` is a DOF-graph stamp, not an analysis-stage marker — key anything "per stage" on the node/element set (WP-165)
- **Bites:** the Ladruno recorder started a new `MODEL_STAGE` (full model copy, every source rebuilt, envelopes
  and energy reset) on every stamp move. The stamp moves for `addSP_Constraint` into a pattern (`sp`,
  `imposedMotion`), `addElementalLoad` (`eleLoad`), `addLoadPattern` of a pattern that already holds SPs, every
  `remove*`, and the ADR-60 contact `-reemit` path, which calls `domainChange()` INSIDE `commit()` before the
  recorder loop — so a `-reemit` run wrote one model copy per re-sort, and stage names drifted between ranks.
  It does NOT move for `addNodalLoad`, `loadConst`, `setTime`.
- **Workaround/status:** WP-165: rebuild only when a fingerprint of the node/element/pressure-constraint set
  (tag AND object address — `remove element 5; element ... 5` is a new object) changed. Anything that caches
  `Element*`/`Response*` across a stamp move must use the same test, not the stamp alone.
