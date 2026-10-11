---
wp: WP-183
title: "modalDamping under constraints Lagrange exits: DOF_Group::getEigenvectors: 0 Node Pointer"
date: 2026-10-11
---
### `modalDamping` + `constraints Lagrange` kills the process with `DOF_Group::getEigenvectors: 0 Node Pointer`
- **Bites:** the first transient step after `modalDamping` when the handler is Lagrange (any MP or SP).
- **Why:** `IncrementalIntegrator::setupModal` gathers eigenvectors from every DOF_Group; a `LagrangeDOF_Group` (the multipliers) has no node, and `DOF_Group::getEigenvectors` `exit(-1)`s on a null node. Pre-existing, independent of WP-183's mass fix.
- **Workaround:** use Transformation, Penalty or Auto with modalDamping. Not changed. (Found while writing WP-183's tests; the Lagrange-modal case there is a documented skip.)
