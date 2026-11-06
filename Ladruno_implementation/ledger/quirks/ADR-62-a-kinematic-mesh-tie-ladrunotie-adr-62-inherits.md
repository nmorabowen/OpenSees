---
wp: ADR-62
title: "A kinematic mesh-tie (LadrunoTie, ADR-62) inherits ALL of the projection handler's requirements; the GENERATOR must refuse-and-hand-off where they aren't met"
legacy_seq: 138
---
## A kinematic mesh-tie (LadrunoTie, ADR-62) inherits ALL of the projection handler's requirements; the GENERATOR must refuse-and-hand-off where they aren't met

`LadrunoTie` emits ordinary `EQ_Constraint`s (`u_s = Σ N_i u_{m,i}`) for the **shipped**
`LadrunoProjectionHandler` (ADR-30) to enforce. So the tie is only as usable as that handler:
`system Diagonal` only, explicit only, no MP-chains / double-constraints, lumped mass on every tied
DOF, ICs on the constraint manifold, partition-interior. The generator front-loads the detectable
violations as **named refusals at model-build** (clearer than a mid-analysis singular solve):
- **node-disjoint + one-facet-per-slave (BLOCKER-1).** A node that is both a slave and a master
  facet node, or a slave listed twice, is refused — those are exactly the MP-chain / double-constraint
  topologies the handler refuses at `handle()`. Collocation (one facet per slave) guarantees disjoint
  slave sets, so no slave appears in two constraints.
- **massed tied DOF (BLOCKER-2).** The projection keeps slave DOFs in the equation set (does NOT
  eliminate them), so a massless tied DOF makes `(LᵀML)` singular. **Subtlety: at the generator's
  model-build time the nodal `mass()` is usually still zero — solid nodes get their mass from element
  `-rho`, ASSEMBLED later.** So the generator's mass check must read `Element::getMass()` (formable
  from rho+geometry at build time, the same call `consistentMassGuard` makes at `handle()`) and mark a
  node "massed" if it has nodal mass OR belongs to an element with a nonzero mass diagonal. Consequence:
  **define `LadrunoTie` AFTER the elements/masses** — emitting it before the mass-bearing elements exist
  trips a false BLOCKER-2 refusal.
- **conforming-at-interface ICs (OQ-3).** The displacement tie is trivially on-manifold for a fresh
  model (all `u=0`), so the IC concern is purely geometric: the slave's reference coords must lie on
  the master surface or the projection/weights are meaningless. The generator refuses a slave whose
  closest-point projection lands farther than `-tol * facet-size` off the surface (default 1e-6).
  v1 does **not** snap ICs (`-projectICs` was declined) — non-conforming means different mesh
  *resolutions* on a *shared* surface, where slave nodes are already on a master facet.

Also: the generator drops near-zero shape weights (`|N_i| < 1e-12`) before emitting, so a slave that
projects onto a facet corner/edge ties only to the master nodes that actually carry a share (a corner
collocation degenerates to a clean `u_s = u_{m,corner}`, i.e. an `equalDOF`), avoiding spurious group
connectivity in the handler's connected-component grouping.
