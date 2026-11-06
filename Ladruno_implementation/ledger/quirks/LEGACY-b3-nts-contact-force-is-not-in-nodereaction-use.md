---
wp: LEGACY
title: "B3: NTS contact force is NOT in nodeReaction — use the ladrunoContactForce query"
legacy_seq: 131
---
### B3: NTS contact force is NOT in nodeReaction — use the `ladrunoContactForce` query
- **Bites:** reading per-node contact pressure. The NTS contact traction is assembled by an injected
  `LadrunoContactFE` adapter (an FE_Element with no backing Domain Element), so it does NOT contribute to
  `Node::addReactionForce` ⇒ `ops.nodeReaction(slave, 3)` returns 0 for the contact force (only real
  elements + nodal loads/inertia accumulate into reactions). **Fix shipped (B3):** the SEGMENT adapter
  reports its `tn = kn·<−gap>₊` into a Domain snapshot (`set/getNtsForce`, cleared each handle in
  `frictionGCBegin`); query `ladrunoContactForce slaveNodeTag` returns the Σ over the node's pairs. Pure
  side-channel (no resid/tang effect).
