---
wp: LEGACY
title: "B3 (P2b-2c): the LadrunoContact handler does NOT enforce non-zero SP (imposed displacement) — NOW WARNS"
legacy_seq: 127
---
### B3 (P2b-2c): the LadrunoContact handler does NOT enforce non-zero SP (imposed displacement) — NOW WARNS
- **Bites:** any model that drives a `LadrunoContact` analysis by imposed nodal displacement
  (`ops.sp(node, dof, value)` with a non-zero value in a `pattern`) — e.g. a displacement-controlled
  indenter. **Verified:** a free-x node with `sp(node,1,0.05)` ends at `ux=0.000` under
  `constraints LadrunoContact`, vs `ux=0.050` under `constraints Transformation`. The contact handler is
  Plain-style — it REMOVES the constrained DOF (ID=-1) but never applies the non-homogeneous value (only
  Transformation/Penalty/Lagrange do). So imposed-displacement DRIVING silently does nothing.
- **NOT a regression:** stock `PlainHandler` does the IDENTICAL thing (`ux=0.000`) — it is the documented
  PlainHandler limitation. The contact handler MUST be Plain-style (it injects the contact FE adapters in
  `handle()`, which Transformation/Penalty can't do), so it inherits this. **Fix shipped:** the handler now
  emits the same non-homogeneous-SP WARNING as `PlainHandler` (it was SILENT before — worse than stock
  Plain), pointing at `DisplacementControl`. So a displacement-driven contact model fails LOUD, not zero.
- **How to drive instead:** external LOAD (`ops.load`, `LoadControl`), the `DisplacementControl` integrator
  on a FREE DOF (it augments the system with the load factor — works with the Plain/contact handler, no SP
  needed), or a fixed-geometry incompatibility (pre-set node positions, like `block_on_block`'s 1e-8
  pre-penetration). FULL imposed-SP support would need a Transformation-style contact handler (no current
  consumer; deferred). Found while building the Hertz benchmark (capstone B3 gate 2).
