---
wp: WP-158
title: "Probing a tangent from Python: setNodeDisp without -commit RESETS the node's other DOFs, and printA -ret is the TRANSPOSE (WP-158)"
legacy_seq: 541
---
### Probing a tangent from Python: `setNodeDisp` without `-commit` RESETS the node's other DOFs, and `printA -ret` is the TRANSPOSE (WP-158)
- **Bites:** an FD check `setNodeDisp n 1 ux; setNodeDisp n 2 uy; setNodeDisp n 3 uz; printB` evaluates the residual at
  (committed_x, committed_y, uz), not at (ux, uy, uz). `OPS_setNodeDisp` copies `getDisp()` (the COMMITTED vector), sets
  one component and calls `setTrialDisp`, so each call wipes the trial value of the DOFs set before it. R0.7 lost a day
  to a "12 % tangent error" that was this. Separately, `printA -ret` hands back the `Matrix` buffer, which is
  column-major, so `np.array(...).reshape(n, n)` is `Kᵀ`. A symmetric tangent hides it; a non-symmetric one
  (`-consistanttan`, or the parked FD oracle patch) looks wrong by exactly a transpose.
- **Also:** `setNodeDisp` does not `update()` elements, so a zeroLength or brick keeps its last-`update()` force in
  `printB`. Add their exact linear part analytically (or FD only contact DOFs).
- **Workaround/status:** use `setNodeDisp ... -commit` (a node-level commit; the contact Domain path state is untouched)
  and `reshape(n, n).T`. Probes: `contact_prototypes/probe_adr158_mortar_tangent_fd.py`, `probe_adr158_newton.py`.
