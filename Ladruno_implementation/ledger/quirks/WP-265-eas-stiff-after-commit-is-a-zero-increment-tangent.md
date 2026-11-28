---
wp: WP-265
title: "Reading an element tangent after analyze() can return the zero-increment tangent, not the step's"
pr: "#1179"
date: 2026-11-28
status: "open"
---
### `eleResponse 'stiff'` equals `'stiffInitial'` at a plastic state (LadrunoBrick eas, ShellMITC4): the tangent was re-formed from the NEW committed state

- **Bites:** a K-vs-FD tangent check that drives the probe state with an imposed
  step (`analyze`) and then reads the tangent. `LadrunoBrick -formulation eas`
  answered exactly its elastic `'stiffInitial'` at a J2 state 10x past yield,
  so the check failed by 0.62 while std/bbar/ssp/uri passed at 1e-10.
- **Why:** some elements re-strain their material points on EVERY force or
  tangent call, not only in `update()`. Examples: eas (`formEAStrue` re-solves
  alpha and calls `setTrialStrain` at the 8 GPs), and ShellMITC4 (which strains
  its section from the nodal trial displacements). After `analyze` commits, a
  tangent call strains the GPs again at the same strain, but from the
  just-committed history. That is a zero increment, so the material returns the
  elastic or unloading tangent. Elements that set strains only in `update()`
  (std, bbar, ssp, uri, LadrunoBrick20) return the material's last TRIAL
  tangent, which is the step's consistent one. Inside Newton both behave the
  same, because K is always formed after `update()` at the iterate. Only
  post-commit probes differ.
- **Fix (probing rule, WP-265):** for such elements, put the probe in a TRIAL
  state on the committed path. Use `setNodeDisp(tag, dof, value, '-commit')` on
  every component. Per WP-158, without `-commit` each call wipes the others; it
  commits the node, never the materials. Then read `eleForce` and the tangent.
  Shells answer no `'stiff'`. Read them from `printA('-ret')`, which re-forms
  the tangent first; reshape transposed (column-major); take the node-DOF block
  via `nodeDOFs`; and keep the last SP pattern alive so printA's `formTangent`
  does not see deleted constraints. The opposite trap: ASDShellQ4 does NOT
  re-strain on a force call (setNodeDisp leaves its force unchanged, K_fd = 0),
  so it needs the committed-step rig.
- **Seen with:** `tests/test_continuum_tangent_gate.py` (`_probe_trial`) and
  `tests/test_shellmod_tangent_gate.py`.
