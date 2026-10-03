---
wp: LEGACY
title: "rigidLink beam across an ndf mismatch (6-DOF master, 3-DOF slave) WARNS and silently adds NOTHING — the model solves DISCONNECTED"
legacy_seq: 159
---
### `rigidLink beam` across an ndf mismatch (6-DOF master, 3-DOF slave) WARNS and silently adds NOTHING — the model solves DISCONNECTED

- **Bites:** the natural first attempt at a shell↔solid seam — `rigidLink('beam',
  shellNode, solidNode)` per matched node column — prints
  `RigidBeam - mismatch in numDOF between constrained Node ... and Retained node ...`
  to opserr and returns WITHOUT adding any constraint and WITHOUT raising. The
  analysis then runs to completion with the two meshes disconnected: under a
  moment the shell rides freely as a hinge and the solid stays unloaded — easy
  to read as "the connection is just flexible".
- **Why:** vanilla `RigidBeam` builds a square ndf×ndf constraint matrix and
  hard-requires matching DOF counts; the failure path is a warning + early
  return, and the interpreter command does not convert it into an error.
- **Workaround/status (measured 2026-07-07, ADR-66 G10):** the seam recipe is
  `LadrunoTie -shellSolid` (exact plane-section rows, works on non-matched
  grids, dt_cr-neutral under the explicit projector). The disconnect trap is
  pinned by `tests/test_ladrunoSolidShell_seam.py::test_rigidlink_cross_ndf_silently_disconnects`
  so a future rigidLink behavior change surfaces loudly.
