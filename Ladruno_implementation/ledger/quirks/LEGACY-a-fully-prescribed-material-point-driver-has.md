---
wp: LEGACY
title: "A fully-prescribed material-point driver has ZERO free equations — so it cannot see a wrong tangent AT ALL"
legacy_seq: 304
---
### A fully-prescribed material-point driver has ZERO free equations — so it cannot see a wrong tangent AT ALL

- **Bites:** every single-element "material point" test in `tests/` that drives
  all 24 DOFs of a unit cube with `sp` constraints (the `lat=(t,v)` flavour of
  `test_asdplastic_mctc`'s driver, and anything copied from it). It is a great
  way to exercise a constitutive law — the strain path is exact, there is no
  global limit point, `nodeDisp` matches the target bit for bit — and that is
  precisely why it is a trap: with every DOF prescribed the global system has
  **no equations**, the Newton loop converges in 1 iteration by construction,
  and **the tangent the material hands the assembler is never used for
  anything**. A material can return the elastic matrix, a blend, or garbage and
  the whole battery still passes.
- **This is how ADR-84 P0 shipped an 87%-wrong corner tangent** past a battery
  that included a finite-difference tangent test: `test_tangent_fd` checks the
  tangent against differences of the material's own response, which a
  self-consistent-but-wrong operator passes, and nothing else in the module
  could observe the tangent at all.
- **Workaround/status:** if a test is meant to gate TANGENT quality (as opposed
  to stress-path correctness), leave some DOFs free so the global Newton has
  real work, and gate on `ops.testIter()`. `test_adr84_p3_confined_corner.py`
  leaves the z-faces unprescribed (`sigma_zz = 0`, 4 free equations) for exactly
  this reason, and the iteration counts then separate the tangent operators
  cleanly (Continuum 2/step, Secant 7, Elastic 9). Note `ops.printA('-ret')`
  returns an EMPTY list after `analyze()` on this driver, so it is not an
  alternative route to the assembled matrix.
