---
wp: LEGACY
title: "ZeroLengthSection requires -ndf 3 (2D) / -ndf 6 (3D) — silently absent otherwise"
legacy_seq: 55
---
### `ZeroLengthSection` requires `-ndf 3` (2D) / `-ndf 6` (3D) — silently absent otherwise
- **Bites:** building a `zeroLengthSection` in a reduced-DOF model (e.g. an axial
  SDOF on `-ndf 2`) prints *"ZeroLengthSection::setDomain() -- element only works
  for 3 (2d) or 6 (3d) dof per node"* ([ZeroLengthSection.cpp:247]) and then the
  element is **not added** — but `analyze()` still runs, on a model with no spring,
  so the response looks like an undamped/zero-stiffness free body (constant disp,
  no oscillation). Plain `ZeroLength`/`TwoNodeLink` have no such restriction.
- **Why:** ZeroLengthSection maps the full section response set (P, Vy, Mz, …) onto
  the element DOFs and assumes the complete 3-/6-dof node layout.
- **Rule:** use `-ndf 3`/`-ndf 6` and fix the unused DOFs; never `-ndf 2`. Caught
  by `tests/test_spring_damping_claims.py`. Learned 2026-06-02.
