---
wp: WP-240
title: "updateDomain() after setResponse() with unchanged displacements updates every element again"
pr: ""
date: 2026-11-19
---
### An explicit integrator that calls `updateDomain()` after `setResponse()` in `update()` runs every `Element::update()` a second time, with the same displacements

- **Bites:** an explicit integrator whose `newStep()` already sets the new displacements and
  calls `theModel->updateDomain(time, dt)`, and whose `update()` then pushes the final
  velocity and acceleration with `setResponse()` followed by `theModel->updateDomain()`.
  `AnalysisModel::updateDomain()` is `Domain::update()` plus `ConstraintHandler::update()`.
  No handler in the tree overrides the second (the base returns 0), so the call is one more
  `Element::update()` per element. Elements that integrate their materials in `update()`
  (LadrunoBrick, LadrunoQuad, stdBrick) do every material point twice per
  step. `LadrunoDynamicRelaxation` did this from v1 until WP-240; José Abell found it on
  jaabell/OpenSees (`c8b0f3bc4`).
- **Cost, measured 2026-11-19** (`Ladruno_files/testbed/perf/wp240/RESULTS.md`): removing the
  second pass gave 68.6 to 49.1 ms/step on a 16x16x16 LadrunoBrick/LadrunoJ2 block, 52.1 to
  36.0 on stdBrick/J2Plasticity and 23.5 to 17.4 on a 64x64 LadrunoQuad panel, with every
  nodal displacement bit-identical. Abell measured 76 to 13 ms/step with a Mohr-Coulomb
  material.
- **Not a pure no-op for rate-dependent elements.** The second pass saw the full-step velocity
  that `update()` had just pushed, where the first saw the half-step one. An element that
  reads the nodal velocity in `update()` (a truss with a `Viscous` material, dampers) therefore
  committed a state its force in the residual was never computed with. With one pass it
  commits the state of `v_{n+1/2}`. Rate-independent elements are unaffected.
- **Safe to drop because:** `TransformationDOF_Group::setNodeDisp` keeps the trial value of a
  DOF that has an `sp` (it does not zero it), so prescribed displacements enforced by
  `ConstraintHandler::applyLoad()` in `newStep()` survive the later `setResponse()`.
- **Do instead:** in `update()`, call `setResponse()` alone when the displacements have not
  changed since `newStep()`. Before removing such a call from another integrator, check that
  its `update()` really pushes the same displacements (CentralDifferenceLadruno and
  ExplicitBathe were not audited here). Gate: `tests/test_wp240_dr_single_update.py` (S-1
  fails on a two-pass build).
