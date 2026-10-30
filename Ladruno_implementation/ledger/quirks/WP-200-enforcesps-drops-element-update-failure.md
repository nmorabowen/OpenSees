---
wp: WP-200
title: "TransformationConstraintHandler::enforceSPs drops element update() failures; an early-return element then commits the failed iterate"
date: 2026-10-29
---
### Static step with prescribed displacements "converges" in one iteration with Norm 0, the element force unchanged and the section states garbage

- **Bites:** with `constraints Transformation` and an `sp` that moves (a Path or
  Linear series on the drive DOF), a force-based element whose update fails when
  the SP values are applied still converges: `NormDispIncr` prints iteration 1,
  Norm 0, Norm R 0. The element keeps its old force, its sections keep the failed
  iterate, and `commitState` commits those section states. `eleResponse(...,
  'section', i, 'deformation')` then reports strains of 1e4-1e10, also after a
  later failed analyze (the committed state is the garbage one).
- **Why:** `TransformationConstraintHandler::enforceSPs`
  (`SRC/analysis/handler/TransformationConstraintHandler.cpp`, the loop at the
  end) calls `FE_Element::updateElement()` for every element and ignores the
  return code. After a failed update the element force is the old one, so the
  residual is zero and the first Newton correction is exactly zero. The force-based
  `update()` returns early when `|getBasicIncrDeltaDisp()| <= DBL_EPSILON`
  (vanilla `ForceBeamColumn2d::update` too), so the second update reports success
  without solving. `Domain::commit` drops element `commitState` return codes as
  well.
- **Fix/status:** fixed 2026-10-29 in `LadrunoForceBeamColumn2d/3d`: a failed
  update sets `trialFailed`, which disables both early returns (the next update
  re-solves and fails again, so the step fails), puts the sections back on the
  last solved trial state `vs`, and makes `commitState` refuse through
  `ladrunoNoteCommitRefusal()` (WP-99 F7). Revert and a successful update clear
  it. Upstream `ForceBeamColumn` and the handler are unchanged: any element with
  a "nothing moved" early return is exposed. Regression:
  `test_lfbc_failed_update_is_not_committed` (2D/3D, `-iter 1 1e-300` on an
  elastic jump). `revertToLastCommit` itself restores what `eleResponse` reads,
  as vanilla does (`test_failed_analyze_restores_committed_section_deformation`).
