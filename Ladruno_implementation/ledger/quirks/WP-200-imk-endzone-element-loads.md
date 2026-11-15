---
wp: WP-200
title: "LadrunoIMKBeam element loads: dropped before WP-200, and wrong under -endZones until the end-zone load fix"
date: 2026-11-15
---
### `LadrunoIMKBeam` under `eleLoad` or a `LadrunoMembrane` gives wrong support/joint moments with `-endZones`, or carries no element load at all on an old build

- **Bites (old builds):** before WP-200 (34df942c7) `LadrunoIMKBeam(2d)::addLoad` refused
  every element load ("element loads not supported in v1"). `eleLoad` was accepted at
  definition and contributed nothing; a `LadrunoMembrane` whose `-frames` were IMK beams
  lost its slab load without a word (it pushes `Beam3dPartialUniformLoad` /
  `Beam2dPartialUniformLoad` to its frames). Seen by YOO on the Windows build 35dff42ab.
- **Bites (WP-200 to the fix):** with `-endZones dI dJ` and any element load, three errors.
  (1) The face reactions p0 of the clear span were never moved to the nodes with their lever
  arm, so the node moments were short by p0·dI (p0·dJ). Displacements at fixed supports
  still matched, which hid it, but the support moments were wrong (81.7 instead of 116.7 kN·m
  for w = 20 kN/m, L = 8 m, d = 0.5 m) and a free node, such as a beam-column joint, got a
  spurious unbalanced moment. (2) Load lying on a rigid zone was dropped (ΣV = w·Lc,
  not w·L). (3) Partial and point positions were read as fractions of Lc instead of L, so
  the loads `LadrunoMembrane` pushes shrank by Lc / L. Also `hingeMoment` / `hingeSummary` /
  `hingeBiaxialI|J` reported the node basic moments, not the face moments the hinges carry.
- **Why:** the hinge kernel works on the clear span and returns q_node = Tᵀ q_face. Tᵀ adds
  the lever arm of the shear that the end MOMENTS make, but the load's own face reactions
  sat in p0, which the transform adds as node shears only. The element-guide item "FB-3"
  described a fix (`ladrunoIMKClearSpanBeamUniformReactions2d/3d`) that was never written.
- **Fix:** `addLoadEndZones` (both classes) reads positions on the full L and splits each
  load into its i-zone, clear and j-zone parts: a zone part goes to its node as a force
  and a moment; the clear part gives q0 on Lc and face reactions, whose lever arm adds to
  the new node-moment state `dM0` (zeroed in `zeroLoad`, serialized). Axial load stays on
  the full length: the element's axial stiffness is EA/L over L, the zones are rigid in
  bending only. `effectiveBasic()`
  folds `dM0` into the basic forces handed to the transform and into `localForce` / the
  end-section responses; the hinge responses report face moments. The no-zone path is
  unchanged. Oracle: `tests/test_ladruno_imk_endzone_loads.py` (an `elasticBeamColumn`
  member with 1e6-stiff end pieces, same load, 3-D and 2-D).
- **Workaround on an unfixed build:** model the rigid zones as separate stiff elements (or
  `geomTransf -jntOffset` with the hinges at the node) and put the load on the clear-span
  element only.
- **Related trap:** `-pmBal P0c P0t Pbal muBal` takes all four as POSITIVE numbers; `P0t`
  is the tension capacity as a magnitude. A signed `P0t` (e.g. `-2000`) is refused; the
  message now says which rule failed. Learned 2026-11-15 (YOO fork request).
