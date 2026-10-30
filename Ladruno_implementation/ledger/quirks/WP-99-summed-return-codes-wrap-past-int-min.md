---
wp: WP-99
title: "Summed setTrialStrain codes wrap past INT_MIN: 64,907 refusing points turn a refusal into an accepted step"
date: 2026-10-30
---
### Summing `setTrialStrain` codes wraps past INT_MIN: 64,907 refusing points turn a refusal into an accepted step

- **Bites:** a model of FORWARD hosts (`Tri31`, `FourNodeQuad`, `LadrunoQuad`;
  the roster in `WP-99-element-refusal-roster-who-acts-on-a-material-s.md`)
  whose material returns `LADRUNO_MATERIAL_REFUSED`
  (-33086) at every point: with about 16k fully refusing quads the step is
  COMMITTED, not cut, and nothing is printed beyond the per-point material
  warnings. Found in the Opus review of #1038 (apeGmsh #2619).
- **Why:** the hosts summed the codes (`ret += setTrialStrain(eps)`) and
  `Domain::update` summed the element codes (`ok += theEle->update()`), all in
  `int`. 64,907 x -33086 < INT_MIN, so the sum wraps to a POSITIVE value, and
  every integrator refuses only `< 0` (`LoadControl.cpp:149`,
  `DisplacementControl.cpp:289`, `Newmark.cpp:498`, `Linear.cpp:145`). A bare
  -1 needs 2^31 points to wrap, which is why upstream never saw it; the far-
  from-zero sentinel (ADR-86b) made the sum reachable.
- **Fix/status:** fold, never sum: keep the most negative code seen (a
  positive one only when nothing negative was seen). Done in the three hosts
  and in both `Domain::update` loops (`ladrunoKeepWorstCode`, Domain.cpp),
  gated by `tests/test_domain_update_refusal_saturates.py` (69,696 refusing
  points on a 132 x 132 mesh, refused; the pre-fix engine returns 0). Still
  summing, outside that slice: `PartitionedDomain::update` (`res +=
  theSub->update()`), `NineNodeQuad`, `SixNodeTri`, `LadrunoCST`, `LadrunoLST`,
  `BezierTri6`, `BezierTet10`, `TenNodeTetrahedron`, `LadrunoBrick20`,
  `LadrunoUP`, `FourNodeQuadUP`, `BBarFourNodeQuadUP`, `NDFiberSection2d/3d`
  (`res +=`). When a host accumulates return codes, write the fold, not `+=`.
  Learned 2026-10-30 (apeGmsh #2619).
