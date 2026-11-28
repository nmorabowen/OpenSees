---
wp: WP-264
title: "A tangent-only change to the Concrete3D kernel moves the replay's R records' STRESS, not just their M lines"
pr: "#1175"
date: 2026-11-28
status: "open"
---
### Concrete3D replay: a change that only touches the damaged tangent still fails R records on O, N, Z and T values (and can flip a record's outcome)

- **Bites:** anyone changing `damagedTangent` (or anything it calls) who expects the record replay
  (`tests/_testbed/concrete3d_replay.cpp`) to show tangent-only diffs.
- **Why:** an R record replays the plane-stress σ₂₂ Newton of `LadrunoConcrete3D::setTrialStrain`, and that Newton
  steps with the tangent's D₂₂ (`s6[2] -= sig[2]/D22`). A 1e-9 change of D₂₂ moves every later iterate. Most R records
  are refusals: a Newton that does not converge, or that has a negative D₂₂, with chaotic iterates. So the last pass's
  stress, State and even the outcome can change by O(1). The same records already differ in outcome across libms
  (libmmt, glibc and mingw-ucrt disagree on seq 380, 788, 1175, 5671 and 5786).
- **Measured (WP-264, dual driver gradients, mingw-ucrt, `record` of the same inputs on the pre- and post-WP
  headers):** the C and F records' O and N lines are bit-identical, and 6 of 14 C records move their D line by
  ≤ 8.0e-9 relative. 7 of 10 R records move O/N/Z/T/W/M, and two change outcome: seq 380 goes from subInfo 0 to 15,
  and seq 1175 from refused at pass 2 to converged at pass 0.
- **Rule:** for a tangent-only kernel change, "stress byte-exact" holds for the C and F records only. Judge the R
  records by their first pass (Z 0: same strain, D₂₂ moved by the tangent change alone). Any kernel edit also changes
  the header blob, so both fixtures are re-recorded (glibc lane `record`, libmmt lane `concrete3d_replay_libmmt.bat
  record`) in the PR that lands it, once the reviewers have judged the fixture diff.
