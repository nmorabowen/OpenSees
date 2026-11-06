---
wp: WP-159
title: "The mortar friction origin gT0 latches at the FIRST Newton iterate that touches, so a shipped S1 \"passes\" by forgiving the first iterate's slip (WP-159)"
legacy_seq: 543
---
### The mortar friction origin `gT0` latches at the FIRST Newton iterate that touches, so a shipped S1 "passes" by forgiving the first iterate's slip (WP-159)
- **Bites:** `addMortarFriction` captures `gT0` (the stick origin) the first time a node evaluates with p < 0,
  inside the Newton loop, and only `revertToLastStep` undoes it. Under `-adjust` every node starts at p = 0
  (open), so on the R3 pile the whole gravity settlement of the first iterate becomes the stick origin. Engage
  the same nodes from the reference instead (`-gapOffset -1e-6`, or the ADR-159 smoothed law, whose
  first iterate has P(0) = S/4 > 0) and the shipped law FAILS the alpha S1 with growing norms (340 kN). The R3 "S1 passes, axial
  stalls at the slip front" picture is partly this artifact: the stick origin depends on which iterate first
  touched, not on the physics.
- **Also:** a probe that FD-checks friction with `printA`/`printB` at a fresh state must engage the nodes first
  (one `printB` at a zero-slip state), or the first `printB` latches `gT0` at the probe state and every slip
  state reads as stick (a 100 % "tangent error" that is the probe).
- **Workaround/status:** recorded; not changed (ADR-159 §5-§6). Compare runs only at the same engagement
  history.
