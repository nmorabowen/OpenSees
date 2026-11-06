---
wp: WP-139
title: "LadrunoBrick -lumped: one mass model for the inertia residual AND the tangent (WP-139)"
pr: "#879"
status: "ready (owner merges)"
section: "table"
legacy_seq: 15
---
| **LadrunoBrick `-lumped`: one mass model for the inertia residual AND the tangent (WP-139)** ([[139_brick_lumped_inertia]]) — WP-124 gap C8. Under `-lumped` the residual inertia was the CONSISTENT mass (inherited from upstream `Brick.cpp:882-890`) while `getMass` — Newton tangent, αM Rayleigh, ground load — was row-sum LUMPED: implicit dynamics integrated a hybrid and Newton converged only linearly. `formInertiaTerms` now adds `M_L(c,c) a(c)` built in the same Gauss-point order as the mass diagonal (LadrunoBrick20's F-1). **Behaviour change by design** for `-lumped` implicit dynamics AND vanilla `CentralDifference` (which reads the residual at a nonzero trial acceleration); `CentralDifferenceLadruno` byte-identical (Azero residual). Fingerprint: 11/707 series differ, all `Brick/lumped` dynamic; everything else byte-identical. Tests: residual == M_L a node by node, Newton ≤ 3 iterations, consistent mass unchanged; mutation rows L1/L2. Upstream `Brick` left alone. | behaviour fix | 33002 (existing) | `SRC/element/ladrunoBrick/LadrunoBrick.cpp`, `tests/test_ladrunoBrick_lumped_inertia.py`, `Ladruno_implementation/wp139_brick_lumped/*`, LEDGER_quirks, `ladruno-new-element` guide | **ready (owner merges)** | [#879](https://github.com/nmorabowen/OpenSees/pull/879) |
