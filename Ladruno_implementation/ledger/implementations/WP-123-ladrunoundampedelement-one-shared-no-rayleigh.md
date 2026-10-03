---
wp: WP-123
title: "LadrunoUndampedElement — one shared no-Rayleigh base for the four couplings (WP-123)"
pr: "#858"
status: "draft"
section: "table"
legacy_seq: 13
---
| **`LadrunoUndampedElement` — one shared no-Rayleigh base for the four couplings (WP-123)** ([[123_coupling_undamped_base]]) — `LadrunoDistributingCoupling`, `LadrunoKinematicCoupling`, `LadrunoEmbeddedNode`, `LadrunoEmbeddedRebar` each carried the same Rayleigh-ignoring trio + `C0`/`dampF` buffers, and the missing-`getDamp` implicit-transient crash was fixed in them one after another (#219, #220). They now derive from a header-only base (refuses the factors; per-instance zero `getDamp`). The four `getRayleighDampingForces` copies are deleted: not virtual in `Element`, they only shadowed it and were never called. Bit-identical: 7,504/7,504 recorded values (56 series) equal to the pre-change build; 178/178 batteries. New zone_a `tests/test_ladruno_undamped_couplings.py` (damped == undamped bit for bit, `dampingForce` == 0, `dt_cr` independent of `betaK`), mutation-verified. Refactor candidate 1 of WP-120 R3. | refactor (no behaviour change) | 33005 / 33006 / 33011 / 33012 (existing) | `SRC/element/ladrunoEmbeddedRebar/LadrunoUndampedElement.h` (new), the four elements' `.h/.cpp`, `SRC/element/ladrunoEmbeddedRebar/CMakeLists.txt`, `tests/test_ladruno_undamped_couplings.py`, `Ladruno_implementation/wp123_undamped/*`, `ladruno-new-element` guide, LEDGER_quirks | draft | #858 |
