---
wp: LEGACY
title: "UniformExcitation writes to theDof with NO per-node ndf bounds check (mixed-ndf footgun)"
legacy_seq: 72
---
### `UniformExcitation` writes to `theDof` with NO per-node ndf bounds check (mixed-ndf footgun)
- **Why it bites:** `UniformExcitation::applyLoad` calls `theNode->setR(theDof, 0, fact)` for every
  node in the domain, guarding only on `ndm` (`theDof < 1/2/3`), NEVER on the node's actual
  `numberDOF` (`SRC/domain/pattern/UniformExcitation.cpp:303-365`, writes at `:318/:323/:335`).
  In a MIXED-ndf model a single `UniformExcitation` hits ALL nodes, so exciting e.g. `dof 2`
  is fine on the ndf>=3 nodes but writes OUT OF BOUNDS on any ndf=2 (plane-solid) node sharing
  the domain. No warning, no skip. Apply ground motion only to a node set that actually owns
  the dof, or fix the affected nodes out of the excited direction. See
  [[ndf_and_mixed_models_guide]] §6. 2026-06-07.
