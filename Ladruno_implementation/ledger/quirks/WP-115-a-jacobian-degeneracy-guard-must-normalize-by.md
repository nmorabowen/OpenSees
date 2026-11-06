---
wp: WP-115
title: "A Jacobian degeneracy guard must normalize by the LARGEST column — |det|/∏‖cols‖ is blind to axis COLLAPSE, and an absolute |det| threshold is not scale-free (…"
legacy_seq: 473
---
### A Jacobian degeneracy guard must normalize by the LARGEST column — `|det|/∏‖cols‖` is blind to axis COLLAPSE, and an absolute `|det|` threshold is not scale-free (#588, recorded WP-115)
- **Bites:** an element that refuses degenerate geometry (here the `-formulation eas` centroid Jacobian of `LadrunoQuad`/`LadrunoBrick`) keeps running on a collapsed element. `|det J0| / (‖c0‖‖c1‖)` is `|sin θ|`: it catches collinear axes, but a quad squashed onto a line leaves one fp-residue column (~1e-17), `det` and the column product shrink together, the ratio stays ~1, `J0⁻¹ ~ 1e17`, and the enhanced Newton goes NaN with `analyze()` returning 0 (see "`analyze()` returns rc=0 on a NaN-poisoned system"). An absolute `|det|` threshold instead scales as `L^dim` and refuses healthy small elements.
- **Fix:** `|det| / max‖col‖^dim` = `(‖c_min‖/‖c_max‖)·|sin θ|` catches collapse and collinearity, is size-invariant, and still passes a healthy 100:1 element (1e-2 ≫ 1e-10). Fixed in `LadrunoQuad.cpp` / `LadrunoBrick.cpp` by #588.
- **Why it survived:** the guard had ZERO direct tests; both adversarial reviews probed collinearity, never collapse, and reverting to an absolute threshold passed all 13 existing tests. Gate: `tests/test_ladrunoQuad_eas.py`, `tests/test_ladrunoBrick_eas.py` (two-mode refusal battery + 1e-6/1/1e6 scale pins). Any new geometry guard needs the same two-mode refusal test.
