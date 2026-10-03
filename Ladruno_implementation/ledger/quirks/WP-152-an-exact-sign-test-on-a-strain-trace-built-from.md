---
wp: WP-152
title: "An exact sign test on a strain TRACE built from B·u is a coin flip under isochoric deformation (WP-152)"
legacy_seq: 540
---
### An exact sign test on a strain TRACE built from B·u is a coin flip under isochoric deformation (WP-152)
- **Bites:** a material branch on `tr(dε) > 0.0` (or `>= 0`, `< 0`) to tell compression from opening.
  - The element forms dε = B·Δu; under a deformation that is isochoric in exact arithmetic, each Gauss point's trace is round-off (~1e-20 on 1e-4 strains) with either sign.
  - Measured: WP-152's first non-compressing gate (`tr dε > 0.0`) on a homogeneous stdBrick under +3e-4/−3e-4 pure shear: 6 of 8 GPs separated, 2 were held as "compressing" and refused. Under stdBrick (which discards the material's code) the 2 refusals aborted the commit through the WP-99 latch.
- **Rule:** compare a derived trace against a tolerance scaled to the increment (`tr dε > 1e-10·‖dε‖`), and test the branch with a pure-shear increment on EVERY Gauss point of an element, not just GP 1.
- **Workaround/status:** ✅ WP-152 re-review (`LadrunoSANISANDSasME.cpp`, the E1/E2 gate; `tests/test_ladruno_sanisand_tension_cutoff.py::test_E2_separates_under_isochoric_shear`).
