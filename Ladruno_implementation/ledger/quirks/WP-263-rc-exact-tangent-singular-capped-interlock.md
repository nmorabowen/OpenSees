---
wp: WP-263
title: "The exact RC tangent of a capped interlock crack is singular, and ASDShellQ4's EAS update blows up on it"
date: 2026-11-28
status: "fixed-by WP-263"
pr: "#1181"
---
### ASDShellQ4 + LayeredShell RC: "sigma_33 condensation did not converge" with -dualTangent, membrane strains ~1e26 after the crack caps

- **Bites:** an exact (consistent, dual, or FD) tangent of a material with a perfectly plastic cap,
  used inside an element that inverts a block of the material tangent. First seen with
  `LadrunoRCConcrete -interlock -cyclic -dualTangent` in `ASDShellQ4` (EAS on by default) +
  `LayeredShell`, `tests/test_ladrunoRCConcrete_dual.py::test_platefiber_view_runs_cyclic_interlock_with_dual`.
  It ran in tension, then 2 shear steps after the crack capped (load factor 1.16):
  - global iteration 2 had `Norm deltaR 9.4e18`;
  - the GP membrane strains reached `eyy, gxy ~ 1e26`;
  - the material then reported "sigma_33 condensation did not converge (|s33|=4.36, then 37.0)";
  - `analyze` returned -3.

  The same run with `-noeas`, and the default tangent, pass.
- **Why:** on a capped crack `tau_ci = +/- v_ci,max` is pinned, so the exact
  `d tau_ci / d gamma = 0`. The in-plane tangent then has a zero crack-shear stiffness: the
  condensed plane-stress `dsig_xy/dgamma_xy` is exactly 0 on that path (`rc_dual_gpp.cpp`
  `ps_shell_test_path`, SH0 = 0). EAS updates the enhanced modes with `KQQ^-1`,
  `KQQ = int G^T D G`, which is singular. The material's sigma_33 Newton was NOT the cause: at
  the material point, on the prescribed strains, it converges in both tangent modes. Its failure
  was downstream of the absurd strains. This is research downside D1 (an exact tangent can be
  singular or indefinite where the secant is SPD), found in an ELEMENT-local solve. It bears on
  decision D2: a production exact tangent needs the same conditioning floors the secant carries.
- **Fix (WP-263):** the Dual keeps mode 0's documented `-betaSrMin` floor on a capped crack,
  derivative only, so the stress is unchanged:
  - on the cyclic cap it adds `betaSrMin G_int m_eps (x) m_full`, equal to mode 0's to 1e-14;
  - on the 2a clip it adds `betaSrMin m_eps (x) d(tau_sm)`.

  `-betaSrMin 0` gives the exact map; the oracle's FD legs run with it.
  `tests/test_ladrunoRCConcrete_dual_cpp.py::test_capped_interlock_floor_and_conditioning` gates
  both:
  - the exact capped shear stiffness is singular;
  - the production one is `betaSrMin G`.
- **Still open:** other exact-tangent singularities or indefiniteness inside element-local
  solves, for example strain softening at a GP under EAS. A floor covers only the cap.
