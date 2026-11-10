---
wp: WP-210
title: "HS plastic shear strain gamma^p is a PRINCIPAL-strain measure (eps1 - eps2 - eps3), not 2 eps_q: the invariant shortcut hardens the cone twice as fast in extension"
date: 2026-11-09
---
### An invariant (p, q) HS return map that books d(gamma^p) = 2 d(eps_q^p) is exact only in triaxial compression

- **Bites:** in a (p, q) radial return it is natural to take the shear-hardening
  increment as `dgam = 2 (q_tr - q) / (3 G)`, i.e. `2 d(eps_q^p)`. PLAXIS HS (Schanz,
  Vermeer & Bonnier 1999, Eq. 12) defines `gamma^p = eps1^p - eps2^p - eps3^p =
  2 eps1^p - eps_v^p`, with `eps_v^p = -sin(psi_m) gamma^p`. The two agree in TXC with
  `psi_m = 0` only. In TXE the invariant version mobilises the cone with HALF the
  major-axis plastic strain (drained TXE E50 about 20 % stiffer than an independent HS
  code with the cap inactive). Whenever the soil dilates, it also drops the
  `-eps_v^p/3` term. LadrunoHS shipped it until WP-210 round 4.
- **Why the gates missed it:** every closed-form and parity gate ran in TXC (where the
  measures coincide) or compared the kernel with its own Python oracle (same shortcut on
  both sides). Only an independent code in extension exposed it
  (`hardening-soil-validation`, numgeo HS-MN-Bricks incrementalDriver, TXE NC and OCR = 3).
- **Fix:** with deviatoric flow coaxial with the trial deviator, `e1^p = (3/2) k1(x)
  d(eps_q^p)` with `k1 = s1/q` (2/3 TXC, 1/3 TXE), so `d(gamma^p) = c_gam d(eps_q^p)` with
  `c_gam = 3 k1 / (1 - sin(psi_m)/3)` (`gamma_factor`, `k1_and_deriv` in
  `LadrunoHSKernel.h`). The Jacobians need `d c_gam/dx` and `d c_gam/d sin(psi_m)`. The
  K0nc cap calibration uses `d(eps_q^p) = dgam (1 - sin(psi_m)/3)/2`.
- **Generalises to:** any model whose hardening variable is defined on principal plastic
  strains (HS, HSsmall, PH/CYSoil in FLAC, HMC) implemented in invariant space. Test it
  in TXE, or at an intermediate Lode angle, against a code that is not your own oracle.
