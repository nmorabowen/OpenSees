---
wp: WP-144
title: "LadrunoNorSand's C is d(sigma_tensor)/d(eps_tensor) with the shear component INDEPENDENT: an engineering-shear caller HALVES the shear COLUMNS (unlike LadrunoJ…"
date: 2026-11-06
---
### LadrunoNorSand's C is d(sigma_tensor)/d(eps_tensor) with the shear component INDEPENDENT: an engineering-shear caller HALVES the shear COLUMNS (unlike `LadrunoJ2Kernel`) (WP-144)
- **Bites:** `LadrunoJ2Kernel` returns a tangent already expressed for the engineering-shear Voigt vector, so the shell passes it through. The NorSand kernel returns C = d sigma / d eps with the six TENSOR components independent (eps_01 an independent variable, sigma_01 its conjugate). A shell that copies the J2 pattern gets a tangent whose shear COLUMNS are 2x too stiff: the normal block and the stresses are right, so a uniaxial or isotropic test passes, and only a shear-carrying test (or a Newton iteration count) shows it.
- **Rule:** T[a][b] = C[a][b] * w_b, w_b = 1 (b < 3) and 1/2 (b >= 3); rows unchanged (shear stress is the same number in tensor and Voigt form). Verify with a finite difference of the kernel stress against the ENGINEERING strain vector at a plastic, non-coaxial state (`kernel_parity/tools/voigt_check.cpp`). The same applies to a plane-strain view (index set {0,1,3}).
- **Workaround/status:** ✅ `LadrunoNorSand::getTangent` / `getInitialTangent` (WP-144 P1b); pinned by `tools/voigt_check.cpp` and, at G2, by the replay finite-difference test.
