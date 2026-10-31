---
wp: WP-205
title: "ASDShellQ4/T3 -corotational: the tangent had no geometric stiffness (G = 0 under the polar frame) and the EICR sign was the Kratos one"
pr: "#1011"
date: 2026-10-22
---
### `ASDShellQ4 -corotational` gets static P-Delta right but eigen after gravity shows NO softening and Newton needs ~34 iterations

- **What bites.** 30 m cantilever wall strip, P = 0.5 Pcr, gravity -> `loadConst` -> eigen: `-corotational` T = 25.10 s (the unloaded period) vs 35.36 s for ShellNLDKGQ and theory; the lateral step takes 34 Newton iterations vs 7. The static amplification was right (1.98), so the RESIDUAL was exact and only the TANGENT was wrong. Same for ASDShellT3.
- **Why.** `USE_POLAR_DECOMP_ALLIGN` is defined in both corotational transformations, and under it `RotationGradient` sets G = 0 on the claim that the polar-aligned frame has a zero spin-lever. It does not: the polar alignment fixes only the in-plane orientation; the normal follows the nodes. With G = 0, P = Pt and both EICR geometric terms (K_GP = G'Fn'P, K_GR = Fnm G) vanish. The second, latent, bug: the code ADDS both terms with the comment "RHS already has the negative sign" (Kratos convention); OpenSees' RHS here is +internal force. With only G restored the period dropped to 20.57 s (stiffening by exactly 1/sqrt(1.5)).
- **Fix (WP-205).** Restore G (rows 0-1 from the diagonal-cross normal, row 2 = linearised polar spin (dF21 - dF12)/tr(U)) and subtract K_GP, K_GR. Wall strip: T = 35.33 s, 4 lateral iterations; P = 0 unchanged.
- **For linear eigen/RSA about gravity** use the new `-pDelta` (linear kinematics + K_sigma of the current membrane forces), the shell counterpart of a PDelta geomTransf.
