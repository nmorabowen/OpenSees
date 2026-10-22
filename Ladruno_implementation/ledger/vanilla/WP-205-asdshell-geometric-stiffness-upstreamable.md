---
wp: WP-205
title: "ASDShellQ4/T3: consistent corotational tangent (EICR rotation gradient + sign) and -pDelta initial-stress stiffness"
date: 2026-10-22
table: "upstreamable"
files: ["`SRC/element/shell/ASDShellQ4CorotationalTransformation.h`", "`SRC/element/shell/ASDShellT3CorotationalTransformation.h`", "`SRC/element/shell/ASDShellQ4.cpp`", "`SRC/element/shell/ASDShellQ4.h`", "`SRC/element/shell/ASDShellT3.cpp`", "`SRC/element/shell/ASDShellT3.h`"]
---
| `SRC/element/shell/ASDShellQ4CorotationalTransformation.h` | `// Ladruno WP-205`: with `USE_POLAR_DECOMP_ALLIGN` (always defined) `RotationGradient` returned **G = 0**, which zeroes both EICR geometric terms (K_GP, K_GR) and reduces P to the translational projector: the corotational tangent had NO geometric stiffness while the residual (exact frame) did — no eigen softening after gravity + `loadConst`, ~34 Newton iterations on a P-Delta step. G now carries the out-of-plane spin of the diagonal-cross normal (rows 0-1) and the linearised in-plane spin of the polar frame, dtheta = (dF21 - dF12)/tr(U) with reference-centre gradients (row 2). Exposing G also exposed the stock sign: K_GP and K_GR were ADDED ("RHS already has the negative sign" — the Kratos convention); here the RHS is +internal force, so both are now subtracted (EICR). Upstream bug, candidate for upstream. | WP-205 |
| `SRC/element/shell/ASDShellT3CorotationalTransformation.h` | `// Ladruno WP-205`: same two fixes for the triangle (G from the constant deformation gradient; the stock non-polar branch also had node 3's out-of-plane row wrong, x32 for x21 — unused, left as is). | WP-205 |
| `SRC/element/shell/ASDShellQ4.cpp` | `// Ladruno WP-205`: new `-pDelta` option (exclusive with `-corotational`): linear kinematics plus the initial-stress stiffness K_sigma = int G'[N11 N12; N12 N22]G dA of the current membrane forces on the three local translations (w: out-of-plane P-Delta; u/v: in-plane P-Delta of a wall swaying in its plane), and its forces K_sigma*UL in the residual (second-order equilibrium, the PDelta-transformation convention; dK_sigma/du*UL omitted, as PDelta does). Initial stiffness stays material-only. `sendSelf`/`recvSelf` ID 21 -> 22 (pDelta flag). | WP-205 |
| `SRC/element/shell/ASDShellQ4.h` | `// Ladruno WP-205`: `setPDelta(bool)`, `m_pdelta`. | WP-205 |
| `SRC/element/shell/ASDShellT3.cpp` | `// Ladruno WP-205`: `-pDelta` as in ASDShellQ4 (constant-gradient triangle); `sendSelf`/`recvSelf` ID 18 -> 19. | WP-205 |
| `SRC/element/shell/ASDShellT3.h` | `// Ladruno WP-205`: `setPDelta(bool)`, `m_pdelta`. | WP-205 |
