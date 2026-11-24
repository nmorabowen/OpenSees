---
wp: WP-252
title: "Force-based RC fibre column past a shear failure: the secant section tangent cannot follow a negative axial stiffness, and the global Newton diverges at any step size"
pr: "TBD"
date: 2026-11-23
---
### Force-based RC fibre column past a shear failure: the secant section tangent cannot follow a negative axial stiffness, and the global Newton diverges at any step size (WP-252 adds `-tangent numerical`)
- **Bites:** a `LadrunoForceBeamColumn` with `LadrunoNDFiber` sections of `LadrunoRCConcrete` fibres under constant axial load and imposed lateral displacement (Lynn et al. 1996, 3CLH18), on the first negative excursion after the shear failure. Every displacement step fails, down to 1e-4 mm; arc length on the lateral load fails on its first increment; tiny steps and smaller sub-steps do not help.
- **Why:** the column now lengthens as it loses strength: the axial reaction rises with the tip's axial displacement, a negative axial stiffness at constant load. The material's tangent is the damaged secant `W_B C0`, positive by construction, so every Newton correction in the axial DOF points away from equilibrium. The Newton iterates show it directly: the first lands at N ≈ P, then each correction grows by about 1.3x in the same direction. It is not a snap-back, so arc length on the lateral load does not help.
- **Fix/status (WP-252):** `section LadrunoNDFiber … -tangent numerical` gives the element the forward-difference section tangent, which carries the negative slope; plain Newton then crosses the step. KrylovNewton with small sub-steps also crossed it, with the secant, by building its correction from previous iterates.
- See [[LadrunoNDFiberSection3d_guide]] §3b.
