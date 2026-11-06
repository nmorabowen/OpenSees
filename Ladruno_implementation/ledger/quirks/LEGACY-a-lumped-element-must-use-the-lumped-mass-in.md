---
wp: LEGACY
title: "A -lumped element must use the lumped mass in the inertia RESIDUAL too — upstream Brick lumps only the tangent"
legacy_seq: 527
---
### A `-lumped` element must use the lumped mass in the inertia RESIDUAL too — upstream `Brick` lumps only the tangent
- **Bites:** upstream `Brick` (and LadrunoBrick until WP-139) builds the inertia residual from the CONSISTENT mass `ρ Σ N_j N_k dV` for every `massType`, while `-lumped` changes only the tangent branch (`getMass`: Newton tangent, αM Rayleigh, ground load) to the row-sum diagonal. Implicit dynamics then integrates a hybrid that no textbook describes, and Newton's Jacobian is not the derivative of the residual: an ELASTIC `-lumped` Newmark step converges only linearly (no 1e-10 in 40 iterations). Vanilla `CentralDifference` reads the residual at a nonzero trial acceleration and integrated the same hybrid. It hid because every total-based check (rigid-body translation, Σ forces, the ground-motion probe) sees only row sums — equal for lumped and consistent mass — and `CentralDifferenceLadruno`'s Azero residual skips the inertia pass entirely.
- **Rule:** one mass model per element: whatever `getMass()` returns, the residual's `M·a` must use the same matrix (LadrunoBrick20 F-1; LadrunoBrick WP-139). Test it node by node with a NON-uniform acceleration field, and with Newton iteration counts.
- **Workaround/status:** ✅ FIXED in LadrunoBrick (WP-139, owner option A). Upstream `Brick` still carries it (vanilla, left alone). *2026-09-27.*
