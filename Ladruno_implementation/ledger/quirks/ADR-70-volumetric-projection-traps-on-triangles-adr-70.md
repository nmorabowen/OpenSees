---
wp: ADR-70
title: "Volumetric-projection traps on triangles (ADR 70 P4 spike)"
legacy_seq: 177
---
## Volumetric-projection traps on triangles (ADR 70 P4 spike)

Three traps from the P4 design spike (all numpy-pinned + adversarially
twin-verified in `Ladruno_implementation/adr70_p4_spike/`):

1. **A quadrature L2 projection is the IDENTITY whenever #GP = dim(projection
   space) at unisolvent points.** On the T6's 3-interior-point rule, projecting
   the dilatation onto P1 (3 modes, 3 samples) interpolates — b̃ ≡ b at the GPs
   to ~4e-15, straight OR curved (weights and detJ cancel algebraically). An
   "element-local P1-projected dilatation" formulation flag would silently ship
   `std`. Same trap one level up: patch-P2 over a 2-T6 patch (6 GPs = dim P2)
   is also the identity. Volumetric relief on the triangle REQUIRES coupling
   beyond the element's own quadrature — there is no element-local escape.

2. **Enlarging the constant-mean averaging region can never fix the T6.** The
   P3 conformal modes re-center: u = a(z−z_p)² about the REGION centroid z_p
   has ε_dev ≡ 0 pointwise and zero region-mean dilatation, for any region. So
   dSNPO §15.1.9 F-bar-Patch (one J̄ per patch) cures the T3 pair (exactly 3
   RBM) but leaves the T6 pair with 5 zero-energy modes. The T6 cure must see
   the LINEAR variation of the dilatation (patch-P1 does; rank exactly 3,
   robust distorted + conforming-curved, inf-sup β_h ≈ 0.46 plateau on straight
   structured pair-meshes — unstructured/curved inf-sup unproven).

3. **F-bar-Patch cross-element tangent blocks are STRESS-PROPORTIONAL** (dSNPO
   eqs. 15.37–15.38): K^(es) vanishes at F = I / zero prestress, where the
   macro tangent is symmetric. An FD consistent-tangent gate run at zero
   stress CANNOT catch a wrong/missing/symmetrized cross block — the gate must
   run at finite stress, and the pair element must declare an unsymmetric
   tangent. Also: a shared edge whose mid-nodes don't exactly coincide silently
   cracks the patch and re-opens spurious modes (bit the adversarial twin's own
   curved test — assert conformity).
