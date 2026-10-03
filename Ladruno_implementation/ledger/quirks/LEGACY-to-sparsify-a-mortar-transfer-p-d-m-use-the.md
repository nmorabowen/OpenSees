---
wp: LEGACY
title: "To sparsify a mortar transfer P=D⁻¹M, use the DUAL basis — row-sum LUMPING D breaks the patch"
legacy_seq: 143
---
### To sparsify a mortar transfer `P=D⁻¹M`, use the DUAL basis — row-sum LUMPING D breaks the patch
- **Bites:** anyone trying to make the mortar slave mass `D` diagonal (so `P` is local/sparse instead of
  dense) — the ADR-62 P2.1 `LadrunoTie -mortar -dual`.
- **Why:** the obvious "diagonalize D" = row-sum LUMPING (`D_II ← Σ_J D_IJ`) keeps partition of unity
  (`ΣP=1`) but is NOT biorthogonal, so it **breaks linear completeness** ⇒ FAILS the constant-stress
  patch (oracle measured 0.28 error vs 2e-15 for dual). The CORRECT diagonalizer is the biorthogonal
  DUAL basis (Wohlmuth): replace the slave TEST functions `N_I` with `ψ_I = Aᵉ·N`, `Aᵉ = diag(cᵉ)(Dᵉ)⁻¹`,
  `cᵉ_a = Σ_b Dᵉ_ab = ∫N_a` ⇒ `Dᵉ_dual = AᵉDᵉ = diag(cᵉ)` EXACTLY (any facet) while preserving the patch.
- **How to apply:** the dual transform is a per-slave-FACET `npsS×npsS` solve applied to that facet's own
  `Mᵉ` — built from the SAME `LadrunoMortarKernel::integratePair` `Dᵉ/Mᵉ` ⇒ **no kernel change, no handler
  change** (same trick as P2's global-D⁻¹ dodge). `P(I,k)=Y(I,k)=(Dᵉ⁻¹Mᵉ)(I,k)`; PoU holds exactly because
  `Dᵉ` and `Mᵉ` integrate over the SAME overlap (`Dᵉ·1 = Mᵉ·1 = cᵉ`). Keep the default (standard dense P)
  byte-identical — the dual path is a separate post-guard block, opt-in `-dual`.
