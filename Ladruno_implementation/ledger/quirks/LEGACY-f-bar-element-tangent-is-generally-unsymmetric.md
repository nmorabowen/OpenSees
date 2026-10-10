---
wp: LEGACY
title: "F-bar element tangent is GENERALLY UNSYMMETRIC (needs an unsymmetric solver)"
legacy_seq: 52
---
### F-bar element tangent is GENERALLY UNSYMMETRIC (needs an unsymmetric solver)
- **Bites:** `LadrunoBrick -geom finite -formulation bbar` (F-bar, dSNPO §15.1). The
  consistent tangent (eq 15.10) is `K = ∫Gᵀa G dv + ∫Gᵀq(G₀−G) dv`; the second
  (F-bar coupling) term is **not symmetric** in general — the book says so
  explicitly (after eq 15.10): "the additional stiffness term … is generally
  unsymmetric and, therefore, requires an unsymmetric solver." So `system BandSPD`
  / `ProfileSPD` (symmetric storage) will silently use only the upper triangle and
  **converge to the wrong answer or stall**. Use `system FullGeneral` (or any
  unsymmetric solver). The symmetric storage path will *not* error — it just drops
  the lower triangle, so this is a silent-wrong-answer trap.
- **Why:** F̄ = (J₀/J)^(1/3)F couples every Gauss point's stress to the element
  centroid dilatation; the linearization mixes the GP gradient G with the centroid
  gradient G₀, and `q⊗(G₀−G)` is a non-symmetric outer product. (The plain
  `-formulation std` finite tangent stays symmetric — the coupling term is absent.)
- **Also:** the eq 15.11 coupling tensor is `q = (1/3)a:(I⊗I) − (2/3)σ⊗I` with `a`
  the **full** spatial tangent (= `c − σδ`, the same modulus as the standard term),
  **not** the material modulus `c` alone. A first-principles spatial shortcut
  (`dσ̄ = c̄:sym(L̄)`) drops the `−(2/3)σ⊗I` and gives a subtly wrong tangent that
  still "looks right" under a crude FD check — the element FD-tangent test against an
  *analytic* material tangent is what catches it. Learned 2026-06-02 (F-bar impl).
