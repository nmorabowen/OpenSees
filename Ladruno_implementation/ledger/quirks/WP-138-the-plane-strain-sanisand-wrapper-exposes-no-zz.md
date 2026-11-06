---
wp: WP-138
title: "The plane-strain SANISAND wrapper exposes no σ_zz (getStressZZ is NaN, stress has 3 components) — recover it EXACTLY from psi and the void ratio (WP-138)"
legacy_seq: 509
---
### The plane-strain SANISAND wrapper exposes no σ_zz (`getStressZZ` is NaN, `stress` has 3 components) — recover it EXACTLY from `psi` and the void ratio (WP-138)
- **Bites:** anyone dumping a `LadrunoSANISAND` plane-strain Gauss point for a material-point replay (`ladrunoSANISANDReplay`, WP-134's oracle), which needs all six stress components. `ManzariDafaliasPlaneStrain::getStress()` returns σ_xx, σ_yy, σ_xy; `getStressToRecord()` has σ_zz but no response reaches it; `LadrunoQuad`'s `stressZZ` column reads NaN from the `NDMaterial` default. SANISAND is hypoelastic, so σ_zz cannot be rebuilt from the elastic strain in `state` either.
- **Workaround:** with `-Presidual 0` the `psi` response is `e − (e0 − λc (p/Pa)^ξ)` at the committed p = tr(σ)/3, and `state[24]` is e, so p = Pa·((e0 − e + ψ)/λc)^(1/ξ) and σ_zz = 3p − σ_xx − σ_yy. Check it against `yieldDistance`: recomputing f from the recovered six-vector and `alpha` matched to ≤ 1.5e-12 on 9 720 points (WP-138 footing, `footing_ab.py::derived`). With p_r ≠ 0 subtract p_r first. *2026-09-27.*
