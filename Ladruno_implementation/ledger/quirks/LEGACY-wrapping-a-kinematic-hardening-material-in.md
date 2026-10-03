---
wp: LEGACY
title: "Wrapping a KINEMATIC-hardening material in LogStrain (finite strain) is NOT objective under large rotation — backstress doesn't co-rotate"
legacy_seq: 53
---
### Wrapping a KINEMATIC-hardening material in `LogStrain` (finite strain) is NOT objective under large rotation — backstress doesn't co-rotate
- **Bites:** `nDMaterial LogStrain $t $j2` over a combined-hardening `LadrunoJ2`
  (or any backstress material) gives finite-strain J2 that is **exact for the
  isotropic part but loses frame-indifference for the kinematic part under a large
  rotation**: superpose a finite rotation `Q` on a plastically-loaded state and the
  Cauchy stress does NOT come back as `Q σ Qᵀ` (principal stresses change).
- **Why:** the log-strain (MATISU) wrapper co-rotates only the elastic state it owns,
  `Bᵉᵗʳ = F_Δ Bᵉ_n F_Δᵀ`. Isotropic yield sees only `‖s‖,ε̄ᵖ` (rotation-invariant), so
  it is objective. The **backstress `α` lives inside the UNCHANGED small-strain inner
  in a FIXED frame** — the wrapper never rotates it — so `‖M‖=‖s−α‖` is not
  rotation-invariant once `α≠0`. It is a framework limit, not a bug: the direct
  Box-14.4 chain shows the identical behaviour. This is exactly the
  kinematic-hardening-at-finite-strain case dSNPO defers to **§14.11**.
- **Rule:** the simple `LogStrain`-wrap of a backstress material is correct only for
  **no / small rotation** (or pure isotropic hardening). For exact large-rotation
  combined hardening you need a finite-strain-NATIVE material that co-rotates `α`
  every step (push `α` forward by the incremental rotation) — which is why the J2
  return map was extracted into the OpenSees-free `LadrunoJ2Kernel.h`: a future
  `FiniteStrainNDMaterial`-subclass J2 can reuse the verified map and add the `α`
  push-forward. Pin the boundary with a **strict xfail** so v2 flips it green
  (`tests/test_ladrunoJ2_finite.py`). Learned 2026-06-02.
