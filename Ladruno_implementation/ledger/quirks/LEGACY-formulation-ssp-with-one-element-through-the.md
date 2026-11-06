---
wp: LEGACY
title: "-formulation ssp with ONE element through the bending depth never plastifies — its P–δ curve is identical to the elastic one"
legacy_seq: 372
---
### `-formulation ssp` with ONE element through the bending depth never plastifies — its P–δ curve is identical to the elastic one
- **Bites:** you mesh a wall/beam one element thick in the bending direction, use
  `-formulation ssp` with a perfectly good inelastic material, and the model
  returns an exactly elastic pushover. No warning: the material is present,
  committed, and reports zero plastic strain because it is never strained.
  `std`/`bbar` form a hinge on the same mesh, so an A/B against a "reference"
  formulation is what exposes it.
- **Why:** the single-point forms evaluate the constitutive model **once, at the
  element centroid** (`isSinglePoint()`, the PR #94 cost fix). Bending strain is
  antisymmetric about the centroid, so `ssp`'s mean-dilatation core sees ~zero
  strain and never yields; the entire bending response is carried by the
  artificial `Kstab`, which is elastic. Measured: `Ehg/W_ext = 98%`, and the
  `elastic` and `J2Plasticity` runs agree to every printed digit.
- **Scope — `uri` is NOT the same:** `uri`+`stiffness` strains its centroid from
  the plain centroid `B`, picks up transverse shear, and does yield, at an
  aspect-ratio-dependent point (elastic to `u/L ≈ 4%` on unit cubes; prompt
  hinging on 1.0 × 0.25 × 2.0 elements). Don't widen the rule to "single-point
  formulations" — both cases are pinned in the regression test.
- **No stabilization-degradation scheme can fix this** — a damage model would not
  damage at that centroid either. It is inherent to one-point integration.
- **Rule:** `nd ≥ 4` elements through the bending depth is the threshold that
  avoids the *never-plastifies* pathology — it is **not** a threshold for an
  accurate hinge load. Measured elastic-calibrated load bias vs `bbar` is
  **+77 / +45 / +15% at `nd` = 2 / 4 / 8**. Production RC walls meshed 2–3 thick
  sit squarely in the +45–77% band.
- **Usually change formulation, not mesh:** on the SAME coarse meshes
  `-formulation eas` biases **+7.5 / +1.0 / −1.5%**, because true Simo-Rifai
  re-condenses `K* = Kdd − Kda Kaa⁻¹ Kad` from the CURRENT tangent at 8 live GPs
  every assembly, where `ssp` freezes one condensation of `C(0)`. `eas` and `ssp`
  agree elastically to 3–4 digits on that mesh family, so the whole plastic gap is
  the stabilization treatment. Costs: `eas` is small-strain-only in this fork and
  pays 8 live GPs + an inner Newton vs `ssp`'s single material evaluation — an
  implicit-analysis choice, not an explicit one. `std`/`bbar` remain the
  no-stabilization reference.
- *2026-07-30 (Tier-A `Kstab` scope challenge).*
