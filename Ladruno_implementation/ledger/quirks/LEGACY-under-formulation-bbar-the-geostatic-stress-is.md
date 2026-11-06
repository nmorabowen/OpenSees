---
wp: LEGACY
title: "Under -formulation bbar the geostatic stress is piecewise constant at the ELEMENT CENTROID value — which buys an exact 1-D patch check, and traps anyone compar…"
legacy_seq: 358
---
## Under `-formulation bbar` the geostatic stress is piecewise constant at the ELEMENT CENTROID value — which buys an exact 1-D patch check, and traps anyone comparing Gauss-point stresses to the pointwise solution

**Found 2026-09-05, ADR-90 WP-A2.**

A laterally-rollered, base-fixed box under self weight has the closed-form geostatic state
`sigma_zz = gamma*z`, `sigma_xx = sigma_yy = K0 sigma_zz` with `K0 = nu/(1-nu)` — exact even for a
pressure-dependent material such as SANISAND (whose `G` goes as `sqrt(p)`), because `nu` is
constant, so `K0` is depth-independent and the field satisfies equilibrium and the roller
kinematics identically.

Under `-formulation bbar` that field comes back **piecewise constant over each element at the
centroid value**, not varying between the element's own Gauss points: MEASURED
`max |sigma_zz(gp) - gamma*z_centroid| / |gamma*z_centroid| = 1.1e-13` over every Gauss point of
every element, on all three meshes of a refinement sequence.

- **Use it.** It is the *field* check `00_canonical_testbed` §1d demands alongside the resultant
  identity, and it is exact to round-off, so it has no tolerance to tune. A wrong body-force
  convention or a wrong lateral boundary moves it O(1) while the resultant identity sits at 1e-16.
- **The trap:** compare Gauss-point stresses against the POINTWISE `gamma*z_gp` and the same
  correct model reads an O(h) "error" that is not an error — it is the b-bar volumetric average
  being visible. Compare against the centroid value.
