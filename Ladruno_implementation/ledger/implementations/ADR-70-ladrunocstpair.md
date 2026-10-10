---
wp: ADR-70
title: "LadrunoCSTPair"
pr: "#558"
status: "shipped — P4a"
class_tags: ["ELE_TAG_LadrunoCSTPair"]
section: "table"
legacy_seq: 109
---
| **LadrunoCSTPair** — disjoint 2-triangle **F-bar-Patch macro-element** (dSNPO §15.1.9), the T3 volumetric cure the ADR-70 P4 design spike decided ([[70_ladruno_plane_finite_triangles_adr]] §9): 4-node CCW quad split along the n1-n3 diagonal into T1=(n1,n2,n3)/T2=(n1,n3,n4), one centroid GP each, both driven at F̄_e=(J̄/J_e)^{1/2}F_e with the shared patch dilatation J̄=v_patch/V_patch (eq 15.36). **Exact patch-local consistent tangent incl. the stress-proportional cross blocks (eqs 15.37/15.38) — generally UNSYMMETRIC** — assembled through the shared `LadrunoFiniteStrain2DKernel` in 4-node zero-padded rows: ONE `addFbarCoupling2D` call per triangle with the volume-weighted patch row ḡ=Σ_s(v_s/v_p)g_s substituted for the centroid g₀ (zero kernel changes). Finite-strain + PlaneStrain only (`FiniteStrainND2DMaterial`, e.g. LogStrain2D); parser refuses PlaneStress + non-finite materials; det F≤0 step-cut per triangle; `Jbar` element response for pressure diagnostics; `getInitialStiff` = symmetric reference BᵀD₀B seed (family convention). Oracle-anchored: `tests/cstpair_reference.py` (FD-exact tangent at finite stress, row-form ≡ dSNPO-split identity, reduce-to-2-CSTs under homogeneous F, 3-RBM rank) + `tests/test_cstpair_reference.py`. | Element (macro) | **`ELE_TAG_LadrunoCSTPair` 33021** | `SRC/element/ladrunoPlane/{LadrunoCSTPair.{cpp,h},OPS_LadrunoCSTPair.cpp}`, `tests/{cstpair_reference.py,test_cstpair_reference.py,test_ladrunocstpair.py}` | shipped — **P4a** | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
