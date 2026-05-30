---
title: MPCO_Ladruno — standard-rule QUADRATURE + GLOBAL_GP_COORDS + NDIR plan
project: Ladruno
status: plan
owner: nmora
tags:
  - implementation
  - recorder
  - schema
  - plan
---

# Standard-rule QUADRATURE + belt-and-suspenders GLOBAL_GP_COORDS + NDIR

Closes review findings **(a)** (standard-rule elements emit no QUADRATURE → reader KeyErrors)
and **(f)** (simplex `ndir` conflated with `len(ORDER)`), and implements the architect's
**belt-and-suspenders** geometry decision. Branch: `feature/mpco-quadrature-global-gp` (off `ladruno`,
which already has the BezierTri6 bucket-as-group + `basisInfo` probe via PR #18).

## Goals
1. **Standard-rule QUADRATURE** — for legacy Gauss-Legendre elements (quad/hex/tri/tet/line), derive
   GP natural coords + weights from a built-in C++ table and write `QUADRATURE/{GP_PARAM,GP_WEIGHT}` +
   `NUM_GP` like custom-rule elements already do. ("Parity by derivation" — moves STKO's reader-side
   hardcoded table to the write side.)
2. **`GLOBAL_GP_COORDS[nElem×nGP×ndim]`** — static, C++-computed reference-config global GP coords per
   element group, `x = Σ Rᵢ(ξ_g)·Xᵢ` from a write-side basis evaluator over `node->getCrds()`. Write-time
   oracle, kills two-reader basis-math drift, makes legacy elements consumable with zero reader basis math.
3. **`NDIR`** — explicit int attr per element group = number of `GP_PARAM` columns = parametric dim
   (line 1, quad 2, hex 3, tri 2, tet 3). Validator/reader stop deriving `ndir` from `len(ORDER)`.

## Resolved sub-decisions (from the planning pass)
- GP coords source = **built-in C++ table keyed on `ElementIntegrationRuleType::Enum`** (elements don't
  expose `integrationPoints` except force-based beams, which already go through the custom path; **no**
  OpenSees element exposes global GP coords).
- GLOBAL_GP_COORDS = **computed in C++** (write-side basis evaluator), reference config only (static).
- NDIR = **introduce now**.

## Steps
1. `MPCOL_ElementResults.h`: `getStandardQuadrature(rule, gp_param, gp_weight, num_gp, ndir)` — switch on
   the rule enum returning hardcoded abscissae/weights. v1 coverage: line 1/2/3-pt, quad 1/4/9-pt,
   hex 1/8/27-pt, tri 1/3/4-pt, tet 4-pt.
2. `MPCOL_ElementResults.h`: `computeGlobalGP(...)` — C++ mirror of `ladruno_basis.reconstruct_global`
   for v1 topologies (lagrange line/quad/hex tensor + tri/tet linear barycentric).
3. `MPCORecorderLadruno.cpp::writeModelElements`: resolve `(num_gp,ndir,gp_param,gp_weight)` once (custom →
   standard → none); always write `NDIR`+`NUM_GP`+`QUADRATURE` when a rule resolves; compute+write
   `GLOBAL_GP_COORDS` (skip for topologies the evaluator can't do — see R3). No change to value channels.
4. `ladruno_format.py`: reader reads `NDIR`, `GLOBAL_GP_COORDS` optional; validator uses `NDIR` not
   `len(ORDER)`, accepts `bary` with `ndir∈{2,3}`, checks `GLOBAL_GP_COORDS` shape if present.
5. `make_synthetic.py`: add `NDIR` + `GLOBAL_GP_COORDS` + one standard-rule simplex (tri) group.
6. `ladruno_basis.py`: add `bary` + tri/tet linear basis for the oracle.
7. New gate `standard_quad_model.py` + `standard_quad_check.py`: validate() clean; GP_PARAM/GP_WEIGHT ==
   analytic rule to 1e-12; **round-trip oracle** reconstruct x(ξ) vs GLOBAL_GP_COORDS ≤1e-12; element
   value parity vs frozen mpco unregressed.
8. Build (one-time full, then incremental) + re-run all gates (nodal 80/80, quad 96/96, beam 144/144,
   bezier 72/72).

## Open risks / sub-decisions for the architect
- **R1 (highest) — GP ordering parity.** The standard-table GP row order MUST match the element engine's
  gauss-result index (the `gauss_id` in COLUMN_MAP). OpenSees quad/brick GP ordering is element-specific,
  not guaranteed tensor-x-fastest. **Verify empirically** against the frozen recorder for FourNodeQuad(4),
  a brick(8), NineNodeQuad(9) before locking the table; GP_PARAM order follows the *engine*, and the basis
  evaluator must consume that same order.
- **R3 — higher-order Lagrange for GLOBAL_GP_COORDS.** 9N quad / 20-27N hex need order-2/serendipity shape
  fns in both the C++ evaluator and `ladruno_basis.py`. **Recommended v1 scope: option (b)** — write
  QUADRATURE+NDIR everywhere (closes the KeyError), but emit GLOBAL_GP_COORDS + the oracle for **linear
  topologies first**, defer higher-order global coords to v2.
- **R6 — ORDER vs NDIR for simplices / reconcile with landed BezierTri6 (`ORDER=[2,2]`).** With NDIR
  decoupled, decide whether tri/tet ORDER stays total-degree `[p]` (len 1) while NDIR carries dim. The
  landed `bezier_check.py` asserts `ORDER==(2,2)`; changing the convention touches that expectation.
- R2 tri 4-pt negative centroid weight sign/scaling; R4 Tet10 rule (read TenNodeTetrahedron.cpp, likely
  defer); R5 confirm `MPCOL_Hdf5.h` rank-3 support else store `[nElem×nGP*ndim]` 2-D + documented reshape.
