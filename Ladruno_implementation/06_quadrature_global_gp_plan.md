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

## Progress

- **Step A DONE / merged (PR #23):** `getStandardQuadrature(...)` table in
  `MPCOL_ElementResults.h` (line/quad/tri/tet/hex GL rules; dormant until Step B).
- **Step B DONE (this branch `feature/mpco-step-b-global-gp`):**
  - `computeGlobalGP(...)` added to `MPCOL_ElementResults.h` — write-side basis
    evaluator `x(ξ)=ΣN_i(ξ)X_i` for line2/quad4/quad9/tri3/tet4/hex8, shape
    functions verified against each element source (FourNodeQuad CCW, NineNodeQuad
    biquadratic Lagrange, Brick/`shp3d` trilinear node order, Tri31 `N0=ξ,N1=η,N2=1-ξ-η`,
    FourNodeTet `N0=r,N1=s,N2=t,N3=1-r-s-t`). Returns `false` (graceful skip) for
    topologies whose basis isn't implemented yet.
  - `writeModelElements` rewired: resolves `(num_gp,ndir,gp_param,gp_weight)` via
    custom → `getStandardQuadrature` → none; writes `NDIR`+`NUM_GP`+`QUADRATURE` for
    **all** resolved rules; writes `GLOBAL_GP_COORDS[nElem × (nGP·ndim)]` (2-D, R5
    reshape) computed from `node->getCrds()`. Custom-rule `GP_X`/`CUSTOM_*` back-compat
    preserved.
  - Reader/validator (`ladruno_format.py`) made **QUADRATURE-tolerant** (warn, not
    KeyError, when a group has no QUADRATURE), **NDIR-authoritative** (D4 — not
    `len(ORDER)`; accepts `bary` ndir∈{2,3}), and **GLOBAL_GP_COORDS-aware** (optional,
    shape-checked).
  - New gate `standard_quad_{model,check}.py` — FourNodeQuad/Tri31/stdBrick on
    axis-aligned unit cells + **write-time round-trip oracle** (independent Python basis
    reconstructs `x(GP_PARAM[k])` from MODEL/NODES+CONNECTIVITY, compares to
    `GLOBAL_GP_COORDS[k]`). **ALL PASS** (quad/tri max_err 0, brick 2.8e-17 ≤ 1e-12);
    GP_PARAM confirmed CCW/centroid/lexicographic (R1).
  - **No regression:** nodal 80/80, element quad 96/96 + beam 144/144, bezier 72/72,
    multistage 108/108, pytest 10/10.
- **Step D PARTIAL DONE (branch `feature/mpco-step-de-higher-order`):** the two
  higher-order elements whose GP order + node order are *directly verifiable from
  source* now have full QUADRATURE + GLOBAL_GP_COORDS:
  - **NineNodeQuad (quad9, Quad_GL_3 9-pt)** — rule already tabulated (Step A), shape
    fn already in `computeGlobalGP` (Step B); this step just gates it. GP_PARAM matches
    the 9-pt corners/edges/center rule; round-trip 1.1e-16.
  - **TenNodeTetrahedron (tet10, Tet_GL_2 4-pt)** — NEW `Tet_GL_2` in
    `getStandardQuadrature` (α=(5+3√5)/20, β=(5−√5)/20; element's (α,β,β)→(β,β,β)
    cycle, w=1/24); NEW tet10 quadratic shape fn in `computeGlobalGP` (with the
    element's node-8/9 edge swap, verified vs `TenNodeTetrahedron::shp3d`). Round-trip
    1.1e-16.
  - Gate `standard_quad_{model,check}.py` extended (quad9n + TenNodeTetrahedron unit
    cells; checker `shape()` + EXPECT_GP cover quad9/tet10). **ALL PASS**, both
    CONFORMANT; no regression (80/80·96/96·144/144·72/72·108/108·pytest 10/10).
- **Still deferred:** **hex20 (Twenty_Node_Brick, Hex_GL_3 27pt)** — its 27-pt GP order
  lives inside `brcshl`/`Jacobian3d`/`computeBasis` (UP-ucsd `shp3dv`), not directly
  readable, and the round-trip oracle can't catch a GP↔result pairing error → needs a
  full `brcshl` trace + a frozen-recorder `GP_X` cross-check before it's safe. 27N
  Lagrange hex has no element mapping; 6N Lagrange tri / 3N line have no standard-rule
  element (beams use the custom force-based path). Also pending: importable
  `ladruno_basis.py` higher-order oracle + `make_synthetic.py` NDIR/GLOBAL_GP_COORDS/
  simplex fixture additions. Minor: geom-derived `ORDER` for 9N/10N still reports the
  linear `(1[,1])` (NDIR is authoritative, so non-load-bearing).

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
