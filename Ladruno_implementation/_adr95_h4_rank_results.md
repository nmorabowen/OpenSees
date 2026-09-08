# ADR-95 H4: fully-plastic element-tangent rank deficiency

Build `db19a30efce91a1b2bd8cd6c8be95b91dfc7344a`. Script:
`Ladruno_files/testbed/hypo_bearing/p1_plastic_rank.py`.

## VERDICT

H4 (rank-5-per-plastic-GP, DP non-associated) is **confirmed for the simplex
tets and the H20-uri hourglass count, refuted for LadrunoBrick8-bbar (off by
+1) and LadrunoBrick20-std (which is NOT robustly full-rank once plastic)**.

| element | nGP | nDOF | zero-SV elastic | zero-SV plastic | predicted (H4) |
|---|---|---|---|---|---|
| LadrunoBrick `-bbar` | 8 | 24 | 6 | **7** | 6 (full rank) |
| LadrunoBrick20 `-uri` | 8 | 60 | 12 | **20** | 20 = 6 rigid + 14 hourglass — **MATCH** |
| LadrunoBrick20 std | 27 | 60 | 6 | **10** | 6 (full rank) — **REFUTED (+4)** |
| BezierTet10 | 4 | 30 | 6 | **10** | 10 = 30 − 4·5 — **MATCH** |
| BezierTet10 `-bbar` | 4 | 30 | 9 | **13** | "same or +" — measured +4 vs plain |
| TenNodeTetrahedron | 4 | 30 | 6 | REFUSED | 10 (not measured) |
| H20-uri MESH 2x2x2 (8 ele, 243 dof) | 64 | 243 | 6 | REFUSED | n/a (qualitative) |

Elastic baseline is always exactly 6 (rigid modes only) **except** the two
formulations that are already singular before yielding: LadrunoBrick20 `-uri`
(12 = 6 rigid + 6 non-communicable single-element hourglass modes) and
BezierTet10 `-bbar` (9 = 6 rigid + 3 bbar volumetric-averaging modes).

The H20-uri mesh's elastic count is 6 (matches the memory note: hourglass
modes are not communicable in a *connected* elastic mesh), so the interesting
plastic-propagation question is exactly the one the extraction method could
not answer here (see below).

## Method notes

- Homogeneous field: traceless proportional shear `eps0=(1e-5,-0.5e-5,-0.5e-5)`
  via `sp()` at every node (corners+midside, exact at any order), ramped
  `LoadControl` dLambda=0.05, Newton. **Pure triaxial compression never yields**
  with this UWmaterials DP convention (`f1=norm_eta+rho*I1-...`, `rho=0.2>0`
  makes f1 *more negative* as compression grows); a deviatoric-dominant path
  was used instead. All elements reach `ladrunoBranch[0]==1` at every GP by
  lambda≈1.33 (matches `2*sqrt(2)*G*eps0*lambda = sqrt(2/3)*sigma_y`).
- Extraction: `eleResponse(tag,'stiff')` preferred (LadrunoBrick, Brick20,
  BezierTet10 all expose a direct `getTangentStiff()` dump, no domain solve,
  immune to singularity). TenNodeTetrahedron/H20-mesh have no such response
  and fall back to `printA('-ret')` after removing all constraints — two
  OpenSees gotchas hit there (Plain handler silently zeroes non-homogeneous
  sp values; `Domain::removeSP_Constraint(node,dof,patternTag)` is broken for
  pattern-scoped SPs — only `remove('loadPattern',tag)` actually detaches
  them); see script docstring.
- **Refused**: TenNodeTetrahedron and the H20-uri 2x2x2 mesh extract cleanly
  elastically via `printA`, but plastically the assembled tangent is exactly
  multiply rank-deficient and the LAPACK dense solve inside `analyze(1)` hits
  an exact-zero pivot, poisoning the in-place buffer with NaN before `printA`
  reads it — a fallback-method limitation, not a physical claim. A quick
  `SparseGeneral`/RCM swap also failed (0-length matrix). Not pursued further.
