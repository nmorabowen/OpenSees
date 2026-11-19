# Limit analysis: the P0 interface (WP-237a)

**Status:** **frozen at `03f59c1a4`** (2026-11-19, after the pre-freeze round), with **amendment 1** recorded after
it (section 19.1). A change after the freeze is an amendment logged in section 19, never a silent edit. Section 18
separates accepted owner decisions from orchestrator P0 rulings.

**Who reads it.** The Python reference author, the independent checker author, the external reference
author and later the C++ author. Each implements from this file and `tests/limit_cases.txt`, without
reading the other side's code. `tests/test_limit_contract.py` holds the schema validators and the
convention checks, and `evaluate_case`, the executable meaning of the case lines (section 17.9). It contains
no assembler and no checker. `Ladruno_files/testbed/limit_analysis/certificate_demo.py` contains both, for a
demonstration only (section 18.3): **the 237b, 237c, 237e and 237f authors do not read it.**

**Scope of version 0.** Associated perfect plasticity on straight tetrahedra; von Mises and
Drucker-Prager; second-order cones only; bonded rigid footings; one-shot solves on one mesh. Footing cases
use `z` upward. The file format does not presuppose the C++ engine or the final `.la.h5` layout.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| `E`, `N`, `n_f` | numbers of tetrahedra, mesh nodes, footings |
| `e`, `i` | element index; local vertex `i` in 0..3 |
| `sigma_ei` | stress at vertex `i` of element `e`, Voigt, 6 components |
| `d_ei` | strain rate at vertex `i` of element `e`, Voigt, 6 components |
| `v_a` | velocity of mesh node `a` |
| `V_f = (u_f, w_f)` | rigid velocity of footing `f`: translation, then angular velocity |
| `R_f = (F_f, M_f)` | resultant applied by footing `f` on the soil: force, then moment about `x_ref` |
| `k`, `alpha` | strength data of a region (section 4) |
| `b0`, `b1` | fixed and reference body force per unit volume |
| `t0`, `t1` | fixed and reference traction |
| `F0`, `F1` | fixed and reference footing resultants |
| `lambda` | load multiplier |
| `J2`, `I1` | `0.5 s:s`, `tr(sigma)` with `s` the deviator |
| `theta`, `e` | `tr(D)` and the deviator of the strain-rate tensor `D` |

All indices are 0-based. Arrays are row-major.

## 2. Conventions

| Item | Decision |
|---|---|
| Units | Any consistent set; nothing converts. `units` strings are labels. |
| Sign | **Tension positive.** Compression is negative stress. Dilation is positive `theta`. |
| Voigt order | `[xx, yy, zz, xy, yz, zx]` for stress and strain rate. |
| Shear | Stress carries tensor shear. Strain rate carries **engineering** shear, `d_xy = 2 D_xy`. So `sigma^T d = sigma : D`. |
| Normal | Facet normals are unit and **outward from the element** that owns the facet. |
| Traction | `t = sigma n`, the force per area on the body across its boundary with outward normal `n`. |
| Body force | `b` is a force per unit volume acting on the soil, signed. A unit weight `gamma` acting along `-y` is `(0, -gamma, 0)`. No axis is privileged. |
| Footing resultant | `R_f` is the action **of the footing on the soil**: `F = int t dA`, `M = int (x - x_ref) x t dA`, with `t = sigma n` and `n` outward from the soil. A downward bearing load on a top surface has a negative vertical force. |
| Resultant order | `[Fx, Fy, Fz, Mx, My, Mz]`, conjugate to `[ux, uy, uz, rx, ry, rz]`, with power `R^T V`. |
| Plane strain | Not a separate mode. It is a 3D slab with `symmetry` faces (section 16.4). |

### 2.1 The operator `N(n)`

For a vector `n`,

```text
N(n) = [ n1  0   0   n2  0   n3 ]
       [ 0   n2  0   n1  n3  0  ]
       [ 0   0   n3  0   n2  n1 ]
```

It gives all three maps used here:

- traction: `t = N(n) sigma`;
- divergence of `sigma(x) = sum_i L_i(x) sigma_i`: `div sigma = sum_i N(grad L_i) sigma_i`;
- strain rate of `v(x) = sum_a N_a(x) v_a`: `d = sum_a N(grad N_a)^T v_a`.

## 3. Mesh

### 3.1 Node order

A `tet10` row lists four vertices, then six mid-edge nodes:

| Local node | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| Position | v0 | v1 | v2 | v3 | (0,1) | (1,2) | (2,0) | (0,3) | (1,3) | (2,3) |

A `tet4` is the first four entries. This is the VTK quadratic-tetrahedron order, and it is **not** the
gmsh order. A bridge must permute, and should do so by matching physical coordinates.

**Orientation.** `det[x1-x0, x2-x0, x3-x0] > 0` for every element. A negative element is refused, not
reordered.

**Straightness.** Every mid-edge node lies at its chord midpoint within `TOL_GEOM_STRAIGHT * h_e`, where
`h_e` is the longest edge of the element. The file's mid-edge coordinates are **only checked**. Every
operator uses the four vertex coordinates, so geometry is exactly affine.

### 3.2 Facet order

Local facet `f` is opposite vertex `f`. The order makes the right-hand normal point outward.

| Facet | Vertices | All six P2 nodes: vertices, then edges (a,b), (b,c), (c,a) |
|---|---|---|
| 0 | 1, 2, 3 | 1, 2, 3, 5, 9, 8 |
| 1 | 0, 3, 2 | 0, 3, 2, 7, 9, 6 |
| 2 | 0, 1, 3 | 0, 1, 3, 4, 8, 7 |
| 3 | 0, 2, 1 | 0, 2, 1, 6, 5, 4 |

`n = (xb - xa) x (xc - xa) / (2 A)`, with `A` the facet area.

### 3.3 Geometry formulas

With `J = [x1-x0, x2-x0, x3-x0]` (columns): `V = det J / 6`; rows 0..2 of `inv(J)` are `grad L_1..3`;
`grad L_0 = -sum`. Equivalently `grad L_i = -(A_i / (3 V)) n_i`, with `A_i`, `n_i` of the facet opposite
vertex `i`. A checker may use this second route.

P2 shape functions: `N_i = L_i (2 L_i - 1)`; `N_ab = 4 L_a L_b`.

At vertex `m` (`L_m = 1`, others 0):

```text
grad N_m = 3 grad L_m
grad N_i = - grad L_i                (i != m)
grad N_ab = 4 grad L_b  if a = m;   4 grad L_a  if b = m;   0 otherwise
```

### 3.4 Canonical facet lists

- An **exterior** facet belongs to one element. `facets` lists every exterior facet once as
  `(element, local_facet)`, sorted lexicographically.
- An **interior** facet belongs to two elements. Its key is `(e-, f-)`, the entry with the smaller
  element index, and the list is sorted by that key. Its normal is the outward normal of `e-`. Its
  vertex order is that of `e-`'s facet.
- The mesh must be conforming: shared edges share their mid-edge node and no facet has three owners.

### 3.5 Structured Kuhn mesh

Cases use this generator so all implementations build identical meshes.

- Box `[x0, x1]`, cell counts `(nx, ny, nz)`. Vertex `(i, j, k)` has id `i + (nx+1)(j + (ny+1) k)`.
- Cell `(i, j, k)` has six tetrahedra, one per axis permutation in lexicographic order
  `(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)`. A tetrahedron walks from the cell's minimum
  corner, adding one unit along each axis of the permutation. If its determinant is negative, swap
  local vertices 2 and 3.
- Element id: `6 * cell_id + permutation_index`, `cell_id = i + nx (j + ny k)`.
- Mid-edge node ids start after the vertices, in lexicographic order of the sorted vertex pair. The
  same rule turns a `mesh_file` tet4 mesh into tet10.
- Doubling every count refines every tetrahedron, because the planes `x_i = c` and `x_i - x_j = c` of
  the coarse mesh remain planes of the fine mesh. Doubling only some counts is **not** nested.

### 3.6 Uniform refinement of a tet4 mesh

`refine r` applies this rule `r` times. Edge midpoints are new vertices, numbered after the old ones in
lexicographic order of the sorted vertex pair. With `mab` the midpoint of local edge `(a, b)`, tetrahedron `e`
yields children `8 e + 0..7`:

```text
(v0, m01, m02, m03)   (m01, v1, m12, m13)   (m02, m12, v2, m23)   (m03, m13, m23, v3)
(m01, m02, m03, m13)  (m01, m02, m12, m13)  (m02, m03, m13, m23)  (m02, m12, m13, m23)
```

If a child's determinant is negative, swap its local vertices 2 and 3. Each child has volume `V/8`. Every facet
splits into four, so conformity is kept, and the fine mesh is nested in the coarse one.

## 4. Strength

### 4.1 Data

A region carries `criterion`, `k`, and, for Drucker-Prager, `alpha`.

| Criterion | Yield function, tension positive | Admissible data |
|---|---|---|
| `von_mises` | `sqrt(J2) - k <= 0` | `k > 0`; `alpha` absent and taken as 0 |
| `drucker_prager` | `sqrt(J2) + alpha I1 - k <= 0` | `0 < alpha < 1/sqrt(3)`, `k >= 0` |

`k` is the strength in pure shear. For von Mises, `k = sigma_y / sqrt(3)`. `tresca` and `mohr_coulomb`
are reserved for P5 and are refused in version 0. **In 3D, Drucker-Prager is not Mohr-Coulomb.**

### 4.2 Cone form

Define

```text
Q = [ 1/sqrt2  -1/sqrt2   0        0       0       0     ]
    [-1/sqrt6  -1/sqrt6   2/sqrt6  0       0       0     ]
    [ 0         0         0        sqrt2   0       0     ]
    [ 0         0         0        0       sqrt2   0     ]
    [ 0         0         0        0       0       sqrt2 ]

one = [1, 1, 1, 0, 0, 0]^T
```

`|Q sigma|^2 = 2 J2` and `Q one = 0`. The stress is admissible exactly when

```text
( sqrt2 (k - alpha one^T sigma),  Q sigma ) in SOC6,     SOC6 = { (t, w in R^5) : t >= |w| }.
```

### 4.3 Support function

`pi(d) = sup { sigma^T d : sigma admissible }`, with `theta = one^T d` and `|e|` the Frobenius norm of
the deviator of `D`, `|e|^2 = |e_diag|^2 + (d_xy^2 + d_yz^2 + d_zx^2)/2`.

| Case | Domain | Value |
|---|---|---|
| `alpha = 0` | `theta = 0` | `sqrt2 k |e|` |
| `alpha > 0` | `theta >= 3 sqrt2 alpha |e|` | `k theta / (3 alpha)` |
| outside the domain | | `+infinity` |

Conic form, used by every Upper formulation:

```text
d = sqrt2 alpha r one + Q^T z,   (r, z) in SOC6,   pi(d) = min sqrt2 k r.
```

With `alpha = 0`, the minimum is at `r = |z|`. With `alpha > 0`, `r = theta / (3 sqrt2 alpha)` is forced.

### 4.4 Plane-strain sections (derived)

Take `R = sqrt(((sxx - syy)/2)^2 + sxy^2)`, `m = (sxx + syy)/2`, `delta = szz - m`, no out-of-plane
shear. Then `J2 = R^2 + delta^2/3` and `I1 = 3 m + delta`. Minimising over `delta` gives the section

```text
R sqrt(1 - 3 alpha^2) + 3 alpha m <= k.
```

So von Mises gives Tresca with `cu = k`. Drucker-Prager gives Mohr-Coulomb `R + m sin(phi) <= c cos(phi)`
exactly when

```text
alpha = tan(phi) / sqrt(9 + 12 tan^2 phi),     k = 3 c / sqrt(9 + 12 tan^2 phi).
```

This match holds only for plane strain with associated flow.

## 5. Loads and boundary

**Load form.** Every load is `fixed + lambda * reference`: `b0 + lambda b1` per region; `t0 + lambda t1`
per traction facet; `F0 + lambda F1` per footing. Traction varies linearly on a facet and is stored at
its three vertices. Body force is constant in a region. A region's strength is constant.

**Gravity** is `b0 = gamma g_hat`, per unit volume. A gravity multiplier uses `b1` instead.

| Kind | Lower, at each facet vertex | Upper, at all six facet nodes of the element |
|---|---|---|
| `free` | `N(n) sigma = 0` | none |
| `traction` | `N(n) sigma = t0 + lambda t1` | none; adds work |
| `fixed` | none | `v = 0` |
| `symmetry` | `T1^T N(n_S) sigma = T2^T N(n_S) sigma = 0` | `n_S^T v = 0` |
| `footing` | enters the resultant (section 6) | rigid map (section 6) |
| `interface` | reserved, refused in version 0 | reserved |

**Tangent basis.** `m` is the index of the smallest `|n_m|`, the first if tied. `T1 = (e_m x n)/|e_m x n|`
and `T2 = n x T1`.

**Symmetry sets** must be planar. All their facets share one normal and one plane within
`TOL_GEOM_PLANE`. The set's normal is the normal of its first facet in `facets` order.

**Symmetry rows use the declared set normal (amendment 1, item 3).** In both bounds every symmetry row of set `S` uses
`n_S`, the set normal above, never the facet's own computed normal: Lower rows `T1^T N(n_S) sigma` and `T2^T N(n_S) sigma`,
with `T1`, `T2` the tangent basis of `n_S`; Upper rows `n_S^T v`. The program (15.4) and the certificate projections
(17.3, 17.4) use the same rows. Facet normals of a planar set differ from `n_S` only by rounding, so this fixes one row
per node or vertex and makes independent assemblies produce identical rows.

**Corners.** Conditions apply per facet in the Lower problem and accumulate per node in the Upper problem.

## 6. Rigid footings

A footing has a reference point `x_ref`, one or more `footing` boundary sets, and an `interface_law`.
Version 0 admits `bonded`. `no_tension` and `coulomb` are reserved and refused.

- **Lower.** No extra traction variable. With facet vertices `a` and area `A`,

  ```text
  F = sum_facets sum_a (A/3) t_a
  M = sum_facets sum_a sum_b A w_ab (x_a - x_ref) x t_b,   w_ab = 1/6 if a = b, else 1/12
  ```

  with `t_a = N(n) sigma_(e,a)`. Both are exact for linear traction on a flat facet.
- **Upper.** Every node of a footing facet follows `v_a = u + w x (x_a - x_ref)` exactly. P2 contains the
  linear rigid field, so the whole facet moves rigidly.
- **Power.** `int t . v dA = F . u + M . w = R^T V`.
- **Restrained components.** A query may list components as restrained. Then `V_c = 0`, and `R_c` is an
  unconstrained reaction. A symmetric half model needs this: the half footing of a strip with `z` upward and
  thickness along `y` restrains `ux` and `ry`.

## 7. Queries

Let `K` be the convex set of resultants, on the unrestrained components, that the soil can carry with the
fixed loads `b0`, `t0` applied. All forms take `restrained`, a list of `[footing, component]`.

| Form | Data | Solves | `value` |
|---|---|---|---|
| `ray` | `F0`, `F1`, shape `(n_f, 6)` | maximise `lambda` with `R = F0 + lambda F1` | `lambda` |
| `fixed_components` | `fixed: [[f, c, value], ...]`, `direction` `(n_f, 6)` | a ray with `F0` = the fixed values, `F1` = `direction` | `lambda` |
| `support` | `direction` `d`, `(n_f, 6)` | maximise `d^T R`, with `b1` and `t1` all zero | `d^T R` |

Rules:

- `F0`, `F1` and `direction` are exactly zero on restrained components.
- For `fixed_components`, `direction` is zero on fixed components and non-zero somewhere.
- A `ray` must have a non-zero reference load somewhere in `F1`, `t1` or `b1`.
- A `support` query needs `n_f >= 1` and a non-zero `d`.
- `d` has velocity units: translation, then angular velocity.
- Problem `length_scale` `ell > 0` defines the reported coordinates `R_hat = (F, M/ell)` and
  `d_hat = (d_u, ell d_w)`, so `d^T R = d_hat^T R_hat`.

**Meaning.** Let `lambda_max = sup { lambda : F0 + lambda F1 is carried }`. A Lower solution gives
`lambda_L <= lambda_max`; an Upper solution gives `lambda_max <= lambda_U`. For `support`, a Lower
solution gives a point `R_L` in `K`, and an Upper solution gives a half-space `{ d^T R <= U }` containing
`K`. So `conv(lower points) subset K subset intersection(upper half-spaces)`. A negative `lambda_L` means
the fixed loads are not shown safe.

## 8. Formulations are named and versioned

A **problem** states physics only. A **formulation id** states the discretisation. It is carried by every
result, program, certificate and feasibility record, and is not part of the problem hash. Case semantics
never change when a formulation is added. A new formulation gets a new id and new `formulation` lines in
a case, and old ids keep their meaning.

| Id | Bound | Version 0 status |
|---|---|---|
| `lower.p1d-vertex.v0` | Lower | admitted; section 9 |
| `upper.p2d-jump-bernstein.v0` | Upper | admitted; **the P0 Upper formulation**; section 12 |
| `upper.p2c-vertex.v0` | Upper | admitted as the **control**; section 10. It is not proven for 3D isochoric flow |

A case's `formulation` lines say which formulations its assertions and gates bind. A `control` line records
a formulation for comparison and binds nothing. Any other id is refused with `formulation_not_admitted`.

## 9. Lower problem, `lower.p1d-vertex.v0`

Unknowns: `sigma_ei`, linear and discontinuous per element; `lambda` for ray forms.

```text
maximise   lambda                      (ray, fixed_components)
           d^T R(sigma)                (support)
subject to
  sum_i N(grad L_i) sigma_ei + b0_e + lambda b1_e = 0                each element
  N(n) (sigma_(e-,a) - sigma_(e+,a)) = 0                             each interior facet, 3 vertices
  boundary rows of section 5                                         each exterior facet, 3 vertices
  R_fc(sigma) = F0_fc + lambda F1_fc                                 each unrestrained component
  ( sqrt2 (k_e - alpha_e one^T sigma_ei), Q sigma_ei ) in SOC6       each element, 4 vertices
```

## 10. Upper control, `upper.p2c-vertex.v0`

This continuous P2 problem is kept as a control. On the committed unstructured meshes with fixed walls its
isochoric space does not let the footing move (section 11). Unknowns: continuous P2 `v_a`; `V_f`;
`(r_ei, z_ei)`.

```text
minimise   D - W0
D  = sum_e sum_i (V_e / 4) sqrt2 k_e r_ei
subject to
  B_ei v = sqrt2 alpha_e r_ei one + Q^T z_ei,  (r_ei, z_ei) in SOC6   each element, 4 vertices
  boundary and footing rows of sections 5, 6
  V_fc = 0                                                            restrained components
  W1 = 1                                                              ray forms
  V_fc = d_fc                                                         support, unrestrained components
```

`B_ei v = sum_a N(grad N_a(x_i))^T v_a`, with the gradients of section 3.3.

Work rates, for a form with reference loads:

```text
W0 = sum_e b0_e . int v dV + sum_facets int t0 . v dA + F0^T V      (F0 term absent for support)
W1 = sum_e b1_e . int v dV + sum_facets int t1 . v dA + F1^T V
```

Exact integrals on straight elements:

- `int N dV` is `-V/20` at each vertex node and `V/5` at each edge node.
- `int L_i N_j dA`, with `L_i` the traction basis at facet vertex `i`:

| P2 node `j` on the facet | Value |
|---|---|
| vertex `i` | `A/30` |
| another vertex | `-A/60` |
| edge touching vertex `i` | `2A/15` |
| the opposite edge | `A/15` |

## 11. True-3D Upper gate

Thin-slab success does not show that an Upper formulation works in 3D, and **feasibility alone never passes
this gate.** Each `upper_gate_3d` case runs one of `phi = 0, 10, 20, 30` degrees with one wall type. That is
eight cases.

**Boxes.** A bonded unit-square footing conforms on `z = 0` at `[-0.5, 0.5]^2`. Walls and base are all `fixed` or
all `symmetry` (rollers). Each angle has its own box, sized so that the plane-strain Prandtl extent of a footing of
half-width `b = 0.5` (section 16.4), **taken as a length scale**, fits inside it:

| `phi` | Extent from the centreline, depth | Box | Meshes `a`, `b`, `c`, `d` (tets) | Mesh 5 |
|---|---|---|---|---|
| 0 | 1.50, 0.71 | `[-2, 2]^2 x [-2, 0]` | `upper3d_{a,b,c,d}`: 243, 463, 949, 2,840 | 22,720 |
| 10 | 2.07, 1.03 | `[-2.5, 2.5]^2 x [-2, 0]` | `upper3d_phi10_{a,b,c,d}`: 276, 390, 976, 2,608 | 20,864 |
| 20 | 3.03, 1.54 | `[-3.5, 3.5]^2 x [-2, 0]` | `upper3d_phi20_{a,b,c,d}`: 291, 429, 951, 2,567 | 20,536 |
| 30 | 4.79, 2.48 | `[-5, 5]^2 x [-3, 0]` | `upper3d_phi30_{a,b,c,d}`: 323, 494, 1,053, 2,850 | 22,800 |

The mechanism of a square footing in 3D is expected to be smaller than the plane-strain one, so the rule is a
conservative sizing, not a claim about the 3D mechanism. The gate does **not** intend wall-confined flow. Mesh
indices are, per case: `a, b, c, d`, then `a` refined once (mesh 4) and `d` refined once (mesh 5). The first four are
independent non-structured sizes; the refined meshes form the nested pairs `(0, 4)` and `(3, 5)`. Section 15.2 lists
the files and their digests.

### 11.1 Kinematic test: rank and residual

Let `E` be the equality rows of the formulation: supports, rollers, rigid footing, isochory at `alpha = 0`
vertices, and zero normal jump on every Bernstein coefficient where `tau = 0`. Let `A` be the set of `(v, V)`
satisfying the boundary and footing rows and the query normalisation.

- **Rank, reported only.** `rank` counts singular values of `E` above `TOL_RANK * sigma_max`.
  `nullity = n_unknowns - rank`. `footing_modes` is the rank of the kernel restricted to `V`. A dense factorisation
  may be skipped when too large; then the three are `null`. Nullity does not enter the verdict: the broken space has a
  large nullity on every mesh, and every `alpha > 0` space has no isochory rows.
- **No footing (amendment 1, item 10).** When the problem has no footing (`n_f = 0`) there is no `V`, so `footing_modes`
  is `null` whether or not `rank` is computed, and the `footing_modes` condition of `formulation_feasible` does not
  apply. A number in its place is refused (`feasibility_rank`). With `n_f >= 1`, `rank`, `nullity` and `footing_modes`
  are computed together or are all `null`.
- **Residual.** With `c_ei` the flow violation at element vertices, `c_FB` the jump violation on Bernstein
  coefficients, `w_ei = V_e / 4` and `w_FB = A_F / (6 h_F)`, `h_F` the longest edge of the facet:

```text
c_ei = |theta|                                  alpha = 0
c_ei = max(0, 3 sqrt2 alpha |e| - theta)        alpha > 0
c_FB = max(0, tau |c_t| - n . c_B)              tau > 0;   |n . c_B| when tau = 0

rho(v, V) = sqrt( sum w_ei c_ei^2 + sum w_FB c_FB^2 ) / sqrt( sum w_ei |D_ei|_F^2 + sum w_FB |c_B|^2 )
rho_feas  = min_A sqrt( sum w_ei c_ei^2 + sum w_FB c_FB^2 ) / min_A sqrt( sum w_ei |D_ei|_F^2 + sum w_FB |c_B|^2 )
```

For the control the facet sums are absent. With every `alpha = 0` both minima are sparse linear least-squares
problems, so `rho_feas` does not depend on a conic solver. With `alpha > 0` the minimisation is conic and no method is
fixed, so `rho_feas` is the witness's own `rho`.

- **Witness.** A result document (section 15.3) of the formulation holding a field `(v, V)` in `A`. A certified Upper
  of the same problem is a witness. The **checker**, not the producer, recomputes `rho` at the witness from the problem
  file with its own operators, and writes the record of section 15.5 with the witness's `result_sha256`.
- **Witness field (amendment 1, item 10).** `rho` is evaluated at the **certified projected field** of the witness: the
  field that the checker's section 17.4 procedure produces from the witness result, after the affine projection of step
  2 and, when a repair is made, the blend of step 4. When that procedure certifies the witness it is the field `x_c` of
  section 17.5. The raw arrays of the result document are **not** the witness field: they satisfy the equality rows only
  to solver accuracy. The projection uses the scales and literal zeros of section 17.1, so the field and its residual
  are reproducible from the problem and result documents alone.

| `status` | Condition |
|---|---|
| `formulation_feasible` | `rho_feas <= TOL_FEAS_OK`; the witness has `rho <= TOL_FEAS_OK`; `footing_modes`, where computed and `n_f >= 1`, is at least 1 |
| `formulation_infeasible` | every `alpha = 0`, and `rho_feas >= TOL_FEAS_LOCKED`. Never claimed for `alpha > 0` |
| `inconclusive` | anything else |

A failure is a **formulation result**, recorded with its numbers per mesh, wall and angle. It may not be hidden
by relaxing a tolerance, by reading `AlmostSolved` as success, or by quoting a conic value whose flow rule is
violated. `formulation_infeasible` says this space yields no Upper bound on that mesh. It says nothing against
the mechanics. Only a record from a real checker can satisfy a gate (`gate_admissible`).

### 11.2 Quality test

For each case and the binding Upper formulation, all of these must hold, using **certified** values only:

1. `kinematic_status` is `formulation_feasible` on every listed mesh;
2. a certified Upper exists on every listed mesh (ADR 10.1, "a certified Upper at every size");
3. a certified bracket exists on mesh 5;
4. its gap is at most `MAX_GAP_3D_VM` for `phi = 0`, or `MAX_GAP_3D_DP` otherwise;
5. `U_fine <= U_coarse + TOL_MONOTONE |U_coarse|` along each nested pair;
6. `L_fine >= L_coarse - TOL_MONOTONE |L_coarse|` along each nested pair.

The gap on every mesh is recorded (`record gap mesh i`) and gates nothing: it is the gap trend the ADR asks for.
Section 13 proves the monotonicity for these formulations. A violation therefore signals a defect or an
unconverged candidate, and is never averaged away.

The limits were set late (section 18.3) and the owner accepted them on 2026-11-19 (OD1). They are closed from this
revision: a limit does not move and no gate mesh changes. A miss is a real P0 failure, recorded.

**Mesh-5 conditioning route (amendment 1, item 8).** At this amendment the checker's capacity is **insufficient** to prove
the conditioning bound of section 17.1 on mesh 5 for two groups: the **Lower** formulation in all eight `upper_gate_3d`
cases, and the **jump Upper** formulation in the von Mises cases (`phi = 0`, fixed and roller walls). On those meshes the
checker writes a capacity or conditioning refusal (section 15.8), not a certificate. So `certified_bracket mesh 5`,
`gap_le ... mesh 5` and `lower_nondecrease 3 5` in every 3D case, and `certified_upper mesh 5` and `upper_nonincrease 3 5`
in the two von Mises cases, **cannot pass** until WP-237d identifies an open-source route that proves the bound at that size
(section 17.1: shape-forced zeros, dense singular values on disjoint row blocks, or an explicit null-space basis with a
stated residual bound, and a proved lower bound on the smallest kept singular value). The gates are **not** weakened: no
mesh, limit, tolerance or gate line moves, a refusal never counts as a pass, and under OD2 it is not a missing certificate
(section 11.3).

### 11.3 Verdicts

The verdict is per case, so **per angle and wall type**, as `evaluate_case` returns it, and `evaluate_program` reads the
verdicts for P0. By owner decision OD2, after the box correction above:

- at `phi = 30` degrees, **only two kinds of failure withhold the 3D Drucker-Prager claim at 30 degrees** instead of
  blocking P0: the gate `gap_le MAX_GAP_3D_DP`, and a missing certificate, meaning a `certified_upper`,
  `certified_bracket`, assertion or monotonicity item that fails because no gate-admissible certificate exists on its
  mesh. The case is then recorded with `claim_withheld`, by angle and wall type;
- at **any** angle, a failed `upper_nonincrease` or `lower_nondecrease` between two certified values, a failed or
  missing `kinematic_status`, or an invalid record keeps P0 open, and any contradiction (certified Lower above certified
  Upper beyond the floors) stops the program;
- at **any** angle, an item whose certificate was replaced by a **capacity or conditioning refusal** (section 15.8)
  keeps P0 open. A refusal is **not** a missing certificate under OD2: it says that the checker could not prove the
  bound, not that no candidate was offered (amendment 1, item 7; orchestrator ruling of 2026-11-19);
- a miss at `phi = 0, 10, 20` degrees is a P0 failure that keeps P0 open (ADR 18.3), and no candidate passing
  triggers K2;
- a failure of an exact formulation case (`exact_on_mesh`, `exact_lower_on_mesh`) or of the thin-slab strip still
  stops P0.

`_withholdable` and `evaluate_program` in the contract test are the executable form.

## 12. P0 Upper formulation, `upper.p2d-jump-bernstein.v0`

Unknowns: an independent P2 velocity `v_(e,a)` on each element (10 nodes, 30 values); `V_f`; `(r_ei, z_ei)` as in
section 10; `(rho_FB, y_FB)` on every interior facet `F` and Bernstein coefficient `B`.

### 12.1 The jump

On an interior facet use the order of `e-` (section 3.4). `n` is the outward normal of `e-`, so it points into
`e+`. The jump is `j = v+ - v-`; `j . n > 0` is opening. Its six nodal values, in facet-node order, are `j_a` at
vertices and `j_ab` at edge midpoints. **A P2 jump is quadratic**, so conditions at facet vertices or nodes do not
cover the facet. Its Bernstein coefficients are

```text
c_a = j_a,      c_ab = 2 j_ab - (j_a + j_b) / 2.
```

Then `j(x) = sum_B B_B(x) c_B` with `B_B >= 0`, `sum_B B_B = 1`, `int_F B_B dA = A/6` for each of the six.

### 12.2 Jump cone and dissipation, derived

`t . j = sigma : sym(j (x) n)`, so the jump uses the same support function with `d = N(n)^T j`. For opening
`o = j . n` and slip `s = |j_t|`: `theta = o` and `|e|^2 = s^2/2 + 2 o^2/3`. Section 4.3 then gives:

```text
root = sqrt(1 - 12 alpha^2),   tau = 3 alpha / root,   k_J = k / root,   requires alpha < 1/sqrt(12)

domain:        o >= tau s
pi_J(j, n) =   k s                  when alpha = 0 and o = 0
               k o / (3 alpha)      when alpha > 0
```

For `alpha >= 1/sqrt(12)`, the only admissible jump is zero. This formulation refuses such a problem with
`formulation_not_applicable`. With the plane-strain match of section 4.4, `tau = tan(phi)` and slip with `o = s tan(phi)`
dissipates `c s`, the associated Coulomb surface. **In 3D this remains the Drucker-Prager jump condition**; it
needs no separate formulation.

### 12.3 Discrete problem

```text
minimise   D_vol + D_jump - W0
D_vol  = sum_e sum_i (V_e / 4) sqrt2 k_e r_ei
D_jump = sum_F sum_B (A_F / 6) k_J rho_FB
subject to
  B_ei v_e = sqrt2 alpha_e r_ei one + Q^T z_ei,   (r_ei, z_ei) in SOC6      each element, 4 vertices
  n . c_FB = tau rho_FB                                                     each facet, 6 coefficients
  T1 . c_FB = y1_FB,  T2 . c_FB = y2_FB,   (rho_FB, y1_FB, y2_FB) in SOC3   each facet, 6 coefficients
  boundary and footing rows, per element node                               section 15.4
  V_fc = 0 (restrained);  W1 = 1 (ray);  V_fc = d_fc (support)
```

`k_J`, `tau` come from the region of `e-`. Every admissible stress meets yield on both sides, so either side is
valid.

- **von Mises** has `tau = 0`. The normal jump is **exactly zero** on all six coefficients, which is zero across
  the facet. The tangential jump dissipates at most `sum_B (A/6) k |c_B|`, which bounds `int k |j| dA` by Jensen.
- **Boundary conditions are strong.** On `fixed`, `symmetry` and `footing` facets the element's own six facet
  nodes satisfy the section 5 and 6 rows. There are no jumps against walls or footings in version 0, which only
  shrinks the admissible set.
- Work rates are those of section 10, integrated per element and per boundary facet.

### 12.4 Why the bound holds

Summing integration by parts over elements, with `n` from `e-` to `e+` and continuous traction `t = sigma n`:

```text
sum_e int_e sigma : D dV + sum_F int_F t . j dA = W0 + lambda W1 .
```

`sigma : D <= pi(D)` and `t . j <= pi_J(j, n)`. On each element Jensen gives section 13's vertex sum. On each facet
`j(x)` is a convex combination of the `c_B`, the cone is convex, so `j(x)` lies in it everywhere, and
`int_F pi_J(j) dA <= sum_B (A/6) pi_J(c_B)` because `pi_J` is convex and positively homogeneous. The check is
conservative: it can refuse a jump that is admissible, and never accepts one that is not.

## 13. Conditions for each bound claim

**Definitions.** `sigma` is statically admissible if it satisfies equilibrium, traction continuity, the traction and
footing conditions, and yield everywhere. `(v, V)` is kinematically admissible if it satisfies the velocity
conditions, its strain rate lies in `dom pi` almost everywhere, and any jump lies in `dom pi_J` on the whole facet.

**Lower claim, `lambda_L <= lambda_max`, holds when:** elements are straight; `k`, `alpha`, `b` are constant per
element; tractions are linear per facet; resultants are integrated exactly; equalities hold; yield holds at four
vertices. Equilibrium and continuity then hold everywhere because every residual is a polynomial fixed by the rows
imposed, and yield holds everywhere by convexity.

**Upper claim, `lambda_max <= lambda_U`, holds when:** elements are straight; `k`, `alpha`, `b` are constant per
element; velocity conditions hold exactly; `W1 > 0`; every vertex strain rate lies in `dom pi`; every Bernstein
coefficient of every jump lies in `dom pi_J`. Virtual work gives `W0 + lambda W1 = int sigma : D dV` plus the jump
term of section 12.4. `D(x)` is linear, so the cone holds throughout, and Jensen with `int L_i dV = V/4` gives
`int pi(D) dV <= sum_i (V/4) pi(d_i)`.

**Monotonicity under nested refinement** (section 3.6). Take a parent field on the fine mesh.

- **Lower.** A parent field is admissible on the children, so the fine maximum is not smaller.
- **Upper.** A parent field is admissible on the children. At child vertices `pi(d(x)) <= sum_i L_i(x) pi(d_i)`.
  Summing `V_c / 4` over child vertices reproduces `int L_i dV = V/4`, so the fine volume sum does not exceed the
  parent's. On a child facet the new Bernstein coefficients are convex combinations of the parent's, with total
  weight `A/6` per parent coefficient, so the jump sum does not grow. New interior facets carry no jump.

**Both are claims about the stated idealisation:** associated flow, perfect plasticity, the stated strength and
interface law. For non-associated real soil they are not bounds on the physical problem.

## 14. Domain claims

`domain_claim` belongs to **each certified bound**. A problem declares `domain.target`, `box` or `halfspace`.
With `halfspace` it also declares `ground` `{axis, coordinate}`, the ground plane; `truncation_sets`, the walls that
cut the half-space; `symmetry_sets`, symmetry planes of the problem itself, which must be `symmetry` sets; and
`exterior_load_free`, true when no load acts outside the mesh. A `box` has no truncation sets, no symmetry sets and
`ground` `null`.

| Value | Meaning |
|---|---|
| `box_stated` | Bound for the stated box with its exact wall conditions. The only value for target `box` |
| `unverified` | Target `halfspace`, extension not established. The bound holds for the box and is an estimate for the half-space |
| `halfspace_by_zero_extension` | The admissible field extends by zero |
| `halfspace_by_extension_elements` | Reserved. Version 0 refuses it because no extension data exists |

**Coverage, required for any half-space claim.** Every exterior facet lies on the ground plane (within
`TOL_GEOM_PLANE`), or belongs to a truncation set, or to a declared symmetry set; and no `traction` or `footing` facet
lies off the ground plane. Otherwise the claim is `unverified`. Without this rule a problem with loaded or free side
faces could claim a half-space bound because its base alone was listed.

`halfspace_by_zero_extension` is then supported only as follows.

- **Upper:** every truncation set is `fixed`. Zero velocity costs no dissipation and no work outside.
- **Lower:** every truncation set is `free`, every `b0` and `b1` is zero, and `exterior_load_free` is true.
  Then zero stress is admissible outside. With self-weight this route is closed.

Roller walls allow tangential motion and never support the Upper claim. A fixed wall supports no Lower
extension. **A bracket's claim is the domain both bounds support:** if both certificates carry
`halfspace_by_zero_extension`, the bracket does. Otherwise it is `box_stated`, and it passes no half-space gate.
In version 0 that first branch **cannot fire**: one problem hash cannot have every truncation set `fixed` (Upper) and
`free` (Lower). Every version 0 bracket is therefore `box_stated`. The rule is kept for extension elements.

## 15. Files, version 0

Every document is UTF-8 JSON, with `NaN` and `Infinity` forbidden, and a number that parses to a non-finite value
(`1e999`) refused after parsing. Unknown keys are refused. Arrays live in one NPZ named by the document's `arrays` key,
relative to the JSON file. Arrays are `<f8` or `<i8` only.

**Names.** A case on mesh `m` has the stem `<case>.m<m>`. Its files are:

```text
<case>.m<m>.problem.json                                    and .npz
<case>.m<m>.<query>.<formulation-id>.result.json            and .npz
<case>.m<m>.<query>.<formulation-id>.program.json           and .npz
<case>.m<m>.<query>.<formulation-id>.certificate.json
<case>.m<m>.<query>.<formulation-id>.certificate.method.json         amendment 1, section 15.9
<case>.m<m>.<query>.<formulation-id>.refusal.json                    amendment 1, section 15.8
<case>.m<m>.<query>.<formulation-id>.feasibility.json
<case>.m<m>.<query>.<lower-formulation-id>__<upper-formulation-id>.bracket.json
```

So two Upper formulations, or six meshes, never collide. A fixture that is not a case uses its own name as `<case>`.

### 15.1 Hash

`sha256` runs over `enc(document without its own hash key)`, then, for each array in name order, the
concatenation `b"a" + enc(name) + enc(dtype) + u64(ndim) + u64(dims...) + raw little-endian C-order bytes`.

```text
enc(null)=b"n"   enc(true)=b"t"   enc(false)=b"f"
enc(number) = b"d" + float64-LE          (integers too)
enc(string) = b"s" + u64(byte length) + UTF-8
enc(list)   = b"l" + u64(count) + items
enc(object) = b"m" + u64(count) + for each key in UTF-8 byte order: enc(key) + enc(value)
```

`u64` is unsigned little-endian. The `arrays` file name is part of the document. Hash keys are
`problem_sha256` in a problem and `result_sha256` in a result. A certificate has no hash key of its own: its
`certificate_sha256` (amendment 1, section 15.9) is this hash of the whole certificate document with no arrays. **Negative zero hashes as zero**: numbers are
encoded as `float64(x) + 0.0`, and float arrays are hashed as `array + 0.0`.

The Kuhn generator (section 3.5) computes vertex `(i, j, k)` as `lo + (hi - lo) * [i/nx, j/ny, k/nz]` in IEEE double,
in that order, and refinement midpoints as `0.5 * (x_a + x_b)`. `tests/data/limit/problem_hashes.json` holds the
expected `problem_sha256` of every active case on every listed mesh, so two builders can compare hashes without
running each other.

Worked digest: the document `{"b": [1, null, true], "a": "x"}` without hash key, with the single array
`v = [0.5, -0.0]` of `<f8`, hashes to `bba0e8c06001808c4f198c1956e7d81423a19daf45369b633cfc14b0c9f12a46`.

### 15.2 Problem: `<stem>.problem.json`, `.npz`

| Key | Content |
|---|---|
| `format`, `format_version` | `ladruno-limit-problem`, `0` |
| `arrays`, `name`, `units` | NPZ file name; free name; `{length, force}` labels |
| `conventions` | the literal object of section 15.7 |
| `domain` | `{target, truncation_sets, symmetry_sets, ground, exterior_load_free}` (section 14) |
| `length_scale` | positive number |
| `regions` | list of `{name, criterion, k, [alpha], b0, b1}` |
| `boundary_sets` | list of `{name, kind, [footing]}`; `footing` is an index, present only for kind `footing` |
| `footings` | list of `{name, x_ref, interface_law}` |
| `queries` | list of `{id, form, restrained, ...}` as in section 7 |
| `problem_sha256` | section 15.1 |

| Array | Shape | Content |
|---|---|---|
| `nodes` | `(N, 3)` | coordinates |
| `tets` | `(E, 10)` | section 3.1 |
| `element_region` | `(E,)` | index into `regions` |
| `facets` | `(Fb, 2)` | section 3.4 |
| `facet_set` | `(Fb,)` | index into `boundary_sets` |
| `t0`, `t1` | `(Fb, 3, 3)` | traction at the facet's 3 vertices; exactly zero except on `traction` facets |

Every boundary set is non-empty and every footing has at least one set.

**From a case to a problem.** `tests/limit_cases.txt` defines the file format of a case. A case on mesh `m`
becomes a problem as follows: all elements take region 0; `target` is `box`; `length_scale` is the longest
box side; exterior facets take the first matching `boundary` line; tractions are evaluated at facet vertices;
`F0` and `F1` default to zero. `build_case_problem` in `tests/test_limit_contract.py` is the executable form.

**The problem `name` (amendment 1, item 4).** `name` is part of the hashed document, so a case problem's name is fixed:
`<case>:<label>`, with `<case>` the case name and `<label>` built from the case's `mesh` or `mesh_file` line for that
mesh index:

- `mesh nx ny nz` gives `<nx>x<ny>x<nz>`, the three counts as decimal integers joined by a lower-case `x`
  (`tresca_strip_vm_thin_slab:8x1x4` for its mesh 0);
- `mesh_file <file> [refine r]` gives `<stem>+<r>`, with `<stem>` the file name up to its **first** `.` and `r` the
  refinement count, `0` when `refine` is absent (`upper3d_fixed_phi00:upper3d_a+0`, `upper3d_fixed_phi00:upper3d_d+1`).

The other derived values follow the same executable form: `units` is `{length, force}` from the first two
`;`-separated parts of the `units` line, each split at its first space; `length_scale` is the longest box side as a
double; `regions` copy the case's region lines with `b0`, `b1` defaulting to `[0, 0, 0]` and `alpha` present only for
Drucker-Prager; `restrained` lists `[footing index, component index]` in case order. `problem_hashes.json` is the check.
A `mesh_file` is an NPZ with `nodes (n, 3)` of `<f8` and `tets (E, 4)` of `<i8`, positively oriented, and
sections 3.5 and 3.6 add refinement and mid-edge nodes. The four committed files are unstructured tetrahedral meshes
of the box `[-2, 2] x [-2, 2] x [-2, 0]`, with the square `[-0.5, 0.5]^2` conforming on `z = 0`, made with gmsh 4.15
(Delaunay 3D, one thread, seed 1). `a`, `b` and `c` use uniform target sizes 1.0, 0.8 and 0.58. `d` grades from 0.08
within 0.1 of the footing edge to 1.0 beyond 1.6. The twelve `upper3d_phi*` files hold the per-angle boxes of section 11,
made by `Ladruno_files/testbed/limit_analysis/make_gate_meshes.py` with gmsh 4.15. All sixteen are data. No test or
implementation depends on gmsh. Each NPZ holds only `nodes` and `tets`.

| File | SHA-256 | Tets |
|---|---|---|
| `upper3d_a.mesh.npz` | `51e1c06981d00932b7ca585eda00c84215f66c24479a23ba64d8cba20ae4cbf3` | 243 |
| `upper3d_b.mesh.npz` | `13985a26835df2df74da6e6d202e7790807282ed7b91652281a00b1dffc7bd75` | 463 |
| `upper3d_c.mesh.npz` | `1e5955642fa0ab339137beef53b66d3ea1bb59dc4099b04085e59c3ec8c8e0af` | 949 |
| `upper3d_d.mesh.npz` | `fb4cf3d23b4aa74f7abe4a0d6602eeefac21815a8a73b69a345a77d99bb7d76a` | 2840 |
| `upper3d_phi10_a.mesh.npz` | `567cf6a3cb7200555b2232dbb9228d69d26c65ea2de4a112c87885779358dad6` | 276 |
| `upper3d_phi10_b.mesh.npz` | `123bc0d36dbd9c8a01f4641359eeaab6feff40c22df1133a090bf18542fe6eb8` | 390 |
| `upper3d_phi10_c.mesh.npz` | `868be0dec382f44b9c60edb94893a171e472521b24801c77f396ff7b71086cbf` | 976 |
| `upper3d_phi10_d.mesh.npz` | `d526b7ee0951d424084028fcf3ec5c0cfd13c61fa5e9835cddc2bd0d44186345` | 2608 |
| `upper3d_phi20_a.mesh.npz` | `f1519598133e4b959a14172b11dd679e1ca4bb61b1aef4913357743e5a59bb4b` | 291 |
| `upper3d_phi20_b.mesh.npz` | `abd4a662b434123a061a630fbec8132f373facb0f4e6e1456fc13215d9c7db94` | 429 |
| `upper3d_phi20_c.mesh.npz` | `ee941c6cb422999ef6be683479f5217b1bcea356f4fc4ed02b6887e3acd81297` | 951 |
| `upper3d_phi20_d.mesh.npz` | `a8e1c3c1b28e8b1705dcd349555eb4384c283d88f42f9a85f439bc182fcf7f80` | 2567 |
| `upper3d_phi30_a.mesh.npz` | `c0647dcbae45efb3f275d081f069a27d2036855fc9d9ac141f1bdba1e1bf90b8` | 323 |
| `upper3d_phi30_b.mesh.npz` | `a841fdc95485e79b83cc150791ce5e1d263f938c969df7d2feee1748aed8c697` | 494 |
| `upper3d_phi30_c.mesh.npz` | `908bdd8bdb6ef761f1d13bbd001b10459ef2ec44a2a8d1b16462c3abadbbb380` | 1053 |
| `upper3d_phi30_d.mesh.npz` | `72005dc6f9c73c6d5a94ab30ef20a6bf5c41a6a5986d8af972fe0c05a3643812` | 2850 |

### 15.3 Result: `<stem>.<query>.<formulation-id>.result.json`, `.npz`

Keys: `format` = `ladruno-limit-result`, `format_version`, `arrays`, `problem_sha256`, `query_id`, `bound`,
`formulation`, `producer` `{name, version}`, `solver` `{name, version, settings}`, `solver_status`,
`objective`, `value`, `wall_time_s`, `result_sha256`.

| `bound` | Required arrays | Optional, never trusted |
|---|---|---|
| `lower` | `sigma (E,4,6)`, `resultant (n_f,6)` | `interior_sigma (E,4,6)` with number `interior_value`; `equilibrium_multiplier (E,3)` |
| `upper`, `upper.p2c-vertex.v0` | `velocity (N,3)`, `footing_velocity (n_f,6)` | `interior_velocity (N,3)` and `interior_footing_velocity (n_f,6)`; `strain_rate (E,4,6)`; `dissipation (E,4)` |
| `upper`, `upper.p2d-jump-bernstein.v0` | `element_velocity (E,10,3)`, `footing_velocity (n_f,6)` | `interior_element_velocity (E,10,3)` and `interior_footing_velocity (n_f,6)`; `strain_rate (E,4,6)`; `dissipation (E,4)` |

**`value` and `objective` (amendment 1, item 2).** `objective` is `c^T x`, the objective of the section 15.4 program at
the producer's solution `x`, with that program's sign. `value` is the bound the field claims:

- Lower: `value = objective_sign * objective`, which is `lambda` for ray forms and `d^T R` for `support`;
- Upper: `value` is the **field value of section 17.4 step 6** evaluated on the stored velocity field, without testing
  the domain: `(D - W0) / W1` for ray forms and `D - W0` for `support`, with `D` the vertex and Bernstein sums of that step.
  It is **not** `c^T x`. The two differ when `W1 = 1` holds only to solver accuracy, or when a cone variable `r` or `rho`
  sits above its minimum. The checker compares `value` with its own step 6 value (17.4 step 1, `value_mismatch`).

`interior_value` is a top-level key, present exactly when `interior_sigma` is. With `solver_status` of
`solved` or `almost_solved`, `value` and the arrays are required. Otherwise `value` and `objective` are
`null` and no field arrays are stored: a time limit gives **no number**. A `lower` result carrying velocity
arrays, or the reverse, is refused. `strain_rate` and `dissipation` are never trusted: a finite `dissipation` charge at a
vertex outside the flow cone is not a repair (section 17.4).

### 15.4 Program: `<stem>.<query>.<formulation-id>.program.json`, `.npz`

Canonical form: minimise `c^T x` subject to `A x + s = b`, `s` in `K`. `K` is the zero cone, then the listed
second-order cones in order. A result's `objective` is `c^T x`; for Lower, `value = objective_sign * c^T x`; for Upper,
`value` is the field value of section 15.3, not `c^T x` (amendment 1, item 2).

Keys: `format` = `ladruno-limit-program`, `format_version`, `arrays`, `problem_sha256`, `query_id`, `bound`,
`formulation`, `objective_sign`, `n_variables`, `n_rows`, `n_zero`, and `cones`, a list of `{dim, count}`. Arrays: CSC
`A_indptr`, `A_indices`, `A_data`; `b`; `c`. There is **no scaling and no presolve** in the program. A solver
adapter may scale or presolve its own copy. Rows may be linearly dependent.

**Lower.** `x = [lambda (ray forms only), sigma_(0,0)...]`, stress in `(e, i, component)` order.
`objective_sign = -1`. `c` is `-1` on `lambda`, or `-d^T R` as a row over `sigma` for `support`.
Rows, in order:

1. equilibrium, 3 per element: `[b1_e | sum N(grad L_i)]`, right side `-b0_e`;
2. interior facets in section 3.4 order, 3 vertices, 3 components: `N(n) sigma_- - N(n) sigma_+ = 0`;
3. exterior facets in `facets` order, 3 vertices: `free` and `traction`, 3 rows, `N(n) sigma - lambda t1 = t0`;
   `symmetry`, rows `T1`, `T2` of the set normal `n_S` (section 5); none for others;
4. ray forms: for each footing and unrestrained component, `R - lambda F1 = F0`;
5. cones, per `(e, i)`: rows `[sqrt2 alpha one^T ; -Q]` on `sigma_ei`, right side `(sqrt2 k, 0)`.

**Upper control.** `x = [v (3N), V (6 n_f), (r, z)_(0,0)...]`. `objective_sign = +1`.
`c = (V_e/4) sqrt2 k_e` on each `r`, minus the `W0` coefficients. Rows, in order:

1. flow, per `(e, i)`, 6 rows: `B_ei v - sqrt2 alpha r one - Q^T z = 0`;
2. `fixed`: for each fixed node ascending, 3 rows `v = 0`;
3. `symmetry`: per set in order, nodes ascending, fixed nodes omitted, `n_S^T v = 0` with the set normal;
4. `footing`: per footing, nodes ascending, 3 rows `v - u - w x (x_a - x_ref) = 0`;
5. restrained components, in list order: `V = 0`;
6. ray forms: `W1 = 1`; `support`: `V_fc = d_fc` over unrestrained `(f, c)`;
7. cones, per `(e, i)`: `-(r, z) + s = 0`. `cones = [{dim: 6, count: 4 E}]`.

**Upper jump.** `x = [v_(e,a,c) (30 E), V (6 n_f), (r, z)_(e,i) (24 E), (rho, y1, y2)_(F,B) (18 F)]`, with `F` the
interior facets of section 3.4 and `B` in facet-node order. `objective_sign = +1`. `c` has `(V_e/4) sqrt2 k_e` on `r`,
`(A_F/6) k_J` on `rho`, minus the `W0` coefficients. Rows, in order:

1. flow, per `(e, i)`, 6 rows, as above with `v_e`;
2. jumps, per `(F, B)`, 3 rows: `n . c_B - tau rho = 0`, `T1 . c_B - y1 = 0`, `T2 . c_B - y2 = 0`. `c_B` uses `e+` values
   minus `e-` values, matched by global node id, and the tangent basis of section 5 for `n`;
3. `fixed`: per `(e, a)` ascending, 3 rows;
4. `symmetry`: per set, per `(e, a)` ascending, those in `fixed` omitted, one row `n_S^T v = 0`;
5. `footing`: per footing, per `(e, a)` ascending, 3 rows;
6. restrained components, then `W1 = 1` or `V_fc = d_fc`;
7. cones: `SOC6` per `(e, i)`, then `SOC3` per `(F, B)`. `cones = [{dim: 6, count: 4 E}, {dim: 3, count: 6 F}]`.

### 15.5 Feasibility record: `<stem>.<query>.<formulation-id>.feasibility.json`

Written by the checker. Keys: `format` = `ladruno-limit-feasibility`, `format_version`, `problem_sha256`, `query_id`,
`formulation`, `status`, `rho_feas`, `numerator`, `denominator`, `n_unknowns`, `n_rows`, `rank`, `nullity`,
`footing_modes`, `method`, `witness_result_sha256`, `witness_residual`, `checker` `{name, version}`, `tolerances`.
`rho_feas = numerator / denominator`; `witness_residual` is `rho` at the witness, and both witness keys are `null`
together. For `alpha > 0`, `rho_feas` equals `witness_residual`. `status` follows section 11.1. `witness_residual` is
evaluated at the certified projected field of the witness (11.1). `footing_modes` is `null` when the problem has no
footing, and `rank`, `nullity`, `footing_modes` are `null` together when the factorisation is skipped (amendment 1).

### 15.6 Certificate and bracket

`<stem>.<query>.<formulation-id>.certificate.json`, written by the checker only:

| Key | Content |
|---|---|
| `format`, `format_version` | `ladruno-limit-certificate`, `0` |
| `problem_sha256`, `result_sha256`, `query_id`, `bound`, `formulation` | identity |
| `classification` | `certified_bound` or `estimate` |
| `reasons` | codes of section 17.6; empty exactly when certified |
| `value_candidate`, `value_certified`, `floor` | numbers; the last two are `null` for an estimate |
| `domain_claim`, `criteria`, `flow_rule`, `interface_laws` | section 14; sorted lists; the literal `associated` |
| `checks` | metrics of section 17.7; a `certified_bound` must show the numbers listed there, within their limits |
| `checker` | `{name, version}` |
| `tolerances` | the table of section 17.8 |

`<stem>.<query>.<lower-id>__<upper-id>.bracket.json`: `format` = `ladruno-limit-bracket`, `format_version`,
`problem_sha256`, `query_id`, `lower_formulation`, `upper_formulation`, `lower_result_sha256`, `upper_result_sha256`,
`lower`, `upper`, `floor_lower`, `floor_upper`, `status`, `gap`, `domain_claim`, `lower_domain_claim`,
`upper_domain_claim`. The two formulation keys must be those of the two certificates.

Two checker names are reserved and no gate accepts them: `schema-fixture`, for hand-made fixtures, and
`procedure-demo`, for the demonstration of section 18.3.

**Gate certificates need a method report (amendment 1, item 5).** A certificate written under any other checker name is a
*gate certificate*. It is gate-admissible only with a valid method report (section 15.9) bound to it by
`certificate_sha256`. `gate_admissible`, the bracket validator and `evaluate_case` all refuse a gate certificate without
one; the bracket validator refuses the bracket with `method_report_missing`, or with the report's own code when the report
is invalid.

### 15.7 Literal conventions object

```json
{"sign": "tension_positive", "voigt": ["xx", "yy", "zz", "xy", "yz", "zx"],
 "stress_shear": "tensor", "strain_rate_shear": "engineering", "tet10": "vtk",
 "facets": "opposite_vertex_outward", "footing_resultant": "footing_on_soil"}
```

### 15.8 Capacity refusal: `<stem>.<query>.<formulation-id>.refusal.json` (amendment 1, item 7)

Written by the checker **instead of** a certificate when it cannot finish the procedure of section 17 for lack of
capacity, or cannot prove the conditioning bound of section 17.1. Keys: `format` = `ladruno-limit-refusal`,
`format_version`, `problem_sha256`, `result_sha256`, `query_id`, `bound`, `formulation`, `kind`, `stage`, `detail`,
`resources`, `checker` `{name, version}`, `tolerances`.

| Key | Content |
|---|---|
| `kind` | `capacity`: the computation did not fit in the checker's memory or time. `conditioning`: it ran, and no bound of section 17.1 was proved |
| `stage` | `projection`, `dependent_rows` or `kappa_bound`: the step that could not be completed |
| `detail` | free text naming the route tried and why it failed |
| `resources` | `{n_rows, n_unknowns, nnz}` of the reduced system, non-negative integers |

Identity keys must match the result document, exactly as for a certificate. A refusal carries no value and no
classification. A result with both a certificate and a refusal is an invalid record set (`refusal_conflict`).

**Meaning.** A refusal is **not** an estimate of the bound and **not** a missing certificate. A case item whose mesh and
formulation have a refusal from a gate checker, and no gate-admissible certificate, fails with cause `refused`. Under OD2
that cause is never withholdable, so it keeps P0 open at every angle (section 11.3). `evaluate_program` lists the cases
with refused items.

### 15.9 Method report: `<stem>.<query>.<formulation-id>.certificate.method.json` (amendment 1, item 5)

Written by the checker beside every certificate, and **required** for a gate certificate (15.6). Keys: `format` =
`ladruno-limit-certificate-method`, `format_version`, `problem_sha256`, `result_sha256`, `certificate_sha256`, `query_id`,
`bound`, `formulation`, `checker`, `kappa`, `dependent_rows`, `rows_skipped`, `interior_field`, `tolerances`.

| Key | Content and rule |
|---|---|
| identity | `problem_sha256`, `result_sha256`, `query_id`, `bound`, `formulation` and `checker` equal the certificate's. `certificate_sha256` is the section 15.1 hash of the certificate document |
| `kappa` | `{method, sigma_max_upper, sigma_min_lower, kappa_upper, perturbation_bound, residual}`. `method` is `dense_svd`, `dense_svd_blocks`, `eigen_normal_rows` (`A_hat A_hat^T`) or `eigen_normal_columns` (`A_hat^T A_hat`). `sigma_max_upper` is a proved upper bound on `sigma_max(A_hat)`; `sigma_min_lower > 0` a proved lower bound on the smallest kept singular value; `perturbation_bound >= 0` the stated bound used for them. `residual` is the Ritz or eigenvalue residual bound, required for the two eigenvalue methods and `null` for the two singular-value methods. `kappa_upper` equals the certificate's `checks.kappa_eq` and is at least `sigma_max_upper / sigma_min_lower` |
| `dependent_rows` | `{n_rows, n_unknowns, rank, proofs}` of the reduced system. `proofs` is a list of `{route, count, residual_bound}`: `route` is `shape_forced_zeros`, `dense_svd_blocks` or `null_space_basis`; `count >= 1`. `residual_bound` is `null` for `shape_forced_zeros`, the stated factorisation perturbation bound for `dense_svd_blocks`, and the bound on the null-space basis residual for `null_space_basis`, which must not exceed `TOL_RANK * sigma_max_upper`. `rank = n_rows - rows_skipped` |
| `rows_skipped` | the number of rows whose zero singular values were skipped. It equals the sum of the proof counts: **no row is skipped without a proof** |
| `interior_field` | `{source, sha256}`. `source` is `none` (no repair; requires `checks.repair_mu = 0`), `zero` (Lower only, under 17.3 step 5), or `result` (the result's interior arrays). `sha256` is `null` unless `source` is `result`; then it is the section 15.1 hash of the empty document `{}` with the interior arrays of the bound (`interior_sigma`; or `interior_velocity` or `interior_element_velocity`, with `interior_footing_velocity`) |
| `tolerances` | the table of section 17.8 |

`validate_method_report` in `tests/test_limit_contract.py` is the executable form. A report that does not match its
certificate is invalid, and the certificate is then not gate-admissible.

## 16. Derivations behind the cases

### 16.1 Uniform cube

Uniform stress lies in the Lower space and linear velocity in every Upper space, so both are exact on every mesh.
With `sigma = s e_x e_x`: `sqrt(J2) = |s|/sqrt3`, `I1 = s`, so tension fails at `k / (1/sqrt3 + alpha)` and
compression at `k / (1/sqrt3 - alpha)`. At the apex, `sigma = p I` and `p = k / (3 alpha)`, with Upper field
`v = x`. For von Mises, `v = (x, -y/2, -z/2)` gives `pi = sqrt3 k`, and `W1 = 1`. A dead traction `f sqrt3 k` leaves
`lambda = (1 - f) sqrt3 k`.

### 16.2 Rigid plate in shear

`sigma_xy = lambda` only, with `t1 = (0, +-1, 0)` on the `x` faces. Upper field `v = (y, 0, 0)` with plate velocity
`(1, 0, 0)`: the side traction does no work, `pi = k`, `lambda = k`. With `x_ref` raised by `h` above the plate
centre, `M = (0, 0, h) lambda A`, which is the offset case.

### 16.3 Reference body force on a shear block

`sigma_xy = lambda gamma (H - y)` is linear, so the Lower field is exact on every mesh. Yield at the base gives
`lambda = k / (gamma H)`. The Upper field `v_x = min(y/h, 1)`, with `h` the first layer, gives
`lambda <= k / (gamma (H - h/2))`, the `shear_band` formula.

### 16.4 Strip on a thin slab

The slab uses `x` across the footing, `y` through the thickness and `z` upward. `symmetry` on both `y` faces gives
plane strain: a plane-strain velocity is admissible in the slab, and a plane-strain stress field extended with the
optimal `syy` of section 4.4 is admissible too, so the two collapse loads coincide.

**Lemma (box and half-plane).** Let a bonded strip of width `B` act on a weightless half-plane with collapse load
`lambda_inf`. Let the box have **fixed** walls at `x = +-W` and `z = -D`, with the same ground conditions inside
it. Suppose (a) a statically admissible field exists for the half-plane at `lambda_inf`, and (b) a collapse mechanism
at `lambda_inf` has support in the closed box. Then `lambda_box = lambda_inf`.

*Proof.* The mechanism vanishes on and beyond the walls, so it is admissible for the box and `lambda_box <= lambda_inf`.
A fixed wall imposes no stress condition, so the half-plane field restricted to the box is admissible and
`lambda_box >= lambda_inf`. A symmetric half model follows by averaging a field with its mirror image, which is
admissible by convexity. *End.*

The gate compares box bounds with the half-plane value only under this lemma. Rollers lose the second inequality.
The Lower assertion uses (b) only, and the Upper assertion uses (a) only.

**Premise (a) is unverified for `phi = 20` degrees.** For Tresca, Bishop (1953) extends the Prandtl field to the
half-plane. For the frictional strips no source with a complete admissible stress field has been read; Shield (1955),
J. Mech. Phys. Solids 4, 10-16, is a recalled candidate. `upper_ge_target` on both Drucker-Prager strips therefore
rests on an unread source, and the case files say so. The Tresca strip cites Prandtl (1920) for the mechanism and
Prandtl (1921), the ADR's [Pra21], for the value; neither text has been read.

**Extents,** footing half-width `b = B/2`, `beta = pi/4 + phi/2`, measured from the **centreline**:
`r0 = b / cos(beta)`, `r1 = r0 exp((pi/2) tan(phi))`. The Prandtl mechanism reaches `x = b + 2 r1 cos(pi/4 - phi/2)` and
depth at most `r1`. For `phi = 0` this is `x = 1.5 B` and depth `B / sqrt2`, the frozen minimum. For `phi = 20 deg`,
`b = 1`: `x = 6.06`, depth at most `3.09`. The gate boxes have `W = 4`, `D = 2` for `phi = 0` and `W = 8`, `D = 4` for
`phi = 20 deg`, with `B = 2`.

**Values.** The slip-line invariant along the fan gives `Nq = exp(pi tan(phi)) tan^2(beta)` and
`Nc = (Nq - 1) cot(phi)`, tending to `2 + pi`. A bonded footing gives the same value as a smooth one: the wedge
moves with the footing, and a smooth-footing stress field is admissible for a bonded footing.

**Gate.** On the listed finest mesh, the certified bracket has gap at most `MAX_GAP_STRIP` for the Tresca strip and
`MAX_GAP_PRANDTL` for the Prandtl cases, and refines monotonically along the nested pair.

### 16.5 Gravity

**Rigid-motion work (L0).** For `v = u + w x (x - x_p)`, `int b . v dV = b . (u V + w x S)`, with `S` the first moment
about `x_p`. Case `gravity_work_rigid_motion` checks the volume weights. Every symmetric rule integrates a rigid motion,
so wrong quadrature is caught by the quadratic field of the negative fixture instead.

**Dead-load shear layer, exact.** Layer `0 <= z <= H`, fixed base, dead body force `(gamma, 0, 0)`, top traction
`(lambda, 0, 0)`, roller `y` faces, side tractions `+-(0, 0, lambda + gamma (H - z))`. The field
`sigma_zx = lambda + gamma (H - z)` is in equilibrium and is linear. It first yields at the base, so
`lambda* >= k - gamma H`. The Upper field `v_x = min(z/h, 1)` has no lateral strain, so it joins the rigid layer above
it without a jump. It gives `lambda <= k - gamma (H - h/2)`. Letting `h -> 0`, **`lambda* = k - gamma H`**.

**Weighted Drucker-Prager column, exact.** Fixed top, roller sides, body force `(0, 0, -gamma)`, base traction
`(0, 0, -lambda)`. The hydrostatic field `sigma = (lambda + gamma z) I` is admissible up to the apex at the top. The Upper
band `v_z = -min((H - z)/h, 1)` is one-dimensional, and for `2 sqrt3 alpha <= 1` its strain lies in the flow cone. It gives
`lambda <= k/(3 alpha) - gamma (H - h/2)`. Hence **`lambda* = k/(3 alpha) - gamma H`**.

**Von Mises column: a sandwich, not an exact value.** Free sides, smooth base, top pressure `lambda`, unit weight `gamma`.

```text
sqrt3 k - gamma H   <=   lambda*   <=   sqrt3 k - gamma H / 2
```

- Lower field: `sigma_zz = -lambda - gamma (H - z)` alone. It reaches yield only at the base.
- Upper field: uniform compression, `v = (x/(2H), y/(2H), -z/H)`.

The case asserts all four sides: `L >= sqrt3 k - gamma H` and `U <= sqrt3 k - gamma H/2` from the explicit fields, and
`L <= sqrt3 k - gamma H/2` and `U >= sqrt3 k - gamma H` from the sandwich. The third kills a Lower zero-gravity mutant,
whose value `sqrt3 k = 1.732` passes the first.

**`sqrt3 k - gamma H` is not the collapse load of a finite column.** The material above the base has spare strength, so
self-equilibrated lateral stress can confine the base: horizontal equilibrium needs only `int sigma_xx dz = 0` on each
vertical section. A rigid block above a thin plastic layer is not available either, because the layer's lateral flow
would need a velocity jump that grows as the layer thins. Uncommitted check, 2026-11-18, Clarabel 0.11.1, 768 tets,
`k = 1, gamma = 0.25, H = 2`: `estimate: L=1.360, U=1.3835`, against the static value `1.232`.

**Which case kills the zero-gravity mutant.** `vm_dead_gravity_shear_layer` is the designated killer. With gravity zeroed
its Upper value rises to `k`, above `upper_le_formula`, and its Lower value falls below the target. Uncommitted check on
24 tets: `estimate: L=0.247, U=0.992` against `L=0.700, U=0.725` with gravity. `dp_weighted_apex_column` and
`vm_weighted_column_sandwich` also fail their Upper limits.

**Level-ground von Mises footings are insensitive to gravity.** On a level ground `z = 0` with fixed or roller walls,
add `gamma z I` to any admissible stress: it balances the weight, leaves `J2` unchanged, adds no shear on rollers, and
vanishes on the ground and footing. For an isochoric field, `int v_z dV = sum of boundary terms of z v . n`, which vanish on
the ground (`z = 0`), fixed walls (`v = 0`) and rollers (`v . n = 0`), and across facets with zero normal jump. So weight
changes **neither bound**. `tresca_strip_weighted_insensitive` records this. **No level-ground von Mises footing case can
detect a gravity defect.**

## 17. Certificate

The checker works from the problem and result only, rebuilds every operator, and never reads solver cones or solver
tolerances. It always runs **affine projection, then inequality repair, then a full recheck of the final field.**

### 17.1 Terms

`A x = g` is the affine system of the bound, with the rows of sections 17.3 and 17.4.

- **Scales.** Each unknown `j` has a scale `d_j >= |x_j|`. Lower: `S_ref` (section 17.7) for every stress component,
  and `|value|` for `lambda` (1 if the value is zero). Upper: `v_ref`, the largest absolute nodal or element-node
  velocity and footing translation, for every velocity and translation, and `max(v_ref / ell, max |w|)` for footing
  rotations. Scaled unknowns `x_j / d_j` make every step below independent of units.
- **Literal zeros.** A homogeneous row (`g_i = 0`) with one non-zero entry fixes its unknown to zero. Mark that unknown,
  drop it from the other rows, and repeat until no homogeneous row with one remaining entry is left. The marked unknowns
  are set to the literal `0.0` before projection and stay zero, so those rows hold exactly. Fixed nodes are the common
  case.
- **Reduced system.** The remaining rows and unknowns. `A_hat` is its matrix times `diag(d)`, with each row divided by
  its largest absolute entry.
- **Backward error.** `r = A x - g` and

  ```text
  omega = max_i |r_i| / ( sum_j |a_ij| d_j + |g_i| )
  ```

  over all rows, evaluated in double precision. A row with zero denominator must have `r_i = 0`. So each row holds up to
  `omega` times its own natural size: `x` solves exactly a system whose right side `g_i` moves by at most
  `omega (sum_j |a_ij| d_j + |g_i|)`. A componentwise form, with `|x_j|` in place of `d_j`, is **not** used: a homogeneous
  row whose entries are all rounding-small has residual equal to its denominator, and no solver field could pass.
- **Conditioning: a proved bound (amendment 1, item 6; replaces the frozen rule).** With `n_dep` the number of dependent
  rows **proved** below and `r = n_rows - n_dep`, let `sigma_r` be the `r`-th largest singular value of `A_hat`. The
  reported `kappa_eq` must be a **proved upper bound** on `sigma_max(A_hat) / sigma_r`: `kappa_eq >= sigma_max_upper /
  sigma_min_lower`, with `sigma_max_upper >= sigma_max(A_hat)` and `0 < sigma_min_lower <= sigma_r` each proved, and the
  perturbation and residual bounds used stated in the method report (15.9). An estimate, including an unbounded Lanczos or
  Ritz value, is not a bound and certifies nothing.
  - *Dependent rows* count, and their zero singular values are skipped, **only** when proved by one of:
    1. **shape-forced zeros**: after literal zeros a row has no remaining entry. Its right side must be zero; otherwise
       the system is inconsistent and fails at the recheck;
    2. **dense singular values on disjoint row blocks**: the rows are split into blocks whose column supports are
       pairwise disjoint, so the singular values of `A_hat` are the union of the blocks' values. Each block gets a dense
       singular value decomposition with a stated perturbation bound `delta`. A value counts as zero only when it is at
       most `TOL_RANK * sigma_max_upper`, and the smallest kept value minus `delta` exceeds `TOL_RANK * sigma_max_upper +
       delta`, so the factorisation error cannot move a value across the split. One block, the whole matrix, is allowed;
    3. **an explicit null-space basis**: `n_dep` left null vectors `Z` with stated bounds on `|Z^T A_hat|_2` and on the
       departure of `Z` from orthonormality, which together bound the skipped singular values by at most
       `TOL_RANK * sigma_max_upper`.
  - *Eigenvalue routes.* `A_hat A_hat^T` and `A_hat^T A_hat` are both allowed. Their non-zero eigenvalues are the
    `sigma_i^2`. `A_hat A_hat^T` has `n_dep` zero eigenvalues and `A_hat^T A_hat` has `n_unknowns - r`; a skipped count
    must match. A Ritz or eigenvalue residual gives `sigma_min_lower` only together with a proof that no unskipped
    eigenvalue lies below it, such as the full spectrum of a dense block or a counted inertia; the residual is reported.
    A shift used on these singular matrices is **positive**.
  - *No proof, no number.* When no route proves `n_dep`, or no lower bound on `sigma_r` is proved within the checker's
    capacity, the checker writes a **capacity or conditioning refusal** (section 15.8) and no certificate. It never
    reports an estimated `kappa_eq`, and it never skips a row it has not proved dependent.
  237c reports `kappa_eq`, or the refusal, on every gate mesh before any gate is read.
- **Move.** `projection_move = max_j |x_P,j - x_j| / d_j`.

### 17.2 Equality acceptance

```text
tau_eq = TOL_EQ_ROUNDING_FACTOR * u,     u = 2^-53            (7.1e-15)

accept only if   omega <= tau_eq,   kappa_eq <= TOL_KAPPA_CAP   and   kappa_eq * omega <= TOL_FLOOR_MAX
```

The threshold is a rounding level. Evaluating a row of at most about 60 entries in double precision is accurate to a
few dozen `u` of its natural size whatever the conditioning, so `tau_eq` carries **no** `kappa` factor. Conditioning
enters the floor only (section 17.5): `kappa_eq` measures how far a rounding-level residual can move the value. The cap
is a number: `TOL_KAPPA_CAP = 1e8`, and with `omega <= tau_eq` it implies `kappa_eq * omega <= 7.1e-7`. `kappa_eq` is the
proved bound of section 17.1. A proved bound above the cap gives an estimate with `equality_ill_conditioned`; a bound that
cannot be proved gives a refusal (15.8), not an estimate; a system with `omega > tau_eq` gives
`equality_residual_above_floor`. No solver setting enters either test. The demonstration of section 18.3 reached
`omega` between 3e-17 and 7e-17 on solver fields.

### 17.3 Lower procedure

1. Check hashes, query, formulation, straightness and volumes. For ray forms the candidate `lambda` **is** `value`;
   recompute `resultant` from `sigma` with section 6 and compare it with the stored array within
   `TOL_VALUE_MATCH * S_ref * A_f` (moments with an extra `ell`). For `support`, recompute `value = d^T R` and compare
   within `TOL_VALUE_MATCH * |value|`. A mismatch is `value_mismatch`.
2. **Project.** `x = (lambda, sigma)`, or `sigma` for `support`. Rows: equilibrium, facet continuity, free and traction
   boundaries, symmetry rows with the set normal (section 5), and unrestrained footing resultants. Apply the literal zeros, then the minimum-norm
   correction of the scaled unknowns in `A_hat`, with iterative refinement until `omega` stops decreasing (at most five
   sweeps). Require `projection_move <= TOL_PROJECTION_MOVE`. An inconsistent system keeps its residual and fails at
   step 6.
3. Evaluate `g_ei = sqrt(J2) + alpha I1 - k` from **principal stresses**. If `max g / S_ref > TOL_CONE_ACCEPT`, stop with
   `cone_violation_too_large`.
4. **Repair inequalities**, if any `g_ei > -TOL_CONE_MARGIN S_ref`, with a strictly admissible field `(sigma_int, value_int)`
   satisfying the same equalities and `g_int <= -TOL_INTERIOR_MIN S_ref`:
   `mu = max (g + m) / (g - g_int)`, `m = TOL_CONE_MARGIN S_ref`, `sigma_new = (1 - mu) sigma + mu sigma_int`,
   `value_new = (1 - mu) value + mu value_int`. Require `mu <= TOL_REPAIR_MAX`.
5. `sigma_int = 0`, `value_int = 0` is available only when `b0`, `t0` and `F0` are all zero and every `k > 0`. Otherwise it
   must be supplied as `interior_sigma` and is itself projected by step 2 and checked. If none exists, the result stays
   an estimate with `no_interior_field`; an invalid one gives `interior_field_invalid`.
6. **Recheck the final field from scratch:** every literal-zero unknown is exactly `0.0` (else `literal_zero_violated`),
   `omega`, `kappa_eq`, every traction and resultant row, and `g_ei <= 0` from principal stresses.

### 17.4 Upper procedure

1. As Lower step 1. The value is recomputed from the candidate with the formulas of sections 4.3 and 12.2 **without
   testing the domain**: `sqrt2 k |e|` or `k theta / (3 alpha)` per vertex, `k s` or `k o / (3 alpha)` per Bernstein
   coefficient. The domain is step 3's business. Stored `strain_rate` and `dissipation` arrays are ignored.
2. **Project.** `x = (v, V)`. Rows: fixed nodes, roller nodes with the set normal (section 5), rigid footprint, restrained components, `V = d` for
   `support`, `theta = 0` at every `alpha = 0` element vertex, and `n . c_B = 0` on every coefficient where `tau = 0`.
   Work normalisation is **not** a row (section 18.2, R9). Literal zeros and the scaled minimum-norm correction as in
   Lower step 2. Require `projection_move <= TOL_PROJECTION_MOVE`.
3. Compute `d_ei` and `c_B` with the checker's own operators, and the **signed** conditions

   ```text
   c~_ei = 3 sqrt2 alpha |e_ei| - theta_ei          each vertex with alpha > 0
   c~_FB = tau |c_t| - n . c_B                      each coefficient with tau > 0
   ```

   which are convex and positively homogeneous; admissible means `c~ <= 0`.

   **`D_scale` has one value per kind of condition, never a mixture (amendment 1, item 10).** A volume condition has the
   units of a strain rate (1/time) and a jump condition those of a velocity (length/time), so they are scaled separately:

   ```text
   D_scale_vol  = max over (e, i) of |D_ei|_F          strain-rate units; |D|_F^2 = |d_diag|^2 + (d_xy^2 + d_yz^2 + d_zx^2)/2
   D_scale_jump = max over (F, B) of |c_B|             velocity units; Euclidean norm of the coefficient vector
   rho_flow      = max(0, max c~_ei) / D_scale_vol
   rho_jump_cone = max(0, max c~_FB) / D_scale_jump
   ```

   A ratio whose scale is zero is `0`: every condition of that kind then vanishes. Wherever a `D_scale` appears below, it
   is the scale of the condition's own kind. If either ratio exceeds `TOL_FLOW_ACCEPT`, stop with
   `flow_violation_too_large`.
4. **Repair inequalities** only by an admissible construction. `pi` and `pi_J` are `+infinity` outside their domains, so a
   finite dissipation charge is not a repair. Let `m = TOL_CONE_MARGIN * D_scale` (per kind) and `P` the conditions with
   `c~ > m`.
   - If `P` is empty, nothing is repaired. This is always the case for von Mises: after projection no inequality remains.
   - Otherwise blend with a supplied field `(v_int, V_int)` (`interior_element_velocity` or `interior_velocity`, with
     `interior_footing_velocity`), projected by step 2 with its own scales. It must have `D_scale_int > 0` for each kind
     that has a condition in `P`, satisfy `c~_int <= m_int` at every condition and
     `c~_int <= -TOL_INTERIOR_MIN * D_scale_int` (same kind) at every condition of `P`,
     else `interior_field_invalid`. Then `mu = max over P of c~ / (c~ - c~_int)`, `x_new = (1 - mu) x + mu x_int`,
     `mu <= TOL_REPAIR_MAX`. By convexity each condition of `P` ends at or below 0 and every other one at or below its
     margin.
   - Conditions that the rows force to zero, such as a vertex of an element with two fixed facets or a coefficient
     between two rigidly moving facets, have `c~ = 0` to rounding, are never in `P`, and need no clearance.
   - A producer can obtain an interior field by maximising one bounded slack per condition, then re-raising the slacks
     still at zero while holding the others, until none moves: the result lies in the relative interior.
   - If no interior field is supplied, the result stays an estimate with `no_admissible_correction`.
5. Require `W1 >= TOL_WORK_MIN * W_scale` and `W1 > 0` for ray forms, with `W_scale` the sum of absolute work terms;
   otherwise `work_nonpositive`.
6. `D = sum (V_e/4) pi(d_ei) + sum (A_F/6) pi_J(c_B)`. Ray forms: `U = (D - W0) / W1`, so no exact normalisation is needed.
   `support`: `U = D - W0`.
7. **Recheck the final field from scratch:** every literal-zero unknown exactly `0.0` (else `literal_zero_violated`),
   `omega`, `kappa_eq`, all conditions of step 2,
   `c~ <= TOL_CONE_MARGIN * D_scale` (per kind) at every `alpha > 0` condition, and `U`.

### 17.5 What "certified" means

A **certified bound** means all of the following held for the final field `x_c`, as evaluated in double precision by the
independent checker:

1. every affine row holds with residual at most `omega (sum_j |a_ij| d_j + |g_i|)`, `omega <= tau_eq`, and the rows of
   literal zeros hold exactly;
2. every inequality passes on the final field: Lower `g <= 0`; Upper `c~ <= TOL_CONE_MARGIN * D_scale`, with the
   `D_scale` of the condition's kind (17.4 step 3). The margin is a
   rounding allowance for evaluating the cone, a few thousand `u` of the field's size, and is unrelated to the floor.
   No flow or jump condition is merely charged for;
3. the value is the theorem's value for `x_c`, recomputed by the checker;
4. the reported `floor = max(FLOOR_MIN_REL, kappa_eq * omega) * |value_certified|` estimates, to first order, how far
   moving `x_c` onto the exact affine set could move the value. The estimate omits the norm-equivalence factor between the
   row-wise `omega` and the 2-norm in which `kappa_eq` acts; a rigorous version multiplies it by at most
   `sqrt(n_rows) * max_i nnz_i`.

It is **not** an interval-arithmetic proof. "Rigorous" means rigorous for the stated idealisation, up to this floor.

### 17.6 Outcome codes and negative fixtures

Reason codes: `no_candidate`, `hash_mismatch`, `query_mismatch`, `formulation_not_admitted`, `formulation_not_applicable`,
`curved_element`, `degenerate_element`, `boundary_violation`, `projection_move_too_large`, `equality_residual_above_floor`,
`equality_ill_conditioned`, `cone_violation_too_large`, `no_interior_field`, `interior_field_invalid`, `repair_too_large`,
`flow_violation_too_large`, `no_admissible_correction`, `work_nonpositive`, `value_mismatch`, `unsupported_domain_claim`,
`literal_zero_violated`.

`tests/data/limit/negative_fixtures.json` lists the fixtures for every false-green item of the harness and every row of
ADR 11.5 (the test maps each ADR row to its fixtures). `stage` is `schema` (the format validators refuse it), `checker` (the
full procedure must end as an estimate with the **one** `expect` code) or `final_check`. Each checker fixture is built so
that every earlier step passes: the stored value or resultant is that of the mutated field unless the defect is the value
itself, and the contract test asserts that the named defect is present. A `final_check` fixture goes straight to the
recheck, **without projection**, and must be refused with `equality_residual_above_floor`: one equilibrium row sits at
three times `tau_eq` in the `omega` of section 17.1, with `S_ref = 1`.

The jump fixture uses the ADR 10.2 pattern on a strictly dilating field with no equality rows: one side of an interior
plane opens by `delta n` at element vertices only, so every nodal jump is admissible and the edge Bernstein coefficients are
`-delta n`. A checker that tests jumps only at facet nodes passes it.

Fixtures added by amendment 1 (item 9):

- `normal_jump_isochoric` replaces `normal_jump_one_node`, which also broke isochory. It moves one mid-edge node of one
  element along `grad L_a x grad L_b`, the direction orthogonal to both vertex gradients of its edge `(a, b)`. Then
  `theta` stays zero at all four vertices (section 3.3), no boundary row is touched, and the only violated equality is the
  normal jump `n . c_B = 0` on the two interior facets through that edge;
- `vm_jump_opening` translates one side of an interior plane of the von Mises cube rigidly, so every Bernstein coefficient on
  that plane **opens** (`n . c_B > 0`) with no strain and no boundary row touched. A checker that reads the von Mises jump
  condition as the cone `n . c_B >= tau |c_t|` with `tau = 0`, instead of the equality of section 12.3, passes it;
- `restrained_footing_component` adds `delta y e_z` to the plate-shear field of a problem whose query restrains the
  plate's `uz`. The field stays isochoric, zero on the fixed base and rigid on the footprint with `V_uz = delta`, so the
  only violated row is the restrained component `V_uz = 0`.

Each expects `projection_move_too_large`.

### 17.7 Scaled metrics

`L_ref` is the bounding-box diagonal.
`S_ref = max(max k, max |sigma|, max |t0|, |value| max |t1|, L_ref max(|b0| + |value| |b1|), max_f |F_f| / A_f)`, using the
components of `sigma` for Lower and omitting that term for Upper. For `support`, the `|value|` terms are absent.

`checks` holds these keys, `null` where not applicable. Diagnostics of the candidate before projection: `rho_div`
(`|div sigma + b|_inf h_e / S_ref`), `rho_jump`, `rho_bc`, `rho_resultant` (force over `S_ref A_f`, moment over
`S_ref A_f ell`). Final-field values: `projection_move`, `omega_eq`, `kappa_eq`, `tau_eq`, `rho_cone`, `repair_mu`,
`straightness_max`, `rho_flow`, `rho_jump_cone`, `work_reference`.

A `certified_bound` must show, as finite numbers within their limits: `omega_eq <= tau_eq`; `tau_eq` equal to section
17.2's value; `1 <= kappa_eq <= TOL_KAPPA_CAP` with `kappa_eq * omega_eq <= TOL_FLOOR_MAX`;
`projection_move <= TOL_PROJECTION_MOVE`; `repair_mu <= TOL_REPAIR_MAX`; `straightness_max <= TOL_GEOM_STRAIGHT`. Lower:
`rho_cone = max g / S_ref <= 0`. Upper: `rho_flow <= TOL_CONE_MARGIN`; for the jump formulation
`rho_jump_cone <= TOL_CONE_MARGIN`; for ray forms `work_reference = W1 > 0`. And
`floor == max(FLOOR_MIN_REL, kappa_eq * omega_eq) * |value_certified|`. `certified_check_failures` is the executable form.

### 17.8 Named tolerances

Values are provisional until 237c and 237d calibrate them on the exact cases. A value changes only by amendment, with
evidence, and **never to rescue a run**.

| Name | Value | Rationale |
|---|---|---|
| `TOL_GEOM_STRAIGHT` | 1e-9 | Relative to `h_e`. Admits meshes written with 12 significant digits at `L/h` up to 1e3 |
| `TOL_GEOM_VOLUME` | 1e-10 | `V_e >= tol h_e^3`. A regular tetrahedron has 0.118; below this, gradients amplify rounding by 1e10 |
| `TOL_GEOM_PLANE` | 1e-9 | Normal deviation and offset, relative to `L_ref` |
| `TOL_RANK` | 1e-10 | Separates numerical zero singular values of row-scaled operators from the smallest real ones |
| `TOL_EQ_ROUNDING_FACTOR` | 64 | Multiple of `u` in `tau_eq`. A row has at most about 60 entries, so its evaluation error is about that many `u` of its natural size |
| `TOL_FLOOR_MAX` | 1e-6 | Largest first-order value effect `kappa_eq * omega` a certificate may carry |
| `TOL_KAPPA_CAP` | 1e8 | Largest `kappa_eq` certified. With `omega <= tau_eq` the floor stays below 7.1e-7. Measured 45 to 190 on the demonstration meshes |
| `TOL_PROJECTION_MOVE` | 1e-6 | A certificate speaks about this candidate. A larger move means a different field |
| `TOL_VALUE_MATCH` | 1e-6 | The result's value must be the value of its own field. Catches normalisation and quadrature faults |
| `TOL_CONE_ACCEPT` | 1e-6 | Largest scaled yield violation worth repairing |
| `TOL_CONE_MARGIN` | 1e-12 | Clearance over eigenvalue rounding, a few thousand `u` of `S_ref`, or of `D_scale` for Upper conditions. Unrelated to the floor |
| `TOL_INTERIOR_MIN` | 1e-6 | A strictly admissible field must clear yield by this, so `mu` stays small |
| `TOL_REPAIR_MAX` | 1e-3 | A repair may not move the value by more than about 0.1% |
| `TOL_FLOW_ACCEPT` | 1e-6 | Flow-rule and jump-cone analogue of `TOL_CONE_ACCEPT`, for inequality conditions only |
| `TOL_WORK_MIN` | 1e-12 | `W1` must be meaningfully positive against its own terms |
| `TOL_FEAS_OK` | 1e-9 | Least-squares residual indistinguishable from zero at this scale |
| `TOL_FEAS_LOCKED` | 1e-6 | Three decades above `TOL_FEAS_OK`, so solver accuracy cannot move a verdict across the band |
| `TOL_MONOTONE` | 1e-6 | Allowance for certified floors in the refinement checks of section 13 |
| `FLOOR_MIN_REL` | 1e-12 | Smallest floating-point floor ever reported |
| `TOL_CASE_ONE_SIDED` | 1e-6 | Allowance in one-sided case assertions, for floors only |
| `TOL_CASE_EXACT_GAP` | 1e-6 | Certified gap allowed where both exact fields lie in the discrete spaces |
| `MAX_GAP_STRIP` | 0.05 | Requirement: resolves `(2 + pi) cu` to about 2.5% either side on the listed finest mesh. Owner accepted (OD1) |
| `MAX_GAP_PRANDTL` | 0.10 | Requirement for the frictional strips, whose fans are wider and less resolved. Owner accepted (OD1) |
| `MAX_GAP_3D_VM` | 0.30 | Requirement for true-3D von Mises on mesh 5. A bracket this wide places the collapse load within 15% of its midpoint, the resolution at which a 3D bracket can still check a finite-element capacity. Owner accepted (OD1) |
| `MAX_GAP_3D_DP` | 0.50 | Requirement for true-3D Drucker-Prager on mesh 5 of each angle: within 25% of the midpoint, the coarsest resolution at which a frictional 3D bracket still says whether a capacity is plausible. Owner accepted (OD1); a 30-degree miss follows OD2 |

### 17.9 Status terms

| Term | Values and meaning |
|---|---|
| `solver_status` | `solved`, `almost_solved`, `primal_infeasible`, `dual_infeasible`, `time_limit`, `iteration_limit`, `numerical_error`, `not_run` |
| `classification` | `certified_bound`: passed this checker. `estimate`: everything else, with reasons |
| Inadmissible | Any estimate. It cannot satisfy a gate |
| Bracket `status` | `bracket`: `L <= U`. `order_within_floor`: `U < L <= U + floors`. `contradiction`: `L > U + floors`, which **stops the program** |
| Gap | `(U - L) / ((U + L) / 2)`, `null` unless `L > 0` and the status is `bracket` |
| Feasibility `status` | section 11.1 |

`primal_infeasible` for Lower or Upper, and `dual_infeasible`, are recorded as results. Write `[L, U]` only for a certified
bracket. Otherwise write `estimate: L=..., U=...`. Report `domain_claim`, criterion, `flow_rule` and interface law with
every result.

**Case lines.** `evaluate_case(case, certificates, results, feasibility, methods, *, brackets, refusals)` in
`tests/test_limit_contract.py` is the executable meaning of `assert`, `gate` and `record`. `methods` are the method
reports of section 15.9 and `refusals` the records of section 15.8 (amendment 1); `brackets` and `refusals` are keyword
arguments.

- Records are matched to meshes through `problem_hashes.json`, so a record of another case or mesh counts as missing.
- Every record used is validated against the problem built from the case: each certificate with its result document
  (`validate_certificate`), each feasibility record (`validate_feasibility`), and each bracket record supplied
  (`validate_bracket`). An invalid record is listed, excluded, and keeps the case from passing.
- Only gate-admissible records count; `schema-fixture` and `procedure-demo` never do. A gate certificate counts only with a
  valid method report bound to it (15.9); an invalid report is listed as an invalid record.
- A refusal (15.8) is validated against its result document. When a mesh and formulation have no gate-admissible
  certificate and a refusal from a gate checker, the item fails with cause `refused`, which is never withholdable.
- Bracket status and gap are **recomputed** from the two gate-admissible certificates with section 17.9's rules. A bracket
  record is checked, never trusted, and is not needed.
- Every `formulation` line binds: a Lower item holds for each Lower formulation, an Upper item for each Upper formulation,
  a pair item for each pair, and `control` lines bind nothing. Assertions hold on every listed mesh.
- For `gap_le`, status `order_within_floor` counts as gap 0 **only in an `exact_on_mesh` case**, where both exact fields lie
  in the discrete spaces; elsewhere it fails the item. `contradiction` fails the case and sets `stop_program`.
- Floors are not added to one-sided assertions: `TOL_CASE_ONE_SIDED` is the allowance for them.
- Each failed item carries a cause (`missing`, `estimate`, `refused`, `violated`, `order_within_floor`, `no_gap`,
  `contradiction`). The verdict sets `claim_withheld` when every failure is withholdable under OD2 (section 11.3), and
  `p0_blocking` otherwise. `evaluate_program(verdicts, required)` reports `stop_program`, `p0_open`, the blocking and
  missing cases, the cases with refused items, and the withheld claims.

## 18. Decision record

### 18.1 Accepted owner amendments

The owner accepted A1 to A4 for ADR 237 on 2026-11-18, and decided OD1 and OD2 on 2026-11-19. The `OD` ids are this
contract's; they are not the ADR's working defaults O1 and O2.

| Id | Amendment | In this contract |
|---|---|---|
| A1 | True-3D locking gate. Continuous P2 is a tested control, not an assumed element | Sections 8, 10, 11 |
| A2 | Whole-facet jump checks. P2 jumps are quadratic | Section 12 |
| A3 | `domain_claim` per bound. A bracket claims only what both bounds support | Section 14 |
| A4 | Certificate repair only by an admissible construction | Sections 17.3, 17.4 |
| OD1 | Keep the disclosed gap limits and gate meshes: `MAX_GAP_STRIP` 0.05, `MAX_GAP_PRANDTL` 0.10, `MAX_GAP_3D_VM` 0.30, `MAX_GAP_3D_DP` 0.50. They close from this review-fix revision, with the chronology of section 18.3. A miss is a real P0 failure | Sections 11.2, 16.4, 17.8 |
| OD2 | After the box correction, a 3D Drucker-Prager miss at 30 degrees withholds that 3D claim and is recorded by angle and wall type. Exact-formulation failures still stop P0 | Section 11.3 |

### 18.2 Orchestrator P0 rulings

These are orchestrator rulings for P0. R1 to R8 were recorded on 2026-11-18, R9 to R11 on 2026-11-19. **They are not
owner decisions.**

| Id | Ruling | In this contract |
|---|---|---|
| R1 | Versioned formulation ids, with cases that do not change meaning | Section 8 |
| R2 | `upper.p2d-jump-bernstein.v0` is the P0 Upper formulation; `upper.p2c-vertex.v0` stays a control | Sections 8, 12 |
| R3 | The true-3D set: `phi = 0, 10, 20, 30`, fixed and roller walls, non-structured meshes at three or more sizes, rank, residual and certified quality | Section 11; `upper_gate_3d` cases |
| R4 | Restrained footing components | Sections 6, 7 |
| R5 | Strip gate: maximum gap, fixed `x` and `z` walls, extents of at least `1.5 B` and depth `B / sqrt2`, and the box lemma. Here `z` is the depth axis and `1.5 B` is measured from the centreline | Section 16.4 |
| R6 | Zero extension may also serve a Lower bound when weightless with free walls; zero is an interior Lower field only without fixed loads | Sections 14, 17.3 |
| R7 | Certificate order: affine projection, inequality repair, full recheck; equality acceptance at rounding level | Section 17 |
| R8 | Gravity set: exact shear layer as the zero-gravity killer, exact weighted Drucker-Prager column, von Mises column as a two-sided sandwich only | Section 16.5 |
| R9 | The Upper work normalisation is not a projection row, unlike ADR 11.2 step 3: the value is `(D - W0) / W1`, so no exact normalisation is needed, and a row would only move the field | Section 17.4 |
| R10 | Equality acceptance is row-normwise with literal zeros, `tau_eq` carries no `kappa` factor, and `TOL_KAPPA_CAP` replaces the ADR's `TOL_EQ_CAP` | Sections 17.1, 17.2 |
| R11 | Per-angle Drucker-Prager boxes; a certified Upper at every size with the gap trend recorded; kinematic verdicts written by the checker from a witness | Section 11 |
| R12 | A checker capacity or conditioning failure is not an OD2 miss; it keeps P0 open. Recorded 2026-11-19, entered by amendment 1 | Sections 11.3, 15.8 |

### 18.3 Evidence notes and chronology

All values here are **estimates**, from runs of 2026-11-18 and 2026-11-19 with Clarabel 0.11.1, unless marked as a probe.
None is a certified bound: the demonstration below is not the independent checker.

**Chronology of the limits and gate meshes (ADR 18.2 item 4).** The ADR asks for each limit before the first solve of its
case. That did not happen, and the record says so:

| When | Event |
|---|---|
| `7909baa4e`, 2026-11-18 | First contract. No gap limit and no `gate` line. The strip meshes had depth along `y` |
| 2026-11-18, uncommitted | Estimates on the strips and the 3D gate, listed below with mesh indices of the present cases |
| `32c5a1842`, 2026-11-18 | The four `MAX_GAP_*` limits and every `gate` line first appear, **in the same commit** that records those estimates. The Drucker-Prager strips gained `mesh 64 1 32` as gate mesh 3, while their estimates are on mesh 2 |
| This revision, 2026-11-19 | Review-fix round. The Drucker-Prager 3D gate meshes are replaced by per-angle boxes (section 11) because the old box did not hold the mechanism scale. No solve has been run on any of the new meshes. The owner accepted the limits and gate meshes as disclosed (OD1); they are **closed** from this revision |

**Estimates and where they sit.**

- `Ladruno_files/testbed/limit_analysis/p2_isochory_probe.py` is a tracked **linear-algebra probe, not a bound**. Its output
  file is the record for the control's locking and for the jump space's mobility: the control locked with **fixed** walls
  (`rho_feas` 1.0e-1, 4.1e-2, 2.1e-2 on 243, 463, 949 tets) but mobile with **roller** walls (`rho_feas` near 1e-11). Roller
  cases therefore do not discriminate locking. The jump space is mobile for both wall types (`rho_feas` near 1e-11).
- Half-model strip without restrained components, mesh 0 (192 tets): `estimate: L=1.18, U=2.15`, against 5.14.
- Tresca strip, jump formulation, **mesh 3** (3,072 tets, the gate mesh): `estimate: L=5.005, U=5.186`, gap 0.036 against
  `MAX_GAP_STRIP` 0.05.
- Prandtl strips, `phi = 20 deg`, jump formulation, **mesh 2** (3,072 tets; the gate is mesh 3, 12,288 tets, unrun): cohesion
  `estimate: L=13.39, U=15.40` (gap 0.14) against 14.83; surcharge `estimate: L=5.873, U=6.604` (gap 0.12) against 6.399.
- True 3D, `phi = 0`, fixed walls, jump formulation: **mesh 4** (1,944 tets) `estimate: L=3.348, U=9.244`; **mesh 3** (2,840
  tets) `estimate: L=4.548, U=6.419` (gap 0.34). Mesh 5 (22,720 tets) is unrun.
- True 3D, `phi = 30 deg`, fixed walls, 1,944 tets of the **superseded** `[-2, 2]^2 x [-2, 0]` box: `estimate: L=9.65,
  U=168.9`. It says nothing about the new box.

**Certificate demonstration, 2026-11-19.** `Ladruno_files/testbed/limit_analysis/certificate_demo.py` assembles the programs of
section 15.4, solves them with Clarabel 0.11.1 at default settings, and runs sections 17.3 and 17.4 on the solver fields;
`certificate_demo_record.json` holds the documents. Producer and checker share one author, so the checker name is
`procedure-demo`, which no gate accepts, and the values are estimates for the program. It shows that the acceptance rule is
reachable on solver fields:

| Case, mesh | Bound | Solver status | Procedure | Value | `kappa_eq` | `omega_eq` | `repair_mu` |
|---|---|---|---|---|---|---|---|
| `tresca_strip_vm_thin_slab`, mesh 0 (192) | Lower | solved | passes | 4.684449 | 52.5 | 6.6e-17 | 0 |
| same | Upper, jump | solved | passes | 5.328414 | 142 | 6.7e-17 | 0 |
| `upper3d_fixed_phi00`, mesh 0 (243) | Lower | almost_solved | passes | 2.620562 | 190 | 4.6e-17 | 0 |
| same | Upper, jump | solved | passes | 12.310834 | 54.9 | 3.1e-17 | 0 |
| `dp_prandtl_cohesion_thin_slab`, unlisted Kuhn `8 1 4` (192) | Lower | almost_solved | passes | 10.078265 | 45.1 | 5.6e-17 | 0 |
| same | Upper, jump | solved | passes after repair | 16.751907 | 3.61 | 2.8e-17 | 1.0e-7 |

As estimates: strip `estimate: L=4.684, U=5.328` against 5.142; 3D `phi = 0` `estimate: L=2.621, U=12.311`; Drucker-Prager strip
`estimate: L=10.08, U=16.75` against 14.83. The Drucker-Prager repair used a relative-interior field (section 17.4 step 4);
its worst clearance on the repaired conditions was -3.9e-4 of its size. The review's componentwise `omega` had stayed at 1
on the same kind of fields.

**Rerun under amendment 1, 2026-11-19** (same driver, Clarabel 0.11.1, numpy 2.3.5, scipy 1.16.3). The demonstration now
writes Upper `objective = c^T x`, builds Lower symmetry rows from the set normal, and scales flow and jump conditions by
`D_scale_vol` and `D_scale_jump` separately. Every value, `kappa_eq`, `omega_eq` and `repair_mu` in the table is unchanged to
the digits shown, and the Drucker-Prager clearance is still -3.9e-4. On the strip the jump Upper has `objective` 5.3284148
against the field `value` 5.3284138, the gap that item 2 of amendment 1 makes explicit. Its `kappa_eq` is the ratio of
computed dense singular values, **not** a proved bound of section 17.1, and it writes no method report: the values remain
estimates and the reserved checker name keeps them from every gate.

## 19. Amendment log

The file is not frozen, so the 2026-11-19 review-fix round is recorded here rather than as amendments.

| Finding | Resolution |
|---|---|
| B1 | Row-normwise `omega` with scales; `tau_eq` without `kappa`; literal zeros; certification demonstrated on strip mesh 0 and `upper3d_a` (sections 17.1, 17.2, 18.3) |
| B2 | Signed conditions, the set `P`, the interior-field requirement and `mu` stated; rounding margin `TOL_CONE_MARGIN` on the recheck; forced-zero conditions explained; repair demonstrated on a Drucker-Prager strip (17.4) |
| B3 | `domain.ground`, `domain.symmetry_sets` and the coverage rule (14); the half-space fixture rebuilt with every side fixed and an admissible compact field; negative fixture `unlisted_side_face` |
| M1 | File names carry case, mesh and formulation; brackets name both formulations and the validator checks them (15) |
| M2 | Witness result hash in the feasibility record, written by the checker; nullity out of the verdict; `alpha > 0` from the witness alone and never `formulation_infeasible`; `gate_admissible` covers feasibility records (11.1, 15.5) |
| M3 | A `certified_bound` must show its numbers within limits, with malformed-certificate tests (17.7) |
| M4 | The six missing ADR fixtures added; the jump fixture rebuilt on a strictly dilating field; one expected code per fixture; coverage tested against the ADR rows (17.6) |
| M5 | `TOL_KAPPA_CAP` = 1e8 stated as a number; method named; 237c reports `kappa_eq` before any gate (17.1, 17.2) |
| M6 | `certified_upper` on meshes 0 to 5 and `record gap` per mesh in all eight 3D cases (11.2) |
| M7 | `evaluate_case` with synthetic-record tests (17.9) |
| M8 | Chronology in 18.3; owner decision OD1 closes the limits and gate meshes |
| M9 | Per-angle boxes and meshes, verdicts per angle and wall type, owner decision OD2 (11) |
| m1 | `lower_le_formula` and `upper_ge_formula` added and asserted on the column (16.5) |
| m2 | Candidate `lambda` is `value`; `resultant` recomputed; scaled unknowns; fixture `wrong_lower_value` (17.3) |
| m3 | Premise (a) at 20 degrees marked unverified in 16.4 and in both case files; Prandtl 1921 cited beside 1920 |
| m4 | The bracket rule's first branch cannot fire in version 0; said in 14 |
| m5 | Worked digest, `-0.0` canonicalised, non-finite refused after parsing, Kuhn arithmetic stated, `problem_hashes.json` (15.1) |
| m6 | Tautology replaced by `k / root == c`; docstring corrected; the probe test renamed to say it reads a record |
| m7 | R9 records the query normalisation; 17.5 states the norm-equivalence factor and the margin's role |
| m8 | Mesh digests in 15.2; the public-file test opens every NPZ |

**Pre-freeze round, 2026-11-19.**

| Item | Resolution |
|---|---|
| N1 | OD2's scope stated in 11.3: only `gap_le MAX_GAP_3D_DP` and missing certificates at 30 degrees withhold; executable in `_withholdable`, `evaluate_case` and `evaluate_program`, with tests |
| M7 residuals | The evaluator recomputes status and gap from certificates, validates every certificate, result-linked certificate, feasibility and bracket record, and limits the zero-gap rule to `exact_on_mesh`; tests include a malformed feasibility record and a forged bracket gap |
| `D_scale_int` | Must be positive (17.4 step 4) |
| Shift-invert | The shift must be positive on the singular normal matrix; skipped eigenvalues must number `n_rows - rank` (17.1) |
| Literal zeros | `literal_zero_violated` when a literal-zero unknown is not exactly zero in the final field (17.3, 17.4, 17.6) |
| Decision ids | O1 and O2 renamed OD1 and OD2 so they cannot be read as the ADR's working defaults |
| Demonstration | Its invalid-interior code is `interior_field_invalid`, it returns `no_interior_field` when zero stress is not interior, it requires `D_scale_int > 0`, and it checks literal zeros on the final field |

**ADR 237 changes this contract required in the review-fix round** (the ADR is not edited here). Items 1 to 10 were applied
on `origin/wp/237-limit-analysis` at `1375c2c46`. Still owed after the pre-freeze round:

- 19.2: the contract's decision labels are now OD1 and OD2, not O1 and O2.
- 10.1 verdicts: only `gap_le MAX_GAP_3D_DP` and a missing certificate withhold the 30-degree claim; monotonicity failures,
  kinematic failures and invalid records at any angle keep P0 open, and a contradiction stops it.
- 11.3 or 11.4: `literal_zero_violated`, the positive shift for shift-invert, and `D_scale_int > 0` for the interior field.

Review-fix list, applied:

1. 11.3: replace the componentwise scaled residual with the row-normwise `omega` and scales of section 17.1, drop the `kappa`
   factor from `tau_eq`, and replace `TOL_EQ_CAP` by `TOL_KAPPA_CAP` on `kappa_eq` (R10). Add literal zeros.
2. 11.2 step 3: the Upper projection does not contain the query normalisation (R9).
3. 11.2 step 4 and 10.4: state the signed Upper conditions, the repaired set `P`, the relative-interior requirement and
   the rounding margin of section 17.4.
4. 10.1 gate table: per-angle boxes sized by the plane-strain extent; a certified Upper at every size with the gap trend
   recorded; verdicts per angle and wall type; the kinematic verdict from a checker-recomputed witness, with nullity
   reported only and no infeasibility verdict for `alpha > 0`.
5. 10.1 and 18.3: record OD2. A 30-degree Drucker-Prager miss withholds that 3D claim instead of triggering K2; exact
   formulation failures still stop P0.
6. 18.2 item 4: record OD1 and the chronology of section 18.3, since the limits were not set before the first solve.
7. 12: the half-space coverage rule with `ground` and `symmetry_sets`; the bracket's half-space branch cannot fire in
   version 0.
8. 11.5: no table change; add that each row is discharged by a named fixture with one expected code.
9. 16.1: the column assertions `L <= S - gamma H/2` and `U >= S - gamma H` are now in the case.
10. 5 and 6 (file formats): the record names of section 15, the bracket's formulation keys, and the reserved checker name
    `procedure-demo`.

### 19.1 Amendment 1, 2026-11-19

**Freeze record.** The orchestrator froze this file, `tests/limit_cases.txt` and `tests/test_limit_contract.py` at
`03f59c1a4` (WP-237a pre-freeze round). Everything above this subsection, except the passages marked "amendment 1", is
that frozen text. Amendment 1 is the first change after the freeze. The orchestrator requested it, and it was written
without reading the reference or checker code.

| Item | Change | Sections | Executable form |
|---|---|---|---|
| 1 | Status: frozen at `03f59c1a4`, amendment 1 recorded after it | header, 19.1 | `test_interface_records_the_freeze_and_amendment_1` |
| 2 | Upper `value` is the field value of 17.4 step 6; `objective` is `c^T x` for both bounds; Lower `value = objective_sign * objective` | 15.3, 15.4 | fixture Upper results carry `objective = D - W0` |
| 3 | Symmetry rows use the declared set normal `n_S` in both bounds, in the program and in both projections | 5, 15.4, 17.3, 17.4 | text; the demonstration's Lower rows now use `n_S` |
| 4 | Problem `name` of a case problem: `<case>:<nx>x<ny>x<nz>` or `<case>:<stem>+<r>`, with the other derived values | 15.2 | `build_case_problem`, `test_case_problem_name_rule` |
| 5 | Method report `<stem>.<query>.<formulation-id>.certificate.method.json`, bound by `certificate_sha256`, recording the kappa method, dependent-row proofs, rows skipped, Ritz or eigenvalue residual and interior-field identity; required for every gate certificate | 15.1, 15.6, 15.9, 17.9 | `validate_method_report`; `gate_admissible(..., method=)`, `validate_bracket(..., lower_method=, upper_method=)`, `evaluate_case(..., methods)` |
| 6 | Conditioning: dependent rows only by shape-forced zeros, dense singular values on disjoint row blocks, or an explicit null-space basis with a residual bound; `kappa_eq` a proved upper bound; `A_hat^T A_hat` allowed; no proof gives a refusal. Replaces the frozen rule, which accepted an estimate within a factor 2 and skipped eigenvalues below a threshold without proof | 17.1, 17.2 | method-report rules |
| 7 | Capacity-refusal record; cause `refused`; not a missing certificate under OD2; keeps P0 open; ruling R12 | 11.3, 15.8, 17.9, 18.2 | `validate_refusal`, `evaluate_case(..., refusals=)`, `evaluate_program` |
| 8 | Mesh-5 conditioning route: present capacity insufficient for mesh-5 Lower (all 3D cases) and von Mises jump Upper; WP-237d must identify an open, proved route; gates unchanged | 11.2 | none; gates unchanged |
| 9 | Fixtures: `normal_jump_one_node` replaced by `normal_jump_isochoric`; `vm_jump_opening` and `restrained_footing_component` added, with the fixture problem `vm_plate_shear_restrained` | 17.6 | `negative_fixtures.json`, `apply_negative` |
| 10 | Feasibility witness field is the certified projected field; `footing_modes` is `null` without a footing; `D_scale` split into `D_scale_vol` (1/time) and `D_scale_jump` (velocity) | 11.1, 15.5, 17.4, 17.5 | `validate_feasibility`; text |

**What changes downstream.** A producer writes Upper `value` from the field and `objective = c^T x`. The checker writes a
method report beside every certificate, a refusal instead of a certificate when no bound is proved, rows with `n_S`, the
witness residual at the projected field, and two flow scales. A caller of `evaluate_case` passes method reports as the
fifth positional argument and brackets and refusals by keyword. Every fixture certificate and result was regenerated; no
case value, limit, mesh, tolerance or gate line changed, and every `problem_sha256` is unchanged.
