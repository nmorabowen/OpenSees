# Limit analysis: the frozen interface for P0 (WP-237a)

**Status:** proposed for freeze. The orchestrator freezes this file. After the freeze a change is an
amendment logged in section 19, never a silent edit. Items in section 18 are **proposed, pending the
owner's decision**; they are written here so that P0 can test them, and they do not record a decision.

**Who reads it.** The Python reference author, the independent checker author, the external reference
author and later the C++ author. Each implements from this file and `tests/limit_cases.txt`, without
reading the other side's code. `tests/test_limit_contract.py` holds the schema validators and the
convention checks. It contains no assembler and no checker.

**Scope of version 0.** Associated perfect plasticity on straight tetrahedra; von Mises and
Drucker-Prager; second-order cones only; bonded rigid footings; one-shot solves on one mesh. The file
format does not presuppose the C++ engine or the final `.la.h5` layout.

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

| Kind | Lower, at each facet vertex | Upper, at all six facet nodes |
|---|---|---|
| `free` | `N(n) sigma = 0` | none |
| `traction` | `N(n) sigma = t0 + lambda t1` | none; adds work |
| `fixed` | none | `v = 0` |
| `symmetry` | `T1^T N(n) sigma = T2^T N(n) sigma = 0` | `n^T v = 0` |
| `footing` | enters the resultant (section 6) | rigid map (section 6) |
| `interface` | reserved, refused in version 0 | reserved |

**Tangent basis.** `m` is the index of the smallest `|n_m|`, the first if tied. `T1 = (e_m x n)/|e_m x n|`
and `T2 = n x T1`.

**Symmetry sets** must be planar. All their facets share one normal and one plane within
`TOL_GEOM_PLANE`. The set's normal is the normal of its first facet in `facets` order.

**Corners.** Conditions apply per facet in the Lower problem and accumulate per node in the Upper problem.

## 6. Rigid footings

A footing has a reference point `x_ref`, one or more `footing` boundary sets, and an `interface_law`.
Version 0 admits `bonded`. `no_tension` and `coulomb` are reserved and refused (section 18, P-3).

- **Lower.** No extra traction variable. With facet vertices `a` and area `A`,

  ```text
  F = sum_facets sum_a (A/3) t_a
  M = sum_facets sum_a sum_b A w_ab (x_a - x_ref) x t_b,   w_ab = 1/6 if a = b, else 1/12
  ```

  with `t_a = N(n) sigma_(e,a)`. Both are exact for linear traction on a flat facet.
- **Upper.** Every node of a footing facet follows `v_a = u + w x (x_a - x_ref)` exactly. P2 contains the
  linear rigid field, so the whole facet moves rigidly.
- **Power.** `int t . v dA = F . u + M . w = R^T V`.
- **Restrained components** (section 18, P-6). A query may list components as restrained. Then
  `V_c = 0`, and `R_c` is an unconstrained reaction. A symmetric half model needs this: the half footing
  of a strip restrains `ux` and `rz`.

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
| `upper.p2c-vertex.v0` | Upper | admitted **as the baseline under test**; section 10. It is not proven for 3D isochoric flow |
| `upper.p2d-jump-bernstein.v0` | Upper | **proposed, refused** by version 0 validators; section 12 |

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

## 10. Upper problem, `upper.p2c-vertex.v0`

Unknowns: continuous P2 `v_a`; `V_f`; `(r_ei, z_ei)`.

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

## 11. True-3D Upper feasibility gate (proposed)

Thin-slab success does not show that an Upper formulation works in 3D. Each `upper_feasibility_3d` case
runs on the committed unstructured meshes for `phi = 0` and for the planned Drucker-Prager angles. Those
angles are provisional (section 18, P-1).

Let `A` be the set of `(v, V)` that satisfies every linear boundary and footing condition, and the query
normalisation. At vertex `(e, i)` let

```text
c_ei(v) = |theta|                                alpha_e = 0
c_ei(v) = max(0, 3 sqrt2 alpha_e |e| - theta)    alpha_e > 0
```

With `w_ei = V_e / 4`, `|c|_w^2 = sum w_ei c_ei^2` and `|D|_w^2 = sum w_ei |D_ei|_F^2`,

```text
rho_feas = min_A |c(v)|_w / min_A |D(v)|_w .
```

It is dimensionless. With every `alpha = 0`, both minima are sparse linear least-squares problems, so the
result does not depend on a conic solver.

| `status` | Condition |
|---|---|
| `formulation_feasible` | `rho_feas <= TOL_FEAS_OK`, **and** a witness `(v, V)` in `A` passed `c_ei` recomputed independently of the producer |
| `formulation_infeasible` | `rho_feas >= TOL_FEAS_LOCKED` |
| `inconclusive` | anything else |

- A failure is a **formulation result**. It is recorded with its `rho_feas`, per mesh and per angle.
- It may not be hidden by relaxing a tolerance, reading `AlmostSolved` as success, or quoting a value from
  a conic run whose flow rule is violated.
- `formulation_infeasible` says no Upper bound comes from this formulation on that mesh. It does not say
  the mechanics has no collapse load.
- Under P0, a failure means selecting and testing another Upper formulation before C++.

## 12. Velocity jumps (proposed, not admitted in version 0)

P2 velocity jumps are **quadratic** on a facet, so conditions at the three facet vertices do not cover the
facet. This section defines the whole-facet condition for `upper.p2d-jump-bernstein.v0`.

Let `j = v+ - v-`, with `n` the outward normal of `e-`, so `j` is the velocity of the `e+` side relative
to the `e-` side along `n`. Its six nodal values are `j_a` at vertices and `j_ab` at edge midpoints.

**Bernstein coefficients.** `c_a = j_a` and `c_ab = 2 j_ab - (j_a + j_b)/2`. Then
`j(x) = sum_B B_B(x) c_B`, with `B_B >= 0`, `sum_B B_B = 1` and `int B_B dA = A/6` for each of the six.

**Jump dissipation.** `t . j = sigma : sym(j (x) n)`, so the jump uses the same support function with
`d = N(n)^T j`:

```text
N(n)^T c_B = sqrt2 alpha r_B one + Q^T z_B,   (r_B, z_B) in SOC6,   each of the 6 coefficients
D_facet    = sum_B (A/6) sqrt2 k r_B
```

**Proof.** The admissible set is a convex cone and `j(x)` is a convex combination of the `c_B`, so the
flow condition holds at every point. `pi` is convex and positively homogeneous, so
`int pi(N^T j) dA <= sum_B (A/6) pi(N^T c_B)`. The result is conservative, as an Upper bound must be.

For a slip `s` and opening `o` this gives `o = 3 alpha s / sqrt(1 - 12 alpha^2)` on the cone boundary. With
the plane-strain match of section 4.4 that is `o = s tan(phi)`, the associated dilation of a Coulomb surface.

On a facet between two regions, use the strength of `e-`. Every admissible stress field meets the yield
condition on both sides, so either side is valid.

## 13. Conditions for each bound claim

**Definitions.** `sigma` is statically admissible if it satisfies equilibrium, traction continuity, the
traction and footing conditions, and yield everywhere. `(v, V)` is kinematically admissible if it is
continuous, or has jumps charged by section 12; it satisfies the velocity conditions, and its strain rate
lies in `dom pi` almost everywhere.

**Lower claim, `lambda_L <= lambda_max`, holds when:** elements are straight; `k`, `alpha`, `b` are
constant per element; tractions are linear per facet; resultants are integrated exactly; equalities hold;
yield holds at four vertices. Then equilibrium and continuity hold everywhere because all residuals are
polynomials determined by the rows imposed, and yield holds everywhere by convexity.

**Upper claim, `lambda_max <= lambda_U`, holds when:** elements are straight; `k`, `alpha`, `b` are
constant per element; velocity conditions hold exactly; `W1 > 0`; every vertex strain rate lies in
`dom pi`. Virtual work gives `W0 + lambda W1 = int sigma : D dV` for any admissible pair, so
`lambda W1 <= int pi(D) dV - W0`. `D(x)` is linear, so the cone holds throughout, and Jensen's inequality
with `int L_i dV = V/4` gives `int pi(D) dV <= sum_i (V/4) pi(d_i)`.

**Both are claims about the stated idealisation:** associated flow, perfect plasticity, the stated
strength and interface law. For non-associated real soil they are not bounds on the physical problem.

## 14. Domain claims

`domain_claim` belongs to **each certified bound**. A problem declares `domain.target`, `box` or
`halfspace`. With `halfspace` it also lists `truncation_sets`, and `exterior_load_free`, which is true when
no load acts outside the mesh.

| Value | Meaning |
|---|---|
| `box_stated` | Bound for the stated box with its exact wall conditions. The only value for target `box` |
| `unverified` | Target `halfspace`, extension not established. The bound holds for the box and is an estimate for the half-space |
| `halfspace_by_zero_extension` | The admissible field extends by zero |
| `halfspace_by_extension_elements` | Reserved. Version 0 refuses it because no extension data exists |

`halfspace_by_zero_extension` is supported only as follows.

- **Upper:** every truncation set is `fixed`. Zero velocity costs no dissipation and no work outside.
- **Lower:** every truncation set is `free`, every `b0` and `b1` is zero, and `exterior_load_free` is true.
  Then zero stress is admissible outside. With self-weight this route is closed.

Roller walls allow tangential motion and never support the Upper claim. A fixed wall supports no Lower
extension. **A bracket's claim is the domain both bounds support:** if both certificates carry
`halfspace_by_zero_extension`, the bracket does. Otherwise it is `box_stated`, and it passes no half-space
gate.

## 15. Files, version 0

Every document is UTF-8 JSON, with `NaN` and `Infinity` forbidden. Unknown keys are refused. Arrays live in
one NPZ named by the document's `arrays` key, relative to the JSON file. Arrays are `<f8` or `<i8` only.

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
`problem_sha256` in a problem and `result_sha256` in a result.

### 15.2 Problem: `<stem>.problem.json`, `.npz`

| Key | Content |
|---|---|
| `format`, `format_version` | `ladruno-limit-problem`, `0` |
| `arrays`, `name`, `units` | NPZ file name; free name; `{length, force}` labels |
| `conventions` | the literal object of section 15.7 |
| `domain` | `{target, truncation_sets, exterior_load_free}` |
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
A `mesh_file` is an NPZ with `nodes (n, 3)` of `<f8` and `tets (E, 4)` of `<i8`, positively oriented, and
section 3.5 adds the mid-edge nodes. The three committed files are unstructured tetrahedral meshes of the box
`[-2, 2] x [-2, 2] x [-2, 0]`, with the square `[-0.5, 0.5]^2` conforming on `z = 0`, made with gmsh 4.15
(Delaunay 3D, one thread, seed 1, target sizes 1.0, 0.8, 0.58). Their 243, 463 and 949 tetrahedra are data.
No test or implementation depends on gmsh.

### 15.3 Result: `<stem>.<query>.<bound>.result.json`, `.npz`

Keys: `format` = `ladruno-limit-result`, `format_version`, `arrays`, `problem_sha256`, `query_id`, `bound`,
`formulation`, `producer` `{name, version}`, `solver` `{name, version, settings}`, `solver_status`,
`objective`, `value`, `wall_time_s`, `result_sha256`.

| `bound` | Required arrays | Optional, never trusted |
|---|---|---|
| `lower` | `sigma (E,4,6)`, `resultant (n_f,6)` | `interior_sigma (E,4,6)` with number `interior_value`; `equilibrium_multiplier (E,3)` |
| `upper` | `velocity (N,3)`, `footing_velocity (n_f,6)` | `interior_velocity (N,3)` and `interior_footing_velocity (n_f,6)`; `strain_rate (E,4,6)`; `dissipation (E,4)` |

`interior_value` is a top-level key, present exactly when `interior_sigma` is. With `solver_status` of
`solved` or `almost_solved`, `value` and the arrays are required. Otherwise `value` and `objective` are
`null` and no field arrays are stored: a time limit gives **no number**. A `lower` result carrying velocity
arrays, or the reverse, is refused.

### 15.4 Program: `<stem>.<query>.<bound>.program.json`, `.npz`

Canonical form: minimise `c^T x` subject to `A x + s = b`, `s` in `K`, with `K` = zero cone, then `n_soc`
copies of `SOC6`. `value = objective_sign * c^T x`.

Keys: `format` = `ladruno-limit-program`, `format_version`, `arrays`, `problem_sha256`, `query_id`, `bound`,
`formulation`, `objective_sign`, `n_variables`, `n_rows`, `n_zero`, `n_soc`, `soc_dim` (= 6). Arrays: CSC
`A_indptr`, `A_indices`, `A_data`; `b`; `c`. There is **no scaling and no presolve** in the program. A solver
adapter may scale or presolve its own copy. Rows may be linearly dependent.

**Lower.** `x = [lambda (ray forms only), sigma_(0,0)...]`, stress in `(e, i, component)` order.
`objective_sign = -1`. `c` is `-1` on `lambda`, or `-d^T R` as a row over `sigma` for `support`.
Rows, in order:

1. equilibrium, 3 per element: `[b1_e | sum N(grad L_i)]`, right side `-b0_e`;
2. interior facets in section 3.4 order, 3 vertices, 3 components: `N(n) sigma_- - N(n) sigma_+ = 0`;
3. exterior facets in `facets` order, 3 vertices: `free` and `traction`, 3 rows, `N(n) sigma - lambda t1 = t0`;
   `symmetry`, rows `T1`, `T2`; none for others;
4. ray forms: for each footing and unrestrained component, `R - lambda F1 = F0`;
5. cones, per `(e, i)`: rows `[sqrt2 alpha one^T ; -Q]` on `sigma_ei`, right side `(sqrt2 k, 0)`.

**Upper.** `x = [v (3N), V (6 n_f), (r, z)_(0,0)...]`. `objective_sign = +1`.
`c = (V_e/4) sqrt2 k_e` on each `r`, minus the `W0` coefficients. Rows, in order:

1. flow, per `(e, i)`, 6 rows: `B_ei v - sqrt2 alpha r one - Q^T z = 0`;
2. `fixed`: for each fixed node ascending, 3 rows `v = 0`;
3. `symmetry`: per set in order, nodes ascending, fixed nodes omitted, `n^T v = 0`;
4. `footing`: per footing, nodes ascending, 3 rows `v - u - w x (x_a - x_ref) = 0`;
5. restrained components, in list order: `V = 0`;
6. ray forms: `W1 = 1`; `support`: `V_fc = d_fc` over unrestrained `(f, c)`;
7. cones, per `(e, i)`: `-(r, z) + s = 0`.

### 15.5 Feasibility record: `<stem>.<query>.feasibility.json`

Keys: `format` = `ladruno-limit-feasibility`, `format_version`, `problem_sha256`, `query_id`, `formulation`,
`status`, `rho_feas`, `numerator`, `denominator`, `method`, `witness_checked`, `producer`, `tolerances`.
`rho_feas = numerator / denominator` and `status` follows section 11.

### 15.6 Certificate and bracket

`<stem>.<query>.<bound>.certificate.json`, written by the checker only:

| Key | Content |
|---|---|
| `format`, `format_version` | `ladruno-limit-certificate`, `0` |
| `problem_sha256`, `result_sha256`, `query_id`, `bound`, `formulation` | identity |
| `classification` | `certified_bound` or `estimate` |
| `reasons` | codes of section 17.3; empty exactly when certified |
| `value_candidate`, `value_certified`, `floor` | numbers; the last two are `null` for an estimate |
| `domain_claim`, `criteria`, `flow_rule`, `interface_laws` | section 14; sorted lists; the literal `associated` |
| `checks` | metrics of section 17.4 |
| `checker` | `{name, version}` |
| `tolerances` | the table of section 17.5 |

`<stem>.<query>.bracket.json`: `format` = `ladruno-limit-bracket`, `format_version`, `problem_sha256`,
`query_id`, `lower_result_sha256`, `upper_result_sha256`, `lower`, `upper`, `floor_lower`, `floor_upper`,
`status`, `gap`, `domain_claim`, `lower_domain_claim`, `upper_domain_claim`.

The checker name `schema-fixture` is reserved for hand-made fixtures. No gate may accept it.

### 15.7 Literal conventions object

```json
{"sign": "tension_positive", "voigt": ["xx", "yy", "zz", "xy", "yz", "zx"],
 "stress_shear": "tensor", "strain_rate_shear": "engineering", "tet10": "vtk",
 "facets": "opposite_vertex_outward", "footing_resultant": "footing_on_soil"}
```

## 16. Derivations behind the cases

### 16.1 Uniform cube

Uniform stress lies in the Lower space and linear velocity in the Upper space, so both are exact on every
mesh. With `sigma = s e_x e_x`: `sqrt(J2) = |s|/sqrt3`, `I1 = s`, so tension fails at
`k / (1/sqrt3 + alpha)` and compression at `k / (1/sqrt3 - alpha)`. At the apex, `sigma = p I` and
`p = k / (3 alpha)`, with Upper field `v = x`. For von Mises, `v = (x, -y/2, -z/2)` gives
`pi = sqrt3 k`, and `W1 = 1`. A dead traction `f sqrt3 k` leaves `lambda = (1 - f) sqrt3 k`.

### 16.2 Rigid plate in shear

`sigma_xy = lambda` only, with `t1 = (0, +-1, 0)` on the `x` faces. Upper field `v = (y, 0, 0)` with plate
velocity `(1, 0, 0)`: the side traction does no work, `pi = k`, `lambda = k`. With `x_ref` raised by `h`
above the plate centre, `M = (0, 0, h) lambda A`, which is the offset case.

### 16.3 Body-force shear block

`sigma_xy = lambda gamma (H - y)` is linear, so the Lower field is exact on every mesh. Yield at the base
gives `lambda = k / (gamma H)`. The Upper field `v_x = min(y/h, 1)`, with `h` the first layer, gives
`lambda <= k / (gamma (H - h/2))`, the `shear_band` formula.

### 16.4 Strip on a thin slab

**Slab equals plane strain.** With `symmetry` on both `z` faces, a plane-strain velocity is admissible in
the slab, and a plane-strain stress field extended with the optimal `szz` of section 4.4 is admissible too.
So the slab and the plane-strain collapse loads coincide.

**The stated box.** The reference is valid for the box because: (a) the walls are **fixed**; (b) the Prandtl
mechanism has compact support and fits inside. Then the mechanism gives `lambda_box <= target`. A
half-space admissible stress field restricted to the box is admissible for it, so `lambda_box >= target`.
With roller walls, (b) still gives the first inequality, but the second is not automatic. **Only the Lower
assertion relies on the mechanism, and only the Upper assertion relies on the stress extension.**

**Mechanism extent,** footing half-width `b`, `beta = pi/4 + phi/2`: `r0 = b / cos(beta)`,
`r1 = r0 exp((pi/2) tan(phi))`. It reaches `x = b + 2 r1 cos(pi/4 - phi/2)` and depth at most `r1`. For
`phi = 0`: `x = 3b`, depth `sqrt2 b`. For `phi = 20 deg`, `b = 1`: `x = 6.06`, depth at most `3.09`.

**Values.** The slip-line invariant along the fan gives `Nq = exp(pi tan(phi)) tan^2(beta)`, and
`Nc = (Nq - 1) cot(phi)`, tending to `2 + pi`. A bonded footing gives the same value as a smooth one here:
the wedge moves with the footing, and a smooth-footing stress field is admissible for a bonded footing.

## 17. Certificate

The checker works from the problem and result only, rebuilds every operator, and never reads solver cones.

### 17.1 Lower procedure

1. Check hashes, query, formulation, straightness and volumes.
2. Compute scaled residuals (section 17.4). If `rho_eq > TOL_EQ_ACCEPT`, the result is an estimate.
3. Correct equalities with a minimum-norm `delta sigma`. Require `rho_eq_final <= TOL_EQ_FINAL`.
4. Evaluate `g_ei = sqrt(J2) + alpha I1 - k` from **principal stresses**. If
   `max g / S_ref > TOL_CONE_ACCEPT`, the result is an estimate.
5. If any `g_ei > -TOL_CONE_MARGIN S_ref`, **repair with a strictly admissible field**
   `(sigma_int, value_int)`, which must satisfy the same equalities and `g_int <= -TOL_INTERIOR_MIN S_ref`. Use
   `mu = max (g + m) / (g - g_int)`, `m = TOL_CONE_MARGIN S_ref`, over the offending vertices, and
   `sigma_new = (1 - mu) sigma + mu sigma_int`. Convexity gives admissibility.
   `value_new = (1 - mu) value + mu value_int`.
6. `sigma_int = 0`, `value_int = 0` is available only when `b0`, `t0` and `F0` are all zero and every `k > 0`.
   Otherwise the field must be supplied as `interior_sigma` and is checked like the candidate. If none
   exists, the result stays an estimate with `no_interior_field`. Require `mu <= TOL_REPAIR_MAX`.

### 17.2 Upper procedure

1. Check identity and geometry as above.
2. Overwrite velocity conditions exactly: fixed nodes to zero, footing nodes from `V`, symmetry nodes by
   removing the normal part. Before overwriting, the change must be within `TOL_EQ_ACCEPT` of the largest
   velocity.
3. Compute `d_ei` by the checker's own operator, and `rho_flow = max c_ei / max |D_ei|_F`. If it exceeds
   `TOL_FLOW_ACCEPT`, the result is an estimate.
4. **Correct to an admissible field.** A finite dissipation charge is not a repair: `pi` is `+infinity`
   outside its domain.
   - `alpha = 0` vertices: the domain is a subspace with empty interior. Apply a minimum-norm correction
     `delta v`, homogeneous in the velocity conditions, with `theta(v + delta v) = 0`. Require
     `rho_iso_final <= TOL_FLOW_FINAL`.
   - `alpha > 0` vertices: blend with a supplied strictly admissible field,
     `v_new = (1 - mu) v + mu v_int`, where `v_int` satisfies the same linear conditions, and the same footing
     velocity for `support`, and has `c < 0` at these vertices and `theta` at the floor at `alpha = 0`
     vertices. `c_ei` is convex and positively homogeneous, so `c(v_new) <= 0` for
     `mu = max c / (c - c_int)`.
   - If no construction exists, the result stays an estimate with `no_admissible_correction`.
5. Require `W1 >= TOL_WORK_MIN * W_scale` for ray forms, with `W_scale` the sum of absolute work terms.
6. `D = sum (V_e/4) pi(d_ei)`, with `pi` from section 4.3. Ray forms: `U = (D - W0) / W1`, so no exact
   normalisation is needed. `support`: `U = D - W0`.

### 17.3 Outcome

Compare the recomputed value before repair with the result's `value`. A difference above
`TOL_EQ_ACCEPT * max(|value|, |recomputed|)` gives `value_mismatch`. That catches a wrong normalisation or
swapped labels.

`floor = max(FLOOR_MIN_REL, final residual) * |value_certified|`. It is a floating-point statement, not a
theorem. "Rigorous" means rigorous for the stated idealisation, up to this floor.

Reason codes: `no_candidate`, `hash_mismatch`, `query_mismatch`, `formulation_not_admitted`,
`curved_element`, `degenerate_element`, `boundary_violation`, `equality_residual_too_large`,
`equality_correction_failed`, `cone_violation_too_large`, `no_interior_field`, `interior_field_invalid`,
`repair_too_large`, `flow_violation_too_large`, `no_admissible_correction`, `work_nonpositive`,
`value_mismatch`, `unsupported_domain_claim`.

### 17.4 Scaled metrics

`L_ref` is the bounding-box diagonal.
`S_ref = max(max k, max |sigma|, max |t0|, |value| max |t1|, L_ref max(|b0| + |value| |b1|), max_f |F_f| / A_f)`,
using the components of `sigma` in the Lower case and omitting that term in the Upper case. For `support`,
the `|value|` terms are absent.

`checks` holds these keys, `null` where not applicable: `rho_div` (`|div sigma + b|_inf h_e / S_ref`),
`rho_jump`, `rho_bc`, `rho_resultant` (force over `S_ref A_f`, moment over `S_ref A_f ell`), `rho_eq_final`,
`rho_cone`, `repair_mu`, `straightness_max`, `rho_flow`, `rho_iso_final`, `work_reference`.

### 17.5 Named tolerances

Values are provisional until 237c and 237d calibrate them on the exact cases. A value changes only by
amendment, with evidence, and **never to rescue a run**.

| Name | Value | Rationale |
|---|---|---|
| `TOL_GEOM_STRAIGHT` | 1e-9 | Relative to `h_e`. Admits meshes written with 12 significant digits at `L/h` up to 1e3 |
| `TOL_GEOM_VOLUME` | 1e-10 | `V_e >= tol h_e^3`. A regular tetrahedron has 0.118; below this, gradients amplify rounding by 1e10 |
| `TOL_GEOM_PLANE` | 1e-9 | Normal deviation and offset, relative to `L_ref` |
| `TOL_EQ_ACCEPT` | 1e-6 | Largest scaled equality residual worth correcting; about 10 to 100 times interior-point feasibility tolerance |
| `TOL_EQ_FINAL` | 1e-11 | After correction; about 5e4 unit roundoffs, allowing for operator conditioning |
| `TOL_CONE_ACCEPT` | 1e-6 | Same scale as `TOL_EQ_ACCEPT` |
| `TOL_CONE_MARGIN` | 1e-12 | Clearance over eigenvalue rounding, a few thousand unit roundoffs of `S_ref` |
| `TOL_INTERIOR_MIN` | 1e-6 | A strictly admissible field must clear yield by this, so `mu` stays small |
| `TOL_REPAIR_MAX` | 1e-3 | A repair may not move the value by more than about 0.1% |
| `TOL_FLOW_ACCEPT` | 1e-6 | Flow-rule analogue of `TOL_CONE_ACCEPT` |
| `TOL_FLOW_FINAL` | 1e-11 | Flow-rule analogue of `TOL_EQ_FINAL` |
| `TOL_WORK_MIN` | 1e-12 | `W1` must be meaningfully positive against its own terms |
| `TOL_FEAS_OK` | 1e-9 | Least-squares residual indistinguishable from zero at this scale |
| `TOL_FEAS_LOCKED` | 1e-6 | Three decades above `TOL_FEAS_OK`, so solver accuracy cannot move a verdict across the band |
| `FLOOR_MIN_REL` | 1e-12 | Smallest floating-point floor ever reported |
| `TOL_CASE_ONE_SIDED` | 1e-6 | Allowance in one-sided case assertions, for floors only |
| `TOL_CASE_EXACT_GAP` | 1e-6 | Certified gap allowed where both exact fields lie in the discrete spaces |

### 17.6 Status terms

| Term | Values and meaning |
|---|---|
| `solver_status` | `solved`, `almost_solved`, `primal_infeasible`, `dual_infeasible`, `time_limit`, `iteration_limit`, `numerical_error`, `not_run` |
| `classification` | `certified_bound`: passed this checker. `estimate`: everything else, with reasons |
| Inadmissible | Any estimate. It cannot satisfy a gate |
| Bracket `status` | `bracket`: `L <= U`. `order_within_floor`: `U < L <= U + floors`. `contradiction`: `L > U + floors`, which **stops the program** |
| Gap | `(U - L) / ((U + L) / 2)`, `null` unless `L > 0` and the status is `bracket` |
| Feasibility `status` | section 11 |

`primal_infeasible` for Lower or Upper, and `dual_infeasible`, are recorded as results. Write `[L, U]` only
for a certified bracket. Otherwise write `estimate: L=..., U=...`. Report `domain_claim`, criterion,
`flow_rule` and interface law with every result.

## 18. Proposed design changes, pending the owner's decision

None of these is an owner decision. Each states what this contract does meanwhile.

| Id | Proposal | Reason | In this contract |
|---|---|---|---|
| P-1 | Treat continuous P2 as a baseline under test, with a true-3D feasibility gate. Confirm the planned Drucker-Prager angles | Exact isochory at tet vertices over-constrains P2 on unstructured tetrahedra. See the note below | Sections 8, 11. Angles 10, 20 and 30 degrees are placeholders |
| P-2 | Name and version formulations | A new Upper element must not change the meaning of existing cases | Section 8 |
| P-3 | Enforce jump cones on Bernstein coefficients | P2 jumps are quadratic, so facet-vertex checks do not cover the facet. The same applies to `no_tension` and `coulomb` on the Upper side | Section 12; refused in version 0 |
| P-4 | `domain_claim` per bound; a bracket claims what both support; the Lower bound can also extend by zero when weightless with free walls | One bracket can mix claims | Section 14 |
| P-5 | Repair only by an admissible construction | `pi` is infinite off the cone, so a charge pays for nothing. A Lower interior field must itself satisfy equilibrium, so zero works only with no fixed load | Sections 17.1, 17.2 |
| P-6 | Add restrained footing components | A fixed force component is not a restrained velocity. Without it a half-model strip is wrongly constrained | Sections 6, 7 |
| P-7 | State that the `(2 + pi) cu` gate needs fixed walls and a box containing the mechanism | With rollers, the box load need not equal the half-space load | Section 16.4; case geometry |

**Measured context, all estimates, none certified.** Uncommitted checks of 2026-11-18, Clarabel 0.11.1,
gmsh 4.15, `upper.p2c-vertex.v0`, square bonded footing with fixed walls: `rho_feas` at `phi = 0` was
`1.0e-1`, `4.1e-2`, `2.1e-2` on the committed meshes of 243, 463 and 949 tetrahedra, and the conic Upper
problem reported `PrimalInfeasible` on the two coarser meshes. The half-model strip without restrained
components gave `estimate: L=1.18, U=2.15` on 192 tetrahedra against a target of 5.14. With them it gave
`estimate: L=4.68, U=6.17`.

## 19. Amendment log

None yet.
