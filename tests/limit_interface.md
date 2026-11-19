# Limit analysis: the P0 interface (WP-237a)

**Status:** proposed for freeze; **not frozen**. The orchestrator freezes this file. After the freeze a
change is an amendment logged in section 19, never a silent edit. Section 18 separates accepted owner
amendments from orchestrator P0 rulings.

**Who reads it.** The Python reference author, the independent checker author, the external reference
author and later the C++ author. Each implements from this file and `tests/limit_cases.txt`, without
reading the other side's code. `tests/test_limit_contract.py` holds the schema validators and the
convention checks. It contains no assembler and no checker.

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

**Frozen set.** The box is `[-2, 2] x [-2, 2] x [-2, 0]`, with a bonded unit-square footing conforming on
`z = 0`. Walls and base are all `fixed` or all `symmetry` (rollers). The meshes, in case order, are
`upper3d_a` (243 tets), `upper3d_b` (463), `upper3d_c` (949), `upper3d_d` (2,840, graded toward the footing edge),
then `upper3d_a` refined once (1,944) and `upper3d_d` refined once (22,720). The first four are independent
non-structured sizes. The refined meshes form the nested pairs `(0, 4)` and `(3, 5)`.

### 11.1 Kinematic test: rank and residual

Let `E` be the equality rows of the formulation: supports, rollers, rigid footing, isochory at `alpha = 0`
vertices, and zero normal jump on every Bernstein coefficient where `tau = 0`. Let `A` be the set of `(v, V)`
satisfying the boundary and footing rows and the query normalisation.

- **Rank.** `rank` counts singular values of `E` above `TOL_RANK * sigma_max`. `nullity = n_unknowns - rank`.
  `footing_modes` is the rank of the kernel restricted to `V`. A dense factorisation may be skipped when too
  large; then `rank`, `nullity` and `footing_modes` are `null`, and the residual test still applies.
- **Residual.** With `c_ei` the flow violation at element vertices, `c_FB` the jump violation on Bernstein
  coefficients, `w_ei = V_e / 4` and `w_FB = A_F / (6 h_F)`, `h_F` the longest edge of the facet:

```text
c_ei = |theta|                                  alpha = 0
c_ei = max(0, 3 sqrt2 alpha |e| - theta)        alpha > 0
c_FB = max(0, tau |c_t| - n . c_B)              tau > 0;   |n . c_B| when tau = 0

rho_feas = min_A sqrt( sum w_ei c_ei^2 + sum w_FB c_FB^2 )
         / min_A sqrt( sum w_ei |D_ei|_F^2 + sum w_FB |c_B|^2 )
```

For the control the facet sums are absent. With every `alpha = 0` both minima are sparse linear least-squares
problems, so the result does not depend on a conic solver.

| `status` | Condition |
|---|---|
| `formulation_feasible` | `rho_feas <= TOL_FEAS_OK`; a witness in `A` passed `c` recomputed independently of the producer; `nullity` and `footing_modes`, where computed, are at least 1 |
| `formulation_infeasible` | `rho_feas >= TOL_FEAS_LOCKED` |
| `inconclusive` | anything else |

A failure is a **formulation result**, recorded with its numbers per mesh, wall and angle. It may not be hidden
by relaxing a tolerance, by reading `AlmostSolved` as success, or by quoting a conic value whose flow rule is
violated. `formulation_infeasible` says this space yields no Upper bound on that mesh. It says nothing against
the mechanics.

### 11.2 Quality test

For each case and the binding Upper formulation, all of these must hold, using **certified** values only:

1. `kinematic_status` is `formulation_feasible` on every listed mesh;
2. a certified bracket exists on mesh 5;
3. its gap is at most `MAX_GAP_3D_VM` for `phi = 0`, or `MAX_GAP_3D_DP` otherwise;
4. `U_fine <= U_coarse + TOL_MONOTONE |U_coarse|` along each nested pair;
5. `L_fine >= L_coarse - TOL_MONOTONE |L_coarse|` along each nested pair.

Section 13 proves the monotonicity for these formulations. A violation therefore signals a defect or an
unconverged candidate, and is never averaged away. The limits are requirements set in advance. If a formulation
misses them, the gate fails and the failure is recorded; the limit does not move.

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
sections 3.5 and 3.6 add refinement and mid-edge nodes. The four committed files are unstructured tetrahedral meshes
of the box `[-2, 2] x [-2, 2] x [-2, 0]`, with the square `[-0.5, 0.5]^2` conforming on `z = 0`, made with gmsh 4.15
(Delaunay 3D, one thread, seed 1). `a`, `b` and `c` use uniform target sizes 1.0, 0.8 and 0.58. `d` grades from 0.08
within 0.1 of the footing edge to 1.0 beyond 1.6. Their 243, 463, 949 and 2,840 tetrahedra are data. No test or
implementation depends on gmsh.

### 15.3 Result: `<stem>.<query>.<bound>.result.json`, `.npz`

Keys: `format` = `ladruno-limit-result`, `format_version`, `arrays`, `problem_sha256`, `query_id`, `bound`,
`formulation`, `producer` `{name, version}`, `solver` `{name, version, settings}`, `solver_status`,
`objective`, `value`, `wall_time_s`, `result_sha256`.

| `bound` | Required arrays | Optional, never trusted |
|---|---|---|
| `lower` | `sigma (E,4,6)`, `resultant (n_f,6)` | `interior_sigma (E,4,6)` with number `interior_value`; `equilibrium_multiplier (E,3)` |
| `upper`, `upper.p2c-vertex.v0` | `velocity (N,3)`, `footing_velocity (n_f,6)` | `interior_velocity (N,3)` and `interior_footing_velocity (n_f,6)`; `strain_rate (E,4,6)`; `dissipation (E,4)` |
| `upper`, `upper.p2d-jump-bernstein.v0` | `element_velocity (E,10,3)`, `footing_velocity (n_f,6)` | `interior_element_velocity (E,10,3)` and `interior_footing_velocity (n_f,6)`; `strain_rate (E,4,6)`; `dissipation (E,4)` |

`interior_value` is a top-level key, present exactly when `interior_sigma` is. With `solver_status` of
`solved` or `almost_solved`, `value` and the arrays are required. Otherwise `value` and `objective` are
`null` and no field arrays are stored: a time limit gives **no number**. A `lower` result carrying velocity
arrays, or the reverse, is refused.

### 15.4 Program: `<stem>.<query>.<bound>.program.json`, `.npz`

Canonical form: minimise `c^T x` subject to `A x + s = b`, `s` in `K`. `K` is the zero cone, then the listed
second-order cones in order. `value = objective_sign * c^T x`.

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
   `symmetry`, rows `T1`, `T2`; none for others;
4. ray forms: for each footing and unrestrained component, `R - lambda F1 = F0`;
5. cones, per `(e, i)`: rows `[sqrt2 alpha one^T ; -Q]` on `sigma_ei`, right side `(sqrt2 k, 0)`.

**Upper control.** `x = [v (3N), V (6 n_f), (r, z)_(0,0)...]`. `objective_sign = +1`.
`c = (V_e/4) sqrt2 k_e` on each `r`, minus the `W0` coefficients. Rows, in order:

1. flow, per `(e, i)`, 6 rows: `B_ei v - sqrt2 alpha r one - Q^T z = 0`;
2. `fixed`: for each fixed node ascending, 3 rows `v = 0`;
3. `symmetry`: per set in order, nodes ascending, fixed nodes omitted, `n^T v = 0`;
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
4. `symmetry`: per set, per `(e, a)` ascending, those in `fixed` omitted, one row;
5. `footing`: per footing, per `(e, a)` ascending, 3 rows;
6. restrained components, then `W1 = 1` or `V_fc = d_fc`;
7. cones: `SOC6` per `(e, i)`, then `SOC3` per `(F, B)`. `cones = [{dim: 6, count: 4 E}, {dim: 3, count: 6 F}]`.

### 15.5 Feasibility record: `<stem>.<query>.feasibility.json`

Keys: `format` = `ladruno-limit-feasibility`, `format_version`, `problem_sha256`, `query_id`, `formulation`,
`status`, `rho_feas`, `numerator`, `denominator`, `n_unknowns`, `n_rows`, `rank`, `nullity`, `footing_modes`,
`method`, `witness_checked`, `producer`, `tolerances`. `rho_feas = numerator / denominator` and `status` follows
section 11.1.

### 15.6 Certificate and bracket

`<stem>.<query>.<bound>.certificate.json`, written by the checker only:

| Key | Content |
|---|---|
| `format`, `format_version` | `ladruno-limit-certificate`, `0` |
| `problem_sha256`, `result_sha256`, `query_id`, `bound`, `formulation` | identity |
| `classification` | `certified_bound` or `estimate` |
| `reasons` | codes of section 17.6; empty exactly when certified |
| `value_candidate`, `value_certified`, `floor` | numbers; the last two are `null` for an estimate |
| `domain_claim`, `criteria`, `flow_rule`, `interface_laws` | section 14; sorted lists; the literal `associated` |
| `checks` | metrics of section 17.7 |
| `checker` | `{name, version}` |
| `tolerances` | the table of section 17.8 |

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

`A x = g` is the affine system of the bound, with the rows of sections 17.3 and 17.4. `A_hat` is `A` with each row divided
by its largest absolute entry.

- **Backward error.** `r = A x - g` and `omega = max_i |r_i| / (|A| |x| + |g|)_i`, evaluated in double precision.
  A row with zero denominator must have `r_i = 0`. By Oettli and Prager, `x` exactly satisfies a system whose entries of
  `A` and `g` differ relatively by at most `omega`.
- **Conditioning.** `kappa_eq = sigma_max(A_hat) / sigma_min+(A_hat)`, the smallest singular value above
  `TOL_RANK sigma_max`, known within a factor 2. Rows may be dependent.
- **Move.** `projection_move = |x_P - x|_inf / max(|x|_inf, S_ref)`, with `S_ref` replaced by the largest velocity for
  Upper fields.

### 17.2 Equality acceptance

```text
tau_eq = TOL_EQ_ROUNDING_FACTOR * u * max(1, kappa_eq),     u = 2^-53

accept only if   omega <= tau_eq   and   kappa_eq * omega <= TOL_FLOOR_MAX
```

The threshold is a rounding level scaled by the conditioning of this system. It is not a solver tolerance and no solver
setting enters it. If `kappa_eq * tau_eq > TOL_FLOOR_MAX`, the system cannot support a claim at this precision and the
result is an estimate with `equality_ill_conditioned`.

### 17.3 Lower procedure

1. Check hashes, query, formulation, straightness and volumes. Recompute the value from the candidate and compare it with
   the result, within `TOL_VALUE_MATCH`.
2. **Project.** `x = (lambda, sigma)`, or `sigma` for `support`. Rows: equilibrium, facet continuity, free and traction
   boundaries, symmetry rows, and unrestrained footing resultants. Use the minimum-norm correction in `A_hat`, with
   iterative refinement. Require `projection_move <= TOL_PROJECTION_MOVE`.
3. Evaluate `g_ei = sqrt(J2) + alpha I1 - k` from **principal stresses**. If `max g / S_ref > TOL_CONE_ACCEPT`, stop.
4. **Repair inequalities**, if any `g_ei > -TOL_CONE_MARGIN S_ref`, with a strictly admissible field `(sigma_int, value_int)`
   satisfying the same equalities and `g_int <= -TOL_INTERIOR_MIN S_ref`:
   `mu = max (g + m) / (g - g_int)`, `m = TOL_CONE_MARGIN S_ref`, `sigma_new = (1 - mu) sigma + mu sigma_int`,
   `value_new = (1 - mu) value + mu value_int`. Require `mu <= TOL_REPAIR_MAX`.
5. `sigma_int = 0`, `value_int = 0` is available only when `b0`, `t0` and `F0` are all zero and every `k > 0`. Otherwise it
   must be supplied as `interior_sigma` and is itself projected and checked. If none exists, the result stays an estimate
   with `no_interior_field`.
6. **Recheck the final field from scratch:** `omega`, every traction and resultant row, and `g_ei <= 0` from principal
   stresses.

### 17.4 Upper procedure

1. As Lower step 1, with the value recomputed from the candidate velocity.
2. **Project.** `x = (v, V)`. Rows: fixed nodes, roller nodes, rigid footprint, restrained components, `V = d` for `support`,
   `theta = 0` at every `alpha = 0` element vertex, and `n . c_B = 0` on every coefficient where `tau = 0`. Work normalisation
   is **not** a row. Require `projection_move <= TOL_PROJECTION_MOVE`.
3. Compute `d_ei` and `c_B` with the checker's own operators. `rho_flow = max c_ei / max |D_ei|_F` and
   `rho_jump_cone = max c_FB / max |c_B|`, over the remaining inequality conditions. If either exceeds `TOL_FLOW_ACCEPT`,
   stop with `flow_violation_too_large`.
4. **Repair inequalities** only by an admissible construction. `pi` and `pi_J` are `+infinity` outside their domains, so a
   finite dissipation charge is not a repair.
   - von Mises: after projection no inequality remains. **No repair is needed.**
   - Drucker-Prager: blend with a supplied field `(v_int, V_int)` satisfying the same equalities, the same `V = d` for
     `support`, and `c < 0` strictly: `v_new = (1 - mu) v + mu v_int`, `mu = max c / (c - c_int)`, `mu <= TOL_REPAIR_MAX`.
     `c_ei` and `c_FB` are convex and positively homogeneous, so the blend is admissible.
   - If no construction exists, the result stays an estimate with `no_admissible_correction`.
5. Require `W1 >= TOL_WORK_MIN * W_scale` for ray forms, with `W_scale` the sum of absolute work terms.
6. `D = sum (V_e/4) pi(d_ei) + sum (A_F/6) pi_J(c_B)`. Ray forms: `U = (D - W0) / W1`, so no exact normalisation is needed.
   `support`: `U = D - W0`.
7. **Recheck the final field from scratch:** `omega`, all conditions of step 2, `c <= 0` at `alpha > 0` conditions, and `U`.

### 17.5 What "certified" means

A **certified bound** means all of the following held for the final field `x_c`, as evaluated in double precision by the
independent checker:

1. `x_c` satisfies every affine condition exactly for a problem whose operator and load entries differ relatively from the
   stated ones by at most `omega <= tau_eq`;
2. every inequality passes on the final field, with the margins above: no flow or jump condition is merely charged for;
3. the value is the theorem's value for `x_c`, recomputed by the checker;
4. the reported `floor = max(FLOOR_MIN_REL, kappa_eq * omega) * |value_certified|` estimates, to first order, the effect of
   that perturbation and of rounding on the value.

It is **not** an interval-arithmetic proof. "Rigorous" means rigorous for the stated idealisation, up to this floor.

### 17.6 Outcome codes and negative fixtures

Reason codes: `no_candidate`, `hash_mismatch`, `query_mismatch`, `formulation_not_admitted`, `formulation_not_applicable`,
`curved_element`, `degenerate_element`, `boundary_violation`, `projection_move_too_large`, `equality_residual_above_floor`,
`equality_ill_conditioned`, `cone_violation_too_large`, `no_interior_field`, `interior_field_invalid`, `repair_too_large`,
`flow_violation_too_large`, `no_admissible_correction`, `work_nonpositive`, `value_mismatch`, `unsupported_domain_claim`.

`tests/data/limit/negative_fixtures.json` lists one fixture per false-green item of the harness. `stage` is `schema` (the
format validators refuse it), `checker` (the full procedure must end as an estimate with one of the `expect` codes) or
`final_check`. A `final_check` fixture goes straight to the recheck, **without projection**, and must be refused with
`equality_residual_above_floor`. It sits at three times `tau_eq`, using the recorded `kappa_reference`.

### 17.7 Scaled metrics

`L_ref` is the bounding-box diagonal.
`S_ref = max(max k, max |sigma|, max |t0|, |value| max |t1|, L_ref max(|b0| + |value| |b1|), max_f |F_f| / A_f)`, using the
components of `sigma` for Lower and omitting that term for Upper. For `support`, the `|value|` terms are absent.

`checks` holds these keys, `null` where not applicable. Diagnostics of the candidate before projection: `rho_div`
(`|div sigma + b|_inf h_e / S_ref`), `rho_jump`, `rho_bc`, `rho_resultant` (force over `S_ref A_f`, moment over
`S_ref A_f ell`). Final-field values: `projection_move`, `omega_eq`, `kappa_eq`, `tau_eq`, `rho_cone`, `repair_mu`,
`straightness_max`, `rho_flow`, `rho_jump_cone`, `work_reference`.

### 17.8 Named tolerances

Values are provisional until 237c and 237d calibrate them on the exact cases. A value changes only by amendment, with
evidence, and **never to rescue a run**.

| Name | Value | Rationale |
|---|---|---|
| `TOL_GEOM_STRAIGHT` | 1e-9 | Relative to `h_e`. Admits meshes written with 12 significant digits at `L/h` up to 1e3 |
| `TOL_GEOM_VOLUME` | 1e-10 | `V_e >= tol h_e^3`. A regular tetrahedron has 0.118; below this, gradients amplify rounding by 1e10 |
| `TOL_GEOM_PLANE` | 1e-9 | Normal deviation and offset, relative to `L_ref` |
| `TOL_RANK` | 1e-10 | Separates numerical zero singular values of row-scaled operators from the smallest real ones |
| `TOL_EQ_ROUNDING_FACTOR` | 64 | Multiple of `u` in `tau_eq`. A row has at most about 60 entries, so its evaluation error is about that many `u` |
| `TOL_FLOOR_MAX` | 1e-6 | Largest first-order value effect `kappa_eq * omega` a certificate may carry |
| `TOL_PROJECTION_MOVE` | 1e-6 | A certificate speaks about this candidate. A larger move means a different field |
| `TOL_VALUE_MATCH` | 1e-6 | The result's value must be the value of its own field. Catches normalisation and quadrature faults |
| `TOL_CONE_ACCEPT` | 1e-6 | Largest scaled yield violation worth repairing |
| `TOL_CONE_MARGIN` | 1e-12 | Clearance over eigenvalue rounding, a few thousand `u` of `S_ref` |
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
| `MAX_GAP_STRIP` | 0.05 | Requirement: resolves `(2 + pi) cu` to about 2.5% either side on the listed finest mesh |
| `MAX_GAP_PRANDTL` | 0.10 | Requirement for the frictional strips, whose fans are wider and less resolved |
| `MAX_GAP_3D_VM` | 0.30 | Requirement for true-3D von Mises on the 22,720-tet mesh; useful to judge a 3D result |
| `MAX_GAP_3D_DP` | 0.50 | Requirement for true-3D Drucker-Prager on the 22,720-tet mesh |

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

## 18. Decision record

### 18.1 Accepted owner amendments

The owner accepted these for ADR 237 on 2026-11-18.

| Id | Amendment | In this contract |
|---|---|---|
| A1 | True-3D locking gate. Continuous P2 is a tested control, not an assumed element | Sections 8, 10, 11 |
| A2 | Whole-facet jump checks. P2 jumps are quadratic | Section 12 |
| A3 | `domain_claim` per bound. A bracket claims only what both bounds support | Section 14 |
| A4 | Certificate repair only by an admissible construction | Sections 17.3, 17.4 |

### 18.2 Orchestrator P0 rulings

These are orchestrator rulings for P0, recorded on 2026-11-18. **They are not owner decisions.**

| Id | Ruling | In this contract |
|---|---|---|
| R1 | Versioned formulation ids, with cases that do not change meaning | Section 8 |
| R2 | `upper.p2d-jump-bernstein.v0` is the P0 Upper formulation; `upper.p2c-vertex.v0` stays a control | Sections 8, 12 |
| R3 | The frozen true-3D set: `phi = 0, 10, 20, 30`, fixed and roller walls, non-structured meshes at three or more sizes, rank, residual and certified quality | Section 11; `upper_gate_3d` cases |
| R4 | Restrained footing components | Sections 6, 7 |
| R5 | Strip gate: maximum gap, fixed `x` and `z` walls, extents of at least `1.5 B` and depth `B / sqrt2`, and the box lemma. Here `z` is the depth axis and `1.5 B` is measured from the centreline | Section 16.4 |
| R6 | Zero extension may also serve a Lower bound when weightless with free walls; zero is an interior Lower field only without fixed loads | Sections 14, 17.3 |
| R7 | Certificate order: affine projection, inequality repair, full recheck; equality acceptance is rounding-level and scaled by conditioning | Section 17 |
| R8 | Gravity set: exact shear layer as the zero-gravity killer, exact weighted Drucker-Prager column, von Mises column as a two-sided sandwich only | Section 16.5 |

### 18.3 Evidence notes

All values here are **estimates**, from uncommitted checks of 2026-11-18 with Clarabel 0.11.1, unless marked as a probe.
None is a certified bound.

- `Ladruno_files/testbed/limit_analysis/p2_isochory_probe.py` is a tracked **linear-algebra probe, not a bound**. Its output
  file is the record for the control's locking and for the jump space's mobility.
- That probe finds the control locked with **fixed** walls (`rho_feas` 1.0e-1, 4.1e-2, 2.1e-2 on 243, 463, 949 tets) but
  mobile with **roller** walls (`rho_feas` near 1e-11). Roller cases therefore do not discriminate locking. The jump
  space is mobile for both wall types (`rho_feas` near 1e-11).
- Half-model strip without restrained components: `estimate: L=1.18, U=2.15` on 192 tets, against 5.14.
- Tresca strip, 3,072 tets, jump formulation: `estimate: L=5.005, U=5.186`.
- Prandtl strips, `phi = 20 deg`, 3,072 tets, jump formulation: cohesion `estimate: L=13.39, U=15.40` against 14.83;
  surcharge `estimate: L=5.873, U=6.604` against 6.399.
- True 3D, `phi = 0`, fixed walls, jump formulation: `estimate: L=3.348, U=9.244` on 1,944 uniform-refined tets;
  `estimate: L=4.548, U=6.419` on the 2,840-tet graded mesh. No value on the 22,720-tet gate mesh exists yet.
- True 3D, `phi = 30 deg`, fixed walls, 1,944 tets: `estimate: L=9.65, U=168.9`. Whether `MAX_GAP_3D_DP` is reachable is
  **not established**.

## 19. Amendment log

None yet.
