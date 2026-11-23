# Lower-bound stall study and the `presolve` / `relgap1e-3` presets

Driver `Ladruno_files/testbed/limit_analysis/solver_study.py` (this commit; base 17ee1b4fc), reference c7b047674, checker 94593291e, clarabel 0.11.1 faer, qoco 0.3.2, date 2026-11-19, workstation (other workers were running; walls are indicative). Data: `limit_study_237d_presolve` (after), `limit_study_237d`, `...\diag`, `limit_study_237d_prandtl` (before). Every value below is a checker-certified bound unless the classification says otherwise.

## Rank deficiency of the Lower equality block (`rankstudy`)

Dependent rows = eigenvalues of A_eq A_eq^T below sigma, counted by Sylvester inertia of an unpivoted SuperLU factorisation of A A^T - sigma I at sigma = 1e-10 and 1e-9 of max|A A^T| (two counts). The count agrees with the checker's independent rank on strip m0 (2732) and cohesion m0 (11324). Spread = max/min equality-row 2-norm.

| case | m | tets | eq rows | dependent | fixed cols | empty rows dropped | dependent after presolve | row-norm spread |
|---|---|---|---|---|---|---|---|---|
| upper3d_roller_phi00 | 0 | 243 | 5478 | [122, 122] | 1129 | 122 | [0, 0] | 60 |
| upper3d_fixed_phi00 | 0 | 243 | 4686 | [42, 42] | 414 | 42 | [0, 0] | 60 |
| tresca_strip_vm_thin_slab | 0 | 192 | 4168 | [272, 272] | 1164 | 168 | [104, 104] | 1881 |
| dp_prandtl_cohesion_thin_slab | 0 | 768 | 16816 | [1064, 1064] | 4428 | 600 | [464, 464] | 333 |
| upper3d_fixed_phi00 | 1 | 463 | 9072 | [110, 110] | 831 | 110 | [0, 0] | 149 |
| upper3d_roller_phi00 | 1 | 463 | 10380 | [214, 214] | 2037 | 214 | [0, 0] | 149 |
| upper3d_fixed_phi00 | 2 | 949 | 18702 | [122, 122] | 1106 | 122 | [0, 0] | 221 |
| upper3d_roller_phi00 | 2 | 949 | 20946 | [266, 268] | 3207 | 266 | [0, 2] | 221 |
| upper3d_fixed_phi00 | 3 | 2840 | 60069 | [790, 790] | 5059 | 790 | [0, 0] | 7274 |
| upper3d_roller_phi00 | 3 | 2840 | 61029 | [878, 882] | 5932 | 878 | [0, 2] | 7274 |
| tresca_strip_vm_thin_slab | 3 | 3072 | 67348 | [4160, 4160] | 16956 | 2208 | [1952, 1952] | 9216 |
| dp_prandtl_cohesion_thin_slab | 2 | 3072 | 67420 | [4176, 4176] | 17052 | 2224 | [1952, 1952] | 815 |
| dp_prandtl_cohesion_thin_slab | 1 | 6144 | 131764 | [7240, 7240] | 18088 | 2312 | [4928, 4928] | 1006 |
| dp_prandtl_cohesion_thin_slab | 3 | 12288 | 270004 | [16544, 16544] | 66876 | 8544 | [8000, 8000] | 2206 |

The full matrix (equality and cone rows) has no column null space in any case measured (m0 to m2): no free stress mode exists.

## Before and after (clarabel-faer, lower.p1d-vertex.v0)

| case | m | preset | data | status | it | solve s | peak GB | rel gap | certificate | value |
|---|---|---|---|---|---|---|---|---|---|---|
| dp_prandtl_cohesion_thin_slab | 1 | default | limit_study_237d_prandtl | numerical_error | 27 | 104.133 | 0.702 | 2.04e-04 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 1 | base | diag | numerical_error | 27 | 149.548 | 0.701 | 2.04e-04 | - | - |
| dp_prandtl_cohesion_thin_slab | 1 | nodyn | limit_study_237d_prandtl | numerical_error | 33 | 152.448 | 0.703 | 1.53e-04 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 2 | default | limit_study_237d_prandtl | numerical_error | 28 | 44.794 | 0.316 | 6.57e-05 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 2 | nodyn | limit_study_237d_prandtl | almost_solved | 30 | 81.993 | 0.324 | 3.18e-05 | certified_bound | 13.3896 |
| dp_prandtl_cohesion_thin_slab | 3 | default | limit_study_237d_prandtl | numerical_error | 26 | 177.131 | 1.065 | 6.72e-04 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 3 | nodyn | limit_study_237d_prandtl | numerical_error | 44 | 458.172 | 1.065 | 8.36e-05 | no_value | - |
| tresca_strip_vm_thin_slab | 3 | default | limit_study_237d | numerical_error | 23 | 32.557 | 0.279 | 5.39e-05 | no_value | - |
| tresca_strip_vm_thin_slab | 3 | base | diag | numerical_error | 23 | 32.585 | 0.28 | 5.39e-05 | - | - |
| tresca_strip_vm_thin_slab | 3 | nodyn | limit_study_237d | almost_solved | 36 | 43.216 | 0.279 | 4.56e-06 | certified_bound | 5.00525 |
| tresca_strip_vm_thin_slab | 3 | nodyn | diag | almost_solved | 36 | 53.156 | 0.281 | 4.56e-06 | - | - |
| upper3d_fixed_phi00 | 1 | default | limit_study_237d | numerical_error | 17 | 5.396 | 0.147 | 1.42e-04 | no_value | - |
| upper3d_fixed_phi00 | 1 | base | diag | numerical_error | 17 | 6.027 | 0.149 | 1.42e-04 | - | - |
| upper3d_fixed_phi00 | 1 | nodyn | limit_study_237d | almost_solved | 34 | 9.502 | 0.149 | 3.36e-05 | certified_bound | 3.95219 |
| upper3d_fixed_phi00 | 1 | nodyn | diag | almost_solved | 34 | 12.207 | 0.148 | 3.36e-05 | - | - |
| upper3d_fixed_phi00 | 2 | default | limit_study_237d | numerical_error | 16 | 12.66 | 0.254 | 1.06e-03 | no_value | - |
| upper3d_fixed_phi00 | 2 | nodyn | limit_study_237d | numerical_error | 18 | 13.073 | 0.255 | 5.31e-04 | no_value | - |
| upper3d_fixed_phi00 | 3 | nodyn | limit_study_237d | numerical_error | 19 | 64.239 | 0.921 | 9.67e-05 | no_value | - |
| upper3d_roller_phi00 | 1 | default | limit_study_237d | almost_solved | 24 | 7.825 | 0.149 | 7.12e-06 | certified_bound | 3.95212 |
| upper3d_roller_phi00 | 1 | base | diag | almost_solved | 24 | 9.093 | 0.149 | 7.12e-06 | - | - |
| upper3d_roller_phi00 | 2 | default | limit_study_237d | numerical_error | 20 | 16.686 | 0.256 | 8.62e-05 | no_value | - |
| upper3d_roller_phi00 | 2 | base | diag | numerical_error | 20 | 19.534 | 0.256 | 8.62e-05 | - | - |
| upper3d_roller_phi00 | 2 | nodyn | limit_study_237d | almost_solved | 25 | 22.956 | 0.257 | 4.03e-05 | certified_bound | 3.11273 |
| upper3d_roller_phi00 | 3 | base | diag | almost_solved | 19 | 100.257 | 0.922 | 4.14e-07 | - | - |
| dp_prandtl_cohesion_thin_slab | 1 | presolve | presolve | numerical_error | 26 | 119.382 | 0.599 | 1.53e-04 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 1 | relgap1e-3 | presolve | almost_solved | 27 | 155.763 | 0.784 | 2.04e-04 | certified_bound | 13.6397 |
| dp_prandtl_cohesion_thin_slab | 2 | presolve | presolve | numerical_error | 27 | 34.083 | 0.273 | 9.61e-05 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 2 | relgap1e-3 | presolve | almost_solved | 28 | 57.197 | 0.318 | 6.57e-05 | certified_bound | 13.3894 |
| dp_prandtl_cohesion_thin_slab | 3 | presolve | presolve | numerical_error | 26 | 167.094 | 0.813 | 3.83e-04 | no_value | - |
| dp_prandtl_cohesion_thin_slab | 3 | relgap1e-3 | presolve | almost_solved | 26 | 223.575 | 1.365 | 6.72e-04 | certified_bound | 13.9588 |
| tresca_strip_vm_thin_slab | 3 | presolve | presolve | numerical_error | 20 | 26.071 | 0.269 | 1.99e-04 | no_value | - |
| tresca_strip_vm_thin_slab | 3 | relgap1e-3 | presolve | almost_solved | 23 | 36.754 | 0.283 | 5.39e-05 | certified_bound | 5.00521 |
| upper3d_fixed_phi00 | 1 | presolve | presolve | numerical_error | 18 | 5.649 | 0.157 | 5.81e-05 | no_value | - |
| upper3d_fixed_phi00 | 1 | presolve_relgap | presolve | almost_solved | 18 | 6.953 | 0.156 | 5.81e-05 | certified_bound | 3.95218 |
| upper3d_fixed_phi00 | 1 | presolve_relgap1e-3 | presolve | almost_solved | 18 | 5.477 | 0.157 | 5.81e-05 | certified_bound | 3.95218 |
| upper3d_fixed_phi00 | 1 | relgap | presolve | almost_solved | 17 | 6.926 | 0.15 | 1.42e-04 | certified_bound | 3.95215 |
| upper3d_fixed_phi00 | 1 | relgap1e-3 | presolve | almost_solved | 17 | 5.739 | 0.15 | 1.42e-04 | certified_bound | 3.95215 |
| upper3d_fixed_phi00 | 2 | presolve | presolve | numerical_error | 21 | 16.375 | 0.267 | 1.11e-04 | no_value | - |
| upper3d_fixed_phi00 | 2 | presolve_relgap1e-3 | presolve | almost_solved | 21 | 27.623 | 0.268 | 1.11e-04 | certified_bound | 3.11272 |
| upper3d_fixed_phi00 | 2 | relgap1e-3 | presolve | numerical_error | 16 | 16.388 | 0.257 | 1.06e-03 | no_value | - |
| upper3d_fixed_phi00 | 3 | presolve | presolve | numerical_error | 15 | 67.862 | 0.897 | 1.67e-04 | no_value | - |
| upper3d_fixed_phi00 | 3 | presolve_relgap1e-3 | presolve | almost_solved | 15 | 56.677 | 0.897 | 1.67e-04 | certified_bound | 4.54794 |
| upper3d_fixed_phi00 | 3 | relgap1e-3 | presolve | almost_solved | 15 | 67.194 | 0.925 | 1.30e-04 | certified_bound | 4.54798 |
| upper3d_fixed_phi00 | 4 | nodyn_relgap1e-3 | presolve | almost_solved | 23 | 55.849 | 0.51 | 1.91e-04 | certified_bound | 3.34843 |
| upper3d_fixed_phi00 | 4 | presolve | presolve | numerical_error | 14 | 27.666 | 0.495 | 3.59e-03 | no_value | - |
| upper3d_fixed_phi00 | 4 | presolve_relgap1e-3 | presolve | numerical_error | 14 | 27.274 | 0.494 | 3.59e-03 | no_value | - |
| upper3d_fixed_phi00 | 4 | relgap1e-3 | presolve | numerical_error | 17 | 37.901 | 0.51 | 1.13e-03 | no_value | - |
| upper3d_roller_phi00 | 1 | presolve | presolve | almost_solved | 21 | 6.984 | 0.145 | 9.59e-06 | certified_bound | 3.95212 |
| upper3d_roller_phi00 | 1 | presolve_relgap1e-3 | presolve | almost_solved | 21 | 6.147 | 0.145 | 9.59e-06 | certified_bound | 3.95212 |
| upper3d_roller_phi00 | 1 | relgap1e-3 | presolve | almost_solved | 24 | 9.615 | 0.152 | 7.12e-06 | certified_bound | 3.95212 |
| upper3d_roller_phi00 | 2 | presolve | presolve | numerical_error | 19 | 16.215 | 0.262 | 1.55e-04 | no_value | - |
| upper3d_roller_phi00 | 2 | presolve_relgap1e-3 | presolve | almost_solved | 19 | 14.968 | 0.262 | 1.55e-04 | certified_bound | 3.11271 |
| upper3d_roller_phi00 | 2 | relgap1e-3 | presolve | almost_solved | 20 | 18.503 | 0.259 | 8.62e-05 | certified_bound | 3.11272 |
| upper3d_roller_phi00 | 3 | presolve | presolve | solved | 21 | 81.745 | 0.904 | 5.09e-09 | certified_bound | 4.54801 |
| upper3d_roller_phi00 | 3 | presolve_relgap1e-3 | presolve | solved | 21 | 85.775 | 0.97 | 5.09e-09 | certified_bound | 4.54801 |
| upper3d_roller_phi00 | 3 | relgap1e-3 | presolve | almost_solved | 19 | 89.244 | 0.925 | 4.14e-07 | certified_bound | 4.54801 |
| upper3d_roller_phi00 | 4 | presolve_relgap1e-3 | presolve | almost_solved | 17 | 30.69 | 0.48 | 5.39e-04 | certified_bound | 3.34839 |
| upper3d_roller_phi00 | 4 | relgap1e-3 | presolve | almost_solved | 19 | 39.679 | 0.59 | 3.65e-04 | certified_bound | 3.34841 |

## QOCO 0.3.2 (adapter fixed to `setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, **settings)`)

| case | m | preset | data | status | it | solve s | peak GB | rel gap | certificate | value |
|---|---|---|---|---|---|---|---|---|---|---|
| dp_prandtl_cohesion_thin_slab | 1 | p1d-vertex | prandtl | numerical_error | 31 | 307.127 | 0.506 | 3.92e-03 | no_value | - |
| tresca_strip_vm_thin_slab | 0 | p1d-vertex | prandtl | solved | 22 | 0.158 | 0.52 | 3.24e-07 | certified_bound | 4.68445 |
| tresca_strip_vm_thin_slab | 0 | p2d-jump-bernstein | prandtl | solved | 17 | 2.093 | 0.495 | 3.84e-07 | certified_bound | 5.32841 |
| tresca_strip_vm_thin_slab | 2 | p1d-vertex | prandtl | almost_solved | 33 | 1.764 | 0.107 | 7.39e-07 | certified_bound | 4.90825 |
| tresca_strip_vm_thin_slab | 2 | p2d-jump-bernstein | prandtl | solved | 25 | 26.639 | 0.214 | 5.86e-07 | certified_bound | 5.233 |

## Reading

- The Lower equality block is rank deficient (42 to 16,544 rows). In 3D every dependent row is an equality row emptied by singleton-row elimination; the presolve removes all of them (0 left, 2 at the larger shift on roller m2/m3). The thin slabs keep further dependent combinations (104 to 8,000).
- Removing them and Ruiz scaling does not stop the failure: 9 of 11 presolve-only runs still end in numerical_error, with primal and dual residuals at 1e-9 or below. The interior-point method breaks down with a relative gap of 5e-5 to 4e-3, above Clarabel's default reduced (AlmostSolved) gap of 5e-5. QOCO stalls the same way on cohesion m1, so the solver is not the cause.
- `relgap1e-3` (reduced gap 1e-3) certifies every case whose stall gap is below 1e-3: 3D phi00 fixed and roller m1 to m3, roller m4, strip m3, cohesion m1 to m3. `presolve_relgap1e-3` also certifies fixed m2 (stall gap 1.06e-3 without presolve, 1.1e-4 with). Fixed m4 needs `nodyn_relgap1e-3`.
- The hypothesis that dependent rows make the KKT system singular does not hold as the cause: the rows exist, but the breakdown persists without them.
