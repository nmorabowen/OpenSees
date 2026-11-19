# WP-237d solver study, rechecked under the amended checker

Base: reference c7b047674 + checker 94593291e (branch head 0a1ed41a0). Driver: `Ladruno_files/testbed/limit_analysis/solver_study.py` (`table`). Data: `C:\Users\nmora\limit_study_237d_a1`. Date: 2026-11-19. Rows 141; failed jobs 5; failed rechecks 1; gate-admissible certificates 96; classification changes 1.

Every value below is a certified bound read through `read_gate_certificate`, unless it is in the estimates column. Gap = (U - L) / ((U + L) / 2), computed only for a certified pair on one mesh.

## Frozen gap gates (contract `evaluate_case`)

| case | gate | limit | gate mesh run | meshes run | gap gate | case passed | withheld | P0 blocking |
|---|---|---|---|---|---|---|---|---|
| dp_prandtl_cohesion_thin_slab | MAX_GAP_PRANDTL m3 | 0.1 | False | [0] | fail: missing | False | False | True |
| dp_prandtl_surcharge_thin_slab | MAX_GAP_PRANDTL m3 | 0.1 | False | [0] | fail: missing | False | False | True |
| tresca_strip_vm_thin_slab | MAX_GAP_STRIP m3 | 0.05 | True | [0, 1, 2, 3] | fail: violated | False | False | True |
| tresca_strip_weighted_insensitive | MAX_GAP_STRIP m3 | 0.05 | True | [0, 1, 2, 3] | fail: missing | False | False | True |
| upper3d_fixed_phi00 | MAX_GAP_3D_VM m5 | 0.3 | False | [0, 1, 2, 3] | fail: missing | False | False | True |
| upper3d_fixed_phi10 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1] | fail: missing | False | False | True |
| upper3d_fixed_phi20 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1] | fail: missing | False | False | True |
| upper3d_fixed_phi30 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1, 2, 3] | fail: missing | False | False | True |
| upper3d_roller_phi00 | MAX_GAP_3D_VM m5 | 0.3 | False | [0, 1, 2] | fail: missing | False | False | True |
| upper3d_roller_phi10 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1] | fail: missing | False | False | True |
| upper3d_roller_phi20 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1] | fail: missing | False | False | True |
| upper3d_roller_phi30 | MAX_GAP_3D_DP m5 | 0.5 | False | [0, 1, 2] | fail: missing | False | False | True |

## Best certified bounds per mesh

| case | m | Lower | Upper jump | gap jump | Upper P2c | gap P2c | estimates (no certificate) |
|---|---|---|---|---|---|---|---|
| dp_prandtl_cohesion_thin_slab | 0 | 12.3166 | 15.9171 | 0.255048 | 18.1116 | 0.380895 |  |
| dp_prandtl_surcharge_thin_slab | 0 | 5.48284 | 6.79336 | 0.213506 | 7.59208 | 0.32264 |  |
| tresca_strip_vm_thin_slab | 0 | 4.68445 | 5.32841 | 0.128627 | 6.16882 | 0.273534 |  |
| tresca_strip_vm_thin_slab | 1 | 4.93982 | 5.23422 | 0.0578727 | 5.67984 | 0.139367 |  |
| tresca_strip_vm_thin_slab | 2 | 4.90825 | 5.23301 | 0.0640474 | 5.67342 | 0.144622 |  |
| tresca_strip_vm_thin_slab | 3 | 5.00525 | 5.18571 | 0.0354148 | 5.41681 | 0.0789786 |  |
| tresca_strip_weighted_insensitive | 0 | 4.68445 | 5.32841 | 0.128627 | 6.16882 | 0.273535 |  |
| tresca_strip_weighted_insensitive | 1 | 4.9398 | 5.23422 | 0.057876 | 5.67984 | 0.139371 |  |
| tresca_strip_weighted_insensitive | 2 | 4.90825 | 5.23301 | 0.0640474 | 5.67342 | 0.144622 |  |
| tresca_strip_weighted_insensitive | 3 | - | 5.18571 | - | 5.41681 | - |  |
| upper3d_fixed_phi00 | 0 | 2.62056 | 12.3108 | 1.29797 | - | - |  |
| upper3d_fixed_phi00 | 1 | 3.95219 | 7.40493 | 0.608032 | - | - |  |
| upper3d_fixed_phi00 | 2 | - | 7.80215 | - | - | - |  |
| upper3d_fixed_phi00 | 3 | - | - | - | - | - | upper_jump=6.41885 |
| upper3d_fixed_phi10 | 0 | 4.49309 | - | - | 50.7729 | 1.6748 | upper_jump=18.5401 |
| upper3d_fixed_phi10 | 1 | 2.20871 | 14.3829 | 1.46751 | 31.2242 | 1.73574 |  |
| upper3d_fixed_phi20 | 0 | 3.41074 | - | - | 69.1953 | 1.8121 | upper_jump=30.9899 |
| upper3d_fixed_phi20 | 1 | 4.4873 | - | - | 56.0415 | 1.70346 | upper_jump=28.1323 |
| upper3d_fixed_phi30 | 0 | 5.55602 | - | - | 278.194 | 1.92168 | upper_jump=103.392 |
| upper3d_fixed_phi30 | 1 | 8.32297 | - | - | 156.729 | 1.79829 | upper_jump=77.4871 |
| upper3d_fixed_phi30 | 2 | 9.28004 | - | - | - | - | upper_jump=66.4383; upper_control=108.095 |
| upper3d_fixed_phi30 | 3 | 16.2399 | - | - | - | - |  |
| upper3d_roller_phi00 | 0 | 2.62056 | 11.9463 | 1.2804 | 26.9961 | 1.64607 |  |
| upper3d_roller_phi00 | 1 | 3.95212 | 7.3766 | 0.604565 | 20.0405 | 1.34111 |  |
| upper3d_roller_phi00 | 2 | 3.11273 | 7.7883 | 0.857822 | 19.3479 | 1.44566 |  |
| upper3d_roller_phi10 | 0 | 4.4928 | 17.4767 | 1.18199 | 32.2713 | 1.51118 |  |
| upper3d_roller_phi10 | 1 | 2.20372 | 14.3392 | 1.46715 | 25.6691 | 1.68375 |  |
| upper3d_roller_phi20 | 0 | 3.38748 | 28.164 | 1.57055 | 46.4855 | 1.72831 |  |
| upper3d_roller_phi20 | 1 | 4.48346 | - | - | 42.006 | 1.61424 | upper_jump=27.6663 |
| upper3d_roller_phi30 | 0 | - | - | - | 116.39 | - | lower=5.53594; upper_jump=82.1974 |
| upper3d_roller_phi30 | 1 | 8.31374 | - | - | 93.3954 | 1.67304 | upper_jump=71.4125 |
| upper3d_roller_phi30 | 2 | 9.2797 | - | - | 82.4031 | 1.59514 | upper_jump=63.5534 |

## Classification changes under the amended checker

- `upper3d_fixed_phi00.m3.q.upper.p2d-jump-bernstein.v0.clarabel-faer.nodyn`: no_value -> no_check, gate admissible False, reasons None  recheck memory_guard exit 1

## Failed jobs (error, crash, incomplete, memory guard, timeout)

- `upper3d_fixed_phi00.m1.q.upper.p2d-jump-bernstein.v0.scs: incomplete`
- `upper3d_fixed_phi00.m3.q.lower.p1d-vertex.v0.clarabel-qdldl: incomplete`
- `upper3d_fixed_phi00.m3.q.upper.p2d-jump-bernstein.v0.clarabel-faer.nodyn: memory_guard`
- `upper3d_fixed_phi30.m2.q.upper.p2d-jump-bernstein.v0.clarabel-faer: incomplete`
- `upper3d_fixed_phi30.m3.q.upper.p2d-jump-bernstein.v0.clarabel-faer.nodyn: incomplete`
- `upper3d_fixed_phi00.m3.q.upper.p2d-jump-bernstein.v0.clarabel-faer.nodyn: recheck memory_guard`

## Rows

| case | m | formulation | solver | outcome | status | it | solve s | peak GB | old cert | cert | value | kappa_eq |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dp_prandtl_cohesion_thin_slab | 0 | p1d-vertex | clarabel-faer | completed | almost_solved | 23 | 8.181 | 0.137 | certified_bound | certified_bound | 12.3166 | 134.12805310653678 |
| dp_prandtl_cohesion_thin_slab | 0 | p1d-vertex | clarabel-qdldl | completed | almost_solved | 24 | 2.031 | 0.125 | certified_bound | certified_bound | 12.3166 | 134.12567407048275 |
| dp_prandtl_cohesion_thin_slab | 0 | p1d-vertex | scs | completed | solved | 72050 | 111.615 | 0.125 | certified_bound | certified_bound | 12.3165 | 134.12815593348034 |
| dp_prandtl_cohesion_thin_slab | 0 | p2c-vertex | clarabel-faer | completed | solved | 22 | 6.877 | 0.187 | certified_bound | certified_bound | 18.1116 | 4.0000000000006075 |
| dp_prandtl_cohesion_thin_slab | 0 | p2c-vertex | clarabel-qdldl | completed | solved | 22 | 5.045 | 0.182 | certified_bound | certified_bound | 18.1116 | 4.0000000000006075 |
| dp_prandtl_cohesion_thin_slab | 0 | p2c-vertex | scs | completed | solved | 79625 | 175.738 | 0.188 | estimate | estimate | None | 4.000000000000607 |
| dp_prandtl_cohesion_thin_slab | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 26 | 19.492 | 0.316 | certified_bound | certified_bound | 15.9171 | 5.000000000001389 |
| dp_prandtl_cohesion_thin_slab | 0 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 26 | 25.868 | 0.313 | certified_bound | certified_bound | 15.9171 | 5.000000000001388 |
| dp_prandtl_cohesion_thin_slab | 0 | p2d-jump-bernstein | scs | completed | time_limit | 47400 | 600.557 | 0.263 | no_value | no_value | None | None |
| dp_prandtl_surcharge_thin_slab | 0 | p1d-vertex | clarabel-faer | completed | solved | 20 | 7.975 | 0.144 | certified_bound | certified_bound | 5.48284 | 136.83628142922657 |
| dp_prandtl_surcharge_thin_slab | 0 | p2c-vertex | clarabel-faer | completed | solved | 21 | 5.726 | 0.187 | certified_bound | certified_bound | 7.59208 | 4.000000000000607 |
| dp_prandtl_surcharge_thin_slab | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 26 | 21.537 | 0.315 | certified_bound | certified_bound | 6.79336 | 5.000000000001388 |
| tresca_strip_vm_thin_slab | 0 | p1d-vertex | clarabel-faer | completed | solved | 21 | 0.13 | 0.52 | certified_bound | certified_bound | 4.68445 | 52.528272420757645 |
| tresca_strip_vm_thin_slab | 0 | p1d-vertex | clarabel-qdldl | completed | solved | 21 | 0.123 | 0.519 | certified_bound | certified_bound | 4.68445 | 52.52827242077357 |
| tresca_strip_vm_thin_slab | 0 | p1d-vertex | scs | completed | solved | 38050 | 13.544 | 0.522 | certified_bound | certified_bound | 4.68444 | 52.528272409134686 |
| tresca_strip_vm_thin_slab | 0 | p2c-vertex | clarabel-faer | completed | solved | 10 | 0.809 | 0.097 | certified_bound | certified_bound | 6.16882 | 80.76769597976754 |
| tresca_strip_vm_thin_slab | 0 | p2c-vertex | clarabel-qdldl | completed | solved | 10 | 0.229 | 0.11 | certified_bound | certified_bound | 6.16882 | 80.76769597976754 |
| tresca_strip_vm_thin_slab | 0 | p2c-vertex | scs | completed | solved | 850 | 0.369 | 0.09 | certified_bound | certified_bound | 6.16882 | 80.76769597976744 |
| tresca_strip_vm_thin_slab | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 19 | 3.696 | 0.51 | certified_bound | certified_bound | 5.32841 | 141.83873472229686 |
| tresca_strip_vm_thin_slab | 0 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 19 | 1.572 | 0.494 | certified_bound | certified_bound | 5.32841 | 141.83873472229672 |
| tresca_strip_vm_thin_slab | 0 | p2d-jump-bernstein | scs | completed | solved | 2175 | 4.776 | 0.498 | certified_bound | certified_bound | 5.32843 | 141.83873472229644 |
| tresca_strip_vm_thin_slab | 1 | p1d-vertex | clarabel-faer | completed | almost_solved | 19 | 15.563 | 0.205 | certified_bound | certified_bound | 4.9398 | 382.08407987331054 |
| tresca_strip_vm_thin_slab | 1 | p1d-vertex | clarabel-qdldl | completed | almost_solved | 20 | 8.92 | 0.189 | certified_bound | certified_bound | 4.93982 | 382.08410760911437 |
| tresca_strip_vm_thin_slab | 1 | p2c-vertex | clarabel-faer | completed | solved | 13 | 11.22 | 0.612 | certified_bound | certified_bound | 5.67984 | 260.21088431812285 |
| tresca_strip_vm_thin_slab | 1 | p2c-vertex | clarabel-qdldl | completed | solved | 13 | 16.336 | 0.59 | certified_bound | certified_bound | 5.67984 | 260.2051714126353 |
| tresca_strip_vm_thin_slab | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 33 | 82.179 | 0.584 | certified_bound | certified_bound | 5.23422 | 625.0670168364265 |
| tresca_strip_vm_thin_slab | 1 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 33 | 473.522 | 0.684 | certified_bound | certified_bound | 5.23422 | 625.0670168429543 |
| tresca_strip_vm_thin_slab | 2 | p1d-vertex | clarabel-faer | completed | almost_solved | 27 | 9.696 | 0.134 | certified_bound | certified_bound | 4.90825 | 179.73184169678288 |
| tresca_strip_vm_thin_slab | 2 | p2c-vertex | clarabel-faer | completed | solved | 13 | 4.553 | 0.59 | certified_bound | certified_bound | 5.67342 | 110.85832814467378 |
| tresca_strip_vm_thin_slab | 2 | p2d-jump-bernstein | clarabel-faer | completed | solved | 31 | 26.838 | 0.253 | certified_bound | certified_bound | 5.23301 | 258.08704384933753 |
| tresca_strip_vm_thin_slab | 3 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 36 | 43.216 | 0.279 | certified_bound | certified_bound | 5.00525 | 465.1382970747314 |
| tresca_strip_vm_thin_slab | 3 | p1d-vertex | clarabel-faer | completed | numerical_error | 23 | 32.557 | 0.279 | no_value | no_value | None | None |
| tresca_strip_vm_thin_slab | 3 | p2c-vertex | clarabel-faer | completed | solved | 16 | 25.878 | 1.354 | certified_bound | certified_bound | 5.41681 | 342.27341099392555 |
| tresca_strip_vm_thin_slab | 3 | p2d-jump-bernstein | clarabel-faer | completed | solved | 45 | 195.982 | 0.784 | certified_bound | certified_bound | 5.18571 | 731.3158800923799 |
| tresca_strip_weighted_insensitive | 0 | p1d-vertex | clarabel-faer | completed | almost_solved | 21 | 0.118 | 0.523 | certified_bound | certified_bound | 4.68445 | 52.711020272431 |
| tresca_strip_weighted_insensitive | 0 | p2c-vertex | clarabel-faer | completed | solved | 10 | 0.807 | 0.119 | certified_bound | certified_bound | 6.16882 | 80.76769597976742 |
| tresca_strip_weighted_insensitive | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 19 | 3.491 | 0.51 | certified_bound | certified_bound | 5.32841 | 141.83873472229686 |
| tresca_strip_weighted_insensitive | 1 | p1d-vertex | clarabel-faer | completed | almost_solved | 19 | 16.311 | 0.206 | certified_bound | certified_bound | 4.9398 | 382.6768465965732 |
| tresca_strip_weighted_insensitive | 1 | p2c-vertex | clarabel-faer | completed | solved | 13 | 10.979 | 0.615 | certified_bound | certified_bound | 5.67984 | 260.2051714126353 |
| tresca_strip_weighted_insensitive | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 33 | 83.038 | 0.593 | certified_bound | certified_bound | 5.23422 | 625.0670168342538 |
| tresca_strip_weighted_insensitive | 2 | p1d-vertex | clarabel-faer | completed | solved | 29 | 10.49 | 0.137 | certified_bound | certified_bound | 4.90825 | 180.4294728218199 |
| tresca_strip_weighted_insensitive | 2 | p2c-vertex | clarabel-faer | completed | solved | 13 | 5.089 | 0.591 | certified_bound | certified_bound | 5.67342 | 110.85832814467378 |
| tresca_strip_weighted_insensitive | 2 | p2d-jump-bernstein | clarabel-faer | completed | solved | 31 | 31.922 | 0.254 | certified_bound | certified_bound | 5.23301 | 258.08704384913983 |
| tresca_strip_weighted_insensitive | 3 | p1d-vertex | clarabel-faer.nodyn | completed | numerical_error | 29 | 47.038 | 0.28 | no_value | no_value | None | None |
| tresca_strip_weighted_insensitive | 3 | p2c-vertex | clarabel-faer | completed | solved | 16 | 27.757 | 1.4 | certified_bound | certified_bound | 5.41681 | 342.27341099392555 |
| tresca_strip_weighted_insensitive | 3 | p2d-jump-bernstein | clarabel-faer | completed | solved | 45 | 195.35 | 0.796 | certified_bound | certified_bound | 5.18571 | 731.3158800923799 |
| upper3d_fixed_phi00 | 0 | p1d-vertex | clarabel-faer | completed | almost_solved | 16 | 3.359 | 1.114 | certified_bound | certified_bound | 2.62056 | 190.3055758162368 |
| upper3d_fixed_phi00 | 0 | p1d-vertex | clarabel-qdldl | completed | almost_solved | 16 | 1.282 | 1.096 | certified_bound | certified_bound | 2.62056 | 190.30557581648875 |
| upper3d_fixed_phi00 | 0 | p1d-vertex | scs | completed | solved | 41375 | 40.148 | 1.1 | certified_bound | certified_bound | 2.62056 | 190.305575700258 |
| upper3d_fixed_phi00 | 0 | p2c-vertex | clarabel-faer | completed | primal_infeasible | 8 | 0.365 | 0.096 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 0 | p2c-vertex | clarabel-qdldl | completed | primal_infeasible | 7 | 0.132 | 0.005 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 26 | 6.765 | 0.742 | certified_bound | certified_bound | 12.3108 | 54.77141396451653 |
| upper3d_fixed_phi00 | 0 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 26 | 4.563 | 0.721 | certified_bound | certified_bound | 12.3108 | 54.77141396451635 |
| upper3d_fixed_phi00 | 0 | p2d-jump-bernstein | scs | completed | solved | 13775 | 87.722 | 0.725 | certified_bound | certified_bound | 12.3109 | 54.771413999861366 |
| upper3d_fixed_phi00 | 1 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 34 | 9.502 | 0.149 | certified_bound | certified_bound | 3.95219 | 326.93009290752434 |
| upper3d_fixed_phi00 | 1 | p1d-vertex | clarabel-faer | completed | numerical_error | 17 | 5.396 | 0.147 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 1 | p1d-vertex | clarabel-qdldl | completed | numerical_error | 18 | 9.549 | 0.129 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 1 | p1d-vertex | scs | completed | iteration_limit | 100000 | 456.772 | 0.123 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 1 | p2c-vertex | clarabel-faer | completed | primal_infeasible | 7 | 0.723 | 0.122 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 22 | 11.668 | 0.213 | certified_bound | certified_bound | 7.40493 | 130.08493954232972 |
| upper3d_fixed_phi00 | 1 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 22 | 23.187 | 0.198 | certified_bound | certified_bound | 7.40493 | 130.09370751007742 |
| upper3d_fixed_phi00 | 1 | p2d-jump-bernstein | scs | incomplete | None | None | None | None | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 2 | p1d-vertex | clarabel-faer.nodyn | completed | numerical_error | 18 | 13.073 | 0.255 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 2 | p1d-vertex | clarabel-faer | completed | numerical_error | 16 | 12.66 | 0.254 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 2 | p1d-vertex | clarabel-qdldl | completed | numerical_error | 19 | 82.571 | 0.247 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 2 | p2c-vertex | clarabel-faer | completed | primal_infeasible | 7 | 1.477 | 0.174 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 2 | p2d-jump-bernstein | clarabel-faer | completed | solved | 25 | 30.665 | 0.414 | certified_bound | certified_bound | 7.80215 | 155.28106686751107 |
| upper3d_fixed_phi00 | 2 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 25 | 197.492 | 0.423 | certified_bound | certified_bound | 7.80215 | 155.28206911184523 |
| upper3d_fixed_phi00 | 3 | p1d-vertex | clarabel-faer.nodyn | completed | numerical_error | 19 | 64.239 | 0.921 | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 3 | p1d-vertex | clarabel-qdldl | incomplete | None | None | None | None | no_value | no_value | None | None |
| upper3d_fixed_phi00 | 3 | p2d-jump-bernstein | clarabel-faer.nodyn | memory_guard | solved | 38 | 219.862 | 2.867 | no_value | no_check | None | None |
| upper3d_fixed_phi10 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 21 | 4.795 | 1.365 | certified_bound | certified_bound | 4.49309 | 201.49417062549557 |
| upper3d_fixed_phi10 | 0 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 15 | 1.142 | 0.112 | certified_bound | certified_bound | 50.7729 | 5.478907440133874 |
| upper3d_fixed_phi10 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 20 | 7.101 | 0.189 | estimate | estimate | None | 7.8164013936733125 |
| upper3d_fixed_phi10 | 1 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 18 | 10.018 | 0.146 | certified_bound | certified_bound | 2.20871 | 431.27893110042606 |
| upper3d_fixed_phi10 | 1 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 18 | 1.678 | 0.134 | certified_bound | certified_bound | 31.2242 | 8.124038404662437 |
| upper3d_fixed_phi10 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 19 | 10.35 | 0.235 | certified_bound | certified_bound | 14.3829 | 12.529964086284531 |
| upper3d_fixed_phi20 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 23 | 4.543 | 1.475 | certified_bound | certified_bound | 3.41074 | 94.36028087990235 |
| upper3d_fixed_phi20 | 0 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 15 | 1.228 | 0.117 | certified_bound | certified_bound | 69.1953 | 5.83112033148995 |
| upper3d_fixed_phi20 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 19 | 10.14 | 0.193 | estimate | estimate | None | 8.544477161709233 |
| upper3d_fixed_phi20 | 1 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 19 | 8.671 | 0.159 | certified_bound | certified_bound | 4.4873 | 796.128162026684 |
| upper3d_fixed_phi20 | 1 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 17 | 2.661 | 0.14 | certified_bound | certified_bound | 56.0415 | 8.124038404662437 |
| upper3d_fixed_phi20 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 21 | 14.627 | 0.264 | estimate | estimate | None | 12.529964086284462 |
| upper3d_fixed_phi30 | 0 | p1d-vertex | clarabel-faer | completed | almost_solved | 21 | 5.086 | 0.132 | certified_bound | certified_bound | 5.55602 | 380.0135742737935 |
| upper3d_fixed_phi30 | 0 | p1d-vertex | clarabel-qdldl | completed | almost_solved | 22 | 5.095 | 0.112 | certified_bound | certified_bound | 5.55602 | 380.0135741428746 |
| upper3d_fixed_phi30 | 0 | p1d-vertex | scs | completed | solved | 34975 | 65.288 | 0.108 | certified_bound | certified_bound | 5.55595 | 380.01357423371667 |
| upper3d_fixed_phi30 | 0 | p2c-vertex | clarabel-faer | completed | solved | 16 | 1.272 | 0.121 | certified_bound | certified_bound | 278.194 | 5.831034236608605 |
| upper3d_fixed_phi30 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | almost_solved | 17 | 5.67 | 0.205 | estimate | estimate | None | 8.544513282861221 |
| upper3d_fixed_phi30 | 0 | p2d-jump-bernstein | clarabel-faer | completed | almost_solved | 17 | 5.874 | 0.206 | estimate | estimate | None | 8.544513282861221 |
| upper3d_fixed_phi30 | 0 | p2d-jump-bernstein | clarabel-qdldl | completed | almost_solved | 17 | 6.263 | 0.203 | estimate | estimate | None | 8.54451328286127 |
| upper3d_fixed_phi30 | 0 | p2d-jump-bernstein | scs | completed | solved | 5425 | 36.404 | 0.209 | estimate | estimate | None | 8.544513278153204 |
| upper3d_fixed_phi30 | 1 | p1d-vertex | clarabel-faer | completed | almost_solved | 20 | 8.073 | 0.168 | certified_bound | certified_bound | 8.32293 | 564.1915935246618 |
| upper3d_fixed_phi30 | 1 | p1d-vertex | clarabel-qdldl | completed | almost_solved | 19 | 20.014 | 0.149 | certified_bound | certified_bound | 8.32297 | 564.191593488409 |
| upper3d_fixed_phi30 | 1 | p2c-vertex | clarabel-faer | completed | solved | 21 | 2.783 | 0.147 | certified_bound | certified_bound | 156.729 | 8.124038404662446 |
| upper3d_fixed_phi30 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 18 | 11.066 | 0.294 | estimate | estimate | None | 12.529964086284453 |
| upper3d_fixed_phi30 | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 18 | 11.479 | 0.291 | estimate | estimate | None | 12.529964086284453 |
| upper3d_fixed_phi30 | 1 | p2d-jump-bernstein | clarabel-qdldl | completed | solved | 18 | 36.75 | 0.289 | estimate | estimate | None | 12.529964086284444 |
| upper3d_fixed_phi30 | 2 | p1d-vertex | clarabel-faer | completed | almost_solved | 21 | 21.569 | 0.339 | certified_bound | certified_bound | 9.28004 | 4783.963925892136 |
| upper3d_fixed_phi30 | 2 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 22 | 7.346 | 0.223 | estimate | estimate | None | 9.899495276735477 |
| upper3d_fixed_phi30 | 2 | p2c-vertex | clarabel-faer | completed | solved | 22 | 7.455 | 0.223 | estimate | estimate | None | 9.899495276735477 |
| upper3d_fixed_phi30 | 2 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 23 | 35.568 | 0.586 | estimate | estimate | None | 15.52417570965536 |
| upper3d_fixed_phi30 | 2 | p2d-jump-bernstein | clarabel-faer | incomplete | solved | 23 | 35.918 | None | no_value | no_value | None | None |
| upper3d_fixed_phi30 | 3 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 21 | 83.01 | 1.061 | certified_bound | certified_bound | 16.2399 | 2406.9092550158366 |
| upper3d_fixed_phi30 | 3 | p2d-jump-bernstein | clarabel-faer.nodyn | incomplete | almost_solved | 35 | 222.153 | None | no_value | no_value | None | None |
| upper3d_roller_phi00 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 21 | 2.647 | 1.046 | certified_bound | certified_bound | 2.62056 | 211.51634180055686 |
| upper3d_roller_phi00 | 0 | p1d-vertex | clarabel-faer | completed | numerical_error | 19 | 2.635 | 0.111 | no_value | no_value | None | None |
| upper3d_roller_phi00 | 0 | p2c-vertex | clarabel-faer | completed | solved | 10 | 0.486 | 0.119 | certified_bound | certified_bound | 26.9961 | 531.7359399492699 |
| upper3d_roller_phi00 | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 15 | 3.73 | 1.032 | certified_bound | certified_bound | 11.9463 | 53.699412313180716 |
| upper3d_roller_phi00 | 1 | p1d-vertex | clarabel-faer | completed | almost_solved | 24 | 7.825 | 0.149 | certified_bound | certified_bound | 3.95212 | 461.9739764158804 |
| upper3d_roller_phi00 | 1 | p2c-vertex | clarabel-faer | completed | solved | 10 | 1.083 | 0.306 | certified_bound | certified_bound | 20.0405 | 588.8366173986342 |
| upper3d_roller_phi00 | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 24 | 12.147 | 0.212 | certified_bound | certified_bound | 7.3766 | 129.40948661452805 |
| upper3d_roller_phi00 | 2 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 25 | 22.956 | 0.257 | certified_bound | certified_bound | 3.11273 | 2233.114080421804 |
| upper3d_roller_phi00 | 2 | p1d-vertex | clarabel-faer | completed | numerical_error | 20 | 16.686 | 0.256 | no_value | no_value | None | None |
| upper3d_roller_phi00 | 2 | p2c-vertex | clarabel-faer | completed | solved | 10 | 3.14 | 0.917 | certified_bound | certified_bound | 19.3479 | 2282.1659550191366 |
| upper3d_roller_phi00 | 2 | p2d-jump-bernstein | clarabel-faer | completed | solved | 25 | 33.463 | 0.404 | certified_bound | certified_bound | 7.7883 | 153.72102963117212 |
| upper3d_roller_phi10 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 18 | 5.044 | 1.3 | certified_bound | certified_bound | 4.4928 | 335.65995463926595 |
| upper3d_roller_phi10 | 0 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 11 | 0.814 | 0.111 | certified_bound | certified_bound | 32.2713 | 5.477727937536614 |
| upper3d_roller_phi10 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 16 | 6.558 | 0.188 | certified_bound | certified_bound | 17.4767 | 7.813764650927314 |
| upper3d_roller_phi10 | 1 | p1d-vertex | clarabel-faer.nodyn | completed | solved | 21 | 7.588 | 0.146 | certified_bound | certified_bound | 2.20372 | 697.1140515974967 |
| upper3d_roller_phi10 | 1 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 13 | 1.362 | 0.133 | certified_bound | certified_bound | 25.6691 | 8.12403840466244 |
| upper3d_roller_phi10 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 20 | 10.972 | 0.234 | certified_bound | certified_bound | 14.3392 | 12.529964086284478 |
| upper3d_roller_phi20 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 29 | 8.517 | 1.417 | certified_bound | certified_bound | 3.38748 | 123.1318269960571 |
| upper3d_roller_phi20 | 0 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 12 | 0.933 | 0.113 | certified_bound | certified_bound | 46.4855 | 5.831120331489953 |
| upper3d_roller_phi20 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 17 | 6.403 | 0.192 | certified_bound | certified_bound | 28.164 | 8.544477161709237 |
| upper3d_roller_phi20 | 1 | p1d-vertex | clarabel-faer.nodyn | completed | solved | 23 | 9.755 | 0.161 | certified_bound | certified_bound | 4.48346 | 797.7804398308838 |
| upper3d_roller_phi20 | 1 | p2c-vertex | clarabel-faer.nodyn | completed | solved | 12 | 1.606 | 0.136 | certified_bound | certified_bound | 42.006 | 8.124038404662441 |
| upper3d_roller_phi20 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 18 | 11.32 | 0.264 | estimate | estimate | None | 12.529964086284501 |
| upper3d_roller_phi30 | 0 | p1d-vertex | clarabel-faer.nodyn | completed | solved | 22 | 5.827 | 1.747 | estimate | estimate | None | 263.3267492486826 |
| upper3d_roller_phi30 | 0 | p1d-vertex | clarabel-faer | completed | almost_solved | 21 | 6.552 | 1.748 | estimate | estimate | None | 263.32674924751814 |
| upper3d_roller_phi30 | 0 | p2c-vertex | clarabel-faer | completed | solved | 13 | 1.331 | 0.121 | certified_bound | certified_bound | 116.39 | 5.8311336849857724 |
| upper3d_roller_phi30 | 0 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 16 | 5.954 | 0.206 | estimate | estimate | None | 8.545003035713972 |
| upper3d_roller_phi30 | 0 | p2d-jump-bernstein | clarabel-faer | completed | solved | 16 | 8.167 | 0.205 | estimate | estimate | None | 8.545003035713979 |
| upper3d_roller_phi30 | 1 | p1d-vertex | clarabel-faer | completed | almost_solved | 19 | 11.003 | 0.17 | certified_bound | certified_bound | 8.31374 | 607.3427601873001 |
| upper3d_roller_phi30 | 1 | p2c-vertex | clarabel-faer | completed | solved | 14 | 1.962 | 0.143 | certified_bound | certified_bound | 93.3954 | 8.124038404662448 |
| upper3d_roller_phi30 | 1 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 18 | 11.236 | 0.292 | estimate | estimate | None | 12.529964086284465 |
| upper3d_roller_phi30 | 1 | p2d-jump-bernstein | clarabel-faer | completed | solved | 18 | 15.108 | 0.292 | estimate | estimate | None | 12.529964086284444 |
| upper3d_roller_phi30 | 2 | p1d-vertex | clarabel-faer.nodyn | completed | almost_solved | 30 | 32.354 | 0.338 | certified_bound | certified_bound | 9.2797 | 4791.361019348518 |
| upper3d_roller_phi30 | 2 | p1d-vertex | clarabel-faer | completed | numerical_error | 32 | 38.304 | 0.341 | no_value | no_value | None | None |
| upper3d_roller_phi30 | 2 | p2c-vertex | clarabel-faer | completed | solved | 14 | 4.905 | 0.218 | certified_bound | certified_bound | 82.4031 | 9.899495770573674 |
| upper3d_roller_phi30 | 2 | p2d-jump-bernstein | clarabel-faer.nodyn | completed | solved | 21 | 37.455 | 0.586 | estimate | estimate | None | 15.52417599729926 |
| upper3d_roller_phi30 | 2 | p2d-jump-bernstein | clarabel-faer | completed | solved | 21 | 34.403 | 0.585 | estimate | estimate | None | 15.52417599729926 |
