
### Mode F, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (10368) | 0.543 (2%) | 0.858 (17%)* | 0.951 (3%) | 1.75 | 1841 |
| bspline2 | elastic | 10000 (10368) | 0.971 (16%)* | 2.99 (2%) | 4.53 (12%)* | 4.66 | 1030 |
| bspline2 | sanisand | 10000 (10368) | 0.0439 (2%) | 0.0457 (2%) | 0.0463 (1%) | 1.05 | 22778 |
| gimp | dp | 10000 (10368) | 0.639 (3%) | 0.944 (1%) | 0.949 (7%) | 1.48 | 1564 |
| gimp | elastic | 10000 (10368) | 1.37 (1%) | 3.88 (19%)* | 5.39 (10%)* | 3.93 | 728 |
| bspline2 | dp | 100000 (98496) | 0.581 (2%) | 0.912 (4%) | 0.97 (5%) | 1.67 | 1722 |
| bspline2 | elastic | 100000 (98496) | 1.04 (2%) | 2.83 (3%) | 3.99 (8%) | 3.84 | 964 |
| bspline2 | sanisand | 100000 (98496) | 0.0992 (16%)* | 0.106 (35%)* | 0.106 (9%) | 1.07 (within noise) | 10085 |
| gimp | dp | 100000 (98496) | 0.683 (2%) | 1.01 (1%) | 1.07 (2%) | 1.57 | 1464 |
| gimp | elastic | 100000 (98496) | 1.38 (4%) | 3.87 (3%) | 5.59 (14%)* | 4.06 | 727 |

### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (9800) | 0.89 (7%) | 1.15 (3%) | 1.15 (32%)* | 1.29 (within noise) | 1124 |
| bspline2 | elastic | 10000 (9800) | 2.67 (1%) | 6.57 (9%) | 6.69 (4%) | 2.51 | 375 |
| bspline2 | sanisand | 10000 (9800) | 0.201 (2%) | 0.214 (1%) | 0.213 (6%) | 1.06 (within noise) | 4965 |
| gimp | dp | 10000 (9800) | 0.928 (3%) | 1.18 (6%) | 1.16 (10%)* | 1.25 | 1077 |
| gimp | elastic | 10000 (9800) | 3.02 (3%) | 7.62 (2%) | 7.82 (31%)* | 2.59 | 331 |
| bspline2 | dp | 100000 (100352) | 0.908 (1%) | 1.14 (2%) | 1.16 (1%) | 1.28 | 1101 |
| bspline2 | elastic | 100000 (100352) | 2.62 (1%) | 6.46 (1%) | 8.28 (6%) | 3.17 | 382 |
| bspline2 | sanisand | 100000 (100352) | 0.341 (67%)* | 0.372 (1%) | 0.374 (1%) | 1.10 (within noise) | 2934 |
| gimp | dp | 100000 (100352) | 0.939 (2%) | 1.15 (1%) | 1.19 (5%) | 1.27 | 1065 |
| gimp | elastic | 100000 (100352) | 3.11 (2%) | 7.62 (8%) | 9.84 (4%) | 3.17 | 322 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3D bspline2 dp N=10000 | 2.5 | 0.2 | 8.5 | 12.6 | 0.4 | 6.3 | 3.4 | 62.3 | 3.5 | 0.0 | 97.2 | 0.4 | 42.1 |
| 3D bspline2 elastic N=10000 | 1.2 | 0.3 | 15.1 | 22.9 | 0.7 | 11.4 | 6.1 | 40.2 | 1.7 | 0.0 | 98.4 | 0.7 | 6.8 |
| 3D bspline2 sanisand N=10000 | 0.2 | 0.0 | 0.7 | 1.0 | 0.0 | 0.5 | 0.3 | 94.3 | 2.9 | 0.0 | 99.8 | 0.0 | 94.7 |
| 3D gimp dp N=10000 | 3.0 | 0.2 | 6.5 | 8.5 | 0.4 | 4.7 | 2.4 | 69.5 | 4.2 | 0.0 | 96.7 | 0.4 | 50.2 |
| 3D gimp elastic N=10000 | 1.5 | 0.5 | 13.8 | 17.7 | 0.7 | 10.0 | 5.2 | 47.7 | 2.1 | 0.1 | 97.9 | 0.7 | 9.3 |
| 3D bspline2 dp N=100000 | 2.7 | 0.2 | 9.0 | 13.4 | 0.3 | 6.7 | 3.6 | 60.3 | 3.7 | 0.0 | 97.1 | 0.2 | 41.4 |
| 3D bspline2 elastic N=100000 | 1.6 | 0.3 | 15.9 | 24.0 | 0.6 | 11.9 | 6.4 | 37.1 | 1.9 | 0.0 | 98.1 | 0.4 | 7.0 |
| 3D bspline2 sanisand N=100000 | 0.5 | 0.0 | 1.5 | 2.3 | 0.1 | 1.1 | 0.6 | 87.4 | 6.4 | 0.0 | 99.5 | 0.0 | 88.8 |
| 3D gimp dp N=100000 | 3.2 | 0.3 | 7.7 | 9.5 | 0.3 | 5.3 | 2.7 | 66.4 | 4.3 | 0.0 | 96.5 | 0.2 | 48.8 |
| 3D gimp elastic N=100000 | 2.1 | 0.5 | 15.7 | 19.5 | 0.7 | 10.7 | 5.5 | 42.2 | 2.6 | 0.0 | 97.3 | 0.4 | 9.5 |
| 2D bspline2 dp N=10000 | 5.2 | 0.2 | 4.8 | 5.6 | 0.5 | 3.0 | 1.6 | 73.0 | 5.4 | 0.0 | 94.4 | 0.5 | 59.3 |
| 2D bspline2 elastic N=10000 | 1.7 | 0.5 | 14.2 | 16.6 | 1.6 | 8.8 | 4.7 | 48.3 | 2.3 | 0.1 | 97.5 | 1.1 | 12.1 |
| 2D bspline2 sanisand N=10000 | 1.5 | 0.1 | 1.1 | 1.3 | 0.1 | 0.7 | 0.4 | 81.5 | 13.2 | 0.0 | 98.4 | 0.1 | 89.1 |
| 2D gimp dp N=10000 | 5.5 | 0.2 | 4.7 | 4.2 | 0.5 | 2.7 | 1.4 | 74.6 | 5.7 | 0.0 | 94.2 | 0.5 | 62.1 |
| 2D gimp elastic N=10000 | 1.9 | 0.7 | 15.1 | 13.0 | 1.6 | 8.6 | 4.3 | 50.8 | 2.6 | 0.1 | 97.1 | 1.2 | 12.5 |
| 2D bspline2 dp N=100000 | 5.6 | 0.2 | 4.9 | 5.7 | 0.5 | 3.1 | 1.6 | 72.3 | 5.6 | 0.0 | 94.2 | 0.3 | 60.6 |
| 2D bspline2 elastic N=100000 | 3.1 | 0.6 | 14.2 | 16.9 | 1.6 | 9.2 | 4.9 | 45.1 | 3.7 | 0.0 | 96.3 | 0.9 | 12.7 |
| 2D bspline2 sanisand N=100000 | 2.5 | 0.1 | 1.8 | 2.2 | 0.2 | 1.2 | 0.6 | 69.3 | 22.0 | 0.0 | 97.4 | 0.1 | 83.2 |
| 2D gimp dp N=100000 | 6.0 | 0.2 | 4.5 | 3.7 | 0.5 | 2.6 | 1.3 | 74.8 | 6.0 | 0.0 | 93.8 | 0.3 | 63.9 |
| 2D gimp elastic N=100000 | 3.5 | 0.8 | 14.8 | 11.9 | 1.8 | 8.6 | 4.2 | 48.8 | 4.6 | 0.0 | 95.6 | 1.0 | 14.4 |

### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)

| row | threads | T ratio | verdict | G(n) |
|---|---|---|---|---|
| F 3D bspline2 dp N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 3D bspline2 elastic N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 sanisand N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp dp N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp elastic N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 dp N=100000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 elastic N=100000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 sanisand N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp dp N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp elastic N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 2D bspline2 dp N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 elastic N=10000 | 1 | 1.04 | within noise | 1.04 |
| F 2D bspline2 sanisand N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 2D gimp dp N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 2D gimp elastic N=10000 | 1 | 1.03 | within noise | 1.03 |
| F 2D bspline2 dp N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 elastic N=100000 | 1 | 1.01 | within noise | 1.01 |
| F 2D bspline2 sanisand N=100000 | 1 | 1.04 | within noise | 1.04 |
| F 2D gimp dp N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 2D gimp elastic N=100000 | 1 | 1.02 | within noise | 1.02 |

### Determinism (particle-array digests against 1 thread)

| row | threads | identical to t=1 |
|---|---|---|
| F 2D bspline2 dp N=10000 | 1 | - |
| F 2D bspline2 dp N=10000 | 4 | yes |
| F 2D bspline2 dp N=10000 | 8 | yes |
| F 2D bspline2 dp N=100000 | 1 | - |
| F 2D bspline2 dp N=100000 | 4 | yes |
| F 2D bspline2 dp N=100000 | 8 | yes |
| F 2D bspline2 elastic N=10000 | 1 | - |
| F 2D bspline2 elastic N=10000 | 4 | yes |
| F 2D bspline2 elastic N=10000 | 8 | yes |
| F 2D bspline2 elastic N=100000 | 1 | - |
| F 2D bspline2 elastic N=100000 | 4 | yes |
| F 2D bspline2 elastic N=100000 | 8 | yes |
| F 2D bspline2 sanisand N=10000 | 1 | - |
| F 2D bspline2 sanisand N=10000 | 4 | yes |
| F 2D bspline2 sanisand N=10000 | 8 | yes |
| F 2D bspline2 sanisand N=100000 | 1 | - |
| F 2D bspline2 sanisand N=100000 | 4 | yes |
| F 2D bspline2 sanisand N=100000 | 8 | yes |
| F 2D gimp dp N=10000 | 1 | - |
| F 2D gimp dp N=10000 | 4 | yes |
| F 2D gimp dp N=10000 | 8 | yes |
| F 2D gimp dp N=100000 | 1 | - |
| F 2D gimp dp N=100000 | 4 | yes |
| F 2D gimp dp N=100000 | 8 | yes |
| F 2D gimp elastic N=10000 | 1 | - |
| F 2D gimp elastic N=10000 | 4 | yes |
| F 2D gimp elastic N=10000 | 8 | yes |
| F 2D gimp elastic N=100000 | 1 | - |
| F 2D gimp elastic N=100000 | 4 | yes |
| F 2D gimp elastic N=100000 | 8 | yes |
| F 3D bspline2 dp N=10000 | 1 | - |
| F 3D bspline2 dp N=10000 | 4 | yes |
| F 3D bspline2 dp N=10000 | 8 | yes |
| F 3D bspline2 dp N=100000 | 1 | - |
| F 3D bspline2 dp N=100000 | 4 | yes |
| F 3D bspline2 dp N=100000 | 8 | yes |
| F 3D bspline2 elastic N=10000 | 1 | - |
| F 3D bspline2 elastic N=10000 | 4 | yes |
| F 3D bspline2 elastic N=10000 | 8 | yes |
| F 3D bspline2 elastic N=100000 | 1 | - |
| F 3D bspline2 elastic N=100000 | 4 | yes |
| F 3D bspline2 elastic N=100000 | 8 | yes |
| F 3D bspline2 sanisand N=10000 | 1 | - |
| F 3D bspline2 sanisand N=10000 | 4 | yes |
| F 3D bspline2 sanisand N=10000 | 8 | yes |
| F 3D bspline2 sanisand N=100000 | 1 | - |
| F 3D bspline2 sanisand N=100000 | 4 | yes |
| F 3D bspline2 sanisand N=100000 | 8 | yes |
| F 3D gimp dp N=10000 | 1 | - |
| F 3D gimp dp N=10000 | 4 | yes |
| F 3D gimp dp N=10000 | 8 | yes |
| F 3D gimp dp N=100000 | 1 | - |
| F 3D gimp dp N=100000 | 4 | yes |
| F 3D gimp dp N=100000 | 8 | yes |
| F 3D gimp elastic N=10000 | 1 | - |
| F 3D gimp elastic N=10000 | 4 | yes |
| F 3D gimp elastic N=10000 | 8 | yes |
| F 3D gimp elastic N=100000 | 1 | - |
| F 3D gimp elastic N=100000 | 4 | yes |
| F 3D gimp elastic N=100000 | 8 | yes |

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('ab255f0ecbc43db7500030dbcc66add15f6fa6ce', 'ab255f0ecbc43db7500030dbcc66add15f6fa6ce')]

12 row(s) with spread > 10 % owed a 5-repeat re-run: F/2D/bspline2/dp/N10000/t8, F/2D/bspline2/sanisand/N100000/t1, F/2D/gimp/dp/N10000/t8, F/2D/gimp/elastic/N10000/t8, F/3D/bspline2/dp/N10000/t4, F/3D/bspline2/elastic/N10000/t1, F/3D/bspline2/elastic/N10000/t8, F/3D/bspline2/sanisand/N100000/t1, F/3D/bspline2/sanisand/N100000/t4, F/3D/gimp/elastic/N10000/t4, F/3D/gimp/elastic/N10000/t8, F/3D/gimp/elastic/N100000/t8
