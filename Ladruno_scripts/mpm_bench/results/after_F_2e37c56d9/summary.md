
### Mode F, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (10368) | 0.537 (1%) | 0.839 (4%) | 0.775 (7%) | 1.44 | 1861 |
| bspline2 | elastic | 10000 (10368) | 0.952 (7%) | 2.96 (7%) | 4.64 (15%)* | 4.88 | 1051 |
| bspline2 | sanisand | 10000 (10368) | 0.0438 (1%) | 0.045 (3%) | 0.0425 (3%) | 0.97 (within noise) | 22839 |
| gimp | dp | 10000 (10368) | 0.646 (2%) | 0.937 (4%) | 0.938 (19%)* | 1.45 | 1548 |
| gimp | elastic | 10000 (10368) | 1.34 (0%) | 4.12 (30%)* | 6.39 (13%)* | 4.77 | 747 |
| bspline2 | dp | 100000 (98496) | 0.568 (2%) | 0.914 (0%) | 0.904 (4%) | 1.59 | 1762 |
| bspline2 | elastic | 100000 (98496) | 0.998 (0%) | 2.86 (7%) | 4.36 (1%) | 4.37 | 1002 |
| bspline2 | sanisand | 100000 (98496) | 0.098 (2%) | 0.1 (6%) | 0.0923 (5%) | 0.94 (within noise) | 10202 |
| gimp | dp | 100000 (98496) | 0.683 (1%) | 0.993 (3%) | 1.01 (11%)* | 1.47 | 1464 |
| gimp | elastic | 100000 (98496) | 1.36 (1%) | 4.22 (2%) | 6.62 (9%) | 4.85 | 733 |

### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (9800) | 0.902 (4%) | 1.17 (5%) | 1.07 (16%)* | 1.19 (within noise) | 1109 |
| bspline2 | elastic | 10000 (9800) | 2.59 (4%) | 7.11 (8%) | 8.86 (17%)* | 3.42 | 386 |
| bspline2 | sanisand | 10000 (9800) | 0.204 (1%) | 0.21 (14%)* | 0.202 (4%) | 0.99 (within noise) | 4900 |
| gimp | dp | 10000 (9800) | 0.944 (2%) | 1.17 (3%) | 1.21 (23%)* | 1.28 | 1060 |
| gimp | elastic | 10000 (9800) | 2.99 (0%) | 7.05 (16%)* | 9.67 (19%)* | 3.23 | 334 |
| bspline2 | dp | 100000 (100352) | 0.885 (4%) | 1.15 (5%) | 1.09 (4%) | 1.23 | 1130 |
| bspline2 | elastic | 100000 (100352) | 2.55 (2%) | 6.42 (10%) | 8.18 (4%) | 3.21 | 392 |
| bspline2 | sanisand | 100000 (100352) | 0.339 (1%) | 0.331 (3%) | 0.301 (2%) | 0.89 | 2946 |
| gimp | dp | 100000 (100352) | 0.945 (1%) | 1.15 (4%) | 1.11 (7%) | 1.18 | 1058 |
| gimp | elastic | 100000 (100352) | 3.03 (1%) | 7.2 (6%) | 9.27 (7%) | 3.06 | 331 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3D bspline2 dp N=10000 | 2.2 | 0.2 | 8.3 | 13.7 | 0.4 | 6.3 | 3.4 | 61.7 | 3.3 | 0.0 | 97.5 | 0.4 | 41.0 |
| 3D bspline2 elastic N=10000 | 1.0 | 0.3 | 14.6 | 24.2 | 0.7 | 11.0 | 6.1 | 39.8 | 1.5 | 0.0 | 98.6 | 0.6 | 6.5 |
| 3D bspline2 sanisand N=10000 | 0.2 | 0.0 | 0.7 | 1.1 | 0.0 | 0.5 | 0.3 | 94.2 | 2.9 | 0.0 | 99.8 | 0.0 | 94.6 |
| 3D gimp dp N=10000 | 2.7 | 0.2 | 6.6 | 9.1 | 0.4 | 4.7 | 2.7 | 69.4 | 3.8 | 0.0 | 97.0 | 0.4 | 48.6 |
| 3D gimp elastic N=10000 | 1.4 | 0.5 | 13.4 | 18.9 | 0.7 | 9.7 | 5.4 | 47.4 | 1.9 | 0.1 | 98.0 | 0.7 | 8.5 |
| 3D bspline2 dp N=100000 | 2.7 | 0.2 | 8.8 | 14.5 | 0.3 | 6.6 | 3.6 | 59.4 | 3.5 | 0.0 | 97.1 | 0.2 | 40.7 |
| 3D bspline2 elastic N=100000 | 1.5 | 0.3 | 15.5 | 25.5 | 0.5 | 11.6 | 6.4 | 36.6 | 1.6 | 0.0 | 98.1 | 0.3 | 6.6 |
| 3D bspline2 sanisand N=100000 | 0.5 | 0.0 | 1.5 | 2.5 | 0.1 | 1.2 | 0.6 | 86.9 | 6.6 | 0.0 | 99.5 | 0.0 | 88.5 |
| 3D gimp dp N=100000 | 3.1 | 0.3 | 7.7 | 10.4 | 0.3 | 5.2 | 2.9 | 65.6 | 4.2 | 0.0 | 96.6 | 0.2 | 46.8 |
| 3D gimp elastic N=100000 | 2.0 | 0.5 | 15.2 | 20.6 | 0.6 | 10.5 | 5.8 | 42.1 | 2.2 | 0.0 | 97.5 | 0.4 | 9.1 |
| 2D bspline2 dp N=10000 | 4.9 | 0.2 | 4.9 | 5.8 | 0.5 | 3.1 | 1.8 | 73.1 | 5.1 | 0.0 | 94.8 | 0.5 | 59.3 |
| 2D bspline2 elastic N=10000 | 2.0 | 0.6 | 14.0 | 16.6 | 1.6 | 8.9 | 5.2 | 47.9 | 2.1 | 0.1 | 97.2 | 1.1 | 10.8 |
| 2D bspline2 sanisand N=10000 | 1.2 | 0.0 | 1.1 | 1.3 | 0.1 | 0.7 | 0.4 | 81.9 | 13.1 | 0.0 | 98.7 | 0.1 | 89.2 |
| 2D gimp dp N=10000 | 4.8 | 0.2 | 4.8 | 4.5 | 0.5 | 2.8 | 1.5 | 75.4 | 4.8 | 0.0 | 94.9 | 0.4 | 61.1 |
| 2D gimp elastic N=10000 | 2.1 | 0.7 | 15.0 | 14.0 | 1.6 | 8.6 | 4.8 | 49.9 | 2.0 | 0.1 | 97.0 | 1.2 | 11.8 |
| 2D bspline2 dp N=100000 | 5.6 | 0.2 | 4.8 | 5.8 | 0.5 | 3.1 | 1.8 | 72.3 | 5.5 | 0.0 | 94.2 | 0.3 | 58.5 |
| 2D bspline2 elastic N=100000 | 3.8 | 0.6 | 14.0 | 16.9 | 1.6 | 9.1 | 5.2 | 44.6 | 3.2 | 0.0 | 95.5 | 0.9 | 12.0 |
| 2D bspline2 sanisand N=100000 | 2.1 | 0.1 | 1.9 | 2.2 | 0.2 | 1.2 | 0.7 | 69.4 | 22.1 | 0.0 | 97.8 | 0.1 | 93.9 |
| 2D gimp dp N=100000 | 5.8 | 0.2 | 4.5 | 4.1 | 0.5 | 2.6 | 1.4 | 74.8 | 5.8 | 0.0 | 94.0 | 0.3 | 61.4 |
| 2D gimp elastic N=100000 | 4.4 | 0.8 | 14.6 | 13.2 | 1.7 | 8.3 | 4.3 | 48.1 | 3.8 | 0.0 | 94.8 | 1.0 | 13.4 |

### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)

| row | threads | T ratio | verdict | G(n) |
|---|---|---|---|---|
| F 3D bspline2 dp N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 3D bspline2 elastic N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 3D bspline2 sanisand N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp dp N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 3D gimp elastic N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 3D bspline2 dp N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 3D bspline2 elastic N=100000 | 1 | 0.97 | slower | 0.97 |
| F 3D bspline2 sanisand N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 3D gimp dp N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 3D gimp elastic N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 dp N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 2D bspline2 elastic N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 2D bspline2 sanisand N=10000 | 1 | 1.00 | within noise | 1.00 |
| F 2D gimp dp N=10000 | 1 | 1.02 | within noise | 1.02 |
| F 2D gimp elastic N=10000 | 1 | 1.02 | faster | 1.02 |
| F 2D bspline2 dp N=100000 | 1 | 0.97 | within noise | 0.97 |
| F 2D bspline2 elastic N=100000 | 1 | 0.98 | within noise | 0.98 |
| F 2D bspline2 sanisand N=100000 | 1 | 1.03 | within noise | 1.03 |
| F 2D gimp dp N=100000 | 1 | 1.00 | within noise | 1.00 |
| F 2D gimp elastic N=100000 | 1 | 1.00 | within noise | 1.00 |

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

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('2e37c56d9a8feac8926575b7e0ac63f9812a04cf', '2e37c56d9a8feac8926575b7e0ac63f9812a04cf')]

11 row(s) with spread > 10 % owed a 5-repeat re-run: F/2D/bspline2/dp/N10000/t8, F/2D/bspline2/elastic/N10000/t8, F/2D/bspline2/sanisand/N10000/t4, F/2D/gimp/dp/N10000/t8, F/2D/gimp/elastic/N10000/t4, F/2D/gimp/elastic/N10000/t8, F/3D/bspline2/elastic/N10000/t8, F/3D/gimp/dp/N10000/t8, F/3D/gimp/dp/N100000/t8, F/3D/gimp/elastic/N10000/t4, F/3D/gimp/elastic/N10000/t8
