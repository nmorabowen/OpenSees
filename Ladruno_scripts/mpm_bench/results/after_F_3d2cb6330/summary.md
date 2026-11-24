
### Mode F, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (10368) | 0.566 (1%) | 0.868 (1%) | 0.972 (6%) | 1.72 | 1766 |
| bspline2 | elastic | 10000 (10368) | 0.996 (1%) | 3.16 (3%) | 5.04 (5%) | 5.06 | 1004 |
| bspline2 | sanisand | 10000 (10368) | 0.0449 (1%) | 0.0458 (2%) | 0.046 (7%) | 1.02 (within noise) | 22258 |
| gimp | dp | 10000 (10368) | 0.673 (2%) | 0.987 (1%) | 1.06 (9%) | 1.57 | 1486 |
| gimp | elastic | 10000 (10368) | 1.4 (2%) | 4.34 (1%) | 6.85 (7%) | 4.91 | 716 |
| bspline2 | dp | 100000 (98496) | 0.603 (3%) | 0.939 (2%) | 1.03 (2%) | 1.71 | 1660 |
| bspline2 | elastic | 100000 (98496) | 1.06 (0%) | 3 (1%) | 4.8 (4%) | 4.54 | 947 |
| bspline2 | sanisand | 100000 (98496) | 0.101 (2%) | 0.109 (2%) | 0.109 (2%) | 1.08 | 9870 |
| gimp | dp | 100000 (98496) | 0.714 (1%) | 1.04 (2%) | 1.1 (3%) | 1.54 | 1401 |
| gimp | elastic | 100000 (98496) | 1.44 (2%) | 4.28 (3%) | 7.18 (3%) | 5.00 | 697 |

### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (9800) | 0.927 (1%) | 1.18 (3%) | 1.22 (5%) | 1.31 | 1079 |
| bspline2 | elastic | 10000 (9800) | 2.67 (1%) | 6.7 (5%) | 7.44 (23%) | 2.78 | 374 |
| bspline2 | sanisand | 10000 (9800) | 0.181 (45%) | 0.206 (41%) | 0.194 (29%) | 1.07 (within noise) | 5514 |
| gimp | dp | 10000 (9800) | 0.972 (1%) | 1.21 (3%) | 1.28 (2%) | 1.32 | 1029 |
| gimp | elastic | 10000 (9800) | 3 (2%) | 7.17 (1%) | 7.57 (30%) | 2.53 | 334 |
| bspline2 | dp | 100000 (100352) | 0.909 (2%) | 1.19 (0%) | 1.19 (3%) | 1.31 | 1100 |
| bspline2 | elastic | 100000 (100352) | 2.67 (1%) | 6.79 (3%) | 9.47 (9%) | 3.55 | 375 |
| bspline2 | sanisand | 100000 (100352) | 0.347 (68%) | 0.37 (2%) | 0.384 (4%) | 1.11 (within noise) | 2884 |
| gimp | dp | 100000 (100352) | 0.973 (2%) | 1.2 (1%) | 1.24 (2%) | 1.27 | 1027 |
| gimp | elastic | 100000 (100352) | 3.16 (1%) | 7.81 (2%) | 11.1 (6%) | 3.50 | 316 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3D bspline2 dp N=10000 | 2.3 | 0.2 | 8.8 | 13.0 | 0.4 | 6.4 | 3.4 | 62.1 | 3.2 | 0.0 | 97.5 | 0.4 | 42.7 |
| 3D bspline2 elastic N=10000 | 1.1 | 0.3 | 15.4 | 22.7 | 0.7 | 11.4 | 6.0 | 40.5 | 1.5 | 0.0 | 98.5 | 0.6 | 6.2 |
| 3D bspline2 sanisand N=10000 | 0.2 | 0.0 | 0.7 | 1.0 | 0.0 | 0.5 | 0.3 | 94.3 | 2.9 | 0.0 | 99.8 | 0.0 | 95.1 |
| 3D gimp dp N=10000 | 2.7 | 0.3 | 6.8 | 8.7 | 0.4 | 4.8 | 2.5 | 69.4 | 3.9 | 0.0 | 96.9 | 0.4 | 50.7 |
| 3D gimp elastic N=10000 | 1.5 | 0.5 | 14.0 | 17.7 | 0.7 | 10.0 | 5.2 | 47.8 | 1.8 | 0.0 | 97.9 | 0.6 | 8.6 |
| 3D bspline2 dp N=100000 | 2.7 | 0.2 | 9.1 | 13.6 | 0.3 | 6.9 | 3.6 | 59.7 | 3.6 | 0.0 | 97.1 | 0.2 | 42.2 |
| 3D bspline2 elastic N=100000 | 1.5 | 0.3 | 15.9 | 23.8 | 0.6 | 12.1 | 6.3 | 37.3 | 1.7 | 0.0 | 98.1 | 0.3 | 6.8 |
| 3D bspline2 sanisand N=100000 | 0.4 | 0.0 | 1.5 | 2.3 | 0.1 | 1.2 | 0.6 | 87.3 | 6.5 | 0.0 | 99.5 | 0.0 | 89.5 |
| 3D gimp dp N=100000 | 3.2 | 0.3 | 8.0 | 9.7 | 0.3 | 5.4 | 2.8 | 65.9 | 4.3 | 0.0 | 96.5 | 0.2 | 49.9 |
| 3D gimp elastic N=100000 | 2.0 | 0.5 | 15.9 | 19.3 | 0.7 | 10.8 | 5.5 | 42.6 | 2.4 | 0.0 | 97.4 | 0.4 | 9.0 |
| 2D bspline2 dp N=10000 | 4.5 | 0.2 | 5.0 | 5.9 | 0.6 | 3.1 | 1.7 | 73.9 | 4.6 | 0.0 | 95.2 | 0.4 | 60.9 |
| 2D bspline2 elastic N=10000 | 1.7 | 0.5 | 14.3 | 16.7 | 1.6 | 8.9 | 4.8 | 48.4 | 1.9 | 0.1 | 97.5 | 1.1 | 11.3 |
| 2D bspline2 sanisand N=10000 | 1.0 | 0.0 | 1.0 | 1.1 | 0.1 | 0.6 | 0.3 | 84.2 | 11.5 | 0.0 | 98.9 | 0.1 | 90.5 |
| 2D gimp dp N=10000 | 4.5 | 0.2 | 4.9 | 4.4 | 0.5 | 2.8 | 1.4 | 75.9 | 4.5 | 0.0 | 95.1 | 0.4 | 62.7 |
| 2D gimp elastic N=10000 | 1.9 | 0.7 | 15.1 | 13.3 | 1.6 | 8.6 | 4.3 | 50.9 | 2.2 | 0.1 | 97.1 | 1.2 | 11.9 |
| 2D bspline2 dp N=100000 | 5.5 | 0.2 | 5.0 | 5.8 | 0.5 | 3.2 | 1.7 | 72.2 | 5.6 | 0.0 | 94.3 | 0.3 | 62.1 |
| 2D bspline2 elastic N=100000 | 2.8 | 0.6 | 14.6 | 16.7 | 1.6 | 9.2 | 4.8 | 45.6 | 3.2 | 0.0 | 96.6 | 0.9 | 12.5 |
| 2D bspline2 sanisand N=100000 | 2.0 | 0.1 | 1.9 | 2.2 | 0.2 | 1.2 | 0.6 | 68.8 | 23.0 | 0.0 | 97.9 | 0.1 | 84.6 |
| 2D gimp dp N=100000 | 5.7 | 0.2 | 4.6 | 3.7 | 0.5 | 2.6 | 1.3 | 75.1 | 5.8 | 0.0 | 94.1 | 0.3 | 65.2 |
| 2D gimp elastic N=100000 | 3.4 | 0.8 | 15.0 | 12.2 | 1.7 | 8.5 | 4.3 | 49.4 | 3.7 | 0.0 | 95.8 | 1.0 | 14.4 |

### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)

| row | threads | T ratio | verdict | G(n) |
|---|---|---|---|---|
| F 3D bspline2 dp N=10000 | 1 | 1.04 | faster | 1.04 |
| F 3D bspline2 elastic N=10000 | 1 | 1.04 | faster | 1.04 |
| F 3D bspline2 sanisand N=10000 | 1 | 1.03 | faster | 1.03 |
| F 3D gimp dp N=10000 | 1 | 1.05 | faster | 1.05 |
| F 3D gimp elastic N=10000 | 1 | 1.03 | within noise | 1.03 |
| F 3D bspline2 dp N=100000 | 1 | 1.05 | faster | 1.05 |
| F 3D bspline2 elastic N=100000 | 1 | 1.03 | faster | 1.03 |
| F 3D bspline2 sanisand N=100000 | 1 | 1.02 | within noise | 1.02 |
| F 3D gimp dp N=100000 | 1 | 1.05 | faster | 1.05 |
| F 3D gimp elastic N=100000 | 1 | 1.04 | within noise | 1.04 |
| F 2D bspline2 dp N=10000 | 1 | 1.03 | within noise | 1.03 |
| F 2D bspline2 elastic N=10000 | 1 | 1.04 | within noise | 1.04 |
| F 2D bspline2 sanisand N=10000 | 1 | 0.89 | within noise | 0.89 |
| F 2D gimp dp N=10000 | 1 | 1.05 | faster | 1.05 |
| F 2D gimp elastic N=10000 | 1 | 1.02 | within noise | 1.02 |
| F 2D bspline2 dp N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 elastic N=100000 | 1 | 1.02 | faster | 1.02 |
| F 2D bspline2 sanisand N=100000 | 1 | 1.06 | within noise | 1.06 |
| F 2D gimp dp N=100000 | 1 | 1.03 | faster | 1.03 |
| F 2D gimp elastic N=100000 | 1 | 1.04 | faster | 1.04 |

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

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('3d2cb6330be3206ae9ab5e3c29b802f9b98a1e96', '3d2cb6330be3206ae9ab5e3c29b802f9b98a1e96')]
