
### Mode F, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (10368) | 0.548 (1%) | 1.75 (6%) | 2.88 (8%) | 5.25 | 1825 |
| bspline2 | sanisand | 10000 (10368) | 0.0444 (2%) | 0.139 (9%) | 0.209 (8%) | 4.70 | 22501 |
| gimp | dp | 10000 (10368) | 0.657 (1%) | 2.05 (2%) | 3.48 (11%) | 5.29 | 1522 |
| bspline2 | dp | 100000 (98496) | 0.581 (2%) | 1.78 (2%) | 2.89 (1%) | 4.97 | 1723 |
| bspline2 | sanisand | 100000 (98496) | 0.102 (0%) | 0.288 (8%) | 0.456 (4%) | 4.47 | 9798 |
| gimp | dp | 100000 (98496) | 0.691 (1%) | 2.19 (3%) | 3.65 (5%) | 5.29 | 1448 |

### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=4 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (9800) | 0.893 (1%) | 3.11 (3%) | 5.29 (16%) | 5.93 | 1120 |
| bspline2 | sanisand | 10000 (9800) | 0.207 (0%) | 0.466 (12%) | 0.669 (4%) | 3.23 | 4823 |
| gimp | dp | 10000 (9800) | 0.958 (2%) | 3.01 (8%) | 5.11 (17%) | 5.33 | 1044 |
| bspline2 | dp | 100000 (100352) | 0.908 (1%) | 2.89 (9%) | 4.72 (9%) | 5.20 | 1102 |
| bspline2 | sanisand | 100000 (100352) | 0.344 (3%) | 0.839 (9%) | 0.907 (3%) | 2.63 | 2903 |
| gimp | dp | 100000 (100352) | 0.964 (2%) | 3.06 (3%) | 4.98 (9%) | 5.17 | 1037 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3D bspline2 dp N=10000 | 2.2 | 0.2 | 8.5 | 12.7 | 0.4 | 6.2 | 3.4 | 63.0 | 3.1 | 0.0 | 97.6 | 0.3 | 43.3 |
| 3D bspline2 sanisand N=10000 | 0.2 | 0.0 | 0.7 | 1.0 | 0.0 | 0.5 | 0.3 | 94.3 | 2.9 | 0.0 | 99.8 | 0.0 | 94.6 |
| 3D gimp dp N=10000 | 2.5 | 0.2 | 6.7 | 8.4 | 0.4 | 4.7 | 2.5 | 70.5 | 3.7 | 0.0 | 97.2 | 0.4 | 51.1 |
| 3D bspline2 dp N=100000 | 2.6 | 0.2 | 8.9 | 13.2 | 0.3 | 6.7 | 3.5 | 60.8 | 3.6 | 0.0 | 97.2 | 0.2 | 42.8 |
| 3D bspline2 sanisand N=100000 | 0.5 | 0.0 | 1.5 | 2.3 | 0.1 | 1.2 | 0.6 | 87.3 | 6.4 | 0.0 | 99.5 | 0.0 | 88.9 |
| 3D gimp dp N=100000 | 3.1 | 0.3 | 7.7 | 9.4 | 0.3 | 5.3 | 2.7 | 66.8 | 4.2 | 0.0 | 96.6 | 0.2 | 50.2 |
| 2D bspline2 dp N=10000 | 4.9 | 0.2 | 4.9 | 5.7 | 0.6 | 3.1 | 1.6 | 73.8 | 4.6 | 0.0 | 94.7 | 0.4 | 60.7 |
| 2D bspline2 sanisand N=10000 | 1.2 | 0.0 | 1.1 | 1.3 | 0.1 | 0.7 | 0.4 | 81.9 | 13.1 | 0.0 | 98.7 | 0.1 | 89.3 |
| 2D gimp dp N=10000 | 4.9 | 0.3 | 4.9 | 4.3 | 0.5 | 2.8 | 1.4 | 75.8 | 4.6 | 0.0 | 94.7 | 0.5 | 63.1 |
| 2D bspline2 dp N=100000 | 5.6 | 0.2 | 4.9 | 5.7 | 0.5 | 3.1 | 1.6 | 72.4 | 5.5 | 0.0 | 94.2 | 0.3 | 61.3 |
| 2D bspline2 sanisand N=100000 | 2.1 | 0.1 | 1.9 | 2.2 | 0.2 | 1.2 | 0.6 | 69.5 | 22.2 | 0.0 | 97.8 | 0.1 | 83.2 |
| 2D gimp dp N=100000 | 6.0 | 0.2 | 4.5 | 3.7 | 0.5 | 2.6 | 1.3 | 74.9 | 5.8 | 0.0 | 93.7 | 0.3 | 64.3 |

### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)

| row | threads | T ratio | verdict | G(n) |
|---|---|---|---|---|
| F 3D bspline2 dp N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 sanisand N=10000 | 1 | 1.01 | within noise | 1.01 |
| F 3D gimp dp N=10000 | 1 | 1.03 | faster | 1.03 |
| F 3D bspline2 dp N=100000 | 1 | 1.01 | within noise | 1.01 |
| F 3D bspline2 sanisand N=100000 | 1 | 1.03 | faster | 1.03 |
| F 3D gimp dp N=100000 | 1 | 1.01 | within noise | 1.01 |
| F 2D bspline2 dp N=10000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 sanisand N=10000 | 1 | 1.02 | faster | 1.02 |
| F 2D gimp dp N=10000 | 1 | 1.03 | within noise | 1.03 |
| F 2D bspline2 dp N=100000 | 1 | 0.99 | within noise | 0.99 |
| F 2D bspline2 sanisand N=100000 | 1 | 1.05 | within noise | 1.05 |
| F 2D gimp dp N=100000 | 1 | 1.02 | faster | 1.02 |

### Determinism (particle-array digests against 1 thread)

| row | threads | identical to t=1 |
|---|---|---|
| F 2D bspline2 dp N=10000 | 1 | - |
| F 2D bspline2 dp N=10000 | 4 | yes |
| F 2D bspline2 dp N=10000 | 8 | yes |
| F 2D bspline2 dp N=100000 | 1 | - |
| F 2D bspline2 dp N=100000 | 4 | yes |
| F 2D bspline2 dp N=100000 | 8 | yes |
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
| F 3D bspline2 dp N=10000 | 1 | - |
| F 3D bspline2 dp N=10000 | 4 | yes |
| F 3D bspline2 dp N=10000 | 8 | yes |
| F 3D bspline2 dp N=100000 | 1 | - |
| F 3D bspline2 dp N=100000 | 4 | yes |
| F 3D bspline2 dp N=100000 | 8 | yes |
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

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('efc4113f1e116e459df26be7efeaa675257882f9', 'efc4113f1e116e459df26be7efeaa675257882f9')]
