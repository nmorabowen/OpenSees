
### Mode A, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=8 | S(8) | ns/p-step (t=1) | k (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (10648) | 0.000329 (0%) | 0.000248 (41%)* | 0.75 (within noise) | 3043232 | 14.3 |
| bspline2 | elastic | 10000 (10648) | 0.00131 (0%) | 0.00137 (0%) | 1.05 | 764947 | 5.8 |
| gimp | dp | 10000 (10648) | 0.00256 (0%) | 0.00301 (5%) | 1.18 | 391377 | 9.2 |
| gimp | elastic | 10000 (10648) | 0.00632 (1%) | 0.00681 (1%) | 1.08 | 158238 | 6.0 |
| gimp | sanisand | 10000 (10648) | 0.0124 (1%) | - | - | 80388 | 2.0 |
| gimp | dp | 100000 (-) | fail | fail | - | - | 0.0 |
| gimp | elastic | 100000 (97336) | 0.00205 (0%) | 0.00322 (1%) | 1.57 | 487446 | 6.3 |

### Mode A, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=8 | S(8) | ns/p-step (t=1) | k (t=1) |
|---|---|---|---|---|---|---|---|
| bspline2 | dp | 10000 (-) | fail | fail | - | - | 0.0 |
| bspline2 | elastic | 10000 (10000) | 0.02 (0%) | 0.022 (1%) | 1.10 | 50025 | 6.2 |
| gimp | dp | 10000 (-) | fail | fail | - | - | 0.0 |
| gimp | elastic | 10000 (10000) | 0.033 (0%) | 0.0354 (1%) | 1.07 | 30294 | 6.2 |

### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)

| row | threads | T ratio | verdict | G(n) |
|---|---|---|---|---|
| A 3D bspline2 dp N=10000 | 1 | 1.02 | within noise | 1.02 |
| A 3D bspline2 elastic N=10000 | 1 | 1.02 | within noise | 1.02 |
| A 3D bspline2 elastic N=10000 | 8 | 1.02 | within noise | 1.07 |
| A 3D gimp dp N=10000 | 1 | 1.02 | faster | 1.02 |
| A 3D gimp dp N=10000 | 8 | 1.15 | faster | 1.21 |
| A 3D gimp elastic N=10000 | 1 | 1.02 | faster | 1.02 |
| A 3D gimp elastic N=10000 | 8 | 1.09 | faster | 1.10 |

### Determinism (particle-array digests against 1 thread)

| row | threads | identical to t=1 |
|---|---|---|
| A 2D bspline2 elastic N=10000 | 1 | - |
| A 2D bspline2 elastic N=10000 | 8 | NO |
| A 2D gimp elastic N=10000 | 1 | - |
| A 2D gimp elastic N=10000 | 8 | NO |
| A 3D bspline2 dp N=10000 | 1 | - |
| A 3D bspline2 dp N=10000 | 8 | NO |
| A 3D bspline2 elastic N=10000 | 1 | - |
| A 3D bspline2 elastic N=10000 | 8 | NO |
| A 3D gimp dp N=10000 | 1 | - |
| A 3D gimp dp N=10000 | 8 | NO |
| A 3D gimp elastic N=10000 | 1 | - |
| A 3D gimp elastic N=10000 | 8 | NO |
| A 3D gimp elastic N=100000 | 1 | - |
| A 3D gimp elastic N=100000 | 8 | NO |

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('3d2cb6330be3206ae9ab5e3c29b802f9b98a1e96', '3d2cb6330be3206ae9ab5e3c29b802f9b98a1e96')]

1 row(s) with spread > 10 % owed a 5-repeat re-run: A/3D/bspline2/dp/N10000/t8

Failures (data, not dropped): A/2D/gimp/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/gimp/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/gimp/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/gimp/dp/N10000/t8: OpenSeesError('See stderr output'); A/2D/gimp/dp/N10000/t8: OpenSeesError('See stderr output'); A/2D/gimp/dp/N10000/t8: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t1: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t8: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t8: OpenSeesError('See stderr output'); A/2D/bspline2/dp/N10000/t8: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t1: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t1: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t1: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t8: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t8: OpenSeesError('See stderr output'); A/3D/gimp/dp/N100000/t8: OpenSeesError('See stderr output')
