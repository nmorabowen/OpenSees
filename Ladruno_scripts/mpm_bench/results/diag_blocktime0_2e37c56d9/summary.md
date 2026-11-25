
### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|
| bspline2 | dp | 100000 (100352) | 0.891 (1%) | 1.18 (1%) | 1.33 | 1122 |
| bspline2 | sanisand | 100000 (100352) | 0.329 (68%)* | 0.372 (7%) | 1.13 (within noise) | 3039 |
| gimp | dp | 100000 (100352) | 0.943 (3%) | 1.23 (2%) | 1.30 | 1060 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2D bspline2 dp N=100000 | 5.8 | 0.2 | 4.8 | 5.7 | 0.5 | 3.1 | 1.8 | 71.8 | 5.8 | 0.0 | 94.0 | 0.3 | - |
| 2D bspline2 sanisand N=100000 | 2.1 | 0.1 | 1.8 | 2.2 | 0.2 | 1.2 | 0.7 | 69.3 | 22.3 | 0.0 | 97.8 | 0.1 | - |
| 2D gimp dp N=100000 | 5.9 | 0.2 | 4.5 | 4.0 | 0.5 | 2.6 | 1.4 | 74.5 | 6.0 | 0.0 | 93.9 | 0.3 | - |

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('2e37c56d9a8feac8926575b7e0ac63f9812a04cf', '2e37c56d9a8feac8926575b7e0ac63f9812a04cf')]

1 row(s) with spread > 10 % owed a 5-repeat re-run: F/2D/bspline2/sanisand/N100000/t1
