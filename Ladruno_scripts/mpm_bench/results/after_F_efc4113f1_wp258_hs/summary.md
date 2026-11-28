
### Mode F, 3D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|
| bspline2 | hs | 100000 (98496) | 0.599 (1%) | 2.71 (5%) | 4.52 | 1669 |

### Mode F, 2D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, rerun owed

| basis | material | N (actual) | t=1 | t=8 | S(8) | ns/p-step (t=1) |
|---|---|---|---|---|---|---|
| bspline2 | hs | 100000 (100352) | 1.07 (0%) | 5.09 (8%) | 4.74 | 931 |

### Mode F phase split at 1 thread: share of the `step` wall [%]

| row | mpm.f.dtcrit | mpm.f.admissible | mpm.f.support | mpm.f.p2g | mpm.f.grid | mpm.f.g2p | mpm.f.project2 | mpm.f.update | mpm.f.commit | commit | mpm.advance | other | material share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3D bspline2 hs N=100000 | 12.0 | 0.2 | 9.1 | 13.7 | 0.3 | 6.9 | 3.6 | 48.2 | 5.8 | 0.0 | 87.8 | 0.2 | 33.7 |
| 2D bspline2 hs N=100000 | 13.6 | 0.2 | 5.8 | 6.7 | 0.6 | 3.7 | 1.9 | 56.7 | 10.4 | 0.0 | 86.1 | 0.4 | 49.6 |

### Determinism (particle-array digests against 1 thread)

| row | threads | identical to t=1 |
|---|---|---|
| F 2D bspline2 hs N=100000 | 1 | - |
| F 2D bspline2 hs N=100000 | 8 | yes |
| F 3D bspline2 hs N=100000 | 1 | - |
| F 3D bspline2 hs N=100000 | 8 | yes |

Provenance: 13th Gen Intel(R) Core(TM) i7-13700KF | 24 logical | Windows-11-10.0.26100-SP0 | MKL mkl_core.3.dll 2026.1.0.0 | builds [('efc4113f1e116e459df26be7efeaa675257882f9', 'efc4113f1e116e459df26be7efeaa675257882f9')]
