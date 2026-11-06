---
wp: PR-39
title: "39 -- 9 vanilla row(s)"
pr: "#39"
files: ["`SRC/element/elasticBeamColumn/ElasticBeam2d.cpp`", "`SRC/element/dispBeamColumn/DispBeamColumn2d.cpp`", "`SRC/element/dispBeamColumn/DispBeamColumn3d.cpp`", "`SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp`", "`SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp`", "`SRC/element/mixedBeamColumn/MixedBeamColumn3d.cpp`", "`SRC/element/mixedBeamColumn/MixedBeamColumn2d.cpp`", "`SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.cpp`", "`SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn2d.cpp`"]
table: "main"
legacy_seq: [77, 78, 79, 80, 81, 82, 83, 84, 85]
---
| `SRC/element/elasticBeamColumn/ElasticBeam2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) → `Vector(9)` dir cosines from `theCoordTransf->getLocalAxes` (2D transf returns full 3D cosines, vz=(0,0,1)). Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/mixedBeamColumn/MixedBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive — completes the remaining-beam localAxes coverage. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/mixedBeamColumn/MixedBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes` (added to the `switch`-based `getResponse` as `case 30`). Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes` (`case 30`). Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
