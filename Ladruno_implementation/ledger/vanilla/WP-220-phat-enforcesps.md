---
wp: WP-220
title: "Path-control reference load without the prior residual; element failures under imposed displacements are reported"
date: 2026-11-02
table: "upstreamable"
files: ["SRC/analysis/integrator/ArcLength.cpp", "SRC/analysis/integrator/ArcLength1.cpp", "SRC/analysis/integrator/DisplacementControl.cpp", "SRC/analysis/integrator/DistributedDisplacementControl.cpp", "SRC/analysis/integrator/EQPath.cpp", "SRC/analysis/integrator/HSConstraint.cpp", "SRC/analysis/integrator/MinUnbalDispNorm.cpp", "SRC/analysis/handler/TransformationConstraintHandler.cpp"]
---
| `SRC/analysis/integrator/ArcLength.cpp` | WP-220: `domainChanged` phat = R(lambda+1) - R(lambda), via `LadrunoReferenceLoad.h` | TBD |
| `SRC/analysis/integrator/ArcLength1.cpp` | WP-220: as ArcLength | TBD |
| `SRC/analysis/integrator/DisplacementControl.cpp` | WP-220: as ArcLength | TBD |
| `SRC/analysis/integrator/DistributedDisplacementControl.cpp` | WP-220: as ArcLength | TBD |
| `SRC/analysis/integrator/EQPath.cpp` | WP-220: as ArcLength (its reference load is `q`) | TBD |
| `SRC/analysis/integrator/HSConstraint.cpp` | WP-220: as ArcLength | TBD |
| `SRC/analysis/integrator/MinUnbalDispNorm.cpp` | WP-220: as ArcLength | TBD |
| `SRC/analysis/handler/TransformationConstraintHandler.cpp` | WP-220: `enforceSPs` returns -1 if an element failed to update | TBD |
