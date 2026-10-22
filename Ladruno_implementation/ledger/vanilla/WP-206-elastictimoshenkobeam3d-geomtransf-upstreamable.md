---
wp: WP-206
title: "ElasticTimoshenkoBeam3d honours its geomTransf (PDelta, Corotational, -jntOffset)"
date: 2026-10-22
table: "upstreamable"
files: ["`SRC/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.cpp`", "`SRC/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.h`"]
---
| `SRC/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.cpp` | `// Ladruno WP-206`: the element built its own rotation-only Tgl from the transformation's axes and took its geometric flag only from `-geomNonlinear`, so a `PDelta` transf gave no P-Delta, a `Corotational` transf was linear, and `-jntOffset` rigid arms were silently dropped (`recvSelf` meanwhile re-derived the flag from the transf type, disagreeing with the constructor). Now, unless the legacy combination `-geomNonlinear` + `Linear` transf is requested (kept byte-for-byte, with a warning if it would drop offsets), it works in the 6-dof basic system: Timoshenko basic stiffness EI/(L(1+phi))[[4+phi, 2-phi],[2-phi, 4+phi]], basic P-delta matrix (the deformational part of the stock klgeo) under `-geomNonlinear`, element loads split into basic q0 + CrdTransf reactions p0, and the CrdTransf maps everything to global (offsets, PDelta chord stiffness, corotational). New `basicForce` response. Consistent mass still uses rotation-only axes (mass at the nodes, as ElasticBeam3d). `sendSelf`/`recvSelf` data 19 -> 20 (flag). Upstream bug, candidate for upstream. Test `tests/test_wp206_timoshenko_geomtransf.py`. | WP-206 |
| `SRC/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.h` | `// Ladruno WP-206`: basic-system members (`useBasic`, `geomFlag`, `kb`, `kgb`, `qb`, `q0b`, `p0b`) and helpers. | WP-206 |
