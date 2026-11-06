---
wp: PR-38
title: "38 -- 1 vanilla row(s)"
pr: "#38"
files: ["`SRC/element/elasticBeamColumn/ElasticBeam3d.cpp`"]
table: "main"
legacy_seq: [76]
---
| `SRC/element/elasticBeamColumn/ElasticBeam3d.cpp` | `// Ladruno`: add a `"localAxes"` element response (id 30) returning the 9 packed direction cosines from `theCoordTransf->getLocalAxes`, so the Ladruno recorder can write `MODEL/LOCAL_AXES` instead of a silent identity-quaternion fallback (apeGmsh beam-orientation gap). Additive — no existing response touched; pattern to replicate on the other beams. | [#38](https://github.com/nmorabowen/OpenSees/pull/38) |
