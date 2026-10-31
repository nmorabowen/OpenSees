---
wp: ADR-39
title: "Contact soft penalty under any explicit integrator: Domain::commit passes the step size"
date: 2026-10-31
pr: "#1048"
table: "main"
files: ["`SRC/domain/domain/Domain.cpp`"]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno`: `Domain::commit()` passes `dT` to the contact engine, `theContactDomain->commit(contactAugmenting, dT)`. The engine keeps it as the last committed step, the `-soft` / `-edgeSoft` step size of an integrator that forms the residual at the committed time (stock `CentralDifference`, where `getDT()` is 0 during the step). Without a soft contact the value is unused, so results are byte-identical. | [#1048](https://github.com/nmorabowen/OpenSees/pull/1048) |
