---
wp: PR-584
title: "584 -- upstreamable-table row(s)"
pr: "#584"
files: ["`SRC/domain/domain/Domain.cpp`"]
table: "upstreamable"
legacy_seq: [323]
---
| `SRC/domain/domain/Domain.cpp` | `// Ladruno` (ADR-69/ADR-72 P3): `clearAll()` now calls `Ladruno::EnergyChannelRegistry::instance().resetOnWipe()` (+include). The ADR-69 energy-channel registry kept PROCESS-lifetime totals; a diverged explicit run with LNVD active published non-finite work and poisoned RES for every later model in the process (NaN−NaN survives the recorder's baseline subtraction — Linux Zone-A, PR #584 r1), and a huge-but-finite total (~1e300 pre-overflow) ABSORBS a later model's small increments (`total + dE == total` in double precision), silently zeroing its channel delta. wipe() destroys every producer/consumer, so `clearAll` is the semantic zero point. Strictly additive; pairs with the `addEnergy` finiteness guard in fork-owned `LadrunoEnergyChannels.h`. Regression: `tests/test_ladrunoBrick20_dynamics.py::test_energy_registry_survives_prior_diverged_lnvd_run`. | [#584](https://github.com/nmorabowen/OpenSees/pull/584) |
