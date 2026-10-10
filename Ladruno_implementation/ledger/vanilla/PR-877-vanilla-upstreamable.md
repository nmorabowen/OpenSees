---
wp: PR-877
title: "PR-877 vanilla rows (upstreamable)"
pr: "#877"
date: 2026-10-10
files: ["`SRC/material/nD/NDMaterial.cpp`", "`SRC/interpreter/OpenSeesNDTestCommands.cpp`"]
table: "upstreamable"
---
| `SRC/material/nD/NDMaterial.cpp` | `// Ladruno WP-99` (#877): `OPS_clearAllNDMaterial()` — the nD-material wipe hook (the same hook as the #841 `ladrunoSanisandResetImplexGlobals()` call, a separate addition) — also calls `ladrunoConcrete3dResetWarningBudget()` (fork-owned `LadrunoConcrete3D.cpp`, local extern; resets the process-wide return-map warning budget) and `ladrunoClearCommitRefusals()` (header-only, `LadrunoMaterialStatus.h`) so a refusal declared by an NDTest/probe commit, or by a model that died mid-commit, cannot fail the NEXT model's first commit. One call each per wipe; stock decks byte-identical. | [#877](https://github.com/nmorabowen/OpenSees/pull/877) |
| `SRC/interpreter/OpenSeesNDTestCommands.cpp` | `// Ladruno WP-99`: `OPS_NDCommitState` commits the material inside a `LadrunoProbeCommitScope` - NDTest commits OUTSIDE `Domain::commit()`, the only place the commit-refusal counter is checked and cleared, so a refusal declared here would otherwise leak into the next real commit. | [#877](https://github.com/nmorabowen/OpenSees/pull/877) |
