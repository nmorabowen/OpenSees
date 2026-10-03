---
wp: ADR-74
title: "ADR-74 N1 PR -- upstreamable-table row(s)"
files: ["`SRC/analysis/numberer/ParallelNumberer.{h,cpp}`", "`SRC/classTags.h`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.{h,cpp}`", "`SRC/analysis/numberer/CMakeLists.txt`"]
table: "upstreamable"
legacy_seq: [325, 326, 327, 328, 329]
---
| `SRC/analysis/numberer/ParallelNumberer.{h,cpp}` | `// Ladruno` (ADR-74 N1, ledgered promotion): private→protected for `theNumberer`/`processID`/`numChannels`/`theChannels` + two tag-forwarding protected ctors (`(int classTag, GraphNumberer&)`, `(int classTag)`), HHT/GeneralizedAlpha pattern. Visibility/ctor-surface only — zero behavior change; enables `LadrunoParallelNumberer`=33000 | ADR-74 N1 PR |
| `SRC/classTags.h` | `// Ladruno` (ADR-74 N1): register `NUMBERER_TAG_LadrunoParallelNumberer`=33000 (numberer registry band; numerically-equal tags across registries are blessed) | ADR-74 N1 PR |
| `SRC/tcl/commands.cpp` | `// Ladruno` (ADR-74 N1): `LadrunoParallelRCM`/`LadrunoParallelPlain` verbs in the `_PARALLEL_INTERPRETERS` numberer chain (mirror the stock wiring: `setProcessID(OPS_rank)` + `setChannels`) + include | ADR-74 N1 PR |
| `SRC/interpreter/OpenSeesCommands.{h,cpp}` | `// Ladruno` (ADR-74 N1): `OPS_LadrunoParallelRCM`/`OPS_LadrunoParallelPlain` factory twins of `OPS_ParallelRCM` + dispatch entries + decls + include | ADR-74 N1 PR |
| `SRC/analysis/numberer/CMakeLists.txt` | Add `LadrunoParallelNumberer.cpp` to `OPS_Analysis` sources + header | ADR-74 N1 PR |
