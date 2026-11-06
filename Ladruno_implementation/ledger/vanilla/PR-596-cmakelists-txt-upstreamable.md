---
wp: PR-596
title: "596 -- upstreamable-table row(s)"
pr: "#596"
files: ["`CMakeLists.txt` (root)", "`SRC/system_of_eqn/CMakeLists.txt`", "`SRC/classTags.h`", "`SRC/analysis/analysis/{StaticAnalysis,DirectIntegrationAnalysis}.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`Ladruno_scripts/stamp_headers.py`"]
table: "upstreamable"
legacy_seq: [334, 335, 336, 337, 338, 339, 340]
---
| `CMakeLists.txt` (root) | `# Ladruno ADR 1000`: add the OFF-by-default `LADRUNO_CMS` option, isolated `OPS_LadrunoCMS` target and MP/PyMP-only linkage. `LADRUNO_CMS_BUILD_TESTS` is also OFF by default and only materializes standalone numerical checks. OFF does not add CMS sources or definitions. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `SRC/system_of_eqn/CMakeLists.txt` | `# Ladruno ADR 1000`: enter the new `ladrunoCMS` subdirectory only when `LADRUNO_CMS=ON`. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `SRC/classTags.h` | `// Ladruno ADR 1000`: reserve/register the independent pair `EigenSOE_TAGS_LadrunoCMS=33025` and `EigenSOLVER_TAGS_LadrunoCMS=33026`. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `SRC/analysis/analysis/{StaticAnalysis,DirectIntegrationAnalysis}.cpp` | `// Ladruno ADR 1000`: `setEigenSOE` compares `EigenSOE` **object identity** (was: classTag) before deleting the analysis-owned SOE. UNGATED but behavior-neutral for stock flows — merge review traced every in-tree caller: (i) different-classTag swap deletes in both semantics; (ii) same-object re-attach is guarded by `!= &theNewSOE`; (iii) repeated same-type `eigen` never reaches `setEigenSOE` (both interpreters reuse the cached SOE). It fixes the documented FEAST-era trap where a fresh same-classTag SOE was silently dropped (stale solve + leak). Empirical: classic-Tcl repeated/switched/transient-path eigen bit-identical; the OPS-path second-`eigen` failure is pre-existing on `ladruno` (identical with/without this edit). | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `SRC/tcl/commands.cpp` | `// Ladruno ADR 1000`: add the opt-in Tcl `eigen -ladrunoCMS` parser and construct the independent SOE/solver pair only in CMS-enabled MP builds. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno ADR 1000`: mirror the same `eigen -ladrunoCMS` option contract in the Python/shared interpreter path. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
| `Ladruno_scripts/stamp_headers.py` | `# Ladruno ADR 1000`: include the new `SRC/system_of_eqn/ladrunoCMS` directory in the fork header-stamping inventory. | [#596](https://github.com/nmorabowen/OpenSees/pull/596) |
