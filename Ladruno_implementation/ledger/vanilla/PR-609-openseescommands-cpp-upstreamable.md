---
wp: PR-609
title: "609 -- upstreamable-table row(s)"
pr: "#609"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "upstreamable"
legacy_seq: [350]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` **upstream bug fix**: every SECOND consecutive `eigen` in openseespy with NO analysis object failed with "no EigenSOE has been set". Mechanism: `OpenSeesCommands::eigen` builds an EPHEMERAL `DirectIntegrationAnalysis` per call and deletes it after the solve (no-op dtor — the cached `theEigenSOE` survives); the next call builds a fresh analysis whose internal EigenSOE is null, but a same-classTag cached SOE leaves `eigenSOEUpdated=false`, so `setEigenSOE` is never re-invoked. Fix: attach also when `newanalysis` is true (`(eigenSOEUpdated \|\| newanalysis)`). Classic Tcl re-attaches and never had the bug. Upstream `OpenSees/OpenSees` master carries the same defect — the ephemeral-analysis build, the `newanalysis` flag, and the `eigenSOEUpdated`-gated attach block all match upstream — so it is an **upstream candidate**, but as a **re-port not a cherry-pick**: the fork's surrounding lines carry the Ladruno FEAST/CMS `providedEigenSOE` seam (ADR-43) that upstream lacks, so the patch will not apply verbatim. Regression: `tests/test_eigen_repeat.py` (zone_a). | [#609](https://github.com/nmorabowen/OpenSees/pull/609) |
