---
wp: PR-539
title: "539 -- 2 vanilla row(s)"
pr: "#539"
files: ["`SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}`", "`SRC/analysis/analysis/CMakeLists.txt`"]
table: "main"
legacy_seq: [251, 252]
---
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.{h,cpp}` | `// Ladruno` ADR44 P1b: additive opt-in `-combine {SRSS\|CQC\|ABS\|TenPercent}` (+ `-damp`/`-modalDamp` for CQC ζ) stage on the interpreter `responseSpectrumAnalysis` — parses the flags, and when present runs `analyzeCombined()` (combines the per-mode peak nodal displacements via `LadrunoModalCombination::combine` and commits ONE field) instead of the per-mode states. `-combine` absent ⇒ byte-identical to stock (Petracca per-mode path untouched). In-class member init (`m_combine=-1`, `m_xi`) keeps the ctor untouched; `#include <LadrunoModalCombination.h>`. | [#539](https://github.com/nmorabowen/OpenSees/pull/539) |
| `SRC/analysis/analysis/CMakeLists.txt` | ADR44 P1b: add `LadrunoModalCombination.h` (header-only combination kernel) to `OPS_Analysis` PUBLIC headers. | [#539](https://github.com/nmorabowen/OpenSees/pull/539) |
