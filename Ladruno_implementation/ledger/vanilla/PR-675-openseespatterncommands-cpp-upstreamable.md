---
wp: PR-675
title: "675 -- upstreamable-table row(s)"
pr: "#675"
files: ["`SRC/interpreter/OpenSeesPatternCommands.cpp`"]
table: "upstreamable"
legacy_seq: [401]
---
| `SRC/interpreter/OpenSeesPatternCommands.cpp` | `// Ladruno`: **`sp ... -subtractInit` was a NO-OP on the openseespy path.** `OPS_SP()` handled the flag with `retZeroInitValue = true` — already the initialiser at `:1123` — so the branch did nothing and incremental SP was unreachable from Python, while **both** Tcl parsers set it to `false` correctly (`SRC/modelbuilder/tcl/TclModelBuilder.cpp:3808-3810`, `SRC/runtime/commands/modeling/constraint.cpp:392-394`). A Python-vs-Tcl behaviour split, not a missing feature. Consequence: a staged deck (gravity, then an imposed displacement) silently **yanked the node to the absolute value** instead of moving it incrementally from its staged state — wrong answer, no warning. Measured A/B on the same deck: `-subtractInit` gave **3.0 before / 5.5 after**, where Tcl already gave 5.5 on the unfixed binary. One-token change (`true` → `false`) plus the comment that records why. Only handlers reading `SP_Constraint::getInitialValue()` honour the flag (`TransformationDOF_Group:1068`, `PenaltySP_FE:125`, `LagrangeSP_FE:143`); under `constraints Plain` it is inert by design, fixed or not. Gated by `tests/test_sp_subtract_init.py` (verified to FAIL on the pre-fix binary). | [#675](https://github.com/nmorabowen/OpenSees/pull/675) |
