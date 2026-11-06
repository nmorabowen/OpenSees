---
wp: PR-346
title: "346 -- 2 vanilla row(s)"
pr: "#346"
files: ["`SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`", "`SRC/analysis/fe_ele/FE_Element.cpp`"]
table: "main"
legacy_seq: [52, 68]
---
| `SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-39 P2a: add the `contactPlane tag slaveSurf nx ny nz px py pz kn` command (rigid analytical plane) — OPS_ body + decl + Py + Tcl wrappers | [#346](https://github.com/nmorabowen/OpenSees/pull/346) |
| `SRC/analysis/fe_ele/FE_Element.cpp` | `// Ladruno` ADR-39 P2a: in the subtype ctor `FE_Element(tag, numDOF_Group, ndof)`, move `numFEs++` BELOW the `if (numFEs == 0)` class-wide-scratch allocation guard (was incremented first ⇒ guard dead). A model whose ONLY FE_Elements are subtype adapters (a contact-only model with no Domain Elements) never allocated `theMatrices`/`theVectors` ⇒ null-deref in `~FE_Element` when the last FE was destroyed (teardown segfault). Zero behaviour change for any model with an element-backed FE (numFEs already >0 there); fixes all FE subtypes (PenaltySP_FE etc.), not just contact. Found by the P2a code-gate (B1 BLOCKER). | [#346](https://github.com/nmorabowen/OpenSees/pull/346) |
