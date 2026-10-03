---
wp: ADR-87
title: "Adding an ADR-87 D2 mutation family: the reporting tables are DUPLICATED, their loop bound was hardcoded, and omitting them makes a real mutant answer none"
legacy_seq: 375
---
## Adding an ADR-87 D2 mutation family: the reporting tables are DUPLICATED, their loop bound was hardcoded, and omitting them makes a real mutant answer `none`

`SRC/Ladruno_mutation.h`'s own "ADDING A FAMILY" checklist used to list four
steps: the `LADRUNO_MUTATE_<FAMILY>` default, the two accessor macros, the
`mutation_gate.py` registration, and the manifest score. **Following exactly
those four steps produces a mutant binary that lies about itself.**

Two sites are missing from that list, and they fail in opposite ways:

- **`_ladruno_valid_families` + the cache `STRINGS` property in the root
  `CMakeLists.txt`** (~lines 247/251). An unlisted family is a `FATAL_ERROR`, so
  the mutant configure dies immediately. Loud, harmless, easy.
- **The `ladrunoMutation` family table — which exists TWICE**, once in
  `SRC/interpreter/OpenSeesOutputCommands.cpp` (`OPS_LadrunoMutation`) and again
  in its classic-Tcl twin in `SRC/tcl/commands.cpp` (~line 11549). This one is
  silent and it is the dangerous one.

**Why silent is dangerous here.** Omit the table entry and everything else works
perfectly: CMake applies `add_compile_definitions`, every TU compiles with the
mutation active, the accessors really do zero the tangent and the resultant —
and then `ladrunoMutation()` reports `none`, because it only knows the families
in its table. `mutation_gate.py --expect` catches it and refuses. **Had it not
asked, the harness would have scored a mutant as a baseline and concluded that
the suite cannot detect deleted physics — the exact inversion of the truth, with
no error anywhere.** Measured on ADR 91 (SHELLMOD): mutant built clean, md5
differed from the baseline, `IDENTITY: none`.

**Both tables also iterated with a hardcoded `for (int i = 0; i < 5; i++)`**, so
even adding the sixth entry would still have reported `none`. Both bounds are now
`sizeof(fams)/sizeof(fams[0])`; a table entry can no longer be silently ignored.
The header's checklist now names every site and says which omissions fail loudly
and which fail silently.

**Corollary for gating a family with no element of its own.** ADR 91's feature is
a *section* decorator, so there is no Element API accessor to instrument. Gate the
section accessors instead (`getStressResultant` / `getSectionTangent` /
`getInitialTangent`) — that is the whole physics surface. Expect parse-time
refusal tests and any pure-numpy companion assertions to SURVIVE `ZERO`: they
never run an analysis, so they cannot notice deleted physics, and that is the
correct answer rather than a coverage gap. Audit each survivor individually
before accepting the score.
