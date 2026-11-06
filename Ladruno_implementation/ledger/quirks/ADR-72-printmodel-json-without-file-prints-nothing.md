---
wp: ADR-72
title: "printModel('-JSON') without -file prints NOTHING (upstream interpreter quirk, ADR 72 P2)"
legacy_seq: 180
---
## `printModel('-JSON')` without `-file` prints NOTHING (upstream interpreter quirk, ADR 72 P2)

`OPS_printModel` (SRC/interpreter/OpenSeesCommands.cpp:2927-2930) treats
`-JSON` as flag-only: the argument loop sets `flag = OPS_PRINT_PRINTMODEL_JSON`
and moves on, but `theDomain->Print(...)` is reached ONLY inside the
filename branch (`print <file>` / `-file <file>`). Bare `ops.printModel()`
prints the whole domain (classic format, opserr); bare
`ops.printModel("-JSON")` silently emits nothing at all — a capfd-based
"assert JSON contains X" test then fails with empty capture (bit the ADR-72
P2 S8 acceptance test). **Rule: to assert on model JSON, always use
`ops.printModel("-JSON", "-file", tmpfile)` and read the file** (OVERWRITE
mode, so the file is complete + fresh). Classic-format assertions via capfd
remain fine. Upstream behavior — fix-on-touch only.
