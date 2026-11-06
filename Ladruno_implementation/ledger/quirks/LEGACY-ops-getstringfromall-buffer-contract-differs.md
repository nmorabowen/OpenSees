---
wp: LEGACY
title: "OPS_GetStringFromAll buffer contract differs between interpreters (classic-Tcl reads garbage)"
legacy_seq: 169
---
## `OPS_GetStringFromAll` buffer contract differs between interpreters (classic-Tcl reads garbage)

**Symptom** (ADR-71 P1, 2026-07-10): `element LadrunoUP ...` in classic Tcl fails
with "got 0 leading integers" while the identical command works in openseespy.

**Cause**: the two backends implement `OPS_GetStringFromAll(buf, len)` with
DIFFERENT buffer semantics. The modern interpreter
(`OpenSeesCommands.cpp:1185`) copies the token into `buf` and returns it; the
classic-Tcl shim (`elementAPI_TCL.cpp:494`) just does `return OPS_GetString();`
— **it never touches `buf`**. Any parser that calls the function and then reads
`buf` (a common family idiom) parses garbage under classic Tcl and silently
"works" everywhere else.

**Rule**: always consume the RETURN VALUE, never the buffer. If you need the
token in a local buffer, normalize:
`const char* s = OPS_GetStringFromAll(buf,len); if (s && s != buf) strncpy(...)`
(see `upGetTok()` in `SRC/element/ladrunoUP/OPS_LadrunoUP.cpp`). Fix the vanilla
shim only via an ADR — several upstream parsers may depend on the current
behavior.
