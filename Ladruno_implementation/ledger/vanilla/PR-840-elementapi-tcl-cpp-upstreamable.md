---
wp: PR-840
title: "840 -- upstreamable-table row(s)"
pr: "#840"
files: ["`SRC/api/elementAPI_TCL.cpp`"]
table: "upstreamable"
legacy_seq: [618]
---
| `SRC/api/elementAPI_TCL.cpp` | `OPS_GetStringFromAll(char* buffer, int len)` ignores `buffer` in the classic-Tcl build while the openseespy build fills it, so the documented (`elementAPI.h:209`) *"does a strcpy"* contract holds on only one of the two interpreters. Vanilla callers that read the buffer -- `SRC/recorder/ElementRecorder.cpp:269`, `SRC/recorder/EnvelopeElementRecorder.cpp:262`, `SRC/interpreter/OpenSeesOutputCommands.cpp:1886` -- read uninitialised heap under Tcl (all three are unreachable from the Tcl command set today, so the bug is latent upstream rather than active). Fix is four lines and preserves the out-of-args return. | [#840](https://github.com/nmorabowen/OpenSees/pull/840) |
