---
wp: PR-718
title: "ladrunoBuild cmd — machine-readable build-stamp query (git hash the binary was compiled from; the constant CMake stamps into OPENSEES_VERSION, same as the bann…"
pr: "#718"
status: "shipped"
section: "table"
legacy_seq: 131
---
| `ladrunoBuild` cmd — machine-readable build-stamp query (git hash the binary was compiled from; the constant CMake stamps into `OPENSEES_VERSION`, same as the banner's "Ladruno OpenSees build:" line). Born from the TIMs T1 incident (2026-08-10): banner-scraping is the only provenance path and it needs the banner un-suppressed, breaks cp1252 text-mode capture (UTF-8 box glyphs), and is not queryable in-process. `ops.ladrunoBuild() -> str` / Tcl `ladrunoBuild`, registered in BOTH Tcl paths (classic `commands.cpp` + `TclWrapper`) and Python; works under `LADRUNO_OPENSEES_QUIET=1`; plain ASCII. No banner line (info command, not a feature — `LadrunoMassCache` precedent). | runtime info command | — | `OPS_LadrunoBuild()` in `SRC/interpreter/OpenSeesOutputCommands.cpp` + decl `OpenSeesCommands.h`; wrappers `SRC/interpreter/{Python,Tcl}Wrapper.cpp`; classic `SRC/tcl/commands.cpp`; `tests/test_ladruno_build_stamp.py` | shipped | [#718](https://github.com/nmorabowen/OpenSees/pull/718) |
