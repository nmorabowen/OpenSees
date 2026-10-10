---
wp: PR-62
title: "62 -- 1 vanilla row(s)"
pr: "#62"
files: ["`SRC/recorder/TclRecorderCommands.cpp`"]
table: "main"
legacy_seq: [19]
---
| `SRC/recorder/TclRecorderCommands.cpp` | `// Ladruno`: classic-Tcl wiring — add the `recorder ladruno` branch (extern `OPS_LadrunoRecorder` + a `strcmp(argv[1],"ladruno")` block mirroring the `mpco` one), so `recorder ladruno` also works in the classic `OpenSees.exe` Tcl shell (`commands.cpp`→`TclAddRecorder`), not only the Python/interpreter path. Additive. | [#62](https://github.com/nmorabowen/OpenSees/pull/62) |
