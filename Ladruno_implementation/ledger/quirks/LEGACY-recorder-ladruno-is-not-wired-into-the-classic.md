---
wp: LEGACY
title: "recorder ladruno is NOT wired into the classic Tcl OpenSees.exe (was; now fixed)"
legacy_seq: 43
---
### `recorder ladruno` is NOT wired into the classic Tcl `OpenSees.exe` (was; now fixed)
- **Bites:** `recorder ladruno ...` works from OpenSeesPy/openseesmp and the
  interpreter-based Tcl (`TclWrapper`→`OPS_Recorder`, the shared map in
  `OpenSeesOutputCommands.cpp`), but the **classic** Tcl `OpenSees.exe`
  (`commands.cpp`→`addRecorder`→`TclAddRecorder` in `TclRecorderCommands.cpp`)
  hardcodes its recorder dispatch in a *separate* file that the rename PR never
  touched — it had `mpco`/`vtkhdf`/`gmsh`/`EnergyBalance` but no `ladruno`. So
  `recorder ladruno` raised "recorder type ladruno is unknown" only in `OpenSees.exe`.
- **Fix:** added the `else if (strcmp(argv[1],"ladruno")==0)` branch (+ extern
  `OPS_LadrunoRecorder`) mirroring the `mpco` block in `TclRecorderCommands.cpp`.
  Lesson: there are TWO recorder-dispatch tables (shared `OPS_Recorder` map for
  Py/interpreter-Tcl; hardcoded `TclAddRecorder` for classic Tcl) — wire new
  recorders into BOTH. Learned 2026-05-31.
