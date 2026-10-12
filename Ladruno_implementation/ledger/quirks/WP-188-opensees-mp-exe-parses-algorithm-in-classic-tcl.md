---
wp: WP-188
title: "`algorithm KrylovNewton` builds AcceleratedNewton in every binary, and OpenSeesMP.exe parses it in SRC/tcl/commands.cpp"
date: 2026-10-11
---
### `algorithm KrylovNewton` builds `AcceleratedNewton` in every binary, and `OpenSeesMP.exe` parses it in `SRC/tcl/commands.cpp`

- **What bites.** A change made to `KrylovNewton.cpp` has no effect on `algorithm KrylovNewton`:
  both the classic Tcl parser (`SRC/tcl/commands.cpp`) and the interpreter parser
  (`OPS_KrylovNewton`, `SRC/interpreter/OpenSeesCommands.cpp`) construct
  `AcceleratedNewton` + `KrylovAccelerator`. The PROFILE_4D prototype lost a probe (N3) to
  this. Second trap: `OpenSeesMP.exe` (and `OpenSees.exe`, `OpenSeesSP.exe`) compile
  `SRC/tcl/commands.cpp`; only openseespy / openseesmp use `OPS_KrylovNewton`. An algorithm
  option added to one parser only is silently ignored by the other binaries (the vanilla
  KrylovNewton loops skip unknown flags).
- **Workaround.** Patch `AcceleratedNewton`, and both parsers (WP-188 does; its
  `test_classic_parser_fails_closed` and `test_openseespy_parser_and_linear_reuse` cover
  one each). Check with `grep -o "OpenSeesMP.dir[^ ]*commands.cpp.obj" build.ninja`.
- **Related.** OpenSeesMP + MUMPS is not byte-reproducible run to run even on a 300-DOF frame
  at 2 ranks (three reruns of one binary, three md5s, 1e-16 relative): an MP output-md5 gate
  cannot work; compare at round-off or gate on the serial build (WP-188 §6.1).
- Learned 2026-10-11 (WP-188, [#PRNUM](https://github.com/nmorabowen/OpenSees/pull/PRNUM)).
