---
wp: LEGACY
title: "A failed OPS_GetIntInput/OPS_GetDoubleInput consumes the arg in openseespy but NOT in the classic Tcl exe — blind OPS_ResetCurrentInputArg(-1) overshoots in Tcl"
legacy_seq: 11
---
### A failed `OPS_GetIntInput`/`OPS_GetDoubleInput` consumes the arg in openseespy but NOT in the classic Tcl exe — blind `OPS_ResetCurrentInputArg(-1)` overshoots in Tcl
- **Bites:** any greedy "read ints until the next flag" parse scan that un-gets the
  flag with `OPS_ResetCurrentInputArg(-1)` after a failed numeric read. Symptom:
  the command parses fine from openseespy but is a parse error from the classic
  Tcl exe **only when the greedy option is followed by another flag** — e.g.
  `recorder ladruno f -G energy 1 -T nsteps 1` died with *"option -T nsteps
  requires a number-of-steps argument"* while `… -T nsteps 1 -G energy 1` ran.
- **Why:** the two API stacks disagree on cursor consumption at parse failure.
  openseespy (`PythonModule.cpp::getInt`) calls `incrCurrentArg()` *before* the
  type check, so the failed token is consumed and `-1` restores it — net zero.
  The classic Tcl exe (`SRC/api/elementAPI_TCL.cpp::OPS_GetIntInput`, ~line 224)
  returns `-1` *without* advancing, so the `-1` rewind steps back onto the
  previous already-consumed token and the whole tail parse shifts by one
  (a local `numdata` countdown then starves the line's last token). String reads
  (`OPS_GetString`/`OPS_GetStringFromAll`) consume in both stacks — only the
  numeric readers diverge.
- **Workaround/status (2026-06-11):** never blind-rewind — capture
  `int oldn = OPS_GetNumRemainingInputArgs();` before the numeric read and rewind
  only `if (OPS_GetNumRemainingInputArgs() < oldn)`. Applied to the
  LadrunoRecorder `-G` scan, both LadrunoMonitorRecorder scans, and the
  LadrunoBrick `-hourglass` coeff probe (LadrunoRCConcrete's `readList` already
  did it right). Guarded by the `TCL FLAG ORDER` regression gate
  (`flag_order_model.tcl`); the openseespy ordering was already covered by
  `energy_model.py`.
