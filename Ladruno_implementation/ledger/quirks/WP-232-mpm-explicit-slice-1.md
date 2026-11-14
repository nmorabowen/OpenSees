---
wp: WP-232
title: "Classic Tcl: OPS_GetNDM() after OPS_ResetInputNoBuilder dereferences a null model builder"
date: 2026-11-13
---
### Classic Tcl: `OPS_GetNDM()` after `OPS_ResetInputNoBuilder` dereferences a null model builder

- **Bites:** a fork command bridged into classic Tcl with `OPS_ResetInputNoBuilder` (the
  WP-127 replay pattern) that calls `OPS_GetNDM()` / `OPS_GetNDF()` before any `model`
  command has run: `elementAPI_TCL.cpp` answers `theModelBuilder->getNDM()`, and
  `theModelBuilder` is set only by the `TclModelBuilder` constructor (and reset to 0 by its
  destructor), so the exe crashes instead of erroring.
- **Workaround:** the classic bridge checks `commands.cpp`'s `theBuilder == 0` and refuses
  before calling the shared entry (`ladrunoMPM`, WP-232). The Python engine is unaffected.
- Learned 2026-11-13 (WP-232c).
