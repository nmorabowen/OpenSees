---
wp: LEGACY
title: "printModel('-ele', <tag>) fails from openseespy — a Tcl-shaped OPS_ResetCurrentInputArg(2) rewind that assumes argv[0] is the command name"
legacy_seq: 330
---
### `printModel('-ele', <tag>)` fails from openseespy — a Tcl-shaped `OPS_ResetCurrentInputArg(2)` rewind that assumes `argv[0]` is the command name
- **Bites:** `ops.printModel('-ele', 1)` aborts with `WARNING print ele failed to get integer:` and prints nothing, while `ops.printModel('-ele')` (all elements) works fine. Reads like a bad tag.
- **Why:** `printElement()` (`SRC/interpreter/OpenSeesOutputCommands.cpp:2012`) rewinds the argument cursor with `OPS_ResetCurrentInputArg(2)`. That constant is correct for the classic-Tcl `argv`, where index 0 is the command name `print` and index 1 is `-ele`. In the Python backend arg 0 is already `-ele`, so rewinding by 2 lands *past* the tag and the subsequent int read consumes the wrong token. Upstream, not Ladruno.
- **Workaround:** print all elements (`printModel('-ele')`), or print to a file (`printModel('-file', path, '-ele')`) and filter. Both reach `NDMaterial::Print` for the Gauss-point material via `Brick::Print` → `materialPointers[0]->Print`.
- **Learned:** 2026-08-26, writing the `Print`-override fingerprint for [[86_ladruno_sanisand_adr]] PR-1.
