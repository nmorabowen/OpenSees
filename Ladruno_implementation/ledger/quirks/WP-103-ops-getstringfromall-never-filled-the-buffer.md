---
wp: WP-103
title: "OPS_GetStringFromAll never filled the buffer under classic Tcl (WP-103, #840)"
legacy_seq: 430
---
## `OPS_GetStringFromAll` never filled the buffer under classic Tcl (WP-103, #840)

**Until WP-103, `char tok[64]; OPS_GetStringFromAll(tok, sizeof(tok));` read
UNINITIALISED MEMORY on `OpenSees.exe`.** openseespy was fine, which is why no
deck and no test ever noticed.

The classic-Tcl definition (`SRC/api/elementAPI_TCL.cpp:493`) was:

```cpp
extern "C" const char* OPS_GetStringFromAll(char *buffer, int len)
{ return OPS_GetString(); }          // `buffer` is never touched
```

That is a correct RETURN value — everything really is a string in Tcl — but
`elementAPI.h:209` declares the function with the comment **`// does a strcpy`**,
and the openseespy backend (`PythonModule::getStringFromAll`) does exactly that:
it stringifies ints/floats/str into `buffer` and returns `buffer`. So the two
interpreters disagreed on half the contract.

- **The idiom that breaks.** The whole Ladruno tie/coupling family parses
  `-k <number>|auto` and greedy integer lists as

  ```cpp
  char kTok[64];
  OPS_GetStringFromAll(kTok, sizeof(kTok));
  if (strcmp(kTok, "auto") == 0) ...
  ```

  Under Tcl `kTok` is whatever was on the stack. Measured on the `9c2f964`
  release build, the garbage is even printable:

  ```
  WARNING LadrunoKinematicCoupling: -k wants a number or 'auto', got '<garbage bytes>'
  WARNING LadrunoEmbeddedNode: nHost must be >= 1 (or use -host eleTag); got '<garbage bytes>'
  ```

  `-k`, `-kr`, `-kt`, `-dof` and — worse — `LadrunoEmbeddedNode`/
  `LadrunoEmbeddedRebar`'s **host spec** all went through it. The explicit
  `<nHost> h1..hN` form of those two elements was therefore **unusable from a
  `.tcl` deck at all**, while the identical openseespy call worked.

- **Why it hid for so long.** The failure is classic-Tcl-ONLY, and the fork's
  batteries are openseespy. 111/111 openseespy cases for the four affected
  elements passed on the *pre-fix* build. The only Tcl deck in the repo that
  builds these elements (`Ladruno_scripts/verify_classic_tcl_parity.tcl`)
  passes **no flags**, so it never touched the code path. Two files had already
  hit this and worked around it *privately* — `upGetTok` in `OPS_LadrunoUP.cpp`
  and `ovGetTok` in `OPS_LadrunoPorousOverlay.cpp`, both with a comment saying
  the classic-Tcl version "NEVER touches `buf`" — without anyone fixing the
  function itself.

- **The rule.** `OPS_GetStringFromAll` is the *only* way to read a token that
  may arrive as a typed Python int/float (openseespy's `OPS_GetString` answers
  `"Invalid String Input!"` for those). Since WP-103 both backends fill the
  buffer AND return it, so **either** is now safe. Prefer the return value: it
  is the one form that was always correct on both.

- **Out-of-args is still `0` in the Tcl build, and deliberately so.** The
  "Python path" has two layers that disagree: `PythonModule::getStringFromAll`
  returns `0`, and the `OPS_GetStringFromAll` wrapper in
  `OpenSeesCommands.cpp:1213` maps that `0` to the sentinel string
  `"Invalid String Input!"`. Mirroring the *sentinel* in Tcl would have made
  `OPS_GetStringFromAll` disagree with `OPS_GetString` **in the same
  translation unit** (which returns `0`), and would have silently killed the
  null checks in `upGetTok`/`ovGetTok`. So the Tcl build returns `0` and
  additionally sets `buffer[0] = '\0'` — a caller that ignores the return now
  reads `""`, not garbage.

- **`SRC/runtime/parsing/InterpreterAPI.cpp:140` carries the identical defect
  and was NOT fixed.** `SRC/runtime` is absent from `SRC/CMakeLists.txt`'s
  `add_subdirectory` list, so `OPS_Runtime`/`OpenSeesRT` is **never built in
  this fork** — the code is dead and a fix there would be unverifiable. If that
  target is ever enabled, fix it the same way first.

- **Gate:** `tests/test_wp103_getstringfromall_tcl.py` +
  `tests/tcl/wp103_getstringfromall.tcl`. It shells out to `dist/bin/OpenSees.exe`
  because nothing else can see the bug.
