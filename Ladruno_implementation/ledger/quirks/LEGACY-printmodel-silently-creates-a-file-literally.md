---
wp: LEGACY
title: "printModel silently creates a file literally named Invalid String Input! when handed a non-string where it expects a filename"
legacy_seq: 331
---
### `printModel` silently creates a file literally named `Invalid String Input!` when handed a non-string where it expects a filename
- **Bites:** `ops.printModel('-material', 1)` (or any variant that lands in the filename branch with an int) does not warn — it creates a file called `Invalid String Input!` in the current working directory and writes the model into it. In a test run or a batch job this litters the repo with a junk file whose name is the error message that was never raised.
- **Why:** the `else` branch calls `OPS_GetString()` on the integer; the Python backend consumes the argument and returns the literal string `"Invalid String Input!"` rather than a null, and the caller opens that as the output path. Same root cause as the `getString()`-consumes-then-fails behaviour that fork parsers must rewind around with `OPS_ResetCurrentInputArg(-1)`.
- **Workaround:** always pass an explicit `'-file', <path>` when you want file output, and check for a stray `Invalid String Input!` after a printing test.
- **Learned:** 2026-08-26, same session as the row above.
