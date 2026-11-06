---
wp: LEGACY
title: "OPS_GetString returns the literal \"Invalid String Input!\" for a Python NUMERIC arg — never use it to peek an optional numeric token"
legacy_seq: 19
---
### `OPS_GetString` returns the literal `"Invalid String Input!"` for a Python NUMERIC arg — never use it to peek an optional numeric token
- **Bites:** a factory that peeks an optional trailing numeric arg by calling
  `OPS_GetString()` and parsing it (e.g. `strtod`) works under Tcl (args arrive as
  strings, `"0.05"` parses) but **silently drops the value under openseespy**. For a
  Python number (`ops.element(..., '-hourglass','viscous', 0.05)`), `0.05` is a
  `PyFloat`, not a `PyUnicode`. `PythonModule::getString()` (`PyUnicode_Check` fails)
  returns NULL → `OPS_GetString` (`OpenSeesCommands.cpp`) substitutes the literal
  string `"Invalid String Input!"`. Your `strtod` then fails, you push the token
  back, the option loop re-reads it (NULL again → unknown-option WARNING), and the
  coefficient is **consumed-and-discarded** — the element silently falls back to the
  default coeff. Tcl path is unaffected, so it hides from Tcl-only testing.
- **Why:** `OPS_GetString` is for *string* options; `OPS_GetDoubleInput` /
  `OPS_GetIntInput` are the number-aware readers (their Python backends accept
  `PyFloat`/`PyLong`/`PyBool`). 
- **Fix / idiom for an OPTIONAL trailing numeric arg that may instead be the next
  flag:** read it with `OPS_GetDoubleInput(&n1,&tmp)`; on success use it, on failure
  `OPS_ResetCurrentInputArg(-1)` to un-get it for the option loop. `OPS_GetDoubleInput`
  advances the arg cursor by one even when the conversion fails, so the single
  `-1` reset un-gets exactly that token — works on **both** Tcl and Python.
- **How it surfaced:** adversarial review of LadrunoBrick PR [#75](https://github.com/nmorabowen/OpenSees/pull/75)
  (eas + viscous). The first `-hourglass <type> -lumped` fix used `OPS_GetString`+`strtod`;
  caught before merge by adding `test_hourglass_coefficient_reaches_kernel` (a numeric
  Python-float coeff must change the response). Date learned: 2026-06-01.
