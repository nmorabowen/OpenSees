---
wp: WP-170
date: 2026-10-31
title: "ExplicitBathe ignored its documented <computeCriticalTimestep> positional"
---
### ExplicitBathe ignored its documented <computeCriticalTimestep> positional
- **Bites:** `integrator ExplicitBathe 0.54 1` (the documented grammar, and the upstream one) printed "unknown option 1 (ignored)" and left the critical time step off, so `criticalTimeStep` returned -1. The deprecated `ExplicitBatheLNVD <p> <alpha> 1` did the same. `-cfl` worked.
- **Why:** the unified parser read p with the typed getter and then sent every remaining token to the option loop, which knows only `-` flags. The positional flag had no branch.
- **Fix (2026-10-31):** `ebReadCtsPositional()` reads one number after the leading positionals (`ebTok` + `ebIsNumber`, so it works in Tcl and OpenSeesPy) and sets the flag (nonzero = on). A non-number is un-read for the option loop. Found while checking the fork against the jaabell round-6 ExplicitBathe port (up/57), whose test uses the positional form. Gate: `tests/test_explicitBathe_positional_cts.py`.
- **Rule:** a positional named in a command's usage comment needs a parser branch and a test; the option loop silently ignores what it does not know.
