---
wp: WP-232
title: "Classic Tcl prints OPS_SetDoubleOutput values as %35.20f, fixed-point"
date: 2026-11-13
---
### Classic Tcl prints `OPS_SetDoubleOutput` values as `%35.20f`, fixed-point

- **Bites:** a test that reads numbers back from `OpenSees.exe` (`puts [cmd]`) gets 20
  decimals, not 17 significant digits: a value below 1e-20 reads as 0 and one near 1e-12
  keeps only 8 digits, while Python returns the full double. A Python-vs-Tcl parity check
  with a purely relative tolerance then fails on small values.
- **Workaround:** compare with an absolute floor (>= 1e-15) plus a relative term, as
  `tests/test_ladruno_mpm_engine.py::test_mpm_engine_parity` does.
- Learned 2026-11-13 (WP-232c).
