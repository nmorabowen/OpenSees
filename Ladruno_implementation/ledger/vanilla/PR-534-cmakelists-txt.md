---
wp: PR-534
title: "534 -- 1 vanilla row(s)"
pr: "#534"
files: ["`CMakeLists.txt` (root)"]
table: "main"
legacy_seq: [268]
---
| `CMakeLists.txt` (root) | `# Ladruno` ADR43 L2-profile: opt-in `option(LADRUNO_MKL_FEAST_LINUX OFF)` — enables MKL FEAST on a non-Windows host (esmeralda). When ON+non-WIN32, find the MKL 3-layer (`mkl_intel_lp64`/`mkl_sequential`/`mkl_core`, HINT `MKL_RT_HINT`), define `_LADRUNO_MKL_FEAST`, and attach macro (PRIVATE) + libs (PUBLIC) to `OPS_SysOfEqn` (the sole compiler of the FEAST sources) so all consumers link MKL in lockstep. **Default OFF ⇒ Windows/oneAPI + Ubuntu Zone-A byte-identical** (verified). | [#534](https://github.com/nmorabowen/OpenSees/pull/534) |
