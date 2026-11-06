---
wp: PR-568
title: "568 -- 1 vanilla row(s)"
pr: "#568"
files: ["`SRC/{interpreter/OpenSeesMiscCommands.cpp,tcl/commands.cpp,runtime/commands/domain/commands.cpp}`"]
table: "main"
legacy_seq: [300]
---
| `SRC/{interpreter/OpenSeesMiscCommands.cpp,tcl/commands.cpp,runtime/commands/domain/commands.cpp}` | `// Ladruno`: backport upstream [`191c67c2d`](https://github.com/OpenSees/OpenSees/commit/191c67c2d) — `InitialStateAnalysis on\|off` deleted the `InitialStateParameter` right after `Domain::addParameter` stored the pointer → double free at `wipe()`/exit (0xC0000374/0xC0000005 heap corruption, found by the ADR-71 P4 init battery), and the fixed tag-0 collision froze `ops_InitialStateAnalysis` on repeat calls (off-after-on silently never flipped the flag → post-ISA `ops.reset()` kept PM4Sand-family state / compounded `InitialStateAnalysisWrapper` ε₀). Fix = stack `InitialStateParameter` + `setDomain` only (nothing registered; adversarial panel 4/4 — SP forwarding it drops was a stuck-on switch that could never propagate `off`). Regression: `tests/test_initial_state_analysis_lifetime.py`. The `runtime/` copy is unbuilt — patched for upstream-sync parity only. | [#568](https://github.com/nmorabowen/OpenSees/pull/568) |
