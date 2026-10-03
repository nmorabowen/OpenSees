---
wp: ADR-74
title: "ADR-74 flush PR -- upstreamable-table row(s)"
files: ["`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "upstreamable"
legacy_seq: [345, 346]
---
| `SRC/tcl/commands.cpp` | `// Ladruno` (ADR-74 pre-G3 hardening): `profiler checkpoint <file> [-run id]` subcommand in `TclCommand_profiler` — mid-run snapshot (same lowering as `report`, via `mergedRollup(quietLive)`), OVERWRITES via `<file>.tmp` + rename so a walltime-killed/crashed run yields the last checkpoint instead of nothing. Kill-gated: a run killed at step 9/20 left all 8 ranks' checkpoints readable with full dc.* attribution + the 8-step per-step series. | ADR-74 flush PR |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` (ADR-74 pre-G3 hardening): the same `checkpoint` subcommand in `OPS_profiler` (interpreter twin) + `<cstdio>`/`<string>` includes. | ADR-74 flush PR |
