---
wp: WP-104
title: "841 -- upstreamable-table row(s)"
pr: "#841"
files: ["`SRC/material/nD/NDMaterial.cpp`"]
table: "upstreamable"
legacy_seq: [626]
---
| `SRC/material/nD/NDMaterial.cpp` | `// Ladruno WP-104`: `OPS_clearAllNDMaterial()` — the nD-material wipe hook every interpreter and `OpenSees.exe` go through — now also calls `ladrunoSanisandResetImplexGlobals()` (defined in fork-owned `LadrunoSANISAND.cpp`, reached through a local `extern`, vanilla's own idiom for the `OPS_clearAll*` hooks in `commands.cpp:55-63`). Zeroes `LadrunoSANISAND`'s process-wide IMPL-EX ledger (`implexRefusals` slots 0-3 + 5, `implexGuards`, `avgImplexError` + commit-round marker) on `wipe`. **The defect:** a FRESH material in a NEW model inherited the previous model's refusal totals (measured `[9,0,0,9,0,9]` on a new tag after `ops.wipe()`, apeGmsh live test 2026-09-15), so any end-of-run `implexRefusals[3] == 0` check depended on which model ran first in the process. Same rule as the ADR-69/72 energy-channel reset in `Domain::clearAll()`; hooked here, not there, because `Domain::clearAll()` also runs from `recvSelf()` on an MP rank, which is not a wipe. One call per wipe; stock decks byte-identical. | [#841](https://github.com/nmorabowen/OpenSees/pull/841) |
