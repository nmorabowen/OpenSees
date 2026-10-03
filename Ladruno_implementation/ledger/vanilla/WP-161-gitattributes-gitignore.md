---
wp: WP-161
title: "Ledger union lines removed; the ledger build dir ignored"
pr: "#915"
date: 2026-10-03
files: ["`.gitattributes`", "`.gitignore`"]
---
| `.gitattributes` | `# --- Ladruno: append-only bookkeeping files ---` block (WP-161): the three `LEDGER_*.md merge=union` lines removed (the ledgers are per-WP fragments now); `Ladruno_scripts/banner_features.txt merge=union` kept, comment rewritten. | [#915](https://github.com/nmorabowen/OpenSees/pull/915) |
| `.gitignore` | WP-161: `Ladruno_implementation/ledger/_build/` (the generated LEDGER_*.md) appended at the end. | [#915](https://github.com/nmorabowen/OpenSees/pull/915) |
