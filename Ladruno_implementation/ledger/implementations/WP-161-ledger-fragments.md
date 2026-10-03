---
wp: WP-161
title: "Per-WP ledger fragments; the LEDGER_*.md are generated"
date: 2026-10-03
status: "draft"
---
| **WP-161 — per-WP ledger fragments; the three LEDGER_*.md are generated** ([[161_ledger_fragments]]; format [[ledger/README]]). One file per ledger entry under `Ladruno_implementation/ledger/{implementations,quirks,vanilla}/`, so no two PRs write the same ledger file (55 of 64 conflicted merge-ups over the 60 PRs before it were ledger-only, and GitHub's server-side merge ignored `merge=union`). `ci/ledger.py` builds the ledgers into the gitignored `ledger/_build/` (CI artifact `ledgers`), splits the old files once, migrates a pre-WP-161 branch's ledger edits; gate `ci/check_ledger_fragments.py` (schema, name, duplicates, class tags, banner, stubs); quirk lint L3 and viewer gate V1 read the fragments; the three ledger `merge=union` lines are gone. Round trip build(split(ledgers)) == ledgers proven in `ci/test_ledger.py`. | CI tooling + docs | — (no class tag) | `ci/ledger.py`, `ci/check_ledger_fragments.py`, `ci/test_ledger.py`, `ci/check_quirk_patterns.py` (L3), `ci/check_viewer_ledger.py` (V1), `.github/workflows/ladruno.yml`, `Ladruno_implementation/ledger/` | draft | TBD |
