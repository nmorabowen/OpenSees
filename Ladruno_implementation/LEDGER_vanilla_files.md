---
title: Ledger — vanilla files we touched (stub)
project: Ladruno
tags:
  - ledger
---

# Ledger — vanilla files we touched

This file is a stub. Since WP-161 the ledger is one fragment per entry in
[`ledger/vanilla/`](ledger/vanilla/), named `WP-<nnn>-<slug>.md`.
**Never edit this file**; write a fragment (format: [`ledger/README.md`](ledger/README.md)).

- Read: `rg <pattern> Ladruno_implementation/ledger/vanilla/`
- Full ledger: `python ci/ledger.py build` writes `Ladruno_implementation/ledger/_build/LEDGER_vanilla_files.md`
  (gitignored; CI uploads it as the `ledgers` artifact of the static-gates job).

<!-- ledger-stub kind=vanilla split-source=dd38dde987c64271ab23276454fd7d98c3d79f7a -->
