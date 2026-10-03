---
title: Ledger — OpenSees quirks & gotchas
project: Ladruno
tags:
  - ledger
  - quirks
  - gotchas
---

# Ledger — OpenSees quirks & gotchas we learned

Surprising, undocumented, or bug-prone behaviours of upstream OpenSees that
cost us time. Recording them here so we (and future us) stop re-discovering
them. This is observation-only — fixes we actually applied are tracked in
[[LEDGER_vanilla_files]] / [[LEDGER_implementations]].

## Conventions

- **One section per quirk.** Title = the symptom you'd search for.
- State: *what bites*, *why*, *workaround/status*, and the *date* learned.
- If a quirk drove a code change, cross-link the ledger row / PR.
- Deep build/toolchain quirks may live in
  [[../Ladruno_internal/01_compilation_journal]]; link rather than duplicate.
- **If the quirk names a greppable pattern, enforce it**: add a rule to
  `ci/check_quirk_patterns.py` (with a self-test case) instead of relying on
  someone re-reading this file. If it bites element or material authors, add a
  one-line pointer to the matching `.claude/skills/ladruno-new-*` guide; the
  gate's L3 check fails if a pointer stops matching a heading here (WP-115).

## Quirks

<!-- ledger:entries -->
