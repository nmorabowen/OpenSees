---
title: WP-249 — MPM hardening before phase 7 (#1095, #1110, #1103, #1106)
project: Ladruno
tags:
  - mpm
  - wp-249
---

# WP-249 — MPM hardening before phase 7

Program #1075. A fix/hardening WP queued by the owner on 2026-11-22, to land before
phase 7 (232m, multi-body contact). It is not a §9 phase of ADR 232.

Scope, in priority order:

1. **#1095** — the mode-A crash on a step-cut retry (`Matrix::operator*(Vector):
   incompatible sizes`, then 0xC0000005). Root cause, fix, and a deterministic test that
   forces cuts across a connectivity change.
2. **#1110** — check whether 232l's `mpm traction` is the single-phase surcharge carrier
   #1110 asks for, in both modes; demonstrate it on 232k's S2 column, or record the gap.
3. **#1103, #1106** — the small, local items of the 232i and 232j code-review follow-ups.

Out of scope: #1107, #1124–#1127, #1141, contact.

## Status

Draft. Findings and results are added below as they land.
