---
wp: LEGACY
title: "Under system MPIDiagonal the global equation numbering exists only TRANSIENTLY — setSize rewrites every DOF id to rank-local 0..n−1"
legacy_seq: 206
---
### Under `system MPIDiagonal` the global equation numbering exists only TRANSIENTLY — setSize rewrites every DOF id to rank-local 0..n−1
- **Bites:** any oracle/debug dump of equation ids taken AFTER analysis setup shows rank-LOCAL ids (shared boundary nodes legitimately disagree across ranks); naive cross-rank identity checks fail on correct runs. MUMPS SOEs do NOT localize — dumps keep globals.
- **Why:** `MPIDiagonalSOE::setSize` builds its shared-DOF exchange from the globals, then deliberately compacts per rank ("renumber DOFs 0 through size") and has FE elements re-cache.
- **Workaround/status:** by design, no physics defect. Numbering oracles must be two-deck: `system Mumps` for strict global identity, MPIDiagonal for end-state identity (`tests/test_adr74_numberer_1.py`). *2026-07-22 (ADR-74 N0's first catch).*
