---
wp: LEGACY
title: "Shared append-point files conflict on stale branches — append at END + merge=union"
legacy_seq: 61
---
## Shared append-point files conflict on stale branches — append at END + `merge=union`

Every new fork feature touches the same hotspot files (`SRC/classTags.h`, the
broker `FEM_ObjectBrokerAllClasses.cpp`, the `OpenSees*Commands.cpp` registries,
`*/CMakeLists.txt`, the banner `banner_features.txt`/`tclMain.cpp`/`PythonModule.cpp`,
and the `LEDGER_*.md` / testbed `manifest.yaml` bookkeeping). When a feature branch
falls behind `ladruno` (e.g. #127 was 48 commits stale), these are exactly where
merge conflicts land — git auto-merges *additions at different lines*, but two edits
to **adjacent** lines (a new row inserted next to a row another PR also edited) do
NOT auto-merge.

Prevention (all three, not just one):
- **Reconcile with latest `ladruno` right before merging** (rebase or merge-from-base
  — equivalent under squash-merge; rebase = linear + force-push, merge = no force-push).
  This fixes staleness, but NOT contention.
- **Append new entries at the END of a list/section, never interleaved.** A `classTags.h`
  tag appended after the last one auto-merges; a `LEDGER_*.md` row inserted *mid-table*
  next to a row another PR edits will conflict (the #127 case — its LadrunoJ2Finite row
  sat above the UniaxialJ2 row that the Lemaitre PR was editing).
- **`merge=union` driver** (set in `.gitattributes`) for the append-only logs
  (`LEDGER_*.md`, `banner_features.txt`) — git keeps BOTH sides instead of conflicting.
  NEVER apply `union` to source code (it interleaves → garbage). Learned 2026-06-02
  (the #127 rebase past the Lemaitre-damage merge; [[16_finite_native_j2_adr]]).
