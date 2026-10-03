---
wp: LEGACY
title: "remove loadPattern does NOT fire domainChange (unless the pattern owned SP_Constraints) — cached pattern pointers dangle silently"
legacy_seq: 197
---
### `remove loadPattern` does NOT fire domainChange (unless the pattern owned SP_Constraints) — cached pattern pointers dangle silently
- **Bites:** any object caching a `LoadPattern*` across steps (a recorder result source, an engine seam, a driver) keeps a freed pointer after `remove loadPattern $tag`: the interpreter command DELETES the pattern object (`OpenSeesMiscCommands.cpp` remove path), but `Domain::removeLoadPattern` calls `domainChange()` only when the pattern carried SP_Constraints — a LadrunoPorousOverlay owns none, so NO domain-change stamp bumps and NO recorder/model rebuild fires. The ADR-73 P4 panel measured the consequence: an early `-overlay` recorder build cached the pattern pointer and wrote subnormal garbage (9.9e-312) into every post-removal row — silent use-after-free, allocator-dependent whether it corrupts or segfaults.
- **Why:** element/node removal goes through paths that mark the domain changed; load-pattern removal is only conditionally marked. Anything keyed off `hasDomainChanged()` (LadrunoRecorder writeModel rebuild, MODEL_STAGE rollover) will NOT observe a pattern removal.
- **Workaround/status (2026-07-18, ADR-73 P4):** never cache a `LoadPattern*` across commits — re-resolve by tag each use (`domain->getLoadPattern(tag)` + classTag check) and act loudly/zero-fill when absent. `OverlayPressureSource::evaluate` and the Monitor overlay path both do this now (battery gate (i) pins it: post-removal rows identically 0.0, file readable). Audit any future pattern-consuming seam for the same hole.
