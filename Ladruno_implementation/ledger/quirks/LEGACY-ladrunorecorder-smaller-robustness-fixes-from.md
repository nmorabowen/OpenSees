---
wp: LEGACY
title: "LadrunoRecorder smaller robustness fixes from the adversarial review (2026-06-03)"
date: 2026-06-03
legacy_seq: 49
---
### LadrunoRecorder smaller robustness fixes from the adversarial review (2026-06-03)
- **Mixed-dimension node OOB (FIXED):** `writeModelNodes` latched `ndim` from the *first*
  node then read `crds[1]`/`crds[2]` unchecked for every node — a node carrying fewer coords
  than `ndim` read past its `Vector`. Upstream MPCO guards this; the port had dropped it.
  Fix: clamp each read to `crds.Size()` (pad 0.0), mirroring the GLOBAL_GP path.
- **`eo_response` leak (FIXED):** in `initElementSources`, a `CompositeResponse` built but
  then rejected because `eo_stream.error_code != OK` (e.g. `ERROR_CODE_GENERIC`) was owned by
  nobody. Fix: `else if (eo_response) delete eo_response;` after the registration block.
- **Recorder `exit(-1)` aborting the analysis (FIXED):** the `OutputDescriptorStream`
  tag/attr parser and `mapElements` in `Ladruno_ElementResults.h` called `exit(-1)` on an
  element output-tag nesting they didn't expect (invalid parent for SectionOutput/FiberOutput,
  invalid tag at level, empty item-list at a walked level) or on two same-classTag elements
  with differing `getNumExternalNodes()` — a *recorder* killing the whole run. **Fix:** the 7
  stream sites now `error_code = ERROR_CODE_GENERIC; return -1;` (the offending element's bucket
  is dropped by the existing `error_code != OK` gate — same mechanism `ensureItemsOfUniformType`
  already used at line ~1036); `mapElements` now `continue;`s past the inconsistent element.
  No live `exit(-1)` remains (two pre-existing commented-out ones at ~770/~1033 untouched).
  Happy path unchanged (full recorder regression green). A runtime trigger needs a custom
  element that emits an unsupported tag nesting (not reachable from standard openseespy), so
  this is verified by the no-regression run + reuse of the already-proven GENERIC drop-path.
- **Still open:** `domainChanged()`/`restart()`/`setDomain()` are no-ops, so the only
  source-rebuild trigger is the `hasDomainChanged()` stamp inside `record()` (cached
  `Element*`/`Response*` can dangle if a model edit doesn't move the stamp across a record).
