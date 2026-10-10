---
wp: LEGACY
title: "Recorder-token aliases must never change the emitted ResponseType labels"
legacy_seq: 246
---
### Recorder-token aliases must never change the emitted `ResponseType` labels
- **Bites:** the obvious way to add an alias is to canonicalize the token and let one branch serve several spellings — which is right — but if the alias is allowed to influence what the branch *emits* (the `output.tag("ResponseType", ...)` labels, the response ID, or the vector width), then a recorder's column headers start depending on how the user spelled the token in their script. MPCO/STKO readers key on those labels; `stress` and `stresses` producing different metadata for the same numbers is worse than `stress` not working at all.
- **Rule:** `LadrunoResp::is()` (`SRC/element/LadrunoResponseTokens.h`) decides only WHICH branch runs. Everything the branch emits stays keyed to the branch. Matching is symmetric, so an element cannot accidentally test against a non-canonical spelling, and an unregistered token falls back to an exact `strcmp` — never worse than upstream.
- **Corollary for new elements:** canonicalization is GLOBAL across the table. Before giving a new element a short token (`k`, `dir`, `L`, `gap`, …), check the table — those are already claimed by the penalty/tie and rigid-body families. *2026-07-28 (recorder-token consistency sweep).*
