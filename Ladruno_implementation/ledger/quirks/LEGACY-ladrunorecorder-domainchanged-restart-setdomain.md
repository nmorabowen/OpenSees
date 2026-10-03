---
wp: LEGACY
title: "LadrunoRecorder domainChanged()/restart()/setDomain() being no-ops is INTENTIONAL — do not \"fix\""
legacy_seq: 66
---
### LadrunoRecorder `domainChanged()`/`restart()`/`setDomain()` being no-ops is INTENTIONAL — do not "fix"
- **Why it looks wrong:** an adversarial review flagged that these lifecycle hooks are inert,
  so cached `Element*`/`Response*` could dangle after a model edit. **Verified NON-issue:**
  the *only* source-rebuild trigger is the `domain->hasDomainChanged()` stamp checked inside
  `record()` (the `rebuild_model` block) — and **every** structural edit (`addElement`/
  `removeElement`/etc.) bumps the domain's geometry tag, so the stamp moves and the rebuild
  fires, re-acquiring fresh pointers and (re)writing the MODEL_STAGE. This is the exact frozen
  `MPCORecorder::record()` pattern (the code comment says so). Forcing a rebuild in
  `domainChanged()` would be redundant with the stamp check and risk breaking the multi-stage
  logic. **Leave them as no-ops.** (The only genuinely-real lifecycle gap is minor: the `-T`
  frequency gate can stall if `commitTag` regresses across a second `analyze()` after
  `wipeAnalysis` — a defensive guard, not yet added.) 2026-06-03.
