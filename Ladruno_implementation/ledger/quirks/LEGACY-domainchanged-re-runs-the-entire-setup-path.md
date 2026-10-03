---
wp: LEGACY
title: "domainChanged() re-runs the ENTIRE setup path (handle + numberDOF + setSize) on EVERY domain change — it is a K× cost, not a one-time one — and eigen() routes…"
legacy_seq: 209
---
### `domainChanged()` re-runs the ENTIRE setup path (handle + numberDOF + setSize) on EVERY domain change — it is a K× cost, not a one-time one — and `eigen()` routes through it too
- **Bites:** treating "setup" as paid once at step 1. It is re-paid on every `domainChange`: apeGmsh emits one per **stage** (gravity → dynamic, staged construction, SSI decks — recorder MODEL_STAGE splitting, apeGmsh PR #633), ADR-51 element removal bumps the stamp per removal event, ADR-55 contact re-discovery likewise, and `DirectIntegrationAnalysis::eigen()` (`:323`) calls `domainChanged()` — so modal/FEAST runs pay it too. A 20-event removal history at 19 M re-crosses the hour line **even post-fix** if any setup term were still super-linear; a progressive-collapse (AEM) run at 10 M+ with hundreds of events multiplies whatever residual remains.
- **Why:** `domainChanged` unconditionally re-forms DOF groups, re-numbers, and re-sizes the SOE — there is no "incremental re-setup" path; a domain that changed by one element pays the full O(model) again.
- **Workaround/status:** the reason the ADR-74 fixes had to make EVERY setup term linear, not just fast-at-K=1: the numberer (T0/T1), setSize (#593), handle (#595), and the MP/EQ `-4` sweeps (#598) are all now ~linear, so K× is bounded. The per-run benchmark hides this — a single-`domainChange` deck (the G3 plane wave) is the K=1 best case and must be flagged as such; staged/removal/contact decks are where K bites. `LadrunoParallelPlain` (no RCM pass) further cuts the per-event cost on the explicit + MUMPS lanes. *2026-07-22.*
