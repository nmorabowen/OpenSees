---
wp: WP-232
title: "Stock recorders cache Node and Element pointers, and Domain never calls Recorder::domainChanged(): removing a watched entity is a use-after-free"
date: 2026-11-13
---
### Stock recorders cache `Node*` / `Element*` / `Response*`, and `Domain` never calls `Recorder::domainChanged()` — removing a watched entity is a use-after-free

- **Bites:** any runtime topology change — `remove element`, `remove node`, element-removal workflows (ADR-51), PFEM remeshing, an MPM grid that retires cells or particle elements (WP-232) — while a stock recorder (`NodeRecorder`, `ElementRecorder`, `EnvelopeElementRecorder`, `RemoveRecorder`) points at a removed entity. Re-adding an entity under the same tag does not help: the cached pointer is still the old one.
- **Why:** `Domain.cpp` calls only `record`, `restart` and `flush` on its recorders; nothing calls `Recorder::domainChanged()`. `NodeRecorder` caches `Node*` on its first record and would re-initialise only in its own `domainChanged()` (`NodeRecorder.cpp:1133`), which is never called; the element recorders cache `Response*` / `Element*`.
- **`LadrunoRecorder` is safer, but not a general escape.** Since WP-165 it opens a new model stage only when a topology fingerprint of tags **and object addresses** changes (`LadrunoRecorder.cpp:664-694`). A removal or re-addition creates or destroys objects, so the fingerprint moves and the recorder rebuilds — at the cost of one model stage per change. An element whose nodes change in place (the same object re-pointed, as WP-232's implicit particle elements are) leaves the fingerprint unchanged, so its connectivity in the output goes stale. See `WP-165-domain-hasdomainchanged-is-a-dof-graph-stamp.md`, `WP-163-a-recorder-that-rebuilds-on-the-domain-change.md` and `LEGACY-ladrunorecorder-domainchanged-restart-setdomain.md`.
- **Workaround / status (2026-11-13, WP-232 — disclosed, not changed):** never point a stock recorder at an entity that may be removed. Use `recorder ladruno` when the set of entities changes at runtime, accepting a model stage per change. For entities that move or re-point every step, write a dedicated output (WP-232's particle group in `.ladruno`).
