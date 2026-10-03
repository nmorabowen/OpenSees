---
wp: LEGACY
title: "ops_TheActiveElement is a mutable GLOBAL written per-element inside Domain::update()'s loop and read by materials — a static-only re-entrancy audit misses it e…"
legacy_seq: 215
---
### `ops_TheActiveElement` is a mutable GLOBAL written per-element inside `Domain::update()`'s loop and read by materials — a `static`-only re-entrancy audit misses it entirely
- **Bites:** you audit element/material kernels for `static` scratch, conclude the `update` loop is clean, thread it, and get plausible-but-wrong regularized softening. The hazard is not a `static` — it is a file-scope global: `Element *ops_TheActiveElement` (`SRC/element/Element.cpp:47`, `extern` in `SRC/G3Globals.h:45`), assigned **per element inside the loop** at `Domain.cpp:2401` (and in `Element`'s ctor `Element.cpp:65`, `Domain.cpp:461`, `OpenSeesCommands.cpp:2865`, `LadrunoDispBeamColumn2d.cpp:528`, `3d.cpp:651`).
- **Why:** it is a deliberate fork idiom — the "Phase-3b lch latch". Materials read it on first `setTrialStrain` to fetch a regularization characteristic length: `LadrunoJ2.cpp:353`, `LadrunoConcrete3D.cpp:347`, `ASDConcrete3DMaterial.cpp:1614`, `LadrunoRCConcrete.cpp:330`, `LadrunoRCFiniteStrain.cpp`, plus documented reliance in `BezierTet10.cpp`, `BezierTri6.cpp`, `LadrunoUP.cpp`.
- **Failure mode:** an element's material latches **another element's** characteristic length ⇒ a converged, plausible, wrong softening response. Silent.
- **Workaround/status:** must become `thread_local` (or be plumbed down the call chain) before `Domain::update()` is threaded — a prerequisite of [[75b_ladruno_threaded_assembly_adr]] L3-1, and it touches vanilla files (`Element.cpp`, `G3Globals.h`) so it owes [[LEDGER_vanilla_files]] rows when it lands. Note there are also several unrelated *definitions* of the same symbol in non-linked TUs (`NeesDataTest.cpp:45`, `TestDataOutput{Database,File,Stream}Handler.cpp:49`) — don't mistake those for the live one. *2026-07-25 (ADR-75b L3-0).*
