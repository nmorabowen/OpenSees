---
wp: LEGACY
title: "Staggered-overlay twin models: plain quad b1 b2 = body FORCE/VOLUME, LadrunoUP -body = ACCELERATION — copying the same number silently unloads (or double-loads…"
legacy_seq: 190
---
### Staggered-overlay twin models: plain quad `b1 b2` = body FORCE/VOLUME, `LadrunoUP -body` = ACCELERATION — copying the same number silently unloads (or double-loads) the solid
- **Bites:** building the ADR-73 staggered twin of a monolithic `LadrunoUP` model (plain `quad` + `LadrunoPorousOverlay`), you copy the monolithic `-body 0 $bY` value into the quad's trailing `b1 b2` args. The staggered solid then carries `bY` per unit VOLUME instead of `rho_mix*bY`, and with a hydrostatic overlay `+Q·p` force the net solid load can cancel to ~zero — the P2 battery measured settlement 1e-9 vs 3.5e-4 (u-trace rel diff exactly 1.0) before the fix. Nothing errors; the fluid side looks perfect.
- **Why:** `FourNodeQuad` applies `b` directly (`P(ia) -= dvol*shp*b`, FourNodeQuad.cpp:900 — force density, rho only builds mass), while `LadrunoUP`/upstream `quadUP` scale `-body` by the mixture density (acceleration semantics). Same flag name, different units.
- **Workaround/status (2026-07-18, ADR-73 P2):** staggered twin recipe — quad gets `b2 = rho_mix * bY_accel` (full mixture weight as force/volume; the overlay's `+Q·p` supplies the pore-pressure part of effective stress), overlay `-fluidBody` keeps the acceleration form (`f_seep` scales by `rhoF` internally, matching `LadrunoUP -fluidBody`). Pinned by battery gate (d)(ii) (`tests/test_ladruno_overlay_driver.py`); the P4 guide inherits the recipe.
