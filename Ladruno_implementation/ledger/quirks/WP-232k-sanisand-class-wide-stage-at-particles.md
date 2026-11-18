---
wp: WP-232k
title: "LadrunoSANISAND's stage is one class-wide static: a flip that reaches one point flips all, and the rest skip Elastic2Plastic"
date: 2026-11-18
---
### SANISAND's stage flag is class-wide; only the points `updateParameter` reaches get `Elastic2Plastic`

- **Bites:** a SANISAND FE model and SANISAND MPM particles in one process, or any mesh where
  `updateMaterialStage` stops at the first element that answers (`MaterialStageParameter.cpp:79`).
  `mElastFlag` is a static on `ManzariDafalias` (`ManzariDafalias.h:557`) shared by every instance
  of every tag. Points reached by `updateParameter` run `Elastic2Plastic()` (the p < Pmin reset,
  `:5787-5795`). The others flip lazily in `ladrunoTrialUpdate` (`LadrunoSANISAND.cpp:3397-3406`),
  and that path skips the reset. A second trap is the double path:
  `parameter ... updateMaterialStage` + `updateParameter` reads id 1's `info.theInt`, which is
  indeterminate there (`LadrunoSANISAND.cpp:5253-5257`); use `materialState` (id 5, theDouble).
- **At MPM particles (WP-232k):** `mpm stage <tag> <s>` walks every particle through
  `materialState`, and the `parameter ... loadPattern` and mode-A `updateMaterialStage` routes
  enter the same walk. `mpm analyze` is refused when the class-wide flag moved by a route that
  missed the particles, and particles cannot be added after a non-zero walk (add all bodies
  first).
