---
wp: WP-232j
title: "NDMaterial: getTangent() and getInitialTangent() can return the SAME member buffer -- hold a copy, never two references"
date: 2026-11-17
---
### `getTangent()` and `getInitialTangent()` may share one buffer: copy the first before calling the second

- **Bites:** any code that holds `const Matrix &Dt = mat->getTangent();` and then calls
  `mat->getInitialTangent()` (or the reverse) before reading `Dt`. `DruckerPragerPlaneStrain`
  fills its single member `tangent` in both getters (`DruckerPragerPlaneStrain.cpp:129-142` from
  `mCep`, `:152-164` from `mCe`) and returns a reference to it, so the held reference silently
  reads the second getter's matrix. Other plane-strain / condensing wrappers may do the same.
- **Seen (WP-232j, L4 r3 diagnosis):** `mpm debug failDump`'s "tangent norm / initial norm" was
  exactly 1.0 for every particle on build `4fc274d51`, plastic or not, so the apex flag never fired.
  The element's own tangent was not affected: `LadrunoMPMPoint::tangent` copies `Dm` into a local
  array before any other material call.
- **Fix:** copy the first matrix (`const Matrix Dt(mat->getTangent());`) before calling the
  second getter (`LadrunoMPMPoint::computeTrial`, the `diagnose` block).
