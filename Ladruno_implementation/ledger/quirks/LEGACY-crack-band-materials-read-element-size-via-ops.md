---
wp: LEGACY
title: "Crack-band materials read element size via ops_TheActiveElement->getCharacteristicLength() — and the base default is wrong for high-order elements"
legacy_seq: 13
---
### Crack-band materials read element size via `ops_TheActiveElement->getCharacteristicLength()` — and the base default is wrong for high-order elements
- **Bites:** `ASDConcrete3DMaterial` (and any crack-band/smeared-crack material) auto-
  regularizes its softening branch by the element characteristic length `lch`,
  read **once** on the first `setTrialStrain` via the global
  `ops_TheActiveElement->getCharacteristicLength()`
  (`ASDConcrete3DMaterial.cpp:1614`, guarded by `regularization_done`,
  `ASDConcrete3DMaterial.h:386`). If your element doesn't override
  `getCharacteristicLength()`, it inherits `Element::getCharacteristicLength()`
  (`SRC/element/Element.cpp:682`), which returns the **minimum inter-node
  distance**. On a quadratic element (BezierTri6, BezierTet10, TenNodeTetrahedron,
  the quad/brick "...N" siblings) the closest node pair is corner-to-mid-edge,
  i.e. ≈ **½** the true edge length → `lch` under-estimated ~2× → fracture energy
  smeared over too small a band → **over-softening / spurious snap-back**.
- **Why:** the global is set by the framework, not the element — `Domain::addElement`
  (`Domain.cpp:447`) and, crucially, the `Domain::update()` loop
  (`Domain.cpp:2263`) sets `ops_TheActiveElement = theEle` immediately before
  `theEle->update()`. So as long as the element pushes strain **only inside
  `update()`** (both Bézier elements do), the active pointer is correct when the
  one-time regularization fires — no wrong-element window. The value is just
  geometrically poor for the min-distance default on high-order elements.
- **Also:** `getCharacteristicLength()` returns a **single element scalar**, read
  once per material instance — so a true per-Gauss-point `lch = sqrt(detJ·w)` is
  *not* expressible through this seam; every GP material on the element gets the
  same value. Pick one representative element size.
- **Workaround/status (2026-06-01):** override `getCharacteristicLength()` with an
  element-size equivalent from the integrated area/volume — BezierTri6 returns
  `sqrt(2·A)` (right-isosceles-triangle edge), BezierTet10 returns `cbrt(6·V)`
  (right-tetrahedron leg); both recover the true edge length on right-simplex
  meshes and err in the safe (under-estimating) direction on curved Bézier edges.
  LadrunoBrick is exempt — its 8 corner nodes give the correct min-edge already,
  and its only strain-push site is `update()` (verified read-only assembly).
  Done in `BezierTri6.cpp` / `BezierTet10.cpp` `getCharacteristicLength()`.
