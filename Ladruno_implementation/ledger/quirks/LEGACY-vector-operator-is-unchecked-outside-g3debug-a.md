---
wp: LEGACY
title: "Vector::operator() is UNCHECKED outside _G3DEBUG — a mis-sized static Vector in getResponse is a silent heap overrun, not a caught index error"
legacy_seq: 260
---
### `Vector::operator()` is UNCHECKED outside `_G3DEBUG` — a mis-sized `static Vector` in `getResponse` is a silent heap overrun, not a caught index error
- **Bites:** any `getResponse` that fills a fixed-size scratch `Vector` in a
  Gauss-point loop. Found in vanilla `TenNodeTetrahedron::getResponse`, where
  `static Vector stresses(6)` — copy-pasted from `FourNodeTetrahedron`, which
  has ONE Gauss point so 6 is right there — is written by BOTH the `stresses`
  and `strains` branches looping over all `NumGaussPoints=4` points × 6
  components = **24 doubles into a 6-double block**. 18 doubles / 144 bytes past
  the end, on every recorder step. `Brick` gets it right (`stresses(48)` for
  8 GP), so the pattern is fine — only the size was left behind when the loop
  bound was edited.
- **Why it is invisible:** `Vector::operator()` bounds-checks ONLY under
  `_G3DEBUG` (`SRC/matrix/Vector.h`), which release builds do not define — the
  checked accessor is `operator[]`, which nobody uses in element code. And the
  buffer is `static`, so it is heap-allocated once on first call and the SAME
  144 bytes past it are stomped forever after. The crash therefore surfaces at
  whatever the allocator placed next — often a later `free`/alloc far from the
  write — so the backtrace does not point at the element.
- **Second-order damage even without a crash:** `Information::setVector` does
  `*theVector = newVector`, and `Vector::operator=` REALLOCATES on size
  mismatch. So the `ElementResponse`'s advertised `Vector(6*nGP)` gets silently
  shrunk to the scratch size, while `ElementRecorder` already sized its columns
  from the advertised size at setup (`ElementRecorder.cpp`) and then copies
  `eleData.Size()` per element — the column layout desynchronises across
  elements. Garbage output, no diagnostic.
- **Measured A/B** (one tet10, uniform-strain patch, `eleResponse(ele,'stresses')`):
  pre-fix build returns **6** values, post-fix returns **24** — so the cheap,
  deterministic tell is a response list that is a WHOLE FRACTION of the
  advertised length, one GP block instead of nGP. Both builds report
  `σxx = 10000.0 = E·ε` exactly, i.e. the physics was never wrong — only the
  buffer. Note the 1-element control did NOT crash: the overrun is real on
  every call but whether it segfaults depends on what the allocator put after
  the block, so a small repro proving "no crash" proves nothing. Trust the
  length, not the absence of a crash.
- **Rule:** size the scratch from the same expression `setResponse` advertises
  (`6*NumGaussPoints`, never a literal), and treat any `static Vector`/`Matrix`
  in a response path as a place to check the loop bound against the declared
  size. When a new element is derived by copy-paste from one with a different
  Gauss-point count, audit every fixed size in the file, not just the loops.
  Cross-ref [[ladruno-adr79-geom-hypo]] (tet10 was the recorder in use).
  *2026-08-04 (tet10 recorder segfault).*
