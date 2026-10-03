---
wp: LEGACY
title: "printA('-ret') (dense) is EMPTY for every SOE except FullGeneral"
legacy_seq: 394
---
### `printA('-ret')` (dense) is EMPTY for every SOE except `FullGeneral`
- **Bites:** `ops.printA('-ret')` returns nothing under `UmfPack` and the other sparse solvers — looks like "no tangent" rather than "wrong API" (`getA()` is null for non-dense SOEs, `OpenSeesCommands.cpp:2718`). `FullGeneral` itself crashes on fully-prescribed drivers (N = 0).
- **Rule:** `printA('-sparse', '-ret')` works with `UmfPack` and returns `{rowIndices, colIndices, values}`. It calls `formTangent()` itself, so it reads the tangent after the last `update()` — exactly the window the ASDP static-tangent defect lives in.
