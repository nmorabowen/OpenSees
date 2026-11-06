---
wp: PR-901
title: "901 -- upstreamable-table row(s)"
pr: "#901"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`"]
table: "upstreamable"
legacy_seq: [676]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno WP-158` — **`ForwardEuler` (IntScheme 5; also 4, 7–9 on small increments and the WP-130 `-cppmStart` walk): three defects, numerics of every other scheme untouched.** (1) `Vector r = GetDevPart(CurStress) / p;` inside `if (p > small)` declared a NEW `r` that shadowed the outer one, so `r` stayed zero and both `(n:r)` terms of the plastic multiplier vanished (upstream master has the same line) → `r = ...`. (2) tangent `temp2 = 2G n - (n:r) I` → `2G n - K (n:r) I` (the multiplier's numerator is `2G n:de - K de_v (n:r)`). (3) tangent `temp1 = 2G mIIdevMix + K mIIvol` (2G on the shear diagonal) → `aC` = `GetStiffness(K, G)`. Gated by `tests/test_manzari_forward_euler_r.py` (one-step consistency order, drained triaxial vs IntScheme 1, tangent vs FD; all three fail on d63f49750). WP-129 byte-identity: only `ls3d_s5` moved, re-pinned. LEDGER_quirks: "has a shadowed `Vector r`". | [#901](https://github.com/nmorabowen/OpenSees/pull/901) |
