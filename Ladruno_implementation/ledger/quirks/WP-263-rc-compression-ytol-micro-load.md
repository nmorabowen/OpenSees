---
wp: WP-263
title: "LadrunoRCConcrete: a principal stress of -1e-12 loads the compression backbone where q is clamped, and the tangent there is garbage"
date: 2026-11-28
status: "open"
---
### RC exact tangent d22 ~ 0 (or x1000) at a plane-stress root: the compression backbone "loads" at x < YTOL/E

- **Bites:** anyone differentiating `LadrunoRCKernel.h` exactly (`-dualTangent`, an FD of the map,
  a future consistent tangent), and, latently, the default mode-0 tangent. Seen first as a
  PlaneStress σ33 Newton that cycles (`tests/_testbed/rc_dual_gpp.cpp`, `ps_interlock_cyclic`
  step 1): at σ33 = -1e-12 the exact d22 was -7e-12, so the Newton step `σ33/d22` jumped and the
  condensed tangent `D - D:2 D2: / d22` would divide by ~0.
- **Why:** the tensile measure has a `STRESS_TOL = 1e-9` gate (`equivTensile`), the compressive
  one has none (`equivCompressive`). So any slightly negative principal stress gives
  `xc_trial > pc.x` with committed `xc = 0`, and `evaluateAt(hc, x)` runs at x ~ 3e-17, where
  `E x < YTOL = 1e-12`: y and q are clamped to YTOL and
  `dc_plastic = 1 - q / seff_eq_c = 1 - YTOL/(E x)` is noise (very negative as x -> 0). The
  stress `beta (1 - dc) SC` stays ~YTOL, so no stress test sees it. Its derivative is a
  micro-plateau (flat, d22 ~ 0) at E x ~ YTOL and is multiplied by YTOL/(E x) below that. The
  mode-0 secant has the same factor through `wc = 1 - beta (1 - dc)`, so it is amplified there too.
  It is rarely hit only because the plane-stress Newton usually lands on σ33 = +1e-16.
- **Fix (WP-263, Dual path only):** in `returnMap3DT`, when `value(pc.q) <= YTOL` the compressive
  term takes the elastic-side derivative `beta dSC` (the limit from x > YTOL/E, where dc = 0
  exactly). The value is untouched: `fix = cterm - beta SC; zero_value(fix); cterm -= fix`
  (subtract, so a -0.0 stays -0.0). The double path is unchanged (bit-identity).
- **Not fixed:** the double-path value of `dc_bar` in that micro range, and the mode-0 tangent
  amplification. Both need a value change, so a `STRESS_TOL` gate on `equivCompressive` or a
  clamp of `dc_plastic` at 0, behind its own bit-identity decision. Do not "fix" the Dual by
  clamping the derivative of dc only: the value `(1 - dc)` still multiplies `dSC`.
