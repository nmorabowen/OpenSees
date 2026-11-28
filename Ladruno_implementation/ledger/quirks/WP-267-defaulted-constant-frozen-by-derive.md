---
wp: WP-267
title: "A default the constructor writes back into a member stops following its driver: the DDM of the driver then misses the FD of a rebuilt model"
pr: "#1176"
date: 2026-11-28
status: "fixed-by WP-267"
---
### DDM `sensNodeDisp` for one parameter is off by tens of percent while the kernel check passes: a defaulted constant (s0 = 0.1·s1) was frozen at construction

- **Bites:** a material whose derive step fills a default INTO the member it was given
  in, e.g. `LadrunoBondSlip::setDerived`: `if (s0 <= 0) s0 = 0.1*s1`. After
  construction `s0 > 0`, so every later re-derive (`updateParameter`, the DDM pass,
  `recvSelf`) keeps it fixed. The DDM wrt `s1` is then d/ds1 at fixed s0. An FD that
  rebuilds the model at s1 ± h, as the end-to-end tests do, moves s0 with s1. On the
  MSVC build of `9327199f5`, `test_ladrunoBondSlip_ddm.py[s1]` read DDM vs FD
  4.429e-01, while the kernel check passed. The check's FD used the same frozen
  semantics as the DDM, so it could not see the gap. With an explicit `-s0 0.1` the
  same leg gave 2.0e-9.
- **Fix (WP-267):** the class keeps a `s0User` flag (s0 > 0 at construction or via
  `updateParameter(s0)`). `setDerived` and the DDM `seededParams` reset a defaulted
  s0 to 0 before deriving, so it is re-derived as 0.1·s1 with its derivative. At
  construction s0 is already 0 there, so no existing value changes. `getCopy` copies
  the flag. `recvSelf` sets it (the flag is not serialized), so after a send a
  defaulted s0 is frozen again. The kernel check's FD now emulates construction at
  theta ± h for the derived defaults.
- **Check for it:** when porting the DDM recipe, list every constant that the derive
  step defaults from another one. `LadrunoCohesiveHinge` already does this right: its
  `userKpen` flag re-derives the default `Kpen = penaltyRatio·Mc²/(2Gf)` from Mc and
  Gf. Make the kernel FD rebuild from the user's inputs rather than re-derive from
  stored ones, or it shares the DDM's blind spot.
