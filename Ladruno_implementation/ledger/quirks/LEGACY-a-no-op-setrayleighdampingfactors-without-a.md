---
wp: LEGACY
title: "A no-op setRayleighDampingFactors WITHOUT a getDamp override ⇒ hard crash in implicit transient"
legacy_seq: 74
---
### A no-op `setRayleighDampingFactors` WITHOUT a `getDamp` override ⇒ hard crash in implicit transient
- **Why it bites:** a pure-coupling/penalty element often overrides
  `setRayleighDampingFactors(...)` to a no-op `return 0` to refuse Rayleigh damping (so a
  `betaK` can't spuriously shrink its explicit `dt_cr`). But the base `Element::getDamp()`
  (`SRC/element/Element.cpp:211`) does `if (index==-1) this->setRayleighDampingFactors(...)`
  then `theMatrices[index]->Zero()` — and it is the BASE `setRayleighDampingFactors` that
  lazily allocates `theMatrices[index]` and sets `index>=0`. A no-op override never
  allocates, so `index` stays at its ctor default −1 and `theMatrices[-1]` is an
  out-of-bounds dereference → **hard crash**. `FE_Element::addCtoTang(fact)` calls
  `getDamp()` whenever `fact!=0`, and the Newmark/HHT velocity coefficient
  `c2=γ/(βΔt)` is ALWAYS nonzero ⇒ getDamp fires in EVERY implicit transient step. The
  residual damping-force path `addD_Force` calls it too. `getMass`,
  `getResistingForceIncInertia`, and `getRayleighDampingForces` share the same
  `theMatrices[index]` landmine — the first two are usually already overridden (so safe);
  **`getDamp` (and `getRayleighDampingForces`) are the ones people forget.**
- **Why quasi-static tests miss it:** `LoadControl`/`Static` never form a C-tangent;
  `CentralDifference` (explicit) dodges it when the model has no Rayleigh. Only an implicit
  transient (Newmark/HHT) — or any transient with Rayleigh — triggers it.
- **Fix:** override `getDamp()` and `getRayleighDampingForces()` to return an element-owned
  ZEROED `Matrix`/`Vector` (sized `nDOF`, allocated alongside the mass matrix), bypassing the
  base index path. `D≡0` is physically correct for a pure coupling; mass/inertia still come
  from `getMass` + bipenalty. Confirmed: `LadrunoDistributingCoupling` (RBE3, 33011) crashed
  a Newmark transient (exit 5) before the override, passes after (regression test added).
  **`LadrunoEmbeddedNode` (33006) + `LadrunoEmbeddedRebar` (33005) had the SAME latent bug**
  (no-op setRayleigh, no getDamp override) — **FIXED 2026-06-09** with the identical
  element-owned zeroed `C0`/`dampF` pattern + a `Newmark 0.5 0.25` `test_transient_newmark_smoke`
  regression in each Zone-A battery (empirically reproduced: pre-fix the smoke test segfaults
  `0xC0000005`, post-fix 77/77 pass). See [[LEDGER_implementations]] rows 33005/33006, PR #220.
  2026-06-07 (RBE3) / 2026-06-09 (embedded).
- **Status (WP-123, 2026-09-25):** the four copies are now ONE base,
  `SRC/element/ladrunoEmbeddedRebar/LadrunoUndampedElement.h` — a new element that ignores
  Rayleigh derives from it. **Correction:** the `getRayleighDampingForces` half of the fix above
  never did anything: `Element::getRayleighDampingForces` is NOT virtual, so the four copies only
  shadowed it and nothing called them (`eleResponse … dampingForce` reaches the base, which answers
  zero because the refused factors stay 0). Only `setRayleighDampingFactors` and `getDamp` are live.
  And since the 2026-07-28 `Element.cpp` fix (next entry, "makes 11 `Element` methods") the base
  `getDamp` no longer crashes either; the override is now defence-in-depth against that vanilla edit
  being lost in an upstream sync. `tests/test_ladruno_undamped_couplings.py` pins the contract.
- **Test-design trap (WP-123 mutation row C):** for an element whose residual carries no D·v, a
  spurious nonzero `getDamp` reaches only the TANGENT, and Newton iterates it away: same converged
  answer, just slower. A "no damping" assertion run under `algorithm Newton` therefore cannot see it.
  In WP-123 it caught 1 of 4 elements, and only because Newton stalled on a tiny mass. Use
  `algorithm Linear` (one solve with the element's own tangent), so the polluted tangent shows up as a
  wrong response. Likewise, comparing damped vs undamped runs in the SAME build cannot see a C that does
  not depend on the Rayleigh factors; use an absolute oracle (energy conservation under Newmark γ=½, β=¼).
