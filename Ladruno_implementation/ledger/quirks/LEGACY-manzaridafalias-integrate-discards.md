---
wp: LEGACY
title: "ManzariDafalias::integrate() discards BackwardEuler_CPPM's return value, and the CPPM's own ladder can never fail anyway — a scheme-2 non-convergence is invisi…"
legacy_seq: 462
---
### `ManzariDafalias::integrate()` discards `BackwardEuler_CPPM`'s return value, and the CPPM's own ladder can never fail anyway — a scheme-2 non-convergence is invisible in every channel
- **Bites:** a CPPM step whose Newton diverged, whose Jacobian was singular, or which recursed
  through up to 512 half-increments and then gave up looks EXACTLY like a clean implicit return —
  no return code, no `opserr` line, no response. Measured cost of one such invisible failure on a
  single-element drained triaxial at `dEz = 1e-4`: **133.75 s for one step against 30 ms for its
  neighbours** (WP-105 / F12).
- **Why (verified on `634824e1f`):** `ManzariDafalias.cpp:1023-1027` calls
  `BackwardEuler_CPPM(...)` for `mScheme == INT_BackwardEuler` with no assignment — the return
  value is simply discarded. It would not matter even if it were kept: the `while(errFlag != 1)`
  ladder always terminates by falling through to `explicit_integrator` and setting
  `errFlag = 1` (`:2584-2588`), and the low-`p` branch does the same (`errFlag = 0` then explicit
  then `errFlag = 1`, on this checkout at `:2418`/`:2431`/`:2436` — the guide's §3 citation of
  `:2264` for the low-`p` branch has moved on this build; WP-105 recorded the new lines rather than
  editing the guide's prose, which is about the mechanism, not the line number). The recursion
  itself increments `implicitLevel` at `:2538` and is capped at `implicitLevel > mMaxSubStep = 10`
  (`:2352-2357`), i.e. up to `2^9` halvings, each its own 19-unknown Newton of up to 30 iterations
  with a 19x19 solve. With `ManzariDafalias::debugFlag` a compile-time `const bool = false`
  (`:57`), none of this prints. The only observable is the shipped `substeps` response
  (`LadrunoSANISAND::setResponse`, `LadrunoSANISAND.cpp:4044-4058`), whose first component
  `mSubstepsTakenInME` is non-zero **only if** `ModifiedEuler` ran — an exact detector of the CPPM
  falling back on scheme 2, but silent whenever the ladder succeeds by substepping rather than by
  falling through. The fork's only trial-time refusal for this material,
  `LadrunoSANISAND::ladrunoUpdateStatus()` (`LadrunoSANISAND.cpp:3993-3996`), returns
  `LADRUNO_MATERIAL_REFUSED` only via `mSubstepCapHitInME`, set at exactly one site inside
  `ModifiedEuler()` (`ManzariDafalias.cpp:1578-1600`) and gated on `mMaxSubstepsInME` — so on scheme
  2 a refusal can arise only AFTER the CPPM has already given up, and only if `-maxSubsteps > 0`
  (next entry). The F7 element-refusal roster (this file, "element refusal roster", PR #838) is
  irrelevant to a bare CPPM failure: there is nothing for the element to forward.
- **Workaround/status (2026-09-16, WP-105 / F12, no code changed):** not fixed. A fix would need
  (a) `integrate()` to check `BackwardEuler_CPPM`'s return and (b) a way to surface it that does not
  depend on `debugFlag` (compile-time off) or on the ladder actually falling through to
  `ModifiedEuler`. Until then, treat any scheme-2 run through a global Newton as unauditable for
  silent quality loss — measure cost and stall rate (this entry's numbers), not correctness,
  because correctness has no channel to fail loudly through.
- **Status (WP-130, #868):** now VISIBLE and, on request, a REFUSAL. `substepStats` columns
  17-27 count every CPPM call, local-Newton failure, half-increment, **silent explicit fallback**
  (`cppmExplicitFail`, vanilla's path, still taken at the defaults) and low-p explicit branch, per
  integration point. `-cppmOnFail refuse` turns the terminal explicit branch into
  `LADRUNO_MATERIAL_REFUSED` (`BackwardEuler_CPPM` returns -5, `ladrunoUpdateStatus()` ORs
  `mLadrunoCPPMRefused`), so a forwarding element cuts the step. `integrate()` still discards the
  return value; the refusal travels in the flag. Default byte-identical.
