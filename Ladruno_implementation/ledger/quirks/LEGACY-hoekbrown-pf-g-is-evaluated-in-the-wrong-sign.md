---
wp: LEGACY
title: "HoekBrown_PF::g is evaluated in the WRONG SIGN FRAME and is a Tresca potential"
legacy_seq: 410
---
### `HoekBrown_PF::g` is evaluated in the WRONG SIGN FRAME and is a Tresca potential

`HoekBrown_YF` negates to the geomechanics frame (`sigma_geo = -sigma`) before
calling `principalStresses()`. `HoekBrown_PF::g` does NOT, then destructures the
ascending tuple as `[sigma3, sigma2, sigma1]` and feeds `sigma3` — the tree's most
COMPRESSIVE principal, not the geo-frame minor — into
`arg = mb_psi*sigma3/sigma_ci + s`. On any compressive state that `arg` is negative,
so `g` always takes its `else` branch, `sigma1 - sigma3 - sigma_ci*s`, which is a
**Tresca** potential. Measured by the ADR-97 P0 oracle (`adr97_oracle/cppm_hb.py`):

* the header's own central-difference `dg/dsigma` at `[-2000,-6000,-25000]` is
  `[1,0,-1,0,0,0]` for `HB_mb_psi = mb`, `mb/2` **and** `0` — the parameter has
  **no effect at all** and the flow is exactly non-dilatant (trace 0);
* it is **32.86 deg** off the frame-consistent normal at `mb_psi == mb`, i.e. exactly
  where the deck is asking for ASSOCIATED flow;
* at the apex all six of its flow directions have negative trace, so the hydrostatic
  direction is not in its return cone and a trial pushed past the tensile corner has
  **no return to the apex at all** — this is the mechanism behind the residual
  recorded in `tests/test_adr94_hlist_hb.py` ("the drive still fails to converge on
  the step that would push strain past the tensile corner").

Fixing `g` changes `Backward_Euler`, which ADR-97 D1 keeps byte-identical, so
ADR-97 P3 did NOT fix it: `Closest_Point` builds a frame-consistent Hoek-Brown
potential in its own code path, and the difference is pinned in BOTH directions by
`tests/test_adr97_p3_hoekbrown.py::test_gate4_cp_and_be_disagree_by_the_measured_potential_gap`
(plastic volumetric strain **+2.208e-05** under the intended potential vs
**+2.1e-13** under the shipped one; a 5.9e-02 relative stress gap, 26.2 % of the
strength scale). That test turns RED the day `g` is fixed, which is its purpose.
