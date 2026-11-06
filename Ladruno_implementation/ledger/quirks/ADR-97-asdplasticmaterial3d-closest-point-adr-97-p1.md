---
wp: ADR-97
title: "ASDPlasticMaterial3D — Closest_Point (ADR-97 P1)"
legacy_seq: 408
---
## ASDPlasticMaterial3D — `Closest_Point` (ADR-97 P1)

**Two integrators, two apex answers for the same yield function — on purpose.**
`DruckerPrager_YF::check_apex_region` is a EUCLIDEAN normal-cone test
(`p − p_apex >= eta*q`) and says so in its own comment: the exact condition is in
the ELASTIC metric, `p_tr − p_apex >= (K*etabar/G) * q_tr`, and the yield
function's signature cannot see K, G or the dilatancy. `Backward_Euler` cannot
fix that; `Closest_Point` can, because classification happens in the integrator
where `E` is in scope — so `Closest_Point` does its own region test and **never
calls `check_apex_region`**. The ADR-97 P0 oracle quantified the cost of the
Euclidean test in BOTH directions: at `etabar = 0.2` (exact slope 0.333 < the
header's 0.4) a trial at `(p−p_apex)/q = 0.36` is classified CONE and the cone
return then gives `sqrt(J2)_{n+1} = −0.47`, an inadmissible negative deviatoric
norm; at `etabar = eta = 0.4` (exact slope 0.667) trials at 0.45 and 0.60 are
classified APEX although the correct return is to the cone. Expect the two
integrators to disagree in that band, and do not "fix" one to match the other.

**A sign-based region test degenerates on the hydrostatic axis.** ADR-97's
elastic-metric apex test is `dot(dev_ret, dev_tr) < 0` on the LINEARISED cone
step. For a trial state sitting exactly on the hydrostatic axis, `dev_tr == 0`
and that reads `0 < 0` — i.e. CONE — after which the cone Newton has no flow
direction at all and exhausts its iterations. That degenerate state is not
exotic: it is ADR-94 B4's hydrostatic-tension reproducer, the deck that used to
commit NaN. Any deviator-direction test needs an explicit
`||dev_tr|| <= tol` short circuit, and `tol` must be the YIELD tolerance so the
comparison stays unit consistent (ADR-94 M5).

**A `Path` time series returns ZERO outside its defined range — including at the
last time point of a multi-step run.** A prescribed-strain driver that ends
exactly on the final `-time` entry unloads the whole path to zero in ONE step,
and a plasticity material then reports a perfectly plausible ON-SURFACE stress
that is simply the wrong point on the surface (it took us one confused debugging
round to see it, because the yield residual was 1e-15 the whole way). Always
give a `Path` series one time point past the end of the analysis.

**`utuple_concat_unique_type` de-duplicates internal variables by TYPE, not by
name.** `VonMises_YF<BackStress<TensorLinearHardening>, …>` paired with
`VonMises_PF<BackStress<NullHardeningTensor>>` gives the storage **two**
`BackStress` entries — the yield function reads one and the flow direction reads
the other, and they evolve independently. A deck that wants one shared back
stress must name the SAME hardening law on both sides (the ADR-97 gate-1
Armstrong–Frederick deck does). This is why ADR-97's `df_dq` / `dm_dq` select on
`std::is_same<IVType, AlphaHardeningType>` inside the YF/PF rather than on the
variable's name: each functor differentiates with respect to the variable IT
reads, and returns zero for the other one — which is the correct derivative.

**`tangent_type Algorithmic` existed upstream with no dispatch case anywhere.**
It was dead, which is the only reason nobody was silently getting it. ADR-97 D2
gives it exactly one meaning — the consistent tangent of the `Closest_Point`
map — and REFUSES it with every other integrator: a consistent tangent is
defined only relative to a specific committed map, and offering it on the
cutting-plane `Backward_Euler` would ship a fourth almost-right tangent, which
is the class of defect ADR-94 M3 found.

**`C_alg` is unsymmetric whenever `m != n`.** Non-associated Drucker–Prager
(`etabar != eta`), Mohr–Coulomb (`psi != phi`), Hoek–Brown (`mb_psi != mb`).
Models running `tangent_type Algorithmic` need `system UmfPack`; `ProfileSPD` is
wrong, PARDISO's symmetric `-matrixType` (ADR-75 P1d) must not be selected, and
`FullGeneral` crashes a fully prescribed material-point rig (`FullGenLinSOE`
N = 0).
