"""ADR-97 WP-b (P1) -- re-runnable ledger / ADR-log updates.

Idempotent: each block is skipped when its marker is already present.

usage: python3.12 apply_p1_docs.py [<worktree root>]
"""
import io
import os
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\asdplastic-review-plan-62585c"
LI = os.path.join(ROOT, "Ladruno_implementation")

VANILLA_ROWS = """| `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` | `// Ladruno (ADR-97 wp/97b)`: the CLOSEST-POINT return map and its consistent tangent. New private members `Closest_Point(const VoigtVector&)`, `cp_assemble()` (residual + Jacobian at one point), `cp_apex_region()` (ELASTIC-metric apex classification -- `Closest_Point` deliberately does NOT call the yield function's Euclidean `check_apex_region`, see [[LEDGER_quirks]]), `cp_apex_return()` (the reduced apex system) and `cp_n_iv()`; a `Closest_Point` case in the `setTrialStrainIncr` dispatch switch; `Closest_Point` added to `set_constitutive_integration_method`'s accept list; a compile-time `ladruno_cp_supported` flag (YF trait AND PF trait AND a fold over every IV's hardening law) plus the public `supportsClosestPoint()` the parser refuses on; `ASDP_CP_MAXN` = 26 with the `cp_matrix_t`/`cp_vector_t` fixed-max Eigen typedefs (no heap traffic in the Gauss-point loop); a per-instance `cp_last_iterations` (NSDMI) exposed as `getCPIterations()` / `setResponse` token `cp_iterations` / `getResponse` id 9. `Backward_Euler` is UNTOUCHED (ADR-97 D1) and its histories on 23 decks are pinned in `Ladruno_implementation/adr97_oracle/baselines/`. `Algorithmic` reuses the ONE Jacobian factorization already in hand at convergence; every failure path (non-convergence, NaN, singular Jacobian, negative converged multiplier, over-cap system size) returns `LADRUNO_MATERIAL_REFUSED` and every commit path passes through the existing `ladruno_strict_rejects` gate. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3DGlobals.h` | `// Ladruno (ADR-97 wp/97b)`: `+1` enum value `Closest_Point` on `ASDPlasticMaterial3D_Constitutive_Integration_Method` (appended last, so no existing value's integer changes). `ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic` already existed upstream with no dispatch case anywhere -- ADR-97 D2 gives it one and refuses it with every other integrator. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctionBase.h` | `// Ladruno (ADR-97 wp/97b)`: additive `yf_has_cp_derivatives` trait (defaults false) + the `YIELD_FUNCTION_IV_DERIVATIVE` macro declaring `df_dq(iv, sigma, ivs, params, double* out)` -- the UNCONTRACTED df/dq the closest-point Jacobian's `J_fq` block needs (the shipped `YIELD_FUNCTION_HARDENING` only ever exposes the scalar contraction `-(df/dq).h`). Base default writes zeros, so an unconverted YF compiles unchanged and is refused at parse time. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowBase.h` | `// Ladruno (ADR-97 wp/97b)`: additive `pf_has_cp_derivatives` trait + the `PLASTIC_FLOW_STRESS_DERIVATIVE` (`dm_dsigma`, 6x6) and `PLASTIC_FLOW_IV_DERIVATIVE` (`dm_dq`, 6 x iv.size()) macros, a protected per-instance `dm_dsigma_buffer`, a ZERO default for `dm_dq` and a ONE-TIME-WARNED CENTRAL DIFFERENCE default for `dm_dsigma` (no family ADR-97 P1 ships reaches it; it exists so a future PF runs loudly while its analytic block is written). `dl * dm/dsigma` is the entire content of ADR-94 M3: setting `dl = 0` turns the algorithmic elastic modulus `Xi = (E^-1 + dl dm/ds)^-1` back into `E` and recovers the shipped `Continuum` operator exactly. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/HardeningFunction.h` | `// Ladruno (ADR-97 wp/97b)`: `+#include <type_traits>`; the `HARDENING_FUNCTION_IV_DERIVATIVE` (`dh_dq`) and `HARDENING_FUNCTION_M_DERIVATIVE` (`dh_dm`) macros; `hardening_policy_has_cp_derivatives` (keyed on the policy) and `hardening_has_cp_derivatives` (keyed on the `HardeningFunction<EVT,Policy>` an internal variable actually stores); and forwarders on `HardeningFunction` that zero the block and call the policy only under `if constexpr` on the trait. `dh/dsigma` is deliberately NOT declared -- no policy in this tree reads `sigma`. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/InternalVariableType.h` | `// Ladruno (ADR-97 wp/97b)`: `+#include "HardeningFunction.h"` (the file had no includes at all) and three members mirroring the existing `hardening_function`: `hardening_dh_dq`, `hardening_dh_dm` (both evaluated at the TRIAL value, i.e. at `q_{n+1}` -- this is what makes the ADR-97 D4 hardening update implicit) and the `static constexpr hardening_supports_cp()` the material folds over the whole IV tuple. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/ElasticityBase.h` | `// Ladruno (ADR-97 wp/97b, D6)`: `+#include <type_traits>`; the `el_is_stress_dependent` trait (false by default -- `LinearIsotropic3D_EL`, the only elasticity registered in all 43 non-StiffSoil specializations, is stress independent, so `dE/dsigma` is identically zero) and the `ELASTICITY_STRESS_DERIVATIVE` macro (`dE_dsigma_contract`, contracted with `m` on the way out so no rank-3 object is ever formed) with a zero default. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/AllASDHardeningFunctions.h` | `// Ladruno (ADR-97 wp/97b)`: analytic `dh_dq`/`dh_dm` for the five policies P1 solves implicitly -- `LinearHardeningForTensor` (dh/dm = H*I_dev), `LinearHardeningForScalar` (dh/dm = H*(2/3)*W_e m / sqrt(2/3<m,m>_e)), both `NullHardening*` (zero), and `ArmstrongFrederick` (dh/dalpha = -c_r |dev m|_eq I_dev; dh/dm = h_a I_dev - c_r dev(alpha) (x) (2/3) W_e dev(m)/|dev m|_eq) -- plus the five `hardening_policy_has_cp_derivatives` specializations. AF is differentiated as the branch `f` ACTUALLY TAKES, saturation branch included: `f` is shared with `Backward_Euler` (D1), so a CP-only `f` would be a second hardening law with the same name; the P0 oracle mirrors the same branch. Every formula was cross-checked against a central difference before the C++ was written (8.8e-10 worst). | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/VonMises_YF.h` | `// Ladruno (ADR-97 wp/97b)`: `df_dq` (df/dk = -SQRT_2_over_3, df/dalpha = -(W_s r)/||r||_s, i.e. the shipped `yf_hardening`'s operand, in the same VOIGT convention) + the `yf_has_cp_derivatives` specialization. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/DruckerPrager_YF.h` | `// Ladruno (ADR-97 wp/97b)`: `df_dq` (df/dalpha = -(r/sqrt(J2)) with the NORMAL slots halved; **df/d(cohesion IV) = ZERO**) + the `yf_has_cp_derivatives` specialization. The zero is ADR-97 P0 header finding 2: this YF's cohesion internal variable is commented out of its own `f` (line 26) while `yf_hardening` still contributes `df/dk = -1` times its rate, so a Drucker-Prager with cohesion hardening is PERFECTLY PLASTIC under `Closest_Point` and hardening under `Backward_Euler`. Pinned by `tests/test_adr97_p1_smooth.py::test_gate1_dp_cohesion_hardening_is_perfectly_plastic_under_cp`; fixing the YF changes BE, which D1 forbids here. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/VonMises_PF.h` | `// Ladruno (ADR-97 wp/97b)`: analytic `dm_dsigma = (diag(W_s) I_dev - m (x) m)/||r||_s` and `dm_dq` (= `-(diag(W_s) - m (x) m)/||r||_s` for the back stress, zero otherwise) + the `pf_has_cp_derivatives` specialization. Also `this->`-qualifies `dm_dsigma_buffer`: it lives in the DEPENDENT base `PlasticFlowBase<T>`, so an unqualified reference is ill-formed on GCC and accepted by MSVC (the trap that cost the ADR-94 wave a red CI). | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/DruckerPrager_PF.h` | `// Ladruno (ADR-97 wp/97b)`: analytic `dm_dsigma = diag(W_s) I_dev/(2q) - md (x) md/q` (`md` the deviatoric part of the shipped flow direction; the `etabar*p` term has zero Hessian) and the matching `dm_dq`, both returning ZERO at `q -> 0` because that is the apex, which `Closest_Point` classifies and returns on its own reduced system. Same `this->` qualification. + the `pf_has_cp_derivatives` specialization. | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/ElasticityModels/{StiffSoil_EL,DuncanChang_EL}.h` | `// Ladruno (ADR-97 wp/97b, D6)`: `el_is_stress_dependent` specialized to `true_type`. It drives a ONE-TIME `opserr` warning that `tangent_type Algorithmic` is missing the `E,sigma : (sigma - sigma_tr)` term; the committed stress is unaffected because `Closest_Point` evaluates `E(sigma_{n+1})` inside the residual regardless. ADR-97 P5 supplies the block. Neither model is reachable under `Closest_Point` today (StiffSoil's YF/PF are refused at parse time, `DuncanChang_EL` is commented out of the generator). | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
| `SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp` | `// Ladruno (ADR-97 wp/97b)`: `integration_method Closest_Point` (refused, naming the YF/PF/IV and the ADR phase, when `instance->supportsClosestPoint()` is false -- MohrCoulomb/MCTC are P2, HoekBrown P3, StiffSoil P5) and `tangent_type Algorithmic`; the D2 cross-refusal after the option loop (`Algorithmic` with any integrator other than `Closest_Point` is rejected naming both tokens); both `Valid values:` lists and `print_usage()` updated. Unknown tokens are still rejected (the ADR-94 wp/94a contract). | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
"""

IMPL_ROW = """| **ADR-97 P1 — `integration_method Closest_Point` + `tangent_type Algorithmic`** ([[97_ladruno_asdp_closest_point_adr]], report [[reviews/adr97_p1_report]]) — the fully implicit CLOSEST-POINT return map for `ASDPlasticMaterial3D`, solved as ONE coupled Newton in `x = (σ_{n+1}, q_{n+1}, dλ)` with `E(σ_{n+1})` inside the residual (D6), plus the EXACT consistent tangent of that map taken from the converged Jacobian. Scope: the smooth families — VonMises and Drucker–Prager including the apex — with Null / LinearScalar / LinearTensor / **ArmstrongFrederick** hardening solved implicitly at n+1: **20 of the 46 registered specializations**. Answers ADR-94 M3: no shipped tangent was the tangent of the committed map (`Continuum` 57 %, `Secant` — the DEFAULT — 80 %, `Elastic` 103 %, `Numerical_Algorithmic` 31 % against a central difference of the material's own response), and `Backward_Euler` is an Ortiz–Simo CUTTING PLANE whose fixed point is `σ_tr − Σ_k dλ_k E m(σ^k)`. **Measured:** `Algorithmic` vs a free-DOF finite difference of the binary's own assembled residual **1.51e-11** (4 DOF) / **2.15e-10** (12-DOF sheared), against the pinned `Backward_Euler`/`Continuum` negative control 0.573447 reproduced to the digit; the ADR-94 two-cube model costs **4, 4, 4, 4** Newton iterations against `Continuum`'s **6, 41, 36, 30** (7.1× fewer); committed stress vs the numpy CPPM oracles 7.4e-14 … 9.2e-9 over 9 von Mises and 3 Drucker–Prager cases; local Newton 1 iteration (perfect / linear), 3 (Armstrong–Frederick); Armstrong–Frederick step-refinement error 4.07× smaller than the cutting plane's (the P0 oracle predicted 4.3×), while with LINEAR hardening the two maps agree to 1.45e-16 — which is exactly why ADR-94 H6 could not see the defect. `Backward_Euler` is **byte-identical** (D1), pinned on 23 decks / 282 committed-stress rows in fresh subprocesses. Six new interface members, each with an INERT base default plus an opt-in trait, so an unconverted family compiles unchanged and is REFUSED at parse time rather than silently approximated. `Algorithmic` is refused with any other integrator (D2 — a consistent tangent is defined only relative to a specific committed map). **Not** the default (a separate PR, gated on the P6 measurement). | constitutive integrator + consistent tangent (member functions on the existing template; oracles + `zone_a` gates) | — (no class tag: ASDP is one templated class, no new class, no new file, no CMake target) | `SRC/material/nD/ASDPlasticMaterial3D/*` (14 files, see [[LEDGER_vanilla_files]]); `Ladruno_implementation/97_ladruno_asdp_closest_point_adr.md`, `adr97_oracle/` (5 numpy oracles + `reference_output.txt`), `adr97_oracle/baselines/` (the D1 byte-identity baseline + its dumper), `adr97_scripts/` (the re-runnable edit scripts); `tests/test_adr97_p1_smooth.py`, `tests/test_adr97_p4_inertness.py`, `tests/test_adr97_p6_failloud.py` | **P1 shipped (opt-in, not the default)**; P2 MohrCoulomb/MCTC, P3 HoekBrown, P4 `Numerical_Algorithmic_*` re-point, P5 explicit-integrator gate + StiffSoil, P6 the default-flip measurement, P7 close-out + banner line still open. No banner line yet (P7). | [#819](https://github.com/nmorabowen/OpenSees/pull/819) |
"""

QUIRKS = """
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
"""

LOG_ENTRY = """- **2026-09-07 — P1 (`wp/97b-cp-smooth`, PR [#819](https://github.com/nmorabowen/OpenSees/pull/819), build `dd05b60a6` then `ec6091c4f`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for VonMises and Drucker–Prager (flank + apex) with Null / LinearScalar / LinearTensor / ArmstrongFrederick hardening — 20 of 46 specializations. Six new interface members with inert base defaults + opt-in traits; the parser refuses `Closest_Point` for unconverted families and `Algorithmic` for every other integrator (D2). **Measured:** gate 1 — committed stress vs the P0 oracles 7.4e-14 … 2.4e-12 (VM perfect/linear, 3 paths), 1.3e-13 / 9.9e-13 / **9.2e-9** (VM+AF triaxial / simple shear / rotating normal; the 9.2e-9 is the yield tolerance accumulated over 20 steps, `f_absolute_tol` 1e-6 on a ~45 stress), DP cone associated 4.9e-13 and non-associated 5.1e-13, DP apex exact; worst committed `|f|` over a whole path 7.1e-15; local Newton **1** iteration (perfect, linear) and **3** (AF), gate `<= 5`. gate 2 — `Algorithmic` vs a central difference of the binary's own assembled residual **1.51e-11** (4 free DOFs) and **2.15e-10** (12-DOF sheared rig) against the pinned `Backward_Euler`/`Continuum` **0.573447** reproduced to the digit; ADR-94 two-cube `testIter` **4,4,4,4** vs `Continuum`'s **6,41,36,30** (7.1× fewer; the ADR's `<= 3` was an estimate, 4 is what an exact tangent costs on that rig). gate 3 — AF step-refinement error at N=10 against each map's own N=160 limit: CP 5.68e-3 vs BE 2.31e-2, **4.07×** (P0 predicted 4.3×); with LINEAR hardening the two maps agree to **1.45e-16**, the reason ADR-94 H6 was blind to this. gate 4 — 23 decks / 282 rows byte-identical in fresh subprocesses; CP≡BE at 0.0 / 1.6e-33 / 3.8e-17 on non-rotating perfectly plastic decks and 2.3e-2 apart on AF. gate 6 — 13/13 loud. **Two P0 findings honoured, not fixed** (both would change `Backward_Euler`, which D1 forbids): DP's cohesion IV is out of its own `f`, so DP-with-cohesion-hardening is perfectly plastic under CP and hardening under BE (pinned as a test); and AF's saturation branch is KEPT, with CP differentiating the branch `f` actually takes, because `f` is shared with BE. **Found while testing:** the elastic-metric apex test is a SIGN test on `dot(dev_ret, dev_tr)` and degenerates to "CONE" for a trial state exactly on the hydrostatic axis — ADR-94 B4's own reproducer; fixed with a `||dev_tr|| <= tol_f` short circuit. Also `dm_dsigma_buffer` needed `this->` (dependent base, GCC-only error) and a `Path` time series returns 0 at its last time point.
"""


def append_once(path, marker, text, before=None):
    src = io.open(path, encoding="utf-8", errors="surrogateescape", newline="").read()
    crlf = "\r\n" in src
    mk = marker.replace("\n", "\r\n") if crlf else marker
    if mk in src:
        print("  skip     %s" % os.path.basename(path))
        return 0
    body = text.replace("\n", "\r\n") if crlf else text
    if before is None:
        if not src.endswith("\n" if not crlf else "\r\n"):
            src += "\r\n" if crlf else "\n"
        src += body
    else:
        bf = before.replace("\n", "\r\n") if crlf else before
        assert src.count(bf) == 1, "anchor not unique in %s" % path
        src = src.replace(bf, body + bf)
    io.open(path, "w", encoding="utf-8", errors="surrogateescape",
            newline="").write(src)
    print("  appended %s" % os.path.basename(path))
    return 1


def main():
    n = 0
    n += append_once(os.path.join(LI, "LEDGER_vanilla_files.md"),
                     "Ladruno (ADR-97 wp/97b)", VANILLA_ROWS)
    n += append_once(os.path.join(LI, "LEDGER_implementations.md"),
                     "ADR-97 P1 — `integration_method Closest_Point`", IMPL_ROW)
    n += append_once(os.path.join(LI, "LEDGER_quirks.md"),
                     "ASDPlasticMaterial3D — `Closest_Point` (ADR-97 P1)", QUIRKS)
    n += append_once(os.path.join(LI, "97_ladruno_asdp_closest_point_adr.md"),
                     "P1 (`wp/97b-cp-smooth`, PR", LOG_ENTRY,
                     before="## See also")
    print("apply_p1_docs: %d block(s) appended" % n)


if __name__ == "__main__":
    main()
