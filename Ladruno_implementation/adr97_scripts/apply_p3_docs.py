#!/usr/bin/env python3.12
"""ADR-97 wp/97d -- idempotent ledger / ADR-log edits.

Same shape as apply_p2_docs.py: exact-string anchors, LOUD on a miss, and a
marker check so the script can be re-run after a partial failure.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # Ladruno_implementation
PR = "[#825](https://github.com/nmorabowen/OpenSees/pull/825)"


def edit(path, anchor, new, marker):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        txt = fh.read()
    if "\r\n" in txt:
        anchor = anchor.replace("\n", "\r\n")
        new = new.replace("\n", "\r\n")
        marker = marker.replace("\n", "\r\n")
    if marker in txt:
        print("  SKIP (already applied): %s" % os.path.basename(path))
        return
    if txt.count(anchor) != 1:
        raise SystemExit("ANCHOR MISS (%d hits) in %s:\n%s"
                         % (txt.count(anchor), path, anchor[:300]))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(txt.replace(anchor, new, 1))
    print("  applied: %s" % os.path.basename(path))


# ---------------------------------------------------------------------------
# LEDGER_vanilla_files.md -- one row per touched upstream file
# ---------------------------------------------------------------------------
VF = os.path.join(ROOT, "LEDGER_vanilla_files.md")

VF_ANCHOR = ("| `SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp` | "
             "`// Ladruno (ADR-97 wp/97c)`")

ROWS_97D = [
    ("`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`",
     "`// Ladruno (ADR-97 wp/97d)`: the PRINCIPAL-STRESS-SPACE closest-point return for "
     "the **Hoek-Brown** family (Clausen & Damkilde 2008) -- a CURVED surface, so P2's "
     "closed-form projection does not apply and P1's smooth 6D Newton cannot be used "
     "either (`HoekBrown_YF` has no analytic gradient at all: it central-differences the "
     "COMPOSITE `max(f_shear, f_tension)` over the six raw Voigt slots). New private "
     "members: the `hb_consts_t` scratch struct; `hb_arg_floor`, `hb_a_ij`, `hb_m_ij`, "
     "`hb_dm_ii` (the CURVATURE term; `dm/dy` has exactly one non-zero entry), "
     "`hb_f_composite`, `hb_f_floor` (the GRADIENT-scaled yield tolerance -- `|df/dy1| = "
     "1 + a mb arg^(a-1)` diverges at the apex, so an absolute gate is unattainable "
     "there), `hb_layout`, `hb_y_of` (the edge constraint built into the "
     "parametrization), `hb_assemble` (residual + ANALYTIC Jacobian), `hb_start`, "
     "`hb_solve`, `hb_in_apex_cone` (exact elastic-metric cone test by the DUAL/facet "
     "form), `hb_classify` and `cp_hb_return`; a `ladruno_cp_principal_family == 3` "
     "branch in `Closest_Point`'s dispatch, taken BEFORE the `yf_has_apex` block. The "
     "Newton runs in the surface's OWN variable `arg = (w^2)^(1/a)` (with `y1` as the "
     "unknown the map diverges near the apex -- the first step overshoots `arg < 0` and "
     "the next Jacobian is singular) with a NORMALIZED flow direction (`|m| ~ "
     "arg^(a-1)` blows up where the near-apex returns land; 5 iterations normalized "
     "against 6 un-normalized). P2's spectral decomposition, eigenprojection "
     "back-transform, `cp_apply_tangent_policy`, plastic-strain convention and "
     "strict-convergence contract are REUSED, not duplicated. The FRAME-CONSISTENT "
     "Hoek-Brown potential lives here, in `Closest_Point`'s own path: the shipped "
     "`HoekBrown_PF::g` is a Tresca potential (see [[LEDGER_quirks]]) and is left "
     "untouched so `Backward_Euler` stays byte-identical (D1). Refuses loudly, with "
     "`LADRUNO_MATERIAL_REFUSED`, on a YF/PF parameter mismatch, out-of-range "
     "`HB_sigci`/`HB_mb`/`HB_s`/`HB_a`/`HB_mb_psi`, `HB_mb_psi > HB_mb` (the potential's "
     "own `arg` goes negative on part of its own yield surface), a non-isotropic elastic "
     "tangent, a failed eigen-decomposition, a singular converged Jacobian, a broken "
     "principal ordering, NaN, and a return outside the surface (checked EXACTLY in "
     "principal space AND against the header's own composite `f`, both with the "
     "gradient-scaled floor)."),
    ("`SRC/material/nD/ASDPlasticMaterial3D/YieldFunctionBase.h`",
     "`// Ladruno (ADR-97 wp/97d)`: the `CP_PRINCIPAL_HB_FACE_PARAMS` macro declaring "
     "`cp_hb_face_params(ivs, params, sigma_ci, mb, s, a)` and its base default, which "
     "returns FALSE so a yield function that is not a Hoek-Brown surface refuses the "
     "step rather than returning to a fabricated one. Additive; nothing existing is "
     "touched."),
    ("`SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowBase.h`",
     "`// Ladruno (ADR-97 wp/97d)`: the twin `CP_PRINCIPAL_HB_FLOW_PARAMS` macro "
     "(`cp_hb_flow_params(ivs, params, sigma_ci, mb_psi, s, a)`) with a false base "
     "default. `Closest_Point` builds its own frame-consistent potential from these and "
     "does NOT call the functor's `g`; `sigma_ci`/`s`/`a` come back so the material can "
     "cross-check them against the yield function's."),
    ("`SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/HoekBrown_YF.h`",
     "`// Ladruno (ADR-97 wp/97d)`: `cp_hb_face_params` (the four surface constants, "
     "read from the SAME `HB_sigci`/`HB_mb`/`HB_s`/`HB_a` this functor's own `f` uses; "
     "`HB_ds` is deliberately NOT exported, because `Closest_Point` never uses the "
     "central-difference gradient) + the `yf_cp_principal_family = 3` specialization. "
     "Nothing existing is touched -- in particular the EUCLIDEAN `CHECK_APEX_REGION` "
     "stays as it is and `Closest_Point` never calls it (it over-claims the apex "
     "region: 32 of the P0 oracle's 400 scanned trials disagree, and on one of them "
     "`APEX_STRESS` would discard 2.29 % of the strength scale)."),
    ("`SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/HoekBrown_PF.h`",
     "`// Ladruno (ADR-97 wp/97d)`: `cp_hb_flow_params` + the "
     "`pf_cp_principal_family = 3` specialization. `g` and the central-difference "
     "`PLASTIC_FLOW_DIRECTION` are left EXACTLY as shipped (D1), even though `g` is "
     "evaluated in the wrong sign frame and collapses to a Tresca potential -- see "
     "[[LEDGER_quirks]] and [[reviews/adr97_p3_report]] Sec. 3; fixing it changes "
     "`Backward_Euler` and is the owner's separate PR."),
    ("`SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp`",
     "`// Ladruno (ADR-97 wp/97d)`: the `Closest_Point` refusal message now names the "
     "Hoek-Brown family as SHIPPED under the matched-pair rule (P3) instead of pending; "
     "`StiffSoil` remains P5. The refusal itself is unchanged -- it is driven by "
     "`instance->supportsClosestPoint()`, which the compile-time family markers widen."),
]


def vf_rows():
    return "".join("| %s | %s | %s |\n" % (f, w, PR) for f, w in ROWS_97D)


# ---------------------------------------------------------------------------
# LEDGER_implementations.md
# ---------------------------------------------------------------------------
IMP = os.path.join(ROOT, "LEDGER_implementations.md")
IMP_ANCHOR = "| **ADR-97 P2 - `Closest_Point` for the Mohr-Coulomb family"

IMP_ROW = (
    "| **ADR-97 P3 - `Closest_Point` for the Hoek-Brown family (principal-space return "
    "to a CURVED surface + curvature-consistent tangent)** "
    "([[97_ladruno_asdp_closest_point_adr]], report [[reviews/adr97_p3_report]], "
    "mutation [[reviews/adr97_p3_mutation]]) - `HoekBrown_YF` has NO analytic gradient "
    "(it central-differences the COMPOSITE `max(f_shear, f_tension)` over the six raw "
    "Voigt slots, which picks up spurious gradient across the branch switch), and the "
    "Hoek-Brown meridian is a CURVE, so neither P1's smooth 6D Newton nor P2's "
    "closed-form projection applies. P3 keeps P2's principal space and its "
    "eigenprojection back-transform verbatim and replaces the projection with a 4x4 "
    "Newton on the face or either curved edge, plus a closed-form vertex. Three things "
    "the P0 oracle measured as REQUIREMENTS rather than choices: the Newton runs in the "
    "surface's own variable `arg = (w^2)^(1/a)` (with `y1` as the unknown the first step "
    "overshoots `arg < 0` and the next Jacobian is singular), the flow direction is "
    "NORMALIZED in the residual (`|m| ~ arg^(a-1)` blows up exactly where near-apex "
    "returns land; 5 Newton iterations against 6), and the admissibility tolerance is "
    "GRADIENT-scaled (`|df/dy1|` diverges at the apex, so an absolute 1e-10 gate is "
    "unattainable within ~1e-2 kPa of the vertex -- the Hoek-Brown instance of ADR-94's "
    "`f_relative_tol` lesson). Region selection uses no boundary planes (a curved "
    "surface has none): the apex by an EXACT elastic-metric cone test in dual/facet "
    "form, then the face, then the edges ordered by the FACE return's own `y1-y2` / "
    "`y2-y3` margins, which are the exact signed boundary functions. The composite's "
    "tension branch is INERT on the surface (`f_shear <= 0` already forces `y1 <= T`; "
    "0 of 4000 sampled surface points have `f_tension` winning), so the shear/tension "
    "corner IS the apex and no Rankine face return exists. **Measured:** all SEVEN "
    "oracle regions to **1.4e-15 / 1.7e-15 / 7.3e-16 / 1.3e-15 / 0.0 / 0.0 / 1.9e-13** "
    "relative with `cp_iterations` **4/4/3/3/1/1/5** (the ADR's `<= 5` face gate, met); "
    "worst committed `f` **7.3e-11** over triaxial / shear / rotating paths; the "
    "consistent tangent vs a central difference of the binary's own assembled internal "
    "force **3.3e-09 / 3.9e-09** (free-node FACE rig, curvature AND rotation terms live) "
    "and **1.5e-10** (degenerate-eigenvalue edge rig); `Backward_Euler` with BOTH "
    "`Continuum` and `Secant` FAILS TO CONVERGE on the Hoek-Brown oedometric deck where "
    "`Closest_Point`/`Algorithmic` converges in 16 global iterations; the ADR-94 uniaxial "
    "tension plateau lands on the EXACT closed-form limit **244.4854419 kPa** "
    "(1.2e-16 relative) rather than on the apex `T` = 245.0152 that ADR-94's `rel=5e-2` "
    "pin accepts, and a trial 2x/20x past the tensile corner now commits the finite "
    "apex instead of stalling. **Support 22 -> 23 of 46**: `HoekBrown_YF` is registered "
    "in 7 specializations and `HoekBrown_PF` in 6, and in exactly ONE are both of the "
    "family; the other ELEVEN Hoek-Brown pairings stay refused (all eleven pinned by "
    "`static_assert` in the g++ pre-flight, the seven deck-reachable ones also at run "
    "time in BOTH directions). `Backward_Euler` byte-identical (D1), re-pinned on the "
    "same 23 decks. **Finding, recorded not fixed:** `HoekBrown_PF::g` is evaluated in "
    "the un-negated frame and collapses to a TRESCA potential -- `HB_mb_psi` inert, "
    "exactly zero dilatancy, 32.86 deg off the normal at `mb_psi == mb`, and NO return "
    "to the apex past the tensile corner (the mechanism behind the ADR-94 HB residual). "
    "`Closest_Point` uses a frame-consistent potential in its own path and the gap is "
    "PINNED in both directions (plastic `eps_vol` **+2.208e-05** vs BE's **+2.1e-13**, "
    "a 5.9e-02 stress gap = 26.2 %% of the strength scale), so the test turns red the "
    "day the shipped `g` is fixed -- the owner's separate PR. **Also fixed:** "
    "`HB_sigma_ci` is not a parameter (`HB_sigci` is), so P1's and P2's Hoek-Brown "
    "refusal assertions were passing on a MISSING PARAMETER and never exercised the "
    "family gate. | constitutive integrator + consistent tangent (member functions on "
    "the existing template) | - (no class tag) | "
    "`SRC/material/nD/ASDPlasticMaterial3D/*` (6 files, see [[LEDGER_vanilla_files]]); "
    "`tests/test_adr97_p3_hoekbrown.py` (32 tests); "
    "`Ladruno_implementation/adr97_scripts/apply_p3_docs.py`; "
    "`reviews/adr97_p3_report.md` | **P3 shipped (opt-in, not the default)**; P4 "
    "`Numerical_Algorithmic_*` re-point, P5 explicit gate + StiffSoil, P6 the "
    "default-flip measurement, P7 close-out + banner still open | %s |\n" % PR)


# ---------------------------------------------------------------------------
# LEDGER_quirks.md
# ---------------------------------------------------------------------------
QK = os.path.join(ROOT, "LEDGER_quirks.md")
QK_ANCHOR = "## ASDPlasticMaterial3D — Mohr-Coulomb principal-space return (ADR-97 P2)"

QK_BLOCK = """## ASDPlasticMaterial3D — Hoek-Brown (ADR-97 P3)

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

### The Hoek-Brown yield tolerance must be scaled by the GRADIENT, not by sigma_ci

`|df/dy1| = 1 + a*mb*arg^(a-1)` **diverges** at the apex (`arg -> 0`), so a last-ulp
error in the returned `y1` carries `eps*|y1|*|df/dy1|` into `f`: an ABSOLUTE 1e-10
admissibility gate is unattainable within ~1e-2 kPa of the vertex. The P0 oracle
measures `max|f| = 1.08e-10` against its own round-off floor of `1.16e-08` there.
This is the Hoek-Brown instance of ADR-94 M5's `f_relative_tol` lesson, and it is
sharper: scaling by `sigma_ci` (or by `strength_scale`) is not enough on its own,
because the offending factor is the gradient's own conditioning.

### The Newton must run in the surface's OWN variable, not in `y1`

The Hoek-Brown surface exists only for `arg = s - mb*y1/sigma_ci >= 0`. With `y1` as
the Newton unknown the first step from the elastic predictor OVERSHOOTS (measured
`arg = -2.2456e-03` at iteration 1 on the oracle's own near-apex trial) and the next
Jacobian is singular. Substituting `arg = (w^2)^(1/a)` — so
`y1 = T - (sigma_ci/mb)(w^2)^(1/a)` and `f = y1 - y3 - sigma_ci*w^2` — is
polynomial-smooth and feasible for ANY real `w`: no clipping, no line search, no
feasibility guard. Write `(w*w)^(1/a)` rather than `w^(2/a)`: the latter is NaN for a
negative iterate, and the two agree for `w > 0`.

Related: NORMALIZE the flow direction in the residual (multiplier rescaled by `|m|`,
the returned stress is invariant). `|m| ~ arg^(a-1)` blows up exactly where the
near-apex returns land — 460.6 there against 3.9 on an ordinary face point — and the
worst-case Newton count over the oracle's 400-trial scan is **6** un-normalized and
**5** normalized.

### `HB_sigma_ci` is not a parameter; it is `HB_sigci`

ADR-97 P1 and P2 both wrote `HB_sigma_ci` in their Hoek-Brown refusal tests. Under
the ADR-94 contract a missing model parameter is an ERROR, so those decks were
rejected for the wrong reason and the refusal assertions never exercised the family
gate at all. A refusal test that is not ALSO checked in the positive direction (the
same deck must CONSTRUCT under an integrator that does support it) cannot tell the
two apart. Every gate-6 row in `tests/test_adr97_p3_hoekbrown.py` is checked both
ways for exactly this reason.

### `AllASDInternalVariableTypes.h` and `AllASDHardeningFunctions.h` have NO include guard

Including either directly in a translation unit that also includes
`ASDPlasticMaterial3D.h` (which pulls both in) is a redefinition storm. Relevant to
any standalone syntax-check / pre-flight translation unit.

### `#define private public` breaks GCC 15's libstdc++ if it precedes `<sstream>`

`std::basic_stringbuf::__xfer_bufptrs` is declared `private` and re-declared later;
flipping the keyword makes the second declaration disagree with the first, which
GCC 15 reports as a hard `-Wtemplate-body` ERROR (not a warning). The ADR-97 P2
pre-flight idiom still works — include the standard library and Eigen FIRST, then
`#define private public`, then the project's own headers.

"""


# ---------------------------------------------------------------------------
# the ADR's implementation log + the phase table
# ---------------------------------------------------------------------------
ADR = os.path.join(ROOT, "97_ladruno_asdp_closest_point_adr.md")
ADR_ANCHOR = "## See also"

ADR_ENTRY = """- **2026-09-07 — P3 (`wp/97d-cp-hoekbrown`, PR [#825](https://github.com/nmorabowen/OpenSees/pull/825), build `d0de2d7c8`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for the **Hoek-Brown family** as a principal-stress-space return to a **CURVED** surface (Clausen & Damkilde 2008). P2's principal machinery — spectral decomposition, eigenprojection back-transform with the rotation term and its l'Hôpital limit, analytic `Rs⁻¹`, tangent policy, plastic-strain convention, strict-convergence contract — is REUSED verbatim; only the projection is replaced, by a 4×4 Newton on the face or either curved edge plus a closed-form vertex. Three things the P0 oracle measured as REQUIREMENTS: (a) the Newton runs in the surface's own variable `arg = (w²)^(1/a)` — with `y1` as the unknown the first step overshoots to `arg = −2.2e-03` and the next Jacobian is singular; (b) the flow direction is NORMALIZED in the residual, because `|m| ~ arg^(a−1)` blows up exactly where the near-apex returns land (460.6 vs 3.9), which is what holds the worst case at **5** Newton iterations instead of 6; (c) the admissibility tolerance is GRADIENT-scaled, because `|df/dy1| = 1 + a·mb·arg^(a−1)` diverges at the apex and an absolute 1e-10 gate is unattainable within ~1e-2 kPa of the vertex. Region selection uses NO boundary planes (a curved surface has none): the apex by an exact elastic-metric cone test in dual/facet form, then the face, then the edges ordered by the FACE return's own `y1−y2` / `y2−y3` margins, which ARE the exact signed boundary functions. The composite's tension branch is INERT on the surface (`f_shear ≤ 0` already forces `y1 ≤ T`; 0 of 4000 sampled surface points have `f_tension` winning, `min(f_shear − f_tension) = 1.004e-04`), so the shear/tension corner IS the apex and no Rankine face return exists — a Rankine return from the oracle's `[645,265,255]` trial leaves `f_shear = +123.34` kPa, inadmissible. **Measured:** gate 1a — all SEVEN oracle regions to **1.448e-15 / 1.735e-15 / 7.259e-16 / 1.268e-15 / 0.0 / 0.0 / 1.851e-13** relative with `cp_iterations` **4/4/3/3/1/1/5** and committed `|f| ≤ 1.8e-11`; the last row is the trial the header's EUCLIDEAN `CHECK_APEX_REGION` calls an apex, where `APEX_STRESS` would discard **122.6 kPa = 2.29 %** of the strength scale (the header always over-claims: `D3·(octant)` is a strict subset of the octant, and 32 of the oracle's 400 scanned trials disagree). gate 1b — worst committed `f` **0.0 / +7.3e-12 / +7.3e-11** on triaxial / simple-shear / rotating paths (rotation checked at **16.80°**). gate 1c — the ADR-94 uniaxial tension plateau lands on the **exact closed-form limit 244.4854419 kPa (1.16e-16 relative)**, NOT on the apex `T` = 245.0152 that ADR-94's `rel = 5e-2` pin accepts (`T` is 0.216 % above the limit, because at `y1 = T` the clamp leaves `f_shear = y1 > 0`); and a trial 2×/20× past the tensile corner now commits the finite apex instead of stalling, because the frame-consistent potential HAS a return to the vertex where the shipped Tresca `g` has none. gate 2 — `Algorithmic` vs a central difference of the binary's own assembled internal force **3.342e-09 / 3.869e-09** (free-node FACE rig, curvature AND rotation terms live) and **1.485e-10** (degenerate-eigenvalue edge rig); the apex tangent is rank 0 by construction so there is no FD gate on it; the iteration contrast is an OUTCOME — `Backward_Euler` with `Continuum` AND with `Secant` fails to converge (`analyze -> -3`) on the HB oedometric deck where `Closest_Point`/`Algorithmic` converges in **16**. gate 3 — nothing to gate: the family is registered perfectly plastic and the all-IVs-inert requirement is folded into the family marker at compile time. gate 4 — the 23-deck / 282-row `Backward_Euler` baseline still byte-identical; step refinement replaces P2's exact step independence, because a CURVED surface cannot be step independent — N = 1 IS the oracle (1.448e-15) and the sequence converges monotonically (**2.594e-3 → 1.413e-3 → 5.885e-4 → 1.260e-4** against N = 160) to a limit 0.26 % away. gate 6 — all seven deck-reachable MIXED Hoek-Brown pairings refused, each checked in BOTH directions; `strict_convergence` byte-inert (0.0); hydrostatic and near-degenerate tension states commit a finite vertex, not NaN. **Support 22 → 23 of 46:** `HoekBrown_YF` is registered in 7 specializations and `HoekBrown_PF` in 6, and in exactly ONE are both of the family; the other eleven stay refused, all pinned by `static_assert` in the g++ pre-flight. **The P3 decision, implemented and recorded:** `HoekBrown_PF::g` is evaluated in the un-negated frame and collapses to a **Tresca** potential (`HB_mb_psi` inert, exactly zero dilatancy, 32.86° off the normal at `mb_psi == mb`, no apex return past the tensile corner). `Closest_Point` uses a frame-consistent Hoek-Brown potential in its OWN path; the shipped `g` is untouched (D1) and the gap is PINNED in both directions — plastic `eps_vol` **+2.208111e-05** (the oracle's own +2.208e-05) against BE's **+2.09e-13**, a **5.910e-02** relative stress gap = **26.20 %** of the strength scale. Fixing `g` is the owner's separate PR; this WP is its warrant. **Also found:** `HB_sigma_ci` is not a parameter — the name is `HB_sigci` — so P1's and P2's Hoek-Brown refusal assertions were passing on a MISSING PARAMETER and never exercised the family gate; both corrected and inverted here. Full battery **204 passed, 2 skipped (pre-existing), 0 failed**. Full report: [[reviews/adr97_p3_report]].
"""

ADR_PHASE_ANCHOR = ("| **P2** | `wp/97c-cp-principal` |")
ADR_PHASE_NOTE = ("| **P2** | `wp/97c-cp-principal` |")


def main():
    print("ADR-97 wp/97d doc edits")
    with open(VF, "r", encoding="utf-8", newline="") as fh:
        vf_txt = fh.read()
    nl = "\r\n" if "\r\n" in vf_txt else "\n"
    anchor_line = None
    for line in vf_txt.split(nl):
        if line.startswith(VF_ANCHOR):
            anchor_line = line
            break
    if anchor_line is None:
        raise SystemExit("ANCHOR MISS: no wp/97c OPS_All... row in LEDGER_vanilla_files.md")
    edit(VF, anchor_line + "\n", anchor_line + "\n" + vf_rows(),
         "the PRINCIPAL-STRESS-SPACE closest-point return for the **Hoek-Brown**")

    with open(IMP, "r", encoding="utf-8", newline="") as fh:
        imp_txt = fh.read()
    nl = "\r\n" if "\r\n" in imp_txt else "\n"
    p2_line = None
    for line in imp_txt.split(nl):
        if line.startswith(IMP_ANCHOR):
            p2_line = line
            break
    if p2_line is None:
        raise SystemExit("ANCHOR MISS: no ADR-97 P2 row in LEDGER_implementations.md")
    edit(IMP, p2_line + "\n", p2_line + "\n" + IMP_ROW,
         "ADR-97 P3 - `Closest_Point` for the Hoek-Brown family")

    edit(QK, QK_ANCHOR, QK_BLOCK + QK_ANCHOR,
         "## ASDPlasticMaterial3D — Hoek-Brown (ADR-97 P3)")

    edit(ADR, ADR_ANCHOR, ADR_ENTRY + ADR_ANCHOR,
         "P3 (`wp/97d-cp-hoekbrown`, PR [#825]")
    print("done")


if __name__ == "__main__":
    sys.exit(main())
