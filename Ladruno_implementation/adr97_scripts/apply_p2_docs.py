"""ADR-97 P2 (wp/97c) -- idempotent doc/ledger staging.

Appends the P2 rows to the three build-control ledgers and the P2 entry to the
ADR's implementation log.  Re-runnable: every block is keyed on a marker string.

    python3.12 Ladruno_implementation/adr97_scripts/apply_p2_docs.py
"""
import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMPL = os.path.join(ROOT, "Ladruno_implementation")

PR = "[#824](https://github.com/nmorabowen/OpenSees/pull/824)"

VANILLA = [
    ("`cp_principal_return()` (spectral ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: the PRINCIPAL-STRESS-SPACE multi-surface "
     "closest-point return for the Mohr-Coulomb family (Clausen, Damkilde & "
     "Andersen 2006/2007). New private members `cp_principal_return()` (spectral "
     "decomposition of the trial stress, Clausen boundary-plane region test, "
     "closed-form face / edge / edge / apex returns, the Koiter tangent per "
     "region, and the 6D eigenprojection back-transform incl. the rotation term "
     "and its l'Hopital limit) and `cp_apply_tangent_policy()`; a dispatch block "
     "in `Closest_Point` that (a) for MohrCoulombTensionCutoff offers the trial "
     "to ADR-84's `special_return` FIRST -- mapping `SR_QUALITY_EXACT`'s raw "
     "Koiter tangent to `Algorithmic` and refusing `SR_QUALITY_FALLBACK` under "
     "`strict_convergence`, exactly as `Backward_Euler` does -- and (b) otherwise "
     "runs the principal return, so the Mohr-Coulomb family NEVER reaches the "
     "smooth 6D Newton or the Euclidean apex test; a compile-time "
     "`asdp_all_ivs_are_inert` fold and the `ladruno_cp_principal_family` "
     "constant (non-zero only when the YF's and PF's family markers are EQUAL and "
     "non-zero AND every IV is inert), with `ladruno_cp_supported` widened to "
     "`smooth || principal`; `+#include <Eigen/Eigenvalues>`. Refuses loudly, with "
     "`LADRUNO_MATERIAL_REFUSED`, on `MC_phi == 0` (Tresca: the apex is at "
     "infinity and Clausen's boundary planes pass through it), on a non-isotropic "
     "elastic tangent, on a singular face/edge system, on NaN, and on a returned "
     "state outside the Mohr-Coulomb cone (checked EXACTLY in principal space "
     "against all three surfaces) or outside the header's composite `f` (checked "
     "at a tolerance RELATIVE to `max(strength_scale, |sigma_ret|)` -- "
     "`yf_tolerance()` alone is the ADR-94 M5 unit trap for a stress reassembled "
     "from a spectral decomposition). `Backward_Euler` UNTOUCHED (ADR-97 D1), "
     "re-pinned byte-identical on the same 23 decks. | " + PR + " |"),
    ("additive `yf_cp_principal_family` trait ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctionBase.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: additive `yf_cp_principal_family` trait "
     "(`integral_constant<int,0>` by default; 1 = Mohr-Coulomb, 2 = MC + Rankine "
     "cutoff) + the `CP_PRINCIPAL_MC_FACE_PARAMS` macro declaring "
     "`cp_mc_face_params(ivs, params, sin_phi, k_coh)` and its base default, which "
     "returns FALSE so a yield function that is not a Mohr-Coulomb surface refuses "
     "the step rather than returning to a fabricated one. | " + PR + " |"),
    ("the twin `pf_cp_principal_family` trait ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowBase.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: the twin `pf_cp_principal_family` trait + the "
     "`CP_PRINCIPAL_MC_FLOW_PARAMS` macro (`cp_mc_flow_params(ivs, params, "
     "sin_phi, sin_psi)`) with a false base default. The material takes the "
     "principal path only when the two family markers are EQUAL and non-zero, "
     "which is what keeps the generator's MIXED pairings (MohrCoulomb_YF x "
     "VonMises_PF, VonMises_YF x MohrCoulomb_PF, ...) refused. | " + PR + " |"),
    ("additive `hardening_policy_is_inert` trait ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/HardeningFunction.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: additive `hardening_policy_is_inert` trait "
     "(keyed on the policy) and its `hardening_is_inert<HardeningFunction<EVT,"
     "Policy>>` lift -- h identically zero, so the internal variable never moves. "
     "The principal-space return is a projection onto a FIXED surface and carries "
     "no `q`-row, so it is only valid for a perfectly plastic specialization; "
     "folding this at compile time is what keeps a hypothetical "
     "MohrCoulomb_YF<ArmstrongFrederick..> out of that path. | " + PR + " |"),
    ("`hardening_policy_is_inert` specialized to ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/AllASDHardeningFunctions.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `hardening_policy_is_inert` specialized to "
     "`true_type` for `NullHardeningScalarPolicy` and `NullHardeningTensorPolicy` "
     "only. Nothing else in the file qualifies -- both Linear laws and "
     "ArmstrongFrederick have a live `h`. | " + PR + " |"),
    ("hardening_is_perfectly_plastic()`, the per-IV accessor ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/InternalVariableType.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `static constexpr "
     "hardening_is_perfectly_plastic()`, the per-IV accessor the material folds "
     "over the whole tuple. | " + PR + " |"),
    ("`cp_mc_face_params` (returns `sin(phi)` and ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/MohrCoulomb_YF.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `cp_mc_face_params` (returns `sin(phi)` and "
     "`k = c cos(phi)` read from the SAME MC_phi / MC_c this functor's own `f` "
     "uses, so the closest-point map returns to exactly the surface `f` measures) "
     "+ the `yf_cp_principal_family = 1` specialization. Nothing existing is "
     "touched; in particular the EUCLIDEAN `check_apex_region` stays as it is and "
     "`Closest_Point` never calls it. | " + PR + " |"),
    ("`cp_mc_flow_params` (`sin(phi)`, `sin(psi)`) ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/MohrCoulomb_PF.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `cp_mc_flow_params` (`sin(phi)`, `sin(psi)`) "
     "+ the `pf_cp_principal_family = 1` specialization. The closest-point map "
     "never reads `c` from the flow potential, so this functor's `MC_c * M_PI/180` "
     "units wart cannot reach it. | " + PR + " |"),
    ("`cp_mc_face_params` (the Mohr-Coulomb HALF of ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/MohrCoulombTensionCutoff_YF.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `cp_mc_face_params` (the Mohr-Coulomb HALF of "
     "the composite -- reached only after `special_return` has declined) + the "
     "`yf_cp_principal_family = 2` specialization. `special_return` itself is "
     "UNCHANGED and is reused verbatim by `Closest_Point`. | " + PR + " |"),
    ("`cp_mc_flow_params` + the ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/MohrCoulombTensionCutoff_PF.h` | "
     "`// Ladruno (ADR-97 wp/97c)`: `cp_mc_flow_params` + the "
     "`pf_cp_principal_family = 2` specialization. | " + PR + " |"),
    ("the parse-time `Closest_Point` refusal message ",
     "| `SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp` | "
     "`// Ladruno (ADR-97 wp/97c)`: the parse-time `Closest_Point` refusal message "
     "updated for P2 -- it now names the MIXED pairings explicitly (MohrCoulomb_YF "
     "with VonMises_PF / DruckerPrager_PF / HoekBrown_PF, VonMises_YF / "
     "DruckerPrager_YF / HoekBrown_YF with MohrCoulomb_PF) and says that no oracle "
     "covers them, instead of listing MohrCoulomb wholesale as 'P2'. The mechanism "
     "(`instance->supportsClosestPoint()`) is unchanged. | " + PR + " |"),
]

IMPLEMENTATIONS = (
    "ADR-97 P2 - `Closest_Point` for the Mohr-Coulomb family (principal-space ",
    "| **ADR-97 P2 - `Closest_Point` for the Mohr-Coulomb family (principal-space "
    "multi-surface CPPM + Koiter tangent)** ([[97_ladruno_asdp_closest_point_adr]], "
    "report [[reviews/adr97_p2_report]], mutation [[reviews/adr97_p2_mutation]]) - "
    "the shipped 6D Mohr-Coulomb yield function is the exact Lode-angle form, whose "
    "gradient carries a `1/cos(3θ)` the header dodges with a `|θ| >= 29°` "
    "Drucker-Prager substitution and, by default, a central difference over the six "
    "raw Voigt slots; a coupled Newton on that gradient cannot converge "
    "quadratically at a corner, which is where MC models live. In PRINCIPAL STRESS "
    "space the same surface is six PLANES, so every return - face, either corner "
    "LINE, or the vertex - is a closed-form linear projection in the elastic metric "
    "(no Newton at all; `cp_iterations` reads 1 on a plastic step, 0 on an elastic "
    "one) and the Koiter tangent is constant per region. Region selection is "
    "Clausen's BOUNDARY-PLANE test, never a `dΛ >= 0` active-set search (at the apex "
    "with ψ<φ the three multipliers are not all positive); the plane signs are taken "
    "ANALYTICALLY (`sgn_L = sign(n_L·ℓ_other)`), verified against the P0 oracle's "
    "calibration over 60 (φ,ψ,c) combinations and 13094 trial states with 0 "
    "mismatches. The 6D back-transform is `C = Rs T Rs⁻¹ E` with the eigenprojection "
    "ROTATION term `(y_i−y_j)/(x_i−x_j)` and its l'Hôpital limit on a degenerate "
    "trial eigenvalue, thresholded RELATIVE to `strength_scale()` (ADR-94 M5). "
    "**MohrCoulombTensionCutoff reuses ADR-84's `special_return` VERBATIM** (cutoff "
    "face / Rankine edge / MC∩TC corner / compound corner / apex), mapping "
    "`SR_QUALITY_EXACT`'s raw Koiter `stiffness_return` to `Algorithmic` and refusing "
    "`SR_QUALITY_FALLBACK` under `strict_convergence`; only when the hook declines "
    "does the plain-MC principal return take over, and the COMPOSITE `f` is "
    "re-checked. **Measured:** all six P0 oracle trial states (face, sheared face, "
    "both corner lines, two apex states) reproduced to **1e-16** relative with "
    "committed `|f| <= 5.9e-14`; step-size independent (N = 1/4/10/40 all 1e-16); "
    "`Algorithmic` vs a central difference of the binary's own assembled internal "
    "force **2.88e-11** (degenerate edge) and **8.4e-9 / 1.1e-8** (face, on a new "
    "free-node rig - neither `fd_tangent_driver` rig can reach a face state), where "
    "the shipped `Backward_Euler` tangents do not converge on the same rig at all; "
    "global Newton **16 vs 62** iterations (MC oedometric) and **22 vs 74** (the "
    "ADR-84 MCTC deck); MCTC hydrostatic-tension and Rankine-face returns "
    "BIT-IDENTICAL to `Backward_Euler` (gap 0.0), which is the check that the "
    "ADR-84 geometry really is reused rather than re-derived. **Support 20 → 22 of "
    "46**, NOT the plan's optimistic 31: `MohrCoulomb_YF` appears in 7 "
    "specializations and `MohrCoulomb_PF` in 6, and in only ONE are both of the "
    "family; the six MIXED pairings are covered by no oracle and stay refused at "
    "parse time, enforced by matching compile-time family markers plus an "
    "all-IVs-inert fold. `Backward_Euler` byte-identical (D1). **Finding, recorded "
    "not fixed:** `Backward_Euler` reproduces the exact return only through its OWN "
    "finite difference (`MC_ds > 0`, error 2.1e-14) - the shipped ANALYTIC "
    "Lode-angle `c1/c2/c3` branch (`MC_ds = 0`, which every deck in this repo uses) "
    "is 2.9e-1 away and step-size dependent. | constitutive integrator + consistent "
    "tangent (member functions on the existing template) | - (no class tag) | "
    "`SRC/material/nD/ASDPlasticMaterial3D/*` (11 files, see [[LEDGER_vanilla_files]]); "
    "`tests/test_adr97_p2_principal.py` (38 tests); "
    "`Ladruno_implementation/adr97_scripts/apply_p2_{cpp,core,fixes,docs}.py`, "
    "`mutate_p2.py`; `reviews/adr97_p2_report.md` | **P2 shipped (opt-in, not the "
    "default)**; P3 HoekBrown, P4 `Numerical_Algorithmic_*` re-point, P5 explicit "
    "gate + StiffSoil, P6 the default-flip measurement, P7 close-out + banner still "
    "open | " + PR + " |")

QUIRKS = ("## ASDPlasticMaterial3D — Mohr-Coulomb principal-space return (ADR-97 P2)",
          """
## ASDPlasticMaterial3D — Mohr-Coulomb principal-space return (ADR-97 P2)

**`Backward_Euler` reproduces the exact Mohr-Coulomb return ONLY through the
header's own finite difference.** `MohrCoulomb_YF::df_dsigma_ij` and
`MohrCoulomb_PF` both branch on `MC_ds`: `> 0` central-differences their own
`f` / `g` over the six Voigt slots, `== 0` uses an ANALYTIC Lode-angle
`c1/c2/c3` expression. On a proportional face-return deck (E 30000, nu 0.25,
phi 30, psi 10, c 10; trial `sigma = [-10,-40,-100]`) `Backward_Euler` lands on
the exact closest point to **2.1e-14** with `MC_ds = 1e-4`, 1.3e-12 with 1e-6,
1.0e-10 with 1e-8 — and **2.9e-1 away, step-size dependent**, with `MC_ds = 0`.
Mohr-Coulomb's flow direction is CONSTANT inside a sextant (the surface is
piecewise linear), so the cutting plane and the closest point are provably the
same point there; `f` is exactly linear in principal stress, so a central
difference of it is the EXACT 6D gradient and its accuracy *improves* with a
larger step. Two independent references therefore agree against the analytic
coefficients. **Every Mohr-Coulomb deck in `tests/` passes `MC_ds 0.0`**, i.e.
runs on the analytic branch; the ADR-84 MCTC battery is largely insulated
because `special_return` does not use that gradient. Pinned in both directions
by `tests/test_adr97_p2_principal.py::
test_gate4_backward_euler_agrees_only_through_its_own_finite_difference`.
Fixing it changes `Backward_Euler`, so ADR-97 D1 defers it to its own PR.

**An oedometric Mohr-Coulomb deck at `nu = 0.25` with `phi = 30` NEVER YIELDS.**
The elastic lateral-stress ratio `K0 = nu/(1-nu) = 1/3` coincides EXACTLY with
the Mohr-Coulomb compression meridian `(1-sin phi)/(1+sin phi) = 1/3`, so the
uniaxial-strain stress path runs PARALLEL to the yield surface and `f` is
identically `-c cos(phi)` at every stress level, however hard you press. A
tangent or iteration gate built on that rig is silently vacuous: it measures the
ELASTIC operator and passes. Assert the state is plastic before measuring
anything, and pick `nu < 0.25` for `phi = 30` (the P2 gates use 0.15). The same
coincidence exists for any `nu = (1-sin phi)/2`.

**`yf_tolerance()` is not a usable admissibility tolerance for a stress
reassembled from a spectral decomposition.** It defaults to the ABSOLUTE
`f_absolute_tol = 1e-6` (`f_relative_tol` defaults to 0, ADR-94 M5). Recomputing
the header's `f` from `Q diag(y) Q^T` on the ADR-84 MCTC deck (kPa,
`|sigma| ~ 5.4e3`) gives 3.4e-6 — 6e-10 RELATIVE, i.e. round-off — and a check
written against `yf_tolerance()` refuses the step. The round-off is amplified
because an EDGE return lands exactly on a corner, where the Lode angle is ill
conditioned (`dtheta/dJ3 ~ 1/cos(3 theta)` diverges and `dA/dtheta` is not
stationary). Check admissibility where the return was COMPUTED — in principal
space, against all three surfaces of the sextant, which is exact and perfectly
conditioned — and keep the invariant-form check only as a loose guard, at a
tolerance relative to `max(strength_scale, |sigma|)`.

**A load-driven cube rig cannot reach a Mohr-Coulomb FACE state.** Both rigs in
`adr97_oracle/fd_tangent_driver.py` land on `s1 == s2`: `uniaxial` fixes x and y
on every node, so the state is axisymmetric by construction, and `full` (top face
free) passes its limit point under any lateral load before it yields in a
three-distinct-principal state. Measuring a tangent in the face region needs a
kinematically over-determined rig — ADR-97 P2 uses a FREE-NODE rig: the
homogeneous strain field prescribed on seven of the eight nodes, node 7 left
free. Seven prescribed nodes mean no limit point at any stress level, so any
state can be reached, and the 3x3 assembled block at the free node is compared
with a central difference of its own reaction (the same two-rig scheme, since
`setNodeDisp` does not trigger `Domain::update`).
""")

QUIRKS2 = ("**A `template class` explicit instantiation is the WRONG shape for a",
           """
**A `template class` explicit instantiation is the WRONG shape for a
`g++ -fsyntax-only` pre-flight of `ASDPlasticMaterial3D`.** It instantiates
EVERY member of the specialization, including members the real build never
touches because their only call site sits under an `if constexpr` -- e.g.
`cp_apex_return`, which calls `yf.apex_stress()` and therefore fails to compile
for any yield function without an apex (`VonMises_YF`,
`MohrCoulombTensionCutoff_YF`). The result is a page of errors about code that is
correct and unreachable. Instantiate the MEMBERS instead:

    #define private public
    #include ".../AllASDPlasticMaterial3Ds.h"
    typedef ASDPlasticMaterial3D<LinearIsotropic3D_EL, ...> MCMC_t;
    template int MCMC_t::Closest_Point(const VoigtVector&);
    static_assert(MCMC_t::supportsClosestPoint(), "...");

`if constexpr` then discards the unreachable branches exactly as it does in the
real build, and the `static_assert`s turn the support matrix itself into a
compile-time gate (P2 pins all six mixed pairings that way, so a widened family
trait fails at pre-flight instead of at run time). The include set comes from
`adr97_scripts/mk_incs.py` and the build tree's `build/build/Release/build.ninja`
-- note the doubled `build/build`, which `mk_incs.py`'s own usage line does not
say.
"""
)

LOG = ("**2026-09-07 — P2 (`wp/97c-cp-principal`",
       """
- **2026-09-07 — P2 (`wp/97c-cp-principal`, PR [#824](https://github.com/nmorabowen/OpenSees/pull/824), build `1072c27ae`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for the **Mohr-Coulomb family** through a PRINCIPAL-STRESS-SPACE multi-surface closest-point return (Clausen, Damkilde & Andersen 2006/2007): on the sorted sextant the surface is one PLANE, two corner LINES and a vertex, so every return is a closed-form linear projection in the elastic metric — no Newton at all, and `cp_iterations` reads 1 on a plastic step and 0 on an elastic one. Region selection is Clausen's BOUNDARY-PLANE test, whose signs are taken analytically (`sgn_L = sign(n_L·ℓ_other)`) rather than from the oracle's dimensional reference point — verified identical over 60 (φ,ψ,c) combinations and 13094 trial states, **0 mismatches**; the edge directions come from a cross product instead of an SVD (2.2e-16) and `Rs⁻¹` is built as the Voigt image of `Qᵀ·Q` instead of a numerical inverse (2.4e-15). **MohrCoulombTensionCutoff reuses ADR-84's `special_return` verbatim** and only falls back to the plain-MC return when that hook declines, re-checking the COMPOSITE `f`. **Measured:** gate 1a — all six P0 oracle trial states (face, sheared face, both corner lines, two apex states) to **1.5e-16 / 9.1e-16 / 0.0 / 8.0e-16 / 0.0 / 0.0** relative with committed `|f| <= 5.9e-14`; gate 1b — worst committed `f` **1.2e-14 / 1.8e-15 / 2.0e-14** on triaxial, simple-shear and rotating-principal-direction paths (the rotation checked at **20.70°**); gate 1c — MCTC hydrostatic tension `24.7·I` and the Rankine face return both **BIT-IDENTICAL to `Backward_Euler` (gap 0.0)**, which is the check that ADR-84's geometry is reused rather than re-derived, and the confined-compression fall-through admissible on BOTH branches (`f_MC` 2.8e-14 on a scale of 94, `f_TC` −217); gate 2 — `Algorithmic` vs a central difference of the binary's own assembled internal force **2.88e-11** (degenerate edge, l'Hôpital branch) and **8.40e-9 / 1.08e-8** (face, separations 0.19 / 0.27) on a NEW free-node rig, with the shipped `Backward_Euler`/`Continuum` and `/Secant` failing to converge on that same rig (`analyze -> -3`), and global-Newton **16 vs 62 / 50** (MC oedometric) and **22 vs 74** (ADR-84 MCTC); gate 4 — the 23-deck / 282-row `Backward_Euler` baseline still byte-identical, CP step-size independent (N = 1/4/10/40 all 1e-16) and CP ≡ BE at the apex to **2.1e-16**; gate 6 — 12 refusal tests, including all six MIXED YF/PF pairings, `MC_phi == 0`, and three hydrostatic-plus-vanishing-deviator trials that commit a finite vertex instead of NaN. **Support 20 → 22 of 46, not the plan's 31:** `MohrCoulomb_YF` is registered in 7 specializations and `MohrCoulomb_PF` in 6, and in only ONE are both of the family; the principal map assumes BOTH the surface and the potential are piecewise linear and the smooth 6D map cannot use MC's Lode-angle gradient, so the six mixed pairings stay refused, enforced by matching compile-time family markers plus an all-IVs-inert fold. **Found while implementing:** (1) `Backward_Euler` reproduces the exact return ONLY through its own finite difference — `MC_ds = 1e-4` gives **2.1e-14**, `MC_ds = 0` (the shipped ANALYTIC Lode-angle branch, which every deck in this repo uses) gives **2.9e-1** and is step-size dependent; pinned in both directions, NOT fixed, because it changes `Backward_Euler` (D1); (2) an oedometric MC deck at `nu = 0.25`, `phi = 30` NEVER yields, because `K0 = nu/(1-nu) = 1/3` coincides exactly with the compression meridian — a silently vacuous tangent gate; (3) `yf_tolerance()` is not a usable admissibility tolerance for a stress reassembled from a spectral decomposition (round-off 6e-10 relative refused a valid MCTC step, amplified because an edge return lands on a Lode-angle corner) — replaced by an EXACT principal-space check plus a stress-relative composite-`f` guard; (4) neither `fd_tangent_driver` rig can reach a Mohr-Coulomb FACE state. All four are in [[LEDGER_quirks]]. **gate 5 (mutation, [[reviews/adr97_p2_mutation]]).** Full report: [[reviews/adr97_p2_report]].""")


def append_block(path, marker, block):
    """Idempotent append.  `marker` MUST be a literal substring of `block`, so a
    second run detects its own previous output."""
    txt = io.open(path, encoding="utf-8", newline="").read()
    assert marker in block, marker
    if marker in txt:
        return False
    crlf = "\r\n" in txt
    b = block.replace("\n", "\r\n") if crlf else block
    if not txt.endswith("\r\n" if crlf else "\n"):
        txt += "\r\n" if crlf else "\n"
    io.open(path, "w", encoding="utf-8", newline="").write(txt + b.lstrip("\r\n"))
    return True


def main():
    n = 0
    vpath = os.path.join(IMPL, "LEDGER_vanilla_files.md")
    for marker, row in VANILLA:
        if append_block(vpath, marker, row + "\n"):
            n += 1
    ipath = os.path.join(IMPL, "LEDGER_implementations.md")
    if append_block(ipath, IMPLEMENTATIONS[0], IMPLEMENTATIONS[1] + "\n"):
        n += 1
    qpath = os.path.join(IMPL, "LEDGER_quirks.md")
    if append_block(qpath, QUIRKS[0], QUIRKS[1]):
        n += 1
    if append_block(qpath, QUIRKS2[0], QUIRKS2[1]):
        n += 1

    # the ADR implementation log: insert before "## See also"
    apath = os.path.join(IMPL, "97_ladruno_asdp_closest_point_adr.md")
    txt = io.open(apath, encoding="utf-8", newline="").read()
    if LOG[0] not in txt:
        crlf = "\r\n" in txt
        blk = LOG[1].replace("\n", "\r\n") if crlf else LOG[1]
        anchor = ("\r\n## See also" if crlf else "\n## See also")
        assert txt.count(anchor) == 1, txt.count(anchor)
        txt = txt.replace(anchor, blk + anchor, 1)
        io.open(apath, "w", encoding="utf-8", newline="").write(txt)
        n += 1
    print("blocks written: %d" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
