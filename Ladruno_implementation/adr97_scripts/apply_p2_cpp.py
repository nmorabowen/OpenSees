"""ADR-97 P2 (wp/97c) -- re-runnable, idempotent C++ staging script.

Adds the PRINCIPAL-STRESS-SPACE multi-surface closest-point return map (Clausen,
Damkilde & Andersen 2006/2007) and its Koiter consistent tangent to
`integration_method Closest_Point` / `tangent_type Algorithmic`, for the
Mohr-Coulomb family:

  * `MohrCoulomb_YF` x `MohrCoulomb_PF`
  * `MohrCoulombTensionCutoff_YF` x `MohrCoulombTensionCutoff_PF`, which first
    calls ADR-84's `special_return` hook (cutoff face / Rankine edge / MC-cutoff
    corner / compound corner / apex, all closed form, with the raw Koiter tangent
    in `stiffness_return`) and only falls back to the plain-MC principal return
    when that hook declines -- no ADR-84 geometry is re-derived here.

`Backward_Euler` and every YF/PF code path it executes are UNTOUCHED (ADR-97 D1).

Run from the worktree root:
    python3.12 Ladruno_implementation/adr97_scripts/apply_p2_cpp.py
Every edit is an exact-string anchor replacement; a missing or already-applied
anchor is reported and, for a missing one, aborts.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

EDITS = []


def edit(relpath, anchor, new, tag):
    EDITS.append((relpath, anchor, new, tag))


# ===========================================================================
# 1. HardeningFunction.h -- an "inert hardening" trait (h == 0 identically)
# ===========================================================================
edit(
    "HardeningFunction.h",
    """// Opt-in trait, keyed on the POLICY (specialized in AllASDHardeningFunctions.h).
template <typename Policy>
struct hardening_policy_has_cp_derivatives : std::false_type {};
""",
    """// Opt-in trait, keyed on the POLICY (specialized in AllASDHardeningFunctions.h).
template <typename Policy>
struct hardening_policy_has_cp_derivatives : std::false_type {};

// Ladruno (ADR-97 wp/97c): is this hardening policy INERT -- h identically zero,
// so the internal variable never moves?  The principal-stress-space Mohr-Coulomb
// return map added in P2 is a closed-form projection onto a FIXED surface: it has
// no q-row, so it is only valid for a perfectly plastic specialization.  Folded
// over the IV tuple at compile time, it is what keeps a hypothetical
// MohrCoulomb_YF<ArmstrongFrederick...> out of the principal path and refused at
// parse time (ADR-97 P2) instead of silently returning to the wrong surface.
template <typename Policy>
struct hardening_policy_is_inert : std::false_type {};
""",
    "hardening_policy_is_inert trait",
)

edit(
    "HardeningFunction.h",
    """struct hardening_has_cp_derivatives<HardeningFunction<EVT, Policy>>
    : hardening_policy_has_cp_derivatives<Policy> {};""",
    """struct hardening_has_cp_derivatives<HardeningFunction<EVT, Policy>>
    : hardening_policy_has_cp_derivatives<Policy> {};

// Ladruno (ADR-97 wp/97c): the same lift, for the inert-hardening trait.
template <typename T>
struct hardening_is_inert : std::false_type {};

template <typename EVT, typename Policy>
struct hardening_is_inert<HardeningFunction<EVT, Policy>>
    : hardening_policy_is_inert<Policy> {};""",
    "hardening_is_inert lift",
)

# ===========================================================================
# 2. AllASDHardeningFunctions.h -- the two Null policies are inert
# ===========================================================================
edit(
    "AllASDHardeningFunctions.h",
    """template <> struct hardening_policy_has_cp_derivatives<ArmstrongFrederickPolicy>       : std::true_type {};""",
    """template <> struct hardening_policy_has_cp_derivatives<ArmstrongFrederickPolicy>       : std::true_type {};

// Ladruno (ADR-97 wp/97c): the two Null policies return h == 0 identically, so
// the internal variable never moves.  Nothing else in this file qualifies: both
// Linear laws and ArmstrongFrederick have a live h.
template <> struct hardening_policy_is_inert<NullHardeningScalarPolicy> : std::true_type {};
template <> struct hardening_policy_is_inert<NullHardeningTensorPolicy> : std::true_type {};""",
    "Null policies inert",
)

# ===========================================================================
# 3. InternalVariableType.h -- per-IV accessor for the inert trait
# ===========================================================================
edit(
    "InternalVariableType.h",
    """    static constexpr bool hardening_supports_cp()
    {
        return hardening_has_cp_derivatives<HardeningType>::value;
    }""",
    """    static constexpr bool hardening_supports_cp()
    {
        return hardening_has_cp_derivatives<HardeningType>::value;
    }

    // Ladruno (ADR-97 wp/97c): does this internal variable's hardening law leave
    // it FIXED (h == 0)?  Required by the principal-space Mohr-Coulomb return,
    // which projects onto a surface it assumes does not move.
    static constexpr bool hardening_is_perfectly_plastic()
    {
        return hardening_is_inert<HardeningType>::value;
    }""",
    "IV hardening_is_perfectly_plastic",
)

# ===========================================================================
# 4. YieldFunctionBase.h -- principal-family marker + the MC face parameters
# ===========================================================================
edit(
    "YieldFunctionBase.h",
    """template <typename T>
struct yf_has_cp_derivatives : std::false_type {};
""",
    """template <typename T>
struct yf_has_cp_derivatives : std::false_type {};

// Ladruno (ADR-97 wp/97c): which PRINCIPAL-STRESS-SPACE closest-point family this
// yield function belongs to (0 = none, 1 = Mohr-Coulomb, 2 = Mohr-Coulomb with a
// Rankine tension cutoff).  `pf_cp_principal_family` in PlasticFlowBase.h is the
// twin; ASDPlasticMaterial3D takes the principal path ONLY when the two markers
// are equal and non-zero.  That equality test is the point of the design: the
// generator registers cross pairings such as MohrCoulomb_YF x VonMises_PF and
// VonMises_YF x MohrCoulomb_PF, for which NEITHER map is verified -- the
// principal return assumes both the surface AND the flow potential are the
// piecewise-linear Mohr-Coulomb ones, and the smooth 6D map of P1 cannot use
// MohrCoulomb's Lode-angle gradient (it swaps in a Drucker-Prager substitution
// for |theta| >= 29 deg and otherwise central-differences f).  Those pairings
// stay REFUSED at parse time, naming ADR-97 P2.
template <typename T>
struct yf_cp_principal_family : std::integral_constant<int, 0> {};

// Ladruno (ADR-97 wp/97c): the Mohr-Coulomb face constants this yield function
// contributes to the principal-space return:
//     f = a . s - k_coh,   a_ij = [ (1+sin phi)/2 , 0 , -(1-sin phi)/2 ]
//     k_coh = c cos(phi),  apex = (k_coh / sin phi) * [1,1,1]
// with s the principal stresses sorted DESCENDING and tension positive.  The P0
// oracle (`adr97_oracle/cppm_mc.py`) verifies numerically, over 2000 random
// states, that this IS the header's own invariant expression
// A(theta) sqrt(J2) + I1 sin(phi)/3 - c cos(phi) to 1e-14 relative.  Returning
// false means "not a Mohr-Coulomb surface" and is the base default.
#define CP_PRINCIPAL_MC_FACE_PARAMS template <typename IVStorageType, typename ParameterStorageType> \\
    bool cp_mc_face_params(const IVStorageType& internal_variables_storage, \\
        const ParameterStorageType& parameters_storage, \\
        double& sin_phi, double& k_coh) const
""",
    "yf_cp_principal_family + face-params macro",
)

edit(
    "YieldFunctionBase.h",
    """    YIELD_FUNCTION_IV_DERIVATIVE
    {
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        const int n = iv.size();
        for (int i = 0; i < n; ++i) out[i] = 0.0;
    }""",
    """    YIELD_FUNCTION_IV_DERIVATIVE
    {
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        const int n = iv.size();
        for (int i = 0; i < n; ++i) out[i] = 0.0;
    }

    // Ladruno (ADR-97 wp/97c): default -- this yield function is not a
    // Mohr-Coulomb surface.  Unreachable in practice (the material only calls it
    // when yf_cp_principal_family != 0), but a false here refuses the step
    // loudly rather than returning to a fabricated surface.
    CP_PRINCIPAL_MC_FACE_PARAMS
    {
        (void) internal_variables_storage;
        (void) parameters_storage;
        sin_phi = 0.0;
        k_coh   = 0.0;
        return false;
    }""",
    "YF base cp_mc_face_params default",
)

# ===========================================================================
# 5. PlasticFlowBase.h -- principal-family marker + the MC dilatancy parameter
# ===========================================================================
edit(
    "PlasticFlowBase.h",
    """template <typename T>
struct pf_has_cp_derivatives : std::false_type {};
""",
    """template <typename T>
struct pf_has_cp_derivatives : std::false_type {};

// Ladruno (ADR-97 wp/97c): twin of `yf_cp_principal_family` (YieldFunctionBase.h).
// The material takes the principal-stress-space closest-point path only when the
// yield function's marker and this one are EQUAL and non-zero.
template <typename T>
struct pf_cp_principal_family : std::integral_constant<int, 0> {};

// Ladruno (ADR-97 wp/97c): the dilatancy this flow direction contributes to the
// principal-space return.  The header's `m` is
//     m = deviator( dg/dsigma evaluated with PHI ) + sin(psi)/3 * delta
// -- the DEVIATORIC shape of the phi-surface plus a psi-controlled volumetric
// part, NOT the textbook non-associated gradient (which would use psi in the
// deviatoric shape too).  In principal space that is exactly
//     m_ij = a_ij - (sin phi)/3 * [1,1,1] + (sin psi)/3 * [1,1,1]
// which reduces to a_ij when psi == phi, as it must.  `sin_phi` is read back for
// a consistency check against the yield function's own (they share ONE MC_phi
// parameter object -- utuple_storage de-duplicates parameters by type).
#define CP_PRINCIPAL_MC_FLOW_PARAMS template <typename StorageType, typename ParameterStorageType> \\
    bool cp_mc_flow_params(const StorageType& internal_variables_storage, \\
        const ParameterStorageType& parameters_storage, \\
        double& sin_phi, double& sin_psi) const
""",
    "pf_cp_principal_family + flow-params macro",
)

edit(
    "PlasticFlowBase.h",
    """    PLASTIC_FLOW_IV_DERIVATIVE
    {
        (void) iv;
        (void) depsilon;
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        out.setZero();
    }""",
    """    PLASTIC_FLOW_IV_DERIVATIVE
    {
        (void) iv;
        (void) depsilon;
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        out.setZero();
    }

    // Ladruno (ADR-97 wp/97c): default -- this flow direction is not a
    // Mohr-Coulomb potential.
    CP_PRINCIPAL_MC_FLOW_PARAMS
    {
        (void) internal_variables_storage;
        (void) parameters_storage;
        sin_phi = 0.0;
        sin_psi = 0.0;
        return false;
    }""",
    "PF base cp_mc_flow_params default",
)

# ===========================================================================
# 6. MohrCoulomb_YF.h -- family 1
# ===========================================================================
edit(
    "YieldFunctions/MohrCoulomb_YF.h",
    """    APEX_STRESS
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double p_apex = c / tan(phi);""",
    """    // Ladruno (ADR-97 wp/97c): the principal-space face constants.  Both are
    // read from the SAME MC_phi / MC_c parameters this functor's own `f` uses,
    // so the closest-point map returns to exactly the surface `f` measures.
    CP_PRINCIPAL_MC_FACE_PARAMS
    {
        (void) internal_variables_storage;
        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double c   = GET_PARAMETER_VALUE(MC_c);
        sin_phi = std::sin(phi);
        k_coh   = c * std::cos(phi);
        return true;
    }

    APEX_STRESS
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double p_apex = c / tan(phi);""",
    "MC_YF cp_mc_face_params",
)

edit(
    "YieldFunctions/MohrCoulomb_YF.h",
    """//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<MohrCoulomb_YF<NO_HARDENING>> : std::true_type {};""",
    """//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<MohrCoulomb_YF<NO_HARDENING>> : std::true_type {};

// Ladruno (ADR-97 wp/97c): principal-stress-space closest-point family 1
// (plain Mohr-Coulomb).  Paired only with MohrCoulomb_PF, which carries the
// same marker.
template<class NO_HARDENING>
struct yf_cp_principal_family<MohrCoulomb_YF<NO_HARDENING>>
    : std::integral_constant<int, 1> {};""",
    "MC_YF family marker",
)

# ===========================================================================
# 7. MohrCoulomb_PF.h -- family 1
# ===========================================================================
edit(
    "PlasticFlowDirections/MohrCoulomb_PF.h",
    """    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi>;""",
    """    // Ladruno (ADR-97 wp/97c): the principal-space dilatancy constants.  NOTE
    // that this functor's `c` is scaled by M_PI/180 a few lines above (a units
    // wart recorded in LEDGER_quirks -- benign there because c is additive in g
    // and dies under differentiation); the closest-point map never reads c from
    // the flow potential, only sin(psi), so the wart cannot reach it.
    CP_PRINCIPAL_MC_FLOW_PARAMS
    {
        (void) internal_variables_storage;
        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double psi = GET_PARAMETER_VALUE(MC_psi)*M_PI/180;
        sin_phi = std::sin(phi);
        sin_psi = std::sin(psi);
        return true;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi>;""",
    "MC_PF cp_mc_flow_params",
)

edit(
    "PlasticFlowDirections/MohrCoulomb_PF.h",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

#endif""",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97c): principal-stress-space closest-point family 1.
template<class NO_HARDENING>
struct pf_cp_principal_family<MohrCoulomb_PF<NO_HARDENING>>
    : std::integral_constant<int, 1> {};

#endif""",
    "MC_PF family marker",
)

# ===========================================================================
# 8. MohrCoulombTensionCutoff_YF.h -- family 2
# ===========================================================================
edit(
    "YieldFunctions/MohrCoulombTensionCutoff_YF.h",
    """    using internal_variables_t = std::tuple<NO_HARDENING>;

    // Ladruno (ADR-94 wp/94c, M5): composite max(f_MC, f_TC).""",
    """    // Ladruno (ADR-97 wp/97c): the Mohr-Coulomb HALF of the composite, for the
    // principal-space closest-point return.  The map only reaches it after
    // `special_return` above has declined -- i.e. when the cutoff is inactive at
    // the trial, or when the trial is MC-dominant and no exact cutoff feature
    // validated (Stage 3c) -- so the surface being returned to really is the MC
    // one.  The integrator re-checks the COMPOSITE f at the returned state and
    // refuses loudly if the cutoff turns out to be violated.
    CP_PRINCIPAL_MC_FACE_PARAMS
    {
        (void) internal_variables_storage;
        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double c   = GET_PARAMETER_VALUE(MC_c);
        sin_phi = std::sin(phi);
        k_coh   = c * std::cos(phi);
        return true;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    // Ladruno (ADR-94 wp/94c, M5): composite max(f_MC, f_TC).""",
    "MCTC_YF cp_mc_face_params",
)

# ===========================================================================
# 9. MohrCoulombTensionCutoff_PF.h -- family 2
# ===========================================================================
edit(
    "PlasticFlowDirections/MohrCoulombTensionCutoff_PF.h",
    """    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi, TC_min_stress>;""",
    """    // Ladruno (ADR-97 wp/97c): principal-space dilatancy constants (see
    // MohrCoulomb_PF for the derivation; this functor's MC branch is verbatim
    // that one).
    CP_PRINCIPAL_MC_FLOW_PARAMS
    {
        (void) internal_variables_storage;
        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double psi = GET_PARAMETER_VALUE(MC_psi)*M_PI/180;
        sin_phi = std::sin(phi);
        sin_psi = std::sin(psi);
        return true;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi, TC_min_stress>;""",
    "MCTC_PF cp_mc_flow_params",
)

edit(
    "PlasticFlowDirections/MohrCoulombTensionCutoff_PF.h",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

#endif""",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97c): principal-stress-space closest-point family 2
// (Mohr-Coulomb + Rankine tension cutoff).
template<class NO_HARDENING>
struct pf_cp_principal_family<MohrCoulombTensionCutoff_PF<NO_HARDENING>>
    : std::integral_constant<int, 2> {};

#endif""",
    "MCTC_PF family marker",
)

# the MCTC YF marker has to live where the class is complete
edit(
    "YieldFunctions/MohrCoulombTensionCutoff_YF.h",
    """template<class NO_HARDENING>
struct yf_has_special_return<MohrCoulombTensionCutoff_YF<NO_HARDENING>> : std::true_type {};""",
    """template<class NO_HARDENING>
struct yf_has_special_return<MohrCoulombTensionCutoff_YF<NO_HARDENING>> : std::true_type {};

// Ladruno (ADR-97 wp/97c): principal-stress-space closest-point family 2.
template<class NO_HARDENING>
struct yf_cp_principal_family<MohrCoulombTensionCutoff_YF<NO_HARDENING>>
    : std::integral_constant<int, 2> {};""",
    "MCTC_YF family marker",
)


def main():
    applied, skipped = [], []
    for relpath, anchor, new, tag in EDITS:
        path = os.path.join(SRC, relpath)
        if not os.path.isfile(path):
            print("MISSING FILE: %s" % path)
            return 2
        with open(path, "r", encoding="utf-8", newline="") as fh:
            txt = fh.read()
        # the sources in this tree are CRLF; the LF-written anchors are
        # converted to match the file before any matching happens (P1 lesson).
        crlf = "\r\n" in txt
        an = anchor.replace("\n", "\r\n") if crlf else anchor
        nw = new.replace("\n", "\r\n") if crlf else new
        if nw in txt:
            skipped.append(tag)
            continue
        n = txt.count(an)
        if n != 1:
            print("ANCHOR MISS (%d matches) in %s for %r" % (n, relpath, tag))
            return 3
        txt = txt.replace(an, nw, 1)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(txt)
        applied.append(tag)
    print("applied  : %d" % len(applied))
    for t in applied:
        print("   + %s" % t)
    print("already  : %d" % len(skipped))
    for t in skipped:
        print("   = %s" % t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
