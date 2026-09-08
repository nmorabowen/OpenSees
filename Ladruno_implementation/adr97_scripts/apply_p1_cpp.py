"""ADR-97 WP-b (P1) -- re-runnable C++ edit script.

Applies the `Closest_Point` closest-point return map + `Algorithmic` consistent
tangent to ASDPlasticMaterial3D for the SMOOTH families (VonMises,
DruckerPrager incl. apex) with Null / LinearScalar / LinearTensor /
ArmstrongFrederick hardening solved implicitly at n+1.

Idempotent: every hunk is skipped when its marker is already present, and the
script FAILS LOUD when an anchor is missing (upstream drifted).

usage:  python3.12 apply_p1_cpp.py [<worktree root>]
"""
import os
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\asdplastic-review-plan-62585c"
ASD = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

EDITS = []          # (relpath, marker, anchor, replacement, last)


def edit(rel, marker, anchor, replacement, last=False):
    EDITS.append((rel, marker, anchor, replacement, last))


# ===========================================================================
# 1. YieldFunctionBase.h -- trait + df/dq (uncontracted) with an inert default
# ===========================================================================
edit(
    "YieldFunctionBase.h",
    "yf_has_cp_derivatives",
    """template <typename T>
struct yf_has_special_return : std::false_type {};
""",
    """template <typename T>
struct yf_has_special_return : std::false_type {};

// Ladruno (ADR-97 wp/97b): opt-in trait for yield functions that supply the
// UNCONTRACTED closest-point Jacobian block df/dq (YIELD_FUNCTION_IV_DERIVATIVE
// below).  `integration_method Closest_Point` is REFUSED at parse time for a YF
// that does not specialize this to true_type -- never silently approximated.
// See ADR-97 D3 (family-by-family) and the base default in YieldFunctionBase.
template <typename T>
struct yf_has_cp_derivatives : std::false_type {};
""")

edit(
    "YieldFunctionBase.h",
    "#define YIELD_FUNCTION_IV_DERIVATIVE",
    """#define GET_INTERNAL_VARIABLE_HARDENING(type) \\""",
    """// Ladruno (ADR-97 wp/97b): df/dq -- the derivative of f with respect to ONE
// internal variable, UNCONTRACTED (the shipped YIELD_FUNCTION_HARDENING only
// ever exposes the scalar contraction -(df/dq).h, which the closest-point
// Jacobian's J_fq block cannot use).  `iv` selects the variable; the YF matches
// it against its own template parameters with `if constexpr` and writes
// `iv.size()` doubles into `out`.  Convention: the derivative with respect to the
// STORED slot (shear slots doubled relative to the tensor derivative), matching
// df_dsigma_ij and the framework's plain-dot contractions (ADR-94 wp/94c B5).
// It must be the TRUE derivative of THIS YF's f: DruckerPrager_YF has its
// cohesion IV commented out of f, so its df/dq for that IV is ZERO under
// Closest_Point even though yf_hardening still contributes -1 times its rate
// (ADR-97 P0 header finding 2 -- recorded, not fixed here: fixing the YF would
// change Backward_Euler, which D1 keeps byte-identical).
#define YIELD_FUNCTION_IV_DERIVATIVE template <typename IVType, typename IVStorageType, typename ParameterStorageType> \\
    void df_dq(const IVType& iv, const VoigtVector& sigma, \\
        const IVStorageType& internal_variables_storage, \\
        const ParameterStorageType& parameters_storage, \\
        double* out) const

#define GET_INTERNAL_VARIABLE_HARDENING(type) \\""")

edit(
    "YieldFunctionBase.h",
    "ADR-97 wp/97b): default -- this YF",
    """    inline const char* getName() const { return static_cast<T*>(this)->NAME; }

};

#endif""",
    """    // Ladruno (ADR-97 wp/97b): default -- this YF declares no closest-point
    // IV derivative.  Inert (zero), and unreachable in practice: the parser
    // refuses `integration_method Closest_Point` unless yf_has_cp_derivatives.
    YIELD_FUNCTION_IV_DERIVATIVE
    {
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        const int n = iv.size();
        for (int i = 0; i < n; ++i) out[i] = 0.0;
    }

    inline const char* getName() const { return static_cast<T*>(this)->NAME; }

};

#endif""")

# ===========================================================================
# 2. PlasticFlowBase.h -- trait + dm/dsigma (FD default) + dm/dq (zero default)
# ===========================================================================
edit(
    "PlasticFlowBase.h",
    "pf_has_cp_derivatives",
    """#define PLASTIC_FLOW_DIRECTION template <typename StorageType, typename ParameterStorageType> \\""",
    """// Ladruno (ADR-97 wp/97b): opt-in trait for plastic-flow directions that supply
// ANALYTIC closest-point Jacobian blocks (dm/dsigma, dm/dq).  A PF without it
// still compiles -- it inherits the one-time-warned central-difference dm/dsigma
// below -- but `integration_method Closest_Point` is refused at parse time for
// it, so no family this ADR has not converted can silently run on a finite
// difference.  See ADR-97 D3.
template <typename T>
struct pf_has_cp_derivatives : std::false_type {};

// Ladruno (ADR-97 wp/97b): dm/dsigma, the 6x6 derivative of the flow direction
// with respect to the STORED stress slots -- the `dl * dm/ds` term of the
// algorithmic elastic modulus Xi = (E^-1 + dl*dm/ds)^-1.  ADR-94 M3 measured
// exactly this term: setting dl = 0 gives Xi = E and recovers the shipped
// `Continuum` operator, whose error against a central difference of the
// material's own committed response was 57 %.
#define PLASTIC_FLOW_STRESS_DERIVATIVE template <typename StorageType, typename ParameterStorageType> \\
    const VoigtMatrix& dm_dsigma( \\
        const VoigtVector &depsilon, \\
        const VoigtVector& sigma, \\
        const StorageType& internal_variables_storage,  \\
        const ParameterStorageType& parameters_storage) const

// Ladruno (ADR-97 wp/97b): dm/dq for ONE internal variable -- a 6 x iv.size()
// block written into the leading columns of `out` (a 6x6 buffer; both internal
// variable kinds in this tree are size 1 or 6).  `out` is zeroed by the caller.
#define PLASTIC_FLOW_IV_DERIVATIVE template <typename IVType, typename StorageType, typename ParameterStorageType> \\
    void dm_dq(const IVType& iv, \\
        const VoigtVector &depsilon, \\
        const VoigtVector& sigma, \\
        const StorageType& internal_variables_storage,  \\
        const ParameterStorageType& parameters_storage, \\
        VoigtMatrix& out) const

#define PLASTIC_FLOW_DIRECTION template <typename StorageType, typename ParameterStorageType> \\""")

edit(
    "PlasticFlowBase.h",
    "ADR-97 wp/97b): default dm/dsigma",
    """    inline const char* getName() const { return static_cast<T*>(this)->NAME; }
};

#endif""",
    """    // Ladruno (ADR-97 wp/97b): default dm/dsigma -- a central difference of this
    // PF's own flow direction, with a ONE-TIME warning naming the class.  No
    // family shipped by ADR-97 P1 reaches it (they are analytic, and an
    // unconverted family is refused at parse time); it exists so a future PF
    // compiles and runs while its analytic block is being written, loudly.
    PLASTIC_FLOW_STRESS_DERIVATIVE
    {
        static bool warned_fd_dm_dsigma = false;
        if (!warned_fd_dm_dsigma)
        {
            opserr << "ASDPlasticMaterial3D (ADR-97) - plastic flow direction '"
                   << static_cast<const T*>(this)->NAME
                   << "' has no analytic dm/dsigma; Closest_Point is using a"
                   << " CENTRAL DIFFERENCE for the J_ss Jacobian block."
                   << " The converged stress is still exact (inexact Newton),"
                   << " but tangent_type Algorithmic is NOT." << endln;
            warned_fd_dm_dsigma = true;
        }
        const T* self = static_cast<const T*>(this);
        double scale = sigma.maxAbs();
        if (!(scale > 0.0)) scale = 1.0;
        const double h = 1e-8 * scale;
        VoigtVector sp, sm, mp, mm;
        for (int j = 0; j < 6; ++j)
        {
            sp = sigma; sm = sigma;
            sp(j) += h;  sm(j) -= h;
            // copy out immediately: operator() returns a reference to the PF's
            // single mutable buffer, which the second call overwrites.
            mp = self->operator()(depsilon, sp, internal_variables_storage, parameters_storage);
            mm = self->operator()(depsilon, sm, internal_variables_storage, parameters_storage);
            for (int i = 0; i < 6; ++i)
                dm_dsigma_buffer(i, j) = (mp(i) - mm(i)) / (2.0 * h);
        }
        return dm_dsigma_buffer;
    }

    // Ladruno (ADR-97 wp/97b): default dm/dq -- zero.  Correct for every flow
    // direction whose `m` does not read an internal variable; a PF that does
    // read one must override this (VonMises_PF and DruckerPrager_PF do, for
    // their back stress).
    PLASTIC_FLOW_IV_DERIVATIVE
    {
        (void) iv;
        (void) depsilon;
        (void) sigma;
        (void) internal_variables_storage;
        (void) parameters_storage;
        out.setZero();
    }

    inline const char* getName() const { return static_cast<T*>(this)->NAME; }

protected:

    // Ladruno (ADR-97 wp/97b): per-instance return buffer for the finite-
    // difference dm/dsigma default (never static -- ADR-94 wp/94b F2).
    mutable VoigtMatrix dm_dsigma_buffer = VoigtMatrix::Zero();
};

#endif""")

# ===========================================================================
# 3. HardeningFunction.h -- traits + dh/dq, dh/dm forwarders
# ===========================================================================
edit(
    "HardeningFunction.h",
    "hardening_policy_has_cp_derivatives",
    """// Function wrapper base class""",
    """// Ladruno (ADR-97 wp/97b): the closest-point Jacobian needs dh/dq and dh/dm
// (see ADR-97 "Residuals and Jacobians": J_qs = -dl*(dh/dm)(dm/ds),
// J_qq = I - dl*(dh/dq + (dh/dm)(dm/dq))).  `dh/ds` is deliberately NOT
// declared: no hardening policy in this tree reads `sigma` (they read
// `current_value`, `m` and `depsilon`), and one that later does must declare it.
//
// The blocks are written into the leading rows/columns of a 6x6 buffer:
//   dh/dq is n_iv x n_iv, dh/dm is n_iv x 6, with n_iv = 1 (scalar IV) or 6.
// A policy without the derivatives keeps a ZERO default and is REFUSED at parse
// time under Closest_Point, so the zero is never silently believed.
#define HARDENING_FUNCTION_IV_DERIVATIVE template <class EVT, class ParameterStorageType> \
    static void dh_dq( \
        const EVT& current_value, \
        const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma, \
        const ParameterStorageType& parameters_storage, \
        VoigtMatrix& out)

#define HARDENING_FUNCTION_M_DERIVATIVE template <class EVT, class ParameterStorageType> \
    static void dh_dm( \
        const EVT& current_value, \
        const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma, \
        const ParameterStorageType& parameters_storage, \
        VoigtMatrix& out)

// Opt-in trait, keyed on the POLICY (specialized in AllASDHardeningFunctions.h).
template <typename Policy>
struct hardening_policy_has_cp_derivatives : std::false_type {};

// Function wrapper base class""")

edit(
    "HardeningFunction.h",
    "hardening_has_cp_derivatives",
    """    static constexpr const char* NAME = HardeningPolicy::NAME;
    using parameters_t = typename HardeningPolicy::parameters_t;
};""",
    """    // Ladruno (ADR-97 wp/97b): forwarders with an inert zero default.
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        out.setZero();
        if constexpr (hardening_policy_has_cp_derivatives<HardeningPolicy>::value)
            HardeningPolicy::dh_dq(current_value, depsilon, m, sigma, parameters_storage, out);
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        out.setZero();
        if constexpr (hardening_policy_has_cp_derivatives<HardeningPolicy>::value)
            HardeningPolicy::dh_dm(current_value, depsilon, m, sigma, parameters_storage, out);
    }

    static constexpr const char* NAME = HardeningPolicy::NAME;
    using parameters_t = typename HardeningPolicy::parameters_t;
};

// Ladruno (ADR-97 wp/97b): the same question asked of a HardeningFunction type
// (which is what an InternalVariableType stores), not of the policy.
template <typename T>
struct hardening_has_cp_derivatives : std::false_type {};

template <typename EVT, class Policy>
struct hardening_has_cp_derivatives<HardeningFunction<EVT, Policy>>
    : hardening_policy_has_cp_derivatives<Policy> {};""")

edit(
    "HardeningFunction.h",
    "#include <type_traits>",
    """#ifndef HardeningFunctionBase_H
#define HardeningFunctionBase_H
""",
    """#ifndef HardeningFunctionBase_H
#define HardeningFunctionBase_H

#include <type_traits>   // Ladruno (ADR-97 wp/97b): std::false_type for the CP traits
""")

# ===========================================================================
# 4. InternalVariableType.h -- per-IV access to the hardening derivatives
# ===========================================================================
edit(
    "InternalVariableType.h",
    "ADR-97 wp/97b): closest-point",
    """    using parameters_t = typename HardeningType::parameters_t;


};""",
    """    // Ladruno (ADR-97 wp/97b): closest-point Jacobian blocks for THIS internal
    // variable, evaluated (like hardening_function above) at the TRIAL value --
    // i.e. at q_{n+1}, which is what makes the ADR-97 hardening update implicit.
    // `out` is a 6x6 buffer; only rows 0..size()-1 (and, for dh_dq, the same
    // number of columns) are meaningful.
    template <class ParameterStorageType>
    void hardening_dh_dq(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters,
                VoigtMatrix& out) const
    {
        HardeningType::dh_dq(trial_value, depsilon, m, sigma, parameters, out);
    }

    template <class ParameterStorageType>
    void hardening_dh_dm(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters,
                VoigtMatrix& out) const
    {
        HardeningType::dh_dm(trial_value, depsilon, m, sigma, parameters, out);
    }

    // Ladruno (ADR-97 wp/97b): does this IV's hardening law have analytic
    // closest-point derivatives?  Folded over the whole IV tuple at compile time
    // by ASDPlasticMaterial3D so the parser can refuse Closest_Point for a
    // specialization carrying an unconverted hardening law.
    static constexpr bool hardening_supports_cp()
    {
        return hardening_has_cp_derivatives<HardeningType>::value;
    }

    using parameters_t = typename HardeningType::parameters_t;


};""")

edit(
    "InternalVariableType.h",
    'ADR-97 wp/97b): the CP hardening',
    """//Base template struct for all internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
struct InternalVariableType {""",
    """// Ladruno (ADR-97 wp/97b): the CP hardening-derivative forwarders and the
// hardening_has_cp_derivatives trait live in HardeningFunction.h.
#include "HardeningFunction.h"

//Base template struct for all internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
struct InternalVariableType {""")

# ===========================================================================
# 5. ElasticityBase.h -- dE/dsigma contraction (D6), inert zero default
# ===========================================================================
edit(
    "ElasticityBase.h",
    "el_is_stress_dependent",
    """#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value
""",
    """#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

// Ladruno (ADR-97 wp/97b, D6): is this elasticity model stress dependent?
// `LinearIsotropic3D_EL` -- the only elasticity registered in all 43
// non-StiffSoil specializations -- is not, so `dE/dsigma` is identically zero
// and the algorithmic tangent is exact without it.  A model that specializes
// this to true_type gets a one-time warning under `tangent_type Algorithmic`
// until it supplies ELASTICITY_STRESS_DERIVATIVE (ADR-97 P5).
template <typename T>
struct el_is_stress_dependent : std::false_type {};

// Ladruno (ADR-97 wp/97b, D6): the J_ss contribution dl * (dE/dsigma : m).
// Contracted with `m` on the way out because only that contraction ever enters
// the Jacobian -- a rank-3 object never has to be formed.
#define ELASTICITY_STRESS_DERIVATIVE template<class ParameterStorageType> \\
    void dE_dsigma_contract(const VoigtVector& stress, \\
        const VoigtVector& m, \\
        const ParameterStorageType& parameters_storage, \\
        VoigtMatrix& out) const
""")

edit(
    "ElasticityBase.h",
    "ADR-97 wp/97b, D6): default -- stress-independent",
    """protected:

    // Ladruno (ADR-94 wp/94b, M1/F2): was `static`""",
    """    // Ladruno (ADR-97 wp/97b, D6): default -- stress-independent elasticity, so
    // the term is exactly zero.  Dropping a NONZERO one from the Jacobian would
    // still converge to the exact residual (inexact Newton), but it would make
    // `tangent_type Algorithmic` inexact; that is why el_is_stress_dependent
    // exists and is warned about rather than silently tolerated.
    ELASTICITY_STRESS_DERIVATIVE
    {
        (void) stress;
        (void) m;
        (void) parameters_storage;
        out.setZero();
    }

protected:

    // Ladruno (ADR-94 wp/94b, M1/F2): was `static`""")

edit(
    "ElasticityBase.h",
    "#include <type_traits>  // Ladruno (ADR-97",
    """#include "EigenAPI.h"
#include <Channel.h>
""",
    """#include "EigenAPI.h"
#include <Channel.h>
#include <type_traits>  // Ladruno (ADR-97 wp/97b): el_is_stress_dependent
""")

# ===========================================================================
# 6. AllASDHardeningFunctions.h -- dh/dq, dh/dm for the five P1 policies
# ===========================================================================
edit(
    "AllASDHardeningFunctions.h",
    "ADR-97 wp/97b): h = H*dev(m)",
    """        VoigtVector h = H * m.deviator();  // best not to use 'auto' here
        return h;
    }
    using parameters_t = tuple<TensorLinearHardeningParameter>;""",
    """        VoigtVector h = H * m.deviator();  // best not to use 'auto' here
        return h;
    }

    // Ladruno (ADR-97 wp/97b): h = H*dev(m) is independent of q and LINEAR in m,
    // so dh/dq = 0 and dh/dm = H * I_dev (the deviatoric projector on stored
    // Voigt slots: d(dev v)_i/d v_j = delta_ij - 1/3 for i,j < 3, delta_ij else).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        double H = GET_PARAMETER_VALUE(TensorLinearHardeningParameter);
        out.setZero();
        for (int i = 0; i < 6; ++i) out(i, i) = H;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= H / 3.0;
    }

    using parameters_t = tuple<TensorLinearHardeningParameter>;""")

edit(
    "AllASDHardeningFunctions.h",
    "ADR-97 wp/97b): h = H*sqrt(2/3",
    """        double h = H * sqrt((2 * tensor_dot_engineering_strain_like(m, m)) / 3);
        return h;
    }
    using parameters_t = tuple<ScalarLinearHardeningParameter>;""",
    """        double h = H * sqrt((2 * tensor_dot_engineering_strain_like(m, m)) / 3);
        return h;
    }

    // Ladruno (ADR-97 wp/97b): h = H*sqrt(2/3 <m,m>_e) with <.,.>_e the
    // ENGINEERING-strain contraction (shear weight 1/2).  Independent of q, so
    //     dh/dq   = 0
    //     dh/dm_j = H * (2/3) * (W_e m)_j / sqrt(2/3 <m,m>_e)
    // (one row).  At <m,m>_e = 0 the rate itself is zero and non-differentiable;
    // the derivative is taken as zero there, which is the correct one-sided limit
    // for a flow direction that has not yet been established.
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) sigma;
        double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);
        out.setZero();
        const double mm = tensor_dot_engineering_strain_like(m, m);
        const double eq = sqrt((2.0 / 3.0) * mm);
        if (eq > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
        {
            const double c = H * (2.0 / 3.0) / eq;
            out(0, 0) = c * m(0);
            out(0, 1) = c * m(1);
            out(0, 2) = c * m(2);
            out(0, 3) = c * 0.5 * m(3);
            out(0, 4) = c * 0.5 * m(4);
            out(0, 5) = c * 0.5 * m(5);
        }
    }

    using parameters_t = tuple<ScalarLinearHardeningParameter>;""")

edit(
    "AllASDHardeningFunctions.h",
    "ADR-97 wp/97b): h == 0, so both derivatives are zero (scalar)",
    """        double zero=0;
        return zero;
    }
    using parameters_t = tuple<>;""",
    """        double zero=0;
        return zero;
    }

    // Ladruno (ADR-97 wp/97b): h == 0, so both derivatives are zero (scalar).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    using parameters_t = tuple<>;""")

edit(
    "AllASDHardeningFunctions.h",
    "ADR-97 wp/97b): h == 0, so both derivatives are zero (tensor)",
    """        VoigtVector zero = VoigtVector::Zero();
        return zero;
    }
    using parameters_t = tuple<>;""",
    """        VoigtVector zero = VoigtVector::Zero();
        return zero;
    }

    // Ladruno (ADR-97 wp/97b): h == 0, so both derivatives are zero (tensor).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    using parameters_t = tuple<>;""")

edit(
    "AllASDHardeningFunctions.h",
    "ADR-97 wp/97b): Armstrong-Frederick",
    """        return derivative;

    }
    using parameters_t = tuple<AF_ha, AF_cr>;
};""",
    """        return derivative;

    }

    // Ladruno (ADR-97 wp/97b): Armstrong-Frederick, differentiated EXACTLY as
    // the branch `f` above actually took -- including the saturation branch.
    //   unsaturated:  h = ha*dev(m) - cr*|dev m|_eq * dev(alpha),
    //                 |v|_eq = sqrt(2/3 <v,v>_e)
    //     dh/dalpha = -cr * |dev m|_eq * I_dev
    //     dh/dm     =  ha * I_dev
    //                 - cr * dev(alpha) (x) [ (2/3) * (W_e dev(m)) / |dev m|_eq ]
    //   saturated (|dev alpha|_eq >= ha/cr):  h == 0, so both blocks are zero.
    //
    // ADR-97 plan-note: the plan suggested Closest_Point DROP the saturation
    // branch (it is non-differentiable, and the implicit AF update is already a
    // contraction toward ||alpha|| = ha/cr so it cannot overshoot).  It is kept
    // here instead, because `f` is SHARED with Backward_Euler (D1: byte
    // identical) and a CP-only `f` would be a second hardening law with the same
    // name.  Differentiating the branch actually taken keeps CP's Jacobian
    // exactly consistent with CP's own residual, which is what both the Newton
    // and the consistent tangent require; the P0 oracle (cppm_vm.py) mirrors the
    // same branch, so the pinned reference values are unaffected.
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) depsilon; (void) sigma;
        out.setZero();
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);
        if (!(cr > 0.0))
            return;                       // no recovery term -> h is q-independent
        VoigtVector alpha_dev = current_value.deviator();
        const double alpha_norm =
            sqrt((2. / 3.) * tensor_dot_stress_like(alpha_dev, alpha_dev));
        if (alpha_norm >= ha / cr)
            return;                       // saturated branch: h == 0
        VoigtVector mdev = m.deviator();
        const double mdev_eq =
            sqrt((2. / 3.) * tensor_dot_engineering_strain_like(mdev, mdev));
        const double c = -cr * mdev_eq;
        for (int i = 0; i < 6; ++i) out(i, i) = c;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= c / 3.0;
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) depsilon; (void) sigma;
        out.setZero();
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);
        VoigtVector alpha_dev = current_value.deviator();
        if (cr > 0.0)
        {
            const double alpha_norm =
                sqrt((2. / 3.) * tensor_dot_stress_like(alpha_dev, alpha_dev));
            if (alpha_norm >= ha / cr)
                return;                   // saturated branch: h == 0
        }
        // ha * I_dev
        for (int i = 0; i < 6; ++i) out(i, i) = ha;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= ha / 3.0;
        if (!(cr > 0.0))
            return;
        VoigtVector mdev = m.deviator();
        const double mdev_eq =
            sqrt((2. / 3.) * tensor_dot_engineering_strain_like(mdev, mdev));
        if (!(mdev_eq > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return;
        // row vector d|dev m|_eq / dm  =  (2/3) * W_e dev(m) / |dev m|_eq
        double drow[6];
        const double s = (2.0 / 3.0) / mdev_eq;
        drow[0] = s * mdev(0);
        drow[1] = s * mdev(1);
        drow[2] = s * mdev(2);
        drow[3] = s * 0.5 * mdev(3);
        drow[4] = s * 0.5 * mdev(4);
        drow[5] = s * 0.5 * mdev(5);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                out(i, j) -= cr * alpha_dev(i) * drow[j];
    }

    using parameters_t = tuple<AF_ha, AF_cr>;
};

// Ladruno (ADR-97 wp/97b): the five hardening policies P1 solves implicitly at
// n+1.  Every other policy in the tree (the StiffSoil pair) keeps the inert zero
// default and is REFUSED at parse time under `integration_method Closest_Point`.
template <> struct hardening_policy_has_cp_derivatives<LinearHardeningForTensorPolicy> : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<LinearHardeningForScalarPolicy> : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<NullHardeningScalarPolicy>      : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<NullHardeningTensorPolicy>      : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<ArmstrongFrederickPolicy>       : std::true_type {};""")

# ===========================================================================
# 7. VonMises_YF.h -- df/dq + trait
# ===========================================================================
edit(
    "YieldFunctions/VonMises_YF.h",
    "ADR-97 wp/97b): df/dq, UNCONTRACTED",
    """    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(r:r) - sqrt(2/3)*k, so the term that""",
    """    // Ladruno (ADR-97 wp/97b): df/dq, UNCONTRACTED, in the same VOIGT convention
    // as df_dsigma_ij and yf_hardening's df_dalpha (shear slots doubled).
    //   f = ||r||_s - sqrt(2/3)*k ,  r = dev(sigma) - alpha, ||r||_s^2 = r' W_s r
    //   df/dk      = -sqrt(2/3)                       (the truncated literal)
    //   df/dalpha  = -(W_s r)/||r||_s  ==  -(r/||r||_s) with the shear slots x2
    // Any other internal variable (e.g. the plastic-flow direction's own back
    // stress when it carries a different hardening law, which makes it a
    // DISTINCT entry of the IV storage) does not appear in f: zero.
    YIELD_FUNCTION_IV_DERIVATIVE
    {
        (void) parameters_storage;
        const int nq = iv.size();
        for (int i = 0; i < nq; ++i) out[i] = 0.0;

        if constexpr (std::is_same<IVType, KHardeningType>::value)
        {
            out[0] = -SQRT_2_over_3;
        }
        else if constexpr (std::is_same<IVType, AlphaHardeningType>::value)
        {
            auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
            VoigtVector r = sigma.deviator() - alpha;
            const double den = sqrt(tensor_dot_stress_like(r, r));
            if (den > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            {
                VoigtVector d = -r / den;
                d(3) *= 2.0;
                d(4) *= 2.0;
                d(5) *= 2.0;
                for (int i = 0; i < 6; ++i) out[i] = d(i);
            }
        }
    }

    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(r:r) - sqrt(2/3)*k, so the term that""")

edit(
    "YieldFunctions/VonMises_YF.h",
    "yf_has_cp_derivatives<VonMises_YF",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.


#endif""",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): von Mises supplies the analytic closest-point df/dq.
template<class AlphaHardeningType, class KHardeningType>
struct yf_has_cp_derivatives<VonMises_YF<AlphaHardeningType, KHardeningType>> : std::true_type {};


#endif""")

# ===========================================================================
# 8. DruckerPrager_YF.h -- df/dq + trait
# ===========================================================================
edit(
    "YieldFunctions/DruckerPrager_YF.h",
    "ADR-97 wp/97b): df/dq, UNCONTRACTED",
    """    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(J2) + eta*p - xi_c, so xi_c (the""",
    """    // Ladruno (ADR-97 wp/97b): df/dq, UNCONTRACTED, in the VOIGT convention
    // (normal slots halved, exactly as df_dsigma_ij and yf_hardening's
    // df_dalpha).  f = sqrt(J2(r)) + eta*p - xi_c , r = dev(sigma) - alpha:
    //   df/dalpha = -(r/sqrt(J2)) with the three NORMAL slots halved
    //   df/d(cohesion IV) = 0
    //
    // That zero is deliberate and is ADR-97 P0 header finding 2: this yield
    // function's cohesion internal variable is COMMENTED OUT of its own `f`
    // (line 26) while `yf_hardening` still contributes df/dk = -1 times its
    // rate.  `Closest_Point` must use the TRUE derivative of the f it solves, so
    // a Drucker-Prager with cohesion hardening is PERFECTLY PLASTIC under
    // Closest_Point (the cohesion IV still evolves, it just does not move the
    // surface) and hardening under Backward_Euler.  Restoring `k` in `f` -- or
    // deleting the hardening term -- changes Backward_Euler, which ADR-97 D1
    // keeps byte-identical; it is a separate PR.
    YIELD_FUNCTION_IV_DERIVATIVE
    {
        (void) parameters_storage;
        const int nq = iv.size();
        for (int i = 0; i < nq; ++i) out[i] = 0.0;

        if constexpr (std::is_same<IVType, AlphaHardeningType>::value)
        {
            auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
            VoigtVector r = sigma.deviator() - alpha;
            const double den = sqrt(0.5 * tensor_dot_stress_like(r, r));
            if (den > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            {
                VoigtVector d = -r / den;
                d(0) *= 0.5;
                d(1) *= 0.5;
                d(2) *= 0.5;
                for (int i = 0; i < 6; ++i) out[i] = d(i);
            }
        }
    }

    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(J2) + eta*p - xi_c, so xi_c (the""")

edit(
    "YieldFunctions/DruckerPrager_YF.h",
    "yf_has_cp_derivatives<DruckerPrager_YF",
    """// Declares this YF as featuring an apex
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};""",
    """// Declares this YF as featuring an apex
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};

// Ladruno (ADR-97 wp/97b): Drucker-Prager supplies the analytic closest-point
// df/dq.  NOTE that `check_apex_region` above stays EUCLIDEAN and is used only
// by Backward_Euler: `Closest_Point` classifies the apex region in the ELASTIC
// metric, inside the integrator where K, G and etabar are in scope (ADR-97
// "Drucker-Prager apex").  Two integrators, two answers for the same YF -- a
// documentation obligation (LEDGER_quirks), not a bug.
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_cp_derivatives<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};""")

# ===========================================================================
# 9. VonMises_PF.h -- dm/dsigma, dm/dq + trait
# ===========================================================================
edit(
    "PlasticFlowDirections/VonMises_PF.h",
    "ADR-97 wp/97b): dm/dsigma and dm/dq",
    """    using internal_variables_t = std::tuple<AlphaHardeningType>;
    using parameters_t = std::tuple<>;""",
    """    // Ladruno (ADR-97 wp/97b): dm/dsigma and dm/dq for the closest-point map.
    //
    // With r = dev(sigma) - alpha, N = ||r||_s = sqrt(r' W_s r) and the shipped
    // (Voigt) flow direction m = W_s r / N:
    //
    //     dm/dsigma = ( diag(W_s) * I_dev - m (x) m ) / N
    //     dm/dalpha = -( diag(W_s)         - m (x) m ) / N
    //
    // (I_dev is the deviatoric projector on STORED slots; the second term is the
    // same in both because I_dev * m == m for a deviatoric r.)  Verified against
    // a central difference of this very expression at 2.6e-11 / 5.7e-11 before
    // the C++ was written.  Associated flow, so dm/dsigma is the yield Hessian
    // and C_alg is symmetric for von Mises.
    PLASTIC_FLOW_STRESS_DERIVATIVE
    {
        (void) depsilon;
        (void) parameters_storage;
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        VoigtVector r = sigma.deviator() - alpha;
        const double N = sqrt(tensor_dot_stress_like(r, r));
        dm_dsigma_buffer.setZero();
        if (!(N > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return dm_dsigma_buffer;      // degenerate: m is not defined here
        VoigtVector mv = r / N;
        mv(3) *= 2.0; mv(4) *= 2.0; mv(5) *= 2.0;      // m = W_s r / N
        // diag(W_s) * I_dev
        for (int i = 0; i < 6; ++i)
            dm_dsigma_buffer(i, i) = (i < 3) ? 1.0 : 2.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) dm_dsigma_buffer(i, j) -= 1.0 / 3.0;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                dm_dsigma_buffer(i, j) = (dm_dsigma_buffer(i, j) - mv(i) * mv(j)) / N;
        return dm_dsigma_buffer;
    }

    PLASTIC_FLOW_IV_DERIVATIVE
    {
        (void) depsilon;
        (void) parameters_storage;
        out.setZero();
        if constexpr (std::is_same<IVType, AlphaHardeningType>::value)
        {
            (void) iv;
            auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
            VoigtVector r = sigma.deviator() - alpha;
            const double N = sqrt(tensor_dot_stress_like(r, r));
            if (!(N > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
                return;
            VoigtVector mv = r / N;
            mv(3) *= 2.0; mv(4) *= 2.0; mv(5) *= 2.0;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j)
                    out(i, j) = ((i == j ? ((i < 3) ? 1.0 : 2.0) : 0.0)
                                 - mv(i) * mv(j)) / (-N);
        }
    }

    using internal_variables_t = std::tuple<AlphaHardeningType>;
    using parameters_t = std::tuple<>;""")

edit(
    "PlasticFlowDirections/VonMises_PF.h",
    "pf_has_cp_derivatives<VonMises_PF",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.""",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): analytic closest-point derivatives available.
template<class AlphaHardeningType>
struct pf_has_cp_derivatives<VonMises_PF<AlphaHardeningType>> : std::true_type {};""")

# ===========================================================================
# 10. DruckerPrager_PF.h -- dm/dsigma, dm/dq + trait
# ===========================================================================
edit(
    "PlasticFlowDirections/DruckerPrager_PF.h",
    "ADR-97 wp/97b): dm/dsigma and dm/dq",
    """    using internal_variables_t = std::tuple<AlphaHardeningType, EtaHardeningType>;
    using parameters_t = std::tuple<DP_etabar>;""",
    """    // Ladruno (ADR-97 wp/97b): dm/dsigma and dm/dq for the closest-point map.
    //
    // With r = dev(sigma) - alpha, q = sqrt(J2) = sqrt(0.5 r' W_s r) and the
    // shipped flow direction m = (W_s/2) r / q + (etabar/3) delta, whose
    // deviatoric part is md = (W_s/2) r / q:
    //
    //     dm/dsigma = diag(W_s) * I_dev / (2q)  -  md (x) md / q
    //     dm/dalpha = -[ diag(W_s) / (2q)       -  md (x) md / q ]
    //
    // (the etabar*p term has zero Hessian).  Verified against a central
    // difference at 1.3e-11 / 3.7e-11 before the C++ was written.  Non-associated
    // whenever etabar != eta, so C_alg is UNSYMMETRIC -- tests run on UmfPack.
    // At q -> 0 the expression diverges as 1/q: that is the apex, which
    // `Closest_Point` classifies in the elastic metric and returns with its own
    // reduced system, so this returns zero rather than an infinity.
    PLASTIC_FLOW_STRESS_DERIVATIVE
    {
        (void) depsilon;
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        double etabar = GET_PARAMETER_VALUE(DP_etabar);
        (void) etabar;
        VoigtVector r = sigma.deviator() - alpha;
        const double q = sqrt(0.5 * tensor_dot_stress_like(r, r));
        dm_dsigma_buffer.setZero();
        if (!(q > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return dm_dsigma_buffer;      // apex: handled by the integrator
        VoigtVector md = r / (2.0 * q);
        md(3) *= 2.0; md(4) *= 2.0; md(5) *= 2.0;      // md = (W_s/2) r / q
        for (int i = 0; i < 6; ++i)
            dm_dsigma_buffer(i, i) = ((i < 3) ? 1.0 : 2.0) / (2.0 * q);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                dm_dsigma_buffer(i, j) -= 1.0 / (3.0 * 2.0 * q);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                dm_dsigma_buffer(i, j) -= md(i) * md(j) / q;
        return dm_dsigma_buffer;
    }

    PLASTIC_FLOW_IV_DERIVATIVE
    {
        (void) depsilon;
        out.setZero();
        if constexpr (std::is_same<IVType, AlphaHardeningType>::value)
        {
            (void) iv;
            auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
            VoigtVector r = sigma.deviator() - alpha;
            const double q = sqrt(0.5 * tensor_dot_stress_like(r, r));
            if (!(q > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
                return;
            VoigtVector md = r / (2.0 * q);
            md(3) *= 2.0; md(4) *= 2.0; md(5) *= 2.0;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j)
                    out(i, j) = -(((i == j) ? (((i < 3) ? 1.0 : 2.0) / (2.0 * q)) : 0.0)
                                  - md(i) * md(j) / q);
        }
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, EtaHardeningType>;
    using parameters_t = std::tuple<DP_etabar>;""")

edit(
    "PlasticFlowDirections/DruckerPrager_PF.h",
    "pf_has_cp_derivatives<DruckerPrager_PF",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

#endif""",
    """// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): analytic closest-point derivatives available.
template<class AlphaHardeningType, class EtaHardeningType>
struct pf_has_cp_derivatives<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType>> : std::true_type {};

#endif""")

# ===========================================================================
# 11. Stress-dependent elasticity models -> el_is_stress_dependent
# ===========================================================================
edit(
    "ElasticityModels/StiffSoil_EL.h",
    "el_is_stress_dependent<StiffSoil_EL>",
    """    using parameters_t = std::tuple<SS_Eur_ref, PoissonsRatio, SS_pref, SS_m, MC_phi, MC_c>;""",
    """    using parameters_t = std::tuple<SS_Eur_ref, PoissonsRatio, SS_pref, SS_m, MC_phi, MC_c>;

    // Ladruno (ADR-97 wp/97b, D6): declared stress dependent -- see the
    // specialization after this class.""")

edit(
    "ElasticityModels/DuncanChang_EL.h",
    "el_is_stress_dependent<DuncanChang_EL>",
    """    using parameters_t = std::tuple<ReferenceYoungsModulus,PoissonsRatio,ReferencePressure,DuncanChang_MaxSigma3,DuncanChang_n>;""",
    """    using parameters_t = std::tuple<ReferenceYoungsModulus,PoissonsRatio,ReferencePressure,DuncanChang_MaxSigma3,DuncanChang_n>;

    // Ladruno (ADR-97 wp/97b, D6): declared stress dependent -- see the
    // specialization after this class.""")

edit(
    "ElasticityModels/StiffSoil_EL.h",
    "struct el_is_stress_dependent<StiffSoil_EL>",
    "#endif",
    """// Ladruno (ADR-97 wp/97b, D6): E depends on sigma, so the exact algorithmic
// tangent carries an `E,sigma : (sigma - sigma_tr)` term that ELASTICITY_STRESS_
// DERIVATIVE does not yet supply.  Closest_Point evaluates E(sigma_{n+1}) inside
// the residual regardless (the converged state IS hyperelastically consistent);
// this trait only drives the one-time `tangent_type Algorithmic` warning.
// ADR-97 P5 supplies the block.
template <>
struct el_is_stress_dependent<StiffSoil_EL> : std::true_type {};

#endif""", last=True)

edit(
    "ElasticityModels/DuncanChang_EL.h",
    "struct el_is_stress_dependent<DuncanChang_EL>",
    "#endif",
    """// Ladruno (ADR-97 wp/97b, D6): see the StiffSoil_EL specialization.
template <>
struct el_is_stress_dependent<DuncanChang_EL> : std::true_type {};

#endif""", last=True)

# ===========================================================================
# 12. ASDPlasticMaterial3DGlobals.h -- the new integration-method enum value
# ===========================================================================
edit(
    "ASDPlasticMaterial3DGlobals.h",
    "Closest_Point,",
    """    Forward_Euler_Subincrement,
    Backward_Euler_LineSearch,
};""",
    """    Forward_Euler_Subincrement,
    Backward_Euler_LineSearch,
    Closest_Point,   // Ladruno (ADR-97 wp/97b): fully implicit closest-point return map
};""")


# ===========================================================================
# apply
# ===========================================================================
def apply_edits(edits, base):
    """Apply (rel, marker, anchor, replacement, last) edits under `base`.

    Line endings are preserved: the sources in this tree are CRLF, so the
    LF-written anchors/replacements are converted to match the file before any
    matching happens.
    """
    n_applied = 0
    n_skipped = 0
    for rel, marker, anchor, replacement, last in edits:
        path = os.path.join(base, rel.replace("/", os.sep))
        if not os.path.exists(path):
            raise SystemExit("MISSING FILE: %s" % path)
        with open(path, "r", encoding="utf-8", errors="surrogateescape",
                  newline="") as f:
            src = f.read()
        crlf = "\r\n" in src
        mk = marker.replace("\n", "\r\n") if crlf else marker
        an = anchor.replace("\n", "\r\n") if crlf else anchor
        rp = replacement.replace("\n", "\r\n") if crlf else replacement
        if mk in src:
            n_skipped += 1
            continue
        cnt = src.count(an)
        if cnt == 0 or (cnt != 1 and not last):
            raise SystemExit(
                "ANCHOR NOT UNIQUE (%d matches) in %s for marker %r:\n---\n%s\n---"
                % (cnt, rel, marker, anchor[:400]))
        if cnt == 1:
            src = src.replace(an, rp)
        else:
            k = src.rfind(an)
            src = src[:k] + rp + src[k + len(an):]
        with open(path, "w", encoding="utf-8", errors="surrogateescape",
                  newline="") as f:
            f.write(src)
        n_applied += 1
        print("  applied  %-46s %s" % (rel, marker[:52]))
    return n_applied, n_skipped


def main():
    a, s = apply_edits(EDITS, ASD)
    print("apply_p1_cpp: %d applied, %d already present" % (a, s))


if __name__ == "__main__":
    main()
