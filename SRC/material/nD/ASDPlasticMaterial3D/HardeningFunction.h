// Ladruno (HB/StiffSoil integration, ledger row 337): adds include guard
#ifndef HardeningFunctionBase_H
#define HardeningFunctionBase_H

#include <type_traits>   // Ladruno (ADR-97 wp/97b): std::false_type for the CP traits


#define HARDENING_FUNCTION_DEFINITION template <class EVT, class ParameterStorageType> \
    static auto f( \
        const EVT& current_value, \
        const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma, \
        const ParameterStorageType& parameters_storage) 


// Ladruno (ADR-97 wp/97b): the closest-point Jacobian needs dh/dq and dh/dm
// (see ADR-97 "Residuals and Jacobians": J_qs = -dl*(dh/dm)(dm/ds),
// J_qq = I - dl*(dh/dq + (dh/dm)(dm/dq))).  `dh/ds` is deliberately NOT
// declared: no hardening policy in this tree reads `sigma` (they read
// `current_value`, `m` and `depsilon`), and one that later does must declare it.
//
// The blocks are written into the leading rows/columns of a 6x6 buffer:
//   dh/dq is n_iv x n_iv, dh/dm is n_iv x 6, with n_iv = 1 (scalar IV) or 6.
// A policy without the derivatives keeps a ZERO default and is REFUSED at parse
// time under Closest_Point, so the zero is never silently believed.
#define HARDENING_FUNCTION_IV_DERIVATIVE template <class EVT, class ParameterStorageType>     static void dh_dq(         const EVT& current_value,         const VoigtVector& depsilon,         const VoigtVector& m,         const VoigtVector& sigma,         const ParameterStorageType& parameters_storage,         VoigtMatrix& out)

#define HARDENING_FUNCTION_M_DERIVATIVE template <class EVT, class ParameterStorageType>     static void dh_dm(         const EVT& current_value,         const VoigtVector& depsilon,         const VoigtVector& m,         const VoigtVector& sigma,         const ParameterStorageType& parameters_storage,         VoigtMatrix& out)

// Opt-in trait, keyed on the POLICY (specialized in AllASDHardeningFunctions.h).
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

// Function wrapper base class
template <typename EvolvingVariableType, class HardeningPolicy>
struct HardeningFunction {
	HARDENING_FUNCTION_DEFINITION 
    {
        return HardeningPolicy::f(current_value, depsilon, m, sigma, parameters_storage);
    }
    // Ladruno (ADR-97 wp/97b): forwarders with an inert zero default.
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
    : hardening_policy_has_cp_derivatives<Policy> {};

// Ladruno (ADR-97 wp/97c): the same lift, for the inert-hardening trait.
template <typename T>
struct hardening_is_inert : std::false_type {};

template <typename EVT, typename Policy>
struct hardening_is_inert<HardeningFunction<EVT, Policy>>
    : hardening_policy_is_inert<Policy> {};


template <typename EvolvingVariableType, class HardeningPolicy>
std::ostream& operator<<(std::ostream& os, const HardeningFunction<EvolvingVariableType, HardeningPolicy>& obj) {
    os << "HardeningFunction<" << typeid(EvolvingVariableType).name() << ", " << typeid(HardeningPolicy).name() << ">";
    return os;
}



#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

#endif  //HardeningFunctionBase_H


