/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial3D
//
// Fully general templated material class for plasticity modeling

#ifndef YieldFunctionBase_H
#define YieldFunctionBase_H

#include "EigenAPI.h"
#include <Channel.h>
#include <type_traits>


// Trait for enabling special algorithms in YF that have apexes

// Base template assuming no apex member exists
template <typename T>
struct yf_has_apex : std::false_type {};

// Ladruno (ADR-84 P0): opt-in trait for YFs that implement their own exact
// return for stress states where the generic scalar-Newton return map is
// invalid (multi-surface corners, cutoff planes, apex cones). A YF that
// specializes this to true_type must define the SPECIAL_RETURN member below.
template <typename T>
struct yf_has_special_return : std::false_type {};

// Ladruno (ADR-97 wp/97b): opt-in trait for yield functions that supply the
// UNCONTRACTED closest-point Jacobian block df/dq (YIELD_FUNCTION_IV_DERIVATIVE
// below).  `integration_method Closest_Point` is REFUSED at parse time for a YF
// that does not specialize this to true_type -- never silently approximated.
// See ADR-97 D3 (family-by-family) and the base default in YieldFunctionBase.
template <typename T>
struct yf_has_cp_derivatives : std::false_type {};



// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct yf_has_parameters_t : std::false_type {};

template <typename T>
struct yf_has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

// Helper template to check if a class has a internal_variables_t type alias
template <typename T, typename = void>
struct yf_has_internal_variables_t : std::false_type {};

template <typename T>
struct yf_has_internal_variables_t<T, typename std::enable_if<!std::is_same<typename T::internal_variables_t, void>::value>::type> : std::true_type {};


#define YIELD_FUNCTION template <typename IVStorageType, typename ParameterStorageType> \
    double operator()( const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage) const 

#define YIELD_FUNCTION_STRESS_DERIVATIVE template <typename IVStorageType, typename ParameterStorageType> \
    const VoigtVector& df_dsigma_ij(const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage) const

#define YIELD_FUNCTION_HARDENING template <typename IVStorageType, typename ParameterStorageType> \
    double hardening(const VoigtVector& depsilon, \
        const VoigtVector& m, \
        const VoigtVector& sigma,\
        const IVStorageType& internal_variables_storage,\
        const ParameterStorageType& parameters_storage) const

#define CHECK_APEX_REGION template <typename IVStorageType, typename ParameterStorageType> \
    bool check_apex_region( const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage) const 

#define APEX_STRESS template <typename IVStorageType, typename ParameterStorageType> \
    const VoigtVector& apex_stress(const IVStorageType& internal_variables_storage, \
                        const ParameterStorageType& parameters_storage) const

// Ladruno (ADR-94 wp/94c, M5): the YF's own strength scale -- the stress magnitude
// that normalises f -- so a RELATIVE yield tolerance can be formed:
//     tol = max(f_absolute_tol, f_relative_tol * strength_scale)
// ADR-94 M5 measured the same Mohr-Coulomb problem passing 20/20 in kPa and being
// refused on step 1 in Pa at the default absolute 1e-6: the unit system decided
// pass/fail.  A YF that does not define one inherits the base's 0.0, which makes
// f_relative_tol inert for it (the absolute tolerance still applies).
#define YF_STRENGTH_SCALE template <typename IVStorageType, typename ParameterStorageType> \
    double strength_scale(const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage) const

// Ladruno (ADR-84 P0): signature for the opt-in special return (see
// yf_has_special_return above). Called by the constitutive integrator after
// the elastic check and before the scalar-Newton plastic correction, with the
// trial (elastic-predictor) stress and the state-frozen elastic tangent.
// Returns false to fall through to the generic return map, or true after
// writing the returned stress, the plastic-strain increment, and the tangent
// to use. The argument names match the other macros so GET_PARAMETER_VALUE /
// YF(...) work inside the body.
// Ladruno (ADR-84 P3): `stiffness_return` is the RAW active-set (Koiter)
// consistent tangent -- the YF does NOT apply a tangent-operator policy. The
// integrator applies the material's configured `tangent_type` to it, exactly as
// it does on the generic path. (P0 secant-blended inside the YF, which silently
// overrode `tangent_type` and fabricated stiffness in the degenerate corner
// direction -- see the ADR's §9 measurement.)
// `return_quality` reports how the state was resolved: SR_QUALITY_EXACT for a
// closed-form face/edge/corner/apex-cone return, SR_QUALITY_FALLBACK for the
// conservative terminal vertex projection (a large, unphysical stress drop that
// the caller may refuse under strict_convergence).
#define SR_QUALITY_EXACT    1
#define SR_QUALITY_FALLBACK 2

#define SPECIAL_RETURN template <typename IVStorageType, typename ParameterStorageType> \
    bool special_return(const VoigtVector& sigma_trial, \
        const VoigtMatrix& Eelastic, \
        const double tol_yf, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage, \
        VoigtVector& sigma_return, \
        VoigtVector& plastic_strain_incr, \
        VoigtMatrix& stiffness_return, \
        int& return_quality) const

// Ladruno (ADR-97 wp/97b): df/dq -- the derivative of f with respect to ONE
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
#define YIELD_FUNCTION_IV_DERIVATIVE template <typename IVType, typename IVStorageType, typename ParameterStorageType> \
    void df_dq(const IVType& iv, const VoigtVector& sigma, \
        const IVStorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage, \
        double* out) const

#define GET_INTERNAL_VARIABLE_HARDENING(type) \
    internal_variables_storage.template get<type> ().hardening_function(depsilon, m, sigma, parameters_storage)

#define GET_TRIAL_INTERNAL_VARIABLE(type) \
    ((internal_variables_storage).template get<type>().trial_value)

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

#define YF(SIGMA) (this)->operator()(SIGMA, internal_variables_storage, parameters_storage)

template <class T>
class YieldFunctionBase
{
public:
    YieldFunctionBase() { 
        static_assert(yf_has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
        static_assert(yf_has_internal_variables_t<T>::value, "Derived class must have a 'internal_variables_t' type alias.");
    }

    YIELD_FUNCTION
    {
        return static_cast<T*>(this)->operator()(sigma, internal_variables_storage, parameters_storage);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        return static_cast<T*>(this)->df_dsigma_ij(sigma, internal_variables_storage, parameters_storage);
    }

    YIELD_FUNCTION_HARDENING
    {
        return static_cast<T*>(this)->df_dxi_star_h_star(depsilon, m , sigma, internal_variables_storage, parameters_storage);
    }

    // Ladruno (ADR-94 wp/94c, M5): default -- this YF declares no strength scale,
    // so `f_relative_tol` contributes nothing and only `f_absolute_tol` binds.
    YF_STRENGTH_SCALE
    {
        (void) internal_variables_storage;
        (void) parameters_storage;
        return 0.0;
    }

    // Ladruno (ADR-97 wp/97b): default -- this YF declares no closest-point
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

#endif
