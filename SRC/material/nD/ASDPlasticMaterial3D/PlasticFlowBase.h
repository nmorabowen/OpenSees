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

#ifndef PlasticFlowBase_H
#define PlasticFlowBase_H

#include <Channel.h>

#include "EigenAPI.h"


// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct pf_has_parameters_t : std::false_type {};

template <typename T>
struct pf_has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

// Helper template to check if a class has a internal_variables_t type alias
template <typename T, typename = void>
struct pf_has_internal_variables_t : std::false_type {};

template <typename T>
struct pf_has_internal_variables_t<T, typename std::enable_if<!std::is_same<typename T::internal_variables_t, void>::value>::type> : std::true_type {};


// Ladruno (ADR-97 wp/97b): opt-in trait for plastic-flow directions that supply
// ANALYTIC closest-point Jacobian blocks (dm/dsigma, dm/dq).  A PF without it
// still compiles -- it inherits the one-time-warned central-difference dm/dsigma
// below -- but `integration_method Closest_Point` is refused at parse time for
// it, so no family this ADR has not converted can silently run on a finite
// difference.  See ADR-97 D3.
template <typename T>
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
#define CP_PRINCIPAL_MC_FLOW_PARAMS template <typename StorageType, typename ParameterStorageType> \
    bool cp_mc_flow_params(const StorageType& internal_variables_storage, \
        const ParameterStorageType& parameters_storage, \
        double& sin_phi, double& sin_psi) const

// Ladruno (ADR-97 wp/97b): dm/dsigma, the 6x6 derivative of the flow direction
// with respect to the STORED stress slots -- the `dl * dm/ds` term of the
// algorithmic elastic modulus Xi = (E^-1 + dl*dm/ds)^-1.  ADR-94 M3 measured
// exactly this term: setting dl = 0 gives Xi = E and recovers the shipped
// `Continuum` operator, whose error against a central difference of the
// material's own committed response was 57 %.
#define PLASTIC_FLOW_STRESS_DERIVATIVE template <typename StorageType, typename ParameterStorageType> \
    const VoigtMatrix& dm_dsigma( \
        const VoigtVector &depsilon, \
        const VoigtVector& sigma, \
        const StorageType& internal_variables_storage,  \
        const ParameterStorageType& parameters_storage) const

// Ladruno (ADR-97 wp/97b): dm/dq for ONE internal variable -- a 6 x iv.size()
// block written into the leading columns of `out` (a 6x6 buffer; both internal
// variable kinds in this tree are size 1 or 6).  `out` is zeroed by the caller.
#define PLASTIC_FLOW_IV_DERIVATIVE template <typename IVType, typename StorageType, typename ParameterStorageType> \
    void dm_dq(const IVType& iv, \
        const VoigtVector &depsilon, \
        const VoigtVector& sigma, \
        const StorageType& internal_variables_storage,  \
        const ParameterStorageType& parameters_storage, \
        VoigtMatrix& out) const

#define PLASTIC_FLOW_DIRECTION template <typename StorageType, typename ParameterStorageType> \
    const VoigtVector& operator()( \
        const VoigtVector &depsilon, \
        const VoigtVector& sigma, \
        const StorageType& internal_variables_storage,  \
        const ParameterStorageType& parameters_storage) const

#define GET_TRIAL_INTERNAL_VARIABLE(type) \
    ((internal_variables_storage).template get<type>().trial_value)

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value


template <class T>
class PlasticFlowBase
{
public:
    PlasticFlowBase() { 
        static_assert(pf_has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
        static_assert(pf_has_internal_variables_t<T>::value, "Derived class must have a 'internal_variables_t' type alias.");
    }

    PLASTIC_FLOW_DIRECTION
    {
        return static_cast<T*>(this)->operator()( depsilon,  sigma, internal_variables_storage, parameters_storage);
    }

    // Ladruno (ADR-97 wp/97b): default dm/dsigma -- a central difference of this
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

    // Ladruno (ADR-97 wp/97c): default -- this flow direction is not a
    // Mohr-Coulomb potential.
    CP_PRINCIPAL_MC_FLOW_PARAMS
    {
        (void) internal_variables_storage;
        (void) parameters_storage;
        sin_phi = 0.0;
        sin_psi = 0.0;
        return false;
    }

    inline const char* getName() const { return static_cast<T*>(this)->NAME; }

protected:

    // Ladruno (ADR-97 wp/97b): per-instance return buffer for the finite-
    // difference dm/dsigma default (never static -- ADR-94 wp/94b F2).
    mutable VoigtMatrix dm_dsigma_buffer = VoigtMatrix::Zero();
};

#endif
