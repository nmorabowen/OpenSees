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

#ifndef ElasticityBase_H
#define ElasticityBase_H

#include "EigenAPI.h"
#include <Channel.h>
#include <type_traits>  // Ladruno (ADR-97 wp/97b): el_is_stress_dependent

// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct has_parameters_t : std::false_type {};

template <typename T>
struct has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

#define ELASTICITY_MATRIX template<class ParameterStorageType> \
    const VoigtMatrix& operator()(const VoigtVector& stress, \
        const ParameterStorageType& parameters_storage) const 

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value

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
#define ELASTICITY_STRESS_DERIVATIVE template<class ParameterStorageType> \
    void dE_dsigma_contract(const VoigtVector& stress, \
        const VoigtVector& m, \
        const ParameterStorageType& parameters_storage, \
        VoigtMatrix& out) const


template <class T>
class ElasticityBase
{
public:

    ElasticityBase() {
        static_assert(has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
        EE_MATRIX.setZero();  // Ladruno (ADR-94 wp/94b): per-instance now, so zero it here
    }
    
    ELASTICITY_MATRIX
    {
        return static_cast<T*>(this)->operator()(stress, parameters_storage);
    }

    // Ladruno (ADR-97 wp/97b, D6): default -- stress-independent elasticity, so
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

    // Ladruno (ADR-94 wp/94b, M1/F2): was `static` -- one elastic-tangent buffer shared
    // by every material instance using this elasticity model (the comment in
    // LinearIsotropic3D_EL::operator() admits as much: "It may have values from another
    // instance with different parameters"). `mutable` because operator() is const.
    mutable VoigtMatrix EE_MATRIX;
};



#endif
