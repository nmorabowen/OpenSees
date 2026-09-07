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

#ifndef VonMises_PF_H
#define VonMises_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"

#include <cmath>
#include <typeinfo>


using namespace std;

template<class AlphaHardeningType>
class VonMises_PF : public PlasticFlowBase<VonMises_PF<AlphaHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "VonMises_PF";

    VonMises_PF():
        PlasticFlowBase<VonMises_PF<AlphaHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);

        auto s = sigma.deviator();
        vv_out = s - alpha;

        double den = sqrt(tensor_dot_stress_like(vv_out, vv_out));
        if (abs(den) > sqrt(tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            vv_out = vv_out / den;

        // Ladruno (ADR-94 wp/94c, B5): VOIGT/engineering convention, matching
        // VonMises_YF::df_dsigma_ij (associated flow).  d(eps^p)_12 = dLambda *
        // dg/d sigma_12, and the Voigt slot stores gamma^p_12 = 2*(eps^p)_12, so
        // the shear slots of m carry the factor 2.  `TrialPlastic_Strain +=
        // dLambda*m` and `Eelastic*m` (whose shear diagonal is mu, not 2*mu) both
        // require exactly this.  See ADR-94 B5.
        vv_out(3) *= 2.0;   // gamma12
        vv_out(4) *= 2.0;   // gamma23
        vv_out(5) *= 2.0;   // gamma13

        return vv_out;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType>;
    using parameters_t = std::tuple<>;

private:

    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.


#endif