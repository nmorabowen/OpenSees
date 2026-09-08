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

    // Ladruno (ADR-97 wp/97b): dm/dsigma and dm/dq for the closest-point map.
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
        this->dm_dsigma_buffer.setZero();
        if (!(N > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return this->dm_dsigma_buffer;      // degenerate: m is not defined here
        VoigtVector mv = r / N;
        mv(3) *= 2.0; mv(4) *= 2.0; mv(5) *= 2.0;      // m = W_s r / N
        // diag(W_s) * I_dev
        for (int i = 0; i < 6; ++i)
            this->dm_dsigma_buffer(i, i) = (i < 3) ? 1.0 : 2.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) this->dm_dsigma_buffer(i, j) -= 1.0 / 3.0;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                this->dm_dsigma_buffer(i, j) = (this->dm_dsigma_buffer(i, j) - mv(i) * mv(j)) / N;
        return this->dm_dsigma_buffer;
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
    using parameters_t = std::tuple<>;

private:

    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): analytic closest-point derivatives available.
template<class AlphaHardeningType>
struct pf_has_cp_derivatives<VonMises_PF<AlphaHardeningType>> : std::true_type {};


#endif