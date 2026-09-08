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

#ifndef DruckerPrager_PF_H
#define DruckerPrager_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
using namespace ASDPlasticMaterial3DGlobals;

#include <cmath>
#include <typeinfo>


template<class AlphaHardeningType, class EtaHardeningType>
class DruckerPrager_PF : public PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_PF";


    DruckerPrager_PF( ):
        PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        auto s = sigma.deviator();
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        
        // Get dilation parameter (etabar) from parameters
        double etabar = GET_PARAMETER_VALUE(DP_etabar);  // dilation parameter (controls plastic volume change)
        
        // Compute deviatoric part
        VoigtVector dev_part = s - alpha;
        double den = sqrt(0.5*tensor_dot_stress_like(dev_part, dev_part));
        
        if (abs(den) > sqrt(0.5*tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            dev_part = dev_part / den;
        else
            dev_part.setZero();  // Ladruno (ADR-94 wp/94a): was `*= 0.0`; NaN*0 == NaN

        // Ladruno (ADR-94 wp/94c, B4/B5): identical correction to
        // DruckerPrager_YF::df_dsigma_ij -- d sqrt(J2)/d v is r/(2 sqrt(J2)) on the
        // three normal slots and r/sqrt(J2) on the three shear slots.  The flat
        // `r/den` was 2x too large on the normal slots, so the plastic flow
        // direction (and with etabar == eta the associativity check that compares
        // m against n) was wrong for every non-pure-shear stress state.
        dev_part(0) *= 0.5;
        dev_part(1) *= 0.5;
        dev_part(2) *= 0.5;
            
        // Add pressure-dependent part: etabar * dp/dsigma = etabar/3 * I
        // Ladruno (ADR-94 wp/94a): ADR-94 B4 -- `VoigtVector x; x *= 0.0;` multiplies
        // UNINITIALISED Eigen storage by zero, which does not clear heap garbage that
        // decodes as NaN. This is where the Drucker-Prager hydrostatic-tension NaN was
        // born (and then committed with analyze() == 0). Same fix as the ADR-84 P0
        // constructor precedent in ASDPlasticMaterial3D.h.
        VoigtVector pressure_part = VoigtVector::Zero();
        pressure_part(0) = etabar / 3.0;  // sigma_xx component
        pressure_part(1) = etabar / 3.0;  // sigma_yy component  
        pressure_part(2) = etabar / 3.0;  // sigma_zz component
        // shear components remain zero
        
        vv_out = dev_part + pressure_part;
        
        return vv_out;
    }



    // Ladruno (ADR-97 wp/97b): dm/dsigma and dm/dq for the closest-point map.
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
        this->dm_dsigma_buffer.setZero();
        if (!(q > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return this->dm_dsigma_buffer;      // apex: handled by the integrator
        VoigtVector md = r / (2.0 * q);
        md(3) *= 2.0; md(4) *= 2.0; md(5) *= 2.0;      // md = (W_s/2) r / q
        for (int i = 0; i < 6; ++i)
            this->dm_dsigma_buffer(i, i) = ((i < 3) ? 1.0 : 2.0) / (2.0 * q);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                this->dm_dsigma_buffer(i, j) -= 1.0 / (3.0 * 2.0 * q);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                this->dm_dsigma_buffer(i, j) -= md(i) * md(j) / q;
        return this->dm_dsigma_buffer;
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
    using parameters_t = std::tuple<DP_etabar>;

private:

    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type

};


// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): analytic closest-point derivatives available.
template<class AlphaHardeningType, class EtaHardeningType>
struct pf_has_cp_derivatives<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType>> : std::true_type {};

#endif
