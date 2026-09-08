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

#ifndef VonMises_YF_H
#define VonMises_YF_H

// #include "../EvolvingVariable.h"
#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>

#include "../ASDPlasticMaterial3DGlobals.h"
using namespace ASDPlasticMaterial3DGlobals;




template<class AlphaHardeningType, class KHardeningType>
class VonMises_YF : public YieldFunctionBase<VonMises_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "VonMises_YF";

    VonMises_YF( ):
        YieldFunctionBase<VonMises_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase() // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                {}

    YIELD_FUNCTION
    {
        auto s = sigma.deviator();

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        double tmp = tensor_dot_stress_like(s - alpha, s - alpha);  // J2 = 0.5 || s : s ||
        tmp = tmp > 0 ? tmp : 0;
        return std::sqrt( tmp ) - SQRT_2_over_3 * k.value() ;  // This one assumes p positive in tension
        //        sqrt(   3  / 2 * s:s  )   -    sigma_y 
        //        sqrt(    s:s  )   -   SQRT_2_over_3 * sigma_y 
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        result = sigma.deviator() - alpha;

        double den = sqrt(tensor_dot_stress_like(result, result));
        if (abs(den) > 100*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            result = result / den;

        // Ladruno (ADR-94 wp/94c, B5): VOIGT convention.
        //   f = sqrt(r_ij r_ij) with r held in raw (undoubled) Voigt storage, so
        //     df/d sigma_12 = r_12/den                (tensor derivative)
        //     df/d v_12     = 2*r_12/den              (Voigt derivative)
        // because the stored slot v_12 drives BOTH sigma_12 and sigma_21.  The
        // normal slots are identical under both conventions.
        // Every consumption site is a plain Voigt contraction (n^T E m), and both
        // `TrialPlastic_Strain += dLambda*m` and `Eelastic*m` require an
        // ENGINEERING-shear operand, so Voigt is the convention the framework
        // actually needs -- and the one MC/HB/MCTC/StiffSoil already return.
        // Returning the bare tensor derivative here (as this did) under-counted
        // shear by 2 in the yield normal, the flow direction and the plastic
        // strain update.  See ADR-94 B5.
        result(3) *= 2.0;   // v12
        result(4) *= 2.0;   // v23
        result(5) *= 2.0;   // v13

        return result;
    }
    
    YIELD_FUNCTION_HARDENING
    {
        double dbl_result = 0.0;

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        //Zero the stress deviator
        auto s = sigma.deviator();

        // This is for the hardening of k
        double df_dk = -SQRT_2_over_3;
        dbl_result +=  (df_dk * GET_INTERNAL_VARIABLE_HARDENING(KHardeningType)).value();

        //This is for the hardening of alpha
        double den = sqrt(tensor_dot_stress_like(s - alpha, s - alpha));

        if (abs(den) < sqrt(tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
        {
            return dbl_result;
        }

        // Ladruno (ADR-94 wp/94c, B5): df/dalpha in the SAME Voigt convention as
        // df_dsigma_ij above (shear slots doubled), contracted with the hardening
        // rate by a plain dot.  Numerically identical to the old tensor-convention
        // `tensor_dot_stress_like(df_dalpha, hh)` -- the factor 2 just moves from
        // the contraction into the operand, and scaling by 2 is exact in IEEE --
        // but it keeps ONE convention across the file.
        VoigtVector df_dalpha = -(s - alpha) / den;
        df_dalpha(3) *= 2.0;
        df_dalpha(4) *= 2.0;
        df_dalpha(5) *= 2.0;
        VoigtVector hh = GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType);
        dbl_result +=  df_dalpha.dot(hh);

        return dbl_result;
    }

    // Ladruno (ADR-97 wp/97b): df/dq, UNCONTRACTED, in the same VOIGT convention
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

    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(r:r) - sqrt(2/3)*k, so the term that
    // sets f's scale is sqrt(2/3)*sigma_y, with sigma_y the CURRENT yield stress
    // internal variable (not a model parameter for this YF).
    YF_STRENGTH_SCALE
    {
        (void) parameters_storage;
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);
        double kv = k.value();
        return SQRT_2_over_3 * (kv < 0 ? -kv : kv);
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;

    using parameters_t = std::tuple<>;


private:

    mutable VoigtVector result = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type

};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97b): von Mises supplies the analytic closest-point df/dq.
template<class AlphaHardeningType, class KHardeningType>
struct yf_has_cp_derivatives<VonMises_YF<AlphaHardeningType, KHardeningType>> : std::true_type {};


#endif
