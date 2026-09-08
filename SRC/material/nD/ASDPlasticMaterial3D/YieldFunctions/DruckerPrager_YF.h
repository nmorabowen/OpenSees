#ifndef DruckerPrager_YF_H
#define DruckerPrager_YF_H

#include "../YieldFunctionBase.h"
#include <cmath>
#include <iostream>
#include "../AllASDModelParameterTypes.h"

template<class AlphaHardeningType, class CohesionHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_YF";

    DruckerPrager_YF():
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>>::YieldFunctionBase()
    {}

    YIELD_FUNCTION
    {
        auto s = sigma.deviator();
        double p = sigma.meanStress();  // mean stress (positive in compression in geomechanics convention)
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        // auto eta = GET_TRIAL_INTERNAL_VARIABLE(CohesionHardeningType);
        
        // Get friction parameter (eta) and cohesion (c) from parameters
        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);      // cohesion
        double eta = GET_PARAMETER_VALUE(DP_eta);      // friction slope
        
        double tmp = tensor_dot_stress_like(s - alpha, s - alpha);
        tmp = tmp > 0 ? tmp : 0;
        double sqrt_J2 = std::sqrt(0.5*tmp);
        
        // Drucker-Prager yield function: sqrt(J2) + eta * p - (c + k)
        // where:
        // - 0.5*sqrt(J2) is the deviatoric stress magnitude
        // - eta is the friction parameter 
        // - p is the mean stress (positive in compression)
        // - xi_c is the adjusted cohesion
        return sqrt_J2 + eta * p - (xi_c);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        auto s = sigma.deviator();
        double p = sigma.meanStress();
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);

        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);      // cohesion
        double eta = GET_PARAMETER_VALUE(DP_eta);      // friction slope

        // Derivative with respect to deviatoric stress
        VoigtVector dev_part = s - alpha;
        double den = sqrt(0.5*tensor_dot_stress_like(dev_part, dev_part));
        
        if (abs(den) > 100*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            dev_part = dev_part / den;
        else
            dev_part.setZero();  // Ladruno (ADR-94 wp/94a): was `*= 0.0`; NaN*0 == NaN

        // Ladruno (ADR-94 wp/94c, B4/B5): sqrt(J2) = sqrt(0.5 * r_ij r_ij), so the
        // per-slot chain rule is NOT flat:
        //     d sqrt(J2) / d sigma_ij = r_ij / (2 sqrt(J2))      (tensor)
        //     d sqrt(J2) / d v_ii     = r_ii / (2 sqrt(J2))      (Voigt, normal)
        //     d sqrt(J2) / d v_ij     = r_ij /    sqrt(J2)       (Voigt, shear)
        // `r/den` (what the line above produces) is the SHEAR answer applied to all
        // six slots: exactly 2x too large on the three normal slots.  That is a
        // real, convention-independent gradient bug -- ADR-94 B4 measured 0.971
        // relative error against a central difference at a pure-normal point, and
        // it is the reason Drucker-Prager matched NEITHER the tensor nor the Voigt
        // finite difference.  Halving the normal slots leaves the gradient in the
        // Voigt convention used everywhere else in the framework.
        dev_part(0) *= 0.5;
        dev_part(1) *= 0.5;
        dev_part(2) *= 0.5;
            
        // Add pressure-dependent part: eta * dp/dsigma = eta/3 * I
        // Ladruno (ADR-94 wp/94a): ADR-94 B4 -- `VoigtVector x; x *= 0.0;` multiplies
        // UNINITIALISED Eigen storage by zero, which does not clear heap garbage that
        // decodes as NaN. This is where the Drucker-Prager hydrostatic-tension NaN was
        // born (and then committed with analyze() == 0). Same fix as the ADR-84 P0
        // constructor precedent in ASDPlasticMaterial3D.h.
        VoigtVector pressure_part = VoigtVector::Zero();
        pressure_part(0) = eta / 3.0;  // sigma_xx component
        pressure_part(1) = eta / 3.0;  // sigma_yy component  
        pressure_part(2) = eta / 3.0;  // sigma_zz component
        // shear components remain zero
        
        vv_out = dev_part + pressure_part;
        
        return vv_out;
    }
    
    YIELD_FUNCTION_HARDENING
    {
        double dbl_result = 0.0;
        
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto eta = GET_TRIAL_INTERNAL_VARIABLE(CohesionHardeningType);
        
        auto s = sigma.deviator();
        
        // Hardening contribution from eta (isotropic hardening)
        double df_deta = -1.0;  // derivative of f with respect to eta
        dbl_result += (df_deta * GET_INTERNAL_VARIABLE_HARDENING(CohesionHardeningType)).value();
        
        // Hardening contribution from alpha (kinematic hardening)
        double den = sqrt(0.5*tensor_dot_stress_like(s - alpha, s - alpha));
        
        if (abs(den) < sqrt(0.5*tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
        {
            return dbl_result;
        }
        
        // Ladruno (ADR-94 wp/94c, B4/B5): same two corrections as df_dsigma_ij --
        // d sqrt(J2)/d alpha carries the 1/2 on the normal slots (the old
        // `-(s-alpha)/den` was 2x too large there, so the kinematic-hardening term
        // of H was 2x too large for every Drucker-Prager with a tensor IV), and the
        // result is expressed in the Voigt convention so the contraction with the
        // hardening rate is a plain dot.
        VoigtVector df_dalpha = -(s - alpha) / den;
        df_dalpha(0) *= 0.5;
        df_dalpha(1) *= 0.5;
        df_dalpha(2) *= 0.5;
        VoigtVector hh = GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType);
        dbl_result += df_dalpha.dot(hh);
        
        return dbl_result;
    }


    // Ladruno (ADR-94 wp/94c, B4): implemented.  Both methods were `Implement!!!`
    // stubs (check returned false, apex returned the zero stress), which is why
    // Drucker-Prager hydrostatic tension ran the flank return map past sqrt(J2) = 0
    // and committed NaN.
    CHECK_APEX_REGION
    {
        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);
        double eta  = GET_PARAMETER_VALUE(DP_eta);

        if (!(eta > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return false;                 // no cone -> no apex

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        VoigtVector r = sigma.deviator() - alpha;
        double tmp = tensor_dot_stress_like(r, r);
        tmp = tmp > 0 ? tmp : 0;
        const double q = std::sqrt(0.5 * tmp);      // sqrt(J2)
        const double p = sigma.meanStress();        // TENSION-POSITIVE (trace/3)
        const double p_apex = xi_c / eta;

        // Cone geometry in the (p, q) half-plane, q >= 0:
        //   f = q + eta*p - xi_c = 0  =>  the yield surface is the straight line
        //   q = xi_c - eta*p, which terminates at the apex (p_apex, 0) with
        //   p_apex = xi_c/eta.  Beyond that point the flank return would have to
        //   drive sqrt(J2) NEGATIVE, which is what produced the NaN.
        //   The unit tangent at the apex pointing AWAY from the flank is
        //   t ~ (1, -eta), so a trial state projects onto the apex rather than
        //   onto the flank exactly when (P - apex) . t >= 0, i.e.
        //         p - p_apex >= eta * q.
        // HONEST CAVEAT: the exact condition in the ELASTIC metric is
        //         p - p_apex >= (K * etabar / G) * q,
        // which needs the bulk and shear moduli and the plastic potential's
        // dilatancy -- none of which this signature can see.  The Euclidean test
        // coincides with it when K*etabar/G == eta.  The integrator independently
        // validates f(sigma_apex) against the yield tolerance before committing an
        // apex return, so a misclassification cannot commit an inadmissible stress;
        // it only falls back to the generic return map.
        return (p - p_apex) >= eta * q;
    }

    APEX_STRESS
    {
        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);
        double eta  = GET_PARAMETER_VALUE(DP_eta);

        // Hydrostatic apex, TENSION-POSITIVE because meanStress() == trace/3:
        // q = 0 and f = eta*p - xi_c = 0  =>  p = xi_c/eta, giving f(sigma_apex)
        // identically 0.
        const double p_apex = (eta != 0.0) ? xi_c / eta : 0.0;
        vv_out = VoigtVector(p_apex, p_apex, p_apex, 0, 0, 0);

        return vv_out;
    }

    // Ladruno (ADR-97 wp/97b): df/dq, UNCONTRACTED, in the VOIGT convention
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

    // Ladruno (ADR-94 wp/94c, M5): f = sqrt(J2) + eta*p - xi_c, so xi_c (the
    // adjusted cohesion) is the term that sets f's scale.
    YF_STRENGTH_SCALE
    {
        (void) internal_variables_storage;
        double xi_c = GET_PARAMETER_VALUE(DP_xi_c);
        return xi_c < 0 ? -xi_c : xi_c;
    }


    using internal_variables_t = std::tuple<AlphaHardeningType, CohesionHardeningType>;
    using parameters_t         = std::tuple<DP_xi_c, DP_eta>;

private:
    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Static member definition
// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Declares this YF as featuring an apex
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};

// Ladruno (ADR-97 wp/97b): Drucker-Prager supplies the analytic closest-point
// df/dq.  NOTE that `check_apex_region` above stays EUCLIDEAN and is used only
// by Backward_Euler: `Closest_Point` classifies the apex region in the ELASTIC
// metric, inside the integrator where K, G and etabar are in scope (ADR-97
// "Drucker-Prager apex").  Two integrators, two answers for the same YF -- a
// documentation obligation (LEDGER_quirks), not a bug.
template<class AlphaHardeningType, class CohesionHardeningType>
struct yf_has_cp_derivatives<DruckerPrager_YF<AlphaHardeningType, CohesionHardeningType>> : std::true_type {};

#endif
