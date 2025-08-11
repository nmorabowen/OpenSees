#ifndef DruckerPrager_YF_H
#define DruckerPrager_YF_H

#include "../YieldFunctionBase.h"
#include <cmath>
#include <iostream>

template<class AlphaHardeningType, class KHardeningType>
class DruckerPrager_YF : public YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_YF";

    DruckerPrager_YF():
        YieldFunctionBase<DruckerPrager_YF<AlphaHardeningType, KHardeningType>>::YieldFunctionBase()
    {}

    // Yield function
    YIELD_FUNCTION
    {
        // p is positive in compression (soil mechanics sign convention)
        double p = -sigma.meanStress();
        auto s = sigma.deviator();

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        // Numerical guards
        constexpr double eps_p = 1e-14;
        constexpr double eps_k = 1e-14;

        double p_safe = (std::abs(p) < eps_p) ? ( (p >= 0.0) ? eps_p : -eps_p ) : p;
        double k_safe = (std::abs(k.value()) < eps_k) ? eps_k : k.value();

        // Deviatoric measure with kinematic shift
        VoigtVector s_shift = s - p_safe * alpha;

        // Guard against small negative inside sqrt due to round-off
        double tmp = s_shift.dot(s_shift);
        if (tmp < 0.0 && tmp > -1e-16) tmp = 0.0;

        double yf = std::sqrt(tmp) - (SQRT_2_over_3 * k_safe * p_safe);

        return yf;
    }

    // df/dsigma
    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        double p = -sigma.meanStress();
        auto s = sigma.deviator();

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);

        constexpr double eps_p = 1e-14;
        constexpr double eps_k = 1e-14;

        double p_safe = (std::abs(p) < eps_p) ? ( (p >= 0.0) ? eps_p : -eps_p ) : p;
        double k_safe = (std::abs(k.value()) < eps_k) ? eps_k : k.value();

        // Stress ratio r = s/p
        VoigtVector r = s / p_safe;

        // Denominator for normalization
        double den = (SQRT_2_over_3 * k_safe);

        // n = (r - alpha) / den
        auto n = (r - alpha) / den;

        double nr = n.dot(r);

        vv_out = n - nr * kronecker_delta() / 3.0;

        return vv_out;
    }

    // Hardening term df/dq * dq/d(lambda)
    YIELD_FUNCTION_HARDENING
    {
        double dbl_result = 0.0;

        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        auto k = GET_TRIAL_INTERNAL_VARIABLE(KHardeningType);
        
        double p = -sigma.meanStress();
        auto s = sigma.deviator();

        constexpr double eps_p = 1e-14;
        constexpr double eps_k = 1e-14;

        double p_safe = (std::abs(p) < eps_p) ? ( (p >= 0.0) ? eps_p : -eps_p ) : p;
        double k_safe = (std::abs(k.value()) < eps_k) ? eps_k : k.value();

        double den = (SQRT_2_over_3 * k_safe);

        // Isotropic hardening part (k)
        double df_dk = -SQRT_2_over_3 * p_safe;
        dbl_result += (df_dk * GET_INTERNAL_VARIABLE_HARDENING(KHardeningType)).value();

        // Kinematic hardening part (alpha)
        auto df_dalpha = (p_safe * alpha - s) / den;
        dbl_result += df_dalpha.dot(GET_INTERNAL_VARIABLE_HARDENING(AlphaHardeningType));

        return dbl_result;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;
    using parameters_t         = std::tuple<>;

private:
    static VoigtVector vv_out; // For returning VoigtVectors
};

// Static member definition
template <class AlphaHardeningType, class KHardeningType>
VoigtVector DruckerPrager_YF<AlphaHardeningType, KHardeningType>::vv_out;

// Declares this YF as featuring an apex
template<class AlphaHardeningType, class KHardeningType>
struct yf_has_apex<DruckerPrager_YF<AlphaHardeningType, KHardeningType>> : std::true_type {};

#endif
