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

#ifndef MohrCoulomb_YF_H
#define MohrCoulomb_YF_H


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"


template<class NO_HARDENING>
class MohrCoulomb_YF : public YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_YF";


    MohrCoulomb_YF( ):
        YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>>::YieldFunctionBase() 
        {}

    YIELD_FUNCTION 
    {
        using namespace std;

        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);

        double rThresh = c * cos(phi);
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        double yf = rEquivalentStress - rThresh;

        return yf;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);

        using namespace std;

        double sigma_norm = sigma.norm();

        // Perturbation to smooth the YF
        ds = std::max(ds, ds*sigma_norm);

        // Helper lambda for numerical differentiation
        auto computeNumericalDerivative = [this, &internal_variables_storage, &parameters_storage](const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double yf1 = YF(SIG1);
                double yf2 = YF(SIG2);

                result(i) = (yf1 - yf2) / (2*perturbation);
            }
            return result;
        };

        // If perturbation is set, use numerical differentiation
        if (ds > 0) {
            vv_out = computeNumericalDerivative(sigma, ds);
        } 
        else {
            // Try analytical solution with fallback to numerical
            bool useNumerical = false;
            
            try {
                double J2 = sigma.getJ2();
                
                // Check for numerical issues - use simplified approach for hydrostatic states
                if (J2 < 1e-15) {
                    vv_out = std::sin(phi) / 3.0 * calculate_first_vector();
                } else {
                    VoigtVector first_vector = calculate_first_vector();
                    VoigtVector second_vector = calculate_second_vector(sigma);
                    VoigtVector third_vector = calculate_third_vector(sigma);

                    double lode_angle = sigma.lodeAngle();
                    double c1, c2, c3;
                    double checker = std::abs(lode_angle * 180.0 / M_PI);

                    if (std::abs(checker) < 29.0) { // Regular case
                        c1 = std::sin(phi) / 3.0;
                        
                        double denominator = 2.0 * J2 * std::cos(3.0 * lode_angle);
                        if (std::abs(denominator) < 1e-15) {
                            useNumerical = true; // Division by zero risk
                        } else {
                            c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(phi) * std::cos(lode_angle)) / denominator;
                            c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                                std::sin(phi) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
                        }
                    } else { // Edge smoothing with Drucker-Prager
                        c1 = 3.0 * (2.0 * std::sin(phi) / (std::sqrt(3.0) * (3.0 - std::sin(phi))));
                        c2 = 1.0;
                        c3 = 0.0;
                    }

                    if (!useNumerical) {
                        vv_out = c1 * first_vector + c2 * second_vector + c3 * third_vector;

                        // Validate result - check for NaN/Inf
                        for (int i = 0; i < 6; ++i) {
                            if (!std::isfinite(vv_out(i))) {
                                useNumerical = true;
                                break;
                            }
                        }
                    }
                }
            }
            catch (...) {
                useNumerical = true;
            }

            // Fallback to numerical if analytical failed
            if (useNumerical) {
                vv_out = computeNumericalDerivative(sigma, 1e-6);
            }
        }

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // This model does not support hardening 
        return 0.0;
    }


    // Ladruno (ADR-94 wp/94c, B4): this test was written when the integrator's
    // apex call site was dead code, so its eagerness never mattered.  Now that the
    // call site is live it decides which trial states are projected onto the apex,
    // so it is the proper NORMAL-CONE test rather than the bare `p beyond p_apex`.
    CHECK_APEX_REGION
    {
        using namespace std;

        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c   = GET_PARAMETER_VALUE(MC_c);

        double sphi = sin(phi);
        if (!(sphi > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return false;         // phi == 0 is Tresca: a cylinder, with no apex

        // f = A(theta)*sqrt(J2) + sin(phi)*p - c*cos(phi), with p TENSION-POSITIVE
        // (p = I1/3) and A(theta) = cos(theta) - sin(theta)*sin(phi)/sqrt(3) > 0 for
        // every admissible Lode angle.  In the (p, r = sqrt(J2)) half-plane the
        // surface is the line r = (c*cos(phi) - sin(phi)*p)/A, which terminates at
        // the apex (p_apex, 0), p_apex = c*cos(phi)/sin(phi) = c/tan(phi).  A trial
        // state projects onto the apex rather than onto the flank exactly when
        //         p - p_apex >= (sin(phi)/A) * r
        // (the tangent-perpendicular test at the endpoint of the flank).  The old
        // form kept only the p term, which classified deviatorically loaded states
        // near the apex as apex states.  Same elastic-metric caveat as
        // DruckerPrager_YF::check_apex_region; the integrator validates
        // f(sigma_apex) before committing.
        double A = cos(sigma.lodeAngle()) - sin(sigma.lodeAngle()) * sphi / sqrt(3.0);
        if (!(A > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return false;

        double p = sigma.getI1() / 3.0;
        double r = sqrt(sigma.getJ2() > 0 ? sigma.getJ2() : 0.0);
        double p_apex = c * cos(phi) / sphi;

        return (p - p_apex) >= (sphi / A) * r;
    }

    // Ladruno (ADR-94 wp/94c, M5): f = A*sqrt(J2) + sin(phi)*p - c*cos(phi), so
    // c*cos(phi) is the term that sets f's scale (it degenerates to c for phi = 0).
    YF_STRENGTH_SCALE
    {
        (void) internal_variables_storage;
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c   = GET_PARAMETER_VALUE(MC_c);
        double sc  = c * std::cos(phi);
        return sc < 0 ? -sc : sc;
    }

    // Ladruno (ADR-97 wp/97c): the principal-space face constants.  Both are
    // read from the SAME MC_phi / MC_c parameters this functor's own `f` uses,
    // so the closest-point map returns to exactly the surface `f` measures.
    CP_PRINCIPAL_MC_FACE_PARAMS
    {
        (void) internal_variables_storage;
        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double c   = GET_PARAMETER_VALUE(MC_c);
        sin_phi = std::sin(phi);
        k_coh   = c * std::cos(phi);
        return true;
    }

    APEX_STRESS
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double p_apex = c / tan(phi);

        vv_out = VoigtVector(p_apex,p_apex,p_apex,0,0,0);

        return vv_out;
    }

  
    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds>;

private:


    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<MohrCoulomb_YF<NO_HARDENING>> : std::true_type {};

// Ladruno (ADR-97 wp/97c): principal-stress-space closest-point family 1
// (plain Mohr-Coulomb).  Paired only with MohrCoulomb_PF, which carries the
// same marker.
template<class NO_HARDENING>
struct yf_cp_principal_family<MohrCoulomb_YF<NO_HARDENING>>
    : std::integral_constant<int, 1> {};

#endif
