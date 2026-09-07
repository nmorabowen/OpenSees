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

// Implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial3D - Hoek-Brown Yield Function
//
// Implements the Hoek-Brown failure criterion for rock masses
// Based on: Hoek, E., Brown, E.T. (2018). The Hoek-Brown failure criterion 
// and GSI – 2018 edition. Journal of Rock Mechanics and Geotechnical Engineering.
//
// The generalized Hoek-Brown criterion is:
//   σ₁ = σ₃ + σci * (mb * σ₃/σci + s)^a
//
// Where:
//   σ₁, σ₃ = major and minor principal stresses (compression positive)
//   σci = unconfined compressive strength of intact rock
//   mb = reduced material constant for rock mass = mi * exp((GSI-100)/(28-14D))
//   s = rock mass parameter = exp((GSI-100)/(9-3D))
//   a = rock mass parameter = 0.5 + (1/6)*(exp(-GSI/15) - exp(-20/3))
//   GSI = Geological Strength Index (0-100)
//   D = Disturbance factor (0-1)
//   mi = material constant for intact rock
//
// Sign convention: Compression is POSITIVE (geomechanics convention)

#ifndef HoekBrown_YF_H
#define HoekBrown_YF_H

#include "../YieldFunctionBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include <cmath>
#include <iostream>
#include <algorithm>

template<class NO_HARDENING>
class HoekBrown_YF : public YieldFunctionBase<HoekBrown_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "HoekBrown_YF";

    HoekBrown_YF():
        YieldFunctionBase<HoekBrown_YF<NO_HARDENING>>::YieldFunctionBase()
    {}

    // YIELD_FUNCTION
    // {
    //     using namespace std;
    //
    //     // Get Hoek-Brown parameters
    //     double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);  // Unconfined compressive strength
    //     double mb = GET_PARAMETER_VALUE(HB_mb);           // Reduced material constant
    //     double s = GET_PARAMETER_VALUE(HB_s);             // Rock mass parameter s
    //     double a = GET_PARAMETER_VALUE(HB_a);             // Rock mass parameter a
    //
    //     // Get principal stresses (sorted: sigma1 >= sigma2 >= sigma3)
    //     // Note: In this code, compression is positive
    //     VoigtVector sigma_geo = -sigma;
    //     auto [sigma3, sigma2, sigma1] = sigma_geo.principalStresses();
    //
    //     // The Hoek-Brown criterion in compression is:
    //     // f = σ₁ - σ₃ - σci * (mb * σ₃/σci + s)^a
    //     // 
    //     // We need to handle the case when mb * σ₃/σci + s < 0
    //     // which can happen in tension
    //
    //     double arg = mb * sigma3 / sigma_ci + s;
    //
    //     double yf;
    //     if (arg > 0) {
    //         // Standard Hoek-Brown
    //         yf = sigma1 - sigma3 - sigma_ci * pow(arg, a);
    //     } else {
    //         // In tension region - use tension cutoff behavior
    //         // The tensile strength is: σt = -s * σci / mb (for a = 0.5)
    //         // For general a, we use the limit as arg -> 0
    //         // When arg <= 0, we're beyond the tensile cutoff
    //         // Return a large positive value to indicate yielding
    //         double sigma_t = -s * sigma_ci / mb;  // Tensile strength (negative)
    //         yf = sigma1 - sigma3 - sigma_ci * s;  // Essentially tension cutoff
    //         // Alternative: linear extrapolation beyond tension cutoff
    //         // yf = sigma1 - sigma3;  // Pure shear criterion
    //     }
    //
    //     return yf;
    //

    // Ladruno (ADR-94 wp/94d): ported jaabell/ASDP 60d9b9b23 composite HB
    // yield function. Replaces the discontinuous if/else above (dead here,
    // preserved as a historical comment) with a continuous
    // max(f_shear, f_tension) composite -- see
    // Ladruno_implementation/_adr94_hb_drift.md for the derivation and the
    // measured tension-plateau defect this fixes (H10a / M4).
    YIELD_FUNCTION
    {
        using namespace std;

        // Hoek–Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb       = GET_PARAMETER_VALUE(HB_mb);
        double s        = GET_PARAMETER_VALUE(HB_s);
        double a        = GET_PARAMETER_VALUE(HB_a);

        // Compression-positive convention; σ1 ≥ σ2 ≥ σ3
        VoigtVector sigma_geo = -sigma;
        auto [sigma3, sigma2, sigma1] = sigma_geo.principalStresses();

        // Rock-mass tensile strength (negative in compression-positive sense)
        double sigma_t = -s * sigma_ci / mb;

        // arg = mb·σ3/σci + s.  arg > 0 ⇔ σ3 > σt (shear regime).
        // Clamp to 0 before pow() so the negative-arg case doesn't NaN.
        double arg      = mb * sigma3 / sigma_ci + s;
        double arg_safe = max(arg, 0.0);

        // Hoek–Brown shear surface
        double f_shear   = sigma1 - sigma3 - sigma_ci * pow(arg_safe, a);

        // Tension cut-off surface (Rankine-style on σ3)
        double f_tension = sigma_t - sigma3;

        // Composite: material yields if EITHER surface is violated.
        // Continuous at the apex: σ3 = σt ⇒ f_shear = σ1 − σt, f_tension = 0,
        // and on the HB envelope at the apex σ1 = σt so the composite is 0.
        return max(f_shear, f_tension);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);
        double a = GET_PARAMETER_VALUE(HB_a);
        double ds = GET_PARAMETER_VALUE(HB_ds);  // Perturbation for numerical derivative

        double sigma_norm = sigma.norm();
        ds = std::max(ds, ds * sigma_norm);

        // Use numerical differentiation for robustness
        // The analytical derivative involves principal stress directions which 
        // can be numerically sensitive
        auto computeNumericalDerivative = [this, &internal_variables_storage, &parameters_storage]
            (const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double yf1 = YF(SIG1);
                double yf2 = YF(SIG2);

                result(i) = (yf1 - yf2) / (2.0 * perturbation);
            }
            return result;
        };

        if (ds > 0) {
            vv_out = computeNumericalDerivative(sigma, ds);
        } else {
            // Default to numerical with small perturbation
            vv_out = computeNumericalDerivative(sigma, 1e-8 * std::max(1.0, sigma_norm));
        }

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // This model does not support hardening (perfectly plastic)
        return 0.0;
    }

    // Ladruno (ADR-94 wp/94d): kept live from our pre-port tree -- jaabell's
    // 60d9b9b23 carries this apex plumbing as a dead comment (below,
    // preserved), but the ASDPlasticMaterial3D.h apex call site is dead
    // code on both trees today (see _adr94_hb_drift.md (a)), so keeping it
    // live here is cosmetic-only and costs nothing.
    CHECK_APEX_REGION
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);
        double a = GET_PARAMETER_VALUE(HB_a);

        // The apex of the Hoek-Brown criterion in principal stress space
        // occurs at the tensile strength point
        // σt = -s * σci / mb (for a = 0.5)
        // For general a, it's approximately the same

        (void) a;

        // Ladruno (ADR-94 wp/94c, B4): this test was INVERTED and its sign
        // convention was wrong, which did not matter while the integrator's apex
        // call site was dead code.  It read `p = I1/3` (which is TENSION-positive
        // here, not compression-positive as its comment claimed) and fired when
        // p < -s*sigma_ci/mb, i.e. deep in hydrostatic COMPRESSION -- the opposite
        // half-axis from the Hoek-Brown tensile apex -- and APEX_STRESS then
        // returned an interior point.  Now that the call site is live, that pair
        // would have replaced a legitimate compressive state with a near-zero
        // stress on every compression deck.
        //
        // Near the apex the composite yield function is its Rankine tension cut-off
        // branch, f_tension = sigma_1(tension-positive) - T, with rock-mass tensile
        // strength T = s*sigma_ci/mb > 0.  That is a three-plane corner whose vertex
        // is sigma = T*I; the normal cone at a Rankine vertex is the positive octant
        // in principal space, so the trial state projects onto the vertex exactly
        // when ALL THREE tension-positive principal stresses are >= T, i.e. when the
        // LARGEST compression-positive principal is <= -T.
        if (!(mb > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return false;

        double T = s * sigma_ci / mb;          // tensile strength, TENSION-positive

        VoigtVector sigma_geo = -sigma;        // compression-positive
        auto [pr3, pr2, pr1] = sigma_geo.principalStresses();   // ascending
        (void) pr3;
        (void) pr2;

        return pr1 <= -T;
    }

    APEX_STRESS
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);

        // Ladruno (ADR-94 wp/94c, B4): sign fixed.  The rock-mass tensile strength
        // is T = s*sigma_ci/mb, and this class stores stress TENSION-POSITIVE, so
        // the hydrostatic-tension apex is +T on the three normal slots.  The old
        // -s*sigma_ci/mb put the "apex" in hydrostatic COMPRESSION, where the
        // composite yield function evaluates to max(-sigma_ci*(2s)^a, -2T) < 0 --
        // an interior point, not the apex.  With +T both branches vanish:
        // f_tension = T - T = 0 and f_shear = 0 - sigma_ci*0^a = 0.
        double sigma_t = s * sigma_ci / mb;

        // Return hydrostatic stress state at tensile strength
        // (all normal stresses equal, no shear)
        vv_out = VoigtVector(sigma_t, sigma_t, sigma_t, 0, 0, 0);

        return vv_out;
    }

    // jaabell/ASDP 60d9b9b23 carries the block above as a dead comment; not
    // reproduced verbatim here since we keep it live instead (see marker).

    using internal_variables_t = std::tuple<NO_HARDENING>;

    // Ladruno (ADR-94 wp/94c, M5): f has the units of stress and its natural
    // magnitude is the rock-mass uniaxial compressive strength sigma_ci * s^a
    // (5350 kPa for the ADR-94 sigma_ci = 50 MPa fixture -- nine orders above the
    // 1e-6 absolute default, which is exactly the M5 defect).
    YF_STRENGTH_SCALE
    {
        (void) internal_variables_storage;
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double s        = GET_PARAMETER_VALUE(HB_s);
        double a        = GET_PARAMETER_VALUE(HB_a);
        double sc = sigma_ci * std::pow(s > 0 ? s : 0.0, a);
        return sc < 0 ? -sc : sc;
    }

    using parameters_t = std::tuple<HB_sigci, HB_mb, HB_s, HB_a, HB_ds>;

private:
    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Declares this YF as featuring an apex
// Ladruno (ADR-94 wp/94d): kept live (see CHECK_APEX_REGION marker above).
template<class NO_HARDENING>
struct yf_has_apex<HoekBrown_YF<NO_HARDENING>> : std::true_type {};

#endif
