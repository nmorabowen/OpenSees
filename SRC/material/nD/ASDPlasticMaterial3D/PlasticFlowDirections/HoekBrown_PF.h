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
// ASDPlasticMaterial3D - Hoek-Brown Plastic Flow Direction
//
// Implements a non-associated plastic flow rule for the Hoek-Brown criterion.
// The plastic potential uses a similar form to the yield function but with
// a different (typically lower) dilation parameter to control volumetric
// plastic strain.
//
// The plastic potential is:
//   g = σ₁ - σ₃ - σci * (mb_psi * σ₃/σci + s)^a
//
// Where mb_psi = mi_psi * exp((GSI-100)/(28-14D)) with mi_psi typically
// corresponding to a dilation angle ψ instead of friction angle φ.
//
// For rock masses, non-associated flow is essential because associated
// flow would predict unrealistically large dilation.
//
// Sign convention: Compression is POSITIVE (geomechanics convention)

#ifndef HoekBrown_PF_H
#define HoekBrown_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include <cmath>
#include <typeinfo>

template<class NO_HARDENING>
class HoekBrown_PF : public PlasticFlowBase<HoekBrown_PF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "HoekBrown_PF";

    HoekBrown_PF():
        PlasticFlowBase<HoekBrown_PF<NO_HARDENING>>::PlasticFlowBase()
    {}

    // Plastic potential function g (similar to yield function but with dilation parameter)
    double g(const VoigtVector& sig, double sigma_ci, double mb_psi, double s, double a) const
    {
        using namespace std;
        
        // Get principal stresses (sorted: sigma1 >= sigma2 >= sigma3)
        auto [sigma3, sigma2, sigma1] = sig.principalStresses();
        
        double arg = mb_psi * sigma3 / sigma_ci + s;
        
        double gval;
        if (arg > 0) {
            gval = sigma1 - sigma3 - sigma_ci * pow(arg, a);
        } else {
            // In tension region - use tension cutoff behavior
            gval = sigma1 - sigma3 - sigma_ci * s;
        }
        
        return gval;
    }

    PLASTIC_FLOW_DIRECTION
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb_psi = GET_PARAMETER_VALUE(HB_mb_psi);  // Dilation parameter (reduced mb)
        double s = GET_PARAMETER_VALUE(HB_s);
        double a = GET_PARAMETER_VALUE(HB_a);
        double ds = GET_PARAMETER_VALUE(HB_ds);

        double sigma_norm = sigma.norm();
        ds = std::max(ds, ds * sigma_norm);

        // Use numerical differentiation of the plastic potential
        auto computeNumericalFlowDirection = [this, sigma_ci, mb_psi, s, a, &internal_variables_storage, &parameters_storage]
            (const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double g1 = this->g(SIG1, sigma_ci, mb_psi, s, a);
                double g2 = this->g(SIG2, sigma_ci, mb_psi, s, a);

                result(i) = (g1 - g2) / (2.0 * perturbation);
            }
            return result;
        };

        if (ds > 0) {
            vv_out = computeNumericalFlowDirection(sigma, ds);
        } else {
            // Default to numerical with small perturbation
            vv_out = computeNumericalFlowDirection(sigma, 1e-8 * std::max(1.0, sigma_norm));
        }

        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    // Ladruno (ADR-97 wp/97d): the principal-space potential constants.  See the
    // macro's comment in PlasticFlowBase.h: `Closest_Point` builds a
    // FRAME-CONSISTENT Hoek-Brown potential from these and does NOT call `g`
    // above, whose `arg` is built from the tree's most COMPRESSIVE principal and
    // is therefore negative on every compressive state -- making `g` a Tresca
    // potential with HB_mb_psi inert.  `g` is left exactly as shipped so
    // Backward_Euler stays byte-identical (ADR-97 D1).
    CP_PRINCIPAL_HB_FLOW_PARAMS
    {
        (void) internal_variables_storage;
        sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        mb_psi   = GET_PARAMETER_VALUE(HB_mb_psi);
        s_hb     = GET_PARAMETER_VALUE(HB_s);
        a_hb     = GET_PARAMETER_VALUE(HB_a);
        return true;
    }

    using parameters_t = std::tuple<HB_sigci, HB_mb_psi, HB_s, HB_a, HB_ds>;

private:
    mutable VoigtVector vv_out = VoigtVector(0., 0., 0., 0., 0., 0.);  // Ladruno (ADR-94 wp/94b, F2): was a class-static return buffer, shared by every material that reuses this functor type
};

// Ladruno (ADR-94 wp/94b, F2): out-of-class static definition removed;
// the return buffer is a per-instance member now.

// Ladruno (ADR-97 wp/97d): principal-stress-space closest-point family 3
// (Hoek-Brown).  Matched against yf_cp_principal_family in HoekBrown_YF.h.
template<class NO_HARDENING>
struct pf_cp_principal_family<HoekBrown_PF<NO_HARDENING>>
    : std::integral_constant<int, 3> {};

#endif
