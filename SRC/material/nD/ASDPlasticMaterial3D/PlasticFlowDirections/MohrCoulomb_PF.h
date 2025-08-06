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

#ifndef MohrCoulomb_PF_H
#define MohrCoulomb_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include <cmath>
#include <typeinfo>

template<class NO_HARDENING>
class MohrCoulomb_PF : public PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>// CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_PF";

    MohrCoulomb_PF( ):
        PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>::PlasticFlowBase()
                { }

    double g(const VoigtVector& sigma, double phi) const
    {
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        return rEquivalentStress;
    }

    PLASTIC_FLOW_DIRECTION
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double ds = GET_PARAMETER_VALUE(MC_ds);
        double D = GET_PARAMETER_VALUE(Dilatancy);

        using namespace std;

        double sigma_norm = sigma.norm();

        // Perturbation to smooth the PF
        ds = std::max(ds, ds*sigma_norm);

        // If the perturbation is set greater than zero, use numerical differentiation to get plastic flow direction
        if (ds > 0)
        {
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sigma;
                VoigtVector SIG2 = sigma;
                
                // Increment SIG at index i by a small amount for numerical differentiation
                SIG1(i) += ds;
                SIG2(i) -= ds;

                // Compute the plastic potential function at the perturbed state
                double g1 = g(SIG1, phi);
                double g2 = g(SIG2, phi);

                // Calculate the derivative
                vv_out(i) = (g1 - g2) / (2*ds);
            }
        } 
        else // Use analytic solution (consistent with YF implementation)
        {
            VoigtVector first_vector = calculate_first_vector();
            VoigtVector second_vector = calculate_second_vector(sigma);
            VoigtVector third_vector = calculate_third_vector(sigma);

            double J2 = sigma.getJ2();
            double J3 = sigma.getJ3();
            double lode_angle = sigma.lodeAngle();

            double c1, c3, c2;
            double checker = std::abs(lode_angle * 180.0 / M_PI);

            if (std::abs(checker) < 29.0) { // If it is not the edge
                c1 = std::sin(phi) / 3.0;
                c3 = (std::sqrt(3.0) * std::sin(lode_angle) + std::sin(phi) * std::cos(lode_angle)) /
                    (2.0 * J2 * std::cos(3.0 * lode_angle));
                c2 = 0.5 * std::cos(lode_angle)*(1.0 + std::tan(lode_angle) * std::sin(3.0 * lode_angle) +
                    std::sin(phi) * (std::tan(3.0 * lode_angle) - std::tan(lode_angle)) / std::sqrt(3.0));
            } else { // smoothing with drucker-prager
                c1 = 3.0 * (2.0 * std::sin(phi) / (std::sqrt(3.0) * (3.0 - std::sin(phi))));
                c2 = 1.0;
                c3 = 0.0;
            }

            vv_out = c1 * first_vector + c2 * second_vector + c3 * third_vector;
        }

        // Apply dilatancy correction: modify the volumetric component
        // For non-associated flow, we modify the mean stress contribution
        // The flow direction becomes: dg/dsigma - D * delta_ij/3
        // where D is the dilatancy parameter
        
        // Alternative implementation: modify only the volumetric part
        // Extract deviatoric and volumetric components
        VoigtVector volumetric_part = kronecker_delta() / 3.0;
        
        // Apply dilatancy modification to the volumetric component
        // For D=0: associated flow, D<sin(phi): non-associated flow with reduced dilation
        double volumetric_multiplier = std::sin(phi) - D;
        vv_out = vv_out.deviator() + volumetric_multiplier * volumetric_part;

        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds,Dilatancy>;

private:

    static VoigtVector vv_out; 
};

template<class NO_HARDENING>
VoigtVector MohrCoulomb_PF<NO_HARDENING>::vv_out;

#endif