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
//
// Ladruno (ADR-84 P0): composite MohrCoulomb + tension-cutoff plastic flow direction (fork addition).

#ifndef MohrCoulombTensionCutoff_PF_H
#define MohrCoulombTensionCutoff_PF_H


#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include <cmath>
#include <typeinfo>
#include <Eigen/Eigenvalues>   // Defensive: SelfAdjointEigenSolver for the Rankine branch

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

template<class NO_HARDENING>
class MohrCoulombTensionCutoff_PF : public PlasticFlowBase<MohrCoulombTensionCutoff_PF<NO_HARDENING>>// CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulombTensionCutoff_PF";

    MohrCoulombTensionCutoff_PF( ):
        PlasticFlowBase<MohrCoulombTensionCutoff_PF<NO_HARDENING>>::PlasticFlowBase()
                { }

    // Mohr-Coulomb plastic potential (same arithmetic as MohrCoulomb_YF's yield function)
    double g(const VoigtVector& sigma, double phi, double c) const
    {
        using namespace std;

        double rThresh = c * cos(phi);
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        double yf = rEquivalentStress - rThresh;

        return yf;
    }

    // Tension-cutoff (Rankine) potential: largest principal stress minus the cutoff
    // Same arithmetic as TensionCutoff_YF's yield function.
    double g_tc(const VoigtVector& sigma, double min_stress) const
    {
        double s1, s2, s3;
        std::tie(s1, s2, s3) = sigma.principalStresses();  //Ordered (ascending), s3 is the largest

        return s3 - min_stress;
    }

    PLASTIC_FLOW_DIRECTION
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double ds = GET_PARAMETER_VALUE(MC_ds);
        double psi = GET_PARAMETER_VALUE(MC_psi)*M_PI/180;
        // NOTE: MohrCoulomb_PF.h scales the cohesion by M_PI/180 (a units bug, benign there
        // because c is additive in g and cancels in the finite differences). Deliberately not
        // replicated here: c is used as-is so the branch test against f_MC is dimensionally right.
        double c = GET_PARAMETER_VALUE(MC_c);
        double min_stress = GET_PARAMETER_VALUE(TC_min_stress);

        using namespace std;

        // ---------------------------------------------------------------------
        // Branch selection: whichever surface is more violated (or less far from
        // being violated) at the current stress governs the flow direction.
        // ---------------------------------------------------------------------
        double f_MC = g(sigma, phi, c);              // MohrCoulomb_YF arithmetic
        double f_TC = g_tc(sigma, min_stress);       // TensionCutoff_YF arithmetic

        if (f_TC >= f_MC)
        {
            // -----------------------------------------------------------------
            // Tension-cutoff branch: associated (Rankine) flow, m = v3 (x) v3
            // with v3 the eigenvector of the largest principal stress.
            // -----------------------------------------------------------------
            double s1, s2, s3;
            std::tie(s1, s2, s3) = sigma.principalStresses();

            bool degenerate = (s3 - s2) < 1e-8 * std::max(1.0, std::abs(s3));

            if (!degenerate)
            {
                Eigen::Matrix3d stress_tensor;
                stress_tensor << sigma.v11(), sigma.v12(), sigma.v13(),
                                 sigma.v12(), sigma.v22(), sigma.v23(),
                                 sigma.v13(), sigma.v23(), sigma.v33();

                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_tensor);

                // Eigenvalues (and their eigenvectors) come out in ascending order,
                // so column 2 belongs to the largest principal stress.
                Eigen::Vector3d v = eigen_solver.eigenvectors().col(2);

                double v0 = v(0);
                double v1 = v(1);
                double v2 = v(2);

                // Voigt assembly [11,22,33,12,23,13] with factor-2 shears (strain-like),
                // matching the convention of calculate_second_vector().
                vv_out(0) = v0 * v0;
                vv_out(1) = v1 * v1;
                vv_out(2) = v2 * v2;
                vv_out(3) = 2.0 * v0 * v1;
                vv_out(4) = 2.0 * v1 * v2;
                vv_out(5) = 2.0 * v0 * v2;

                // Validate result - check for NaN/Inf
                for (int i = 0; i < 6; ++i) {
                    if (!std::isfinite(vv_out(i))) {
                        degenerate = true;
                        break;
                    }
                }
            }

            if (degenerate)
            {
                // Repeated maximum principal stress: the eigenvector is not unique.
                // Fall back to central finite differences of the Rankine potential,
                // mirroring TensionCutoff_YF's stress-derivative loop.
                double ds_tc = ds == 0 ? 1e-8 : ds;

                for (int i = 0; i < 6; ++i) {
                    VoigtVector SIG1 = sigma;
                    VoigtVector SIG2 = sigma;

                    SIG1(i) += ds_tc;
                    SIG2(i) -= ds_tc;

                    double g1 = g_tc(SIG1, min_stress);
                    double g2 = g_tc(SIG2, min_stress);

                    vv_out(i) = (g1 - g2) / (2*ds_tc);
                }
            }

            // Associated flow on the cutoff: no dilatancy correction here.
            return vv_out;
        }

        // ---------------------------------------------------------------------
        // Mohr-Coulomb branch (verbatim MohrCoulomb_PF logic)
        // ---------------------------------------------------------------------

        double sigma_norm = sigma.norm();

        // Perturbation to smooth the PF
        ds = std::max(ds, ds*sigma_norm);

        // Helper lambda for numerical differentiation of plastic potential
        auto computeNumericalFlowDirection = [this, phi, c, &internal_variables_storage, &parameters_storage](const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;

                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double g1 = g(SIG1, phi, c);
                double g2 = g(SIG2, phi, c);

                if (i < 3)
                {
                    result(i) = (g1 - g2) / (2*perturbation);
                }
                else
                {
                    // result(i) = 2 * (g1 - g2) / (2*perturbation);   // Engineering strain
                    result(i) = (g1 - g2) / (2*perturbation);   // Engineering strain
                }
            }
            return result;
        };

        // If perturbation is set, use numerical differentiation
        if (ds > 0) {
            vv_out = computeNumericalFlowDirection(sigma, ds);
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
                vv_out = computeNumericalFlowDirection(sigma, 1e-6);
            }
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
        double volumetric_multiplier = std::sin(psi);
        vv_out = vv_out.deviator() + volumetric_multiplier * volumetric_part;


        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi, TC_min_stress>;

private:

    static VoigtVector vv_out;
};

template<class NO_HARDENING>
VoigtVector MohrCoulombTensionCutoff_PF<NO_HARDENING>::vv_out;

#endif
