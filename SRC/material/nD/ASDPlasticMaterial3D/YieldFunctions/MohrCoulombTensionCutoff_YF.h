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

// Ladruno (ADR-84 P0): composite Mohr-Coulomb + tension-cutoff yield surface
// (fork addition). f = max(f_MC, f_TC) with f_MC verbatim MohrCoulomb_YF and
// f_TC = sigma_max - TC_min_stress (Rankine, tension positive). Implements
// SPECIAL_RETURN with exact principal-space returns on the cutoff face, the
// Rankine edge, the MC∩TC corner (Koiter, non-associated MC flow retained)
// and the apex point T_eff*I, T_eff = min(TC_min_stress, c*cot(phi)).
// Every trial state with the cutoff active resolves INSIDE the hook to a
// candidate validated against both yield values -- it is never handed to the
// generic scalar-Newton return map, whose exhaustion-accept path is what let
// f > 0 states through for plain Mohr-Coulomb near the apex.
//
// Ladruno (ADR-84 P3): the hook returns the RAW active-set (Koiter) consistent
// tangent. It used to secant-blend (E+D)/2 in here, which (a) silently
// overrode the material's configured `tangent_type` and (b) at the MC-cutoff
// corner fabricated ~2e6 kPa of stiffness in the direction whose true tangent
// is ~0, giving an 87%-wrong tangent -- measured as the cause of the Cerro
// Lindo M5 confined-shear Newton stall. Blending is now the integrator's
// decision, so `tangent_type Secant` (the default) is unchanged and
// `tangent_type Continuum` finally does what it says.
// See Ladruno_implementation/84_ladruno_mc_tension_cutoff_adr.md §9.

#ifndef MohrCoulombTensionCutoff_YF_H
#define MohrCoulombTensionCutoff_YF_H


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include <limits>
#include <Eigen/Eigenvalues>
#include "../AllASDModelParameterTypes.h"


template<class NO_HARDENING>
class MohrCoulombTensionCutoff_YF : public YieldFunctionBase<MohrCoulombTensionCutoff_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulombTensionCutoff_YF";


    MohrCoulombTensionCutoff_YF( ):
        YieldFunctionBase<MohrCoulombTensionCutoff_YF<NO_HARDENING>>::YieldFunctionBase()
        {}

    // f_MC: byte-for-byte the arithmetic of MohrCoulomb_YF::operator(), so that
    // stress paths on which the cutoff never competes are bit-identical to the
    // plain MohrCoulomb_YF material.
    static double yf_mc_value(const VoigtVector& sigma, double phi, double c)
    {
        using namespace std;

        double rThresh = c * cos(phi);
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        return rEquivalentStress - rThresh;
    }

    // f_TC: byte-for-byte the arithmetic of TensionCutoff_YF::operator().
    // principalStresses() returns ascending values, so s3 is the largest
    // (most tensile) principal stress.
    static double yf_tc_value(const VoigtVector& sigma, double min_stress)
    {
        double s1, s2, s3;
        std::tie(s1, s2, s3) = sigma.principalStresses();

        return s3 - min_stress;
    }

    YIELD_FUNCTION
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double min_stress = GET_PARAMETER_VALUE(TC_min_stress);

        double yf_mc = yf_mc_value(sigma, phi, c);
        double yf_tc = yf_tc_value(sigma, min_stress);

        return std::max(yf_mc, yf_tc);
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);
        double min_stress = GET_PARAMETER_VALUE(TC_min_stress);

        using namespace std;

        // Branch on the active surface at this stress state. Corner states
        // (both surfaces active) are resolved by special_return, never by the
        // scalar return map, so a hard branch here cannot chatter.
        double yf_mc = yf_mc_value(sigma, phi, c);
        double yf_tc = yf_tc_value(sigma, min_stress);

        if (yf_tc >= yf_mc)
        {
            // ---- Tension-cutoff (Rankine) branch: n = v3 (x) v3, the max
            //      principal direction, with factor-2 shear components in
            //      Voigt (same convention as the finite-difference gradients
            //      used elsewhere in this framework).
            Eigen::Matrix3d sig3;
            sig3 << sigma.v11(), sigma.v12(), sigma.v13(),
                    sigma.v12(), sigma.v22(), sigma.v23(),
                    sigma.v13(), sigma.v23(), sigma.v33();
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig3);

            const double s2 = es.eigenvalues()(1);
            const double s3 = es.eigenvalues()(2);

            bool degenerate = (s3 - s2) < 1e-8 * std::max(1.0, std::abs(s3));

            if (!degenerate)
            {
                Eigen::Vector3d v3 = es.eigenvectors().col(2);
                vv_out = direction_voigt(v3);

                for (int i = 0; i < 6; ++i)
                {
                    if (!std::isfinite(vv_out(i)))
                    {
                        degenerate = true;
                        break;
                    }
                }
            }

            if (degenerate)
            {
                // Coalescing principal stresses: the max-principal direction is
                // not unique. Fall back to central differences of f_TC (the
                // active branch), mirroring TensionCutoff_YF.
                double ds_tc = ds == 0 ? 1e-8 : ds;
                for (int i = 0; i < 6; ++i)
                {
                    VoigtVector SIG1 = sigma;
                    VoigtVector SIG2 = sigma;

                    SIG1(i) += ds_tc;
                    SIG2(i) -= ds_tc;

                    double yf1 = yf_tc_value(SIG1, min_stress);
                    double yf2 = yf_tc_value(SIG2, min_stress);

                    vv_out(i) = (yf1 - yf2) / (2 * ds_tc);
                }
            }

            return vv_out;
        }

        // ---- Mohr-Coulomb branch: verbatim MohrCoulomb_YF::df_dsigma_ij.
        //      (The YF(...) probes inside the numerical path evaluate the
        //      composite max, which equals f_MC wherever this branch is taken.)
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

    SPECIAL_RETURN
    {
        using namespace std;

        // Every rung below is a closed-form return except the terminal Stage-4
        // vertex projection, which downgrades this to SR_QUALITY_FALLBACK.
        return_quality = SR_QUALITY_EXACT;

        const double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        const double c   = GET_PARAMETER_VALUE(MC_c);
        const double psi = GET_PARAMETER_VALUE(MC_psi)*M_PI/180;
        const double T   = GET_PARAMETER_VALUE(TC_min_stress);

        // MC apex (hydrostatic tension) and effective cutoff. When T >= p_apex
        // the cutoff plane sits at/beyond the MC apex and cannot truncate it:
        // the composite degenerates to plain MC and the only special region
        // left is the MC apex cone itself.
        const double tan_phi = std::tan(phi);
        const bool  has_apex = tan_phi > 1e-12;
        const double p_apex  = has_apex ? c / tan_phi : std::numeric_limits<double>::infinity();
        const double T_eff   = std::min(T, p_apex);

        if (T >= p_apex)
        {
            static bool warned_T_above_apex = false;
            if (!warned_T_above_apex)
            {
                std::cout << "MohrCoulombTensionCutoff_YF: TC_min_stress = " << T
                          << " >= c/tan(phi) = " << p_apex
                          << " -- the cutoff sits at/beyond the MC apex and is clamped to it."
                          << " The material behaves as plain Mohr-Coulomb with an apex return." << std::endl;
                warned_T_above_apex = true;
            }
        }

        // Finite-trial guard: a NaN/Inf-contaminated trial must fail loudly in
        // the generic return map (which has an explicit NaN check), NOT be
        // silently projected onto the surface by the escalation chain below --
        // that would convert garbage into a plausible-looking admissible state.
        for (int i = 0; i < 6; ++i)
        {
            if (!std::isfinite(sigma_trial(i)))
                return false;
        }

        // Isotropy guard: the closed-form returns below assume an isotropic
        // state-frozen elastic tangent (the registered pairing uses
        // LinearIsotropic3D_EL). If a future pairing violates this, fall
        // through to the generic return map rather than returning wrong stress.
        const double lambda_L = Eelastic(0, 1);
        const double G        = (Eelastic(0, 0) - Eelastic(0, 1)) / 2.0;
        const double M_oed    = lambda_L + 2.0 * G;
        {
            const double scale = std::abs(Eelastic(0, 0)) + std::abs(G);
            if (G <= 0.0 ||
                std::abs(Eelastic(0, 2) - lambda_L) > 1e-8 * scale ||
                std::abs(Eelastic(3, 3) - G)        > 1e-8 * scale)
            {
                return false;
            }
        }

        // Spectral decomposition of the trial stress (ascending eigenvalues).
        Eigen::Matrix3d sig3;
        sig3 << sigma_trial.v11(), sigma_trial.v12(), sigma_trial.v13(),
                sigma_trial.v12(), sigma_trial.v22(), sigma_trial.v23(),
                sigma_trial.v13(), sigma_trial.v23(), sigma_trial.v33();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig3);
        const Eigen::Vector3d sp = es.eigenvalues();
        const Eigen::Matrix3d V  = es.eigenvectors();
        const double s1 = sp(0), s2 = sp(1), s3 = sp(2);

        if (s3 - T <= tol_yf)
        {
            // Cutoff inactive at the trial state. With T_eff < p_apex the
            // composite surface has no reachable MC apex, and the MC return
            // direction strictly decreases the max principal stress, so the
            // scalar-Newton return cannot cross the cutoff: fall through
            // (bit-identical to plain MohrCoulomb_YF).
            if (T >= p_apex && has_apex)
            {
                // Residual plain-MC apex cone (mirrors MohrCoulomb_YF::check_apex_region).
                const double p_trial = sigma_trial.getI1() / 3.0;
                if (p_trial > p_apex)
                {
                    return apex_projection(p_apex, sigma_trial, Eelastic,
                                           sigma_return, plastic_strain_incr, stiffness_return);
                }
            }
            return false;
        }

        // ------------------------------------------------------------------
        // Cutoff active (f_TC > tol at the trial). Resolve fully in-hook:
        // face -> edge -> corner -> frozen-direction MC -> apex. Every accepted
        // candidate is validated against both yield values; the chain ends at
        // the always-feasible apex point, so no f > 0 state can be accepted.
        // ------------------------------------------------------------------

        // ---- Stage 1: Rankine face return, s3 -> T (exact).
        bool face_needs_corner = false;
        {
            const double dl_face = (s3 - T) / M_oed;
            const double s1f = s1 - dl_face * lambda_L;
            const double s2f = s2 - dl_face * lambda_L;

            if (s2f <= T)
            {
                VoigtVector sigma_cand = assemble_from_principal(V, s1f, s2f, T);
                if (yf_mc_value(sigma_cand, phi, c) <= tol_yf)
                {
                    VoigtVector m3 = direction_voigt(V.col(2));
                    sigma_return        = sigma_cand;
                    plastic_strain_incr = dl_face * m3;
                    stiffness_return    = rank1_tangent(Eelastic, m3, m3);
                    return_quality      = SR_QUALITY_EXACT;
                    return true;
                }
                face_needs_corner = true; // MC also active at the face candidate
            }
        }

        // ---- Stage 2: Rankine edge return, s2 -> T and s3 -> T (exact 2x2).
        //      Only reached when the face candidate had s2f > T.
        if (!face_needs_corner)
        {
            const double det = M_oed * M_oed - lambda_L * lambda_L; // > 0 since G > 0
            const double la  = (M_oed * (s3 - T) - lambda_L * (s2 - T)) / det;
            const double lb  = (M_oed * (s2 - T) - lambda_L * (s3 - T)) / det;
            const double s1e = s1 - (la + lb) * lambda_L;

            if (la >= 0.0 && lb >= 0.0 && s1e <= T)
            {
                VoigtVector sigma_cand = assemble_from_principal(V, s1e, T, T);
                if (yf_mc_value(sigma_cand, phi, c) <= tol_yf)
                {
                    VoigtVector m3 = direction_voigt(V.col(2));
                    VoigtVector m2 = direction_voigt(V.col(1));
                    sigma_return        = sigma_cand;
                    plastic_strain_incr = la * m3 + lb * m2;
                    stiffness_return    = rank2_tangent(Eelastic, m3, m2, m3, m2);
                    return_quality      = SR_QUALITY_EXACT;
                    return true;
                }
                // Edge candidate violates MC: fall through to the corner.
            }
            else if (s1e > T)
            {
                // All three principals beyond the cutoff: composite vertex.
                return apex_projection(T_eff, sigma_trial, Eelastic,
                                       sigma_return, plastic_strain_incr, stiffness_return);
            }
        }

        // ---- Stage 3: MC∩TC(face) corner, Koiter two-multiplier return.
        {
            CornerResult res = corner_return(sigma_trial, Eelastic, tol_yf, V,
                                             phi, c, psi, T,
                                             sigma_return, plastic_strain_incr, stiffness_return);
            if (res == CornerResult::ACCEPTED)
                return true;

            if (res == CornerResult::MC_ONLY)
            {
                // The cutoff multiplier came out negative: the true return is
                // on the MC face alone. Resolve it here with a frozen-direction
                // 1D Newton rather than handing a cutoff-active trial to the
                // generic return map (whose branch switching could chatter and
                // whose loop accepts on exhaustion).
                if (frozen_mc_return(sigma_trial, Eelastic, tol_yf, V, phi, c, psi, T,
                                     sigma_return, plastic_strain_incr, stiffness_return))
                    return true;
            }
        }

        // ---- Stage 3b: compound corner MC ∩ TC(s3) ∩ TC(s2) — the
        //      triaxial-extension corner line of the composite. The attractor
        //      of symmetric lateral-tension drives (s2 = s3): the 2-surface
        //      corner leaves s2 unmanaged and would collapse these states to
        //      the vertex, losing the entire MC confinement.
        if (compound_corner_return(sigma_trial, Eelastic, tol_yf, V, phi, c, psi, T,
                                   sigma_return, plastic_strain_incr, stiffness_return))
            return true;

        // ---- Stage 3c: every exact cutoff-feature return rejected this
        //      trial. For an MC-DOMINANT trial (f_MC >= f_TC) that is the
        //      geometrically consistent outcome: its plastic increment is
        //      infeasible for every cutoff-feature cone because the true
        //      return lies on the MC surface BELOW the cutoff (e.g. the
        //      triaxial-extension ridge reached from a slightly-tensile
        //      lateral overshoot). Hand it to the scalar Newton, whose
        //      MC-branch gradient treats it exactly as plain MC would; if the
        //      Newton converges, max(f_MC, f_TC) <= tol holds by construction.
        if (yf_mc_value(sigma_trial, phi, c) >= yf_tc_value(sigma_trial, T))
            return false;

        // ---- Stage 4: apex of the composite surface (always feasible:
        //      f_TC(T_eff*I) <= 0 by clamp, f_MC(T_eff*I) <= 0 since T_eff <= p_apex).
        //      Reached only by TC-dominant trials none of whose exact returns
        //      validated. This is the ONE inexact rung: from a confined trial it
        //      replaces the whole deviatoric state with T_eff*I -- a large,
        //      unphysical stress drop that reads downstream as a modelling
        //      problem. Flagged so the integrator can refuse it loudly under
        //      strict_convergence instead of committing it silently.
        return_quality = SR_QUALITY_FALLBACK;
        return apex_projection(T_eff, sigma_trial, Eelastic,
                               sigma_return, plastic_strain_incr, stiffness_return);
    }


    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi, MC_c, MC_ds, MC_psi, TC_min_stress>;

private:

    enum class CornerResult { ACCEPTED, MC_ONLY, ESCALATE };

    // Assemble a Voigt stress from principal values and the trial eigenbasis.
    static VoigtVector assemble_from_principal(const Eigen::Matrix3d& V,
                                               double p1, double p2, double p3)
    {
        Eigen::Matrix3d S = p1 * V.col(0) * V.col(0).transpose()
                          + p2 * V.col(1) * V.col(1).transpose()
                          + p3 * V.col(2) * V.col(2).transpose();
        return VoigtVector(S(0, 0), S(1, 1), S(2, 2), S(0, 1), S(1, 2), S(0, 2));
    }

    // Rank-one direction v (x) v in Voigt with factor-2 shear components
    // (strain-like convention, consistent with E * m products and with
    // TrialPlastic_Strain updates in the integrator).
    static VoigtVector direction_voigt(const Eigen::Vector3d& v)
    {
        return VoigtVector(v(0) * v(0), v(1) * v(1), v(2) * v(2),
                           2.0 * v(0) * v(1), 2.0 * v(1) * v(2), 2.0 * v(0) * v(2));
    }

    // Assemble a principal-frame direction diag(d) into the trial eigenbasis,
    // Voigt with factor-2 shears.
    static VoigtVector principal_direction_voigt(const Eigen::Matrix3d& V,
                                                 const Eigen::Vector3d& d)
    {
        Eigen::Matrix3d Mten = d(0) * V.col(0) * V.col(0).transpose()
                             + d(1) * V.col(1) * V.col(1).transpose()
                             + d(2) * V.col(2) * V.col(2).transpose();
        return VoigtVector(Mten(0, 0), Mten(1, 1), Mten(2, 2),
                           2.0 * Mten(0, 1), 2.0 * Mten(1, 2), 2.0 * Mten(0, 2));
    }

    // Tensor contraction n : (E m). Both directions carry factor-2 shear
    // slots, E*m is single-count stress-like Voigt, so the plain dot product
    // is the correct double contraction.
    static double nEm_dot(const VoigtVector& n, const VoigtVector& Em)
    {
        return n.dot(Em);
    }

    // Consistent tangent for a single plane surface with frozen direction:
    // D = E - (E m)(n^T E) / (n : E m). For planes the consistent and
    // continuum tangents coincide.
    static VoigtMatrix rank1_tangent(const VoigtMatrix& E,
                                     const VoigtVector& n, const VoigtVector& m)
    {
        VoigtVector Em = E * m;
        VoigtVector En = E * n; // E symmetric: (n^T E)^T
        const double denom = nEm_dot(n, Em);
        VoigtMatrix D = E;
        if (denom > 0.0)
            D = E - (Em * En.transpose()) / denom;
        return D;
    }

    // Koiter two-surface consistent tangent with frozen directions:
    // D = E - [E m_a, E m_b] A^{-1} [E n_a, E n_b]^T,  A_ij = n_i : E m_j.
    static VoigtMatrix rank2_tangent(const VoigtMatrix& E,
                                     const VoigtVector& n_a, const VoigtVector& n_b,
                                     const VoigtVector& m_a, const VoigtVector& m_b)
    {
        VoigtVector Ema = E * m_a;
        VoigtVector Emb = E * m_b;
        VoigtVector Ena = E * n_a;
        VoigtVector Enb = E * n_b;

        Eigen::Matrix2d A;
        A << nEm_dot(n_a, Ema), nEm_dot(n_a, Emb),
             nEm_dot(n_b, Ema), nEm_dot(n_b, Emb);

        VoigtMatrix D = E;
        const double detA = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
        if (std::abs(detA) > 1e-20)
        {
            Eigen::Matrix2d Ainv;
            Ainv <<  A(1, 1) / detA, -A(0, 1) / detA,
                    -A(1, 0) / detA,  A(0, 0) / detA;

            Eigen::Matrix<double, 6, 2> EM;
            EM.col(0) = Ema; EM.col(1) = Emb;
            Eigen::Matrix<double, 6, 2> EN;
            EN.col(0) = Ena; EN.col(1) = Enb;

            D = E - EM * Ainv * EN.transpose();
        }
        return D;
    }

    // One-shot projection to the hydrostatic point p * I. Always feasible for
    // p <= min(TC_min_stress, c*cot(phi)). Terminal fallback of the escalation
    // chain: worst case a small conservative stress drop, never an accepted
    // f > 0 state.
    static bool apex_projection(double p, const VoigtVector& sigma_trial,
                                const VoigtMatrix& Eelastic,
                                VoigtVector& sigma_return,
                                VoigtVector& plastic_strain_incr,
                                VoigtMatrix& stiffness_return)
    {
        static int apex_event_count = 0;

        VoigtVector sigma_apex = p * kronecker_delta();
        Eigen::Matrix<double, 6, 1> rhs = sigma_trial - sigma_apex;
        Eigen::Matrix<double, 6, 1> dep = Eelastic.ldlt().solve(rhs);

        sigma_return        = sigma_apex;
        plastic_strain_incr = VoigtVector(dep);

        // Vertex tangent: rank-3 Koiter over the three (degenerate) Rankine
        // planes -- projects out ALL normal stiffness and keeps the shear
        // block. Full Eelastic here makes the global Newton limit-cycle when a
        // transient iterate overshoots into the vertex regime (measured: norm
        // stalls at ~2e-5 while plain MC recovers on the same path); the true
        // consistent tangent at a vertex is ~0.
        {
            VoigtVector m1(1, 0, 0, 0, 0, 0);
            VoigtVector m2(0, 1, 0, 0, 0, 0);
            VoigtVector m3(0, 0, 1, 0, 0, 0);
            VoigtVector Em1 = Eelastic * m1;
            VoigtVector Em2 = Eelastic * m2;
            VoigtVector Em3 = Eelastic * m3;

            Eigen::Matrix3d A;
            A << nEm_dot(m1, Em1), nEm_dot(m1, Em2), nEm_dot(m1, Em3),
                 nEm_dot(m2, Em1), nEm_dot(m2, Em2), nEm_dot(m2, Em3),
                 nEm_dot(m3, Em1), nEm_dot(m3, Em2), nEm_dot(m3, Em3);

            stiffness_return = Eelastic;
            if (std::abs(A.determinant()) > 1e-20)
            {
                Eigen::Matrix<double, 6, 3> EM;
                EM.col(0) = Em1; EM.col(1) = Em2; EM.col(2) = Em3;
                stiffness_return = VoigtMatrix(Eelastic - EM * A.inverse() * EM.transpose());
            }
        }

        if (apex_event_count < 5)
        {
            std::cout << "MohrCoulombTensionCutoff_YF: apex/vertex return to p = " << p
                      << " (event " << (apex_event_count + 1) << "; further events silent)" << std::endl;
        }
        ++apex_event_count;

        return true;
    }

    // MC∩TC(face) corner: Koiter two-multiplier return with directions frozen
    // at the trial state (planes in principal space; the residuals use the
    // framework's smoothed f_MC so the accepted point is admissible w.r.t.
    // the surface the material actually evaluates). Non-associated MC flow
    // (psi) is retained; the cutoff flow is associated Rankine.
    CornerResult corner_return(const VoigtVector& sigma_trial,
                               const VoigtMatrix& Eelastic,
                               double tol_yf,
                               const Eigen::Matrix3d& V,
                               double phi, double c, double psi, double T,
                               VoigtVector& sigma_return,
                               VoigtVector& plastic_strain_incr,
                               VoigtMatrix& stiffness_return) const
    {
        using namespace std;

        const double sin_phi = std::sin(phi);
        const double sin_psi = std::sin(psi);

        // Principal-frame plane normals/directions in the ordered sextant
        // (s1 <= s2 <= s3):
        //   MC yield normal   dF/ds = (-(1-sin phi)/2, 0, (1+sin phi)/2)
        //   MC flow direction dG/ds = (-(1-sin psi)/2, 0, (1+sin psi)/2)
        //   TC normal = flow  = (0, 0, 1)
        Eigen::Vector3d nmc_p(-(1.0 - sin_phi) / 2.0, 0.0, (1.0 + sin_phi) / 2.0);
        Eigen::Vector3d mmc_p(-(1.0 - sin_psi) / 2.0, 0.0, (1.0 + sin_psi) / 2.0);

        VoigtVector n_mc = principal_direction_voigt(V, nmc_p);
        VoigtVector m_mc = principal_direction_voigt(V, mmc_p);
        VoigtVector n_tc = direction_voigt(V.col(2));
        const VoigtVector& m_tc = n_tc; // associated on the cutoff

        VoigtVector Em_mc = Eelastic * m_mc;
        VoigtVector Em_tc = Eelastic * m_tc;

        Eigen::Matrix2d J;
        J << -nEm_dot(n_mc, Em_mc), -nEm_dot(n_mc, Em_tc),
             -nEm_dot(n_tc, Em_mc), -nEm_dot(n_tc, Em_tc);

        const double detJ = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
        if (std::abs(detJ) < 1e-20)
            return CornerResult::ESCALATE;

        Eigen::Matrix2d Jinv;
        Jinv <<  J(1, 1) / detJ, -J(0, 1) / detJ,
                -J(1, 0) / detJ,  J(0, 0) / detJ;

        double l_mc = 0.0, l_tc = 0.0;
        VoigtVector sigma_k = sigma_trial;

        const int max_corner_iter = 10;
        bool converged = false;
        for (int it = 0; it < max_corner_iter; ++it)
        {
            const double r0 = yf_mc_value(sigma_k, phi, c);
            const double r1 = yf_tc_value(sigma_k, T);

            if (std::abs(r0) < tol_yf && std::abs(r1) < tol_yf)
            {
                converged = true;
                break;
            }

            const double dl0 = -(Jinv(0, 0) * r0 + Jinv(0, 1) * r1);
            const double dl1 = -(Jinv(1, 0) * r0 + Jinv(1, 1) * r1);

            l_mc += dl0;
            l_tc += dl1;

            if (l_tc < 0.0)
                return CornerResult::MC_ONLY;
            if (l_mc < 0.0)
                return CornerResult::ESCALATE;

            sigma_k = sigma_trial - l_mc * Em_mc - l_tc * Em_tc;
        }

        if (!converged)
            return CornerResult::ESCALATE;

        // Reject if the middle principal stress crossed the cutoff (compound
        // corner not representable by these two surfaces): caller falls to apex.
        {
            double q1, q2, q3;
            std::tie(q1, q2, q3) = sigma_k.principalStresses();
            (void)q1; (void)q3;
            if (q2 > T + tol_yf)
                return CornerResult::ESCALATE;
        }

        sigma_return        = sigma_k;
        plastic_strain_incr = l_mc * m_mc + l_tc * m_tc;
        stiffness_return    = rank2_tangent(Eelastic, n_mc, n_tc, m_mc, m_tc);
        return CornerResult::ACCEPTED;
    }

    // Compound corner MC ∩ TC(s3) ∩ TC(s2): the target point is CLOSED FORM in
    // the trial principal frame — (s1*, T, T) with
    //   s1* = (T (1+sin phi) − 2 c cos phi) / (1 − sin phi),
    // the triaxial-extension corner of the composite (the framework's smoothed
    // f_MC collapses to the textbook form at theta = +30°, so s1* is exact).
    // Multipliers come from a cone-membership check: the coaxial plastic
    // strain increment must decompose non-negatively over {m_MC, e2⊗e2, e3⊗e3}.
    bool compound_corner_return(const VoigtVector& sigma_trial,
                                const VoigtMatrix& Eelastic,
                                double tol_yf,
                                const Eigen::Matrix3d& V,
                                double phi, double c, double psi, double T,
                                VoigtVector& sigma_return,
                                VoigtVector& plastic_strain_incr,
                                VoigtMatrix& stiffness_return) const
    {
        using namespace std;

        const double sin_phi = std::sin(phi);
        const double cos_phi = std::cos(phi);
        const double sin_psi = std::sin(psi);

        const double denom_phi = 1.0 - sin_phi;
        if (denom_phi < 1e-12)
            return false;
        const double s1_star = (T * (1.0 + sin_phi) - 2.0 * c * cos_phi) / denom_phi;

        VoigtVector sigma_cand = assemble_from_principal(V, s1_star, T, T);

        // Validate against the framework yield values.
        if (yf_mc_value(sigma_cand, phi, c) > tol_yf)
            return false;
        if (yf_tc_value(sigma_cand, T) > tol_yf)
            return false;

        // Coaxial plastic strain increment: dep = E^{-1} (sigma_trial - sigma_cand).
        Eigen::Matrix<double, 6, 1> rhs = sigma_trial - sigma_cand;
        Eigen::Matrix<double, 6, 1> dep = Eelastic.ldlt().solve(rhs);

        // Cone membership in the trial principal frame. At this corner BOTH
        // MC faces meeting at the triaxial-extension ridge are active
        // (s2 = s3 = T), so the flow cone is spanned by FOUR directions:
        //   m_MC(1,3) = (-a, 0, b),  m_MC(1,2) = (-a, b, 0),
        //   e2 = (0,1,0),  e3 = (0,0,1),   a = (1-sin psi)/2, b = (1+sin psi)/2.
        // dep decomposes non-negatively over these iff
        //   dep_1 <= 0,  dep_2 >= 0,  dep_3 >= 0,  and
        //   b * (-dep_1 / a) <= dep_2 + dep_3
        // (mu_total = -dep_1/a split over the two MC faces; the split exists
        // exactly under the last inequality). Testing only ONE MC face makes
        // symmetric lateral drives infeasible and collapses them to the
        // vertex, losing the entire MC confinement.
        {
            VoigtVector dep_v = VoigtVector(dep);
            Eigen::Vector3d dep_p;
            for (int i = 0; i < 3; ++i)
            {
                Eigen::Vector3d vi = V.col(i);
                VoigtVector ni = direction_voigt(vi);
                // strain-like Voigt dot with a unit rank-one direction gives
                // the normal component v_i . eps . v_i (shear slots carry the
                // factor 2 on exactly one side).
                dep_p(i) = dep_v.v11() * ni.v11() + dep_v.v22() * ni.v22() + dep_v.v33() * ni.v33()
                         + (dep_v.v12() * ni.v12() + dep_v.v23() * ni.v23() + dep_v.v13() * ni.v13()) / 2.0;
            }

            const double a = (1.0 - sin_psi) / 2.0;
            const double b = (1.0 + sin_psi) / 2.0;
            if (a <= 1e-12)
                return false;

            const double scale = std::abs(dep_p(0)) + std::abs(dep_p(1)) + std::abs(dep_p(2)) + 1e-300;
            const double slack = 1e-10 * scale;

            const double mu_total = -dep_p(0) / a;
            if (dep_p(0) > slack)               return false; // needs axial extension: not this corner
            if (dep_p(1) < -slack)              return false;
            if (dep_p(2) < -slack)              return false;
            if (b * mu_total > dep_p(1) + dep_p(2) + slack) return false;
        }

        // Rank-3 Koiter consistent tangent over {MC, TC(s3), TC(s2)} with
        // frozen directions (phi normals, psi flow on MC; associated Rankine).
        {
            Eigen::Vector3d nmc_p(-(1.0 - sin_phi) / 2.0, 0.0, (1.0 + sin_phi) / 2.0);
            Eigen::Vector3d mmc_p(-(1.0 - sin_psi) / 2.0, 0.0, (1.0 + sin_psi) / 2.0);
            VoigtVector n1 = principal_direction_voigt(V, nmc_p);
            VoigtVector m1 = principal_direction_voigt(V, mmc_p);
            VoigtVector n2 = direction_voigt(V.col(1));
            VoigtVector n3 = direction_voigt(V.col(2));

            VoigtVector Em1 = Eelastic * m1;
            VoigtVector Em2 = Eelastic * n2; // associated on the cutoffs
            VoigtVector Em3 = Eelastic * n3;
            VoigtVector En1 = Eelastic * n1;
            VoigtVector En2 = Em2;
            VoigtVector En3 = Em3;

            Eigen::Matrix3d A;
            A << nEm_dot(n1, Em1), nEm_dot(n1, Em2), nEm_dot(n1, Em3),
                 nEm_dot(n2, Em1), nEm_dot(n2, Em2), nEm_dot(n2, Em3),
                 nEm_dot(n3, Em1), nEm_dot(n3, Em2), nEm_dot(n3, Em3);

            stiffness_return = Eelastic;
            if (std::abs(A.determinant()) > 1e-20)
            {
                Eigen::Matrix<double, 6, 3> EM;
                EM.col(0) = Em1; EM.col(1) = Em2; EM.col(2) = Em3;
                Eigen::Matrix<double, 6, 3> EN;
                EN.col(0) = En1; EN.col(1) = En2; EN.col(2) = En3;
                stiffness_return = VoigtMatrix(Eelastic - EM * A.inverse() * EN.transpose());
            }
        }

        sigma_return        = sigma_cand;
        plastic_strain_incr = VoigtVector(dep);
        return true;
    }

    // Frozen-direction single-surface MC return: 1D Newton on the framework's
    // smoothed f_MC along the fixed direction E*m_MC from the trial principal
    // frame. Used when the corner solve reports the cutoff inactive, so that a
    // cutoff-active trial is never handed to the generic scalar return map.
    bool frozen_mc_return(const VoigtVector& sigma_trial,
                          const VoigtMatrix& Eelastic,
                          double tol_yf,
                          const Eigen::Matrix3d& V,
                          double phi, double c, double psi, double T,
                          VoigtVector& sigma_return,
                          VoigtVector& plastic_strain_incr,
                          VoigtMatrix& stiffness_return) const
    {
        using namespace std;

        const double sin_phi = std::sin(phi);
        const double sin_psi = std::sin(psi);

        Eigen::Vector3d nmc_p(-(1.0 - sin_phi) / 2.0, 0.0, (1.0 + sin_phi) / 2.0);
        Eigen::Vector3d mmc_p(-(1.0 - sin_psi) / 2.0, 0.0, (1.0 + sin_psi) / 2.0);

        VoigtVector n_mc = principal_direction_voigt(V, nmc_p);
        VoigtVector m_mc = principal_direction_voigt(V, mmc_p);
        VoigtVector Em_mc = Eelastic * m_mc;

        const double dPhi_dl = -nEm_dot(n_mc, Em_mc);
        if (std::abs(dPhi_dl) < 1e-20)
            return false;

        double l_mc = 0.0;
        VoigtVector sigma_k = sigma_trial;

        const int max_iter = 15;
        bool converged = false;
        for (int it = 0; it < max_iter; ++it)
        {
            const double r = yf_mc_value(sigma_k, phi, c);
            if (std::abs(r) < tol_yf)
            {
                converged = true;
                break;
            }
            l_mc += -r / dPhi_dl;
            if (l_mc < 0.0)
                return false;
            sigma_k = sigma_trial - l_mc * Em_mc;
        }

        if (!converged)
            return false;

        // The corner solve claimed the cutoff is inactive here -- verify.
        if (yf_tc_value(sigma_k, T) > tol_yf)
            return false;

        sigma_return        = sigma_k;
        plastic_strain_incr = l_mc * m_mc;
        stiffness_return    = rank1_tangent(Eelastic, n_mc, m_mc);
        return true;
    }

    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class NO_HARDENING>
VoigtVector MohrCoulombTensionCutoff_YF<NO_HARDENING>::vv_out;

// Opts in to the exact special-region return in the constitutive integrator.
// Deliberately does NOT declare yf_has_apex: the legacy apex hook takes no
// trial stress and cannot express the cutoff-corner targets.
template<class NO_HARDENING>
struct yf_has_special_return<MohrCoulombTensionCutoff_YF<NO_HARDENING>> : std::true_type {};

#endif
