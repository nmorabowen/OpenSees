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


#ifndef _AllASDHardeningFunctions
#define _AllASDHardeningFunctions

#include "ASDPlasticMaterial3DGlobals.h" 
#include "AllASDModelParameterTypes.h" 
#include "HardeningFunction.h" 


// Hardening policies
struct LinearHardeningForTensorPolicy {
    static constexpr const char* NAME = "TensorLinearHardeningFunction";
    HARDENING_FUNCTION_DEFINITION
    {
        double H = GET_PARAMETER_VALUE(TensorLinearHardeningParameter);
        VoigtVector h = H * m.deviator();  // best not to use 'auto' here
        return h;
    }

    // Ladruno (ADR-97 wp/97b): h = H*dev(m) is independent of q and LINEAR in m,
    // so dh/dq = 0 and dh/dm = H * I_dev (the deviatoric projector on stored
    // Voigt slots: d(dev v)_i/d v_j = delta_ij - 1/3 for i,j < 3, delta_ij else).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        double H = GET_PARAMETER_VALUE(TensorLinearHardeningParameter);
        out.setZero();
        for (int i = 0; i < 6; ++i) out(i, i) = H;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= H / 3.0;
    }

    using parameters_t = tuple<TensorLinearHardeningParameter>;
};

struct LinearHardeningForScalarPolicy {
    static constexpr const char* NAME = "ScalarLinearHardeningFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);
        // Ladruno (ADR-94 wp/94c, B5): the equivalent plastic strain rate is
        // sqrt(2/3 * m_ij m_ij), and m is an ENGINEERING-shear Voigt vector, so the
        // shear slots enter with weight 1/2, not 1.  A plain `m.dot(m)` mis-weighted
        // them for every yield function.  Identical on any shear-free path; for the
        // (now Voigt) von Mises flow direction it restores h == H*sqrt(2/3) exactly,
        // on every stress path rather than only on the pure-normal axis.
        double h = H * sqrt((2 * tensor_dot_engineering_strain_like(m, m)) / 3);
        return h;
    }

    // Ladruno (ADR-97 wp/97b): h = H*sqrt(2/3 <m,m>_e) with <.,.>_e the
    // ENGINEERING-strain contraction (shear weight 1/2).  Independent of q, so
    //     dh/dq   = 0
    //     dh/dm_j = H * (2/3) * (W_e m)_j / sqrt(2/3 <m,m>_e)
    // (one row).  At <m,m>_e = 0 the rate itself is zero and non-differentiable;
    // the derivative is taken as zero there, which is the correct one-sided limit
    // for a flow direction that has not yet been established.
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) sigma;
        double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);
        out.setZero();
        const double mm = tensor_dot_engineering_strain_like(m, m);
        const double eq = sqrt((2.0 / 3.0) * mm);
        if (eq > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
        {
            const double c = H * (2.0 / 3.0) / eq;
            out(0, 0) = c * m(0);
            out(0, 1) = c * m(1);
            out(0, 2) = c * m(2);
            out(0, 3) = c * 0.5 * m(3);
            out(0, 4) = c * 0.5 * m(4);
            out(0, 5) = c * 0.5 * m(5);
        }
    }

    using parameters_t = tuple<ScalarLinearHardeningParameter>;
};

struct NullHardeningScalarPolicy {
    static constexpr const char* NAME = "NullHardeningScalarFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        double zero=0;
        return zero;
    }

    // Ladruno (ADR-97 wp/97b): h == 0, so both derivatives are zero (scalar).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    using parameters_t = tuple<>;
};

struct NullHardeningTensorPolicy {
    static constexpr const char* NAME = "NullHardeningTensorFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        // Ladruno (ADR-94 wp/94a): ADR-94 B4 -- this RETURNED uninitialised Eigen
        // storage. Every consumer of a Null tensor hardening law read whatever the
        // heap happened to hold, NaN included.
        VoigtVector zero = VoigtVector::Zero();
        return zero;
    }

    // Ladruno (ADR-97 wp/97b): h == 0, so both derivatives are zero (tensor).
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) current_value; (void) depsilon; (void) m; (void) sigma;
        (void) parameters_storage;
        out.setZero();
    }

    using parameters_t = tuple<>;
};

// struct ExponentialHardeningForScalarPolicy {
//     static constexpr const char* NAME = "ScalarExponentialLinear";
//     HARDENING_FUNCTION_DEFINITION 
//     {
//         double Sigma0 = GET_PARAMETER_VALUE(ScalarExponentialLinear_Sigma0);
//         double SigmaInf = GET_PARAMETER_VALUE(ScalarExponentialLinear_SigmaInf);
//         double delta = GET_PARAMETER_VALUE(ScalarExponentialLinear_delta);
//         double H = GET_PARAMETER_VALUE(ScalarLinearHardeningParameter);

//         // This is the same as the hardening of the "q(xi)" in terms of the internal variable "xi" in J2Plasticity,
//         // xi is root23 accumulated plastic which we here compute
//         double dxi = QRT_2_over_3 * dLambda;
        
//         // Because of how ASDPlasticMaterial3D works, instead of using the explicit form of 
//         // q(xi) used in J2 plasticity, we here use the incremental form. Thus, this
//         // function returns dq / dxi 
//         // double dq_over_dxi = -delta * (SigmaInf - Sigma0) * exp(-delta*)
        
//     }
//     using parameters_t = tuple<ScalarExponentialLinear_Sigma0,ScalarExponentialLinear_SigmaInf, ScalarExponentialLinear_delta, ScalarLinearHardeningParameter>;
// };

struct ArmstrongFrederickPolicy {
    static constexpr const char* NAME = "ArmstrongFrederickHardeningFunction";
    HARDENING_FUNCTION_DEFINITION
    {
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);

        auto alpha = current_value;
        auto alpha_dev = alpha.deviator();


        // Ladruno (ADR-94 wp/94c, B5): `sqrt(2/3 * v.squaredNorm())` is a PLAIN sum
        // of squares, which mis-weights the three shear slots for BOTH kinds of
        // operand this lambda was used on.  Strain-like Voigt vectors (the flow
        // direction m, depsilon) store engineering shear (gamma = 2 eps), so
        // eps_ij eps_ij = normals + 0.5*shears; the back stress alpha is stress-like
        // and stores tensor shear, so a_ij a_ij = normals + 2*shears (red2 Q3
        // reported the alpha side as "under-weighting shear by 2").  Each operand
        // now uses the matching contraction.  Identical on any shear-free path.
        auto eq_norm_strain = [](const VoigtVector& v){ return sqrt((2. / 3.) * tensor_dot_engineering_strain_like(v, v)); };
        auto eq_norm_stress = [](const VoigtVector& v){ return sqrt((2. / 3.) * tensor_dot_stress_like(v, v)); };

        VoigtVector mdev = m.deviator();
        double mdev_eq = eq_norm_strain(mdev);
        VoigtVector dEPS_dev = depsilon.deviator();
        double dEPS_dev_eq = eq_norm_strain(dEPS_dev);
        VoigtVector alpha_dev_v = alpha_dev;
        double alpha_norm = eq_norm_stress(alpha_dev_v);
        // double alpha_norm = alpha_dev.norm();
        // double alpha_limit = sqrt(2. / 3.) * ha / cr;
        double alpha_limit =  ha / cr;

        cout << "depsilon    = " << depsilon.transpose() << endl;
        cout << "m           = " << m.transpose() << endl;
        cout << "mdev        = " << mdev.transpose() << endl;
        cout << "alpha       = " << alpha.transpose() << endl;
        cout << "alpha_norm  = " << alpha_norm << "  <= alpha_limit = " << alpha_limit <<  endl;
        VoigtVector derivative = VoigtVector::Zero();   // Ladruno (ADR-94 wp/94a): ADR-94 B4

        //Compute the derivative (hardening function)
        if (alpha_norm >= alpha_limit)
        {
            cout << "Saturation!" << endl;
            derivative.setZero();  // Ladruno (ADR-94 wp/94a): was `*= 0` on uninitialised storage
        }
        else
        {
            derivative =   ha * mdev - cr * mdev_eq * alpha_dev;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha_dev;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha;
            // derivative =  (2. / 3.) * ha * mdev - cr * sqrt((2. / 3.) * mdev.dot(mdev)) * alpha;
            // derivative =  ha * dEPS_dev - cr * dEPS_dev_eq * alpha_dev;
        }
        cout << "----> derivative = " << derivative.transpose() << endl;

        return derivative;

    }

    // Ladruno (ADR-97 wp/97b): Armstrong-Frederick, differentiated EXACTLY as
    // the branch `f` above actually took -- including the saturation branch.
    //   unsaturated:  h = ha*dev(m) - cr*|dev m|_eq * dev(alpha),
    //                 |v|_eq = sqrt(2/3 <v,v>_e)
    //     dh/dalpha = -cr * |dev m|_eq * I_dev
    //     dh/dm     =  ha * I_dev
    //                 - cr * dev(alpha) (x) [ (2/3) * (W_e dev(m)) / |dev m|_eq ]
    //   saturated (|dev alpha|_eq >= ha/cr):  h == 0, so both blocks are zero.
    //
    // ADR-97 plan-note: the plan suggested Closest_Point DROP the saturation
    // branch (it is non-differentiable, and the implicit AF update is already a
    // contraction toward ||alpha|| = ha/cr so it cannot overshoot).  It is kept
    // here instead, because `f` is SHARED with Backward_Euler (D1: byte
    // identical) and a CP-only `f` would be a second hardening law with the same
    // name.  Differentiating the branch actually taken keeps CP's Jacobian
    // exactly consistent with CP's own residual, which is what both the Newton
    // and the consistent tangent require; the P0 oracle (cppm_vm.py) mirrors the
    // same branch, so the pinned reference values are unaffected.
    HARDENING_FUNCTION_IV_DERIVATIVE
    {
        (void) depsilon; (void) sigma;
        out.setZero();
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);
        if (!(cr > 0.0))
            return;                       // no recovery term -> h is q-independent
        VoigtVector alpha_dev = current_value.deviator();
        const double alpha_norm =
            sqrt((2. / 3.) * tensor_dot_stress_like(alpha_dev, alpha_dev));
        if (alpha_norm >= ha / cr)
            return;                       // saturated branch: h == 0
        VoigtVector mdev = m.deviator();
        const double mdev_eq =
            sqrt((2. / 3.) * tensor_dot_engineering_strain_like(mdev, mdev));
        const double c = -cr * mdev_eq;
        for (int i = 0; i < 6; ++i) out(i, i) = c;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= c / 3.0;
    }

    HARDENING_FUNCTION_M_DERIVATIVE
    {
        (void) depsilon; (void) sigma;
        out.setZero();
        double ha = GET_PARAMETER_VALUE(AF_ha);
        double cr = GET_PARAMETER_VALUE(AF_cr);
        VoigtVector alpha_dev = current_value.deviator();
        if (cr > 0.0)
        {
            const double alpha_norm =
                sqrt((2. / 3.) * tensor_dot_stress_like(alpha_dev, alpha_dev));
            if (alpha_norm >= ha / cr)
                return;                   // saturated branch: h == 0
        }
        // ha * I_dev
        for (int i = 0; i < 6; ++i) out(i, i) = ha;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out(i, j) -= ha / 3.0;
        if (!(cr > 0.0))
            return;
        VoigtVector mdev = m.deviator();
        const double mdev_eq =
            sqrt((2. / 3.) * tensor_dot_engineering_strain_like(mdev, mdev));
        if (!(mdev_eq > 100 * ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return;
        // row vector d|dev m|_eq / dm  =  (2/3) * W_e dev(m) / |dev m|_eq
        double drow[6];
        const double s = (2.0 / 3.0) / mdev_eq;
        drow[0] = s * mdev(0);
        drow[1] = s * mdev(1);
        drow[2] = s * mdev(2);
        drow[3] = s * 0.5 * mdev(3);
        drow[4] = s * 0.5 * mdev(4);
        drow[5] = s * 0.5 * mdev(5);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                out(i, j) -= cr * alpha_dev(i) * drow[j];
    }

    using parameters_t = tuple<AF_ha, AF_cr>;
};

// Ladruno (ADR-97 wp/97b): the five hardening policies P1 solves implicitly at
// n+1.  Every other policy in the tree (the StiffSoil pair) keeps the inert zero
// default and is REFUSED at parse time under `integration_method Closest_Point`.
template <> struct hardening_policy_has_cp_derivatives<LinearHardeningForTensorPolicy> : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<LinearHardeningForScalarPolicy> : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<NullHardeningScalarPolicy>      : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<NullHardeningTensorPolicy>      : std::true_type {};
template <> struct hardening_policy_has_cp_derivatives<ArmstrongFrederickPolicy>       : std::true_type {};



// Ladruno (HB/StiffSoil integration, ledger row 337): StiffSoil shear/cap hardening IVs
#include "StiffSoil_HardeningFunctions.h"

// Plastic deviatoric strain for shear mechanism
struct EpsQpShearName { static constexpr const char* name = "EpsQpShear"; };
using EpsQpShear = InternalVariableType<VoigtScalar, StiffSoilShearHardeningFunction, EpsQpShearName>;

// Cap pressure for volumetric mechanism  
struct CapPressureName { static constexpr const char* name = "CapPressure"; };
using CapPressure = InternalVariableType<VoigtScalar, StiffSoilCapHardeningFunction, CapPressureName>;

// Cap pressure with linear hardening (alternative)
using CapPressureLinear = InternalVariableType<VoigtScalar, StiffSoilCapLinearHardeningFunction, CapPressureName>;




// Aliases for HardeningFunction with specific hardening policies
using TensorLinearHardeningFunction = HardeningFunction<VoigtVector, LinearHardeningForTensorPolicy>;
using ScalarLinearHardeningFunction = HardeningFunction<VoigtScalar, LinearHardeningForScalarPolicy>;
// using ScalarExponentialHardeningFunction = HardeningFunction<VoigtScalar, ExponentialHardeningForScalarPolicy>;
using ArmstrongFrederickHardeningFunction = HardeningFunction<VoigtVector, ArmstrongFrederickPolicy>;
using NullHardeningScalarFunction = HardeningFunction<VoigtScalar, NullHardeningScalarPolicy>;
using NullHardeningTensorFunction = HardeningFunction<VoigtVector, NullHardeningTensorPolicy>;


#endif //not defined _AllASDHardeningFunctions
