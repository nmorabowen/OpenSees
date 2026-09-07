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
    using parameters_t = tuple<ScalarLinearHardeningParameter>;
};

struct NullHardeningScalarPolicy {
    static constexpr const char* NAME = "NullHardeningScalarFunction";
    HARDENING_FUNCTION_DEFINITION 
    {
        double zero=0;
        return zero;
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
    using parameters_t = tuple<AF_ha, AF_cr>;
};



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
