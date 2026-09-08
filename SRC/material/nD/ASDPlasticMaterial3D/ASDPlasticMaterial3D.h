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



#ifndef ASDPlasticMaterial3D_H
#define ASDPlasticMaterial3D_H


#include "NDMaterial.h"
#include <G3Globals.h>
#include <iostream>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <LadrunoMaterialStatus.h>  // Ladruno (ADR-84 -> ADR-86b): LADRUNO_MATERIAL_REFUSED

#include "ASDPlasticMaterial3DTraits.h"

#include "ASDPlasticMaterial3DGlobals.h"


// #include "MaterialInternalVariables.h"
#include "YieldFunctions/AllYieldFunctions.h"
#include "PlasticFlowDirections/AllPlasticFlowDirections.h"
#include "ElasticityModels/AllElasticityModels.h"
#include "AllASDModelParameterTypes.h"
#include "AllASDInternalVariableTypes.h"
#include "AllASDHardeningFunctions.h"
#include <Eigen/Eigenvalues>   // Ladruno (ADR-97 wp/97c): SelfAdjointEigenSolver

#include "utuple_storage.h"

// for print p, q, theta
#include <tuple>
#include <utility> // For std::pair
#include <map> // For std::pair
#include <cstdio> // Ladruno: snprintf for internal-variable ResponseType names
#include <limits>
#include <type_traits>

#include "std_tuple_concat.h"

// for debugging printing
#include <fstream>

#define ASDPlasticMaterial3D_MAXITER_BRENT 50


using namespace ASDPlasticMaterial3DGlobals;

#define ASDP_TAG this->getTag()

// Ladruno (ADR-97 wp/97b): does EVERY internal variable of this specialization
// carry a hardening law with analytic closest-point derivatives (dh/dq, dh/dm)?
// Folded over the concatenated IV tuple at compile time so the parser can refuse
// `integration_method Closest_Point` for a specialization this ADR has not
// converted, instead of silently believing the inert zero defaults.
template <typename Tuple>
struct asdp_all_ivs_support_cp;

template <typename... Ts>
struct asdp_all_ivs_support_cp<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_supports_cp() && ... && true);
};

// Ladruno (ADR-97 wp/97c): is EVERY internal variable of this specialization
// FIXED (h == 0)?  The principal-stress-space Mohr-Coulomb return of P2 is a
// closed-form projection onto a surface it assumes does not move: it carries no
// q-row, so a live hardening law would silently be ignored.  Folding this at
// compile time is what keeps a hypothetical MohrCoulomb_YF<ArmstrongFrederick..>
// out of that path and REFUSED at parse time instead.
template <typename Tuple>
struct asdp_all_ivs_are_inert;

template <typename... Ts>
struct asdp_all_ivs_are_inert<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_is_perfectly_plastic() && ... && true);
};

template <
    class ElasticityType,
    class YieldFunctionType,
    class PlasticFlowType,
    int thisClassTag >
class ASDPlasticMaterial3D : public NDMaterial
{

public:

    // Concatenate the internal varibles into the storage
    using iv_concat_types = utuple_concat_unique_type <
                            typename YieldFunctionType::internal_variables_t,
                            typename PlasticFlowType::internal_variables_t >;
    using iv_storage_t = utuple_storage<iv_concat_types>;

    // Concatenate the model parameters into the parameters storage
    using extracted_parameters_t = utuple_concat_unique_type <
                                   ExtractNestedParameterTypes_t<typename YieldFunctionType::internal_variables_t>,
                                   ExtractNestedParameterTypes_t<typename PlasticFlowType::internal_variables_t>
                                   >;

    using parameters_concat_types = utuple_concat_type <
                                    typename YieldFunctionType::parameters_t,
                                    typename PlasticFlowType::parameters_t,
                                    typename ElasticityType::parameters_t,
                                    extracted_parameters_t,
                                    std::tuple<MassDensity>,
                                    std::tuple<InitialP0>
                                    >;
    using parameters_storage_t = utuple_storage<parameters_concat_types>;

    // Ladruno (ADR-97 wp/97b): is `integration_method Closest_Point` available for
    // THIS specialization?  All three legs must opt in: the yield function must
    // supply the uncontracted df/dq, the plastic flow direction the analytic
    // dm/dsigma and dm/dq, and every internal variable's hardening law dh/dq and
    // dh/dm.  P1 ships VonMises + DruckerPrager x {Null, Linear scalar/tensor,
    // ArmstrongFrederick} = 20 of the 46 registered specializations; the rest are
    // refused at parse time naming the ADR phase that will deliver them.
    static constexpr bool ladruno_cp_smooth_supported =
        yf_has_cp_derivatives<YieldFunctionType>::value &&
        pf_has_cp_derivatives<PlasticFlowType>::value &&
        asdp_all_ivs_support_cp<iv_concat_types>::value;

    // Ladruno (ADR-97 wp/97c): the PRINCIPAL-STRESS-SPACE family of this
    // specialization (0 = none).  Non-zero only when the yield function and the
    // plastic flow direction carry the SAME non-zero marker -- so the generator's
    // cross pairings (MohrCoulomb_YF x VonMises_PF, VonMises_YF x MohrCoulomb_PF,
    // HoekBrown_YF x MohrCoulomb_PF, ...) are NOT enabled by P2: the principal
    // return assumes both the surface and the potential are the piecewise-linear
    // Mohr-Coulomb ones, and no oracle covers the mixed maps.  They stay refused.
    // Perfect plasticity is required as well (the map has no q-row).
    static constexpr int ladruno_cp_principal_family =
        (yf_cp_principal_family<YieldFunctionType>::value != 0
         && yf_cp_principal_family<YieldFunctionType>::value
            == pf_cp_principal_family<PlasticFlowType>::value
         && asdp_all_ivs_are_inert<iv_concat_types>::value)
        ? yf_cp_principal_family<YieldFunctionType>::value : 0;

    static constexpr bool ladruno_cp_supported =
        ladruno_cp_smooth_supported || (ladruno_cp_principal_family != 0);

    static constexpr bool supportsClosestPoint() { return ladruno_cp_supported; }

    // Ladruno (ADR-97 wp/97b): Newton system size cap.  6 stress rows + the
    // internal variables + one consistency row; the widest registered
    // specialization has 14 IV components, so 26 leaves headroom and keeps the
    // Jacobian on the stack (Eigen fixed-max dynamic storage, no heap traffic in
    // the Gauss-point loop).  A specialization that exceeds it is REFUSED, loudly.
    static constexpr int ASDP_CP_MAXN = 26;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::ColMajor, ASDP_CP_MAXN, ASDP_CP_MAXN> cp_matrix_t;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1,
                          Eigen::ColMajor, ASDP_CP_MAXN, 1> cp_vector_t;


    //==================================================================================================
    //  Constructors
    //==================================================================================================

    ASDPlasticMaterial3D( )
        : NDMaterial(0, thisClassTag)
    {
        stress_set_externally = false; // Ladruno (HB/StiffSoil integration, ledger row 337): init new flag

        // Ladruno (ADR-94 wp/94b, M1/F1): dsigma / depsilon_elpl / intersection_* /
        // Stiffness used to be class-statics (one copy per <E,Y,P,tag> specialization,
        // zero-initialised once by the loader and then shared by every Gauss point,
        // element and material tag). They are ordinary members now, so each instance
        // must zero its own -- Eigen fixed-size storage is NOT zero-initialised.
        // first_step was never initialised in this constructor at all.
        first_step = true;
        TrialStress.setZero();
        CommitStress.setZero();
        TrialStrain.setZero();
        CommitStrain.setZero();
        TrialPlastic_Strain.setZero();
        CommitPlastic_Strain.setZero();
        dsigma.setZero();
        depsilon_elpl.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();
        Stiffness.setZero();
        initial_iv_captured = false;
    }


    ASDPlasticMaterial3D(int tag)
        : NDMaterial(tag, thisClassTag)
    {


        // Ladruno (ADR-84 P0): these members are uninitialized Eigen storage,
        // and "*= 0" does NOT clear heap garbage that decodes as NaN/Inf
        // (NaN*0 == NaN). On a churned heap the shear slots kept NaN and every
        // yf comparison silently went false. setZero() is unconditional.
        TrialStress.setZero();
        CommitStress.setZero();
        TrialStrain.setZero();
        CommitStrain.setZero();
        TrialPlastic_Strain.setZero();
        CommitPlastic_Strain.setZero();

        first_step = true;
        stress_set_externally = false;

        // Ladruno (ADR-94 wp/94b, M1/F1): per-instance now (were class-statics).
        dsigma.setZero();
        depsilon_elpl.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();
        Stiffness.setZero();
        initial_iv_captured = false;
    }


    ~ASDPlasticMaterial3D(void)
    {

    }

    //==================================================================================================
    // To set internal variables values for the model
    //==================================================================================================
    auto getInternalVariablesNames() const
    {
        return iv_storage.getParameterNames();
    }

    int getInternalVariableSizeByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableSizeByName(iv_name);
    }

    int getInternalVariableIndexByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableIndexByName(iv_name);
    }

    auto setInternalVariableByName(const char * iv_name, int iv_size, double* iv_values)
    {
        cout << "  --->  Setting " << iv_name << " = ";
        for (int i = 0; i < iv_size; ++i)
        {
            cout << iv_values[i] << " ";
        }
        cout << endl;
        return iv_storage.setInternalVariableByName(iv_name, iv_size, iv_values);
    }



    //==================================================================================================
    // To set parameter values for the model
    //==================================================================================================
    auto getParameterNames() const
    {
        return parameters_storage.getParameterNames();
    }

    auto setParameterByName(const char * param_name, double param_value)
    {
        cout << "  --->  Setting " << param_name << " = " << param_value << endl;
        return parameters_storage.setParameterByName(param_name, param_value);
    }


    //==================================================================================================
    //  Class type function
    //==================================================================================================
    const char *getClassType(void) const
    {
        std::string name("ASDPlasticMaterial3D");

        return name.c_str();
    };

    double getRho(void)
    {
        return parameters_storage.template get<MassDensity>().value;
    }

    double getPressure(void)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        return -CommitStress.trace() / 3;
    }

    std::string getYFName() const {return YieldFunctionType::NAME;}
    std::string getPFName() const {return PlasticFlowType::NAME;}
    std::string getELName() const {return ElasticityType::NAME;}
    std::string getIVName() const {return iv_storage.getVariableNamesAndHardeningLaws();}


    //==================================================================================================
    //  Set Trial strain and trial strain increment
    //==================================================================================================
    // For total strain-based elements.
    // Receives the current total strain at a GP.
    // This function then computes the incremental strain (subtracting from the committed one)
    // and sets the increment.
    // Returns a success flag from the call to setTrialStrainIncr

    int setTrialStrain(const Vector &v)
    {
        // Ladruno (ADR-94 wp/94b, M6): snapshot the internal variables exactly once, on
        // the first strain any element hands us. At that point the parser has finished
        // configuring this instance and no integration has run, so this IS the "start"
        // state revertToStart() has to restore. utuple_storage offers commit_all() and
        // revert_all() only -- it has no initial-value facility of its own.
        if (!initial_iv_captured)
        {
            iv_storage_initial = iv_storage;
            initial_iv_captured = true;
        }

        // Ladruno (HB/StiffSoil integration, ledger row 337): skip K0 init if stress was set externally
        if (first_step && !stress_set_externally)
        {
            double p0 = parameters_storage.template get<InitialP0>().value;
            TrialStress(0) = p0;
            TrialStress(1) = p0;
            TrialStress(2) = p0;
            CommitStress(0) = p0;
            CommitStress(1) = p0;
            CommitStress(2) = p0;
        }


        TrialStrain = VoigtVector::fromStrain(v);
        return setTrialStrainIncr( TrialStrain - CommitStrain );
    }

    int setTrialStrainIncr( const VoigtVector &strain_increment )
    {

        int exitflag = -1;

        switch (INT_OPT_constitutive_integration_method[ASDP_TAG])
        {
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Not_Set :
            exitflag = -1;
            cerr << "CEP::setTrialStrainIncr - Integration method not set!\n" ;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler :
            exitflag = this->Forward_Euler(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement :
            exitflag = this->Forward_Euler_Subincrement(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler :
            exitflag = this->Backward_Euler(strain_increment);;
            break;
        // Ladruno (ADR-97 wp/97b): the fully implicit closest-point return map.
        // Backward_Euler above is untouched (ADR-97 D1: byte-identical).
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point :
            exitflag = this->Closest_Point(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch :
            exitflag = this->Backward_Euler_LineSearch(strain_increment);;
            break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement :
        //     exitflag = this->Backward_Euler_ddlambda_Subincrement(strain_increment);;
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Crisfield :
        //     exitflag = this->Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield :
        //     exitflag = this->Multistep_Forward_Euler(strain_increment, true);
        //     break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control :
            exitflag = this->Modified_Euler_Error_Control(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control :
            exitflag = this->Runge_Kutta_45_Error_Control(strain_increment);;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control_old :
            exitflag = this->Runge_Kutta_45_Error_Control_old(strain_increment);;
            break;
        default:
            cerr << "ASDPlasticMaterial3D::setTrialStrainIncr - Integration method not available!\n" ;
            exitflag = -1;
        }

        return exitflag;
    }

//==================================================================================================
//  Getters
//==================================================================================================

    const VoigtMatrix &getTangentTensor( void )
    {

        return Stiffness;
    }

    const Vector &getStress(void)
    {
        static Vector result(6);
        TrialStress.toStress(result);
        return result;
    }



    const Vector &getStrain(void)
    {
        static Vector result(6);
        TrialStrain.toStrain(result);
        return result;
    }

    const Vector &getPstrain(void)
    {
        static Vector result(6);
        TrialPlastic_Strain.toStress(result);
        return result;
    }

    const Vector &getEQPstrain(void)
    {
        static Vector result(1);
        result(0) = std::sqrt(2./3.*TrialPlastic_Strain.squaredNorm());
        // opserr << "JOSE: getEQPstrain(void) result = " << result << endln;
        return result;
    }

    const Vector &getPStress(void)
    {
        static Vector result(1);
        result(0) = TrialStress.meanStress();
        // opserr << "JOSE: getPStress(void) result = " << result << endln;
        return result;
    }

    const Vector &getJ2Stress(void)
    {
        static Vector result(1);
        result(0) = TrialStress.getJ2();
        // opserr << "JOSE: getJ2Stress(void) result = " << result << endln;
        return result;
    }

    const Vector &getVolStrain(void)
    {
        static Vector result(1);
        result(0) = TrialStrain.getI1();
        // opserr << "JOSE: getVolStrain(void) result = " << result << endln;
        return result;
    }

    const Vector &getJ2Strain(void)
    {
        static Vector result(1);
        result(0) = TrialStrain.getJ2();
        // opserr << "JOSE: getJ2Strain(void) result = " << result << endln;
        return result;
    }

    // Ladruno (ADR-97 wp/97b): how many Newton iterations the LAST Closest_Point
    // step took at this Gauss point (0 for an elastic step).  Exposed so ADR-97's
    // gate 1 can assert the quadratic <= 5 without parsing stdout -- the material
    // prints a page per construction and pytest's capfd cannot see the .pyd's
    // cout anyway (LEDGER_quirks).
    const Vector &getCPIterations(void)
    {
        static Vector result(1);
        result(0) = (double) cp_last_iterations;
        return result;
    }

    const Vector &getInternalVariableByPos(int pos)
    {
        static Vector return_vector(6);
        int find_pos = 0;

        // Note by J. Abell on Wed 06 Dec 2023 11:44:24
        //
        // The lambda capture of return_vector which is declared static above triggers
        // the following warning in g++ 11.4.0
        //
        // warning: capture of variable ‘return_vector’ with non-automatic storage duration
        //
        // Explanation:
        //   Because usually the lambda functions are used for threaded applications
        //   this warning is there to help with race conditions on the return_vector.
        //   In this case the warning is benign because opensees is not threaded
        //   at this level.
        //   We want to keep the static allocation of return_vector to avoid
        //   many calls to malloc (new) every time this function is called
        //   for performance reasons, so we have to live with the warning.
        iv_storage.apply([&pos, &find_pos, this](auto & internal_variable)
        {
            if (pos == find_pos)
            {
                auto &iv = internal_variable.trial_value;
                int iv_size = iv.size();
                return_vector.resize(iv_size);
                for (int i = 0; i < iv_size; ++i)
                {
                    return_vector(i) = iv(i);
                }
            }
            find_pos ++;
        });

        return return_vector;
    }

    const VoigtVector &getStressTensor( void )
    {
        return TrialStress;
    }

    const VoigtVector &getStrainTensor( void )
    {
        return TrialStrain;
    }

    const VoigtVector &getPlasticStrainTensor( void )
    {
        return TrialPlastic_Strain;
    }

    const VoigtVector  &getCommittedStressTensor(void)
    {
        return CommitStress;
    }

    const VoigtVector &getCommittedStrainTensor(void)
    {
        return CommitStrain;
    }

    const VoigtVector &getCommittedPlasticStrainTensor(void)
    {
        return CommitPlastic_Strain;
    }

    void ComputeTangentStiffness()
    {
        if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic)
        {
            VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
            Stiffness = Eelastic;
        }
    // ADR97_P4_MARKER:adr97_p4_guard_computetangentstiffness_dispatch
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_FirstOrder)
        {
            // Ladruno (ADR-97 wp/97e): a perturbed sub-call made BY
            // numerical_tangent_of_committed_map() re-enters this exact
            // dispatch through the SAME integrator; suppress_numerical_tangent
            // is set for exactly the duration of that sub-call so this branch
            // does not recurse. The sub-call only needs TrialStress, not
            // Stiffness, so it is safe to leave Stiffness untouched here.
            if (!suppress_numerical_tangent)
                compute_numerical_tangent_firstorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_SecondOrder)
        {
            if (!suppress_numerical_tangent)
                compute_numerical_tangent_secondorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Continuum)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = n.transpose() * Eelastic * m - hardening;

            VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
            Stiffness = Econtinuum;
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Secant)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = n.transpose() * Eelastic * m - hardening;

            VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
            Stiffness = (Econtinuum + Eelastic)/2;
        }
    }


    // Ladruno (ADR-94 wp/94c, M5): THE yield tolerance -- every |f| test in this
    // file goes through here.
    //
    // ADR-94 M5: `f_absolute_tol` is an ABSOLUTE tolerance in stress units, tested
    // against a |Phi| whose natural scale is sigma_y (VM), xi_c (DP), c*cos(phi)
    // (MC/MCTC), sigma_ci*s^a (HB) -- four orders across the catalogue before the
    // user picks a unit system.  Measured: the same MC problem completes 20/20 in
    // kPa and is refused on step 1 in Pa, same physics, same strains.
    //
    // `f_relative_tol` (default 0, i.e. OFF and byte-identical to before) adds a
    // unit-invariant floor:  tol = max(f_absolute_tol, f_relative_tol * scale),
    // with `scale` supplied by the yield function itself (YF_STRENGTH_SCALE).  A YF
    // that declares no scale returns 0 and is unaffected even with the option on.
    double yf_tolerance() const
    {
        const double abs_tol = DBL_OPT_f_absolute_tol[ASDP_TAG]; // KEEP_RAW_ABS_TOL (this IS the accessor)
        const double rel = DBL_OPT_f_relative_tol[ASDP_TAG];
        if (!(rel > 0.0))
            return abs_tol;
        double scale = yf.strength_scale(iv_storage, parameters_storage);
        if (scale < 0) scale = -scale;
        const double rel_tol = rel * scale;
        return (rel_tol > abs_tol) ? rel_tol : abs_tol;
    }

    // Ladruno (ADR-94 wp/94a): ONE admissibility gate for strict_convergence.
    //
    // ADR-94 B2/M2: `strict_convergence` only ever reached Backward_Euler, and the
    // "f merely DECREASED => call it elastic" shortcut exists in EIGHT places. Under
    // the flag, none of them may commit a state that is still outside the surface.
    // Returns true (and says so on opserr) when `stress_state` must be refused.
    // Flag off: constant-false, so every caller is byte-identical to upstream.
    //
    // NaN-safe on purpose: `!(yf_val <= tol)` is TRUE for NaN, so a NaN yield value
    // is a refusal, not an accept (ADR-94 B4).
    bool ladruno_strict_rejects(const char* where, const VoigtVector& stress_state) const
    {
        if (INT_OPT_strict_convergence[ASDP_TAG] == 0)
            return false;
        const double tol_yf = yf_tolerance();
        const double yf_val = yf(stress_state, iv_storage, parameters_storage);
        if (yf_val <= tol_yf)
            return false;
        opserr << "ASDPlasticMaterial3D::" << where << " (tag " << ASDP_TAG
               << ") - refusing to commit an inadmissible state: f = " << yf_val
               << " > yield tolerance = " << tol_yf
               << " -- rejecting step (strict_convergence)" << endln;
        return true;
    }
    // ADR97_P4_MARKER:adr97_p4_compute_local_stress_not_a_map
    // Ladruno (ADR-97 wp/97e, D2): NOT A MAP. Dead code after the P4
    // re-point below -- nothing calls this any more (verify with
    // `grep -n compute_local_stress` before ever wiring a new caller to it).
    // This is the simplified single-shot elastic-predictor / one-step
    // plastic-corrector that ADR-94 M3 and ADR-84 P4 identified as a THIRD
    // return map: it evaluates n/m/H once at the yield-crossing intersection
    // and takes one closed-form dLambda correction, with no resemblance to
    // Backward_Euler's cutting-plane Newton loop or Closest_Point's coupled
    // implicit solve. `Numerical_Algorithmic_FirstOrder/SecondOrder` used to
    // differentiate THIS function (measured 31%/4.6% off the true consistent
    // tangent, ADR-94 H6); they now differentiate the actual committed map
    // via numerical_tangent_of_committed_map(). Retained per ADR-97 D2 for
    // provenance / possible future standalone use -- do not resurrect it as
    // a tangent source without re-reading that decision.
    int compute_local_stress(
        const VoigtVector& local_stress, const VoigtVector& local_strain,
        const VoigtVector& strain_incr, VoigtVector& stress_incr) const
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // Initialize local variables for stress and strain increments
        VoigtVector depsilon = strain_incr;  // Strain increment (perturbation)
        VoigtVector dsigma = VoigtVector();  // Stress increment
        VoigtVector trial_stress = VoigtVector();  // Trial stress

        // Compute elastic stiffness matrix based on the local stress and parameters
        VoigtMatrix Eelastic = et(local_stress, parameters_storage);  // Elasticity tensor
        
        // Compute the elastic stress increment: dsigma = E * depsilon
        dsigma = Eelastic * depsilon;
        
        // Compute the trial stress
        trial_stress = local_stress + dsigma;

        // // Evaluate the yield function for the current stress and trial stress
        double yf_val_start = yf(local_stress, iv_storage, parameters_storage);
        double yf_val_end = yf(trial_stress, iv_storage, parameters_storage);

        // Check if the material response is elastic or plastic
        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) {
            // Elastic response: no plastic correction
            stress_incr = dsigma;  // Stress increment is purely elastic
            // Ladruno (ADR-94 wp/94a): the tangent-probe twin of the eight
            // f-decreasing exits. Its three callers ignore the return code, so
            // `stress_incr` is left assigned above and the numerical tangent is
            // unchanged; the code is reported for contract uniformity.
            if (ladruno_strict_rejects("compute_local_stress", trial_stress))
                return LADRUNO_MATERIAL_REFUSED;
        } else {
            // Plastic response: need to apply plastic correction
            // Compute plastic correction by finding the intersection of the yield surface
            VoigtVector intersection_stress = local_stress;
            VoigtVector intersection_strain = local_strain;
            depsilon_elpl = strain_incr;

            if (yf_val_start < 0) {
                // Find the intersection of the yield surface between the start and trial stress
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing(
                    local_stress, trial_stress, 0.0, 1.0, tol_yf);
                
                // Ensure the intersection factor is within valid bounds [0, 1]
                intersection_factor = std::max(0.0, std::min(1.0, intersection_factor));

                // Compute the intersection stress and strain
                intersection_stress = local_stress * (1 - intersection_factor) +
                                      trial_stress * intersection_factor;
                intersection_strain = local_strain + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * strain_incr;
            }

            // The trial stress is updated based on the intersection
            // stress_incr = intersection_stress - local_stress;
            trial_stress = intersection_stress;

            Eelastic = et(intersection_stress, parameters_storage);
            trial_stress  += Eelastic * depsilon_elpl;

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
            double den = n.transpose() * Eelastic * m - hardening;

            //Compute the plastic multiplier
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("m", m);
                printTensor1("n", n);
                cout << "hardening = " << hardening << endl;
                cout << "den = " << den << endl;
                printTensor1("depsilon_elpl", depsilon_elpl);
                return -1;
            }

            double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
            dLambda /= den;

            if (dLambda <= 0)
            {
                // cout << "CEP - dLambda = " << dLambda << " <= 0\n";
                // printTensor1("m", m);
                // printTensor1("n", n);
                // cout << "hardening = " << hardening << endl;
                // cout << "den = " << den << endl;
                // printTensor1("depsilon_elpl", depsilon_elpl);
                dLambda = 0;
            }

            trial_stress = trial_stress - dLambda * Eelastic * m;
        }

        stress_incr = trial_stress - CommitStress;

        return 0;  // Return success
    }



    // ADR97_P4_MARKER:adr97_p4_numalg_repoint_firstorder
    // Ladruno (ADR-97 wp/97e, D2 / closes ADR-84 P4 and ADR-94/ADR-97 M3's
    // "third map"): the CONSISTENT numerical tangent of the map this instance
    // actually commits -- whichever `integration_method` is configured
    // (`Backward_Euler`, `Closest_Point`, or in principle any other) -- taken
    // as a forward (this function) or central (compute_numerical_tangent_
    // secondorder) finite difference of `setTrialStrainIncr()` ITSELF, the
    // SAME dispatch every host element drives, instead of the simplified
    // single-shot `compute_local_stress()` a real element never calls. A
    // correct FD of this function is BY CONSTRUCTION the consistent tangent
    // up to stencil truncation error: there is no third algorithm here to be
    // wrong.
    //
    // Recursion: `Backward_Euler`/`Closest_Point` both call
    // `ComputeTangentStiffness()` at the end of every successful commit; a
    // perturbed sub-call below would otherwise re-enter that exact path and
    // try to compute ANOTHER numerical tangent of its own perturbed state,
    // unbounded. `suppress_numerical_tangent` is raised for the duration of
    // every sub-call; `ComputeTangentStiffness()` skips the `Numerical_
    // Algorithmic_*` branches while it is set (a sub-call only needs
    // `TrialStress`, never `Stiffness`).
    //
    // State: the integrator mutates `TrialStrain`, `TrialStress`,
    // `TrialPlastic_Strain`, every IV's `trial_value` (all of `iv_storage`),
    // `cp_last_iterations`, and the four `mutable` scratch buffers (`dsigma`,
    // `depsilon_elpl`, `intersection_stress`, `intersection_strain`). Several
    // early-return branches inside `Backward_Euler` (the Drucker-Prager apex
    // switch, the `special_return` hook switch, and the "PLASTIC
    // INCONSISTENCY" elastic-fallback exit) also assign `Stiffness`
    // UNCONDITIONALLY, bypassing the tangent-type dispatch entirely. All of
    // these are snapshotted once before the loop and restored after EVERY
    // perturbation -- both integrators already reset `TrialPlastic_Strain`/
    // `iv_storage` from Commit* at their own top (`revert_all()`), so nothing
    // would actually accumulate across perturbations even without the
    // per-iteration restore, but restoring immediately keeps the instance
    // well-defined if a LATER perturbation in the same loop refuses.
    //
    // The FD is assembled into a LOCAL matrix, never into `tangent_matrix`
    // (== `this->Stiffness`, passed by reference from
    // `ComputeTangentStiffness()`) column-by-column: the unconditional
    // `Stiffness = ...` branches named above would otherwise clobber columns
    // already written by an earlier perturbation. `tangent_matrix` is
    // assigned exactly once, after the loop and after the final restore.
    //
    // Refusals: a nonzero return from any perturbed `setTrialStrainIncr()`
    // call aborts the loop immediately, restores the snapshot, and
    // propagates that SAME code outward -- never assembles a partial tangent
    // from whatever columns happened to finish.
    //
    // Cost: 6 (first order) or 12 (second order) EXTRA full integrations per
    // call (the forward-difference baseline is the real, already-computed
    // `TrialStress` sitting in the snapshot -- no extra unperturbed call is
    // needed, unlike the old compute_local_stress()-based version).
    // `Backward_Euler`'s singular-tangent/NaN/plastic-inconsistency
    // diagnostics are plain `cout <<`, not `opserr`-gated, so a perturbation
    // landing near one of those guards multiplies that chatter 6x-12x.
    int numerical_tangent_of_committed_map(const VoigtVector& strain_incr,
                                            VoigtMatrix& tangent_matrix,
                                            bool second_order,
                                            double epsilon_ref = 1e-8,
                                            double delta_min = 1e-12)
    {
        // -------- snapshot the real (unperturbed) converged state ----------
        const VoigtVector snap_TrialStrain          = TrialStrain;
        const VoigtVector snap_TrialStress          = TrialStress;
        const VoigtVector snap_TrialPlastic_Strain  = TrialPlastic_Strain;
        const iv_storage_t snap_iv_storage          = iv_storage;
        const int          snap_cp_last_iterations  = cp_last_iterations;
        const VoigtVector snap_dsigma               = dsigma;
        const VoigtVector snap_depsilon_elpl        = depsilon_elpl;
        const VoigtVector snap_intersection_stress  = intersection_stress;
        const VoigtVector snap_intersection_strain  = intersection_strain;
        const VoigtMatrix snap_Stiffness            = Stiffness;

        auto restore_snapshot = [&]()
        {
            TrialStrain          = snap_TrialStrain;
            TrialStress          = snap_TrialStress;
            TrialPlastic_Strain  = snap_TrialPlastic_Strain;
            iv_storage           = snap_iv_storage;
            cp_last_iterations   = snap_cp_last_iterations;
            dsigma               = snap_dsigma;
            depsilon_elpl        = snap_depsilon_elpl;
            intersection_stress  = snap_intersection_stress;
            intersection_strain  = snap_intersection_strain;
            Stiffness            = snap_Stiffness;
        };

        const int n = 6;
        VoigtMatrix local_tangent = VoigtMatrix::Zero();
        const double delta = std::max(epsilon_ref * strain_incr.norm(), delta_min);
        const VoigtVector unperturbed_stress = snap_TrialStress;

        suppress_numerical_tangent = true;
        int rc = 0;
        for (int i = 0; i < n && rc == 0; ++i)
        {
            VoigtVector strain_incr_p1 = strain_incr;
            strain_incr_p1(i) += delta;
            rc = this->setTrialStrainIncr(strain_incr_p1);
            if (rc != 0) break;
            const VoigtVector perturbed_stress_p1 = TrialStress;
            restore_snapshot();

            if (second_order)
            {
                VoigtVector strain_incr_p2 = strain_incr;
                strain_incr_p2(i) -= delta;
                rc = this->setTrialStrainIncr(strain_incr_p2);
                if (rc != 0) break;
                const VoigtVector perturbed_stress_p2 = TrialStress;
                restore_snapshot();

                for (int j = 0; j < n; ++j)
                    local_tangent(j, i) = (perturbed_stress_p1(j) - perturbed_stress_p2(j)) / (2.0 * delta);
            }
            else
            {
                for (int j = 0; j < n; ++j)
                    local_tangent(j, i) = (perturbed_stress_p1(j) - unperturbed_stress(j)) / delta;
            }
        }
        suppress_numerical_tangent = false;
        restore_snapshot();

        if (rc != 0)
            return rc;

        tangent_matrix = local_tangent;
        return 0; // Return success
    }

    // Ladruno (ADR-97 wp/97e): thin wrapper kept for source compatibility --
    // now differentiates the ACTUAL committed map (forward difference) via
    // numerical_tangent_of_committed_map(), instead of compute_local_stress().
    int compute_numerical_tangent_firstorder(
        const VoigtVector& strain_incr, VoigtMatrix& tangent_matrix, double epsilon_ref = 1e-8, double delta_min = 1e-12)
    {
        return numerical_tangent_of_committed_map(strain_incr, tangent_matrix, false, epsilon_ref, delta_min);
    }



    // ADR97_P4_MARKER:adr97_p4_numalg_repoint_secondorder
    // Ladruno (ADR-97 wp/97e): thin wrapper kept for source compatibility --
    // now differentiates the ACTUAL committed map (central difference) via
    // numerical_tangent_of_committed_map(), instead of compute_local_stress().
    // See that function's doc comment for the full design (recursion guard,
    // state snapshot/restore, refusal propagation, cost).
    int compute_numerical_tangent_secondorder(
        const VoigtVector& strain_incr, VoigtMatrix& tangent_matrix, double epsilon_ref = 1e-8, double delta_min = 1e-12)
    {
        return numerical_tangent_of_committed_map(strain_incr, tangent_matrix, true, epsilon_ref, delta_min);
    }

    const Matrix& getTangent()
    {
        static Matrix return_matrix(6, 6);

        // Ladruno (ADR-94 wp/94b, M1): Stiffness is per-instance now. An instance whose
        // integrator has never run would hand the assembler an exactly-zero (singular)
        // block, where the old class-static happened to carry whatever the last
        // instance to integrate had left in it. Fall back to the elastic tangent at the
        // committed stress -- the only defensible tangent for a state nobody has
        // integrated yet.
        if (Stiffness.isZero(0.0))
        {
            Stiffness = et(CommitStress, parameters_storage);
        }

        copyToMatrixReference(Stiffness, return_matrix);

        return return_matrix;
    }


    const Matrix& getInitialTangent()
    {
        static Matrix return_matrix(6, 6);

        // Ladruno (ADR-94 wp/94b, F1): this was `Stiffness = Eelastic;` -- a "getter"
        // that clobbered the shared tangent every other instance's later getTangent()
        // would read. Compute into a local and copy that to the return buffer; the
        // member state is left alone.
        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);

        copyToMatrixReference(Eelastic, return_matrix);

        return return_matrix;
    }



//==================================================================================================
//  State commiting and reversion
//==================================================================================================


    int commitState(void)
    {

        CommitStress = TrialStress;
        CommitStrain = TrialStrain;
        CommitPlastic_Strain = TrialPlastic_Strain;

        iv_storage.commit_all();

        if (first_step)
        {
            first_step = false;
        }

        if (GLOBAL_INT_max_iter[ASDP_TAG] > 0 || GLOBAL_DBL_max_error[ASDP_TAG] > 0.)
        {
            cout << "  () ASDP Integration Info. Tag = " << ASDP_TAG << " max_iter = " << GLOBAL_INT_max_iter[ASDP_TAG] << " max_error = " << GLOBAL_DBL_max_error[ASDP_TAG] << endl;
            GLOBAL_INT_max_iter[ASDP_TAG] = 0;
            GLOBAL_DBL_max_error[ASDP_TAG] = 0.;
        }

        return 0;
    }

    //Reverts the commited variables to the trials and calls revert on BET Classes.
    int revertToLastCommit(void)
    {
        // Ladruno (ADR-94 wp/94b, M6/H4): every statement in this body used to be
        // commented out, so a step that failed to converge left the dirty trial state
        // behind while the method reported success. Domain::revertToLastCommit() calls
        // this on every element of a failed step; revert means revert.
        TrialStress = CommitStress;
        TrialStrain = CommitStrain;
        TrialPlastic_Strain = CommitPlastic_Strain;

        iv_storage.revert_all();

        // No committed tangent is stored anywhere, so the only tangent consistent with
        // the committed state is the elastic one evaluated at the committed stress.
        Stiffness = et(CommitStress, parameters_storage);

        dsigma.setZero();
        depsilon_elpl.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

        if (GLOBAL_INT_max_iter[ASDP_TAG] > 0 || GLOBAL_DBL_max_error[ASDP_TAG] > 0.)
        {
            GLOBAL_INT_max_iter[ASDP_TAG] = 0;
            GLOBAL_DBL_max_error[ASDP_TAG] = 0.;
        }

        return 0;
    }

    int revertToStart(void)
    {
        // Ladruno (ADR-94 wp/94b, M6/H4): was `cerr << "not implemented"; return -1;`,
        // and Domain::revertToStart() / OPS_resetModel() both discard that -1, so
        // ops.reset() left the material's committed state alive underneath a
        // zeroed geometry (contract doc S4: "a third, inconsistent number").
        // Restore exactly the freshly-constructed, freshly-parsed state.
        TrialStress.setZero();
        CommitStress.setZero();
        TrialStrain.setZero();
        CommitStrain.setZero();
        TrialPlastic_Strain.setZero();
        CommitPlastic_Strain.setZero();

        if (initial_iv_captured)
        {
            iv_storage = iv_storage_initial;
        }
        iv_storage.revert_all();

        dsigma.setZero();
        depsilon_elpl.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();
        Stiffness.setZero();

        // first_step = true re-arms the InitialP0 geostatic seed in setTrialStrain(),
        // which is what a fresh instance would do on its own first strain.
        first_step = true;
        stress_set_externally = false;

        GLOBAL_INT_max_iter[ASDP_TAG] = 0;
        GLOBAL_DBL_max_error[ASDP_TAG] = 0.;

        return 0;
    }

    NDMaterial *getCopy(void)
    {

        ASDPlasticMaterial3D <ElasticityType,
                           YieldFunctionType,
                           PlasticFlowType,
                           thisClassTag> *newmaterial = new ASDPlasticMaterial3D<ElasticityType,
        YieldFunctionType,
        PlasticFlowType,
        thisClassTag>(ASDP_TAG);
        newmaterial->TrialStrain = this->TrialStrain;
        newmaterial->TrialStress = this->TrialStress;
        newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
        newmaterial->CommitStress = this->CommitStress;
        newmaterial->CommitStrain = this->CommitStrain;
        newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
        newmaterial->iv_storage = this->iv_storage;
        newmaterial->parameters_storage = this->parameters_storage;
        newmaterial->stress_set_externally = this->stress_set_externally; // Ladruno (HB/StiffSoil integration, ledger row 337)
        // Ladruno (ADR-94 wp/94b, H14/F6): first_step was never copied, so a copy made
        // from an already-advanced instance silently re-armed the InitialP0 seed.
        newmaterial->first_step = this->first_step;
        newmaterial->iv_storage_initial = this->iv_storage_initial;
        newmaterial->initial_iv_captured = this->initial_iv_captured;

        return newmaterial;
    }


    NDMaterial *getCopy(const char *type) {
        if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
            ASDPlasticMaterial3D <ElasticityType,
                               YieldFunctionType,
                               PlasticFlowType,
                               thisClassTag> *newmaterial = new ASDPlasticMaterial3D<ElasticityType,
            YieldFunctionType,
            PlasticFlowType,
            thisClassTag>(ASDP_TAG);
            newmaterial->TrialStrain = this->TrialStrain;
            newmaterial->TrialStress = this->TrialStress;
            newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
            newmaterial->CommitStress = this->CommitStress;
            newmaterial->CommitStrain = this->CommitStrain;
            newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
            newmaterial->iv_storage = this->iv_storage;
            newmaterial->parameters_storage = this->parameters_storage;
            newmaterial->stress_set_externally = this->stress_set_externally; // Ladruno (HB/StiffSoil integration, ledger row 337)
        // Ladruno (ADR-94 wp/94b, H14/F6): first_step was never copied, so a copy made
        // from an already-advanced instance silently re-armed the InitialP0 seed.
        newmaterial->first_step = this->first_step;
        newmaterial->iv_storage_initial = this->iv_storage_initial;
        newmaterial->initial_iv_captured = this->initial_iv_captured;

            return newmaterial;
        } else
        {
            cout << "ASDPlasticMaterial3D::getCopy(const char *type) - Only 3D is currently supported. Use 3D elements!" << endl;
        }
        return 0;
    }


    int setParameter(const char **argv, int argc, Parameter &param)
    {

        cout << "ASDPlasticMaterial3D::setParameter  argv = " << *argv << endl;

        // if (argc < 2)
        //     return -1;
        
        // int theMaterialTag;
        // theMaterialTag = atoi(argv[1]);
        
        // if (theMaterialTag == this->getTag()) {
        if (true) {
            
            // State variables (use specific response IDs)
            if (strcmp(argv[0], "stress") == 0) {
                return param.addObject(1, this);
            }
            else if (strcmp(argv[0], "strain") == 0) {
                return param.addObject(2, this);
            }
            else if (strcmp(argv[0], "plasticStrain") == 0) {
                return param.addObject(3, this);
            }
            else if (strcmp(argv[0], "trialStress") == 0) {
                return param.addObject(4, this);
            }
            else if (strcmp(argv[0], "trialStrain") == 0) {
                return param.addObject(5, this);
            }
            else if (strcmp(argv[0], "trialPlasticStrain") == 0) {
                return param.addObject(6, this);
            }
            else if (strcmp(argv[0], "K02D") == 0) {
                cout << "       ---->  K02D" << endl;
                return param.addObject(7, this);
            }
            else if (strcmp(argv[0], "K03D") == 0) {
                cout << "       ---->  K03D" << endl;
                return param.addObject(8, this);
            }
            // Ladruno (HB/StiffSoil integration, ledger row 337): register stress-increment setParameter tokens
            else if (strcmp(argv[0], "trialStressIncrement") == 0) {
                return param.addObject(9, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXX") == 0) {
                return param.addObject(10, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementYY") == 0) {
                return param.addObject(11, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementZZ") == 0) {
                return param.addObject(12, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXY") == 0) {
                return param.addObject(13, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementYZ") == 0) {
                return param.addObject(14, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXZ") == 0) {
                return param.addObject(15, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXX") == 0) {
                return param.addObject(16, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementYY") == 0) {
                return param.addObject(17, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementZZ") == 0) {
                return param.addObject(18, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXY") == 0) {
                return param.addObject(19, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementYZ") == 0) {
                return param.addObject(20, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXZ") == 0) {
                return param.addObject(21, this);
            }
            else {
                // For all other parameter names, use the parameter system to pass the name
                // Store the parameter name in the Parameter object (if supported)
                // or use a generic response ID
                current_parameter_name = argv[0]; // Store for use in updateParameter
                return param.addObject(1000, this);
            }
        }
        
        return -1;
    }

    int updateParameter(int responseID, Information &info)
    {

        cout << "ASDPlasticMaterial3D::updateParameter  responseID = " << responseID << endl;

        // Ladruno (HB/StiffSoil integration, ledger row 337): debug-print the Information payload
        opserr << " info = "; // << info << endln;
        info.Print(opserr);

        // State variables (committed values)
        if (responseID == 1) { // stress
            if (info.theType == VectorType) {
                const Vector& newStress = *(info.theVector);
                CommitStress = VoigtVector::fromStress(newStress);
                TrialStress = CommitStress;
                stress_set_externally = true; // Ladruno (HB/StiffSoil integration, ledger row 337)
            }
            return 0;
        }
        else if (responseID == 2) { // strain
            if (info.theType == VectorType) {
                const Vector& newStrain = *(info.theVector);
                CommitStrain = VoigtVector::fromStrain(newStrain);
                TrialStrain = CommitStrain;
            }
            return 0;
        }
        else if (responseID == 3) { // plasticStrain
            if (info.theType == VectorType) {
                const Vector& newPlasticStrain = *(info.theVector);
                CommitPlastic_Strain = VoigtVector::fromStrain(newPlasticStrain);
                TrialPlastic_Strain = CommitPlastic_Strain;
            }
            return 0;
        }
        // Trial state variables
        else if (responseID == 4) { // trialStress
            if (info.theType == VectorType) {
                const Vector& newTrialStress = *(info.theVector);
                TrialStress = VoigtVector::fromStress(newTrialStress);
                stress_set_externally = true; // Ladruno (HB/StiffSoil integration, ledger row 337)
            }
            return 0;
        }
        else if (responseID == 5) { // trialStrain
            if (info.theType == VectorType) {
                const Vector& newTrialStrain = *(info.theVector);
                TrialStrain = VoigtVector::fromStrain(newTrialStrain);
            }
            return 0;
        }
        else if (responseID == 6) { // trialPlasticStrain
            if (info.theType == VectorType) {
                const Vector& newTrialPlasticStrain = *(info.theVector);
                TrialPlastic_Strain = VoigtVector::fromStrain(newTrialPlasticStrain);
            }
            return 0;
        }
        else if (responseID == 7) { // K02D - Use Sigma_Y
            // cout << "responseID == 7 !! info.theType = " << info.theType << " DoubleType = " << DoubleType << endl;
            // if (info.theType == DoubleType) {
                const double& K02D = info.theDouble;
                cout << "ASDPL @ tag = " << this->getTag() << " K02D  K0 = " << K02D << endl;
                CommitStress(0) = K02D * CommitStress(1);
                CommitStress(2) = K02D * CommitStress(1);
                stress_set_externally = true; // Ladruno (HB/StiffSoil integration, ledger row 337)
            // }
            return 0;
        }
        else if (responseID == 8) { // K03D - Use Sigma_Z
            // if (info.theType == DoubleType) {
                const double& K03D = info.theDouble;
                cout << "ASDPL @ tag = " << this->getTag() << " K03D  K0 = " << K03D << endl;
                CommitStress(0) = K03D * CommitStress(2);
                CommitStress(1) = K03D * CommitStress(2);
                stress_set_externally = true; // Ladruno (HB/StiffSoil integration, ledger row 337)
            // }
            return 0;
        }
        // Ladruno (HB/StiffSoil integration, ledger row 337): responseID 9-21, trial/commit stress-increment updateParameter handlers
        else if (responseID == 9) { // trialStressIncrement
            if (info.theType == VectorType) {
                const Vector& newTrialStress = *(info.theVector);
                opserr << "ASDPL @ tag = " << this->getTag() << "  newTrialStress = " << newTrialStress   << endln;
                TrialStress += VoigtVector::fromStress(newTrialStress);
                stress_set_externally = true;
            }
            return 0;
        }
        else if (responseID == 10) { // trialStressIncrementXX
            TrialStress(0) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 11) { // trialStressIncrementYY
            TrialStress(1) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 12) { // trialStressIncrementZZ
            TrialStress(2) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 13) { // trialStressIncrementXY
            TrialStress(3) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 14) { // trialStressIncrementYZ
            TrialStress(4) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 15) { // trialStressIncrementXZ
            TrialStress(5) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 16) { // commitStressIncrementXX
            CommitStress(0) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 17) { // commitStressIncrementYY
            CommitStress(1) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 18) { // commitStressIncrementZZ
            CommitStress(2) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 19) { // commitStressIncrementXY
            CommitStress(3) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 20) { // commitStressIncrementYZ
            CommitStress(4) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 21) { // commitStressIncrementXZ
            CommitStress(5) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        // Generic parameter update (model parameters and internal variables)
        else if (responseID == 1000) {
            // Use the stored parameter name from the most recent setParameter call
            const char* param_name = current_parameter_name.c_str();
            double param_value = info.theDouble;
            
            // Try to set as model parameter first
            // utuple_storage::setParameterByName handles non-existent parameters gracefully (does nothing)
            parameters_storage.setParameterByName(param_name, param_value);
            
            // Also try to set as internal variable (for scalar internal variables)
            // utuple_storage::setInternalVariableByName also handles non-existent variables gracefully
            iv_storage.setInternalVariableByName(param_name, 1, &param_value);
            
            return 0;
        }
        
        return -1;
    }



    Response *setResponse (const char **argv, int argc,
                           OPS_Stream & s)
    {

        static Vector return_vector(6);

        if (argc < 1)
            return 0;

        // Ladruno: emit NdMaterialOutput + ResponseType XML so recorders
        // (MPCO / LadrunoRecorder) label the per-Gauss-point material columns
        // with real component names instead of the generic C1..CN fallback.
        // Requested via `-E material.<token>` (the `material.` prefix makes the
        // recorder iterate Gauss points and forward here); a bare `-E <token>`
        // never reaches the material and records nothing.
        s.tag("NdMaterialOutput");
        s.attr("matType", this->getClassType());
        s.attr("matTag", this->getTag());

        static const char *TENSOR_STRESS[6] =
            {"sigma11", "sigma22", "sigma33", "sigma12", "sigma23", "sigma13"};
        static const char *TENSOR_STRAIN[6] =
            {"eps11", "eps22", "eps33", "eps12", "eps23", "eps13"};
        static const char *TENSOR_PSTRAIN[6] =
            {"epsP11", "epsP22", "epsP33", "epsP12", "epsP23", "epsP13"};

        if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0) {
            for (int i = 0; i < 6; ++i) s.tag("ResponseType", TENSOR_STRESS[i]);
            return new MaterialResponse(this, 1, this->getStress());
        }
        else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
            for (int i = 0; i < 6; ++i) s.tag("ResponseType", TENSOR_STRAIN[i]);
            return new MaterialResponse(this, 2, this->getStrain());
        }
        else if (strcmp(argv[0], "pstrain") == 0 || strcmp(argv[0], "pstrains") == 0) {
            for (int i = 0; i < 6; ++i) s.tag("ResponseType", TENSOR_PSTRAIN[i]);
            return new MaterialResponse(this, 3, this->getPstrain());
        }
        else if (strcmp(argv[0], "eqpstrain") == 0 ) {
            s.tag("ResponseType", "eqpstrain");
            return new MaterialResponse(this, 4, this->getEQPstrain());
        }
        else if (strcmp(argv[0], "PStress") == 0 ) {
            s.tag("ResponseType", "p");                  // mean (hydrostatic) stress
            return new MaterialResponse(this, 5, this->getPStress());
        }
        else if (strcmp(argv[0], "J2Stress") == 0 ) {
            s.tag("ResponseType", "J2stress");           // 2nd deviatoric stress invariant
            return new MaterialResponse(this, 6, this->getJ2Stress());
        }
        else if (strcmp(argv[0], "VolStrain") == 0 ) {
            s.tag("ResponseType", "epsVol");             // volumetric strain (I1 of strain)
            return new MaterialResponse(this, 7, this->getVolStrain());
        }
        else if (strcmp(argv[0], "J2Strain") == 0 ) {
            s.tag("ResponseType", "J2strain");           // 2nd deviatoric strain invariant
            return new MaterialResponse(this, 8, this->getJ2Strain());
        }
        else if (strcmp(argv[0], "cp_iterations") == 0 ) {
            // Ladruno (ADR-97 wp/97b)
            s.tag("ResponseType", "cp_iterations");
            return new MaterialResponse(this, 9, this->getCPIterations());
        }
        else
        {
            const char *iv_name = argv[0];

            int iv_size = this->getInternalVariableSizeByName(iv_name);
            int pos = this->getInternalVariableIndexByName(iv_name);

            // Ladruno: an unrecognized token (pos < 0) previously fell through to
            // MaterialResponse(1000+pos, Vector(iv_size)) with iv_size == -1 -- a
            // malformed response. Return 0 so the recorder records nothing for the
            // bad token (e.g. `material.plasticStrain`, which is spelled `pstrain`
            // here) instead of a garbage bucket.
            if (pos < 0 || iv_size <= 0)
                return 0;

            // label each internal-variable component <ivName> or <ivName>_<i>
            if (iv_size == 1) {
                s.tag("ResponseType", iv_name);
            } else {
                char buf[64];
                for (int i = 0; i < iv_size; ++i) {
                    snprintf(buf, sizeof(buf), "%s_%d", iv_name, i + 1);
                    s.tag("ResponseType", buf);
                }
            }
            return new MaterialResponse(this, 1000 + pos, Vector(iv_size));
        }

        return 0;
    }


    int getResponse (int responseID, Information & matInformation)
    {


        if (matInformation.theVector == 0)
            return 0;

        if (responseID == -1)
        {
            return -1;
        }
        else if (responseID == 1)
        {
            *(matInformation.theVector) = getStress();
        }
        else if (responseID == 2)
            *(matInformation.theVector) = getStrain();
        else if (responseID == 3)
            *(matInformation.theVector) = getPstrain();
        else if (responseID == 4)
            *(matInformation.theVector) = getEQPstrain();        
        else if (responseID == 5)
            *(matInformation.theVector) = getPStress();        
        else if (responseID == 6)
            *(matInformation.theVector) = getJ2Stress();        
        else if (responseID == 7)
            *(matInformation.theVector) = getVolStrain();        
        else if (responseID == 8)
            *(matInformation.theVector) = getJ2Strain();
        else if (responseID == 9)   // Ladruno (ADR-97 wp/97b)
            *(matInformation.theVector) = getCPIterations();
        else if (responseID >= 1000)
        {
            int pos = responseID - 1000;
            *(matInformation.theVector) = getInternalVariableByPos(pos);
        }

        return 0;
    }

    const char *getType(void) const {return "ThreeDimensional";}

    int sendSelf(int commitTag, Channel & theChannel)
    {
        cerr << "ASDPlasticMaterial3D::sendSelf - not implemented!!!\n" ;


        return 0;
    }

    int recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
    {
        cerr << "ASDPlasticMaterial3D::recvSelf - not implemented!!!\n" ;


        return 0;
    }

    void Print(OPS_Stream & s, int flag = 0) {
        s << "ASDPlasticMaterial3D" << endln;
        s << "  Yield Function          : " << yf.NAME << endln;
        s << "  Plastic flow direction  : " << pf.NAME << endln;
        s << "  Elasticity              : " << et.NAME << endln;
        s << "  # of Internal variables : " << (int) iv_storage.size() <<  endln;
        iv_storage.print_components();
        s << "  # of Parameters         : " << (int) parameters_storage.size() <<  endln;
        parameters_storage.print_components();
    }

    void Print(ostream & s, int flag = 0)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        s << "ASDPlasticMaterial3D" << endl;
        s << "  Yield Function          : " << yf.NAME << endl;
        s << "  Plastic flow direction  : " << pf.NAME << endl;
        s << "  Elasticity              : " << et.NAME << endl;
        s << "  # of Internal variables : " << iv_storage.size() <<  endl;
        iv_storage.print_components();
        s << "  # of Parameters         : " << parameters_storage.size() <<  endl;
        parameters_storage.print_components();

    }

    int getObjectSize()
    {
        int size = 0;

        // 6 3x3 VoigtVectors and 1 VoigtMatrix (3x3x3x3)
        size += (3 * 3 * 6 + 3 * 3 * 3 * 3) * sizeof(double);

        //Four pointers
        size += 4 * sizeof(YieldFunctionType*);

        //Whatever the base components size is
        size += sizeof(yf);//yf->getObjectSize();
        size += sizeof(et);//et->getObjectSize();
        size += sizeof(pf);//pf->getObjectSize();
        // size += sizeof(internal_variables);//internal_variables->getObjectSize();
        size += sizeof(NDMaterial);

        // size += static_cast<T*>(this)->getObjectSize();

        return size;
    }

    void setTrialStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        TrialStress = stress;
    }

    bool set_constitutive_integration_method(int method, int tangent, double f_absolute_tol, double stress_absolute_tol, int n_max_iterations, int return_to_yield_surface, int rk45_niter_max, double rk45_dT_min, int strict_convergence = 0, double f_relative_tol = 0.0) // Ladruno (ADR-84 P2a): opt-in strict_convergence flag, default 0 preserves upstream behavior. Ladruno (ADR-94 wp/94c, M5): opt-in f_relative_tol, default 0 = off
    {
        // Ladruno (ADR-94 wp/94a): ADR-94 M8 -- Forward_Euler_Crisfield,
        // Multistep_Forward_Euler, Multistep_Forward_Euler_Crisfield and
        // Full_Backward_Euler were ACCEPTED here but have NO case in the
        // setTrialStrain dispatch switch, so selecting one silently fell through to
        // the switch default. Only values that actually dispatch are accepted.
        if ( method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Not_Set
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control_old
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point   // Ladruno (ADR-97 wp/97b)
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch)
        {
            INT_OPT_constitutive_integration_method[ASDP_TAG] = (ASDPlasticMaterial3D_Constitutive_Integration_Method) method ;
            INT_OPT_tangent_operator_type[ASDP_TAG] = (ASDPlasticMaterial3D_Tangent_Operator_Type) tangent ;
            DBL_OPT_f_absolute_tol[ASDP_TAG] = f_absolute_tol ;
            DBL_OPT_stress_absolute_tol[ASDP_TAG] = stress_absolute_tol ;
            INT_OPT_n_max_iterations[ASDP_TAG] = n_max_iterations ;
            INT_OPT_return_to_yield_surface[ASDP_TAG] = return_to_yield_surface ;
            DBL_OPT_RK45_dT_min[ASDP_TAG] = rk45_dT_min ;
            INT_OPT_RK45_niter_max[ASDP_TAG] = rk45_niter_max ;
            INT_OPT_strict_convergence[ASDP_TAG] = strict_convergence ; // Ladruno (ADR-84 P2a)
            DBL_OPT_f_relative_tol[ASDP_TAG] = f_relative_tol ; // Ladruno (ADR-94 wp/94c, M5)

            GLOBAL_INT_max_iter[ASDP_TAG] = 0;
            GLOBAL_DBL_max_error[ASDP_TAG] = 0.;

            cout << "set_constitutive_integration_method tag = " << ASDP_TAG << " ::: " << endl;
            cout << "   method = " << method << endl;
            cout << "   tanget_type = " << tangent << endl;
            cout << "   f_absolute_tol = " << f_absolute_tol << endl;
            cout << "   stress_absolute_tol = " << stress_absolute_tol << endl;
            cout << "   n_max_iterations = " << n_max_iterations << endl;
            cout << "   return_to_yield_surface = " << return_to_yield_surface << endl;
            cout << "   rk45_niter_max = " << rk45_niter_max << endl;
            cout << "   rk45_dT_min = " << rk45_dT_min << endl;
            cout << "   strict_convergence = " << strict_convergence << endl; // Ladruno (ADR-84 P2a)
            cout << "   f_relative_tol = " << f_relative_tol << endl; // Ladruno (ADR-94 wp/94c, M5)

            return true;
        }
        else
        {
            // Ladruno (ADR-94 wp/94a): cerr is invisible under most OpenSees front ends.
            opserr << "ASDPlasticMaterial3D::set_constitutive_integration_method - refusing "
                   << "constitutive_integration_method " << method
                   << " (unknown, or an enum value with no dispatch case -- ADR-94 M8)" << endln;
            return false;
        }
    }

protected:

    void setTrialPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        TrialPlastic_Strain = strain;
    }

    void setCommitStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStress = stress;
    }

    void setCommitStrain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStrain = strain;
    }

    void setCommitPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitPlastic_Strain = strain;
    }

    void setStiffness(const VoigtMatrix & stiff)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        Stiffness = stiff;
    }


    void setStressTensor(VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStress = stress;
        TrialStress = stress;
        return;
    }


private:


    int Forward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;



        int errorcode = -1;

        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        dsigma = Eelastic * depsilon;

        TrialStress = sigma + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(sigma, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Forward_Euler", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = epsilon  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;

            Eelastic = et(intersection_stress, parameters_storage);
            TrialStress  += Eelastic * depsilon_elpl;

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
            double den = n.transpose() * Eelastic * m - hardening;

            //Compute the plastic multiplier
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("m", m);
                printTensor1("n", n);
                cout << "hardening = " << hardening << endl;
                cout << "den = " << den << endl;
                printTensor1("depsilon_elpl", depsilon_elpl);
                return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }

            double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
            dLambda /= den;

            if (dLambda <= 0)
            {
                // cout << "CEP - dLambda = " << dLambda << " <= 0\n";
                // printTensor1("m", m);
                // printTensor1("n", n);
                // cout << "hardening = " << hardening << endl;
                // cout << "den = " << den << endl;
                // printTensor1("depsilon_elpl", depsilon_elpl);
                dLambda = 0;
            }

            // Update the trial plastic strain.
            TrialPlastic_Strain += dLambda * m;

            // This code iterates internal variables and updates the trial values
            iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                internal_variable.trial_value += dLambda * h;
            });


            // Deal with the APEX if needed
            // if constexpr (yf_has_apex<YieldFunctionType>::value) {
            // {
            //    // TODO
            // }

            //Correct the trial stress
            TrialStress = TrialStress - dLambda * Eelastic * m;


            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     n_corr.transpose() * Eelastic * m_corr - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * m_corr;
                TrialPlastic_Strain += dLambda_corr * m_corr;
            }
            else
            {
                // Do nothing
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                // cout << "Numeric error!\n";
                // printTensor1("TrialStress = " , TrialStress);
                // printTensor1("CommitStress = " , CommitStress);
                // printTensor1("depsilon = " , depsilon);
                // printTensor1("dsigma   = " , dsigma);
                // printTensor1("intersection_stress = " , intersection_stress);
                // printTensor2("Eelastic = " , Eelastic);
                // printTensor2("Stiffness = " , Stiffness);
                // cout << "yf_val_start = " << yf_val_start << endl;
                // cout << "yf_val_end = " << yf_val_end << endl;
                // printTensor1("n = " , n );
                // printTensor1("m = " , m );
                // cout << "hardening  = " << hardening << endl;
                // cout << "den = " << den << endl;
                // cout << "dLambda = " << dLambda << endl;

                return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }
            else
            {
                // Ladruno (ADR-94 wp/94a)
                if (ladruno_strict_rejects("Forward_Euler (post return-to-yield)", TrialStress))
                    return LADRUNO_MATERIAL_REFUSED;
                return 0;
            }

            ComputeTangentStiffness();

        }

        return errorcode;
    }



    // int Forward_Euler_Subincrement(const VoigtVector &strain_incr, bool const& with_return2yield_surface)
    int Forward_Euler_Subincrement(const VoigtVector & strain_incr)
    {
       using namespace ASDPlasticMaterial3DGlobals;



        int errorcode = -1;

        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        dsigma = Eelastic * depsilon;

        TrialStress = sigma + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(sigma, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Forward_Euler_Subincrement", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = epsilon  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;

            int Nsubsteps = INT_OPT_n_max_iterations[ASDP_TAG];
            depsilon_elpl = depsilon_elpl/Nsubsteps;

            for (int substep = 0; substep < Nsubsteps; ++substep)
            {
          
                Eelastic = et(intersection_stress, parameters_storage);
                TrialStress  += Eelastic * depsilon_elpl;

                //Compute normal to YF (n) and Plastic Flow direction (m)
                const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

                double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
                double den = n.transpose() * Eelastic * m - hardening;

                //Compute the plastic multiplier
                if (abs(den) < MACHINE_EPSILON)
                {
                    cout << "CEP - den = 0\n";
                    cout << "yf_val_start = " << yf_val_start << endl;
                    cout << "yf_val_end = " << yf_val_end << endl;
                    printTensor1("m", m);
                    printTensor1("n", n);
                    cout << "hardening = " << hardening << endl;
                    cout << "den = " << den << endl;
                    printTensor1("depsilon_elpl", depsilon_elpl);
                    return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
                }

                double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
                dLambda /= den;

                if (dLambda <= 0)
                {
                    dLambda = 0;
                }

                // Update the trial plastic strain.
                TrialPlastic_Strain += dLambda * m;

                // This code iterates internal variables and updates the trial values
                iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
                {
                    auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                    internal_variable.trial_value += dLambda * h;
                });

                //Correct the trial stress
                TrialStress = TrialStress - dLambda * Eelastic * m;
            }

            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m_corr,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     n_corr.transpose() * Eelastic * m_corr - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * m_corr;
                TrialPlastic_Strain += dLambda_corr * m_corr;
            }
            else
            {
                // Do nothing
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                // cout << "Numeric error!\n";
                // printTensor1("TrialStress = " , TrialStress);
                // printTensor1("CommitStress = " , CommitStress);
                // printTensor1("depsilon = " , depsilon);
                // printTensor1("dsigma   = " , dsigma);
                // printTensor1("intersection_stress = " , intersection_stress);
                // printTensor2("Eelastic = " , Eelastic);
                // printTensor2("Stiffness = " , Stiffness);
                // cout << "yf_val_start = " << yf_val_start << endl;
                // cout << "yf_val_end = " << yf_val_end << endl;
                // printTensor1("n = " , n );
                // printTensor1("m = " , m );13
                // cout << "hardening  = " << hardening << endl;
                // cout << "den = " << den << endl;
                // cout << "dLambda = " << dLambda << endl;

                return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }
            else
            {
                // Ladruno (ADR-94 wp/94a)
                if (ladruno_strict_rejects("Forward_Euler_Subincrement (post return-to-yield)", TrialStress))
                    return LADRUNO_MATERIAL_REFUSED;
                return 0;
            }

            ComputeTangentStiffness();

        }

        return errorcode;
    }



    // int Backward_Euler(const VoigtVector & strain_incr)
    // {
    //     using namespace ASDPlasticMaterial3DGlobals;

    //     int errorcode = -1;

    //     static VoigtVector depsilon;
    //     depsilon *= 0;
    //     depsilon = strain_incr;

    //     const VoigtVector& sigma = CommitStress;
    //     const VoigtVector& epsilon = CommitStrain;

    //     iv_storage.revert_all();

    //     dsigma *= 0;
    //     intersection_stress *= 0;
    //     intersection_strain *= 0;

    //     VoigtMatrix Eelastic = et(sigma, parameters_storage);

    //     // Initial elastic predictor
    //     dsigma = Eelastic * depsilon;
    //     TrialStress = sigma + dsigma;
    //     TrialStrain = CommitStrain + depsilon;
    //     TrialPlastic_Strain = CommitPlastic_Strain;

    //     double yf_val_start = yf(sigma, iv_storage, parameters_storage);
    //     double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

    //     VoigtVector start_stress = CommitStress;
    //     VoigtVector end_stress = TrialStress;

    //     intersection_stress = start_stress;

    //     if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //     {
    //         Stiffness = Eelastic;
    //         return 0;
    //     }
    //     else  //Plasticity - Backward Euler Implementation
    //     {
    //         depsilon_elpl = depsilon;
            
    //         // Initialize backward Euler iteration variables
    //         VoigtVector sigma_trial = TrialStress;
    //         VoigtVector sigma_n_plus_1 = sigma_trial + Eelastic * depsilon_elpl; // Initial guess
    //         double dLambda = 0.0;
            
    //         // Newton-Raphson iteration for backward Euler
    //         int max_iterations = INT_OPT_n_max_iterations[ASDP_TAG];
    //         double f_tol = DBL_OPT_f_absolute_tol[ASDP_TAG];
    //         double stress_tol = DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
    //         for (int iter = 0; iter < max_iterations; iter++)
    //         {
    //             // Evaluate yield function at current stress state
    //             double f_val = yf(sigma_n_plus_1, iv_storage, parameters_storage);
                
    //             // Check convergence on yield function
    //             if (abs(f_val) <= f_tol)
    //             {
    //                 TrialStress = sigma_n_plus_1;
    //                 TrialPlastic_Strain = CommitPlastic_Strain + dLambda * pf(depsilon_elpl, sigma_n_plus_1, iv_storage, parameters_storage);
                    
    //                 // Check for NaN
    //                 double norm_trial_stress = TrialStress.transpose() * TrialStress;
    //                 if (norm_trial_stress != norm_trial_stress)
    //                 {
    //                     return -1;
    //                 }
                    
    //                 ComputeTangentStiffness();
    //                 return 0;
    //             }
                
    //             // Compute derivatives for Newton-Raphson
    //             const VoigtVector& n = yf.df_dsigma_ij(sigma_n_plus_1, iv_storage, parameters_storage);
    //             const VoigtVector& m = pf(depsilon_elpl, sigma_n_plus_1, iv_storage, parameters_storage);
    //             double hardening = yf.hardening(depsilon_elpl, m, sigma_n_plus_1, iv_storage, parameters_storage);
                
    //             // Build residual vector and Jacobian for backward Euler
    //             // Residual: R1 = sigma_n+1 - sigma_trial - E * (depsilon - dLambda * m)
    //             //          R2 = f(sigma_n+1, alpha_n+1)
                
    //             VoigtVector R1 = sigma_n_plus_1 - sigma_trial - Eelastic * (depsilon_elpl - dLambda * m);
    //             double R2 = f_val;
                
    //             // Check stress residual convergence
    //             double stress_residual_norm = R1.norm();
    //             if (stress_residual_norm <= stress_tol && abs(R2) <= f_tol)
    //             {
    //                 TrialStress = sigma_n_plus_1;
    //                 TrialPlastic_Strain = CommitPlastic_Strain + dLambda * m;
                    
    //                 // Check for NaN
    //                 double norm_trial_stress = TrialStress.transpose() * TrialStress;
    //                 if (norm_trial_stress != norm_trial_stress)
    //                 {
    //                     return -1;
    //                 }
                    
    //                 ComputeTangentStiffness();
    //                 return 0;
    //             }
                
    //             // Build Jacobian matrix for Newton-Raphson
    //             // J11 = I + dLambda * E * dm/dsigma
    //             // J12 = E * m
    //             // J21 = df/dsigma
    //             // J22 = -hardening
                
    //             VoigtMatrix I = VoigtMatrix::Identity();
    //             VoigtMatrix J11 = I; // Simplified - could add dLambda * E * dm/dsigma for better convergence
    //             VoigtVector J12 = Eelastic * m;
    //             VoigtVector J21 = n;
    //             double J22 = -hardening;
                
    //             // Solve Newton system using block elimination
    //             // [J11  J12] [Δσ   ]   [R1]
    //             // [J21  J22] [ΔdL  ] = [R2]
                
    //             // Eliminate to get: (J22 - J21 * J11^-1 * J12) * ΔdL = R2 - J21 * J11^-1 * R1
    //             double den = J22 - J21.transpose() * J12;
                
    //             if (abs(den) < MACHINE_EPSILON)
    //             {
    //                 cout << "Backward Euler - Singular Jacobian, den = " << den << endl;
    //                 return -1;
    //             }
                
    //             double delta_dLambda = (R2 - J21.transpose() * R1) / den;
    //             VoigtVector delta_sigma = -R1 - delta_dLambda * J12;
                
    //             // Update variables
    //             sigma_n_plus_1 += delta_sigma;
    //             dLambda += delta_dLambda;
                
    //             // Ensure plastic multiplier is non-negative
    //             if (dLambda < 0)
    //             {
    //                 dLambda = 0;
    //             }
                
    //             // Update internal variables based on current plastic multiplier
    //             iv_storage.apply([&m, &dLambda, &sigma_n_plus_1,  this](auto & iv)
    //             {
    //                 auto h = iv.hardening_function(depsilon_elpl, m, sigma_n_plus_1, parameters_storage);
    //                 iv.trial_value = iv.committed_value + dLambda * h;
    //             });
    //         }
            
    //         // If we reach here, Newton-Raphson did not converge
    //         cout << "Backward Euler - Newton-Raphson did not converge after " << max_iterations << " iterations" << endl;
    //         return -1;
    //     }

    //     return errorcode;
    // }


// int Backward_Euler(const VoigtVector & strain_incr)
// {
//     using namespace ASDPlasticMaterial3DGlobals;

//     int errorcode = -1;

//     static VoigtVector depsilon;
//     depsilon *= 0;
//     depsilon = strain_incr;

//     const VoigtVector& sigma   = CommitStress;
//     const VoigtVector& epsilon = CommitStrain;

//     iv_storage.revert_all();

//     dsigma *= 0;
//     intersection_stress *= 0;
//     intersection_strain *= 0;

//     // Elastic stiffness at committed state
//     VoigtMatrix Eelastic = et(sigma, parameters_storage);

//     // Elastic predictor
//     dsigma      = Eelastic * depsilon;
//     TrialStress = sigma + dsigma;
//     TrialStrain = CommitStrain + depsilon;
//     TrialPlastic_Strain = CommitPlastic_Strain;

//     // Simple elastic check at committed internal vars
//     const double f_tol     = DBL_OPT_f_absolute_tol[ASDP_TAG];
//     const double stress_tol= DBL_OPT_stress_absolute_tol[ASDP_TAG];
//     int max_iterations     = INT_OPT_n_max_iterations[ASDP_TAG];

//     double f_trial = yf(TrialStress, iv_storage, parameters_storage);

//     // Keep some of your existing bookkeeping
//     intersection_stress = CommitStress;

//     if (f_trial <= f_tol) {
//         // Elastic step
//         Stiffness = Eelastic;
//         return 0;
//     }

//     // --- Plastic corrector (Backward Euler / return mapping) ---
//     // Keep your depsilon_elpl usage to preserve pf()/hardening() signatures
//     depsilon_elpl = depsilon;

//     // Unknowns
//     VoigtVector sigma_n1 = TrialStress; // start from trial
//     double dLambda = 0.0;

//     for (int iter = 0; iter < max_iterations; ++iter)
//     {
//         // Gradients at (sigma_n1, committed internal vars)
//         const VoigtVector n = yf.df_dsigma_ij(sigma_n1, iv_storage, parameters_storage);               // ∂f/∂σ
//         const VoigtVector m = pf(depsilon_elpl, sigma_n1, iv_storage, parameters_storage);             // ∂g/∂σ (flow)
//         const double      H = yf.hardening(depsilon_elpl, m, sigma_n1, iv_storage, parameters_storage);// effective hardening slope (scalar)

//         // Residuals:
//         // R1 = σ_{n+1} - σ_trial + E * (dΛ m) = 0
//         // R2 = f(σ_{n+1}, α_n) + H dΛ = 0
//         const VoigtVector R1 = sigma_n1 - TrialStress + Eelastic * (dLambda * m);
//         const double      f_sigma = yf(sigma_n1, iv_storage, parameters_storage);
//         const double      R2 = f_sigma + H * dLambda;

//         // Convergence check
//         if (R1.norm() <= stress_tol && std::fabs(R2) <= f_tol) {
//             // Finalize state
//             const VoigtVector m_fin = pf(depsilon_elpl, sigma_n1, iv_storage, parameters_storage);
//             const VoigtVector n_fin = yf.df_dsigma_ij(sigma_n1, iv_storage, parameters_storage);
//             const double      H_fin = yf.hardening(depsilon_elpl, m_fin, sigma_n1, iv_storage, parameters_storage);

//             TrialStress = sigma_n1;
//             TrialStrain = CommitStrain + depsilon;
//             TrialPlastic_Strain = CommitPlastic_Strain + dLambda * m_fin;

//             // NaN guard (keep same style as your code)
//             double norm_trial_stress = TrialStress.transpose() * TrialStress;
//             if (norm_trial_stress != norm_trial_stress) {
//                 return -1;
//             }

//             // Consistent algorithmic tangent (non-associative allowed)
//             const VoigtVector Em = Eelastic * m_fin;
//             const VoigtVector En = Eelastic * n_fin;
//             double denom_tan = H_fin + n_fin.transpose() * Em;
//             if (std::fabs(denom_tan) < MACHINE_EPSILON) {
//                 Stiffness = Eelastic; // fallback
//             } else {
//                 // rank-1 update: C_ep = E - (E m) ⊗ (E n) / (H + n^T E m)
//                 Stiffness = Eelastic - (Em * En.transpose()) / denom_tan;
//             }

//             // Update internal vars to trial values: α_{n+1} = α_n + dΛ h
//             iv_storage.apply([&](auto & iv)
//             {
//                 auto h = iv.hardening_function(depsilon_elpl, m_fin, sigma_n1, parameters_storage);
//                 iv.trial_value = iv.committed_value + dLambda * h;
//             });

//             // Maintain dsigma for downstream use
//             dsigma = TrialStress - CommitStress;
//             return 0;
//         }

//         // Newton step with J11 = I (neglecting dm/dσ for robustness)
//         const VoigtVector Em = Eelastic * m;
//         const double denom = H + n.transpose() * Em;  // Schur complement denominator

//         if (std::fabs(denom) < MACHINE_EPSILON) {
//             // singular / ill-conditioned -> fail gracefully
//             return -1;
//         }

//         // ΔdΛ = (n^T R1 - R2) / (H + n^T E m)
//         const double delta_dLambda = (n.transpose() * R1 - R2) / denom;

//         // Project to dΛ >= 0
//         double d_dLambda = delta_dLambda;
//         if (dLambda + d_dLambda < 0.0) d_dLambda = -dLambda;

//         // Δσ = -(R1 + E m ΔdΛ)
//         const VoigtVector delta_sigma = -(R1 + Em * d_dLambda);

//         // Update unknowns
//         sigma_n1 += delta_sigma;
//         dLambda  += d_dLambda;

//         // Basic sanity
//         double norm_s = sigma_n1.transpose() * sigma_n1;
//         if (norm_s != norm_s) {
//             return -1;
//         }
//     }

//     // If we reach here, Newton failed to converge
//     std::cout << "Backward Euler - Newton-Raphson did not converge after " << max_iterations << " iterations\n";
//     return -1;
// }




    //==================================================================================================
    // Ladruno (ADR-97 wp/97b): CLOSEST-POINT return map (CPPM) + consistent tangent
    //==================================================================================================
    //
    // Solves, as ONE coupled Newton system in x = (sigma_{n+1}, q_{n+1}, dlambda):
    //
    //     R_sigma = s - s_tr + dl * E(s) * m(s, q)                     (6 rows)
    //     R_q     = q - q_n  - dl * h(s, q, m(s,q))                    (n_q rows)
    //     R_f     = f(s, q)                                            (1 row)
    //
    // with s_tr = sigma_n + E(sigma_n)*depsilon.  This is NOT what
    // `Backward_Euler` does: that is an Ortiz-Simo CUTTING PLANE which
    // re-evaluates n, m and H at the running iterate and accumulates
    // sigma -= deltaLambda*(E*m), so its fixed point is
    // s_tr - sum_k dl_k E m(sigma^k), and its internal-variable update is a
    // Newton-path quadrature exact only for constant h.  The two coincide exactly
    // when the flow direction does not rotate over the step (proportional loading)
    // and diverge otherwise -- measurably so for Armstrong-Frederick, which 22 of
    // the 46 registered specializations carry (ADR-94 M3; ADR-97 P0 measured a
    // 9.15 stress spread over four equally valid cutting-plane iterate paths on a
    // stress of ~46, against 1.1e-12 for this map).
    //
    // `Backward_Euler` is untouched (ADR-97 D1).  Every result the fork has
    // published on this material was obtained on it.
    //
    // Conventions (ADR-94 wp/94c, do not re-derive): tension positive,
    // p = trace/3, Voigt [11 22 33 12 23 13], strain-like slots ENGINEERING, every
    // YF/PF gradient taken w.r.t. the STORED slot, every contraction a plain dot.
    //
    // E(sigma_{n+1}) is evaluated INSIDE the residual (ADR-97 D6), which makes the
    // converged state hyperelastically consistent; the dE/dsigma JACOBIAN block is
    // optional (identically zero for LinearIsotropic3D_EL, the only elasticity
    // registered in all 43 non-StiffSoil specializations) -- dropping it costs
    // iterations, never accuracy, because the converged x still satisfies the
    // exact residual.
    //
    // There is NO intersection / elastic-fraction split, and none is missing:
    // `Backward_Euler` zeroes `intersection_stress`/`intersection_strain` at setup
    // and never reads either again (the yield-crossing bisection lives in
    // Forward_Euler and compute_local_stress only).
    //
    // No line search: for convex f with associated flow the CPPM residual is the
    // stationarity system of a strictly convex closest-point projection in the
    // E-metric, so Newton from the elastic predictor is globally convergent on the
    // yielding branch.  A line search that "succeeds" on a non-converged system is
    // ADR-94 M7's Backward_Euler_LineSearch failure mode (2/20 vs BE's 20/20).

    int cp_n_iv()
    {
        int n = 0;
        iv_storage.apply([&n](auto & internal_variable)
        {
            n += internal_variable.size();
        });
        return n;
    }

    // Scatter x into (TrialStress, IV trial values) and form the residual; when
    // `J` is non-null, form the Jacobian at the SAME point.  Returns false on a
    // non-finite residual.
    bool cp_assemble(const cp_vector_t& x, int N,
                     const VoigtVector& sigma_tr, const VoigtVector& depsilon,
                     VoigtVector& m_out, cp_vector_t& R, cp_matrix_t* J)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        const double dl = x(N - 1);
        for (int i = 0; i < 6; ++i) TrialStress(i) = x(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i) internal_variable.trial_value(i) = x(off + i);
                off += nq;
            });
        }

        // GCC: never bind an Eigen product (or a functor's mutable return buffer,
        // which the NEXT call overwrites) to a reference -- named value locals only.
        VoigtMatrix Ecur = et(TrialStress, parameters_storage);          // ADR-97 D6
        VoigtVector m    = pf(depsilon, TrialStress, iv_storage, parameters_storage);
        m_out = m;
        VoigtVector Em   = Ecur * m;
        const double f_val = yf(TrialStress, iv_storage, parameters_storage);

        for (int i = 0; i < 6; ++i)
            R(i) = x(i) - sigma_tr(i) + dl * Em(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                // hardening_function reads the TRIAL value, i.e. q_{n+1}: this is
                // what makes the internal-variable update implicit (ADR-97 D4).
                auto h = internal_variable.hardening_function(depsilon, m, TrialStress, parameters_storage);
                for (int i = 0; i < nq; ++i)
                    R(off + i) = x(off + i) - internal_variable.committed_value(i) - dl * h(i);
                off += nq;
            });
        }
        R(N - 1) = f_val;

        for (int i = 0; i < N; ++i)
            if (!(R(i) == R(i))) return false;      // NaN

        if (J == 0) return true;

        VoigtVector nvec = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
        VoigtMatrix dmds = pf.dm_dsigma(depsilon, TrialStress, iv_storage, parameters_storage);
        VoigtMatrix Edmds = Ecur * dmds;

        cp_matrix_t& Jm = *J;
        Jm.setZero();

        //   J_ss = I + dl * E * dm/ds     J_sl = E*m     J_fs = (df/ds)^T
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
                Jm(i, j) = ((i == j) ? 1.0 : 0.0) + dl * Edmds(i, j);
            Jm(i, N - 1) = Em(i);
            Jm(N - 1, i) = nvec(i);
        }
        //   + dl * (dE/ds : m)  -- ADR-97 D6; identically zero for LinearIsotropic3D_EL
        if constexpr (el_is_stress_dependent<ElasticityType>::value)
        {
            VoigtMatrix dEm = VoigtMatrix::Zero();
            et.dE_dsigma_contract(TrialStress, m, parameters_storage, dEm);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Jm(i, j) += dl * dEm(i, j);
        }

        //   J_sq = dl * E * dm/dq        J_fq = (df/dq)^T
        {
            int offb = 6;
            iv_storage.apply([&](auto & iv_b)
            {
                const int nb = iv_b.size();
                VoigtMatrix dmdq = VoigtMatrix::Zero();
                pf.dm_dq(iv_b, depsilon, TrialStress, iv_storage, parameters_storage, dmdq);
                VoigtMatrix Edmdq = Ecur * dmdq;
                double dfdq[6] = {0., 0., 0., 0., 0., 0.};
                yf.df_dq(iv_b, TrialStress, iv_storage, parameters_storage, dfdq);
                for (int c = 0; c < nb; ++c)
                {
                    for (int i = 0; i < 6; ++i) Jm(i, offb + c) = dl * Edmdq(i, c);
                    Jm(N - 1, offb + c) = dfdq[c];
                }
                offb += nb;
            });
        }

        //   J_ql = -h
        //   J_qs = -dl * (dh/dm)(dm/ds)                 [dh/ds is zero by construction:
        //                                                no policy in this tree reads sigma]
        //   J_qq = I - dl * (dh/dq + (dh/dm)(dm/dq))    [the OFF-DIAGONAL IV blocks are
        //                                                NOT zero: every h reads m, and m
        //                                                reads the flow direction's own
        //                                                back stress]
        {
            int offa = 6;
            iv_storage.apply([&](auto & iv_a)
            {
                const int na = iv_a.size();
                auto h = iv_a.hardening_function(depsilon, m, TrialStress, parameters_storage);
                VoigtMatrix dhdm = VoigtMatrix::Zero();
                VoigtMatrix dhdq = VoigtMatrix::Zero();
                iv_a.hardening_dh_dm(depsilon, m, TrialStress, parameters_storage, dhdm);
                iv_a.hardening_dh_dq(depsilon, m, TrialStress, parameters_storage, dhdq);
                for (int rr = 0; rr < na; ++rr)
                {
                    Jm(offa + rr, N - 1) = -h(rr);
                    for (int j = 0; j < 6; ++j)
                    {
                        double v = 0.0;
                        for (int k = 0; k < 6; ++k) v += dhdm(rr, k) * dmds(k, j);
                        Jm(offa + rr, j) = -dl * v;
                    }
                }
                int offb = 6;
                iv_storage.apply([&](auto & iv_b)
                {
                    const int nb = iv_b.size();
                    VoigtMatrix dmdq_b = VoigtMatrix::Zero();
                    pf.dm_dq(iv_b, depsilon, TrialStress, iv_storage, parameters_storage, dmdq_b);
                    const bool same_iv = (offa == offb);
                    for (int rr = 0; rr < na; ++rr)
                        for (int c = 0; c < nb; ++c)
                        {
                            double v = same_iv ? dhdq(rr, c) : 0.0;
                            for (int k = 0; k < 6; ++k) v += dhdm(rr, k) * dmdq_b(k, c);
                            Jm(offa + rr, offb + c) =
                                ((same_iv && rr == c) ? 1.0 : 0.0) - dl * v;
                        }
                    offb += nb;
                });
                offa += na;
            });
        }
        return true;
    }

    // ELASTIC-METRIC apex classification.  `check_apex_region` in the yield
    // function is EUCLIDEAN (p - p_apex >= eta*q) and says so in its own comment:
    // the exact condition needs K, G and the dilatancy, none of which the YF's
    // signature can see.  `Closest_Point` classifies HERE, where E is in scope, and
    // does not call `check_apex_region` at all.  ADR-97's P0 oracle quantified the
    // cost of the Euclidean test in BOTH directions: at etabar = 0.2 (exact slope
    // 0.333 < the header's 0.4) a trial at (p-p_apex)/q = 0.36 is classified CONE
    // and the cone return then gives sqrt(J2)_{n+1} = -0.47, an INADMISSIBLE
    // negative deviatoric norm; at etabar = eta = 0.4 (exact slope 0.667) trials at
    // 0.45 and 0.60 are classified APEX although the correct return is to the cone.
    //
    // The test itself is family-agnostic: take the LINEARISED cone step at the
    // trial state, dl0 = f_tr/(n:E:m - H), and ask whether the deviatoric part of
    // the returned stress has FLIPPED sign.  For Drucker-Prager E*m has the
    // deviatoric part (G/q)*r exactly, so this reduces to q_tr - G*dl0 < 0 -- the
    // oracle's exact test, with K*etabar and G entering through E.
    bool cp_apex_region(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double f_tr, double tol_f)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        VoigtVector m0 = pf(depsilon, sigma_tr, iv_storage, parameters_storage);
        VoigtVector n0 = yf.df_dsigma_ij(sigma_tr, iv_storage, parameters_storage);
        VoigtVector Em0 = Eelastic * m0;
        // The plastic modulus formed from CLOSEST_POINT's own df/dq, not from
        // yf.hardening(): for Drucker-Prager the two differ, because the shipped
        // yf_hardening carries a df/dk = -1 term for a cohesion IV that its own f
        // does not contain (ADR-97 P0 header finding 2).
        double H0 = 0.0;
        iv_storage.apply([&](auto & internal_variable)
        {
            double d[6] = {0., 0., 0., 0., 0., 0.};
            yf.df_dq(internal_variable, sigma_tr, iv_storage, parameters_storage, d);
            auto h = internal_variable.hardening_function(depsilon, m0, sigma_tr, parameters_storage);
            const int nq = internal_variable.size();
            for (int i = 0; i < nq; ++i) H0 += d[i] * h(i);
        });
        VoigtVector dev_tr = sigma_tr.deviator();
        double q2_tr = tensor_dot_stress_like(dev_tr, dev_tr);
        if (!(q2_tr > 0.0)) q2_tr = 0.0;
        const double q_tr = std::sqrt(q2_tr);
        // A trial state ON the hydrostatic axis and outside the surface can only
        // return to the vertex.  The flip test below is a SIGN test on
        // dot(dev_ret, dev_tr), which degenerates to 0 < 0 -- i.e. says CONE --
        // when the trial deviator vanishes; the cone Newton then has no flow
        // direction at all and exhausts its iterations.  That degenerate state
        // is not exotic: it is exactly ADR-94 B4's hydrostatic-tension
        // reproducer, the one that used to commit NaN.  `tol_f` is the yield
        // tolerance, so this comparison is in stress units and unit consistent
        // (ADR-94 M5).
        if (q_tr <= tol_f) return true;
        const double den = n0.dot(Em0) - H0;
        if (!(den > MACHINE_EPSILON)) return false;
        const double dl0 = f_tr / den;
        VoigtVector dev_Em0 = Em0.deviator();
        VoigtVector dev_ret = dev_tr - dl0 * dev_Em0;
        return tensor_dot_stress_like(dev_ret, dev_tr) < 0.0;
    }

    // Reduced apex system.  At the vertex dm/ds blows up as 1/sqrt(J2) and the
    // plastic flow direction is a whole subdifferential cone, so the 6-row
    // residual cannot produce it: sigma_{n+1} is PINNED at the yield function's
    // own apex, and the only unknowns are the internal variables and dlambda.
    //     sigma_{n+1} = sigma_apex(q_{n+1})
    //     d eps^p     = E^{-1}(sigma_tr - sigma_apex)
    //     dl          = <m_apex, d eps^p>_e / <m_apex, m_apex>_e   (>= 0)
    //     q_{n+1}     = q_n + dl * h(sigma_apex, q_{n+1}, m_apex)
    // solved by fixed point (ONE pass for perfect plasticity, where h == 0).
    // Returns false when the yield function's own apex is not on its own surface,
    // in which case the caller falls through to the generic Newton -- the same
    // admissibility gate `Backward_Euler` uses, kept verbatim.
    bool cp_apex_return(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double tol_f, double tol_s,
                        int max_iter, int& rc)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        Eigen::Matrix<double, 6, 6> E_eig;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) E_eig(i, j) = Eelastic(i, j);
        Eigen::FullPivLU< Eigen::Matrix<double, 6, 6> > elu(E_eig);
        if (!elu.isInvertible())
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the elastic tangent is singular at the apex return"
                   << " -- rejecting step" << endln;
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }

        VoigtVector sigma_apex = VoigtVector::Zero();
        VoigtVector m_apex     = VoigtVector::Zero();
        VoigtVector dep        = VoigtVector::Zero();
        double dl = 0.0;
        bool ok = false;

        for (int it = 0; it < max_iter; ++it)
        {
            // Copy out of the YF's return buffer immediately: apex_stress and
            // df_dsigma_ij share one mutable member per functor.
            VoigtVector sa = yf.apex_stress(iv_storage, parameters_storage);
            sigma_apex = sa;
            const double f_apex = yf(sigma_apex, iv_storage, parameters_storage);
            const double f_apex_abs = (f_apex < 0) ? -f_apex : f_apex;
            if (!(f_apex_abs <= tol_f))
            {
                // The apex this YF names is not on its own surface: fall through
                // to the generic map rather than commit a fabricated stress.
                iv_storage.revert_all();
                return false;
            }
            VoigtVector ma = pf(depsilon, sigma_apex, iv_storage, parameters_storage);
            m_apex = ma;

            Eigen::Matrix<double, 6, 1> rhs_e;
            for (int i = 0; i < 6; ++i) rhs_e(i) = sigma_tr(i) - sigma_apex(i);
            Eigen::Matrix<double, 6, 1> dep_e = elu.solve(rhs_e);
            for (int i = 0; i < 6; ++i) dep(i) = dep_e(i);

            // Both operands are ENGINEERING-strain-like Voigt vectors, so the
            // projection uses the engineering contraction (ADR-94 wp/94c B5).
            const double mm = tensor_dot_engineering_strain_like(m_apex, m_apex);
            double dl_new = 0.0;
            if (mm > MACHINE_EPSILON)
            {
                dl_new = tensor_dot_engineering_strain_like(m_apex, dep) / mm;
                if (dl_new < 0.0) dl_new = 0.0;
            }

            double dq_max = 0.0;
            iv_storage.apply([&](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon, m_apex, sigma_apex, parameters_storage);
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i)
                {
                    const double nv = internal_variable.committed_value(i) + dl_new * h(i);
                    const double d = nv - internal_variable.trial_value(i);
                    const double ad = (d < 0) ? -d : d;
                    if (ad > dq_max) dq_max = ad;
                    internal_variable.trial_value(i) = nv;
                }
            });
            dl = dl_new;
            cp_last_iterations = it + 1;
            if (dq_max <= tol_s) { ok = true; break; }
        }

        if (!ok)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the apex internal-variable fixed point did not converge in "
                   << max_iter << " iterations -- rejecting step" << endln;
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }

        {
            VoigtVector sa = yf.apex_stress(iv_storage, parameters_storage);
            sigma_apex = sa;
            Eigen::Matrix<double, 6, 1> rhs_e;
            for (int i = 0; i < 6; ++i) rhs_e(i) = sigma_tr(i) - sigma_apex(i);
            Eigen::Matrix<double, 6, 1> dep_e = elu.solve(rhs_e);
            for (int i = 0; i < 6; ++i) dep(i) = dep_e(i);
        }
        TrialStress         = sigma_apex;
        TrialPlastic_Strain = CommitPlastic_Strain + dep;
        (void) dl;

        // Algorithmic apex tangent.  sigma_{n+1} = sigma_apex(q_{n+1}), so
        // C_alg = (d sigma_apex/dq)(dq/d eps).  For every yield function this ADR
        // ships, apex_stress() is built from MODEL PARAMETERS only
        // (DruckerPrager_YF: p_apex = xi_c/eta), so d sigma_apex/dq is identically
        // zero and the consistent tangent is EXACTLY the zero matrix -- rank 0,
        // which is what the P0 oracle pins (cppm_dp.py, "apex, associated,
        // perfect": rank C = 0, FD error 0).  That is rank deficient by
        // construction and WILL make an element whose every Gauss point sits at
        // the apex singular, which is the true state of affairs.  The finite
        // difference below MEASURES d sigma_apex/dq instead of assuming it; a
        // yield function whose apex moves with an internal variable gets a
        // one-time warning and the same zero, which ADR-97 P5 replaces.
        VoigtMatrix apex_stiff = VoigtMatrix::Zero();
        {
            double a_max = 0.0;
            VoigtVector s0 = yf.apex_stress(iv_storage, parameters_storage);
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i)
                {
                    const double v0 = internal_variable.trial_value(i);
                    const double av = (v0 < 0) ? -v0 : v0;
                    const double dv = 1e-8 * (av + 1.0);
                    internal_variable.trial_value(i) = v0 + dv;
                    VoigtVector sp = yf.apex_stress(iv_storage, parameters_storage);
                    internal_variable.trial_value(i) = v0;
                    for (int r = 0; r < 6; ++r)
                    {
                        const double d = (sp(r) - s0(r)) / dv;
                        const double ad = (d < 0) ? -d : d;
                        if (ad > a_max) a_max = ad;
                    }
                }
            });
            if (a_max > 0.0
                && INT_OPT_tangent_operator_type[ASDP_TAG]
                   == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
            {
                static bool warned_apex_tangent = false;
                if (!warned_apex_tangent)
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - this yield function's apex MOVES with an internal"
                           << " variable (|d sigma_apex/dq| = " << a_max << "), so the"
                           << " zero apex tangent reported under tangent_type"
                           << " Algorithmic is approximate (ADR-97 P5)." << endln;
                    warned_apex_tangent = true;
                }
            }
        }

        {
            using TOT = ASDPlasticMaterial3D_Tangent_Operator_Type;
            switch (INT_OPT_tangent_operator_type[ASDP_TAG])
            {
            case TOT::Elastic:
                Stiffness = Eelastic;
                break;
            case TOT::Continuum:
            case TOT::Algorithmic:
            case TOT::Numerical_Algorithmic_FirstOrder:
            case TOT::Numerical_Algorithmic_SecondOrder:
                Stiffness = apex_stiff;
                break;
            case TOT::Secant:
            default:
                Stiffness = VoigtMatrix((apex_stiff + Eelastic) / 2.0);
                break;
            }
        }

        if (ladruno_strict_rejects("Closest_Point (apex)", TrialStress))
        {
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }
        rc = 0;
        return true;
    }


    //==================================================================================================
    // Ladruno (ADR-97 wp/97c): PRINCIPAL-STRESS-SPACE multi-surface closest-point
    // return for the Mohr-Coulomb family (Clausen, Damkilde & Andersen 2006/2007)
    //==================================================================================================
    //
    // WHY A SECOND ALGORITHM (ADR-97 D3).  The shipped 6D Mohr-Coulomb yield
    // function is the exact Lode-angle form
    //     f = A(theta) sqrt(J2) + I1 sin(phi)/3 - c cos(phi),
    //     A(theta) = cos(theta) - sin(theta) sin(phi)/sqrt(3),
    // whose gradient carries a 1/cos(3 theta) that BLOWS UP on the corners; the
    // header dodges it by swapping in a Drucker-Prager gradient for
    // |theta| >= 29 deg and otherwise central-differencing f over the six raw
    // Voigt slots.  A coupled Newton built on that gradient cannot converge
    // quadratically at a corner -- which is precisely where Mohr-Coulomb models
    // spend their time.  In PRINCIPAL STRESS SPACE the same surface is six
    // PLANES: on the sorted sextant s1 >= s2 >= s3 (tension positive) it is the
    // single plane
    //     f = a_13 . s - k ,  a_ij = [ (1+sin phi)/2, 0, -(1-sin phi)/2 ], k = c cos(phi),
    // its two corners are LINES and its vertex is the point (k/sin phi)*[1,1,1].
    // Every return is then a LINEAR projection in the elastic metric -- closed
    // form, no Newton at all -- and the tangent is constant on each region.
    // The P0 oracle `adr97_oracle/cppm_mc.py` verifies numerically, over 2000
    // random states, that this principal-stress f IS the header's own invariant
    // expression to 1e-14 relative; every number below is pinned against it.
    //
    // FLOW DIRECTION.  The header's `m` is deviator(dg/dsigma with PHI) plus
    // sin(psi)/3 * delta -- the deviatoric shape of the phi-surface with a
    // psi-controlled volumetric part, NOT the textbook non-associated gradient.
    // In principal space that is m_ij = a_ij - (sin phi)/3 * 1 + (sin psi)/3 * 1,
    // which collapses to a_ij when psi == phi.  Mirrored here, not corrected.
    //
    // REGION SELECTION IS BY BOUNDARY PLANES, NEVER BY AN ACTIVE-SET SEARCH.
    // Each edge line L has a boundary plane spanned by the line direction ell_L
    // and the face return direction rp = D m_13, passing through the apex; the
    // side of it the trial principal point falls on decides face vs edge, and the
    // sign of the edge-line parameter of the edge return decides edge vs apex.
    // This is exact and branch-free, and it REPLACES the |theta| < 29 deg guard.
    // It is also not optional: at the apex with psi < phi the three Koiter
    // multipliers are NOT all positive (the cone of return directions no longer
    // contains the hydrostatic direction), so a "grow the active set while
    // dLambda >= 0" search classifies those states wrongly -- ADR-97 P0 header
    // finding 6, and the reason the oracle prints the multipliers with that note.
    //
    // The boundary-plane SIGNS are taken analytically here, not from the oracle's
    // dimensional reference point: any point strictly inside the open face is
    // apex + b1*ell_1 + b2*ell_2 with b1, b2 > 0, and n_1 . ell_1 == 0 by
    // construction, so sgn_1 = sign(n_1 . ell_2) and symmetrically.  Verified
    // identical to the oracle's calibration over 60 (phi, psi, c) combinations
    // and 13094 trial states, mismatches 0 (see the P2 report).
    //
    // BACK-TRANSFORM.  sigma_ret is an isotropic tensor function of sigma_tr, so
    //     d sigma / d sigma_tr = Rs * T * Rs^{-1}
    // with Rs the Voigt image of E -> Q E Q^T (Q the TRIAL eigenvectors),
    // T's normal block = dy/dx (the principal 3x3 Koiter block) and T's shear
    // slots = (y_i - y_j)/(x_i - x_j) -- the eigenprojection ROTATION term.  That
    // ratio is 0/0 on every axisymmetric path (most triaxial decks, so this is
    // the common case, not the corner case); the l'Hopital limit is
    // dy_i/dx_i - dy_i/dx_j.  The switch threshold is RELATIVE to the yield
    // function's strength scale (ADR-94 M5: an absolute one reintroduces the
    // unit dependence that decided pass/fail in kPa vs Pa).  Rs^{-1} is built
    // directly as the Voigt image of E -> Q^T E Q, not inverted numerically.
    //
    // NO NEWTON, so `cp_iterations` reads 1 on a plastic step here and 0 on an
    // elastic one -- it is a region count, not an iteration count.
    //
    // `Backward_Euler` is untouched (ADR-97 D1): it still runs the scalar Newton
    // on the smoothed 6D gradient, and for MohrCoulombTensionCutoff still calls
    // the same ADR-84 `special_return` hook this map calls.

    // Apply the material's configured tangent_type to a raw active-set (Koiter)
    // tangent.  Mirrors the policy ADR-84 P3 set for `special_return` and the one
    // P1 set for the Drucker-Prager apex: `Continuum` maps to the same raw
    // operator, because the continuum elastoplastic operator of a multi-surface
    // corner is not defined by a single n : E : m.
    void cp_apply_tangent_policy(const VoigtMatrix& raw, const VoigtMatrix& Eelastic)
    {
        using TOT = ASDPlasticMaterial3D_Tangent_Operator_Type;
        switch (INT_OPT_tangent_operator_type[ASDP_TAG])
        {
        case TOT::Elastic:
            Stiffness = Eelastic;
            break;
        case TOT::Continuum:
        case TOT::Algorithmic:
        case TOT::Numerical_Algorithmic_FirstOrder:
        case TOT::Numerical_Algorithmic_SecondOrder:
            Stiffness = raw;
            break;
        case TOT::Secant:
        default:
            Stiffness = VoigtMatrix((raw + Eelastic) / 2.0);
            break;
        }
    }

    int cp_principal_return(const VoigtVector& depsilon,
                            const VoigtVector& sigma_tr,
                            const VoigtMatrix& Eelastic,
                            double tol_f)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        (void) depsilon;

        // ---- surface constants ------------------------------------------
        double sin_phi = 0.0, k_coh = 0.0, sin_phi_pf = 0.0, sin_psi = 0.0;
        if (!yf.cp_mc_face_params(iv_storage, parameters_storage, sin_phi, k_coh)
                || !pf.cp_mc_flow_params(iv_storage, parameters_storage,
                                         sin_phi_pf, sin_psi))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - this yield function / plastic flow pair is marked as a"
                   << " principal-space Mohr-Coulomb family but does not supply the"
                   << " face parameters (ADR-97 P2) -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        {
            const double d = sin_phi - sin_phi_pf;
            if ((d > 1e-12) || (d < -1e-12))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the yield function and the plastic flow direction"
                       << " disagree on MC_phi (" << sin_phi << " vs " << sin_phi_pf
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        if (!(sin_phi > 100 * MACHINE_EPSILON))
        {
            // phi == 0 is Tresca: the Mohr-Coulomb apex is at infinity, and
            // Clausen's boundary planes are defined THROUGH it.  Rather than
            // silently reclassify, refuse and name the alternative.
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - integration_method Closest_Point needs MC_phi > 0"
                   << " (ADR-97 P2): at phi = 0 the Mohr-Coulomb surface is a"
                   << " Tresca prism whose apex is at infinity, so the"
                   << " boundary-plane region test is undefined."
                   << " Use Backward_Euler." << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        // ---- isotropy of the elastic tangent ----------------------------
        // The principal-space projection needs D in the principal frame, which
        // exists only for an isotropic E.  Same guard ADR-84's special_return
        // uses; LinearIsotropic3D_EL is the only elasticity registered with this
        // family, so this is a future-proofing refusal, not a live branch.
        const double lam = Eelastic(0, 1);
        const double Gmu = (Eelastic(0, 0) - Eelastic(0, 1)) / 2.0;
        {
            const double sE = ((Eelastic(0, 0) < 0) ? -Eelastic(0, 0) : Eelastic(0, 0))
                              + ((Gmu < 0) ? -Gmu : Gmu);
            const double t = 1e-8 * sE;
            const double d02 = Eelastic(0, 2) - lam;
            const double d12 = Eelastic(1, 2) - lam;
            const double d33 = Eelastic(3, 3) - Gmu;
            const double d11 = Eelastic(1, 1) - Eelastic(0, 0);
            if (!(Gmu > 0.0)
                    || (d02 > t) || (d02 < -t) || (d12 > t) || (d12 < -t)
                    || (d33 > t) || (d33 < -t) || (d11 > t) || (d11 < -t))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space Mohr-Coulomb return needs an"
                       << " ISOTROPIC elastic tangent -- rejecting step (ADR-97 P2)"
                       << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        Eigen::Matrix3d D3;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                D3(i, j) = lam + ((i == j) ? 2.0 * Gmu : 0.0);

        // ---- spectral decomposition of the TRIAL stress, DESCENDING ------
        Eigen::Matrix3d st;
        st << sigma_tr.v11(), sigma_tr.v12(), sigma_tr.v13(),
              sigma_tr.v12(), sigma_tr.v22(), sigma_tr.v23(),
              sigma_tr.v13(), sigma_tr.v23(), sigma_tr.v33();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(st);
        if (es.info() != Eigen::Success)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the trial stress eigen-decomposition failed"
                   << " -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        Eigen::Vector3d xv;
        Eigen::Matrix3d Q;
        for (int i = 0; i < 3; ++i)          // eigenvalues come out ASCENDING
        {
            xv(i) = es.eigenvalues()(2 - i);
            Q.col(i) = es.eigenvectors().col(2 - i);
        }

        // ---- the three surfaces of the sorted sextant --------------------
        // index 0 = (1,3), 1 = (2,3), 2 = (1,2), in 0-based principal slots.
        static const int PIJ[3][2] = {{0, 2}, {1, 2}, {0, 1}};
        const double Aco = 0.5 * (1.0 + sin_phi);
        const double Bco = 0.5 * (1.0 - sin_phi);
        Eigen::Vector3d av[3], mv[3];
        for (int t = 0; t < 3; ++t)
        {
            av[t].setZero();
            av[t](PIJ[t][0]) =  Aco;
            av[t](PIJ[t][1]) = -Bco;
            mv[t] = av[t];
            const double vol = (sin_psi - sin_phi) / 3.0;
            mv[t](0) += vol; mv[t](1) += vol; mv[t](2) += vol;
        }

        const Eigen::Vector3d apexv = Eigen::Vector3d::Constant(k_coh / sin_phi);
        Eigen::Vector3d Dm13 = D3 * mv[0];

        // edge directions: null(A) of the two active gradients = their cross
        // product, oriented (like the oracle) so ell(0) >= ell(2).
        Eigen::Vector3d ell1 = av[0].cross(av[1]);      // edge s1 == s2
        if (ell1(0) < ell1(2)) ell1 = -ell1;
        ell1.normalize();
        Eigen::Vector3d ell2 = av[2].cross(av[0]);      // edge s2 == s3
        if (ell2(0) < ell2(2)) ell2 = -ell2;
        ell2.normalize();
        Eigen::Vector3d nb1 = ell1.cross(Dm13);
        Eigen::Vector3d nb2 = ell2.cross(Dm13);
        const double dot12 = nb1.dot(ell2);
        const double dot21 = nb2.dot(ell1);
        if (!(dot12 > 0.0 || dot12 < 0.0) || !(dot21 > 0.0 || dot21 < 0.0))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - degenerate Mohr-Coulomb boundary-plane geometry"
                   << " -- rejecting step (ADR-97 P2)" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        const double sgn1 = (dot12 > 0.0) ? 1.0 : -1.0;
        const double sgn2 = (dot21 > 0.0) ? 1.0 : -1.0;

        // ---- region classification --------------------------------------
        const Eigen::Vector3d dtr = xv - apexv;
        const double pl1 = sgn1 * nb1.dot(dtr);
        const double pl2 = sgn2 * nb2.dot(dtr);

        Eigen::Vector3d yv;
        Eigen::Matrix3d dydx = Eigen::Matrix3d::Zero();
        const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
        int region = 0;                      // 0 face, 1 edge s1==s2, 2 edge s2==s3, 3 apex

        if (pl1 >= 0.0 && pl2 >= 0.0)
        {
            const double den = av[0].dot(Dm13);
            if (!(den > MACHINE_EPSILON))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular Mohr-Coulomb face return (a . E m = " << den
                       << ") -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            const double dl = (av[0].dot(xv) - k_coh) / den;
            if (dl < 0.0)
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the Mohr-Coulomb face return produced a NEGATIVE"
                       << " plastic multiplier (" << dl << ") from a trial state"
                       << " outside the surface -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            yv = xv - dl * Dm13;
            Eigen::Matrix3d outer = Dm13 * av[0].transpose();
            dydx = I3 - outer / den;
            region = 0;
        }
        else
        {
            const bool use1 = (pl1 < pl2);
            const int t0 = use1 ? 0 : 2;     // LINE1 = (13,23); LINE2 = (12,13)
            const int t1 = use1 ? 1 : 0;
            Eigen::Matrix<double, 2, 3> Am;
            Am.row(0) = av[t0].transpose();
            Am.row(1) = av[t1].transpose();
            Eigen::Matrix<double, 3, 2> Mm;
            Mm.col(0) = mv[t0];
            Mm.col(1) = mv[t1];
            Eigen::Matrix<double, 3, 2> DM = D3 * Mm;
            Eigen::Matrix2d ADM = Am * DM;
            const double det = ADM(0, 0) * ADM(1, 1) - ADM(0, 1) * ADM(1, 0);
            const double dscale = (ADM.cwiseAbs().maxCoeff() + MACHINE_EPSILON);
            const double adet = (det < 0) ? -det : det;
            if (!(adet > MACHINE_EPSILON * dscale * dscale))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular Mohr-Coulomb edge return (det A E M = "
                       << det << ") -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            Eigen::Matrix2d ADMi;
            ADMi(0, 0) =  ADM(1, 1) / det;  ADMi(0, 1) = -ADM(0, 1) / det;
            ADMi(1, 0) = -ADM(1, 0) / det;  ADMi(1, 1) =  ADM(0, 0) / det;
            Eigen::Vector2d rhs;
            rhs(0) = av[t0].dot(xv) - k_coh;
            rhs(1) = av[t1].dot(xv) - k_coh;
            Eigen::Vector2d dl2 = ADMi * rhs;
            Eigen::Vector3d yline = xv - DM * dl2;
            const double tpar = (yline - apexv).dot(use1 ? ell1 : ell2);
            if (tpar < 0.0)
            {
                // Beyond the apex along the edge: the vertex is the return.  The
                // three Koiter multipliers there need NOT be positive for psi <
                // phi, which is exactly why this is decided geometrically.
                yv = apexv;
                dydx.setZero();
                region = 3;
            }
            else
            {
                yv = yline;
                Eigen::Matrix<double, 3, 2> DMi = DM * ADMi;
                Eigen::Matrix3d proj = DMi * Am;
                dydx = I3 - proj;
                region = use1 ? 1 : 2;
            }
        }

        // ---- back to 6D --------------------------------------------------
        Eigen::Matrix3d Sret = Q * yv.asDiagonal() * Q.transpose();
        VoigtVector sigma_ret(Sret(0, 0), Sret(1, 1), Sret(2, 2),
                              Sret(0, 1), Sret(1, 2), Sret(0, 2));
        for (int i = 0; i < 6; ++i)
        {
            if (!(sigma_ret(i) == sigma_ret(i)))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - NaN in the principal-space Mohr-Coulomb return"
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // Admissibility, in TWO places, for two different reasons.
        //
        // (a) EXACTLY, in principal space, against ALL THREE surfaces of the
        //     sextant.  This is where the return was computed and it is perfectly
        //     conditioned, so it is the check that catches a coding error -- a
        //     misclassified region, or an edge return that overshoots past a
        //     third surface.
        //
        // (b) LOOSELY, against the header's own COMPOSITE f.  For plain
        //     Mohr-Coulomb (b) is redundant with (a); for
        //     MohrCoulombTensionCutoff it is the real gate, because this path is
        //     only reached after `special_return` declined and the plain-MC
        //     return it then performs must still respect the cutoff plane.  Its
        //     tolerance is RELATIVE to the returned stress magnitude, NOT the
        //     bare `yf_tolerance()`: that accessor defaults to the ABSOLUTE
        //     `f_absolute_tol = 1e-6` (`f_relative_tol` defaults to 0), and on
        //     the ADR-84 MCTC deck (kPa, |sigma| ~ 5.4e3, strength scale 94) the
        //     header's f recomputed from the reassembled `Q diag(y) Q^T` is
        //     3.4e-6 -- 6e-10 relative, i.e. round-off, amplified because an EDGE
        //     return lands exactly on a corner where the Lode angle is ill
        //     conditioned (dtheta/dJ3 ~ 1/cos(3 theta)).  Refusing on that would
        //     be ADR-94 M5 again: the same model in Pa and in kPa disagreeing.
        //     A genuine fall-through error is O(|sigma|), which 1e-8 relative
        //     still catches by eight orders of magnitude.
        double sig_max = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            const double a = (sigma_ret(i) < 0) ? -sigma_ret(i) : sigma_ret(i);
            if (a > sig_max) sig_max = a;
        }
        double scale_ref = yf.strength_scale(iv_storage, parameters_storage);
        if (scale_ref < 0) scale_ref = -scale_ref;
        const double ref_mag = (sig_max > scale_ref) ? sig_max : scale_ref;
        {
            double f_princ = -1e300;
            for (int t = 0; t < 3; ++t)
            {
                const double v = av[t].dot(yv) - k_coh;
                if (v > f_princ) f_princ = v;
            }
            if (!(f_princ <= 1e-10 * ref_mag))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE the"
                       << " Mohr-Coulomb cone (max_k a_k . y - k = " << f_princ
                       << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        {
            const double f_ret = yf(sigma_ret, iv_storage, parameters_storage);
            const double tol_adm = (tol_f > 1e-8 * ref_mag) ? tol_f : 1e-8 * ref_mag;
            if (!(f_ret <= tol_adm))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE this yield"
                       << " function's own surface (f = " << f_ret << " > tol = "
                       << tol_adm << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // ---- the Koiter tangent, rotated to the global Voigt frame -------
        VoigtMatrix Tp = VoigtMatrix::Zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Tp(i, j) = dydx(i, j);
        {
            double xmax = 0.0;
            for (int i = 0; i < 3; ++i)
            {
                const double a = (xv(i) < 0) ? -xv(i) : xv(i);
                if (a > xmax) xmax = a;
            }
            const double eps_deg = 1e-9 * ((xmax > scale_ref) ? xmax : scale_ref);
            static const int SIJ[3][2] = {{0, 1}, {1, 2}, {0, 2}};   // slots 3,4,5
            for (int s = 0; s < 3; ++s)
            {
                const int i = SIJ[s][0], j = SIJ[s][1];
                const double dx = xv(i) - xv(j);
                const double adx = (dx < 0) ? -dx : dx;
                Tp(3 + s, 3 + s) = (adx > eps_deg)
                                   ? (yv(i) - yv(j)) / dx
                                   : (dydx(i, i) - dydx(i, j));
            }
        }

        VoigtMatrix Rs = VoigtMatrix::Zero();
        VoigtMatrix Rsi = VoigtMatrix::Zero();
        {
            static const int BI[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}};
            for (int k = 0; k < 6; ++k)
            {
                Eigen::Matrix3d Ek = Eigen::Matrix3d::Zero();
                Ek(BI[k][0], BI[k][1]) = 1.0;
                Ek(BI[k][1], BI[k][0]) = 1.0;
                Eigen::Matrix3d F  = Q * Ek * Q.transpose();
                Eigen::Matrix3d Fi = Q.transpose() * Ek * Q;
                for (int r = 0; r < 6; ++r)
                {
                    Rs(r, k)  = F(BI[r][0], BI[r][1]);
                    Rsi(r, k) = Fi(BI[r][0], BI[r][1]);
                }
            }
        }
        VoigtMatrix RT   = Rs * Tp;
        VoigtMatrix RTR  = RT * Rsi;
        VoigtMatrix Calg = RTR * Eelastic;

        // ---- plastic strain increment ------------------------------------
        // d eps^p = E^{-1} (sigma_tr - sigma_ret): convention-safe (it never
        // touches the engineering-shear factor of the principal flow vectors)
        // and exact for every region including the vertex.
        VoigtVector dep = VoigtVector::Zero();
        {
            Eigen::Matrix<double, 6, 6> Ee;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Ee(i, j) = Eelastic(i, j);
            Eigen::FullPivLU< Eigen::Matrix<double, 6, 6> > elu(Ee);
            if (!elu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the elastic tangent is singular"
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            Eigen::Matrix<double, 6, 1> rhs6;
            for (int i = 0; i < 6; ++i) rhs6(i) = sigma_tr(i) - sigma_ret(i);
            Eigen::Matrix<double, 6, 1> d6 = elu.solve(rhs6);
            for (int i = 0; i < 6; ++i) dep(i) = d6(i);
        }

        TrialStress         = sigma_ret;
        TrialPlastic_Strain = CommitPlastic_Strain + dep;
        cp_last_iterations  = 1;      // closed form: no Newton on any region
        cp_apply_tangent_policy(Calg, Eelastic);

        if (ladruno_strict_rejects("Closest_Point (principal)", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;
        return 0;
    }

    //==================================================================================================
    // Ladruno (ADR-97 wp/97d): PRINCIPAL-STRESS-SPACE closest-point return for the
    // HOEK-BROWN family (Clausen & Damkilde 2008) -- a CURVED surface
    //==================================================================================================
    //
    // WHY A THIRD ALGORITHM.  P1's smooth 6D Newton cannot be used: `HoekBrown_YF`
    // has NO analytic gradient at all -- `YIELD_FUNCTION_STRESS_DERIVATIVE` central
    // differences the COMPOSITE max(f_shear, f_tension) over the six raw Voigt
    // slots, which picks up spurious gradient across the branch switch (ADR-94 M5)
    // and across every principal-stress crossing.  P2's principal-space map cannot
    // be used either: it projects onto PLANES in closed form, and the Hoek-Brown
    // meridian is the curve
    //     f_ij = y_i - y_j - sigma_ci * (s - mb*y_i/sigma_ci)^a ,   0 < a < 1
    // (tension positive, y1 >= y2 >= y3; the header stores tension positive and
    // negates internally, so `sigma_geo = -sigma` turns its (sigma_1, sigma_3)
    // compression-positive pair into this (y_i, y_j) tension-positive one).
    // So P3 keeps P2's PRINCIPAL SPACE and its back-transform, and replaces the
    // closed-form projection with a small Newton on the curved surface.
    //
    // THE NATURAL VARIABLE (a requirement, not taste -- ADR-97 P0 measured it).
    // With y1 as the Newton unknown the map DIVERGES near the apex: the surface
    // exists only for arg = s - mb*y1/sigma_ci >= 0, the first step from the
    // elastic predictor overshoots (arg = -2.2e-03 at iteration 1 on the oracle's
    // own near-apex trial) and the next Jacobian is singular.  Substituting
    //     arg = (w^2)^(1/a)   <=>   sigma_ci * arg^a = sigma_ci * w^2
    // -- the surface's OWN variable, since sigma_ci*arg^a IS the strength term --
    // gives  y1 = T - (sigma_ci/mb) (w^2)^(1/a)  and  f_13 = y1 - y3 - sigma_ci w^2:
    // polynomial-smooth (2/a = 3.98 for the ADR-94 fixture) and feasible for ANY
    // real w.  No clipping, no line search, no feasibility guard.  `(w*w)^(1/a)`
    // rather than `w^(2/a)` because the latter is NaN for a negative iterate and
    // the two agree for w > 0; the residual is then even in w, so a Newton that
    // crosses zero converges to -w* and returns the SAME stress.
    //
    // NORMALIZED FLOW DIRECTION.  |m| ~ arg^(a-1) blows up exactly where the
    // near-apex returns land (460.6 there against 3.9 on an ordinary face point).
    // The residual therefore carries m/|m| with the multiplier rescaled by |m|;
    // the returned stress is invariant, and the worst-case Newton count over the
    // oracle's 400-trial scan plus every finite-difference stencil point drops
    // from 6 to 5.
    //
    // THE POTENTIAL (ADR-97 P3 decision -- implemented here, NOT in the shipped
    // `HoekBrown_PF::g`).  `HoekBrown_PF::g` is evaluated in the un-negated
    // (tension-positive) frame while `HoekBrown_YF` negates first, and then
    // destructures the ascending tuple as [sigma3, sigma2, sigma1] and feeds the
    // tree's most COMPRESSIVE principal into arg = mb_psi*sigma3/sigma_ci + s.
    // On any compressive state that arg is negative, so `g` always takes its
    // `else` branch and collapses to a TRESCA potential: `HB_mb_psi` has NO
    // effect, the flow is exactly non-dilatant, it sits 32.86 deg off the surface
    // normal even at mb_psi == mb (where the deck is asking for ASSOCIATED flow),
    // and at the apex all six of its flow directions have negative trace, so the
    // hydrostatic direction is not in its return cone and a trial pushed past the
    // tensile corner has NO return to the apex at all -- the mechanism behind the
    // residual recorded in tests/test_adr94_hlist_hb.py.  `Closest_Point` uses the
    // FRAME-CONSISTENT potential below in its own code path; the shipped `g` is
    // untouched, so `Backward_Euler` stays byte-identical (ADR-97 D1) and the
    // difference is PINNED by a test rather than hidden.  Fixing the shipped `g`
    // is the owner's separate PR (see LEDGER_quirks and the ADR implementation
    // log).  When HB_mb_psi == HB_mb the potential below IS the yield function's
    // own gradient, i.e. associated flow is exact under Closest_Point.
    //
    // THE TENSION BRANCH IS INERT ON THE SURFACE, so there is no tension-plane
    // return and none is missing: f_shear <= 0 already forces y1 <= T (above the
    // apex the clamp leaves f_shear = y1 - y3 > 0 for every non-hydrostatic
    // state), so the composite's zero set is the HB shear surface plus its own
    // natural vertex T*[1,1,1] -- 0 of 4000 surface points sampled by the P0
    // oracle have f_tension winning, with min(f_shear - f_tension) = 1.0e-04.
    // Outside the domain the tension branch wins exactly on {y1 > T and y3 > T},
    // which IS the header's CHECK_APEX_REGION set.  A Rankine return onto the
    // plane y1 = T from the oracle's [645, 265, 255] trial leaves
    // f_shear = +123.34 kPa (2.31 % of the strength scale) -- INADMISSIBLE.  The
    // shear/tension corner IS the apex.
    //
    // REGION SELECTION.  Not boundary PLANES (P2's construction needs a
    // piecewise-linear surface).  APEX: the elastic domain's tangent cone at the
    // vertex is exactly the negative octant (the meridian meets the hydrostatic
    // axis vertically, df/dy1 -> inf), so the apex region is exactly
    // apex + cone{D3 m_ij(apex)} -- for associated flow that is apex + D3*(the
    // positive octant), i.e. C3 (x_tr - apex) >= 0.  The header's
    // CHECK_APEX_REGION is the EUCLIDEAN octant x_tr >= T, and D3*(octant) is a
    // STRICT subset of the octant, so the header always OVER-claims: 32 of the
    // oracle's 400 scanned trials disagree, and on [645, 265, 255] APEX_STRESS
    // would commit T*[1,1,1] instead of the correct face return -- 122.6 kPa
    // (2.29 % of the strength scale) of silently lost strength.  `Closest_Point`
    // classifies here, where E is in scope, and does not call CHECK_APEX_REGION.
    // FACE vs EDGE: a curved surface has no precomputed boundary planes, but
    // their defining property survives -- the boundary is the ruled surface where
    // the FACE return lands on the edge, so the face return's OWN margins
    // y1 - y2 and y2 - y3 are the exact signed boundary functions.
    //
    // The elastic domain is CONVEX (-sigma_ci*arg^a is convex in y1: arg affine,
    // 0 < a < 1; a max of such functions over the six permutations stays convex),
    // so the elastic-metric closest point is unique and KKT is SUFFICIENT.  The
    // oracle's 400-trial scan closes at a max return-direction residual of
    // 1.672e-14 with every Koiter multiplier >= 0, which is a PROOF that no trial
    // is misclassified, not a spot check.
    //
    // `Backward_Euler` is untouched (ADR-97 D1).

    // Every constant the HB return needs, read ONCE per call from the yield
    // function and the plastic flow direction (never re-read inside the Newton).
    struct hb_consts_t
    {
        double sigci, mb, s, a, mb_psi, T, scale;
        Eigen::Matrix3d D3;
    };

    // arg = s - mb_psi*y_i/sigma_ci, floored so pow(0, a-1) cannot become inf.
    // The floor is RELATIVE to s (ADR-94 M5: an absolute one reintroduces the
    // unit dependence) and is nine orders below anything a converged return can
    // produce, so it can only ever fire on a diverging iterate.
    double hb_arg_floor(const hb_consts_t& C) const
    {
        const double sc = (C.s > 0.0) ? C.s : 1.0;
        return 1e-30 * sc;
    }

    // df/dy for the pair (i, j) -- the TRUE gradient of the surface this map
    // returns to.  Used only for the start guess (the residual carries the yield
    // row in the natural variable instead), and NEVER the header's central
    // difference of the composite max().
    void hb_a_ij(const hb_consts_t& C, const Eigen::Vector3d& y, int i, int j,
                 Eigen::Vector3d& out) const
    {
        double arg = C.s - C.mb * y(i) / C.sigci;
        const double fl = hb_arg_floor(C);
        if (!(arg > fl)) arg = fl;
        out.setZero();
        out(i) = 1.0 + C.a * C.mb * std::pow(arg, C.a - 1.0);
        out(j) = -1.0;
    }

    // The FRAME-CONSISTENT Hoek-Brown potential's gradient (see the block comment):
    //     g_ij = y_i - y_j - sigma_ci (s - mb_psi y_i/sigma_ci)^a
    // which is a_ij exactly when mb_psi == mb (associated flow).
    void hb_m_ij(const hb_consts_t& C, const Eigen::Vector3d& y, int i, int j,
                 Eigen::Vector3d& out) const
    {
        double argp = C.s - C.mb_psi * y(i) / C.sigci;
        const double fl = hb_arg_floor(C);
        if (!(argp > fl)) argp = fl;
        out.setZero();
        out(i) = 1.0 + C.a * C.mb_psi * std::pow(argp, C.a - 1.0);
        out(j) = -1.0;
    }

    // dm/dy has ONE non-zero entry, at (i, i).  This IS the curvature term that
    // the algorithmic tangent needs and that the P3 mutation gate removes.
    double hb_dm_ii(const hb_consts_t& C, const Eigen::Vector3d& y, int i) const
    {
        double argp = C.s - C.mb_psi * y(i) / C.sigci;
        const double fl = hb_arg_floor(C);
        if (!(argp > fl)) argp = fl;
        return -C.a * (C.a - 1.0) * C.mb_psi * C.mb_psi / C.sigci
               * std::pow(argp, C.a - 2.0);
    }

    // The header's OWN composite, mirrored in principal space (y DESCENDING).
    double hb_f_composite(const hb_consts_t& C, const Eigen::Vector3d& yd,
                          double* fs_out, double* ft_out) const
    {
        double arg = C.s - C.mb * yd(0) / C.sigci;
        if (!(arg > 0.0)) arg = 0.0;
        const double fs = yd(0) - yd(2) - C.sigci * std::pow(arg, C.a);
        const double ft = yd(0) - C.T;
        if (fs_out) *fs_out = fs;
        if (ft_out) *ft_out = ft;
        return (fs > ft) ? fs : ft;
    }

    // Round-off floor of f at y -- NOT a fudge factor.  |df/dy1| = 1 + a mb
    // arg^(a-1) DIVERGES at the apex, so a last-ulp error in y1 carries
    // eps*|y1|*|df/dy1| into f and an ABSOLUTE gate is unattainable within
    // ~1e-2 kPa of the vertex (the oracle measures max|f| = 1.08e-10 against its
    // own floor of 1.16e-08 there).  This is the Hoek-Brown instance of ADR-94's
    // f_relative_tol lesson: scale the yield tolerance by the GRADIENT, not by
    // sigma_ci.  64 ulp rather than the oracle's 4, because the C++ path
    // reassembles through a spectral decomposition the oracle does not.
    double hb_f_floor(const hb_consts_t& C, const Eigen::Vector3d& yd) const
    {
        double arg = C.s - C.mb * yd(0) / C.sigci;
        const double fl = hb_arg_floor(C);
        if (!(arg > fl)) arg = fl;
        const double cond = 1.0 + C.a * C.mb * std::pow(arg, C.a - 1.0);
        const double y0 = (yd(0) < 0) ? -yd(0) : yd(0);
        const double mag = (y0 > C.scale) ? y0 : C.scale;
        return 64.0 * MACHINE_EPSILON * mag * cond;
    }

    // Region layout.  nshape + ndl == 4 on every region, so the Newton system is
    // 4x4 everywhere.  key[] indexes the sextant's surfaces in the SAME order
    // P2 uses: 0 = (1,3), 1 = (2,3), 2 = (1,2), in 0-based principal slots.
    static void hb_layout(int region, int& ndl, int& nshape, int* key)
    {
        if (region == 0)      { ndl = 1; nshape = 3; key[0] = 0; key[1] = -1; }
        else if (region == 1) { ndl = 2; nshape = 2; key[0] = 0; key[1] = 1; }
        else                  { ndl = 2; nshape = 2; key[0] = 2; key[1] = 0; }
    }

    // Unknowns -> principal stresses, with the EDGE CONSTRAINT BUILT IN.  One
    // surface row then suffices on an edge: f_13 and f_23 share a root iff
    // y1 == y2 (y - sigma_ci arg(y)^a is strictly increasing), and likewise
    // f_12 / f_13 iff y2 == y3.
    void hb_y_of(const hb_consts_t& C, const double* z, int region,
                 Eigen::Vector3d& y) const
    {
        const double y1 = C.T - (C.sigci / C.mb) * std::pow(z[0] * z[0], 1.0 / C.a);
        if (region == 0)      { y(0) = y1; y(1) = z[1]; y(2) = z[2]; }
        else if (region == 1) { y(0) = y1; y(1) = y1;   y(2) = z[1]; }
        else                  { y(0) = y1; y(1) = z[1]; y(2) = z[1]; }
    }

    // Residual (and, when J != 0, the ANALYTIC Jacobian at the same point) of
    //     R_y = y(z) - x + sum_k dl_k D3 mhat_k(y)          (3 rows)
    //     R_f = (y_i - y_j) - sigma_ci w^2                  (1 row)
    bool hb_assemble(const hb_consts_t& C, const Eigen::Vector3d& x, int region,
                     const double* z, Eigen::Vector4d& R,
                     Eigen::Matrix4d* J, Eigen::Vector3d& y_out,
                     Eigen::Matrix3d& Ydz_out) const
    {
        static const int PIJ[3][2] = {{0, 2}, {1, 2}, {0, 1}};
        int ndl = 0, nshape = 0, key[2] = {0, 0};
        hb_layout(region, ndl, nshape, key);

        Eigen::Vector3d y;
        hb_y_of(C, z, region, y);
        y_out = y;

        const double w = z[0];
        const double dy1dw = -(C.sigci / C.mb) * (2.0 / C.a) * w
                             * std::pow(w * w, 1.0 / C.a - 1.0);

        Eigen::Matrix3d Ydz = Eigen::Matrix3d::Zero();
        if (region == 0)
        {
            Ydz(0, 0) = dy1dw;
            Ydz(1, 1) = 1.0;
            Ydz(2, 2) = 1.0;
        }
        else if (region == 1)
        {
            Ydz(0, 0) = dy1dw;  Ydz(1, 0) = dy1dw;
            Ydz(2, 1) = 1.0;
        }
        else
        {
            Ydz(0, 0) = dy1dw;
            Ydz(1, 1) = 1.0;    Ydz(2, 1) = 1.0;
        }
        Ydz_out = Ydz;

        Eigen::Vector3d corr = Eigen::Vector3d::Zero();
        Eigen::Vector3d mh[2], Dmh[2];
        double dscale[2] = {0.0, 0.0};
        int    dslot[2]  = {0, 0};
        for (int k = 0; k < ndl; ++k)
        {
            const int t = key[k], i = PIJ[t][0], jj = PIJ[t][1];
            Eigen::Vector3d mv;
            hb_m_ij(C, y, i, jj, mv);
            const double nm = mv.norm();
            if (!(nm > 0.0)) return false;
            mh[k]  = mv / nm;
            Dmh[k] = C.D3 * mh[k];
            corr  += z[nshape + k] * Dmh[k];
            dscale[k] = hb_dm_ii(C, y, i) / nm;
            dslot[k]  = i;
        }

        const int fi = PIJ[key[0]][0], fj = PIJ[key[0]][1];
        R(0) = y(0) - x(0) + corr(0);
        R(1) = y(1) - x(1) + corr(1);
        R(2) = y(2) - x(2) + corr(2);
        R(3) = (y(fi) - y(fj)) - C.sigci * w * w;
        for (int i = 0; i < 4; ++i)
            if (!(R(i) == R(i))) return false;

        if (J == 0) return true;

        // A = I + sum_k dl_k D3 (dmhat_k/dy).  mhat = m/|m| so
        // dmhat/dy = (1/|m|)(I - mhat mhat^T) dm/dy, and dm/dy has ONE non-zero
        // column, so each k contributes exactly one rank-one column update.
        // Two surfaces of an edge CAN share the slot (LINE2 is (1,2) and (1,3),
        // both differentiating y1), hence the accumulation.
        Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
        for (int k = 0; k < ndl; ++k)
        {
            const int i = dslot[k];
            Eigen::Vector3d col = Eigen::Vector3d::Zero();
            col(i) = 1.0;
            col -= mh[k](i) * mh[k];
            Eigen::Vector3d contrib = C.D3 * col;
            A.col(i) += (z[nshape + k] * dscale[k]) * contrib;
        }

        Eigen::Matrix3d AY = A * Ydz;
        Eigen::Matrix4d& Jm = *J;
        Jm.setZero();
        for (int c = 0; c < nshape; ++c)
            for (int r = 0; r < 3; ++r) Jm(r, c) = AY(r, c);
        for (int k = 0; k < ndl; ++k)
            for (int r = 0; r < 3; ++r) Jm(r, nshape + k) = Dmh[k](r);

        Eigen::Vector3d gf = Eigen::Vector3d::Zero();
        gf(fi) = 1.0;
        gf(fj) = -1.0;
        for (int c = 0; c < nshape; ++c)
            Jm(3, c) = gf.dot(Ydz.col(c)) - ((c == 0) ? 2.0 * C.sigci * w : 0.0);
        return true;
    }

    // Elastic predictor in the natural variable (seed sigma_ci w^2 with the
    // trial's own spread, keep the start sorted) PLUS the first-order
    // cutting-plane multiplier f/(a . D3 mhat).  Measured by the P0 oracle over
    // 400+ random trials: with dl = 0 the worst case is 6 Newton iterations,
    // with this seed it is 5.
    void hb_start(const hb_consts_t& C, const Eigen::Vector3d& x, int region,
                  double* z) const
    {
        int ndl = 0, nshape = 0, key[2] = {0, 0};
        hb_layout(region, ndl, nshape, key);
        double sp = x(0) - x(2);
        if (!(sp > 1e-12)) sp = 1e-12;
        z[0] = std::sqrt(sp / C.sigci);
        z[1] = z[2] = z[3] = 0.0;
        if (region == 0)      { z[1] = x(1); z[2] = x(2); }
        else if (region == 1) { z[1] = x(2); }
        else                  { z[1] = x(1); }
        const double y1 = C.T - (C.sigci / C.mb) * std::pow(z[0] * z[0], 1.0 / C.a);
        if (z[1] > y1) z[1] = y1;
        if (region == 0 && z[2] > y1) z[2] = y1;

        Eigen::Vector3d y0;
        hb_y_of(C, z, region, y0);
        Eigen::Vector3d yd = y0;
        for (int p = 0; p < 2; ++p)
            for (int q = 0; q < 2 - p; ++q)
                if (yd(q) < yd(q + 1)) { const double t = yd(q); yd(q) = yd(q + 1); yd(q + 1) = t; }

        Eigen::Vector3d av, mv;
        hb_a_ij(C, y0, 0, 2, av);
        hb_m_ij(C, y0, 0, 2, mv);
        const double nm = mv.norm();
        double f0 = hb_f_composite(C, yd, 0, 0);
        if (f0 < 0.0) f0 = 0.0;
        double den = 0.0;
        if (nm > 0.0)
        {
            Eigen::Vector3d Dm = C.D3 * (mv / nm);
            den = av.dot(Dm);
        }
        z[nshape] = (den > 0.0) ? (f0 / den) : 0.0;
    }

    // Newton on one region.  Converged when the residual reaches the yield
    // tolerance floor OR the correction stagnates at round-off -- the second
    // exit matters near the apex, where the residual floor is set by the
    // diverging gradient rather than by the tolerance.
    bool hb_solve(const hb_consts_t& C, const Eigen::Vector3d& x, int region,
                  double* z, int max_iter, double tol_R, int& iters,
                  Eigen::Vector3d& y, Eigen::Matrix3d& Ydz,
                  Eigen::Matrix4d& J) const
    {
        hb_start(C, x, region, z);
        iters = 0;
        Eigen::Vector4d R;
        for (int it = 0; it < max_iter; ++it)
        {
            if (!hb_assemble(C, x, region, z, R, &J, y, Ydz)) return false;
            double rn = 0.0;
            for (int i = 0; i < 4; ++i)
            {
                const double a = (R(i) < 0) ? -R(i) : R(i);
                if (a > rn) rn = a;
            }
            if (rn <= tol_R) return true;

            Eigen::FullPivLU<Eigen::Matrix4d> lu(J);
            if (!lu.isInvertible()) return false;
            Eigen::Vector4d dz = lu.solve(R);
            double zn = 1.0, dn = 0.0;
            for (int i = 0; i < 4; ++i)
            {
                if (!(dz(i) == dz(i))) return false;
                const double az = (z[i] < 0) ? -z[i] : z[i];
                if (az > zn) zn = az;
                const double ad = (dz(i) < 0) ? -dz(i) : dz(i);
                if (ad > dn) dn = ad;
            }
            for (int i = 0; i < 4; ++i) z[i] -= dz(i);
            ++iters;
            if (dn <= 1e-13 * zn)
            {
                // at the round-off floor: re-form at the accepted point so `J`,
                // `y` and `Ydz` belong to the state that is returned.
                if (!hb_assemble(C, x, region, z, R, &J, y, Ydz)) return false;
                return true;
            }
        }
        return false;
    }

    // Is `d` in the cone of return directions at the apex?  The generators are
    // the six limiting directions D3 m_ij(apex), taken at apex*(1 - 1e-9) exactly
    // as the P0 oracle does (m_ij diverges AT the vertex).  For associated flow
    // they collapse to D3*(the positive octant), i.e. C3 (x - apex) >= 0, which
    // is what the oracle's README states as the exact apex region.
    //
    // Membership is decided by the DUAL (facet) test: for a POINTED,
    // FULL-DIMENSIONAL polyhedral cone in R^3 the extreme rays of the dual are
    // among the pairwise cross products of the generators, so d is in the cone
    // iff every supporting normal has n . d >= 0.  Pointedness is guaranteed
    // here, not assumed: every m_ij has component sum M - 1 > 0 (M = 1 + a mb_psi
    // arg^(a-1) > 1) and D3 is positive definite with
    // trace(D3 m) = (3 lambda + 2 mu) sum(m) > 0, so all six generators lie
    // strictly inside the half space trace > 0.  Near-duplicate generators (the
    // associated case pairs them) are removed first, so a cross product of two
    // copies of the same ray cannot masquerade as a facet.
    bool hb_in_apex_cone(const hb_consts_t& C, const Eigen::Vector3d& d) const
    {
        static const int OP[6][2] = {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}};
        Eigen::Vector3d G[6];
        int nG = 0;
        Eigen::Vector3d ye = Eigen::Vector3d::Constant(C.T - 1e-9 * C.T);
        for (int k = 0; k < 6; ++k)
        {
            Eigen::Vector3d mv;
            hb_m_ij(C, ye, OP[k][0], OP[k][1], mv);
            Eigen::Vector3d gv = C.D3 * mv;
            const double n = gv.norm();
            if (!(n > 0.0)) continue;
            gv /= n;
            bool dup = false;
            for (int q = 0; q < nG; ++q)
                if ((G[q] - gv).norm() <= 1e-4) { dup = true; break; }
            if (!dup) G[nG++] = gv;
        }
        const double nd = d.norm();
        if (!(nd > 0.0)) return true;
        if (nG < 3) return false;

        bool any_facet = false;
        for (int p = 0; p < nG; ++p)
            for (int q = p + 1; q < nG; ++q)
            {
                Eigen::Vector3d nv = G[p].cross(G[q]);
                const double nn = nv.norm();
                if (!(nn > 1e-3)) continue;
                nv /= nn;
                int npos = 0, nneg = 0;
                for (int k = 0; k < nG; ++k)
                {
                    const double v = nv.dot(G[k]);
                    if (v >  1e-9) ++npos;
                    if (v < -1e-9) ++nneg;
                }
                double sgn = 0.0;
                if (nneg == 0)      sgn =  1.0;
                else if (npos == 0) sgn = -1.0;
                else continue;
                any_facet = true;
                if (sgn * nv.dot(d) < -1e-12 * nd) return false;
            }
        return any_facet;
    }

    // apex (exact elastic-metric cone test) -> face -> edge (ordered by the FACE
    // return's own margins, which ARE the exact signed boundary functions of a
    // curved surface) -> apex as the fallback.  Returns the region, and leaves
    // the ACCEPTED region's converged Newton state in z / y / Ydz / J so the
    // tangent reuses the factorization-point data (no second solve).
    int hb_classify(const hb_consts_t& C, const Eigen::Vector3d& x,
                    int max_iter, double tol_R, double* z, int& iters,
                    Eigen::Vector3d& y, Eigen::Matrix3d& Ydz,
                    Eigen::Matrix4d& J) const
    {
        if (hb_in_apex_cone(C, Eigen::Vector3d(x - Eigen::Vector3d::Constant(C.T))))
            return 3;

        int order[2] = {1, 2};
        double zf[4] = {0., 0., 0., 0.};
        int itf = 0;
        Eigen::Vector3d yf;
        Eigen::Matrix3d Yf;
        Eigen::Matrix4d Jf;
        if (hb_solve(C, x, 0, zf, max_iter, tol_R, itf, yf, Yf, Jf))
        {
            const double m0 = yf(0) - yf(1);
            const double m1 = yf(1) - yf(2);
            if (m0 >= 0.0 && m1 >= 0.0)
            {
                for (int i = 0; i < 4; ++i) z[i] = zf[i];
                iters = itf; y = yf; Ydz = Yf; J = Jf;
                return 0;
            }
            if (m0 >= m1) { order[0] = 2; order[1] = 1; }
        }

        for (int q = 0; q < 2; ++q)
        {
            const int reg = order[q];
            double ze[4] = {0., 0., 0., 0.};
            int ite = 0;
            Eigen::Vector3d ye;
            Eigen::Matrix3d Ye;
            Eigen::Matrix4d Je;
            if (!hb_solve(C, x, reg, ze, max_iter, tol_R, ite, ye, Ye, Je))
                continue;
            int ndl = 0, nshape = 0, key[2] = {0, 0};
            hb_layout(reg, ndl, nshape, key);
            bool ok = true;
            for (int k = 0; k < ndl; ++k)
                if (!(ze[nshape + k] >= -1e-12)) ok = false;
            Eigen::Vector3d yd = ye;
            for (int p = 0; p < 2; ++p)
                for (int r = 0; r < 2 - p; ++r)
                    if (yd(r) < yd(r + 1)) { const double t = yd(r); yd(r) = yd(r + 1); yd(r + 1) = t; }
            if (ok && hb_f_composite(C, yd, 0, 0) <= 1e-8 * C.scale)
            {
                for (int i = 0; i < 4; ++i) z[i] = ze[i];
                iters = ite; y = ye; Ydz = Ye; J = Je;
                return reg;
            }
        }
        return 3;
    }

    int cp_hb_return(const VoigtVector& depsilon,
                     const VoigtVector& sigma_tr,
                     const VoigtMatrix& Eelastic,
                     double tol_f, int max_iter)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        (void) depsilon;

        // ---- surface / potential constants -------------------------------
        hb_consts_t C;
        double pf_sigci = 0.0, pf_s = 0.0, pf_a = 0.0;
        if (!yf.cp_hb_face_params(iv_storage, parameters_storage,
                                  C.sigci, C.mb, C.s, C.a)
                || !pf.cp_hb_flow_params(iv_storage, parameters_storage,
                                         pf_sigci, C.mb_psi, pf_s, pf_a))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - this yield function / plastic flow pair is marked as a"
                   << " principal-space Hoek-Brown family but does not supply the"
                   << " surface parameters (ADR-97 P3) -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        {
            // The two functors share ONE HB_sigci / HB_s / HB_a parameter object
            // (utuple_storage de-duplicates parameters by type), so a mismatch
            // here can only mean the traits were widened to a pair that does not
            // actually share them.  Refuse rather than return to a surface the
            // potential does not know about.
            const double d1 = pf_sigci - C.sigci, d2 = pf_s - C.s, d3 = pf_a - C.a;
            if ((d1 > 1e-12 * C.sigci) || (d1 < -1e-12 * C.sigci)
                    || (d2 > 1e-12) || (d2 < -1e-12)
                    || (d3 > 1e-12) || (d3 < -1e-12))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the Hoek-Brown yield function and plastic flow"
                       << " direction disagree on (sigma_ci, s, a)"
                       << " -- rejecting step (ADR-97 P3)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        if (!(C.sigci > 0.0) || !(C.mb > 100 * MACHINE_EPSILON) || !(C.s > 0.0)
                || !(C.a > 0.0) || !(C.a < 1.0) || !(C.mb_psi > 0.0))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - integration_method Closest_Point needs HB_sigci > 0,"
                   << " HB_mb > 0, HB_s > 0, HB_mb_psi > 0 and 0 < HB_a < 1"
                   << " (got sigci = " << C.sigci << ", mb = " << C.mb
                   << ", s = " << C.s << ", a = " << C.a
                   << ", mb_psi = " << C.mb_psi << ")"
                   << " -- rejecting step (ADR-97 P3). Use Backward_Euler."
                   << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        if (C.mb_psi > C.mb * (1.0 + 1e-12))
        {
            // mb_psi > mb is dilation in excess of friction: the potential's own
            // arg goes NEGATIVE on the surface, where the flow direction does not
            // exist.  No registered deck does this; refuse rather than clamp.
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - HB_mb_psi (" << C.mb_psi << ") exceeds HB_mb (" << C.mb
                   << "): the Hoek-Brown flow potential is then undefined on part"
                   << " of its own yield surface -- rejecting step (ADR-97 P3)"
                   << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        C.T = C.s * C.sigci / C.mb;
        C.scale = C.sigci * std::pow(C.s, C.a);

        // ---- isotropy of the elastic tangent -----------------------------
        const double lam = Eelastic(0, 1);
        const double Gmu = (Eelastic(0, 0) - Eelastic(0, 1)) / 2.0;
        {
            const double sE = ((Eelastic(0, 0) < 0) ? -Eelastic(0, 0) : Eelastic(0, 0))
                              + ((Gmu < 0) ? -Gmu : Gmu);
            const double t = 1e-8 * sE;
            const double d02 = Eelastic(0, 2) - lam;
            const double d12 = Eelastic(1, 2) - lam;
            const double d33 = Eelastic(3, 3) - Gmu;
            const double d11 = Eelastic(1, 1) - Eelastic(0, 0);
            if (!(Gmu > 0.0)
                    || (d02 > t) || (d02 < -t) || (d12 > t) || (d12 < -t)
                    || (d33 > t) || (d33 < -t) || (d11 > t) || (d11 < -t))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space Hoek-Brown return needs an"
                       << " ISOTROPIC elastic tangent -- rejecting step (ADR-97 P3)"
                       << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                C.D3(i, j) = lam + ((i == j) ? 2.0 * Gmu : 0.0);

        // ---- spectral decomposition of the TRIAL stress, DESCENDING ------
        Eigen::Matrix3d st;
        st << sigma_tr.v11(), sigma_tr.v12(), sigma_tr.v13(),
              sigma_tr.v12(), sigma_tr.v22(), sigma_tr.v23(),
              sigma_tr.v13(), sigma_tr.v23(), sigma_tr.v33();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(st);
        if (es.info() != Eigen::Success)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the trial stress eigen-decomposition failed"
                   << " -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        Eigen::Vector3d xv;
        Eigen::Matrix3d Q;
        for (int i = 0; i < 3; ++i)          // eigenvalues come out ASCENDING
        {
            xv(i) = es.eigenvalues()(2 - i);
            Q.col(i) = es.eigenvectors().col(2 - i);
        }

        double xmax = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            const double a = (xv(i) < 0) ? -xv(i) : xv(i);
            if (a > xmax) xmax = a;
        }
        const double scale_x = (xmax > C.scale) ? xmax : C.scale;
        double tol_R = 1e-14 * scale_x;
        {
            const double floor_R = 64.0 * MACHINE_EPSILON * scale_x;
            if (floor_R > tol_R) tol_R = floor_R;
        }

        // ---- region + return ---------------------------------------------
        double z[4] = {0., 0., 0., 0.};
        int iters = 0;
        Eigen::Vector3d yv = Eigen::Vector3d::Zero();
        Eigen::Matrix3d Ydz = Eigen::Matrix3d::Zero();
        Eigen::Matrix4d J = Eigen::Matrix4d::Zero();
        const int region = hb_classify(C, xv, max_iter, tol_R, z, iters,
                                       yv, Ydz, J);

        Eigen::Matrix3d dydx = Eigen::Matrix3d::Zero();
        if (region == 3)
        {
            // The vertex does not move, so the consistent tangent is EXACTLY the
            // zero matrix -- rank 0, as for the Drucker-Prager (P1) and
            // Mohr-Coulomb (P2) apices.  That is rank deficient by construction
            // and WILL make an element whose every Gauss point sits at the apex
            // singular, which is the true state of affairs.
            yv = Eigen::Vector3d::Constant(C.T);
            iters = 1;
        }
        else
        {
            int ndl = 0, nshape = 0, key[2] = {0, 0};
            hb_layout(region, ndl, nshape, key);
            Eigen::FullPivLU<Eigen::Matrix4d> lu(J);
            if (!lu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the converged Hoek-Brown Jacobian is singular, so"
                       << " no algorithmic tangent exists (region " << region
                       << ") -- rejecting step (ADR-97 P3)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            // dz/dx from the converged system (dR/dx = -I on the three stress
            // rows), then dy/dx = (dy/dz)(dz/dx).  The curvature term dm/dy is
            // ALREADY inside J -- it is the `dl * D3 dmhat/dy` block of A.
            Eigen::Matrix<double, 4, 3> rhs = Eigen::Matrix<double, 4, 3>::Zero();
            for (int i = 0; i < 3; ++i) rhs(i, i) = 1.0;
            Eigen::Matrix<double, 4, 3> Z = lu.solve(rhs);
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double v = 0.0;
                    for (int k = 0; k < nshape; ++k) v += Ydz(r, k) * Z(k, c);
                    dydx(r, c) = v;
                }
            if (!(yv(0) >= yv(1) - 1e-9 * scale_x)
                    || !(yv(1) >= yv(2) - 1e-9 * scale_x))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the Hoek-Brown return broke the principal ordering"
                       << " (region " << region << ", y = " << yv(0) << " "
                       << yv(1) << " " << yv(2) << ") -- rejecting step (ADR-97 P3)"
                       << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // ---- admissibility, in TWO places, for two different reasons ------
        // (a) EXACTLY, in principal space, against the composite this map
        //     mirrors.  Perfectly conditioned; this is the check that catches a
        //     coding error.  Its tolerance carries the GRADIENT-scaled floor
        //     (hb_f_floor): |df/dy1| diverges at the apex, so an absolute gate is
        //     unattainable within ~1e-2 kPa of the vertex.
        // (b) LOOSELY, against the header's own 6D f recomputed from the
        //     reassembled Q diag(y) Q^T.  Redundant with (a) up to the spectral
        //     round-trip, which is exactly what it is there to bound.
        double sig_max = 0.0;
        Eigen::Matrix3d Sret = Q * yv.asDiagonal() * Q.transpose();
        VoigtVector sigma_ret(Sret(0, 0), Sret(1, 1), Sret(2, 2),
                              Sret(0, 1), Sret(1, 2), Sret(0, 2));
        for (int i = 0; i < 6; ++i)
        {
            if (!(sigma_ret(i) == sigma_ret(i)))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - NaN in the principal-space Hoek-Brown return"
                       << " (region " << region << ") -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            const double av = (sigma_ret(i) < 0) ? -sigma_ret(i) : sigma_ret(i);
            if (av > sig_max) sig_max = av;
        }
        double scale_ref = yf.strength_scale(iv_storage, parameters_storage);
        if (scale_ref < 0) scale_ref = -scale_ref;
        const double ref_mag = (sig_max > scale_ref) ? sig_max : scale_ref;
        const double f_floor = hb_f_floor(C, yv);
        {
            double gate = 1e-10 * ref_mag;
            if (tol_f > gate)   gate = tol_f;
            if (f_floor > gate) gate = f_floor;
            const double f_pr = hb_f_composite(C, yv, 0, 0);
            if (!(f_pr <= gate))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE the"
                       << " Hoek-Brown surface (f = " << f_pr << " > tol = " << gate
                       << ", region " << region << ") -- rejecting step (ADR-97 P3)"
                       << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        {
            double gate = 1e-8 * ref_mag;
            if (tol_f > gate)   gate = tol_f;
            if (f_floor > gate) gate = f_floor;
            const double f_ret = yf(sigma_ret, iv_storage, parameters_storage);
            if (!(f_ret <= gate))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the reassembled Hoek-Brown return is outside this"
                       << " yield function's own surface (f = " << f_ret
                       << " > tol = " << gate << ", region " << region
                       << ") -- rejecting step (ADR-97 P3)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // ---- the tangent, rotated to the global Voigt frame ---------------
        // Identical construction to P2 (cp_principal_return): the normal block is
        // dy/dx, the shear slots the eigenprojection ROTATION term
        // (y_i - y_j)/(x_i - x_j) with the l'Hopital limit dy_i/dx_i - dy_i/dx_j
        // on a degenerate trial eigenvalue, and Rs^-1 built as the Voigt image of
        // E -> Q^T E Q rather than inverted numerically.
        VoigtMatrix Tp = VoigtMatrix::Zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Tp(i, j) = dydx(i, j);
        {
            const double eps_deg = 1e-9 * ((xmax > scale_ref) ? xmax : scale_ref);
            static const int SIJ[3][2] = {{0, 1}, {1, 2}, {0, 2}};   // slots 3,4,5
            for (int sslot = 0; sslot < 3; ++sslot)
            {
                const int i = SIJ[sslot][0], j = SIJ[sslot][1];
                const double dx = xv(i) - xv(j);
                const double adx = (dx < 0) ? -dx : dx;
                Tp(3 + sslot, 3 + sslot) = (adx > eps_deg)
                                           ? (yv(i) - yv(j)) / dx
                                           : (dydx(i, i) - dydx(i, j));
            }
        }
        VoigtMatrix Rs = VoigtMatrix::Zero();
        VoigtMatrix Rsi = VoigtMatrix::Zero();
        {
            static const int BI[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}};
            for (int k = 0; k < 6; ++k)
            {
                Eigen::Matrix3d Ek = Eigen::Matrix3d::Zero();
                Ek(BI[k][0], BI[k][1]) = 1.0;
                Ek(BI[k][1], BI[k][0]) = 1.0;
                Eigen::Matrix3d F  = Q * Ek * Q.transpose();
                Eigen::Matrix3d Fi = Q.transpose() * Ek * Q;
                for (int r = 0; r < 6; ++r)
                {
                    Rs(r, k)  = F(BI[r][0], BI[r][1]);
                    Rsi(r, k) = Fi(BI[r][0], BI[r][1]);
                }
            }
        }
        VoigtMatrix RT   = Rs * Tp;
        VoigtMatrix RTR  = RT * Rsi;
        VoigtMatrix Calg = RTR * Eelastic;

        // ---- plastic strain increment ------------------------------------
        VoigtVector dep = VoigtVector::Zero();
        {
            Eigen::Matrix<double, 6, 6> Ee;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Ee(i, j) = Eelastic(i, j);
            Eigen::FullPivLU< Eigen::Matrix<double, 6, 6> > elu(Ee);
            if (!elu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the elastic tangent is singular"
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            Eigen::Matrix<double, 6, 1> rhs6;
            for (int i = 0; i < 6; ++i) rhs6(i) = sigma_tr(i) - sigma_ret(i);
            Eigen::Matrix<double, 6, 1> d6 = elu.solve(rhs6);
            for (int i = 0; i < 6; ++i) dep(i) = d6(i);
        }

        TrialStress         = sigma_ret;
        TrialPlastic_Strain = CommitPlastic_Strain + dep;
        cp_last_iterations  = iters;
        cp_apply_tangent_policy(Calg, Eelastic);

        if (ladruno_strict_rejects("Closest_Point (Hoek-Brown)", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;
        return 0;
    }

    int Closest_Point(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        if (!ladruno_cp_supported)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - integration_method Closest_Point is not implemented for"
                   << " this yield function / plastic flow / hardening family"
                   << " (ADR-97 D3). Supported today: VonMises and DruckerPrager"
                   << " with Null / Linear (scalar or tensor) / ArmstrongFrederick"
                   << " hardening (P1), and the perfectly plastic MohrCoulomb and"
                   << " MohrCoulombTensionCutoff families where BOTH the yield"
                   << " function and the plastic flow direction are of that family"
                   << " (P2). MIXED pairings such as MohrCoulomb_YF x VonMises_PF or"
                   << " VonMises_YF x MohrCoulomb_PF are NOT supported: the"
                   << " principal-space map assumes both are piecewise linear, and"
                   << " the smooth 6D map cannot use MohrCoulomb's Lode-angle"
                   << " gradient. The Hoek-Brown family (P3) is supported for"
                   << " HoekBrown_YF x HoekBrown_PF only, on the same"
                   << " matched-pair rule. StiffSoil is P5."   // Ladruno (ADR-97 wp/97d)
                   << " Use Backward_Euler." << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        const int    max_iter = INT_OPT_n_max_iterations[ASDP_TAG];
        const double tol_f    = yf_tolerance();

        VoigtVector depsilon = strain_incr;
        const VoigtVector& sigma_n = CommitStress;

        iv_storage.revert_all();
        dsigma.setZero();
        // The non-Algorithmic tangent types read `depsilon_elpl` as the `depsilon`
        // argument of pf()/hardening(); Backward_Euler leaves it stale (it never
        // sets it), Closest_Point does not.
        depsilon_elpl = depsilon;

        VoigtMatrix Eelastic = et(sigma_n, parameters_storage);
        VoigtVector sigma_tr = sigma_n + Eelastic * depsilon;

        TrialStrain         = CommitStrain + depsilon;
        TrialStress         = sigma_tr;
        TrialPlastic_Strain = CommitPlastic_Strain;
        cp_last_iterations  = 0;

        // The stress- and IV-row tolerances get ADR-94 M5's relative floor too, so
        // the same problem in Pa and in kPa converges identically.  Every internal
        // variable in the P1 families is stress dimensioned (a yield stress, a
        // cohesion, a back stress), so one scaled tolerance serves all of them.
        double scale = yf.strength_scale(iv_storage, parameters_storage);
        if (scale < 0) scale = -scale;
        double tol_s = DBL_OPT_stress_absolute_tol[ASDP_TAG];
        {
            const double rel_s = DBL_OPT_f_relative_tol[ASDP_TAG] * scale;
            if (rel_s > tol_s) tol_s = rel_s;
        }

        const double f_tr = yf(sigma_tr, iv_storage, parameters_storage);
        if (!(f_tr == f_tr))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the elastic predictor's yield value is NaN"
                   << " -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        if (f_tr <= tol_f)
        {
            // Genuinely elastic.  Note this is NOT Backward_Euler's second
            // disjunct ("f merely DECREASED by more than tol => call it elastic"),
            // which accepts an end state that is still outside the surface; ADR-84
            // P2a had to gate that with strict_convergence.  A new integrator has
            // no compatibility reason to inherit it.
            Stiffness = Eelastic;
            return 0;
        }

        // Ladruno (ADR-97 wp/97c): the Mohr-Coulomb family leaves here and never
        // reaches the smooth 6D Newton below -- its yield function has no usable
        // 6D gradient (Drucker-Prager substitution for |theta| >= 29 deg,
        // central differences otherwise).  MohrCoulombTensionCutoff first offers
        // the trial to ADR-84's `special_return`, whose closed-form cutoff
        // face / Rankine edge / MC-cutoff corner / compound-corner / apex returns
        // and RAW Koiter tangent are reused verbatim -- no geometry is
        // re-derived; when that hook declines (cutoff inactive at the trial, or
        // an MC-dominant trial none of whose cutoff features validated) the plain
        // Mohr-Coulomb principal return takes over and re-checks the COMPOSITE f.
        // Ladruno (ADR-97 wp/97d): the Hoek-Brown family leaves here.  Like the
        // Mohr-Coulomb family below it never reaches the smooth 6D Newton -- its
        // yield function has NO analytic gradient at all (a central difference
        // of the COMPOSITE max() over the six raw Voigt slots).  Unlike it, the
        // surface is CURVED, so the return is a small Newton in the surface's own
        // variable rather than a closed-form projection.  Note this branch is
        // taken BEFORE the yf_has_apex block: HoekBrown_YF declares an apex, but
        // its CHECK_APEX_REGION is the Euclidean octant and over-claims the apex
        // region (P0 finding 7), so the elastic-metric cone test inside
        // cp_hb_return is used instead.
        if constexpr (ladruno_cp_principal_family == 3)
        {
            return cp_hb_return(depsilon, sigma_tr, Eelastic, tol_f, max_iter);
        }
        else if constexpr (ladruno_cp_principal_family != 0)
        {
            if constexpr (yf_has_special_return<YieldFunctionType>::value)
            {
                VoigtVector sigma_sr, dep_sr;
                VoigtMatrix stiff_sr;
                int sr_quality = SR_QUALITY_EXACT;
                if (yf.special_return(sigma_tr, Eelastic, tol_f,
                                      iv_storage, parameters_storage,
                                      sigma_sr, dep_sr, stiff_sr, sr_quality))
                {
                    if (sr_quality == SR_QUALITY_FALLBACK
                            && INT_OPT_strict_convergence[ASDP_TAG] != 0)
                    {
                        opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                               << ") - special_return could not land this trial on any"
                               << " exact cutoff feature (face/edge/corner/compound);"
                               << " the only remaining candidate is the conservative"
                               << " vertex projection, which would discard the"
                               << " deviatoric state -- rejecting step"
                               << " (strict_convergence)" << endln;
                        return LADRUNO_MATERIAL_REFUSED;
                    }
                    TrialStress         = sigma_sr;
                    TrialPlastic_Strain = CommitPlastic_Strain + dep_sr;
                    cp_last_iterations  = 1;
                    // ADR-84 P3: stiff_sr is the RAW active-set (Koiter) tangent,
                    // i.e. the consistent tangent of THAT exact return -- which is
                    // precisely what `Algorithmic` means here.
                    cp_apply_tangent_policy(stiff_sr, Eelastic);
                    if (ladruno_strict_rejects("Closest_Point (special_return)", TrialStress))
                        return LADRUNO_MATERIAL_REFUSED;
                    return 0;
                }
            }
            return cp_principal_return(depsilon, sigma_tr, Eelastic, tol_f);
        }

        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr, tol_f))
            {
                int rc = -1;
                if (cp_apex_return(depsilon, sigma_tr, Eelastic, tol_f, tol_s, max_iter, rc))
                    return rc;
                // else: the YF's apex failed its own admissibility gate; fall
                // through to the generic Newton below.
                TrialStress = sigma_tr;
            }
        }

        const int n_iv = cp_n_iv();
        const int N = 6 + n_iv + 1;
        if (N > ASDP_CP_MAXN)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - this specialization needs a " << N << "x" << N
                   << " Newton system, over the ASDP_CP_MAXN = " << ASDP_CP_MAXN
                   << " cap -- rejecting step (raise the cap in "
                   << "ASDPlasticMaterial3D.h)" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        cp_vector_t x(N), R(N), dx(N);
        cp_matrix_t J(N, N);

        // Elastic-predictor start: s = s_tr, q = q_n, dl = 0.
        for (int i = 0; i < 6; ++i) x(i) = sigma_tr(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i) x(off + i) = internal_variable.trial_value(i);
                off += nq;
            });
        }
        x(N - 1) = 0.0;

        VoigtVector m_conv = VoigtVector::Zero();
        bool converged = false;
        int iter = 0;
        double rs = 0.0, rq = 0.0, rf = 0.0;

        for (iter = 0; iter < max_iter; ++iter)
        {
            if (!cp_assemble(x, N, sigma_tr, depsilon, m_conv, R, &J))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - NaN in the closest-point residual at iteration "
                       << iter << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }

            rs = 0.0; rq = 0.0;
            for (int i = 0; i < 6; ++i)
            {
                const double a = (R(i) < 0) ? -R(i) : R(i);
                if (a > rs) rs = a;
            }
            for (int i = 6; i < N - 1; ++i)
            {
                const double a = (R(i) < 0) ? -R(i) : R(i);
                if (a > rq) rq = a;
            }
            rf = (R(N - 1) < 0) ? -R(N - 1) : R(N - 1);

            if (rf <= tol_f && rs <= tol_s && rq <= tol_s)
            {
                converged = true;
                break;      // J is the Jacobian AT the converged point: reuse it
            }

            Eigen::FullPivLU<cp_matrix_t> lu(J);
            if (!lu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular closest-point Jacobian at iteration " << iter
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            dx = lu.solve(R);
            for (int i = 0; i < N; ++i)
                if (!(dx(i) == dx(i)))
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - NaN in the Newton correction at iteration " << iter
                           << " -- rejecting step" << endln;
                    return LADRUNO_MATERIAL_REFUSED;
                }
            x -= dx;
        }

        if (!converged)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the closest-point Newton exhausted " << max_iter
                   << " iterations: |f| = " << rf << " (tol " << tol_f
                   << "), |R_sigma| = " << rs << ", |R_q| = " << rq
                   << " (tol " << tol_s << ") -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        cp_last_iterations = iter;

        const double dl = x(N - 1);
        if (dl < 0.0)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - converged with a NEGATIVE plastic multiplier (" << dl
                   << "), which is inadmissible -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        TrialPlastic_Strain = CommitPlastic_Strain + dl * m_conv;

        if (INT_OPT_tangent_operator_type[ASDP_TAG]
                == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
        {
            // The converged system depends on epsilon_{n+1} only through s_tr, and
            // d(s_tr)/d(eps) = E(sigma_n).  Differentiating R(x(eps)) = 0:
            //     J * [dsigma ; dq ; ddl] = [E*deps ; 0 ; 0]
            //     C_alg = (J^{-1})_{sigma,sigma} * E
            // Solved with the ONE factorization already in hand -- the Schur form
            // in the ADR is for the doc, not for the code.  C_alg is UNSYMMETRIC
            // whenever m != n (non-associated Drucker-Prager, etabar != eta), so
            // models using it must run on an unsymmetric solver (UmfPack; NOT
            // ProfileSPD, and NOT PARDISO's symmetric -matrixType).
            if constexpr (el_is_stress_dependent<ElasticityType>::value)
            {
                static bool warned_dEds = false;
                if (!warned_dEds)
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - this elasticity model is stress dependent and has"
                           << " no ELASTICITY_STRESS_DERIVATIVE, so tangent_type"
                           << " Algorithmic is missing the E,sigma : (sigma - sigma_tr)"
                           << " term (ADR-97 D6/P5). The committed stress is"
                           << " unaffected." << endln;
                    warned_dEds = true;
                }
            }
            Eigen::FullPivLU<cp_matrix_t> lu(J);
            if (!lu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the converged closest-point Jacobian is singular, so"
                       << " no algorithmic tangent exists -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            cp_matrix_t rhs(N, 6);
            rhs.setZero();
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) rhs(i, j) = Eelastic(i, j);
            cp_matrix_t Z = lu.solve(rhs);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Stiffness(i, j) = Z(i, j);
        }
        else
        {
            // Secant / Continuum / Elastic / Numerical_Algorithmic_* are evaluated
            // at the CONVERGED closest-point state, but they are still the same
            // operators ADR-94 M3 measured at 80 / 57 / 103 / 5 % against a central
            // difference of the material's own response: none of them is the
            // tangent of THIS map.  Only `Algorithmic` is.
            ComputeTangentStiffness();
        }

        if (ladruno_strict_rejects("Closest_Point", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;

        return 0;
    }

    // Ladruno (ADR-94 wp/94f): the one-shot APEX projection, factored out of
    // `Backward_Euler` so it can be reached from TWO places -- the predictor
    // classification (as before) and the flank-failure fallback below.  The body
    // is the wp/94c projection VERBATIM; only the trial stress it projects and
    // the diagnostic prefix are parameters now.
    //
    //     sigma   = sigma_apex
    //     d eps^p = inv(E) : (sigma_trial - sigma_apex)
    //     dlambda = <m_apex, d eps^p>_e / <m_apex, m_apex>_e   (>= 0)
    //
    // Returns  1 : apex return taken; TrialStress / TrialPlastic_Strain / the
    //              internal variables / Stiffness are all set, caller returns 0.
    //          0 : the yield function's own apex is not on its own surface, so
    //              nothing was written; caller falls through.
    //         -1 : same as 0, but strict_convergence says refuse loudly.
    //
    // The caller owns the state it hands in: at the fallback site TrialStress,
    // TrialPlastic_Strain and the IVs must be reset to the elastic predictor
    // BEFORE calling, because the failed flank iteration has already moved them.
    int be_apex_project(const VoigtVector & depsilon,
                        const VoigtVector & sigma_trial,
                        const VoigtMatrix & Eelastic,
                        double tol_yf, bool be_strict, const char* origin)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // Copy out of the YF's return buffer immediately: `apex_stress`
        // and `df_dsigma_ij` share one mutable member per functor.
        const VoigtVector sigma_apex = yf.apex_stress(iv_storage, parameters_storage);

        // The apex the YF names must actually BE on its own surface.  This
        // is what makes the classification safe: a yield function whose
        // apex geometry is wrong (or whose region test misfires) falls
        // through to the generic return map instead of committing a
        // fabricated stress.  NaN-safe: !(x <= tol) is true for NaN.
        const double f_apex = yf(sigma_apex, iv_storage, parameters_storage);
        const double f_apex_abs = (f_apex < 0) ? -f_apex : f_apex;
        if (!(f_apex_abs <= tol_yf))
        {
            if (be_strict)
            {
                opserr << "ASDPlasticMaterial3D::Backward_Euler (tag " << ASDP_TAG
                       << ") - " << origin << ": the trial stress classifies as an APEX state but the "
                       << "yield function's own apex is not on its surface: |f(sigma_apex)| = "
                       << f_apex_abs << " > tol = " << tol_yf
                       << " -- rejecting step (strict_convergence)" << endln;
                return -1;
            }
            return 0;   // fall through to the generic return map
        }

        // d eps^p = inv(E) * (sigma_trial - sigma_apex)
        const VoigtVector rhs = sigma_trial - sigma_apex;
        VoigtVector dep = VoigtVector::Zero();
        {
            Eigen::Matrix<double, 6, 6> E_eig;
            Eigen::Matrix<double, 6, 1> rhs_eig;
            for (int i = 0; i < 6; ++i)
            {
                rhs_eig(i) = rhs(i);
                for (int j = 0; j < 6; ++j) E_eig(i, j) = Eelastic(i, j);
            }
            Eigen::Matrix<double, 6, 1> dep_eig;
            auto chol = E_eig.selfadjointView<Eigen::Lower>().llt();
            if (chol.info() == Eigen::Success) dep_eig = chol.solve(rhs_eig);
            else                               dep_eig = E_eig.ldlt().solve(rhs_eig);
            for (int i = 0; i < 6; ++i) dep(i) = dep_eig(i);
        }

        // Plastic multiplier: least-squares projection of dep onto the
        // apex flow direction.  Both are ENGINEERING-strain-like Voigt
        // vectors, so the inner product is the engineering contraction
        // (ADR-94 wp/94c, B5) -- not the stress-like one the dead code
        // used.
        const VoigtVector m_apex = pf(depsilon, sigma_apex, iv_storage, parameters_storage);
        const double m_dot_m = tensor_dot_engineering_strain_like(m_apex, m_apex);
        double dLambda_apex = 0.0;
        if (m_dot_m > MACHINE_EPSILON)
        {
            dLambda_apex = tensor_dot_engineering_strain_like(m_apex, dep) / m_dot_m;
            if (dLambda_apex < 0.0) dLambda_apex = 0.0;   // keep lambda >= 0
        }

        TrialStress          = sigma_apex;
        TrialPlastic_Strain  = TrialPlastic_Strain + dep;

        // Internal variables: ONE hardening evaluation, at the apex.
        // Every yield function that opts into `yf_has_apex` today is
        // perfectly plastic (Null hardening), so this term is exactly
        // zero for them and the IVs are unchanged; it is written this
        // way so a hardening Drucker-Prager does not silently freeze.
        iv_storage.apply([&](auto & internal_variable)
        {
            auto h = internal_variable.hardening_function(depsilon, m_apex, TrialStress, parameters_storage);
            internal_variable.trial_value += dLambda_apex * h;
        });

        // Tangent.  The honest continuum operator at a perfectly plastic
        // apex is ZERO: the stress is pinned at sigma_apex, so no strain
        // increment that stays in the apex region changes it.  That is
        // rank-deficient by construction and will make an element whose
        // every Gauss point sits at the apex singular -- which is the
        // true state of affairs, and why it is only produced for the
        // tangent types the user opts into.  Secant (the DEFAULT) blends
        // it with the elastic operator, exactly as the special_return
        // path above does, and stays invertible.
        {
            VoigtMatrix apex_stiff = VoigtMatrix::Zero();
            using TOT = ASDPlasticMaterial3D_Tangent_Operator_Type;
            switch (INT_OPT_tangent_operator_type[ASDP_TAG])
            {
            case TOT::Elastic:
                Stiffness = Eelastic;
                break;
            case TOT::Continuum:
            case TOT::Algorithmic:
            case TOT::Numerical_Algorithmic_FirstOrder:
            case TOT::Numerical_Algorithmic_SecondOrder:
                Stiffness = apex_stiff;
                break;
            case TOT::Secant:
            default:
                Stiffness = VoigtMatrix((apex_stiff + Eelastic) / 2.0);
                break;
            }
        }

        return 1;
    }

    int Backward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;

        // Ladruno (HB/StiffSoil integration, ledger row 337): hoisted earlier for the elastic-exit checks below
        int    max_iter = INT_OPT_n_max_iterations[ASDP_TAG];
        double tol_yf   = yf_tolerance(); 

        // -------- setup
        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon = strain_incr;

        const VoigtVector& sigma   = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        // -------- elastic predictor
        dsigma      = Eelastic * depsilon;
        TrialStress = sigma + dsigma;
        TrialStrain = epsilon + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        // cout << "BE - TrialStress = " << TrialStress.transpose() <<  endl;
        // cout << "BE - CommitStress = " << CommitStress.transpose() <<  endl;
        
        const double yf_val_start = yf(sigma,        iv_storage, parameters_storage);
        // cout << "----->  yf_val_start = " << yf_val_start <<  endl;
        
        const double yf_val_end   = yf(TrialStress,  iv_storage, parameters_storage);
        // cout << "----->  yf_val_end = " << yf_val_end << " tol_yf = " << tol_yf <<  endl;

        // Ladruno (ADR-84 P2a): strict_convergence gate on the f-decreasing early-exit.
        // Upstream, the second disjunct accepts a step as "elastic" whenever f merely
        // DECREASED by more than tol_yf, even if the end state is still outside the
        // surface (yf_val_end > 0) -- from an inadmissible committed state this
        // perpetuates the violation instead of correcting it. With strict_convergence
        // the exit additionally requires yf_val_end <= tol_yf (elastic unloading from
        // a plastic state has yf_end < 0, so it still exits here). Flag off: identical
        // to upstream.
        const bool be_strict = INT_OPT_strict_convergence[ASDP_TAG] != 0; // Ladruno (ADR-84 P2a)

        // purely elastic or moving deeper inside the surface
        if ( (yf_val_start <= 0.0 && yf_val_end <= 0.0)
                || ( (yf_val_start - yf_val_end > tol_yf)
                     && (!be_strict || yf_val_end <= tol_yf) ) ) { // Ladruno (ADR-84 P2a)
            // cout << "BE - ELASTIC!" << endl << endl;
            Stiffness = Eelastic;
            return 0;
        }


        // Ladruno (ADR-94 wp/94f): the elastic predictor is needed AGAIN at the
        // flank-failure fallback below, and the Newton loop overwrites TrialStress
        // in place, so keep a copy.
        const VoigtVector be_sigma_trial_elastic = TrialStress;

        // Deal with the APEX if needed
        // -------- APEX return (one-shot projection) -------------------------------
        // Ladruno (ADR-94 wp/94c, B4): REVIVED.  The whole body below was commented
        // out, so `check_apex_region` was called and its answer discarded for every
        // yield function that declares `yf_has_apex` -- Drucker-Prager then ran the
        // flank return map past sqrt(J2) = 0 and committed NaN on hydrostatic
        // tension.  The projection itself now lives in `be_apex_project` above.
        //
        // Ladruno (ADR-94 wp/94f): LAYER (a) -- ELASTIC-METRIC classification.
        // `check_apex_region` is EUCLIDEAN, `(p - p_apex) >= eta*q`, and says so in
        // its own comment; the exact test is `(p - p_apex) >= (K*etabar/G)*q`, which
        // at etabar = 0 (zero dilatancy -- the ADR-95 Prandtl footing deck) collapses
        // to `p >= p_apex`, because a non-dilatant flank return cannot change the
        // mean stress AT ALL.  The Euclidean test is then strictly too narrow: every
        // over-apex state with q > eta*(p - p_apex) is routed to the flank map, whose
        // scalar Newton cannot close f and exhausts (measured: 435 refusals, floor at
        // s/B 0.01122).  `cp_apex_region` -- written for ADR-97's `Closest_Point`,
        // reused verbatim here -- is the same test done family-agnostically inside the
        // integrator, where E (hence K and G) and the plastic flow direction (hence
        // etabar) are both in scope: take the linearised cone step and ask whether the
        // returned deviator has FLIPPED sign.  It is applied ONLY to yield functions
        // that opt in via `yf_apex_elastic_metric` (today: Drucker-Prager, whose apex
        // is a cone vertex in the (p, sqrt(J2)) half-plane); every other yield
        // function keeps its own Euclidean answer and is byte-identical here.
        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            bool be_in_apex = yf.check_apex_region(TrialStress, iv_storage, parameters_storage);
            if constexpr (yf_apex_elastic_metric<YieldFunctionType>::value)
            {
                if (!be_in_apex)
                    be_in_apex = cp_apex_region(depsilon, TrialStress, Eelastic, yf_val_end, tol_yf);
            }
            if (be_in_apex)
            {
                const int apex_rc = be_apex_project(depsilon, TrialStress, Eelastic,
                                                    tol_yf, be_strict, "predictor-classified apex");
                if (apex_rc > 0)  return 0;
                if (apex_rc < 0)  return LADRUNO_MATERIAL_REFUSED;
                // apex_rc == 0: fall through to the generic return map
            }
        }

        // Ladruno (ADR-84 P0): predictor-classified exact special-region return
        // (tension-cutoff face/edge, MC∩TC corner, apex cone) for yield functions
        // that opt in via yf_has_special_return. The YF validates every candidate
        // against its own yield values before accepting, so no f > 0 state can be
        // committed through this path. For every existing YF the trait is false
        // and this block compiles away -- behavior is bit-identical.
        if constexpr (yf_has_special_return<YieldFunctionType>::value)
        {
            VoigtVector sigma_sr, dep_sr;
            VoigtMatrix stiff_sr;
            int sr_quality = SR_QUALITY_EXACT;
            if (yf.special_return(TrialStress, Eelastic, tol_yf,
                                  iv_storage, parameters_storage,
                                  sigma_sr, dep_sr, stiff_sr, sr_quality))
            {
                // Ladruno (ADR-84 P3): a conservative terminal vertex projection
                // replaces a confined deviatoric state with a hydrostatic one.
                // It keeps the "no f > 0 commit" guarantee, but it is NOT the
                // right answer, and committing it silently is what makes a
                // material failure read downstream as a modelling problem.
                // Under strict_convergence, refuse loudly instead.
                if (sr_quality == SR_QUALITY_FALLBACK && be_strict)
                {
                    opserr << "ASDPlasticMaterial3D::Backward_Euler (tag " << ASDP_TAG
                           << ") - special_return could not land this trial on any exact "
                           << "cutoff feature (face/edge/corner/compound); the only "
                           << "remaining candidate is the conservative vertex projection, "
                           << "which would discard the deviatoric state"
                           << " -- rejecting step (strict_convergence)" << endln;
                    // Ladruno (ADR-84 -> ADR-86b): was a bare -1. LadrunoBrick's
                    // update() paths (updateHypo/formEAStrue) now check ONLY the
                    // LADRUNO_MATERIAL_REFUSED sentinel, not a blanket `< 0` -- a
                    // bare -1 here was silently SWALLOWED under -geom hypo/EAS,
                    // which is strictly worse than the pre-sentinel behaviour this
                    // rejection was written to replace. See LadrunoMaterialStatus.h.
                    return LADRUNO_MATERIAL_REFUSED;
                }

                TrialStress          = sigma_sr;
                TrialPlastic_Strain += dep_sr;
                // Internal variables: the opted-in YFs are perfectly plastic
                // (Null hardening), so no IV update is required here.

                // Ladruno (ADR-84 P3): stiff_sr is the RAW active-set (Koiter)
                // tangent. Apply the SAME tangent-operator policy the generic
                // path applies, instead of the hardcoded secant blend P0 used --
                // that blend ignored `tangent_type` and, at the MC-cutoff corner,
                // reported ~2e6 kPa of stiffness where the true tangent is ~0
                // (87% error), stalling the global Newton. Default is Secant, so
                // this is byte-identical unless the user asks for another type.
                using TOT = ASDPlasticMaterial3D_Tangent_Operator_Type;
                switch (INT_OPT_tangent_operator_type[ASDP_TAG])
                {
                case TOT::Elastic:
                    Stiffness = Eelastic;
                    break;
                case TOT::Continuum:
                case TOT::Algorithmic:
                    Stiffness = stiff_sr;
                    break;
                // Numerical_Algorithmic_* differentiate compute_local_stress(),
                // a simplified map that does NOT call this hook -- it would be
                // strictly worse here than the exact active-set tangent, so the
                // numerical types intentionally fall back to it.
                case TOT::Numerical_Algorithmic_FirstOrder:
                case TOT::Numerical_Algorithmic_SecondOrder:
                    Stiffness = stiff_sr;
                    break;
                case TOT::Secant:
                default:
                    Stiffness = VoigtMatrix((stiff_sr + Eelastic) / 2.0);
                    break;
                }
                return 0;
            }
        }

        // -------- plastic correction (Backward Euler)

        double dLambda = 0.0;

        bool be_converged = false; // Ladruno (ADR-84 P2a): only read when strict_convergence is on

        // Ladruno (ADR-94 wp/94f): a HARD failure inside the flank Newton (singular
        // local tangent, NaN, or -- under strict_convergence -- the plastic
        // inconsistency branch) no longer returns on the spot.  It breaks out so the
        // apex fallback below gets its chance; if the fallback declines, the very
        // same refusal is issued, so the flag-off / no-apex behaviour is unchanged.
        bool be_flank_failed = false;
        const char* be_flank_reason = nullptr;

        // cout << "BE - Plastic! Begin iterations----------" << endl << endl;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            // directions at current (trial) end state
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);

            // use total step depsilon as proxy argument; no derivatives of m required
            const VoigtVector& m = pf(depsilon, TrialStress, iv_storage, parameters_storage);

            // hardening at end state
            const double H = yf.hardening(depsilon, m, TrialStress, iv_storage, parameters_storage);

            // consistency residual Φ(σ_{n+1}) = 0
            const double Phi = yf(TrialStress, iv_storage, parameters_storage);

            if (std::abs(Phi) < tol_yf) {
                GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], iter);
                GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], std::abs(Phi));

                // cout << "  =>  n = " << n.transpose() << endl;
                // cout << "  =>  m = " << m.transpose() << endl;
                // cout << "  =>  H = " << H << endl;
                be_converged = true; // Ladruno (ADR-84 P2a)
                break; // converged
            }

            // dΦ/dλ ≈ - n^T E m + H
            // const double dPhi_dLambda = - (n.transpose() * Eelastic * m) + H;
            // Ladruno (ADR-94 wp/94c, B5): ONE convention, everywhere.  With a
            // Voigt-convention `n` (shear slots already carrying the symmetric
            // partner) the plastic modulus is a PLAIN contraction n_i (E m)_i --
            // the factor 2 lives in n, not in the dot product.  This line was the
            // ONLY site of six that doubled it a second time; it happened to be
            // right for the old tensor-convention von Mises and wrong for the five
            // Voigt-convention families (MC/HB/MCTC/StiffSoil), and no single
            // choice here could be right for both.  Now it matches
            // ComputeTangentStiffness, compute_local_stress, Forward_Euler,
            // Forward_Euler_Subincrement and Backward_Euler_LineSearch.
            // GCC: bind the Eigen product to a local before contracting.
            const VoigtVector Em = Eelastic * m;
            const double nEm = n.dot(Em);
            const double dPhi_dLambda = H - nEm;   // == - n^T E m + H

            if (std::abs(dPhi_dLambda) < MACHINE_EPSILON) {
                // singular local tangent
                cout << " SINGULAR LOCAL TANGENT - FAILING!" << endl;
                cout << "  =>  n = " << n.transpose() << endl;
                cout << "  =>  m = " << m.transpose() << endl;
                cout << "  =>  H = " << H << endl;
                be_flank_failed = true;  // Ladruno (ADR-94 wp/94f): apex fallback first, then refuse
                be_flank_reason = "singular local tangent (|H - n:E:m| < eps)";
                break;                   // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }

            const double deltaLambda = - Phi / dPhi_dLambda;

            // cout << " ---- iter = " << iter << " / " << max_iter << endl;
            // cout << "  =>  n = " << n.transpose() << endl;
            // cout << "  =>  m = " << m.transpose() << endl;
            // cout << "  =>  H = " << H << endl;
            // cout << "  =>  nEm = " << nEm << endl;
            // cout << "  =>  Phi = " << Phi << endl;
            // cout << "  =>  dPhi_dLambda = H - nEm = " << dPhi_dLambda << endl;
            // cout << "  =>  deltaLambda = -Phi / dPhi_dLambda = " << deltaLambda << endl;
            // cout << "  =>  dLambda = " << dLambda << endl ;

            // keep λ >= 0
            if (dLambda + deltaLambda < 0.0) {
                cout << " PLASTIC INCONSISTENCY - ELASTIC STEP! (dLambda + deltaLambda < 0.0)" << endl << endl;
                // step cannot be plastic; fall back to elastic in this rare case
                // Ladruno (ADR-94 wp/94a): ADR-94 B3/H7 -- this branch commits the
                // ELASTIC PREDICTOR exactly and returns success, and it fires
                // before the be_strict exhaustion check below, so it was an
                // eighth silent-accept site inside the DEFAULT integrator.
                if (be_strict) {
                    // Ladruno (ADR-94 wp/94f): message deferred to the refusal site
                    // below, so a rescued step does not print a refusal it did not make.
                    be_flank_failed = true;
                    be_flank_reason = "plastic inconsistency (dLambda + deltaLambda < 0): the"
                                      " elastic predictor would be committed uncorrected";
                    break;
                }
                Stiffness = Eelastic;
                return 0;
            }

            // incremental updates (use delta to avoid re-summing from commit each iter)
            dLambda += deltaLambda;
            TrialStress          = TrialStress - deltaLambda * (Eelastic * m);
            TrialPlastic_Strain  = TrialPlastic_Strain + deltaLambda * m;
            // cout << "  =>  dLambda + deltaLambda = " << dLambda << endl;
            // cout << "  =>  CommitStress = " << CommitStress.transpose() << endl;
            // cout << "  =>  TrialStress = " << TrialStress.transpose() << endl;
            // cout << "  =>  TrialPlastic_Strain = " << TrialPlastic_Strain.transpose() << endl;
            iv_storage.apply([&](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon, m, TrialStress, parameters_storage);
                internal_variable.trial_value += deltaLambda * h;
                // cout << "  => " <<  internal_variable << endl;
            });

            // NaN guard
            const double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (!(norm_trial_stress == norm_trial_stress)) { // NaN chec
               cout << "NaN!" << endl;
                be_flank_failed = true;  // Ladruno (ADR-94 wp/94f): apex fallback first, then refuse
                be_flank_reason = "the flank return produced a NaN stress";
                break;                   // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }
        }
        // cout << "BE - END iterations----------" << endl << endl;

        // Ladruno (ADR-84 P2a): strict_convergence gate on the exhaustion-accept.
        // Upstream, falling out of the loop after max_iter iterations falls through
        // to ComputeTangentStiffness()/return 0, silently committing a non-converged
        // state with f > tol_yf as success. With strict_convergence, exhaustion with
        // |Phi| >= tol_yf fails loud so the element reports the failure upward.
        // (The last Newton update may have converged without being checked -- the
        // check runs at the top of the next iteration -- so re-evaluate yf here.)
        double be_Phi_final = 0.0;
        bool   be_exhausted = false;   // Ladruno (ADR-94 wp/94f)
        if (be_strict && !be_converged && !be_flank_failed)
        {
            be_Phi_final = std::abs(yf(TrialStress, iv_storage, parameters_storage));
            be_exhausted = (be_Phi_final >= tol_yf);
        }

        // Ladruno (ADR-94 wp/94f): LAYER (b) -- FLANK-FIRST APEX FALLBACK.
        // This block fires at EXACTLY the sites that would otherwise return
        // LADRUNO_MATERIAL_REFUSED (Newton exhaustion under strict_convergence, a
        // singular local tangent, a NaN stress, the strict plastic-inconsistency
        // branch) and only for yield functions that declare an apex.  A trial state
        // whose MEAN STRESS is already beyond the apex has no admissible flank
        // return -- with zero dilatancy the flank map cannot move p at all -- so the
        // apex projection, not a refusal, is the right answer.  Every other outcome
        // is unchanged: if the trial is not beyond the apex, or the yield function's
        // own apex fails the |f(sigma_apex)| <= tol_yf guard, the same refusal is
        // issued with the same sentinel.  Yield functions that opt in through
        // `yf_has_apex` today: DruckerPrager, MohrCoulomb<NO_HARDENING>,
        // HoekBrown<NO_HARDENING>, TensionCutoff<NO_HARDENING>; for the latter three
        // this can only convert a refusal into an admissible vertex state.
        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            if (be_flank_failed || be_exhausted)
            {
                // Rewind to the elastic predictor: the failed iteration has already
                // moved the stress, the plastic strain and the internal variables.
                iv_storage.revert_all();
                TrialStress         = be_sigma_trial_elastic;
                TrialPlastic_Strain = CommitPlastic_Strain;

                const VoigtVector sigma_apex_probe = yf.apex_stress(iv_storage, parameters_storage);
                const double p_apex  = sigma_apex_probe.meanStress();
                const double p_trial = be_sigma_trial_elastic.meanStress();

                if (p_trial > p_apex)
                {
                    const int apex_rc = be_apex_project(depsilon, be_sigma_trial_elastic, Eelastic,
                                                        tol_yf, be_strict, "flank-first apex fallback");
                    if (apex_rc > 0)
                    {
                        if (ladruno_strict_rejects("Backward_Euler (apex fallback)", TrialStress))
                            return LADRUNO_MATERIAL_REFUSED;
                        return 0;
                    }
                    // apex_rc <= 0: the yield function's own apex is not on its own
                    // surface -- refuse below, exactly as before the fallback existed.
                }
            }
        }

        if (be_flank_failed)
        {
            opserr << "ASDPlasticMaterial3D::Backward_Euler (tag " << ASDP_TAG
                   << ") - " << (be_flank_reason ? be_flank_reason : "flank return map failed")
                   << " -- rejecting step" << endln;
            // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            return LADRUNO_MATERIAL_REFUSED;
        }

        if (be_exhausted)
        {
            opserr << "ASDPlasticMaterial3D::Backward_Euler (tag " << ASDP_TAG
                   << ") - scalar Newton exhausted " << max_iter
                   << " iterations without converging: |Phi| = " << be_Phi_final
                   << " >= tol_yf = " << tol_yf
                   << " -- rejecting step (strict_convergence)" << endln;
            // Ladruno (ADR-84 -> ADR-86b): was a bare -1; same rationale as the
            // special_return fallback above. See LadrunoMaterialStatus.h.
            return LADRUNO_MATERIAL_REFUSED;
        }

        // Ladruno (ADR-94 wp/94a)
        if (ladruno_strict_rejects("Backward_Euler (post return-to-yield)", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;
        ComputeTangentStiffness();

        return 0;
    }

    int Backward_Euler_LineSearch(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;

        // ------------------ setup ------------------
        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;

        const VoigtVector& sigma   = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        // tolerancias
        const double tol_abs = (yf_tolerance() > 0.0) ? yf_tolerance() : 1e-10;
        const double tol_rel = 1e-8;
        const int    max_iter = 30;
        // Ladruno (ADR-94 wp/94a): refusal flag out of the solve lambda (which
        // can only say true/false, and whose `false` means "halve the step").
        bool ls_strict_refused = false;

        auto converged = [&](double Phi, double Phi_scale)->bool {
            double tol = std::max(tol_abs, tol_rel * std::max(1.0, Phi_scale));
            return std::abs(Phi) <= tol;
        };

        // ------------------ solver (un sub-paso Δε) ------------------
        auto solve_increment = [&](const VoigtVector& dEps)->bool
        {
            // Predictor elástico desde el estado commit
            TrialStrain = epsilon + dEps;
            TrialPlastic_Strain = CommitPlastic_Strain;
            TrialStress = sigma + Eelastic * dEps;

            // Estado interno trial = commit
            iv_storage.revert_all();

            const double yf_start = yf(sigma,       iv_storage, parameters_storage);
            const double yf_end   = yf(TrialStress, iv_storage, parameters_storage);

            // Misma lógica que usabas: puramente elástico o moviéndose "hacia adentro"
            if ( (yf_start <= 0.0 && yf_end <= 0.0) || (yf_start > yf_end) ) {
                Stiffness = Eelastic;
                // Ladruno (ADR-94 wp/94a)
                if (ladruno_strict_rejects("Backward_Euler_LineSearch", TrialStress)) {
                    ls_strict_refused = true;
                    return false;
                }
                return true;
            }

            // ---------- Newton escalar con backtracking ----------
            double dLambda = 0.0;

            // Arranque seguro (si el trial está fuera)
            {
                const VoigtVector& n0 = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m0 = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H0 = yf.hardening(dEps, m0,       TrialStress, iv_storage, parameters_storage);
                const double       Phi0 = yf(TrialStress, iv_storage, parameters_storage);
                const double       nEm0 = n0.dot(Eelastic * m0);
                const double       denom0 = nEm0 - H0; // == -(H0 - nEm0) == -dPhi/dλ

                if (denom0 > MACHINE_EPSILON && Phi0 > 0.0) {
                    // dλ inicial ≥ 0
                    double dL0 = Phi0 / denom0;
                    // límite físico simple
                    double dLmax = dEps.norm() / std::max(m0.norm(), 1e-16);
                    dLambda = std::min(dL0, dLmax);

                    // aplicar arranque
                    TrialStress         = TrialStress - dLambda * (Eelastic * m0);
                    TrialPlastic_Strain = TrialPlastic_Strain + dLambda * m0;
                    iv_storage.apply([&](auto & iv){
                        auto h = iv.hardening_function(dEps, m0, TrialStress, parameters_storage);
                        iv.trial_value += dLambda * h;
                    });
                }
            }

            bool newton_ok = false;

            for (int iter = 0; iter < max_iter; ++iter)
            {
                // (opcional) si E depende fuerte de σ, descomenta:
                // Eelastic = et(TrialStress, parameters_storage);

                const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H = yf.hardening(dEps, m,        TrialStress, iv_storage, parameters_storage);

                const double Phi = yf(TrialStress, iv_storage, parameters_storage);
                const double Phi_scale = std::abs(n.dot(TrialStress)) + 1.0;

                if (converged(Phi, Phi_scale)) { newton_ok = true; break; }

                const double nEm = n.dot(Eelastic * m);
                const double dPhi_dLambda = H - nEm;  // total deriv. aprox

                if (!std::isfinite(dPhi_dLambda) || std::abs(dPhi_dLambda) < 1e-20) {
                    // Singular → deja que el substepping maneje
                    newton_ok = false;
                    break;
                }

                // Paso de Newton propuesto
                double deltaLambda = - Phi / dPhi_dLambda;

                // Limitar tamaño (previene runaway)
                const double dLmax = dEps.norm() / std::max(m.norm(), 1e-16);
                if (std::abs(deltaLambda) > dLmax)
                    deltaLambda = (deltaLambda > 0.0 ? 1.0 : -1.0) * dLmax;

                // Enforce λ ≥ 0
                if (dLambda + deltaLambda < 0.0)
                    deltaLambda = -dLambda * 0.5; // reduce para no cruzar a negativo

                // -------- line search (backtracking) con Φ linealizado --------
                double alpha = 1.0;
                const double c = 1e-4;
                double dl_accepted = 0.0;

                for (int ls = 0; ls < 8; ++ls) {
                    const double dl = alpha * deltaLambda;

                    // Predicción lineal de Φ
                    const double Phi_pred = Phi + dPhi_dLambda * dl;

                    if (std::abs(Phi_pred) <= (1.0 - c*alpha) * std::abs(Phi)) {
                        dl_accepted = dl;
                        break;
                    }
                    alpha *= 0.5;
                }

                if (dl_accepted == 0.0) {
                    // No se pudo aceptar un paso útil → dejar a substepping
                    newton_ok = false;
                    break;
                }

                // Aplicar actualización aceptada
                const VoigtVector Em = Eelastic * m;

                TrialStress         = TrialStress - dl_accepted * Em;
                TrialPlastic_Strain = TrialPlastic_Strain + dl_accepted * m;

                iv_storage.apply([&](auto & iv){
                    auto h = iv.hardening_function(dEps, m, TrialStress, parameters_storage);
                    iv.trial_value += dl_accepted * h;
                });

                dLambda += dl_accepted;

                // Guard NaN
                if (!std::isfinite(TrialStress.squaredNorm())) {
                    newton_ok = false;
                    break;
                }
            }

            if (!newton_ok) return false;

            // -------- Optional one-shot polish: return-to-surface --------
            if (INT_OPT_return_to_yield_surface[ASDP_TAG]) {
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H_corr = yf.hardening(dEps, m_corr,   TrialStress, iv_storage, parameters_storage);
                const double       Phi_corr = yf(TrialStress, iv_storage, parameters_storage);

                const VoigtVector Em_corr = Eelastic * m_corr;
                const double nEm_corr = n_corr.dot(Em_corr);
                const double denom_corr = nEm_corr - H_corr;

                if (std::abs(Phi_corr) > tol_abs && std::abs(denom_corr) > MACHINE_EPSILON) {
                    const double dLambda_corr =  Phi_corr / denom_corr;
                    TrialStress         = TrialStress - dLambda_corr * Em_corr;
                    TrialPlastic_Strain = TrialPlastic_Strain + dLambda_corr * m_corr;

                    iv_storage.apply([&](auto & iv){
                        auto h = iv.hardening_function(dEps, m_corr, TrialStress, parameters_storage);
                        iv.trial_value += dLambda_corr * h;
                    });
                    dLambda += dLambda_corr;
                }
            }

            // -------- Tangente algorítmica consistente --------
            {
                const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H = yf.hardening(dEps, m,        TrialStress, iv_storage, parameters_storage);

                const VoigtVector Em = Eelastic * m;
                const double nEm = n.dot(Em);
                const double denom = nEm - H; // ojo: este es el de la fórmula de C_alg

                if (std::abs(denom) > MACHINE_EPSILON) {
                    // (6x6) = E - (6x1)*(1x6)/denom
                    const auto row_nE = (n.transpose() * Eelastic); // 1x6
                    Stiffness = Eelastic - (Em * row_nE) / denom;
                } else {
                    Stiffness = Eelastic; // fallback
                }
            }

            return true;
        };

        // ------------------ substepping automático ------------------
        VoigtVector dEps = depsilon;

        // intenta paso completo; si falla, corta a la mitad repetidamente (hasta 1/32)
        bool ok = false;
        for (int split = 0; split <= 5; ++split)   // 0..5 → 1, 2, 4, 8, 16, 32 subpasos
        {
            // resetear a commit antes de intentar este tamaño de paso
            TrialStress = sigma + Eelastic * dEps;   // predictor para flags de arriba
            TrialStrain = epsilon + dEps;
            TrialPlastic_Strain = CommitPlastic_Strain;
            iv_storage.revert_all();

            if (solve_increment(dEps)) { ok = true; break; }
            // Ladruno (ADR-94 wp/94a): a strict refusal is not a "try a smaller step".
            if (ls_strict_refused) break;

            // reducir paso y reintentar
            dEps *= 0.5;
        }

        // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
        if (ls_strict_refused) return LADRUNO_MATERIAL_REFUSED;
        if (!ok) return LADRUNO_MATERIAL_REFUSED;

        // Ladruno (ADR-94 wp/94a)
        if (ladruno_strict_rejects("Backward_Euler_LineSearch (post return-to-yield)", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;
        return 0;
    }



    std::pair<double, VoigtVector> CalculateLambdaM(
        const VoigtVector & thisSigma,
        const VoigtVector & depsilon_elpl,
        const parameters_storage_t &this_parameters_storage,
        const iv_storage_t &this_iv_storage)
    {
        VoigtMatrix Eelastic = et(thisSigma, this_parameters_storage);
        const VoigtVector& n = yf.df_dsigma_ij(thisSigma, this_iv_storage, this_parameters_storage);
        VoigtVector m = pf(depsilon_elpl, thisSigma, this_iv_storage, this_parameters_storage);
        double hardening = yf.hardening(depsilon_elpl, m, thisSigma, this_iv_storage, this_parameters_storage);

        double den = n.transpose() * Eelastic * m - hardening;
        double dLambda = (den != 0) ? (n.transpose() * Eelastic * depsilon_elpl).value() / den : 0;

        if (dLambda != dLambda)
        {
            cout << "CalculateLambdaM error" << endl;
            cout << "yf = " << yf(thisSigma, this_iv_storage, this_parameters_storage) << endl;
            cout << "thisSigma = " << thisSigma.transpose() << endl;
            cout << "depsilon_elpl = " << depsilon_elpl.transpose() << endl;
            cout << "n = " << n.transpose() << endl;
            cout << "m = " << m.transpose() << endl;
            cout << "hardening = " << hardening << endl;
            cout << "den = " << den << endl;
            cout << "Eelastic = " << Eelastic << endl;

            cout << "IVSTORAGE" << endln;
            this_iv_storage.print_components();

            cout << "PARAMSTORAGE" << endln;
            this_parameters_storage.print_components();
        }


        return std::make_pair(dLambda, m);
    }

    int Runge_Kutta_45_Error_Control_old(const VoigtVector & strain_incr)
    {
        // cout << "Runge_Kutta_45_Error_Control" << endl;

        int errorcode = -1;

        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;


        iv_storage.revert_all();


        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);

        dsigma = Eelastic * depsilon;


        TrialStress = CommitStress + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(CommitStress, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Runge_Kutta_45_Error_Control_old", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];


            VoigtVector next_Sigma = TrialStress;
            iv_storage_t next_iv_storage = iv_storage;
            VoigtVector next_EpsilonPl = CommitPlastic_Strain;

            VoigtVector prev_Sigma = TrialStress;
            VoigtVector prev_EpsilonPl = CommitPlastic_Strain;
            iv_storage_t prev_iv_storage = iv_storage;

            VoigtVector dSigma, dSigma1, dSigma2, dSigma3, dSigma4, dSigma5, dSigma6;
            VoigtVector dEpsilonPl, dEpsilonPl1, dEpsilonPl2, dEpsilonPl3, dEpsilonPl4, dEpsilonPl5, dEpsilonPl6;
            iv_storage_t iv_storage1 = iv_storage;
            iv_storage_t iv_storage2 = iv_storage;
            iv_storage_t iv_storage3 = iv_storage;
            iv_storage_t iv_storage4 = iv_storage;
            iv_storage_t iv_storage5 = iv_storage;


            int niter = 0;
            double maxStepError = 0;
            while (T < 1.0)
            {
                niter ++;

                VoigtVector dEPS = dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                //Delta 1
                std::tie(dLambda, m) = CalculateLambdaM(prev_Sigma, dEPS, parameters_storage, prev_iv_storage);
                dSigma1  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl1 = dLambda * m;
                iv_storage1.apply([&m, &dLambda, &prev_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, prev_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + 0.5 * dLambda * h;
                });
                next_Sigma =  prev_Sigma + 0.5 * dSigma1;

                //Delta 2
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage1);
                dSigma2  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl2 = dLambda * m;
                iv_storage2.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, this](auto & iv2)
                {
                    using VT1 = std::decay_t<decltype(iv2)>;
                    const VT1 &iv1 = iv_storage1.template get<VT1>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;

                    auto h2 = iv2.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH2 = dLambda * h2;

                    iv2.trial_value = iv2.committed_value + 0.25 * (DH1 + DH2);
                });
                next_Sigma =  prev_Sigma + 0.25 * (dSigma1 + dSigma2);

                //Delta 3
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage2);
                dSigma3  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl3 = dLambda * m;
                iv_storage3.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, this](auto & iv3)
                {
                    using VT = std::decay_t<decltype(iv3)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;

                    auto h3 = iv3.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH3 = dLambda * h3;

                    iv3.trial_value = iv3.committed_value +  -DH2 + 2 * DH3;
                });
                next_Sigma =  prev_Sigma - dSigma2 + 2 * dSigma3;

                //Delta 4
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage3);
                dSigma4  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl4 = dLambda * m;
                iv_storage4.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, this](auto & iv4)
                {
                    using VT = std::decay_t<decltype(iv4)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;

                    auto h4 = iv4.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH4 = dLambda * h4;

                    iv4.trial_value = iv4.committed_value +  (7 * DH1 + 10 * DH2 + DH4) / 27;
                });
                next_Sigma =  prev_Sigma + (7 * dSigma1 + 10 * dSigma2 + dSigma4) / 27;

                //Delta 5
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage4);
                dSigma5  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl5 = dLambda * m;
                iv_storage5.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, this](auto & iv5)
                {
                    using VT = std::decay_t<decltype(iv5)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;

                    auto h5 = iv5.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH5 = dLambda * h5;

                    iv5.trial_value = iv5.committed_value +  (28 * DH1 - 125 * DH2 + 546 * DH3 + 54 * DH4 - 378 * DH5) / 625;
                });
                next_Sigma =  prev_Sigma + (28 * dSigma1 - 125 * dSigma2 + 546 * dSigma3 + 54 * dSigma4 - 378 * dSigma5) / 625;

                //Delta 6 - The final predictor
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage5);
                dSigma6  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl6 = dLambda * m;
                next_iv_storage.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, &iv_storage5, this](auto & niv)
                {
                    using VT = std::decay_t<decltype(niv)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;
                    const VT &iv5 = iv_storage5.template get<VT>();
                    auto DH5 = iv5.trial_value - iv5.committed_value;

                    auto h6 = niv.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH6 = dLambda * h6;

                    niv.trial_value = niv.committed_value +  ( DH1 +  4 * DH3 + DH4) / 6;
                });
                dSigma =  (dSigma1 + 4 * dSigma3 + dSigma4) / 6;
                dEpsilonPl =  (dEpsilonPl1 + 4 * dEpsilonPl3 + dEpsilonPl4) / 6;
                next_Sigma =  prev_Sigma + dSigma;
                next_EpsilonPl = prev_EpsilonPl + dEpsilonPl;


                if (dSigma != dSigma)
                {
                    cout << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control Integration error" << endl;
                    cout << "T = " << T << endl;
                    cout << "dT = " << dT << endl;
                    cout << "dSigma1 = " << dSigma1.transpose() << endl;
                    cout << "dSigma2 = " << dSigma2.transpose() << endl;
                    cout << "dSigma3 = " << dSigma3.transpose() << endl;
                    cout << "dSigma4 = " << dSigma4.transpose() << endl;
                    cout << "dSigma5 = " << dSigma5.transpose() << endl;
                    cout << "dSigma6 = " << dSigma6.transpose() << endl;
                    cout << "m = " << m.transpose() << endl;
                    cout << "dEpsilonPl = " << dEpsilonPl.transpose() << endl;
                    cout << "TrialStress = " << TrialStress.transpose() << endl;
                    cout << "dEPS = " << dEPS.transpose() << endl;
                    exit(-1);
                }




                //Stress norm and stress error
                double stressNorm = next_Sigma.norm();
                double curStepError1 = (-42 * dSigma1 - 224 * dSigma3 - 21 * dSigma4 + 162 * dSigma5 + 125 * dSigma6).norm() / 336;
                if (stressNorm >= 0.5) { curStepError1 /= (2 * stressNorm); }
                //Internal variables norm and internal variables error (TODO)

                // double curStepError = curStepError1; //fmax(curStepError1, curStepError2);
                // double curStepError = curStepError1/stressNorm; //fmax(curStepError1, curStepError2);
                double curStepError = curStepError1; //fmax(curStepError1, curStepError2);

                //Check convergence and adjust integration timestep
                if (curStepError > TolE)
                {
                    double q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.1);

                    if (dT == dT_min) {

                        prev_Sigma = next_Sigma;
                        prev_EpsilonPl = next_EpsilonPl;
                        prev_iv_storage = next_iv_storage;

                        T += dT;
                    }
                    dT = fmax(q * dT, dT_min);
                }
                else {

                    prev_Sigma = next_Sigma;
                    prev_EpsilonPl = next_EpsilonPl;
                    prev_iv_storage = next_iv_storage;

                    double q = fmin(0.8 * pow(TolE / curStepError, 0.2), 2.0);
                    T += dT;
                    dT = fmax(q * dT, dT_min);
                    dT = fmin(dT, 1 - T);

                    maxStepError = max(maxStepError, curStepError);
                }

                if (niter > this->INT_OPT_RK45_niter_max[ASDP_TAG])
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " exceeded number of iterations. niter = " << niter << " niter_max = " <<this->INT_OPT_RK45_niter_max[ASDP_TAG] << " T= " << T << " dT = " << dT << endl;
                    // throw std::runtime_error("ASDPLasticMaterial3D - Unable to find a valid bracket in compute_yf_crossing");
                    return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
                }

            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = next_Sigma;
            TrialPlastic_Strain = next_EpsilonPl;
            iv_storage = next_iv_storage;

           //Return to Yield
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double yf_val_after_corrector;
                int iter =0;
                // double TOL = 10*this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
                // double NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                // do
                {
                    yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    // In the function below, depsilon_elpl is actually not used at all in hardening
                    // double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double dLambda_after_corrector = yf_val_after_corrector / (
                                                         n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector
                                                     );
                    TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                    TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;

                    // iter++;
                }
                // while(yf_val_after_corrector > TOL && iter < NITER);
            }


           //Return to Yield with bisection
           else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with iterations
           {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double y0  = yf(TrialStress, iv_storage, parameters_storage) ;
                int iter = 0;
                double TOL = this->yf_tolerance();
                double NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                // do
                if(y0 > 0 && iter < NITER)
                {
                    y0 = yf(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    // In the function below, depsilon_elpl is actually not used at all in hardening
                    // double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double dL = y0 / (
                                                         n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector
                                                     );

                    
                    VoigtVector TS = TrialStress - dL * Eelastic * m_after_corrector;

                    double y1 = yf(TS, iv_storage, parameters_storage);

                    //Try to bracket solution
                    // cout << "   y0 = " << y0 << "  dL =" << dL << "   y1 =" << y1 << endl;
                    while( y1 > 0 && iter < NITER)
                    {
                        dL = dL * 1.1;
                        TS = TrialStress - dL * Eelastic * m_after_corrector;
                        y1 = yf(TS, iv_storage, parameters_storage);
                        iter ++;

                        // cout << "   iter = " << iter << "  dL =" << dL << "   y1 =" << y1 << endl;
                    }

                    iter = 0;

                    //Once solution is bracketed, use bisection to get to YS
                    if (y1 < 0)
                    {
                        double dL_min = 0;
                        double dL_max = dL;
                        double dL_mid = dL / 2;

                        VoigtVector TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                        double y_mid = yf(TS2, iv_storage, parameters_storage);

                        // cout << "   y_mid = " << y_mid << "  dL_min =" << dL_min << "   dL_mid =" << dL_mid << "   dL_max =" << dL_max << endl;
                        
                        while(abs(y_mid) > TOL && iter < NITER)
                        {
                            // cout << "   iter = " << iter << "   y_mid = " << y_mid << "  dL_min =" << dL_min << "   dL_mid =" << dL_mid << "   dL_max =" << dL_max << endl;
                            if (y_mid  > 0)
                            {
                                dL_min = dL_mid;
                            } else
                            {
                                dL_max = dL_mid;
                            }
                            dL_mid = 0.5*(dL_min + dL_max);
                            TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                            y_mid = yf(TS2, iv_storage, parameters_storage);
                            iter++;
                        }
                        dL = dL_mid;
                    }
                    

                    TrialStress = TrialStress - dL * Eelastic * m_after_corrector;
                    TrialPlastic_Strain += dL * m_after_corrector;

                    // iter++;
                }
            }

            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control  Numeric error!\n";
                printTensor1("TrialStress = " , TrialStress);
                printTensor1("CommitStress = " , CommitStress);
                printTensor1("depsilon = " , depsilon);
                printTensor1("dsigma   = " , dsigma);
                printTensor1("intersection_stress = " , intersection_stress);
                printTensor2("Eelastic = " , Eelastic);
                printTensor2("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("n = " , yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage) );
                printTensor1("m = " , pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage) );
                cout << "hardening  = " << yf.hardening( depsilon_elpl, pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage),  TrialStress, iv_storage, parameters_storage) << endl;

                errorcode = -1;

                exit(-1);
            }

            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Runge_Kutta_45_Error_Control_old (post return-to-yield)", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            ComputeTangentStiffness();
        }

        return 0;
    }

    // Modified Euler with Error Control - Following RK45 coding style
    int Modified_Euler_Error_Control(const VoigtVector & strain_incr)
    {
        // Modified Euler coefficients (Heun's method)
        constexpr double a21 = 1.0;  // For predictor step
        constexpr double b1 = 0.5;   // Weight for corrector (average)
        constexpr double b2 = 0.5;   // Weight for corrector (average)

        int errorcode = -1;

        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;

        iv_storage.revert_all();

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
        dsigma = Eelastic * depsilon;

        TrialStress = CommitStress + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(CommitStress, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Modified_Euler_Error_Control", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
            // Adaptive step control parameters
            const double safety_factor = 0.9;
            const double min_scale = 0.25;
            const double max_scale = 4.0;
            const double beta = 0.04;  // PI controller parameter
            double previous_error = TolE;

            VoigtVector current_Sigma = TrialStress;
            iv_storage_t current_iv_storage = iv_storage;
            VoigtVector current_EpsilonPl = CommitPlastic_Strain;

            // Storage for k values - stress derivatives
            VoigtVector k1_sigma, k2_sigma;
            // Storage for k values - plastic strain derivatives  
            VoigtVector k1_pstrain, k2_pstrain;
            // Storage for intermediate iv storages
            iv_storage_t iv_k1 = iv_storage, iv_k2 = iv_storage;

            // Storage for predictor and corrector solutions
            VoigtVector predictor_sigma, corrector_sigma;
            VoigtVector predictor_pstrain, corrector_pstrain;
            iv_storage_t predictor_iv = iv_storage, corrector_iv = iv_storage;

            int niter = 0;
            double maxStepError = 0;
            int max_iterations = this->INT_OPT_RK45_niter_max[ASDP_TAG];
            
            while (T < 1.0)
            {
                niter++;
                
                double effective_dT = std::min(dT, 1.0 - T);
                VoigtVector dEPS = effective_dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                // Update elasticity matrix for current state
                Eelastic = et(current_Sigma, parameters_storage);

                // PREDICTOR STEP (Forward Euler)
                // k1 = f(t, y)
                std::tie(dLambda, m) = CalculateLambdaM(current_Sigma, dEPS, parameters_storage, current_iv_storage);
                k1_sigma = Eelastic * (dEPS - dLambda * m);
                k1_pstrain = dLambda * m;
                iv_k1 = current_iv_storage;
                iv_k1.apply([&m, &dLambda, &current_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, current_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + dLambda * h;
                });

                // Predictor solution: y_pred = y_n + h*k1
                predictor_sigma = current_Sigma + k1_sigma;
                predictor_pstrain = current_EpsilonPl + k1_pstrain;
                predictor_iv = current_iv_storage;
                predictor_iv.apply([&iv_k1, &current_iv_storage](auto & pred_var)
                {
                    using VT = std::decay_t<decltype(pred_var)>;
                    const VT &iv1_var = iv_k1.template get<VT>();
                    const VT &current_var = current_iv_storage.template get<VT>();
                    auto dk1 = iv1_var.trial_value - current_var.committed_value;
                    pred_var.trial_value = current_var.committed_value + dk1;
                });

                // CORRECTOR STEP (Modified Euler)
                // k2 = f(t + h, y_pred)
                Eelastic = et(predictor_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(predictor_sigma, dEPS, parameters_storage, predictor_iv);
                k2_sigma = Eelastic * (dEPS - dLambda * m);
                k2_pstrain = dLambda * m;
                iv_k2 = current_iv_storage;
                iv_k2.apply([&m, &dLambda, &predictor_sigma, &dEPS, this](auto & iv2)
                {
                    auto h = iv2.hardening_function(dEPS, m, predictor_sigma, parameters_storage);
                    iv2.trial_value = iv2.committed_value + dLambda * h;
                });

                // Corrector solution: y_corr = y_n + h/2*(k1 + k2)
                corrector_sigma = current_Sigma + b1 * k1_sigma + b2 * k2_sigma;
                corrector_pstrain = current_EpsilonPl + b1 * k1_pstrain + b2 * k2_pstrain;
                corrector_iv = current_iv_storage;
                corrector_iv.apply([&iv_k1, &iv_k2, &current_iv_storage, b1, b2](auto & corr_var)
                {
                    using VT = std::decay_t<decltype(corr_var)>;
                    const VT &iv1_var = iv_k1.template get<VT>();
                    const VT &iv2_var = iv_k2.template get<VT>();
                    const VT &current_var = current_iv_storage.template get<VT>();
                    auto dk1 = iv1_var.trial_value - current_var.committed_value;
                    auto dk2 = iv2_var.trial_value - current_var.committed_value;
                    corr_var.trial_value = current_var.committed_value + b1 * dk1 + b2 * dk2;
                });

                // Error estimation: difference between predictor and corrector
                VoigtVector sigma_error = corrector_sigma - predictor_sigma;
                VoigtVector pstrain_error = corrector_pstrain - predictor_pstrain;
                
                double step_error = sigma_error.norm() + pstrain_error.norm();
                
                // Normalize error by solution magnitude
                double solution_norm = corrector_sigma.norm() + corrector_pstrain.norm();
                if (solution_norm > 0.1) {
                    step_error /= solution_norm;
                }

                // Check for NaN
                if (std::isnan(step_error) || std::isnan(corrector_sigma.norm()) || std::isnan(corrector_pstrain.norm()))
                {
                    cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control - NaN encountered, reducing step size" << endl;
                    dT *= 0.5;
                    if (dT < dT_min) {
                        cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control - Minimum step size reached with NaN" << endl;
                        return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
                    }
                    continue;
                }

                // Step size control with PI controller
                double error_ratio = TolE / std::max(step_error, 1e-15);
                double scale_factor = safety_factor * std::pow(error_ratio, 0.5) * std::pow(previous_error / step_error, beta);
                scale_factor = std::max(min_scale, std::min(max_scale, scale_factor));

                // Accept or reject step
                if (step_error <= TolE || effective_dT <= dT_min) {
                    // Ladruno (ADR-94 wp/94a): ADR-94 M2(c) -- at dT_min the error control degenerates to
                    // accept-everything. Under strict_convergence an out-of-tolerance
                    // substep is a refusal, not an acceptance.
                    if (step_error > TolE && INT_OPT_strict_convergence[ASDP_TAG] != 0) {
                        opserr << "ASDPlasticMaterial3D::Modified_Euler_Error_Control (tag " << ASDP_TAG
                               << ") - substep accepted at dT_min with step_error = " << step_error
                               << " > stress_absolute_tol = " << TolE
                               << " -- rejecting step (strict_convergence)" << endln;
                        return LADRUNO_MATERIAL_REFUSED;
                    }
                    // Accept step - use corrector solution
                    current_Sigma = corrector_sigma;
                    current_EpsilonPl = corrector_pstrain;
                    current_iv_storage = corrector_iv;
                    
                    T += effective_dT;
                    maxStepError = std::max(maxStepError, step_error);
                    previous_error = step_error;
                    
                    // Validate yield function drift
                    double yf_val = yf(current_Sigma, current_iv_storage, parameters_storage);
                    if (yf_val > 10 * yf_tolerance()) {
                        // cout << "Warning: Yield function drift detected: f = " << yf_val << endl;
                    }
                }

                // Update step size for next iteration
                double new_dT = scale_factor * effective_dT;
                dT = std::max(dT_min, std::min(new_dT, 1.0 - T));

                if (niter > max_iterations)
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " Modified Euler exceeded number of iterations. niter = " << niter << " niter_max = " << max_iterations << " T= " << T << " dT = " << dT << endl;
                    return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
                }
            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = current_Sigma;
            TrialPlastic_Strain = current_EpsilonPl;
            iv_storage = current_iv_storage;

            //Return to Yield
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                if (yf_val_after_corrector > yf_tolerance()) {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dLambda_after_corrector = yf_val_after_corrector / denominator;
                        TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                        TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dLambda_after_corrector, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dLambda_after_corrector * h;
                        });
                    }
                }
            }
            else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with bisection
            {
                double y0 = yf(TrialStress, iv_storage, parameters_storage);
                int iter = 0;
                double TOL = this->yf_tolerance();
                int NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                
                if(y0 > TOL && iter < NITER)
                {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dL = y0 / denominator;
                        VoigtVector TS = TrialStress - dL * Eelastic * m_after_corrector;
                        double y1 = yf(TS, iv_storage, parameters_storage);

                        // Try to bracket solution
                        while( y1 > 0 && iter < NITER)
                        {
                            dL = dL * 1.1;
                            TS = TrialStress - dL * Eelastic * m_after_corrector;
                            y1 = yf(TS, iv_storage, parameters_storage);
                            iter++;
                        }

                        iter = 0;

                        // Once solution is bracketed, use bisection to get to YS
                        if (y1 < 0)
                        {
                            double dL_min = 0;
                            double dL_max = dL;
                            double dL_mid = dL / 2;

                            VoigtVector TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                            double y_mid = yf(TS2, iv_storage, parameters_storage);
                            
                            while(std::abs(y_mid) > TOL && iter < NITER)
                            {
                                if (y_mid > 0) {
                                    dL_min = dL_mid;
                                } else {
                                    dL_max = dL_mid;
                                }
                                dL_mid = 0.5*(dL_min + dL_max);
                                TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                                y_mid = yf(TS2, iv_storage, parameters_storage);
                                iter++;
                            }
                            dL = dL_mid;
                        }
                        
                        TrialStress = TrialStress - dL * Eelastic * m_after_corrector;
                        TrialPlastic_Strain += dL * m_after_corrector;
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                    }
                }
            }

            // Final validation
            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control  Numeric error!\n";
                printTensor1("TrialStress = " , TrialStress);
                printTensor1("CommitStress = " , CommitStress);
                printTensor1("depsilon = " , depsilon);
                printTensor1("dsigma   = " , dsigma);
                printTensor1("intersection_stress = " , intersection_stress);
                printTensor2("Eelastic = " , Eelastic);
                printTensor2("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("n = " , yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage) );
                printTensor1("m = " , pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage) );
                cout << "hardening  = " << yf.hardening( depsilon_elpl, pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage),  TrialStress, iv_storage, parameters_storage) << endl;

                return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }

            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Modified_Euler_Error_Control (post return-to-yield)", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            ComputeTangentStiffness();
        }

        return 0;
    }

    // Complete corrected RK45 implementation with proper internal variable handling
    int Runge_Kutta_45_Error_Control(const VoigtVector & strain_incr)
    {
        int RK45_EXIT_FLAG = 0;

        // Dormand-Prince RK45 coefficients
        constexpr double a21 = 1.0/5.0;
        constexpr double a31 = 3.0/40.0, a32 = 9.0/40.0;
        constexpr double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
        constexpr double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
        constexpr double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
        
        // 5th order solution coefficients
        constexpr double b1 = 35.0/384.0, b2 = 0.0, b3 = 500.0/1113.0, b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0;
        
        // 4th order solution coefficients for error estimation
        constexpr double bhat1 = 5179.0/57600.0, bhat2 = 0.0, bhat3 = 7571.0/16695.0, bhat4 = 393.0/640.0;
        constexpr double bhat5 = -92097.0/339200.0, bhat6 = 187.0/2100.0, bhat7 = 1.0/40.0;

        int errorcode = -1;

        VoigtVector depsilon;  // Ladruno (ADR-94 wp/94b, M1): was `static` -- a function-local buffer shared across every instance of this specialization
        depsilon.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        depsilon = strain_incr;

        // CRITICAL FIX: Initialize trial values from committed values at start
        iv_storage.apply([](auto & iv) {
            iv.trial_value = iv.committed_value;
        });

        dsigma.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_stress.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN
        intersection_strain.setZero();  // Ladruno (ADR-94 wp/94a): *= 0 does not clear NaN

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
        dsigma = Eelastic * depsilon;

        TrialStress = CommitStress + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(CommitStress, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Runge_Kutta_45_Error_Control", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = yf_tolerance();
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
            // Adaptive step control parameters
            const double safety_factor = 0.9;
            const double min_scale = 0.2;
            const double max_scale = 5.0;
            const double beta = 0.04;  // PI controller parameter
            double previous_error = TolE;

            VoigtVector current_Sigma = TrialStress;
            iv_storage_t current_iv_storage = iv_storage;
            VoigtVector current_EpsilonPl = CommitPlastic_Strain;

            // Storage for k values - stress derivatives
            VoigtVector k1_sigma, k2_sigma, k3_sigma, k4_sigma, k5_sigma, k6_sigma;
            // Storage for k values - plastic strain derivatives  
            VoigtVector k1_pstrain, k2_pstrain, k3_pstrain, k4_pstrain, k5_pstrain, k6_pstrain;
            
            // CORRECTED: Storage for internal variable derivatives (not absolute states)
            iv_storage_t iv_k1_derivatives, iv_k2_derivatives, iv_k3_derivatives;
            iv_storage_t iv_k4_derivatives, iv_k5_derivatives, iv_k6_derivatives;

            int niter = 0;
            double maxStepError = 0;
            
            while (T < 1.0)
            {
                niter++;
                
                double effective_dT = std::min(dT, 1.0 - T);
                VoigtVector dEPS = effective_dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                // Update elasticity matrix for current state
                Eelastic = et(current_Sigma, parameters_storage);

                // k1 = f(t, y) - compute derivatives at current state
                std::tie(dLambda, m) = CalculateLambdaM(current_Sigma, dEPS, parameters_storage, current_iv_storage);
                k1_sigma = Eelastic * (dEPS - dLambda * m);
                k1_pstrain = dLambda * m;
                
                // CORRECTED: Store pure derivatives for internal variables
                iv_k1_derivatives = current_iv_storage;  // Copy structure
                iv_k1_derivatives.apply([&m, &dLambda, &current_Sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, current_Sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;  // Pure derivative
                });

                // k2 = f(t + c2*h, y + h*(a21*k1))
                VoigtVector y2_sigma = current_Sigma + a21 * k1_sigma;
                VoigtVector y2_pstrain = current_EpsilonPl + a21 * k1_pstrain;
                
                // CORRECTED: Properly propagate internal variable state using RK formula
                iv_storage_t iv2_state = current_iv_storage;
                iv2_state.apply([&iv_k1_derivatives, &current_iv_storage, a21](auto & iv2_var)
                {
                    using VT = std::decay_t<decltype(iv2_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    iv2_var.trial_value = current_var.trial_value + a21 * k1_deriv.trial_value;
                });
                
                Eelastic = et(y2_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y2_sigma, dEPS, parameters_storage, iv2_state);
                k2_sigma = Eelastic * (dEPS - dLambda * m);
                k2_pstrain = dLambda * m;
                
                iv_k2_derivatives = current_iv_storage;
                iv_k2_derivatives.apply([&m, &dLambda, &y2_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y2_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k3 = f(t + c3*h, y + h*(a31*k1 + a32*k2))
                VoigtVector y3_sigma = current_Sigma + a31 * k1_sigma + a32 * k2_sigma;
                VoigtVector y3_pstrain = current_EpsilonPl + a31 * k1_pstrain + a32 * k2_pstrain;
                
                iv_storage_t iv3_state = current_iv_storage;
                iv3_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &current_iv_storage, a31, a32](auto & iv3_var)
                {
                    using VT = std::decay_t<decltype(iv3_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    iv3_var.trial_value = current_var.trial_value + a31 * k1_deriv.trial_value + a32 * k2_deriv.trial_value;
                });
                
                Eelastic = et(y3_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y3_sigma, dEPS, parameters_storage, iv3_state);
                k3_sigma = Eelastic * (dEPS - dLambda * m);
                k3_pstrain = dLambda * m;
                
                iv_k3_derivatives = current_iv_storage;
                iv_k3_derivatives.apply([&m, &dLambda, &y3_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y3_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k4 = f(t + c4*h, y + h*(a41*k1 + a42*k2 + a43*k3))
                VoigtVector y4_sigma = current_Sigma + a41 * k1_sigma + a42 * k2_sigma + a43 * k3_sigma;
                VoigtVector y4_pstrain = current_EpsilonPl + a41 * k1_pstrain + a42 * k2_pstrain + a43 * k3_pstrain;
                
                iv_storage_t iv4_state = current_iv_storage;
                iv4_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &current_iv_storage, 
                               a41, a42, a43](auto & iv4_var)
                {
                    using VT = std::decay_t<decltype(iv4_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    iv4_var.trial_value = current_var.trial_value + a41 * k1_deriv.trial_value + 
                                         a42 * k2_deriv.trial_value + a43 * k3_deriv.trial_value;
                });
                
                Eelastic = et(y4_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y4_sigma, dEPS, parameters_storage, iv4_state);
                k4_sigma = Eelastic * (dEPS - dLambda * m);
                k4_pstrain = dLambda * m;
                
                iv_k4_derivatives = current_iv_storage;
                iv_k4_derivatives.apply([&m, &dLambda, &y4_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y4_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k5 = f(t + c5*h, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
                VoigtVector y5_sigma = current_Sigma + a51 * k1_sigma + a52 * k2_sigma + a53 * k3_sigma + a54 * k4_sigma;
                VoigtVector y5_pstrain = current_EpsilonPl + a51 * k1_pstrain + a52 * k2_pstrain + 
                                        a53 * k3_pstrain + a54 * k4_pstrain;
                
                iv_storage_t iv5_state = current_iv_storage;
                iv5_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                               &current_iv_storage, a51, a52, a53, a54](auto & iv5_var)
                {
                    using VT = std::decay_t<decltype(iv5_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    iv5_var.trial_value = current_var.trial_value + a51 * k1_deriv.trial_value + 
                                         a52 * k2_deriv.trial_value + a53 * k3_deriv.trial_value + 
                                         a54 * k4_deriv.trial_value;
                });
                
                Eelastic = et(y5_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y5_sigma, dEPS, parameters_storage, iv5_state);
                k5_sigma = Eelastic * (dEPS - dLambda * m);
                k5_pstrain = dLambda * m;
                
                iv_k5_derivatives = current_iv_storage;
                iv_k5_derivatives.apply([&m, &dLambda, &y5_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y5_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k6 = f(t + h, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
                VoigtVector y6_sigma = current_Sigma + a61 * k1_sigma + a62 * k2_sigma + a63 * k3_sigma + 
                                      a64 * k4_sigma + a65 * k5_sigma;
                VoigtVector y6_pstrain = current_EpsilonPl + a61 * k1_pstrain + a62 * k2_pstrain + 
                                        a63 * k3_pstrain + a64 * k4_pstrain + a65 * k5_pstrain;
                
                iv_storage_t iv6_state = current_iv_storage;
                iv6_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                               &iv_k5_derivatives, &current_iv_storage, a61, a62, a63, a64, a65](auto & iv6_var)
                {
                    using VT = std::decay_t<decltype(iv6_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    iv6_var.trial_value = current_var.trial_value + a61 * k1_deriv.trial_value + 
                                         a62 * k2_deriv.trial_value + a63 * k3_deriv.trial_value + 
                                         a64 * k4_deriv.trial_value + a65 * k5_deriv.trial_value;
                });
                
                Eelastic = et(y6_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y6_sigma, dEPS, parameters_storage, iv6_state);
                k6_sigma = Eelastic * (dEPS - dLambda * m);
                k6_pstrain = dLambda * m;
                
                iv_k6_derivatives = current_iv_storage;
                iv_k6_derivatives.apply([&m, &dLambda, &y6_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y6_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // 5th order solution
                VoigtVector next_Sigma_5th = current_Sigma + (b1 * k1_sigma + b2 * k2_sigma + b3 * k3_sigma + 
                                                             b4 * k4_sigma + b5 * k5_sigma + b6 * k6_sigma);
                VoigtVector next_EpsilonPl_5th = current_EpsilonPl + (b1 * k1_pstrain + b2 * k2_pstrain + 
                                                                     b3 * k3_pstrain + b4 * k4_pstrain + 
                                                                     b5 * k5_pstrain + b6 * k6_pstrain);

                // CORRECTED: 5th order solution for internal variables including ALL terms
                iv_storage_t next_iv_5th = current_iv_storage;
                next_iv_5th.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                                 &iv_k5_derivatives, &iv_k6_derivatives, &current_iv_storage, 
                                 b1, b2, b3, b4, b5, b6](auto & next_var)
                {
                    using VT = std::decay_t<decltype(next_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    const VT &k6_deriv = iv_k6_derivatives.template get<VT>();
                    
                    next_var.trial_value = current_var.trial_value + 
                        b1 * k1_deriv.trial_value + b2 * k2_deriv.trial_value + b3 * k3_deriv.trial_value + 
                        b4 * k4_deriv.trial_value + b5 * k5_deriv.trial_value + b6 * k6_deriv.trial_value;
                });

                // 4th order solution for error estimation
                VoigtVector next_Sigma_4th = current_Sigma + (bhat1 * k1_sigma + bhat2 * k2_sigma + bhat3 * k3_sigma + 
                                                             bhat4 * k4_sigma + bhat5 * k5_sigma + bhat6 * k6_sigma);
                VoigtVector next_EpsilonPl_4th = current_EpsilonPl + (bhat1 * k1_pstrain + bhat2 * k2_pstrain + 
                                                                     bhat3 * k3_pstrain + bhat4 * k4_pstrain + 
                                                                     bhat5 * k5_pstrain + bhat6 * k6_pstrain);

                // 4th order solution for internal variables (for error estimation)
                iv_storage_t next_iv_4th = current_iv_storage;
                next_iv_4th.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                                 &iv_k5_derivatives, &iv_k6_derivatives, &current_iv_storage, 
                                 bhat1, bhat2, bhat3, bhat4, bhat5, bhat6](auto & next_var_4th)
                {
                    using VT = std::decay_t<decltype(next_var_4th)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    const VT &k6_deriv = iv_k6_derivatives.template get<VT>();
                    
                    next_var_4th.trial_value = current_var.trial_value + 
                        bhat1 * k1_deriv.trial_value + bhat2 * k2_deriv.trial_value + bhat3 * k3_deriv.trial_value + 
                        bhat4 * k4_deriv.trial_value + bhat5 * k5_deriv.trial_value + bhat6 * k6_deriv.trial_value;
                });

                // Error estimation including internal variables
                VoigtVector sigma_error = next_Sigma_5th - next_Sigma_4th;
                VoigtVector pstrain_error = next_EpsilonPl_5th - next_EpsilonPl_4th;
                
                // ADDED: Include internal variable errors in step control
                double iv_error = 0.0;
                next_iv_5th.apply([&next_iv_4th, &iv_error](const auto & iv_5th)
                {
                    using VT = std::decay_t<decltype(iv_5th)>;
                    const VT &iv_4th = next_iv_4th.template get<VT>();
                    auto diff = iv_5th.trial_value - iv_4th.trial_value;
                    iv_error += diff.norm();
                });
                
                double step_error = sigma_error.norm() + pstrain_error.norm() + iv_error;
                
                // Normalize error by solution magnitude
                double solution_norm = next_Sigma_5th.norm() + next_EpsilonPl_5th.norm();
                double iv_norm = 0.0;
                next_iv_5th.apply([&iv_norm](const auto & iv)
                {
                    iv_norm += iv.trial_value.norm();
                });
                solution_norm += iv_norm;
                
                if (solution_norm > 0.1) {
                    step_error /= solution_norm;
                }

                // Check for NaN
                if (std::isnan(step_error) || std::isnan(next_Sigma_5th.norm()) || std::isnan(next_EpsilonPl_5th.norm()))
                {
                    cout << "ASDPlasticMaterial3D::RK45 - NaN encountered, reducing step size" << endl;
                    dT *= 0.5;
                    if (dT < dT_min) {
                        cout << "ASDPlasticMaterial3D::RK45 - Minimum step size reached with NaN" << endl;
                        return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
                    }
                    continue;
                }

                // Step size control with PI controller
                double error_ratio = TolE / std::max(step_error, 1e-15);
                double scale_factor = safety_factor * std::pow(error_ratio, 0.2) * std::pow(previous_error / step_error, beta);
                scale_factor = std::max(min_scale, std::min(max_scale, scale_factor));

                // Accept or reject step
                if (step_error <= TolE || effective_dT <= dT_min) {
                    // Ladruno (ADR-94 wp/94a): ADR-94 M2(c) -- at dT_min the error control degenerates to
                    // accept-everything. Under strict_convergence an out-of-tolerance
                    // substep is a refusal, not an acceptance.
                    if (step_error > TolE && INT_OPT_strict_convergence[ASDP_TAG] != 0) {
                        opserr << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control (tag " << ASDP_TAG
                               << ") - substep accepted at dT_min with step_error = " << step_error
                               << " > stress_absolute_tol = " << TolE
                               << " -- rejecting step (strict_convergence)" << endln;
                        return LADRUNO_MATERIAL_REFUSED;
                    }
                    // Accept step - use 5th order solution
                    current_Sigma = next_Sigma_5th;
                    current_EpsilonPl = next_EpsilonPl_5th;
                    current_iv_storage = next_iv_5th;
                    
                    T += effective_dT;
                    maxStepError = std::max(maxStepError, step_error);
                    previous_error = step_error;
                    
                    // Validate yield function drift
                    double yf_val = yf(current_Sigma, current_iv_storage, parameters_storage);
                    if (yf_val > 10 * TolE) {
                        // Optional: Add drift correction here
                        // cout << "Warning: Yield function drift detected: f = " << yf_val << endl;
                    }
                }

                // Update step size for next iteration
                double new_dT = scale_factor * effective_dT;
                dT = std::max(dT_min, std::min(new_dT, 1.0 - T));

                if (niter > this->INT_OPT_RK45_niter_max[ASDP_TAG])
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " RK45 exceeded number of iterations. niter = " << niter << " niter_max = " <<this->INT_OPT_RK45_niter_max[ASDP_TAG] << " T= " << T << " dT = " << dT << endl;
                    RK45_EXIT_FLAG = -1;
                    break;
                }
            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = current_Sigma;
            TrialPlastic_Strain = current_EpsilonPl;
            // iv_storage = current_iv_storage;
            iv_storage.updateTrialValueFromOther(current_iv_storage);


            //Return to Yield Surface
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                if (yf_val_after_corrector > yf_tolerance()) {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dLambda_after_corrector = yf_val_after_corrector / denominator;
                        TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                        TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dLambda_after_corrector, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dLambda_after_corrector * h;
                        });
                    }
                }
            }
            else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with bisection
            {
                double y0 = yf(TrialStress, iv_storage, parameters_storage);
                int iter = 0;
                double TOL = this->yf_tolerance();
                int NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                
                if(y0 > TOL && iter < NITER)
                {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dL = y0 / denominator;
                        VoigtVector TS = TrialStress - dL * Eelastic * m_after_corrector;
                        
                        // Create temporary iv storage for testing
                        iv_storage_t temp_iv = iv_storage;
                        temp_iv.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                        
                        double y1 = yf(TS, temp_iv, parameters_storage);

                        // Try to bracket solution
                        while( y1 > 0 && iter < NITER)
                        {
                            dL = dL * 1.1;
                            TS = TrialStress - dL * Eelastic * m_after_corrector;
                            
                            temp_iv = iv_storage;
                            temp_iv.apply([&m_after_corrector, &dL, this](auto& iv) {
                                auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                iv.trial_value += dL * h;
                            });
                            
                            y1 = yf(TS, temp_iv, parameters_storage);
                            iter++;
                        }

                        iter = 0;

                        // Once solution is bracketed, use bisection to get to YS
                        if (y1 < 0)
                        {
                            double dL_min = 0;
                            double dL_max = dL;
                            double dL_mid = dL / 2;

                            VoigtVector TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                            
                            iv_storage_t mid_iv = iv_storage;
                            mid_iv.apply([&m_after_corrector, &dL_mid, this](auto& iv) {
                                auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                iv.trial_value += dL_mid * h;
                            });
                            
                            double y_mid = yf(TS2, mid_iv, parameters_storage);
                            
                            while(std::abs(y_mid) > TOL && iter < NITER)
                            {
                                if (y_mid > 0) {
                                    dL_min = dL_mid;
                                } else {
                                    dL_max = dL_mid;
                                }
                                dL_mid = 0.5*(dL_min + dL_max);
                                TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                                
                                mid_iv = iv_storage;
                                mid_iv.apply([&m_after_corrector, &dL_mid, this](auto& iv) {
                                    auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                    iv.trial_value += dL_mid * h;
                                });
                                
                                y_mid = yf(TS2, mid_iv, parameters_storage);
                                iter++;
                            }
                            dL = dL_mid;
                        }
                        
                        TrialStress = TrialStress - dL * Eelastic * m_after_corrector;
                        TrialPlastic_Strain += dL * m_after_corrector;
                        
                        // Update internal variables with final correction
                        iv_storage.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                    }
                }
            }

            // Final validation
            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control  Numeric error!\n";
                printTensor1("TrialStress = " , TrialStress);
                printTensor1("CommitStress = " , CommitStress);
                printTensor1("depsilon = " , depsilon);
                printTensor1("dsigma   = " , dsigma);
                printTensor1("intersection_stress = " , intersection_stress);
                printTensor2("Eelastic = " , Eelastic);
                printTensor2("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("n = " , yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage) );
                printTensor1("m = " , pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage) );
                cout << "hardening  = " << yf.hardening( depsilon_elpl, pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage),  TrialStress, iv_storage, parameters_storage) << endl;

                return LADRUNO_MATERIAL_REFUSED;  // Ladruno (ADR-94 wp/94a): was a bare -1, dropped by every hex host
            }

            // ADDED: Final consistency check for internal variables
            #ifdef DEBUG_INTERNAL_VARIABLES
            cout << "=== Final Internal Variable State ===" << endl;
            iv_storage.apply([](const auto & iv) {
                cout << iv.getName() << ": trial=" << iv.trial_value.transpose() 
                     << ", committed=" << iv.committed_value.transpose() << endl;
            });
            
            double final_yf = yf(TrialStress, iv_storage, parameters_storage);
            cout << "Final yield function value: " << final_yf << endl;
            #endif

            // Ladruno (ADR-94 wp/94a)
            if (ladruno_strict_rejects("Runge_Kutta_45_Error_Control (post return-to-yield)", TrialStress))
                return LADRUNO_MATERIAL_REFUSED;
            ComputeTangentStiffness();
        }

        // Ladruno (ADR-94 wp/94a): RK45_EXIT_FLAG is set to -1 when the substep
        // iteration budget is exhausted; that bare -1 was dropped by every hex host.
        return RK45_EXIT_FLAG == 0 ? 0 : LADRUNO_MATERIAL_REFUSED;
    }


    //Robust Brent algorithm
    double compute_yf_crossing(const VoigtVector & start_stress, const VoigtVector & end_stress, double x1, double x2, double tol) const
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // Constants
        const double MAX_BRACKET_EXPANSION = 100.0; // Maximum factor to expand the bracket
        const double BRACKET_EXPANSION_FACTOR = 1.6; // Factor to expand the bracket each iteration
        const double EPS_MULTIPLE = 2.0;
        const double TOL_MULTIPLE = 0.5;
        const double EPS = std::numeric_limits<double>::epsilon();

        // Lambda functions
        auto calculateSigma = [](const VoigtVector & start, const VoigtVector & end, double multiplier) {
            return start * (1 - multiplier) + end * multiplier;
        };

        auto calculateYf = [this](const VoigtVector & sigma) {
            return yf(sigma, iv_storage, parameters_storage);
        };

        // Bracket Adjustment: Expand the bracket if necessary
        double a = x1;
        double b = x2;
        double c = x2;
        double d = 0;
        double e = 0.0;
        double fa = calculateYf(calculateSigma(start_stress, end_stress, a));
        double fb = calculateYf(calculateSigma(start_stress, end_stress, b));
        double fc = fb;


        for (double factor = 1.0; factor <= MAX_BRACKET_EXPANSION; factor *= BRACKET_EXPANSION_FACTOR) {
            if ((fb * fa) <= 0.0) {
                break; // Valid bracket found
            }
            if (std::abs(fa) < std::abs(fb)) {
                a -= (b - a) * factor;
                fa = calculateYf(calculateSigma(start_stress, end_stress, a));
            } else {
                b += (b - a) * factor;
                fb = calculateYf(calculateSigma(start_stress, end_stress, b));
            }
        }

        // Brent's Method Implementation
        if ((fb * fa) > 0.0) {
            throw std::runtime_error("ASDPLasticMaterial3D - Unable to find a valid bracket in compute_yf_crossing");
        }

        for (int iter = 1; iter <= ASDPlasticMaterial3D_MAXITER_BRENT; iter++) {
            if ((fb * fc) > 0.0) {
                c = a;   // Rename a, b, c and adjust bounding interval d
                fc = fa;
                e = d = b - a;
            }

            if (std::abs(fc) < std::abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            double tol1 = EPS_MULTIPLE * EPS * std::abs(b) + TOL_MULTIPLE * tol;
            double xm = 0.5 * (c - b);

            if (std::abs(xm) <= tol1 || fb == 0.0) {
                return b;  // Convergence
            }

            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
                // Attempt inverse quadratic interpolation
                double s = fb / fa;
                double p, q;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) q = -q;
                p = std::abs(p);
                double min1 = 3.0 * xm * q - std::abs(tol1 * q);
                double min2 = std::abs(e * q);

                if (2.0 * p < std::min(min1, min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;  // Bounds decreasing too slowly, use bisection
                e = d;
            }

            a = b;
            fa = fb;

            if (std::abs(d) > tol1) {
                b += d;
            } else {
                b += (xm > 0.0 ? std::abs(tol1) : -std::abs(tol1));
            }

            fb = calculateYf(calculateSigma(start_stress, end_stress, b));

            // Check for NaN (Not a Number) values
            if (std::isnan(fb)) {
                throw std::runtime_error("compute_yf_crossing: NaN encountered in function evaluation.");
            }
        }

        throw std::runtime_error("Maximum iterations reached without convergence in compute_yf_crossing");
    }


protected:

    VoigtVector TrialStrain;
    VoigtVector TrialStress;
    VoigtVector TrialPlastic_Strain;

    VoigtVector CommitStress;
    VoigtVector CommitStrain;
    VoigtVector CommitPlastic_Strain;


    YieldFunctionType yf;
    ElasticityType    et;
    PlasticFlowType   pf;

    iv_storage_t iv_storage;
    parameters_storage_t parameters_storage;

    // Ladruno (ADR-94 wp/94b, M6): snapshot of iv_storage as configured at model-build
    // time so revertToStart() can restore the internal variables' initial values.
    iv_storage_t iv_storage_initial;
    bool initial_iv_captured;

    std::string current_parameter_name; // Stores the most recent parameter name from setParameter

protected:

    static std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> INT_OPT_constitutive_integration_method;     //
    static std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> INT_OPT_tangent_operator_type;     //
    static std::map<int, double> DBL_OPT_f_absolute_tol;
    static std::map<int, double> DBL_OPT_stress_absolute_tol;
    static std::map<int, int> INT_OPT_n_max_iterations;
    static std::map<int, int> INT_OPT_return_to_yield_surface;
    static std::map<int, double> DBL_OPT_RK45_dT_min;
    static std::map<int, int> INT_OPT_RK45_niter_max;
    static std::map<int, int> INT_OPT_strict_convergence; // Ladruno (ADR-84 P2a): 1 = fail loud on BE non-convergence instead of silently accepting
    static std::map<int, double> DBL_OPT_f_relative_tol; // Ladruno (ADR-94 wp/94c, M5): 0 = off (absolute tolerance only)

    static std::map<int, int> GLOBAL_INT_max_iter; 
    static std::map<int, double> GLOBAL_DBL_max_error; 

    // ADR97_P4_MARKER:adr97_p4_suppress_numerical_tangent_member
    bool first_step;
    bool stress_set_externally;

    // Ladruno (ADR-97 wp/97e): recursion guard for
    // numerical_tangent_of_committed_map(). Backward_Euler and Closest_Point
    // both call ComputeTangentStiffness() at the end of every successful
    // commit; a perturbed sub-call made BY the numerical-tangent helper would
    // otherwise re-enter that same call and try to compute another numerical
    // tangent of ITS OWN perturbed state, unbounded. Per-instance (not one of
    // the static per-tag option maps): each material instance's own
    // in-flight tangent call must suppress recursion independently.
    bool suppress_numerical_tangent = false;

    // Ladruno (ADR-97 wp/97b): per-instance Closest_Point Newton iteration count
    // (NSDMI, so both constructors get it without touching their bodies).
    int cp_last_iterations = 0;

    // Ladruno (ADR-94 wp/94b, M1/F1): these five were `static` -- ONE copy shared by
    // every Gauss point, element and material tag of a given <E,Y,P,tag> specialization,
    // so the whole model was assembled with the tangent of the last GP integrated
    // (ADR-94 H1) and no threaded element loop (ADR-75b) could ever be deterministic.
    // The per-tag INT_OPT_*/GLOBAL_* maps above stay static: they are keyed by material
    // tag and shared by design.
    // The four scratch buffers are `mutable` because `compute_local_stress()` --
    // the const helper the numerical-tangent probe calls -- writes them. It used to
    // write the class-STATIC copies, i.e. it scribbled on every other instance's
    // scratch state; `mutable` keeps that behaviour byte-identical while confining
    // the damage to the probing instance. `Stiffness` is deliberately NOT mutable:
    // no const method may set the tangent.
    mutable VoigtVector dsigma;
    mutable VoigtVector depsilon_elpl;    //Elastoplastic strain increment : For a strain increment that causes first yield, the step is divided into an elastic one (until yield) and an elastoplastic one.
    mutable VoigtVector intersection_stress;
    mutable VoigtVector intersection_strain;
    VoigtMatrix Stiffness;


};

template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_constitutive_integration_method;
template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_tangent_operator_type;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_f_absolute_tol;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_stress_absolute_tol;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_n_max_iterations;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_return_to_yield_surface;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_RK45_dT_min;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_RK45_niter_max;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_strict_convergence; // Ladruno (ADR-84 P2a)
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_f_relative_tol; // Ladruno (ADR-94 wp/94c, M5)

template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::GLOBAL_DBL_max_error;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::GLOBAL_INT_max_iter; 

// Ladruno (ADR-94 wp/94b, M1/F1): the out-of-class definitions of dsigma,
// depsilon_elpl, intersection_stress, intersection_strain and Stiffness were here.
// They are ordinary per-instance members now -- see the declarations above.


#endif
