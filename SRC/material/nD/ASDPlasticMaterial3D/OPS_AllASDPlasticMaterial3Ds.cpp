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

#ifdef _EIGEN3

#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <list>
#include <set>       // Ladruno (ADR-94 wp/94a): required-parameter bookkeeping
#include <string>    // Ladruno (ADR-94 wp/94a)

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "AllASDPlasticMaterial3Ds.h"

// Ladruno (ADR-94 wp/94a): was void. ADR-94 B1/H13 -- an unrecognised token inside
// Begin_Integration_Options, an unknown model-parameter name, and a model parameter
// never given a value were ALL accepted silently (the last one running at 0.0, so a
// typo'd MC_phi ran the model at phi = 0). Parsing now FAILS and the material is not
// created. `asdp_parse_rejected` distinguishes "your options were wrong" from the
// pre-existing "no such YF/PF/EL/IV combination" report downstream.
template <typename T>
bool populate_ASDPlasticMaterial3D(T* instance);

static bool asdp_parse_rejected = false;   // Ladruno (ADR-94 wp/94a)

typedef std::tuple<std::string, std::string, std::string, std::string> model_spec_t;


NDMaterial*  ASDPlasticMaterial3DFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models);


template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial3D(int instance_tag, const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, std::list<NDMaterial*> &instance_pointers, std::list<model_spec_t> &available_models);


void print_usage(void)
{
    opserr <<
       "SYNTAX:\n"
       "nDMaterial ASDPlasticMaterial3D $tag \\ \n"
       "    YF_type \\ \n"
       "    PF_type \\ \n"
       "    EL_type \\ \n"
       "    IV_type \\ \n"
       "Begin_Internal_Variables \\ \n"
       "    IV1 (initial value(s)) \\ \n"
       "    IV2 (initial value(s)) \\ \n"
       "    ....\n"
       "End_Internal_Variables \\ \n"
       "Begin_Model_Parameters \\ \n"
       "    PARAM1 (value) \\\n"
       "    PARAM2 (value) \\\n"
       "    ....\n"
       "End_Model_Parameters \\ \n"
       "Begin_Integration_Options \\ \n"
       "    f_absolute_tol (double value)\\ \n"
       "    f_relative_tol (double value, default 0 = off) : tolerance becomes\n"
       "        max(f_absolute_tol, f_relative_tol * yield-function strength scale),\n"
       "        which makes convergence independent of the unit system\\ \n" // Ladruno (ADR-94 wp/94c, M5)
       "    stress_absolute_tol (double value)\\ \n"
       "    n_max_iterations (int value)\\ \n"
       "    return_to_yield_surface (0 or 1)\\ \n"
       "    strict_convergence (0 or 1, default 0) : 1 = Backward_Euler fails loud on non-convergence\\ \n" // Ladruno (ADR-84 P2a)
       "    integration_method (string) : Backward_Euler | Closest_Point (ADR-97) |\n"
       "        Forward_Euler | Forward_Euler_Subincrement |\n"
       "        Modified_Euler_Error_Control | Runge_Kutta_45_Error_Control\n"
       "    (Closest_Point is the fully implicit closest-point return map; it is the\n"
       "     ONLY integrator for which tangent_type Algorithmic -- the exact\n"
       "     consistent tangent -- is defined, and it is currently available for\n"
       "     VonMises and DruckerPrager only.)\n"
       "    (ADR-94: Backward_Euler_LineSearch and Runge_Kutta_45_Error_Control_old are REFUSED;\n"
       "     every model parameter except MassDensity and InitialP0 is REQUIRED; and any\n"
       "     unrecognised token here is an ERROR -- the material is not created.)\n"
       "    method (string) : Forward_Euler | Runge_Kutta_45_Error_Control\\ \n"
       "    tangent_type (string) : Secant (default) | Continuum | Elastic |\n"
       "        Algorithmic (Closest_Point only) | Numerical_Algorithmic_FirstOrder |\n"
       "        Numerical_Algorithmic_SecondOrder\\ \n"
       "End_Integration_Options \\ \n"
       "\n";
}

void *OPS_AllASDPlasticMaterial3Ds(void)
{
    // some kudos
    static bool first_done = false;
    if (!first_done) {
        opserr << "\n\nUsing ASDPlasticMaterial3D - Developed by: Jose Abell (UANDES), Massimo Petracca and Guido Camata (ASDEA Software Technology)\n\n";
        first_done = true;
    }

    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    // Ladruno (HB/StiffSoil integration, ledger row 337): allow 2-arg calls (yf_type only)
    if (numArgs < 2) {
        "nDMaterial ASDPlasticMaterial3D Error: Few arguments \n";
        opserr << "    numArgs = " << numArgs << endln << endln;
        print_usage();
        return nullptr;
    }

    int numData;

    int tag = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0)  {
        opserr << "nDMaterial ASDPlasticMaterial3D Error: invalid 'tag'.\n";
        return nullptr;
    }


    const char *yf_type = nullptr;
    const char *pf_type = nullptr;
    const char *el_type = nullptr;
    const char *iv_type = nullptr;

    // Now use conditional checks based on numArgs to assign values
    yf_type = numArgs >= 2 ? OPS_GetString() : " X ";
    pf_type = numArgs >= 3 ? OPS_GetString() : " X ";
    el_type = numArgs >= 4 ? OPS_GetString() : " X ";
    iv_type = numArgs >= 5 ? OPS_GetString() : " X ";
    
    // Ladruno (HB/StiffSoil integration, ledger row 337): debug-print numArgs for wildcard search
    opserr << "    numArgs = " << numArgs << endln << endln;


    cout << "Searching for instance with:\n";


    cout << "yf_type = " << yf_type << "\n";
    cout << "pf_type = " << pf_type << "\n";
    cout << "el_type = " << el_type << "\n";
    cout << "iv_type = " << iv_type << "\n";

    std::list<model_spec_t> available_models;

    asdp_parse_rejected = false;   // Ladruno (ADR-94 wp/94a)

    NDMaterial* instance = ASDPlasticMaterial3DFactory(tag, yf_type, pf_type, el_type, iv_type, available_models);

    if(std::strcmp(yf_type, "list")==0)
    {
        cout << "Available models for YF = " << pf_type << endln;
        for(model_spec_t &model : available_models)
        {
            std::string model_yf_type = std::get<0>(model);
            std::string model_pf_type = std::get<1>(model);
            std::string model_el_type = std::get<2>(model);
            std::string model_iv_type = std::get<3>(model);
            
            // Ladruno (HB/StiffSoil integration, ledger row 337): wildcard placeholder search
            if (std::strcmp(pf_type, model_yf_type.c_str())==0 || std::strcmp(pf_type, " X ")==0 )
            {
                cout << "  YF = " << model_yf_type << endl;
                cout << "  PF = " << model_pf_type << endl;
                cout << "  EL = " << model_el_type << endl;
                cout << "  IV = " << model_iv_type << endl << endln;
            }
        }
    }

    // Ladruno (ADR-94 wp/94a): the YF/PF/EL/IV combination WAS found; its options or
    // parameters were rejected. Printing the "material not found" table here would send
    // the user hunting for the wrong bug.
    if(instance==nullptr && asdp_parse_rejected)
    {
        opserr << "nDMaterial ASDPlasticMaterial3D " << tag
               << " - REJECTED: see the ASDPlasticMaterial3D error(s) above."
               << " The material was NOT created (ADR-94)." << endln;
        return nullptr;
    }

    if(instance==nullptr)
    {

        bool matches_yf = false;
        bool matches_pf = false;
        bool matches_el = false;


        for(model_spec_t &model : available_models)
        {
            std::string model_yf_type = std::get<0>(model);
            std::string model_pf_type = std::get<1>(model);
            std::string model_el_type = std::get<2>(model);
            // std::string model_iv_type = std::get<3>(model);
            if (std::strcmp(yf_type, model_yf_type.c_str())==0)
                matches_yf = true;
            if (std::strcmp(pf_type, model_pf_type.c_str())==0)
                matches_pf = true;
            if (std::strcmp(el_type, model_el_type.c_str())==0)
                matches_el = true;
        }


        cout << "\n";
        cout << "ASDPlasticMaterial3D -- ERROR! Material not found for input specification:\n";
        cout << "yf_type = " << yf_type << (matches_yf ? " :) " : " ") <<"\n";
        cout << "pf_type = " << pf_type << (matches_pf ? " :) " : " ") <<"\n";
        cout << "el_type = " << el_type << (matches_el ? " :) " : " ") <<"\n";
        cout << "iv_type = " << iv_type  <<"\n";
        cout << "\n";

        cout << endl;
        print_usage();
        cout << endl;

    }


    return instance;

}



NDMaterial*  ASDPlasticMaterial3DFactory(int instance_tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models)
{

    std::list<NDMaterial*> instance_pointers;


    #include "ASD_material_definitions.cpp"    

    //Search for the valid pointer and return that one
    for(auto instance : instance_pointers)
    {
        if(instance != nullptr)
        {
            return instance;
        }
    }

    return nullptr;
}






// Ladruno (HB/StiffSoil integration, ledger row 337): reformatted signature, no functional change
template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial3D(int instance_tag, 
        const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, 
        std::list<NDMaterial*> &instance_pointers, std::list<model_spec_t> &available_models) 
{
    auto instance = new ASDPlasticMaterial3D<EL, YF, PF, ND_TAG_ASDPlasticMaterial3D>(instance_tag);


    available_models.push_back(std::make_tuple(
        instance->getYFName(),
        instance->getPFName(),
        instance->getELName(),
        instance->getIVName())
    );

    if(
        std::strcmp(yf_type,instance->getYFName().c_str())==0 &&
        std::strcmp(pf_type,instance->getPFName().c_str())==0 &&
        std::strcmp(el_type,instance->getELName().c_str())==0 &&
        std::strcmp(iv_type,instance->getIVName().c_str())==0 
        )
    {
        // Ladruno (ADR-94 wp/94a): a rejected deck must not produce a material.
        if (!populate_ASDPlasticMaterial3D(instance))
        {
            delete instance;
            instance_pointers.push_back(nullptr);
            return nullptr;
        }
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial3D \n\n\n";
        instance_pointers.push_back(static_cast<NDMaterial*>(instance));
        return static_cast<NDMaterial*>(instance);
    }  
    else
    {
        delete instance;
        instance_pointers.push_back(nullptr);        
        return nullptr;
    }
}








template <typename T>
bool populate_ASDPlasticMaterial3D(T* instance)   // Ladruno (ADR-94 wp/94a): was void
{

    int get_one_value = 1;

    // Ladruno (ADR-94 wp/94a): names the deck actually assigned, for the
    // required-parameter check at the end of this function.
    std::set<std::string> assigned_parameters;

    // Ladruno (ADR-94 wp/94a): the valid Begin_Integration_Options tokens, printed
    // verbatim when one is not recognised so the user can see the spelling.
    static const char* const ASDP_VALID_INTEGRATION_OPTIONS =
        "f_absolute_tol, f_relative_tol, stress_absolute_tol, n_max_iterations, strict_convergence, "
        "rk45_dT_min, rk45_niter_max, return_to_yield_surface, integration_method, "
        "tangent_type, End_Integration_Options";

    cout << "\n\nDefined internal variables: \n";
    auto iv_names = instance->getInternalVariablesNames();
    for_each_in_tuple(iv_names, [instance](auto & name)
    {
        std::cout << "   "  << name << " size = " << instance->getInternalVariableSizeByName(name) << std::endl;
    });

    cout << "\n\nDefined model parameters\n";
    auto parameter_names = instance->getParameterNames();
    for_each_in_tuple(parameter_names, [](auto & name)
    {
        std::cout << "   "  <<  name << std::endl;
    });

    // Default integration options
    // Ladruno (HB/StiffSoil integration, ledger row 337): default Backward_Euler/Secant
    int method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;
    int tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Secant;
    double f_absolute_tol = 1e-6;
    double f_relative_tol = 0.0; // Ladruno (ADR-94 wp/94c, M5): 0 = off, i.e. absolute tolerance only (byte-identical to before)
    double stress_absolute_tol = 1e-6; 
    int n_max_iterations = 100;
    int return_to_yield_surface = 1;
    int strict_convergence = 0; // Ladruno (ADR-84 P2a): opt-in fail-loud on Backward_Euler non-convergence
    double rk45_dT_min = 1e-2;
    int rk45_niter_max = 110;

    // Loop over input arguments
    while (OPS_GetNumRemainingInputArgs() > 0) {

        //Current command
        const char *cmd = OPS_GetString();


        // Specifying internal variables
        if (std::strcmp(cmd, "Begin_Internal_Variables") == 0)
        {
            cout << "\n\nReading internal variables from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double iv_values[6];
                const char *iv_name = OPS_GetString();

                if (std::strcmp(iv_name, "End_Internal_Variables") == 0)
                {
                    cout << "\n\n Done reading internal variables from input\n";
                    break;
                }

                int iv_size = instance->getInternalVariableSizeByName(iv_name);
                // Ladruno (ADR-94 wp/94a): an unknown IV name returns -1 here, which was
                // then handed straight to OPS_GetDouble as a count. Fail instead.
                if (iv_size < 0 || iv_size > 6)
                {
                    opserr << "nDMaterial ASDPlasticMaterial3D - unknown internal variable '"
                           << iv_name << "' inside Begin_Internal_Variables." << endln;
                    opserr << "   Defined internal variables for this model:" << endln;
                    auto iv_names_tuple = instance->getInternalVariablesNames(); // Ladruno (ADR-94 wp/94a): GCC cannot bind a temporary to for_each_in_tuple's non-const reference
                    for_each_in_tuple(iv_names_tuple, [](auto & nm)
                    {
                        opserr << "      " << nm << endln;
                    });
                    asdp_parse_rejected = true;
                    return false;
                }
                OPS_GetDouble(&iv_size, iv_values);
                cout << iv_name << " = ";
                for (int i = 0; i < iv_size; ++i)
                {
                    cout << iv_values[i] << " ";
                }
                cout << endl;
                if (!instance->setInternalVariableByName(iv_name, iv_size, &iv_values[0]))
                {
                    // Ladruno (ADR-94 wp/94a)
                    opserr << "nDMaterial ASDPlasticMaterial3D - internal variable '"
                           << iv_name << "' was not accepted by any component." << endln;
                    asdp_parse_rejected = true;
                    return false;
                }
            }

        }


        // Specifying model parameters
        if (std::strcmp(cmd, "Begin_Model_Parameters") == 0)
        {
            cout << "\n\nReading parameters from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double param_value;
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Model_Parameters") == 0)
                {
                    cout << "\n\n Done reading parameters from input\n";
                    break;
                }

                OPS_GetDouble(&get_one_value, &param_value);
                cout << param_name << " = " << param_value << endl;
                // Ladruno (ADR-94 wp/94a): ADR-94 B1/H13 -- a misspelled parameter name
                // used to be dropped by utuple_storage's base case, leaving the intended
                // parameter at its 0.0 default (a typo'd MC_phi ran the model at phi = 0).
                if (!instance->setParameterByName(param_name, param_value))
                {
                    opserr << "nDMaterial ASDPlasticMaterial3D - unknown model parameter '"
                           << param_name << "' inside Begin_Model_Parameters." << endln;
                    opserr << "   Valid parameters for this model:" << endln;
                    auto param_names_tuple = instance->getParameterNames(); // Ladruno (ADR-94 wp/94a): GCC cannot bind a temporary to for_each_in_tuple's non-const reference
                    for_each_in_tuple(param_names_tuple, [](auto & name)
                    {
                        opserr << "      " << name << endln;
                    });
                    asdp_parse_rejected = true;
                    return false;
                }
                assigned_parameters.insert(std::string(param_name));
            }
        }


        // set_constitutive_integration_method(int method, double f_absolute_tol, double stress_absolute_tol, int n_max_iterations)
        if (std::strcmp(cmd, "Begin_Integration_Options") == 0)
        {
            cout << "\n\nReading Integration Options\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Integration_Options") == 0)
                {
                    cout << "\n\nDone reading Integration Options\n";
                    break;
                }

                // Ladruno (ADR-94 wp/94a): the chain below is a run of independent `if`s
                // with no `else`, so an unrecognised token (and its value) was silently
                // dropped -- ADR-94 B1/H13 measured `strict_convergance` (typo) producing
                // byte-identical results to omitting the option entirely.
                bool option_recognised = false;

                if (std::strcmp(param_name, "f_absolute_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &f_absolute_tol);
                    cout << "   Setting f_absolute_tol = " << f_absolute_tol << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }

                if (std::strcmp(param_name, "f_relative_tol") == 0) // Ladruno (ADR-94 wp/94c, M5)
                {
                    OPS_GetDouble(&get_one_value, &f_relative_tol);
                    cout << "   Setting f_relative_tol = " << f_relative_tol << endl;
                    option_recognised = true;
                }

                if (std::strcmp(param_name, "stress_absolute_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &stress_absolute_tol);
                    cout << "   Setting stress_absolute_tol = " << stress_absolute_tol << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }

                if (std::strcmp(param_name, "n_max_iterations") == 0)
                {
                    OPS_GetInt(&get_one_value, &n_max_iterations);
                    cout << "   Setting n_max_iterations = " << n_max_iterations << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }
                if (std::strcmp(param_name, "strict_convergence") == 0) // Ladruno (ADR-84 P2a)
                {
                    OPS_GetInt(&get_one_value, &strict_convergence);
                    cout << "   Setting strict_convergence = " << strict_convergence << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }
                if (std::strcmp(param_name, "rk45_dT_min") == 0)
                {
                    OPS_GetDouble(&get_one_value, &rk45_dT_min);
                    cout << "   Setting rk45_dT_min = " << rk45_dT_min << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }

                if (std::strcmp(param_name, "rk45_niter_max") == 0)
                {
                    OPS_GetInt(&get_one_value, &rk45_niter_max);
                    cout << "   Setting rk45_niter_max = " << rk45_niter_max << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }
                if (std::strcmp(param_name, "return_to_yield_surface") == 0)
                {
                    // OPS_GetInt(&get_one_value, &return_to_yield_surface);
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Disabled") == 0)
                        return_to_yield_surface = 0;
                    else if (std::strcmp(method_name, "One_Step_Return") == 0)
                        return_to_yield_surface = 1;
                    else if (std::strcmp(method_name, "Iterative_Return") == 0)
                        return_to_yield_surface = 2;
                    else
                    {
                        // Ladruno (ADR-94 wp/94a): was a silent default to One_Step_Return.
                        opserr << "nDMaterial ASDPlasticMaterial3D - unknown "
                               << "return_to_yield_surface '" << method_name << "'." << endln;
                        opserr << "   Valid values: Disabled, One_Step_Return, Iterative_Return"
                               << endln;
                        asdp_parse_rejected = true;
                        return false;
                    }

                    cout << "   Setting return_to_yield_surface = " << return_to_yield_surface << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }

                if (std::strcmp(param_name, "integration_method") == 0)
                {
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Forward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler;
                    else if (std::strcmp(method_name, "Forward_Euler_Subincrement") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement;
                    else if (std::strcmp(method_name, "Modified_Euler_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control_old") == 0)
                    {
                        // Ladruno (ADR-94 wp/94a): ADR-94 M8/H9 -- this integrator's drift
                        // check is an empty `if` body and its NaN guard calls exit(-1) on
                        // the whole process. The code is kept (it is the reference the
                        // non-_old RK45 was derived from) but it may not be selected.
                        opserr << "nDMaterial ASDPlasticMaterial3D - integration_method "
                               << "'Runge_Kutta_45_Error_Control_old' is REFUSED (ADR-94 M8):"
                               << " its yield-drift check is dead code and its NaN guard"
                               << " calls exit() on the process."
                               << " Use Runge_Kutta_45_Error_Control or Backward_Euler."
                               << endln;
                        asdp_parse_rejected = true;
                        return false;
                    }
                    else if (std::strcmp(method_name, "Backward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;
                    else if (std::strcmp(method_name, "Closest_Point") == 0)
                    {
                        // Ladruno (ADR-97 wp/97b): the closest-point return map.
                        // Available family by family (ADR-97 D3): P1 ships
                        // VonMises and DruckerPrager (flank + apex) with Null,
                        // Linear scalar/tensor and ArmstrongFrederick hardening.
                        // A specialization whose YF, PF or hardening law has not
                        // opted in is refused HERE rather than silently
                        // approximated with the inert zero/finite-difference
                        // defaults -- the class of defect ADR-94 M3 found.
                        if (!instance->supportsClosestPoint())
                        {
                            opserr << "nDMaterial ASDPlasticMaterial3D - "
                                   << "integration_method 'Closest_Point' is not"
                                   << " available for this model (YF "
                                   << instance->getYFName().c_str() << ", PF "
                                   << instance->getPFName().c_str() << ", IV "
                                   << instance->getIVName().c_str() << ")." << endln;
                            opserr << "   ADR-97 D3: the closest-point map is added"
                                   << " family by family. P1 covers VonMises and"
                                   << " DruckerPrager with Null / Linear / "
                                   << "ArmstrongFrederick hardening; MohrCoulomb and"
                                   << " MohrCoulombTensionCutoff are P2, HoekBrown"
                                   << " P3, StiffSoil P5. Use Backward_Euler."
                                   << endln;
                            asdp_parse_rejected = true;
                            return false;
                        }
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point;
                    }
                    else if (std::strcmp(method_name, "Backward_Euler_LineSearch") == 0)
                    {
                        // Ladruno (ADR-94 wp/94a): ADR-94 M7/H8 -- measured 2/20 steps
                        // where plain Backward_Euler does 20/20. n_max_iterations and
                        // strict_convergence are both inert inside it, its line search
                        // accepts alpha = 1 unconditionally, and its "substepping" reports
                        // success for a strain the element never asked for. Code kept,
                        // selection refused.
                        opserr << "nDMaterial ASDPlasticMaterial3D - integration_method "
                               << "'Backward_Euler_LineSearch' is REFUSED (ADR-94 M7):"
                               << " it ignores n_max_iterations, its line search cannot cut"
                               << " the step, and its substepping returns success for a"
                               << " strain increment the element never asked for."
                               << " Use Backward_Euler." << endln;
                        asdp_parse_rejected = true;
                        return false;
                    }
                    else
                    {
                        // Ladruno (ADR-94 wp/94a): was a silent default to Modified_Euler.
                        opserr << "nDMaterial ASDPlasticMaterial3D - unknown "
                               << "integration_method '" << method_name << "'." << endln;
                        opserr << "   Valid values: Forward_Euler, Forward_Euler_Subincrement, "
                               << "Modified_Euler_Error_Control, Runge_Kutta_45_Error_Control, "
                               << "Backward_Euler, Closest_Point" << endln;  // Ladruno (ADR-97 wp/97b)
                        asdp_parse_rejected = true;
                        return false;
                    }
                    cout << "   Setting integration method = " << method_name << " method_int = " << method << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)
                }

                if (std::strcmp(param_name, "tangent_type") == 0)
                {
                    const char *tangent_type_name = OPS_GetString();
                    if (std::strcmp(tangent_type_name, "Elastic") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
                    else if (std::strcmp(tangent_type_name, "Continuum") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Continuum;
                    else if (std::strcmp(tangent_type_name, "Secant") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Secant;
                    else if (std::strcmp(tangent_type_name, "Algorithmic") == 0)
                        // Ladruno (ADR-97 wp/97b, D2): the EXACT consistent tangent
                        // of the Closest_Point map.  The enum value has existed
                        // since upstream with no dispatch case anywhere -- it was
                        // dead, which is why nobody has been silently getting it.
                        // The cross-check against integration_method is below,
                        // after the whole option loop (the two tokens may arrive in
                        // either order).
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic;
                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_FirstOrder") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_FirstOrder;
                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_SecondOrder") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_SecondOrder;
                    else
                    {
                        // Ladruno (ADR-94 wp/94a): was a silent default to Elastic.
                        opserr << "nDMaterial ASDPlasticMaterial3D - unknown tangent_type '"
                               << tangent_type_name << "'." << endln;
                        opserr << "   Valid values: Elastic, Continuum, Secant, Algorithmic, "  // Ladruno (ADR-97 wp/97b)
                               << "Numerical_Algorithmic_FirstOrder, "
                               << "Numerical_Algorithmic_SecondOrder" << endln;
                        asdp_parse_rejected = true;
                        return false;
                    }
                    cout << "   Setting tangent type = " << tangent_type_name << " tangent int = " << tangent << endl;
                    option_recognised = true;   // Ladruno (ADR-94 wp/94a)

                }

                // Ladruno (ADR-94 wp/94a): the `else` the if-chain above never had.
                if (!option_recognised)
                {
                    opserr << "nDMaterial ASDPlasticMaterial3D - unknown option '"
                           << param_name << "' inside Begin_Integration_Options." << endln;
                    opserr << "   Valid options: " << ASDP_VALID_INTEGRATION_OPTIONS << endln;
                    asdp_parse_rejected = true;
                    return false;
                }

            }
        }
    }

    // Ladruno (ADR-94 wp/94a): ADR-94 B1 / ADR-84 P2(e) -- every model parameter this
    // specialization declares must be given a value. Unset ones silently ran at 0.0,
    // which for a friction angle or a cohesion is a DIFFERENT, weaker material that
    // still converges. MassDensity and InitialP0 are genuinely optional (0 = no mass,
    // no geostatic seed).
    {
        std::string missing;
        int n_missing = 0;
        auto required_names_tuple = instance->getParameterNames(); // Ladruno (ADR-94 wp/94a): GCC cannot bind a temporary to for_each_in_tuple's non-const reference
        for_each_in_tuple(required_names_tuple,
            [&missing, &n_missing, &assigned_parameters](auto & name)
        {
            std::string n(name);
            if (n == "MassDensity" || n == "InitialP0")
                return;
            if (assigned_parameters.find(n) == assigned_parameters.end())
            {
                missing += (n_missing ? ", " : "");
                missing += n;
                ++n_missing;
            }
        });

        if (n_missing > 0)
        {
            opserr << "nDMaterial ASDPlasticMaterial3D - " << n_missing
                   << " required model parameter(s) were never given a value: "
                   << missing.c_str() << endln;
            opserr << "   An unset parameter defaults to 0.0, which is a different"
                   << " material that still converges (ADR-94 B1). Set them inside"
                   << " Begin_Model_Parameters ... End_Model_Parameters." << endln;
            asdp_parse_rejected = true;
            return false;
        }
    }

    // Ladruno (ADR-97 wp/97b, D2): a consistent tangent is defined only relative
    // to a specific committed map.  Offering `Algorithmic` on the cutting-plane
    // Backward_Euler would ship a FOURTH almost-right tangent, which is exactly
    // the class of defect ADR-94 M3 found (Continuum 57 %, Secant 80 %, Elastic
    // 103 %, Numerical_Algorithmic 31 % against a central difference of the
    // material's own committed response).  Refuse, naming both tokens.
    if (tangent == (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic
            && method != (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point)
    {
        opserr << "nDMaterial ASDPlasticMaterial3D - tangent_type 'Algorithmic' is"
               << " the exact consistent tangent of the 'Closest_Point' return map"
               << " and is REFUSED with any other integration_method (ADR-97 D2)."
               << endln;
        opserr << "   Set 'integration_method Closest_Point', or pick a tangent_type"
               << " that the chosen integrator defines: Secant, Continuum, Elastic,"
               << " Numerical_Algorithmic_FirstOrder, Numerical_Algorithmic_SecondOrder."
               << endln;
        asdp_parse_rejected = true;
        return false;
    }

    instance->set_constitutive_integration_method(method, tangent, f_absolute_tol, stress_absolute_tol, n_max_iterations, return_to_yield_surface, rk45_niter_max, rk45_dT_min, strict_convergence, f_relative_tol); // Ladruno (ADR-84 P2a); Ladruno (ADR-94 wp/94c, M5)

    return true;   // Ladruno (ADR-94 wp/94a)
}


#endif // _EIGEN3
