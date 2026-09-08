
// Ladruno (ADR-97 wp/97b): the CP hardening-derivative forwarders and the
// hardening_has_cp_derivatives trait live in HardeningFunction.h.
#include "HardeningFunction.h"

//Base template struct for all internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
struct InternalVariableType {
    static constexpr const char* NAME = NAMER::name;
    EvolvingVariableType trial_value;
    EvolvingVariableType committed_value;
    InternalVariableType() = default;
    InternalVariableType(EvolvingVariableType x) : trial_value(x), committed_value(x) {}
    void commit() {committed_value = trial_value;};
    void revert() {trial_value = committed_value;};
    const char *  getName() const { return NAME;}//(std::string(NAME) + std::string("(") + std::string(HardeningType::NAME) + std::string(")")).c_str(); }
    std::string getFullName() const 
    { 
        std::string fullname = std::string(NAME) + std::string("(") + std::string(HardeningType::NAME) + std::string(")");
        // cout << "fullname = " << fullname << endl;
        return fullname; 
    }
    int size() const {return trial_value.size();}

    template <class ParameterStorageType>
    EvolvingVariableType hardening_function(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters) const
    {
        return HardeningType::f(trial_value, depsilon, m, sigma, parameters);
    }

    // Assignment operator
    InternalVariableType& operator=(const InternalVariableType& other) {
        if (this != &other) { // Check for self-assignment
            trial_value = other.trial_value;
            committed_value = other.committed_value;
        }
        return *this;
    }

    void updateTrialValueFromOther(const InternalVariableType& other)
    {
        if (this != &other) { // Check for self-assignment
            trial_value = other.trial_value;
        }
        return;
    }

    // Ladruno (ADR-97 wp/97b): closest-point Jacobian blocks for THIS internal
    // variable, evaluated (like hardening_function above) at the TRIAL value --
    // i.e. at q_{n+1}, which is what makes the ADR-97 hardening update implicit.
    // `out` is a 6x6 buffer; only rows 0..size()-1 (and, for dh_dq, the same
    // number of columns) are meaningful.
    template <class ParameterStorageType>
    void hardening_dh_dq(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters,
                VoigtMatrix& out) const
    {
        HardeningType::dh_dq(trial_value, depsilon, m, sigma, parameters, out);
    }

    template <class ParameterStorageType>
    void hardening_dh_dm(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters,
                VoigtMatrix& out) const
    {
        HardeningType::dh_dm(trial_value, depsilon, m, sigma, parameters, out);
    }

    // Ladruno (ADR-97 wp/97b): does this IV's hardening law have analytic
    // closest-point derivatives?  Folded over the whole IV tuple at compile time
    // by ASDPlasticMaterial3D so the parser can refuse Closest_Point for a
    // specialization carrying an unconverted hardening law.
    static constexpr bool hardening_supports_cp()
    {
        return hardening_has_cp_derivatives<HardeningType>::value;
    }

    // Ladruno (ADR-97 wp/97c): does this internal variable's hardening law leave
    // it FIXED (h == 0)?  Required by the principal-space Mohr-Coulomb return,
    // which projects onto a surface it assumes does not move.
    static constexpr bool hardening_is_perfectly_plastic()
    {
        return hardening_is_inert<HardeningType>::value;
    }

    using parameters_t = typename HardeningType::parameters_t;


};

//Stream operator for internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
std::ostream& operator<<(std::ostream& os, const InternalVariableType<EvolvingVariableType, HardeningType, NAMER>& param) {
    os << param.getFullName() << 
    "\n               : [Trial: " << param.trial_value.transpose() << ", Committed: " << param.committed_value.transpose() << "]";
    os << endl;
    // os << "   HardeningType::parameters --> " << typeid(typename InternalVariableType<EvolvingVariableType, HardeningType, NAMER>::parameters_t).name() << endl;
    return os;
}
