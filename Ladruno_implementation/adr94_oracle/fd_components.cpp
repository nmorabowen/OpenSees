// ADR-94 R3a -- component derivative harness (standalone, header-only route).
//
// Instantiates the SAME YF/PF/EL template specializations the ASDP registry
// uses (copied verbatim from SRC/material/nD/ASDPlasticMaterial3D/
// ASD_material_definitions.cpp, one canonical combo per registered YF
// family) and checks:
//   - YF: analytic df_dsigma_ij vs. central-difference of operator() over a
//     stress cloud (random + targeted Lode-edge/apex/J2->0 points).
//   - PF: m finite everywhere; for the associated VM/DP pairs, m equals the
//     YF normal direction.
//   - EL: E(sigma) symmetric + positive-definite over the cloud.
//
// Build (MinGW g++, header-only -- no OpenSees link required, see report):
//   g++ -std=c++17 -O2 -I <eigen3> @incflags.txt -I <YieldFunctions dir>
//       -I <PlasticFlowDirections dir> -I <ElasticityModels dir>
//       fd_components.cpp -o fd_components.exe
//
// No SRC/ edits. Read-only against the ASDPlasticMaterial3D headers.

#include <iostream>
#include <iomanip>
#include <tuple>
using namespace std;

#include "VonMises_YF.h"
#include "DruckerPrager_YF.h"
#include "MohrCoulomb_YF.h"
#include "HoekBrown_YF.h"
#include "StiffSoilCap_YF.h"
#include "StiffSoilShear_YF.h"
#include "MohrCoulombTensionCutoff_YF.h"

#include "VonMises_PF.h"
#include "DruckerPrager_PF.h"
#include "MohrCoulomb_PF.h"
#include "HoekBrown_PF.h"
#include "StiffSoilCap_PF.h"
#include "StiffSoilShear_PF.h"
#include "MohrCoulombTensionCutoff_PF.h"

#include "LinearIsotropic3D_EL.h"
#include "StiffSoil_EL.h"

#include "../utuple_storage.h"
#include "../AllASDInternalVariableTypes.h"
#include "../AllASDHardeningFunctions.h"
#include "../AllASDModelParameterTypes.h"

#include <vector>
#include <string>
#include <random>
#include <cmath>

// ---------------------------------------------------------------------
// Stress cloud: random + targeted (Lode edges, apex/hydrostatic, J2->0)
// ---------------------------------------------------------------------
// Category tags: "smooth" points (random + Lode-edge, where the yield
// surfaces are expected to be differentiable) vs. "singular" points (the
// exact hydrostatic/apex axis and J2->0, where sqrt(J2)-type surfaces have a
// genuine cone corner and a central-difference gradient is NOT expected to
// match any single analytic sub-gradient -- see the report).
struct TaggedSigma { VoigtVector sigma; std::string category; };

static std::vector<TaggedSigma> build_stress_cloud()
{
    std::vector<TaggedSigma> cloud;
    std::mt19937 rng(94010203u); // fixed seed -- reproducible
    std::uniform_real_distribution<double> normal_d(-60.0, 60.0);
    std::uniform_real_distribution<double> shear_d(-25.0, 25.0);

    // 1) ~160 random general (non-diagonal) stress states
    for (int i = 0; i < 160; ++i)
    {
        double sxx = normal_d(rng), syy = normal_d(rng), szz = normal_d(rng);
        double sxy = shear_d(rng), syz = shear_d(rng), sxz = shear_d(rng);
        cloud.push_back({VoigtVector(sxx, syy, szz, sxy, syz, sxz), "smooth"});
    }

    // 2) Lode-edge points: theta = +-30 deg (triaxial compression/extension),
    //    for several (p, r=sqrt(2 J2)) pairs
    std::vector<double> ps = {1.0, 5.0, 20.0, 50.0};
    std::vector<double> rs = {2.0, 10.0, 30.0};
    for (double p : ps)
        for (double r : rs)
            for (double theta_deg : {30.0, -30.0})
            {
                double theta = theta_deg * M_PI / 180.0;
                double s1 = p + std::sqrt(2.0 / 3.0) * r * std::cos(theta);
                double s2 = p + std::sqrt(2.0 / 3.0) * r * std::cos(theta - 2.0 * M_PI / 3.0);
                double s3 = p + std::sqrt(2.0 / 3.0) * r * std::cos(theta + 2.0 * M_PI / 3.0);
                cloud.push_back({VoigtVector(s1, s2, s3, 0, 0, 0), "smooth"});
            }

    // 3) Apex / hydrostatic axis: pure hydrostatic states, deviator = 0
    //    (a genuine corner of any J2-based surface -- expected non-smooth)
    for (double p : {0.0, 1.0, 5.0, 20.0, 50.0, 100.0})
        cloud.push_back({VoigtVector(p, p, p, 0, 0, 0), "singular"});

    // 4) J2 -> 0: hydrostatic plus a tiny deviatoric perturbation
    //    (still within FD step size of the corner -- expected non-smooth)
    for (double p : {1.0, 5.0, 20.0, 50.0})
        cloud.push_back({VoigtVector(p + 1e-7, p - 0.5e-7, p - 0.5e-7, 1e-8, 0, 0), "singular"});

    return cloud; // ~ 160 + 24 + 6 + 4 = 194 points
}

// ---------------------------------------------------------------------
// FD checker for a YF instance
// ---------------------------------------------------------------------
struct YFResult
{
    std::string name;
    int n_points = 0;
    int n_nonfinite_f = 0;
    int n_nonfinite_df = 0;
    // "voigt" = analytic vs the RAW central difference (perturbing the six
    // independent Voigt slots this object actually stores -- unambiguous
    // ground truth for d f(v)/dv_i).
    double max_rel_err_smooth = 0.0;
    VoigtVector worst_point_smooth;
    double max_rel_err_singular = 0.0;
    VoigtVector worst_point_singular;
    // "tensor" = analytic vs the SAME central difference with its three
    // shear slots (indices 3,4,5 = xy,yz,xz) halved. This is not a separate
    // measurement -- it is the fixed identity d/dv12 = d/ds12 + d/ds21 =
    // 2*(tensor derivative) for any scalar function of a symmetric tensor,
    // so tensor_fd = voigt_fd with shear halved, always. Comparing analytic
    // against BOTH tells us which convention df_dsigma_ij actually returns
    // (coordinator-requested convention adjudication, ADR-94 R3a followup).
    double max_rel_err_smooth_tensor = 0.0;
    VoigtVector worst_point_smooth_tensor;
    double h_used = 0.0;
};

template <class YF, class IVS, class PS>
YFResult check_yf(const std::string& name, const YF& yf, const IVS& ivs, const PS& ps,
                   const std::vector<TaggedSigma>& cloud, double h = 1e-6)
{
    YFResult r;
    r.name = name;
    r.h_used = h;
    for (const auto& ts : cloud)
    {
        const VoigtVector& sigma = ts.sigma;
        r.n_points++;
        double f0 = yf(sigma, ivs, ps);
        if (!std::isfinite(f0)) { r.n_nonfinite_f++; continue; }

        VoigtVector analytic = yf.df_dsigma_ij(sigma, ivs, ps); // copy out of the cached ref
        VoigtVector fd;
        bool df_ok = true;
        for (int i = 0; i < 6; ++i)
        {
            VoigtVector sp = sigma, sm = sigma;
            sp(i) += h; sm(i) -= h;
            double fp = yf(sp, ivs, ps);
            double fm = yf(sm, ivs, ps);
            if (!std::isfinite(fp) || !std::isfinite(fm)) { df_ok = false; break; }
            fd(i) = (fp - fm) / (2.0 * h);
        }
        if (!df_ok || !analytic.allFinite())
        {
            r.n_nonfinite_df++;
            continue;
        }

        VoigtVector fd_tensor = fd;
        fd_tensor(3) *= 0.5; fd_tensor(4) *= 0.5; fd_tensor(5) *= 0.5;

        double num = (analytic - fd).norm();
        double den = std::max(fd.norm(), 1e-8);
        double rel = num / den;

        double num_t = (analytic - fd_tensor).norm();
        double den_t = std::max(fd_tensor.norm(), 1e-8);
        double rel_t = num_t / den_t;

        if (ts.category == "smooth")
        {
            if (rel > r.max_rel_err_smooth) { r.max_rel_err_smooth = rel; r.worst_point_smooth = sigma; }
            if (rel_t > r.max_rel_err_smooth_tensor) { r.max_rel_err_smooth_tensor = rel_t; r.worst_point_smooth_tensor = sigma; }
        }
        else
        {
            if (rel > r.max_rel_err_singular) { r.max_rel_err_singular = rel; r.worst_point_singular = sigma; }
        }
    }
    return r;
}


// ---------------------------------------------------------------------
// PF checker: finiteness (+ optional equality to a reference normal)
// ---------------------------------------------------------------------
struct PFResult
{
    std::string name;
    int n_points = 0;
    int n_nonfinite = 0;
    double max_assoc_err = -1.0; // -1 => not checked
};

template <class PF, class IVS, class PS>
PFResult check_pf(const std::string& name, const PF& pf, const IVS& ivs, const PS& ps,
                   const std::vector<VoigtVector>& cloud,
                   const std::vector<VoigtVector>* reference_normals = nullptr)
{
    PFResult r;
    r.name = name;
    VoigtVector depsilon; // unused by these PFs, pass zero
    for (size_t idx = 0; idx < cloud.size(); ++idx)
    {
        r.n_points++;
        VoigtVector m = pf(depsilon, cloud[idx], ivs, ps);
        if (!m.allFinite()) { r.n_nonfinite++; continue; }
        if (reference_normals)
        {
            double err = (m - (*reference_normals)[idx]).norm();
            r.max_assoc_err = std::max(r.max_assoc_err, err);
        }
    }
    return r;
}

// ---------------------------------------------------------------------
// EL checker: symmetry + SPD
// ---------------------------------------------------------------------
struct ELResult
{
    std::string name;
    int n_points = 0;
    int n_asym = 0;
    int n_not_spd = 0;
    double max_asym = 0.0;
    double min_eig = 1e300;
};

template <class EL, class PS>
ELResult check_el(const std::string& name, const EL& el, const PS& ps,
                   const std::vector<VoigtVector>& cloud)
{
    ELResult r;
    r.name = name;
    for (const auto& sigma : cloud)
    {
        r.n_points++;
        VoigtMatrix E = el(sigma, ps);
        double asym = (E - E.transpose()).norm();
        r.max_asym = std::max(r.max_asym, asym);
        if (asym > 1e-6 * std::max(1.0, E.norm())) r.n_asym++;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> es(0.5 * (E + E.transpose()));
        double lmin = es.eigenvalues().minCoeff();
        r.min_eig = std::min(r.min_eig, lmin);
        if (lmin <= 0.0) r.n_not_spd++;
    }
    return r;
}

static void print_yf(const YFResult& r)
{
    cout << left << setw(30) << r.name
         << " max_rel_err(smooth)=" << scientific << setprecision(3) << r.max_rel_err_smooth
         << "  vs_tensor(smooth)=" << scientific << setprecision(3) << r.max_rel_err_smooth_tensor
         << "  max_rel_err(singular)=" << r.max_rel_err_singular
         << "  nonfinite(f)=" << r.n_nonfinite_f
         << "  nonfinite(df)=" << r.n_nonfinite_df
         << "  n=" << r.n_points
         << "  worst_smooth=[" << fixed << setprecision(3) << r.worst_point_smooth.transpose() << "]"
         << "\n";
}

static void print_pf(const PFResult& r)
{
    cout << left << setw(30) << r.name
         << " nonfinite=" << r.n_nonfinite << "/" << r.n_points;
    if (r.max_assoc_err >= 0.0)
        cout << "  max_assoc_err=" << scientific << setprecision(3) << r.max_assoc_err;
    cout << "\n";
}

static void print_el(const ELResult& r)
{
    cout << left << setw(30) << r.name
         << " max_asym=" << scientific << setprecision(3) << r.max_asym
         << "  min_eig=" << r.min_eig
         << "  not_spd=" << r.n_not_spd << "/" << r.n_points
         << "\n";
}

int main()
{
    cout << std::boolalpha;
    auto cloud = build_stress_cloud();
    std::vector<VoigtVector> plain_cloud;
    plain_cloud.reserve(cloud.size());
    for (auto& ts : cloud) plain_cloud.push_back(ts.sigma);
    cout << "Stress cloud size: " << cloud.size() << "\n\n";

    // ---------------- VonMises_YF / VonMises_PF ----------------
    {
        using AlphaT = BackStress<TensorLinearHardeningFunction>;
        using KT = YieldStress<ScalarLinearHardeningFunction>;
        using YF = VonMises_YF<AlphaT, KT>;
        using PFAssoc = VonMises_PF<BackStress<NullHardeningTensorFunction>>; // registry pair #1

        YF yf;
        PFAssoc pf;

        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs;
        ivs.set(AlphaT(VoigtVector(0,0,0,0,0,0)));
        ivs.set(KT(10.0));
        using ps_t = utuple_storage<YF::parameters_t>; // empty tuple
        ps_t ps;

        using iv_pf_t = utuple_storage<PFAssoc::internal_variables_t>;
        iv_pf_t iv_pf;
        iv_pf.set(BackStress<NullHardeningTensorFunction>(VoigtVector(0,0,0,0,0,0)));
        using ps_pf_t = utuple_storage<PFAssoc::parameters_t>;
        ps_pf_t ps_pf;

        auto r = check_yf("VonMises_YF", yf, ivs, ps, cloud);
        print_yf(r);

        std::vector<VoigtVector> normals;
        normals.reserve(cloud.size());
        for (auto& s : plain_cloud) normals.push_back(yf.df_dsigma_ij(s, ivs, ps));
        auto rpf = check_pf("VonMises_PF (assoc)", pf, iv_pf, ps_pf, plain_cloud, &normals);
        print_pf(rpf);
    }

    // ---------------- DruckerPrager_YF / DruckerPrager_PF ----------------
    {
        using AlphaT = BackStress<TensorLinearHardeningFunction>;
        using CohT = DP_cohesion<ScalarLinearHardeningFunction>;
        using YF = DruckerPrager_YF<AlphaT, CohT>;
        using EtaT = DP_cohesion<ScalarLinearHardeningFunction>; // PF's 2nd IV slot (registry pair)
        using PFAssoc = DruckerPrager_PF<AlphaT, EtaT>;

        YF yf;
        PFAssoc pf;

        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs;
        ivs.set(AlphaT(VoigtVector(0,0,0,0,0,0)));
        ivs.set(CohT(5.0));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(DP_xi_c(5.0));
        ps.set(DP_eta(0.3));

        using iv_pf_t = utuple_storage<PFAssoc::internal_variables_t>;
        iv_pf_t iv_pf;
        iv_pf.set(AlphaT(VoigtVector(0,0,0,0,0,0)));
        iv_pf.set(EtaT(5.0));
        using ps_pf_t = utuple_storage<PFAssoc::parameters_t>;
        ps_pf_t ps_pf;
        ps_pf.set(DP_etabar(0.3)); // == eta -> associated flow

        auto r = check_yf("DruckerPrager_YF", yf, ivs, ps, cloud);
        print_yf(r);

        std::vector<VoigtVector> normals;
        normals.reserve(cloud.size());
        for (auto& s : plain_cloud) normals.push_back(yf.df_dsigma_ij(s, ivs, ps));
        auto rpf = check_pf("DruckerPrager_PF (etabar=eta)", pf, iv_pf, ps_pf, plain_cloud, &normals);
        print_pf(rpf);
    }

    // ---------------- MohrCoulomb_YF / MohrCoulomb_PF ----------------
    {
        using IVT = BackStress<NullHardeningTensorFunction>;
        using YF = MohrCoulomb_YF<IVT>;
        using PF = MohrCoulomb_PF<IVT>;
        YF yf; PF pf;

        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs; ivs.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(MC_phi(30.0));
        ps.set(MC_c(10.0));
        ps.set(MC_ds(1e-4));

        auto r = check_yf("MohrCoulomb_YF", yf, ivs, ps, cloud);
        print_yf(r);

        using iv_pf_t = utuple_storage<PF::internal_variables_t>;
        iv_pf_t iv_pf; iv_pf.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_pf_t = utuple_storage<PF::parameters_t>; // <MC_phi,MC_c,MC_ds,MC_psi>
        ps_pf_t ps_pf;
        ps_pf.set(MC_phi(30.0));
        ps_pf.set(MC_c(10.0));
        ps_pf.set(MC_ds(1e-4));
        ps_pf.set(MC_psi(10.0)); // dilatancy != phi -> non-associated by design
        auto rpf = check_pf("MohrCoulomb_PF", pf, iv_pf, ps_pf, plain_cloud, nullptr);
        print_pf(rpf);
    }

    // ---------------- HoekBrown_YF / HoekBrown_PF ----------------
    {
        using IVT = BackStress<NullHardeningTensorFunction>;
        using YF = HoekBrown_YF<IVT>;
        using PF = HoekBrown_PF<IVT>;
        YF yf; PF pf;

        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs; ivs.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(HB_sigci(30.0));
        ps.set(HB_mb(2.0));
        ps.set(HB_s(0.01));
        ps.set(HB_a(0.5));
        ps.set(HB_ds(1e-4));

        auto r = check_yf("HoekBrown_YF", yf, ivs, ps, cloud);
        print_yf(r);

        using iv_pf_t = utuple_storage<PF::internal_variables_t>;
        iv_pf_t iv_pf; iv_pf.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_pf_t = utuple_storage<PF::parameters_t>; // <HB_sigci,HB_mb_psi,HB_s,HB_a,HB_ds>
        ps_pf_t ps_pf;
        ps_pf.set(HB_sigci(30.0));
        ps_pf.set(HB_mb_psi(1.0)); // reduced (dilation) mb -> non-associated by design
        ps_pf.set(HB_s(0.01));
        ps_pf.set(HB_a(0.5));
        ps_pf.set(HB_ds(1e-4));
        auto rpf = check_pf("HoekBrown_PF", pf, iv_pf, ps_pf, plain_cloud, nullptr);
        print_pf(rpf);
    }

    // ---------------- StiffSoilCap_YF ----------------
    {
        using IVT = CapPressure;
        using YF = StiffSoilCap_YF<IVT>;
        YF yf;
        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs; ivs.set(IVT(20.0));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(MC_phi(30.0));
        ps.set(MC_ds(1e-4));
        ps.set(SS_alpha(0.5));
        ps.set(SS_pref(100.0));
        ps.set(SS_m(0.5));
        ps.set(SS_beta(1.0));

        auto r = check_yf("StiffSoilCap_YF", yf, ivs, ps, cloud);
        print_yf(r);

        using PF = StiffSoilCap_PF<IVT>;
        PF pf;
        using iv_pf_t = utuple_storage<PF::internal_variables_t>;
        iv_pf_t iv_pf; iv_pf.set(IVT(20.0));
        using ps_pf_t = utuple_storage<PF::parameters_t>; // <MC_phi,MC_ds,SS_alpha>
        ps_pf_t ps_pf;
        ps_pf.set(MC_phi(30.0));
        ps_pf.set(MC_ds(1e-4));
        ps_pf.set(SS_alpha(0.5));
        auto rpf = check_pf("StiffSoilCap_PF", pf, iv_pf, ps_pf, plain_cloud, nullptr);
        print_pf(rpf);
    }

    // ---------------- StiffSoilShear_YF ----------------
    {
        using IVT = EpsQpShear;
        using YF = StiffSoilShear_YF<IVT>;
        YF yf;
        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs; ivs.set(IVT(0.001));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(MC_phi(30.0));
        ps.set(MC_c(10.0));
        ps.set(MC_ds(1e-4));
        ps.set(SS_E50_ref(20000.0));
        ps.set(SS_Eur_ref(60000.0));
        ps.set(SS_Rf(0.9));
        ps.set(SS_m(0.5));
        ps.set(SS_pref(100.0));

        auto r = check_yf("StiffSoilShear_YF", yf, ivs, ps, cloud);
        print_yf(r);

        using PF = StiffSoilShear_PF<IVT>;
        PF pf;
        using iv_pf_t = utuple_storage<PF::internal_variables_t>;
        iv_pf_t iv_pf; iv_pf.set(IVT(0.001));
        using ps_pf_t = utuple_storage<PF::parameters_t>; // <MC_phi,MC_psi,MC_c,MC_ds>
        ps_pf_t ps_pf;
        ps_pf.set(MC_phi(30.0));
        ps_pf.set(MC_psi(10.0));
        ps_pf.set(MC_c(10.0));
        ps_pf.set(MC_ds(1e-4));
        auto rpf = check_pf("StiffSoilShear_PF", pf, iv_pf, ps_pf, plain_cloud, nullptr);
        print_pf(rpf);
    }

    // ---------------- MohrCoulombTensionCutoff_YF ----------------
    {
        using IVT = BackStress<NullHardeningTensorFunction>;
        using YF = MohrCoulombTensionCutoff_YF<IVT>;
        YF yf;
        using iv_t = utuple_storage<YF::internal_variables_t>;
        iv_t ivs; ivs.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_t = utuple_storage<YF::parameters_t>;
        ps_t ps;
        ps.set(MC_phi(30.0));
        ps.set(MC_c(10.0));
        ps.set(MC_ds(1e-4));
        ps.set(MC_psi(10.0));
        ps.set(TC_min_stress(-5.0));

        auto r = check_yf("MohrCoulombTensionCutoff_YF", yf, ivs, ps, cloud);
        print_yf(r);

        using PF = MohrCoulombTensionCutoff_PF<IVT>;
        PF pf;
        using iv_pf_t = utuple_storage<PF::internal_variables_t>;
        iv_pf_t iv_pf; iv_pf.set(IVT(VoigtVector(0,0,0,0,0,0)));
        using ps_pf_t = utuple_storage<PF::parameters_t>; // <MC_phi,MC_c,MC_ds,MC_psi,TC_min_stress>
        ps_pf_t ps_pf;
        ps_pf.set(MC_phi(30.0));
        ps_pf.set(MC_c(10.0));
        ps_pf.set(MC_ds(1e-4));
        ps_pf.set(MC_psi(10.0));
        ps_pf.set(TC_min_stress(-5.0));
        auto rpf = check_pf("MohrCoulombTensionCutoff_PF", pf, iv_pf, ps_pf, plain_cloud, nullptr);
        print_pf(rpf);
    }

    cout << "\n-- Elasticity SPD checks --\n";
    // ---------------- LinearIsotropic3D_EL ----------------
    {
        LinearIsotropic3D_EL el;
        using ps_t = utuple_storage<LinearIsotropic3D_EL::parameters_t>;
        ps_t ps;
        ps.set(YoungsModulus(30000.0));
        ps.set(PoissonsRatio(0.3));
        auto r = check_el("LinearIsotropic3D_EL", el, ps, plain_cloud);
        print_el(r);
    }
    // ---------------- StiffSoil_EL ----------------
    {
        StiffSoil_EL el;
        using ps_t = utuple_storage<StiffSoil_EL::parameters_t>;
        ps_t ps;
        ps.set(SS_Eur_ref(60000.0));
        ps.set(PoissonsRatio(0.3));
        ps.set(SS_pref(100.0));
        ps.set(SS_m(0.5));
        ps.set(MC_phi(30.0));
        ps.set(MC_c(10.0));
        // StiffSoil_EL is stress-dependent (power-law) -- restrict to the
        // positive-mean-stress subset of the cloud (its intended domain).
        std::vector<VoigtVector> admissible;
        for (auto& s : plain_cloud) if (s.meanStress() > 0.5) admissible.push_back(s);
        auto r = check_el("StiffSoil_EL", el, ps, admissible);
        print_el(r);
    }

    return 0;
}
