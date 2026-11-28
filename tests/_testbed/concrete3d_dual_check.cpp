// concrete3d_dual_check.cpp -- WP-264 (dual-number program, #1174): the damage-driver gradients of
// LadrunoConcrete3DKernel.h by dual numbers, against the finite differences they replace.
//
// damagedTangent differentiates four scalar damage drivers f(sig_bar) (eps_tilde, x_s, alpha_c, eps_tilde_c) with
// respect to the effective stress. Before WP-264 it took a per-Voigt-slot central FD (dscalarDsig: 12 driver calls, 12
// eigendecompositions per gradient, step DRIVER_FD_REL * S). WP-264 replaces it off the hydrostatic axis by
// dscalarDsigDual: a Dual<3> pass in the principal values and the spectral chain rule, with the Lode angle carried
// through cos 3 theta (lodeRT). The axis keeps its one-sided FD (dscalarDsigAxis). This check proves, with g++ only:
//   B  BIT-IDENTITY. The templated double entry points (invariants, equivStrainGeneral, alphaCompression, tcPhi,
//      damageDrivers, compressiveEquivStrain, scalarDriver, dscalarDsig) against a verbatim transcription of the
//      pre-WP-264 code (namespace old below), bit for bit, on random, degenerate and record states; and the value of
//      every Dual pass == the double result bit for bit.
//   G  GRADIENT ORACLE. dscalarDsigDual vs a step-plateau central FD (h = S 10^-k, k = 2..9, the k where successive
//      steps agree best) and vs the shipped FD (DRIVER_FD_REL S), in MPa, kPa and Pa, per driver, on random states,
//      on exactly repeated principal pairs (uniaxial, equal-biaxial, confined), on near-repeated pairs, and on the
//      replay record states (tests/_testbed/concrete3d_replay_records.txt: the committed and converged effective
//      stresses of every pass). Plus two exact cross-checks the FD cannot give: the gradient is independent of the
//      eigenbasis eig3sym picks inside a repeated pair, and for eps_tilde and x_s it equals a Dual<6> pass through the
//      invariants of the FULL tensor (no eigendecomposition at all).
//   M  BREAK-THE-GATE. Four mutants the oracle must flag: the shear factor 2 dropped, the Lode derivative dropped, the
//      eigenvalue/eigenvector pairing scrambled, and the naive dual through acos (the meridian 0 * inf).
//   T  DAMAGED TANGENT. returnMap's damaged tangent with dual drivers vs FD drivers (driverGradFD(), same build) on the
//      replay record states and on three damage-advancing paths in MPa, kPa and Pa; both against the FD of the damaged
//      map; the stress and returned State must be bit-identical between the two modes.
//   P  TIMING. Per-call time of the driver gradient (dual vs FD) and of returnMap with the tangent.
//
// Build + run from the repo root (no -march=native):
//   g++ -std=c++17 -O2 -I. tests/_testbed/concrete3d_dual_check.cpp -o c3ddual && ./c3ddual
//   (optional argv[1]: the replay fixture; default tests/_testbed/concrete3d_replay_records.txt). Exit 0 = all pass.
//
// The replay testbed is included whole for its fixture loader (fillParams, fillState, the record inputs); its main is
// renamed away. It is compiled under this file's own libm tag, which only labels output: nothing here compares a value
// with the fixture's recorded outputs.

#define C3D_REPLAY_LIBM "wp264-dual-check"
#define main c3d_replay_main_unused
#include "concrete3d_replay.cpp"
#undef main

#include <array>
#include <random>

using ladruno_dual::Dual;

static int gFails = 0;
static void gate(bool ok, const std::string& name)
{
    std::printf("  %s: %s\n", ok ? "ok  " : "FAIL", name.c_str());
    if (!ok) ++gFails;
}
static bool bitEq(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// ===========================================================================
// old: the pre-WP-264 code, VERBATIM from origin/ladruno 94dbd6a85b95d3603f1cf449802fbe91ab45d6b3 (the kernel blob
// the replay fixtures were recorded on). lodeR and eig3sym are unchanged by WP-264 and are taken from the kernel.
// ===========================================================================
namespace old {
using namespace Ladruno::Concrete3D;
inline void invariants(const double sig[6], double& xi, double& rho, double& theta)
{
    const double I1 = sig[0] + sig[1] + sig[2];
    const double p = I1 / 3.0;
    const double d0 = sig[0] - p, d1 = sig[1] - p, d2 = sig[2] - p;
    const double s01 = sig[3], s12 = sig[4], s02 = sig[5];
    const double J2 = 0.5 * (d0 * d0 + d1 * d1 + d2 * d2) + s01 * s01 + s12 * s12 + s02 * s02;
    const double J3 = d0 * (d1 * d2 - s12 * s12)
                    - s01 * (s01 * d2 - s12 * s02)
                    + s02 * (s01 * s12 - d1 * s02);
    xi = I1 / SQRT3;
    rho = (J2 > 0.0) ? std::sqrt(2.0 * J2) : 0.0;
    if (J2 <= 1.0e-300) {
        theta = 0.0;
    } else {
        double c3 = (3.0 * SQRT3 / 2.0) * J3 / std::pow(J2, 1.5);
        if (c3 > 1.0) c3 = 1.0; else if (c3 < -1.0) c3 = -1.0;
        theta = std::acos(c3) / 3.0;
    }
}
inline double equivStrainGeneral(const double sig_pr[3], const Params& mp)
{
    const double eps0 = mp.ft / mp.E;
    const double sv[6] = { sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0 };
    double xi, rho, theta; old::invariants(sv, xi, rho, theta);
    const double sigV = xi / SQRT3;
    const double A = rho * lodeR(theta, mp.e) / (SQRT6 * mp.fc) + sigV / mp.fc;
    double rad = (eps0 * eps0 * mp.m0 * mp.m0 / 4.0) * A * A
               + 3.0 * eps0 * eps0 * rho * rho / (2.0 * mp.fc * mp.fc);
    if (rad < 0.0) rad = 0.0;
    return (eps0 * mp.m0 / 2.0) * A + std::sqrt(rad);
}
inline double alphaCompression(const double sig_pr[3])
{
    double nrm2 = 0.0, num = 0.0;
    for (int i = 0; i < 3; ++i) {
        nrm2 += sig_pr[i] * sig_pr[i];
        const double spc = sig_pr[i] < 0.0 ? sig_pr[i] : 0.0;   // negative part
        num += spc * sig_pr[i];
    }
    if (nrm2 <= 1.0e-300) return 0.0;
    return num / nrm2;
}
inline double tcPhi(const Params& mp, double s)
{
    const double v = 1.0 - (s - TC_DEAD * mp.ft) / ((TC_BAND - TC_DEAD) * mp.ft);
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}
inline void damageDrivers(const double sig_pr[3], const Params& mp, double& et, double& ac, double& xs)
{
    et = old::equivStrainGeneral(sig_pr, mp);
    ac = old::alphaCompression(sig_pr);
    const double sv[6] = { sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0 };
    double xi, rho, theta; old::invariants(sv, xi, rho, theta);
    const double sigV = xi / SQRT3;
    const double Rs = (sigV <= 0.0 && rho > 1.0e-12) ? (-SQRT6 * sigV / rho) : 0.0;   // Eq.57
    xs = 1.0 + (mp.As - 1.0) * Rs;                                                     // Eq.56
}
inline double compressiveEquivStrain(const Params& mp, const double w_stress[3], double et)
{
    if (mp.tcTemper != 2) return et;
    double mx = w_stress[0]; for (int a = 1; a < 3; ++a) if (w_stress[a] > mx) mx = w_stress[a];
    if (mx <= TC_DEAD * mp.ft) return et;
    double wc[3];
    for (int a = 0; a < 3; ++a) wc[a] = w_stress[a] < 0.0 ? w_stress[a] : old::tcPhi(mp, w_stress[a]) * w_stress[a];
    return old::equivStrainGeneral(wc, mp);
}
inline double scalarDriver(int which, const double sig6[6], const Params& mp)
{
    double A[3][3], w[3], V[3][3]; voigtToMat(sig6, A); eig3sym(A, w, V);
    if (which == 0) return old::equivStrainGeneral(w, mp);
    if (which == 1) { double et, ac, xs; old::damageDrivers(w, mp, et, ac, xs); return xs; }
    if (which == 3) return old::compressiveEquivStrain(mp, w, old::equivStrainGeneral(w, mp));   // PV20 tcTemper
    return old::alphaCompression(w);
}
inline void dscalarDsig(int which, const double sig6[6], const Params& mp, double g[6], double h)
{
    for (int k = 0; k < 6; ++k) {
        double dp[6], dm[6];
        for (int i = 0; i < 6; ++i) { dp[i] = sig6[i]; dm[i] = sig6[i]; }
        dp[k] += h; dm[k] -= h;
        g[k] = (old::scalarDriver(which, dp, mp) - old::scalarDriver(which, dm, mp)) / (2.0 * h);
    }
}
} // namespace old

// ===========================================================================
// states and parameters
// ===========================================================================
static K::Params baseParams(double s, int tc)
{
    K::Params q; q.E = 30000 * s; q.nu = 0.2; q.fc = 30 * s; q.ft = 3 * s;
    q.e = K::eccentricityFromKupfer(30 * s, 3 * s, 1.16); q.m0 = K::m0Of(30 * s, 3 * s, q.e);
    q.As = 2.0; q.tcTemper = tc;
    return q;
}
static K::Params scaleParams(const K::Params& p, double s)   // stresses and moduli x s (e, m0, As are dimensionless)
{
    K::Params q = p; q.E *= s; q.fc *= s; q.ft *= s; q.Gf *= s; q.Gc *= s; return q;
}
static void randRot(std::mt19937_64& rng, double R[3][3])   // a uniformly random rotation (normalized quaternion)
{
    std::normal_distribution<double> n(0.0, 1.0);
    double q0 = n(rng), q1 = n(rng), q2 = n(rng), q3 = n(rng);
    const double m = std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3); q0 /= m; q1 /= m; q2 /= m; q3 /= m;
    R[0][0] = 1 - 2 * (q2 * q2 + q3 * q3); R[0][1] = 2 * (q1 * q2 - q0 * q3); R[0][2] = 2 * (q1 * q3 + q0 * q2);
    R[1][0] = 2 * (q1 * q2 + q0 * q3); R[1][1] = 1 - 2 * (q1 * q1 + q3 * q3); R[1][2] = 2 * (q2 * q3 - q0 * q1);
    R[2][0] = 2 * (q1 * q3 - q0 * q2); R[2][1] = 2 * (q2 * q3 + q0 * q1); R[2][2] = 1 - 2 * (q1 * q1 + q2 * q2);
}
static void rotated(const double w[3], const double R[3][3], double s6[6])   // R diag(w) R^T in Voigt
{
    double M[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        double v = 0.0; for (int a = 0; a < 3; ++a) v += R[i][a] * w[a] * R[j][a]; M[i][j] = v;
    }
    K::matToVoigt(M, s6);
}
struct St { double s6[6]; int fam; };   // fam: 0 random, 1 exactly repeated pair, 2 near-repeated, 3 record
static const char* famName[4] = {"random", "repeated pair", "near-repeated", "record"};

// synthetic states in MPa-like units (scaled by s afterwards); every one off the hydrostatic axis
static std::vector<St> synthStates(std::mt19937_64& rng, int nRand)
{
    std::vector<St> v; std::uniform_real_distribution<double> u(-1.0, 1.0);
    double R[3][3];
    for (int i = 0; i < nRand; ++i) {   // principal values spanning tension (to ~ft) and compression (to ~3 fc)
        double w[3] = {3.0 * u(rng), 10.0 * u(rng) - 5.0, 45.0 * u(rng) - 45.0};
        randRot(rng, R); St s; rotated(w, R, s.s6); s.fam = 0; v.push_back(s);
    }
    // exactly repeated pairs (the meridians): uniaxial tension/compression, equal-biaxial, confined; rotated, so
    // eig3sym returns an arbitrary basis of the repeated eigenspace
    const double rep[][3] = {{2.0, 0, 0}, {-30.0, 0, 0}, {1.5, 1.5, 0}, {-20.0, -20.0, 0}, {-40.0, -5.0, -5.0},
                             {-10.0, -10.0, -60.0}, {1.0, 1.0, -30.0}, {2.5, -1.0, -1.0}, {-3.0, -3.0, 1.0},
                             {0.5, -25.0, -25.0}};
    for (const auto& w : rep) for (int r = 0; r < 4; ++r) {
        St s; if (r == 0) { for (int k = 0; k < 6; ++k) s.s6[k] = (k < 3) ? w[k] : 0.0; }
        else { randRot(rng, R); rotated(w, R, s.s6); }
        s.fam = 1; v.push_back(s);
    }
    // near-repeated pairs: w_a - w_b = delta * |w|, delta = 1e-4 .. 1e-12
    for (const auto& w0 : rep) for (double dl : {1e-4, 1e-8, 1e-12}) {
        double w[3] = {w0[0], w0[1], w0[2]};
        const double mag = std::fmax(std::fabs(w[0]), std::fmax(std::fabs(w[1]), std::fabs(w[2])));
        if (w[0] == w[1]) w[1] += dl * mag; else w[2] += dl * mag;
        randRot(rng, R); St s; rotated(w, R, s.s6); s.fam = 2; v.push_back(s);
    }
    return v;
}

static bool onAxis(const double s6[6], const K::Params& mp)
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    double xi, rho, th; const double sv[6] = {w[0], w[1], w[2], 0, 0, 0}; K::invariants(sv, xi, rho, th);
    return rho <= K::AXIS_RHO_REL * K::driverStressScale(w, mp);
}
static double normInf(const double g[6]) { double m = 0; for (int k = 0; k < 6; ++k) m = std::fmax(m, std::fabs(g[k])); return m; }
static double diffInf(const double a[6], const double b[6]) { double m = 0; for (int k = 0; k < 6; ++k) m = std::fmax(m, std::fabs(a[k] - b[k])); return m; }

static void dualGrad(int which, const double s6[6], const K::Params& mp, double g[6])
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    K::dscalarDsigDual(which, w, V, mp, g);
}
// the gradient's natural scale (the floor of the relative error): |d eps_tilde/d sig| ~ 1/E, the degree-0 drivers ~ 1/S
static double gradScale(int which, const double s6[6], const K::Params& mp)
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    return (which == 0 || which == 3) ? 1.0 / mp.E : 1.0 / K::driverStressScale(w, mp);
}
// step-plateau central FD: h = S 10^-k, k = 2..9; the plateau is the k whose g(h_k), g(h_k+1) agree best
static bool plateauFD(int which, const double s6[6], const K::Params& mp, double ref[6], double& quality)
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    const double S = K::driverStressScale(w, mp), sc = gradScale(which, s6, mp);
    double g[9][6];
    for (int k = 2; k <= 10; ++k) K::dscalarDsig(which, s6, mp, g[k - 2], S * std::pow(10.0, -k));
    int best = 0; quality = 1e300;
    for (int k = 0; k < 8; ++k) {
        const double q = diffInf(g[k], g[k + 1]) / std::fmax(normInf(g[k + 1]), 1e-6 * sc);
        if (q < quality) { quality = q; best = k; }
    }
    for (int i = 0; i < 6; ++i) ref[i] = g[best + 1][i];
    return quality < 1e-6;
}
static double relErr(int which, const double a[6], const double ref[6], const double s6[6], const K::Params& mp)
{
    return diffInf(a, ref) / std::fmax(normInf(ref), 1e-6 * gradScale(which, s6, mp));
}

struct Dist { std::vector<double> v; void add(double x) { v.push_back(x); }
    double q(double f) { if (v.empty()) return 0; std::vector<double> s = v; std::sort(s.begin(), s.end());
        return s[(size_t)std::min<double>(s.size() - 1, std::floor(f * (s.size() - 1) + 0.5))]; }
    double mx() { double m = 0; for (double x : v) m = std::fmax(m, x); return m; } };

// ===========================================================================
// record states: every record's committed effective stress X and, run through returnMap at every pass strain, the
// converged effective stress; with each record's own Params (MPa). And the inputs of every pass, for part T.
// ===========================================================================
struct RecPass { K::Params mp; K::State in; double strain[6]; double dt; };
static void recordInputs(const Fixture& fx, std::vector<RecPass>& passes, std::vector<std::pair<St, K::Params>>& states)
{
    for (const Rec& r : fx.recs) {
        K::Params mp; double dt = 0; fillParams(r, mp, dt); K::State in; fillState(r, in);
        St x; for (int i = 0; i < 6; ++i) x.s6[i] = r.sigEff[i]; x.fam = 3; states.push_back({x, mp});
        std::vector<double> e22;
        if (r.kind == 'R') { e22.push_back(r.c[1]); for (const Pass& p : r.z) e22.push_back(p.eps22); }
        else e22.push_back(r.strain[2]);
        for (double e : e22) {
            RecPass p; p.mp = mp; p.in = in; p.dt = dt; for (int i = 0; i < 6; ++i) p.strain[i] = r.strain[i]; p.strain[2] = e;
            passes.push_back(p);
            K::State o; double sg[6], se[6], D[6][6];
            K::returnMap(mp, p.strain, in, o, sg, se, D, false, dt, true);
            St y; for (int i = 0; i < 6; ++i) y.s6[i] = se[i]; y.fam = 3; states.push_back({y, mp});
        }
    }
}

// ===========================================================================
// B: bit-identity
// ===========================================================================
static void runBitIdentity(const std::vector<std::pair<St, K::Params>>& all)
{
    std::printf("[B] bit-identity of the templated double entry points vs the pre-WP-264 code (verbatim)\n");
    long nCmp = 0, nBad = 0, nDualBad = 0;
    auto cmp = [&](double a, double b) { ++nCmp; if (!bitEq(a, b)) ++nBad; };
    for (const auto& sp : all) {
        const double* s6 = sp.first.s6; const K::Params& mp = sp.second;
        { double a1, a2, a3, b1, b2, b3; K::invariants(s6, a1, a2, a3); old::invariants(s6, b1, b2, b3); cmp(a1, b1); cmp(a2, b2); cmp(a3, b3); }
        double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
        cmp(K::equivStrainGeneral(w, mp), old::equivStrainGeneral(w, mp));
        cmp(K::alphaCompression(w), old::alphaCompression(w));
        for (int a = 0; a < 3; ++a) cmp(K::tcPhi(mp, w[a]), old::tcPhi(mp, w[a]));
        { double a1, a2, a3, b1, b2, b3; K::damageDrivers(w, mp, a1, a2, a3); old::damageDrivers(w, mp, b1, b2, b3); cmp(a1, b1); cmp(a2, b2); cmp(a3, b3); }
        const double et = K::equivStrainGeneral(w, mp);
        cmp(K::compressiveEquivStrain(mp, w, et), old::compressiveEquivStrain(mp, w, et));
        const double h = K::DRIVER_FD_REL * K::driverStressScale(w, mp);
        for (int which = 0; which < 4; ++which) {
            cmp(K::scalarDriver(which, s6, mp), old::scalarDriver(which, s6, mp));
            double g1[6], g2[6]; K::dscalarDsig(which, s6, mp, g1, h); old::dscalarDsig(which, s6, mp, g2, h);
            for (int k = 0; k < 6; ++k) cmp(g1[k], g2[k]);
            // the value of the Dual pass is the double result, bit for bit
            const Dual<3> wd[3] = {Dual<3>::variable(w[0], 0), Dual<3>::variable(w[1], 1), Dual<3>::variable(w[2], 2)};
            const Dual<3> f = K::scalarDriverOfPrincipals(which, wd, mp);
            if (!bitEq(f.v, K::scalarDriverOfPrincipals(which, w, mp))) ++nDualBad;
        }
    }
    std::printf("       (%zu states x 3 unit systems x tcTemper 0/2 folded in; %ld double compares, %ld differ; dual value mismatches %ld)\n",
                all.size(), nCmp, nBad, nDualBad);
    gate(nBad == 0 && nCmp > 10000, "B1 the double entry points are bit-identical to the pre-WP-264 code");
    gate(nDualBad == 0, "B2 value(f(Dual<3>)) == f(double) bit for bit for every driver");
}

// ===========================================================================
// G: gradient oracle
// ===========================================================================
// an exact, eigen-free cross-check for eps_tilde (which 0) and x_s (which 1): Dual<6> seeded on the full tensor
template <class T>
static T eigenFreeDriver(int which, const T s6[6], const K::Params& mp)
{
    using std::sqrt;
    T xi, rho, c3; double theta; K::invariantsT(s6, xi, rho, c3, theta);
    const T sigV = xi / K::SQRT3;
    if (which == 1) {
        const T Rs = (sigV <= 0.0 && rho > 1.0e-12) ? T(-K::SQRT6 * sigV / rho) : T(0.0);
        return 1.0 + (mp.As - 1.0) * Rs;
    }
    const double eps0 = mp.ft / mp.E;
    const T A = rho * K::lodeRT(c3, theta, mp.e) / (K::SQRT6 * mp.fc) + sigV / mp.fc;
    T rad = (eps0 * eps0 * mp.m0 * mp.m0 / 4.0) * A * A + 3.0 * eps0 * eps0 * rho * rho / (2.0 * mp.fc * mp.fc);
    if (rad < 0.0) rad = T(0.0);
    return (eps0 * mp.m0 / 2.0) * A + sqrt(rad);
}
static void eigenFreeGrad(int which, const double s6[6], const K::Params& mp, double g[6])
{
    Dual<6> x[6]; for (int k = 0; k < 6; ++k) x[k] = Dual<6>::variable(s6[k], k);
    const Dual<6> f = eigenFreeDriver(which, x, mp);
    for (int k = 0; k < 6; ++k) g[k] = f.d[k];   // per Voigt component: d/d sig6[k] moves both shear entries already
}

// The distance (relative to S) from the nearest KINK of driver `which` off the axis, where the derivative jumps and a
// central FD whose stencil straddles it averages two one-sided slopes: alpha_c at w_a = 0 (the negative part); x_s at
// sigV = 0 (Rs switches on); eps_tilde_c at the two tcPhi corners (TC_DEAD ft, TC_BAND ft) and at the VERTEX of its
// argument wc (rho(wc) == 0: e.g. every tensile principal past TC_BAND ft and the others zero, plane-stress tension).
// eps_tilde has none off the axis. A kink-adjacent state (distance < 1e-3) is reported, not gated against the FD.
static double kinkDistance(int which, const double s6[6], const K::Params& mp)
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    const double S = K::driverStressScale(w, mp);
    double d = 1e300;
    if (which == 1) d = std::fabs(w[0] + w[1] + w[2]) / 3.0;
    if (which == 2) for (int a = 0; a < 3; ++a) d = std::fmin(d, std::fabs(w[a]));
    if (which == 3) {
        for (int a = 0; a < 3; ++a)
            d = std::fmin(d, std::fmin(std::fabs(w[a] - K::TC_DEAD * mp.ft), std::fabs(w[a] - K::TC_BAND * mp.ft)));
        double wc[3];
        for (int a = 0; a < 3; ++a) wc[a] = w[a] < 0.0 ? w[a] : K::tcPhi(mp, w[a]) * w[a];
        double xi, rho, th; const double sv[6] = {wc[0], wc[1], wc[2], 0, 0, 0}; K::invariants(sv, xi, rho, th);
        d = std::fmin(d, rho);
    }
    return d / S;
}

struct OracleRow { Dist plateau, shipped, kinkPlateau; int n = 0, skipped = 0, kink = 0; };
static void runGradientOracle(const std::vector<std::pair<St, K::Params>>& all, const char* unitTag)
{
    OracleRow rows[4][4];   // [family][driver]
    int nonFinite = 0;
    Dist basisDist, eigenFree;
    std::mt19937_64 rng(264);
    for (const auto& sp : all) {
        const double* s6 = sp.first.s6; const K::Params& mp = sp.second; const int fam = sp.first.fam;
        if (onAxis(s6, mp)) continue;   // the axis keeps its dedicated FD (dscalarDsigAxis)
        for (int which = 0; which < 4; ++which) {
            if (which == 3 && mp.tcTemper != 2) continue;
            double gd[6], gs[6], gp[6], q;
            dualGrad(which, s6, mp, gd);
            bool fin = true; for (int k = 0; k < 6; ++k) fin = fin && std::isfinite(gd[k]);
            auto dbg = [&](const char* what, double e) {
                if (!std::getenv("LADRUNO_DUAL_DEBUG")) return;
                double A2[3][3], w2[3], V2[3][3]; K::voigtToMat(s6, A2); K::eig3sym(A2, w2, V2);
                std::printf("DBG %s fam=%d which=%d tc=%d err=%.3e w=(%.17g %.17g %.17g) ft=%g\n", what, fam, which, mp.tcTemper, e,
                            w2[0], w2[1], w2[2], mp.ft);
            };
            if (!fin) { dbg("nonfinite", 0); ++nonFinite; continue; }
            double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
            K::dscalarDsig(which, s6, mp, gs, K::DRIVER_FD_REL * K::driverStressScale(w, mp));
            OracleRow& r = rows[fam][which];
            const bool kink = kinkDistance(which, s6, mp) < 1e-3;
            if (kink) {
                ++r.kink;
                if (plateauFD(which, s6, mp, gp, q)) r.kinkPlateau.add(relErr(which, gd, gp, s6, mp));
            } else {
                r.shipped.add(relErr(which, gd, gs, s6, mp));
                if (plateauFD(which, s6, mp, gp, q)) {
                    const double e = relErr(which, gd, gp, s6, mp); if (e > 1e-6) dbg("plateau", e); r.plateau.add(e); ++r.n;
                } else { dbg("noplateau", q); ++r.skipped; }
            }
            // exact cross-checks (kink-adjacent states included)
            if (which <= 1) {
                double ge[6]; eigenFreeGrad(which, s6, mp, ge);
                const double e = relErr(which, gd, ge, s6, mp); if (e > 1e-9) dbg("eigenfree", e); eigenFree.add(e);
            }
            if (fam == 1 && !kink) {   // rotate the eigenbasis inside the repeated pair: the gradient must not move
                int a = -1, b = -1;
                for (int i = 0; i < 3 && a < 0; ++i) for (int j = i + 1; j < 3; ++j)
                    if (std::fabs(w[i] - w[j]) <= 1e-12 * K::driverStressScale(w, mp)) { a = i; b = j; break; }
                if (a >= 0) {
                    const double phi = std::uniform_real_distribution<double>(0.1, 3.0)(rng), c = std::cos(phi), s = std::sin(phi);
                    double V2[3][3]; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) V2[i][j] = V[i][j];
                    for (int i = 0; i < 3; ++i) { V2[i][a] = c * V[i][a] - s * V[i][b]; V2[i][b] = s * V[i][a] + c * V[i][b]; }
                    double g2[6]; K::dscalarDsigDual(which, w, V2, mp, g2);
                    const double e = relErr(which, g2, gd, s6, mp); if (e > 1e-9) dbg("basis", e); basisDist.add(e);
                }
            }
        }
    }
    static const char* dn[4] = {"eps_tilde", "x_s", "alpha_c", "eps_tilde_c"};
    std::printf("       units %s: dual vs plateau FD (p50 / p90 / max, n, no-plateau) | vs shipped FD 1e-6 S (p50 / p90 / max)"
                " | kink-adjacent (n, max vs plateau)\n", unitTag);
    double worst = 0; int nTot = 0, skipTot = 0;
    for (int fam = 0; fam < 4; ++fam) for (int which = 0; which < 4; ++which) {
        OracleRow& r = rows[fam][which];
        if (r.shipped.v.empty() && r.kink == 0) continue;
        std::printf("       %-13s %-11s %.1e / %.1e / %.1e  n=%3d skip=%d | %.1e / %.1e / %.1e | kink n=%3d max %.1e\n",
                    famName[fam], dn[which], r.plateau.q(0.5), r.plateau.q(0.9), r.plateau.mx(), r.n, r.skipped,
                    r.shipped.q(0.5), r.shipped.q(0.9), r.shipped.mx(), r.kink, r.kinkPlateau.mx());
        worst = std::fmax(worst, r.plateau.mx()); nTot += r.n; skipTot += r.skipped;
    }
    std::printf("       exact cross-checks: eigen-free Dual<6> (eps_tilde, x_s) max %.1e (n %zu); basis rotation inside a repeated"
                " pair max %.1e (n %zu)\n", eigenFree.mx(), eigenFree.v.size(), basisDist.mx(), basisDist.v.size());
    gate(nonFinite == 0, std::string("G1 [") + unitTag + "] the dual gradient is finite on every off-axis state, repeated pairs and kinks included");
    gate(nTot > 200 && worst < 1e-6 && skipTot * 50 <= nTot,
         std::string("G2 [") + unitTag + "] dual == plateau FD to 1e-6 off the kinks (no-plateau <= 2 %)");
    gate(eigenFree.mx() < 1e-9 && basisDist.mx() < 1e-9 && basisDist.v.size() >= 40,
         std::string("G3 [") + unitTag + "] dual == eigen-free Dual<6> and basis-independent inside repeated pairs to 1e-9");
}

// ===========================================================================
// M: break-the-gate mutants (each must be flagged on the same oracle set)
// ===========================================================================
struct DualAcos {   // the naive rule: theta = acos(c3)/3 pushed through a dual (d acos = -dx / sqrt(1 - x^2))
    template <int N> static Dual<N> acosD(const Dual<N>& x)
    {
        Dual<N> r(std::acos(x.v)); const double k = -1.0 / std::sqrt(1.0 - x.v * x.v);
        for (int i = 0; i < N; ++i) r.d[i] = k * x.d[i];
        return r;
    }
};
template <int MUT>
static Dual<3> mutEquivStrain(const Dual<3> sp[3], const K::Params& mp)
{
    typedef Dual<3> T; using std::sqrt;
    const double eps0 = mp.ft / mp.E;
    const T sv[6] = {sp[0], sp[1], sp[2], T(0.0), T(0.0), T(0.0)};
    T xi, rho, c3; double theta; K::invariantsT(sv, xi, rho, c3, theta);
    T r;
    if (MUT == 2) r = T(K::lodeR(theta, mp.e));                              // M2: the Lode derivative dropped
    else {                                                                    // M4: d theta through acos (0 * inf)
        const T th = DualAcos::acosD(c3) / 3.0;
        r = T(K::lodeR(theta, mp.e));
        const double h = 1e-7, drdth = (K::lodeR(theta + h, mp.e) - K::lodeR(theta - h, mp.e)) / (2 * h);
        for (int i = 0; i < 3; ++i) r.d[i] = drdth * th.d[i];
    }
    const T sigV = xi / K::SQRT3;
    const T A = rho * r / (K::SQRT6 * mp.fc) + sigV / mp.fc;
    T rad = (eps0 * eps0 * mp.m0 * mp.m0 / 4.0) * A * A + 3.0 * eps0 * eps0 * rho * rho / (2.0 * mp.fc * mp.fc);
    if (rad < 0.0) rad = T(0.0);
    return (eps0 * mp.m0 / 2.0) * A + sqrt(rad);
}
static void mutantGrad(int mut, const double s6[6], const K::Params& mp, double g[6])
{
    double A[3][3], w[3], V[3][3]; K::voigtToMat(s6, A); K::eig3sym(A, w, V);
    const Dual<3> wd[3] = {Dual<3>::variable(w[0], 0), Dual<3>::variable(w[1], 1), Dual<3>::variable(w[2], 2)};
    Dual<3> f;
    if (mut == 2) f = mutEquivStrain<2>(wd, mp);
    else if (mut == 4) f = mutEquivStrain<4>(wd, mp);
    else f = K::scalarDriverOfPrincipals(0, wd, mp);
    static const int I6[6] = {0, 1, 2, 0, 1, 0}, J6[6] = {0, 1, 2, 1, 2, 2};
    for (int k = 0; k < 6; ++k) {
        double v = 0.0;
        for (int a = 0; a < 3; ++a) {
            const int b = (mut == 3) ? (a + 1) % 3 : a;                       // M3: pairing scrambled
            v += f.d[a] * V[I6[k]][b] * V[J6[k]][b];
        }
        g[k] = (k < 3 || mut == 1) ? v : 2.0 * v;                             // M1: the shear factor 2 dropped
    }
}
static void runMutants(const std::vector<std::pair<St, K::Params>>& all)
{
    std::printf("[M] break-the-gate: each mutant of the eps_tilde gradient must fail the plateau-FD gate (1e-6)\n");
    static const char* mn[5] = {"", "M1 shear factor 2 dropped", "M2 Lode derivative dropped",
                                "M3 eigenvalue/eigenvector pairing scrambled", "M4 naive dual through acos (meridians)"};
    for (int mut = 1; mut <= 4; ++mut) {
        double worst = 0; int nBad = 0, nNonFinite = 0, n = 0;
        for (const auto& sp : all) {
            const double* s6 = sp.first.s6; const K::Params& mp = sp.second;
            if (onAxis(s6, mp)) continue;
            if (mut == 4 && sp.first.fam != 1) continue;                      // M4 is the meridian hazard: repeated pairs
            double gm[6], gp[6], q; mutantGrad(mut, s6, mp, gm);
            if (!plateauFD(0, s6, mp, gp, q)) continue;
            ++n;
            bool fin = true; for (int k = 0; k < 6; ++k) fin = fin && std::isfinite(gm[k]);
            if (!fin) { ++nNonFinite; ++nBad; continue; }
            const double e = relErr(0, gm, gp, s6, mp);
            worst = std::fmax(worst, e); if (e > 1e-6) ++nBad;
        }
        std::printf("       %-44s flagged on %d of %d states (non-finite %d), worst finite %.2e\n", mn[mut], nBad, n, nNonFinite, worst);
        gate(nBad > 0, std::string(mn[mut]) + ": flagged");
    }
}

// ===========================================================================
// T: the damaged tangent with dual vs FD drivers
// ===========================================================================
static double relF6(const double A[6][6], const double R[6][6])
{
    double n = 0, d = 0;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { n += (A[i][j] - R[i][j]) * (A[i][j] - R[i][j]); d += R[i][j] * R[i][j]; }
    return d > 0 ? std::sqrt(n / d) : std::sqrt(n);
}
struct TanEval { int st; K::State o; double sg[6], se[6], D[6][6]; };
static TanEval evalMode(const K::Params& mp, const K::State& in, const double e[6], double dt, bool fd)
{
    K::driverGradFD() = fd;
    TanEval t; t.st = K::returnMap(mp, e, in, t.o, t.sg, t.se, t.D, true, dt, true);
    K::driverGradFD() = false;
    return t;
}
static bool sameStress(const TanEval& a, const TanEval& b)
{
    if (a.st != b.st) return false;
    double na[N_STATE], nb[N_STATE]; stateOut(a.o, na); stateOut(b.o, nb);   // field by field (State has padding)
    for (int i = 0; i < N_STATE; ++i) if (!bitEq(na[i], nb[i])) return false;
    for (int i = 0; i < 6; ++i) if (!bitEq(a.sg[i], b.sg[i]) || !bitEq(a.se[i], b.se[i])) return false;
    return true;
}
// The central FD of the damaged map at e (h = 1e-8, 1e-9, 1e-10; every stencil map converged on the centre's piece
// schedule), best of h, for two tangents at once. false when no h keeps every stencil point on the schedule.
static bool mapFDBest(const K::Params& q, const K::State& in, const double e[6], double dt, int sub,
                      const double DA[6][6], const double DB[6][6], double& bestA, double& bestB)
{
    bestA = bestB = 1e300;
    for (double h : {1e-8, 1e-9, 1e-10}) {
        double F[6][6]; bool ok = true;
        for (int j = 0; j < 6 && ok; ++j) {
            double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
            K::State u, v; double a1[6], b1[6], c1[6], d1[6], J[6][6];
            ok = K::returnMap(q, ep, in, u, a1, b1, J, false, dt, true) == 0 && K::returnMap(q, em, in, v, c1, d1, J, false, dt, true) == 0
                 && u.subInfo == sub && v.subInfo == sub;
            for (int i = 0; i < 6; ++i) F[i][j] = (a1[i] - c1[i]) / (2 * h);
        }
        if (ok) { bestA = std::fmin(bestA, relF6(DA, F)); bestB = std::fmin(bestB, relF6(DB, F)); }
    }
    return bestA < 1e300;
}
static void runTangentRecords(const std::vector<RecPass>& passes)
{
    std::printf("[T1] damaged tangent, dual vs FD drivers, at every replay record pass (MPa, the records' own Params)\n");
    Dist d, eA, eB; int n = 0, nSame = 0, nStressBad = 0, nRef = 0, noRef = 0, worse = 0;
    for (const RecPass& p : passes) {
        const TanEval a = evalMode(p.mp, p.in, p.strain, p.dt, false), b = evalMode(p.mp, p.in, p.strain, p.dt, true);
        if (!sameStress(a, b)) ++nStressBad;
        ++n;
        const double r = relF6(a.D, b.D);
        if (r == 0.0) { ++nSame; continue; }
        d.add(r);
        double ba, bb;   // which is right where they differ: both against the FD of the damaged map
        if (a.st == 0 && mapFDBest(p.mp, p.in, p.strain, p.dt, a.o.subInfo, a.D, b.D, ba, bb)) {
            ++nRef; eA.add(ba); eB.add(bb);
            if (ba > std::fmax(2.0 * bb, 1e-6)) ++worse;
        } else ++noRef;
    }
    std::printf("       passes %d: tangent bit-identical %d (no interior omega, or on the axis), differing %zu: p50 %.2e p90 %.2e max %.2e;"
                " stress/State differ %d\n", n, nSame, d.v.size(), d.q(0.5), d.q(0.9), d.mx(), nStressBad);
    std::printf("       vs the FD of the damaged map (%d with a reference, %d without: refused or off-schedule stencil):"
                " dual p50 %.2e p90 %.2e max %.2e | FD drivers p50 %.2e p90 %.2e max %.2e\n",
                nRef, noRef, eA.q(0.5), eA.q(0.9), eA.mx(), eB.q(0.5), eB.q(0.9), eB.mx());
    gate(nStressBad == 0, "T1a the stress and the returned State are bit-identical between the two driver modes");
    gate(d.mx() < 1e-5 && worse == 0, "T1b the two tangents differ by < 1e-5, and the dual one is nowhere worse vs the map FD "
                                      "(> max(2x the FD-driver error, 1e-6))");
}
static void runTangentPaths(double s, double L, const char* tag)
{
    K::Params q; q = baseParams(s, 2);
    q.tensionLaw = 1; q.compDrive = 1; q.Gf = 0.1 * s * L; q.Gc = 5.0 * s * L; q.lch = 50.0 * L; q.epsFc = 1e-3; q.As = 2.0;
    q.flowPotential = 1; q.maxSubIncr = 10;
    const char* names[3] = {"far compressive", "compression + lateral expansion", "tension + shear"};
    std::printf("[T2] damaged tangent on damage-advancing paths, units %s (tcTemper proj): vs the FD of the damaged map (best h)\n", tag);
    for (int path = 0; path < 3; ++path) {
        K::State in; Dist dd, df, dx; int skipped = 0, nStressBad = 0;
        for (int k = 1; k <= 120; ++k) {
            const double t = k / 120.0; double e[6] = {0, 0, 0, 0, 0, 0};
            if (path == 0)      { const double M = 1.5e-1; e[0] = -M * t; e[3] = 0.04 * M * t; e[4] = 0.02 * M * t; }
            else if (path == 1) { const double M = 3e-2; e[0] = -M * t; e[1] = 0.5 * M * t; e[2] = 0.4 * M * t; e[3] = 0.1 * M * t; }
            else                { e[0] = 1.5e-3 * t; e[1] = -3e-4 * t; e[3] = 6e-4 * t; e[5] = 2e-4 * t; }
            const TanEval a = evalMode(q, in, e, -1.0, false), b = evalMode(q, in, e, -1.0, true);
            if (!sameStress(a, b)) ++nStressBad;
            if (a.st == 0 && ((a.o.wc > in.wc && in.wc > 0.0) || (a.o.wt > in.wt && in.wt > 0.0))) {
                double bestA, bestB;
                if (mapFDBest(q, in, e, -1.0, a.o.subInfo, a.D, b.D, bestA, bestB)) { dd.add(bestA); df.add(bestB); dx.add(relF6(a.D, b.D)); }
                else ++skipped;
            }
            in = a.o;
        }
        std::printf("       %-32s n=%3zu skip=%d | dual: p50 %.2e p90 %.2e max %.2e | FD drivers: p50 %.2e p90 %.2e max %.2e"
                    " | dual vs FD drivers: p50 %.2e p90 %.2e max %.2e\n",
                    names[path], dd.v.size(), skipped, dd.q(0.5), dd.q(0.9), dd.mx(), df.q(0.5), df.q(0.9), df.mx(), dx.q(0.5), dx.q(0.9), dx.mx());
        gate(nStressBad == 0 && skipped == 0 && dd.v.size() >= 20 && dd.mx() < 1e-6,
             std::string("T2 ") + names[path] + " [" + tag + "]: dual-driver damaged tangent == FD of the damaged map to 1e-6, stress bit-identical");
    }
}

// ===========================================================================
// P: timing (g++ -O2, no -march=native)
// ===========================================================================
static volatile double gSink = 0.0;
static void runTiming(const std::vector<std::pair<St, K::Params>>& all, const std::vector<RecPass>& passes)
{
    std::printf("[P] per-call time, g++ -O2 (no -march=native), off-axis states of the oracle set\n");
    std::vector<const std::pair<St, K::Params>*> st;
    for (const auto& sp : all) if (!onAxis(sp.first.s6, sp.second)) st.push_back(&sp);
    typedef std::chrono::steady_clock C;
    auto perCall = [&](auto&& fn, long reps) {
        const auto t0 = C::now();
        for (long r = 0; r < reps; ++r) for (const auto* sp : st) fn(*sp);
        return std::chrono::duration<double, std::nano>(C::now() - t0).count() / (double)(reps * (long)st.size());
    };
    static const char* dn[4] = {"eps_tilde", "x_s", "alpha_c", "eps_tilde_c"};
    std::printf("       %-12s %12s %12s %14s %8s\n", "driver", "FD [ns]", "dual [ns]", "dual+eig [ns]", "FD/dual");
    double sumFD = 0, sumD = 0;
    for (int which = 0; which < 4; ++which) {
        const long reps = 20;
        auto fdFn = [&](const std::pair<St, K::Params>& sp) {
            double A[3][3], w[3], V[3][3]; K::voigtToMat(sp.first.s6, A); K::eig3sym(A, w, V);   // driverStressScale needs w
            double g[6]; K::dscalarDsig(which, sp.first.s6, sp.second, g, K::DRIVER_FD_REL * K::driverStressScale(w, sp.second)); gSink += g[0];
        };
        // damagedTangent already holds w and V (its own eig3sym of sig_bar): the dual path is the call alone
        std::vector<std::array<double, 12>> wv; wv.reserve(st.size());
        for (const auto* sp : st) { double A[3][3], w[3], V[3][3]; K::voigtToMat(sp->first.s6, A); K::eig3sym(A, w, V);
            std::array<double, 12> x; for (int a = 0; a < 3; ++a) x[a] = w[a]; for (int i = 0; i < 9; ++i) x[3 + i] = V[i / 3][i % 3]; wv.push_back(x); }
        size_t idx = 0;
        auto dualFn = [&](const std::pair<St, K::Params>& sp) {
            const std::array<double, 12>& x = wv[idx++ % wv.size()];
            const double w[3] = {x[0], x[1], x[2]}; double V[3][3]; for (int i = 0; i < 9; ++i) V[i / 3][i % 3] = x[3 + i];
            double g[6]; K::dscalarDsigDual(which, w, V, sp.second, g); gSink += g[0];
        };
        auto dualEigFn = [&](const std::pair<St, K::Params>& sp) { double g[6]; dualGrad(which, sp.first.s6, sp.second, g); gSink += g[0]; };
        const double tFD = perCall(fdFn, reps), tEig = perCall(dualEigFn, reps);
        idx = 0; const double tD = perCall(dualFn, reps);
        std::printf("       %-12s %12.1f %12.1f %14.1f %8.1fx\n", dn[which], tFD, tD, tEig, tFD / tD);
        if (which < 3) { sumFD += tFD; sumD += tD; }
    }
    std::printf("       the three gradients every needDw call takes (eps_tilde, x_s, alpha_c): FD %.0f ns, dual %.0f ns (%.1fx)\n",
                sumFD, sumD, sumFD / sumD);
    // returnMap with and without the tangent, in both modes: the end-to-end cost of one material call. Two sets: the
    // damage-advancing states of the three [T2] paths in MPa (a typical damaged step: the gradients run every call) and
    // the record passes (hard states, many sub-increments, the return map dominates). Modes interleaved, best of 5.
    std::vector<RecPass> dmg;
    {
        K::Params q = baseParams(1.0, 2);
        q.tensionLaw = 1; q.compDrive = 1; q.Gf = 0.1; q.Gc = 5.0; q.lch = 50.0; q.epsFc = 1e-3; q.As = 2.0; q.flowPotential = 1; q.maxSubIncr = 10;
        for (int path = 0; path < 3; ++path) {
            K::State in;
            for (int k = 1; k <= 120; ++k) {
                const double t = k / 120.0; double e[6] = {0, 0, 0, 0, 0, 0};
                if (path == 0)      { const double M = 1.5e-1; e[0] = -M * t; e[3] = 0.04 * M * t; e[4] = 0.02 * M * t; }
                else if (path == 1) { const double M = 3e-2; e[0] = -M * t; e[1] = 0.5 * M * t; e[2] = 0.4 * M * t; e[3] = 0.1 * M * t; }
                else                { e[0] = 1.5e-3 * t; e[1] = -3e-4 * t; e[3] = 6e-4 * t; e[5] = 2e-4 * t; }
                K::State o; double sg[6], se[6], D[6][6];
                const int st = K::returnMap(q, e, in, o, sg, se, D, false, -1.0, true);
                if (st == 0 && ((o.wc > in.wc && in.wc > 0.0) || (o.wt > in.wt && in.wt > 0.0))) {
                    RecPass p; p.mp = q; p.in = in; p.dt = -1.0; for (int i = 0; i < 6; ++i) p.strain[i] = e[i]; dmg.push_back(p);
                }
                in = o;
            }
        }
    }
    auto rmTime = [&](const std::vector<RecPass>& set, bool fd, bool tan, long reps) {
        K::driverGradFD() = fd;
        const auto t0 = C::now();
        for (long r = 0; r < reps; ++r) for (const RecPass& p : set) {
            K::State o; double sg[6], se[6], D[6][6]; K::returnMap(p.mp, p.strain, p.in, o, sg, se, D, tan, p.dt, true); gSink += D[0][0] * tan + sg[0];
        }
        K::driverGradFD() = false;
        return std::chrono::duration<double, std::micro>(C::now() - t0).count() / (double)(reps * (long)set.size());
    };
    for (int s = 0; s < 2; ++s) {
        const std::vector<RecPass>& set = s == 0 ? dmg : passes;
        const long reps = s == 0 ? 40 : 6;
        double tNo = 1e300, tFD = 1e300, tD = 1e300;
        for (int k = 0; k < 5; ++k) {
            tNo = std::fmin(tNo, rmTime(set, false, false, reps));
            tFD = std::fmin(tFD, rmTime(set, true, true, reps));
            tD = std::fmin(tD, rmTime(set, false, true, reps));
        }
        std::printf("       returnMap per call, %-28s (%3zu): no tangent %8.2f us | tangent, FD drivers %8.2f us | dual drivers %8.2f us"
                    " | the tangent's cost FD %.2f -> dual %.2f us\n", s == 0 ? "damage-advancing [T2] states" : "record passes",
                    set.size(), tNo, tFD, tD, tFD - tNo, tD - tNo);
    }
}

int main(int argc, char** argv)
{
    const char* fixPath = argc > 1 ? argv[1] : "tests/_testbed/concrete3d_replay_records.txt";
    const Fixture fx = load(fixPath);
    std::vector<RecPass> passes; std::vector<std::pair<St, K::Params>> recStates;
    recordInputs(fx, passes, recStates);
    std::printf("WP-264 dual-driver check: %zu records, %zu passes, %zu record states (%s)\n",
                fx.recs.size(), passes.size(), recStates.size(), fixPath);

    std::mt19937_64 rng(1175);
    const std::vector<St> syn = synthStates(rng, 300);
    const double units[3] = {1.0, 1.0e3, 1.0e6};
    const char* tags[3] = {"MPa", "kPa", "Pa"};
    std::vector<std::pair<St, K::Params>> perUnit[3], every;
    for (int u = 0; u < 3; ++u) {
        for (const St& s0 : syn) for (int tc = 0; tc <= 2; tc += 2) {
            St s = s0; for (int k = 0; k < 6; ++k) s.s6[k] *= units[u];
            perUnit[u].push_back({s, baseParams(units[u], tc)});
        }
        for (const auto& rs : recStates) {
            St s = rs.first; for (int k = 0; k < 6; ++k) s.s6[k] *= units[u];
            perUnit[u].push_back({s, scaleParams(rs.second, units[u])});
        }
        every.insert(every.end(), perUnit[u].begin(), perUnit[u].end());
    }
    // the bit-identity set also takes on-axis and zero states (they are not differentiated, but they are evaluated)
    { K::Params p = baseParams(1.0, 2); St z; for (int k = 0; k < 6; ++k) z.s6[k] = 0.0; z.fam = 0; every.push_back({z, p});
      for (double h : {2.0, -10.0, 1e-9}) { St a; for (int k = 0; k < 6; ++k) a.s6[k] = (k < 3) ? h : 0.0; a.fam = 1; every.push_back({a, p}); } }

    runBitIdentity(every);
    std::printf("[G] gradient oracle: dscalarDsigDual vs FD, per family and driver\n");
    for (int u = 0; u < 3; ++u) runGradientOracle(perUnit[u], tags[u]);
    runMutants(perUnit[0]);
    runTangentRecords(passes);
    runTangentPaths(1.0, 1.0, "MPa");
    runTangentPaths(1.0e3, 1.0e-3, "kPa-m");
    runTangentPaths(1.0e6, 1.0e-3, "Pa-m");
    if (!std::getenv("LADRUNO_SKIP_TIMING")) runTiming(perUnit[0], passes);
    std::printf(gFails ? "\nDUAL CHECK: %d FAIL\n" : "\nDUAL CHECK: ALL PASS\n", gFails);
    return gFails ? 1 : 0;
}
