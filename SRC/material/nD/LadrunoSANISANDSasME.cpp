/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno WP-129 (TIMs F18(a), F20(b)/(c); amended by the WP-128 verdict and
// the WP-134 oracle): SAS-ME, IntScheme 129 -- a Sloan-Abbo-Sheng-style
// explicit modified Euler for ManzariDafalias / LadrunoSANISAND.
//
// These are MEMBER functions of the vanilla class ManzariDafalias (declared in
// UWmaterials/ManzariDafalias.h behind `// Ladruno WP-129`), defined here so the
// vanilla footprint is a state block, the declarations and one dispatch branch
// in integrate(). The scheme is reachable only when mLadrunoSas.allowed is
// true, i.e. from a LadrunoSANISAND instance, whose wrappers forward the
// refusal to the element (LADRUNO_MATERIAL_REFUSED) and so to analyze().
//
// THE TARGET is WP-134's `uw_model` oracle (DM04 rate equations with the UW
// constitutive additions U1-U5, the PAPER's alpha_in rule, continuous moduli,
// integrated exactly): SAS-ME is an explicit integrator of THOSE equations.
//
// WHAT IT CLOSES, relative to ModifiedEuler (IntScheme 1), item by item:
//   U9 every stage evaluates EVERYTHING at its own state: K, G (GetElasticModuli
//      at the stage stress and the stage's void ratio), n, b, d, h, D, B, C.
//      ModifiedEuler freezes K, G at the committed state over the whole
//      increment (WP-134: 0.6 / 6 / 24 % of the stress increment at
//      1e-5 / 1e-4 / 1e-3, invisible to its error test). The elastic predictor
//      and the elastic part of an intersected increment use Heun on the moduli.
//   E  the substep error measures stress, back-stress alpha AND fabric z:
//        err = max( |dS2-dS1| / max(2|S|, sigma_ref),
//                   |dA2-dA1| / max(2|A|, 1),
//                   |dZ2-dZ1| / max(2|Z|, 1),  DBL_EPSILON )
//      norms at the substep START (as ModifiedEuler), accept iff err <= TolR
//      (the deck's TolR, always). sigma_ref defaults to P_atm/101 = 1 kPa at
//      P_atm 101 -- exactly ModifiedEuler's scale (WP-128 sec. 5.1: its
//      "0.5 kPa switch" is the continuous max(2|S|, 1 kPa)); `-errFloor`
//      sets it. alpha and z are dimensionless and use the unit reference: the
//      alpha error IS a stress error in units of p (d(p alpha)/p = d alpha).
//   F/U10 each stage is classified from the elastic trial with the TRUE yield
//      gradient: N = df/dsigma : C : d_eps, df/dsigma = n - (1/3)(n:alpha +
//      sqrt(2/3) m) I (ModifiedEuler's n:r differs off the surface; the
//      predictor's unload test used n:dsigma). N <= 0: elastic stage, alpha and
//      z UNCHANGED (never re-derived from the stress ratio). N > 0 with a
//      non-positive denominator H = Kp + df/dsigma:C:R: at stage 1 (a property
//      of the accepted start state, which no cut can change) a REFUSAL; at
//      stage 2 a cut, a refusal at dT_min. Never "elastic" (WP-134: that path,
//      with the uncapped q, drives err to exactly 0 and is behind all 25
//      campaign-ME ring escapes from admissible starts).
//   q  step factor clamp(0.9 sqrt(TolR/err), 0.1, 1.1), err floored at
//      DBL_EPSILON (no q = inf jump), no growth on the substep after a
//      rejection.
//   C  refuse instead of force-accept: at dT_min, when the drift correction
//      cannot reach TolF, and when the committed START is inadmissible
//      (f > TolF; alpha/alpha^b > 1 + kappa; non-deviatoric alpha or z;
//      p + p_r <= 0; non-finite). Each refusal has a named code (sasStats).
//   drift  own consistent (Potts-Gens: sigma, alpha AND z along one lambda)
//      -then-normal correction; if neither reduces |f| it FAILS (WP-128 Q2:
//      Stress_Correction's give-up returned f > 0 as success). Stress_Correction
//      itself is not called and not changed.
//   G  alpha_in follows the PAPER's rule inside the increment (default
//      `-sasAlphaIn reseat`, = WP-134's oracle): wherever a stage finds
//      (alpha - alpha_in):n < 0 a new loading process starts there and
//      alpha_in := alpha at THAT stage; an accepted substep ending with it
//      negative re-seats at its end. integrate()'s once-per-increment test on
//      the elastic trial (UW's rule) is undone. WP-134: with this rule
//      h < 0 is impossible (0/960 runs); UW's once-per-increment rule gives it
//      in 146/480 exact ring runs. `bracket` keeps the stage sentinel only;
//      `stale` reproduces ModifiedEuler's rule (attribution only).
//   alpha  after every accepted substep alpha is checked against the bounding
//      surface: rho_alpha = sqrt(3/2)|alpha| / alpha^b(theta_alpha, psi), the
//      Lode angle of alpha ITSELF (WP-134's geometric inside/outside test,
//      n-independent). Default: reject and cut (refuse at dT_min);
//      `-alphaProject 1`: radial projection alpha -> s alpha with the
//      deviatoric stress translated by p*(dalpha), which keeps f and n exactly,
//      counted.
//   tangent  rate-form stages (dS = C:d_eps - lambda C:R directly, no 6x6 per
//      stage, TIMs F18(b)); ONE continuum tangent at the end state for TanType
//      1 AND 2 (SAS practice). ModifiedEuler's TanType-2 chain is not
//      reproduced: it accumulates `T` where the recurrence needs `dT` (quirks).
//
// Ladruno WP-151 (R1, an OPT-IN DM04 variant; every flag OFF by default and the
// code path then byte-identical). The exact DM04 rate equations have a Zeno
// accumulation of alpha_in re-seats near the peak: after a re-seat h = inf and
// alpha slides along b; with b nearly normal to n that slide rotates n on the
// thin yield cone (radius sqrt(2/3) m), (alpha - alpha_in):n turns negative
// again, and the re-seats accumulate in finite pseudo-time while b:n -> 0+ and
// |d alpha| -> inf (Ladruno_implementation/151_sanisand_reseat_singularity.md).
// Its discrete image is `loadingWithNonPositiveDenominator` (one-to-one on the
// WP-138 wall states). Three independent options, each counted in sasStats:
//   -sasHFloor c_A     h = b0 / max((alpha - alpha_in):n, c_A sqrt(2/3) m):
//                      bounded everywhere (the 1e10 sentinel is not used);
//   -sasReseatHyst c_rev  alpha_in re-seats only when (alpha - alpha_in):n <
//                      -c_rev sqrt(2/3) m (a FINITE reversal): two re-seats then
//                      need a finite alpha travel, so they cannot accumulate;
//   -sasSoftCap kappa  where b:n < 0, h <= (1 - kappa) X / ((2/3) p |b:n|), X the
//                      elastic part of the loading denominator, so H >= kappa X.
// The floor and the hysteresis are needed TOGETHER (oracle: each alone leaves
// 97-102 of 320 wall trials singular); the cap closes deep softening (b:n << 0
// at low p), e.g. the inadmissible ring point b8 1950/3.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "UWmaterials/ManzariDafalias.h"
#include <profiler/ProfilerMacros.h>
#include <OPS_Globals.h>

#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>

namespace {
const double kSasDTmin     = 1.0e-6;   // as ModifiedEuler
const int    kSasDriftIter = 20;
// substep trace codes (the WP-127 replay trace; 0/1/4/5/7 as ModifiedEuler)
enum { TR_ACCEPT = 0, TR_REJ_ERR = 1, TR_REJ_LOWP1 = 4, TR_REJ_LOWP2 = 5, TR_CAP = 7,
       TR_REJ_NONPOS_H = 8, TR_REJ_DRIFT = 9, TR_REJ_ALPHA = 10, TR_REFUSED = 11,
       TR_ACCEPT_PROJECTED = 12, TR_REJ_REVERSAL = 13 };
// refusal codes (sasStats LAST_REFUSE_CODE)
enum { RC_START_F = 1, RC_START_ALPHA = 2, RC_START_OTHER = 3, RC_DTMIN = 4,
       RC_NONPOS_H = 5, RC_LOWP = 6, RC_DRIFT = 7, RC_ALPHA = 8, RC_CAP = 9 };
// stage kinds
enum { ST_ELASTIC = 0, ST_PLASTIC = 1, ST_NONPOS_H = -1, ST_TENSION = -2 };

bool finite6(const Vector& v)
{
    for (int i = 0; i < v.Size(); i++)
        if (!std::isfinite(v(i)))
            return false;
    return true;
}

const char* refuseName(int c)
{
    switch (c) {
    case RC_START_F:     return "startOutsideYield";
    case RC_START_ALPHA: return "startAlphaOutsideBounding";
    case RC_START_OTHER: return "startInadmissible(trace/tension/non-finite)";
    case RC_DTMIN:       return "errorAtDTmin";
    case RC_NONPOS_H:    return "loadingWithNonPositiveDenominator";
    case RC_LOWP:        return "tensionAtDTmin";
    case RC_DRIFT:       return "driftCorrectionFailed";
    case RC_ALPHA:       return "alphaOutsideBoundingAtDTmin";
    case RC_CAP:         return "maxSubsteps";
    default:             return "?";
    }
}
}  // namespace

void
ManzariDafalias::ladrunoResetSasStats(void)
{
    for (int i = 0; i < LSAS_COUNT; i++)
        mLadrunoSas.stats[i] = 0.0;
    mLadrunoSas.refused = false;
}

// Ladruno WP-152: the separation state is STATE, not census, so it is reset
// apart from ladrunoResetSasStats: by revertToStart only outside
// InitialStateAnalysis (which keeps the stress, so it must keep the state that
// produced it -- review #5), and by the replay command, which loads a NORMAL state.
void
ManzariDafalias::ladrunoResetSasSep(void)
{
    mLadrunoSas.sep = mLadrunoSas.sep_n = false;
    mLadrunoSas.sepTr = mLadrunoSas.sepTr_n = 0.0;
    mLadrunoSas.sepEvent = 0;
    mLadrunoSas.sepCode = 0;
    mLadrunoSas.sepP0 = 0.0;
    mLadrunoSas.stats[LSAS_SEP_ACTIVE] = 0.0;
}

// Ladruno WP-152 (tension cutoff): put the TRIAL on an isotropic state of model
// pressure pModel (sigma = (pModel - p_r) I) with alpha = alpha_in = 0 -- the only
// alpha consistent with an isotropic stress on DM04's thin cone -- the committed
// fabric kept (vanilla's low-p reset keeps it too), the elastic strain unchanged
// (the separation strain is not elastic), and the elastic stiffness AT that state
// as the tangent (at p_min this is the model's own moduli floor: a declared Newton
// regularisation, since a separated point's stress does not depend on the strain).
// mVoidRatio must already hold the end-of-increment value.
void
ManzariDafalias::ladrunoSasSetIsotropic(double pModel)
{
    mSigma = mI1;
    mSigma *= (pModel - m_Presidual);
    mAlpha.Zero();
    mAlpha_in.Zero();
    mFabric = mFabric_n;
    mEpsilonE = mEpsilonE_n;
    mDGamma = 0.0;
    double K, G;
    GetElasticModuli(mSigma, mVoidRatio, K, G);
    mCe = GetStiffness(K, G); mCep = mCe; mCep_Consistent = mCe;
    mLadrunoSas.stats[LSAS_LAST_RATIO_B] = 0.0;
    mLadrunoSas.stats[LSAS_LAST_F] = GetF(mSigma, mAlpha);
}

// h with the alpha_in rule applied: a stage whose (alpha - alpha_in):n is not
// positive takes the model's own sentinel for (alpha - alpha_in):n = 0
// (alpha_in re-seated at that stage).
// Ladruno WP-151: with -sasHFloor c_A > 0, h = b0 / max(x, c_A sqrt(2/3) m) for
// every x (inside a -sasReseatHyst band x < 0 too): bounded, never negative,
// and never the 1e10 sentinel. With the floor OFF the lines below it are the
// WP-129 ones, unchanged.
double
ManzariDafalias::ladrunoSasBracketH(const Vector& a, const Vector& ain, const Vector& n, double h,
    double b0)
{
    if (mLadrunoSas.opt.alphaInMode == 2)   // stale: ModifiedEuler's behaviour
        return h;
    Vector t(a);
    t -= ain;
    const double x = DoubleDot2_2_Contr(t, n);
    if (mLadrunoSas.opt.hFloor > 0.0) {                               // Ladruno WP-151
        const double eps = mLadrunoSas.opt.hFloor * root23 * m_m;
        return b0 / (x > eps ? x : eps);
    }
    if (x < small)
        return 1.0e10;
    return h;
}

// Ladruno WP-151 (-sasSoftCap kappa): where b:n < 0 (softening), cap h so that
// Kp = (2/3) p h b:n >= -(1 - kappa) X, i.e. H = Kp + X >= kappa X > 0. X is the
// elastic part of the loading denominator at the same state. h itself is capped
// (not Kp alone) so the consistency condition keeps holding: d alpha =
// (2/3) lambda h b uses the same h. OFF (kappa <= 0): h unchanged.
double
ManzariDafalias::ladrunoSasSoftCapH(double h, double bn, double p, double X)
{
    const double kappa = mLadrunoSas.opt.softCap;
    if (!(kappa > 0.0) || !(bn < 0.0) || !(X > 0.0) || !(p > 0.0))
        return h;
    const double hCap = (1.0 - kappa) * X / (two3 * p * (-bn));
    return (h > hCap) ? hCap : h;
}

// Ladruno WP-151 (-sasReseatHyst c_rev): the re-seat threshold delta = c_rev
// sqrt(2/3) m; alpha_in re-seats only where (alpha - alpha_in):n < -delta. OFF:
// 0.0, and `x < -0.0` is `x < 0.0` in IEEE arithmetic, so the WP-129 tests are
// unchanged bit for bit.
double
ManzariDafalias::ladrunoSasReseatDelta(void) const
{
    const double c = mLadrunoSas.opt.reseatHyst;
    return (c > 0.0) ? c * root23 * m_m : 0.0;
}

// rho_alpha: alpha / alpha^b with the Lode angle of alpha itself and psi at (e, p(S)).
double
ManzariDafalias::ladrunoSasAlphaRatio(const Vector& a, const Vector& s, double e)
{
    const double na = GetNorm_Contr(a);
    if (!(na > 1.0e-14))
        return (na == 0.0 || std::isfinite(na)) ? 0.0 : std::numeric_limits<double>::infinity();
    Vector u(a);
    u /= na;
    const double c3 = GetLodeAngle(u);
    double p = one3 * GetTrace(s) + m_Presidual;
    p = (p < small) ? small : p;
    const double psi = GetPSI(e, p);
    const double aB = g(c3, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
    if (!(aB > 0.0))
        return std::numeric_limits<double>::infinity();
    return sqrt(1.5) * na / aB;
}

// -alphaProject 1: alpha -> s*alpha onto the bounding surface (the Lode angle
// of alpha is unchanged by the scaling, so rho_alpha == 1 exactly), deviatoric
// stress translated by p*(alpha' - alpha) so dev(S) - p alpha, hence f and n,
// are unchanged. Elastic strain follows the stress change.
void
ManzariDafalias::ladrunoSasProject(Vector& S, Vector& A, Vector& Ee, double e)
{
    const double r = ladrunoSasAlphaRatio(A, S, e);
    if (!(r > 1.0) || !std::isfinite(r))
        return;
    double K, G;
    GetElasticModuli(S, e, K, G);
    const double p = one3 * GetTrace(S) + m_Presidual;
    Vector dA(A);
    dA *= (1.0 / r - 1.0);
    A += dA;
    Vector dS(dA);
    dS *= p;
    S += dS;
    Ee += DoubleDot4_2(GetCompliance(K, G), dS);
    mLadrunoSas.stats[LSAS_ALPHA_PROJECTED] += 1.0;
}

// Elastic stress increment for d_eps along a straight strain path, EXACT
// (review of #871, numerics item 1: one Heun step was 4-32 % off the oracle,
// independent of TolR). With the calibrated form G = g*sqrt(max(p + pRe, p_min)),
// K = c*G (c fixed by nu), dp = K dv has the closed form
//     sqrt(x) = sqrt(x0) + c*g*t*dv/2   (x = p + pRe > p_min)
//     x linear in t at the floor modulus  (x <= p_min),
// and the deviatoric part is ds = 2 G(t) de_dev dt with G linear in t on the
// first branch and constant on the second, so int_0^1 G dt is exact too.
// Only the `mUseCurrentVoidRatioInG` seam (G through the CURRENT e) breaks the
// closed form; that case falls back to 64 Heun sub-steps.
Vector
ManzariDafalias::ladrunoSasElastic(const Vector& S, const Vector& dEps, double e0, double e1)
{
    if (mUseCurrentVoidRatioInG) {
        const int nsub = 64;
        Vector Sx(S), dE(dEps), d1(6), d2(6), S1(6);
        dE /= (double)nsub;
        double K, G;
        for (int k = 0; k < nsub; k++) {
            const double ea = e0 + (e1 - e0) * (double)k / nsub;
            const double eb = e0 + (e1 - e0) * (double)(k + 1) / nsub;
            GetElasticModuli(Sx, ea, K, G);
            d1 = DoubleDot4_2(GetStiffness(K, G), dE);
            S1 = Sx; S1 += d1;
            GetElasticModuli(S1, eb, K, G);
            d2 = DoubleDot4_2(GetStiffness(K, G), dE);
            d1 += d2; d1 *= 0.5;
            Sx += d1;
        }
        Sx -= S;
        return Sx;
    }
    double K0, G0;
    GetElasticModuli(S, e0, K0, G0);
    const double x0 = one3 * GetTrace(S) + m_PreElastic;
    const double xe0 = (x0 <= m_Pmin) ? m_Pmin : x0;
    const double g = G0 / sqrt(xe0);
    const double c = K0 / G0;
    const double dv = GetTrace(dEps);
    const double spm = sqrt(m_Pmin);
    double t = 0.0, x = x0, I = 0.0;   // I = int_0^1 G(t) dt
    for (int seg = 0; seg < 4 && t < 1.0; seg++) {
        const double rest = 1.0 - t;
        if (x > m_Pmin || (x == m_Pmin && dv > 0.0)) {
            // sqrt branch
            const double u0 = sqrt(x);
            double dt = rest;
            if (dv < 0.0) {
                const double tHit = 2.0 * (spm - u0) / (c * g * dv);   // >= 0
                if (tHit < rest)
                    dt = tHit;
            }
            const double u1 = u0 + 0.5 * c * g * dt * dv;
            I += g * 0.5 * (u0 + u1) * dt;
            x = (dt < rest) ? m_Pmin : u1 * u1;
            if (dt < rest && dv < 0.0)
                x = m_Pmin - 1.0e-300;   // continue on the floor branch
            t += dt;
        } else {
            // floor branch: G = g*sqrt(p_min), K = c*G, both constant
            const double kf = c * g * spm;
            double dt = rest;
            if (dv > 0.0) {
                const double tHit = (m_Pmin - x) / (kf * dv);
                if (tHit < rest)
                    dt = tHit;
            }
            I += g * spm * dt;
            x = (dt < rest) ? m_Pmin : x + kf * dv * dt;
            t += dt;
        }
    }
    Vector dS = ToContraviant(GetDevPart(dEps));
    dS *= (2.0 * I);
    Vector dp(mI1);
    dp *= (x - x0);
    dS += dp;
    return dS;
}

// Where the EXACT elastic path from S along dEps meets the yield surface of A:
// Pegasus on [lo, hi] with f(lo) < 0 < f(hi) (the same path the predictor uses,
// so the plastic portion starts ON the surface -- review numerics item 3).
double
ManzariDafalias::ladrunoSasIntersect(const Vector& S, const Vector& A, const Vector& dEps,
    double e0, double lo, double hi)
{
    Vector dE(dEps), St(6);
    auto fAt = [&](double a) {
        dE = dEps; dE *= a;
        St = S; St += ladrunoSasElastic(S, dE, e0, e0);
        return GetF(St, A);
    };
    double f0 = fAt(lo), f1 = fAt(hi);
    if (!(f0 < 0.0 && f1 > 0.0))
        return -1.0;
    double a0 = lo, a1 = hi, a = hi;
    for (int it = 0; it < 60; it++) {
        a = a1 - f1 * (a1 - a0) / (f1 - f0);
        if (!(a > a0 && a < a1))
            a = 0.5 * (a0 + a1);
        const double f = fAt(a);
        if (fabs(f) < mTolF)
            return a;
        if (f * f1 < 0.0) {
            a0 = a1; f0 = f1;
        } else {
            f0 = f0 * f1 / (f1 + f);   // Pegasus
        }
        a1 = a; f1 = f;
        if (fabs(a1 - a0) < 1.0e-15)
            return a;
    }
    return a;
}

// One rate-form stage, everything evaluated at (s, a, z, e).
// Returns ST_ELASTIC / ST_PLASTIC / ST_NONPOS_H / ST_TENSION.
int
ManzariDafalias::ladrunoSasStage(const Vector& s, const Vector& a, const Vector& z, double e,
    const Vector& ain, double dv, const Vector& ddev,
    Vector& ds, Vector& da, Vector& dz, Vector& dep, double& lam)
{
    const double p = one3 * GetTrace(s) + m_Presidual;
    if (!(p > 0.0))
        return ST_TENSION;

    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, A, D, B, C, K, G;
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.stateDependent");
        GetStateDependent(s, a, z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, A, D, B, C, R);
        GetElasticModuli(s, e, K, G);                      // U9: the stage's own moduli
    }
    OPS_PROFILE_SCOPE("sanisand.sasME.stageArithmetic");
    {
        const double hIn = h;
        h = ladrunoSasBracketH(a, ain, n, h, b0);
        if (mLadrunoSas.opt.hFloor > 0.0) {                               // Ladruno WP-151
            Vector t(a);
            t -= ain;
            if (DoubleDot2_2_Contr(t, n) < mLadrunoSas.opt.hFloor * root23 * m_m)
                mLadrunoSas.stats[LSAS_H_FLOORED] += 1.0;
        } else if (h != hIn)   // (alpha - alpha_in):n <= -small: a reversal the stale alpha_in hides
            mLadrunoSas.stats[LSAS_H_BRACKETS] += 1.0;
    }

    // U10: the TRUE yield gradient Q = n - (1/3)(n:alpha + sqrt(2/3) m) I, so
    // Q : C : d_eps = 2G n:de_dev - K dv (n:alpha + sqrt(2/3) m)
    const double qv = DoubleDot2_2_Contr(n, a) + root23 * m_m;
    if (mLadrunoSas.opt.softCap > 0.0) {                                  // Ladruno WP-151
        const double X = 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n)))) - K * D * qv;
        const double hIn = h;
        h = ladrunoSasSoftCapH(h, DoubleDot2_2_Contr(b, n), p, X);
        if (h != hIn)
            mLadrunoSas.stats[LSAS_H_SOFTCAPPED] += 1.0;
    }
    const double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
    const double H  = Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n)))) - K * D * qv;
    const double N  = 2.0 * G * DoubleDot2_2_Mixed(n, ddev) - K * dv * qv;

    // elastic trial part, common to both kinds
    ds = ToContraviant(ddev);
    ds *= (2.0 * G);
    Vector tmp(mI1);
    tmp *= (K * dv);
    ds += tmp;

    if (!(N > 0.0)) {                       // unloading / neutral: elastic stage
        da.Zero();
        dz.Zero();
        dep.Zero();
        lam = 0.0;
        mLadrunoSas.stats[LSAS_ELASTIC_STAGES] += 1.0;
        return ST_ELASTIC;
    }
    if (!(H > small) || !std::isfinite(H))   // loading with H <= 0: no plastic solution
        return ST_NONPOS_H;

    lam = N / H;
    // C:R = 2G (B n - C (n.n - I/3)) + K D I
    Vector nn = SingleDot(n, n);
    Vector cr(n);
    cr *= B;
    Vector t2(mI1);
    t2 *= (-one3);
    t2 += nn;
    t2 *= C;
    cr -= t2;
    cr *= (2.0 * G);
    Vector t3(mI1);
    t3 *= (K * D);
    cr += t3;
    cr *= lam;
    ds -= cr;

    da = b;
    da *= (lam * two3 * h);

    dz = n;
    dz *= m_z_max;
    dz += z;
    dz *= (-1.0 * lam * m_cz * Macauley(-1.0 * D));

    dep = ToCovariant(R);
    dep *= lam;
    return ST_PLASTIC;
}

// Consistent (Potts-Gens) then normal drift correction, moduli at the current
// state. true = |f| <= TolF (or f <= TolF when !bothSides). false = neither
// direction reduces |f|, the iteration limit, or tension -- the caller turns
// that into a cut / refusal, never into an accepted f > TolF.
bool
ManzariDafalias::ladrunoSasDrift(Vector& S, Vector& A, Vector& Z, Vector& Ee, double e,
    const Vector& ain, bool bothSides)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.drift");
    bool corrected = false;
    for (int it = 0; it <= kSasDriftIter; it++) {
        const double f = GetF(S, A);
        if (!std::isfinite(f))
            return false;
        if ((bothSides ? fabs(f) : f) <= mTolF) {
            if (corrected)
                mLadrunoSas.stats[LSAS_DRIFT_CORRECTIONS] += 1.0;
            return true;
        }
        if (it == kSasDriftIter)
            return false;
        const double p = one3 * GetTrace(S) + m_Presidual;
        if (!(p > 0.0))
            return false;

        Vector n(6), d(6), b(6), R(6);
        double cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C, K, G;
        GetStateDependent(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C, R);
        GetElasticModuli(S, e, K, G);
        h = ladrunoSasBracketH(A, ain, n, h, b0);
        const Matrix aC = GetStiffness(K, G);

        // exact df/dsigma = n - (1/3)(n:alpha + sqrt(2/3) m) I ; df/dalpha = -p n
        Vector Q(mI1);
        Q *= (-one3 * (DoubleDot2_2_Contr(n, A) + root23 * m_m));
        Q += n;
        Vector Rc = ToCovariant(R);
        Vector dSp = DoubleDot4_2(aC, Rc);
        if (mLadrunoSas.opt.softCap > 0.0)                                // Ladruno WP-151: X = Q:C:R here
            h = ladrunoSasSoftCapH(h, DoubleDot2_2_Contr(b, n), p, DoubleDot2_2_Contr(Q, dSp));
        Vector aBar(b);
        aBar *= (two3 * h);
        Vector zBar(n);
        zBar *= m_z_max;
        zBar += Z;
        zBar *= (-1.0 * m_cz * Macauley(-1.0 * D));

        bool moved = false;
        const double den = DoubleDot2_2_Contr(Q, dSp) + p * DoubleDot2_2_Contr(n, aBar);
        if (den > small && std::isfinite(den)) {
            const double lam = f / den;
            Vector S1(S), A1(A);
            Vector t(dSp); t *= lam; S1 -= t;
            Vector ta(aBar); ta *= lam; A1 += ta;
            const double f1 = GetF(S1, A1);
            if (std::isfinite(f1) && fabs(f1) < fabs(f) && one3 * GetTrace(S1) + m_Presidual > 0.0) {
                S = S1;
                A = A1;
                Vector tz(zBar); tz *= lam; Z += tz;
                Vector te(Rc); te *= lam; Ee -= te;
                moved = true;
            }
        }
        if (!moved) {
            const double qq = DoubleDot2_2_Contr(Q, Q);
            if (!(qq > 0.0))
                return false;
            const double lam = f / qq;
            Vector S1(S);
            Vector t(Q); t *= lam; S1 -= t;
            const double f1 = GetF(S1, A);
            if (!(std::isfinite(f1) && fabs(f1) < fabs(f) && one3 * GetTrace(S1) + m_Presidual > 0.0))
                return false;    // the give-up: FAIL, never hand back f > TolF
            S = S1;
            Ee -= DoubleDot4_2(GetCompliance(K, G), t);
        }
        corrected = true;
    }
    return false;
}

// The continuum elastoplastic tangent at (S, A, Z), plastic loading assumed,
// moduli at that state.
void
ManzariDafalias::ladrunoSasContinuumTangent(const Vector& S, const Vector& A, const Vector& Z,
    const Vector& ain, double e, Matrix& Cep)
{
    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, Af, D, B, C, K, G;
    GetStateDependent(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, Af, D, B, C, R);
    GetElasticModuli(S, e, K, G);
    h = ladrunoSasBracketH(A, ain, n, h, b0);
    if (mLadrunoSas.opt.softCap > 0.0) {                                  // Ladruno WP-151: the stage's X
        const double p = one3 * GetTrace(S) + m_Presidual;
        const double qv = DoubleDot2_2_Contr(n, A) + root23 * m_m;
        const double X = 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n)))) - K * D * qv;
        h = ladrunoSasSoftCapH(h, DoubleDot2_2_Contr(b, n), p, X);
    }
    Vector dummy(6);
    Cep = GetElastoPlasticTangent(S, 1.0, dummy, dummy, G, K, B, C, D, h, n, d, b);
}

// The substep loop over the plastic portion curStrain -> nextStrain. Returns 0
// or a refusal code (RC_*). On 0, (S, Ee, A, Z, ain) hold the end state.
// `onset`: the portion starts at a plastic onset (paper alpha_in rule).
int
ManzariDafalias::ladrunoSasSubsteps(Vector& S, Vector& Ee, Vector& A, Vector& Z, Vector& ain,
    const Vector& curStrain, const Vector& nextStrain, bool onset,
    double& lamSum, bool& lastPlastic)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.substeps");
    double* st = mLadrunoSas.stats;
    const LadrunoSasOptions& o = mLadrunoSas.opt;
    const bool paperRule = (o.alphaInMode == 0);
    const double TolE   = mTolR;
    const double sigRef = (o.errFloor < 0.0) ? m_P_atm / 101.0 : o.errFloor;
    const double kappa  = o.alphaBoundTol;
    const double dRev   = ladrunoSasReseatDelta();   // Ladruno WP-151: 0.0 = the paper rule

    Vector dStrain(nextStrain);
    dStrain -= curStrain;
    const double trD = GetTrace(dStrain);
    const Vector devD = GetDevPart(dStrain);

    Vector ds1(6), da1(6), dz1(6), dep1(6), ds2(6), da2(6), dz2(6), dep2(6);
    Vector S1(6), A1(6), Z1(6), nS(6), nA(6), nZ(6), nEe(6), tmp(6), ddev(6);
    Vector ain1(6), ain2(6), tmpN(6);
    double T = 0.0, dT = 1.0;
    bool lastRejected = false;
    (void)onset;   // the paper rule below needs no onset flag: DM04 re-seats at an
                   // onset only when (alpha - alpha_in):n < 0 there (WP-134 oracle)
    lamSum = 0.0;
    lastPlastic = false;

    while (T < 1.0) {
        if (mSubstepsTakenInME < INT_MAX)
            ++mSubstepsTakenInME;
        st[LSAS_SUBSTEPS] += 1.0;
        st[LSAS_LAST_SUBSTEPS] += 1.0;
        if (st[LSAS_LAST_SUBSTEPS] > st[LSAS_MAX_ONE_UPDATE])
            st[LSAS_MAX_ONE_UPDATE] = st[LSAS_LAST_SUBSTEPS];
        const bool atMin = (dT <= kSasDTmin);
        if (mMaxSubstepsInME > 0 && mSubstepsTakenInME > mMaxSubstepsInME) {
            mSubstepCapHitInME = true;
            ladrunoTraceSubstep(T, dT, std::numeric_limits<double>::quiet_NaN(), TR_CAP, atMin);
            return RC_CAP;
        }

        tmp = dStrain; tmp *= T; tmp += curStrain;
        const double e0 = m_e_init - (1 + m_e_init) * GetTrace(tmp);
        tmp = dStrain; tmp *= (T + dT); tmp += curStrain;
        const double e1 = m_e_init - (1 + m_e_init) * GetTrace(tmp);
        const double dv = dT * trD;
        ddev = devD; ddev *= dT;

        double lam1 = 0.0, lam2 = 0.0;
        int k1, k2 = ST_ELASTIC;
        // paper rule (DM04, = WP-134's decide()): at a stage where
        // (alpha - alpha_in):n < 0 a new loading process starts there:
        // alpha_in := alpha at THAT stage's state.
        // Only ON the surface (the oracle decides modes, and re-seats, only
        // there -- review numerics item 7).
        bool reseat1 = false, reseat2 = false;
        const bool onSurface = (fabs(GetF(S, A)) <= mTolF);
        ain1 = ain;
        if (paperRule && onSurface) {
            tmpN = GetNormalToYield(S, A);
            tmp = A; tmp -= ain;
            const double x1 = DoubleDot2_2_Contr(tmp, tmpN);
            if (x1 < -dRev) { ain1 = A; reseat1 = true; }
            else if (x1 < 0.0) st[LSAS_RESEAT_HELD] += 1.0;   // Ladruno WP-151: a sub-threshold reversal
        }
        {
            OPS_PROFILE_SCOPE("sanisand.sasME.stages");
            k1 = ladrunoSasStage(S, A, Z, e0, ain1, dv, ddev, ds1, da1, dz1, dep1, lam1);
        }
        if (k1 == ST_NONPOS_H) {
            // a property of the accepted start of this substep: no cut can change it
            ladrunoTraceSubstep(T, dT, std::numeric_limits<double>::quiet_NaN(), TR_REFUSED, atMin);
            return RC_NONPOS_H;
        }
        int rej = 0;   // 0 none, > 0 trace code of the rejection, -1 accepted after projection
        if (k1 == ST_TENSION) {
            rej = TR_REJ_LOWP1;
        } else {
            S1 = S; S1 += ds1;
            A1 = A; A1 += da1;
            Z1 = Z; Z1 += dz1;
            if (!(one3 * GetTrace(S1) + m_Presidual > 0.0)) {
                rej = TR_REJ_LOWP1;
            } else {
                ain2 = ain1;
                if (paperRule && k1 == ST_PLASTIC) {
                    // a reversal INSIDE a plastic substep: locate it at a substep
                    // start (cut) rather than re-seating on the Euler predictor;
                    // at dT_min re-seat at the stage state.
                    tmpN = GetNormalToYield(S1, A1);
                    tmp = A1; tmp -= ain1;
                    if (DoubleDot2_2_Contr(tmp, tmpN) < -dRev) {   // Ladruno WP-151: -dRev (0 = paper)
                        if (atMin) { ain2 = A1; reseat2 = true; }
                        else rej = TR_REJ_REVERSAL;
                    }
                }
                if (rej == 0) {
                    OPS_PROFILE_SCOPE("sanisand.sasME.stages");
                    k2 = ladrunoSasStage(S1, A1, Z1, e1, ain2, dv, ddev, ds2, da2, dz2, dep2, lam2);
                }
                if (k2 == ST_TENSION)
                    rej = TR_REJ_LOWP2;
                else if (k2 == ST_NONPOS_H)
                    rej = TR_REJ_NONPOS_H;
            }
        }
        double err = std::numeric_limits<double>::quiet_NaN();
        if (rej == 0) {
            nS = ds1; nS += ds2; nS *= 0.5; nS += S;
            nA = da1; nA += da2; nA *= 0.5; nA += A;
            nZ = dz1; nZ += dz2; nZ *= 0.5; nZ += Z;
            if (!(one3 * GetTrace(nS) + m_Presidual > 0.0))
                rej = TR_REJ_LOWP2;
            else {
                tmp = ds2; tmp -= ds1;
                err = GetNorm_Contr(tmp) / fmax(2.0 * GetNorm_Contr(S), sigRef);
                if (o.errorVars == 0) {
                    tmp = da2; tmp -= da1;
                    err = fmax(err, GetNorm_Contr(tmp) / fmax(2.0 * GetNorm_Contr(A), 1.0));
                    tmp = dz2; tmp -= dz1;
                    err = fmax(err, GetNorm_Contr(tmp) / fmax(2.0 * GetNorm_Contr(Z), 1.0));
                }
                if (!std::isfinite(err) || !finite6(nS) || !finite6(nA) || !finite6(nZ))
                    err = std::numeric_limits<double>::infinity();
                err = fmax(err, DBL_EPSILON);
                if (err > TolE)
                    rej = TR_REJ_ERR;
            }
        }

        // the alpha_in this substep ends with (paper rule: the re-seat point)
        Vector ainEnd(reseat2 ? ain2 : (reseat1 ? ain1 : ain));
        if (rej == 0) {
            // candidate accepted by the error test: drift, then alpha
            nEe = dStrain; nEe *= dT; nEe += Ee;
            tmp = dep1; tmp += dep2; tmp *= 0.5; nEe -= tmp;
            const bool anyElastic = (k1 == ST_ELASTIC || k2 == ST_ELASTIC);
            const bool anyPlastic = (k1 == ST_PLASTIC || k2 == ST_PLASTIC);
            // the plastic portion belongs ON the surface; an unloading stage may leave it inside
            // ... and only when the substep STARTED on it: a start inside (f < -TolF)
            // is corrected only if it ends outside (review numerics item 3)
            const bool bothSides = anyPlastic && !anyElastic && onSurface;
            if (!ladrunoSasDrift(nS, nA, nZ, nEe, e1, ainEnd, bothSides)) {
                rej = TR_REJ_DRIFT;
            } else {
                OPS_PROFILE_SCOPE("sanisand.sasME.alphaCheck");
                // Rejected only when PLASTIC FLOW carried alpha outward beyond
                // 1 + kappa -- measured against the END surface with the start
                // alpha, so psi moving the surface is never the reason (review
                // numerics item 2: that was a dead end).
                const double ratio = ladrunoSasAlphaRatio(nA, nS, e1);
                const double ratioFixed = ladrunoSasAlphaRatio(A, nS, e1);
                if (!(ratio <= 1.0 + kappa) && !(ratio <= ratioFixed)) {
                    if (o.alphaProject != 0 && std::isfinite(ratio)) {
                        ladrunoSasProject(nS, nA, nEe, e1);
                        rej = -1;   // accepted, projected
                    } else {
                        rej = TR_REJ_ALPHA;
                    }
                }
                if (rej <= 0) {
                    const double r2 = ladrunoSasAlphaRatio(nA, nS, e1);
                    if (r2 > st[LSAS_MAX_RATIO_B])
                        st[LSAS_MAX_RATIO_B] = r2;
                }
            }
        }

        if (rej > 0) {
            // rejected: cut, or refuse at dT_min
            double q = 0.5;
            int rc = RC_DTMIN;
            switch (rej) {
            case TR_REJ_ERR:      q = fmax(0.9 * sqrt(TolE / err), 0.1); rc = RC_DTMIN;
                                  st[LSAS_REJ_ERR] += 1.0; break;
            case TR_REJ_LOWP1:
            case TR_REJ_LOWP2:    q = 0.1; rc = RC_LOWP; st[LSAS_REJ_LOWP] += 1.0; break;
            case TR_REJ_NONPOS_H: q = 0.5; rc = RC_NONPOS_H; st[LSAS_REJ_NONPOS_H] += 1.0; break;
            case TR_REJ_DRIFT:    q = 0.5; rc = RC_DRIFT; st[LSAS_REJ_DRIFT] += 1.0; break;
            case TR_REJ_ALPHA:    q = 0.5; rc = RC_ALPHA; st[LSAS_REJ_ALPHA] += 1.0; break;
            case TR_REJ_REVERSAL: q = 0.5; rc = RC_DTMIN; st[LSAS_REJ_REVERSAL] += 1.0; break;
            }
            if (!std::isfinite(q)) q = 0.1;
            ladrunoTraceSubstep(T, dT, err, atMin ? TR_REFUSED : rej, atMin);
            if (atMin)
                return rc;
            dT = fmax(q * dT, kSasDTmin);
            lastRejected = true;
            continue;
        }

        // accept
        st[LSAS_ACCEPTED] += 1.0;
        ladrunoTraceSubstep(T, dT, err, rej < 0 ? TR_ACCEPT_PROJECTED : TR_ACCEPT, atMin);
        S = nS; A = nA; Z = nZ; Ee = nEe;
        lamSum += 0.5 * (lam1 + lam2);
        lastPlastic = (k2 == ST_PLASTIC);   // the loading index nearest the end state
        if (reseat1 || reseat2) {
            ain = ainEnd;
            st[LSAS_ALPHA_IN_RESEATS] += 1.0;
        }
        // (alpha - alpha_in):n reaching 0 re-seats alpha_in (the paper's rule),
        // on the surface only
        if (paperRule && fabs(GetF(S, A)) <= mTolF) {
            Vector nEnd = GetNormalToYield(S, A);
            tmp = A; tmp -= ain;
            const double xEnd = DoubleDot2_2_Contr(tmp, nEnd);
            if (xEnd < -dRev) {                                          // Ladruno WP-151: -dRev (0 = paper)
                ain = A;
                st[LSAS_ALPHA_IN_RESEATS] += 1.0;
            } else if (xEnd < 0.0)
                st[LSAS_RESEAT_HELD] += 1.0;                             // Ladruno WP-151
        }
        T += dT;
        double q = fmin(fmax(0.9 * sqrt(TolE / err), 0.1), 1.1);
        if (lastRejected)
            q = fmin(q, 1.0);
        lastRejected = false;
        dT = fmax(q * dT, kSasDTmin);
        dT = fmin(dT, 1.0 - T);
    }
    return 0;
}

// The SAS-ME update: mSigma_n.. -> mSigma.. for the strain mEpsilon.
void
ManzariDafalias::ladrunoSasIntegrate(void)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.update");
    double* st = mLadrunoSas.stats;
    const LadrunoSasOptions& o = mLadrunoSas.opt;
    st[LSAS_UPDATES] += 1.0;
    st[LSAS_LAST_SUBSTEPS] = 0.0;
    st[LSAS_LAST_REFUSE_CODE] = 0.0;
    mLadrunoSas.refused = false;

    const Vector& CurStrain = mEpsilon_n;
    const Vector& NextStrain = mEpsilon;

    Vector dStrain(NextStrain);
    dStrain -= CurStrain;
    const double eN = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
    mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);

    // ---- Ladruno WP-152: the tension cutoff (separation) ----------------------
    // The trial starts as the committed state. A SEPARATED point carries no
    // tension and no shear (model p = p_min) and absorbs the strain; it
    // re-contacts once the volumetric opening since entry has closed with an
    // overlap that gives p_contact: g = tr(eps) - tr(eps_entry) (compression
    // positive) >= g_c = (p_contact - p_min) / K(p_contact). Then
    // p_re = p_min + K(p_contact) g >= p_contact, alpha = alpha_in = 0, and SAS-ME
    // resumes at the next update. Ladruno_implementation/152_sanisand_tension_cutoff.md.
    mLadrunoSas.sep = mLadrunoSas.sep_n;
    mLadrunoSas.sepTr = mLadrunoSas.sepTr_n;
    mLadrunoSas.sepEvent = 0;   // an element may update several times per step: the census
                                // counts the COMMITTED transition (LadrunoSANISAND::commitState)
    mLadrunoSas.sepCode = 0;
    mLadrunoSas.sepP0 = 0.0;
    const bool tcOn = (o.tcPcontact > 0.0);
    if (tcOn && mLadrunoSas.sep_n) {
        const double g = GetTrace(NextStrain) - mLadrunoSas.sepTr_n;
        Vector Sc(mI1);
        Sc *= (o.tcPcontact - m_Presidual);
        double Kc, Gc;
        GetElasticModuli(Sc, mVoidRatio, Kc, Gc);
        const double gc = (o.tcPcontact - m_Pmin) / Kc;
        if (g >= gc) {
            ladrunoSasSetIsotropic(m_Pmin + Kc * g);
            mLadrunoSas.sep = false;
            mLadrunoSas.sepEvent = 3;
            mLadrunoLastPath = 9;
        } else {
            ladrunoSasSetIsotropic(m_Pmin);
            mLadrunoLastPath = 8;
        }
        return;
    }

    // Paper alpha_in rule (the default): alpha_in changes ONLY at a plastic
    // onset or where (alpha - alpha_in):n reaches 0 -- both decided inside this
    // update. integrate()'s once-per-increment test on the elastic trial
    // ((alpha_n - alpha_in_n):Ce:d_eps < 0, UW's rule U6) is therefore UNDONE
    // here: start from the committed alpha_in. (WP-134: with UW's rule 146/480
    // exact ring runs reach h < 0; with the paper's, 0/960.) bracket / stale
    // keep UW's pre-decision.
    const bool paperRule = (o.alphaInMode == 0);
    Vector S(mSigma_n), A(mAlpha_n), Z(mFabric_n), Ee(mEpsilonE_n),
           ain(paperRule ? mAlpha_in_n : mAlpha_in);
    int code = 0;
    const double p0c = one3 * GetTrace(S) + m_Presidual;   // Ladruno WP-152: committed p (the E2 test)
    bool startTension = false;                             // Ladruno WP-152: code 3 BECAUSE p0 <= 0 (E1)

    // ---- 0. entry: the committed state must be admissible ------------------
    {
        const double p0 = one3 * GetTrace(S) + m_Presidual;
        const double ta = GetTrace(A), tz = GetTrace(Z);
        if (!finite6(S) || !finite6(A) || !finite6(Z) || !finite6(ain) || !(p0 > 0.0)
            || fabs(ta) > 1.0e-6 * fmax(GetNorm_Contr(A), m_m)
            || fabs(tz) > 1.0e-6 * fmax(GetNorm_Contr(Z), m_m)) {
            code = RC_START_OTHER;
            // Ladruno WP-152: only "finite, traces clean, p0 <= 0" is a tension
            // start; a non-finite value or a trace violation stays a loud refusal.
            startTension = finite6(S) && finite6(A) && finite6(Z) && finite6(ain) && !(p0 > 0.0)
                && !(fabs(ta) > 1.0e-6 * fmax(GetNorm_Contr(A), m_m))
                && !(fabs(tz) > 1.0e-6 * fmax(GetNorm_Contr(Z), m_m));
        }
        else if (GetF(S, A) > mTolF)
            code = RC_START_F;
        else {
            // ENTRY threshold 1 + kappa_entry (review numerics item 2): the
            // continuum itself carries alpha past the bounding surface when psi
            // shrinks it around a fixed alpha (elastic compression), so a start
            // between 1 + kappa and 1 + kappa_entry is ADMISSIBLE (counted, not
            // refused -- refusing it was a dead end no cut could leave). Beyond
            // it is a state no path of the model reaches (b8 1950/2-3: 6.8/7.3).
            const double r0 = ladrunoSasAlphaRatio(A, S, eN);
            if (!(r0 <= 1.0 + o.alphaEntryTol)) {
                if (o.alphaProject != 0 && std::isfinite(r0))
                    ladrunoSasProject(S, A, Ee, eN);
                else
                    code = RC_START_ALPHA;
            } else if (r0 > 1.0 + o.alphaBoundTol) {
                st[LSAS_ENTRY_OVER_KAPPA] += 1.0;
            }
        }
    }

    double K, G;
    double lamSum = 0.0;
    bool lastPlastic = false;
    if (code == 0) {
        // ---- 1. elastic predictor / intersection ----------------------------
        double a = 0.0;
        bool elastic = false, onset = false;
        {
            OPS_PROFILE_SCOPE("sanisand.sasME.predictor");
            Vector dSe = ladrunoSasElastic(S, dStrain, eN, mVoidRatio);
            Vector St(S);
            St += dSe;
            const double ft = GetF(St, A);
            const double pt = one3 * GetTrace(St) + m_Presidual;
            if (!(pt < m_Presidual) && ft <= mTolF) {
                elastic = true;
                S = St;
                Ee += dStrain;
                mLadrunoLastPath = 0;
            } else {
                const double f0 = GetF(S, A);
                Vector nY = GetNormalToYield(S, A);
                // Paper alpha_in rule at the increment START (review numerics
                // item 4) with WP-134's convention exactly: a start is ON the
                // surface unless f0 < -1e-8 of the cone radius sqrt(2/3) m p, and
                // on the surface (alpha - alpha_in):n < 0 re-seats alpha_in before
                // the mode is chosen. A start further inside is interior (the
                // onset, if any, is at the far side, decided there). Measured: this
                // is what separates b16 5496 (re-seat) from b16 5471 / b8 1961
                // (no re-seat) on the ring, all within 4e-8 of the cone radius.
                if (paperRule) {
                    const double pS = one3 * GetTrace(S) + m_Presidual;
                    if (f0 >= -1.0e-8 * root23 * m_m * pS) {
                        Vector t0(A); t0 -= ain;
                        const double x0 = DoubleDot2_2_Contr(t0, nY);
                        if (x0 < -ladrunoSasReseatDelta()) {                // Ladruno WP-151: 0 = paper
                            ain = A;
                            st[LSAS_ALPHA_IN_RESEATS] += 1.0;
                        } else if (x0 < 0.0)
                            st[LSAS_RESEAT_HELD] += 1.0;                     // Ladruno WP-151
                    }
                }
                // U10: the loading test on the TRUE gradient
                Vector Q(mI1);
                Q *= (-one3 * (DoubleDot2_2_Contr(nY, A) + root23 * m_m));
                Q += nY;
                const double qn = GetNorm_Contr(Q), nd = GetNorm_Contr(dSe);
                if (f0 < -mTolF) {
                    a = ladrunoSasIntersect(S, A, dStrain, eN, 0.0, 1.0);
                    if (a < 0.0) { st[LSAS_INTERSECT_FAIL] += 1.0; a = 0.0; }
                    mLadrunoLastPath = 2;
                    onset = true;
                } else if (DoubleDot2_2_Contr(Q, dSe) / ((qn == 0 ? 1.0 : qn) * (nd == 0 ? 1.0 : nd))
                           > (-sqrt(mTolF))) {
                    a = 0.0;
                    mLadrunoLastPath = 3;
                } else {
                    // unload then reload: find the first sampled point inside,
                    // then the first one back outside, and Pegasus between them
                    a = -1.0;
                    const int ns = 64;
                    Vector dE(6), St2(6);
                    double prevA = -1.0;
                    bool inside = false;
                    for (int k = 1; k <= ns; k++) {
                        const double ak = (double)k / ns;
                        dE = dStrain; dE *= ak;
                        St2 = S; St2 += ladrunoSasElastic(S, dE, eN, eN);
                        const double fk = GetF(St2, A);
                        if (!inside) {
                            if (fk < -mTolF) inside = true;
                        } else if (fk > 0.0) {
                            a = ladrunoSasIntersect(S, A, dStrain, eN, prevA, ak);
                            break;
                        }
                        prevA = ak;
                    }
                    if (a < 0.0) { st[LSAS_INTERSECT_FAIL] += 1.0; a = 0.0; }
                    mLadrunoLastPath = 4;
                    onset = true;
                }
                mLadrunoLastElasticRatio = a;
                if (a > 0.0) {
                    Vector dEa(dStrain);
                    dEa *= a;
                    const double ea = eN - (1 + m_e_init) * GetTrace(dEa);
                    S += ladrunoSasElastic(S, dEa, eN, ea);
                    Ee += dEa;
                }
            }
        }
        if (elastic) {
            st[LSAS_ELASTIC] += 1.0;
            mSigma = S; mAlpha = A; mFabric = Z; mEpsilonE = Ee;
            mAlpha_in = ain;
            mDGamma = 0.0;
            GetElasticModuli(S, mVoidRatio, K, G);
            mCe = GetStiffness(K, G); mCep = mCe; mCep_Consistent = mCe;
            st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
            st[LSAS_LAST_F] = GetF(S, A);
            return;
        }
        // ---- 2. the plastic portion ------------------------------------------
        Vector curP(dStrain);
        curP *= a;
        curP += CurStrain;
        code = ladrunoSasSubsteps(S, Ee, A, Z, ain, curP, NextStrain, onset, lamSum, lastPlastic);
    } else {
        mLadrunoLastPath = 7;
        ladrunoTraceSubstep(0.0, 1.0, std::numeric_limits<double>::quiet_NaN(), TR_REFUSED, false);
    }

    // ---- Ladruno WP-152: tension cutoff ENTRY, masking ONLY low-p / tension -----
    // E1: tension -- code 6 (a stage or the predictor at p <= 0 at dT_min), or
    //     code 3 because the committed p0 <= 0 (finite, traces clean) -- AND the
    //     committed p0 <= p0max (review #2: one Newton iterate can carry a 10-20 kPa
    //     point through p = 0; that is a step to cut, not a separation).
    // E2: a low-confinement accuracy or cost failure -- code 4 or 9 while the
    //     committed p0 < p_sep (p_sep = 0 disables E2: a pure tension cutoff) --
    //     AND a volumetric increment that does not compress (tr d_eps <= 0,
    //     compression positive; review #1: the B/8 top row sits at p' 0.2-0.35 in
    //     situ, below p_sep, so the committed p0 alone would let an accuracy
    //     failure under COMPRESSION separate). A non-compressing increment also
    //     keeps the elastic predictor's p at or below p0 < p_sep.
    // A qualifying refusal the bound or the gate holds back is counted
    // (sepHeldHighP, sepHeldCompressing) and still REFUSES with its own code.
    // Everything else (code 5 loadingNonPosH = the alpha_in singularity, code 2,
    // code 3 non-finite / trace, codes 7, 8, and code 4/9 at p0 >= p_sep) still
    // refuses below, at any p.
    if (tcOn && code != 0) {
        const bool e1q = (code == RC_LOWP) || (code == RC_START_OTHER && startTension);
        const bool e2q = (code == RC_DTMIN || code == RC_CAP) && (p0c < o.tcPsep);
        const bool e1 = e1q && (p0c <= o.tcP0Max);
        const bool e2 = e2q && !(GetTrace(dStrain) > 0.0);
        if (e1q && !e1) st[LSAS_SEP_HELD_HIGHP] += 1.0;
        if (e2q && !e2) st[LSAS_SEP_HELD_COMPRESSING] += 1.0;
        if (e1 || e2) {
            mLadrunoSas.sepEvent = e1 ? 1 : 2;
            mLadrunoSas.sepCode = code;   // the masked refusal (review #3), recorded at commit
            mLadrunoSas.sepP0 = p0c;
            mSubstepCapHitInME = false;   // a code 9 set it; the cap was reached by a separating point
            mLadrunoSas.sep = true;
            mLadrunoSas.sepTr = GetTrace(NextStrain);
            ladrunoSasSetIsotropic(m_Pmin);
            mLadrunoLastPath = 8;
            return;
        }
    }

    if (code != 0) {
        // REFUSED: the trial is put back on the committed state -- a refusing
        // material must not also be an inventing one -- and the wrapper returns
        // LADRUNO_MATERIAL_REFUSED (LadrunoSANISAND::ladrunoUpdateStatus).
        mLadrunoSas.refused = true;
        st[LSAS_REFUSALS] += 1.0;
        st[LSAS_LAST_REFUSE_CODE] = (double)code;
        switch (code) {
        case RC_START_F:     st[LSAS_REF_START_F] += 1.0; break;
        case RC_START_ALPHA: st[LSAS_REF_START_ALPHA] += 1.0; break;
        case RC_START_OTHER: st[LSAS_REF_START_OTHER] += 1.0; break;
        case RC_DTMIN:       st[LSAS_REF_DTMIN] += 1.0; break;
        case RC_NONPOS_H:    st[LSAS_REF_NONPOS_H] += 1.0; break;
        case RC_LOWP:        st[LSAS_REF_LOWP] += 1.0; break;
        case RC_DRIFT:       st[LSAS_REF_DRIFT] += 1.0; break;
        case RC_ALPHA:       st[LSAS_REF_ALPHA] += 1.0; break;
        case RC_CAP:         st[LSAS_REF_CAP] += 1.0; break;
        }
        mSigma = mSigma_n; mAlpha = mAlpha_n; mFabric = mFabric_n; mEpsilonE = mEpsilonE_n;
        if (paperRule)
            mAlpha_in = mAlpha_in_n;
        mDGamma = 0.0;
        mVoidRatio = eN;          // the refused trial carries nothing forward
        GetElasticModuli(mSigma_n, eN, K, G);
        mCe = GetStiffness(K, G); mCep = mCe; mCep_Consistent = mCe;
        // "last" columns describe no valid end state after a refusal
        st[LSAS_LAST_RATIO_B] = std::numeric_limits<double>::quiet_NaN();
        st[LSAS_LAST_F] = std::numeric_limits<double>::quiet_NaN();
        // PER-INSTANCE warn-once (review #871 item 2: no process-wide state,
        // nothing to reset on wipe); every refusal is counted in sasStats.
        if (!mLadrunoSas.warned) {
            mLadrunoSas.warned = true;
            opserr << "WARNING LadrunoSANISAND (SAS-ME, IntScheme 129) material tag " << this->getTag()
                   << ": update REFUSED (" << refuseName(code) << ", code " << code
                   << "); the trial is left on the committed state and the step must be cut."
                   << " Once per integration point; every refusal is counted in `sasStats`."
                   << endln;
        }
        return;
    }

    // ---- 3. accepted: state and ONE end-state tangent -----------------------
    mSigma = S; mAlpha = A; mFabric = Z; mEpsilonE = Ee;
    mAlpha_in = ain;
    mDGamma = lamSum;
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.tangent");
        GetElasticModuli(S, mVoidRatio, K, G);
        mCe = GetStiffness(K, G);
        if (lastPlastic)
            ladrunoSasContinuumTangent(S, A, Z, ain, mVoidRatio, mCep);
        else
            mCep = mCe;
        mCep_Consistent = mCep;
    }
    st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
    st[LSAS_LAST_F] = GetF(S, A);
}
