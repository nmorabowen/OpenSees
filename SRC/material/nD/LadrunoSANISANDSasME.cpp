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

// Ladruno WP-129 (TIMs F18(a), F20(b)/(c); amended by the WP-128 verdict):
// SAS-ME, IntScheme 129 -- a Sloan-Abbo-Sheng-style explicit modified Euler
// for ManzariDafalias / LadrunoSANISAND.
//
// These are MEMBER functions of the vanilla class ManzariDafalias (declared in
// UWmaterials/ManzariDafalias.h behind `// Ladruno WP-129`), defined here so the
// vanilla footprint is a state block, the declarations and one dispatch branch
// in integrate(). The scheme is reachable only when mLadrunoSas.allowed is
// true, i.e. from a LadrunoSANISAND instance, whose wrappers forward the
// refusal to the element (LADRUNO_MATERIAL_REFUSED) and so to analyze().
//
// WHAT IT CLOSES, relative to ModifiedEuler (IntScheme 1), item by item:
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
//   F  each stage is classified from the elastic trial -- the loading-index
//      numerator N = n_f : C : d_eps -- never from the sign of the full plastic
//      multiplier. N <= 0: elastic stage, alpha and z UNCHANGED (never
//      re-derived from the stress ratio). N > 0 with a non-positive
//      denominator H = Kp + n_f:C:R: at stage 1 (a property of the accepted
//      start state, which no cut can change) a REFUSAL; at stage 2 a cut, a
//      refusal at dT_min.
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
//   G  alpha_in is kept consistent INSIDE the increment (WP-128: the trigger of
//      the alpha escape). Default `reseat`: a stage with (alpha - alpha_in):n
//      <= 0 uses the model's own h = 1e10 sentinel (the value GetStateDependent
//      gives for (alpha - alpha_in):n = 0, i.e. alpha_in re-seated at that
//      stage), and after every ACCEPTED substep alpha_in := alpha where
//      (alpha - alpha_in):n < 0 -- Dafalias-Manzari's alpha_in is the
//      back-stress at the last reversal, and (alpha - alpha_in):n < 0 is a
//      reversal the once-per-increment test in integrate() cannot see.
//      `bracket` keeps the stage rule only (h = b0/<(alpha-alpha_in):n>);
//      `stale` reproduces ModifiedEuler's defect (attribution only).
//   alpha  after every accepted substep alpha is checked against the bounding
//      surface: ratio = sqrt(3/2)|alpha| / alpha^b(theta_alpha, psi) with the
//      Lode angle of alpha ITSELF (an alpha-space radius, n-independent, so an
//      elastic rotation of n cannot make a state inadmissible). Default:
//      reject and cut (refuse at dT_min); `-alphaProject 1`: radial projection
//      alpha -> s alpha with the deviatoric stress translated by p*(dalpha),
//      which keeps f and n exactly, counted.
//   tangent  rate-form stages (dS = C:d_eps - lambda C:R directly, no 6x6 per
//      stage, TIMs F18(b)); ONE continuum tangent at the end state for TanType
//      1 AND 2 (SAS practice). ModifiedEuler's TanType-2 chain is not
//      reproduced: it accumulates `T` where the recurrence needs `dT`, so with
//      N equal substeps and elastic stages it returns (N+1)/2 x Ce (quirks).
//
// MODULI: K, G are those of the committed state (as every other explicit
// scheme and the elastic predictor), frozen over the increment, and recomputed
// here rather than read from mK/mG (which can be stale-by-stage after the
// elastic->plastic flip). e is taken at the substep start for both stages.
// Hence SAS-ME and ModifiedEuler integrate the SAME rate equations; where none
// of the fixes binds, they agree to integration tolerance.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "UWmaterials/ManzariDafalias.h"
#include <profiler/ProfilerMacros.h>
#include <OPS_Globals.h>

#include <atomic>
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>

namespace {
const double kSasDTmin     = 1.0e-6;   // as ModifiedEuler
const int    kSasDriftIter = 20;
// substep trace codes (the WP-127 replay trace; 0/1/4/5 as ModifiedEuler)
enum { TR_ACCEPT = 0, TR_REJ_ERR = 1, TR_REJ_LOWP1 = 4, TR_REJ_LOWP2 = 5, TR_CAP = 7,
       TR_REJ_NONPOS_H = 8, TR_REJ_DRIFT = 9, TR_REJ_ALPHA = 10, TR_REFUSED = 11,
       TR_ACCEPT_PROJECTED = 12 };
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

// h with the alpha_in rule applied: a stage whose (alpha - alpha_in):n is not
// positive takes the model's own sentinel for (alpha - alpha_in):n = 0.
double
ManzariDafalias::ladrunoSasBracketH(const Vector& a, const Vector& ain, const Vector& n, double h)
{
    if (mLadrunoSas.opt.alphaInMode == 2)   // stale: ModifiedEuler's behaviour
        return h;
    Vector t(a);
    t -= ain;
    const double x = DoubleDot2_2_Contr(t, n);
    if (x < small)
        return 1.0e10;
    return h;
}

// alpha / alpha^b with the Lode angle of alpha itself and psi at (e, p(S)).
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
// of alpha is unchanged by the scaling, so ratio == 1 exactly), deviatoric
// stress translated by p*(alpha' - alpha) so dev(S) - p alpha, hence f and n,
// are unchanged. Elastic strain follows the stress change.
void
ManzariDafalias::ladrunoSasProject(Vector& S, Vector& A, Vector& Ee, double e, double K, double G)
{
    const double r = ladrunoSasAlphaRatio(A, S, e);
    if (!(r > 1.0) || !std::isfinite(r))
        return;
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

// One rate-form stage. Returns ST_ELASTIC / ST_PLASTIC / ST_NONPOS_H / ST_TENSION.
int
ManzariDafalias::ladrunoSasStage(const Vector& s, const Vector& a, const Vector& z, double e,
    const Vector& ain, double dv, const Vector& ddev, double K, double G,
    Vector& ds, Vector& da, Vector& dz, Vector& dep, double& lam)
{
    const double p = one3 * GetTrace(s) + m_Presidual;
    if (!(p > 0.0))
        return ST_TENSION;

    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, A, D, B, C;
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.stateDependent");
        GetStateDependent(s, a, z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, A, D, B, C, R);
    }
    OPS_PROFILE_SCOPE("sanisand.sasME.stageArithmetic");
    {
        const double hIn = h;
        h = ladrunoSasBracketH(a, ain, n, h);
        if (h != hIn)   // (alpha - alpha_in):n <= -small: a reversal the stale alpha_in hides
            mLadrunoSas.stats[LSAS_H_BRACKETS] += 1.0;
    }

    Vector r = GetDevPart(s);
    r /= p;
    const double nr = DoubleDot2_2_Contr(n, r);
    const double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
    const double H  = Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n)))) - K * D * nr;
    const double N  = 2.0 * G * DoubleDot2_2_Mixed(n, ddev) - K * dv * nr;

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
    if (!(H > small))                        // loading with H <= 0: no plastic solution
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

// Consistent (Potts-Gens) then normal drift correction. true = |f| <= TolF
// (or f <= TolF when !bothSides). false = neither direction reduces |f|, the
// iteration limit, or tension -- the caller turns that into a cut / refusal.
bool
ManzariDafalias::ladrunoSasDrift(Vector& S, Vector& A, Vector& Z, Vector& Ee, double e,
    const Vector& ain, double K, double G, bool bothSides)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.drift");
    const Matrix aC = GetStiffness(K, G);
    const Matrix aD = GetCompliance(K, G);
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
        double cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C;
        GetStateDependent(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C, R);
        h = ladrunoSasBracketH(A, ain, n, h);

        // exact df/dsigma = n - (1/3)(n:alpha + sqrt(2/3) m) I ; df/dalpha = -p n
        Vector Q(mI1);
        Q *= (-one3 * (DoubleDot2_2_Contr(n, A) + root23 * m_m));
        Q += n;
        Vector Rc = ToCovariant(R);
        Vector dSp = DoubleDot4_2(aC, Rc);
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
            Ee -= DoubleDot4_2(aD, t);
        }
        corrected = true;
    }
    return false;
}

// The continuum elastoplastic tangent at (S, A, Z), loading assumed.
void
ManzariDafalias::ladrunoSasContinuumTangent(const Vector& S, const Vector& A, const Vector& Z,
    const Vector& ain, double e, double K, double G, Matrix& Cep)
{
    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, Af, D, B, C;
    GetStateDependent(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, Af, D, B, C, R);
    h = ladrunoSasBracketH(A, ain, n, h);
    Vector dummy(6);
    Cep = GetElastoPlasticTangent(S, 1.0, dummy, dummy, G, K, B, C, D, h, n, d, b);
}

// The substep loop over the plastic portion curStrain -> nextStrain. Returns 0
// or a refusal code (RC_*). On 0, (S, Ee, A, Z, ain) hold the end state.
int
ManzariDafalias::ladrunoSasSubsteps(Vector& S, Vector& Ee, Vector& A, Vector& Z, Vector& ain,
    const Vector& curStrain, const Vector& nextStrain, double K, double G,
    double& lamSum, bool& lastPlastic)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.substeps");
    double* st = mLadrunoSas.stats;
    const LadrunoSasOptions& o = mLadrunoSas.opt;
    const double TolE   = mTolR;
    const double sigRef = (o.errFloor < 0.0) ? m_P_atm / 101.0 : o.errFloor;
    const double kappa  = o.alphaBoundTol;

    Vector dStrain(nextStrain);
    dStrain -= curStrain;
    const double trD = GetTrace(dStrain);
    const Vector devD = GetDevPart(dStrain);

    Vector ds1(6), da1(6), dz1(6), dep1(6), ds2(6), da2(6), dz2(6), dep2(6);
    Vector S1(6), A1(6), Z1(6), nS(6), nA(6), nZ(6), nEe(6), tmp(6), ddev(6);
    double T = 0.0, dT = 1.0;
    bool lastRejected = false;
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
        const double e = m_e_init - (1 + m_e_init) * GetTrace(tmp);
        const double dv = dT * trD;
        ddev = devD; ddev *= dT;

        double lam1 = 0.0, lam2 = 0.0;
        int k1, k2 = ST_ELASTIC;
        {
            OPS_PROFILE_SCOPE("sanisand.sasME.stages");
            k1 = ladrunoSasStage(S, A, Z, e, ain, dv, ddev, K, G, ds1, da1, dz1, dep1, lam1);
        }
        if (k1 == ST_NONPOS_H) {
            // a property of the accepted start of this substep: no cut can change it
            ladrunoTraceSubstep(T, dT, std::numeric_limits<double>::quiet_NaN(), TR_REFUSED, atMin);
            return RC_NONPOS_H;
        }
        int rej = 0;   // 0 none, else trace code
        if (k1 == ST_TENSION) {
            rej = TR_REJ_LOWP1;
        } else {
            S1 = S; S1 += ds1;
            A1 = A; A1 += da1;
            Z1 = Z; Z1 += dz1;
            if (!(one3 * GetTrace(S1) + m_Presidual > 0.0)) {
                rej = TR_REJ_LOWP1;
            } else {
                {
                    OPS_PROFILE_SCOPE("sanisand.sasME.stages");
                    k2 = ladrunoSasStage(S1, A1, Z1, e, ain, dv, ddev, K, G, ds2, da2, dz2, dep2, lam2);
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

        if (rej == 0) {
            // candidate accepted by the error test: drift, then alpha
            nEe = dStrain; nEe *= dT; nEe += Ee;
            tmp = dep1; tmp += dep2; tmp *= 0.5; nEe -= tmp;
            const bool anyElastic = (k1 == ST_ELASTIC || k2 == ST_ELASTIC);
            const bool anyPlastic = (k1 == ST_PLASTIC || k2 == ST_PLASTIC);
            // the plastic portion belongs ON the surface; an unloading stage may leave it inside
            const bool bothSides = anyPlastic && !anyElastic;
            if (!ladrunoSasDrift(nS, nA, nZ, nEe, e, ain, K, G, bothSides)) {
                rej = TR_REJ_DRIFT;
            } else {
                OPS_PROFILE_SCOPE("sanisand.sasME.alphaCheck");
                tmp = dStrain; tmp *= (T + dT); tmp += curStrain;
                const double eEnd = m_e_init - (1 + m_e_init) * GetTrace(tmp);
                const double ratio = ladrunoSasAlphaRatio(nA, nS, eEnd);
                if (!(ratio <= 1.0 + kappa)) {
                    if (o.alphaProject != 0 && std::isfinite(ratio)) {
                        ladrunoSasProject(nS, nA, nEe, eEnd, K, G);
                        rej = -1;   // accepted, projected
                    } else {
                        rej = TR_REJ_ALPHA;
                    }
                }
                if (rej <= 0) {
                    const double r2 = ladrunoSasAlphaRatio(nA, nS, eEnd);
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
        // G: a reversal inside the increment re-seats alpha_in (Dafalias-Manzari's definition)
        if (o.alphaInMode == 0) {
            Vector nEnd = GetNormalToYield(S, A);
            tmp = A; tmp -= ain;
            if (DoubleDot2_2_Contr(tmp, nEnd) < 0.0) {
                ain = A;
                st[LSAS_ALPHA_IN_RESEATS] += 1.0;
            }
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

    const Vector& CurStress = mSigma_n;
    const Vector& CurStrain = mEpsilon_n;
    const Vector& NextStrain = mEpsilon;

    Vector dStrain(NextStrain);
    dStrain -= CurStrain;
    const double eN = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
    mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);

    double K, G;
    GetElasticModuli(CurStress, eN, K, G);
    const Matrix aC = GetStiffness(K, G);

    Vector S(mSigma_n), A(mAlpha_n), Z(mFabric_n), Ee(mEpsilonE_n), ain(mAlpha_in);
    int code = 0;

    // ---- 0. entry: the committed state must be admissible ------------------
    {
        const double p0 = one3 * GetTrace(S) + m_Presidual;
        const double ta = GetTrace(A), tz = GetTrace(Z);
        if (!finite6(S) || !finite6(A) || !finite6(Z) || !finite6(ain) || !(p0 > 0.0)
            || fabs(ta) > 1.0e-6 * fmax(GetNorm_Contr(A), m_m)
            || fabs(tz) > 1.0e-6 * fmax(GetNorm_Contr(Z), m_m))
            code = RC_START_OTHER;
        else if (GetF(S, A) > mTolF)
            code = RC_START_F;
        else {
            const double r0 = ladrunoSasAlphaRatio(A, S, eN);
            if (!(r0 <= 1.0 + o.alphaBoundTol)) {
                if (o.alphaProject != 0 && std::isfinite(r0))
                    ladrunoSasProject(S, A, Ee, eN, K, G);
                else
                    code = RC_START_ALPHA;
            }
        }
    }

    double lamSum = 0.0;
    bool lastPlastic = false;
    if (code == 0) {
        // ---- 1. elastic predictor / intersection ----------------------------
        double a = 0.0;
        bool elastic = false;
        {
            OPS_PROFILE_SCOPE("sanisand.sasME.predictor");
            Vector dSe = DoubleDot4_2(aC, dStrain);
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
                const double nd = GetNorm_Contr(dSe);
                if (f0 < -mTolF) {
                    a = IntersectionFactor(S, CurStrain, NextStrain, A, 0.0, 1.0);
                    mLadrunoLastPath = 2;
                } else if (DoubleDot2_2_Contr(GetNormalToYield(S, A), dSe) / (nd == 0 ? 1.0 : nd)
                           > (-sqrt(mTolF))) {
                    a = 0.0;
                    mLadrunoLastPath = 3;
                } else {
                    a = IntersectionFactor_Unloading(S, CurStrain, NextStrain, A);
                    mLadrunoLastPath = 4;
                }
                mLadrunoLastElasticRatio = a;
                if (a > 0.0) {
                    Vector dEa(dStrain);
                    dEa *= a;
                    S += DoubleDot4_2(aC, dEa);
                    Ee += dEa;
                    if (fabs(GetF(S, A)) > mTolF)
                        st[LSAS_INTERSECT_FAIL] += 1.0;
                }
            }
        }
        if (elastic) {
            st[LSAS_ELASTIC] += 1.0;
            mSigma = S; mAlpha = A; mFabric = Z; mEpsilonE = Ee;
            mDGamma = 0.0;
            mCe = aC; mCep = aC; mCep_Consistent = aC;
            st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
            st[LSAS_LAST_F] = GetF(S, A);
            return;
        }
        // ---- 2. the plastic portion ------------------------------------------
        Vector curP(dStrain);
        curP *= a;
        curP += CurStrain;
        code = ladrunoSasSubsteps(S, Ee, A, Z, ain, curP, NextStrain, K, G, lamSum, lastPlastic);
    } else {
        mLadrunoLastPath = 7;
        ladrunoTraceSubstep(0.0, 1.0, std::numeric_limits<double>::quiet_NaN(), TR_REFUSED, false);
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
        mDGamma = 0.0;
        mCe = aC; mCep = aC; mCep_Consistent = aC;
        // PROCESS-WIDE warning budget (every Gauss point is an instance)
        static std::atomic<int> ladrunoSasWarnCount(0);   // Ladruno WP-129 (diagnostic budget)
        if (ladrunoSasWarnCount.load() < 10) {
            opserr << "WARNING LadrunoSANISAND (SAS-ME, IntScheme 129) material tag " << this->getTag()
                   << ": update REFUSED (" << refuseName(code) << ", code " << code
                   << "); the trial is left on the committed state and the step must be cut."
                   << " Per-point counts: the `sasStats` response." << endln;
            if (ladrunoSasWarnCount.fetch_add(1) + 1 == 10)
                opserr << "WARNING LadrunoSANISAND SAS-ME: further refusal warnings suppressed"
                          " (budget 10 per process)." << endln;
        }
        return;
    }

    // ---- 3. accepted: state and ONE end-state tangent -----------------------
    mSigma = S; mAlpha = A; mFabric = Z; mEpsilonE = Ee;
    mAlpha_in = ain;
    mDGamma = lamSum;
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.tangent");
        mCe = aC;
        if (lastPlastic)
            ladrunoSasContinuumTangent(S, A, Z, ain, mVoidRatio, K, G, mCep);
        else
            mCep = aC;
        mCep_Consistent = mCep;
    }
    st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
    st[LSAS_LAST_F] = GetF(S, A);
}
