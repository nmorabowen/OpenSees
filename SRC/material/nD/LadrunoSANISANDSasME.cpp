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
// Ladruno band scaling (an OPT-IN variant; -sasBandScale S or -sasBandWidth w c_b,
// OFF by default and the code path then byte-identical). The model softens without
// a length: a shear band in a mesh is as thick as the elements. With the option a
// point that has passed its peak (the LATCH, below) stands for a band of sand of
// thickness w in series with elastically unloading sand, S = w / (c_b h) being the
// share of the element that is band: given by the driver (-sasBandScale S, one
// material per element size), or computed once per point from the element's
// characteristic length h (-sasBandWidth w c_b; see ladrunoSasIntegrate).
// After the latch, with Kp, X, R, D the model's own quantities at the state:
//     lambda = Q:C:d_eps / (Kp + S X)        the plastic multiplier of the BAND
//     d_eps  = C^-1 : d_sigma + S lambda R   the strain of the element
//     alpha, z advance with lambda, as in DM04
//     de     = -(1 + e_init) (d_eps_v + (1 - S) lambda D)   the volume change of
//              the band: the void ratio is e(strain) + an integrated offset
//     tangent  C - S (C:R)(x)(Q:C) / (Kp + S X)   (= the DM04 tangent with h / S)
//     -sasSoftCap kappa becomes H >= kappa S X
// This is exactly a mixture at ONE stress of a fraction S of DM04 and 1 - S of its
// own elastic law (the moduli depend on p and e_init only, so both parts share
// them; the `mUseCurrentVoidRatioInG` seam would break that and no constructor
// sets it). S = 1 is DM04. S > 1 (an element thinner than the band) slows the
// softening. Kp + S X <= 0 is the snap-back of the point itself under full strain
// control: a refusal (loadingWithNonPositiveDenominator), as for DM04.
// THE LATCH is permanent and committed with the state:
//   -sasBandTrigger peak (DEFAULT)  an accepted plastic substep over which the
//        major/minor principal stress ratio fell, ending with b:n <= 0. On a
//        uniform specimen this is the peak of the specimen itself, so first
//        loading and the peak are untouched.
//   -sasBandTrigger bn   an accepted plastic substep ending with b:n <= 0 (Kp <= 0
//        at the point; it acts 0.1-0.5 % of strain before the peak of a
//        plane-strain specimen and moves the softening branch with it).
// Form, latch and limits were chosen on an exact Python integration of DM04 (the
// rate equations to rtol 1e-8, the latch as an event) before this was written:
// on a plane-strain point under constant side stress, against a series model of
// the unmodified band and elastically unloading blocks, the latch `peak` keeps R
// within 3 % for S from 5 to 0.1 at sigma3 5-300 kPa (0.79 % at S 0.2, 1.76 % at
// 0.1, 3.9 % at 0.05); the latch `bn` and a latch on loss of ellipticity miss by
// 5-7 % already at S 0.2 (they move the softening branch to before the peak); a
// cheaper form that scales the void ratio only misses by 8 % at S 0.2. The scale
// at which a constant-side-stress specimen snaps back is S 0.01-0.03.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "UWmaterials/ManzariDafalias.h"
#include <profiler/ProfilerMacros.h>
#include <OPS_Globals.h>
#include <Element.h>                          // Ladruno band scaling: ops_TheActiveElement->getCharacteristicLength()

#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>

// Ladruno sasme-speed (item 1): the update runs on fixed-size stack arrays, not
// on `Vector` / `Matrix` temporaries. The sampling profile of a footing leg put
// heap allocation and release at 20 % of the SAS-ME update and `Vector`
// construction, copies and operators at a further 32 %. Every helper below
// repeats, operation for operation and in the same order, the arithmetic of the
// ManzariDafalias / Vector / Matrix function it stands for (named at each one),
// so the update is BIT-IDENTICAL to the Vector form on an IEEE build without
// -ffast-math and without FMA contraction (the fork's x86-64 builds do not pass
// -march, so there is no FMA to contract into). Signed zeros are kept too: the
// identity's zero shear terms are multiplied, not skipped.

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

const double kOne3 = 1.0 / 3.0;                              // ManzariDafalias::one3, the same expression
const double kI1[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };      // ManzariDafalias::mI1

inline bool finite6(const double* v)
{
    for (int i = 0; i < 6; i++)
        if (!std::isfinite(v[i]))
            return false;
    return true;
}
inline void cp6(double* d, const double* s)
{
    for (int i = 0; i < 6; i++)
        d[i] = s[i];
}
inline double tr6(const double* v)                            // GetTrace
{
    return (v[0] + v[1] + v[2]);
}
inline double ddc6(const double* a, const double* b)          // DoubleDot2_2_Contr
{
    double r = 0.0;
    for (int i = 0; i < 6; i++)
        r += a[i] * b[i] + (i > 2) * a[i] * b[i];
    return r;
}
inline double ddm6(const double* a, const double* b)          // DoubleDot2_2_Mixed
{
    double r = 0.0;
    for (int i = 0; i < 6; i++)
        r += a[i] * b[i];
    return r;
}
inline double norm6(const double* v)                          // GetNorm_Contr
{
    return sqrt(ddc6(v, v));
}
inline void dev6(double* r, const double* v)                  // GetDevPart
{
    const double p = tr6(v);
    cp6(r, v);
    r[0] -= kOne3 * p;
    r[1] -= kOne3 * p;
    r[2] -= kOne3 * p;
}
inline void sdot6(double* r, const double* v1, const double* v2)   // SingleDot
{
    r[0] = v1[0]*v2[0] + v1[3]*v2[3] + v1[5]*v2[5];
    r[1] = v1[3]*v2[3] + v1[1]*v2[1] + v1[4]*v2[4];
    r[2] = v1[5]*v2[5] + v1[4]*v2[4] + v1[2]*v2[2];
    r[3] = 0.5*(v1[0]*v2[3] + v1[3]*v2[0] + v1[3]*v2[1] + v1[1]*v2[3] + v1[5]*v2[4] + v1[4]*v2[5]);
    r[4] = 0.5*(v1[3]*v2[5] + v1[5]*v2[3] + v1[1]*v2[4] + v1[4]*v2[1] + v1[4]*v2[2] + v1[2]*v2[4]);
    r[5] = 0.5*(v1[0]*v2[5] + v1[5]*v2[0] + v1[3]*v2[4] + v1[4]*v2[3] + v1[5]*v2[2] + v1[2]*v2[5]);
}
inline void div6(double* v, double f)                          // Vector::operator/=
{
    if (f == 0.0) {
        for (int i = 0; i < 6; i++)
            v[i] = VECTOR_VERY_LARGE_VALUE;
    } else {
        for (int i = 0; i < 6; i++)
            v[i] /= f;
    }
}
inline double lode6(const double* n)                           // GetLodeAngle: cos(3 theta)
{
    double t[6], u[6];
    sdot6(t, n, n);
    sdot6(u, n, t);
    double c = sqrt(6.0) * tr6(u);
    c = c > 1 ? 1 : c;
    c = c < -1 ? -1 : c;
    return c;
}
// GetStiffness / GetCompliance into a column-major 6x6 (Matrix's layout), and
// Matrix::operator*(Vector) on it (DoubleDot4_2), with Matrix's loop order.
inline void stiff36(double* M, double K, double G)
{
    for (int i = 0; i < 36; i++)
        M[i] = 0.0;
    const double a = K + 4.0*kOne3 * G;
    const double b = K - 2.0*kOne3 * G;
    M[0] = M[7] = M[14] = a;
    M[21] = M[28] = M[35] = G;
    M[6] = M[12] = M[13] = b;    // (0,1) (0,2) (1,2)
    M[1] = M[2] = M[8] = b;      // (1,0) (2,0) (2,1)
}
inline void compl36(double* M, double K, double G)
{
    for (int i = 0; i < 36; i++)
        M[i] = 0.0;
    const double a = 1 / (9*K) + 1 / (3*G);
    const double b = 1 / (9*K) - 1 / (6*G);
    const double c = 1 / G;
    M[0] = M[7] = M[14] = a;
    M[21] = M[28] = M[35] = c;
    M[6] = M[12] = M[13] = b;
    M[1] = M[2] = M[8] = b;
}
inline void mul36(double* r, const double* M, const double* v)
{
    for (int j = 0; j < 6; j++)
        r[j] = 0.0;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            r[j] += M[i*6 + j] * v[i];
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

// ---- Ladruno sasme-speed: the model's own functions on stack arrays ----------

// GetF
double
ManzariDafalias::ladrunoSasF(const double* S, const double* A)
{
    double s[6];
    dev6(s, S);
    const double p = one3 * tr6(S) + m_Presidual;
    for (int i = 0; i < 6; i++)
        s[i] -= A[i] * p;
    return norm6(s) - root23 * m_m * p;
}

// GetNormalToYield
void
ManzariDafalias::ladrunoSasNormal(const double* S, const double* A, double* n)
{
    const double p = one3 * tr6(S) + m_Presidual;
    if (fabs(p) < small) {
        for (int i = 0; i < 6; i++)
            n[i] = 0.0;
        return;
    }
    double dv[6];
    dev6(dv, S);
    for (int i = 0; i < 6; i++) {
        n[i] = A[i] * (-p);
        n[i] += dv[i];
    }
    double normN = norm6(n);
    normN = (normN < small) ? 1.0 : normN;
    div6(n, normN);
}

// GetElasticModuli (the 4-argument overload, called as it is): a Vector that
// WRAPS the stack array (Vector(double*, int) does not allocate or free).
void
ManzariDafalias::ladrunoSasModuli(const double* S, double e, double& K, double& G)
{
    Vector sv(const_cast<double*>(S), 6);
    GetElasticModuli(sv, e, K, G);
}

// GetStateDependent
void
ManzariDafalias::ladrunoSasStateDep(const double* stress, const double* alpha, const double* fabric,
    double e, const double* alpha_in, double* n, double* d, double* b, double& cos3Theta, double& h,
    double& psi, double& alphaBtheta, double& alphaDtheta, double& b0, double& A, double& D,
    double& B, double& C, double* R)
{
    double tmp0[6], tmp1[6], nn[6];
    double D_factor = 1.0;
    double p = one3 * tr6(stress) + m_Presidual;
    p = (p < small) ? small : p;

    ladrunoSasNormal(stress, alpha, n);

    for (int i = 0; i < 6; i++)
        tmp0[i] = alpha[i] - alpha_in[i];
    const double AlphaAlphaInDotN = ddc6(tmp0, n);

    psi = GetPSI(e, p);
    cos3Theta = lode6(n);
    alphaBtheta = g(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
    alphaDtheta = g(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;
    b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);

    for (int i = 0; i < 6; i++) {
        d[i] = n[i] * (root23 * alphaDtheta);
        d[i] -= alpha[i];
        b[i] = n[i] * (root23 * alphaBtheta);
        b[i] -= alpha[i];
    }
    if (fabs(AlphaAlphaInDotN) < small)
        h = 1.0e10;
    else
        h = b0 / AlphaAlphaInDotN;

    A = m_A0 * (1 + Macauley(ddc6(fabric, n)));
    D = A * ddc6(d, n);
    if (p < 0.05 * m_P_atm)
        D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * 101.0 / m_P_atm * p)));   // as GetStateDependent
    else
        D_factor = 1.0;
    D *= D_factor;

    B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(cos3Theta, m_c) * cos3Theta;
    C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(cos3Theta, m_c);

    sdot6(nn, n, n);
    for (int i = 0; i < 6; i++) {
        tmp0[i] = kI1[i] * (-1.0 * one3);
        tmp0[i] += nn[i];
        tmp0[i] *= C;
        tmp1[i] = kI1[i] * (one3 * D);
        R[i] = n[i] * B;
        R[i] -= tmp0[i];
        R[i] += tmp1[i];
    }
}

void
ManzariDafalias::ladrunoResetSasStats(void)
{
    for (int i = 0; i < LSAS_COUNT; i++)
        mLadrunoSas.stats[i] = 0.0;
    mLadrunoSas.refused = false;
    mLadrunoSas.carryH = mLadrunoSas.carryH_n = 0.0;   // Ladruno sasme-speed: revertToStart / replay
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

// Ladruno band scaling: the latch and the void-ratio offset are STATE (reset as
// the separation state is: by revertToStart outside InitialStateAnalysis, and by
// the replay command); the three counters are census.
void
ManzariDafalias::ladrunoResetSasBand(void)
{
    mLadrunoSas.band = mLadrunoSas.band_n = false;
    mLadrunoSas.bandS = mLadrunoSas.bandS_n = 0.0;
    mLadrunoSas.bandDe = mLadrunoSas.bandDe_n = 0.0;
    mLadrunoSas.bandStages = mLadrunoSas.bandCapped = mLadrunoSas.bandLatches = 0.0;
}

// Ladruno band scaling: major / minor principal stress of s (compression
// positive, the true stress: no p_r), from the invariants. With cos(3 theta) =
// sqrt(6) tr(n^3) of the unit deviator n (1 in triaxial compression), the
// principal deviatoric stresses are sqrt(2/3) |s_dev| cos(theta - 2 pi k / 3).
// Returns +inf when the minor one is not positive (such a state never latches).
// Ladruno sasme-speed: on a stack array; GetTrace, GetDevPart, GetNorm_Contr,
// Vector::operator/= and GetLodeAngle as tr6, dev6, norm6, div6, lode6 (the
// Vector form of the band branch, operation for operation).
double
ManzariDafalias::ladrunoSasPrincipalRatio(const double* s)
{
    const double p = one3 * tr6(s);
    double dv[6];
    dev6(dv, s);
    const double nd = norm6(dv);
    if (!(nd > 1.0e-14 * fabs(p)) || !std::isfinite(nd))
        return 1.0;
    div6(dv, nd);
    double c3 = lode6(dv);
    c3 = (c3 > 1.0) ? 1.0 : ((c3 < -1.0) ? -1.0 : c3);
    const double th = acos(c3) / 3.0;
    const double r = root23 * nd;
    const double sMax = p + r * cos(th);
    const double sMin = p + r * cos(th + 2.0943951023931954923);   // + 2 pi / 3
    if (!(sMin > 0.0))
        return std::numeric_limits<double>::infinity();
    return sMax / sMin;
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
ManzariDafalias::ladrunoSasBracketH(const double* a, const double* ain, const double* n, double h,
    double b0)
{
    if (mLadrunoSas.opt.alphaInMode == 2)   // stale: ModifiedEuler's behaviour
        return h;
    double t[6];
    for (int i = 0; i < 6; i++)
        t[i] = a[i] - ain[i];
    const double x = ddc6(t, n);
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
ManzariDafalias::ladrunoSasAlphaRatio(const double* a, const double* s, double e)
{
    const double na = norm6(a);
    if (!(na > 1.0e-14))
        return (na == 0.0 || std::isfinite(na)) ? 0.0 : std::numeric_limits<double>::infinity();
    double u[6];
    cp6(u, a);
    div6(u, na);
    const double c3 = lode6(u);
    double p = one3 * tr6(s) + m_Presidual;
    p = (p < small) ? small : p;
    const double psi = GetPSI(e, p);
    const double aB = g(c3, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
    if (!(aB > 0.0))
        return std::numeric_limits<double>::infinity();
    return sqrt(1.5) * na / aB;
}

// The Vector form, for the replay command and the census (not on the update path).
double
ManzariDafalias::ladrunoSasAlphaRatio(const Vector& a, const Vector& s, double e)
{
    double av[6], sv[6];
    for (int i = 0; i < 6; i++) {
        av[i] = a(i);
        sv[i] = s(i);
    }
    return ladrunoSasAlphaRatio(av, sv, e);
}

// -alphaProject 1: alpha -> s*alpha onto the bounding surface (the Lode angle
// of alpha is unchanged by the scaling, so rho_alpha == 1 exactly), deviatoric
// stress translated by p*(alpha' - alpha) so dev(S) - p alpha, hence f and n,
// are unchanged. Elastic strain follows the stress change.
void
ManzariDafalias::ladrunoSasProject(double* S, double* A, double* Ee, double e)
{
    const double r = ladrunoSasAlphaRatio(A, S, e);
    if (!(r > 1.0) || !std::isfinite(r))
        return;
    double K, G;
    ladrunoSasModuli(S, e, K, G);
    const double p = one3 * tr6(S) + m_Presidual;
    double dA[6], dS[6], M[36], de[6];
    for (int i = 0; i < 6; i++) {
        dA[i] = A[i] * (1.0 / r - 1.0);
        A[i] += dA[i];
        dS[i] = dA[i] * p;
        S[i] += dS[i];
    }
    compl36(M, K, G);
    mul36(de, M, dS);
    for (int i = 0; i < 6; i++)
        Ee[i] += de[i];
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
void
ManzariDafalias::ladrunoSasElastic(const double* S, const double* dEps, double e0, double e1,
    double* dS)
{
    if (mUseCurrentVoidRatioInG) {
        const int nsub = 64;
        double Sx[6], dE[6], d1[6], d2[6], S1[6], M[36];
        cp6(Sx, S);
        cp6(dE, dEps);
        div6(dE, (double)nsub);
        double K, G;
        for (int k = 0; k < nsub; k++) {
            const double ea = e0 + (e1 - e0) * (double)k / nsub;
            const double eb = e0 + (e1 - e0) * (double)(k + 1) / nsub;
            ladrunoSasModuli(Sx, ea, K, G);
            stiff36(M, K, G);
            mul36(d1, M, dE);
            for (int i = 0; i < 6; i++) { S1[i] = Sx[i]; S1[i] += d1[i]; }
            ladrunoSasModuli(S1, eb, K, G);
            stiff36(M, K, G);
            mul36(d2, M, dE);
            for (int i = 0; i < 6; i++) {
                d1[i] += d2[i];
                d1[i] *= 0.5;
                Sx[i] += d1[i];
            }
        }
        for (int i = 0; i < 6; i++)
            dS[i] = Sx[i] - S[i];
        return;
    }
    double K0, G0;
    ladrunoSasModuli(S, e0, K0, G0);
    const double x0 = one3 * tr6(S) + m_PreElastic;
    const double xe0 = (x0 <= m_Pmin) ? m_Pmin : x0;
    const double g = G0 / sqrt(xe0);
    const double c = K0 / G0;
    const double dv = tr6(dEps);
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
    double dd[6];
    dev6(dd, dEps);
    dd[3] *= 0.5;   // ToContraviant
    dd[4] *= 0.5;
    dd[5] *= 0.5;
    for (int i = 0; i < 6; i++) {
        dS[i] = dd[i] * (2.0 * I);
        dS[i] += kI1[i] * (x - x0);
    }
}

// Where the EXACT elastic path from S along dEps meets the yield surface of A:
// Pegasus on [lo, hi] with f(lo) < 0 < f(hi) (the same path the predictor uses,
// so the plastic portion starts ON the surface -- review numerics item 3).
double
ManzariDafalias::ladrunoSasIntersect(const double* S, const double* A, const double* dEps,
    double e0, double lo, double hi)
{
    double dE[6], St[6], dSe[6];
    auto fAt = [&](double a) {
        for (int i = 0; i < 6; i++) dE[i] = dEps[i] * a;
        ladrunoSasElastic(S, dE, e0, e0, dSe);
        for (int i = 0; i < 6; i++) { St[i] = S[i]; St[i] += dSe[i]; }
        return ladrunoSasF(St, A);
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
ManzariDafalias::ladrunoSasStage(const double* s, const double* a, const double* z, double e,
    const double* ain, double dv, const double* ddev,
    double* ds, double* da, double* dz, double* dep, double& lam,
    double bandS, double* lamD)
{
    // Ladruno band scaling: bandS > 0 = a latched point, a fraction bandS of band in
    // series with elastic sand. lam is then the multiplier of the BAND (alpha, z),
    // the plastic strain of the element is bandS lam R, and *lamD returns lam D,
    // the plastic volume change of the band. bandS == 0: the WP-129 stage, unchanged.
    const bool bandOn = (bandS > 0.0);
    if (lamD != 0)
        *lamD = 0.0;
    const double p = one3 * tr6(s) + m_Presidual;
    if (!(p > 0.0))
        return ST_TENSION;

    double n[6], d[6], b[6], R[6];
    double cos3Theta, h, psi, aB, aD, b0, A, D, B, C, K, G;
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.stateDependent");
        ladrunoSasStateDep(s, a, z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, A, D, B, C, R);
        ladrunoSasModuli(s, e, K, G);                      // U9: the stage's own moduli
    }
    OPS_PROFILE_SCOPE("sanisand.sasME.stageArithmetic");
    {
        const double hIn = h;
        h = ladrunoSasBracketH(a, ain, n, h, b0);
        if (mLadrunoSas.opt.hFloor > 0.0) {                               // Ladruno WP-151
            double t[6];
            for (int i = 0; i < 6; i++)
                t[i] = a[i] - ain[i];
            if (ddc6(t, n) < mLadrunoSas.opt.hFloor * root23 * m_m)
                mLadrunoSas.stats[LSAS_H_FLOORED] += 1.0;
        } else if (h != hIn)   // (alpha - alpha_in):n <= -small: a reversal the stale alpha_in hides
            mLadrunoSas.stats[LSAS_H_BRACKETS] += 1.0;
    }

    // U10: the TRUE yield gradient Q = n - (1/3)(n:alpha + sqrt(2/3) m) I, so
    // Q : C : d_eps = 2G n:de_dev - K dv (n:alpha + sqrt(2/3) m)
    double nn[6], nnn[6];
    sdot6(nn, n, n);
    sdot6(nnn, n, nn);
    const double trn3 = tr6(nnn);                  // GetTrace(SingleDot(n, SingleDot(n, n)))
    const double qv = ddc6(n, a) + root23 * m_m;
    if (mLadrunoSas.opt.softCap > 0.0) {                                  // Ladruno WP-151
        const double X = 2.0 * G * (B - C * trn3) - K * D * qv;
        const double hIn = h;
        h = ladrunoSasSoftCapH(h, ddc6(b, n), p, bandOn ? bandS * X : X);   // band: H >= kappa S X
        if (h != hIn) {
            mLadrunoSas.stats[LSAS_H_SOFTCAPPED] += 1.0;
            if (bandOn)
                mLadrunoSas.bandCapped += 1.0;                            // Ladruno band scaling
        }
    }
    const double Kp = two3 * p * h * ddc6(b, n);
    double H  = Kp + 2.0 * G * (B - C * trn3) - K * D * qv;
    if (bandOn)   // Ladruno band scaling: H = Kp + S X (the elastic part belongs to the whole element)
        H = Kp + bandS * (2.0 * G * (B - C * trn3) - K * D * qv);
    const double N  = 2.0 * G * ddm6(n, ddev) - K * dv * qv;

    // elastic trial part, common to both kinds: ToContraviant(ddev) 2G + K dv I
    for (int i = 0; i < 6; i++) {
        double c = ddev[i];
        if (i > 2)
            c *= 0.5;
        ds[i] = c * (2.0 * G);
        ds[i] += kI1[i] * (K * dv);
    }

    if (!(N > 0.0)) {                       // unloading / neutral: elastic stage
        for (int i = 0; i < 6; i++)
            da[i] = dz[i] = dep[i] = 0.0;
        lam = 0.0;
        mLadrunoSas.stats[LSAS_ELASTIC_STAGES] += 1.0;
        return ST_ELASTIC;
    }
    if (!(H > small) || !std::isfinite(H))   // loading with H <= 0: no plastic solution
        return ST_NONPOS_H;

    lam = N / H;
    // Ladruno band scaling: the plastic strain of the ELEMENT is S lam R (OFF: lam)
    const double lamE = (bandOn ? bandS * lam : lam);
    // C:R = 2G (B n - C (n.n - I/3)) + K D I
    const double fz = (-1.0 * lam * m_cz * Macauley(-1.0 * D));
    for (int i = 0; i < 6; i++) {
        double cr = n[i] * B;
        double t2 = kI1[i] * (-one3);
        t2 += nn[i];
        t2 *= C;
        cr -= t2;
        cr *= (2.0 * G);
        cr += kI1[i] * (K * D);
        cr *= lamE;
        ds[i] -= cr;

        da[i] = b[i] * (lam * two3 * h);

        dz[i] = n[i] * m_z_max;
        dz[i] += z[i];
        dz[i] *= fz;

        double r = R[i];
        if (i > 2)
            r *= 2.0;                       // ToCovariant
        dep[i] = r * lamE;
    }
    if (bandOn)
        mLadrunoSas.bandStages += 1.0;      // Ladruno band scaling
    if (lamD != 0)
        *lamD = lam * D;                    // Ladruno band scaling: the plastic volume change of the band
    return ST_PLASTIC;
}

// Consistent (Potts-Gens) then normal drift correction, moduli at the current
// state. true = |f| <= TolF (or f <= TolF when !bothSides). false = neither
// direction reduces |f|, the iteration limit, or tension -- the caller turns
// that into a cut / refusal, never into an accepted f > TolF.
bool
ManzariDafalias::ladrunoSasDrift(double* S, double* A, double* Z, double* Ee, double e,
    const double* ain, bool bothSides, double bandS)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.drift");
    // Ladruno band scaling: on a latched point the correction moves the element by
    // S lam R and the alpha, z of the band by lam: lam = f / (S Q:C:R + p n:abar).
    // The void-ratio offset is not corrected (lam is of the order of TolF / H).
    // sB = 1.0 when OFF: every product below is then exact, the WP-129 arithmetic.
    const double sB = (bandS > 0.0) ? bandS : 1.0;
    bool corrected = false;
    for (int it = 0; it <= kSasDriftIter; it++) {
        const double f = ladrunoSasF(S, A);
        if (!std::isfinite(f))
            return false;
        if ((bothSides ? fabs(f) : f) <= mTolF) {
            if (corrected)
                mLadrunoSas.stats[LSAS_DRIFT_CORRECTIONS] += 1.0;
            return true;
        }
        if (it == kSasDriftIter)
            return false;
        const double p = one3 * tr6(S) + m_Presidual;
        if (!(p > 0.0))
            return false;

        double n[6], d[6], b[6], R[6], M[36];
        double cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C, K, G;
        ladrunoSasStateDep(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD2, b0, Af, D, B, C, R);
        ladrunoSasModuli(S, e, K, G);
        h = ladrunoSasBracketH(A, ain, n, h, b0);
        stiff36(M, K, G);

        // exact df/dsigma = n - (1/3)(n:alpha + sqrt(2/3) m) I ; df/dalpha = -p n
        double Q[6], Rc[6], dSp[6], aBar[6], zBar[6];
        const double qs = (-one3 * (ddc6(n, A) + root23 * m_m));
        for (int i = 0; i < 6; i++) {
            Q[i] = kI1[i] * qs;
            Q[i] += n[i];
            Rc[i] = R[i];
        }
        Rc[3] *= 2.0;   // ToCovariant
        Rc[4] *= 2.0;
        Rc[5] *= 2.0;
        mul36(dSp, M, Rc);
        if (mLadrunoSas.opt.softCap > 0.0)                                // Ladruno WP-151: X = Q:C:R here
            h = ladrunoSasSoftCapH(h, ddc6(b, n), p, sB * ddc6(Q, dSp));
        const double fz = (-1.0 * m_cz * Macauley(-1.0 * D));
        for (int i = 0; i < 6; i++) {
            aBar[i] = b[i] * (two3 * h);
            zBar[i] = n[i] * m_z_max;
            zBar[i] += Z[i];
            zBar[i] *= fz;
        }

        bool moved = false;
        const double den = sB * ddc6(Q, dSp) + p * ddc6(n, aBar);
        if (den > small && std::isfinite(den)) {
            const double lam = f / den;
            double S1[6], A1[6];
            for (int i = 0; i < 6; i++) {
                S1[i] = S[i];
                S1[i] -= dSp[i] * (sB * lam);
                A1[i] = A[i];
                A1[i] += aBar[i] * lam;
            }
            const double f1 = ladrunoSasF(S1, A1);
            if (std::isfinite(f1) && fabs(f1) < fabs(f) && one3 * tr6(S1) + m_Presidual > 0.0) {
                for (int i = 0; i < 6; i++) {
                    S[i] = S1[i];
                    A[i] = A1[i];
                    Z[i] += zBar[i] * lam;
                    Ee[i] -= Rc[i] * (sB * lam);
                }
                moved = true;
            }
        }
        if (!moved) {
            const double qq = ddc6(Q, Q);
            if (!(qq > 0.0))
                return false;
            const double lam = f / qq;
            double S1[6], t[6], de[6];
            for (int i = 0; i < 6; i++) {
                t[i] = Q[i] * lam;
                S1[i] = S[i];
                S1[i] -= t[i];
            }
            const double f1 = ladrunoSasF(S1, A);
            if (!(std::isfinite(f1) && fabs(f1) < fabs(f) && one3 * tr6(S1) + m_Presidual > 0.0))
                return false;    // the give-up: FAIL, never hand back f > TolF
            cp6(S, S1);
            compl36(M, K, G);
            mul36(de, M, t);
            for (int i = 0; i < 6; i++)
                Ee[i] -= de[i];
        }
        corrected = true;
    }
    return false;
}

// The continuum elastoplastic tangent at (S, A, Z), plastic loading assumed,
// moduli at that state. Once per accepted update (not in the substep loop), so it
// keeps the Vector form.
void
ManzariDafalias::ladrunoSasContinuumTangent(const Vector& S, const Vector& A, const Vector& Z,
    const Vector& ain, double e, Matrix& Cep, double bandS)
{
    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, Af, D, B, C, K, G;
    GetStateDependent(S, A, Z, e, ain, n, d, b, cos3Theta, h, psi, aB, aD, b0, Af, D, B, C, R);
    GetElasticModuli(S, e, K, G);
    {
        double av[6], ainv[6], nv[6];
        for (int i = 0; i < 6; i++) { av[i] = A(i); ainv[i] = ain(i); nv[i] = n(i); }
        h = ladrunoSasBracketH(av, ainv, nv, h, b0);
    }
    if (mLadrunoSas.opt.softCap > 0.0) {                                  // Ladruno WP-151: the stage's X
        const double p = one3 * GetTrace(S) + m_Presidual;
        const double qv = DoubleDot2_2_Contr(n, A) + root23 * m_m;
        const double X = 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n)))) - K * D * qv;
        h = ladrunoSasSoftCapH(h, DoubleDot2_2_Contr(b, n), p, (bandS > 0.0) ? bandS * X : X);
    }
    // Ladruno band scaling: C - S (C:R)(x)(Q:C)/(Kp + S X) is the DM04 tangent with
    // Kp / S, i.e. with h / S (h enters GetElastoPlasticTangent through Kp only).
    if (bandS > 0.0)
        h /= bandS;
    Vector dummy(6);
    Cep = GetElastoPlasticTangent(S, 1.0, dummy, dummy, G, K, B, C, D, h, n, d, b);
}

// The substep loop over the plastic portion curStrain -> nextStrain. Returns 0
// or a refusal code (RC_*). On 0, (S, Ee, A, Z, ain) hold the end state.
// `onset`: the portion starts at a plastic onset (paper alpha_in rule).
int
ManzariDafalias::ladrunoSasSubsteps(double* S, double* Ee, double* A, double* Z, double* ain,
    const double* curStrain, const double* nextStrain, bool onset,
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
    const bool bandOpt  = o.bandOn();                // Ladruno band scaling: OFF = every value below is the WP-129 one

    double dStrain[6], devD[6];
    for (int i = 0; i < 6; i++) {
        dStrain[i] = nextStrain[i];
        dStrain[i] -= curStrain[i];
    }
    const double trD = tr6(dStrain);
    dev6(devD, dStrain);

    double ds1[6], da1[6], dz1[6], dep1[6], ds2[6], da2[6], dz2[6], dep2[6];
    double S1[6], A1[6], Z1[6], nS[6], nA[6], nZ[6], nEe[6], tmp[6], ddev[6];
    double ain1[6], ain2[6], tmpN[6], ainEnd[6];
    double T = 0.0, dT = 1.0;
    bool lastRejected = false;
    // Ladruno sasme-speed (-sasCarryDT 1): the first trial is the substep the
    // controller proposed at the end of this point's last COMMITTED plastic
    // update, carried as a strain norm (the portions of two updates differ in
    // size). OFF: dn = 0 and dT = 1, every line below as before. dTprop is the
    // controller's proposal before the cut to 1 - T; it is bookkeeping only and
    // never changes dT.
    const bool carry = (o.carryDT != 0);
    const double dn = carry ? sqrt(ddm6(dStrain, dStrain)) : 0.0;
    if (carry && mLadrunoSas.carryH_n > 0.0 && dn > 0.0 && std::isfinite(dn))
        dT = fmax(fmin(1.0, mLadrunoSas.carryH_n / dn), kSasDTmin);
    double dTprop = dT;
    // Ladruno sasme-speed (TanType 3): a normal update logs its accepted dT; the
    // perturbed replays of ladrunoSasFDTangent take exactly those dT, accept
    // every substep the error test would judge (the controller is FROZEN, so the
    // map is smooth in the strain), and fail on any other rejection.
    const bool replay = (mLadrunoSas.fdMode == 2);
    const bool logDT  = (!replay && (int)mTangType == 3);
    std::vector<double>& dtLog = mLadrunoSas.dtLog;
    size_t kRep = 0;
    if (logDT)
        dtLog.clear();
    if (replay) {
        if (dtLog.empty())
            return RC_DTMIN;
        dT = dtLog[0];
    }
    int nAcc = 0;  // Ladruno band scaling x TanType 3: accepted substeps of this update
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

        for (int i = 0; i < 6; i++) { tmp[i] = dStrain[i] * T; tmp[i] += curStrain[i]; }
        const double e0s = m_e_init - (1 + m_e_init) * tr6(tmp);
        for (int i = 0; i < 6; i++) { tmp[i] = dStrain[i] * (T + dT); tmp[i] += curStrain[i]; }
        const double e1s = m_e_init - (1 + m_e_init) * tr6(tmp);
        // Ladruno band scaling: on a latched point the void ratio is that of the
        // band: e(strain) + the offset integrated so far. Stage 1 sees it at the
        // substep start; stage 2 at the end as stage 1 predicts it (Euler); the
        // accepted value is their mean (Heun), like every other variable.
        const bool   bandOn = mLadrunoSas.band;
        const double bS = bandOn ? mLadrunoSas.bandS : 0.0;
        const double cE = bandOn ? -(1 + m_e_init) * (1.0 - bS) : 0.0;
        const double e0 = bandOn ? e0s + mLadrunoSas.bandDe : e0s;
        double e1 = bandOn ? e1s + mLadrunoSas.bandDe : e1s;
        double ld1 = 0.0, ld2 = 0.0;   // lam D of the two stages (plastic volume change of the band)
        const double rat0 = (bandOpt && !bandOn && o.bandTrigger == 0) ? ladrunoSasPrincipalRatio(S) : 0.0;
        const double dv = dT * trD;
        for (int i = 0; i < 6; i++)
            ddev[i] = devD[i] * dT;

        double lam1 = 0.0, lam2 = 0.0;
        int k1, k2 = ST_ELASTIC;
        // paper rule (DM04, = WP-134's decide()): at a stage where
        // (alpha - alpha_in):n < 0 a new loading process starts there:
        // alpha_in := alpha at THAT stage's state.
        // Only ON the surface (the oracle decides modes, and re-seats, only
        // there -- review numerics item 7).
        bool reseat1 = false, reseat2 = false;
        const bool onSurface = (fabs(ladrunoSasF(S, A)) <= mTolF);
        cp6(ain1, ain);
        if (paperRule && onSurface) {
            ladrunoSasNormal(S, A, tmpN);
            for (int i = 0; i < 6; i++) tmp[i] = A[i] - ain[i];
            const double x1 = ddc6(tmp, tmpN);
            if (x1 < -dRev) { cp6(ain1, A); reseat1 = true; }
            else if (x1 < 0.0) st[LSAS_RESEAT_HELD] += 1.0;   // Ladruno WP-151: a sub-threshold reversal
        }
        {
            OPS_PROFILE_SCOPE("sanisand.sasME.stages");
            k1 = ladrunoSasStage(S, A, Z, e0, ain1, dv, ddev, ds1, da1, dz1, dep1, lam1, bS, &ld1);
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
            for (int i = 0; i < 6; i++) {
                S1[i] = S[i]; S1[i] += ds1[i];
                A1[i] = A[i]; A1[i] += da1[i];
                Z1[i] = Z[i]; Z1[i] += dz1[i];
            }
            if (!(one3 * tr6(S1) + m_Presidual > 0.0)) {
                rej = TR_REJ_LOWP1;
            } else {
                cp6(ain2, ain1);
                if (paperRule && k1 == ST_PLASTIC) {
                    // a reversal INSIDE a plastic substep: locate it at a substep
                    // start (cut) rather than re-seating on the Euler predictor;
                    // at dT_min re-seat at the stage state.
                    ladrunoSasNormal(S1, A1, tmpN);
                    for (int i = 0; i < 6; i++) tmp[i] = A1[i] - ain1[i];
                    if (ddc6(tmp, tmpN) < -dRev) {   // Ladruno WP-151: -dRev (0 = paper)
                        if (atMin) { cp6(ain2, A1); reseat2 = true; }
                        else rej = TR_REJ_REVERSAL;
                    }
                }
                if (rej == 0) {
                    OPS_PROFILE_SCOPE("sanisand.sasME.stages");
                    if (bandOn)
                        e1 += cE * ld1;   // Ladruno band scaling: the void ratio of the band at the predictor
                    k2 = ladrunoSasStage(S1, A1, Z1, e1, ain2, dv, ddev, ds2, da2, dz2, dep2, lam2, bS, &ld2);
                }
                if (k2 == ST_TENSION)
                    rej = TR_REJ_LOWP2;
                else if (k2 == ST_NONPOS_H)
                    rej = TR_REJ_NONPOS_H;
            }
        }
        double err = std::numeric_limits<double>::quiet_NaN();
        if (rej == 0) {
            for (int i = 0; i < 6; i++) {
                nS[i] = ds1[i]; nS[i] += ds2[i]; nS[i] *= 0.5; nS[i] += S[i];
                nA[i] = da1[i]; nA[i] += da2[i]; nA[i] *= 0.5; nA[i] += A[i];
                nZ[i] = dz1[i]; nZ[i] += dz2[i]; nZ[i] *= 0.5; nZ[i] += Z[i];
            }
            if (!(one3 * tr6(nS) + m_Presidual > 0.0))
                rej = TR_REJ_LOWP2;
            else {
                for (int i = 0; i < 6; i++) tmp[i] = ds2[i] - ds1[i];
                err = norm6(tmp) / fmax(2.0 * norm6(S), sigRef);
                if (o.errorVars == 0) {
                    for (int i = 0; i < 6; i++) tmp[i] = da2[i] - da1[i];
                    err = fmax(err, norm6(tmp) / fmax(2.0 * norm6(A), 1.0));
                    for (int i = 0; i < 6; i++) tmp[i] = dz2[i] - dz1[i];
                    err = fmax(err, norm6(tmp) / fmax(2.0 * norm6(Z), 1.0));
                }
                if (bandOn)   // Ladruno band scaling: the void-ratio offset is an integrated variable too
                    err = fmax(err, fabs(cE * (ld2 - ld1)) / fmax(2.0 * fabs(e0), 1.0));
                if (!std::isfinite(err) || !finite6(nS) || !finite6(nA) || !finite6(nZ))
                    err = std::numeric_limits<double>::infinity();
                err = fmax(err, DBL_EPSILON);
                if (err > TolE && !replay)
                    rej = TR_REJ_ERR;
            }
        }

        // the alpha_in this substep ends with (paper rule: the re-seat point)
        cp6(ainEnd, reseat2 ? ain2 : (reseat1 ? ain1 : ain));
        // Ladruno band scaling: the void ratio this substep ends with (OFF: e1, unchanged)
        const double eEnd = bandOn ? e1s + mLadrunoSas.bandDe + 0.5 * cE * (ld1 + ld2) : e1;
        if (rej == 0) {
            // candidate accepted by the error test: drift, then alpha
            for (int i = 0; i < 6; i++) {
                nEe[i] = dStrain[i] * dT; nEe[i] += Ee[i];
                tmp[i] = dep1[i]; tmp[i] += dep2[i]; tmp[i] *= 0.5; nEe[i] -= tmp[i];
            }
            const bool anyElastic = (k1 == ST_ELASTIC || k2 == ST_ELASTIC);
            const bool anyPlastic = (k1 == ST_PLASTIC || k2 == ST_PLASTIC);
            // the plastic portion belongs ON the surface; an unloading stage may leave it inside
            // ... and only when the substep STARTED on it: a start inside (f < -TolF)
            // is corrected only if it ends outside (review numerics item 3)
            const bool bothSides = anyPlastic && !anyElastic && onSurface;
            if (!ladrunoSasDrift(nS, nA, nZ, nEe, eEnd, ainEnd, bothSides, bS)) {
                rej = TR_REJ_DRIFT;
            } else {
                OPS_PROFILE_SCOPE("sanisand.sasME.alphaCheck");
                // Rejected only when PLASTIC FLOW carried alpha outward beyond
                // 1 + kappa -- measured against the END surface with the start
                // alpha, so psi moving the surface is never the reason (review
                // numerics item 2: that was a dead end).
                const double ratio = ladrunoSasAlphaRatio(nA, nS, eEnd);
                const double ratioFixed = ladrunoSasAlphaRatio(A, nS, eEnd);
                if (!(ratio <= 1.0 + kappa) && !(ratio <= ratioFixed)) {
                    if (o.alphaProject != 0 && std::isfinite(ratio)) {
                        ladrunoSasProject(nS, nA, nEe, eEnd);
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

        if (rej > 0 && replay)
            return RC_DTMIN;   // Ladruno sasme-speed: a frozen replay cannot cut
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
            dTprop = dT;
            lastRejected = true;
            continue;
        }

        // accept
        st[LSAS_ACCEPTED] += 1.0;
        ++nAcc;
        ladrunoTraceSubstep(T, dT, err, rej < 0 ? TR_ACCEPT_PROJECTED : TR_ACCEPT, atMin);
        cp6(S, nS); cp6(A, nA); cp6(Z, nZ); cp6(Ee, nEe);
        if (bandOn)
            mLadrunoSas.bandDe += 0.5 * cE * (ld1 + ld2);   // Ladruno band scaling: the volume change of the band
        lamSum += 0.5 * (lam1 + lam2);
        lastPlastic = (k2 == ST_PLASTIC);   // the loading index nearest the end state
        if (reseat1 || reseat2) {
            cp6(ain, ainEnd);
            st[LSAS_ALPHA_IN_RESEATS] += 1.0;
        }
        // (alpha - alpha_in):n reaching 0 re-seats alpha_in (the paper's rule),
        // on the surface only
        if (paperRule && fabs(ladrunoSasF(S, A)) <= mTolF) {
            ladrunoSasNormal(S, A, tmpN);
            for (int i = 0; i < 6; i++) tmp[i] = A[i] - ain[i];
            const double xEnd = ddc6(tmp, tmpN);
            if (xEnd < -dRev) {                                          // Ladruno WP-151: -dRev (0 = paper)
                cp6(ain, A);
                st[LSAS_ALPHA_IN_RESEATS] += 1.0;
            } else if (xEnd < 0.0)
                st[LSAS_RESEAT_HELD] += 1.0;                             // Ladruno WP-151
        }
        // Ladruno band scaling: THE LATCH, tested on an accepted, fully plastic
        // substep of a point not yet latched (see the note at the top of the file).
        // Later substeps of this update already use the scale; the latch becomes
        // permanent at commit and is undone with a refused or reverted trial.
        const double sLatch = (o.bandScale > 0.0) ? o.bandScale : mLadrunoSas.bandSElem;
        if (bandOpt && !bandOn && sLatch > 0.0 && k1 == ST_PLASTIC && k2 == ST_PLASTIC) {
            bool trig = true;
            if (o.bandTrigger == 0)
                trig = (ladrunoSasPrincipalRatio(S) < rat0);
            if (trig) {
                double nL[6], dL[6], bL[6], RL[6];
                double c3L, hL, psiL, aBL, aDL, b0L, AL, DL, BL, CL;
                ladrunoSasStateDep(S, A, Z, eEnd, ain, nL, dL, bL, c3L, hL, psiL, aBL, aDL, b0L, AL, DL, BL, CL, RL);
                if (ddc6(bL, nL) <= 0.0) {
                    mLadrunoSas.band = true;
                    mLadrunoSas.bandS = sLatch;   // -sasBandWidth: w/(c_b h), read at the first update
                    mLadrunoSas.bandLatches += 1.0;
                    mLadrunoSas.bandLatchAt = nAcc;   // TanType 3: the replays must latch here too
                }
            }
        }
        T += dT;
        double q = fmin(fmax(0.9 * sqrt(TolE / err), 0.1), 1.1);
        if (lastRejected)
            q = fmin(q, 1.0);
        lastRejected = false;
        const bool cut = (dT < dTprop);   // this substep was the remainder 1 - T, not a proposal
        if (logDT)
            dtLog.push_back(dT);
        dT = fmax(q * dT, kSasDTmin);
        dTprop = cut ? fmax(dT, dTprop) : dT;
        dT = fmin(dT, 1.0 - T);
        if (replay) {                       // Ladruno sasme-speed: the next logged dT, or done
            if (++kRep < dtLog.size())
                dT = dtLog[kRep];
            else
                break;
        }
    }
    if (carry)
        mLadrunoSas.carryH = dTprop * dn;   // Ladruno sasme-speed: the proposal, in strain norm
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
    mLadrunoSas.carryH = 0.0;   // Ladruno sasme-speed: set again only by an accepted plastic update
    const double reseats0 = st[LSAS_ALPHA_IN_RESEATS];
    // Ladruno band scaling: the trial starts from the committed latch and offset
    mLadrunoSas.band = mLadrunoSas.band_n;
    mLadrunoSas.bandS = mLadrunoSas.bandS_n;
    mLadrunoSas.bandDe = mLadrunoSas.bandDe_n;
    mLadrunoSas.bandLatchAt = -1;

    double CurStrain[6], NextStrain[6], dStrain[6];
    for (int i = 0; i < 6; i++) {
        CurStrain[i] = mEpsilon_n(i);
        NextStrain[i] = mEpsilon(i);
        dStrain[i] = NextStrain[i] - CurStrain[i];
    }
    double eN = m_e_init - (1 + m_e_init) * tr6(CurStrain);
    mVoidRatio = m_e_init - (1 + m_e_init) * tr6(NextStrain);
    if (mLadrunoSas.band_n) {   // Ladruno band scaling: a latched point carries the void ratio of the band
        eN += mLadrunoSas.bandDe_n;
        mVoidRatio += mLadrunoSas.bandDe_n;   // the end value of an elastic update; the plastic portion updates it below
    }

    // ---- Ladruno WP-152: the tension cutoff (separation) ----------------------
    // The trial starts as the committed state. A SEPARATED point carries no
    // tension and no shear. With g = tr(eps) - tr(eps_entry) (compression
    // positive) its model p is p_c(g) = p_min + K(p_contact) max(g, 0): OPEN
    // (g <= 0) it sits at p_min and absorbs the strain; CLOSING (g > 0) it
    // reloads isotropically and elastically -- continuous in g (review #9: the
    // first build jumped p_min -> p_contact at g_c, which leaves a band of no
    // equilibrium wherever the point has compliance around it; a free-node Newton
    // column failed there). It re-contacts at g >= g_c = (p_contact - p_min) /
    // K(p_contact), with p_re = p_c(g) >= p_contact, alpha = alpha_in = 0, and
    // SAS-ME resumes at the next update. Re-contact is VOLUMETRIC only: isochoric
    // shear never closes the gap (review #6, documented).
    // Ladruno_implementation/152_sanisand_tension_cutoff.md.
    mLadrunoSas.sep = mLadrunoSas.sep_n;
    mLadrunoSas.sepTr = mLadrunoSas.sepTr_n;
    mLadrunoSas.sepEvent = 0;   // an element may update several times per step: the census
                                // counts the COMMITTED transition (LadrunoSANISAND::commitState)
    mLadrunoSas.sepCode = 0;
    mLadrunoSas.sepP0 = 0.0;
    const bool tcOn = (o.tcPcontact > 0.0);
    if (tcOn && mLadrunoSas.sep_n) {
        const double g = tr6(NextStrain) - mLadrunoSas.sepTr_n;
        Vector Sc(mI1);
        Sc *= (o.tcPcontact - m_Presidual);
        double Kc, Gc;
        GetElasticModuli(Sc, mVoidRatio, Kc, Gc);
        const double gc = (o.tcPcontact - m_Pmin) / Kc;
        ladrunoSasSetIsotropic(m_Pmin + Kc * fmax(g, 0.0));
        if (g > 0.0) {
            // closing, and the exit iterate: within THIS update sigma = (p_min + Kc g) I,
            // so the consistent bulk modulus is K(p_contact) (re-review MED-3). The shear
            // stiffness G(p) is a declared regularisation (the stress carries no shear);
            // SAS-ME's next update starts from C_e at p_re on its own.
            double Kp, Gp;
            GetElasticModuli(mSigma, mVoidRatio, Kp, Gp);
            mCe = GetStiffness(Kc, Gp); mCep = mCe; mCep_Consistent = mCe;
        }
        if (g >= gc) {
            mLadrunoSas.sep = false;
            mLadrunoSas.sepEvent = 3;
            mLadrunoLastPath = 9;
        } else {
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
    double S[6], A[6], Z[6], Ee[6], ain[6];
    for (int i = 0; i < 6; i++) {
        S[i] = mSigma_n(i);
        A[i] = mAlpha_n(i);
        Z[i] = mFabric_n(i);
        Ee[i] = mEpsilonE_n(i);
        ain[i] = paperRule ? mAlpha_in_n(i) : mAlpha_in(i);
    }
    int code = 0;
    const double p0c = one3 * tr6(S) + m_Presidual;   // Ladruno WP-152: committed p (the E2 test)
    bool startTension = false;                        // Ladruno WP-152: code 3 BECAUSE p0 <= 0 (E1)

    // ---- 0. entry: the committed state must be admissible ------------------
    {
        const double p0 = one3 * tr6(S) + m_Presidual;
        const double ta = tr6(A), tz = tr6(Z);
        if (!finite6(S) || !finite6(A) || !finite6(Z) || !finite6(ain) || !(p0 > 0.0)
            || fabs(ta) > 1.0e-6 * fmax(norm6(A), m_m)
            || fabs(tz) > 1.0e-6 * fmax(norm6(Z), m_m)) {
            code = RC_START_OTHER;
            // Ladruno WP-152: only "finite, traces clean, p0 <= 0" is a tension
            // start; a non-finite value or a trace violation stays a loud refusal.
            startTension = finite6(S) && finite6(A) && finite6(Z) && finite6(ain) && !(p0 > 0.0)
                && !(fabs(ta) > 1.0e-6 * fmax(norm6(A), m_m))
                && !(fabs(tz) > 1.0e-6 * fmax(norm6(Z), m_m));
        }
        else if (ladrunoSasF(S, A) > mTolF)
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

    // Ladruno band scaling (-sasBandWidth w c_b): S = w/(c_b h), h the
    // characteristic length of the element that drives this point, read ONCE, on
    // the first SAS-ME update, from the element under update (Domain::update sets
    // ops_TheActiveElement before Element::update; LadrunoQuad hands sqrt(area),
    // LadrunoBrick cbrt(volume); the crack-band materials read it the same way).
    // No usable length (a material-point driver, a degenerate element) REFUSES:
    // the option must not fall back to DM04 silently.
    if (code == 0 && !(o.bandScale > 0.0) && o.bandOn() && !(mLadrunoSas.bandSElem > 0.0)) {
        const double h = (ops_TheActiveElement != 0) ? ops_TheActiveElement->getCharacteristicLength() : 0.0;
        const double sE = (std::isfinite(h) && h > 0.0) ? o.bandWidth / (o.bandCb * h) : 0.0;
        if (std::isfinite(sE) && sE > 0.0) {
            mLadrunoSas.bandSElem = sE;
        } else {
            code = RC_START_OTHER;
            if (!mLadrunoSas.warned) {
                mLadrunoSas.warned = true;
                opserr << "WARNING LadrunoSANISAND (SAS-ME) material tag " << this->getTag()
                       << ": -sasBandWidth needs the element's characteristic length and got "
                       << h << " (no element under update, or a degenerate one); update REFUSED."
                          " Drive the point from an element, or give -sasBandScale S." << endln;
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
            double dSe[6], St[6];
            ladrunoSasElastic(S, dStrain, eN, mVoidRatio, dSe);
            for (int i = 0; i < 6; i++) { St[i] = S[i]; St[i] += dSe[i]; }
            const double ft = ladrunoSasF(St, A);
            const double pt = one3 * tr6(St) + m_Presidual;
            if (!(pt < m_Presidual) && ft <= mTolF) {
                elastic = true;
                cp6(S, St);
                for (int i = 0; i < 6; i++)
                    Ee[i] += dStrain[i];
                mLadrunoLastPath = 0;
            } else {
                const double f0 = ladrunoSasF(S, A);
                double nY[6];
                ladrunoSasNormal(S, A, nY);
                // Paper alpha_in rule at the increment START (review numerics
                // item 4) with WP-134's convention exactly: a start is ON the
                // surface unless f0 < -1e-8 of the cone radius sqrt(2/3) m p, and
                // on the surface (alpha - alpha_in):n < 0 re-seats alpha_in before
                // the mode is chosen. A start further inside is interior (the
                // onset, if any, is at the far side, decided there). Measured: this
                // is what separates b16 5496 (re-seat) from b16 5471 / b8 1961
                // (no re-seat) on the ring, all within 4e-8 of the cone radius.
                if (paperRule) {
                    const double pS = one3 * tr6(S) + m_Presidual;
                    if (f0 >= -1.0e-8 * root23 * m_m * pS) {
                        double t0[6];
                        for (int i = 0; i < 6; i++) t0[i] = A[i] - ain[i];
                        const double x0 = ddc6(t0, nY);
                        if (x0 < -ladrunoSasReseatDelta()) {                // Ladruno WP-151: 0 = paper
                            cp6(ain, A);
                            st[LSAS_ALPHA_IN_RESEATS] += 1.0;
                        } else if (x0 < 0.0)
                            st[LSAS_RESEAT_HELD] += 1.0;                     // Ladruno WP-151
                    }
                }
                // U10: the loading test on the TRUE gradient
                double Q[6];
                const double qs = (-one3 * (ddc6(nY, A) + root23 * m_m));
                for (int i = 0; i < 6; i++) {
                    Q[i] = kI1[i] * qs;
                    Q[i] += nY[i];
                }
                const double qn = norm6(Q), nd = norm6(dSe);
                if (f0 < -mTolF) {
                    a = ladrunoSasIntersect(S, A, dStrain, eN, 0.0, 1.0);
                    if (a < 0.0) { st[LSAS_INTERSECT_FAIL] += 1.0; a = 0.0; }
                    mLadrunoLastPath = 2;
                    onset = true;
                } else if (ddc6(Q, dSe) / ((qn == 0 ? 1.0 : qn) * (nd == 0 ? 1.0 : nd))
                           > (-sqrt(mTolF))) {
                    a = 0.0;
                    mLadrunoLastPath = 3;
                } else {
                    // unload then reload: find the first sampled point inside,
                    // then the first one back outside, and Pegasus between them
                    a = -1.0;
                    const int ns = 64;
                    double dE[6], St2[6], dS2[6];
                    double prevA = -1.0;
                    bool inside = false;
                    for (int k = 1; k <= ns; k++) {
                        const double ak = (double)k / ns;
                        for (int i = 0; i < 6; i++) dE[i] = dStrain[i] * ak;
                        ladrunoSasElastic(S, dE, eN, eN, dS2);
                        for (int i = 0; i < 6; i++) { St2[i] = S[i]; St2[i] += dS2[i]; }
                        const double fk = ladrunoSasF(St2, A);
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
                    double dEa[6], dSa[6];
                    for (int i = 0; i < 6; i++) dEa[i] = dStrain[i] * a;
                    const double ea = eN - (1 + m_e_init) * tr6(dEa);
                    ladrunoSasElastic(S, dEa, eN, ea, dSa);
                    for (int i = 0; i < 6; i++) {
                        S[i] += dSa[i];
                        Ee[i] += dEa[i];
                    }
                }
            }
        }
        if (elastic) {
            st[LSAS_ELASTIC] += 1.0;
            for (int i = 0; i < 6; i++) {
                mSigma(i) = S[i]; mAlpha(i) = A[i]; mFabric(i) = Z[i]; mEpsilonE(i) = Ee[i];
                mAlpha_in(i) = ain[i];
            }
            mDGamma = 0.0;
            ladrunoSasModuli(S, mVoidRatio, K, G);
            mCe = GetStiffness(K, G); mCep = mCe; mCep_Consistent = mCe;
            st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
            st[LSAS_LAST_F] = ladrunoSasF(S, A);
            return;
        }
        // ---- 2. the plastic portion ------------------------------------------
        double curP[6];
        for (int i = 0; i < 6; i++) {
            curP[i] = dStrain[i] * a;
            curP[i] += CurStrain[i];
        }
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
        // re-review MED-1: E1 takes the non-compressing gate too -- a compressing
        // increment (e.g. gravity onto a zero-stress start, a deck error) must refuse,
        // not be swallowed by a separation that re-contacts only after a further g_c.
        // "compressing" with a tolerance scaled to the increment: a trace built from
        // B u is round-off, not zero, under isochoric shear -- an exact `> 0.0` split one
        // homogeneous stdBrick's 8 GPs into 6 separating and 2 refusing (re-review test).
        const bool compressing = (tr6(dStrain) > 1.0e-10 * norm6(dStrain));
        const bool e1 = e1q && (p0c <= o.tcP0Max) && !compressing;
        const bool e2 = e2q && !compressing;
        if (e1q && !(p0c <= o.tcP0Max)) st[LSAS_SEP_HELD_HIGHP] += 1.0;
        else if ((e1q && !e1) || (e2q && !e2)) st[LSAS_SEP_HELD_COMPRESSING] += 1.0;
        if (e1 || e2) {
            mLadrunoSas.sepEvent = e1 ? 1 : 2;
            mLadrunoSas.sepCode = code;   // the masked refusal (review #3), recorded at commit
            mLadrunoSas.sepP0 = p0c;
            mSubstepCapHitInME = false;   // a code 9 set it; the cap was reached by a separating point
            mLadrunoSas.sep = true;
            mLadrunoSas.sepTr = tr6(NextStrain);
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
        mLadrunoSas.band = mLadrunoSas.band_n;       // Ladruno band scaling: nor a latch, nor an offset
        mLadrunoSas.bandS = mLadrunoSas.bandS_n;
        mLadrunoSas.bandDe = mLadrunoSas.bandDe_n;
        GetElasticModuli(mSigma_n, eN, K, G);
        mCe = GetStiffness(K, G); mCep = mCe; mCep_Consistent = mCe;
        // "last" columns describe no valid end state after a refusal
        st[LSAS_LAST_RATIO_B] = std::numeric_limits<double>::quiet_NaN();
        st[LSAS_LAST_F] = std::numeric_limits<double>::quiet_NaN();
        // PER-INSTANCE warn-once (review #871 item 2: no process-wide state,
        // nothing to reset on wipe); every refusal is counted in sasStats.
        if (!mLadrunoSas.warned && mLadrunoSas.fdMode == 0) {
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
    for (int i = 0; i < 6; i++) {
        mSigma(i) = S[i]; mAlpha(i) = A[i]; mFabric(i) = Z[i]; mEpsilonE(i) = Ee[i];
        mAlpha_in(i) = ain[i];
    }
    mDGamma = lamSum;
    if (mLadrunoSas.band)   // Ladruno band scaling: the void ratio of the band at the end of this update
        mVoidRatio = m_e_init - (1 + m_e_init) * tr6(NextStrain) + mLadrunoSas.bandDe;
    // Ladruno sasme-speed: a reversal (an alpha_in re-seat) in this update starts
    // a new loading process; its next update starts at dT = 1 again.
    if (st[LSAS_ALPHA_IN_RESEATS] != reseats0)
        mLadrunoSas.carryH = 0.0;
    if (mLadrunoSas.fdMode == 2)
        return;   // Ladruno sasme-speed: a perturbed replay needs the stress only
    {
        OPS_PROFILE_SCOPE("sanisand.sasME.tangent");
        ladrunoSasModuli(S, mVoidRatio, K, G);
        mCe = GetStiffness(K, G);
        if (lastPlastic)
            ladrunoSasContinuumTangent(mSigma, mAlpha, mFabric, mAlpha_in, mVoidRatio, mCep,
                                       mLadrunoSas.band ? mLadrunoSas.bandS : 0.0);   // Ladruno band scaling
        else
            mCep = mCe;
        mCep_Consistent = mCep;
    }
    st[LSAS_LAST_RATIO_B] = ladrunoSasAlphaRatio(A, S, mVoidRatio);
    st[LSAS_LAST_F] = ladrunoSasF(S, A);
    if ((int)mTangType == 3 && lastPlastic && mLadrunoSas.fdMode == 0)
        ladrunoSasFDTangent();   // Ladruno sasme-speed
}

// Ladruno sasme-speed (TanType 3, a PROTOTYPE of the consistent tangent): the
// derivative of this update's map eps -> sigma with its substep sequence FROZEN
// (the accepted dT of the update just made, replayed without error control), by
// forward differences, one replay per strain column. With the controller frozen
// the map is smooth, and its derivative is what the analytic chain over the
// accepted substeps gives (TANGENT_sasme.md of the TIMs speed study). Plane
// strain (eps_zz, gamma_yz, gamma_zx zero in the committed and the trial
// strain): columns 0, 1, 3 only; the rest keep the continuum tangent. A column
// whose replay takes another path (elastic, a rejection, a refusal) keeps the
// continuum column (counted in fdFallbacks). Every trial quantity the replays
// touch is restored. Cost: 3 (plane strain) or 6 replays of the accepted
// substeps per plastic update, on top of the update itself.
void
ManzariDafalias::ladrunoSasFDTangent(void)
{
    OPS_PROFILE_SCOPE("sanisand.sasME.fdTangent");
    double s0[6], a0[6], z0[6], ee0[6], ain0[6], eps0[6];
    for (int i = 0; i < 6; i++) {
        s0[i] = mSigma(i); a0[i] = mAlpha(i); z0[i] = mFabric(i); ee0[i] = mEpsilonE(i);
        ain0[i] = mAlpha_in(i); eps0[i] = mEpsilon(i);
    }
    double stats0[LSAS_COUNT];
    for (int i = 0; i < LSAS_COUNT; i++)
        stats0[i] = mLadrunoSas.stats[i];
    const double dg0 = mDGamma, vr0 = mVoidRatio, carry0 = mLadrunoSas.carryH;
    const bool ref0 = mLadrunoSas.refused, warned0 = mLadrunoSas.warned, cap0 = mSubstepCapHitInME;
    const bool sep0 = mLadrunoSas.sep;
    const double sepTr0 = mLadrunoSas.sepTr, sepP00 = mLadrunoSas.sepP0;
    const int sepEvent0 = mLadrunoSas.sepEvent, sepCode0 = mLadrunoSas.sepCode;
    const int path0 = mLadrunoLastPath, subs0 = mSubstepsTakenInME;
    // Ladruno band scaling x TanType 3: the band's trial state and census are
    // trial quantities too; a replay must end latched as the update did, at the
    // same accepted substep (a latch elsewhere is a jump in the map, not its
    // derivative), or the column keeps the continuum one.
    const bool band0 = mLadrunoSas.band;
    const double bandS0 = mLadrunoSas.bandS, bandDe0 = mLadrunoSas.bandDe;
    const double bandStages0 = mLadrunoSas.bandStages, bandCapped0 = mLadrunoSas.bandCapped;
    const double bandLatches0 = mLadrunoSas.bandLatches, bandSElem0 = mLadrunoSas.bandSElem;
    const int bandLatchAt0 = mLadrunoSas.bandLatchAt;
    const double er0 = mLadrunoLastElasticRatio;
    LadrunoMESubstepTrace* trace0 = mLadrunoTrace.p;
    mLadrunoTrace.p = 0;
    Matrix Ce0(mCe), Cep0(mCep);
    mCep_Consistent = mCep;

    const bool ps = (eps0[2] == 0.0 && eps0[4] == 0.0 && eps0[5] == 0.0
                     && mEpsilon_n(2) == 0.0 && mEpsilon_n(4) == 0.0 && mEpsilon_n(5) == 0.0);
    const int colsPS[3] = {0, 1, 3}, colsAll[6] = {0, 1, 2, 3, 4, 5};
    const int* cols = ps ? colsPS : colsAll;
    const int nc = ps ? 3 : 6;
    const double h = 1.0e-8;
    mLadrunoSas.fdMode = 2;
    for (int k = 0; k < nc; k++) {
        const int c = cols[k];
        mEpsilon(c) = eps0[c] + h;
        mSubstepsTakenInME = 0;
        ladrunoSasIntegrate();
        bool ok = !mLadrunoSas.refused && mLadrunoLastPath == path0
                  && mLadrunoSas.band == band0 && mLadrunoSas.bandLatchAt == bandLatchAt0;   // band scaling
        for (int r = 0; r < 6 && ok; r++)
            ok = std::isfinite(mSigma(r));
        if (ok) {
            const double dc = mEpsilon(c) - eps0[c];   // the perturbation as represented
            for (int r = 0; r < 6; r++)
                mCep_Consistent(r, c) = (mSigma(r) - s0[r]) / dc;
        } else {
            mLadrunoSas.fdFallbacks += 1.0;
        }
        mEpsilon(c) = eps0[c];
    }
    mLadrunoSas.fdMode = 0;

    for (int i = 0; i < 6; i++) {
        mSigma(i) = s0[i]; mAlpha(i) = a0[i]; mFabric(i) = z0[i]; mEpsilonE(i) = ee0[i];
        mAlpha_in(i) = ain0[i];
    }
    for (int i = 0; i < LSAS_COUNT; i++)
        mLadrunoSas.stats[i] = stats0[i];
    mDGamma = dg0; mVoidRatio = vr0; mLadrunoSas.carryH = carry0;
    mLadrunoSas.refused = ref0; mLadrunoSas.warned = warned0; mSubstepCapHitInME = cap0;
    mLadrunoSas.sep = sep0; mLadrunoSas.sepTr = sepTr0; mLadrunoSas.sepEvent = sepEvent0;
    mLadrunoSas.sepCode = sepCode0; mLadrunoSas.sepP0 = sepP00;
    mLadrunoLastPath = path0; mSubstepsTakenInME = subs0; mLadrunoLastElasticRatio = er0;
    mLadrunoSas.band = band0; mLadrunoSas.bandS = bandS0; mLadrunoSas.bandDe = bandDe0;   // band scaling
    mLadrunoSas.bandStages = bandStages0; mLadrunoSas.bandCapped = bandCapped0;
    mLadrunoSas.bandLatches = bandLatches0; mLadrunoSas.bandSElem = bandSElem0;
    mLadrunoSas.bandLatchAt = bandLatchAt0;
    mLadrunoTrace.p = trace0;
    mCe = Ce0; mCep = Cep0;
}
