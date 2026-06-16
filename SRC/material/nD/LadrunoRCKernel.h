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

#ifndef LadrunoRCKernel_h
#define LadrunoRCKernel_h

// ============================================================================
//  LadrunoRCKernel.h  --  header-only, OpenSees-free RC plastic-damage kernel
//
//  Phase-1 of the Ladruno nonlinear RC shell stack (ADR 19_ladruno_rc_shell_adr.md).
//  Faithful clone of ASDConcrete3DMaterial::compute() (SRC/material/nD/
//  ASDConcrete3DMaterial.cpp:2322-2523) -- spectral effective-stress split,
//  dual dt/dc scalar damage, Lubliner-Lee-Fenves biaxial envelope, crack-band
//  hardening backbones -- PLUS the one Phase-1 physics addition (Derivation B,
//  gate-verified): MCFT compression softening applied to the STRENGTH axis:
//
//      sigma = (1 - dt_bar) ST  +  beta * (1 - dc_bar) SC
//      beta  = clamp( 1/(0.8 + 170 eps_1), betaFloor, 1 )
//      eps_1 = in-plane principal TENSILE strain (transverse to the strut)
//
//  beta multiplies ONLY the assembled compressive cone -- never the strain
//  abscissa, the equivalent-strain measure, the effective stress, dc_plastic,
//  or dc_bar.  With betaOn=false the kernel is byte-faithful to ASDConcrete3D.
//
//  Conventions (match ASDConcrete3D + the OpenSees 3D nDMaterial interface):
//    - Voigt order {00,11,22,01,12,02} (indices 0..5).
//    - strain uses ENGINEERING shear (gamma = 2*eps_tensor); C0 shear rows = mu.
//    - stress Voigt is true stress (no factors).
//
//  Pattern cloned from LadrunoJ2Kernel.h ("one core, many views").
//  Verified against tests/_testbed/rc_shell_ref.py (A1 gate PASS).
// ============================================================================

#include <math.h>

namespace ladruno_rc_kernel {

enum { STATUS_OK = 0, STATUS_SINGULAR = 1, STATUS_NO_CONVERGE = 2 };

const int    MAXPTS    = 32;       // hardening backbone point cap
const double YTOL      = 1.0e-12;  // HardeningLaw y/q floor (m_ytolerance)
const double XTOL      = 1.0e-12;  // HardeningLaw x tolerance (m_xtolerance)
const double STRESS_TOL = 1.0e-9;  // tensile-measure gate (ht.stressTolerance)

// ---------------------------------------------------------------------------
// Hardening backbone: (x = equiv strain, y = nominal stress, q = effective stress)
// d = 1 - y/q.  Clone of HardeningLaw::evaluateAt (cpp:1010-1054).
// ---------------------------------------------------------------------------
struct Backbone {
    int    n;
    double x[MAXPTS], y[MAXPTS], q[MAXPTS];
};

struct HPoint { double x, y, d, q; };

inline HPoint evaluateAt(const Backbone& b, double x)
{
    double x1, x2, y1, y2, q1, q2;
    bool found = false;
    for (int i = 1; i < b.n; ++i) {
        if (x <= b.x[i] + XTOL) {
            x1 = b.x[i - 1]; x2 = b.x[i];
            y1 = b.y[i - 1]; y2 = b.y[i];
            q1 = b.q[i - 1]; q2 = b.q[i];
            found = true;
            break;
        }
    }
    if (!found) {
        int last = b.n - 1;
        x1 = b.x[last];
        x2 = x;
        double span = x - x1;
        y1 = b.y[last];
        double tan_y = (b.y[last] - b.y[last - 1]) / (b.x[last] - b.x[last - 1]);
        y2 = tan_y > 0.0 ? y1 + span * tan_y : y1;
        q1 = b.q[last];
        double tan_q = (b.q[last] - b.q[last - 1]) / (b.x[last] - b.x[last - 1]);
        q2 = tan_q > 0.0 ? q1 + span * tan_q : q1;
    }
    double xspan  = x2 - x1;
    double xratio = xspan > 0.0 ? (x - x1) / xspan : 0.0;
    double y = y1 + (y2 - y1) * xratio; if (y < YTOL) y = YTOL;
    double q = q1 + (q2 - q1) * xratio; if (q < YTOL) q = YTOL;
    HPoint p; p.x = x; p.y = y; p.q = q; p.d = 1.0 - y / q;
    return p;
}

inline double plasticStrain(const HPoint& p, double E) { return p.x - p.q / E; }

inline double backboneMaxStress(const Backbone& b)
{
    double m = 0.0;
    for (int i = 0; i < b.n; ++i) if (b.y[i] > m) m = b.y[i];
    return m;
}

// Build the (x,y,q) backbone from user (strain, nominal-stress, damage) EXACTLY as
// ASDConcrete3D::HardeningLaw c-tor + adjust() (cpp:869-939, 1134-1180): prepend (0,0,0)
// if needed, force d[0]=0, force the first segment elastic (y[1]=E*x[1]), cap every
// secant slope at E, enforce monotone plastic strain + non-decreasing damage, then
// q = y/(1-d) on the ADJUSTED points. This is what makes the elastic branch E-consistent
// (dc_plastic=0 until yield) and is required for byte-faithful reduce-to-ASDConcrete3D.
inline void buildBackbone(Backbone& b, double E,
                          const double* xin, const double* yin, const double* din, int nin)
{
    double x[MAXPTS], y[MAXPTS], d[MAXPTS];
    int n = nin > MAXPTS ? MAXPTS : nin;
    double xmax = 0.0, ymax = 0.0;
    for (int i = 0; i < n; ++i) {
        x[i] = fabs(xin[i]); y[i] = fabs(yin[i]);
        d[i] = din[i] < 0.0 ? 0.0 : (din[i] > 1.0 ? 1.0 : din[i]);
        if (x[i] > xmax) xmax = x[i];
        if (y[i] > ymax) ymax = y[i];
    }
    // first-point handling
    if (x[0] > 0.0) {
        if (y[0] == 0.0) {
            y[0] = E * x[0];
        } else if (n < MAXPTS) {                 // insert (0,0,0) at front
            for (int i = n; i >= 1; --i) { x[i] = x[i-1]; y[i] = y[i-1]; d[i] = d[i-1]; }
            x[0] = 0.0; y[0] = 0.0; d[0] = 0.0; ++n;
        }
    } else {
        y[0] = 0.0;
    }
    d[0] = 0.0;
    if (n > 1) y[1] = E * x[1];                  // force elastic first segment
    // tolerances
    double dxmin = xmax, dymin = ymax;
    for (int i = 1; i < n; ++i) {
        double dx = fabs(x[i] - x[i-1]);
        double dy = fabs(y[i] - y[i-1]);
        if (dx > 0.0 && dx < dxmin) dxmin = dx;
        if (dy > 0.0 && dy < dymin) dymin = dy;
    }
    double xtol = 1.0e-6 * dxmin, ytol = 1.0e-6 * dymin;
    // adjust() — compute q on corrected points
    double Eloc = (n > 1 && x[1] > 0.0) ? y[1] / x[1] : E;
    b.x[0] = x[0]; b.y[0] = y[0];
    b.q[0] = (d[0] < 1.0) ? y[0] / (1.0 - d[0]) : y[0];
    for (int i = 1; i < n; ++i) {
        double xi = x[i], xold = x[i-1], yi = y[i], yold = y[i-1], di = d[i], dold = d[i-1];
        if (xi <= xold) xi += xtol;
        if (yi < ytol) yi = ytol;
        double Ei = (yi - yold) / (xi - xold);
        if (Ei > Eloc) yi = yold + (xi - xold) * Eloc;
        double eepd_old = dold < 1.0 ? xold - yold / ((1.0 - dold) * Eloc) : xold;
        double eepd     = di   < 1.0 ? xi   - yi   / ((1.0 - di)   * Eloc) : -1.0;
        if (eepd < eepd_old) {
            double v = 1.0 - yi / Eloc / (xi - eepd_old);
            di = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
        }
        if (di < dold) di = dold;
        Ei = (yi - yold) / (xi - xold);
        double Ed = (1.0 - di) * Eloc;
        if (Ei > Ed) yi = Ed * (xi - eepd_old);
        x[i] = xi; y[i] = yi; d[i] = di;
        b.x[i] = xi; b.y[i] = yi; b.q[i] = yi / (1.0 - di);
    }
    b.n = n;
}

inline double macauley(double x) { return x > 0.0 ? x : 0.0; }

// ---------------------------------------------------------------------------
// Lubliner-Lee-Fenves criterion (clone of cpp:2483-2497).
// Local shape parameter renamed beta_lub to reserve `beta` for MCFT softening.
// ---------------------------------------------------------------------------
inline double lublinerCriterion(double s1, double s2, double s3,
                                double ft, double fc, double k1, double scale, double Kc)
{
    double fb    = 1.16 * fc;
    double gamma = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0);
    double alpha = (fb - fc) / (2.0 * fb - fc);
    double I1    = s1 + s2 + s3;
    double p     = I1 / 3.0;
    double d1 = s1 - p, d2 = s2 - p, d3 = s3 - p;
    double J2 = 0.5 * (d1 * d1 + d2 * d2 + d3 * d3);
    double beta_lub = fc / ft * (1.0 - alpha) - (1.0 + alpha);
    double smax = macauley(s1);
    double smin = macauley(-s1);
    return (1.0 / (1.0 - alpha) *
            (alpha * I1 + sqrt(3.0 * J2) + k1 * beta_lub * smax - gamma * smin)) * scale;
}

// ---------------------------------------------------------------------------
// Symmetric 3x3 eigen-decomposition (cyclic Jacobi). Eigenvalues sorted
// DESCENDING (so Si[0] is the largest principal, matching the spine's smax use).
// V[r][c] = component r of eigenvector c.
// ---------------------------------------------------------------------------
inline void eig3x3sym(const double A[3][3], double Si[3], double V[3][3])
{
    double a[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a[i][j] = A[i][j];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) V[i][j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
        if (off < 1.0e-20) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (fabs(a[p][q]) < 1.0e-300) continue;
                double app = a[p][p], aqq = a[q][q], apq = a[p][q];
                double phi = 0.5 * atan2(2.0 * apq, aqq - app);
                double c = cos(phi), s = sin(phi);
                // rotate A
                for (int k = 0; k < 3; ++k) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                // accumulate eigenvectors
                for (int k = 0; k < 3; ++k) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    Si[0] = a[0][0]; Si[1] = a[1][1]; Si[2] = a[2][2];
    // sort descending (selection sort on 3, carry eigenvectors)
    for (int i = 0; i < 2; ++i) {
        int m = i;
        for (int j = i + 1; j < 3; ++j) if (Si[j] > Si[m]) m = j;
        if (m != i) {
            double t = Si[i]; Si[i] = Si[m]; Si[m] = t;
            for (int k = 0; k < 3; ++k) { double tv = V[k][i]; V[k][i] = V[k][m]; V[k][m] = tv; }
        }
    }
}

// 6x6 eigenprojector pjj for eigenvector column j  (clone of computePjj,
// cpp:255-300, verified compact form: pjj[a][b] = E[a]*(b<3?1:2)*E[b]).
inline void computePjj(const double V[3][3], int j, double pjj[6][6])
{
    double vx = V[0][j], vy = V[1][j], vz = V[2][j];
    double E[6] = { vx*vx, vy*vy, vz*vz, vx*vy, vy*vz, vx*vz };
    for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b)
            pjj[a][b] = E[a] * (b < 3 ? 1.0 : 2.0) * E[b];
}

inline double heavyside(double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? 0.0 : 0.5); }

// ---------------------------------------------------------------------------
// Spectral stress decomposition (clone of StressDecomposition::compute,
// cpp:781-851). Builds PT, PC (6x6), ST, SC (6, = PT:S, PC:S), R.
// ---------------------------------------------------------------------------
struct Decomp {
    double Si[3], V[3][3];
    double PT[6][6], PC[6][6];
    double ST[6], SC[6];
    double R;
};

inline void decompose(const double S[6], double cdf, Decomp& D)
{
    double M[3][3] = {
        { S[0], S[3], S[5] },
        { S[3], S[1], S[4] },
        { S[5], S[4], S[2] }
    };
    eig3x3sym(M, D.Si, D.V);

    for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) { D.PT[a][b] = 0.0; D.PC[a][b] = 0.0; }
    double pjj[6][6];
    for (int j = 0; j < 3; ++j) {
        computePjj(D.V, j, pjj);
        double hp = heavyside(D.Si[j]);
        double hn = heavyside(-D.Si[j]);
        if (hp > 0.0) for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) D.PT[a][b] += hp * pjj[a][b];
        if (hn > 0.0) for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) D.PC[a][b] += hn * pjj[a][b];
    }
    // PO = I - PT - PC; PT += PO/2; PC += PO/2
    for (int a = 0; a < 6; ++a) {
        for (int b = 0; b < 6; ++b) {
            double po = (a == b ? 1.0 : 0.0) - D.PT[a][b] - D.PC[a][b];
            D.PT[a][b] += 0.5 * po;
            D.PC[a][b] += 0.5 * po;
        }
    }
    for (int a = 0; a < 6; ++a) {
        double st = 0.0, sc = 0.0;
        for (int b = 0; b < 6; ++b) { st += D.PT[a][b] * S[b]; sc += D.PC[a][b] * S[b]; }
        D.ST[a] = st; D.SC[a] = sc;
    }
    // R factor
    double Rnum = 0.0, Rden = 0.0;
    for (int j = 0; j < 3; ++j) {
        double Sj = D.Si[j];
        if (Sj > 0.0) Sj *= cdf;
        if (Sj > 0.0) Rnum += Sj;
        Rden += fabs(Sj);
    }
    D.R = Rden > 0.0 ? Rnum / Rden : 0.0;
}

// ---------------------------------------------------------------------------
// MCFT compression-softening factor (applied to the STRENGTH axis only).
// ---------------------------------------------------------------------------
inline double betaCompr(double e1, double floorv)
{
    double b = 1.0 / (0.8 + 170.0 * e1);
    if (b > 1.0) b = 1.0;
    if (b < floorv) b = floorv;
    return b;
}
inline double dBetaCompr(double e1, double floorv)   // 0 when clamped
{
    double b = 1.0 / (0.8 + 170.0 * e1);
    if (b >= 1.0 || b <= floorv) return 0.0;
    return -170.0 * b * b;
}

// Largest eigenvalue of the in-plane 2x2 TOTAL strain tensor (clamped >=0).
// gxy is ENGINEERING shear (tensor = gxy/2). p1[2] = in-plane unit eigenvector;
// degen<0 flags an equibiaxial/degenerate state (caller blends the tangent term).
inline double eps1FromMembrane(double exx, double eyy, double gxy, double p1[2], int* degen)
{
    double c = 0.5 * (exx + eyy);
    double half_diff = 0.5 * (exx - eyy);
    double half_gxy  = 0.5 * gxy;
    double r = sqrt(half_diff * half_diff + half_gxy * half_gxy);
    double scale = fabs(exx) + fabs(eyy) + fabs(gxy) + 1.0e-30;
    if (r < 1.0e-9 * scale) {            // degenerate / equibiaxial
        if (degen) *degen = 1;
        p1[0] = 1.0; p1[1] = 0.0;        // arbitrary; tangent term blended by caller
        return macauley(c);
    }
    if (degen) *degen = 0;
    double lam1 = c + r;
    // eigenvector of [[exx, gxy/2],[gxy/2, eyy]] for lam1
    double vx = half_gxy;
    double vy = lam1 - exx;
    double nrm = sqrt(vx * vx + vy * vy);
    if (nrm < 1.0e-300) { p1[0] = 1.0; p1[1] = 0.0; }
    else { p1[0] = vx / nrm; p1[1] = vy / nrm; }
    return macauley(lam1);
}

// isotropic elastic 6x6 (clone of getInitialTangent, cpp:1733-1743). Engineering shear.
inline void elasticTangent(double E, double nu, double C0[6][6])
{
    for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) C0[a][b] = 0.0;
    double mu2 = E / (1.0 + nu);
    double lam = nu * mu2 / (1.0 - 2.0 * nu);
    double mu  = 0.5 * mu2;
    mu2 += lam;
    C0[0][0] = C0[1][1] = C0[2][2] = mu2;
    C0[0][1] = C0[1][0] = C0[0][2] = C0[2][0] = C0[1][2] = C0[2][1] = lam;
    C0[3][3] = C0[4][4] = C0[5][5] = mu;
}

// ---------------------------------------------------------------------------
// Parameters + committed history
// ---------------------------------------------------------------------------
struct Params {
    double E, nu, Kc;
    double fcft_ratio;                 // = max(5, fcmax/ftmax), set by setupParams
    double betaFloor;                  // MCFT clamp lower bound
    double cdf;                        // cross-damage factor (default 0)
    double eta;                        // viscous regularization (default 0)
    bool   betaOn;                     // false => beta=1 (identity gate)
    bool   lublinerTCReduced;          // true  => k1=0 in compressive measure
    int    tangentMode;                // 0 algorithmic, 1 secant, 2 numerical
    // --- Phase 2a: fixed-crack aggregate-interlock shear retention (membrane) ---
    bool   interlockOn;                // false => no interlock (Phase-1 baseline)
    int    shearRetMode;              // 0 = MCFT v_ci cap (default); reserved 1=dsfm
    double aggSize;                   // max aggregate size a_g (SI: mm)
    double crackStrain;              // eps_cr crack-formation strain (<=0 => ftmax/E)
    double crackSpacing;            // s_theta crack spacing (<=0 => lch, else 1)
    double lch;                     // characteristic length (s_theta fallback)
    double betaSrMin;              // residual shear-retention floor (tangent safety)
    double sqrtFc;                // sqrt(fc') cache, set by setupParams (SI: sqrt(MPa))
    Backbone ht, hc;                   // tension / compression backbones
};

inline void setupParams(Params& P)    // call after ht/hc are filled
{
    double fcmax = backboneMaxStress(P.hc);
    double ftmax = backboneMaxStress(P.ht);
    P.fcft_ratio = (ftmax > 0.0) ? (fcmax / ftmax > 5.0 ? fcmax / ftmax : 5.0) : 1.0;
    P.sqrtFc = fcmax > 0.0 ? sqrt(fcmax) : 0.0;
    if (P.crackStrain <= 0.0) P.crackStrain = (P.E > 0.0) ? ftmax / P.E : 0.0;
}

struct RCHist {
    double stress_eff[6];   // committed effective stress
    double strain[6];       // committed total strain (engineering shear)
    double xt, xc;          // committed tensile/compressive equiv strain (svt/svc)
    double dt_bar, dc_bar;  // last nominal damage scalars (output/diagnostic)
    double beta;            // last MCFT factor (diagnostic)
    double eps1;            // last in-plane principal tensile strain (diagnostic)
    // --- Phase 2a: fixed-crack aggregate-interlock state ---
    double cracked;         // 0/1 fixed-crack flag (double for serialization)
    double crackC, crackS;  // in-plane crack-NORMAL direction (cos,sin), frozen at capture
    double wmax;            // max crack width reached (monotonic interlock degradation)
    double betaSr;          // last realized shear-retention secant (diagnostic)
};

inline void histZero(RCHist& h)
{
    for (int i = 0; i < 6; ++i) { h.stress_eff[i] = 0.0; h.strain[i] = 0.0; }
    h.xt = h.xc = 0.0; h.dt_bar = h.dc_bar = 0.0; h.beta = 1.0; h.eps1 = 0.0;
    h.cracked = 0.0; h.crackC = 1.0; h.crackS = 0.0; h.wmax = 0.0; h.betaSr = 1.0;
}

inline double equivTensile(const double Si[3], double fcft, double Kc)
{
    if (Si[0] < STRESS_TOL) return 0.0;
    return lublinerCriterion(Si[0], Si[1], Si[2], 1.0 / fcft, 1.0, 1.0, 1.0 / fcft, Kc);
}
inline double equivCompressive(const double Si[3], double fcft, double Kc, double k1)
{
    double s0 = Si[0] < 0.0 ? Si[0] : 0.0;
    double s1 = Si[1] < 0.0 ? Si[1] : 0.0;
    double s2 = Si[2] < 0.0 ? Si[2] : 0.0;
    return lublinerCriterion(s0, s1, s2, 1.0 / fcft, 1.0, k1, 1.0, Kc);
}

// ---------------------------------------------------------------------------
// The 3D return/evaluation. Pure function of (P, eps6, committed in) ->
// (sig6, Dtan6, trial out). No mutation of `in`. Phase-1: implex off, eta-rate
// optional. Mirrors compute() (cpp:2322-2481) step for step.
// ---------------------------------------------------------------------------
inline int returnMap3D(const Params& P, const double eps6[6], const RCHist& in,
                       double sig6[6], double Dtan6[6][6], RCHist& out,
                       bool do_tangent = true, double* residual = 0)
{
    const double E = P.E;
    double C0[6][6]; elasticTangent(E, P.nu, C0);

    // elastic predictor: SEFFn = SEFF_commit + C0:(En - En-1)
    double seff[6];
    for (int a = 0; a < 6; ++a) {
        double acc = in.stress_eff[a];
        for (int b = 0; b < 6; ++b) acc += C0[a][b] * (eps6[b] - in.strain[b]);
        seff[a] = acc;
    }

    Decomp D;
    decompose(seff, P.cdf, D);

    // MCFT beta from in-plane membrane principal tensile strain (gxy = eps6[3] engineering)
    double p1[2]; int degen = 0;
    double e1 = eps1FromMembrane(eps6[0], eps6[1], eps6[3], p1, &degen);
    double beta = P.betaOn ? betaCompr(e1, P.betaFloor) : 1.0;
    double k1c  = (P.betaOn && P.lublinerTCReduced) ? 0.0 : 1.0;

    // committed hardening
    HPoint pt = evaluateAt(P.ht, in.xt);
    HPoint pc = evaluateAt(P.hc, in.xc);
    double xt_pl = plasticStrain(pt, E);
    double xc_pl = plasticStrain(pc, E);

    // rate coefficients (eta=0 -> rc1=0, rc2=1)
    double rc1 = (P.eta > 0.0) ? P.eta / (P.eta + 1.0) : 0.0;
    double rc2 = 1.0 - rc1;

    // /E converts the Lubliner STRESS measure to the strain abscissa (cpp:2509/2522)
    double xt_meas = equivTensile(D.Si, P.fcft_ratio, P.Kc) / E;
    double xc_meas = equivCompressive(D.Si, P.fcft_ratio, P.Kc, k1c) / E;
    double xt_trial = xt_meas + xt_pl;
    double xc_trial = xc_meas + xc_pl;

    double new_xt = in.xt, new_xc = in.xc;
    if (xt_trial > pt.x) { double xn = rc1 * pt.x + rc2 * xt_trial; pt = evaluateAt(P.ht, xn); new_xt = pt.x; }
    if (xc_trial > pc.x) { double xn = rc1 * pc.x + rc2 * xc_trial; pc = evaluateAt(P.hc, xn); new_xc = pc.x; }

    double seff_eq_t = (pt.x - xt_pl) * E;
    double dt_plastic = seff_eq_t > 0.0 ? 1.0 - pt.q / seff_eq_t : 0.0;
    double seff_eq_c = (pc.x - xc_pl) * E;
    double dc_plastic = seff_eq_c > 0.0 ? 1.0 - pc.q / seff_eq_c : 0.0;

    // mix damages (cdf=0 -> R=0 -> identity)
    HPoint auxp = evaluateAt(P.hc, pt.x * P.fcft_ratio);
    double aux_pl = plasticStrain(auxp, E);
    double aux_seff = (auxp.x - aux_pl) * E;
    double aux_d = aux_seff > 0.0 ? 1.0 - auxp.q / aux_seff : 0.0;
    dc_plastic = 1.0 - (1.0 - dc_plastic) * (1.0 - D.R * aux_d);
    double pc_d = 1.0 - (1.0 - pc.d) * (1.0 - D.R * auxp.d);

    double dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
    double dc_bar = pc_d + dc_plastic - pc_d * dc_plastic;

    // nominal stress -- THE Phase-1 insertion: beta multiplies the assembled compressive cone
    for (int a = 0; a < 6; ++a)
        sig6[a] = (1.0 - dt_bar) * D.ST[a] + beta * (1.0 - dc_bar) * D.SC[a];

    // ---- Phase 2a: fixed-crack aggregate-interlock shear retention (membrane) ----
    // Replaces the smeared (isotropic-damage) membrane shear tau_xy on a FROZEN crack
    // plane with tau_ci = clamp(G*gamma_nt, +/- v_ci,max), v_ci,max degrading with the
    // (monotonically grown) crack width w = macauley(eps_n)*s_theta (Vecchio-Collins).
    // SI units assumed (fc in MPa = N/mm^2, w & a_g in mm). Default-off => no change.
    // Phase 2a = a v_ci,max BOUND on the already-assembled (damage-reduced) crack-plane
    // shear, NOT a substitution with bare-elastic shear. So below the cap the stress is
    // UNCHANGED (continuous across cracking, no double-count), and the cap clips the
    // smeared crack-shear to the Vecchio-Collins limit. This keeps stress+tangent a
    // matched secant (see the tangent block; m_sigma . m_eps == 1).
    double crk_c = in.crackC, crk_s = in.crackS, cracked = in.cracked, wmax = in.wmax;
    double betaSr = 1.0;
    bool   il_active = false, il_capped = false;
    double il_m[3] = { 0.0, 0.0, 0.0 };             // m_eps: stress back-rotation of a shear change
    if (P.interlockOn) {
        if (cracked < 0.5 && e1 >= P.crackStrain && !degen) {   // freeze crack normal
            cracked = 1.0; crk_c = p1[0]; crk_s = p1[1];
        }
        if (cracked >= 0.5) {
            il_active = true;
            double exx = eps6[0], eyy = eps6[1], gxy = eps6[3];
            double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
            double en = exx*c2 + eyy*s2 + gxy*cs;               // strain normal to the FROZEN crack
            double sth = (P.crackSpacing > 0.0) ? P.crackSpacing : (P.lch > 0.0 ? P.lch : 1.0);
            double w = macauley(en) * sth;
            if (w > wmax) wmax = w;                              // irreversible (monotonic)
            double denom = 0.31 + 24.0 * wmax / (P.aggSize + 16.0);
            double vcimax = (denom > 0.0) ? 0.18 * P.sqrtFc / denom : 0.18 * P.sqrtFc;
            // tau_sm = smeared (damaged) crack-plane shear stress = m_sigma . sig_ip
            double sxx = sig6[0], syy = sig6[1], sxy = sig6[3];
            double tau_sm = cs * (syy - sxx) + (c2 - s2) * sxy;
            double tau_ci = tau_sm;
            if (tau_ci >  vcimax) { tau_ci =  vcimax; il_capped = true; }
            if (tau_ci < -vcimax) { tau_ci = -vcimax; il_capped = true; }
            double dtau = tau_ci - tau_sm;                       // <=0 in magnitude (a clip)
            betaSr = (fabs(tau_sm) > 1.0e-30) ? tau_ci / tau_sm : 1.0;
            // back-rotate the shear change: global stress deviation = dtau * m_eps
            il_m[0] = -2.0*cs; il_m[1] = 2.0*cs; il_m[2] = c2 - s2;
            sig6[0] += dtau * il_m[0];
            sig6[1] += dtau * il_m[1];
            sig6[3] += dtau * il_m[2];
        }
    }

    // effective stress reassembly for next-step commit (beta NOT applied)
    for (int a = 0; a < 6; ++a)
        out.stress_eff[a] = (1.0 - dt_plastic) * D.ST[a] + (1.0 - dc_plastic) * D.SC[a];
    for (int a = 0; a < 6; ++a) out.strain[a] = eps6[a];
    out.xt = new_xt; out.xc = new_xc;
    out.dt_bar = dt_bar; out.dc_bar = dc_bar; out.beta = beta; out.eps1 = e1;
    out.cracked = cracked; out.crackC = crk_c; out.crackS = crk_s;
    out.wmax = wmax; out.betaSr = betaSr;

    // tangent
    if (do_tangent) {
        if (P.tangentMode == 2) {
            // numerical (forward difference) -- reserved diagnostic
            double base[6]; for (int a = 0; a < 6; ++a) base[a] = sig6[a];
            for (int b = 0; b < 6; ++b) {
                double epert[6]; for (int a = 0; a < 6; ++a) epert[a] = eps6[a];
                double h = 1.0e-8 * (fabs(eps6[b]) > 1.0 ? fabs(eps6[b]) : 1.0);
                epert[b] += h;
                double sp[6]; RCHist dump;
                returnMap3D(P, epert, in, sp, Dtan6 /*ignored*/, dump, false, 0);
                for (int a = 0; a < 6; ++a) Dtan6[a][b] = (sp[a] - base[a]) / h;
            }
        } else {
            // W_B = I - dt_bar PT - (1 - beta(1-dc_bar)) PC ;  C_sec = W_B C0
            double wc = 1.0 - beta * (1.0 - dc_bar);
            double WB[6][6];
            for (int a = 0; a < 6; ++a)
                for (int c = 0; c < 6; ++c)
                    WB[a][c] = (a == c ? 1.0 : 0.0) - dt_bar * D.PT[a][c] - wc * D.PC[a][c];
            for (int a = 0; a < 6; ++a) {
                for (int c = 0; c < 6; ++c) {
                    double acc = 0.0;
                    for (int k = 0; k < 6; ++k) acc += WB[a][k] * C0[k][c];
                    Dtan6[a][c] = acc;
                }
            }
            // beta cross-term (algorithmic only, tangentMode 0):
            //   d sigma/d eps += (1-dc_bar) SC (x) [dBeta/deps1 * deps1/deps]
            //   deps1/deps (engineering dual) = [p1x^2, p1y^2, 0, p1x p1y, 0, 0]
            if (P.tangentMode == 0 && P.betaOn) {
                double dbeta = dBetaCompr(e1, P.betaFloor);
                if (dbeta != 0.0) {
                    double P1[6] = { 0,0,0,0,0,0 };
                    if (degen) { P1[0] = 0.5; P1[1] = 0.5; }   // equibiaxial blend
                    else { P1[0] = p1[0]*p1[0]; P1[1] = p1[1]*p1[1]; P1[3] = p1[0]*p1[1]; }
                    double coef = (1.0 - dc_bar) * dbeta;
                    for (int a = 0; a < 6; ++a)
                        for (int c = 0; c < 6; ++c)
                            Dtan6[a][c] += coef * D.SC[a] * P1[c];
                }
            }
            // Phase-2a interlock CONSISTENT tangent. The stress clips the crack-shear
            // tau_sm = m_sigma.sig_ip to +/- v_ci,max via sig_ip += dtau*m_eps. Below the
            // cap dtau==0 => stress unchanged => tangent unchanged (no cross-term). When
            // capped, tau_ci is pinned, so d(tau_ci)/deps=0 and the consistent secant is
            //   Dtan_ip -= (1-betaSrMin) * m_eps (x) (m_sigma . Dtan_ip).
            // Since m_sigma . m_eps == (c^2+s^2)^2 == 1, this removes the crack-shear
            // stiffness exactly (the full removal at betaSrMin=0; betaSrMin keeps a small
            // residual for conditioning — a documented stabilization, stress stays fully
            // capped regardless). NOTE: the dependence of v_ci,max on the crack-opening
            // strain (d v_ci,max/d eps_n, active only while w is at its monotone max) is a
            // 2nd-order normal->shear coupling that is intentionally OMITTED here, matching
            // the secant character of the Phase-1 W_B tangent (which is itself W_B*C0, not
            // the fully-consistent damage tangent). Use tangentMode 2 (numerical) if an
            // exact cracked tangent is ever needed.
            if (il_capped) {
                double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
                const int vidx[3] = { 0, 1, 3 };
                double msig[3] = { -cs, cs, c2 - s2 };       // m_sigma (tensor shear projector)
                double row[3];
                for (int j = 0; j < 3; ++j) {
                    double acc = 0.0;
                    for (int k = 0; k < 3; ++k) acc += msig[k] * Dtan6[vidx[k]][vidx[j]];
                    row[j] = acc;
                }
                double fac = 1.0 - P.betaSrMin;
                for (int a = 0; a < 3; ++a)
                    for (int c3 = 0; c3 < 3; ++c3)
                        Dtan6[vidx[a]][vidx[c3]] -= fac * il_m[a] * row[c3];
            }
        }
    }

    if (residual) {
        double sn = 0.0; for (int a = 0; a < 6; ++a) sn += sig6[a] * sig6[a];
        *residual = sqrt(sn);
    }
    return STATUS_OK;
}

} // namespace ladruno_rc_kernel

#endif // LadrunoRCKernel_h
