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

// ---------------------------------------------------------------------------
// Phase 3b: crack-band (Bazant-Oh) fracture-energy regularization. Faithful clone of
// ASDConcrete3D HardeningLaw::computeFractureEnergy + regularize. fractureEnergy returns the
// specific (per-length) area under the softening branch; regularize() rescales the post-peak
// strain abscissa in place so g_reg = G_f0*(lch_ref/lch) => the dissipated energy density * lch
// (the physical fracture energy) is mesh-objective. Self-contained (regularize recomputes g0),
// so the Backbone struct + its serialization are unchanged. Default off => never called.
// ---------------------------------------------------------------------------
inline double fractureEnergy(const Backbone& b, double E, int* bounded, int* pos1o, int* pos2o)
{
    *bounded = 0; *pos1o = 0; *pos2o = 0;
    const int n = b.n;
    int pos1 = 0; bool f1 = false;
    for (int i = 1; i < n; ++i)
        if ((b.y[i] - b.y[i-1]) / (b.x[i] - b.x[i-1]) < 0.0) { pos1 = i-1; f1 = true; break; }
    if (!f1) return 0.0;                                   // no softening => unbounded energy
    int pos2 = 0; bool f2 = false;
    for (int i = pos1+1; i < n; ++i)
        if ((b.y[i] - b.y[i-1]) / (b.x[i] - b.x[i-1]) >= 0.0) { pos2 = i-1; f2 = true; break; }
    double g_add = 0.0;
    if (!f2) {                                             // extend the last segment to zero stress
        if (b.y[n-1] > 0.0) {
            double k = (b.y[n-1] - b.y[n-2]) / (b.x[n-1] - b.x[n-2]);
            double x3 = b.x[n-1] - b.y[n-1] / k;
            g_add = b.y[n-1] * (x3 - b.x[n-1]) / 2.0;
        }
        pos2 = n-1;
    }
    double d1 = b.q[pos1] > 0.0 ? 1.0 - b.y[pos1] / b.q[pos1] : 0.0;
    double Ed = (1.0 - d1) * E;
    double g = Ed > 0.0 ? b.y[pos1] * b.y[pos1] / Ed / 2.0 : 0.0;   // unloading triangle at the peak
    for (int i = pos1+1; i <= pos2; ++i) g += (b.x[i] - b.x[i-1]) * (b.y[i-1] + b.y[i]) / 2.0;
    g += g_add;
    *bounded = 1; *pos1o = pos1; *pos2o = pos2;
    return g;
}

inline void regularize(Backbone& b, double lch, double lch_ref, double E)
{
    int bounded, pos1, pos2;
    double g0 = fractureEnergy(b, E, &bounded, &pos1, &pos2);
    double lch_scale = lch > 0.0 ? lch_ref / lch : 0.0;
    if (!bounded || lch_scale <= 0.0 || lch_scale == 1.0) return;
    double gnew = g0 * lch_scale;
    double gmin = (b.y[pos1] * b.x[pos1] / 2.0) * 1.01;    // floor (lch too large): keep monotone x
    if (gnew < gmin) gnew = gmin;
    double tol = 1.0e-3 * gnew;
    double x0 = b.x[pos1];
    double g = g0, dscale = (g0 > 0.0) ? gnew / g0 : 1.0;
    for (int iter = 0; iter < 10; ++iter) {
        for (int i = pos1+1; i < b.n; ++i) {
            double xi_inel = b.x[i] - b.y[i] / E; if (xi_inel < 0.0) xi_inel = 0.0;
            double xi_pl = b.x[i] - b.q[i] / E;
            double xi_ratio = xi_inel > 0.0 ? xi_pl / xi_inel : 0.0;
            b.x[i] = x0 + (b.x[i] - x0) * dscale;          // stretch the post-peak strain
            xi_inel = b.x[i] - b.y[i] / E; if (xi_inel < 0.0) xi_inel = 0.0;
            xi_pl = xi_inel * xi_ratio;                    // keep the plastic-to-inelastic ratio
            b.q[i] = E * (b.x[i] - xi_pl);
        }
        int bd, p1, p2;
        g = fractureEnergy(b, E, &bd, &p1, &p2);
        if (fabs(g - gnew) < tol) break;
        dscale = (g > 0.0) ? gnew / g : 1.0;
    }
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

// ---------------------------------------------------------------------------
// Phase 3: tension-stiffening average tensile stress carried between cracks by
// bonded reinforcement (a stress FLOOR above the bare fracture-energy softening).
//   mode 1 vc (MCFT / Bentz):     sigma_ts = ft / (1 + sqrt(c   * e1))
//   mode 2 cm (Collins-Mitchell): sigma_ts = alpha * ft / (1 + sqrt(500 * e1))
// e1 = COMPOSITE (reinforced) in-plane principal tensile strain (same membrane e1
// the MCFT beta uses; perfect-bond strain compatibility). Writes *dsig = d sigma_ts/d e1.
// ---------------------------------------------------------------------------
inline double tensStiff(double e1, int mode, double ft, double c, double alpha, double* dsig)
{
    double cc, ft_eff;
    if      (mode == 1) { cc = c;     ft_eff = ft; }
    else if (mode == 2) { cc = 500.0; ft_eff = alpha * ft; }
    else { if (dsig) *dsig = 0.0; return 0.0; }
    double u = (e1 > 0.0) ? sqrt(cc * e1) : 0.0;
    double denom = 1.0 + u;
    double sig = ft_eff / denom;
    if (dsig) *dsig = (u > 0.0) ? (-ft_eff * (cc / (2.0 * u)) / (denom * denom)) : 0.0;
    return sig;
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
    // SCALE-RELATIVE fallback: when the eigenvector magnitude is numerical dust (tiny gxy +
    // a near-cancelling vy that rounds to a sub-ULP nonzero), normalizing it would freeze a
    // garbage crack normal (e.g. (0,-1) instead of (1,0)). Fall back to the major-strain axis.
    if (nrm < 1.0e-9 * (fabs(half_diff) + fabs(half_gxy) + 1.0e-30)) {
        p1[0] = (exx >= eyy) ? 1.0 : 0.0;
        p1[1] = (exx >= eyy) ? 0.0 : 1.0;
    }
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
    bool   interlockCyclic;          // Phase 2b: incremental friction-slip + reversible w
    bool   xcrackOn;                 // Phase 2b.2b: 2nd orthogonal crack + slip-driven wear
    double degKappa;                 // wear knockdown slope (default 0.5)
    double degSlipRef;               // reference peak slip for full wear (<=0 => no wear)
    double degMin;                   // residual cap fraction floor (default 0.1)
    // --- Phase 2b.2c: -shearRetention {mcft|const|dsfm|rots} crack-shear retention curve ---
    // Selects the CYCLIC crack-plane slip stiffness G_slip used in the friction predictor
    // tau_cr = clamp(tau_cr_c + G_slip*dgamma_nt, +/- v_ci,max(w)). The v_ci,max(w) CAP is
    // unchanged in every mode (keeps the ADR bound). Modes (all reduce to mcft):
    //   0 mcft  (default): G_slip = G = E/2(1+nu)  -- the shipped full-elastic slip stiffness.
    //   1 const          : G_slip = shearRetFactor*G -- classic constant shear retention
    //                      (Rots/Cervenka). shearRetFactor in (0,1]; ==1 -> mcft.
    //   2 dsfm           : G_slip = G*(0.31/denom)  -- DSFM-flavored width-degraded slip
    //                      stiffness (softens as the crack opens). At w=0 -> mcft.
    //   3 rots           : rotating-coaxial -- NO independent crack-shear; the membrane shear
    //                      stays smeared/spectral (the interlock shear block is skipped). The
    //                      ADR's monotonic-only rotating choice; lets a user A/B vs fixed-crack.
    // const/dsfm parametrize the slip stiffness, so they only bite on the -cyclic path; on the
    // 2a monotone-clip path they are inert (same v_ci,max plateau). The parser implies the
    // needed flags (const/dsfm -> -interlock -cyclic; rots -> -interlock).
    int    shearRetMode;             // 0 mcft (default) | 1 const | 2 dsfm | 3 rots
    double shearRetFactor;           // const-mode retention factor mu in (0,1] (default 0.4)
    double aggSize;                   // max aggregate size a_g (SI: mm)
    double crackStrain;              // eps_cr crack-formation strain (<=0 => ftmax/E)
    double crackSpacing;            // s_theta crack spacing (<=0 => lch, else 1)
    double lch;                     // characteristic length (s_theta fallback)
    double betaSrMin;              // residual shear-retention floor (tangent safety)
    double sqrtFc;                // sqrt(fc') cache, set by setupParams (SI: sqrt(MPa))
    // --- Phase 4 (pulled forward): IMPL-EX (implicit-explicit) ---
    bool   implex;                 // false => fully implicit (default, baseline-identical)
    double implexAlpha;            // extrapolation factor multiplier (default 1.0)
    bool   implexControl;          // adaptive time-step error control (advisory)
    double implexErrTol;           // error tolerance for -implexControl (default 0.05)
    double implexTimeRedLim;       // min dt fraction for -implexControl (default 0.01)
    // --- Phase 3: tension stiffening (default off => baseline-identical) ---
    int    tensStiffMode;          // 0 off | 1 vc (Bentz) | 2 cm (Collins-Mitchell)
    double tensStiffC;             // vc-mode sqrt coefficient c (default 500)
    double tensStiffAlpha;         // cm-mode alpha1*alpha2 (default 1)
    double ftPeak;                 // tension-backbone peak ft, cached by setupParams
    // --- Phase 3b: crack-band (Bazant-Oh) regularization (default off => baseline-identical) ---
    bool   autoReg;                // -autoRegularization: scale softening so G_f is mesh-objective
    double lchRef;                 // reference characteristic length the backbone was authored at
    Backbone ht, hc;                   // tension / compression backbones
};

inline void setupParams(Params& P)    // call after ht/hc are filled
{
    double fcmax = backboneMaxStress(P.hc);
    double ftmax = backboneMaxStress(P.ht);
    P.fcft_ratio = (ftmax > 0.0) ? (fcmax / ftmax > 5.0 ? fcmax / ftmax : 5.0) : 1.0;
    P.sqrtFc = fcmax > 0.0 ? sqrt(fcmax) : 0.0;
    P.ftPeak = ftmax;
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
    // --- Phase 2b: cyclic crack-shear friction-slip state ---
    double tauCr;           // committed crack-plane shear stress (incremental friction)
    double gammaCr;         // committed crack-plane engineering shear strain (slip)
    // --- Phase 2b.2b: second orthogonal crack + interlock-surface wear ---
    double cracked2;        // 0/1 second (orthogonal) crack flag
    double slipCum;         // cumulative |crack slip| (abrasion/Archard wear driver, irreversible)
    // --- Phase 4: IMPL-EX previous-committed (n-1) generation for explicit extrapolation ---
    double xt_old, xc_old;  // committed tensile/compressive equiv strain at step n-1
    double eps1_old;        // committed in-plane principal tensile strain at step n-1
};

inline void histZero(RCHist& h)
{
    for (int i = 0; i < 6; ++i) { h.stress_eff[i] = 0.0; h.strain[i] = 0.0; }
    h.xt = h.xc = 0.0; h.dt_bar = h.dc_bar = 0.0; h.beta = 1.0; h.eps1 = 0.0;
    h.cracked = 0.0; h.crackC = 1.0; h.crackS = 0.0; h.wmax = 0.0; h.betaSr = 1.0;
    h.tauCr = 0.0; h.gammaCr = 0.0;
    h.cracked2 = 0.0; h.slipCum = 0.0;
    h.xt_old = h.xc_old = 0.0; h.eps1_old = 0.0;
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
// do_implex: when true (and P.implex), the damage thresholds (xt,xc) and the MCFT
//   beta are taken from an EXPLICIT extrapolation of the committed history
//   (x_ext = x_n + tf*(x_n - x_{n-1}), tf = time_factor), FROZEN over the step. The
//   returned tangent is then the secant W_B*C0 with frozen dt_bar/dc_bar/beta (no
//   softening-rate term, no beta cross-term) -> removes the dominant indefinite
//   contribution (the damage rate) and is robust for cyclic softening. The TRUE
//   implicit thresholds are recomputed at commit (do_implex=false) to advance the
//   state + measure the implex error. Mirrors ASDConcrete3DMaterial::compute(do_implex).
//   SCOPE (honest): unlike the reference, the spectral projectors PT/PC are NOT frozen
//   here -- they are recomputed from the live effective stress every call, consistent
//   with this material's fixed-projector-secant philosophy (dP/deps is omitted in ALL
//   tangents, implicit included). So the tangent is NOT perfectly strain-constant under
//   a rotating principal axis; the damage (the softening source) IS frozen, which is the
//   robustness that matters. The interlock/friction block (below) likewise runs its
//   normal implicit logic in the explicit pass (it is SPD/secant-character). Freezing
//   PT_commit for a fully strain-constant tangent under rotating principals is a scoped
//   follow-up (see ADR 19 Phase 4).
inline int returnMap3D(const Params& P, const double eps6[6], const RCHist& in,
                       double sig6[6], double Dtan6[6][6], RCHist& out,
                       bool do_tangent = true, double* residual = 0,
                       bool do_implex = false, double time_factor = 1.0)
{
    const bool implexExplicit = (P.implex && do_implex);
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
    // IMPL-EX: freeze beta over the step from the extrapolated (committed) eps1, so the
    // returned tangent carries no beta cross-term (constant secant). eps1 itself is a
    // pure function of strain, so this is the only path-frozen quantity here.
    double e1_beta = e1;
    if (implexExplicit && P.betaOn) {
        e1_beta = in.eps1 + time_factor * (in.eps1 - in.eps1_old);
        if (e1_beta < 0.0) e1_beta = 0.0;
    }
    double beta = P.betaOn ? betaCompr(e1_beta, P.betaFloor) : 1.0;
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
    if (implexExplicit) {
        // explicit extrapolation of the damage thresholds (frozen over the step)
        double xt_ext = in.xt + time_factor * (in.xt - in.xt_old);
        double xc_ext = in.xc + time_factor * (in.xc - in.xc_old);
        if (xt_ext < in.xt) xt_ext = in.xt;     // threshold is monotone non-decreasing
        if (xc_ext < in.xc) xc_ext = in.xc;
        pt = evaluateAt(P.ht, xt_ext); new_xt = pt.x;
        pc = evaluateAt(P.hc, xc_ext); new_xc = pc.x;
    } else {
        if (xt_trial > pt.x) { double xn = rc1 * pt.x + rc2 * xt_trial; pt = evaluateAt(P.ht, xn); new_xt = pt.x; }
        if (xc_trial > pc.x) { double xn = rc1 * pc.x + rc2 * xc_trial; pc = evaluateAt(P.hc, xn); new_xc = pc.x; }
    }

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

    // ---- Phase 3: tension stiffening (rank-1 FLOOR along the live principal tensile axis) ----
    // Between cracks, bonded reinforcement holds the average concrete tension above the bare
    // softening curve: inject delta = sigma_ts(e1) - n^T sigma n along p1 (only when delta>0),
    // active ONLY post-crack (e1 >= eps_cr) so the elastic pre-crack branch is untouched (the
    // gate is mandatory: pre-crack the bare stress ~E*e1 < sigma_ts(e1) would inject a spurious
    // floor). Degenerate (equibiaxial) -> isotropic in-plane. Default off => baseline-identical.
    // Inserted BEFORE the interlock block so the interlock shear clip reads the floored normal
    // stress (a normal stress on the ~p1 plane contributes ~0 crack-plane shear). IMPL-EX: e1 is
    // frozen from the extrapolated committed eps1 (like beta), and the tangent cross-term is
    // omitted below (constant secant).
    // SCOPE / LIMITATIONS (v1, documented boundaries -- see LadrunoRCConcrete_guide + LEDGER_quirks):
    //  (1) MONOTONIC backbone floor: sigma_ts is a pure function of the LIVE e1 (no eps1max
    //      memory). sigma_ts DECREASES with e1, so on UNLOADING (e1 decreasing) the floor
    //      RE-INFLATES (tracks sigma_ts(live e1) back up). This is correct on the monotone
    //      loading branch (the slab / distributed-reinforcement pushover use case TS is scoped
    //      for) but is NOT a hysteretic cyclic-tension model -- use -tensStiff for MONOTONIC /
    //      pushover analyses. (A future eps1max-envelope + secant-unload is the cyclic upgrade.)
    //  (2) TS uses the LIVE principal axis p1 (rotating-crack MCFT view). When combined with the
    //      FIXED-crack -interlock block (which uses the frozen crack normal), TS's normal-stress
    //      injection adds a small shear on the FROZEN crack plane once principal axes rotate
    //      (delta*sin(th_rel)cos(th_rel)) which the interlock then bounds. Combined TS+interlock
    //      is validated for NON-rotating (proportional) loading; rotating-axis combined use is a
    //      documented boundary, not a guaranteed-correct superposition.
    double e1_ts = e1;
    if (implexExplicit && P.tensStiffMode != 0) {
        e1_ts = in.eps1 + time_factor * (in.eps1 - in.eps1_old);
        if (e1_ts < 0.0) e1_ts = 0.0;
    }
    // Three in-plane vectors over Voigt indices {0,1,3}, set per degenerate state:
    //   ts_inj  = injection weights      (sig_i += ts_inj_i * delta)
    //   ts_meas = the measured-quantity  (q = ts_meas . sig_ip), pinned to sigma_ts
    //   ts_P1eps= de1/deps strain dual (6-vec, membrane only)
    // NON-degenerate: q = n^T sig n on p1 (rank-1) -> ts_inj=(a,b,ab), ts_meas=(a,b,2ab),
    //   self-consistency ts_meas.ts_inj = a^2+b^2+2(ab)^2 = (p1x^2+p1y^2)^2 = 1.
    // DEGENERATE (equibiaxial, p1 ill-defined): q = in-plane MEAN normal 0.5(s0+s1);
    //   inject g=sigma_ts-q to BOTH normals -> ts_inj=(1,1,0), ts_meas=(0.5,0.5,0),
    //   ts_meas.ts_inj = 1 (each in-plane normal reaches sigma_ts; the rank-1 (0.5,0.5,0)
    //   reuse would only reach the half-way point -- self-consistency coeff 0.5).
    bool   ts_active = false;
    double ts_inj[3] = {0,0,0}, ts_meas[3] = {0,0,0}, ts_P1eps[6] = {0,0,0,0,0,0}, ts_dsig = 0.0;
    if (P.tensStiffMode != 0 && e1_ts >= P.crackStrain) {
        if (degen) {
            ts_inj[0] = 1.0;  ts_inj[1] = 1.0;
            ts_meas[0] = 0.5; ts_meas[1] = 0.5;
            ts_P1eps[0] = 0.5; ts_P1eps[1] = 0.5;         // de1/d(exx)=de1/d(eyy)=0.5 (mean)
        } else {
            double a = p1[0]*p1[0], b = p1[1]*p1[1], ab = p1[0]*p1[1];
            ts_inj[0] = a;  ts_inj[1] = b;  ts_inj[2] = ab;
            ts_meas[0] = a; ts_meas[1] = b; ts_meas[2] = 2.0*ab;
            ts_P1eps[0] = a; ts_P1eps[1] = b; ts_P1eps[3] = ab;   // de1/deps (membrane)
        }
        double q      = ts_meas[0]*sig6[0] + ts_meas[1]*sig6[1] + ts_meas[2]*sig6[3];
        double sig_ts = tensStiff(e1_ts, P.tensStiffMode, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, &ts_dsig);
        double delta  = sig_ts - q;
        if (delta > 0.0) {
            sig6[0] += delta * ts_inj[0];
            sig6[1] += delta * ts_inj[1];
            sig6[3] += delta * ts_inj[2];
            ts_active = true;
        }
    }

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
    // Phase 2a (default): a v_ci,max BOUND that CLIPS the already-assembled (damage-
    // reduced) crack-plane shear tau_sm = m_sigma.sig_ip to +/- v_ci,max -- NOT a bare-
    // elastic substitution; below the cap the stress is unchanged (continuous). w is
    // monotone (irreversible). m_sigma . m_eps == 1 => matched secant tangent.
    // Phase 2b (-cyclic): the crack-plane shear becomes an INDEPENDENT incremental
    // friction-slip state (FSAM-style): tau_cr = clamp(tau_cr_committed + G*dgamma_nt,
    // +/- v_ci,max), with v_ci,max(w) now REVERSIBLE (crack closure raises the cap).
    // Cyclic reversal unloads along G and re-caps in the opposite sense; the cap
    // collapsing as the crack re-opens (w grows) gives pinching. Reduces to a monotone
    // ramp-then-cap (same plateau as 2a) under monotonic loading.
    double crk_c = in.crackC, crk_s = in.crackS, cracked = in.cracked, wmax = in.wmax;
    double tauCr = in.tauCr, gammaCr = in.gammaCr;   // 2b cyclic crack-shear state
    double cracked2 = in.cracked2, slipCum = in.slipCum;   // 2b.2b X-crack + wear state
    double betaSr = 1.0;
    bool   il_active = false, il_capped = false, il_cyclic = false;
    double il_m[3] = { 0.0, 0.0, 0.0 };             // m_eps: stress back-rotation of a shear change
    double il_Gint = 0.0;                            // interlock shear modulus (cyclic tangent)
    // shearRetMode 3 (rots) = rotating-coaxial: skip the fixed-crack shear entirely, so the
    // membrane shear stays smeared/spectral (no capture, no clip, no friction). il_active
    // stays false => the tangent carries no crack cross-term either.
    if (P.interlockOn && P.shearRetMode != 3) {
        if (cracked < 0.5 && e1 >= P.crackStrain && !degen) {   // freeze crack normal
            cracked = 1.0; crk_c = p1[0]; crk_s = p1[1];
        }
        bool just_captured = (in.cracked < 0.5 && cracked >= 0.5);
        if (cracked >= 0.5) {
            il_active = true;
            double exx = eps6[0], eyy = eps6[1], gxy = eps6[3];
            double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
            double en = exx*c2 + eyy*s2 + gxy*cs;               // strain normal to the FROZEN crack 1
            double sth = (P.crackSpacing > 0.0) ? P.crackSpacing : (P.lch > 0.0 ? P.lch : 1.0);
            // ---- Phase 2b.2b: second orthogonal crack (X-cracking) ----
            // Crack 2 normal = (-s,c) (orthogonal to crack 1), frozen when the strain normal to
            // THAT plane reaches eps_cr (FSAM CepsA2>=eunpA2 + orthogonality lock). The cap then
            // uses the GOVERNING (most-open) crack's opening, so the REVERSE shear direction is
            // interlock-capped too (2a/2b.1 single crack only capped the crack-1 direction).
            double en_gov = en;
            if (P.xcrackOn && P.interlockCyclic) {
                double c2n = -s, s2n = c;                       // crack-2 (orthogonal) normal
                double en2 = exx*c2n*c2n + eyy*s2n*s2n + gxy*c2n*s2n;
                if (cracked2 < 0.5 && en2 >= P.crackStrain) cracked2 = 1.0;
                if (cracked2 >= 0.5 && en2 > en_gov) en_gov = en2;   // most-open crack governs
            }
            double w = macauley(en_gov) * sth;
            if (w > wmax) wmax = w;                              // tracked max (diagnostic / 2a)
            double w_use = P.interlockCyclic ? w : wmax;        // reversible (cyclic) vs monotone (2a)
            double denom = 0.31 + 24.0 * w_use / (P.aggSize + 16.0);
            double vcimax = (denom > 0.0) ? 0.18 * P.sqrtFc / denom : 0.18 * P.sqrtFc;
            // ---- Phase 2b.2b: slip-driven irreversible interlock-surface wear (degradation) ----
            // Aggregate-interlock abrasion (Archard) scales with the CUMULATIVE crack sliding
            // distance; knock the cap down by kn = max(degMin, 1 - degKappa*clamp(slipCum/degSlipRef,
            // 0,1)). Driven by the COMMITTED slipCum (in.slipCum) so it is tangent-neutral. Unlike a
            // peak-slip law (which saturates in cycle 1 under constant amplitude), cumulative sliding
            // keeps wearing every reversal -> cyclic strength decay (and the low-cap pinched waist).
            if (P.xcrackOn && P.degSlipRef > 0.0) {
                double kn = 1.0 - P.degKappa * (in.slipCum / P.degSlipRef);
                if (kn < P.degMin) kn = P.degMin;
                if (kn > 1.0) kn = 1.0;
                vcimax *= kn;
            }
            double sxx = sig6[0], syy = sig6[1], sxy = sig6[3];
            double tau_sm = cs * (syy - sxx) + (c2 - s2) * sxy; // smeared crack-plane shear
            double tau_ci;
            if (P.interlockCyclic) {
                il_cyclic = true;
                // crack-shear slip stiffness G_slip per -shearRetention mode (cap unchanged).
                double Gfull = 0.5 * E / (1.0 + P.nu);          // full elastic shear modulus G
                il_Gint = Gfull;                                 // mcft (default): full G
                if (P.shearRetMode == 1)      il_Gint = P.shearRetFactor * Gfull;  // const: mu*G
                else if (P.shearRetMode == 2) il_Gint = Gfull * (0.31 / denom);    // dsfm: width-degraded
                double g_nt = 2.0*(eyy - exx)*cs + gxy*(c2 - s2);
                // 2b.2b wear driver: accumulate the sliding distance (skip the capture step,
                // where in.gammaCr is still the uncracked 0 and would inject a spurious jump).
                slipCum = in.slipCum + (just_captured ? 0.0 : fabs(g_nt - in.gammaCr));
                if (just_captured) {                            // seed for stress CONTINUITY:
                    gammaCr = g_nt;                             // slip origin = capture slip,
                    tauCr = (tau_sm >  vcimax) ?  vcimax        // tau origin = the clipped smeared
                          : (tau_sm < -vcimax) ? -vcimax : tau_sm;  // shear (== the 2a value) so the
                }                                               // capture step is continuous (no jump)
                tau_ci = tauCr + il_Gint * (g_nt - gammaCr);    // incremental friction predictor
                if (tau_ci >  vcimax) { tau_ci =  vcimax; il_capped = true; }
                if (tau_ci < -vcimax) { tau_ci = -vcimax; il_capped = true; }
                gammaCr = g_nt; tauCr = tau_ci;                 // advance committed slip/stress
                // betaSr in cyclic = cap UTILIZATION tau_ci/v_ci,max (in [-1,1]); the 2a
                // retention-secant tau_ci/tau_sm is non-physical here (tau_ci is independent).
                betaSr = (vcimax > 1.0e-30) ? tau_ci / vcimax : 0.0;
            } else {
                tau_ci = tau_sm;                                // 2a: clip the smeared shear
                if (tau_ci >  vcimax) { tau_ci =  vcimax; il_capped = true; }
                if (tau_ci < -vcimax) { tau_ci = -vcimax; il_capped = true; }
                betaSr = (fabs(tau_sm) > 1.0e-30) ? tau_ci / tau_sm : 1.0;
            }
            double dtau = tau_ci - tau_sm;
            il_m[0] = -2.0*cs; il_m[1] = 2.0*cs; il_m[2] = c2 - s2;   // m_eps
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
    out.tauCr = tauCr; out.gammaCr = gammaCr;
    out.cracked2 = cracked2; out.slipCum = slipCum;
    // IMPL-EX (n-1) generation: carried through verbatim; the rolling n->n-1 is done
    // by the host material at commit (so a mid-step explicit pass never corrupts it).
    out.xt_old = in.xt_old; out.xc_old = in.xc_old; out.eps1_old = in.eps1_old;

    // tangent
    if (do_tangent) {
        // Under IMPL-EX the numerical (forward-difference) tangent is INVALID: the base
        // point is the explicit stress while a perturbed re-call would re-enter implicitly
        // (do_implex defaults false), so the quotient mixes branches (O(implexError/h),
        // huge). IMPL-EX always uses the frozen-damage secant below — matching the
        // reference, which disables the numerical tangent whenever implex is on.
        if (P.tangentMode == 2 && !implexExplicit) {
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
            if (P.tangentMode == 0 && P.betaOn && !implexExplicit) {
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
            // Phase-3 tension-stiffening CONSISTENT tangent. The floor adds, to the in-plane
            // stress rows i in {0,1,3} with weights ts_inj=(a,b,ab):
            //   sig_i += ts_inj_i * (sigma_ts(e1) - n^T sig n),   n^T sig n = a*s0 + b*s1 + 2ab*s3.
            //   d sig_i/d eps_c += ts_inj_i * [ d sigma_ts/de1 * de1/d eps_c  -  d(n^T sig n)/d eps_c ].
            // de1/d eps is the membrane dual P1eps = (a,b,ab) over {0,1,3}, ZERO out of plane.
            // d(n^T sig n)/d eps_c = ts_meas : Dtan[.][c] spans ALL 6 columns (the bare in-plane
            // stresses depend on out-of-plane strains too, e.g. eps_zz) -- pinning sig to sigma_ts(e1)
            // must remove that full sensitivity, else the eps_zz column is left at the (now-pinned-
            // away) bare value. dp1/d eps (projector rotation) is OMITTED, matching the fixed-
            // projector-secant character of the W_B/beta tangents. Skipped under IMPL-EX (e1 frozen
            // => constant secant). Placed before the interlock tangent so the interlock cross-term
            // reads the floored Dtan (consistent with the floored stress order).
            if (P.tangentMode == 0 && ts_active && !implexExplicit) {
                const int vidx[3] = { 0, 1, 3 };
                double row[6];                                              // d(q)/deps, all cols
                for (int c = 0; c < 6; ++c) {
                    double acc = 0.0;
                    for (int k = 0; k < 3; ++k) acc += ts_meas[k] * Dtan6[vidx[k]][c];
                    row[c] = acc;
                }
                for (int a = 0; a < 3; ++a)
                    for (int c = 0; c < 6; ++c)
                        Dtan6[vidx[a]][c] += ts_inj[a] * (ts_dsig * ts_P1eps[c] - row[c]);
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
            if (il_cyclic) {
                // Phase-2b cyclic consistent tangent. Crack-shear is INDEPENDENT
                // (tau_cr = tau_cr_c + G*dgamma_nt, clamped): remove the smeared crack-shear
                // sensitivity m_eps (x) (m_sigma.Dtan_ip), then add the friction stiffness
                // G*m_eps (x) m_eps (full when sliding-elastic; floored to betaSrMin*G when
                // capped/plastic for conditioning). v_ci,max(w) coupling omitted (2nd order).
                double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
                const int vidx[3] = { 0, 1, 3 };
                double msig[3] = { -cs, cs, c2 - s2 };
                double row[3];
                for (int j = 0; j < 3; ++j) {
                    double acc = 0.0;
                    for (int k = 0; k < 3; ++k) acc += msig[k] * Dtan6[vidx[k]][vidx[j]];
                    row[j] = acc;
                }
                double gcoef = (il_capped ? P.betaSrMin : 1.0) * il_Gint;
                for (int a = 0; a < 3; ++a)
                    for (int c3 = 0; c3 < 3; ++c3)
                        Dtan6[vidx[a]][vidx[c3]] += -il_m[a] * row[c3] + gcoef * il_m[a] * il_m[c3];
            } else if (il_capped) {
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
