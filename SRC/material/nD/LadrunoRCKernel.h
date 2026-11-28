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
//      beta  = clamp( 1/(0.8 + C eps_1), betaFloor, 1 ),  C = P.betaC (default 170,
//              -betaC; Vecchio-Collins 1986 is 0.34/|eps'c|, e.g. 189 for PV20)
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
//
//  Ladruno WP-263 (dual tangent): returnMap3D is written ONCE as returnMap3DT<T>. T = double is
//  the production path, arithmetic unchanged (bit-identical to the pre-WP kernel). T =
//  ladruno_dual::Dual<6> with the six engineering strains seeded gives d(sigma)/d(eps) of the map
//  actually executed, in one pass: tangentMode 3 (-dualTangent). The value part of the Dual run
//  equals the double run bit for bit. Three helpers do not push the Dual through their own
//  arithmetic, because the eigen-solver and the membrane eigenvector are kinked or ill-conditioned
//  exactly where RC states sit (repeated principals, gxy = 0): decompose, eps1FromMembrane and
//  eps1Spatial take their VALUES from the double function and their DERIVATIVES from closed-form
//  isotropic-function / eigenprojector formulas (Dual overloads below each double function).
// ============================================================================

#include <math.h>
#include <type_traits>
#include "LadrunoDual.h"   // Ladruno WP-263: dual numbers (forward-mode AD)

namespace ladruno_rc_kernel {

using ladruno_dual::value;   // WP-263: value(x) is the identity on double, the value part on a Dual

// WP-263: makes a parameter a non-deduced context, so tensStiff(e1, ..., 0) still deduces T from e1.
template <class T> struct rc_ident { typedef T type; };

// WP-263: sqrt at a cone apex (J2 = 0, the hydrostatic axis). The value is ::sqrt; a Dual at a
// non-positive argument gets a ZERO derivative (the subgradient at the apex) instead of 0.5/0.
inline double sqrtKink(double x) { return sqrt(x); }
template <int N, class S>
inline ladruno_dual::Dual<N, S> sqrtKink(const ladruno_dual::Dual<N, S>& x)
{
    if (!(value(x) > 0.0)) return ladruno_dual::Dual<N, S>(sqrtKink(x.v));
    return sqrt(x);
}

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

// WP-263: templated on the strain scalar (T = double is the original arithmetic, term for term).
template <class T> struct HPointT { T x, y, d, q; };
typedef HPointT<double> HPoint;

template <class T>
inline HPointT<T> evaluateAt(const Backbone& b, T x)
{
    double x1, y1, q1;
    T x2, y2, q2;
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
        T span = x - x1;
        y1 = b.y[last];
        double tan_y = (b.y[last] - b.y[last - 1]) / (b.x[last] - b.x[last - 1]);
        y2 = tan_y > 0.0 ? T(y1 + span * tan_y) : T(y1);
        q1 = b.q[last];
        double tan_q = (b.q[last] - b.q[last - 1]) / (b.x[last] - b.x[last - 1]);
        q2 = tan_q > 0.0 ? T(q1 + span * tan_q) : T(q1);
    }
    T xspan  = x2 - x1;
    T xratio = xspan > 0.0 ? T((x - x1) / xspan) : T(0.0);
    T y = y1 + (y2 - y1) * xratio; if (y < YTOL) y = YTOL;
    T q = q1 + (q2 - q1) * xratio; if (q < YTOL) q = YTOL;
    HPointT<T> p; p.x = x; p.y = y; p.q = q; p.d = 1.0 - y / q;
    return p;
}

template <class T>
inline T plasticStrain(const HPointT<T>& p, double E) { return p.x - p.q / E; }

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

// Re-enforce the HardeningLaw invariants on an already-built backbone (the clone of
// ASDConcrete3D HardeningLaw::adjust used inside regularize): cap each secant slope at E,
// keep the plastic strain (x - q/E) monotone non-decreasing, keep damage non-decreasing, and
// recompute q = y/(1-d). d is derived from the stored q (d = 1 - y/q). Mirrors the inline
// adjust block in buildBackbone; xtol/ytol are inert post-scaling (points stay strictly
// increasing) but kept for parity. Without this, scaling can drive the post-peak plastic
// strain backward and corrupt q/damage on steep-softening backbones (energy stays correct).
inline void adjustBackbone(Backbone& b, double E, double xtol, double ytol)
{
    const int n = b.n;
    double Eloc = (n > 1 && b.x[1] > 0.0) ? b.y[1] / b.x[1] : E;
    for (int i = 1; i < n; ++i) {
        double xi = b.x[i], xold = b.x[i-1], yi = b.y[i], yold = b.y[i-1];
        double di   = b.q[i]   > 0.0 ? 1.0 - yi   / b.q[i]   : 0.0;
        double dold = b.q[i-1] > 0.0 ? 1.0 - yold / b.q[i-1] : 0.0;
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
        if (di > 1.0 - 1.0e-12) di = 1.0 - 1.0e-12;        // guard q against 1/(1-d) blow-up
        Ei = (yi - yold) / (xi - xold);
        double Ed = (1.0 - di) * Eloc;
        if (Ei > Ed) yi = Ed * (xi - eepd_old);
        b.x[i] = xi; b.y[i] = yi; b.q[i] = yi / (1.0 - di);
    }
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
        adjustBackbone(b, E, 0.0, 0.0);                    // re-enforce monotonicity (ASD parity)
        int bd, p1, p2;
        g = fractureEnergy(b, E, &bd, &p1, &p2);
        if (fabs(g - gnew) < tol) break;
        dscale = (g > 0.0) ? gnew / g : 1.0;
    }
    adjustBackbone(b, E, 0.0, 0.0);                        // final re-adjust (ASDConcrete3D:1000)
}

template <class T>
inline T macauley(const T& x) { return x > 0.0 ? x : T(0.0); }

// ---------------------------------------------------------------------------
// Lubliner-Lee-Fenves criterion (clone of cpp:2483-2497).
// Local shape parameter renamed beta_lub to reserve `beta` for MCFT softening.
// WP-263: templated on the principal stresses; sqrt(3 J2) has a zero derivative at J2 = 0.
// ---------------------------------------------------------------------------
template <class T>
inline T lublinerCriterion(const T& s1, const T& s2, const T& s3,
                           double ft, double fc, double k1, double scale, double Kc)
{
    double fb    = 1.16 * fc;
    double gamma = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0);
    double alpha = (fb - fc) / (2.0 * fb - fc);
    T I1    = s1 + s2 + s3;
    T p     = I1 / 3.0;
    T d1 = s1 - p, d2 = s2 - p, d3 = s3 - p;
    T J2 = 0.5 * (d1 * d1 + d2 * d2 + d3 * d3);
    double beta_lub = fc / ft * (1.0 - alpha) - (1.0 + alpha);
    T smax = macauley(s1);
    T smin = macauley(T(-s1));
    return (1.0 / (1.0 - alpha) *
            (alpha * I1 + sqrtKink(3.0 * J2) + k1 * beta_lub * smax - gamma * smin)) * scale;
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
// WP-263: Si, ST, SC and R carry the stress scalar T; the eigenvectors and the 6x6 projectors stay
// double (the secant tangents of modes 0/1 read them on the double path only).
template <class T>
struct DecompT {
    T Si[3]; double V[3][3];
    double PT[6][6], PC[6][6];
    T ST[6], SC[6];
    T R;
};
typedef DecompT<double> Decomp;

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

// WP-263: decompose for a Dual stress. Values: the double decompose of the value part (so every
// value is bit-identical to the double run). Derivatives, never through the Jacobi solver (rule 6):
//   * ST = sum_j <lambda_j>_+ v_j v_j^T is an isotropic tensor function of S, so in the eigenbasis
//       d(ST)~_ij = G_ij dS~_ij,  G_ii = H(lambda_i),
//       G_ij = (<lambda_i>_+ - <lambda_j>_+) / (lambda_i - lambda_j)   (lambda_i != lambda_j),
//       G_ij = (H(lambda_i) + H(lambda_j)) / 2                          (lambda_i == lambda_j),
//     with H(0) = 1/2 as heavyside(); this is the exact derivative of PT(S):S, projector rotation
//     included, and it has no 1/(lambda_i - lambda_j) singularity. SC = S - ST, so dSC = dS - dST.
//   * Si: d lambda_i = dS~_ii, averaged over a cluster of equal principals (|gap| < 1e-9 of the
//     stress scale): the derivative of the symmetric functions the measures use, independent of
//     the arbitrary basis Jacobi returns inside a repeated eigenspace.
//   * R: the double formula, templated, on the Dual Si.
template <int N>
inline void decompose(const ladruno_dual::Dual<N> Sd[6], double cdf, DecompT<ladruno_dual::Dual<N> >& D)
{
    typedef ladruno_dual::Dual<N> DN;
    static const int vr[6] = {0, 1, 2, 0, 1, 0}, vc[6] = {0, 1, 2, 1, 2, 2};
    double sv[6];
    for (int a = 0; a < 6; ++a) sv[a] = Sd[a].v;
    Decomp Dv;
    decompose(sv, cdf, Dv);
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) D.V[i][j] = Dv.V[i][j];
    for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) { D.PT[a][b] = Dv.PT[a][b]; D.PC[a][b] = Dv.PC[a][b]; }
    for (int a = 0; a < 6; ++a) { D.ST[a] = DN(Dv.ST[a]); D.SC[a] = DN(Dv.SC[a]); }
    for (int i = 0; i < 3; ++i) D.Si[i] = DN(Dv.Si[i]);

    const double* l = Dv.Si;
    double scale = 0.0;
    for (int a = 0; a < 6; ++a) scale += fabs(sv[a]);
    scale += 1.0e-30;
    const double ctol = 1.0e-9 * scale;
    const bool g01 = (l[0] - l[1] < ctol), g12 = (l[1] - l[2] < ctol);
    double G[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j)              G[i][j] = heavyside(l[i]);
            else if (l[i] == l[j])   G[i][j] = 0.5 * (heavyside(l[i]) + heavyside(l[j]));
            else                     G[i][j] = (macauley(l[i]) - macauley(l[j])) / (l[i] - l[j]);
        }
    for (int k = 0; k < N; ++k) {
        const double M[3][3] = {
            { Sd[0].d[k], Sd[3].d[k], Sd[5].d[k] },
            { Sd[3].d[k], Sd[1].d[k], Sd[4].d[k] },
            { Sd[5].d[k], Sd[4].d[k], Sd[2].d[k] }
        };
        double At[3][3];                                  // dS~ = V^T dS V
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double acc = 0.0;
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) acc += Dv.V[r][i] * M[r][c] * Dv.V[c][j];
                At[i][j] = acc;
            }
        double dl[3] = { At[0][0], At[1][1], At[2][2] };
        if (g01 && g12) { const double m = (dl[0] + dl[1] + dl[2]) / 3.0; dl[0] = dl[1] = dl[2] = m; }
        else if (g01)   { const double m = 0.5 * (dl[0] + dl[1]); dl[0] = dl[1] = m; }
        else if (g12)   { const double m = 0.5 * (dl[1] + dl[2]); dl[1] = dl[2] = m; }
        for (int i = 0; i < 3; ++i) D.Si[i].d[k] = dl[i];
        for (int a = 0; a < 6; ++a) {
            double acc = 0.0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    acc += Dv.V[vr[a]][i] * Dv.V[vc[a]][j] * G[i][j] * At[i][j];
            D.ST[a].d[k] = acc;
            D.SC[a].d[k] = Sd[a].d[k] - acc;
        }
    }
    DN Rnum(0.0), Rden(0.0);
    for (int j = 0; j < 3; ++j) {
        DN Sj = D.Si[j];
        if (Sj > 0.0) Sj *= cdf;
        if (Sj > 0.0) Rnum += Sj;
        Rden += fabs(Sj);
    }
    D.R = Rden > 0.0 ? DN(Rnum / Rden) : DN(0.0);
}

// ---------------------------------------------------------------------------
// MCFT compression-softening factor (applied to the STRENGTH axis only).
// ---------------------------------------------------------------------------
// Ladruno (C2): the coefficient C (hard-wired 170 = 0.34/0.002 before C2) is P.betaC;
// the default 170.0 keeps the arithmetic bit-identical to the pre-C2 kernel.
template <class T>
inline T betaCompr(const T& e1, double floorv, double C = 170.0)
{
    T b = 1.0 / (0.8 + C * e1);
    if (b > 1.0) b = 1.0;
    if (b < floorv) b = floorv;
    return b;
}
inline double dBetaCompr(double e1, double floorv, double C = 170.0)   // 0 when clamped
{
    double b = 1.0 / (0.8 + C * e1);
    if (b >= 1.0 || b <= floorv) return 0.0;
    return -C * b * b;
}

// ---------------------------------------------------------------------------
// Phase 3: tension-stiffening average tensile stress carried between cracks by
// bonded reinforcement (a stress FLOOR above the bare fracture-energy softening).
//   mode 1 vc (MCFT / Bentz):     sigma_ts = ft / (1 + sqrt(c   * e1))
//   mode 2 cm (Collins-Mitchell): sigma_ts = alpha * ft / (1 + sqrt(500 * e1))
// e1 = COMPOSITE (reinforced) in-plane principal tensile strain (same membrane e1
// the MCFT beta uses; perfect-bond strain compatibility). Writes *dsig = d sigma_ts/d e1.
// ---------------------------------------------------------------------------
template <class T>
inline T tensStiff(const T& e1, int mode, double ft, double c, double alpha, typename rc_ident<T>::type* dsig)
{
    double cc, ft_eff;
    if      (mode == 1) { cc = c;     ft_eff = ft; }
    else if (mode == 2) { cc = 500.0; ft_eff = alpha * ft; }
    else { if (dsig) *dsig = 0.0; return T(0.0); }
    T u = (e1 > 0.0) ? T(sqrt(cc * e1)) : T(0.0);
    T denom = 1.0 + u;
    T sig = ft_eff / denom;
    if (dsig) *dsig = (u > 0.0) ? T(-ft_eff * (cc / (2.0 * u)) / (denom * denom)) : T(0.0);
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

// WP-263: eps1FromMembrane for a Dual strain. Values (e1, p1, degen) are the double function's.
// e1 = <c + r>_+ is differentiated through the same expressions (r > 0 off the degenerate branch;
// the degenerate branch is <c>_+, the equibiaxial blend). p1 is NOT differentiated through
// (gxy/2, lam1 - exx)/nrm: at gxy = 0 that vector cancels to zero (the fallback branch) although
// p1 rotates. Its derivative is the rotation of the major axis, 2 theta = atan2(gxy, exx - eyy):
//     d theta = ((exx - eyy) d gxy - gxy d(exx - eyy)) / (2 ((exx - eyy)^2 + gxy^2)),
//     d p1 = d theta (-p1y, p1x)   (either sign of p1). Zero on the degenerate branch.
template <int N>
inline ladruno_dual::Dual<N> eps1FromMembrane(const ladruno_dual::Dual<N>& exx, const ladruno_dual::Dual<N>& eyy,
                                              const ladruno_dual::Dual<N>& gxy, ladruno_dual::Dual<N> p1[2], int* degen)
{
    typedef ladruno_dual::Dual<N> DN;
    double pv[2]; int dg = 0;
    eps1FromMembrane(exx.v, eyy.v, gxy.v, pv, &dg);
    if (degen) *degen = dg;
    p1[0] = DN(pv[0]); p1[1] = DN(pv[1]);
    DN c = 0.5 * (exx + eyy);
    if (dg) return macauley(c);
    DN half_diff = 0.5 * (exx - eyy);
    DN half_gxy  = 0.5 * gxy;
    DN r = sqrt(half_diff * half_diff + half_gxy * half_gxy);
    DN lam1 = c + r;
    const double Dv = exx.v - eyy.v, gv = gxy.v, den = 2.0 * (Dv * Dv + gv * gv);
    for (int k = 0; k < N; ++k) {
        const double dth = (Dv * gxy.d[k] - gv * (exx.d[k] - eyy.d[k])) / den;
        p1[0].d[k] = -pv[1] * dth;
        p1[1].d[k] =  pv[0] * dth;
    }
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
    int    tangentMode;                // 0 algorithmic, 1 secant, 2 numerical, 3 dual (WP-263, -dualTangent)
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
    double tensStiffC;             // vc-mode sqrt coefficient c (parser default 200 = Vecchio-Collins
                                   // 1986 since C2; 500 before, = Collins-Mitchell 1991)
    double tensStiffAlpha;         // cm-mode alpha1*alpha2 (default 1)
    double ftPeak;                 // tension-backbone peak ft, cached by setupParams
    // --- Phase 3b: crack-band (Bazant-Oh) regularization (default off => baseline-identical) ---
    bool   autoReg;                // -autoRegularization: scale softening so G_f is mesh-objective
    double lchRef;                 // reference characteristic length the backbone was authored at
    Backbone ht, hc;                   // tension / compression backbones
    // --- Ladruno C2 (defaults = pre-C2 behaviour, so every existing Params builder is unchanged) ---
    double betaC = 170.0;          // MCFT softening coefficient C in beta = 1/(0.8 + C eps1) (-betaC)
    double crackedNu = -1.0;       // <0 => off. >=0: Poisson ratio of C0 once eps1 >= crackStrain
                                   // (irreversible latch RCHist::nuCracked) (-crackedNu)
    // --- Ladruno WP-241 (defaults = pre-WP-241 behaviour, byte-identical) ---
    // addonPlane: where the MCFT add-ons (beta, tension stiffening, cracked-nu latch) read the
    //   principal tensile strain eps1 and its direction p1.
    //   0 membrane (default): the x-y block {exx, eyy, gxy}, 2D p1 (shell / plane views).
    //   1 spatial: the full 3x3 strain tensor, 3D p1 (frame-covariant; for beam fibres, where
    //     the shear sits in gxy AND gxz). The fixed-crack interlock is x-y only and is refused
    //     with spatial (parser); the kernel skips it defensively.
    int    addonPlane = 0;         // -addonPlane {membrane|spatial}
    // crackGate: what opens the cracked-nu latch and the tension-stiffening floor.
    //   0 strain (default): eps1 >= crackStrain. eps1 is the TOTAL principal strain, so Poisson
    //     expansion under compression opens it with no crack (WP-241 probe: at ~0.42 fc in
    //     uniaxial compression of a plane-stress point).
    //   1 damage: tensile damage has started, i.e. the tensile threshold xt is past the peak of
    //     the tension backbone (latch: committed xt; floor: trial xt), and for the floor also
    //     eps1 >= crackStrain. The Lubliner tensile measure is zero without a tensile principal
    //     stress, so compression alone never opens it.
    int    crackGate = 0;          // -crackGate {strain|damage}
    // --- Ladruno WP-243 BeamFiber view: smeared transverse ties (read only by driveBeamFiber) ---
    // Tie ratios per lateral direction of the fibre frame (y = index 1, z = index 2), yield stress,
    // modulus and the bilinear kinematic hardening ratio b (Et = b Es once yielded). rho = 0 in a
    // direction leaves that lateral normal stress free (sigma = 0), so a cover fibre is the
    // default. Inert in every other view.
    double tieRhoY = 0.0, tieRhoZ = 0.0;   // -ties rhoY rhoZ fyt Es
    double tieFy = 0.0, tieEs = 0.0;
    double tieB = 0.0;                     // -tieHardening b, 0 <= b < 1
    // Model C crack check (BeamFiber view only, read by driveBeamFiber): 0 none | 1 sum.
    //   sum: the principal tension q = n1.sig.n1 <= v_ci,max tan(theta) + sum_i rho_i <fyt - f_si> w_i
    // with sin(theta) = |n1_x| (theta = strut angle from the fibre axis), crack width
    // w = e1 / (sin(theta)/smx + cos(theta)/smt), v_ci,max = 0.18 sqrt(fc) / (0.31 + 24 w/(a_g + 16))
    // (SI: MPa, mm) and tie weights w_i = n1_i^2 / (n1_y^2 + n1_z^2).
    int    crackCheck = 0;                 // -crackCheck sum smx smt (absent = 0)
    double ccSmx = 0.0, ccSmt = 0.0;       // crack-spacing parameters along x and transverse
    // Ladruno WP-247: what the BeamFiber view does when the lateral balance has no root (read by
    // driveBeamFiber): 0 refuse the trial (default) | 1 modelC, keep Model C's bisection end point.
    int    lateralNoRoot = 0;              // -lateralNoRoot {refuse|modelC}
};

// Strain abscissa of the tension backbone's maximum stress (the crack-onset threshold used by
// crackGate = damage). The post-peak regularization stretches only the softening branch, so
// this does not move under -autoRegularization.
inline double backbonePeakAbscissa(const Backbone& b)
{
    double ymax = -1.0e300, xp = 0.0;
    for (int i = 0; i < b.n; ++i) if (b.y[i] > ymax) { ymax = b.y[i]; xp = b.x[i]; }
    return xp;
}

// Ladruno WP-241: principal tensile strain of the FULL 3x3 strain tensor (engineering shear in,
// tensor shear gamma/2 used). Returns e1 = <lambda_max> and three Voigt 6-vectors over
// {00,11,22,01,12,02}, built from the m eigenvectors v_k that share lambda_max (multiplicity
// m = 1, 2 or 3, relative tolerance 1e-9 of the strain scale, as eps1FromMembrane):
//   P1[]   = d e1 / d eps (engineering dual) = (1/m) sum_k (vx^2, vy^2, vz^2, vx vy, vy vz, vx vz)
//   inj[]  = TS injection weights            =       sum_k (vx^2, vy^2, vz^2, vx vy, vy vz, vx vz)
//   meas[] = TS measured quantity weights    = (1/m) sum_k (vx^2, vy^2, vz^2, 2vx vy, 2vy vz, 2vx vz)
// so meas . inj = 1 for every m (each of the m normals is pinned to sigma_ts). For an x-y state
// whose largest principal lies in x-y these equal the membrane vectors exactly: m = 1 gives
// (a, b, ab) / (a, b, 2ab); the in-plane equibiaxial m = 2 gives (1,1,0) / (0.5,0.5,0).
// WP-263: the eigen-solve and the multiplicity of eps1Spatial, split out so its Dual overload reads
// the same eigenbasis and the same m (arithmetic unchanged).
inline void strainEig3(const double eps6[6], double S[3], double V[3][3])
{
    double A[3][3] = {
        { eps6[0],       0.5 * eps6[3], 0.5 * eps6[5] },
        { 0.5 * eps6[3], eps6[1],       0.5 * eps6[4] },
        { 0.5 * eps6[5], 0.5 * eps6[4], eps6[2]       }
    };
    eig3x3sym(A, S, V);
}
inline int eps1SpatialMult(const double eps6[6], const double S[3])
{
    double scale = 0.0;
    for (int a = 0; a < 6; ++a) scale += fabs(eps6[a]);
    scale += 1.0e-30;
    int m = 1;
    if (S[0] - S[1] < 1.0e-9 * scale) m = 2;
    if (S[0] - S[2] < 1.0e-9 * scale) m = 3;
    return m;
}

inline double eps1Spatial(const double eps6[6], double P1[6], double inj[6], double meas[6], int* mult)
{
    double S[3], V[3][3];
    strainEig3(eps6, S, V);
    int m = eps1SpatialMult(eps6, S);
    for (int a = 0; a < 6; ++a) { P1[a] = 0.0; inj[a] = 0.0; meas[a] = 0.0; }
    for (int k = 0; k < m; ++k) {
        double vx = V[0][k], vy = V[1][k], vz = V[2][k];
        double w[6] = { vx*vx, vy*vy, vz*vz, vx*vy, vy*vz, vx*vz };
        for (int a = 0; a < 6; ++a) {
            inj[a]  += w[a];
            P1[a]   += w[a] / m;
            meas[a] += (a < 3 ? 1.0 : 2.0) * w[a] / m;
        }
    }
    if (mult) *mult = m;
    double lam = S[0];
    if (m > 1) { lam = 0.0; for (int k = 0; k < m; ++k) lam += S[k]; lam /= m; }
    return macauley(lam);
}

// WP-263: eps1Spatial for a Dual strain. Values are the double function's. Derivatives from the
// eigenprojector of the top cluster C (the m eigenvectors sharing lambda_max), which is smooth
// while C stays separated from the rest even though its individual eigenvectors are not:
//     d lambda_C = (1/m) tr(dA~_CC),
//     d P_C = sum_{c in C, j not in C} dA~_jc / (lambda_c - lambda_j) (v_j v_c^T + v_c v_j^T),
// dA~ = V^T dA V, A the tensor strain (shear = gamma/2). inj = P_C, P1 = P_C/m, meas = the
// Voigt-weighted P_C/m, as in the double function. m = 3 (isotropic strain): P_C = I, d P_C = 0.
template <int N>
inline ladruno_dual::Dual<N> eps1Spatial(const ladruno_dual::Dual<N> eps6[6], ladruno_dual::Dual<N> P1[6],
                                         ladruno_dual::Dual<N> inj[6], ladruno_dual::Dual<N> meas[6], int* mult)
{
    typedef ladruno_dual::Dual<N> DN;
    static const int vr[6] = {0, 1, 2, 0, 1, 0}, vc[6] = {0, 1, 2, 1, 2, 2};
    double ev[6], P1v[6], injv[6], measv[6];
    for (int a = 0; a < 6; ++a) ev[a] = eps6[a].v;
    int m = 1;
    const double e1v = eps1Spatial(ev, P1v, injv, measv, &m);
    if (mult) *mult = m;
    double S[3], V[3][3];
    strainEig3(ev, S, V);
    DN e1(e1v);
    for (int a = 0; a < 6; ++a) { P1[a] = DN(P1v[a]); inj[a] = DN(injv[a]); meas[a] = DN(measv[a]); }
    for (int k = 0; k < N; ++k) {
        const double dA[3][3] = {
            { eps6[0].d[k],       0.5 * eps6[3].d[k], 0.5 * eps6[5].d[k] },
            { 0.5 * eps6[3].d[k], eps6[1].d[k],       0.5 * eps6[4].d[k] },
            { 0.5 * eps6[5].d[k], 0.5 * eps6[4].d[k], eps6[2].d[k]       }
        };
        double At[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double acc = 0.0;
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) acc += V[r][i] * dA[r][c] * V[c][j];
                At[i][j] = acc;
            }
        double dlam = 0.0;
        for (int c = 0; c < m; ++c) dlam += At[c][c];
        dlam /= m;
        e1.d[k] = (e1v > 0.0) ? dlam : 0.0;
        double dP[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
        for (int c = 0; c < m; ++c)
            for (int j = m; j < 3; ++j) {
                const double f = At[j][c] / (S[c] - S[j]);
                for (int r = 0; r < 3; ++r)
                    for (int q = 0; q < 3; ++q) dP[r][q] += f * (V[r][j] * V[q][c] + V[r][c] * V[q][j]);
            }
        for (int a = 0; a < 6; ++a) {
            const double dw = dP[vr[a]][vc[a]];
            inj[a].d[k]  = dw;
            P1[a].d[k]   = dw / m;
            meas[a].d[k] = (a < 3 ? 1.0 : 2.0) * dw / m;
        }
    }
    return e1;
}

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
    // --- Ladruno C2: cracked-Poisson latch (0/1), set once eps1 >= crackStrain with crackedNu on ---
    double nuCracked;
    // --- Ladruno WP-241: crackGate = damage tension-stiffening onset scale (latch) ---
    // < 0: the floor has not opened. >= 0: s = min(1, q_on / sigma_ts(e1_on)), frozen on the step
    // the damage-gated floor first opens; the floor is then s * sigma_ts(e1). Unused (stays -1)
    // under crackGate = strain.
    double tsScale;
};

inline void histZero(RCHist& h)
{
    for (int i = 0; i < 6; ++i) { h.stress_eff[i] = 0.0; h.strain[i] = 0.0; }
    h.xt = h.xc = 0.0; h.dt_bar = h.dc_bar = 0.0; h.beta = 1.0; h.eps1 = 0.0;
    h.cracked = 0.0; h.crackC = 1.0; h.crackS = 0.0; h.wmax = 0.0; h.betaSr = 1.0;
    h.tauCr = 0.0; h.gammaCr = 0.0;
    h.cracked2 = 0.0; h.slipCum = 0.0;
    h.xt_old = h.xc_old = 0.0; h.eps1_old = 0.0;
    h.nuCracked = 0.0;
    h.tsScale = -1.0;
}

template <class T>
inline T equivTensile(const T Si[3], double fcft, double Kc)
{
    if (Si[0] < STRESS_TOL) return T(0.0);
    return lublinerCriterion(Si[0], Si[1], Si[2], 1.0 / fcft, 1.0, 1.0, 1.0 / fcft, Kc);
}
template <class T>
inline T equivCompressive(const T Si[3], double fcft, double Kc, double k1)
{
    T s0 = Si[0] < 0.0 ? Si[0] : T(0.0);
    T s1 = Si[1] < 0.0 ? Si[1] : T(0.0);
    T s2 = Si[2] < 0.0 ? Si[2] : T(0.0);
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
// Ladruno WP-263: returnMap3D (signature unchanged, every caller unchanged) is the double entry
// point. tangentMode 3 (-dualTangent) routes a tangent request through returnMap3DDual: one
// Dual<6> pass of returnMap3DT, whose value part is the double run bit for bit and whose
// derivative part is d(sigma)/d(eps_eng) of the map executed -- the implicit map, or under
// IMPL-EX (do_implex) the extrapolated map with the committed history frozen. It therefore
// carries every term the mode-0 tangent omits: the projector rotation d(PT:S)/deps, the damage
// rate d(dt_bar, dc_bar)/deps, the tension-stiffening direction dp1/deps, dv_ci,max/de_n.
inline int returnMap3D(const Params& P, const double eps6[6], const RCHist& in,
                       double sig6[6], double Dtan6[6][6], RCHist& out,
                       bool do_tangent = true, double* residual = 0,
                       bool do_implex = false, double time_factor = 1.0);

template <class T>
inline int returnMap3DT(const Params& P, const T eps6[6], const RCHist& in,
                        T sig6[6], double Dtan6[6][6], RCHist& out,
                        bool do_tangent, double* residual,
                        bool do_implex, double time_factor)
{
    const bool implexExplicit = (P.implex && do_implex);
    const double E = P.E;
    // Ladruno (C2) -crackedNu: once the in-plane principal tensile strain reaches the cracking
    // strain, the elastic operator of the effective-stress predictor (and of every tangent built
    // on it) switches to nu = crackedNu (irreversible latch). The MCFT/DSFM convention is nu = 0
    // once cracked: with nu retained, the transverse tensile strain unloads the strut effective
    // stress (sig2_eff ~ eps2 + nu eps1). The effective stress is incremental, so the switch is
    // stress-continuous; the step that cracks uses the cracked operator for its whole increment.
    // Off (crackedNu < 0, the default) => C0(E, nu) exactly as before.
    // Ladruno WP-241: the add-on principal strain, membrane (default) or spatial, computed once.
    const bool spatial = (P.addonPlane == 1);
    T sp_P1[6], sp_inj[6], sp_meas[6]; int sp_mult = 1;
    const T e1_sp = spatial ? eps1Spatial(eps6, sp_P1, sp_inj, sp_meas, &sp_mult) : T(0.0);
    // crackGate = damage: crack onset = tensile threshold past the backbone peak
    const double xtPeak = (P.crackGate == 1) ? backbonePeakAbscissa(P.ht) : 0.0;
    double nuCr = in.nuCracked;
    if (P.crackedNu >= 0.0 && nuCr < 0.5) {
        if (P.crackGate == 1) {
            // WP-241: the COMMITTED threshold (the predictor needs nu before this step's
            // damage is known), so the latch lags the crack by one step; stress-continuous.
            if (in.xt > xtPeak) nuCr = 1.0;
        } else {
            T p1c[2]; int dgc = 0;
            T e1c = spatial ? e1_sp : eps1FromMembrane(eps6[0], eps6[1], eps6[3], p1c, &dgc);
            if (e1c >= P.crackStrain && P.crackStrain > 0.0)
                nuCr = 1.0;
        }
    }
    const double nuAct = (P.crackedNu >= 0.0 && nuCr >= 0.5) ? P.crackedNu : P.nu;
    double C0[6][6]; elasticTangent(E, nuAct, C0);

    // elastic predictor: SEFFn = SEFF_commit + C0:(En - En-1)
    T seff[6];
    for (int a = 0; a < 6; ++a) {
        T acc = in.stress_eff[a];
        for (int b = 0; b < 6; ++b) acc += C0[a][b] * (eps6[b] - in.strain[b]);
        seff[a] = acc;
    }

    DecompT<T> D;
    decompose(seff, P.cdf, D);

    // MCFT beta from in-plane membrane principal tensile strain (gxy = eps6[3] engineering)
    T p1[2] = { T(1.0), T(0.0) }; int degen = 0;
    T e1 = spatial ? e1_sp : eps1FromMembrane(eps6[0], eps6[1], eps6[3], p1, &degen);
    // IMPL-EX: freeze beta over the step from the extrapolated (committed) eps1, so the
    // returned tangent carries no beta cross-term (constant secant). eps1 itself is a
    // pure function of strain, so this is the only path-frozen quantity here.
    T e1_beta = e1;
    if (implexExplicit && P.betaOn) {
        e1_beta = in.eps1 + time_factor * (in.eps1 - in.eps1_old);
        if (e1_beta < 0.0) e1_beta = 0.0;
    }
    T beta = P.betaOn ? betaCompr(e1_beta, P.betaFloor, P.betaC) : T(1.0);
    double k1c  = (P.betaOn && P.lublinerTCReduced) ? 0.0 : 1.0;

    // committed hardening
    HPointT<T> pt = evaluateAt(P.ht, T(in.xt));
    HPointT<T> pc = evaluateAt(P.hc, T(in.xc));
    T xt_pl = plasticStrain(pt, E);
    T xc_pl = plasticStrain(pc, E);

    // rate coefficients (eta=0 -> rc1=0, rc2=1)
    double rc1 = (P.eta > 0.0) ? P.eta / (P.eta + 1.0) : 0.0;
    double rc2 = 1.0 - rc1;

    // /E converts the Lubliner STRESS measure to the strain abscissa (cpp:2509/2522)
    T xt_meas = equivTensile(D.Si, P.fcft_ratio, P.Kc) / E;
    T xc_meas = equivCompressive(D.Si, P.fcft_ratio, P.Kc, k1c) / E;
    T xt_trial = xt_meas + xt_pl;
    T xc_trial = xc_meas + xc_pl;
    // Ladruno WP-246: committed-surface consistency shift. The return rescales the tensile and the
    // compressive parts of the effective stress by different factors (1 - dt_plastic, 1 - dc_plastic,
    // below), and the Lubliner measure of a MIXED-sign stress is not homogeneous in either part. So
    // the committed effective stress can sit OUTSIDE its committed surface: measure + plastic strain
    // > committed threshold. Re-entering at the committed strain then advances the threshold again
    // and the stress relaxes on every hold step (inherited from ASDConcrete3D::compute). Measure the
    // committed overshoot and take it off the trial measure, so the threshold advances by the
    // INCREMENT of the measure from the committed state: a zero increment returns the committed
    // stress exactly. Only an overshoot is removed (a committed state inside its surface is untouched),
    // and the shift is a constant of the step (no tangent term).
    if (!implexExplicit && (in.xt > 0.0 || in.xc > 0.0)) {
        Decomp Dn;
        decompose(in.stress_eff, P.cdf, Dn);
        const double ot = (in.xt - value(xt_pl)) - equivTensile(Dn.Si, P.fcft_ratio, P.Kc) / E;
        const double oc = (in.xc - value(xc_pl)) - equivCompressive(Dn.Si, P.fcft_ratio, P.Kc, k1c) / E;
        if (ot < 0.0) xt_trial += ot;
        if (oc < 0.0) xc_trial += oc;
    }

    T new_xt = in.xt, new_xc = in.xc;
    if (implexExplicit) {
        // explicit extrapolation of the damage thresholds (frozen over the step)
        double xt_ext = in.xt + time_factor * (in.xt - in.xt_old);
        double xc_ext = in.xc + time_factor * (in.xc - in.xc_old);
        if (xt_ext < in.xt) xt_ext = in.xt;     // threshold is monotone non-decreasing
        if (xc_ext < in.xc) xc_ext = in.xc;
        pt = evaluateAt(P.ht, T(xt_ext)); new_xt = pt.x;
        pc = evaluateAt(P.hc, T(xc_ext)); new_xc = pc.x;
    } else {
        if (xt_trial > pt.x) { T xn = rc1 * pt.x + rc2 * xt_trial; pt = evaluateAt(P.ht, xn); new_xt = pt.x; }
        if (xc_trial > pc.x) { T xn = rc1 * pc.x + rc2 * xc_trial; pc = evaluateAt(P.hc, xn); new_xc = pc.x; }
    }

    T seff_eq_t = (pt.x - xt_pl) * E;
    T dt_plastic = seff_eq_t > 0.0 ? T(1.0 - pt.q / seff_eq_t) : T(0.0);
    T seff_eq_c = (pc.x - xc_pl) * E;
    T dc_plastic = seff_eq_c > 0.0 ? T(1.0 - pc.q / seff_eq_c) : T(0.0);

    // mix damages (cdf=0 -> R=0 -> identity)
    HPointT<T> auxp = evaluateAt(P.hc, T(pt.x * P.fcft_ratio));
    T aux_pl = plasticStrain(auxp, E);
    T aux_seff = (auxp.x - aux_pl) * E;
    T aux_d = aux_seff > 0.0 ? T(1.0 - auxp.q / aux_seff) : T(0.0);
    dc_plastic = 1.0 - (1.0 - dc_plastic) * (1.0 - D.R * aux_d);
    T pc_d = 1.0 - (1.0 - pc.d) * (1.0 - D.R * auxp.d);

    T dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
    T dc_bar = pc_d + dc_plastic - pc_d * dc_plastic;

    // nominal stress -- THE Phase-1 insertion: beta multiplies the assembled compressive cone
    // WP-263 (Dual only): the compressive measure has no STRESS_TOL gate (the tensile one does), so a
    // principal stress of ~-1e-12 "loads" the compression backbone at x < YTOL/E, where evaluateAt
    // clamps q to YTOL and dc_plastic = 1 - YTOL/(E x) is noise (<< 0 as x -> 0). The stress there
    // is ~YTOL, harmless, but its exact derivative is not: d22 ~ 0 on the PlaneStress sigma_33 root
    // (a micro-plateau of the map), which a condensation divides by. On that clamp the Dual takes
    // the elastic-side derivative beta dSC (the limit from x > YTOL/E, where dc = 0 exactly); the
    // value is untouched (fix has a zero value part). Quirk WP-263-rc-compression-ytol-micro-load.
    const bool microC = (value(pc.q) <= YTOL); (void)microC;
    for (int a = 0; a < 6; ++a) {
        T cterm = beta * (1.0 - dc_bar) * D.SC[a];
        if constexpr (ladruno_dual::is_dual<T>::value) {
            if (microC) { T fix = cterm - beta * D.SC[a]; ladruno_dual::zero_value(fix); cterm -= fix; }
        }
        sig6[a] = (1.0 - dt_bar) * D.ST[a] + cterm;
    }

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
    T e1_ts = e1;
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
    // Ladruno WP-241: the vectors live on an index list -- membrane {0,1,3} (the original
    // arithmetic, term for term), spatial all six (eps1Spatial's multiplicity-aware vectors).
    bool   ts_active = false;
    T ts_inj[6], ts_meas[6], ts_P1eps[6], ts_dsig = 0.0;
    for (int k = 0; k < 6; ++k) { ts_inj[k] = 0.0; ts_meas[k] = 0.0; ts_P1eps[k] = 0.0; }
    const int ts_idx_m[3] = { 0, 1, 3 }, ts_idx_s[6] = { 0, 1, 2, 3, 4, 5 };
    const int* ts_idx = spatial ? ts_idx_s : ts_idx_m;
    const int  ts_n   = spatial ? 6 : 3;
    const bool ts_gate = (P.crackGate == 1) ? (new_xt > xtPeak) : true;   // WP-241 crackGate
    // WP-241 crackGate = damage: the gate opens when tensile damage starts, and under tension-
    // compression (Lubliner) that is BELOW ft, where the bare principal stress q is already under
    // sigma_ts(e1). Injecting the full floor there is a stress JUMP (no root for a lateral balance).
    // So the floor is scaled by s = min(1, q / sigma_ts(e1)) taken on the step it first opens
    // (delta = 0 there: continuous), and frozen at commit (RCHist::tsScale). The gate stays shut
    // while q <= 0. crackGate = strain: s = 1 always (byte-identical).
    T tsScale = in.tsScale;
    if (P.tensStiffMode != 0 && e1_ts >= P.crackStrain && ts_gate) {
        if (spatial) {
            for (int k = 0; k < 6; ++k) { ts_inj[k] = sp_inj[k]; ts_meas[k] = sp_meas[k]; ts_P1eps[k] = sp_P1[k]; }
        } else if (degen) {
            ts_inj[0] = 1.0;  ts_inj[1] = 1.0;
            ts_meas[0] = 0.5; ts_meas[1] = 0.5;
            ts_P1eps[0] = 0.5; ts_P1eps[1] = 0.5;         // de1/d(exx)=de1/d(eyy)=0.5 (mean)
        } else {
            T a = p1[0]*p1[0], b = p1[1]*p1[1], ab = p1[0]*p1[1];
            ts_inj[0] = a;  ts_inj[1] = b;  ts_inj[3] = ab;
            ts_meas[0] = a; ts_meas[1] = b; ts_meas[3] = 2.0*ab;
            ts_P1eps[0] = a; ts_P1eps[1] = b; ts_P1eps[3] = ab;   // de1/deps (membrane)
        }
        T q = 0.0;
        for (int k = 0; k < ts_n; ++k) q += ts_meas[ts_idx[k]] * sig6[ts_idx[k]];
        T sig_ts = tensStiff(e1_ts, P.tensStiffMode, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, &ts_dsig);
        if (P.crackGate == 1) {
            if (tsScale < 0.0 && q > 0.0 && sig_ts > 0.0)      // onset step: open at q
                tsScale = (q < sig_ts) ? T(q / sig_ts) : T(1.0);
            const T s = (tsScale < 0.0) ? T(0.0) : tsScale;      // q <= 0: still shut
            sig_ts *= s; ts_dsig *= s;
        }
        T delta  = sig_ts - q;
        if (delta > 0.0) {
            for (int k = 0; k < ts_n; ++k) sig6[ts_idx[k]] += delta * ts_inj[ts_idx[k]];
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
    T crk_c = in.crackC, crk_s = in.crackS, wmax = in.wmax;
    double cracked = in.cracked;
    T tauCr = in.tauCr, gammaCr = in.gammaCr;        // 2b cyclic crack-shear state
    double cracked2 = in.cracked2;                   // 2b.2b X-crack + wear state
    T slipCum = in.slipCum;
    T betaSr = 1.0;
    bool   il_active = false, il_capped = false, il_cyclic = false;
    T il_m[3] = { T(0.0), T(0.0), T(0.0) };          // m_eps: stress back-rotation of a shear change
    T il_Gint = 0.0;                                 // interlock shear modulus (cyclic tangent)
    // shearRetMode 3 (rots) = rotating-coaxial: skip the fixed-crack shear entirely, so the
    // membrane shear stays smeared/spectral (no capture, no clip, no friction). il_active
    // stays false => the tangent carries no crack cross-term either.
    // WP-241: the fixed crack is an x-y object (2D normal, scalar slip); the spatial add-on
    // plane refuses -interlock in the parser and the kernel skips it here as a backstop.
    if (P.interlockOn && P.shearRetMode != 3 && !spatial) {
        if (cracked < 0.5 && e1 >= P.crackStrain && !degen) {   // freeze crack normal
            cracked = 1.0; crk_c = p1[0]; crk_s = p1[1];
        }
        bool just_captured = (in.cracked < 0.5 && cracked >= 0.5);
        if (cracked >= 0.5) {
            il_active = true;
            T exx = eps6[0], eyy = eps6[1], gxy = eps6[3];
            T c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
            T en = exx*c2 + eyy*s2 + gxy*cs;                    // strain normal to the FROZEN crack 1
            double sth = (P.crackSpacing > 0.0) ? P.crackSpacing : (P.lch > 0.0 ? P.lch : 1.0);
            // ---- Phase 2b.2b: second orthogonal crack (X-cracking) ----
            // Crack 2 normal = (-s,c) (orthogonal to crack 1), frozen when the strain normal to
            // THAT plane reaches eps_cr (FSAM CepsA2>=eunpA2 + orthogonality lock). The cap then
            // uses the GOVERNING (most-open) crack's opening, so the REVERSE shear direction is
            // interlock-capped too (2a/2b.1 single crack only capped the crack-1 direction).
            T en_gov = en;
            if (P.xcrackOn && P.interlockCyclic) {
                T c2n = -s, s2n = c;                            // crack-2 (orthogonal) normal
                T en2 = exx*c2n*c2n + eyy*s2n*s2n + gxy*c2n*s2n;
                if (cracked2 < 0.5 && en2 >= P.crackStrain) cracked2 = 1.0;
                if (cracked2 >= 0.5 && en2 > en_gov) en_gov = en2;   // most-open crack governs
            }
            T w = macauley(en_gov) * sth;
            if (w > wmax) wmax = w;                              // tracked max (diagnostic / 2a)
            T w_use = P.interlockCyclic ? w : wmax;             // reversible (cyclic) vs monotone (2a)
            T denom = 0.31 + 24.0 * w_use / (P.aggSize + 16.0);
            T vcimax = (denom > 0.0) ? T(0.18 * P.sqrtFc / denom) : T(0.18 * P.sqrtFc);
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
            T sxx = sig6[0], syy = sig6[1], sxy = sig6[3];
            T tau_sm = cs * (syy - sxx) + (c2 - s2) * sxy;      // smeared crack-plane shear
            T tau_ci;
            if (P.interlockCyclic) {
                il_cyclic = true;
                // crack-shear slip stiffness G_slip per -shearRetention mode (cap unchanged).
                double Gfull = 0.5 * E / (1.0 + nuAct);         // full elastic shear modulus G (C2: live nu)
                il_Gint = Gfull;                                 // mcft (default): full G
                if (P.shearRetMode == 1)      il_Gint = P.shearRetFactor * Gfull;  // const: mu*G
                else if (P.shearRetMode == 2) il_Gint = Gfull * (0.31 / denom);    // dsfm: width-degraded
                T g_nt = 2.0*(eyy - exx)*cs + gxy*(c2 - s2);
                // 2b.2b wear driver: accumulate the sliding distance (skip the capture step,
                // where in.gammaCr is still the uncracked 0 and would inject a spurious jump).
                slipCum = in.slipCum + (just_captured ? T(0.0) : T(fabs(g_nt - in.gammaCr)));
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
                betaSr = (vcimax > 1.0e-30) ? T(tau_ci / vcimax) : T(0.0);
            } else {
                tau_ci = tau_sm;                                // 2a: clip the smeared shear
                if (tau_ci >  vcimax) { tau_ci =  vcimax; il_capped = true; }
                if (tau_ci < -vcimax) { tau_ci = -vcimax; il_capped = true; }
                betaSr = (fabs(tau_sm) > 1.0e-30) ? T(tau_ci / tau_sm) : T(1.0);
            }
            T dtau = tau_ci - tau_sm;
            il_m[0] = -2.0*cs; il_m[1] = 2.0*cs; il_m[2] = c2 - s2;   // m_eps
            sig6[0] += dtau * il_m[0];
            sig6[1] += dtau * il_m[1];
            sig6[3] += dtau * il_m[2];
        }
    }

    // effective stress reassembly for next-step commit (beta NOT applied)
    for (int a = 0; a < 6; ++a)
        out.stress_eff[a] = value((1.0 - dt_plastic) * D.ST[a] + (1.0 - dc_plastic) * D.SC[a]);
    for (int a = 0; a < 6; ++a) out.strain[a] = value(eps6[a]);
    out.xt = value(new_xt); out.xc = value(new_xc);
    out.dt_bar = value(dt_bar); out.dc_bar = value(dc_bar); out.beta = value(beta); out.eps1 = value(e1);
    out.cracked = cracked; out.crackC = value(crk_c); out.crackS = value(crk_s);
    out.wmax = value(wmax); out.betaSr = value(betaSr);
    out.tauCr = value(tauCr); out.gammaCr = value(gammaCr);
    out.cracked2 = cracked2; out.slipCum = value(slipCum);
    // IMPL-EX (n-1) generation: carried through verbatim; the rolling n->n-1 is done
    // by the host material at commit (so a mid-step explicit pass never corrupts it).
    out.xt_old = in.xt_old; out.xc_old = in.xc_old; out.eps1_old = in.eps1_old;
    out.nuCracked = nuCr;   // Ladruno (C2)
    out.tsScale = value(tsScale);  // Ladruno WP-241 (damage-gated TS onset scale)

    // tangent. WP-263: built on the double path only. A Dual run (tangentMode 3, via
    // returnMap3DDual) reads the tangent off the derivative part of sig6 instead.
    if constexpr (!ladruno_dual::is_dual<T>::value) if (do_tangent) {
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
                double dbeta = dBetaCompr(e1, P.betaFloor, P.betaC);
                if (dbeta != 0.0) {
                    double P1[6] = { 0,0,0,0,0,0 };
                    if (spatial) { for (int k = 0; k < 6; ++k) P1[k] = sp_P1[k]; }   // WP-241: 3D dual
                    else if (degen) { P1[0] = 0.5; P1[1] = 0.5; }   // equibiaxial blend
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
                double row[6];                                              // d(q)/deps, all cols
                for (int c = 0; c < 6; ++c) {
                    double acc = 0.0;
                    for (int k = 0; k < ts_n; ++k) acc += ts_meas[ts_idx[k]] * Dtan6[ts_idx[k]][c];
                    row[c] = acc;
                }
                for (int a = 0; a < ts_n; ++a)
                    for (int c = 0; c < 6; ++c)
                        Dtan6[ts_idx[a]][c] += ts_inj[ts_idx[a]] * (ts_dsig * ts_P1eps[c] - row[c]);
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
                // Ladruno WP-241: the smeared crack shear depends on ALL six strains (eps_zz and
                // the transverse shears enter the in-plane stresses through C0 and the spectral
                // spine), so its removal spans all six columns -- the same fix the TS pinning
                // tangent got (quirk LEGACY-tension-stiffening-pinning-tangent-must-span). With
                // columns {0,1,3} only, a condensed view (eps_zz Newton in PlaneStress /
                // PlateFiber) kept the pinned-away eps_zz sensitivity. The friction stiffness
                // G m_eps (x) m_eps is a function of the in-plane slip only (columns {0,1,3}).
                double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
                const int vidx[3] = { 0, 1, 3 };
                double msig[3] = { -cs, cs, c2 - s2 };
                double row[6];
                for (int j = 0; j < 6; ++j) {
                    double acc = 0.0;
                    for (int k = 0; k < 3; ++k) acc += msig[k] * Dtan6[vidx[k]][j];
                    row[j] = acc;
                }
                double gcoef = (il_capped ? P.betaSrMin : 1.0) * il_Gint;
                double mfull[6] = { il_m[0], il_m[1], 0.0, il_m[2], 0.0, 0.0 };
                for (int a = 0; a < 3; ++a)
                    for (int j = 0; j < 6; ++j)
                        Dtan6[vidx[a]][j] += -il_m[a] * row[j] + gcoef * il_m[a] * mfull[j];
            } else if (il_capped) {
                double c = crk_c, s = crk_s, c2 = c*c, s2 = s*s, cs = c*s;
                const int vidx[3] = { 0, 1, 3 };
                double msig[3] = { -cs, cs, c2 - s2 };       // m_sigma (tensor shear projector)
                double row[6];                                // WP-241: all six columns (above)
                for (int j = 0; j < 6; ++j) {
                    double acc = 0.0;
                    for (int k = 0; k < 3; ++k) acc += msig[k] * Dtan6[vidx[k]][j];
                    row[j] = acc;
                }
                double fac = 1.0 - P.betaSrMin;
                for (int a = 0; a < 3; ++a)
                    for (int j = 0; j < 6; ++j)
                        Dtan6[vidx[a]][j] -= fac * il_m[a] * row[j];
            }
        }
    }

    if (residual) {
        double sn = 0.0; for (int a = 0; a < 6; ++a) sn += value(sig6[a]) * value(sig6[a]);
        *residual = sqrt(sn);
    }
    return STATUS_OK;
}

// Ladruno WP-263: the dual tangent. Seeds the six engineering strains, runs returnMap3DT once in
// Dual<6>, and returns sig6 / out from the value part (bit-identical to the double run) and
// Dtan6[a][b] = d sig_a / d eps_b from the derivative part. There is no local Newton in the map, so
// no implicit-function-theorem step is needed here (the condensing views -- sigma_33 = 0, the
// BeamFiber lateral balance -- condense this Dtan6 at their converged point, which is that step).
inline int returnMap3DDual(const Params& P, const double eps6[6], const RCHist& in,
                           double sig6[6], double Dtan6[6][6], RCHist& out,
                           double* residual, bool do_implex, double time_factor)
{
    typedef ladruno_dual::Dual<6> D6;
    D6 e[6], s[6];
    for (int i = 0; i < 6; ++i) e[i] = D6::variable(eps6[i], i);
    const int st = returnMap3DT<D6>(P, e, in, s, Dtan6, out, true, residual, do_implex, time_factor);
    for (int a = 0; a < 6; ++a) {
        sig6[a] = s[a].v;
        for (int b = 0; b < 6; ++b) Dtan6[a][b] = s[a].d[b];
    }
    return st;
}

inline int returnMap3D(const Params& P, const double eps6[6], const RCHist& in,
                       double sig6[6], double Dtan6[6][6], RCHist& out,
                       bool do_tangent, double* residual,
                       bool do_implex, double time_factor)
{
    if (do_tangent && P.tangentMode == 3)
        return returnMap3DDual(P, eps6, in, sig6, Dtan6, out, residual, do_implex, time_factor);
    return returnMap3DT<double>(P, eps6, in, sig6, Dtan6, out, do_tangent, residual, do_implex, time_factor);
}

// ===========================================================================
// Ladruno WP-243: the BeamFiber view.
//
// A beam fibre imposes the retained components {0, 3, 5} = (eps_xx, gamma_xy, gamma_xz). The
// lateral block {1, 2, 4} = (eps_yy, eps_zz, gamma_yz) is solved so that the concrete balances the
// smeared ties of its own direction:
//   r = ( sig_yy + rhoY f_s(eps_yy),  sig_zz + rhoZ f_s(eps_zz),  sig_yz ) = 0.
// A direction with rho = 0 is stress-free, so a fibre with no ties is the uniaxial-stress fibre
// with two shears. The retained tangent is the static condensation of the same Jacobian:
//   D_RR - D_RL (D_LL + T)^-1 D_LR,   T = diag(rhoY Et_y, rhoZ Et_z, 0).
// Root selection: Newton starts from the COMMITTED lateral strain (in.strain), with a backtracking
// line search on |r|, so it follows the branch the fibre is on through the cracking fold (where
// rho fyt < f_cr, the lateral residual is not monotonic in eps_yy).
// The caller supplies the view's Params (addonPlane = spatial). Returns STATUS_OK, or
// STATUS_NO_CONVERGE when the lateral balance is not met: the caller must refuse the trial, never
// report a stress (stock BeamFiberMaterial returns 0 there).
// ===========================================================================

// Bilinear kinematic-hardening bar, one direction. Committed plastic strain epN in, trial epNew and
// tangent Et out. b = 0 is elastic-perfectly plastic; b > 0 gives Et = b Es past yield.
inline double tieUpdate(double eps, double epN, double Es, double fy, double b, double& Et, double& epNew)
{
    const double H   = (b > 0.0) ? b * Es / (1.0 - b) : 0.0;   // linear kinematic modulus
    const double sTr = Es * (eps - epN);
    const double xi  = sTr - H * epN;                          // relative to the back stress H epN
    const double f   = fabs(xi) - fy;
    if (f <= 0.0) { Et = Es; epNew = epN; return sTr; }
    const double dg = f / (Es + H);
    epNew = epN + (xi > 0.0 ? dg : -dg);
    Et = Es * H / (Es + H);
    return Es * (eps - epNew);
}

inline bool solve3x3(const double A[3][3], const double b[3], double x[3])
{
    double M[3][4];
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) M[i][j] = A[i][j]; M[i][3] = b[i]; }
    for (int c = 0; c < 3; ++c) {
        int piv = c; double best = fabs(M[c][c]);
        for (int r = c + 1; r < 3; ++r) if (fabs(M[r][c]) > best) { best = fabs(M[r][c]); piv = r; }
        if (!(best > 1.0e-300)) return false;
        if (piv != c) for (int j = 0; j < 4; ++j) { double t = M[piv][j]; M[piv][j] = M[c][j]; M[c][j] = t; }
        for (int r = 0; r < 3; ++r) {
            if (r == c) continue;
            const double f = M[r][c] / M[c][c];
            for (int j = c; j < 4; ++j) M[r][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < 3; ++i) x[i] = M[i][3] / M[i][i];
    return true;
}

// Model C `sum` crack-check cap at a full strain (Params::crackCheck). Pure function of the strain
// and the committed tie plastic strains (the tie stresses are re-evaluated here), so its strain
// gradient, ties included, can be taken by central differences. Writes the principal data the
// caller injects along (as eps1Spatial). Returns the cap; *e1out = e1.
inline double crackCheckCap(const Params& P, const double eps6[6], const double tieEpN[2],
                            double P1[6], double inj[6], double meas[6], double* e1out)
{
    int m = 1;
    const double e1 = eps1Spatial(eps6, P1, inj, meas, &m);
    if (e1out) *e1out = e1;
    double s2 = P1[0];                                   // n1_x^2 (multiplicity-averaged)
    if (s2 < 0.0) s2 = 0.0;
    if (s2 > 1.0) s2 = 1.0;
    const double st = sqrt(s2), ct = sqrt(1.0 - s2);
    const double sth = 1.0 / (st / P.ccSmx + ct / P.ccSmt + 1.0e-12);
    const double w = e1 * sth;
    const double vci = 0.18 * P.sqrtFc / (0.31 + 24.0 * w / (P.aggSize + 16.0));
    double th = atan2(st, ct);
    const double thLo = 1.0e-3, thHi = 1.5707963267948966 - 1.0e-3;   // pi/2 (no M_PI on MSVC)
    if (th < thLo) th = thLo;
    if (th > thHi) th = thHi;
    double cap = vci * tan(th);
    const double rho[2] = {P.tieRhoY, P.tieRhoZ};
    const double ty = P1[1], tz = P1[2], tt = ty + tz;
    const double wt[2] = {tt > 1.0e-30 ? ty / tt : 0.5, tt > 1.0e-30 ? tz / tt : 0.5};
    for (int d = 0; d < 2; ++d) {
        if (!(rho[d] > 0.0)) continue;
        double Et, epNew;
        const double fs = tieUpdate(eps6[1 + d], tieEpN[d], P.tieEs, P.tieFy, P.tieB, Et, epNew);
        const double res = P.tieFy - fs;
        if (res > 0.0) cap += rho[d] * res * wt[d];
    }
    return cap;
}

// Apply the cap to (sig6, Dtan6) in place: when q = meas.sig > cap, inject (cap - q) along inj
// (the rank-1 machinery of the TS floor, pushing down instead of up). Tangent:
// Dtan += inj (x) (d cap/d eps - meas.Dtan), d cap/d eps by central differences (h = 1e-9),
// the direction derivative of n1 omitted as in the TS floor. Returns true when the cap binds.
inline bool applyCrackCheck(const Params& P, const double eps6[6], const double tieEpN[2],
                            double sig6[6], double Dtan6[6][6])
{
    double P1[6], inj[6], meas[6], e1 = 0.0;
    const double cap = crackCheckCap(P, eps6, tieEpN, P1, inj, meas, &e1);
    if (!(e1 > 0.0)) return false;
    double q = 0.0;
    for (int a = 0; a < 6; ++a) q += meas[a] * sig6[a];
    if (q <= cap) return false;
    double g[6];
    for (int k = 0; k < 6; ++k) {
        const double h = 1.0e-9;
        double ep[6], em[6], d1[6], d2[6], d3[6];
        for (int a = 0; a < 6; ++a) { ep[a] = eps6[a]; em[a] = eps6[a]; }
        ep[k] += h; em[k] -= h;
        g[k] = (crackCheckCap(P, ep, tieEpN, d1, d2, d3, 0) - crackCheckCap(P, em, tieEpN, d1, d2, d3, 0)) / (2.0 * h);
    }
    double mD[6];
    for (int k = 0; k < 6; ++k) { mD[k] = 0.0; for (int a = 0; a < 6; ++a) mD[k] += meas[a] * Dtan6[a][k]; }
    const double delta = cap - q;
    for (int a = 0; a < 6; ++a) {
        sig6[a] += delta * inj[a];
        for (int k = 0; k < 6; ++k) Dtan6[a][k] += inj[a] * (g[k] - mD[k]);
    }
    return true;
}

// eps6: retained {0,3,5} set by the caller; lateral {1,2,4} overwritten with the converged values.
// tieEpN: committed tie plastic strains (y, z). tieEp/tieSig: trial plastic strains and bar stresses.
// iters (optional): Newton iterations used. rootScan (optional, WP-247): true when the Newton stalled and
// the Model C root scan found the balanced state. noRootAccepted (optional, WP-247): true when
// -lateralNoRoot modelC accepted a state whose r_y is not balanced.
inline int driveBeamFiber(const Params& P, double eps6[6], const RCHist& in, const double tieEpN[2],
                          double sig6[6], double Dtan6[6][6], RCHist& out,
                          double tieEp[2], double tieSig[2], int* iters = 0,
                          bool do_implex = false, double time_factor = 1.0, bool* crackCheckBinds = 0,
                          bool* rootScan = 0, bool* noRootAccepted = 0)
{
    static const int L[3] = {1, 2, 4};
    const double rho[2] = {P.tieRhoY, P.tieRhoZ};
    double Et[2] = {0.0, 0.0};
    bool ccActive = false;
    double r[3], J[3][3];
    for (int a = 0; a < 3; ++a) eps6[L[a]] = in.strain[L[a]];

    auto evaluate = [&]() -> double {
        returnMap3D(P, eps6, in, sig6, Dtan6, out, true, 0, do_implex, time_factor);
        for (int d = 0; d < 2; ++d) {
            Et[d] = 0.0; tieEp[d] = tieEpN[d]; tieSig[d] = 0.0;
            if (rho[d] > 0.0)
                tieSig[d] = tieUpdate(eps6[1 + d], tieEpN[d], P.tieEs, P.tieFy, P.tieB, Et[d], tieEp[d]);
        }
        // Model C crack check: its tie dependence sits in the eps_yy/eps_zz columns of Dtan6
        ccActive = (P.crackCheck == 1) && applyCrackCheck(P, eps6, tieEpN, sig6, Dtan6);
        r[0] = sig6[1] + rho[0] * tieSig[0];
        r[1] = sig6[2] + rho[1] * tieSig[1];
        r[2] = sig6[4];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) J[i][j] = Dtan6[L[i]][L[j]];
        J[0][0] += rho[0] * Et[0];
        J[1][1] += rho[1] * Et[1];
        return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    };
    auto stressScale = [&]() -> double {
        double s = 0.0; for (int a = 0; a < 6; ++a) s += sig6[a] * sig6[a];
        s = sqrt(s);
        return (s > P.ftPeak) ? s : (P.ftPeak > 0.0 ? P.ftPeak : 1.0);
    };

    // Ladruno WP-247: the lateral unknowns live in a physical box, Model C's bracket [-0.01, 0.15] for the
    // normals (and +-0.15 for gamma_yz). A fibre fully cracked in a lateral direction has a nearly flat
    // sigma there, so without the box a near-singular Jacobian sent Newton to eps_yy ~ 3e4 with a small
    // residual (a cover fibre of a high-axial-load column, Sezen-Moehle 2CHD12), and the state never recovered.
    auto box = [](int a, double x) -> double {
        const double lo = (a == 2) ? -0.15 : -0.01, hi = 0.15;
        return x < lo ? lo : (x > hi ? hi : x);
    };

    double rn = evaluate();
    bool conv = false;
    int it = 0;
    // Ladruno WP-247: a Newton that cycles (on the falling branch of a lightly tied or cracked fibre) is left
    // as soon as its best residual has not halved in kStall iterations, instead of running all 50: past the
    // peak of a force-based column these failing Newtons were a third of the fibre cost. The balance tests
    // are unchanged, and a stalled Newton goes to the root scan as before.
    static const int kStall = 8;
    double best = rn; int sinceBest = 0;
    // Ladruno WP-247: Broyden-corrected Jacobian. The material's default tangent is the damaged secant
    // W_B C0, stiffer than the true slope on a softening branch, so plain Newton on the lateral balance
    // contracted by a constant factor (~0.6 per iteration, ~30 iterations to the 1e-10 test). Each accepted
    // step corrects Jb with the residual change it produced (Broyden's rank-one update); a step that
    // fails the line search resets Jb to the fresh tangent. The balance tests, and the tangent returned
    // to the caller (condensed from the material's own Dtan6 and J), are unchanged.
    double Jb[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Jb[i][j] = J[i][j];
    for (; it < 50; ++it) {
        if (rn <= 1.0e-10 * stressScale()) { conv = true; break; }
        if (rn < 0.5 * best) { best = rn; sinceBest = 0; }
        else if (++sinceBest >= kStall) break;
        double dx[3];
        if (!solve3x3(Jb, r, dx)) break;
        const double base[3] = {eps6[1], eps6[2], eps6[4]};
        const double r0[3] = {r[0], r[1], r[2]};
        double lam = 1.0, rnew = rn;
        for (int ls = 0; ls < 8; ++ls) {               // backtracking on |r|; the last try is kept
            for (int a = 0; a < 3; ++a) eps6[L[a]] = box(a, base[a] - lam * dx[a]);
            rnew = evaluate();
            if (rnew < rn) break;
            lam *= 0.5;
        }
        if (rnew < rn) {
            double s[3], ss = 0.0;
            for (int a = 0; a < 3; ++a) { s[a] = eps6[L[a]] - base[a]; ss += s[a] * s[a]; }
            if (ss > 0.0) {
                for (int i = 0; i < 3; ++i) {
                    double Js = 0.0;
                    for (int j = 0; j < 3; ++j) Js += Jb[i][j] * s[j];
                    const double c = (r[i] - r0[i] - Js) / ss;
                    for (int j = 0; j < 3; ++j) Jb[i][j] += c * s[j];
                }
            }
        }
        else {
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Jb[i][j] = J[i][j];
        }
        rn = rnew;
    }
    if (!conv) conv = (rn <= 1.0e-6 * stressScale());  // a stalled but balanced point is accepted

    // Ladruno WP-247: when Newton stalls, Model C's root rule. Where rho fyt < f_cr the lateral residual
    // r_y(eps_yy) = sigma_yy + rho f_s is not monotone past the tension peak (the concrete softens faster
    // than the ties stiffen), and a full Newton step from the uncracked side lands on the falling branch and
    // cycles there, although a root sits next to the committed strain. The rule: along one lateral normal
    // (y, then z), sample the committed value plus geometric offsets +-[1e-8, 1e-1], take the INCREASING
    // sign change of that direction's residual nearest the committed value (the stable root on the branch
    // continuous with the committed state), and bisect it; the other two lateral unknowns are solved by a
    // 2x2 Newton at every sample. The point is accepted on the same balance test as the stalled-Newton
    // branch above, so a bracket around a jump (no root) still refuses. It runs only when Newton failed,
    // so every state Newton reaches is kept (to the balance tolerance).
    bool scanned = false, noRoot = false;
    // Cost: the walk goes outward from the committed value and stops at the first bracket, so a root next
    // to the committed strain costs a few samples; the inner solve is capped (12 iterations, 4 halvings)
    // and warm-started from the previous sample, and the whole scan has a budget of kScanBudget return
    // maps, past which the trial is refused as before. Without the caps, a force-based element whose
    // global Newton was diverging scanned every refused fibre in full on every element iteration, and one
    // analyze step ran for minutes.
    static const int kScanBudget = 3000;
    int nEval = 0;
    bool innerOk = false, warmOk = false;
    double warm[2] = {0.0, 0.0};
    // balance the other two lateral unknowns (o1 and gamma_yz) with eps6[L[ax]] held, from the current eps6
    auto innerSolve = [&](int o1) {
        auto q2 = [&]() { return sqrt(r[o1] * r[o1] + r[2] * r[2]); };
        double qn = q2();
        for (int k = 0; k < 12 && qn > 1.0e-10 * stressScale(); ++k) {
            const double a = J[o1][o1], b = J[o1][2], c = J[2][o1], d = J[2][2];
            const double det = a * d - b * c;
            if (!(fabs(det) > 1.0e-300)) break;
            const double d1 = ( d * r[o1] - b * r[2]) / det;
            const double d2 = (-c * r[o1] + a * r[2]) / det;
            const double b1 = eps6[L[o1]], b2 = eps6[L[2]];
            double lam = 1.0, qt = qn;
            for (int ls = 0; ls < 4; ++ls) {
                eps6[L[o1]] = box(o1, b1 - lam * d1); eps6[L[2]] = box(2, b2 - lam * d2);
                evaluate(); ++nEval; qt = q2();
                if (qt < qn) break;
                lam *= 0.5;
            }
            if (!(qt < qn)) {
                // no descent: a flat direction (a cover fibre cracked in z has sigma_zz ~ 0). Further
                // iterations cannot help, so stop here; if the last point was already balanced, go back to it.
                if (qn <= 1.0e-8 * stressScale()) {
                    eps6[L[o1]] = b1; eps6[L[2]] = b2;
                    evaluate(); ++nEval;
                    qn = q2();
                }
                else qn = qt;
                break;
            }
            qn = qt;
        }
        innerOk = (qn <= 1.0e-8 * stressScale());
    };
    // r[ax] at eps6[L[ax]] = x after balancing the other two: warm start from the previous balanced
    // sample, then (if that fails) from the committed state
    auto innerAt = [&](int ax, double x) -> double {
        const int o1 = (ax == 0) ? 1 : 0;
        for (int pass = (warmOk ? 0 : 1); pass < 2; ++pass) {
            eps6[L[ax]] = box(ax, x);
            eps6[L[o1]] = (pass == 0) ? warm[0] : box(o1, in.strain[L[o1]]);
            eps6[L[2]]  = (pass == 0) ? warm[1] : box(2, in.strain[L[2]]);
            evaluate(); ++nEval;
            innerSolve(o1);
            if (innerOk) break;
        }
        warmOk = innerOk;
        if (innerOk) { warm[0] = eps6[L[o1]]; warm[1] = eps6[L[2]]; }
        return r[ax];
    };
    bool foundY = false; double loY = 0.0, hiY = 0.0;
    if (!conv) {
        for (int ax = 0; ax < 2 && !conv && nEval < kScanBudget; ++ax) {
            warmOk = false;
            // walk outward: interval k spans offsets (off[k-1], off[k]] on each side (off[-1] = 0, the committed
            // value itself); the first interval holding an increasing sign change is the nearest one. Each
            // side warm-starts from its own previous sample.
            const double x0 = in.strain[L[ax]];
            double lo = 0.0, hi = 0.0; bool found = false;
            const double g0 = innerAt(ax, x0); const bool ok0 = innerOk;
            double gR = g0, gL = g0; bool okR = ok0, okL = ok0; double xR = x0, xL = x0;
            double wR[2] = {warm[0], warm[1]}, wL[2] = {warm[0], warm[1]};
            bool wOkR = warmOk, wOkL = warmOk;
            for (int k = 0; k < 41 && !found && nEval < kScanBudget; ++k) {
                const double off = 1.0e-8 * pow(10.0, 7.0 * k / 40.0);
                warm[0] = wR[0]; warm[1] = wR[1]; warmOk = wOkR;
                const double xr = x0 + off, gr = innerAt(ax, xr); const bool okr = innerOk;
                wR[0] = warm[0]; wR[1] = warm[1]; wOkR = warmOk;
                if (okR && okr && gR <= 0.0 && gr >= 0.0) { lo = xR; hi = xr; found = true; break; }
                xR = xr; gR = gr; okR = okr;
                warm[0] = wL[0]; warm[1] = wL[1]; warmOk = wOkL;
                const double xl = x0 - off, gl = innerAt(ax, xl); const bool okl = innerOk;
                wL[0] = warm[0]; wL[1] = warm[1]; wOkL = warmOk;
                if (okl && okL && gl <= 0.0 && gL >= 0.0) { lo = xl; hi = xL; found = true; break; }
                xL = xl; gL = gl; okL = okl;
            }
            if (!found) continue;
            if (ax == 0) { foundY = true; loY = lo; hiY = hi; }
            for (int k = 0; k < 60; ++k) {
                const double m = 0.5 * (lo + hi);
                if (innerAt(ax, m) <= 0.0) lo = m; else hi = m;
            }
            innerAt(ax, 0.5 * (lo + hi));                   // leaves eps6, sig6, Dtan6, J, out at the root
            rn = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            conv = innerOk && (rn <= 1.0e-6 * stressScale());
            scanned = conv;
        }
        // Ladruno WP-247, -lateralNoRoot modelC (opt-in): Model C never refuses a strip. Without a balanced
        // root it bisects r_y over the nearest increasing y-bracket, or over its full range [-0.004, 0.15]
        // when there is none, and keeps the end point: the jump, or a range end where sigma_yy has one sign
        // throughout (a spalled cover strip without ties, after a shear failure). eps_zz and gamma_yz are
        // still balanced, to kNoRootInner x the stress scale: Model C's strip is a membrane with no zz or yz
        // equation, and a cracked fibre's sigma_yz can be flat (no yz shear stiffness left), so a residual of
        // a few 1e-6 MPa there, which the inner solve cannot remove, does not refuse the strip. r_y is left
        // unbalanced.
        static const double kNoRootInner = 1.0e-4;
        if (!conv && P.lateralNoRoot == 1) {
            warmOk = false;
            double lo = foundY ? loY : -0.004, hi = foundY ? hiY : 0.15;
            for (int k = 0; k < 60; ++k) {
                const double m = 0.5 * (lo + hi);
                if (innerAt(0, m) <= 0.0) lo = m; else hi = m;
            }
            innerAt(0, 0.5 * (lo + hi));
            rn = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            const double qInner = sqrt(r[1] * r[1] + r[2] * r[2]);
            conv = innerOk || (qInner <= kNoRootInner * stressScale());
            noRoot = conv;
        }
    }
    if (iters) *iters = it;
    if (rootScan) *rootScan = scanned;
    if (noRootAccepted) *noRootAccepted = noRoot;
    if (crackCheckBinds) *crackCheckBinds = ccActive;
    if (!conv) return STATUS_NO_CONVERGE;

    // condensed tangent at the converged state (J and Dtan6 belong to the last evaluate())
    double T[3][6];
    for (int k = 0; k < 6; ++k) {
        const double rhs[3] = {Dtan6[L[0]][k], Dtan6[L[1]][k], Dtan6[L[2]][k]};
        double c[3];
        if (!solve3x3(J, rhs, c)) {
            // a -lateralNoRoot modelC state can sit on a flat branch (J singular): a small elastic
            // regularization of the lateral block keeps the condensed tangent finite
            if (!noRoot) return STATUS_NO_CONVERGE;
            double Jr[3][3];
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Jr[i][j] = J[i][j] + (i == j ? 1.0e-6 * P.E : 0.0);
            if (!solve3x3(Jr, rhs, c)) return STATUS_NO_CONVERGE;
        }
        for (int i = 0; i < 3; ++i) T[i][k] = c[i];
    }
    double Dc[6][6];
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
        double v = Dtan6[i][j];
        for (int a = 0; a < 3; ++a) v -= Dtan6[i][L[a]] * T[a][j];
        Dc[i][j] = v;
    }
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Dtan6[i][j] = Dc[i][j];
    return STATUS_OK;
}

// Elastic BeamFiber tangent (getInitialTangent, revertToStart): the same condensation on C0 with
// the elastic tie stiffness rho Es.
inline void beamFiberInitialTangent(const Params& P, double D[6][6])
{
    static const int L[3] = {1, 2, 4};
    elasticTangent(P.E, P.nu, D);
    double J[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) J[i][j] = D[L[i]][L[j]];
    J[0][0] += P.tieRhoY * P.tieEs;
    J[1][1] += P.tieRhoZ * P.tieEs;
    double T[3][6];
    for (int k = 0; k < 6; ++k) {
        const double rhs[3] = {D[L[0]][k], D[L[1]][k], D[L[2]][k]};
        double c[3] = {0.0, 0.0, 0.0};
        solve3x3(J, rhs, c);
        for (int i = 0; i < 3; ++i) T[i][k] = c[i];
    }
    double Dc[6][6];
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
        double v = D[i][j];
        for (int a = 0; a < 3; ++a) v -= D[i][L[a]] * T[a][j];
        Dc[i][j] = v;
    }
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) D[i][j] = Dc[i][j];
}

} // namespace ladruno_rc_kernel

#endif // LadrunoRCKernel_h
