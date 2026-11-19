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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 11/2026
//
// LadrunoMPMCore — the PURE numerical core of the explicit MPM engine (WP-232,
// mode F; ADR 232 §5.2 steps 1 and 3–7). Plain C++17 + the header-only hypo
// kernel, NO OpenSees dependency, so it is unit-tested standalone with g++
// against tests/mpm_cases.txt (tests/mpm_engine_check.cpp) without building
// OpenSees — the LadrunoHypoKernel.h / LadrunoUPKernel.h idiom. The OpenSees
// host (LadrunoMPMPattern, a LoadPattern) owns the Domain side: the load factor
// (step 2), Domain::commit and recorders (step 8), the profiler rows (step 9)
// and the NDMaterial copies, which it hands in through MaterialPoints.
//
// The numerics follow tests/mpm_interface.md (frozen by WP-232b; this file was
// written from it and never from the oracle tests/mpm_reference.py):
//   §2 grid and node numbering; §4 uGIMP and linear hat; §5 lattice seeding;
//   §6 one MUSL step with PIC/FLIP, mass cutoff, local damping and the
//   critical time step; §7 grid-face conditions (fix, frictional wall);
//   §8 the incremental kernel entry; §9 sig_hat = sig0_hat + material(eps_hat).
//
// Conventions: tension positive; row-major F[ndm*i+j] = dx_i/dX_j; Voigt
// [xx yy zz xy yz zx] (3D) or [xx yy xy] (plane strain), engineering shear in
// strain, tensor shear in stress. Masses and volumes in 2D are full values
// (thickness 1).

#ifndef LadrunoMPMCore_h
#define LadrunoMPMCore_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../../element/solidTransformation/LadrunoHypoKernel.h"

// ADR-87 D2 mutation family MPM. The OpenSees build sets it through
// SRC/Ladruno_mutation.h's CMake plumbing (-DLADRUNO_MUTATE_MPM=<mode>); the
// default 0 is LADRUNO_MUT_NONE. Mode codes mirror Ladruno_mutation.h (1 ZERO,
// 2 SCALE); LadrunoMPMPattern.cpp static_asserts that they agree.
#ifndef LADRUNO_MUTATE_MPM
#  define LADRUNO_MUTATE_MPM 0
#endif

namespace ladruno_mpm {

const int    MUT_ZERO = 1;         // ZERO : the P2G internal force is dropped
const int    MUT_SCALE = 2;        // SCALE: the G2P position increment is x1.5
const double MUT_SCALE_FACTOR = 1.5;

// WP-232i (tests/mpm_bspline_interface.md B4): node-centred uniform B-splines,
// degree 2 and 3, on the same vertex lattice and ghost layer.
enum Basis { BASIS_GIMP = 0, BASIS_LINEAR = 1, BASIS_BSPLINE2 = 2, BASIS_BSPLINE3 = 3 };

inline bool isBSpline(Basis b) { return b == BASIS_BSPLINE2 || b == BASIS_BSPLINE3; }

inline const char *basisName(Basis b)
{
  switch (b) {
    case BASIS_GIMP: return "gimp";
    case BASIS_LINEAR: return "linear";
    case BASIS_BSPLINE2: return "bspline2";
    case BASIS_BSPLINE3: return "bspline3";
  }
  return "?";
}

// name -> basis; false for an unknown name
inline bool basisFromName(const std::string &s, Basis &b)
{
  if (s == "gimp") { b = BASIS_GIMP; return true; }
  if (s == "linear") { b = BASIS_LINEAR; return true; }
  if (s == "bspline2") { b = BASIS_BSPLINE2; return true; }
  if (s == "bspline3") { b = BASIS_BSPLINE3; return true; }
  return false;
}

// --------------------------------------------------------------------------- //
//  Basis (interface §4)                                                        //
// --------------------------------------------------------------------------- //

// Uniform B-spline 1D factor (tests/mpm_bspline_interface.md B1), r = x_p - x_i,
// degree 2 (support |r| < 3h/2) or 3 (|r| < 2h). dS = d/dx_p.
inline void bsplineFactor(int degree, double r, double h, double &S, double &dS)
{
  const double xi = r / h;
  const double a = std::fabs(xi);
  const double sg = (xi > 0.0) ? 1.0 : ((xi < 0.0) ? -1.0 : 0.0);
  if (degree == 2) {
    if (a < 0.5) {
      S  = 0.75 - xi * xi;
      dS = -2.0 * xi / h;
    } else if (a < 1.5) {
      const double t = 1.5 - a;
      S  = 0.5 * t * t;
      dS = -sg * t / h;
    } else {
      S = 0.0;
      dS = 0.0;
    }
    return;
  }
  if (a < 1.0) {
    S  = 2.0 / 3.0 - xi * xi + 0.5 * a * a * a;
    dS = (-2.0 * xi + 1.5 * xi * a) / h;
  } else if (a < 2.0) {
    const double t = 2.0 - a;
    S  = t * t * t / 6.0;
    dS = -sg * t * t / (2.0 * h);
  } else {
    S = 0.0;
    dS = 0.0;
  }
}

// uGIMP 1D factor, r = x_p - x_i, 0 < lp <= h/2 (interface §4.1).
inline void gimpFactor(double r, double h, double lp, double &S, double &dS)
{
  const double a = std::fabs(r);
  const double sg = (r > 0.0) ? 1.0 : ((r < 0.0) ? -1.0 : 0.0);
  if (a < lp) {
    S  = 1.0 - (r * r + lp * lp) / (2.0 * h * lp);
    dS = -r / (h * lp);
  } else if (a < h - lp) {
    S  = 1.0 - a / h;
    dS = -sg / h;
  } else if (a < h + lp) {
    const double t = h + lp - a;
    S  = t * t / (4.0 * h * lp);
    dS = -sg * t / (2.0 * h * lp);
  } else {
    S = 0.0;
    dS = 0.0;
  }
}

// --------------------------------------------------------------------------- //
//  Grid (interface §2)                                                         //
// --------------------------------------------------------------------------- //

// GHOST LAYERS (interface Amendment 5, replacing A2's truncation). The engine
// grid carries `ghost` node layers beyond every face, enough to hold every
// support node of an admissible particle: centre in the real grid and
// lp <= h/2 give |r| < h + lp <= 1.5 h, so one layer suffices. Nothing is
// dropped, so the partition of unity and Σ∇S = 0 hold for every admissible
// particle (truncation broke Σ∇S along a face: spurious L_tn and drag). A face's
// grid condition acts on its face nodes and on its ghost layers alike; a face
// with no condition just has free ghost nodes.
//
// Numbering: the engine index (extId) runs over the extended lattice
// i in [-ghost, n+ghost], x fastest; with ghost = 1 it IS the interface §2 id of
// Amendment 5, which the hooks report (publicId).
struct Grid {
  int    ndm = 0;
  double x0[3] = {0.0, 0.0, 0.0};
  double h = 0.0;
  int    n[3] = {1, 1, 1};
  int    ghost = 1;

  int  extCount(int d) const { return n[d] + 1 + 2 * ghost; }
  long numNodes() const       // engine (extended) node count
  {
    long c = 1;
    for (int d = 0; d < ndm; d++) c *= (long)extCount(d);
    return c;
  }
  // engine index of a node, -ghost <= i <= n + ghost per axis
  long extId(const int ijk[3]) const
  {
    const long a = ijk[0] + ghost, b = ijk[1] + ghost;
    if (ndm == 2) return a + (long)extCount(0) * b;
    const long c = ijk[2] + ghost;
    return a + (long)extCount(0) * (b + (long)extCount(1) * c);
  }
  void extIJK(long id, int ijk[3]) const
  {
    ijk[2] = 0;
    ijk[0] = (int)(id % extCount(0)) - ghost;
    long rest = id / extCount(0);
    if (ndm == 2) { ijk[1] = (int)rest - ghost; return; }
    ijk[1] = (int)(rest % extCount(1)) - ghost;
    ijk[2] = (int)(rest / extCount(1)) - ghost;
  }
  bool isReal(const int ijk[3]) const
  {
    for (int d = 0; d < ndm; d++)
      if (ijk[d] < 0 || ijk[d] > n[d]) return false;
    return true;
  }
  // The id a hook reports. Amendment 5 numbers the nodes WITH their ghost layer:
  // id = (i+1) + (nx+3)((j+1) + (ny+3)(k+1)), which is exactly the engine index
  // for ghost = 1 (asserted by the g++ core check).
  long publicId(long ext) const { return ext; }
  // cell index of coordinate x on axis d: floor((x - x0)/h), clamped to n-1
  int cellOf(int d, double x) const
  {
    int c = (int)std::floor((x - x0[d]) / h);
    if (c > n[d] - 1) c = n[d] - 1;
    if (c < 0) c = 0;
    return c;
  }
};

// One axis of a particle's support: node indices with S, dS/dx.
struct AxisSupport {
  int    count = 0;
  int    idx[4];
  double S[4];
  double dS[4];
};

// Admissibility on one axis (interface §4.1 as amended by A2 and A5, §4.2): the
// particle CENTRE inside the closed real grid box, and for GIMP 0 < lp <= h/2.
// Its support nodes beyond a face live in the ghost layers (A5): none dropped.
inline bool admissibleAxis(const Grid &g, Basis b, int d, double x, double lp)
{
  const double lo = g.x0[d], hi = g.x0[d] + g.n[d] * g.h;
  if (b == BASIS_LINEAR || isBSpline(b))      // B3: the centre rule; lp not used
    return (lo <= x) && (x <= hi);
  return (lp > 0.0) && (lp <= 0.5 * g.h) && (lo <= x) && (x <= hi);
}

inline void axisSupport(const Grid &g, Basis b, int d, double x, double lp, AxisSupport &a)
{
  a.count = 0;
  const double h = g.h;
  if (b == BASIS_LINEAR) {
    const int c = g.cellOf(d, x);
    const double xi = (x - (g.x0[d] + c * h)) / h;
    a.idx[0] = c;     a.S[0] = 1.0 - xi; a.dS[0] = -1.0 / h;
    a.idx[1] = c + 1; a.S[1] = xi;       a.dS[1] = 1.0 / h;
    a.count = 2;
    return;
  }
  const int c = (int)std::floor((x - g.x0[d]) / h);
  const int deg = (b == BASIS_BSPLINE2) ? 2 : ((b == BASIS_BSPLINE3) ? 3 : 0);
  for (int i = c - 1; i <= c + 2; i++) {
    if (i < -g.ghost || i > g.n[d] + g.ghost) continue;   // never hit when admissible
    double S, dS;
    if (deg) {
      // B2: the same scan. r is measured from the NEARER face, so a node and its
      // image across that face get |r| ordered exactly in floating point (the
      // folded S' = S_image - S_ghost of K4 is then never negative by round-off).
      const double hi = g.x0[d] + g.n[d] * h;
      const double r = (2.0 * (x - g.x0[d]) <= g.n[d] * h) ? (x - g.x0[d]) - i * h
                                                           : (x - hi) - (i - g.n[d]) * h;
      bsplineFactor(deg, r, h, S, dS);
    }
    else gimpFactor(x - (g.x0[d] + i * h), h, lp, S, dS);
    if (S == 0.0 && dS == 0.0) continue;
    a.idx[a.count] = i; a.S[a.count] = S; a.dS[a.count] = dS;
    a.count++;
  }
}

// Tensor-product support (interface §4.3). Appends (engine node index, S,
// dS[ndm]) entries; Grid::publicId turns an index into the reported id.
inline void particleSupport(const Grid &g, Basis b, const double *x, const double *lp,
                            std::vector<long> &node, std::vector<double> &S,
                            std::vector<double> &dS)
{
  const int nd = g.ndm;
  AxisSupport ax[3];
  for (int d = 0; d < nd; d++) axisSupport(g, b, d, x[d], lp[d], ax[d]);
  const int nz = (nd == 3) ? ax[2].count : 1;
  for (int kz = 0; kz < nz; kz++)
    for (int ky = 0; ky < ax[1].count; ky++)
      for (int kx = 0; kx < ax[0].count; kx++) {
        int ijk[3] = {ax[0].idx[kx], ax[1].idx[ky], (nd == 3) ? ax[2].idx[kz] : 0};
        const double sx = ax[0].S[kx], sy = ax[1].S[ky], sz = (nd == 3) ? ax[2].S[kz] : 1.0;
        const double gx = ax[0].dS[kx], gy = ax[1].dS[ky], gz = (nd == 3) ? ax[2].dS[kz] : 0.0;
        node.push_back(g.extId(ijk));
        S.push_back(sx * sy * sz);
        dS.push_back(gx * sy * sz);
        dS.push_back(sx * gy * sz);
        if (nd == 3) dS.push_back(sx * sy * gz);
      }
}

// --------------------------------------------------------------------------- //
//  The material seam                                                           //
// --------------------------------------------------------------------------- //

// The engine's view of its particle materials. The OpenSees host implements it
// over owned NDMaterial copies; the standalone check over the oracle's
// hypoelastic law. setTrialStrain receives the TOTAL accumulated unrotated
// strain (Voigt, engineering shear) and returns the material's unrotated stress
// WITHOUT the particle's initial stress (the core adds sig0_hat, §9).
class MaterialPoints {
 public:
  virtual ~MaterialPoints() {}
  virtual int    setTrialStrain(int p, const double *epsHat, double *sigMat) = 0;
  virtual int    commitState(int p) = 0;
  virtual int    revertToLastCommit(int p) = 0;
  // Number of commit refusals declared since the last call (then cleared).
  virtual int    pendingRefusals() = 0;
  // Discard a stale refusal count before the particle commit loop, so only a
  // refusal declared by THIS commit can abort the run (review F3).
  virtual void   clearRefusals() {}
  // M_p: largest normal diagonal entry of the current tangent (interface §6.1).
  virtual double normalModulus(int p) = 0;
};

// --------------------------------------------------------------------------- //
//  Parameters, grid conditions, particles                                      //
// --------------------------------------------------------------------------- //

struct Params {
  Basis  basis = BASIS_GIMP;
  double flip = 1.0;           // alpha_F, FLIP fraction
  double damping = 0.0;        // alpha_D, local non-viscous damping
  double masscut = 1.0e-10;    // c_m
  double safety = 0.5;         // s
  double gravity[3] = {0.0, 0.0, 0.0};
  // WP-232j (tests/mpm_antilocking_interface.md J1, K3): 0 off (default, the shipped
  // arithmetic bit for bit), 1 fbar (the cell-centre F-bar on the step increment)
  int    antiLocking = 0;
  // WP-232l (tests/mpm_twophase_interface.md P6): 0 hypo (default, the shipped arithmetic
  // bit for bit), 1 small (dEps = sym(dF - I), R = I, V = V0, no geometric stiffness)
  int    smallKin = 0;
  // WP-232l (P5): the pore-pressure G2P blend alpha_p (1 = incremental, FLIP-like)
  double pBlend = 1.0;
  // WP-232l (§7a Q4, ruling c84f71fb3): route S's p-activity cutoff, `mpm pressureCutoff`;
  // `mpm massCutoff` never moves it
  double pCut = 5.0e-2;
};

// --------------------------------------------------------------------------- //
//  Anti-locking (WP-232j, tests/mpm_antilocking_interface.md J1-J5)            //
// --------------------------------------------------------------------------- //

// The background cell holding the centre x (Grid::cellOf per axis, clamped to
// n-1) and its 2^d corners: engine indices ext[k] (bit d of k set = the +h side
// on axis d) and the Q1 shape-function gradients at the cell centre, gN[k][d] =
// +-1/(2h) in 2D, +-1/(4h) in 3D (J1). Returns 2^d. The corners always lie in
// the particle's own support, for every basis (J3).
inline int cellCorners(const Grid &g, const double *x, long ext[8], double gN[8][3])
{
  const int nd = g.ndm, nc = 1 << nd;
  int c[3] = {0, 0, 0};
  for (int d = 0; d < nd; d++) c[d] = g.cellOf(d, x[d]);
  const double a = (nd == 3) ? 0.25 / g.h : 0.5 / g.h;
  for (int k = 0; k < nc; k++) {
    int ijk[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      const int s = (k >> d) & 1;
      if (d < nd) { ijk[d] = c[d] + s; gN[k][d] = s ? a : -a; }
      else gN[k][d] = 0.0;
    }
    ext[k] = g.extId(ijk);
  }
  return nc;
}

inline double detN(int nd, const double *A)
{
  if (nd == 2) return A[0] * A[3] - A[1] * A[2];
  return A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) +
         A[2] * (A[3] * A[7] - A[4] * A[6]);
}

// J1: dFbar = (det dFc / det dF)^(1/d) dF on the step increment (d = ndm; in plane
// strain the 2x2 in-plane block, F_zz = 1 untouched). Returns 0, or -1 when either
// determinant is <= 0 (dFbar untouched). An equal pair gives the factor 1 exactly.
inline int fbarIncrement(int nd, const double *dF, const double *dFc, double *dFbar)
{
  const double J = detN(nd, dF), Jc = detN(nd, dFc);
  if (!(J > 0.0) || !(Jc > 0.0)) return -1;
  const double r = Jc / J;
  const double s = (r == 1.0) ? 1.0 : ((nd == 3) ? std::cbrt(r) : std::sqrt(r));
  for (int k = 0; k < nd * nd; k++) dFbar[k] = s * dF[k];
  return 0;
}

// Faces: 0 xmin, 1 xmax, 2 ymin, 3 ymax, 4 zmin, 5 zmax (interface §2).
struct FaceConditions {
  unsigned fixMask[6] = {0, 0, 0, 0, 0, 0};          // bit d = component d fixed
  double   wallMu[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};   // < 0: no wall
  // WP-232l (tests/mpm_twophase_interface.md P8, P10): `mpm pfix <face> <p0>` and
  // `mpm plate <face> <F>` (a rigid frictionless plate: one normal unknown U)
  bool     pfixOn[6] = {false, false, false, false, false, false};
  double   pfixVal[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  bool     plateOn[6] = {false, false, false, false, false, false};
  double   plateF[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

// WP-232l (P9): a face folds the u components of its ghosts iff it carries a u condition
// (`fix` or `plate`). Without a plate this is 232i's `fixMask != 0` exactly.
inline bool foldsU(const FaceConditions &fc, int f)
{
  return fc.fixMask[f] != 0u || fc.plateOn[f];
}

// WP-232l (P1): the fluid constants of one two-phase body (LadrunoUP's flags and defaults)
struct TwoPhaseProps {
  double Kf = 0.0, n = 0.0, rhoF = 0.0;
  double kbar[3] = {0.0, 0.0, 0.0};   // k_hydraulic / gamma_w per axis
  double alpha = 1.0, Ks = -1.0;      // Biot alpha; grain bulk modulus (<= 0: incompressible)
  int    stabMode = 1;                // 0 off, 1 auto (stabValue = alpha0), 2 manual (stabValue = alpha)
  double stabValue = 0.25;
  bool   fluidBodyGiven = false;      // else the drive is `mpm gravity`
  double fluidBody[3] = {0.0, 0.0, 0.0};
  double oneOverQbar() const
  {
    double v = n / Kf;
    if (Ks > 0.0) v += (alpha - n) / Ks;
    return v;
  }
};

// The id of plate face f's unknown in support lists (never a lattice index)
inline long plateId(int f) { return -(2L + f); }

inline int faceIndex(const std::string &s)
{
  static const char *names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  for (int f = 0; f < 6; f++)
    if (s == names[f]) return f;
  return -1;
}

// --------------------------------------------------------------------------- //
//  K4 — images on `fix` faces of a B-spline grid (tests/mpm_bspline_interface.md //
//  B5, B6, B8). A ghost node behind a face that carries `fix` is slaved to its  //
//  image: per axis, i < 0 behind a fixed min face maps to -i, i > n behind a    //
//  fixed max face to 2n - i, and the velocity reflects with R = prod diag(s),  //
//  s_a = -1 for a component that face fixes (odd), +1 for a free one (even):   //
//  v_g = R v_image. GIMP and linear grids have no images (B5).                 //
// --------------------------------------------------------------------------- //

// The image of extended-lattice node `ext`; false when the node is its own.
// sg[3] receives R's diagonal.
inline bool imageOfNode(const Grid &g, Basis b, const FaceConditions &fc, long ext, long &img,
                        double sg[3])
{
  sg[0] = sg[1] = sg[2] = 1.0;
  if (!isBSpline(b)) return false;
  int ijk[3];
  g.extIJK(ext, ijk);
  bool moved = false;
  for (int d = 0; d < g.ndm; d++) {
    int f = -1;
    if (ijk[d] < 0 && fc.fixMask[2 * d] != 0u) { ijk[d] = -ijk[d]; f = 2 * d; }
    else if (ijk[d] > g.n[d] && fc.fixMask[2 * d + 1] != 0u) { ijk[d] = 2 * g.n[d] - ijk[d]; f = 2 * d + 1; }
    if (f < 0) continue;
    moved = true;
    for (int a = 0; a < g.ndm; a++)
      if (fc.fixMask[f] & (1u << a)) sg[a] = -sg[a];
  }
  if (!moved) return false;
  img = g.extId(ijk);
  return true;
}

// The mass the active-node cutoff tests (critic A1, WP-232i): on a B-spline grid a
// face-layer node of a `fix` face is its own image, so its one-sided mass is 1 / 2^k of the
// body-plus-reflection mass, k = the number of `fix` faces through the node; the cutoff tests
// m * 2^k, as the full symmetric model would. GIMP and linear: m, untouched.
inline double imageCutoffMass(const Grid &g, Basis b, const FaceConditions &fc, long ext, double m)
{
  if (!isBSpline(b)) return m;
  int ijk[3];
  g.extIJK(ext, ijk);
  double f = 1.0;
  for (int d = 0; d < g.ndm; d++) {
    if (ijk[d] == 0 && foldsU(fc, 2 * d)) f *= 2.0;              // WP-232l: plate faces fold too
    if (ijk[d] == g.n[d] && foldsU(fc, 2 * d + 1)) f *= 2.0;
  }
  return m * f;
}

// The folded support of ONE particle (B5, B8), built per axis and then as a
// tensor product, so the image map (separable per axis) and the signs (a product
// of per-axis signs) are exact, and a folded S' is a product of non-negative 1D
// factors. Appends per own node: node id; the unsigned weight W = S + S_ghost and
// grad W (nd); the per-component S'_{.,a} = S + s_a S_ghost (nd values) and
// grad S'_{.,a} (nd x nd, component-major). For GIMP / linear, or with no `fix`
// face in reach, W = S'_a = S. Returns true when any sign is -1.
inline bool foldedSupport(const Grid &g, Basis b, const FaceConditions &fc, const double *x,
                          const double *lp, std::vector<long> &oNode, std::vector<double> &oW,
                          std::vector<double> &odW, std::vector<double> &oSc,
                          std::vector<double> &odSc)
{
  const int nd = g.ndm;
  // per axis: own targets, unsigned and per-component 1D factors
  int nt[3] = {1, 1, 1}, tgt[3][4];
  double w[3][4], dw[3][4], sc[3][4][3], dsc[3][4][3];
  bool odd = false;
  for (int d = 0; d < nd; d++) {
    AxisSupport ax;
    axisSupport(g, b, d, x[d], lp[d], ax);
    nt[d] = 0;
    for (int k = 0; k < ax.count; k++) {
      int i = ax.idx[k], f = -1;
      if (isBSpline(b)) {
        if (i < 0 && fc.fixMask[2 * d] != 0u) { i = -i; f = 2 * d; }
        else if (i > g.n[d] && fc.fixMask[2 * d + 1] != 0u) { i = 2 * g.n[d] - i; f = 2 * d + 1; }
      }
      int j = 0;
      while (j < nt[d] && tgt[d][j] != i) j++;
      if (j == nt[d]) {
        tgt[d][j] = i; w[d][j] = dw[d][j] = 0.0;
        for (int a = 0; a < 3; a++) sc[d][j][a] = dsc[d][j][a] = 0.0;
        nt[d]++;
      }
      w[d][j] += ax.S[k];
      dw[d][j] += ax.dS[k];
      for (int a = 0; a < nd; a++) {
        const double s = (f >= 0 && (fc.fixMask[f] & (1u << a))) ? -1.0 : 1.0;
        odd = odd || s < 0.0;
        sc[d][j][a] += s * ax.S[k];
        dsc[d][j][a] += s * ax.dS[k];
      }
    }
  }
  const int nz = (nd == 3) ? nt[2] : 1;
  for (int kz = 0; kz < nz; kz++)
    for (int ky = 0; ky < nt[1]; ky++)
      for (int kx = 0; kx < nt[0]; kx++) {
        const int kk[3] = {kx, ky, kz};
        int ijk[3] = {tgt[0][kx], tgt[1][ky], (nd == 3) ? tgt[2][kz] : 0};
        oNode.push_back(g.extId(ijk));
        double W = 1.0;
        for (int d = 0; d < nd; d++) W *= w[d][kk[d]];
        oW.push_back(W);
        for (int e = 0; e < nd; e++) {
          double gW = dw[e][kk[e]];
          for (int d = 0; d < nd; d++) if (d != e) gW *= w[d][kk[d]];
          odW.push_back(gW);
        }
        for (int a = 0; a < nd; a++) {
          double S = 1.0;
          for (int d = 0; d < nd; d++) S *= sc[d][kk[d]][a];
          oSc.push_back(S);
        }
        for (int a = 0; a < nd; a++)
          for (int e = 0; e < nd; e++) {
            double gS = dsc[e][kk[e]][a];
            for (int d = 0; d < nd; d++) if (d != e) gS *= sc[d][kk[d]][a];
            odSc.push_back(gS);
          }
      }
  return odd;
}

// --------------------------------------------------------------------------- //
//  WP-232l (tests/mpm_twophase_interface.md P9, P10, §7a): the general support  //
//  of one particle, per component c = 0..nd-1 (u) and c = nd (p), for a         //
//  two-phase epoch or a model with a plate. Per axis, a 1D entry i maps, per     //
//  component, to (target, sign) pairs:                                           //
//   - beyond a face f (B-spline only): u folds onto the mirror iff foldsU(f),    //
//     sign -1 for a fixed component or the plate's normal one (which also adds  //
//     +2 on the plate unknown: u_g = 2U - u_mirror), +1 otherwise; else the      //
//     ghost stays its own node. p folds iff foldsU(f) or pfix(f), odd under      //
//     pfix (drained), even otherwise (impervious).                              //
//   - the face row of a plate face, normal component: the plate unknown, +1.     //
//  The unsigned weight W (the lumped mass) folds ghosts onto mirrors where u    //
//  folds, never onto the plate; the plate entry's W is its face-row W.          //
//  Appends per entry, in rising engine index (x fastest), plate entries last:   //
//  node id (plateId(f) for a plate), W, grad W (nd), S'_c (nc = nd + 1), grad    //
//  S'_c (nc x nd, component-major). Entries with every weight and gradient 0    //
//  are dropped. GIMP / linear without a plate: S'_c = W = S for every c.        //
// --------------------------------------------------------------------------- //
inline void upSupport(const Grid &g, Basis b, const FaceConditions &fc, const double *x,
                      const double *lp, std::vector<long> &oNode, std::vector<double> &oW,
                      std::vector<double> &odW, std::vector<double> &oSc, std::vector<double> &odSc,
                      double *affOut = 0)
{
  const int nd = g.ndm, nc = nd + 1;
  // WP-232l fix 2 (interface 7b, P9 amended): per axis d, the ghost entries folded EVEN for p
  // (an impervious face) carry the affine image's known term; aff1[d] = sum_g S_g (X_g - X_m),
  // daff1[d] its derivative along d
  double aff1[3] = {0, 0, 0}, daff1[3] = {0, 0, 0};
  const int PLATE = -1000000;            // a 1D target meaning "the plate unknown of this axis"
  struct E1 { int t; double s, ds; };
  std::vector<E1> wE[3], cE[3][4];
  int plateFace[3] = {-1, -1, -1};
  double wFace[3] = {0, 0, 0}, dwFace[3] = {0, 0, 0};   // the W of the plate's face row
  auto add = [](std::vector<E1> &L, int t, double s, double ds) {
    for (E1 &e : L) if (e.t == t) { e.s += s; e.ds += ds; return; }
    L.push_back(E1{t, s, ds});
  };
  const bool bs = isBSpline(b);
  for (int d = 0; d < nd; d++) {
    AxisSupport ax;
    axisSupport(g, b, d, x[d], lp[d], ax);
    for (int k = 0; k < ax.count; k++) {
      const int i = ax.idx[k];
      const double S = ax.S[k], dS = ax.dS[k];
      int f = -1, mir = i;
      if (bs) {
        if (i < 0) { f = 2 * d; mir = -i; }
        else if (i > g.n[d]) { f = 2 * d + 1; mir = 2 * g.n[d] - i; }
      }
      const int rowFace = (i == 0) ? 2 * d : ((i == g.n[d]) ? 2 * d + 1 : -1);
      if (f >= 0 && foldsU(fc, f)) add(wE[d], mir, S, dS);
      else add(wE[d], i, S, dS);
      // (Fable audit finding 1, wp-briefs "232l": this `if` once sat between the W fold and its
      // `else`, so a u-fold + pfix face put the ghost W on the mirror AND the ghost)
      if (f >= 0 && (foldsU(fc, f) || fc.pfixOn[f]) && !fc.pfixOn[f]) {
        aff1[d] += S * (i - mir) * g.h;
        daff1[d] += dS * (i - mir) * g.h;
      }
      if (rowFace >= 0 && fc.plateOn[rowFace]) {
        plateFace[d] = rowFace;
        wFace[d] += S;
        dwFace[d] += dS;
      }
      for (int c = 0; c < nc; c++) {
        std::vector<E1> &L = cE[d][c];
        if (f >= 0) {
          if (c < nd) {
            if (foldsU(fc, f)) {
              const bool plateN = fc.plateOn[f] && c == d;
              const double sg = ((fc.fixMask[f] & (1u << c)) || plateN) ? -1.0 : 1.0;
              add(L, mir, sg * S, sg * dS);
              if (plateN) { add(L, PLATE, 2.0 * S, 2.0 * dS); plateFace[d] = f; }
            } else {
              add(L, i, S, dS);
            }
          } else {
            if (foldsU(fc, f) || fc.pfixOn[f]) {
              const double sg = fc.pfixOn[f] ? -1.0 : 1.0;
              add(L, mir, sg * S, sg * dS);
            } else {
              add(L, i, S, dS);
            }
          }
        } else if (c == d && rowFace >= 0 && fc.plateOn[rowFace]) {
          add(L, PLATE, S, dS);
        } else {
          add(L, i, S, dS);
        }
      }
    }
  }
  struct Acc { double W = 0.0, dW[3] = {0, 0, 0}, Sc[4] = {0, 0, 0, 0}, dSc[4][3] = {}; };
  std::vector<std::pair<long, Acc>> acc;
  auto slot = [&](long id) -> Acc & {
    for (auto &pr : acc) if (pr.first == id) return pr.second;
    acc.push_back(std::make_pair(id, Acc()));
    return acc.back().second;
  };
  const bool three = (nd == 3);
  {
    const int cz = three ? (int)wE[2].size() : 1;
    for (int kz = 0; kz < cz; kz++)
      for (int ky = 0; ky < (int)wE[1].size(); ky++)
        for (int kx = 0; kx < (int)wE[0].size(); kx++) {
          const E1 *e[3] = {&wE[0][kx], &wE[1][ky], three ? &wE[2][kz] : 0};
          int ijk[3] = {e[0]->t, e[1]->t, three ? e[2]->t : 0};
          Acc &A = slot(g.extId(ijk));
          double w = 1.0;
          for (int d = 0; d < nd; d++) w *= e[d]->s;
          A.W += w;
          for (int q = 0; q < nd; q++) {
            double gq = e[q]->ds;
            for (int d = 0; d < nd; d++) if (d != q) gq *= e[d]->s;
            A.dW[q] += gq;
          }
        }
  }
  for (int c = 0; c < nc; c++) {
    const int cz = three ? (int)cE[2][c].size() : 1;
    for (int kz = 0; kz < cz; kz++)
      for (int ky = 0; ky < (int)cE[1][c].size(); ky++)
        for (int kx = 0; kx < (int)cE[0][c].size(); kx++) {
          const E1 *e[3] = {&cE[0][c][kx], &cE[1][c][ky], three ? &cE[2][c][kz] : 0};
          int pAxis = -1;
          for (int d = 0; d < nd; d++) if (e[d]->t == PLATE) pAxis = d;
          long id;
          if (pAxis >= 0) {
            id = plateId(plateFace[pAxis]);
          } else {
            int ijk[3] = {e[0]->t, e[1]->t, three ? e[2]->t : 0};
            id = g.extId(ijk);
          }
          Acc &A = slot(id);
          double w = 1.0;
          for (int d = 0; d < nd; d++) w *= e[d]->s;
          A.Sc[c] += w;
          for (int q = 0; q < nd; q++) {
            double gq = e[q]->ds;
            for (int d = 0; d < nd; d++) if (d != q) gq *= e[d]->s;
            A.dSc[c][q] += gq;
          }
        }
  }
  // WP-232l fix 2: a_d = aff1[d] prod_{e != d} (sum of axis e's p weights), its gradient by the
  // product rule; affOut[d (nd + 1)] = a_d, affOut[d (nd + 1) + 1 + q] = d a_d / d x_q
  if (affOut != 0) {
    double sp[3] = {1, 1, 1}, dsp[3] = {0, 0, 0};
    for (int e = 0; e < nd; e++) {
      sp[e] = 0.0;
      for (const E1 &q : cE[e][nd]) { sp[e] += q.s; dsp[e] += q.ds; }
    }
    for (int d = 0; d < nd; d++) {
      double v = aff1[d];
      for (int e = 0; e < nd; e++) if (e != d) v *= sp[e];
      affOut[d * (nd + 1)] = v;
      for (int q = 0; q < nd; q++) {
        double gq = (q == d) ? daff1[d] : aff1[d] * dsp[q];
        for (int e = 0; e < nd; e++) if (e != d && e != q) gq *= sp[e];
        affOut[d * (nd + 1) + 1 + q] = gq;
      }
    }
  }
  // the plate entries' W: the face-row W on the plate axis times the other axes' sums
  for (int d = 0; d < nd; d++) {
    if (plateFace[d] < 0 || (wFace[d] == 0.0 && dwFace[d] == 0.0)) continue;
    Acc &A = slot(plateId(plateFace[d]));
    double sum[3] = {1, 1, 1}, dsum[3] = {0, 0, 0};
    for (int e = 0; e < nd; e++) {
      if (e == d) continue;
      sum[e] = 0.0;
      for (const E1 &q : wE[e]) { sum[e] += q.s; dsum[e] += q.ds; }
    }
    double w = wFace[d];
    for (int e = 0; e < nd; e++) if (e != d) w *= sum[e];
    A.W += w;
    for (int q = 0; q < nd; q++) {
      double gq = (q == d) ? dwFace[d] : wFace[d] * dsum[q];
      for (int e = 0; e < nd; e++) if (e != d && e != q) gq *= sum[e];
      A.dW[q] += gq;
    }
  }
  std::stable_sort(acc.begin(), acc.end(),
                   [](const std::pair<long, Acc> &a, const std::pair<long, Acc> &b2) {
                     const bool pa = a.first < 0, pb = b2.first < 0;
                     if (pa != pb) return !pa;          // lattice entries first
                     return pa ? a.first > b2.first : a.first < b2.first;
                   });
  for (const auto &pr : acc) {
    const Acc &A = pr.second;
    bool any = A.W != 0.0;
    for (int q = 0; q < nd; q++) any = any || A.dW[q] != 0.0;
    for (int c = 0; c < nc; c++) {
      any = any || A.Sc[c] != 0.0;
      for (int q = 0; q < nd; q++) any = any || A.dSc[c][q] != 0.0;
    }
    if (!any) continue;
    oNode.push_back(pr.first);
    oW.push_back(A.W);
    for (int q = 0; q < nd; q++) odW.push_back(A.dW[q]);
    for (int c = 0; c < nc; c++) oSc.push_back(A.Sc[c]);
    for (int c = 0; c < nc; c++)
      for (int q = 0; q < nd; q++) odSc.push_back(A.dSc[c][q]);
  }
}

// WP-232l (P6, `mpm kinematics small`): the small-strain increment of a relative
// deformation gradient dF (row-major nd x nd): deps = sym(dF - I) in the hypo kernel's
// Voigt order, engineering shear ({xx,yy,zz,xy,yz,zx} in 3D, {xx,yy,xy} in plane strain),
// R1 = I. False when det(dF) <= 0 (an inverted increment is refused as in hypo).
inline bool smallIncrement(int nd, const double *dF, double *deps, double *R1)
{
  if (nd == 3) {
    const double det = dF[0] * (dF[4] * dF[8] - dF[5] * dF[7]) -
                       dF[1] * (dF[3] * dF[8] - dF[5] * dF[6]) +
                       dF[2] * (dF[3] * dF[7] - dF[4] * dF[6]);
    if (!(det > 0.0)) return false;
    deps[0] = dF[0] - 1.0; deps[1] = dF[4] - 1.0; deps[2] = dF[8] - 1.0;
    deps[3] = dF[1] + dF[3]; deps[4] = dF[5] + dF[7]; deps[5] = dF[2] + dF[6];
    for (int k = 0; k < 9; k++) R1[k] = (k % 4 == 0) ? 1.0 : 0.0;
  } else {
    if (!(dF[0] * dF[3] - dF[1] * dF[2] > 0.0)) return false;
    deps[0] = dF[0] - 1.0; deps[1] = dF[3] - 1.0; deps[2] = dF[1] + dF[2];
    R1[0] = 1.0; R1[1] = 0.0; R1[2] = 0.0; R1[3] = 1.0;
  }
  return true;
}

const int64_t MAX_PARTICLE_ID = (int64_t(1) << 53) - 1;   // 232a ruling 2

struct Particles {
  int ndm = 0, nv = 0;
  std::vector<int64_t>  id;
  std::vector<int>      body, mat;
  std::vector<unsigned> flags;
  std::vector<double>   m, V0, lp;                 // constant
  std::vector<double>   x, v, V, F, R, eps, sigHat, sig, sig0;   // committed
  std::vector<double>   a;      // particle acceleration: mode A only (WP-232f); mode F leaves it 0
  // WP-232l (P2): pore pressure p (compression positive), and its rate and acceleration
  // (mode A, Newmark's re-seed); 0 for single-phase particles
  std::vector<double>   pw, pwd, pwdd;

  int size() const { return (int)id.size(); }
  void clear()
  {
    id.clear(); body.clear(); mat.clear(); flags.clear();
    m.clear(); V0.clear(); lp.clear();
    x.clear(); v.clear(); V.clear(); F.clear(); R.clear(); a.clear();
    eps.clear(); sigHat.clear(); sig.clear(); sig0.clear();
    pw.clear(); pwd.clear(); pwdd.clear();
  }
};

// --------------------------------------------------------------------------- //
//  The engine                                                                  //
// --------------------------------------------------------------------------- //

class Engine {
 public:
  Grid           grid;
  Params         par;
  FaceConditions fc;
  Particles      P;
  // WP-232i (B11): per-particle body force per unit MASS (acceleration units),
  // ndm per particle in storage order, NOT scaled by the load factor; empty = none.
  std::vector<double> bodyForce;
  // WP-232j (J5, freeze F1): per particle, the anti-locking path of the last
  // ACCEPTED step -- 0 not averaged (off), 1 averaged, 2 fallback (a corner of its
  // cell not active); empty until a step with anti-locking on. Cumulative
  // particle-step counts over accepted steps (both modes; `mpm reset` zeroes them).
  std::vector<unsigned char> alFlag;
  long alAveraged = 0, alFallback = 0;
  // WP-232k (tests/mpm_soil_interface.md D18, #1109): this step's prescribed grid
  // velocities (the host's `mpm gridvel`), keyed by engine node: bit a of mask = component
  // a set to val[a]. Applied in both projections' grid conditions after the face
  // conditions, so the FLIP increment a = (v1 - vn)/dt sees them. Empty = bit for bit
  // the engine without them.
  struct Presc { unsigned mask; double val[3]; };
  std::unordered_map<long, Presc> presc;
  void countAntiLocking(const std::vector<unsigned char> &f)
  {
    alFlag = f;
    for (unsigned char c : f) {
      if (c == 1) alAveraged++;
      else if (c == 2) alFallback++;
    }
  }

  int ndm() const { return grid.ndm; }
  int nv() const { return grid.ndm == 3 ? 6 : 3; }

  // Set up the grid (clears nothing else). Returns false with a message on bad
  // input.
  bool setGrid(int ndm, const double *x0, double h, const int *n, std::string &err)
  {
    if (ndm != 2 && ndm != 3) { err = "ndm must be 2 (plane strain) or 3"; return false; }
    if (!(h > 0.0)) { err = "cell size h must be > 0"; return false; }
    for (int d = 0; d < ndm; d++)
      if (n[d] < 1) { err = "every cell count must be >= 1"; return false; }
    if (P.size() > 0) { err = "the grid cannot change once particles exist"; return false; }
    grid.ndm = ndm;
    grid.h = h;
    for (int d = 0; d < 3; d++) {
      grid.x0[d] = (d < ndm) ? x0[d] : 0.0;
      grid.n[d] = (d < ndm) ? n[d] : 0;
    }
    P.ndm = ndm;
    P.nv = nv();
    return true;
  }

  // Per-axis cell sizes (a grid read from a file, whose format allows unequal
  // sizes): v1 has one h for every axis (interface §2), so unequal components
  // are refused loudly.
  bool setGrid(int ndm, const double *x0, const double *hAxis, const int *n, std::string &err)
  {
    for (int d = 1; d < ndm; d++)
      if (hAxis[d] != hAxis[0]) {
        err = "the grid's cell size differs between axes; v1 has one cell size h for every"
              " axis (interface §2)";
        return false;
      }
    return setGrid(ndm, x0, hAxis[0], n, err);
  }

  // Add one particle in its initial state (F = R = I, eps = 0, sig = sig0).
  // Returns its index, or -1 with a message.
  int addParticle(int64_t pid, int body, int mat, double rho, double V0,
                  const double *x, const double *lp, const double *v,
                  const double *sig0, std::string &err)
  {
    const int nd = grid.ndm, nvv = nv();
    if (nd == 0) { err = "define the grid first"; return -1; }
    if (pid < 0 || pid > MAX_PARTICLE_ID) { err = "particle id out of [0, 2^53-1]"; return -1; }
    if (ids.count(pid)) { err = "particle id " + std::to_string((long long)pid) + " is already used"; return -1; }
    if (!(V0 > 0.0)) { err = "particle volume must be > 0"; return -1; }
    if (!(rho > 0.0)) { err = "particle density must be > 0"; return -1; }
    ids.insert(pid);
    P.id.push_back(pid);
    P.body.push_back(body);
    P.mat.push_back(mat);
    P.flags.push_back(0u);
    P.m.push_back(rho * V0);
    P.V0.push_back(V0);
    P.V.push_back(V0);
    P.pw.push_back(0.0);                // WP-232l (P2)
    P.pwd.push_back(0.0);
    P.pwdd.push_back(0.0);
    for (int d = 0; d < nd; d++) {
      P.x.push_back(x[d]);
      P.v.push_back(v ? v[d] : 0.0);
      P.a.push_back(0.0);
      P.lp.push_back(lp[d]);
    }
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < nd; j++) {
        P.F.push_back(i == j ? 1.0 : 0.0);
        P.R.push_back(i == j ? 1.0 : 0.0);
      }
    for (int k = 0; k < nvv; k++) {
      const double s = sig0 ? sig0[k] : 0.0;
      P.eps.push_back(0.0);
      P.sig0.push_back(s);
      P.sigHat.push_back(s);
      P.sig.push_back(s);
    }
    return P.size() - 1;
  }

  // Lattice seeding (interface §5). shape: 0 box (geom = lo[ndm], hi[ndm]),
  // 1 disc (2D: cx cy R), 2 ball (3D: cx cy cz R). Ids nextId, nextId+1, ...
  // in a_x-fastest order. Returns the number seeded, or -1 with a message.
  int seedLattice(int shape, const double *geom, int ppc, int body, int mat, double rho,
                  int64_t &nextId, std::string &err)
  {
    const int nd = grid.ndm;
    if (nd == 0) { err = "define the grid first"; return -1; }
    if (ppc < 1) { err = "ppc must be >= 1"; return -1; }
    if (shape == 1 && nd != 2) { err = "disc is a 2D body"; return -1; }
    if (shape == 2 && nd != 3) { err = "ball is a 3D body"; return -1; }
    const double hs = grid.h / ppc;
    const double V0 = std::pow(hs, nd);
    double lp[3], zero[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 3; d++) lp[d] = grid.h / (2.0 * ppc);
    double lo[3], hi[3];
    for (int d = 0; d < nd; d++) {
      if (shape == 0) { lo[d] = geom[d]; hi[d] = geom[nd + d]; }
      else { lo[d] = geom[d] - geom[nd]; hi[d] = geom[d] + geom[nd]; }
    }
    int alo[3] = {0, 0, 0}, ahi[3] = {0, 0, 0};
    for (int d = 0; d < nd; d++) {
      const int total = ppc * grid.n[d];
      alo[d] = std::max(0, (int)std::floor((lo[d] - grid.x0[d]) / hs - 0.5) - 1);
      ahi[d] = std::min(total - 1, (int)std::ceil((hi[d] - grid.x0[d]) / hs - 0.5) + 1);
    }
    int count = 0;
    const int az0 = (nd == 3) ? alo[2] : 0, az1 = (nd == 3) ? ahi[2] : 0;
    for (int az = az0; az <= az1; az++)
      for (int ay = alo[1]; ay <= ahi[1]; ay++)
        for (int axx = alo[0]; axx <= ahi[0]; axx++) {
          const int a[3] = {axx, ay, az};
          double x[3];
          for (int d = 0; d < nd; d++) x[d] = grid.x0[d] + (a[d] + 0.5) * grid.h / ppc;
          bool in = true;
          if (shape == 0) {
            for (int d = 0; d < nd; d++) in = in && (lo[d] <= x[d]) && (x[d] <= hi[d]);
          } else {
            double r2 = 0.0;
            for (int d = 0; d < nd; d++) r2 += (x[d] - geom[d]) * (x[d] - geom[d]);
            in = r2 < geom[nd] * geom[nd];
          }
          if (!in) continue;
          if (addParticle(nextId, body, mat, rho, V0, x, lp, zero, 0, err) < 0) return -1;
          nextId++;
          count++;
        }
    return count;
  }

  // Index of the first particle not admissible at its committed position
  // (interface §4.1 as amended by A2 and A5, §4.2), or -1.
  // When several are inadmissible, the one with the LOWEST id (A5 F6).
  int firstInadmissible() const
  {
    const int nd = grid.ndm;
    int worst = -1;
    for (int p = 0; p < P.size(); p++)
      for (int d = 0; d < nd; d++)
        if (!admissibleAxis(grid, par.basis, d, P.x[nd * p + d], P.lp[nd * p + d])) {
          if (worst < 0 || P.id[p] < P.id[worst]) worst = p;
          break;
        }
    return worst;
  }

  // Translate every particle of a body (the cases' `shift`, A5). Before the
  // first step only; the host enforces that. Returns the count.
  int shiftBody(int body, const double *dx)
  {
    const int nd = grid.ndm;
    int count = 0;
    for (int p = 0; p < P.size(); p++) {
      if (P.body[p] != body) continue;
      for (int d = 0; d < nd; d++) P.x[nd * p + d] += dx[d];
      count++;
    }
    return count;
  }

  // Critical time step (interface §6.1) from the committed state; +inf with no
  // particles.
  double criticalTimeStep(MaterialPoints &mats) const
  {
    const int nd = grid.ndm;
    double wmax = 0.0;
    for (int p = 0; p < P.size(); p++) {
      const double rho = P.m[p] / P.V[p];
      const double M = mats.normalModulus(p);
      const double c = (M > 0.0) ? std::sqrt(M / rho) : 0.0;
      double v2 = 0.0;
      for (int d = 0; d < nd; d++) v2 += P.v[nd * p + d] * P.v[nd * p + d];
      const double w = c + std::sqrt(v2);
      if (w > wmax) wmax = w;
    }
    if (!(wmax > 0.0)) return std::numeric_limits<double>::infinity();
    return grid.h / wmax;
  }

  // Sorted reported ids (Grid::publicId) of the nodes that would be active
  // after P2G at the committed positions (interface §6.4) — the unit-level
  // hook of interface Q4.
  void activeNodes(std::vector<long> &out)
  {
    buildSupports();
    p2gMass();
    collectImages();                    // K4: slaved ghosts are never active nodes
    for (size_t s = 0; s < slave.size(); s++) gm[slaveImg[s]] += gm[slave[s]];
    out.clear();
    const double mcut = par.masscut * massRef();
    for (long i : touched)
      if (!gslave[i] && gm[i] > 0.0 &&
          imageCutoffMass(grid, par.basis, fc, i, gm[i]) > mcut)       // critic A1
        out.push_back(grid.publicId(i));
    std::sort(out.begin(), out.end());
    clearGrid();
  }

  // Steps 3–7 of ADR §5.2 for one Δt at load factor lambda (step 2 is the
  // host's). Returns 0 (committed), -1 (the step failed before any material
  // committed: every particle keeps its committed state, materials reverted),
  // -2 (a material refused at commit: the run must be aborted — W17).
  int advance(double dt, double lambda, MaterialPoints &mats, std::string &err)
  {
    const int N = P.size(), nd = grid.ndm, nvv = nv();
    if (N == 0) return 0;
    if (!bodyForce.empty() && (long)bodyForce.size() != (long)nd * N) {
      err = "the particle body force holds " + std::to_string((long long)bodyForce.size()) +
            " values for " + std::to_string(N) + " particles (ndm per particle): set it again;"
            " step refused, nothing changed";
      return -1;
    }
    buildSupports();

    // ---- step 3: P2G -------------------------------------------------------
    p2gMass();
    // WP-232j (J5, interface §11): each particle's anti-locking path, decided ONCE per step
    // from the P2G mass (a corner below the cutoff -> the fallback) and used by both the
    // conjugate force below and the step-5 kinematics, so the two always agree.
    if (par.antiLocking) {
      const double mcut0 = par.masscut * massRef();
      talFlag.assign(N, 1);
      for (int p = 0; p < N; p++) {
        long ce[8];
        double gN[8][3];
        const int nc = cellCorners(grid, &P.x[nd * p], ce, gN);
        for (int k = 0; k < nc; k++)
          if (!(gm[ce[k]] > 0.0 && imageCutoffMass(grid, par.basis, fc, ce[k], gm[ce[k]]) > mcut0)) {
            talFlag[p] = 2;
            break;
          }
      }
    }
    double b[3];
    for (int d = 0; d < 3; d++) b[d] = lambda * par.gravity[d];
    for (int p = 0; p < N; p++) {
      const double mp = P.m[p], Vp = P.V[p];
      const double *bfp = bodyForce.empty() ? nullptr : &bodyForce[nd * p];
      double sg[3][3];
      voigtToTensor(&P.sig[nvv * p], sg);
      // WP-232j (C1, interface §11): the force work-conjugate to the F-bar kinematics,
      // f_a = V [sigma g_a + pbar (ghat_a - g_a)], pbar = tr_d(sigma) / d, ghat the cell
      // centre's Q1 gradients on its corners (dF_c = I at the step start).
      const bool avg = par.antiLocking && talFlag[p] == 1;
      double pbar = 0.0;
      if (avg) {
        for (int a = 0; a < nd; a++) pbar += sg[a][a];
        pbar /= nd;
      }
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        const double S = supS[k];
        const double *g = &supdS[nd * k];
        for (int a = 0; a < nd; a++) {
          gp[nd * i + a]   += S * mp * P.v[nd * p + a];
          gfe[nd * i + a]  += S * mp * (bfp ? b[a] + bfp[a] : b[a]);
          if (LADRUNO_MUTATE_MPM != MUT_ZERO) {
            double s = 0.0;
            for (int c = 0; c < nd; c++) s += sg[a][c] * g[c];
            if (avg) s -= pbar * g[a];
            gfi[nd * i + a] -= Vp * s;
          }
        }
      }
      if (avg && LADRUNO_MUTATE_MPM != MUT_ZERO) {
        long ce[8];
        double gN[8][3];
        const int nc = cellCorners(grid, &P.x[nd * p], ce, gN);
        for (int k = 0; k < nc; k++)
          for (int a = 0; a < nd; a++) gfi[nd * ce[k] + a] -= Vp * pbar * gN[k][a];
      }
    }

    // K4 (B6): fold the slaved ghosts onto their images, mass +, momentum and
    // forces reflected. GIMP and linear grids have no images.
    collectImages();
    for (size_t s = 0; s < slave.size(); s++) {
      gm[slaveImg[s]] += gm[slave[s]];
      foldVec(gp, s);
      foldVec(gfi, s);
      foldVec(gfe, s);
    }
    // A face-layer node lies ON the symmetry plane: its own image is itself, so
    // the components its face fixes carry no momentum and no force (the full
    // body + reflection gives exactly 0 there; G8). Its FLIP increment is then
    // a = 0, where A5's zeroing gives a = -v^n / dt.
    zeroFaceOdd(gp);
    zeroFaceOdd(gfi);
    zeroFaceOdd(gfe);

    // ---- step 4: damping, grid momentum update, grid conditions ------------
    const double mcut = par.masscut * massRef();
    for (long i : touched) {
      if (gslave[i]) continue;          // K4: set from the image below
      const bool act = gm[i] > 0.0 && imageCutoffMass(grid, par.basis, fc, i, gm[i]) > mcut;
      gact[i] = act ? 1 : 0;
      if (!act) {
        for (int a = 0; a < nd; a++) gv1[nd * i + a] = ga[nd * i + a] = 0.0;
        continue;
      }
      double vn[3], v1[3];
      for (int a = 0; a < nd; a++) {
        vn[a] = gp[nd * i + a] / gm[i];
        double f = gfi[nd * i + a] + gfe[nd * i + a];
        const double sg = (vn[a] > 0.0) ? 1.0 : ((vn[a] < 0.0) ? -1.0 : 0.0);
        f -= par.damping * std::fabs(f) * sg;
        v1[a] = vn[a] + dt * f / gm[i];
      }
      applyConditions(i, v1, true);
      for (int a = 0; a < nd; a++) {
        gv1[nd * i + a] = v1[a];
        ga[nd * i + a] = (v1[a] - vn[a]) / dt;
      }
    }
    for (size_t s = 0; s < slave.size(); s++) {
      gact[slave[s]] = gact[slaveImg[s]];
      reflectVec(gv1, s);
      reflectVec(ga, s);
    }

    // ---- step 5: MUSL particle update (trial) ------------------------------
    tx.assign(P.x.begin(), P.x.end());
    tv.assign(P.v.begin(), P.v.end());
    for (int p = 0; p < N; p++) {
      double sa[3] = {0, 0, 0}, sv[3] = {0, 0, 0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        for (int a = 0; a < nd; a++) {
          sa[a] += supS[k] * ga[nd * i + a];
          sv[a] += supS[k] * gv1[nd * i + a];
        }
      }
      for (int a = 0; a < nd; a++) {
        const double vp = P.v[nd * p + a];
        tv[nd * p + a] = par.flip * (vp + dt * sa[a]) + (1.0 - par.flip) * sv[a];
        double dx = dt * sv[a];
        if (LADRUNO_MUTATE_MPM == MUT_SCALE) dx *= MUT_SCALE_FACTOR;
        tx[nd * p + a] = P.x[nd * p + a] + dx;
      }
    }
    // Interface Amendment 7: an accepted step never leaves a centre outside the
    // real grid. Every TRIAL centre is checked before the second projection and
    // before any material is touched; a step that would carry one out is refused
    // with nothing changed, naming the lowest offending id.
    {
      int worst = -1;
      for (int p = 0; p < N; p++)
        for (int a = 0; a < nd; a++)
          if (!admissibleAxis(grid, par.basis, a, tx[nd * p + a], P.lp[nd * p + a])) {
            if (worst < 0 || P.id[p] < P.id[worst]) worst = p;
            break;
          }
      if (worst >= 0) {
        err = "the step would carry particle " + std::to_string((long long)P.id[worst]) +
              " out of the grid (interface A7); step refused, nothing changed";
        clearGrid();
        return -1;
      }
    }
    // second projection of the updated particle momentum
    for (long i : touched)
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = 0.0;
    for (int p = 0; p < N; p++)
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        for (int a = 0; a < nd; a++)
          gvs[nd * i + a] += supS[k] * P.m[p] * tv[nd * p + a];
      }
    for (size_t s = 0; s < slave.size(); s++) foldVec(gvs, s);   // K4 (B6)
    for (long i : touched) {
      if (gslave[i]) continue;
      double vs[3] = {0, 0, 0};
      if (gact[i]) {
        for (int a = 0; a < nd; a++) vs[a] = gvs[nd * i + a] / gm[i];
        applyConditions(i, vs, false);
      }
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = vs[a];
    }
    for (size_t s = 0; s < slave.size(); s++) reflectVec(gvs, s);
    // velocity gradient, F, the kernel, the constitutive update
    const int nn = nd * nd;
    tF.resize(P.F.size()); tR.resize(P.R.size()); tV.resize(N);
    teps.resize(P.eps.size()); tsh.resize(P.sigHat.size()); tsig.resize(P.sig.size());
    for (int p = 0; p < N; p++) {
      const double *Fn = &P.F[nn * p];
      double dF[9], dFu[9], F1[9];
      // WP-232j (C2): the shared step-5 kinematics -- the cell-centre F-bar on the increment
      // (J1) for an averaged particle, after the grid update and the second projection,
      // before the kernel and the material.
      if (stepIncrement(p, gvs, dt, par.antiLocking ? talFlag[p] : 0, dF, dFu) != 0) {
        err = "the anti-locking increment of particle " + std::to_string((long long)P.id[p]) +
              " refused: det(dF) <= 0 or det(dF_cell) <= 0 (inverted or collapsed)";
        revertAll(mats, p);
        clearGrid();
        return -1;
      }
      for (int a = 0; a < nd; a++)
        for (int c = 0; c < nd; c++) {
          double s = 0.0;
          for (int e = 0; e < nd; e++) s += dF[nd * a + e] * Fn[nd * e + c];
          F1[nd * a + c] = s;
        }
      double deps[6], R1[9], J1 = 0.0;
      int st;
      if (par.smallKin) {              // WP-232l (P6): small strain, R = I, V = V0
        st = smallIncrement(nd, dF, deps, R1) ? 0 : -1;
        J1 = 1.0;
      } else {
        st = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1, deps, R1, &J1)
                       : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1, deps, R1, &J1);
      }
      if (st != 0) {
        err = "the kinematic update of particle " + std::to_string((long long)P.id[p]) +
              " refused: det(F_mid) <= 0 or det(F) <= 0 (inverted or collapsed)";
        revertAll(mats, p);
        clearGrid();
        return -1;
      }
      double e1[6], sm[6];
      for (int k = 0; k < nvv; k++) e1[k] = P.eps[nvv * p + k] + deps[k];
      if (mats.setTrialStrain(p, e1, sm) < 0) {
        err = "the material of particle " + std::to_string((long long)P.id[p]) +
              " failed to integrate the trial strain";
        revertAll(mats, p + 1);
        clearGrid();
        return -1;
      }
      for (int k = 0; k < nvv; k++) {
        teps[nvv * p + k] = e1[k];
        tsh[nvv * p + k] = P.sig0[nvv * p + k] + sm[k];
      }
      if (par.smallKin) for (int k = 0; k < nvv; k++) tsig[nvv * p + k] = tsh[nvv * p + k];
      else if (nd == 3) ladruno_hypo::pushStress6(R1, &tsh[nvv * p], &tsig[nvv * p]);
      else ladruno_hypo::pushStressPlaneStrain(R1, &tsh[nvv * p], &tsig[nvv * p]);
      for (int k = 0; k < nn; k++) { tF[nn * p + k] = F1[k]; tR[nn * p + k] = R1[k]; }
      tV[p] = J1 * P.V0[p];
    }
    clearGrid();

    // ---- step 6: commit the particle materials; a refusal aborts (W17) -----
    mats.clearRefusals();
    for (int p = 0; p < N; p++) mats.commitState(p);
    const int nRef = mats.pendingRefusals();
    if (nRef > 0) {
      err = std::to_string(nRef) + " particle material point(s) REFUSED this commit; sibling "
            "points have already committed, so the run is aborted (ADR 232 W17)";
      return -2;
    }

    // ---- step 7: G2P kinematics take their trial values ---------------------
    P.x.swap(tx); P.v.swap(tv); P.F.swap(tF); P.R.swap(tR); P.V.swap(tV);
    P.eps.swap(teps); P.sigHat.swap(tsh); P.sig.swap(tsig);
    if (par.antiLocking) countAntiLocking(talFlag);
    return 0;
  }

  // `mpm debug affineCheck` in mode F (WP-232j L1, C2): every node of the extended lattice
  // carries the grid velocity A x + b (A row-major ndm x ndm), dt = 1, and every particle's
  // increment comes from the SAME helper as Engine::advance (stepIncrement), averaged. The
  // material is evaluated on the result and reverted: nothing is committed. out =
  // [max_p |dFbar - dF|_max / |dF - I|_max, max_p |deps|_max, #averaged, #refused,
  //  max_p |sigHat_trial - sigHat_n|_max (the spurious stress of a rigid rotation)].
  void affineCheck(const double *A, const double *b, MaterialPoints &mats, std::vector<double> &out)
  {
    const int N = P.size(), nd = grid.ndm, nn = nd * nd, nvv = nv();
    buildSupports();
    const long nN = grid.numNodes();
    std::vector<double> gvel((size_t)nd * nN, 0.0);
    for (long i = 0; i < nN; i++) {
      int ijk[3];
      grid.extIJK(i, ijk);
      for (int a = 0; a < nd; a++) {
        double s = b[a];
        for (int c = 0; c < nd; c++) s += A[nd * a + c] * (grid.x0[c] + ijk[c] * grid.h);
        gvel[nd * i + a] = s;
      }
    }
    double worst = 0.0, depsMax = 0.0, spur = 0.0;
    long nAvg = 0, nBad = 0;
    for (int p = 0; p < N; p++) {
      double dF[9], dFu[9];
      if (stepIncrement(p, gvel, 1.0, 1, dF, dFu) != 0) { nBad++; continue; }
      double num = 0.0, den = 0.0;
      for (int k = 0; k < nn; k++) {
        num = std::max(num, std::fabs(dF[k] - dFu[k]));
        den = std::max(den, std::fabs(dFu[k] - ((k % (nd + 1) == 0) ? 1.0 : 0.0)));
      }
      if (den > 0.0) worst = std::max(worst, num / den);
      const double *Fn = &P.F[nn * p];
      double F1[9], deps[6], R1[9], J1;
      for (int a = 0; a < nd; a++)
        for (int c = 0; c < nd; c++) {
          double s = 0.0;
          for (int e = 0; e < nd; e++) s += dF[nd * a + e] * Fn[nd * e + c];
          F1[nd * a + c] = s;
        }
      const int st = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1, deps, R1, &J1)
                               : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1, deps, R1, &J1);
      if (st != 0) { nBad++; continue; }
      nAvg++;
      double e1[6], sm[6];
      for (int k = 0; k < nvv; k++) {
        depsMax = std::max(depsMax, std::fabs(deps[k]));
        e1[k] = P.eps[nvv * p + k] + deps[k];
      }
      if (mats.setTrialStrain(p, e1, sm) == 0)
        for (int k = 0; k < nvv; k++)
          spur = std::max(spur, std::fabs(P.sig0[nvv * p + k] + sm[k] - P.sigHat[nvv * p + k]));
      mats.revertToLastCommit(p);
    }
    out = {worst, depsMax, (double)nAvg, (double)nBad, spur};
  }

  // WP-232j (C2): the step-5 kinematics of particle p, shared by advance() and
  // affineCheck(): dFu = I + dt sum_a v_a (x) gradS_a over its support (the grid velocity
  // array gvel, ndm per engine node), and dF = dFu, or for flag 1 (averaged) the J1
  // dFbar = (det dF_c / det dFu)^(1/d) dFu with dF_c = I + dt sum_corners v_k (x) gradN_k(x_c).
  // Returns 0, or -1 on a non-positive determinant.
  int stepIncrement(int p, const std::vector<double> &gvel, double dt, int flag, double *dF,
                    double *dFu)
  {
    const int nd = grid.ndm;
    double L[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (long k = supStart[p]; k < supStart[p + 1]; k++) {
      const long i = supNode[k];
      const double *g = &supdS[nd * k];
      for (int a = 0; a < nd; a++)
        for (int c = 0; c < nd; c++)
          L[nd * a + c] += gvel[nd * i + a] * g[c];
    }
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++)
        dF[nd * a + c] = dFu[nd * a + c] = (a == c ? 1.0 : 0.0) + dt * L[nd * a + c];
    if (flag != 1) return 0;
    long ce[8];
    double gN[8][3], dFc[9];
    const int nc = cellCorners(grid, &P.x[nd * p], ce, gN);
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++) {
        double s = (a == c) ? 1.0 : 0.0;
        for (int k = 0; k < nc; k++) s += dt * gvel[nd * ce[k] + a] * gN[k][c];
        dFc[nd * a + c] = s;
      }
    return fbarIncrement(nd, dFu, dFc, dF);
  }

  // Forget the id registry (used when the host clears its particles).
  void clearParticles() { P.clear(); ids.clear(); P.ndm = grid.ndm; P.nv = nv(); }
  // Replace the particle set (e.g. `mpm reset` restoring a snapshot) and rebuild
  // the id registry from it, so ids added after the snapshot are free again.
  void setParticles(const Particles &s)
  {
    P = s;
    ids.clear();
    for (int64_t pid : P.id) ids.insert(pid);
  }
  bool hasId(int64_t pid) const { return ids.count(pid) != 0; }

 private:
  std::unordered_set<int64_t> ids;
  // supports at the step-start positions, flattened
  std::vector<long>   supStart, supNode;
  std::vector<double> supS, supdS;
  // grid arrays, dense, reset through the touched list
  std::vector<double> gm, gp, gfi, gfe, gv1, ga, gvs;
  std::vector<char>   gmark, gact, gslave;
  std::vector<long>   touched;
  // K4 (B6): this step's slaved ghost nodes, their images and R's diagonal
  std::vector<long>   slave, slaveImg;
  std::vector<double> slaveSg;
  // trial particle arrays
  std::vector<double> tx, tv, tF, tR, tV, teps, tsh, tsig;
  std::vector<unsigned char> talFlag;     // WP-232j: this step's anti-locking path per particle

  double massRef() const
  {
    double r = 0.0;
    for (double m : P.m) r = std::max(r, m);
    return r;
  }

  void voigtToTensor(const double *s, double t[3][3]) const
  {
    if (grid.ndm == 3) {
      t[0][0] = s[0]; t[1][1] = s[1]; t[2][2] = s[2];
      t[0][1] = t[1][0] = s[3]; t[1][2] = t[2][1] = s[4]; t[2][0] = t[0][2] = s[5];
    } else {
      t[0][0] = s[0]; t[1][1] = s[1]; t[0][1] = t[1][0] = s[2];
      t[0][2] = t[2][0] = t[1][2] = t[2][1] = t[2][2] = 0.0;
    }
  }

  void buildSupports()
  {
    const int N = P.size(), nd = grid.ndm;
    supStart.assign(1, 0);
    supNode.clear(); supS.clear(); supdS.clear();
    for (int p = 0; p < N; p++) {
      particleSupport(grid, par.basis, &P.x[nd * p], &P.lp[nd * p], supNode, supS, supdS);
      supStart.push_back((long)supNode.size());
    }
    const long nN = grid.numNodes();
    if ((long)gm.size() != nN) {
      gm.assign(nN, 0.0);
      gp.assign(nd * nN, 0.0); gfi.assign(nd * nN, 0.0); gfe.assign(nd * nN, 0.0);
      gv1.assign(nd * nN, 0.0); ga.assign(nd * nN, 0.0); gvs.assign(nd * nN, 0.0);
      gmark.assign(nN, 0); gact.assign(nN, 0); gslave.assign(nN, 0);
      touched.clear();
    }
  }

  void p2gMass()
  {
    const int N = P.size();
    for (int p = 0; p < N; p++)
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        if (!gmark[i]) { gmark[i] = 1; touched.push_back(i); }
        gm[i] += supS[k] * P.m[p];
      }
  }

  void clearGrid()
  {
    const int nd = grid.ndm;
    for (long i : touched) {
      gm[i] = 0.0;
      gmark[i] = 0;
      gact[i] = 0;
      gslave[i] = 0;
      for (int a = 0; a < nd; a++)
        gp[nd * i + a] = gfi[nd * i + a] = gfe[nd * i + a] =
        gv1[nd * i + a] = ga[nd * i + a] = gvs[nd * i + a] = 0.0;
    }
    touched.clear();
    slave.clear(); slaveImg.clear(); slaveSg.clear();
  }

  // K4 (B6): find the touched ghosts slaved to an image (B-spline grids with a
  // `fix` face only); an untouched image joins `touched`.
  void collectImages()
  {
    slave.clear(); slaveImg.clear(); slaveSg.clear();
    if (!isBSpline(par.basis)) return;
    const size_t nt = touched.size();
    for (size_t t = 0; t < nt; t++) {
      long img = -1;
      double sg[3];
      if (!imageOfNode(grid, par.basis, fc, touched[t], img, sg)) continue;
      gslave[touched[t]] = 1;
      slave.push_back(touched[t]);
      slaveImg.push_back(img);
      for (int d = 0; d < 3; d++) slaveSg.push_back(sg[d]);
      if (!gmark[img]) { gmark[img] = 1; touched.push_back(img); }
    }
  }
  // arr[image] += R arr[slave]
  void foldVec(std::vector<double> &arr, size_t s)
  {
    const int nd = grid.ndm;
    for (int a = 0; a < nd; a++)
      arr[nd * slaveImg[s] + a] += slaveSg[3 * s + a] * arr[nd * slave[s] + a];
  }
  // K4 (B6): on a B-spline grid, zero the fixed components of arr on the face
  // layer of every `fix` face (touched nodes)
  void zeroFaceOdd(std::vector<double> &arr)
  {
    if (!isBSpline(par.basis)) return;
    const int nd = grid.ndm;
    for (long i : touched) {
      int ijk[3];
      grid.extIJK(i, ijk);
      for (int d = 0; d < nd; d++) {
        unsigned m = 0u;
        if (ijk[d] == 0) m |= fc.fixMask[2 * d];
        if (ijk[d] == grid.n[d]) m |= fc.fixMask[2 * d + 1];
        for (int a = 0; a < nd; a++)
          if (m & (1u << a)) arr[nd * i + a] = 0.0;
      }
    }
  }
  // arr[slave] = R arr[image]
  void reflectVec(std::vector<double> &arr, size_t s)
  {
    const int nd = grid.ndm;
    for (int a = 0; a < nd; a++)
      arr[nd * slave[s] + a] = slaveSg[3 * s + a] * arr[nd * slaveImg[s] + a];
  }

  void revertAll(MaterialPoints &mats, int upTo)
  {
    for (int q = 0; q < upTo; q++) mats.revertToLastCommit(q);
  }

  // Grid-face conditions on one active node (interface §7). full = the
  // operator on v^{n+1} (walls with friction, then fixes); otherwise the
  // second-projection operator (walls' no-penetration part only, then fixes).
  // A face's condition acts on its face nodes AND its ghost layers (A5).
  void applyConditions(long i, double *v, bool full) const
  {
    const int nd = grid.ndm;
    int ijk[3];
    grid.extIJK(i, ijk);
    bool on[6];
    for (int d = 0; d < 3; d++) {
      on[2 * d]     = (d < nd) && ijk[d] <= 0;
      on[2 * d + 1] = (d < nd) && ijk[d] >= grid.n[d];
    }
    for (int f = 0; f < 2 * nd; f++) {
      if (!on[f] || fc.wallMu[f] < 0.0) continue;
      const int d = f / 2;
      const double sgn = (f % 2 == 0) ? 1.0 : -1.0;   // inward normal n = sgn * e_d
      const double vn = sgn * v[d];
      if (!(vn < 0.0)) continue;
      double vt[3];
      for (int a = 0; a < nd; a++) vt[a] = v[a];
      vt[d] -= vn * sgn;                                // v - vn n
      if (!full) {
        for (int a = 0; a < nd; a++) v[a] = vt[a];
        continue;
      }
      double tn = 0.0;
      for (int a = 0; a < nd; a++) tn += vt[a] * vt[a];
      tn = std::sqrt(tn);
      if (tn > 0.0) {
        const double red = std::min(fc.wallMu[f] * std::fabs(vn), tn);
        for (int a = 0; a < nd; a++) v[a] = vt[a] - red * vt[a] / tn;
      } else {
        for (int a = 0; a < nd; a++) v[a] = vt[a];
      }
    }
    for (int f = 0; f < 2 * nd; f++) {
      if (!on[f] || fc.fixMask[f] == 0u) continue;
      for (int a = 0; a < nd; a++)
        if (fc.fixMask[f] & (1u << a)) v[a] = 0.0;
    }
    if (!presc.empty()) {                    // WP-232k (D18): `mpm gridvel` in mode F
      const auto it = presc.find(i);
      if (it != presc.end())
        for (int a = 0; a < nd; a++)
          if (it->second.mask & (1u << a)) v[a] = it->second.val[a];
    }
  }
};

} // namespace ladruno_mpm

#endif
