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

enum Basis { BASIS_GIMP = 0, BASIS_LINEAR = 1 };

// --------------------------------------------------------------------------- //
//  Basis (interface §4)                                                        //
// --------------------------------------------------------------------------- //

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

struct Grid {
  int    ndm = 0;
  double x0[3] = {0.0, 0.0, 0.0};
  double h = 0.0;
  int    n[3] = {1, 1, 1};

  long numNodes() const
  {
    long c = 1;
    for (int d = 0; d < ndm; d++) c *= (long)(n[d] + 1);
    return c;
  }
  long nodeId(const int ijk[3]) const
  {
    if (ndm == 2) return ijk[0] + (long)(n[0] + 1) * ijk[1];
    return ijk[0] + (long)(n[0] + 1) * (ijk[1] + (long)(n[1] + 1) * ijk[2]);
  }
  void nodeIJK(long id, int ijk[3]) const
  {
    ijk[2] = 0;
    ijk[0] = (int)(id % (n[0] + 1));
    long rest = id / (n[0] + 1);
    if (ndm == 2) { ijk[1] = (int)rest; return; }
    ijk[1] = (int)(rest % (n[1] + 1));
    ijk[2] = (int)(rest / (n[1] + 1));
  }
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

// Admissibility on one axis (interface §4.1 as amended by A2, §4.2): the
// particle CENTRE inside the closed grid box, and for GIMP 0 < lp <= h/2.
// Support nodes outside the grid are dropped by axisSupport (truncation).
inline bool admissibleAxis(const Grid &g, Basis b, int d, double x, double lp)
{
  const double lo = g.x0[d], hi = g.x0[d] + g.n[d] * g.h;
  if (b == BASIS_LINEAR)
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
  for (int i = c - 1; i <= c + 2; i++) {
    if (i < 0 || i > g.n[d]) continue;
    double S, dS;
    gimpFactor(x - (g.x0[d] + i * h), h, lp, S, dS);
    if (S == 0.0 && dS == 0.0) continue;
    a.idx[a.count] = i; a.S[a.count] = S; a.dS[a.count] = dS;
    a.count++;
  }
}

// Tensor-product support (interface §4.3). Appends (node, S, dS[ndm]) entries.
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
        node.push_back(g.nodeId(ijk));
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
};

// Faces: 0 xmin, 1 xmax, 2 ymin, 3 ymax, 4 zmin, 5 zmax (interface §2).
struct FaceConditions {
  unsigned fixMask[6] = {0, 0, 0, 0, 0, 0};          // bit d = component d fixed
  double   wallMu[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};   // < 0: no wall
};

inline int faceIndex(const std::string &s)
{
  static const char *names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  for (int f = 0; f < 6; f++)
    if (s == names[f]) return f;
  return -1;
}

const int64_t MAX_PARTICLE_ID = (int64_t(1) << 53) - 1;   // 232a ruling 2

struct Particles {
  int ndm = 0, nv = 0;
  std::vector<int64_t>  id;
  std::vector<int>      body, mat;
  std::vector<unsigned> flags;
  std::vector<double>   m, V0, lp;                 // constant
  std::vector<double>   x, v, V, F, R, eps, sigHat, sig, sig0;   // committed

  int size() const { return (int)id.size(); }
  void clear()
  {
    id.clear(); body.clear(); mat.clear(); flags.clear();
    m.clear(); V0.clear(); lp.clear();
    x.clear(); v.clear(); V.clear(); F.clear(); R.clear();
    eps.clear(); sigHat.clear(); sig.clear(); sig0.clear();
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
    for (int d = 0; d < nd; d++) {
      P.x.push_back(x[d]);
      P.v.push_back(v ? v[d] : 0.0);
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
  // (interface §4.1 as amended by A2, §4.2), or -1.
  int firstInadmissible() const
  {
    const int nd = grid.ndm;
    for (int p = 0; p < P.size(); p++)
      for (int d = 0; d < nd; d++)
        if (!admissibleAxis(grid, par.basis, d, P.x[nd * p + d], P.lp[nd * p + d]))
          return p;
    return -1;
  }

  // Index of the first GIMP particle whose domain reaches within lp of a grid
  // face that carries no grid condition (no fix, no wall), or -1. Such a
  // particle loses partition of unity there (A2): the host warns once per run.
  int firstNearFreeFace() const
  {
    if (par.basis != BASIS_GIMP) return -1;
    const int nd = grid.ndm;
    for (int p = 0; p < P.size(); p++)
      for (int d = 0; d < nd; d++) {
        const double x = P.x[nd * p + d], lp = P.lp[nd * p + d];
        const double lo = grid.x0[d], hi = grid.x0[d] + grid.n[d] * grid.h;
        const bool freeLo = fc.fixMask[2 * d] == 0u && fc.wallMu[2 * d] < 0.0;
        const bool freeHi = fc.fixMask[2 * d + 1] == 0u && fc.wallMu[2 * d + 1] < 0.0;
        if ((freeLo && x - lp < lo) || (freeHi && x + lp > hi)) return p;
      }
    return -1;
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

  // Sorted ids of the nodes that would be active after P2G at the committed
  // positions (interface §6.4) — the unit-level hook of interface Q4.
  void activeNodes(std::vector<long> &out)
  {
    buildSupports();
    p2gMass();
    out.clear();
    const double mcut = par.masscut * massRef();
    for (long i : touched)
      if (gm[i] > 0.0 && gm[i] > mcut) out.push_back(i);
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
    buildSupports();

    // ---- step 3: P2G -------------------------------------------------------
    p2gMass();
    double b[3];
    for (int d = 0; d < 3; d++) b[d] = lambda * par.gravity[d];
    for (int p = 0; p < N; p++) {
      const double mp = P.m[p], Vp = P.V[p];
      double sg[3][3];
      voigtToTensor(&P.sig[nvv * p], sg);
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        const double S = supS[k];
        const double *g = &supdS[nd * k];
        for (int a = 0; a < nd; a++) {
          gp[nd * i + a]   += S * mp * P.v[nd * p + a];
          gfe[nd * i + a]  += S * mp * b[a];
          if (LADRUNO_MUTATE_MPM != MUT_ZERO) {
            double s = 0.0;
            for (int c = 0; c < nd; c++) s += sg[a][c] * g[c];
            gfi[nd * i + a] -= Vp * s;
          }
        }
      }
    }

    // ---- step 4: damping, grid momentum update, grid conditions ------------
    const double mcut = par.masscut * massRef();
    for (long i : touched) {
      const bool act = gm[i] > 0.0 && gm[i] > mcut;
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
    // second projection of the updated particle momentum
    for (long i : touched)
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = 0.0;
    for (int p = 0; p < N; p++)
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        for (int a = 0; a < nd; a++)
          gvs[nd * i + a] += supS[k] * P.m[p] * tv[nd * p + a];
      }
    for (long i : touched) {
      double vs[3] = {0, 0, 0};
      if (gact[i]) {
        for (int a = 0; a < nd; a++) vs[a] = gvs[nd * i + a] / gm[i];
        applyConditions(i, vs, false);
      }
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = vs[a];
    }
    // velocity gradient, F, the kernel, the constitutive update
    const int nn = nd * nd;
    tF.resize(P.F.size()); tR.resize(P.R.size()); tV.resize(N);
    teps.resize(P.eps.size()); tsh.resize(P.sigHat.size()); tsig.resize(P.sig.size());
    for (int p = 0; p < N; p++) {
      double L[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        const double *g = &supdS[nd * k];
        for (int a = 0; a < nd; a++)
          for (int c = 0; c < nd; c++)
            L[nd * a + c] += gvs[nd * i + a] * g[c];
      }
      const double *Fn = &P.F[nn * p];
      double dF[9], F1[9];
      for (int a = 0; a < nd; a++)
        for (int c = 0; c < nd; c++)
          dF[nd * a + c] = (a == c ? 1.0 : 0.0) + dt * L[nd * a + c];
      for (int a = 0; a < nd; a++)
        for (int c = 0; c < nd; c++) {
          double s = 0.0;
          for (int e = 0; e < nd; e++) s += dF[nd * a + e] * Fn[nd * e + c];
          F1[nd * a + c] = s;
        }
      double deps[6], R1[9], J1 = 0.0;
      const int st = (nd == 3)
        ? ladruno_hypo::hypoIncrementFromF(Fn, F1, deps, R1, &J1)
        : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1, deps, R1, &J1);
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
      if (nd == 3) ladruno_hypo::pushStress6(R1, &tsh[nvv * p], &tsig[nvv * p]);
      else ladruno_hypo::pushStressPlaneStrain(R1, &tsh[nvv * p], &tsig[nvv * p]);
      for (int k = 0; k < nn; k++) { tF[nn * p + k] = F1[k]; tR[nn * p + k] = R1[k]; }
      tV[p] = J1 * P.V0[p];
    }
    clearGrid();

    // ---- step 6: commit the particle materials; a refusal aborts (W17) -----
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
    return 0;
  }

  // Forget the id registry (used when the host clears its particles).
  void clearParticles() { P.clear(); ids.clear(); P.ndm = grid.ndm; P.nv = nv(); }
  bool hasId(int64_t pid) const { return ids.count(pid) != 0; }

 private:
  std::unordered_set<int64_t> ids;
  // supports at the step-start positions, flattened
  std::vector<long>   supStart, supNode;
  std::vector<double> supS, supdS;
  // grid arrays, dense, reset through the touched list
  std::vector<double> gm, gp, gfi, gfe, gv1, ga, gvs;
  std::vector<char>   gmark, gact;
  std::vector<long>   touched;
  // trial particle arrays
  std::vector<double> tx, tv, tF, tR, tV, teps, tsh, tsig;

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
      gmark.assign(nN, 0); gact.assign(nN, 0);
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
      for (int a = 0; a < nd; a++)
        gp[nd * i + a] = gfi[nd * i + a] = gfe[nd * i + a] =
        gv1[nd * i + a] = ga[nd * i + a] = gvs[nd * i + a] = 0.0;
    }
    touched.clear();
  }

  void revertAll(MaterialPoints &mats, int upTo)
  {
    for (int q = 0; q < upTo; q++) mats.revertToLastCommit(q);
  }

  // Grid-face conditions on one active node (interface §7). full = the
  // operator on v^{n+1} (walls with friction, then fixes); otherwise the
  // second-projection operator (walls' no-penetration part only, then fixes).
  void applyConditions(long i, double *v, bool full) const
  {
    const int nd = grid.ndm;
    int ijk[3];
    grid.nodeIJK(i, ijk);
    bool on[6];
    for (int d = 0; d < 3; d++) {
      on[2 * d]     = (d < nd) && ijk[d] == 0;
      on[2 * d + 1] = (d < nd) && ijk[d] == grid.n[d];
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
  }
};

} // namespace ladruno_mpm

#endif
