/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// LadrunoCPMLKernel.h — Standalone header-only C++17 C-PML (Complex
// Frequency-Shifted Perfectly Matched Layer) time-domain finite element kernel
// via exact piecewise-linear recursive convolution.
// Zero OpenSees Element/Domain/Matrix/Vector dependencies.

#ifndef LadrunoCPMLKernel_h
#define LadrunoCPMLKernel_h

#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>

#ifndef LadrunoCPMLShapes_h
#  define LadrunoCPML ladruno_cpml_impl
#  include "LadrunoCPMLShapes.h"
#  undef LadrunoCPML
#  define LADRUNO_CPML_KERNEL_NS ladruno_cpml_impl
#else
#  define LADRUNO_CPML_KERNEL_NS LadrunoCPML
#endif

// Mutation hooks for P2 mutation testing (standalone compile-time flags, zero
// OpenSees header dependencies under all preprocessor configurations).
#ifndef LADRUNO_CPML_MUT_ACTIVE
#  define LADRUNO_CPML_MUT_ACTIVE(NAME) (LADRUNO_MUTATION_##NAME != 0)
#endif
#ifndef LADRUNO_MUTATION_CPML_NO_DAMPING
#  define LADRUNO_MUTATION_CPML_NO_DAMPING 0
#endif
#ifndef LADRUNO_MUTATION_CPML_NO_KAPPA
#  define LADRUNO_MUTATION_CPML_NO_KAPPA 0
#endif
#ifndef LADRUNO_MUTATION_CPML_NO_ALPHA
#  define LADRUNO_MUTATION_CPML_NO_ALPHA 0
#endif
#ifndef LADRUNO_MUTATION_CPML_WRONG_SIGN_MEMORY
#  define LADRUNO_MUTATION_CPML_WRONG_SIGN_MEMORY 0
#endif
#ifndef LADRUNO_MUTATION_CPML_NO_CORNER_CROSS
#  define LADRUNO_MUTATION_CPML_NO_CORNER_CROSS 0
#endif

namespace LADRUNO_CPML_KERNEL_NS {

// =========================================================================== //
// 1. Exact Piecewise-Linear Recursive Convolution Filter (RCFilter1D)
//    First-order rational operator H(w) = c0 + c1 / (beta + i*w)
// =========================================================================== //
struct RCFilter1D {
  double c0 = 1.0;
  double c1 = 0.0;
  double beta = 0.0;

  // Committed state at t_n
  double zeta_c = 0.0;
  double x_c = 0.0;
  double y_c = 0.0;
  double D_c = 1.0;

  // Trial state at t_{n+1}
  double zeta_t = 0.0;
  double x_t = 0.0;
  double y_t = 0.0;
  double D_t = 1.0;

  // Cached exact piecewise-linear coefficients for constant-dt stepping
  double cached_dt = -1.0;
  double cached_b = 1.0;
  double cached_phi0 = 0.0;
  double cached_phi1 = 0.0;
  double cached_D = 1.0;

  RCFilter1D() = default;
  RCFilter1D(double c0_, double c1_, double beta_)
      : c0(c0_), c1(c1_), beta(beta_),
        zeta_c(0.0), x_c(0.0), y_c(0.0), D_c(c0_),
        zeta_t(0.0), x_t(0.0), y_t(0.0), D_t(c0_),
        cached_dt(-1.0), cached_b(1.0), cached_phi0(0.0), cached_phi1(0.0), cached_D(c0_) {}

  static void computeCoeffs(double beta_val, double dt,
                            double& b, double& phi0, double& phi1) {
    const double theta = beta_val * dt;
    if (std::fabs(theta) < 1.0e-8) {
      const double t2 = theta * theta;
      const double t3 = t2 * theta;
      const double t4 = t2 * t2;
      b    = 1.0 - theta + 0.5 * t2 - (1.0 / 6.0) * t3 + (1.0 / 24.0) * t4;
      phi0 = dt * (0.5 - (1.0 / 3.0) * theta + 0.125 * t2 - (1.0 / 30.0) * t3);
      phi1 = dt * (0.5 - (1.0 / 6.0) * theta + (1.0 / 24.0) * t2 - (1.0 / 120.0) * t3);
    } else {
      b = std::exp(-theta);
      const double denom = beta_val * beta_val * dt;
      phi0 = (1.0 - b * (1.0 + theta)) / denom;
      phi1 = (theta - 1.0 + b) / denom;
    }
  }

  double instantaneousMultiplier(double dt) const {
    if (c1 == 0.0 || dt <= 0.0) {
      return c0;
    }
    if (dt == cached_dt) {
      return cached_D;
    }
    double b, phi0, phi1;
    computeCoeffs(beta, dt, b, phi0, phi1);
    return c0 + c1 * phi1;
  }

  double evalTrial(double x_new, double dt) {
    x_t = x_new;
    if (c1 == 0.0) {
      zeta_t = 0.0;
      y_t = c0 * x_new;
      D_t = c0;
      return y_t;
    }
    const double sgn = LADRUNO_CPML_MUT_ACTIVE(CPML_WRONG_SIGN_MEMORY) ? -1.0 : 1.0;
    if (dt <= 0.0) {
      zeta_t = zeta_c;
      y_t = c0 * x_new + sgn * c1 * zeta_t;
      D_t = c0;
      return y_t;
    }
    if (dt != cached_dt) {
      computeCoeffs(beta, dt, cached_b, cached_phi0, cached_phi1);
      cached_D = c0 + c1 * cached_phi1;
      cached_dt = dt;
    }
    zeta_t = cached_b * zeta_c + cached_phi0 * x_c + cached_phi1 * x_new;
    y_t = c0 * x_new + sgn * c1 * zeta_t;
    D_t = cached_D;
    return y_t;
  }

  void commit() {
    zeta_c = zeta_t;
    x_c = x_t;
    y_c = y_t;
    D_c = D_t;
  }

  void revertToLastCommit() {
    zeta_t = zeta_c;
    x_t = x_c;
    y_t = y_c;
    D_t = D_c;
  }

  void revertToStart() {
    zeta_c = 0.0;
    x_c = 0.0;
    y_c = 0.0;
    D_c = c0;
    zeta_t = 0.0;
    x_t = 0.0;
    y_t = 0.0;
    D_t = c0;
    cached_dt = -1.0;
  }

  static constexpr int numStateDoubles() { return 8; }

  void packState(double*& ptr) const {
    *ptr++ = zeta_c;
    *ptr++ = x_c;
    *ptr++ = y_c;
    *ptr++ = D_c;
    *ptr++ = zeta_t;
    *ptr++ = x_t;
    *ptr++ = y_t;
    *ptr++ = D_t;
  }

  void unpackState(const double*& ptr) {
    zeta_c = *ptr++;
    x_c    = *ptr++;
    y_c    = *ptr++;
    D_c    = *ptr++;
    zeta_t = *ptr++;
    x_t    = *ptr++;
    y_t    = *ptr++;
    D_t    = *ptr++;
    cached_dt = -1.0;
  }
};

// =========================================================================== //
// 2. Cascade of K Scalar RCFilter1D Stages (RCFilterCascade<MaxStages>)
// =========================================================================== //
template <int MaxStages>
struct RCFilterCascade {
  int nStages = 0;
  std::array<RCFilter1D, MaxStages> stages{};
  double y_c = 0.0;
  double D_c = 1.0;
  double y_t = 0.0;
  double D_t = 1.0;

  void setIdentity() {
    nStages = 0;
    y_c = 0.0;
    D_c = 1.0;
    y_t = 0.0;
    D_t = 1.0;
  }

  void configure(int kStages, const RCFilter1D* stg) {
    assert(kStages >= 0 && kStages <= MaxStages);
    nStages = kStages;
    for (int k = 0; k < nStages; ++k) {
      stages[k] = stg[k];
    }
    const double d0 = instantaneousMultiplier(0.0);
    y_c = 0.0;
    D_c = d0;
    y_t = 0.0;
    D_t = d0;
  }

  void reconfigureCoeffs(int kStages, const RCFilter1D* stg) {
    assert(kStages >= 0 && kStages <= MaxStages);
    nStages = kStages;
    for (int k = 0; k < nStages; ++k) {
      stages[k].c0 = stg[k].c0;
      stages[k].c1 = stg[k].c1;
      stages[k].beta = stg[k].beta;
      stages[k].cached_dt = -1.0;
    }
  }

  double instantaneousMultiplier(double dt) const {
    if (nStages == 0) return 1.0;
    double prod = 1.0;
    for (int k = 0; k < nStages; ++k) {
      prod *= stages[k].instantaneousMultiplier(dt);
    }
    return prod;
  }

  double evalTrial(double x_new, double dt) {
    if (nStages == 0) {
      y_t = x_new;
      D_t = 1.0;
      return y_t;
    }
    double sig = x_new;
    double prod = 1.0;
    for (int k = 0; k < nStages; ++k) {
      sig = stages[k].evalTrial(sig, dt);
      prod *= stages[k].D_t;
    }
    y_t = sig;
    D_t = prod;
    return y_t;
  }

  void commit() {
    for (int k = 0; k < nStages; ++k) stages[k].commit();
    y_c = y_t;
    D_c = D_t;
  }

  void revertToLastCommit() {
    for (int k = 0; k < nStages; ++k) stages[k].revertToLastCommit();
    y_t = y_c;
    D_t = D_c;
  }

  void revertToStart() {
    for (int k = 0; k < nStages; ++k) stages[k].revertToStart();
    const double d0 = instantaneousMultiplier(0.0);
    y_c = 0.0;
    D_c = d0;
    y_t = 0.0;
    D_t = d0;
  }

  int numStateDoubles() const {
    return nStages * RCFilter1D::numStateDoubles() + 4;
  }

  void packState(double*& ptr) const {
    for (int k = 0; k < nStages; ++k) stages[k].packState(ptr);
    *ptr++ = y_c;
    *ptr++ = D_c;
    *ptr++ = y_t;
    *ptr++ = D_t;
  }

  void unpackState(const double*& ptr) {
    for (int k = 0; k < nStages; ++k) stages[k].unpackState(ptr);
    y_c = *ptr++;
    D_c = *ptr++;
    y_t = *ptr++;
    D_t = *ptr++;
  }
};

// =========================================================================== //
// 3. Vector Filter Cascade (RCVectorCascade<MaxDim, MaxStages>)
//    Applies an identical K-stage scalar cascade independently to each
//    component d = 0..dim-1 of a vector signal.
// =========================================================================== //
template <int MaxDim, int MaxStages>
struct RCVectorCascade {
  int dim = 0;
  std::array<RCFilterCascade<MaxStages>, MaxDim> comp{};
  std::array<double, MaxDim> y_c{};
  std::array<double, MaxDim> y_t{};
  double D_c = 1.0;
  double D_t = 1.0;

  void configureIdentity(int d) {
    assert(d >= 0 && d <= MaxDim);
    dim = d;
    for (int c = 0; c < dim; ++c) {
      comp[c].setIdentity();
      y_c[c] = 0.0;
      y_t[c] = 0.0;
    }
    D_c = 1.0;
    D_t = 1.0;
  }

  void configure(int d, int kStages, const RCFilter1D* stg) {
    assert(d >= 0 && d <= MaxDim);
    dim = d;
    for (int c = 0; c < dim; ++c) {
      comp[c].configure(kStages, stg);
      y_c[c] = 0.0;
      y_t[c] = 0.0;
    }
    const double d0 = (dim > 0) ? comp[0].instantaneousMultiplier(0.0) : 1.0;
    D_c = d0;
    D_t = d0;
  }

  void reconfigureCoeffs(int d, int kStages, const RCFilter1D* stg) {
    assert(d >= 0 && d <= MaxDim);
    dim = d;
    for (int c = 0; c < dim; ++c) {
      comp[c].reconfigureCoeffs(kStages, stg);
    }
  }

  double instantaneousMultiplier(double dt) const {
    return (dim > 0) ? comp[0].instantaneousMultiplier(dt) : 1.0;
  }

  void evalTrial(const double* x_new, double dt, double* y_out = nullptr) {
    D_t = 1.0;
    for (int c = 0; c < dim; ++c) {
      y_t[c] = comp[c].evalTrial(x_new[c], dt);
      D_t = comp[c].D_t;
      if (y_out) y_out[c] = y_t[c];
    }
  }

  void commit() {
    for (int c = 0; c < dim; ++c) {
      comp[c].commit();
      y_c[c] = y_t[c];
    }
    D_c = D_t;
  }

  void revertToLastCommit() {
    for (int c = 0; c < dim; ++c) {
      comp[c].revertToLastCommit();
      y_t[c] = y_c[c];
    }
    D_t = D_c;
  }

  void revertToStart() {
    for (int c = 0; c < dim; ++c) {
      comp[c].revertToStart();
      y_c[c] = 0.0;
      y_t[c] = 0.0;
    }
    const double d0 = (dim > 0) ? comp[0].D_t : 1.0;
    D_c = d0;
    D_t = d0;
  }

  int numStateDoubles() const {
    int total = 2; // D_c, D_t
    for (int c = 0; c < dim; ++c) total += comp[c].numStateDoubles();
    return total;
  }

  void packState(double*& ptr) const {
    for (int c = 0; c < dim; ++c) comp[c].packState(ptr);
    *ptr++ = D_c;
    *ptr++ = D_t;
  }

  void unpackState(const double*& ptr) {
    for (int c = 0; c < dim; ++c) {
      comp[c].unpackState(ptr);
      y_c[c] = comp[c].y_c;
      y_t[c] = comp[c].y_t;
    }
    D_c = *ptr++;
    D_t = *ptr++;
  }
};

// =========================================================================== //
// 4. Material, Grading, and Geometry Configurations
// =========================================================================== //
struct CPMLMaterialParams {
  double E = 0.0;
  double nu = 0.0;
  double rho = 0.0;

  double mu = 0.0;
  double lambda = 0.0;
  double Vp = 0.0;
  double C[3][3][3][3]{};

  CPMLMaterialParams() = default;
  CPMLMaterialParams(double E_, double nu_, double rho_) {
    init(E_, nu_, rho_);
  }

  bool isValid() const noexcept {
    return std::isfinite(E) && std::isfinite(nu) && std::isfinite(rho) &&
           (E > 0.0) && (nu > -1.0 && nu < 0.499999) && (rho > 0.0);
  }

  void init(double E_, double nu_, double rho_) {
    E = E_;
    nu = nu_;
    rho = rho_;
    if (!isValid()) {
      mu = 0.0;
      lambda = 0.0;
      Vp = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l)
              C[i][j][k][l] = 0.0;
      return;
    }
    mu = E / (2.0 * (1.0 + nu));
    lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mod_p = std::max(0.0, (lambda + 2.0 * mu) / rho);
    Vp = std::sqrt(mod_p);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l) {
            const double d_ij = (i == j) ? 1.0 : 0.0;
            const double d_kl = (k == l) ? 1.0 : 0.0;
            const double d_ik = (i == k) ? 1.0 : 0.0;
            const double d_jl = (j == l) ? 1.0 : 0.0;
            const double d_il = (i == l) ? 1.0 : 0.0;
            const double d_jk = (j == k) ? 1.0 : 0.0;
            C[i][j][k][l] = lambda * d_ij * d_kl + mu * (d_ik * d_jl + d_il * d_jk);
          }
  }
};

struct CPMLGradingParams {
  double m_deg = 2.0;
  double R0 = 1.0e-6;
  double kappa_max = 1.0;
  double alpha_max = 0.0;

  CPMLGradingParams() = default;
  CPMLGradingParams(double m_deg_, double R0_, double kappa_max_, double alpha_max_)
      : m_deg(m_deg_), R0(R0_), kappa_max(kappa_max_), alpha_max(alpha_max_) {}

  bool isValid() const noexcept {
    return std::isfinite(m_deg) && std::isfinite(R0) &&
           std::isfinite(kappa_max) && std::isfinite(alpha_max) &&
           (m_deg >= 0.0) && (R0 > 0.0 && R0 < 1.0) &&
           (kappa_max >= 1.0) && (alpha_max >= 0.0);
  }
};

enum class CPMLGeomMode {
  Box = 0,
  Plane = 1
};

struct CPMLGeomConfig {
  CPMLGeomMode mode = CPMLGeomMode::Box;
  int ndm = 2;
  bool valid_normal = true;

  // Box mode parameters
  double box_min[3] = {0.0, 0.0, 0.0};
  double box_max[3] = {0.0, 0.0, 0.0};
  double L_box[3]   = {1.0, 1.0, 1.0};

  // Plane mode parameters
  double x0[3] = {0.0, 0.0, 0.0};
  double n_hat[3] = {1.0, 0.0, 0.0};
  double L_pml = 1.0;

  // Local orthonormal basis R in SO(ndm), stored row-major R[row][col]
  double R[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
  };

  bool isValid() const noexcept {
    if (ndm != 2 && ndm != 3) return false;
    if (mode == CPMLGeomMode::Box) {
      for (int m = 0; m < ndm; ++m) {
        if (!std::isfinite(box_min[m]) || !std::isfinite(box_max[m]) || !std::isfinite(L_box[m]))
          return false;
        if (!(L_box[m] > 0.0) || (box_max[m] < box_min[m]))
          return false;
      }
      return true;
    } else {
      if (!valid_normal || !std::isfinite(L_pml) || !(L_pml > 0.0))
        return false;
      for (int m = 0; m < ndm; ++m) {
        if (!std::isfinite(x0[m]) || !std::isfinite(n_hat[m]))
          return false;
      }
      return true;
    }
  }

  static CPMLGeomConfig makeBox2D(double xmin, double xmax,
                                  double ymin, double ymax,
                                  double Lx, double Ly) {
    CPMLGeomConfig cfg;
    cfg.mode = CPMLGeomMode::Box;
    cfg.ndm = 2;
    cfg.box_min[0] = xmin; cfg.box_max[0] = xmax; cfg.L_box[0] = Lx;
    cfg.box_min[1] = ymin; cfg.box_max[1] = ymax; cfg.L_box[1] = Ly;
    return cfg;
  }

  static CPMLGeomConfig makeBox3D(double xmin, double xmax,
                                  double ymin, double ymax,
                                  double zmin, double zmax,
                                  double Lx, double Ly, double Lz) {
    CPMLGeomConfig cfg;
    cfg.mode = CPMLGeomMode::Box;
    cfg.ndm = 3;
    cfg.box_min[0] = xmin; cfg.box_max[0] = xmax; cfg.L_box[0] = Lx;
    cfg.box_min[1] = ymin; cfg.box_max[1] = ymax; cfg.L_box[1] = Ly;
    cfg.box_min[2] = zmin; cfg.box_max[2] = zmax; cfg.L_box[2] = Lz;
    return cfg;
  }

  static CPMLGeomConfig makePlane2D(double x0_, double y0_,
                                    double nx_, double ny_,
                                    double L_pml_) {
    CPMLGeomConfig cfg;
    cfg.mode = CPMLGeomMode::Plane;
    cfg.ndm = 2;
    cfg.x0[0] = x0_;
    cfg.x0[1] = y0_;
    const double norm = std::sqrt(nx_ * nx_ + ny_ * ny_);
    cfg.valid_normal = std::isfinite(norm) && (norm > 1.0e-15);
    const double nx = cfg.valid_normal ? (nx_ / norm) : 1.0;
    const double ny = cfg.valid_normal ? (ny_ / norm) : 0.0;
    cfg.n_hat[0] = nx;
    cfg.n_hat[1] = ny;
    cfg.L_pml = L_pml_;
    // Row 0: (nx, ny), Row 1: (-ny, nx)
    cfg.R[0][0] =  nx; cfg.R[0][1] = ny; cfg.R[0][2] = 0.0;
    cfg.R[1][0] = -ny; cfg.R[1][1] = nx; cfg.R[1][2] = 0.0;
    cfg.R[2][0] = 0.0; cfg.R[2][1] = 0.0; cfg.R[2][2] = 1.0;
    return cfg;
  }

  static CPMLGeomConfig makePlane3D(double x0_, double y0_, double z0_,
                                    double nx_, double ny_, double nz_,
                                    double L_pml_) {
    CPMLGeomConfig cfg;
    cfg.mode = CPMLGeomMode::Plane;
    cfg.ndm = 3;
    cfg.x0[0] = x0_;
    cfg.x0[1] = y0_;
    cfg.x0[2] = z0_;
    const double norm = std::sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_);
    cfg.valid_normal = std::isfinite(norm) && (norm > 1.0e-15);
    const double nx = cfg.valid_normal ? (nx_ / norm) : 1.0;
    const double ny = cfg.valid_normal ? (ny_ / norm) : 0.0;
    const double nz = cfg.valid_normal ? (nz_ / norm) : 0.0;
    cfg.n_hat[0] = nx;
    cfg.n_hat[1] = ny;
    cfg.n_hat[2] = nz;
    cfg.L_pml = L_pml_;

    // Row 0: e1' = n_hat
    cfg.R[0][0] = nx;
    cfg.R[0][1] = ny;
    cfg.R[0][2] = nz;

    // Helper axis a = (1,0,0) if |nx| < 0.9 else (0,1,0)
    const double ax = (std::fabs(nx) < 0.9) ? 1.0 : 0.0;
    const double ay = (std::fabs(nx) < 0.9) ? 0.0 : 1.0;
    const double az = 0.0;

    // Row 1: e2' = (n_hat x a) / ||n_hat x a||
    double c1x = ny * az - nz * ay;
    double c1y = nz * ax - nx * az;
    double c1z = nx * ay - ny * ax;
    const double c1norm = std::sqrt(c1x * c1x + c1y * c1y + c1z * c1z);
    c1x /= c1norm;
    c1y /= c1norm;
    c1z /= c1norm;
    cfg.R[1][0] = c1x;
    cfg.R[1][1] = c1y;
    cfg.R[1][2] = c1z;

    // Row 2: e3' = e1' x e2'
    cfg.R[2][0] = ny * c1z - nz * c1y;
    cfg.R[2][1] = nz * c1x - nx * c1z;
    cfg.R[2][2] = nx * c1y - ny * c1x;
    return cfg;
  }

  void evalDepthAndLength(const double x_gp[3], double xi_depth[3], double L_dir[3]) const {
    for (int m = 0; m < 3; ++m) {
      xi_depth[m] = 0.0;
      L_dir[m] = 1.0;
    }
    if (mode == CPMLGeomMode::Box) {
      for (int m = 0; m < ndm; ++m) {
        const double d_low  = box_min[m] - x_gp[m];
        const double d_high = x_gp[m] - box_max[m];
        xi_depth[m] = std::max(0.0, std::max(d_low, d_high));
        L_dir[m] = L_box[m];
      }
    } else {
      double dot = 0.0;
      for (int m = 0; m < ndm; ++m) {
        dot += n_hat[m] * (x_gp[m] - x0[m]);
        L_dir[m] = L_pml;
      }
      xi_depth[0] = std::max(0.0, dot);
      for (int m = 1; m < ndm; ++m) {
        xi_depth[m] = 0.0;
      }
    }
  }
};

struct CPMLBiotParams {
  bool enabled = false;
  double Kf = 0.0;
  double poro = 0.0;
  double rhoF = 0.0;
  double biotAlpha = 1.0;
  double Ks = 0.0;
  int pOrder = 0;       // 0 = equal, 1 = linear/Taylor-Hood
  int formulation = 0;  // 0 = std, 1 = bbar
  int stabMode = 1;     // 0 = off, 1 = auto, 2 = manual
  double stabValue = 0.25;
  double kbar[3] = {0.0, 0.0, 0.0};

  CPMLBiotParams() = default;

  double oneOverQbar() const noexcept {
    if (!enabled || !(Kf > 0.0)) return 0.0;
    double val = poro / Kf;
    if (Ks > 0.0) {
      val += (biotAlpha - poro) / Ks;
    }
    return val;
  }

  double stabAlpha(double h, double lambda_plus_2mu) const noexcept {
    if (!enabled || pOrder == 1 || stabMode == 0) return 0.0;
    if (stabMode == 1) {
      return (lambda_plus_2mu > 0.0) ? (stabValue * (h * h) / lambda_plus_2mu) : 0.0;
    }
    if (stabMode == 2) {
      return stabValue;
    }
    return 0.0;
  }

  double undrainedVp(const CPMLMaterialParams& mat) const noexcept {
    if (!enabled) return mat.Vp;
    if (!(mat.rho > 0.0)) return 0.0;
    const double invQ = oneOverQbar();
    const double M_u = (mat.lambda + 2.0 * mat.mu) +
                       ((invQ > 0.0) ? ((biotAlpha * biotAlpha) / invQ) : 0.0);
    return std::sqrt(std::max(0.0, M_u / mat.rho));
  }

  bool isValid(int ndm, CPMLTopology topo) const noexcept {
    if (!enabled) return true;
    if (ndm != 2 && ndm != 3) return false;
    if (!std::isfinite(Kf) || !(Kf > 0.0)) return false;
    if (!std::isfinite(poro) || !(poro >= 0.0 && poro < 1.0)) return false;
    if (!std::isfinite(rhoF) || !(rhoF >= 0.0)) return false;
    if (!std::isfinite(biotAlpha) || !(biotAlpha >= poro && biotAlpha <= 1.0)) return false;
    if (!std::isfinite(Ks) || !(Ks >= 0.0)) return false;
    if (pOrder != 0 && pOrder != 1) return false;
    if (pOrder == 1 && shapes::numVertices(topo) >= shapes::topologyNumNodes(topo)) return false;
    if (formulation != 0 && formulation != 1) return false;
    if (stabMode != 0 && stabMode != 1 && stabMode != 2) return false;
    if (!std::isfinite(stabValue) || !(stabValue >= 0.0)) return false;
    for (int m = 0; m < ndm; ++m) {
      if (!std::isfinite(kbar[m]) || kbar[m] < 0.0) return false;
    }
    const double invQ = oneOverQbar();
    return std::isfinite(invQ) && (invQ >= 0.0);
  }
};

// =========================================================================== //
// 5. Per-Gauss-Point State & Filter Cascades (CPMLGaussState)
// =========================================================================== //
struct CPMLGaussState {
  double x_gp[3] = {0.0, 0.0, 0.0};
  double detJ = 0.0;
  double w_vol = 0.0;
  double N[shapes::MAX_NODES]{};
  double dNdx[shapes::MAX_NODES][3]{};
  double Phi_local[shapes::MAX_NODES][3]{};

  // Pressure basis and B-bar / coupling local operators
  int nP = 0;
  double Np[shapes::MAX_NODES]{};
  double dNpdx[shapes::MAX_NODES][3]{};
  double PhiP_local[shapes::MAX_NODES][3]{};
  double Phi_tilde_local[shapes::MAX_NODES][3]{};
  double B_local[3][3][shapes::MAX_NODES][3]{};

  // Stretching parameters (in 2D, index 2 is fixed to d=0, kappa=1, alpha=0)
  double d[3]     = {0.0, 0.0, 0.0};
  double kappa[3] = {1.0, 1.0, 1.0};
  double alpha[3] = {0.0, 0.0, 0.0};
  bool has_d[3]   = {false, false, false};
  bool has_any_d  = false;
  bool has_d01    = false;
  bool has_d02    = false;
  bool has_d12    = false;
  bool has_d012   = false;

  // Stiffness filter cascades F_{jl}(w) applied to local gradient vector k -> G'_{kl}
  RCVectorCascade<3, 3> F_jl[3][3]{};

  // Inertial velocity filters g_0, g_1, g_2 applied to global velocity v_g
  RCVectorCascade<3, 1> g_vel[3]{};

  // Inertial displacement cascades applied to global displacement u_g
  RCVectorCascade<3, 2> g1_g0{};
  RCVectorCascade<3, 2> g2_g0{};
  RCVectorCascade<3, 2> g2_g1{};
  RCVectorCascade<3, 3> phi2_g1_g0{};

  // Biot u-p filter cascades
  RCFilterCascade<3> F_p[3][3]{};  // Darcy permeability cascades F_{jl}^{(p)}
  RCFilterCascade<2> E_p[3]{};     // Coupling pressure cascades E_j^{(p)} on u-rows
  RCFilterCascade<2> E_v[3]{};     // Coupling velocity cascades E_j^{(\dot{u})} on p-rows
  RCFilterCascade<3> J_s{};        // Storage cascade J_s^{(p)} on p-rows

  void refreshActiveFlags() noexcept {
    has_d[0] = (d[0] != 0.0);
    has_d[1] = (d[1] != 0.0);
    has_d[2] = (d[2] != 0.0);
    has_any_d = has_d[0] || has_d[1] || has_d[2];
    has_d01 = has_d[0] && has_d[1];
    has_d02 = has_d[0] && has_d[2];
    has_d12 = has_d[1] && has_d[2];
    has_d012 = has_d01 && has_d[2];
  }

  void init(int ndm, int nNodes,
            const double x_gp_[3], double detJ_, double w_vol_,
            const double* N_, const double* dNdx_,
            const CPMLMaterialParams& mat,
            const CPMLGradingParams& grading,
            const CPMLGeomConfig& geom,
            const CPMLBiotParams& biot = CPMLBiotParams{},
            int nP_ = 0,
            const double* Np_ = nullptr,
            const double* dNpdx_ = nullptr) {
    detJ = detJ_;
    w_vol = w_vol_;
    for (int m = 0; m < 3; ++m) x_gp[m] = (m < ndm) ? x_gp_[m] : 0.0;

    for (int a = 0; a < nNodes; ++a) {
      N[a] = N_[a];
      for (int m = 0; m < 3; ++m) {
        dNdx[a][m] = (m < ndm) ? dNdx_[a * ndm + m] : 0.0;
      }
      for (int j = 0; j < ndm; ++j) {
        double s = 0.0;
        for (int p = 0; p < ndm; ++p) {
          s += geom.R[j][p] * dNdx[a][p];
        }
        Phi_local[a][j] = s;
        Phi_tilde_local[a][j] = s;
      }
      for (int j = ndm; j < 3; ++j) {
        Phi_local[a][j] = 0.0;
        Phi_tilde_local[a][j] = 0.0;
      }
    }

    nP = nP_;
    for (int b = 0; b < shapes::MAX_NODES; ++b) {
      Np[b] = (b < nP && Np_ != nullptr) ? Np_[b] : 0.0;
      for (int m = 0; m < 3; ++m) {
        dNpdx[b][m] = (b < nP && m < ndm && dNpdx_ != nullptr) ? dNpdx_[b * ndm + m] : 0.0;
      }
      for (int j = 0; j < ndm; ++j) {
        double s = 0.0;
        if (b < nP) {
          for (int p = 0; p < ndm; ++p) {
            s += geom.R[j][p] * dNpdx[b][p];
          }
        }
        PhiP_local[b][j] = s;
      }
      for (int j = ndm; j < 3; ++j) {
        PhiP_local[b][j] = 0.0;
      }
    }

    // Evaluate polynomial grading along each local axis m = 0..ndm-1
    double xi_depth[3], L_dir[3];
    geom.evalDepthAndLength(x_gp, xi_depth, L_dir);
    const double Vp_grade = biot.enabled ? biot.undrainedVp(mat) : mat.Vp;

    for (int m = 0; m < 3; ++m) {
      d[m] = 0.0;
      kappa[m] = 1.0;
      alpha[m] = 0.0;
      if (m < ndm && L_dir[m] > 0.0) {
        const double xi_bar = std::min(1.0, std::max(0.0, xi_depth[m] / L_dir[m]));
        if (xi_bar > 0.0) {
          const double p = std::pow(xi_bar, grading.m_deg);
          const double d0 = -((grading.m_deg + 1.0) * Vp_grade * std::log(grading.R0)) / (2.0 * L_dir[m]);
          d[m]     = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_DAMPING) ? 0.0 : (d0 * p);
          kappa[m] = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_KAPPA)   ? 1.0 : (1.0 + (grading.kappa_max - 1.0) * p);
          alpha[m] = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_ALPHA)   ? 0.0 : (grading.alpha_max * (1.0 - xi_bar));
        }
      }
    }
    refreshActiveFlags();

    // Elementary 1D filters for m = 0, 1, 2
    RCFilter1D s_fwd[3], s_inv[3], g_flt[3], phi_flt[3];
    for (int m = 0; m < 3; ++m) {
      s_fwd[m]   = RCFilter1D(kappa[m], d[m], alpha[m]);
      s_inv[m]   = RCFilter1D(1.0 / kappa[m],
                              -d[m] / (kappa[m] * kappa[m]),
                              alpha[m] + d[m] / kappa[m]);
      g_flt[m]   = RCFilter1D(d[m], -alpha[m] * d[m], alpha[m]);
      phi_flt[m] = RCFilter1D(0.0, d[m], alpha[m]);
    }

    // Configure stiffness filter cascades F_{jl}(w) and Biot cascades
    if (ndm == 2) {
      // F_{00} = s_1 o s_0^{-1} (Stage 1 = s_0^{-1}, Stage 2 = s_1)
      const RCFilter1D stg00[2] = {s_inv[0], s_fwd[1]};
      F_jl[0][0].configure(2, 2, stg00);
      F_p[0][0].configure(2, stg00);

      // F_{11} = s_0 o s_1^{-1} (Stage 1 = s_1^{-1}, Stage 2 = s_0)
      const RCFilter1D stg11[2] = {s_inv[1], s_fwd[0]};
      F_jl[1][1].configure(2, 2, stg11);
      F_p[1][1].configure(2, stg11);

      // F_{01} = F_{10} = 1 (0 stages, identity)
      F_jl[0][1].configureIdentity(2);
      F_jl[1][0].configureIdentity(2);
      F_p[0][1].setIdentity();
      F_p[1][0].setIdentity();

      // Coupling E_0 = [s_1], E_1 = [s_0]
      const RCFilter1D stgE0[1] = {s_fwd[1]};
      const RCFilter1D stgE1[1] = {s_fwd[0]};
      E_p[0].configure(1, stgE0);
      E_v[0].configure(1, stgE0);
      E_p[1].configure(1, stgE1);
      E_v[1].configure(1, stgE1);
      E_p[2].setIdentity();
      E_v[2].setIdentity();

      // Storage J_s = [s_0, s_1]
      const RCFilter1D stgJs[2] = {s_fwd[0], s_fwd[1]};
      J_s.configure(2, stgJs);
    } else {
      // 3D stiffness cascades
      // F_{00} = s_2 o s_1 o s_0^{-1}
      const RCFilter1D stg00[3] = {s_inv[0], s_fwd[1], s_fwd[2]};
      F_jl[0][0].configure(3, 3, stg00);
      F_p[0][0].configure(3, stg00);

      // F_{11} = s_2 o s_0 o s_1^{-1}
      const RCFilter1D stg11[3] = {s_inv[1], s_fwd[0], s_fwd[2]};
      F_jl[1][1].configure(3, 3, stg11);
      F_p[1][1].configure(3, stg11);

      // F_{22} = s_1 o s_0 o s_2^{-1}
      const RCFilter1D stg22[3] = {s_inv[2], s_fwd[0], s_fwd[1]};
      F_jl[2][2].configure(3, 3, stg22);
      F_p[2][2].configure(3, stg22);

      // Off-diagonal 1-stage filters
      const RCFilter1D stg01[1] = {s_fwd[2]};
      F_jl[0][1].configure(3, 1, stg01);
      F_jl[1][0].configure(3, 1, stg01);
      F_p[0][1].configure(1, stg01);
      F_p[1][0].configure(1, stg01);

      const RCFilter1D stg02[1] = {s_fwd[1]};
      F_jl[0][2].configure(3, 1, stg02);
      F_jl[2][0].configure(3, 1, stg02);
      F_p[0][2].configure(1, stg02);
      F_p[2][0].configure(1, stg02);

      const RCFilter1D stg12[1] = {s_fwd[0]};
      F_jl[1][2].configure(3, 1, stg12);
      F_jl[2][1].configure(3, 1, stg12);
      F_p[1][2].configure(1, stg12);
      F_p[2][1].configure(1, stg12);

      // Coupling E_0 = [s_1, s_2], E_1 = [s_0, s_2], E_2 = [s_0, s_1]
      const RCFilter1D stgE0[2] = {s_fwd[1], s_fwd[2]};
      const RCFilter1D stgE1[2] = {s_fwd[0], s_fwd[2]};
      const RCFilter1D stgE2[2] = {s_fwd[0], s_fwd[1]};
      E_p[0].configure(2, stgE0);
      E_v[0].configure(2, stgE0);
      E_p[1].configure(2, stgE1);
      E_v[1].configure(2, stgE1);
      E_p[2].configure(2, stgE2);
      E_v[2].configure(2, stgE2);

      // Storage J_s = [s_0, s_1, s_2]
      const RCFilter1D stgJs[3] = {s_fwd[0], s_fwd[1], s_fwd[2]};
      J_s.configure(3, stgJs);
    }

    // Configure velocity filters g_0, g_1, g_2
    for (int m = 0; m < 3; ++m) {
      const RCFilter1D stg_v[1] = {g_flt[m]};
      g_vel[m].configure(ndm, 1, stg_v);
    }

    // Configure displacement inertial cascades
    const RCFilter1D stg_g1_g0[2] = {g_flt[0], g_flt[1]};
    g1_g0.configure(ndm, 2, stg_g1_g0);

    const RCFilter1D stg_g2_g0[2] = {g_flt[0], g_flt[2]};
    g2_g0.configure(ndm, 2, stg_g2_g0);

    const RCFilter1D stg_g2_g1[2] = {g_flt[1], g_flt[2]};
    g2_g1.configure(ndm, 2, stg_g2_g1);

    const RCFilter1D stg_phi2_g1_g0[3] = {g_flt[0], g_flt[1], phi_flt[2]};
    phi2_g1_g0.configure(ndm, 3, stg_phi2_g1_g0);
  }

  void updateMaterialParams(int ndm,
                            const CPMLMaterialParams& mat,
                            const CPMLGradingParams& grading,
                            const CPMLGeomConfig& geom,
                            const CPMLBiotParams& biot = CPMLBiotParams{}) {
    double xi_depth[3], L_dir[3];
    geom.evalDepthAndLength(x_gp, xi_depth, L_dir);
    const double Vp_grade = biot.enabled ? biot.undrainedVp(mat) : mat.Vp;

    for (int m = 0; m < 3; ++m) {
      d[m] = 0.0;
      kappa[m] = 1.0;
      alpha[m] = 0.0;
      if (m < ndm && L_dir[m] > 0.0) {
        const double xi_bar = std::min(1.0, std::max(0.0, xi_depth[m] / L_dir[m]));
        if (xi_bar > 0.0) {
          const double p = std::pow(xi_bar, grading.m_deg);
          const double d0 = -((grading.m_deg + 1.0) * Vp_grade * std::log(grading.R0)) / (2.0 * L_dir[m]);
          d[m]     = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_DAMPING) ? 0.0 : (d0 * p);
          kappa[m] = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_KAPPA)   ? 1.0 : (1.0 + (grading.kappa_max - 1.0) * p);
          alpha[m] = LADRUNO_CPML_MUT_ACTIVE(CPML_NO_ALPHA)   ? 0.0 : (grading.alpha_max * (1.0 - xi_bar));
        }
      }
    }
    refreshActiveFlags();

    RCFilter1D s_fwd[3], s_inv[3], g_flt[3], phi_flt[3];
    for (int m = 0; m < 3; ++m) {
      s_fwd[m]   = RCFilter1D(kappa[m], d[m], alpha[m]);
      s_inv[m]   = RCFilter1D(1.0 / kappa[m],
                              -d[m] / (kappa[m] * kappa[m]),
                              alpha[m] + d[m] / kappa[m]);
      g_flt[m]   = RCFilter1D(d[m], -alpha[m] * d[m], alpha[m]);
      phi_flt[m] = RCFilter1D(0.0, d[m], alpha[m]);
    }

    if (ndm == 2) {
      const RCFilter1D stg00[2] = {s_inv[0], s_fwd[1]};
      F_jl[0][0].reconfigureCoeffs(2, 2, stg00);
      F_p[0][0].reconfigureCoeffs(2, stg00);
      const RCFilter1D stg11[2] = {s_inv[1], s_fwd[0]};
      F_jl[1][1].reconfigureCoeffs(2, 2, stg11);
      F_p[1][1].reconfigureCoeffs(2, stg11);

      const RCFilter1D stgE0[1] = {s_fwd[1]};
      const RCFilter1D stgE1[1] = {s_fwd[0]};
      E_p[0].reconfigureCoeffs(1, stgE0);
      E_v[0].reconfigureCoeffs(1, stgE0);
      E_p[1].reconfigureCoeffs(1, stgE1);
      E_v[1].reconfigureCoeffs(1, stgE1);

      const RCFilter1D stgJs[2] = {s_fwd[0], s_fwd[1]};
      J_s.reconfigureCoeffs(2, stgJs);
    } else {
      const RCFilter1D stg00[3] = {s_inv[0], s_fwd[1], s_fwd[2]};
      F_jl[0][0].reconfigureCoeffs(3, 3, stg00);
      F_p[0][0].reconfigureCoeffs(3, stg00);
      const RCFilter1D stg11[3] = {s_inv[1], s_fwd[0], s_fwd[2]};
      F_jl[1][1].reconfigureCoeffs(3, 3, stg11);
      F_p[1][1].reconfigureCoeffs(3, stg11);
      const RCFilter1D stg22[3] = {s_inv[2], s_fwd[0], s_fwd[1]};
      F_jl[2][2].reconfigureCoeffs(3, 3, stg22);
      F_p[2][2].reconfigureCoeffs(3, stg22);

      const RCFilter1D stg01[1] = {s_fwd[2]};
      F_jl[0][1].reconfigureCoeffs(3, 1, stg01);
      F_jl[1][0].reconfigureCoeffs(3, 1, stg01);
      F_p[0][1].reconfigureCoeffs(1, stg01);
      F_p[1][0].reconfigureCoeffs(1, stg01);
      const RCFilter1D stg02[1] = {s_fwd[1]};
      F_jl[0][2].reconfigureCoeffs(3, 1, stg02);
      F_jl[2][0].reconfigureCoeffs(3, 1, stg02);
      F_p[0][2].reconfigureCoeffs(1, stg02);
      F_p[2][0].reconfigureCoeffs(1, stg02);
      const RCFilter1D stg12[1] = {s_fwd[0]};
      F_jl[1][2].reconfigureCoeffs(3, 1, stg12);
      F_jl[2][1].reconfigureCoeffs(3, 1, stg12);
      F_p[1][2].reconfigureCoeffs(1, stg12);
      F_p[2][1].reconfigureCoeffs(1, stg12);

      const RCFilter1D stgE0[2] = {s_fwd[1], s_fwd[2]};
      const RCFilter1D stgE1[2] = {s_fwd[0], s_fwd[2]};
      const RCFilter1D stgE2[2] = {s_fwd[0], s_fwd[1]};
      E_p[0].reconfigureCoeffs(2, stgE0);
      E_v[0].reconfigureCoeffs(2, stgE0);
      E_p[1].reconfigureCoeffs(2, stgE1);
      E_v[1].reconfigureCoeffs(2, stgE1);
      E_p[2].reconfigureCoeffs(2, stgE2);
      E_v[2].reconfigureCoeffs(2, stgE2);

      const RCFilter1D stgJs[3] = {s_fwd[0], s_fwd[1], s_fwd[2]};
      J_s.reconfigureCoeffs(3, stgJs);
    }

    for (int m = 0; m < 3; ++m) {
      const RCFilter1D stg_v[1] = {g_flt[m]};
      g_vel[m].reconfigureCoeffs(ndm, 1, stg_v);
    }
    const RCFilter1D stg_g1_g0[2] = {g_flt[0], g_flt[1]};
    g1_g0.reconfigureCoeffs(ndm, 2, stg_g1_g0);
    const RCFilter1D stg_g2_g0[2] = {g_flt[0], g_flt[2]};
    g2_g0.reconfigureCoeffs(ndm, 2, stg_g2_g0);
    const RCFilter1D stg_g2_g1[2] = {g_flt[1], g_flt[2]};
    g2_g1.reconfigureCoeffs(ndm, 2, stg_g2_g1);
    const RCFilter1D stg_phi2_g1_g0[3] = {g_flt[0], g_flt[1], phi_flt[2]};
    phi2_g1_g0.reconfigureCoeffs(ndm, 3, stg_phi2_g1_g0);
  }

  void commit(int ndm, bool biot_enabled = false) {
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        F_jl[j][l].commit();
    if (has_any_d) {
      for (int m = 0; m < 3; ++m) {
        if (has_d[m]) g_vel[m].commit();
      }
      if (has_d01)  g1_g0.commit();
      if (has_d02)  g2_g0.commit();
      if (has_d12)  g2_g1.commit();
      if (has_d012) phi2_g1_g0.commit();
    }
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) F_p[j][l].commit();
        E_p[j].commit();
        E_v[j].commit();
      }
      J_s.commit();
    }
  }

  void revertToLastCommit(int ndm, bool biot_enabled = false) {
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        F_jl[j][l].revertToLastCommit();
    if (has_any_d) {
      for (int m = 0; m < 3; ++m) {
        if (has_d[m]) g_vel[m].revertToLastCommit();
      }
      if (has_d01)  g1_g0.revertToLastCommit();
      if (has_d02)  g2_g0.revertToLastCommit();
      if (has_d12)  g2_g1.revertToLastCommit();
      if (has_d012) phi2_g1_g0.revertToLastCommit();
    }
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) F_p[j][l].revertToLastCommit();
        E_p[j].revertToLastCommit();
        E_v[j].revertToLastCommit();
      }
      J_s.revertToLastCommit();
    }
  }

  void revertToStart(int ndm, bool biot_enabled = false) {
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        F_jl[j][l].revertToStart();
    for (int m = 0; m < 3; ++m) g_vel[m].revertToStart();
    g1_g0.revertToStart();
    g2_g0.revertToStart();
    g2_g1.revertToStart();
    phi2_g1_g0.revertToStart();
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) F_p[j][l].revertToStart();
        E_p[j].revertToStart();
        E_v[j].revertToStart();
      }
      J_s.revertToStart();
    }
  }

  int numStateDoubles(int ndm, bool biot_enabled = false) const {
    int total = 0;
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        total += F_jl[j][l].numStateDoubles();
    for (int m = 0; m < 3; ++m) total += g_vel[m].numStateDoubles();
    total += g1_g0.numStateDoubles();
    total += g2_g0.numStateDoubles();
    total += g2_g1.numStateDoubles();
    total += phi2_g1_g0.numStateDoubles();
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) total += F_p[j][l].numStateDoubles();
        total += E_p[j].numStateDoubles();
        total += E_v[j].numStateDoubles();
      }
      total += J_s.numStateDoubles();
    }
    return total;
  }

  void packState(int ndm, double*& ptr, bool biot_enabled = false) const {
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        F_jl[j][l].packState(ptr);
    for (int m = 0; m < 3; ++m) g_vel[m].packState(ptr);
    g1_g0.packState(ptr);
    g2_g0.packState(ptr);
    g2_g1.packState(ptr);
    phi2_g1_g0.packState(ptr);
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) F_p[j][l].packState(ptr);
        E_p[j].packState(ptr);
        E_v[j].packState(ptr);
      }
      J_s.packState(ptr);
    }
  }

  void unpackState(int ndm, const double*& ptr, bool biot_enabled = false) {
    for (int j = 0; j < ndm; ++j)
      for (int l = 0; l < ndm; ++l)
        F_jl[j][l].unpackState(ptr);
    for (int m = 0; m < 3; ++m) g_vel[m].unpackState(ptr);
    g1_g0.unpackState(ptr);
    g2_g0.unpackState(ptr);
    g2_g1.unpackState(ptr);
    phi2_g1_g0.unpackState(ptr);
    if (biot_enabled) {
      for (int j = 0; j < ndm; ++j) {
        for (int l = 0; l < ndm; ++l) F_p[j][l].unpackState(ptr);
        E_p[j].unpackState(ptr);
        E_v[j].unpackState(ptr);
      }
      J_s.unpackState(ptr);
    }
  }
};

// =========================================================================== //
// 6. Complete Element C-PML Kernel (CPMLElementKernel)
// =========================================================================== //
class CPMLElementKernel {
public:
  CPMLElementKernel() = default;

  CPMLElementKernel(CPMLTopology topo, bool lumped, const double* coords,
                    const CPMLMaterialParams& mat,
                    const CPMLGradingParams& grading,
                    const CPMLGeomConfig& geom) {
    init(topo, lumped, coords, mat, grading, geom, CPMLBiotParams{});
  }

  CPMLElementKernel(CPMLTopology topo, bool lumped, const double* coords,
                    const CPMLMaterialParams& mat,
                    const CPMLGradingParams& grading,
                    const CPMLGeomConfig& geom,
                    const CPMLBiotParams& biot) {
    init(topo, lumped, coords, mat, grading, geom, biot);
  }

  bool init(CPMLTopology topo, bool lumped, const double* coords,
            const CPMLMaterialParams& mat,
            const CPMLGradingParams& grading,
            const CPMLGeomConfig& geom) {
    return init(topo, lumped, coords, mat, grading, geom, CPMLBiotParams{});
  }

  bool init(CPMLTopology topo, bool lumped, const double* coords,
            const CPMLMaterialParams& mat,
            const CPMLGradingParams& grading,
            const CPMLGeomConfig& geom,
            const CPMLBiotParams& biot) {
    topo_ = topo;
    lumped_ = lumped;
    static_mode_ = false;
    ndm_ = shapes::topologyNdm(topo);
    nNodes_ = shapes::topologyNumNodes(topo);
    nGP_ = shapes::topologyNumGP(topo);
    mat_ = mat;
    grading_ = grading;
    geom_ = geom;
    biot_ = biot;

    const int nCoordDoubles = nNodes_ * ndm_;
    if (coords == nullptr || nCoordDoubles <= 0) {
      valid_ = false;
      return false;
    }

    // Build per-node interleaved DOF map
    uOff_.assign(nNodes_, 0);
    pOff_.assign(nNodes_, -1);
    pCarrier_.assign(nNodes_, 0);
    if (!biot_.enabled) {
      nP_ = 0;
      for (int a = 0; a < nNodes_; ++a) {
        uOff_[a] = a * ndm_;
        pOff_[a] = -1;
      }
      nDOF_ = nNodes_ * ndm_;
    } else {
      nP_ = (biot_.pOrder == 1) ? shapes::numVertices(topo_) : nNodes_;
      int slot = 0;
      for (int a = 0; a < nNodes_; ++a) {
        uOff_[a] = slot;
        slot += ndm_;
        if (a < nP_) {
          pCarrier_[a] = 1;
          pOff_[a] = slot;
          slot += 1;
        } else {
          pCarrier_[a] = 0;
          pOff_[a] = -1;
        }
      }
      nDOF_ = slot;
    }

    coords_.assign(coords, coords + nCoordDoubles);
    gpStates_.resize(nGP_);

    const bool allowSignedTet10 =
        (topo_ == CPMLTopology::TET10 || topo_ == CPMLTopology::BTET10);
    int nPosDet = 0;
    int nNegDet = 0;
    bool pShapesOK = true;

    for (int g = 0; g < nGP_; ++g) {
      double N_g[shapes::MAX_NODES]{};
      double dNdx_g[shapes::MAX_NODES * 3]{};
      double x_gp[3] = {0.0, 0.0, 0.0};
      double detJ = 0.0;
      double w_vol = 0.0;
      shapes::evalTopologyGP(topo_, g, coords_.data(), N_g, dNdx_g, &detJ, x_gp, &w_vol);
      if (std::isfinite(detJ) && detJ != 0.0) {
        if (detJ > 0.0) nPosDet++;
        else if (detJ < 0.0) nNegDet++;
      }

      double Np_g[shapes::MAX_NODES]{};
      double dNpdx_g[shapes::MAX_NODES * 3]{};
      int nP_g = 0;
      if (biot_.enabled) {
        if (!shapes::evalShapeAndGradP(topo_, g, coords_.data(), biot_.pOrder,
                                       Np_g, dNpdx_g, nP_g) || nP_g != nP_) {
          pShapesOK = false;
        }
      }
      gpStates_[g].init(ndm_, nNodes_, x_gp, detJ, w_vol, N_g, dNdx_g,
                        mat_, grading_, geom_, biot_, nP_g, Np_g, dNpdx_g);
    }

    // Precompute element-mean local solid gradients Phi_bar_local_ and B_local
    double vol_elem = 0.0;
    for (int g = 0; g < nGP_; ++g) {
      vol_elem += gpStates_[g].w_vol;
    }
    for (int a = 0; a < shapes::MAX_NODES; ++a) {
      for (int j = 0; j < 3; ++j) {
        Phi_bar_local_[a][j] = 0.0;
      }
    }
    if (vol_elem > 0.0) {
      for (int g = 0; g < nGP_; ++g) {
        const double w_g = gpStates_[g].w_vol / vol_elem;
        for (int a = 0; a < nNodes_; ++a) {
          for (int j = 0; j < ndm_; ++j) {
            Phi_bar_local_[a][j] += w_g * gpStates_[g].Phi_local[a][j];
          }
        }
      }
    }
    for (int g = 0; g < nGP_; ++g) {
      auto& gs = gpStates_[g];
      for (int a = 0; a < nNodes_; ++a) {
        for (int j = 0; j < ndm_; ++j) {
          gs.Phi_tilde_local[a][j] =
              (biot_.enabled && biot_.formulation == 1) ? Phi_bar_local_[a][j] : gs.Phi_local[a][j];
        }
        for (int k = 0; k < ndm_; ++k) {
          for (int l = 0; l < ndm_; ++l) {
            for (int m = 0; m < ndm_; ++m) {
              double b_val = (k == m) ? gs.Phi_local[a][l] : 0.0;
              if (biot_.enabled && biot_.formulation == 1 && k == l) {
                b_val += (Phi_bar_local_[a][m] - gs.Phi_local[a][m]) / static_cast<double>(ndm_);
              }
              gs.B_local[k][l][a][m] = b_val;
            }
          }
        }
      }
    }

    const bool geomOK = allowSignedTet10
                            ? (nPosDet == nGP_ || nNegDet == nGP_)
                            : (nPosDet == nGP_);

    valid_ = mat_.isValid() && grading_.isValid() && geom_.isValid() &&
             biot_.isValid(ndm_, topo_) && geomOK && pShapesOK;

    rebuildStaticAndMassMatrices();

    R_int_.assign(nDOF_, 0.0);
    R_int_c_.assign(nDOF_, 0.0);
    R_disp_.assign(nDOF_, 0.0);
    R_disp_c_.assign(nDOF_, 0.0);
    R_vel_.assign(nDOF_, 0.0);
    R_vel_c_.assign(nDOF_, 0.0);
    K_eff_ = K_0_;
    K_eff_c_ = K_0_;
    C_eff_ = C_0_;
    C_eff_c_ = C_0_;
    P_ext_ = 0.0;
    P_ext_c_ = 0.0;
    last_dt_ = 0.0;
    last_dt_c_ = 0.0;
    hadPendingStateUnpacked_ = false;

    if (!pendingState_.empty() &&
        static_cast<int>(pendingState_.size()) == numStateDoubles()) {
      unpackState(pendingState_.data());
      pendingState_.clear();
      hadPendingStateUnpacked_ = true;
    }

    return valid_;
  }

  bool updateMaterialParams(const CPMLMaterialParams& mat, bool lumped) {
    if (!mat.isValid() || !valid_) {
      return false;
    }
    mat_ = mat;
    lumped_ = lumped;
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].updateMaterialParams(ndm_, mat_, grading_, geom_, biot_);
    }
    rebuildStaticAndMassMatrices();
    if (static_mode_ || last_dt_ <= 0.0) {
      K_eff_ = K_0_;
      K_eff_c_ = K_0_;
      C_eff_ = C_0_;
      C_eff_c_ = C_0_;
    }
    return true;
  }

  bool updateBiotParams(const CPMLBiotParams& biot) {
    if (!biot.isValid(ndm_, topo_) || !valid_ || biot.enabled != biot_.enabled ||
        biot.pOrder != biot_.pOrder || biot.formulation != biot_.formulation) {
      return false;
    }
    biot_ = biot;
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].updateMaterialParams(ndm_, mat_, grading_, geom_, biot_);
    }
    rebuildStaticAndMassMatrices();
    if (static_mode_ || last_dt_ <= 0.0) {
      K_eff_ = K_0_;
      K_eff_c_ = K_0_;
      C_eff_ = C_0_;
      C_eff_c_ = C_0_;
    }
    return true;
  }

  void setStaticMode(bool s) noexcept {
    static_mode_ = s;
    if (s) {
      tangent_cache_dt_ = -1.0;
    }
  }
  bool isStaticMode() const noexcept { return static_mode_; }
  bool isValid() const noexcept { return valid_; }
  double getLastDt() const noexcept { return last_dt_; }

  bool consumeUnpackedStateFlag() noexcept {
    const bool f = hadPendingStateUnpacked_;
    hadPendingStateUnpacked_ = false;
    return f;
  }

  void invalidateTangentCache() noexcept {
    tangent_cache_dt_ = -1.0;
    tangent_cache_dt_c_ = -1.0;
  }

  void addMassForce(bool useM0, double scale, const double* a, double* out) const noexcept {
    if (a == nullptr || out == nullptr || !(mat_.rho > 0.0)) return;
    if (lumped_) {
      const auto& ml = useM0 ? M0_lumped_ : M_lumped_;
      for (int aNode = 0; aNode < nNodes_; ++aNode) {
        const double ma = scale * ml[aNode];
        const int uBaseA = uOff_[aNode];
        for (int d = 0; d < ndm_; ++d) {
          out[uBaseA + d] += ma * a[uBaseA + d];
        }
      }
    } else {
      const auto& ms = useM0 ? M0_scalar_ : M_scalar_;
      for (int aNode = 0; aNode < nNodes_; ++aNode) {
        const int uBaseA = uOff_[aNode];
        const size_t rowOff = (size_t)aNode * nNodes_;
        for (int d = 0; d < ndm_; ++d) {
          double s = 0.0;
          for (int bNode = 0; bNode < nNodes_; ++bNode) {
            s += ms[rowOff + bNode] * a[uOff_[bNode] + d];
          }
          out[uBaseA + d] += scale * s;
        }
      }
    }
  }

  void setTrial(const double* u, const double* v, const double* a, double dt) {
    // Static / Initial-tangent regime (dt <= 0.0 or static_mode_ == true):
    // Reduce to unstretched static matrices without perturbing any committed
    // recursive convolution memory variables.
    if (dt <= 0.0 || static_mode_) {
      for (int g = 0; g < nGP_; ++g) {
        gpStates_[g].revertToLastCommit(ndm_, biot_.enabled);
      }
      tangent_cache_dt_ = -1.0;
      K_eff_ = K_0_;
      if (!biot_.enabled) {
        std::fill(C_eff_.begin(), C_eff_.end(), 0.0);
        std::fill(R_vel_.begin(), R_vel_.end(), 0.0);
        double p_ext = 0.0;
        for (int i = 0; i < nDOF_; ++i) {
          double r_i = 0.0;
          for (int j = 0; j < nDOF_; ++j) {
            r_i += K_0_[(size_t)i * nDOF_ + j] * u[j];
          }
          R_disp_[i] = r_i;
          R_int_[i] = r_i;
          p_ext += v[i] * r_i;
        }
        P_ext_ = p_ext;
      } else {
        C_eff_ = C_0_;
        for (int i = 0; i < nDOF_; ++i) {
          double r_disp_i = 0.0;
          double r_vel_i  = 0.0;
          for (int j = 0; j < nDOF_; ++j) {
            const size_t idx = (size_t)i * nDOF_ + j;
            r_disp_i += K_0_[idx] * u[j];
            r_vel_i  += C_0_[idx] * v[j];
          }
          R_disp_[i] = r_disp_i;
          R_vel_[i]  = r_vel_i;
          R_int_[i]  = r_disp_i + r_vel_i;
        }
        addMassForce(true, 1.0, a, R_int_.data());
        double p_ext = 0.0;
        for (int i = 0; i < nDOF_; ++i) {
          p_ext += v[i] * R_int_[i];
        }
        P_ext_ = p_ext;
      }
      last_dt_ = 0.0;
      return;
    }

    last_dt_ = dt;
    const bool rebuildTangents = (dt != tangent_cache_dt_);
    const bool isBox = (geom_.mode == CPMLGeomMode::Box);

    std::fill(R_int_.begin(), R_int_.end(), 0.0);
    std::fill(R_disp_.begin(), R_disp_.end(), 0.0);
    std::fill(R_vel_.begin(), R_vel_.end(), 0.0);
    if (rebuildTangents) {
      std::fill(K_eff_.begin(), K_eff_.end(), 0.0);
      std::fill(C_eff_.begin(), C_eff_.end(), 0.0);
    }

    // Local nodal displacements u'_a = R * u_a and velocities v'_a = R * v_a
    double u_local[shapes::MAX_NODES][3]{};
    double v_local[shapes::MAX_NODES][3]{};
    if (isBox) {
      for (int node = 0; node < nNodes_; ++node) {
        const int uBase = uOff_[node];
        for (int i = 0; i < ndm_; ++i) {
          u_local[node][i] = u[uBase + i];
          v_local[node][i] = v[uBase + i];
        }
      }
    } else {
      for (int node = 0; node < nNodes_; ++node) {
        const int uBase = uOff_[node];
        for (int i = 0; i < ndm_; ++i) {
          double su = 0.0;
          double sv = 0.0;
          for (int p = 0; p < ndm_; ++p) {
            su += geom_.R[i][p] * u[uBase + p];
            sv += geom_.R[i][p] * v[uBase + p];
          }
          u_local[node][i] = su;
          v_local[node][i] = sv;
        }
      }
    }

    const bool useBbar = biot_.enabled && (biot_.formulation == 1);

    // Loop over Gauss points
    for (int g = 0; g < nGP_; ++g) {
      auto& gs = gpStates_[g];
      const double w_vol = gs.w_vol;

      // ---- 1. Stiffness part (in local frame x') ------------------------- //
      double G_local[3][3]{};
      if (!useBbar) {
        for (int k = 0; k < ndm_; ++k) {
          for (int l = 0; l < ndm_; ++l) {
            double s = 0.0;
            for (int node = 0; node < nNodes_; ++node) {
              s += u_local[node][k] * gs.Phi_local[node][l];
            }
            G_local[k][l] = s;
          }
        }
      } else {
        for (int k = 0; k < ndm_; ++k) {
          for (int l = 0; l < ndm_; ++l) {
            double s = 0.0;
            for (int node = 0; node < nNodes_; ++node) {
              for (int m = 0; m < ndm_; ++m) {
                s += gs.B_local[k][l][node][m] * u_local[node][m];
              }
            }
            G_local[k][l] = s;
          }
        }
      }

      // Filter G'_{kl} through F_{jl}(w) for each pair (j, l)
      double G_tilde[3][3][3]{}; // G_tilde[j][l][k] = \widetilde{G}'_{kl}((j, l); t_{n+1})
      double D_F[3][3]{};
      for (int l = 0; l < ndm_; ++l) {
        double col_l[3] = {G_local[0][l], G_local[1][l], G_local[2][l]};
        for (int j = 0; j < ndm_; ++j) {
          gs.F_jl[j][l].evalTrial(col_l, dt, G_tilde[j][l]);
          D_F[j][l] = gs.F_jl[j][l].D_t;
        }
      }

      if (!useBbar) {
        // Hoist sigma_local[i][j] outside nodeA (0-ULP bit-identical to original loop)
        double sigma_local[3][3]{};
        for (int i = 0; i < ndm_; ++i) {
          for (int j = 0; j < ndm_; ++j) {
            double sigma_ij = 0.0;
            for (int k = 0; k < ndm_; ++k) {
              for (int l = 0; l < ndm_; ++l) {
                sigma_ij += mat_.C[i][j][k][l] * G_tilde[j][l][k];
              }
            }
            sigma_local[i][j] = sigma_ij;
          }
        }

        for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
          double fK_local[3] = {0.0, 0.0, 0.0};
          for (int i = 0; i < ndm_; ++i) {
            double s_i = 0.0;
            for (int j = 0; j < ndm_; ++j) {
              s_i += gs.Phi_local[nodeA][j] * sigma_local[i][j];
            }
            fK_local[i] = w_vol * s_i;
          }
          const int uBaseA = uOff_[nodeA];
          if (isBox) {
            for (int i = 0; i < ndm_; ++i) {
              R_disp_[uBaseA + i] += fK_local[i];
            }
          } else {
            for (int p = 0; p < ndm_; ++p) {
              double s_p = 0.0;
              for (int i = 0; i < ndm_; ++i) {
                s_p += geom_.R[i][p] * fK_local[i];
              }
              R_disp_[uBaseA + p] += s_p;
            }
          }
        }

        if (rebuildTangents) {
          for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
            const int uBaseA = uOff_[nodeA];
            for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
              const int uBaseB = uOff_[nodeB];
              double K_local_ab[3][3]{};
              for (int i = 0; i < ndm_; ++i) {
                for (int k = 0; k < ndm_; ++k) {
                  double s_ik = 0.0;
                  for (int j = 0; j < ndm_; ++j) {
                    const double phi_aj = gs.Phi_local[nodeA][j];
                    for (int l = 0; l < ndm_; ++l) {
                      s_ik += phi_aj * mat_.C[i][j][k][l] * D_F[j][l] * gs.Phi_local[nodeB][l];
                    }
                  }
                  K_local_ab[i][k] = w_vol * s_ik;
                }
              }
              if (isBox) {
                for (int i = 0; i < ndm_; ++i) {
                  const int row = uBaseA + i;
                  for (int k = 0; k < ndm_; ++k) {
                    const int col = uBaseB + k;
                    K_eff_[(size_t)row * nDOF_ + col] += K_local_ab[i][k];
                  }
                }
              } else {
                for (int p = 0; p < ndm_; ++p) {
                  for (int q = 0; q < ndm_; ++q) {
                    double s_pq = 0.0;
                    for (int i = 0; i < ndm_; ++i) {
                      const double Rip = geom_.R[i][p];
                      for (int k = 0; k < ndm_; ++k) {
                        s_pq += Rip * K_local_ab[i][k] * geom_.R[k][q];
                      }
                    }
                    const int row = uBaseA + p;
                    const int col = uBaseB + q;
                    K_eff_[(size_t)row * nDOF_ + col] += s_pq;
                  }
                }
              }
            }
          }
        }
      } else {
        // B-bar formulation
        double sigma_local[3][3]{};
        for (int i = 0; i < ndm_; ++i) {
          for (int j = 0; j < ndm_; ++j) {
            double s = 0.0;
            for (int k = 0; k < ndm_; ++k) {
              for (int l = 0; l < ndm_; ++l) {
                s += mat_.C[i][j][k][l] * G_tilde[j][l][k];
              }
            }
            sigma_local[i][j] = s;
          }
        }

        for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
          double fK_local[3] = {0.0, 0.0, 0.0};
          for (int r = 0; r < ndm_; ++r) {
            double s_r = 0.0;
            for (int i = 0; i < ndm_; ++i) {
              for (int j = 0; j < ndm_; ++j) {
                s_r += gs.B_local[i][j][nodeA][r] * sigma_local[i][j];
              }
            }
            fK_local[r] = w_vol * s_r;
          }
          const int uBaseA = uOff_[nodeA];
          if (isBox) {
            for (int r = 0; r < ndm_; ++r) {
              R_disp_[uBaseA + r] += fK_local[r];
            }
          } else {
            for (int p = 0; p < ndm_; ++p) {
              double s_p = 0.0;
              for (int r = 0; r < ndm_; ++r) {
                s_p += geom_.R[r][p] * fK_local[r];
              }
              R_disp_[uBaseA + p] += s_p;
            }
          }
        }

        if (rebuildTangents) {
          double CB_local[3][3][shapes::MAX_NODES][3]{};
          for (int i = 0; i < ndm_; ++i) {
            for (int j = 0; j < ndm_; ++j) {
              for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
                for (int m = 0; m < ndm_; ++m) {
                  double s = 0.0;
                  for (int k = 0; k < ndm_; ++k) {
                    for (int l = 0; l < ndm_; ++l) {
                      s += mat_.C[i][j][k][l] * D_F[j][l] * gs.B_local[k][l][nodeB][m];
                    }
                  }
                  CB_local[i][j][nodeB][m] = s;
                }
              }
            }
          }

          for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
            const int uBaseA = uOff_[nodeA];
            for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
              const int uBaseB = uOff_[nodeB];
              double K_local_ab[3][3]{};
              for (int r = 0; r < ndm_; ++r) {
                for (int m = 0; m < ndm_; ++m) {
                  double s_rm = 0.0;
                  for (int i = 0; i < ndm_; ++i) {
                    for (int j = 0; j < ndm_; ++j) {
                      s_rm += gs.B_local[i][j][nodeA][r] * CB_local[i][j][nodeB][m];
                    }
                  }
                  K_local_ab[r][m] = w_vol * s_rm;
                }
              }
              if (isBox) {
                for (int r = 0; r < ndm_; ++r) {
                  const int row = uBaseA + r;
                  for (int m = 0; m < ndm_; ++m) {
                    const int col = uBaseB + m;
                    K_eff_[(size_t)row * nDOF_ + col] += K_local_ab[r][m];
                  }
                }
              } else {
                for (int p = 0; p < ndm_; ++p) {
                  for (int q = 0; q < ndm_; ++q) {
                    double s_pq = 0.0;
                    for (int r = 0; r < ndm_; ++r) {
                      const double Rrp = geom_.R[r][p];
                      for (int m = 0; m < ndm_; ++m) {
                        s_pq += Rrp * K_local_ab[r][m] * geom_.R[m][q];
                      }
                    }
                    const int row = uBaseA + p;
                    const int col = uBaseB + q;
                    K_eff_[(size_t)row * nDOF_ + col] += s_pq;
                  }
                }
              }
            }
          }
        }
      }

      // ---- 2 & 3. Inertial Velocity & Corner Displacement Parts ---------- //
      if (gs.has_any_d) {
        double v_gp[3] = {0.0, 0.0, 0.0};
        double u_gp[3] = {0.0, 0.0, 0.0};
        for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
          const double Nb = gs.N[nodeB];
          const int uBaseB = uOff_[nodeB];
          for (int d = 0; d < ndm_; ++d) {
            v_gp[d] += Nb * v[uBaseB + d];
            u_gp[d] += Nb * u[uBaseB + d];
          }
        }

        double gv_out[3][3]{};
        for (int m = 0; m < 3; ++m) {
          if (gs.has_d[m]) {
            gs.g_vel[m].evalTrial(v_gp, dt, gv_out[m]);
          }
        }
        const double k0 = gs.kappa[0];
        const double k1 = gs.kappa[1];
        const double k2 = gs.kappa[2];

        double A_vel[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < ndm_; ++d) {
          A_vel[d] = (k1 * k2) * gv_out[0][d]
                   + (k0 * k2) * gv_out[1][d]
                   + (k0 * k1) * gv_out[2][d];
        }
        const double D_vel_g = (gs.has_d[0] ? (k1 * k2) * gs.g_vel[0].D_t : 0.0)
                             + (gs.has_d[1] ? (k0 * k2) * gs.g_vel[1].D_t : 0.0)
                             + (gs.has_d[2] ? (k0 * k1) * gs.g_vel[2].D_t : 0.0);

        double y_g1g0[3]{}, y_g2g0[3]{}, y_g2g1[3]{}, y_phi2g1g0[3]{};
        if (gs.has_d01)  gs.g1_g0.evalTrial(u_gp, dt, y_g1g0);
        if (gs.has_d02)  gs.g2_g0.evalTrial(u_gp, dt, y_g2g0);
        if (gs.has_d12)  gs.g2_g1.evalTrial(u_gp, dt, y_g2g1);
        if (gs.has_d012) gs.phi2_g1_g0.evalTrial(u_gp, dt, y_phi2g1g0);

        double A_disp[3] = {0.0, 0.0, 0.0};
        double D_disp_g = 0.0;
        if (!LADRUNO_CPML_MUT_ACTIVE(CPML_NO_CORNER_CROSS) &&
            (gs.has_d01 || gs.has_d02 || gs.has_d12)) {
          for (int d = 0; d < ndm_; ++d) {
            A_disp[d] = k2 * y_g1g0[d]
                      + k1 * y_g2g0[d]
                      + k0 * y_g2g1[d]
                      + y_phi2g1g0[d];
          }
          D_disp_g = (gs.has_d01  ? k2 * gs.g1_g0.D_t : 0.0)
                   + (gs.has_d02  ? k1 * gs.g2_g0.D_t : 0.0)
                   + (gs.has_d12  ? k0 * gs.g2_g1.D_t : 0.0)
                   + (gs.has_d012 ? gs.phi2_g1_g0.D_t : 0.0);
        }

        const double rho_w = mat_.rho * w_vol;
        for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
          const double Na = gs.N[nodeA];
          const int uBaseA = uOff_[nodeA];
          for (int d = 0; d < ndm_; ++d) {
            const int dof = uBaseA + d;
            R_vel_[dof]  += rho_w * Na * A_vel[d];
            R_disp_[dof] += rho_w * Na * A_disp[d];
          }
          if (rebuildTangents) {
            for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
              const int uBaseB = uOff_[nodeB];
              const double Nab = rho_w * Na * gs.N[nodeB];
              const double c_ab = Nab * D_vel_g;
              const double k_ab = Nab * D_disp_g;
              for (int d = 0; d < ndm_; ++d) {
                const int row = uBaseA + d;
                const int col = uBaseB + d;
                C_eff_[(size_t)row * nDOF_ + col] += c_ab;
                K_eff_[(size_t)row * nDOF_ + col] += k_ab;
              }
            }
          }
        }
      }

      // ---- 3b. Biot u-p Cascades (Darcy F_p, Coupling E_p / E_v, Storage J_s) ---- //
      if (biot_.enabled) {
        // (i) Darcy filter F_{jl}^{(p)}
        double g_p_local[3] = {0.0, 0.0, 0.0};
        for (int l = 0; l < ndm_; ++l) {
          double s = 0.0;
          for (int c = 0; c < nP_; ++c) {
            s += u[pOff_[c]] * gs.PhiP_local[c][l];
          }
          g_p_local[l] = s;
        }
        double g_tilde_p[3][3]{};
        double D_Fp[3][3]{};
        for (int l = 0; l < ndm_; ++l) {
          for (int j = 0; j < ndm_; ++j) {
            g_tilde_p[j][l] = gs.F_p[j][l].evalTrial(g_p_local[l], dt);
            D_Fp[j][l]      = gs.F_p[j][l].D_t;
          }
        }
        for (int b = 0; b < nP_; ++b) {
          const int pRow = pOff_[b];
          double fp_H_b = 0.0;
          for (int j = 0; j < ndm_; ++j) {
            double flux_j = 0.0;
            for (int l = 0; l < ndm_; ++l) {
              flux_j += kbar_local_[j][l] * g_tilde_p[j][l];
            }
            fp_H_b += gs.PhiP_local[b][j] * flux_j;
          }
          R_disp_[pRow] += w_vol * fp_H_b;

          if (rebuildTangents) {
            for (int c = 0; c < nP_; ++c) {
              const int pCol = pOff_[c];
              double h_bc = 0.0;
              for (int j = 0; j < ndm_; ++j) {
                for (int l = 0; l < ndm_; ++l) {
                  h_bc += gs.PhiP_local[b][j] * kbar_local_[j][l] * D_Fp[j][l] * gs.PhiP_local[c][l];
                }
              }
              K_eff_[(size_t)pRow * nDOF_ + pCol] += w_vol * h_bc;
            }
          }
        }

        // (ii) Coupling filters E_j^{(p)} (on u-rows) & E_j^{(\dot{u})} (on p-rows)
        double p_g = 0.0;
        for (int b = 0; b < nP_; ++b) {
          p_g += gs.Np[b] * u[pOff_[b]];
        }
        double p_tilde[3] = {0.0, 0.0, 0.0};
        double D_E[3]     = {0.0, 0.0, 0.0};
        double eps_dot_tilde_sum = 0.0;
        for (int j = 0; j < ndm_; ++j) {
          p_tilde[j] = gs.E_p[j].evalTrial(p_g, dt);
          D_E[j]     = gs.E_p[j].D_t;

          double eps_dot_jj = 0.0;
          for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
            eps_dot_jj += v_local[nodeA][j] * gs.Phi_tilde_local[nodeA][j];
          }
          eps_dot_tilde_sum += gs.E_v[j].evalTrial(eps_dot_jj, dt);
        }

        const double alpha_w = w_vol * biot_.biotAlpha;
        for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
          const int uBaseA = uOff_[nodeA];
          double f_up_local[3] = {0.0, 0.0, 0.0};
          double q_local[3]    = {0.0, 0.0, 0.0};
          for (int j = 0; j < ndm_; ++j) {
            const double phi_aj = gs.Phi_tilde_local[nodeA][j];
            f_up_local[j] = -alpha_w * phi_aj * p_tilde[j];
            if (rebuildTangents) {
              q_local[j]  =  alpha_w * phi_aj * D_E[j];
            }
          }
          for (int d = 0; d < ndm_; ++d) {
            double f_up_global_d = 0.0;
            double q_global_d    = 0.0;
            if (isBox) {
              f_up_global_d = f_up_local[d];
              q_global_d    = q_local[d];
            } else {
              for (int j = 0; j < ndm_; ++j) {
                f_up_global_d += geom_.R[j][d] * f_up_local[j];
                if (rebuildTangents) {
                  q_global_d  += geom_.R[j][d] * q_local[j];
                }
              }
            }
            const int uRow = uBaseA + d;
            R_disp_[uRow] += f_up_global_d;
            if (rebuildTangents) {
              for (int b = 0; b < nP_; ++b) {
                const int pCol = pOff_[b];
                const double q_adb = q_global_d * gs.Np[b];
                K_eff_[(size_t)uRow * nDOF_ + pCol] -= q_adb;
                C_eff_[(size_t)pCol * nDOF_ + uRow] += q_adb;
              }
            }
          }
        }

        for (int b = 0; b < nP_; ++b) {
          R_vel_[pOff_[b]] += alpha_w * gs.Np[b] * eps_dot_tilde_sum;
        }

        // (iii) Storage filter J_s^{(p)}
        double p_dot_g = 0.0;
        for (int c = 0; c < nP_; ++c) {
          p_dot_g += gs.Np[c] * v[pOff_[c]];
        }
        const double p_dot_tilde = gs.J_s.evalTrial(p_dot_g, dt);
        const double D_Js        = gs.J_s.D_t;
        const double invQ_w      = w_vol * biot_.oneOverQbar();
        for (int b = 0; b < nP_; ++b) {
          const int pRow = pOff_[b];
          R_vel_[pRow] += invQ_w * gs.Np[b] * p_dot_tilde;
          if (rebuildTangents) {
            for (int c = 0; c < nP_; ++c) {
              const int pCol = pOff_[c];
              C_eff_[(size_t)pRow * nDOF_ + pCol] += invQ_w * D_Js * gs.Np[b] * gs.Np[c];
            }
          }
        }
      }
    }

    // Add unstretched McGann stabilization H_tilde to R_vel_ and C_eff_ in -up mode
    if (biot_.enabled) {
      for (int b = 0; b < nP_; ++b) {
        const int pRow = pOff_[b];
        double f_stab_b = 0.0;
        for (int c = 0; c < nP_; ++c) {
          const int pCol = pOff_[c];
          const double h_bc = H_tilde_[(size_t)b * nP_ + c];
          f_stab_b += h_bc * v[pCol];
          if (rebuildTangents) {
            C_eff_[(size_t)pRow * nDOF_ + pCol] += h_bc;
          }
        }
        R_vel_[pRow] += f_stab_b;
      }
    }

    if (rebuildTangents) {
      tangent_cache_dt_ = dt;
    }

    // ---- 4. Acceleration Part f_acc = M_eff * a -------------------------- //
    for (int i = 0; i < nDOF_; ++i) {
      R_int_[i] = R_disp_[i] + R_vel_[i];
    }
    addMassForce(false, 1.0, a, R_int_.data());

    // ---- 5. External Power P_ext = v^T R_int ----------------------------- //
    double p_ext = 0.0;
    for (int i = 0; i < nDOF_; ++i) {
      p_ext += v[i] * R_int_[i];
    }
    P_ext_ = p_ext;
  }

  void commit() {
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].commit(ndm_, biot_.enabled);
    }
    R_int_c_ = R_int_;
    R_disp_c_ = R_disp_;
    R_vel_c_ = R_vel_;
    K_eff_c_ = K_eff_;
    C_eff_c_ = C_eff_;
    P_ext_c_ = P_ext_;
    last_dt_c_ = last_dt_;
    tangent_cache_dt_c_ = tangent_cache_dt_;
  }

  void revertToLastCommit() {
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].revertToLastCommit(ndm_, biot_.enabled);
    }
    R_int_ = R_int_c_;
    R_disp_ = R_disp_c_;
    R_vel_ = R_vel_c_;
    K_eff_ = K_eff_c_;
    C_eff_ = C_eff_c_;
    P_ext_ = P_ext_c_;
    last_dt_ = last_dt_c_;
    tangent_cache_dt_ = (!static_mode_ && last_dt_c_ > 0.0) ? tangent_cache_dt_c_ : -1.0;
  }

  void revertToStart() {
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].revertToStart(ndm_, biot_.enabled);
    }
    std::fill(R_int_.begin(), R_int_.end(), 0.0);
    std::fill(R_int_c_.begin(), R_int_c_.end(), 0.0);
    std::fill(R_disp_.begin(), R_disp_.end(), 0.0);
    std::fill(R_disp_c_.begin(), R_disp_c_.end(), 0.0);
    std::fill(R_vel_.begin(), R_vel_.end(), 0.0);
    std::fill(R_vel_c_.begin(), R_vel_c_.end(), 0.0);
    K_eff_ = K_0_;
    K_eff_c_ = K_0_;
    C_eff_ = C_0_;
    C_eff_c_ = C_0_;
    P_ext_ = 0.0;
    P_ext_c_ = 0.0;
    last_dt_ = 0.0;
    last_dt_c_ = 0.0;
    invalidateTangentCache();
  }

  // ---- Accessors --------------------------------------------------------- //
  int ndm() const { return ndm_; }
  int numNodes() const { return nNodes_; }
  int numGP() const { return nGP_; }
  int numDOF() const { return nDOF_; }
  bool isLumped() const { return lumped_; }
  CPMLTopology topology() const { return topo_; }

  bool isUPMode() const noexcept { return biot_.enabled; }
  int numPressureDOFs() const noexcept { return nP_; }
  int uOffset(int node) const noexcept { return uOff_[node]; }
  int pOffset(int node) const noexcept { return pOff_[node]; }
  const int* uOffsets() const noexcept { return uOff_.data(); }
  const int* pOffsets() const noexcept { return pOff_.data(); }
  const CPMLBiotParams& getBiotParams() const noexcept { return biot_; }
  const CPMLMaterialParams& getMaterialParams() const noexcept { return mat_; }
  const CPMLGradingParams& getGradingParams() const noexcept { return grading_; }
  const CPMLGeomConfig& getGeomConfig() const noexcept { return geom_; }

  double computeVolume() const {
    double vol = 0.0;
    for (int g = 0; g < nGP_; ++g) {
      vol += gpStates_[g].w_vol;
    }
    return vol;
  }

  void computeNodalVolumes(double* vNode) const {
    for (int a = 0; a < nNodes_; ++a) {
      vNode[a] = 0.0;
    }
    for (int g = 0; g < nGP_; ++g) {
      const double w_vol = gpStates_[g].w_vol;
      for (int a = 0; a < nNodes_; ++a) {
        vNode[a] += w_vol * gpStates_[g].N[a];
      }
    }
  }

  void computePressureGradients(int gp, double* dNp_dx) const {
    if (dNp_dx == nullptr || gp < 0 || gp >= nGP_) return;
    const auto& gs = gpStates_[gp];
    for (int b = 0; b < nP_; ++b) {
      for (int m = 0; m < ndm_; ++m) {
        dNp_dx[b * ndm_ + m] = gs.dNpdx[b][m];
      }
    }
  }

  void computeSeepageLoad(const double* fluidBody, double* fseep_out) const {
    if (fseep_out == nullptr) return;
    for (int b = 0; b < nP_; ++b) {
      fseep_out[b] = 0.0;
    }
    if (!biot_.enabled || fluidBody == nullptr) return;
    for (int g = 0; g < nGP_; ++g) {
      const auto& gs = gpStates_[g];
      const double w_vol = gs.w_vol;
      for (int b = 0; b < nP_; ++b) {
        double s = 0.0;
        for (int i = 0; i < ndm_; ++i) {
          s += gs.dNpdx[b][i] * biot_.kbar[i] * biot_.rhoF * fluidBody[i];
        }
        fseep_out[b] += w_vol * s;
      }
    }
  }

  const std::vector<double>& getR_int() const { return R_int_; }
  const std::vector<double>& getR_disp() const { return R_disp_; }
  const std::vector<double>& getR_vel() const { return R_vel_; }
  const std::vector<double>& getK_eff() const { return K_eff_; }
  const std::vector<double>& getC_eff() const { return C_eff_; }
  const std::vector<double>& getM_eff() const { return M_eff_; }
  const std::vector<double>& getK0() const { return K_0_; }
  const std::vector<double>& getC0() const { return C_0_; }
  const std::vector<double>& getM0() const { return M_0_; }
  const std::vector<double>& getHtilde() const { return H_tilde_; }
  double getP_ext() const { return P_ext_; }

  void getR_int(double* out) const {
    std::copy(R_int_.begin(), R_int_.end(), out);
  }
  void getR_disp(double* out) const {
    std::copy(R_disp_.begin(), R_disp_.end(), out);
  }
  void getR_vel(double* out) const {
    std::copy(R_vel_.begin(), R_vel_.end(), out);
  }
  void getK_eff(double* out) const {
    std::copy(K_eff_.begin(), K_eff_.end(), out);
  }
  void getC_eff(double* out) const {
    std::copy(C_eff_.begin(), C_eff_.end(), out);
  }
  void getM_eff(double* out) const {
    std::copy(M_eff_.begin(), M_eff_.end(), out);
  }
  void getK0(double* out) const {
    std::copy(K_0_.begin(), K_0_.end(), out);
  }
  void getC0(double* out) const {
    std::copy(C_0_.begin(), C_0_.end(), out);
  }
  void getM0(double* out) const {
    std::copy(M_0_.begin(), M_0_.end(), out);
  }

  void getK_T_dyn(double c1, double c2, double c3, double* out) const {
    const size_t n2 = (size_t)nDOF_ * nDOF_;
    for (size_t i = 0; i < n2; ++i) {
      out[i] = c1 * K_eff_[i] + c2 * C_eff_[i] + c3 * M_eff_[i];
    }
  }

  std::vector<double> getK_T_dyn(double c1, double c2, double c3) const {
    std::vector<double> out((size_t)nDOF_ * nDOF_, 0.0);
    getK_T_dyn(c1, c2, c3, out.data());
    return out;
  }

  // ---- State Serialization (for sendSelf / recvSelf) --------------------- //
  int numStateDoubles() const {
    int total = 0;
    for (int g = 0; g < nGP_; ++g) {
      total += gpStates_[g].numStateDoubles(ndm_, biot_.enabled);
    }
    const int n2 = nDOF_ * nDOF_;
    // Also pack trial & committed R_int, R_disp, R_vel, K_eff, C_eff, P_ext, last_dt
    total += 6 * nDOF_ + 4 * n2 + 4;
    return total;
  }

  void packState(double* buf) const {
    double* ptr = buf;
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].packState(ndm_, ptr, biot_.enabled);
    }
    for (double v : R_int_)   *ptr++ = v;
    for (double v : R_int_c_) *ptr++ = v;
    for (double v : R_disp_)   *ptr++ = v;
    for (double v : R_disp_c_) *ptr++ = v;
    for (double v : R_vel_)   *ptr++ = v;
    for (double v : R_vel_c_) *ptr++ = v;
    for (double v : K_eff_)   *ptr++ = v;
    for (double v : K_eff_c_) *ptr++ = v;
    for (double v : C_eff_)   *ptr++ = v;
    for (double v : C_eff_c_) *ptr++ = v;
    *ptr++ = P_ext_;
    *ptr++ = P_ext_c_;
    *ptr++ = last_dt_;
    *ptr++ = last_dt_c_;
  }

  void unpackState(const double* buf) {
    const double* ptr = buf;
    for (int g = 0; g < nGP_; ++g) {
      gpStates_[g].unpackState(ndm_, ptr, biot_.enabled);
    }
    for (double& v : R_int_)   v = *ptr++;
    for (double& v : R_int_c_) v = *ptr++;
    for (double& v : R_disp_)   v = *ptr++;
    for (double& v : R_disp_c_) v = *ptr++;
    for (double& v : R_vel_)   v = *ptr++;
    for (double& v : R_vel_c_) v = *ptr++;
    for (double& v : K_eff_)   v = *ptr++;
    for (double& v : K_eff_c_) v = *ptr++;
    for (double& v : C_eff_)   v = *ptr++;
    for (double& v : C_eff_c_) v = *ptr++;
    P_ext_     = *ptr++;
    P_ext_c_   = *ptr++;
    last_dt_   = *ptr++;
    last_dt_c_ = *ptr++;
    invalidateTangentCache();
  }

  void setPendingState(const double* buf, int n) {
    if (buf != nullptr && n > 0) {
      pendingState_.assign(buf, buf + n);
    } else {
      pendingState_.clear();
    }
  }

private:
  void rebuildStaticAndMassMatrices() {
    invalidateTangentCache();
    K_0_.assign((size_t)nDOF_ * nDOF_, 0.0);
    C_0_.assign((size_t)nDOF_ * nDOF_, 0.0);
    M_0_.assign((size_t)nDOF_ * nDOF_, 0.0);
    M_eff_.assign((size_t)nDOF_ * nDOF_, 0.0);
    H_tilde_.assign((size_t)nP_ * nP_, 0.0);

    for (int j = 0; j < 3; ++j) {
      for (int l = 0; l < 3; ++l) {
        double s = 0.0;
        if (biot_.enabled && j < ndm_ && l < ndm_) {
          for (int m = 0; m < ndm_; ++m) {
            s += geom_.R[j][m] * biot_.kbar[m] * geom_.R[l][m];
          }
        }
        kbar_local_[j][l] = s;
      }
    }

    const double h_elem = biot_.enabled ? shapes::elemSizeH(topo_, coords_.data()) : 0.0;
    const double alpha_stab = biot_.stabAlpha(h_elem, mat_.lambda + 2.0 * mat_.mu);
    const double invQ = biot_.oneOverQbar();
    const bool useBbar = biot_.enabled && (biot_.formulation == 1);

    M0_scalar_.assign((size_t)nNodes_ * nNodes_, 0.0);
    M_scalar_.assign((size_t)nNodes_ * nNodes_, 0.0);
    M0_lumped_.assign(nNodes_, 0.0);
    M_lumped_.assign(nNodes_, 0.0);

    for (int g = 0; g < nGP_; ++g) {
      const auto& gs = gpStates_[g];
      const double w_vol = gs.w_vol;
      const double D_acc = gs.kappa[0] * gs.kappa[1] * gs.kappa[2];
      const double coeff0 = mat_.rho * w_vol;
      const double coeff  = coeff0 * D_acc;

      for (int a = 0; a < nNodes_; ++a) {
        for (int b = 0; b < nNodes_; ++b) {
          M0_scalar_[(size_t)a * nNodes_ + b] += coeff0 * gs.N[a] * gs.N[b];
          M_scalar_[(size_t)a * nNodes_ + b]  += coeff  * gs.N[a] * gs.N[b];
        }
      }

      if (!useBbar) {
        for (int a = 0; a < nNodes_; ++a) {
          const int uBaseA = uOff_[a];
          for (int b = 0; b < nNodes_; ++b) {
            const int uBaseB = uOff_[b];
            for (int i = 0; i < ndm_; ++i) {
              const int row = uBaseA + i;
              for (int k = 0; k < ndm_; ++k) {
                double s_ik = 0.0;
                for (int j = 0; j < ndm_; ++j) {
                  for (int l = 0; l < ndm_; ++l) {
                    s_ik += gs.dNdx[a][j] * mat_.C[i][j][k][l] * gs.dNdx[b][l];
                  }
                }
                const int col = uBaseB + k;
                K_0_[(size_t)row * nDOF_ + col] += w_vol * s_ik;
              }
            }
          }
        }
      } else {
        double CB0_local[3][3][shapes::MAX_NODES][3]{};
        for (int i = 0; i < ndm_; ++i) {
          for (int j = 0; j < ndm_; ++j) {
            for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
              for (int m = 0; m < ndm_; ++m) {
                double s = 0.0;
                for (int k = 0; k < ndm_; ++k) {
                  for (int l = 0; l < ndm_; ++l) {
                    s += mat_.C[i][j][k][l] * gs.B_local[k][l][nodeB][m];
                  }
                }
                CB0_local[i][j][nodeB][m] = s;
              }
            }
          }
        }
        for (int nodeA = 0; nodeA < nNodes_; ++nodeA) {
          const int uBaseA = uOff_[nodeA];
          for (int nodeB = 0; nodeB < nNodes_; ++nodeB) {
            const int uBaseB = uOff_[nodeB];
            double K_local_ab[3][3]{};
            for (int r = 0; r < ndm_; ++r) {
              for (int m = 0; m < ndm_; ++m) {
                double s_rm = 0.0;
                for (int i = 0; i < ndm_; ++i) {
                  for (int j = 0; j < ndm_; ++j) {
                    s_rm += gs.B_local[i][j][nodeA][r] * CB0_local[i][j][nodeB][m];
                  }
                }
                K_local_ab[r][m] = w_vol * s_rm;
              }
            }
            for (int p = 0; p < ndm_; ++p) {
              for (int q = 0; q < ndm_; ++q) {
                double s_pq = 0.0;
                for (int r = 0; r < ndm_; ++r) {
                  const double Rrp = geom_.R[r][p];
                  for (int m = 0; m < ndm_; ++m) {
                    s_pq += Rrp * K_local_ab[r][m] * geom_.R[m][q];
                  }
                }
                const int row = uBaseA + p;
                const int col = uBaseB + q;
                K_0_[(size_t)row * nDOF_ + col] += s_pq;
              }
            }
          }
        }
      }

      if (biot_.enabled) {
        const double alpha_w = w_vol * biot_.biotAlpha;
        for (int a = 0; a < nNodes_; ++a) {
          const int uBaseA = uOff_[a];
          for (int d = 0; d < ndm_; ++d) {
            double phi_tilde_global_d = 0.0;
            for (int j = 0; j < ndm_; ++j) {
              phi_tilde_global_d += geom_.R[j][d] * gs.Phi_tilde_local[a][j];
            }
            const int uRow = uBaseA + d;
            for (int b = 0; b < nP_; ++b) {
              const int pCol = pOff_[b];
              const double q0_adb = alpha_w * phi_tilde_global_d * gs.Np[b];
              K_0_[(size_t)uRow * nDOF_ + pCol] -= q0_adb;
              C_0_[(size_t)pCol * nDOF_ + uRow] += q0_adb;
            }
          }
        }

        const double invQ_w = w_vol * invQ;
        const double stab_w = w_vol * alpha_stab;
        for (int b = 0; b < nP_; ++b) {
          const int pRow = pOff_[b];
          for (int c = 0; c < nP_; ++c) {
            const int pCol = pOff_[c];
            double h0_bc = 0.0;
            for (int j = 0; j < ndm_; ++j) {
              for (int l = 0; l < ndm_; ++l) {
                h0_bc += gs.PhiP_local[b][j] * kbar_local_[j][l] * gs.PhiP_local[c][l];
              }
            }
            K_0_[(size_t)pRow * nDOF_ + pCol] += w_vol * h0_bc;
            C_0_[(size_t)pRow * nDOF_ + pCol] += invQ_w * gs.Np[b] * gs.Np[c];

            double grad_dot = 0.0;
            for (int m = 0; m < ndm_; ++m) {
              grad_dot += gs.dNpdx[b][m] * gs.dNpdx[c][m];
            }
            H_tilde_[(size_t)b * nP_ + c] += stab_w * grad_dot;
          }
        }
      }
    }

    if (biot_.enabled) {
      for (int b = 0; b < nP_; ++b) {
        const int pRow = pOff_[b];
        for (int c = 0; c < nP_; ++c) {
          const int pCol = pOff_[c];
          C_0_[(size_t)pRow * nDOF_ + pCol] += H_tilde_[(size_t)b * nP_ + c];
        }
      }
    }

    auto fillMassFromScalar = [&](const std::vector<double>& Ms,
                                  std::vector<double>& MlumpedOut,
                                  std::vector<double>& Mout) {
      if (!lumped_) {
        for (int a = 0; a < nNodes_; ++a) {
          const int uBaseA = uOff_[a];
          for (int b = 0; b < nNodes_; ++b) {
            const int uBaseB = uOff_[b];
            const double mab = Ms[(size_t)a * nNodes_ + b];
            for (int d = 0; d < ndm_; ++d) {
              const int row = uBaseA + d;
              const int col = uBaseB + d;
              Mout[(size_t)row * nDOF_ + col] = mab;
            }
          }
        }
      } else {
        double M_tot = 0.0;
        double D_sum = 0.0;
        for (int a = 0; a < nNodes_; ++a) {
          D_sum += Ms[(size_t)a * nNodes_ + a];
          for (int b = 0; b < nNodes_; ++b) {
            M_tot += Ms[(size_t)a * nNodes_ + b];
          }
        }
        const double scale = (D_sum > 0.0) ? (M_tot / D_sum) : 1.0;
        for (int a = 0; a < nNodes_; ++a) {
          const int uBaseA = uOff_[a];
          const double ma_hrz = scale * Ms[(size_t)a * nNodes_ + a];
          MlumpedOut[a] = ma_hrz;
          for (int d = 0; d < ndm_; ++d) {
            const int row = uBaseA + d;
            Mout[(size_t)row * nDOF_ + row] = ma_hrz;
          }
        }
      }
    };

    fillMassFromScalar(M0_scalar_, M0_lumped_, M_0_);
    fillMassFromScalar(M_scalar_,  M_lumped_,  M_eff_);
  }

  CPMLTopology topo_ = CPMLTopology::Q4;
  bool lumped_ = false;
  bool static_mode_ = false;
  bool valid_ = false;
  int ndm_ = 2;
  int nNodes_ = 4;
  int nGP_ = 4;
  int nP_ = 0;
  int nDOF_ = 8;

  CPMLMaterialParams mat_{};
  CPMLGradingParams grading_{};
  CPMLGeomConfig geom_{};
  CPMLBiotParams biot_{};

  std::vector<int> uOff_{};
  std::vector<int> pOff_{};
  std::vector<int> pCarrier_{};
  double Phi_bar_local_[shapes::MAX_NODES][3]{};
  double kbar_local_[3][3]{};
  std::vector<double> H_tilde_{};
  std::vector<double> M0_scalar_{}, M_scalar_{};
  std::vector<double> M0_lumped_{}, M_lumped_{};

  std::vector<double> coords_{};
  std::vector<CPMLGaussState> gpStates_{};

  std::vector<double> K_0_{};
  std::vector<double> C_0_{};
  std::vector<double> M_0_{};
  std::vector<double> M_eff_{};
  std::vector<double> R_int_{}, R_int_c_{};
  std::vector<double> R_disp_{}, R_disp_c_{};
  std::vector<double> R_vel_{}, R_vel_c_{};
  std::vector<double> K_eff_{}, K_eff_c_{};
  std::vector<double> C_eff_{}, C_eff_c_{};
  double P_ext_ = 0.0, P_ext_c_ = 0.0;
  double last_dt_ = 0.0, last_dt_c_ = 0.0;
  double tangent_cache_dt_ = -1.0, tangent_cache_dt_c_ = -1.0;
  std::vector<double> pendingState_{};
  bool hadPendingStateUnpacked_ = false;
};

} // namespace LADRUNO_CPML_KERNEL_NS

namespace ladruno_cpml = LADRUNO_CPML_KERNEL_NS;

#endif // LadrunoCPMLKernel_h
