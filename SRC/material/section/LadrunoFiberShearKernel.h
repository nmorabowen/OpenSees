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

// Ladruno WP-245 (program #1131, rcbf-3): the shear map of LadrunoNDFiberSection3d.
//
// Header-only and dependency-free, so tests/_testbed/rc_section_gpp.cpp checks
// it with plain g++ (closed forms, FD, symmetry) without an OpenSees build.
//
// A fibre at (y, z), measured from the section centroid, takes the section
// deformation e = (eps0, kz, ky, gy, gz, tw) -- the stock NDFiberSection3d
// order P, MZ, MY, VY, VZ, T -- through a strain map a (3 x 6):
//
//     eps_xx = eps0 - y kz + z ky
//     gam_xy = phiY(y) gy - z tw
//     gam_xz = phiZ(z) gz + y tw
//
// and returns its stress through a force map b (3 x 6):
//
//     s += A b^T sigma ,    K += A b^T D a
//
// Profiles (the shear-strain shape across the section):
//   uniform   : phi = sqrt(alpha), b = a. This is the stock NDFiberSection3d
//               rule, kept so a uniform Ladruno section reproduces stock.
//   parabolic : phiY(y) = 1.5 (1 - (2y/hy)^2), phiZ(z) = 1.5 (1 - (2z/hz)^2),
//               Petrangeli, Pinto & Ciampi (1999), as Model C uses it: zero at
//               the faces, 1.5 at mid-depth, evaluated at the fibre centroid with
//               no renormalisation (the discrete mean of 40 equal strips is
//               1 + 1/3200).
//
// Resultants (parabolic only; the force map's shear rows):
//   modelC    : b's shear rows carry 1, not phi. V = sum tau A, as Model C
//               computes it (06 §1.4). Not work-conjugate: an elastic rectangle
//               gives K_VV = G A sum(phi A)/A (= G A for the continuum), and the
//               tangent is NOT symmetric once the fibre couples eps and gamma.
//               It needs an unsymmetric solver (UmfPack, SparseGEN, Mumps).
//   conjugate : b = a. Work-conjugate and symmetric whenever D is. An elastic
//               rectangle gives K_VV = G A mean(phi^2) = 6/5 G A in the
//               continuum: a kinematically imposed parabola STIFFENS the section.
//               (Timoshenko's 5/6 needs a parabolic STRESS, which a kinematic
//               fibre section cannot impose.)
//
// Torsion keeps the stock Saint-Venant rows (-z tw, +y tw) in both maps.

#ifndef LadrunoFiberShearKernel_h
#define LadrunoFiberShearKernel_h

#include <math.h>

namespace ladruno_shear {

enum Profile   { UNIFORM = 0, PARABOLIC = 1 };
enum Resultant { MODELC = 0, CONJUGATE = 1 };

struct Spec {
  int    profile;    // Profile
  int    resultant;  // Resultant (parabolic only)
  double hy, hz;     // section extent along y and z (parabolic only)
  double alpha;      // stock shear shape factor (uniform only)
  Spec() : profile(UNIFORM), resultant(MODELC), hy(0.0), hz(0.0), alpha(1.0) {}
};

inline double parabola(double c, double h)
{
  const double u = 2.0 * c / h;
  return 1.5 * (1.0 - u * u);
}

// Is the fibre centroid inside the profile band? (parabolic only; a fibre past
// the face would get a negative phi, i.e. shear strain against the section's.)
inline bool inside(const Spec& S, double y, double z)
{
  if (S.profile != PARABOLIC) return true;
  const double ty = 0.5 * S.hy * (1.0 + 1e-9), tz = 0.5 * S.hz * (1.0 + 1e-9);
  return fabs(y) <= ty && fabs(z) <= tz;
}

// Strain map a and force map b for a fibre at (y, z).
inline void maps(const Spec& S, double y, double z, double a[3][6], double b[3][6])
{
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 6; ++c) a[r][c] = 0.0;
  double py, pz;
  if (S.profile == PARABOLIC) { py = parabola(y, S.hy); pz = parabola(z, S.hz); }
  else                        { py = pz = (S.alpha != 1.0) ? sqrt(S.alpha) : 1.0; }
  a[0][0] = 1.0; a[0][1] = -y; a[0][2] = z;
  a[1][3] = py;  a[1][5] = -z;
  a[2][4] = pz;  a[2][5] = y;
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 6; ++c) b[r][c] = a[r][c];
  if (S.profile == PARABOLIC && S.resultant == MODELC) { b[1][3] = 1.0; b[2][4] = 1.0; }
}

// eps = a e
inline void strain(const double a[3][6], const double e[6], double eps[3])
{
  for (int r = 0; r < 3; ++r) {
    double v = 0.0;
    for (int c = 0; c < 6; ++c) v += a[r][c] * e[c];
    eps[r] = v;
  }
}

// s += A b^T sig ; K += A b^T D a   (s or K may be null)
inline void accumulate(double A, const double a[3][6], const double b[3][6],
                       const double sig[3], const double D[3][3],
                       double* s, double (*K)[6])
{
  if (s != 0)
    for (int c = 0; c < 6; ++c) {
      double v = 0.0;
      for (int r = 0; r < 3; ++r) v += b[r][c] * sig[r];
      s[c] += A * v;
    }
  if (K != 0) {
    double Da[3][6];
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 6; ++c) {
        double v = 0.0;
        for (int k = 0; k < 3; ++k) v += D[r][k] * a[k][c];
        Da[r][c] = A * v;
      }
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) {
        double v = 0.0;
        for (int r = 0; r < 3; ++r) v += b[r][i] * Da[r][j];
        K[i][j] += v;
      }
  }
}

}  // namespace ladruno_shear

#endif
