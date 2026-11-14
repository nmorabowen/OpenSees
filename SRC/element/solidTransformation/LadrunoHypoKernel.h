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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoHypoKernel — the PURE numerical core of the `-geom hypo` rate-form
// updated-Lagrangian geometry method (ADR 79). Plain doubles + <math.h>, NO
// OpenSees dependency, so it is unit-tested standalone against the numpy
// oracle (tests/hypo_reference.py) WITHOUT building OpenSees (the
// LadrunoUPKernel.h / ADR-70 idiom), and reused by every consuming element
// (LadrunoBrick P1, LadrunoUP P2) — each supplies only its reference shape
// gradients and nodal displacements; the objective-increment mechanics live
// here, written once.
//
// THE ALGORITHM (ADR 79 §3 — Hughes–Winget midpoint increment, integrated in
// the UNROTATED / Green–Naghdi material frame):
//
//   F_n   = I + Σ_a u_n,a   ⊗ ∂N_a/∂X          (committed config)
//   F_1   = I + Σ_a u_1,a   ⊗ ∂N_a/∂X          (trial config)
//   F_mid = ½ (F_n + F_1)                       (HW midpoint configuration)
//   ∂N_a/∂x̄ = (∂N_a/∂X) · F_mid⁻¹              (midpoint spatial gradients)
//   G     = Σ_a (u_1,a − u_n,a) ⊗ ∂N_a/∂x̄      (incremental disp gradient)
//   Δε_sp = sym(G)                              (objective strain increment —
//                                                EXACTLY zero for a rigid-
//                                                rotation increment of any
//                                                magnitude: G is then a pure
//                                                skew 2·tan(Δθ/2)·ŵ)
//   R_mid = polar(F_mid)                        (Higham iteration)
//   Δε_mat = R_midᵀ · Δε_sp · R_mid            (de-rotate into the FIXED
//                                                unrotated material frame —
//                                                the material's state never
//                                                needs rotating; ADR 79 §2)
//
// The consuming element accumulates ε_feed = ε_feed_n + Δε_mat (Voigt,
// ENGINEERING shear, order {00,11,22,01,12,20}) and drives the UNCHANGED
// small-strain material via setTrialStrain(ε_feed). At assembly it pushes the
// material (unrotated-frame) stress and 6×6 tangent forward with
// R_1 = polar(F_1) via pushStress6 / pushTangent6 below.
//
// LAYOUTS (all row-major, caller-owned):
//   3×3 tensors: A[3*i+j].
//   dNdX[a*3 + J] = ∂N_a/∂X_J   reference cartesian shape gradients at the GP
//   u[a*3 + i]                  nodal displacement component i of node a
//   Voigt 6-vectors: {00, 11, 22, 01, 12, 20}; STRAIN carries engineering
//   shear (γ = 2ε_offdiag); STRESS carries the plain tensor components.
//   6×6 matrices: D[6*r + c].
//
// The kernel neither allocates nor prints: failures return false and the
// caller (element) owns the diagnostic. Geometry guards (det > 0) ARE
// enforced here — a non-positive det(F_mid) or det(F_1) is a collapsed /
// inverted configuration and must fail loudly at the call site.

#ifndef LadrunoHypoKernel_h
#define LadrunoHypoKernel_h

#include <math.h>

namespace ladruno_hypo {

// --------------------------------------------------------------------------- //
//  3×3 helpers (row-major)                                                    //
// --------------------------------------------------------------------------- //

inline double det3(const double A[9])
{
  return A[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(A[3]*A[8]-A[5]*A[6])
       + A[2]*(A[3]*A[7]-A[4]*A[6]);
}

// Ai = A⁻¹ via adjugate/det. Returns det(A); Ai zero-filled if det underflows.
inline double inv3(const double A[9], double Ai[9])
{
  double det = det3(A);
  double id = (det != 0.0) ? 1.0 / det : 0.0;
  Ai[0] =  (A[4]*A[8]-A[5]*A[7])*id;  Ai[1] = -(A[1]*A[8]-A[2]*A[7])*id;
  Ai[2] =  (A[1]*A[5]-A[2]*A[4])*id;  Ai[3] = -(A[3]*A[8]-A[5]*A[6])*id;
  Ai[4] =  (A[0]*A[8]-A[2]*A[6])*id;  Ai[5] = -(A[0]*A[5]-A[2]*A[3])*id;
  Ai[6] =  (A[3]*A[7]-A[4]*A[6])*id;  Ai[7] = -(A[0]*A[7]-A[1]*A[6])*id;
  Ai[8] =  (A[0]*A[4]-A[1]*A[3])*id;
  return det;
}

inline double frob3(const double A[9])
{
  double s = 0.0; for (int i = 0; i < 9; i++) s += A[i]*A[i]; return sqrt(s);
}

// Displacement gradient H = Σ_a u_a ⊗ ∂N_a/∂X  (H[3i+J] = Σ u_a[i]·dNdX[a][J]).
inline void dispGradient(int nNodes, const double *dNdX, const double *u,
                         double H[9])
{
  for (int i = 0; i < 9; i++) H[i] = 0.0;
  for (int a = 0; a < nNodes; a++)
    for (int i = 0; i < 3; i++)
      for (int J = 0; J < 3; J++)
        H[3*i + J] += u[a*3 + i] * dNdX[a*3 + J];
}

// Polar rotation R = polar(F) by the scaled Higham Newton iteration
//   X_{k+1} = ½ ( γ X + γ⁻¹ X⁻ᵀ ),  γ = (‖X⁻¹‖_F / ‖X‖_F)^½
// (the SolidTransformationCorot algorithm — one source of truth for the
// numerics; duplicated here because the kernel must stay OpenSees-free).
// Returns false if det(F) ≤ 0 (collapsed / inverted).
inline bool polar3(const double F[9], double R[9])
{
  double X[9]; for (int i = 0; i < 9; i++) X[i] = F[i];
  if (det3(X) <= 0.0) return false;

  const int    maxIter = 50;
  const double tol     = 1.0e-15;
  double Xinv[9], Xnew[9];

  for (int it = 0; it < maxIter; it++) {
    if (inv3(X, Xinv) == 0.0) return false;
    double g = sqrt(frob3(Xinv) / frob3(X));
    if (!(g > 0.0)) g = 1.0;
    double diff = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double xn = 0.5 * (g * X[3*i+j] + (1.0/g) * Xinv[3*j+i]);  // X⁻ᵀ
        Xnew[3*i+j] = xn;
      }
    for (int i = 0; i < 9; i++) {
      double d = Xnew[i] - X[i];
      diff += d*d;
      X[i] = Xnew[i];
    }
    if (sqrt(diff) < tol) break;
  }
  for (int i = 0; i < 9; i++) R[i] = X[i];
  return det3(R) > 0.0;
}

// --------------------------------------------------------------------------- //
//  The objective increment (the P0 heart)                                     //
// --------------------------------------------------------------------------- //

// Hughes–Winget midpoint objective strain increment, de-rotated into the
// unrotated (Green–Naghdi) material frame — ADR 79 §3 steps 1–4.
//
//   in : nNodes, dNdX (reference gradients at this GP), un (committed nodal
//        displacements), u1 (trial nodal displacements)
//   out: depsMat6 — Δε_mat, Voigt {00,11,22,01,12,20}, ENGINEERING shear
//        R1[9]    — polar(F_1), the push-forward rotation for assembly
//        J1       — det F_1 (> 0)
//
// Returns false when det(F_mid) ≤ 0, det(F_1) ≤ 0, or a polar fails — a
// collapsed/inverted trial or midpoint configuration the element must reject.
inline bool hypoIncrement(int nNodes, const double *dNdX,
                          const double *un, const double *u1,
                          double depsMat6[6], double R1[9], double &J1)
{
  // F_n, F_1, F_mid (as displacement gradients + I)
  double Hn[9], H1[9];
  dispGradient(nNodes, dNdX, un, Hn);
  dispGradient(nNodes, dNdX, u1, H1);

  double F1[9], Fmid[9];
  for (int i = 0; i < 9; i++) {
    F1[i]   = H1[i];
    Fmid[i] = 0.5 * (Hn[i] + H1[i]);
  }
  F1[0] += 1.0; F1[4] += 1.0; F1[8] += 1.0;
  Fmid[0] += 1.0; Fmid[4] += 1.0; Fmid[8] += 1.0;

  double FmidInv[9];
  if (inv3(Fmid, FmidInv) <= 0.0) return false;   // midpoint collapsed/inverted

  // Incremental displacement gradient on the midpoint configuration:
  //   G_ij = Σ_a Δu_a[i] · ∂N_a/∂x̄_j ,  ∂N_a/∂x̄_j = Σ_K dNdX[a][K]·FmidInv[K][j]
  double G[9] = {0,0,0, 0,0,0, 0,0,0};
  for (int a = 0; a < nNodes; a++) {
    double gax[3];
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int K = 0; K < 3; K++)
        s += dNdX[a*3 + K] * FmidInv[3*K + j];
      gax[j] = s;
    }
    for (int i = 0; i < 3; i++) {
      double du = u1[a*3 + i] - un[a*3 + i];
      for (int j = 0; j < 3; j++)
        G[3*i + j] += du * gax[j];
    }
  }

  // Objective spatial increment Δε_sp = sym(G) (tensor components).
  double e[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      e[3*i+j] = 0.5 * (G[3*i+j] + G[3*j+i]);

  // De-rotate into the unrotated material frame: Δε_mat = R_midᵀ·Δε_sp·R_mid.
  double Rmid[9];
  if (!polar3(Fmid, Rmid)) return false;
  double tmp[9];                                   // tmp = Δε_sp · R_mid
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += e[3*i+k] * Rmid[3*k+j];
      tmp[3*i+j] = s;
    }
  double em[9];                                    // em = R_midᵀ · tmp
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += Rmid[3*k+i] * tmp[3*k+j];
      em[3*i+j] = s;
    }

  // Voigt {00,11,22,01,12,20}, engineering shear.
  depsMat6[0] = em[0];
  depsMat6[1] = em[4];
  depsMat6[2] = em[8];
  depsMat6[3] = em[1] + em[3];    // 2·ε01
  depsMat6[4] = em[5] + em[7];    // 2·ε12
  depsMat6[5] = em[2] + em[6];    // 2·ε20

  // Trial-config rotation + Jacobian for the assembly push-forward.
  J1 = det3(F1);
  if (J1 <= 0.0) return false;
  return polar3(F1, R1);
}

// --------------------------------------------------------------------------- //
//  The incremental entry on (F_n, F_1) — WP-232 (MPM), interface §8           //
// --------------------------------------------------------------------------- //
//
// hypoIncrement() above builds F from reference gradients and committed / trial
// NODAL displacements. A material point carried by the MPM engine has neither:
// the background grid is reset every step, so the particle carries F_n itself
// and the step supplies F_1 = ΔF·F_n. This entry takes the two deformation
// gradients directly and runs the same Hughes–Winget algebra
// (tests/mpm_interface.md §8.2, frozen by WP-232b):
//
//   F_mid = ½ (F_n + F_1)                      refuse if det F_mid ≤ 0
//   G     = (F_1 − F_n) · F_mid⁻¹
//   Δε    = sym(G)
//   R_mid = polar(F_mid)
//   deps  = Voigt(R_midᵀ · Δε · R_mid), engineering shear
//   R1    = polar(F_1),  J1 = det F_1          refuse if J1 ≤ 0
//
// With F_1 = I + Σ u_1 ⊗ ∇N and F_n = I + Σ u_n ⊗ ∇N this is the existing entry
// (G = Σ Δu ⊗ ∇N · F_mid⁻¹ = (F_1 − F_n)·F_mid⁻¹); it is written separately, not
// as a refactor, so the -geom hypo lane's arithmetic stays bit-for-bit as
// verified (ADR-79). tests/hypo_incremental_kernel_check.cpp checks the two
// entries against each other and against tests/hypo_incremental_cases.txt.
//
// Returns 0, or -1 when det(F_mid) ≤ 0, det(F_1) ≤ 0 or a polar fails; on -1
// the outputs are left untouched.
inline int hypoIncrementFromF(const double Fn[9], const double F1[9],
                              double deps[6], double R1[9], double* J1)
{
  double Fmid[9], dF[9];
  for (int i = 0; i < 9; i++) {
    Fmid[i] = 0.5 * (Fn[i] + F1[i]);
    dF[i]   = F1[i] - Fn[i];
  }

  double FmidInv[9];
  if (!(inv3(Fmid, FmidInv) > 0.0)) return -1;     // midpoint collapsed/inverted
  const double j1 = det3(F1);
  if (!(j1 > 0.0)) return -1;                      // trial collapsed/inverted

  // G = (F_1 − F_n) · F_mid⁻¹ ; Δε = sym(G)
  double G[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += dF[3*i+k] * FmidInv[3*k+j];
      G[3*i+j] = s;
    }
  double e[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      e[3*i+j] = 0.5 * (G[3*i+j] + G[3*j+i]);

  // De-rotate: em = R_midᵀ · Δε · R_mid
  double Rmid[9];
  if (!polar3(Fmid, Rmid)) return -1;
  double tmp[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += e[3*i+k] * Rmid[3*k+j];
      tmp[3*i+j] = s;
    }
  double em[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += Rmid[3*k+i] * tmp[3*k+j];
      em[3*i+j] = s;
    }

  double R1loc[9];
  if (!polar3(F1, R1loc)) return -1;

  deps[0] = em[0];
  deps[1] = em[4];
  deps[2] = em[8];
  deps[3] = em[1] + em[3];    // 2·ε01
  deps[4] = em[5] + em[7];    // 2·ε12
  deps[5] = em[2] + em[6];    // 2·ε20
  for (int i = 0; i < 9; i++) R1[i] = R1loc[i];
  *J1 = j1;
  return 0;
}

// Plane strain (interface §8.3): embed each 2×2 F as [[F00 F01 0][F10 F11 0]
// [0 0 1]], run the 3D entry, return deps = [xx yy xy] (engineering shear), R1 =
// the upper-left 2×2 block and J1 = the 3×3 determinant. By the block structure
// the out-of-plane strain components are exactly zero and R1 is a rotation about
// z. Same return codes; outputs untouched on -1.
inline int hypoIncrementFromFPlaneStrain(const double Fn[4], const double F1[4],
                                         double deps[3], double R1[4], double* J1)
{
  const double Fn3[9] = {Fn[0], Fn[1], 0.0,  Fn[2], Fn[3], 0.0,  0.0, 0.0, 1.0};
  const double F13[9] = {F1[0], F1[1], 0.0,  F1[2], F1[3], 0.0,  0.0, 0.0, 1.0};
  double d6[6], R3[9], j1;
  if (hypoIncrementFromF(Fn3, F13, d6, R3, &j1) != 0) return -1;
  deps[0] = d6[0];
  deps[1] = d6[1];
  deps[2] = d6[3];
  R1[0] = R3[0]; R1[1] = R3[1];
  R1[2] = R3[3]; R1[3] = R3[4];
  *J1 = j1;
  return 0;
}

// --------------------------------------------------------------------------- //
//  Voigt push-forwards (assembly side)                                        //
// --------------------------------------------------------------------------- //

// 6×6 stress bond (rotation) matrix M for Voigt order {00,11,22,01,12,20}:
//   σ'_voigt = M · σ_voigt   realizes   σ' = R σ Rᵀ.
// The same M push-forwards the ENGINEERING-shear-convention 6×6 tangent as
//   D' = M D Mᵀ  (energy-consistent: strain transforms by M⁻ᵀ, and D maps
//   eng-strain → stress).
inline void bond6(const double R[9], double M[36])
{
  // index pairs for the Voigt slots: (0,0)(1,1)(2,2)(0,1)(1,2)(2,0)
  const int vi[6] = {0, 1, 2, 0, 1, 2};
  const int vj[6] = {0, 1, 2, 1, 2, 0};
  for (int r = 0; r < 6; r++) {
    const int i = vi[r], j = vj[r];
    for (int c = 0; c < 6; c++) {
      const int p = vi[c], q = vj[c];
      // σ'_ij = R_ip R_jq σ_pq ; off-diagonal source slots carry both (p,q)
      // and (q,p) of the symmetric tensor.
      double m = R[3*i+p] * R[3*j+q];
      if (p != q)
        m += R[3*i+q] * R[3*j+p];
      M[6*r + c] = m;
    }
  }
}

// σ_out = R σ Rᵀ in Voigt (both plain-tensor stress vectors). sOut may NOT
// alias sMat.
inline void pushStress6(const double R[9], const double sMat[6], double sOut[6])
{
  double M[36];
  bond6(R, M);
  for (int r = 0; r < 6; r++) {
    double s = 0.0;
    for (int c = 0; c < 6; c++) s += M[6*r+c] * sMat[c];
    sOut[r] = s;
  }
}

// Plane-strain analogue of pushStress6 on [xx, yy, xy] (tensor shear):
// σ_out = R σ Rᵀ with R row-major 2×2 (WP-232, interface §8.4). Exact for an
// in-plane rotation because the unrotated stress has no xz / yz part. sOut may
// NOT alias sMat.
inline void pushStressPlaneStrain(const double R[4], const double sMat[3], double sOut[3])
{
  const double s00 = sMat[0], s11 = sMat[1], s01 = sMat[2];
  // T = R σ
  const double t00 = R[0]*s00 + R[1]*s01, t01 = R[0]*s01 + R[1]*s11;
  const double t10 = R[2]*s00 + R[3]*s01, t11 = R[2]*s01 + R[3]*s11;
  // σ_out = T Rᵀ
  sOut[0] = t00*R[0] + t01*R[1];
  sOut[1] = t10*R[2] + t11*R[3];
  sOut[2] = t00*R[2] + t01*R[3];
}

// D_out = M D Mᵀ — push the material 6×6 tangent (eng-shear Voigt convention)
// to the current frame. Dout may NOT alias Dmat.
inline void pushTangent6(const double R[9], const double Dmat[36], double Dout[36])
{
  double M[36], T[36];
  bond6(R, M);
  for (int r = 0; r < 6; r++)          // T = M · D
    for (int c = 0; c < 6; c++) {
      double s = 0.0;
      for (int k = 0; k < 6; k++) s += M[6*r+k] * Dmat[6*k+c];
      T[6*r+c] = s;
    }
  for (int r = 0; r < 6; r++)          // Dout = T · Mᵀ
    for (int c = 0; c < 6; c++) {
      double s = 0.0;
      for (int k = 0; k < 6; k++) s += T[6*r+k] * M[6*c+k];
      Dout[6*r+c] = s;
    }
}

// Plane-strain analogue of pushTangent6 on {xx, yy, xy} (engineering shear):
// D_out = M3 D M3ᵀ with M3 the {00, 11, 01} rows and columns of bond6 for the
// in-plane rotation R (row-major 2×2) embedded as [[R 0][0 1]] (WP-232f, mode-A
// particle tangent). Exact: for an in-plane R those three rows of bond6 have no
// entry in the {22, 12, 20} columns. Dout may NOT alias Dmat.
inline void pushTangentPlaneStrain(const double R[4], const double Dmat[9], double Dout[9])
{
  const double R3[9] = {R[0], R[1], 0.0,  R[2], R[3], 0.0,  0.0, 0.0, 1.0};
  double M[36];
  bond6(R3, M);
  const int idx[3] = {0, 1, 3};
  double M3[9], T[9];
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++) M3[3*r + c] = M[6*idx[r] + idx[c]];
  for (int r = 0; r < 3; r++)          // T = M3 · D
    for (int c = 0; c < 3; c++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += M3[3*r+k] * Dmat[3*k+c];
      T[3*r+c] = s;
    }
  for (int r = 0; r < 3; r++)          // Dout = T · M3ᵀ
    for (int c = 0; c < 3; c++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += T[3*r+k] * M3[3*c+k];
      Dout[3*r+c] = s;
    }
}

} // namespace ladruno_hypo

#endif
