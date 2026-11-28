// WP-266 (D-UP1): the LadrunoUP -geom hypo element Jacobian by Dual<32>.
//
// Includes ONLY header-only, OpenSees-free code:
//   SRC/element/ladrunoUP/LadrunoUPHypoKernel.h   (the per-GP residual, templated)
//   SRC/element/ladrunoUP/LadrunoUPShapes.h       (the H8 provider the class uses)
//   SRC/material/nD/LadrunoJ2Kernel.h             (a templated plastic material)
//
// Legs, each printing "ok  : ..." or "FAIL: ...":
//
//  A. BIT-IDENTITY of the double path (worker rule 2, the WP-262 pattern).
//     `namespace pre` holds a VERBATIM copy of the pre-WP-266 functions of
//     LadrunoHypoKernel.h (det3, inv3, frob3, dispGradient, polar3,
//     hypoIncrement, bond6, pushStress6, pushTangent6) and of
//     LadrunoUPKernel.h::oneOverQbar, plus a transcription of the pre-WP-266
//     LadrunoUP class code (updateHypo, hypoPorosity, buildCoreForceHypo,
//     buildSolidKHypo, the hypo f_seep branch, the hypoState KC formula). The
//     templated kernels at T = double, and the new class call sequence, must
//     reproduce every cache, force and tangent entry BIT FOR BIT.
//  B. POLAR DERIVATIVE: polar3<Dual<9>> vs the closed form
//     dR = R·[ω]×, (tr U·I − U)·ω = 2·axial(skew(Rᵀ dF)).
//  C. ELEMENT JACOBIAN: elementResidual<Dual<32>> (24 u + 8 p seeded) vs a
//     central FD of elementResidual<double> on a step plateau, on large-strain,
//     compacted (near and at the Kozeny–Carman floor), rotated, small-strain and
//     J2-plastic states. The Dual run's value part must equal the double run bit
//     for bit. The chainStress seam (any NDMaterial: stress + consistent tangent)
//     must equal the native dual material.
//  D. THE DROPPED-BLOCK TABLE: the full Jacobian vs the shipped tangent
//     (shippedTangent, c2 = 1/Δt), block by block, with each dropped dependency
//     attributed by freezing it (Freeze), plus the Newmark c2 = γ/(βΔt) factor.
//  E. BREAK-THE-GATE: the shipped tangent, a Jacobian without the push-forward
//     rotation term and one without the porosity terms must each FAIL the
//     oracle tolerance; the incremental coupling is exactly zero on a rigid
//     increment where the velocity form Q_curᵀ·Δu/Δt is not (ADR 78).
//  F. COST: ns per element residual, double vs Dual<32> (g++ -O2).
//
// Build: g++ -O2 -std=c++17 -I SRC/element/ladrunoUP -I SRC/element/solidTransformation
//            -I SRC/material/nD tests/ladruno_up_hypo_dual_check.cpp -o uphd
// Exit code 0 and "RESULT: PASS" iff every gate holds.

#include <LadrunoUPHypoKernel.h>
#include <LadrunoUPShapes.h>
#include <LadrunoJ2Kernel.h>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>

typedef ladruno_dual::Dual<32> D32;
typedef ladruno_dual::Dual<9>  D9;
namespace UH = ladruno_up_hypo;
using UH::ElemIn;

static const int NN = 8, NP = 8, NG = 8, NU = 24, NDOF = 32;
static const double PI = 3.14159265358979323846;
static int g_fail = 0, g_checks = 0;

static void gate(bool ok, const char* fmt, ...)
{
  g_checks++;
  if (!ok) g_fail++;
  printf(ok ? "ok  : " : "FAIL: ");
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  printf("\n");
}

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
static bool sameArr(const double* a, const double* b, int n)
{
  for (int i = 0; i < n; i++) if (!sameBits(a[i], b[i])) return false;
  return true;
}

// =========================================================================== //
//  A. the pre-WP-266 code, verbatim                                           //
// =========================================================================== //
namespace pre {

// ---- SRC/element/solidTransformation/LadrunoHypoKernel.h @ origin/ladruno 757c6845d
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

// ---- SRC/element/ladrunoUP/LadrunoUPKernel.h @ origin/ladruno 757c6845d
inline double oneOverQbar(double n, double Kf, double biotAlpha, double Ks) {
  double val = n / Kf;
  if (Ks > 0.0)
    val += (biotAlpha - n) / Ks;
  return val;
}

// ---- SRC/element/ladrunoUP/LadrunoUP.cpp @ origin/ladruno 757c6845d,
//      member state -> plain arrays (nNodes_ = 8, nP_ = 8, nGP_ = 8, ndm_ = 3).
struct Elem {
  // inputs
  const ElemIn* in;
  // per-evaluation caches (the class members)
  double hypoFeed_[NG * 6], hypoR_[NG * 9], hypoJ_[NG], hypoTrDeps_[NG], hypoDvCur_[NG];
  double hypoDNuX_[NG * NN * 3], hypoDNpX_[NG * NP * 3], hypoInvQbar_[NG], hypoKsp_[NG * 9];
  double QrotBlk_[NU * NP], HcurBlk_[NP * NP], ScurBlk_[NP * NP];

  double hypoPorosity(int gp) const
  {
    double n = 1.0 - (1.0 - in->n0) / hypoJ_[gp];
    if (n < 1.0e-6)        n = 1.0e-6;
    else if (n > 1.0 - 1.0e-6) n = 1.0 - 1.0e-6;
    return n;
  }

  // LadrunoUP::updateHypo (the material call replaced by keeping the feed)
  int updateHypo(const double* u1)
  {
    const int nNodes_ = in->nNodes, nP_ = in->nP, nGP_ = in->nGP;
    const double *un = in->un;
    const double biotAlpha_ = in->biotAlpha, Kf_ = in->Kf, Ks_ = in->Ks;
    const bool kozenyCarman_ = in->kozenyCarman;
    const double *kbar_ = in->kbar;
    for (int k = 0; k < NU * NP; k++) QrotBlk_[k] = 0.0;
    for (int k = 0; k < NP * NP; k++) { HcurBlk_[k] = 0.0; ScurBlk_[k] = 0.0; }

    int ret = 0;
    for (int gp = 0; gp < nGP_; gp++) {
      const double *dNdX = &in->dNdX[gp * nNodes_ * 3];   // reference gradients

      double deps[6], R1[9], J1;
      if (!pre::hypoIncrement(nNodes_, dNdX, un, u1, deps, R1, J1))
        return -1;

      for (int c = 0; c < 6; c++)
        hypoFeed_[gp * 6 + c] = in->feedCommit[gp * 6 + c] + deps[c];

      for (int c = 0; c < 9; c++)
        hypoR_[gp * 9 + c] = R1[c];
      hypoJ_[gp] = J1;
      hypoTrDeps_[gp] = deps[0] + deps[1] + deps[2];     // tr is frame-invariant
      double dvCur = J1 * in->dV[gp];
      hypoDvCur_[gp] = dvCur;

      // spatial gradients via F_trial⁻¹ (solid + pressure shape sets)
      double H1[9], F1[9], Fi[9];
      pre::dispGradient(nNodes_, dNdX, u1, H1);
      for (int c = 0; c < 9; c++) F1[c] = H1[c];
      F1[0] += 1.0; F1[4] += 1.0; F1[8] += 1.0;
      pre::inv3(F1, Fi);
      double *gU = &hypoDNuX_[gp * nNodes_ * 3];
      for (int a = 0; a < nNodes_; a++)
        for (int j = 0; j < 3; j++) {
          double s = 0.0;
          for (int K = 0; K < 3; K++)
            s += dNdX[a * 3 + K] * Fi[3 * K + j];
          gU[a * 3 + j] = s;
        }
      const double *dNpRef = &in->dNp[gp * nP_ * 3];
      double *gP = &hypoDNpX_[gp * nP_ * 3];
      for (int bb = 0; bb < nP_; bb++)
        for (int j = 0; j < 3; j++) {
          double s = 0.0;
          for (int K = 0; K < 3; K++)
            s += dNpRef[bb * 3 + K] * Fi[3 * K + j];
          gP[bb * 3 + j] = s;
        }

      // porosity-evolved storage coefficient + spatial permeability
      double n = this->hypoPorosity(gp);
      hypoInvQbar_[gp] = (Kf_ > 0.0)
          ? pre::oneOverQbar(n, Kf_, biotAlpha_, Ks_) : 0.0;
      double kc = 1.0;
      if (kozenyCarman_) {
        double n0 = in->n0;
        kc = (n * n * n / ((1.0 - n) * (1.0 - n)))
           * ((1.0 - n0) * (1.0 - n0) / (n0 * n0 * n0));
        if (kc < 1.0e-6) kc = 1.0e-6;
      }
      double ksp[9];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          double s = 0.0;
          for (int m = 0; m < 3; m++)
            s += R1[3 * i + m] * kbar_[m] * R1[3 * j + m];
          ksp[3 * i + j] = kc * s;
        }
      for (int c = 0; c < 9; c++)
        hypoKsp_[gp * 9 + c] = ksp[c];

      const double *Np = &in->Np[gp * nP_];
      for (int a = 0; a < nNodes_; a++)
        for (int i = 0; i < 3; i++) {
          int row = a * 3 + i;
          double g = gU[row] * biotAlpha_ * dvCur;
          for (int bb = 0; bb < nP_; bb++)
            QrotBlk_[row * nP_ + bb] += g * Np[bb];
        }
      for (int bb = 0; bb < nP_; bb++)
        for (int cc = 0; cc < nP_; cc++) {
          double s = 0.0;
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
              s += gP[bb * 3 + i] * ksp[3 * i + j] * gP[cc * 3 + j];
          HcurBlk_[bb * nP_ + cc] += s * dvCur;
          ScurBlk_[bb * nP_ + cc] += Np[bb] * hypoInvQbar_[gp] * Np[cc] * dvCur;
        }
    }
    return ret;
  }

  // LadrunoUP::buildCoreForceHypo (sMatAll: the material stress per GP)
  void buildCoreForceHypo(const double* sMatAll, const double* pScratch_, double* fu)
  {
    const int nNodes_ = in->nNodes, nP_ = in->nP, nGP_ = in->nGP;
    const double biotAlpha_ = in->biotAlpha;
    for (int k = 0; k < nNodes_ * 3; k++) fu[k] = 0.0;
    for (int gp = 0; gp < nGP_; gp++) {
      double sMat[6], sSp[6];
      for (int c = 0; c < 6; c++) sMat[c] = sMatAll[gp * 6 + c];
      pre::pushStress6(&hypoR_[gp * 9], sMat, sSp);

      double p = 0.0;
      const double *Np = &in->Np[gp * nP_];
      for (int bb = 0; bb < nP_; bb++)
        p += Np[bb] * pScratch_[bb];

      double sig[3][3];
      sig[0][0] = sSp[0] - biotAlpha_ * p;
      sig[1][1] = sSp[1] - biotAlpha_ * p;
      sig[2][2] = sSp[2] - biotAlpha_ * p;
      sig[0][1] = sig[1][0] = sSp[3];
      sig[1][2] = sig[2][1] = sSp[4];
      sig[2][0] = sig[0][2] = sSp[5];

      const double *gU = &hypoDNuX_[gp * nNodes_ * 3];
      double dvCur = hypoDvCur_[gp];
      for (int a = 0; a < nNodes_; a++)
        for (int i = 0; i < 3; i++) {
          double s = 0.0;
          for (int j = 0; j < 3; j++)
            s += sig[i][j] * gU[a * 3 + j];
          fu[a * 3 + i] += s * dvCur;
        }
    }
  }

  // LadrunoUP::buildSolidKHypo (Ks row-major nU×nU; Matrix(i,j) -> Ks[i*nU+j])
  void buildSolidKHypo(const double* sMatAll, const double* DmatAll, const double* pNodal,
                       double* Ks)
  {
    const int nNodes_ = in->nNodes, nP_ = in->nP, nGP_ = in->nGP, nU = 3 * nNodes_;
    const double biotAlpha_ = in->biotAlpha;
    for (int k = 0; k < nU * nU; k++) Ks[k] = 0.0;
    for (int gp = 0; gp < nGP_; gp++) {
      double Dmat[36], Dcur[36];
      for (int r = 0; r < 6; r++)
        for (int c = 0; c < 6; c++)
          Dmat[6 * r + c] = DmatAll[gp * 36 + 6 * r + c];
      pre::pushTangent6(&hypoR_[gp * 9], Dmat, Dcur);

      const double *gU = &hypoDNuX_[gp * nNodes_ * 3];
      double dvCur = hypoDvCur_[gp];

      double Ba[6][3], Bb[6][3], DB[6][3];
      for (int a = 0; a < nNodes_; a++) {
        for (int r = 0; r < 6; r++)
          for (int c = 0; c < 3; c++) Ba[r][c] = 0.0;
        Ba[0][0] = gU[a*3+0]; Ba[1][1] = gU[a*3+1]; Ba[2][2] = gU[a*3+2];
        Ba[3][0] = gU[a*3+1]; Ba[3][1] = gU[a*3+0];
        Ba[4][1] = gU[a*3+2]; Ba[4][2] = gU[a*3+1];
        Ba[5][0] = gU[a*3+2]; Ba[5][2] = gU[a*3+0];

        for (int bn = 0; bn < nNodes_; bn++) {
          for (int r = 0; r < 6; r++)
            for (int c = 0; c < 3; c++) Bb[r][c] = 0.0;
          Bb[0][0] = gU[bn*3+0]; Bb[1][1] = gU[bn*3+1]; Bb[2][2] = gU[bn*3+2];
          Bb[3][0] = gU[bn*3+1]; Bb[3][1] = gU[bn*3+0];
          Bb[4][1] = gU[bn*3+2]; Bb[4][2] = gU[bn*3+1];
          Bb[5][0] = gU[bn*3+2]; Bb[5][2] = gU[bn*3+0];

          for (int r = 0; r < 6; r++)
            for (int c = 0; c < 3; c++) {
              double s = 0.0;
              for (int q = 0; q < 6; q++)
                s += Dcur[6 * r + q] * Bb[q][c];
              DB[r][c] = s;
            }
          for (int i = 0; i < 3; i++)
            for (int k = 0; k < 3; k++) {
              double s = 0.0;
              for (int r = 0; r < 6; r++)
                s += Ba[r][i] * DB[r][k];
              Ks[(3 * a + i) * nU + 3 * bn + k] += s * dvCur;
            }
        }
      }

      double sMat[6], sSp[6];
      for (int c = 0; c < 6; c++) sMat[c] = sMatAll[gp * 6 + c];
      pre::pushStress6(&hypoR_[gp * 9], sMat, sSp);
      double p = 0.0;                                    // gpPressure(gp)
      const double *Np = &in->Np[gp * nP_];
      for (int bb = 0; bb < nP_; bb++)
        p += Np[bb] * pNodal[bb];
      double sig[3][3];
      sig[0][0] = sSp[0] - biotAlpha_ * p;
      sig[1][1] = sSp[1] - biotAlpha_ * p;
      sig[2][2] = sSp[2] - biotAlpha_ * p;
      sig[0][1] = sig[1][0] = sSp[3];
      sig[1][2] = sig[2][1] = sSp[4];
      sig[2][0] = sig[0][2] = sSp[5];

      for (int a = 0; a < nNodes_; a++)
        for (int bn = 0; bn < nNodes_; bn++) {
          double s = 0.0;
          for (int j = 0; j < 3; j++)
            for (int l = 0; l < 3; l++)
              s += gU[a * 3 + j] * sig[j][l] * gU[bn * 3 + l];
          s *= dvCur;
          for (int i = 0; i < 3; i++)
            Ks[(3 * a + i) * nU + 3 * bn + i] += s;
        }
    }
  }

  // the hypo branch of LadrunoUP::getResistingForce's f_seep loop
  void fseepHypo(double* fseep)
  {
    const int nP_ = in->nP, nGP_ = in->nGP;
    const double rhoF_ = in->rhoF;
    for (int bb = 0; bb < nP_; bb++) fseep[bb] = 0.0;
    for (int gp = 0; gp < nGP_; gp++) {
      const double *drive = in->drive;
      const double *gP = &hypoDNpX_[gp * nP_ * 3];
      const double *ksp = &hypoKsp_[gp * 9];
      for (int bb = 0; bb < nP_; bb++) {
        double s = 0.0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            s += gP[bb * 3 + i] * ksp[3 * i + j] * rhoF_ * drive[j];
        fseep[bb] += s * hypoDvCur_[gp];
      }
    }
  }

  // the p-rows of getResistingForce + getResistingForceIncInertia under hypo:
  //   (H·p − f_seep) + [Qᵀ·u̇ + (S+αH̃)·ṗ] + (sNew/Δt − sOld)
  // with M·ü = 0 and C_Ray = 0 (u-rows: fu exactly). vel: block order.
  void residual(const double* fu, const double* p1, const double* vel, double* r)
  {
    const int nP_ = in->nP, nGP_ = in->nGP, nU = 3 * in->nNodes;
    double fseep[NP];
    fseepHypo(fseep);
    for (int k = 0; k < nU + nP_; k++) r[k] = 0.0;
    for (int k = 0; k < nU; k++) r[k] += fu[k];
    for (int bb = 0; bb < nP_; bb++) {
      double s = 0.0;
      for (int cc = 0; cc < nP_; cc++)
        s += HcurBlk_[bb * nP_ + cc] * p1[cc];
      r[nU + bb] += s;
    }
    for (int bb = 0; bb < nP_; bb++)
      r[nU + bb] -= fseep[bb];
    // + getDamp()·v, p-rows: Qᵀ (QrotBlk_) and S + al·H̃
    for (int bb = 0; bb < nP_; bb++) {
      double s = 0.0;
      for (int k = 0; k < nU; k++) s += QrotBlk_[k * nP_ + bb] * vel[k];
      for (int cc = 0; cc < nP_; cc++) {
        double ht = in->HtAl ? in->HtAl[bb * nP_ + cc] : 0.0;
        s += (ScurBlk_[bb * nP_ + cc] + ht) * vel[nU + cc];
      }
      r[nU + bb] += s;
    }
    double dt = in->dt;
    if (dt > 0.0) {
      for (int bb = 0; bb < nP_; bb++) {
        double sOld = 0.0, sNew = 0.0;
        for (int k = 0; k < nU; k++)
          sOld += QrotBlk_[k * nP_ + bb] * vel[k];
        for (int gp = 0; gp < nGP_; gp++)
          sNew += in->Np[gp * nP_ + bb] * in->biotAlpha
                * hypoTrDeps_[gp] * hypoDvCur_[gp];
        for (int cc = 0; cc < nP_; cc++) {
          double ht = in->HtAl ? in->HtAl[bb * nP_ + cc] : 0.0;
          double spc = ScurBlk_[bb * nP_ + cc] + ht;
          sOld += spc * vel[nU + cc];
          sNew += spc * (p1[cc] - in->pn[cc]);
        }
        r[nU + bb] += sNew / dt - sOld;
      }
    }
  }
};

// the hypoState response's KC scale (LadrunoUP::getResponse, responseID 12)
inline double kcResponse(double n, double poro_)
{
  double kc = 1.0;
  double n0 = poro_;
  kc = (n * n * n / ((1.0 - n) * (1.0 - n)))
     * ((1.0 - n0) * (1.0 - n0) / (n0 * n0 * n0));
  if (kc < 1.0e-6) kc = 1.0e-6;
  return kc;
}

} // namespace pre

// =========================================================================== //
//  The NEW class call sequence (LadrunoUP.cpp after WP-266), on arrays         //
// =========================================================================== //
struct NewElem {
  const ElemIn* in;
  double hypoFeed_[NG * 6], hypoR_[NG * 9], hypoJ_[NG], hypoTrDeps_[NG], hypoDvCur_[NG];
  double hypoDNuX_[NG * NN * 3], hypoDNpX_[NG * NP * 3], hypoInvQbar_[NG], hypoKsp_[NG * 9];
  double QrotBlk_[NU * NP], HcurBlk_[NP * NP], ScurBlk_[NP * NP];

  int updateHypo(const double* u1)
  {
    const int nNodes_ = in->nNodes, nP_ = in->nP, nGP_ = in->nGP;
    for (int k = 0; k < NU * NP; k++) QrotBlk_[k] = 0.0;
    for (int k = 0; k < NP * NP; k++) { HcurBlk_[k] = 0.0; ScurBlk_[k] = 0.0; }
    for (int gp = 0; gp < nGP_; gp++) {
      const double *dNdX = &in->dNdX[gp * nNodes_ * 3];
      const double *dNpRef = &in->dNp[gp * nP_ * 3];
      double *gU = &hypoDNuX_[gp * nNodes_ * 3];
      double *gP = &hypoDNpX_[gp * nP_ * 3];
      double deps[6], R1[9], J1, trDeps, dvCur;
      if (!UH::gpKinematics(nNodes_, nP_, dNdX, dNpRef, in->dV[gp], in->un, u1,
                            deps, R1, J1, trDeps, dvCur, gU, gP))
        return -1;
      for (int c = 0; c < 6; c++)
        hypoFeed_[gp * 6 + c] = in->feedCommit[gp * 6 + c] + deps[c];
      for (int c = 0; c < 9; c++) hypoR_[gp * 9 + c] = R1[c];
      hypoJ_[gp] = J1;
      hypoTrDeps_[gp] = trDeps;
      hypoDvCur_[gp] = dvCur;
      double ksp[9];
      UH::gpPoroPerm(J1, R1, in->n0, in->Kf, in->biotAlpha, in->Ks, in->kbar,
                     in->kozenyCarman, hypoInvQbar_[gp], ksp);
      for (int c = 0; c < 9; c++) hypoKsp_[gp * 9 + c] = ksp[c];
      const double *Np = &in->Np[gp * nP_];
      UH::addQcur(nNodes_, nP_, gU, Np, in->biotAlpha, dvCur, QrotBlk_);
      UH::addHScur(nP_, gP, ksp, Np, hypoInvQbar_[gp], dvCur, HcurBlk_, ScurBlk_);
    }
    return 0;
  }
  void buildCoreForceHypo(const double* sMatAll, const double* pNodal, double* fu)
  {
    for (int k = 0; k < in->nNodes * 3; k++) fu[k] = 0.0;
    for (int gp = 0; gp < in->nGP; gp++) {
      double sMat[6];
      for (int c = 0; c < 6; c++) sMat[c] = sMatAll[gp * 6 + c];
      double p = UH::gpPressure(in->nP, &in->Np[gp * in->nP], pNodal);
      UH::addForceUHypo(in->nNodes, &hypoR_[gp * 9], sMat, p,
                        &hypoDNuX_[gp * in->nNodes * 3], hypoDvCur_[gp],
                        in->biotAlpha, fu);
    }
  }
  void buildSolidKHypo(const double* sMatAll, const double* DmatAll, const double* pNodal,
                       double* Ks)
  {
    const int nU = in->nNodes * 3;
    double ksBuf[UH::MAXU * UH::MAXU];
    for (int k = 0; k < nU * nU; k++) ksBuf[k] = 0.0;
    for (int gp = 0; gp < in->nGP; gp++) {
      double Dmat[36], sMat[6];
      for (int c = 0; c < 36; c++) Dmat[c] = DmatAll[gp * 36 + c];
      for (int c = 0; c < 6; c++) sMat[c] = sMatAll[gp * 6 + c];
      double p = UH::gpPressure(in->nP, &in->Np[gp * in->nP], pNodal);
      UH::addSolidKHypoGP(in->nNodes, &hypoR_[gp * 9], Dmat, sMat, p, in->biotAlpha,
                          &hypoDNuX_[gp * in->nNodes * 3], hypoDvCur_[gp], ksBuf, nU);
    }
    for (int i = 0; i < nU; i++)
      for (int j = 0; j < nU; j++) Ks[i * nU + j] = ksBuf[i * nU + j];
  }
  void fseepHypo(double* fseep)
  {
    for (int bb = 0; bb < in->nP; bb++) fseep[bb] = 0.0;
    for (int gp = 0; gp < in->nGP; gp++)
      UH::addFseepCur(in->nP, &hypoDNpX_[gp * in->nP * 3], &hypoKsp_[gp * 9], in->rhoF,
                      in->drive, hypoDvCur_[gp], fseep);
  }
};

// =========================================================================== //
//  Materials                                                                  //
// =========================================================================== //

// Isotropic linear elastic in the unrotated frame with a prestress:
// σ̂ = σ̂₀ + D ε_feed (engineering shear).
struct ElasticMat {
  double D[36], s0[6];
  void init(double E, double nu, const double pre[6])
  {
    const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), mu = E / (2.0 * (1.0 + nu));
    for (int k = 0; k < 36; k++) D[k] = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) D[6 * i + j] = lam;
      D[6 * i + i] += 2.0 * mu;
      D[6 * (i + 3) + i + 3] = mu;
    }
    for (int i = 0; i < 6; i++) s0[i] = pre[i];
  }
  template <class T> void operator()(int, const T feed[6], T sig[6]) const
  {
    for (int i = 0; i < 6; i++) {
      T s = s0[i];
      for (int j = 0; j < 6; j++) s += D[6 * i + j] * feed[j];
      sig[i] = s;
    }
  }
  void tangent(int, const double feed[6], double sig[6], double Dout[36]) const
  {
    (*this)(0, feed, sig);
    for (int k = 0; k < 36; k++) Dout[k] = D[k];
  }
};

// LadrunoJ2 kernel (templated return map, WP-262) from a virgin committed
// state: the NDMaterial's engineering-shear strain -> tensor strain for the
// kernel; its Dtan is already d(stress)/d(engineering strain).
struct J2Mat {
  ladruno_j2_kernel::Params p;
  template <class T> void operator()(int, const T feed[6], T sig[6]) const
  {
    using namespace ladruno_j2_kernel;
    T eps[6], epsPn[6], ebarPn(0.0), alphaN[MAXBACK][6];
    for (int i = 0; i < 6; i++) {
      eps[i] = (i < 3) ? feed[i] : 0.5 * feed[i];
      epsPn[i] = 0.0;
    }
    for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alphaN[k][i] = 0.0;
    T Dt[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG;
    returnMapT<T>(p, eps, epsPn, ebarPn, alphaN, sig, Dt, epsP, ebarP, alpha, dG);
  }
  void tangent(int, const double feed[6], double sig[6], double Dout[36], double* dGout = 0) const
  {
    using namespace ladruno_j2_kernel;
    double eps[6], epsPn[6] = {0, 0, 0, 0, 0, 0}, alphaN[MAXBACK][6];
    for (int i = 0; i < 6; i++) eps[i] = (i < 3) ? feed[i] : 0.5 * feed[i];
    for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alphaN[k][i] = 0.0;
    double Dt[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG;
    returnMap(p, eps, epsPn, 0.0, alphaN, sig, Dt, epsP, ebarP, alpha, dG);
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) Dout[6 * i + j] = Dt[i][j];
    if (dGout) *dGout = dG;
  }
};

// The chainStress seam: any double material (stress + consistent tangent) as a
// T-valued response.
template <class M>
struct ChainMat {
  const M* m;
  template <class T> void operator()(int gp, const T feed[6], T sig[6]) const
  {
    double fv[6], s[6], D[36];
    for (int i = 0; i < 6; i++) fv[i] = ladruno_dual::value(feed[i]);
    m->tangent(gp, fv, s, D);
    UH::chainStress(s, D, feed, sig);
  }
};

// =========================================================================== //
//  Element fixture                                                            //
// =========================================================================== //
static const double HEX[24] = {   // the distorted hex of test_ladruno_up_element_hypo.py
  0.00, 0.00, 0.00,  1.00, 0.10, 0.05,  1.10, 1.00, 0.00,  0.05, 0.95, 0.10,
  0.00, 0.05, 1.00,  1.00, 0.00, 1.05,  1.05, 1.00, 1.10,  0.00, 1.00, 0.95 };

struct Fixture {
  double dNdX[NG * NN * 3], Np[NG * NP], dNp[NG * NP * 3], dV[NG];
  double un[NU], u1[NU], pn[NP], p1[NP], feedCommit[NG * 6], HtAl[NP * NP];
  ElemIn in;

  void shapes()
  {
    using ladruno_up::shapes::H8;
    for (int gp = 0; gp < NG; gp++) {
      double Nu[8], dNu[24], detJ;
      H8::eval(gp, HEX, Nu, dNu, &detJ);
      dV[gp] = H8::gpWeight[gp] * std::fabs(detJ);
      for (int a = 0; a < 8; a++) {
        Np[gp * NP + a] = Nu[a];
        for (int i = 0; i < 3; i++) {
          dNdX[gp * 24 + a * 3 + i] = dNu[a * 3 + i];
          dNp[gp * 24 + a * 3 + i] = dNu[a * 3 + i];
        }
      }
    }
    // the unit-alpha stabilization Laplacian H~ (LadrunoUPKernel addHtilde), x alpha
    double Ht[NP * NP];
    for (int k = 0; k < NP * NP; k++) Ht[k] = 0.0;
    for (int gp = 0; gp < NG; gp++)
      ladruno_up::addHtilde(&dNp[gp * 24], NP, 3, 1.0, dV[gp], Ht, NP);
    for (int k = 0; k < NP * NP; k++) HtAl[k] = 2.0e-6 * Ht[k];
  }
};

// nodal displacement of the map x = c + A X + q·(X·X) e_z, at node a
static void mapDisp(const double A[9], const double c[3], double q, double* u)
{
  for (int a = 0; a < NN; a++) {
    const double* X = &HEX[3 * a];
    const double xx = X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
    for (int i = 0; i < 3; i++) {
      double x = c[i] + A[3 * i] * X[0] + A[3 * i + 1] * X[1] + A[3 * i + 2] * X[2];
      if (i == 2) x += q * xx;
      u[3 * a + i] = x - X[i];
    }
  }
}

static void rotAxis(double ax, double ay, double az, double th, double R[9])
{
  const double n = std::sqrt(ax * ax + ay * ay + az * az);
  ax /= n; ay /= n; az /= n;
  const double c = std::cos(th), s = std::sin(th), C = 1.0 - c;
  R[0] = c + ax*ax*C;      R[1] = ax*ay*C - az*s; R[2] = ax*az*C + ay*s;
  R[3] = ay*ax*C + az*s;   R[4] = c + ay*ay*C;    R[5] = ay*az*C - ax*s;
  R[6] = az*ax*C - ay*s;   R[7] = az*ay*C + ax*s; R[8] = c + az*az*C;
}
static void mat3mul(const double A[9], const double B[9], double C[9])
{
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
    double s = 0.0; for (int k = 0; k < 3; k++) s += A[3*i+k] * B[3*k+j]; C[3*i+j] = s; }
}

struct State {
  const char* name;
  double An[9], A1[9], cn[3], c1[3], qn, q1;
  double pBase, pAmp, dpAmp;
  double feedAmp;
  double kbar[3];
  bool kozeny;
  double dt;
  int material;         // 0 elastic, 1 J2
};

static void setup(Fixture& f, const State& st)
{
  f.shapes();
  mapDisp(st.An, st.cn, st.qn, f.un);
  mapDisp(st.A1, st.c1, st.q1, f.u1);
  for (int b = 0; b < NP; b++) {
    const double* X = &HEX[3 * b];
    f.p1[b] = st.pBase + st.pAmp * (0.3 * X[0] - 0.5 * X[1] + 0.8 * X[2] + 0.1 * b);
    f.pn[b] = f.p1[b] - st.dpAmp * (1.0 + 0.2 * X[2] - 0.1 * X[0]);
  }
  for (int gp = 0; gp < NG; gp++)
    for (int c = 0; c < 6; c++)
      f.feedCommit[gp * 6 + c] = st.feedAmp * ((c < 3) ? (0.3 - 0.2 * c) : 0.1 * (c - 2))
                                 * (1.0 + 0.05 * gp);
  ElemIn& in = f.in;
  in.nNodes = NN; in.nP = NP; in.nGP = NG;
  in.dNdX = f.dNdX; in.Np = f.Np; in.dNp = f.dNp; in.dV = f.dV;
  in.un = f.un; in.pn = f.pn; in.feedCommit = f.feedCommit; in.HtAl = f.HtAl;
  in.biotAlpha = 0.9; in.Kf = 2.2e6; in.Ks = 3.0e7; in.n0 = 0.4;
  for (int i = 0; i < 3; i++) in.kbar[i] = st.kbar[i];
  in.kozenyCarman = st.kozeny;
  in.rhoF = 1.0;
  in.drive[0] = 0.0; in.drive[1] = 0.0; in.drive[2] = -9.81;
  in.dt = st.dt;
}

// =========================================================================== //
//  Jacobian helpers                                                           //
// =========================================================================== //
template <class M>
static bool dualJacobian(const ElemIn& in, const double* u1, const double* p1, M& mat,
                         double J[NDOF * NDOF], double rval[NDOF],
                         UH::ElemParts<D32>* parts = 0, const UH::Freeze& fz = UH::Freeze())
{
  D32 uD[NU], pD[NP], rD[NDOF];
  for (int k = 0; k < NU; k++) uD[k] = D32::variable(u1[k], k);
  for (int k = 0; k < NP; k++) pD[k] = D32::variable(p1[k], NU + k);
  if (!UH::elementResidual(in, uD, pD, mat, rD, parts, fz)) return false;
  for (int i = 0; i < NDOF; i++) {
    rval[i] = rD[i].v;
    for (int j = 0; j < NDOF; j++) J[i * NDOF + j] = rD[i].d[j];
  }
  return true;
}

// Jacobian of one residual part (u-columns only: nRows x 24)
static void partJac(const D32* part, int nRows, double* Jp)
{
  for (int i = 0; i < nRows; i++)
    for (int j = 0; j < NU; j++) Jp[i * NU + j] = part[i].d[j];
}

// block Frobenius norm: rows [r0,r1), cols [c0,c1) of an n x n matrix, or of A-B
static double blkNorm(const double* A, const double* B, int n, int r0, int r1, int c0, int c1)
{
  double s = 0.0;
  for (int i = r0; i < r1; i++)
    for (int j = c0; j < c1; j++) {
      double d = A[i * n + j] - (B ? B[i * n + j] : 0.0);
      s += d * d;
    }
  return std::sqrt(s);
}

// Central FD of the double residual, per column on a step plateau; returns the
// worst over columns of min-over-steps of the block-normalized error against
// Jref. errBlk[0..3] = uu, up, pu, pp. bestH (optional): the plateau step per col.
template <class M>
static void fdCompare(const ElemIn& in, const double* u1, const double* p1, M& mat,
                      const double* Jref, double errBlk[4], double* Jfd = 0)
{
  double nrm[4];
  nrm[0] = blkNorm(Jref, 0, NDOF, 0, NU, 0, NU);
  nrm[1] = blkNorm(Jref, 0, NDOF, 0, NU, NU, NDOF);
  nrm[2] = blkNorm(Jref, 0, NDOF, NU, NDOF, 0, NU);
  nrm[3] = blkNorm(Jref, 0, NDOF, NU, NDOF, NU, NDOF);
  // max-abs block scales (per-entry error normalization)
  double mx[4] = {0, 0, 0, 0};
  for (int i = 0; i < NDOF; i++)
    for (int j = 0; j < NDOF; j++) {
      int b = (i < NU ? 0 : 2) + (j < NU ? 0 : 1);
      mx[b] = std::fmax(mx[b], std::fabs(Jref[i * NDOF + j]));
    }
  (void)nrm;
  for (int b = 0; b < 4; b++) errBlk[b] = 0.0;
  double pScale = 0.0;
  for (int k = 0; k < NP; k++) pScale = std::fmax(pScale, std::fabs(p1[k]));
  if (pScale == 0.0) pScale = 1.0;
  const double hs[7] = {1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6};

  for (int j = 0; j < NDOF; j++) {
    double best[2] = {1e300, 1e300};    // per row block: u-rows, p-rows
    double bestCol[NDOF];
    double bestTot = 1e300;
    for (int ih = 0; ih < 7; ih++) {
      double up[NU], um[NU], pp[NP], pm[NP];
      std::memcpy(up, u1, sizeof up); std::memcpy(um, u1, sizeof um);
      std::memcpy(pp, p1, sizeof pp); std::memcpy(pm, p1, sizeof pm);
      double h;
      if (j < NU) { h = hs[ih]; up[j] += h; um[j] -= h; }
      else        { h = hs[ih] * pScale; pp[j - NU] += h; pm[j - NU] -= h; }
      double rp[NDOF], rm[NDOF];
      if (!UH::elementResidual(in, up, pp, mat, rp) || !UH::elementResidual(in, um, pm, mat, rm))
        continue;
      double e[2] = {0, 0}, col[NDOF];
      for (int i = 0; i < NDOF; i++) {
        col[i] = (rp[i] - rm[i]) / (2.0 * h);
        int b = (i < NU ? 0 : 2) + (j < NU ? 0 : 1);
        double sc = mx[b] > 0.0 ? mx[b] : 1.0;
        double d = std::fabs(col[i] - Jref[i * NDOF + j]) / sc;
        e[i < NU ? 0 : 1] = std::fmax(e[i < NU ? 0 : 1], d);
      }
      double tot = std::fmax(e[0], e[1]);
      if (tot < bestTot) {
        bestTot = tot; best[0] = e[0]; best[1] = e[1];
        std::memcpy(bestCol, col, sizeof bestCol);
      }
    }
    int cb = (j < NU) ? 0 : 1;
    errBlk[0 + cb] = std::fmax(errBlk[0 + cb], best[0]);
    errBlk[2 + cb] = std::fmax(errBlk[2 + cb], best[1]);
    if (Jfd) for (int i = 0; i < NDOF; i++) Jfd[i * NDOF + j] = bestCol[i];
  }
}

// the per-entry block error of a candidate Jacobian Jc against a reference Jr
static void blkErr(const double* Jc, const double* Jr, double err[4])
{
  double mx[4] = {0, 0, 0, 0}, e[4] = {0, 0, 0, 0};
  for (int i = 0; i < NDOF; i++)
    for (int j = 0; j < NDOF; j++) {
      int b = (i < NU ? 0 : 2) + (j < NU ? 0 : 1);
      mx[b] = std::fmax(mx[b], std::fabs(Jr[i * NDOF + j]));
    }
  for (int i = 0; i < NDOF; i++)
    for (int j = 0; j < NDOF; j++) {
      int b = (i < NU ? 0 : 2) + (j < NU ? 0 : 1);
      e[b] = std::fmax(e[b], std::fabs(Jc[i * NDOF + j] - Jr[i * NDOF + j]) /
                                 (mx[b] > 0.0 ? mx[b] : 1.0));
    }
  for (int b = 0; b < 4; b++) err[b] = e[b];
}

template <class M>
static void materialAtTrial(const ElemIn& in, const double* u1, const M& mat,
                            double sMat[NG * 6], double Dmat[NG * 36])
{
  for (int gp = 0; gp < in.nGP; gp++) {
    double deps[6], R1[9], J1, tr, dv, gU[NU], gP[NP * 3];
    UH::gpKinematics(in.nNodes, in.nP, in.dNdX + gp * 24, in.dNp + gp * 24, in.dV[gp],
                     in.un, u1, deps, R1, J1, tr, dv, gU, gP);
    double feed[6];
    for (int c = 0; c < 6; c++) feed[c] = in.feedCommit[gp * 6 + c] + deps[c];
    mat.tangent(gp, feed, &sMat[gp * 6], &Dmat[gp * 36]);
  }
}

// =========================================================================== //
//  A. bit-identity                                                            //
// =========================================================================== //
static unsigned long long g_rng = 0x9E3779B97F4A7C15ull;
static double urand()   // [-1, 1)
{
  g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17;
  return ((g_rng >> 11) * (1.0 / 9007199254740992.0)) * 2.0 - 1.0;
}

static void legA_kernel()
{
  // the shared hypo kernel at T = double vs its verbatim pre-WP copy
  int nCase = 0;
  bool okInc = true, okPol = true, okPush = true, okTan = true, okInv = true;
  for (int c = 0; c < 400; c++) {
    double dNdX[24], un[24], u1[24];
    for (int k = 0; k < 24; k++) {
      dNdX[k] = 0.6 * urand();
      un[k] = 0.15 * urand();
      u1[k] = un[k] + ((c % 4 == 0) ? 0.3 : 0.05) * urand();
    }
    double dA[6], RA[9], JA, dB[6], RB[9], JB;
    bool a = pre::hypoIncrement(8, dNdX, un, u1, dA, RA, JA);
    bool b = ladruno_hypo::hypoIncrement(8, dNdX, un, u1, dB, RB, JB);
    if (a != b) okInc = false;
    if (a && b) {
      nCase++;
      okInc = okInc && sameArr(dA, dB, 6) && sameArr(RA, RB, 9) && sameBits(JA, JB);
    }
    double F[9], Ra[9], Rb[9], Fi1[9], Fi2[9];
    for (int k = 0; k < 9; k++) F[k] = ((k % 4 == 0) ? 1.0 : 0.0) + 0.4 * urand();
    bool pa = pre::polar3(F, Ra), pb = ladruno_hypo::polar3(F, Rb);
    okPol = okPol && (pa == pb) && (!pa || sameArr(Ra, Rb, 9));
    double da = pre::inv3(F, Fi1), db = ladruno_hypo::inv3(F, Fi2);
    okInv = okInv && sameBits(da, db) && sameArr(Fi1, Fi2, 9)
                  && sameBits(pre::det3(F), ladruno_hypo::det3(F))
                  && sameBits(pre::frob3(F), ladruno_hypo::frob3(F));
    double s6[6], o1[6], o2[6], D[36], Do1[36], Do2[36];
    for (int k = 0; k < 6; k++) s6[k] = 100.0 * urand();
    for (int k = 0; k < 36; k++) D[k] = 1000.0 * urand();
    if (pa) {
      pre::pushStress6(Ra, s6, o1); ladruno_hypo::pushStress6(Ra, s6, o2);
      okPush = okPush && sameArr(o1, o2, 6);
      pre::pushTangent6(Ra, D, Do1); ladruno_hypo::pushTangent6(Ra, D, Do2);
      okTan = okTan && sameArr(Do1, Do2, 36);
    }
  }
  gate(okInc && nCase > 300, "A1 hypoIncrement<double> == pre-WP hypoIncrement bit for bit (%d cases)", nCase);
  gate(okPol, "A1 polar3<double> == pre-WP polar3 bit for bit (400 cases)");
  gate(okInv, "A1 det3 / inv3 / frob3 <double> == pre-WP bit for bit");
  gate(okPush, "A1 pushStress6<double> == pre-WP bit for bit");
  gate(okTan, "A1 pushTangent6 (now on bond6<double>) == pre-WP bit for bit");
  bool okQ = true;
  for (int c = 0; c < 200; c++) {
    double n = 0.5 * (urand() + 1.0), Kf = 1e6 * (2.0 + urand()), al = 0.5 + 0.4 * urand();
    double Ks = (c % 3 == 0) ? -1.0 : 1e7 * (1.5 + urand());
    okQ = okQ && sameBits(pre::oneOverQbar(n, Kf, al, Ks), ladruno_up::oneOverQbar(n, Kf, al, Ks));
    double nn = 0.05 + 0.5 * (urand() + 1.0) * 0.4, n0 = 0.3 + 0.1 * urand();
    okQ = okQ && sameBits(pre::kcResponse(nn, n0), UH::kozenyCarmanScale(nn, n0));
  }
  gate(okQ, "A1 oneOverQbar<double> and the KC response formula == pre-WP bit for bit (200 cases)");
}

template <class M>
static void legA_class(const char* tag, Fixture& f, const M& mat)
{
  pre::Elem A; A.in = &f.in;
  NewElem B;   B.in = &f.in;
  int ra = A.updateHypo(f.u1), rb = B.updateHypo(f.u1);
  bool ok = (ra == rb) && ra == 0;
  ok = ok && sameArr(A.hypoFeed_, B.hypoFeed_, NG * 6) && sameArr(A.hypoR_, B.hypoR_, NG * 9)
          && sameArr(A.hypoJ_, B.hypoJ_, NG) && sameArr(A.hypoTrDeps_, B.hypoTrDeps_, NG)
          && sameArr(A.hypoDvCur_, B.hypoDvCur_, NG)
          && sameArr(A.hypoDNuX_, B.hypoDNuX_, NG * NN * 3)
          && sameArr(A.hypoDNpX_, B.hypoDNpX_, NG * NP * 3)
          && sameArr(A.hypoInvQbar_, B.hypoInvQbar_, NG) && sameArr(A.hypoKsp_, B.hypoKsp_, NG * 9)
          && sameArr(A.QrotBlk_, B.QrotBlk_, NU * NP) && sameArr(A.HcurBlk_, B.HcurBlk_, NP * NP)
          && sameArr(A.ScurBlk_, B.ScurBlk_, NP * NP);
  bool okPor = true;
  for (int gp = 0; gp < NG; gp++)
    okPor = okPor && sameBits(A.hypoPorosity(gp), UH::porosity(B.hypoJ_[gp], f.in.n0));
  gate(ok, "A2 [%s] updateHypo caches (feed, R, J, tr, dv, grad u/p, 1/Qbar, k_sp, Q/H/S) bit-identical", tag);
  gate(okPor, "A2 [%s] hypoPorosity bit-identical", tag);

  double sMat[NG * 6], Dmat[NG * 36];
  for (int gp = 0; gp < NG; gp++)
    mat.tangent(gp, &A.hypoFeed_[gp * 6], &sMat[gp * 6], &Dmat[gp * 36]);
  double fuA[NU], fuB[NU], KsA[NU * NU], KsB[NU * NU], fsA[NP], fsB[NP];
  A.buildCoreForceHypo(sMat, f.p1, fuA);
  B.buildCoreForceHypo(sMat, f.p1, fuB);
  A.buildSolidKHypo(sMat, Dmat, f.p1, KsA);
  B.buildSolidKHypo(sMat, Dmat, f.p1, KsB);
  A.fseepHypo(fsA);
  B.fseepHypo(fsB);
  gate(sameArr(fuA, fuB, NU), "A2 [%s] buildCoreForceHypo bit-identical", tag);
  gate(sameArr(KsA, KsB, NU * NU), "A2 [%s] buildSolidKHypo bit-identical", tag);
  gate(sameArr(fsA, fsB, NP), "A2 [%s] hypo f_seep bit-identical", tag);

  // A3: the kernel's elementResidual<double> vs the class residual (velocity
  // form cancelled by sOld) — round-off; shippedTangent vs the class blocks.
  double vel[NDOF];
  for (int k = 0; k < NDOF; k++) vel[k] = 3.0 * urand();
  double rCls[NDOF], rKer[NDOF];
  A.residual(fuA, f.p1, vel, rCls);
  UH::elementResidual(f.in, f.u1, f.p1, mat, rKer);
  double eu = 0, ep = 0, su = 0, sp = 0;
  for (int i = 0; i < NU; i++) { eu = std::fmax(eu, std::fabs(rCls[i] - rKer[i])); su = std::fmax(su, std::fabs(rCls[i])); }
  for (int i = NU; i < NDOF; i++) { ep = std::fmax(ep, std::fabs(rCls[i] - rKer[i])); sp = std::fmax(sp, std::fabs(rCls[i])); }
  bool uBits = sameArr(rCls, rKer, NU);
  gate(uBits && ep <= 1e-12 * sp,
       "A3 [%s] elementResidual<double> == class residual: u-rows bit-identical %d, p-rows %.2e rel (the sOld cancellation)",
       tag, uBits ? 1 : 0, ep / sp);
  double K[NDOF * NDOF];
  const double c2 = 1.7 / f.in.dt;
  UH::shippedTangent(f.in, f.u1, f.p1, sMat, Dmat, c2, K);
  bool okK = true;
  for (int i = 0; i < NU; i++)
    for (int j = 0; j < NU; j++) okK = okK && sameBits(K[i * NDOF + j], KsA[i * NU + j]);
  for (int i = 0; i < NU; i++)
    for (int b = 0; b < NP; b++)
      okK = okK && sameBits(K[i * NDOF + NU + b], -A.QrotBlk_[i * NP + b])
                && sameBits(K[(NU + b) * NDOF + i], c2 * A.QrotBlk_[i * NP + b]);
  for (int b = 0; b < NP; b++)
    for (int c = 0; c < NP; c++)
      okK = okK && sameBits(K[(NU + b) * NDOF + NU + c],
                            A.HcurBlk_[b * NP + c] + c2 * (A.ScurBlk_[b * NP + c] + f.HtAl[b * NP + c]));
  gate(okK, "A3 [%s] shippedTangent == the class's [Ks, -Q; c2 Q', H + c2(S + aH~)] bit-identical", tag);
}

// =========================================================================== //
//  B. polar derivative                                                        //
// =========================================================================== //
static void legB()
{
  double worst = 0.0;
  int maxIt = 0;
  for (int c = 0; c < 60; c++) {
    double Rr[9], S[9], F[9];
    rotAxis(urand(), urand(), urand() + 1.5, 3.0 * urand(), Rr);
    for (int k = 0; k < 9; k++) S[k] = ((k % 4 == 0) ? 1.0 : 0.0);
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++) {
        double v = ((i == j) ? 0.6 : 0.3) * urand();
        S[3 * i + j] += v; if (i != j) S[3 * j + i] += v;
      }
    mat3mul(Rr, S, F);
    D9 FD[9], RD[9];
    for (int k = 0; k < 9; k++) FD[k] = D9::variable(F[k], k);
    if (!ladruno_hypo::polar3(FD, RD)) continue;
    double R[9], U[9];
    for (int k = 0; k < 9; k++) R[k] = RD[k].v;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
      double s = 0.0; for (int k = 0; k < 3; k++) s += R[3*k+i] * F[3*k+j]; U[3*i+j] = s; }
    const double tr = U[0] + U[4] + U[8];
    double M[9];
    for (int k = 0; k < 9; k++) M[k] = -U[k];
    M[0] += tr; M[4] += tr; M[8] += tr;
    double Mi[9];
    ladruno_hypo::inv3(M, Mi);
    double scale = 0.0, err = 0.0;
    for (int k = 0; k < 9; k++) {
      // dF = E_k: A = Rᵀ dF, a = axial(skew A), ω = M⁻¹ (2a), dR = R [ω]×
      const int I = k / 3, Jc = k % 3;
      double A[9] = {0,0,0,0,0,0,0,0,0};
      for (int i = 0; i < 3; i++) A[3 * i + Jc] = R[3 * I + i];
      double a[3] = {0.5 * (A[7] - A[5]), 0.5 * (A[2] - A[6]), 0.5 * (A[3] - A[1])};
      double w[3];
      for (int i = 0; i < 3; i++) w[i] = 2.0 * (Mi[3*i] * a[0] + Mi[3*i+1] * a[1] + Mi[3*i+2] * a[2]);
      double W[9] = {0, -w[2], w[1],  w[2], 0, -w[0],  -w[1], w[0], 0};
      double dR[9];
      mat3mul(R, W, dR);
      for (int m = 0; m < 9; m++) {
        scale = std::fmax(scale, std::fabs(dR[m]));
        err = std::fmax(err, std::fabs(dR[m] - RD[m].d[k]));
      }
    }
    worst = std::fmax(worst, err / scale);
    (void)maxIt;
  }
  gate(worst <= 1e-12, "B polar3<Dual<9>> dR == closed form R[w]x, (trU I - U) w = 2 axial(skew(R'dF)): max rel %.2e (60 F)", worst);
}

// =========================================================================== //
//  C/D/E. the element Jacobian, the dropped blocks, the mutants               //
// =========================================================================== //
// Oracle tolerance on the FD plateau (per entry, normalized by the block's max
// entry). The plateau itself sits at ~1e-11 on every state (central FD
// truncation meets round-off); 1e-8 keeps a 400x margin. A mutant must exceed
// 100·TOL.
static const double TOL = 1e-8;

struct TableRow { const char* name; double v[8]; };
static const int NSTATE_MAX = 8;
static TableRow g_rows[32];
static int g_nrows = 0;
static const char* g_stateNames[NSTATE_MAX];
static int g_nstate = 0;

static void rowSet(int r, const char* name, int s, double v)
{
  g_rows[r].name = name;
  if (r >= g_nrows) g_nrows = r + 1;
  g_rows[r].v[s] = v;
}

template <class M>
static void runState(const State& st, Fixture& f, M& mat, int sIdx)
{
  const ElemIn& in = f.in;
  printf("STATE %s\n", st.name);
  g_stateNames[sIdx] = st.name;

  // kinematic diagnostics
  double Jmin = 1e9, Jmax = -1e9, kcMin = 1e9, kcMax = -1e9;
  for (int gp = 0; gp < NG; gp++) {
    double deps[6], R1[9], J1, tr, dv, gU[NU], gP[24];
    UH::gpKinematics(NN, NP, in.dNdX + gp * 24, in.dNp + gp * 24, in.dV[gp], in.un, f.u1,
                     deps, R1, J1, tr, dv, gU, gP);
    double n = UH::porosity(J1, in.n0), kc = UH::kozenyCarmanScale(n, in.n0);
    Jmin = std::fmin(Jmin, J1); Jmax = std::fmax(Jmax, J1);
    kcMin = std::fmin(kcMin, kc); kcMax = std::fmax(kcMax, kc);
  }
  printf("INFO %s J in [%.4f, %.4f]  KC scale in [%.3e, %.3e]%s\n", st.name, Jmin, Jmax, kcMin,
         kcMax, st.kozeny ? "" : " (KC off)");

  // ---- C: dual Jacobian, value bits, FD plateau ----
  double J[NDOF * NDOF], rv[NDOF], rd[NDOF];
  UH::ElemParts<D32> parts;
  bool ok = dualJacobian(in, f.u1, f.p1, mat, J, rv, &parts);
  UH::elementResidual(in, f.u1, f.p1, mat, rd);
  gate(ok && sameArr(rv, rd, NDOF), "C [%s] Dual<32> residual value part == double residual bit for bit", st.name);

  double eFD[4], Jfd[NDOF * NDOF];
  fdCompare(in, f.u1, f.p1, mat, J, eFD, Jfd);
  gate(std::fmax(std::fmax(eFD[0], eFD[1]), std::fmax(eFD[2], eFD[3])) <= TOL,
       "C [%s] dual Jacobian == central FD (plateau): uu %.2e  up %.2e  pu %.2e  pp %.2e",
       st.name, eFD[0], eFD[1], eFD[2], eFD[3]);

  // ---- D: shipped tangent and the dropped blocks ----
  double sMat[NG * 6], Dmat[NG * 36], K[NDOF * NDOF];
  materialAtTrial(in, f.u1, mat, sMat, Dmat);
  const bool transient = in.dt > 0.0;
  const double c2 = transient ? 1.0 / in.dt : 0.0;
  UH::shippedTangent(in, f.u1, f.p1, sMat, Dmat, c2, K);

  const double nUU = blkNorm(J, 0, NDOF, 0, NU, 0, NU), nUP = blkNorm(J, 0, NDOF, 0, NU, NU, NDOF);
  const double nPU = blkNorm(J, 0, NDOF, NU, NDOF, 0, NU), nPP = blkNorm(J, 0, NDOF, NU, NDOF, NU, NDOF);
  const double nProw = std::sqrt(nPU * nPU + nPP * nPP);

  // attribution by freezing
  double Jrot[NDOF * NDOF], Jkr[NDOF * NDOF], Jpor[NDOF * NDOF], tmp[NDOF];
  UH::ElemParts<D32> pRot, pKr, pPor;
  UH::Freeze fz;
  fz.pushRotation = true;  dualJacobian(in, f.u1, f.p1, mat, Jrot, tmp, &pRot, fz);
  fz = UH::Freeze(); fz.kspRotation = true; dualJacobian(in, f.u1, f.p1, mat, Jkr, tmp, &pKr, fz);
  fz = UH::Freeze(); fz.porosity = true;    dualJacobian(in, f.u1, f.p1, mat, Jpor, tmp, &pPor, fz);

  // pieces (u-columns)
  double dfuP[NU * NU], dH[NP * NU], dHkr[NP * NU], dHpor[NP * NU], dSeep[NP * NU];
  double dS[NP * NU], dCpl[NP * NU];
  partJac(parts.fuP, NU, dfuP);
  partJac(parts.fpH, NP, dH);
  partJac(pKr.fpH, NP, dHkr);
  partJac(pPor.fpH, NP, dHpor);
  partJac(parts.fpSeep, NP, dSeep);
  partJac(parts.fpS, NP, dS);
  partJac(parts.fpCpl, NP, dCpl);

  // the shipped -alpha p geometric term alone (Dmat = 0, sMat = 0 => only -alpha p m)
  double z6[NG * 6], z36[NG * 36], Kp[NDOF * NDOF];
  for (int k = 0; k < NG * 6; k++) z6[k] = 0.0;
  for (int k = 0; k < NG * 36; k++) z36[k] = 0.0;
  UH::shippedTangent(in, f.u1, f.p1, z6, z36, 0.0, Kp);

  double sq;
  // uu
  rowSet(0, "uu  total ||J_uu - K_uu|| / ||J_uu||", sIdx, blkNorm(J, K, NDOF, 0, NU, 0, NU) / nUU);
  rowSet(1, "uu  rate-convection dR1/du in R1 s R1' (frozen R1)", sIdx, blkNorm(J, Jrot, NDOF, 0, NU, 0, NU) / nUU);
  sq = 0.0;
  for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) {
    double d = dfuP[i * NU + j] - Kp[i * NDOF + j]; sq += d * d; }
  rowSet(2, "uu  d(Q(u) p)/du beyond the -alpha p initial-stress term", sIdx, std::sqrt(sq) / nUU);
  // the rest of the uu gap: J_uu|R1 frozen − K_uu − (dfuP/du − Kp_uu): the
  // effective-stress initial-stress form (the shipped δ_ik term vs the exact
  // d(∇N)·σ + ∇N·σ·div δu) and the Hughes–Winget midpoint kinematics of Δε
  sq = 0.0;
  for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) {
    double d = Jrot[i * NDOF + j] - K[i * NDOF + j] - (dfuP[i * NU + j] - Kp[i * NDOF + j]);
    sq += d * d; }
  rowSet(3, "uu  remainder: s' initial-stress form + midpoint kinematics", sIdx, std::sqrt(sq) / nUU);
  // up
  rowSet(4, "up  ||J_up - (-Q_cur)|| / ||J_up||  (exact by design)", sIdx, blkNorm(J, K, NDOF, 0, NU, NU, NDOF) / nUP);
  // pu
  if (transient) {
    rowSet(5, "pu  total ||J_pu - Q_cur'/dt|| / ||J_pu||", sIdx, blkNorm(J, K, NDOF, NU, NDOF, 0, NU) / nPU);
    sq = 0.0;
    for (int b = 0; b < NP; b++) for (int j = 0; j < NU; j++) {
      double d = dCpl[b * NU + j] / in.dt - K[(NU + b) * NDOF + j]; sq += d * d; }
    rowSet(6, "pu  d(alpha tr(deps) dv)/du /dt - Q_cur'/dt  (the Q. coupling)", sIdx, std::sqrt(sq) / nPU);
    sq = 0.0;
    for (int k = 0; k < NP * NU; k++) sq += dS[k] * dS[k];
    rowSet(7, "pu  d(S(u)(p - p_n))/du /dt  (storage 1/Qbar(n(J)), dv, grads)", sIdx, std::sqrt(sq) / in.dt / nPU);
  } else {
    rowSet(5, "pu  total ||J_pu - Q_cur'/dt|| / ||J_pu||", sIdx, blkNorm(J, K, NDOF, NU, NDOF, 0, NU) / (nPU > 0 ? nPU : 1));
    rowSet(6, "pu  d(alpha tr(deps) dv)/du /dt - Q_cur'/dt  (the Q. coupling)", sIdx, 0.0);
    rowSet(7, "pu  d(S(u)(p - p_n))/du /dt  (storage 1/Qbar(n(J)), dv, grads)", sIdx, 0.0);
  }
  sq = 0.0;
  for (int k = 0; k < NP * NU; k++) sq += dH[k] * dH[k];
  rowSet(8, "pu  d(H(u) p)/du  (all of it is dropped)", sIdx, std::sqrt(sq) / nPU);
  sq = 0.0;
  for (int k = 0; k < NP * NU; k++) { double d = dH[k] - dHpor[k]; sq += d * d; }
  rowSet(9, "pu    of which Kozeny-Carman kc(n(J))", sIdx, std::sqrt(sq) / nPU);
  sq = 0.0;
  for (int k = 0; k < NP * NU; k++) { double d = dH[k] - dHkr[k]; sq += d * d; }
  rowSet(10, "pu    of which the rotating axes R kbar R'", sIdx, std::sqrt(sq) / nPU);
  sq = 0.0;
  for (int k = 0; k < NP * NU; k++) sq += dSeep[k] * dSeep[k];
  rowSet(11, "pu  d(f_seep(u))/du  (all of it is dropped)", sIdx, std::sqrt(sq) / nPU);
  // pp
  rowSet(12, "pp  ||J_pp - (H + (S + aH~)/dt)|| / ||J_pp||  (exact by design)", sIdx,
         blkNorm(J, K, NDOF, NU, NDOF, NU, NDOF) / nPP);
  // row blocks
  rowSet(13, "u-rows ||[J_uu J_up] - K|| / ||[J_uu J_up]||", sIdx,
         std::sqrt(std::pow(blkNorm(J, K, NDOF, 0, NU, 0, NU), 2) + std::pow(blkNorm(J, K, NDOF, 0, NU, NU, NDOF), 2))
         / std::sqrt(nUU * nUU + nUP * nUP));
  rowSet(14, "p-rows ||[J_pu J_pp] - K|| / ||[J_pu J_pp]||", sIdx,
         std::sqrt(std::pow(blkNorm(J, K, NDOF, NU, NDOF, 0, NU), 2) + std::pow(blkNorm(J, K, NDOF, NU, NDOF, NU, NDOF), 2))
         / nProw);
  if (transient) {
    // Newmark gamma = 0.6, beta = 0.3025 (the hypo battery): c2 = gamma/(beta dt)
    double KN[NDOF * NDOF];
    UH::shippedTangent(in, f.u1, f.p1, sMat, Dmat, 0.6 / (0.3025 * in.dt), KN);
    rowSet(15, "p-rows Newmark c2 = g/(b dt) vs 1/dt: ||K_N - K|| / ||[J_pu J_pp]||", sIdx,
           std::sqrt(std::pow(blkNorm(KN, K, NDOF, NU, NDOF, 0, NU), 2) + std::pow(blkNorm(KN, K, NDOF, NU, NDOF, NU, NDOF), 2))
           / nProw);
  } else {
    rowSet(15, "p-rows Newmark c2 = g/(b dt) vs 1/dt: ||K_N - K|| / ||[J_pu J_pp]||", sIdx, 0.0);
  }
  gate(blkNorm(J, K, NDOF, 0, NU, NU, NDOF) <= 1e-13 * nUP && blkNorm(J, K, NDOF, NU, NDOF, NU, NDOF) <= 1e-13 * nPP,
       "D [%s] the two blocks the shipped tangent claims exact are exact: up %.2e, pp %.2e",
       st.name, blkNorm(J, K, NDOF, 0, NU, NU, NDOF) / nUP, blkNorm(J, K, NDOF, NU, NDOF, NU, NDOF) / nPP);

  // ---- E: break-the-gate legs (judged against the FD columns) ----
  double eK[4], eR[4], eP[4];
  blkErr(K, Jfd, eK);
  blkErr(Jrot, Jfd, eR);
  blkErr(Jpor, Jfd, eP);
  printf("MUTANT %s shipped(uu %.2e pu %.2e)  noPushRotation(uu %.2e)  noPorosity(pu %.2e)\n",
         st.name, eK[0], eK[2], eR[0], eP[2]);
  if (sIdx == 0 || std::strcmp(st.name, "rotated") == 0 || std::strcmp(st.name, "j2_plastic") == 0)
    gate(std::fmax(eK[0], eK[2]) > 100.0 * TOL, "E [%s] the shipped tangent FAILS the oracle (%.2e > 100 tol)",
         st.name, std::fmax(eK[0], eK[2]));
  if (std::strcmp(st.name, "rotated") == 0)
    gate(eR[0] > 100.0 * TOL, "E [%s] a Jacobian without dR1/du in the stress push FAILS the oracle (uu %.2e)",
         st.name, eR[0]);
  if (std::strcmp(st.name, "compacted_kc") == 0 || std::strcmp(st.name, "fast_drain") == 0)
    gate(eP[2] > 100.0 * TOL, "E [%s] a Jacobian without the porosity terms FAILS the oracle (pu %.2e)",
         st.name, eP[2]);
  if (std::strcmp(st.name, "kc_floor") == 0) {
    // at the floor dkc/dn = 0: freezing porosity changes only the storage path
    double sqd = 0.0;
    for (int k = 0; k < NP * NU; k++) { double d = dH[k] - dHpor[k]; sqd += d * d; }
    gate(std::sqrt(sqd) <= 1e-14 * (nPU > 0 ? nPU : 1),
         "C [%s] at the Kozeny-Carman floor the dual KC slope is exactly 0 (one-sided, %.2e)",
         st.name, std::sqrt(sqd) / nPU);
  }
}

// chainStress seam: dual through an NDMaterial-style (stress, tangent) pair
template <class M>
static void legChain(const char* tag, Fixture& f, M& mat)
{
  ChainMat<M> cm; cm.m = &mat;
  double Jn[NDOF * NDOF], Jc[NDOF * NDOF], r1[NDOF], r2[NDOF];
  dualJacobian(f.in, f.u1, f.p1, mat, Jn, r1);
  dualJacobian(f.in, f.u1, f.p1, cm, Jc, r2);
  double e[4];
  blkErr(Jc, Jn, e);
  gate(std::fmax(std::fmax(e[0], e[1]), std::fmax(e[2], e[3])) <= 1e-11,
       "C [%s] chainStress (stress + consistent tangent) == native dual material: max %.2e",
       tag, std::fmax(std::fmax(e[0], e[1]), std::fmax(e[2], e[3])));
}

// E: the incremental coupling on a rigid increment (ADR 78)
static void legRigid()
{
  Fixture f;
  State st = {"rigid", {0}, {0}, {0.1, -0.2, 0.05}, {0}, 0.0, 0.0, 20.0, 4.0, 2.0, 0.0,
              {1e-4, 5e-5, 1e-5}, true, 0.05, 0};
  double Rn[9], R1[9], S[9] = {1.04, 0.02, 0.0, 0.02, 0.97, 0.01, 0.0, 0.01, 1.01};
  rotAxis(1, 2, 3, 0.4, Rn);
  rotAxis(1, 2, 3, 0.4 + 0.35, R1);   // a 20 degree rigid increment on a stretched state
  mat3mul(Rn, S, st.An);
  mat3mul(R1, S, st.A1);
  for (int i = 0; i < 3; i++) st.c1[i] = st.cn[i] + 0.3;
  setup(f, st);
  ElasticMat em; const double pre0[6] = {0, 0, 0, 0, 0, 0}; em.init(2.0e4, 0.3, pre0);
  UH::ElemParts<double> parts;
  double r[NDOF];
  UH::elementResidual(f.in, f.u1, f.p1, em, r, &parts);
  double cpl = 0.0, vel = 0.0, scale = 0.0;
  NewElem B; B.in = &f.in; B.updateHypo(f.u1);
  for (int b = 0; b < NP; b++) {
    cpl = std::fmax(cpl, std::fabs(parts.fpCpl[b]));
    double s = 0.0;
    for (int k = 0; k < NU; k++) s += B.QrotBlk_[k * NP + b] * (f.u1[k] - f.un[k]);
    vel = std::fmax(vel, std::fabs(s));
    for (int k = 0; k < NU; k++) scale = std::fmax(scale, std::fabs(B.QrotBlk_[k * NP + b]));
  }
  gate(cpl <= 1e-14 && vel > 1e-3,
       "E rigid 20-degree increment: incremental coupling sum|alpha tr(deps) dv| = %.2e (exactly zero), velocity form |Q_cur' du| = %.2e",
       cpl, vel);
}

// F: cost
template <class M>
static void legCost(Fixture& f, M& mat)
{
  const int N = 200;
  double r[NDOF], J[NDOF * NDOF], rv[NDOF], sink = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < N; k++) { UH::elementResidual(f.in, f.u1, f.p1, mat, r); sink += r[k % NDOF]; }
  auto t1 = std::chrono::steady_clock::now();
  for (int k = 0; k < N; k++) { dualJacobian(f.in, f.u1, f.p1, mat, J, rv); sink += J[k % NDOF]; }
  auto t2 = std::chrono::steady_clock::now();
  double sMat[NG * 6], Dmat[NG * 36], K[NDOF * NDOF];
  materialAtTrial(f.in, f.u1, mat, sMat, Dmat);
  for (int k = 0; k < N; k++) { UH::shippedTangent(f.in, f.u1, f.p1, sMat, Dmat, 20.0, K); sink += K[k % NDOF]; }
  auto t3 = std::chrono::steady_clock::now();
  double a = std::chrono::duration<double, std::micro>(t1 - t0).count() / N;
  double b = std::chrono::duration<double, std::micro>(t2 - t1).count() / N;
  double c = std::chrono::duration<double, std::micro>(t3 - t2).count() / N;
  printf("COST us/element: residual<double> %.2f  shippedTangent %.2f  residual<Dual<32>> (full Jacobian) %.2f  "
         "ratio dual/(residual+tangent) %.1f  [g++ -O2, elastic, sink %.3g]\n",
         a, c, b, b / (a + c), sink);
}

int main()
{
  printf("ladruno_up_hypo_dual_check (WP-266)\n");
  legA_kernel();
  legB();

  // ---- the states ----
  const double I9[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  const int NST = 7;
  State S[NST];
  // 1. large strain: stretch + shear, a 4% increment with a quadratic bump
  S[0] = {"large_strain", {1.12, 0.15, 0.0, 0.0, 0.93, 0.05, 0.02, 0.0, 0.85},
          {1.15, 0.17, 0.01, -0.01, 0.91, 0.06, 0.03, 0.02, 0.81}, {0, 0, 0}, {0.01, -0.02, 0.0},
          0.01, 0.03, 40.0, 20.0, 6.0, 0.02, {1e-4, 5e-5, 1e-5}, true, 0.05, 0};
  // 2. compacted near the Kozeny-Carman floor (J ~ 0.615, kc ~ 1e-4)
  S[1] = {"compacted_kc", {1.0, 0.01, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.63},
          {1.0, 0.012, 0.0, 0.0, 1.0, 0.003, 0.0, 0.0, 0.616}, {0, 0, 0}, {0, 0, 0},
          0.0, 0.001, 80.0, 10.0, 8.0, 0.01, {1e-3, 1e-3, 2e-4}, true, 0.05, 0};
  // 3. at the floor (J ~ 0.6025: kc = 1e-6 clamped)
  S[2] = {"kc_floor", {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.61},
          {1.0, 0.004, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.6025}, {0, 0, 0}, {0, 0, 0},
          0.0, 0.0002, 80.0, 10.0, 8.0, 0.01, {1e-3, 1e-3, 2e-4}, true, 0.05, 0};
  // 4. rotated: 50 deg + 12 deg increment about an oblique axis, moderate stretch
  {
    double Rn[9], R1[9], U[9] = {1.05, 0.03, 0.0, 0.03, 0.97, 0.02, 0.0, 0.02, 0.99};
    double U1[9] = {1.06, 0.035, 0.0, 0.035, 0.96, 0.02, 0.0, 0.02, 0.985};
    rotAxis(1.0, 2.0, 3.0, 50.0 * PI / 180.0, Rn);
    rotAxis(1.0, 2.0, 3.0, 62.0 * PI / 180.0, R1);
    S[3] = {"rotated", {0}, {0}, {0.2, 0.1, -0.1}, {0.25, 0.12, -0.08}, 0.0, 0.0,
            30.0, 15.0, 5.0, 0.01, {1e-4, 2e-5, 1e-6}, true, 0.05, 0};
    mat3mul(Rn, U, S[3].An);
    mat3mul(R1, U1, S[3].A1);
  }
  // 5. small strain: everything near the reference configuration
  S[4] = {"small_strain", {0}, {0}, {0, 0, 0}, {0, 0, 0}, 0.0, 0.0, 10.0, 2.0, 1.0, 1e-4,
          {1e-4, 5e-5, 1e-5}, true, 0.05, 0};
  for (int k = 0; k < 9; k++) {
    S[4].An[k] = I9[k] + 1e-4 * ((k % 3) - 1.0);
    S[4].A1[k] = S[4].An[k] + 2e-5 * (1.0 + 0.3 * k) * ((k % 2) ? 1.0 : -1.0);
  }
  // 6. J2 plastic at large strain (the large-strain kinematics, a yielding skeleton)
  S[5] = S[0];
  S[5].name = "j2_plastic";
  S[5].material = 1;
  // 7. fast-draining large strain: the large-strain kinematics with k x 100 and
  //    dt x 10 (c_v dt / h^2 ~ 10^3 larger): the regime where the Darcy terms
  //    H(u)p and f_seep(u) carry weight in the p-rows
  S[6] = S[0];
  S[6].name = "fast_drain";
  S[6].kbar[0] = 1e-2; S[6].kbar[1] = 5e-3; S[6].kbar[2] = 1e-3;
  S[6].dt = 0.5;

  const double prestress[6] = {-150.0, -120.0, -200.0, 15.0, -8.0, 5.0};
  ElasticMat em; em.init(2.0e4, 0.3, prestress);
  J2Mat j2;
  j2.p = ladruno_j2_kernel::Params();
  j2.p.K = 1.6e4; j2.p.G = 7.5e3; j2.p.sig0 = 60.0; j2.p.Qinf = 20.0; j2.p.bIso = 30.0;
  j2.p.Hiso = 300.0; j2.p.nBack = 1; j2.p.C[0] = 2000.0; j2.p.gam[0] = 40.0;

  for (int s = 0; s < NST; s++) {
    Fixture f;
    setup(f, S[s]);
    if (S[s].material == 0) {
      if (s == 0) legA_class(S[s].name, f, em);
      if (s == 3) legA_class(S[s].name, f, em);
      runState(S[s], f, em, s);
      if (s == 0) legChain(S[s].name, f, em);
    } else {
      // confirm the skeleton yields at every GP
      double sMat[NG * 6], Dmat[NG * 36];
      int nPl = 0;
      for (int gp = 0; gp < NG; gp++) {
        double deps[6], R1[9], J1, tr, dv, gU[NU], gP[24], feed[6], dG;
        UH::gpKinematics(NN, NP, f.in.dNdX + gp * 24, f.in.dNp + gp * 24, f.in.dV[gp], f.in.un,
                         f.u1, deps, R1, J1, tr, dv, gU, gP);
        for (int c = 0; c < 6; c++) feed[c] = f.in.feedCommit[gp * 6 + c] + deps[c];
        j2.tangent(gp, feed, &sMat[gp * 6], &Dmat[gp * 36], &dG);
        if (dG > 0.0) nPl++;
      }
      gate(nPl == NG, "C [%s] the skeleton is plastic at %d/8 GPs", S[s].name, nPl);
      legA_class(S[s].name, f, j2);
      runState(S[s], f, j2, s);
      legChain(S[s].name, f, j2);
    }
    g_nstate = s + 1;
  }
  legRigid();
  {
    Fixture f; setup(f, S[0]);
    legCost(f, em);
  }

  // ---- the dropped-block table ----
  printf("\nTABLE dropped blocks: ||part|| / ||block of the full dual Jacobian|| (Frobenius); c2 = 1/dt (dt = 0.05; fast_drain 0.5)\n");
  printf("TABLE %-72s", "term");
  for (int s = 0; s < g_nstate; s++) printf(" | %-12s", g_stateNames[s]);
  printf("\n");
  for (int r = 0; r < g_nrows; r++) {
    printf("TABLE %-72s", g_rows[r].name);
    for (int s = 0; s < g_nstate; s++) printf(" | %12.3e", g_rows[r].v[s]);
    printf("\n");
  }

  printf("\n[checks] %d gates\n", g_checks);
  if (g_fail == 0) printf("RESULT: PASS\n");
  else printf("RESULT: FAIL (%d)\n", g_fail);
  return g_fail == 0 ? 0 : 1;
}
