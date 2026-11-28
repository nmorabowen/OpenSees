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
// Created: 11/2026  (WP-266, dual-number program #1174)
//
// LadrunoUPHypoKernel — the per-Gauss-point residual of the LadrunoUP
// `-geom hypo` lane (ADR 79 P2), written ONCE as `template <class T>`.
//
// LadrunoUP::updateHypo / buildCoreForceHypo / buildSolidKHypo /
// getResistingForce call the T = double instantiations below; their operations
// and their order are the pre-WP-266 class code's, so the class's double path is
// bit-identical (tests/ladruno_up_hypo_dual_check.cpp compares every per-GP
// function bit for bit against a verbatim copy of the pre-WP code).
//
// T = ladruno_dual::Dual<32> seeds the 24 nodal displacements and the 8 nodal
// pressures of the H8 equal-order element. elementResidual<Dual<32>> then
// returns the COMPLETE element Jacobian in one pass: every block the shipped
// "residual exact, tangent dominant-term" tangent drops (shippedTangent below)
// is in it. Today that is an oracle (WP-266); the opt-in production flag is
// wave-2 P4 of the plan (Ladruno_implementation/_wp262_dual_orchestration_plan.md).
//
// THE RESIDUAL that elementResidual returns (block order [u (3·nNodes) ; p (nP)]):
//
//   r_u = Σ_gp ∫ ∇ₓNᵀ (R₁ σ̂(ε_feed) R₁ᵀ − α p m) dv               (u-rows)
//   r_p = H(u) p − f_seep(u)                                        (p-rows)
//       + [ Σ_gp N_p α tr(Δε) dv + (S(u) + αH̃)(p − p_n) ] / Δt      (Δt > 0)
//
// i.e. the INCREMENTAL storage-coupling form of ADR 78 §3.3 / ADR 79 P2:
// α·tr(Δε)·dv is exactly zero on a rigid increment. Its derivative is the
// "Q·" sensitivity. Differentiating the velocity form Q_curᵀ·u̇ instead would
// bring back the rigid-rotation spurious pressure ADR 78 removed, so this kernel
// never offers it. LadrunoUP::getResistingForceIncInertia adds the same
// rate terms after cancelling the integrator's damp·v p-rows. The class result
// therefore equals this one to round-off, not bit for bit; the check measures
// the gap. Left out, all u-independent or exact in the shipped effective tangent:
// the dead body load, nodal loads, the integrator's M·ü and C_Ray·u̇, and
// -dynSeepage's −ü in the seepage drive.
//
// The material enters through a functor  mat(gp, const T feed[6], T sigHat[6])
// (unrotated-frame stress, Voigt {00,11,22,01,12,20}, engineering-shear strain).
// A T-templated material (a kernel such as LadrunoJ2Kernel.h's returnMapT) runs
// natively in Dual. Any NDMaterial can join through chainStress(): the value is
// its stress, the derivative D·dε_feed from its consistent tangent.
//
// Kinks keep their one-sided handling: the porosity clamp n ∈ [1e-6, 1−1e-6]
// and the Kozeny–Carman floor kc ≥ 1e-6 return zero slope on the clamped side
// (correct for the floored function). The polar rotations go through
// ladruno_hypo::polar3 (Higham iteration, differentiated as executed; the
// check proves dR against the closed-form polar derivative).

#ifndef LadrunoUPHypoKernel_h
#define LadrunoUPHypoKernel_h

#include <math.h>
#include <LadrunoDual.h>        // WP-262: forward-mode dual numbers
#include <LadrunoUPKernel.h>    // oneOverQbar (templated on n, WP-266)
#include <LadrunoHypoKernel.h>  // hypoIncrement, inv3, polar3, pushStress6 (templated, WP-266)

namespace ladruno_up_hypo {

using ladruno_dual::value;

// Hard caps (the class's hypo lane is H8 only; sized for its UP_MAX_NODES = 10).
const int MAXN  = 10;          // nodes
const int MAXU  = 3 * MAXN;    // solid DOFs
const int MAXP  = 10;          // pressure carriers
const int MAXGP = 8;           // Gauss points

// frozen(x): the value of x with every derivative dropped (identity on double).
// The oracle uses it to switch off one dependency and attribute a block.
template <class T>
inline T frozen(const T& x) { return T(value(x)); }

// --------------------------------------------------------------------------- //
//  Per-GP pieces — the class calls each of these with T = double              //
// --------------------------------------------------------------------------- //

// Kinematics at one GP: the Hughes–Winget objective increment (de-rotated into
// the unrotated frame), R₁ = polar(F₁), J₁ = det F₁, tr(Δε), dv = J₁·dV, and the
// spatial solid / pressure gradients through F₁⁻¹. Returns false on a
// collapsed / inverted trial or midpoint configuration (gU, gP untouched).
//   dNdX[a*3+K], dNpRef[b*3+K]: reference gradients; un: committed nodal disp;
//   u1: trial nodal disp; gU[a*3+j], gP[b*3+j]: spatial gradients (out).
template <class T>
inline bool gpKinematics(int nNodes, int nP, const double *dNdX,
                         const double *dNpRef, double dV,
                         const double *un, const T *u1,
                         T deps[6], T R1[9], T &J1, T &trDeps, T &dvCur,
                         T *gU, T *gP)
{
  if (!ladruno_hypo::hypoIncrement(nNodes, dNdX, un, u1, deps, R1, J1))
    return false;
  trDeps = deps[0] + deps[1] + deps[2];     // tr is frame-invariant
  dvCur = J1 * dV;

  T H1[9], F1[9], Fi[9];
  ladruno_hypo::dispGradient(nNodes, dNdX, u1, H1);
  for (int c = 0; c < 9; c++) F1[c] = H1[c];
  F1[0] += 1.0; F1[4] += 1.0; F1[8] += 1.0;
  ladruno_hypo::inv3(F1, Fi);
  for (int a = 0; a < nNodes; a++)
    for (int j = 0; j < 3; j++) {
      T s = 0.0;
      for (int K = 0; K < 3; K++)
        s += dNdX[a * 3 + K] * Fi[3 * K + j];
      gU[a * 3 + j] = s;
    }
  for (int bb = 0; bb < nP; bb++)
    for (int j = 0; j < 3; j++) {
      T s = 0.0;
      for (int K = 0; K < 3; K++)
        s += dNpRef[bb * 3 + K] * Fi[3 * K + j];
      gP[bb * 3 + j] = s;
    }
  return true;
}

// n(J) = 1 − (1−n₀)/J — the kinematic porosity update, clamped away from 0/1
// so 1/Q̄ and Kozeny–Carman stay finite (zero slope on a clamped side).
template <class T>
inline T porosity(const T &J, double n0)
{
  T n = 1.0 - (1.0 - n0) / J;
  if (n < 1.0e-6)        n = 1.0e-6;
  else if (n > 1.0 - 1.0e-6) n = 1.0 - 1.0e-6;
  return n;
}

// Kozeny–Carman scale normalized to 1 at n = n₀, floored at 1e-6 so a heavily
// compacted GP degrades to near-impermeable without singularizing H.
template <class T>
inline T kozenyCarmanScale(const T &n, double n0)
{
  T kc = (n * n * n / ((1.0 - n) * (1.0 - n)))
       * ((1.0 - n0) * (1.0 - n0) / (n0 * n0 * n0));
  if (kc < 1.0e-6) kc = 1.0e-6;
  return kc;
}

// Porosity-evolved storage coefficient 1/Q̄(n(J)) and the spatial permeability
// k_sp = kc(n)·R·diag(k̄)·Rᵀ (axes skeleton-attached, ADR 78 §3.4).
template <class T>
inline void gpPoroPerm(const T &J1, const T R1[9], double n0, double Kf,
                       double biotAlpha, double Ks, const double kbar[3],
                       bool kozenyCarman, T &invQbar, T ksp[9])
{
  T n = porosity(J1, n0);
  invQbar = (Kf > 0.0) ? ladruno_up::oneOverQbar(n, Kf, biotAlpha, Ks) : T(0.0);
  T kc = 1.0;
  if (kozenyCarman)
    kc = kozenyCarmanScale(n, n0);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      T s = 0.0;
      for (int m = 0; m < 3; m++)
        s += R1[3 * i + m] * kbar[m] * R1[3 * j + m];
      ksp[3 * i + j] = kc * s;
    }
}

// Current-config coupling block, ACCUMULATED: Q[(a·3+i)·nP + b] += ∇ₓN_a,i α N_b dv.
template <class T>
inline void addQcur(int nNodes, int nP, const T *gU, const double *Np,
                    double biotAlpha, const T &dvCur, T *Q)
{
  for (int a = 0; a < nNodes; a++)
    for (int i = 0; i < 3; i++) {
      int row = a * 3 + i;
      T g = gU[row] * biotAlpha * dvCur;
      for (int bb = 0; bb < nP; bb++)
        Q[row * nP + bb] += g * Np[bb];
    }
}

// Current-config Darcy and storage blocks, ACCUMULATED (row stride nP):
//   H[b,c] += ∇ₓN_bᵀ k_sp ∇ₓN_c dv ,  S[b,c] += N_b (1/Q̄) N_c dv.
template <class T>
inline void addHScur(int nP, const T *gP, const T ksp[9], const double *Np,
                     const T &invQbar, const T &dvCur, T *H, T *S)
{
  for (int bb = 0; bb < nP; bb++)
    for (int cc = 0; cc < nP; cc++) {
      T s = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          s += gP[bb * 3 + i] * ksp[3 * i + j] * gP[cc * 3 + j];
      H[bb * nP + cc] += s * dvCur;
      S[bb * nP + cc] += Np[bb] * invQbar * Np[cc] * dvCur;
    }
}

// Seepage source on the current configuration, ACCUMULATED:
//   fseep[b] += ∇ₓN_bᵀ k_sp ρ_f drive dv   (drive: global frame, built by CALLER).
template <class T>
inline void addFseepCur(int nP, const T *gP, const T ksp[9], double rhoF,
                        const double drive[3], const T &dvCur, T *fseep)
{
  for (int bb = 0; bb < nP; bb++) {
    T s = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        s += gP[bb * 3 + i] * ksp[3 * i + j] * rhoF * drive[j];
    fseep[bb] += s * dvCur;
  }
}

// p at the GP: Σ_b N_b p_b.
template <class T>
inline T gpPressure(int nP, const double *Np, const T *pNodal)
{
  T p = 0.0;
  for (int bb = 0; bb < nP; bb++)
    p += Np[bb] * pNodal[bb];
  return p;
}

// u-row internal force at one GP, ACCUMULATED: fu[a·3+i] += Σ_j σ_ij ∇ₓN_a,j dv
// with the TOTAL Cauchy stress σ = R₁ σ̂ R₁ᵀ − α p m (−α p m against the same
// spatial gradients is exactly −Q_cur·p: K/Q consistency is structural).
template <class T>
inline void addForceUHypo(int nNodes, const T R1[9], const T sMat[6], const T &p,
                          const T *gU, const T &dvCur, double biotAlpha, T *fu)
{
  T sSp[6];
  ladruno_hypo::pushStress6(R1, sMat, sSp);
  T sig[3][3];
  sig[0][0] = sSp[0] - biotAlpha * p;
  sig[1][1] = sSp[1] - biotAlpha * p;
  sig[2][2] = sSp[2] - biotAlpha * p;
  sig[0][1] = sig[1][0] = sSp[3];
  sig[1][2] = sig[2][1] = sSp[4];
  sig[2][0] = sig[0][2] = sSp[5];
  for (int a = 0; a < nNodes; a++)
    for (int i = 0; i < 3; i++) {
      T s = 0.0;
      for (int j = 0; j < 3; j++)
        s += sig[i][j] * gU[a * 3 + j];
      fu[a * 3 + i] += s * dvCur;
    }
}

// The SHIPPED hypo solid tangent at one GP (double only), ACCUMULATED into the
// row-major (3·nNodes)² block Ks (row stride ldK): ∫ Bᵀ c B dv with c = M D Mᵀ
// the pushed-forward material tangent, plus the symmetric initial-stress term
// δ_ik ∫ ∇ₓN_a,j σ_jl ∇ₓN_b,l dv from the TOTAL stress σ̂_sp − α p m. The
// rate-convection terms are omitted (the "residual exact, tangent dominant-
// term" policy). This is LadrunoUP::buildSolidKHypo's per-GP body verbatim.
inline void addSolidKHypoGP(int nNodes, const double R1[9], const double Dmat[36],
                            const double sMat[6], double p, double biotAlpha,
                            const double *gU, double dvCur, double *Ks, int ldK)
{
  double Dcur[36];
  ladruno_hypo::pushTangent6(R1, Dmat, Dcur);

  double Ba[6][3], Bb[6][3], DB[6][3];
  for (int a = 0; a < nNodes; a++) {
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 3; c++) Ba[r][c] = 0.0;
    Ba[0][0] = gU[a*3+0]; Ba[1][1] = gU[a*3+1]; Ba[2][2] = gU[a*3+2];
    Ba[3][0] = gU[a*3+1]; Ba[3][1] = gU[a*3+0];
    Ba[4][1] = gU[a*3+2]; Ba[4][2] = gU[a*3+1];
    Ba[5][0] = gU[a*3+2]; Ba[5][2] = gU[a*3+0];

    for (int bn = 0; bn < nNodes; bn++) {
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
          Ks[(3 * a + i) * ldK + 3 * bn + k] += s * dvCur;
        }
    }
  }

  // geometric / initial-stress term from the TOTAL stress
  double sSp[6];
  ladruno_hypo::pushStress6(R1, sMat, sSp);
  double sig[3][3];
  sig[0][0] = sSp[0] - biotAlpha * p;
  sig[1][1] = sSp[1] - biotAlpha * p;
  sig[2][2] = sSp[2] - biotAlpha * p;
  sig[0][1] = sig[1][0] = sSp[3];
  sig[1][2] = sig[2][1] = sSp[4];
  sig[2][0] = sig[0][2] = sSp[5];

  for (int a = 0; a < nNodes; a++)
    for (int bn = 0; bn < nNodes; bn++) {
      double s = 0.0;
      for (int j = 0; j < 3; j++)
        for (int l = 0; l < 3; l++)
          s += gU[a * 3 + j] * sig[j][l] * gU[bn * 3 + l];
      s *= dvCur;
      for (int i = 0; i < 3; i++)
        Ks[(3 * a + i) * ldK + 3 * bn + i] += s;
    }
}

// chainStress: an NDMaterial's response as a T-valued stress. The value is its
// stress sig; the derivative is D·dε_feed from its 6×6 consistent tangent D
// (engineering-shear strain → stress, row-major). For T = Dual only (on double
// it returns sig + 0).
template <class T>
inline void chainStress(const double sig[6], const double D[36], const T feed[6],
                        T out[6])
{
  T df[6];
  for (int j = 0; j < 6; j++) {
    df[j] = feed[j];
    ladruno_dual::zero_value(df[j]);
  }
  for (int i = 0; i < 6; i++) {
    out[i] = T(sig[i]);
    for (int j = 0; j < 6; j++)
      out[i] += D[6 * i + j] * df[j];
  }
}

// --------------------------------------------------------------------------- //
//  Element level — the oracle's residual and the shipped tangent              //
// --------------------------------------------------------------------------- //

// Fixed (double) element data. Per-GP arrays are GP-major.
struct ElemIn {
  int nNodes, nP, nGP;
  const double *dNdX;        // nGP·nNodes·3  reference solid gradients
  const double *Np;          // nGP·nP        pressure shape values
  const double *dNp;         // nGP·nP·3      reference pressure gradients
  const double *dV;          // nGP           w·detJ (reference)
  const double *un;          // nNodes·3      committed nodal displacements
  const double *pn;          // nP            committed nodal pressures
  const double *feedCommit;  // nGP·6         committed accumulated feed strain
  const double *HtAl;        // nP·nP         α_stab·H̃ (null: none)
  double biotAlpha, Kf, Ks, n0;
  double kbar[3];
  bool   kozenyCarman;
  double rhoF;
  double drive[3];           // seepage drive (fluid body), global frame
  double dt;                 // > 0: transient incremental p-row; <= 0: static
};

// Dependencies the oracle can switch off one at a time (each removes exactly
// its own term of the chain rule, so the attributions add up).
struct Freeze {
  bool pushRotation;   // R₁ in σ = R₁ σ̂ R₁ᵀ   (the rate-convection / spin term)
  bool kspRotation;    // R₁ in k_sp = R k̄ Rᵀ   (rotating permeability axes)
  bool porosity;       // J in n(J) → 1/Q̄, kc   (storage + Kozeny–Carman)
  Freeze() : pushRotation(false), kspRotation(false), porosity(false) {}
};

// Additive pieces of the residual, for block attribution.
template <class T>
struct ElemParts {
  T fuP[MAXU];     // −Q_cur(u)·p   (the pressure part of r_u)
  T fpH[MAXP];     // H(u)·p
  T fpSeep[MAXP];  // f_seep(u)
  T fpCpl[MAXP];   // Σ_gp N_p α tr(Δε) dv        (÷Δt in r_p)
  T fpS[MAXP];     // S(u)·(p − p_n)              (÷Δt in r_p)
  T fpHt[MAXP];    // αH̃·(p − p_n)                (÷Δt in r_p)
};

// The element residual (see the file header) at trial (u1, p1). r has
// 3·nNodes + nP entries in block order. Returns false on a collapsed /
// inverted configuration.
template <class T, class Material>
inline bool elementResidual(const ElemIn &in, const T *u1, const T *p1,
                            Material &mat, T *r, ElemParts<T> *parts = 0,
                            const Freeze &fz = Freeze())
{
  const int nN = in.nNodes, nP = in.nP, nU = 3 * in.nNodes;
  T Q[MAXU * MAXP], H[MAXP * MAXP], S[MAXP * MAXP], fseep[MAXP], fu[MAXU];
  T trDeps[MAXGP], dvCur[MAXGP];
  for (int k = 0; k < nU * nP; k++) Q[k] = 0.0;
  for (int k = 0; k < nP * nP; k++) { H[k] = 0.0; S[k] = 0.0; }
  for (int k = 0; k < nP; k++) fseep[k] = 0.0;
  for (int k = 0; k < nU; k++) fu[k] = 0.0;
  if (parts)
    for (int k = 0; k < nU; k++) parts->fuP[k] = 0.0;

  for (int gp = 0; gp < in.nGP; gp++) {
    const double *dNdX = in.dNdX + gp * nN * 3;
    const double *Np   = in.Np + gp * nP;
    const double *dNp  = in.dNp + gp * nP * 3;
    T deps[6], R1[9], J1, gU[MAXU], gP[MAXP * 3];
    if (!gpKinematics(nN, nP, dNdX, dNp, in.dV[gp], in.un, u1, deps, R1, J1,
                      trDeps[gp], dvCur[gp], gU, gP))
      return false;

    T feed[6], sMat[6];
    for (int c = 0; c < 6; c++)
      feed[c] = in.feedCommit[gp * 6 + c] + deps[c];
    mat(gp, feed, sMat);

    T Rk[9], Rs[9];
    for (int c = 0; c < 9; c++) {
      Rk[c] = fz.kspRotation  ? frozen(R1[c]) : R1[c];
      Rs[c] = fz.pushRotation ? frozen(R1[c]) : R1[c];
    }
    T Jp = fz.porosity ? frozen(J1) : J1;
    T invQbar, ksp[9];
    gpPoroPerm(Jp, Rk, in.n0, in.Kf, in.biotAlpha, in.Ks, in.kbar,
               in.kozenyCarman, invQbar, ksp);
    addQcur(nN, nP, gU, Np, in.biotAlpha, dvCur[gp], Q);
    addHScur(nP, gP, ksp, Np, invQbar, dvCur[gp], H, S);

    T pg = gpPressure(nP, Np, p1);
    addForceUHypo(nN, Rs, sMat, pg, gU, dvCur[gp], in.biotAlpha, fu);
    if (parts) {
      T zero6[6];
      for (int c = 0; c < 6; c++) zero6[c] = 0.0;
      addForceUHypo(nN, Rs, zero6, pg, gU, dvCur[gp], in.biotAlpha, parts->fuP);
    }
    addFseepCur(nP, gP, ksp, in.rhoF, in.drive, dvCur[gp], fseep);
  }

  for (int k = 0; k < nU; k++)
    r[k] = fu[k];

  for (int bb = 0; bb < nP; bb++) {
    T s = 0.0;
    for (int cc = 0; cc < nP; cc++)
      s += H[bb * nP + cc] * p1[cc];
    T rp = s;
    rp -= fseep[bb];

    T cpl = 0.0, sto = 0.0, stab = 0.0;
    if (in.dt > 0.0) {
      T sNew = 0.0;
      for (int gp = 0; gp < in.nGP; gp++)                         // α·Δln J
        sNew += in.Np[gp * nP + bb] * in.biotAlpha * trDeps[gp] * dvCur[gp];
      cpl = sNew;
      for (int cc = 0; cc < nP; cc++) {
        double ht = in.HtAl ? in.HtAl[bb * nP + cc] : 0.0;
        T spc = S[bb * nP + cc] + ht;
        T dp = p1[cc] - in.pn[cc];
        sNew += spc * dp;
        sto += S[bb * nP + cc] * dp;
        stab += ht * dp;
      }
      rp += sNew / in.dt;
    }
    r[nU + bb] = rp;

    if (parts) {
      parts->fpH[bb] = s;
      parts->fpSeep[bb] = fseep[bb];
      parts->fpCpl[bb] = cpl;
      parts->fpS[bb] = sto;
      parts->fpHt[bb] = stab;
    }
  }
  return true;
}

// The SHIPPED element tangent (double), block order, row-major n×n with
// n = 3·nNodes + nP: the class's getTangentStiff + c2·getDamp restricted to
// the terms elementResidual carries (C_Ray and M are left out):
//   [ Ks_hypo          , −Q_cur              ]
//   [ c2·Q_curᵀ        , H_cur + c2(S_cur + αH̃) ]
// sMat / Dmat: the material stress and 6×6 tangent at the trial feed, per GP
// (GP-major, 6 and 36 entries). c2 is the integrator's damping factor
// (γ/(βΔt) for Newmark; the incremental residual differentiates to 1/Δt).
inline bool shippedTangent(const ElemIn &in, const double *u1, const double *p1,
                           const double *sMat, const double *Dmat, double c2,
                           double *K)
{
  const int nN = in.nNodes, nP = in.nP, nU = 3 * in.nNodes, n = nU + nP;
  double Q[MAXU * MAXP], H[MAXP * MAXP], S[MAXP * MAXP], Ks[MAXU * MAXU];
  for (int k = 0; k < nU * nP; k++) Q[k] = 0.0;
  for (int k = 0; k < nP * nP; k++) { H[k] = 0.0; S[k] = 0.0; }
  for (int k = 0; k < nU * nU; k++) Ks[k] = 0.0;
  for (int k = 0; k < n * n; k++) K[k] = 0.0;

  for (int gp = 0; gp < in.nGP; gp++) {
    const double *dNdX = in.dNdX + gp * nN * 3;
    const double *Np   = in.Np + gp * nP;
    const double *dNp  = in.dNp + gp * nP * 3;
    double deps[6], R1[9], J1, trDeps, dvCur, gU[MAXU], gP[MAXP * 3];
    if (!gpKinematics(nN, nP, dNdX, dNp, in.dV[gp], in.un, u1, deps, R1, J1,
                      trDeps, dvCur, gU, gP))
      return false;
    double invQbar, ksp[9];
    gpPoroPerm(J1, R1, in.n0, in.Kf, in.biotAlpha, in.Ks, in.kbar,
               in.kozenyCarman, invQbar, ksp);
    addQcur(nN, nP, gU, Np, in.biotAlpha, dvCur, Q);
    addHScur(nP, gP, ksp, Np, invQbar, dvCur, H, S);
    double pg = gpPressure(nP, Np, p1);
    addSolidKHypoGP(nN, R1, Dmat + 36 * gp, sMat + 6 * gp, pg, in.biotAlpha,
                    gU, dvCur, Ks, nU);
  }

  for (int i = 0; i < nU; i++)
    for (int j = 0; j < nU; j++)
      K[i * n + j] = Ks[i * nU + j];
  for (int i = 0; i < nU; i++)
    for (int bb = 0; bb < nP; bb++) {
      K[i * n + nU + bb] = -Q[i * nP + bb];
      K[(nU + bb) * n + i] = c2 * Q[i * nP + bb];
    }
  for (int bb = 0; bb < nP; bb++)
    for (int cc = 0; cc < nP; cc++) {
      double ht = in.HtAl ? in.HtAl[bb * nP + cc] : 0.0;
      K[(nU + bb) * n + nU + cc] = H[bb * nP + cc] + c2 * (S[bb * nP + cc] + ht);
    }
  return true;
}

} // namespace ladruno_up_hypo

#endif
