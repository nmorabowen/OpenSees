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
// Created: 06/2026
//
// LadrunoJ2Finite — finite-strain-native combined-hardening J2 with a co-rotating
// kinematic backstress (dSNPO §14.11 "v2"). See LadrunoJ2Finite.h for the contract
// and Ladruno_implementation/16_finite_native_j2_adr.md for the derivation. The
// return map is the SHARED LadrunoJ2Kernel.h; the Hencky kinematics + spatial
// tangent are LogStrainKernel.h. Verified against tests/ladrunoj2_finite_native_reference.py
// and tests/test_ladrunoJ2_finite_native_tangent.py.

#include <LadrunoJ2Finite.h>
#include <LadrunoJ2Kernel.h>     // shared combined-hardening return map
#include <LogStrainKernel.h>     // Hencky kinematics + spatial tangent + jacobi3
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Response.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <math.h>

using ladruno_j2_kernel::Params;

// =========================================================================== //
//  small file-local helpers (3×3 row-major; reuse LogStrainKernel primitives)  //
// =========================================================================== //
static const double ZERO6[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

static void fillParams(Params &p, double K, double G, double sig0, double Qinf,
                       double bIso, double Hiso, int nBack,
                       const double *Ckin, const double *gKin)
{
  p.K = K; p.G = G;
  p.sig0 = sig0; p.Qinf = Qinf; p.bIso = bIso; p.Hiso = Hiso;
  p.nBack = nBack;
  for (int k = 0; k < ladruno_j2_kernel::MAXBACK; k++) { p.C[k] = Ckin[k]; p.gam[k] = gKin[k]; }
}

// Rotation part R of the polar decomposition F_Δ = R U (R orthogonal, U SPD), via
// the spectral square root: C = F_Δᵀ F_Δ (SPD), U⁻¹ = Σ_a λ_a^{-1/2} E_a,
// R = F_Δ U⁻¹. This is the incremental MATERIAL rotation that objectively
// transports the backstress (Hughes–Winget style; dSNPO §14.11). Reuses the
// degeneracy-robust cyclic-Jacobi eigensolver from LogStrainKernel.
static void polarRotation(const double Fd[9], double R[9])
{
  using namespace logstrain_kernel;
  double FdT[9]; mat3_transpose(Fd, FdT);
  double C[9];   mat3_mul(FdT, Fd, C);
  double val[3], vec[9];
  jacobi3(C, val, vec);
  // Clamp the squared-stretch eigenvalues to a RELATIVE positive floor before the
  // inverse-sqrt. det F>0 (checked in setTrialF) and det Fn>0 ⇒ val[a]>0 strictly
  // in exact arithmetic, but the channel-B tangent perturbs F internally (bypassing
  // that guard), so a near-singular F_pert could otherwise give 1/sqrt(≈0)=Inf and a
  // rank-deficient (non-orthogonal) U⁻¹. The clamp keeps U⁻¹ full-rank and finite so
  // R stays ~orthogonal; the perturbation is O(1e-7) so the clamp never trips on a
  // well-conditioned step.
  double vmax = fmax(val[0], fmax(val[1], val[2]));
  double vfloor = 1.0e-30 * (vmax > 0.0 ? vmax : 1.0);
  double Uinv[9];
  for (int i = 0; i < 9; i++) Uinv[i] = 0.0;
  for (int a = 0; a < 3; a++) {
    double lam = (val[a] > vfloor) ? val[a] : vfloor;
    double s = 1.0 / sqrt(lam);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Uinv[3*i+j] += s * vec[3*i+a] * vec[3*j+a];
  }
  mat3_mul(Fd, Uinv, R);
}

// out = R · a · Rᵀ for a symmetric tensor a stored as 6 TENSOR components
// {00,11,22,01,12,02} (stress/backstress convention — off-diagonals are true
// tensor, NOT engineering). Result re-symmetrised.
static void rotateSym6(const double R[9], const double a6[6], double out6[6])
{
  using namespace logstrain_kernel;
  double A[9];
  A[0] = a6[0]; A[4] = a6[1]; A[8] = a6[2];
  A[1] = A[3] = a6[3];
  A[5] = A[7] = a6[4];
  A[2] = A[6] = a6[5];
  double RA[9];  mat3_mul(R, A, RA);
  double RT[9];  mat3_transpose(R, RT);
  double RAR[9]; mat3_mul(RA, RT, RAR);
  out6[0] = RAR[0]; out6[1] = RAR[4]; out6[2] = RAR[8];
  out6[3] = 0.5*(RAR[1] + RAR[3]);
  out6[4] = 0.5*(RAR[5] + RAR[7]);
  out6[5] = 0.5*(RAR[2] + RAR[6]);
}

#ifdef LADRUNO_J2FINITE_CHANNELB_NUMERIC
// Cauchy σ (as a 3×3) for a given backstress rotation R, holding the trial elastic
// strain and J fixed — the channel-B (R-only) probe for the NUMERIC tangent fallback.
static void cauchyFromR(const Params &p, const double strain6_t[6], double ebarP_n,
                        const double alpha_n_sp[][6], int nBack,
                        const double R[9], double invJ, double sig3[3][3])
{
  double at[ladruno_j2_kernel::MAXBACK][6];
  for (int k = 0; k < nBack; k++) rotateSym6(R, alpha_n_sp[k], at[k]);
  double s6[6], D[6][6], epsP[6], eb, aout[ladruno_j2_kernel::MAXBACK][6], dG;
  ladruno_j2_kernel::returnMap(p, strain6_t, ZERO6, ebarP_n, at,
                               s6, D, epsP, eb, aout, dG);
  sig3[0][0] = s6[0]*invJ; sig3[1][1] = s6[1]*invJ; sig3[2][2] = s6[2]*invJ;
  sig3[0][1] = sig3[1][0] = s6[3]*invJ;
  sig3[1][2] = sig3[2][1] = s6[4]*invJ;
  sig3[0][2] = sig3[2][0] = s6[5]*invJ;
}
#endif  // LADRUNO_J2FINITE_CHANNELB_NUMERIC

// --------------------------------------------------------------------------- //
//  ANALYTIC channel B helpers (ADR-16 P3). Replace the numeric R-perturbation  //
//  with the exact chain  ∂R/∂F (polar derivative) ∘ ∂α̃/∂R ∘ ∂τ/∂α̃ .          //
//  Oracle: tests/ladrunoj2_finite_channelB_reference.py.                       //
// --------------------------------------------------------------------------- //
static inline double ddot3(const double A[9], const double B[9])  // Frobenius A:B
{ double s = 0.0; for (int i = 0; i < 9; i++) s += A[i]*B[i]; return s; }

// dR for f = R U (polar), perturbation df. Ω skew solves Ω U + U Ω = A − Aᵀ
// (A = Rᵀ df); axial form (tr U·I − U) ω = axial(A−Aᵀ); dR = R Ω. R = polar(f) given.
static void polarRotationDeriv(const double f[9], const double df[9],
                               const double R[9], double dR[9])
{
  using namespace logstrain_kernel;
  double RT[9]; mat3_transpose(R, RT);
  double U[9];  mat3_mul(RT, f, U);                 // right stretch (SPD)
  for (int i = 0; i < 3; i++)                       // symmetrise
    for (int j = i+1; j < 3; j++) { double a = 0.5*(U[3*i+j]+U[3*j+i]); U[3*i+j]=U[3*j+i]=a; }
  double A[9]; mat3_mul(RT, df, A);
  double w[3] = { A[7]-A[5], A[2]-A[6], A[3]-A[1] }; // axial(A − Aᵀ)
  double trU = U[0]+U[4]+U[8];
  double K[9];                                       // trU·I − U  (SPD)
  for (int i = 0; i < 9; i++) K[i] = -U[i];
  K[0]+=trU; K[4]+=trU; K[8]+=trU;
  double Kinv[9]; mat3_inv(K, Kinv);
  double om[3];                                      // ω = K⁻¹ w
  for (int i = 0; i < 3; i++) om[i] = Kinv[3*i+0]*w[0]+Kinv[3*i+1]*w[1]+Kinv[3*i+2]*w[2];
  double Om[9] = {0.0,-om[2],om[1], om[2],0.0,-om[0], -om[1],om[0],0.0};  // skew(ω)
  mat3_mul(R, Om, dR);
}

// Channel-B  dSdF[i][j][k][m] = ∂σ_ij/∂F_km  (through the co-rotated backstress only;
// εᵉᵗʳ and J held fixed). aChi_p are the co-rotated backstresses α̃_p = R α_p Rᵀ as
// 3×3; aRef_p are the committed α_p as 3×3 (so dα̃ = dR α_p Rᵀ + R α_p dRᵀ). dG is the
// converged plastic multiplier (> 0 here). Mirrors ladrunoj2_finite_channelB_reference.
static void channelBAnalytic(const Params &p, const double strain6_t[6], double ebarP_n,
                             const double aRef[][9], int nBack, const double R[9],
                             const double Fd[9], const double Fninv[9], double dG,
                             double invJ, double dSdF[3][3][3][3])
{
  using namespace logstrain_kernel;
  const double G = p.G, root23 = sqrt(2.0/3.0);

  // co-rotated backstresses α̃_p (3×3) + aux at the converged state
  double aChi[ladruno_j2_kernel::MAXBACK][9];
  for (int q = 0; q < nBack; q++) {
    double RA[9]; mat3_mul(R, aRef[q], RA);
    double RT[9]; mat3_transpose(R, RT);
    mat3_mul(RA, RT, aChi[q]);
  }
  double eps[9] = { strain6_t[0],strain6_t[3],strain6_t[5],
                    strain6_t[3],strain6_t[1],strain6_t[4],
                    strain6_t[5],strain6_t[4],strain6_t[2] };
  double trE = eps[0]+eps[4]+eps[8];
  double sTr[9]; for (int i=0;i<9;i++) sTr[i] = 2.0*G*eps[i];
  sTr[0]-=2.0*G*trE/3.0; sTr[4]-=2.0*G*trE/3.0; sTr[8]-=2.0*G*trE/3.0;
  double Dk[ladruno_j2_kernel::MAXBACK];
  double M[9]; for (int i=0;i<9;i++) M[i]=sTr[i];
  double Mp[9] = {0,0,0,0,0,0,0,0,0};
  for (int q = 0; q < nBack; q++) {
    Dk[q] = 1.0 + root23*p.gam[q]*dG;
    for (int i=0;i<9;i++) { M[i]-=aChi[q][i]/Dk[q]; Mp[i]+=aChi[q][i]*root23*p.gam[q]/(Dk[q]*Dk[q]); }
  }
  double normM = sqrt(ddot3(M,M));
  double n[9]; for (int i=0;i<9;i++) n[i]=M[i]/normM;
  double dtheta = 2.0*G;
  for (int q=0;q<nBack;q++) dtheta += (2.0/3.0)*p.C[q]/(Dk[q]*Dk[q]);
  double pbar = ebarP_n + root23*dG;
  double h = dtheta + (2.0/3.0)*ladruno_j2_kernel::yieldSlope(p, pbar) - ddot3(n, Mp);

  for (int k = 0; k < 3; k++)
    for (int mIdx = 0; mIdx < 3; mIdx++) {
      double dFd[9] = {0,0,0,0,0,0,0,0,0};          // dFd = e_{k,mIdx} · Fn⁻¹ : row k = Fninv row mIdx
      dFd[3*k+0]=Fninv[3*mIdx+0]; dFd[3*k+1]=Fninv[3*mIdx+1]; dFd[3*k+2]=Fninv[3*mIdx+2];
      double dR[9]; polarRotationDeriv(Fd, dFd, R, dR);
      double dRT[9]; mat3_transpose(dR, dRT);
      double RT[9];  mat3_transpose(R, RT);
      double ds[9] = {0,0,0,0,0,0,0,0,0};            // Σ_q ∂τ_dev/∂α̃_q : dα̃_q
      for (int q = 0; q < nBack; q++) {
        double t1[9], dA[9], t2[9];                  // dα̃_q = dR A_q Rᵀ + R A_q dRᵀ
        mat3_mul(dR, aRef[q], t1); mat3_mul(t1, RT, dA);
        mat3_mul(R, aRef[q], t1);  mat3_mul(t1, dRT, t2);
        for (int i=0;i<9;i++) dA[i]+=t2[i];
        double dGp = -ddot3(n, dA)/(h*Dk[q]);        // ∂Δγ/∂α̃_q : dα̃_q
        double dM[9]; for (int i=0;i<9;i++) dM[i] = -dA[i]/Dk[q] + Mp[i]*dGp;
        double ndM = ddot3(n, dM);
        for (int i=0;i<9;i++) {
          double dn = (dM[i] - n[i]*ndM)/normM;
          ds[i] += -2.0*G*(dGp*n[i] + dG*dn);
        }
      }
      // dσ_ij = ds_ij/J  (channel-B, deviatoric)
      for (int i=0;i<3;i++) for (int j=0;j<3;j++) dSdF[i][j][k][mIdx] = ds[3*i+j]*invJ;
    }
}

// =========================================================================== //
//  Factory:  nDMaterial LadrunoJ2Finite tag K G -iso voce s0 Qinf b Hiso        //
//                                       -kin N C1 g1 ... [-rho rho]            //
// =========================================================================== //
void *OPS_LadrunoJ2Finite(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoJ2Finite tag? K? G? "
           << "-iso voce s0? Qinf? b? Hiso? -kin N? C1? g1? ... <-rho rho?> <-implex>\n";
    opserr << "  -implex: explicit Delta-gamma-extrapolation tangent (constant SPD). "
           << "ASSUMES UNIFORM TIME STEPS (the extrapolation is Delta-gamma~ = Delta-gamma_n, "
           << "i.e. Delta-t ratio = 1); accuracy degrades under variable Delta-t.\n";
    return 0;
  }

  int tag; int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoJ2Finite tag\n"; return 0;
  }
  double KG[2]; numData = 2;
  if (OPS_GetDoubleInput(&numData, KG) < 0) {
    opserr << "WARNING invalid LadrunoJ2Finite K G\n"; return 0;
  }
  double K = KG[0], G = KG[1];

  double s0 = 0.0, Qinf = 0.0, bIso = 0.0, Hiso = 0.0, rho = 0.0;
  int nBack = 0;
  bool useImplex = false;
  double C[LadrunoJ2Finite::MAXBACK], gam[LadrunoJ2Finite::MAXBACK];
  for (int i = 0; i < LadrunoJ2Finite::MAXBACK; i++) { C[i] = 0.0; gam[i] = 0.0; }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *flag = OPS_GetString();
    if (strcmp(flag, "-iso") == 0) {
      const char *law = OPS_GetString();
      if (strcmp(law, "voce") == 0 || strcmp(law, "Voce") == 0) {
        double d[4]; numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoJ2Finite: -iso voce wants s0 Qinf b Hiso\n"; return 0;
        }
        s0 = d[0]; Qinf = d[1]; bIso = d[2]; Hiso = d[3];
      } else {
        opserr << "WARNING LadrunoJ2Finite: unknown -iso law '" << law
               << "' (only 'voce' in v1)\n"; return 0;
      }
    }
    else if (strcmp(flag, "-kin") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nBack) < 0) {
        opserr << "WARNING LadrunoJ2Finite: -kin wants N\n"; return 0;
      }
      if (nBack < 0 || nBack > LadrunoJ2Finite::MAXBACK) {
        opserr << "WARNING LadrunoJ2Finite: -kin N out of range [0,"
               << LadrunoJ2Finite::MAXBACK << "]\n"; return 0;
      }
      for (int k = 0; k < nBack; k++) {
        double cg[2]; numData = 2;
        if (OPS_GetDoubleInput(&numData, cg) < 0) {
          opserr << "WARNING LadrunoJ2Finite: -kin term " << k+1 << " wants C gamma\n"; return 0;
        }
        C[k] = cg[0]; gam[k] = cg[1];
      }
    }
    else if (strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) {
        opserr << "WARNING LadrunoJ2Finite: -rho wants rho\n"; return 0;
      }
    }
    else if (strcmp(flag, "-implex") == 0) {
      useImplex = true;
    }
    else {
      opserr << "WARNING LadrunoJ2Finite: unknown option '" << flag << "'\n"; return 0;
    }
  }

  // Softening-parameter guard (LEDGER_quirks): the consistent-tangent denominator is
  // h = 2G + (2/3)·σ_y'(p̄) + (kinematic ≥ 0) − dotT(n,Mp). The binding (pure-iso) floor
  // is the MINIMUM isotropic hardening slope over p̄ ≥ 0:
  //   σ_y'(p̄) = Qinf·b·exp(−b·p̄) + Hiso  ⇒  min = Hiso + (Qinf<0 ? Qinf·b : 0).
  // If that drops below 0 the material softens and the implicit consistent tangent can
  // become indefinite (Newton stall); below −3G it loses positive-definiteness outright.
  // We WARN rather than reject: with -implex the reported tangent is the constant SPD
  // elastic operator regardless, and crack-band regularization is a valid use too.
  {
    double slopeMin = Hiso + ((Qinf < 0.0) ? Qinf * bIso : 0.0);
    if (slopeMin < 0.0 && !useImplex) {
      opserr << "WARNING LadrunoJ2Finite (tag " << tag << "): isotropic hardening softens "
             << "(min slope sig_y' = " << slopeMin << " < 0";
      if (slopeMin <= -3.0 * G)
        opserr << " <= -3G = " << (-3.0 * G)
               << "; the consistent tangent LOSES positive-definiteness";
      else
        opserr << "; the consistent tangent may become indefinite";
      opserr << "). Use -implex for a stable constant-SPD tangent, or add crack-band "
             << "regularization. (Hardening params Hiso=" << Hiso << ", Qinf=" << Qinf
             << ", b=" << bIso << ".)\n";
    }
  }

  NDMaterial *mat = new LadrunoJ2Finite(tag, K, G, s0, Qinf, bIso, Hiso, nBack, C, gam, rho,
                                        useImplex);
  if (mat == 0) { opserr << "WARNING LadrunoJ2Finite: allocation failed\n"; return 0; }
  return mat;
}

// =========================================================================== //
//  Construction                                                               //
// =========================================================================== //
void LadrunoJ2Finite::setIdentity(double M[9])
{
  for (int i = 0; i < 9; i++) M[i] = 0.0;
  M[0] = M[4] = M[8] = 1.0;
}

LadrunoJ2Finite::LadrunoJ2Finite(int tag, double K, double G,
                                 double s0, double Qi, double b, double Hi,
                                 int nb, const double *C, const double *gam, double r,
                                 bool implex)
  : FiniteStrainNDMaterial(tag, ND_TAG_LadrunoJ2Finite),
    bulk(K), shear(G), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), rho(r),
    nBack(nb), useImplex(implex), ebarP_n(0.0), dGamma_n(0.0),
    ebarP_trial(0.0), dGammaTrial(0.0),
    sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), K0init(6, 6), Jdet(1.0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->revertToStart();
}

LadrunoJ2Finite::LadrunoJ2Finite()
  : FiniteStrainNDMaterial(0, ND_TAG_LadrunoJ2Finite),
    bulk(0.0), shear(0.0), sig0(0.0), Qinf(0.0), bIso(0.0), Hiso(0.0), rho(0.0),
    nBack(0), useImplex(false), ebarP_n(0.0), dGamma_n(0.0),
    ebarP_trial(0.0), dGammaTrial(0.0),
    sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), K0init(6, 6), Jdet(1.0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->revertToStart();
}

LadrunoJ2Finite::~LadrunoJ2Finite() {}

// =========================================================================== //
//  The finite-strain seam: setTrialF (Box 14.3 + §14.11 backstress push)      //
// =========================================================================== //
int LadrunoJ2Finite::setTrialF(const Matrix &F)
{
  using namespace logstrain_kernel;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Ftrial9[3*i+j] = F(i, j);

  Jdet = mat3_det(Ftrial9);
  if (Jdet <= 0.0) {
    opserr << "LadrunoJ2Finite::setTrialF - non-positive det F (" << Jdet
           << "); element inverted/degenerate\n";
    return -1;
  }

  // f_Δ = F F_n⁻¹, trial bᵉᵗʳ = f_Δ bᵉ_n f_Δᵀ, trial Hencky strain εᵉᵗʳ
  double Fninv[9]; mat3_inv(Fn, Fninv);
  double Fd[9];    mat3_mul(Ftrial9, Fninv, Fd);
  double Betr[9];  trial_Be(Ftrial9, Fn, Be_n, Betr);
  for (int i = 0; i < 9; i++) BeTrial9[i] = Betr[i];

  double epsEng[6];                                  // engineering-shear Voigt
  hencky_voigt(Betr, epsEng);
  for (int k = 0; k < 6; k++) henckyStrain(k) = epsEng[k];
  double strain6_t[6];                               // tensor comps for the kernel
  strain6_t[0] = epsEng[0]; strain6_t[1] = epsEng[1]; strain6_t[2] = epsEng[2];
  strain6_t[3] = 0.5*epsEng[3]; strain6_t[4] = 0.5*epsEng[4]; strain6_t[5] = 0.5*epsEng[5];

  // co-rotate the committed SPATIAL backstresses into the current frame
  double R[9]; polarRotation(Fd, R);
  double alpha_t[MAXBACK][6];
  for (int k = 0; k < nBack; k++) rotateSym6(R, alpha_n[k], alpha_t[k]);

  // the UNCHANGED return map, elastic-trial form (epsP_n = 0; bᵉ_n carries history)
  Params p; fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  double s6[6], epsP[6], resid = 0.0;
  int status = ladruno_j2_kernel::returnMap(p, strain6_t, ZERO6, ebarP_n, alpha_t,
                                            s6, DtanTrial, epsP, ebarP_trial,
                                            alpha_trial, dGammaTrial, &resid);
  if (status == ladruno_j2_kernel::STATUS_SINGULAR)
    opserr << "WARNING LadrunoJ2Finite: singular local Jacobian (tag " << this->getTag() << ")\n";
  else if (status == ladruno_j2_kernel::STATUS_NO_CONVERGE)
    opserr << "WARNING LadrunoJ2Finite: local Newton did not converge (tag "
           << this->getTag() << ", |R|=" << resid << ")\n";

  // Recover the trial log-strain flow direction N = Δεᵖ/Δγ (current frame, unit
  // deviatoric, tensor comps) for the IMPL-EX extrapolation history. epsP is the
  // tensor plastic increment dG·n, so n = epsP/dG. Zero on an elastic step.
  if (dGammaTrial > 0.0) {
    double inv = 1.0 / dGammaTrial;
    for (int k = 0; k < 6; k++) Nflow_trial[k] = epsP[k] * inv;
  } else {
    for (int k = 0; k < 6; k++) Nflow_trial[k] = 0.0;
  }

  // ----- IMPL-EX (Oliver–Huespe–Cante): override the REPORTED stress + tangent ----
  // The implicit return map above stays the committed truth (epsP, alpha_trial,
  // ebarP_trial, Be below are untouched). Only the stress fed to the solver and the
  // reported tangent are replaced: freeze the extrapolated multiplier Δγ̃ = Δγ_n
  // (uniform step) and the committed flow direction N_n co-rotated to the current
  // frame, so the plastic strain εᵖ̃ = Δγ̃·Ñ_n is a pure history quantity and the
  // stress is a linear-elastic evaluation ⇒ the tangent is the constant SPD elastic
  // operator (no plastic h, no backstress co-rotation channel B). One-step lag at
  // first yield (Δγ_n = 0). See Ladruno_implementation/16_finite_native_j2_adr.md and
  // tests/ladrunoj2_finite_implex_reference.py.
  if (useImplex) {
    double Nflow_pf[6]; rotateSym6(R, Nflow_n, Nflow_pf);   // co-rotate committed N_n
    const double lam = bulk - 2.0 * shear / 3.0;
    double ee[6];                                           // εᵉᵗʳ − Δγ̃·Ñ_n (tensor)
    for (int k = 0; k < 6; k++) ee[k] = strain6_t[k] - dGamma_n * Nflow_pf[k];
    double tre = ee[0] + ee[1] + ee[2];
    s6[0] = lam*tre + 2.0*shear*ee[0];
    s6[1] = lam*tre + 2.0*shear*ee[1];
    s6[2] = lam*tre + 2.0*shear*ee[2];
    s6[3] = 2.0*shear*ee[3];
    s6[4] = 2.0*shear*ee[4];
    s6[5] = 2.0*shear*ee[5];
    ladruno_j2_kernel::elasticTangent(p, DtanTrial);        // constant SPD reported tangent
  }

  // Cauchy σ = τ/J (s6 is the Kirchhoff τ — implicit, or the IMPL-EX explicit τ̃)
  double invJ = 1.0 / Jdet;
  for (int k = 0; k < 6; k++) sigmaCauchy(k) = s6[k] * invJ;

  // channel-A spatial constitutive 6×6 (lossy projection): c = (1/2J)[D:L:B]
  double D6flat[36];
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) D6flat[6*I+J] = DtanTrial[I][J];
  double sig6[6], c6[36];
  assemble_material(Betr, D6flat, s6, Jdet, sig6, c6);
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) aTangent(I, J) = c6[6*I+J];

  // updated elastic strain εᵉ_{n+1} = εᵉᵗʳ − Δεᵖ (epsP is the tensor plastic
  // increment dG·n), then bᵉ to commit = exp[2 εᵉ_{n+1}] (engineering feed).
  double epsE_eng[6];
  epsE_eng[0] = epsEng[0] - epsP[0];
  epsE_eng[1] = epsEng[1] - epsP[1];
  epsE_eng[2] = epsEng[2] - epsP[2];
  epsE_eng[3] = epsEng[3] - 2.0*epsP[3];
  epsE_eng[4] = epsEng[4] - 2.0*epsP[4];
  epsE_eng[5] = epsEng[5] - 2.0*epsP[5];
  be_from_hencky_voigt(epsE_eng, Be_trialUpd);
  return 0;
}

// =========================================================================== //
//  Query                                                                      //
// =========================================================================== //
const Vector &LadrunoJ2Finite::getStress(void)  { return sigmaCauchy; }
const Matrix &LadrunoJ2Finite::getTangent(void) { return aTangent; }
const Vector &LadrunoJ2Finite::getStrain(void)  { return henckyStrain; }

const Matrix &LadrunoJ2Finite::getInitialTangent(void)
{
  // per-instance member (NOT a function-static) so concurrent / cached queries on
  // distinct (K,G) materials never alias — mirrors LadrunoJ2::tangentOut.
  double Kt[6][6];
  Params p; fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  ladruno_j2_kernel::elasticTangent(p, Kt);
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) K0init(I, J) = Kt[I][J];
  return K0init;
}

// FULL 4th-order spatial constitutive modulus c_ijkl. Channel A (log-strain,
// analytic) + channel B (backstress co-rotation, plastic steps only). The element
// forms a_ijkl = c_ijkl − σ_il δ_jk. See the header. Channel B is ANALYTIC by default
// (∂R/∂F polar derivative ∘ ∂τ/∂α̃ return-map sensitivity); define
// LADRUNO_J2FINITE_CHANNELB_NUMERIC to fall back to the numeric R-perturbation
// (validation: both agree with tests/ladrunoj2_finite_channelB_reference.py).
int LadrunoJ2Finite::getSpatialTangentTensor(double c[3][3][3][3])
{
  using namespace logstrain_kernel;
  double D6flat[36];
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) D6flat[6*I+J] = DtanTrial[I][J];
  spatial_tangent_full(BeTrial9, D6flat, Jdet, c);           // channel A

  // IMPL-EX reports the constant SPD elastic tangent (DtanTrial set elastic in
  // setTrialF): channel A with the elastic D IS the whole tangent — no co-rotation
  // channel B (the backstress drops out of the explicit stress).
  if (useImplex) return 0;
  if (dGammaTrial <= 0.0) return 0;                          // elastic ⇒ no channel B

  // common setup: εᵉᵗʳ (tensor comps), f_Δ, R = polar(f_Δ), Fn⁻¹.
  Params p; fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  double epsEng[6]; hencky_voigt(BeTrial9, epsEng);
  double strain6_t[6];
  strain6_t[0] = epsEng[0]; strain6_t[1] = epsEng[1]; strain6_t[2] = epsEng[2];
  strain6_t[3] = 0.5*epsEng[3]; strain6_t[4] = 0.5*epsEng[4]; strain6_t[5] = 0.5*epsEng[5];
  double Fninv[9]; mat3_inv(Fn, Fninv);
  double Fd[9];    mat3_mul(Ftrial9, Fninv, Fd);
  double R[9];     polarRotation(Fd, R);
  double invJ = 1.0 / Jdet;
  double dSdF[3][3][3][3];                                   // ∂σ_ij/∂F_km (channel B)

#ifdef LADRUNO_J2FINITE_CHANNELB_NUMERIC
  // numeric fallback: ∂σ/∂F by central-differencing R = polar(f_Δ) through the map.
  const double hFD = 1.0e-7;
  for (int k = 0; k < 3; k++)
    for (int m = 0; m < 3; m++) {
      double Fp[9], Fm[9];
      for (int q = 0; q < 9; q++) { Fp[q] = Ftrial9[q]; Fm[q] = Ftrial9[q]; }
      Fp[3*k+m] += hFD; Fm[3*k+m] -= hFD;
      double Fdp[9]; mat3_mul(Fp, Fninv, Fdp); double Rp[9]; polarRotation(Fdp, Rp);
      double Fdm[9]; mat3_mul(Fm, Fninv, Fdm); double Rm[9]; polarRotation(Fdm, Rm);
      double sp[3][3], sm[3][3];
      cauchyFromR(p, strain6_t, ebarP_n, alpha_n, nBack, Rp, invJ, sp);
      cauchyFromR(p, strain6_t, ebarP_n, alpha_n, nBack, Rm, invJ, sm);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) dSdF[i][j][k][m] = (sp[i][j] - sm[i][j]) / (2.0*hFD);
    }
#else
  // analytic channel B (default): committed backstresses as 3×3, then the exact chain.
  double aRef[ladruno_j2_kernel::MAXBACK][9];
  for (int q = 0; q < nBack; q++) {
    aRef[q][0]=alpha_n[q][0]; aRef[q][4]=alpha_n[q][1]; aRef[q][8]=alpha_n[q][2];
    aRef[q][1]=aRef[q][3]=alpha_n[q][3];
    aRef[q][5]=aRef[q][7]=alpha_n[q][4];
    aRef[q][2]=aRef[q][6]=alpha_n[q][5];
  }
  channelBAnalytic(p, strain6_t, ebarP_n, aRef, nBack, R, Fd, Fninv, dGammaTrial, invJ, dSdF);
#endif

  // map into cmat via Ḟ = L F:  cmatB_ijkl = Σ_m (∂σ_ij/∂F_km) F_lm.
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) {
          double acc = 0.0;
          for (int m = 0; m < 3; m++) acc += dSdF[i][j][k][m] * Ftrial9[3*l+m];
          c[i][j][k][l] += acc;
        }
  return 0;
}

// =========================================================================== //
//  State cycle                                                                //
// =========================================================================== //
int LadrunoJ2Finite::commitState(void)
{
  for (int i = 0; i < 9; i++) { Fn[i] = Ftrial9[i]; Be_n[i] = Be_trialUpd[i]; }
  ebarP_n = ebarP_trial;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha_trial[k][i];
  // IMPL-EX extrapolation history: commit the implicit multiplier increment + the
  // recovered flow direction (becomes N_n, co-rotated next step).
  dGamma_n = dGammaTrial;
  for (int i = 0; i < 6; i++) Nflow_n[i] = Nflow_trial[i];
  return 0;
}

int LadrunoJ2Finite::revertToLastCommit(void)
{
  return 0;   // trial state is recomputed from the committed state on next setTrialF
}

int LadrunoJ2Finite::revertToStart(void)
{
  setIdentity(Fn); setIdentity(Be_n);
  setIdentity(Ftrial9); setIdentity(BeTrial9); setIdentity(Be_trialUpd);
  ebarP_n = 0.0; ebarP_trial = 0.0; dGammaTrial = 0.0; dGamma_n = 0.0; Jdet = 1.0;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) { alpha_n[k][i] = 0.0; alpha_trial[k][i] = 0.0; }
  for (int i = 0; i < 6; i++) { Nflow_n[i] = 0.0; Nflow_trial[i] = 0.0; }
  sigmaCauchy.Zero();
  henckyStrain.Zero();
  aTangent.Zero();
  Params p; fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  ladruno_j2_kernel::elasticTangent(p, DtanTrial);
  return 0;
}

// =========================================================================== //
//  Copies                                                                     //
// =========================================================================== //
NDMaterial *LadrunoJ2Finite::getCopy(void)
{
  LadrunoJ2Finite *c = new LadrunoJ2Finite(this->getTag(), bulk, shear, sig0, Qinf,
                                           bIso, Hiso, nBack, Ckin, gKin, rho, useImplex);
  for (int i = 0; i < 9; i++) { c->Fn[i] = Fn[i]; c->Be_n[i] = Be_n[i]; }
  c->ebarP_n = ebarP_n;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) c->alpha_n[k][i] = alpha_n[k][i];
  c->dGamma_n = dGamma_n;
  for (int i = 0; i < 6; i++) c->Nflow_n[i] = Nflow_n[i];
  return c;
}

NDMaterial *LadrunoJ2Finite::getCopy(const char *type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    return this->getCopy();
  opserr << "LadrunoJ2Finite::getCopy - only ThreeDimensional is supported\n";
  return 0;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int LadrunoJ2Finite::sendSelf(int commitTag, Channel &theChannel)
{
  // tag + 7 params + nBack + 2*MAXBACK + Fn(9) + Be_n(9) + ebarP_n(1) + alpha_n(MAXBACK*6)
  // + useImplex(1) + dGamma_n(1) + Nflow_n(6)
  static Vector data(1 + 7 + 1 + 2*MAXBACK + 9 + 9 + 1 + MAXBACK*6 + 1 + 1 + 6);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = bulk; data(c++) = shear;
  data(c++) = sig0; data(c++) = Qinf; data(c++) = bIso; data(c++) = Hiso;
  data(c++) = rho;
  data(c++) = nBack;
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  for (int i = 0; i < 9; i++) data(c++) = Fn[i];
  for (int i = 0; i < 9; i++) data(c++) = Be_n[i];
  data(c++) = ebarP_n;
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) data(c++) = alpha_n[k][i];
  data(c++) = useImplex ? 1.0 : 0.0;
  data(c++) = dGamma_n;
  for (int i = 0; i < 6; i++) data(c++) = Nflow_n[i];

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2Finite::sendSelf - failed to send vector\n"; return -1;
  }
  return 0;
}

int LadrunoJ2Finite::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(1 + 7 + 1 + 2*MAXBACK + 9 + 9 + 1 + MAXBACK*6 + 1 + 1 + 6);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2Finite::recvSelf - failed to recv vector\n"; return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  bulk = data(c++); shear = data(c++);
  sig0 = data(c++); Qinf = data(c++); bIso = data(c++); Hiso = data(c++);
  rho = data(c++);
  nBack = (int)data(c++);
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  for (int i = 0; i < 9; i++) Fn[i] = data(c++);
  for (int i = 0; i < 9; i++) Be_n[i] = data(c++);
  ebarP_n = data(c++);
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = data(c++);
  useImplex = (data(c++) != 0.0);
  dGamma_n = data(c++);
  for (int i = 0; i < 6; i++) Nflow_n[i] = data(c++);

  // sync trial to committed
  for (int i = 0; i < 9; i++) { Ftrial9[i] = Fn[i]; BeTrial9[i] = Be_n[i]; Be_trialUpd[i] = Be_n[i]; }
  ebarP_trial = ebarP_n; dGammaTrial = 0.0;
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alpha_trial[k][i] = alpha_n[k][i];
  for (int i = 0; i < 6; i++) Nflow_trial[i] = Nflow_n[i];
  Params p; fillParams(p, bulk, shear, sig0, Qinf, bIso, Hiso, nBack, Ckin, gKin);
  ladruno_j2_kernel::elasticTangent(p, DtanTrial);
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void LadrunoJ2Finite::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoJ2Finite\", ";
    s << "\"implex\": " << (useImplex ? "true" : "false");
    s << "}";
  } else {
    s << endln;
    s << "LadrunoJ2Finite (finite-strain combined-hardening J2, co-rotating backstress)" << endln;
    s << "  tag    : " << this->getTag() << endln;
    s << "  K, G   : " << bulk << ", " << shear << endln;
    s << "  iso    : sig0=" << sig0 << " Qinf=" << Qinf << " b=" << bIso << " Hiso=" << Hiso << endln;
    s << "  nBack  : " << nBack << endln;
    for (int k = 0; k < nBack; k++)
      s << "    term " << k+1 << ": C=" << Ckin[k] << " gamma=" << gKin[k] << endln;
    s << "  rho    : " << rho << endln;
    s << "  implex : " << (useImplex ? "on (Delta-gamma extrapolation; assumes uniform dt)"
                                     : "off") << endln;
  }
}

int LadrunoJ2Finite::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "K") == 0)    return param.addObject(1, this);
  if (strcmp(argv[0], "G") == 0 || strcmp(argv[0], "mu") == 0) return param.addObject(2, this);
  if (strcmp(argv[0], "rho") == 0)  return param.addObject(3, this);
  if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "sig0") == 0) return param.addObject(4, this);
  if (strcmp(argv[0], "Hiso") == 0) return param.addObject(5, this);
  if (strcmp(argv[0], "Qinf") == 0) return param.addObject(6, this);
  if (strcmp(argv[0], "b") == 0)    return param.addObject(7, this);
  return -1;
}

int LadrunoJ2Finite::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case 1: bulk  = info.theDouble; return 0;
    case 2: shear = info.theDouble; return 0;
    case 3: rho   = info.theDouble; return 0;
    case 4: sig0  = info.theDouble; return 0;
    case 5: Hiso  = info.theDouble; return 0;
    case 6: Qinf  = info.theDouble; return 0;
    case 7: bIso  = info.theDouble; return 0;
    default: return -1;
  }
}

// =========================================================================== //
//  Recordable responses                                                       //
// =========================================================================== //
Response *LadrunoJ2Finite::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char *a = argv[0];
  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "backStress") == 0 || strcmp(a, "backstress") == 0 || strcmp(a, "alpha") == 0)
    return new MaterialResponse(this, 4, Vector(6));
  if (strcmp(a, "equivalentPlasticStrain") == 0 || strcmp(a, "plasticStrainEq") == 0 ||
      strcmp(a, "ebarP") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoJ2Finite::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
    case 1: if (matInfo.theVector) *(matInfo.theVector) = this->getStress(); return 0;
    case 2: if (matInfo.theVector) *(matInfo.theVector) = this->getStrain(); return 0;
    case 3: if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent(); return 0;
    case 4:                                  // total committed backstress (spatial, stress conv.)
      // tensor-component order {00,11,22,01,12,02}; slot 5 is (0,2) == (2,0) for a
      // symmetric tensor, i.e. interchangeable with the getStress {..,12,20} contract.
      if (matInfo.theVector) {
        Vector &v = *(matInfo.theVector);
        for (int i = 0; i < 6; i++) {
          double tot = 0.0;
          for (int k = 0; k < nBack; k++) tot += alpha_n[k][i];
          v(i) = tot;
        }
      }
      return 0;
    case 6: if (matInfo.theVector) (*(matInfo.theVector))(0) = ebarP_n; return 0;
    default: return -1;
  }
}
