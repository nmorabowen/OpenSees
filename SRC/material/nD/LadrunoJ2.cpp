/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: combined isotropic (Voce+linear) + Chaboche Armstrong-Frederick
// kinematic von Mises (J2) nDMaterial. See LadrunoJ2.h and
// Ladruno_implementation/10_ladruno_j2_plasticity.md for the full derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoJ2.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

// ---- symmetric-tensor helpers (6 tensor components {00,11,22,01,12,02}) ----
// Double contraction A:B with the factor-2 weights on the off-diagonal pairs.
static inline double dotT(const double* a, const double* b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
       + 2.0*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]);
}

// 6x6 deviatoric projector (Belytschko IIdev) in the J2ThreeDimensional mapping:
// normal block (2/3,-1/3,...), shear-diagonal 0.5.
static const double IIDEV6[6][6] = {
  { 2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0 },
  {-1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0 },
  {-1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0 },
  { 0.0, 0.0, 0.0, 0.5, 0.0, 0.0 },
  { 0.0, 0.0, 0.0, 0.0, 0.5, 0.0 },
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5 }
};

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoJ2 tag K G -iso voce s0 Qinf b Hiso -kin N C1 g1 ... [-rho r]
// ===========================================================================
void* OPS_LadrunoJ2(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoJ2 tag? K? G? "
           << "-iso voce s0? Qinf? b? Hiso? -kin N? C1? g1? ... <-rho rho?>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoJ2 tag\n";
    return 0;
  }

  double KG[2];
  numData = 2;
  if (OPS_GetDoubleInput(&numData, KG) < 0) {
    opserr << "WARNING invalid LadrunoJ2 K G\n";
    return 0;
  }
  double K = KG[0], G = KG[1];

  // defaults: perfectly plastic, no kinematic hardening
  double s0 = 0.0, Qinf = 0.0, bIso = 0.0, Hiso = 0.0, rho = 0.0;
  int nBack = 0;
  double C[LadrunoJ2::MAXBACK], gam[LadrunoJ2::MAXBACK];
  for (int i = 0; i < LadrunoJ2::MAXBACK; i++) { C[i] = 0.0; gam[i] = 0.0; }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();

    if (strcmp(flag, "-iso") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "voce") == 0 || strcmp(law, "Voce") == 0) {
        double d[4];
        numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoJ2: -iso voce wants s0 Qinf b Hiso\n";
          return 0;
        }
        s0 = d[0]; Qinf = d[1]; bIso = d[2]; Hiso = d[3];
      } else {
        opserr << "WARNING LadrunoJ2: unknown -iso law '" << law
               << "' (only 'voce' in v1)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-kin") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nBack) < 0) {
        opserr << "WARNING LadrunoJ2: -kin wants N\n";
        return 0;
      }
      if (nBack < 0 || nBack > LadrunoJ2::MAXBACK) {
        opserr << "WARNING LadrunoJ2: -kin N out of range [0," << LadrunoJ2::MAXBACK << "]\n";
        return 0;
      }
      for (int k = 0; k < nBack; k++) {
        double cg[2];
        numData = 2;
        if (OPS_GetDoubleInput(&numData, cg) < 0) {
          opserr << "WARNING LadrunoJ2: -kin term " << k+1 << " wants C gamma\n";
          return 0;
        }
        C[k] = cg[0]; gam[k] = cg[1];
      }
    }
    else if (strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) {
        opserr << "WARNING LadrunoJ2: -rho wants rho\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoJ2: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  NDMaterial* mat = new LadrunoJ2(tag, K, G, s0, Qinf, bIso, Hiso, nBack, C, gam, rho);
  if (mat == 0) {
    opserr << "WARNING LadrunoJ2: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoJ2::LadrunoJ2()
  : NDMaterial(0, ND_TAG_LadrunoJ2),
    bulk(0.0), shear(0.0), sig0(0.0), Qinf(0.0), bIso(0.0), Hiso(0.0), rho(0.0),
    nBack(0), dim(DIM_3D), ncomp(6), condense(false), cEps22(0.0),
    dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->setupDim();
  this->revertToStart();
}

LadrunoJ2::LadrunoJ2(int tag, double K, double G,
                     double s0, double Qi, double b, double Hi,
                     int nb, const double* C, const double* gam, double r, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoJ2),
    bulk(K), shear(G), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), rho(r),
    nBack(nb), dim(dimMode), ncomp(6), condense(false), cEps22(0.0),
    dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->setupDim();
  this->revertToStart();
}

// ---------------------------------------------------------------------------
//  dimensional view setup: reduced-vector order -> full 6-comp index map
// ---------------------------------------------------------------------------
void LadrunoJ2::setupDim(void)
{
  // full ordering: 0:00 1:11 2:22 3:01 4:12 5:02
  switch (dim) {
    case DIM_PSTRAIN:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_AXISYM:     ncomp = 4; { int m[4]={0,1,2,3};     for(int a=0;a<4;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5};   for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

// Static condensation of the out-of-plane (33, index 2) dof so sigma_22 = 0:
//   Dtan[I][J] -= Dtan[I][2] Dtan[2][J] / Dtan[2][2].
void LadrunoJ2::condenseTangent(void)
{
  double d22 = Dtan[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan[i][2]; row[i] = Dtan[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] -= col[I]*row[J]/d22;
}

LadrunoJ2::~LadrunoJ2() {}

// ===========================================================================
//  hardening
// ===========================================================================
double LadrunoJ2::yieldStress(double p) const
{
  // Voce saturation + linear; reduces to perfect (Qinf=Hiso=0) or linear (Qinf=0)
  return sig0 + Qinf*(1.0 - exp(-bIso*p)) + Hiso*p;
}

double LadrunoJ2::yieldSlope(double p) const
{
  return Qinf*bIso*exp(-bIso*p) + Hiso;
}

// ===========================================================================
//  strain interface
// ===========================================================================
int LadrunoJ2::setTrialStrain(const Vector& e)
{
  // map the element's reduced (engineering-shear) vector into the full
  // 6-comp tensor strain. For condensed modes (PlaneStress/PlateFiber) the
  // out-of-plane eps_22 (index 2) is an internal unknown — preserve its
  // running value as the starting guess.
  double eps22 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps22;

  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double val = e(a);
    if (full >= 3) val *= 0.5;     // engineering -> tensor shear
    strain6[full] = val;
  }

  if (!condense) {
    this->integrate();
    return 0;
  }

  // enforce sigma_22 = 0: Newton on eps_22 (= strain6[2]); dSNPO sec 9.2.3.
  const int maxIt = 25;
  for (int it = 0; it < maxIt; it++) {
    this->integrate();
    double d22 = Dtan[2][2];
    double smag = 0.0;
    for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-10 * (smag > 1.0 ? smag : 1.0);
    if (fabs(stress6[2]) <= tol22) break;
    if (fabs(d22) < 1.0e-300) break;
    strain6[2] -= stress6[2] / d22;
    if (it == maxIt - 1)
      opserr << "WARNING LadrunoJ2: sigma_22 condensation did not converge (tag "
             << this->getTag() << ", |s22|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return 0;
}

int LadrunoJ2::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoJ2::setTrialStrainIncr(const Vector& v)
{
  // current reduced strain + increment
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double cur = (full >= 3) ? 2.0*strain6[full] : strain6[full];
    ne(a) = cur + v(a);
  }
  return this->setTrialStrain(ne);
}

int LadrunoJ2::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  the return map (scalar Newton, Kobayashi-Ohno) + consistent tangent
// ===========================================================================
void LadrunoJ2::integrate(void)
{
  const double G = shear, K = bulk;
  const double root23 = sqrt(2.0/3.0);

  // deviatoric strain (tensor comps)
  double tr = strain6[0] + strain6[1] + strain6[2];
  double edev[6];
  for (int i = 0; i < 6; i++) edev[i] = strain6[i];
  for (int i = 0; i < 3; i++) edev[i] -= tr/3.0;

  // trial deviatoric stress: s_tr = 2G (edev - epsP_n)
  double s_tr[6];
  for (int i = 0; i < 6; i++) s_tr[i] = 2.0*G*(edev[i] - epsP_n[i]);

  // trial relative stress M_tr = s_tr - sum_k alpha_k,n  (dG = 0 => denominators 1)
  double M[6];
  for (int i = 0; i < 6; i++) {
    M[i] = s_tr[i];
    for (int k = 0; k < nBack; k++) M[i] -= alpha_n[k][i];
  }
  double normM = sqrt(dotT(M, M));
  double sy0   = yieldStress(ebarP_n);
  double f_tr  = normM - root23*sy0;

  // Tolerances scaled to the operative stress magnitude (unit-aware; never
  // collapses to a fixed absolute when sig0 == 0, e.g. kinematic-only). The
  // normFloor guards the 1/||M|| divisions below — the relative-stress norm
  // ||M(dG)|| can collapse toward zero at a sharp reversal with backstress
  // (J2Plasticity guards the analogous norm_tau division; we restore that).
  const double stressScale = fmax(normM, root23*sy0);
  const double tolY      = 1.0e-12 * fmax(stressScale, 1.0e-300);
  const double normFloor = 1.0e-10 * fmax(stressScale, 1.0e-300);

  if (f_tr <= tolY) {
    // ---- elastic step: state frozen ----
    this->setStateToCommitted();
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    this->buildElasticTangent(Dtan);
    return;
  }

  // ---- plastic corrector: scalar Newton on dG ----
  double dG = 0.0;
  double Dk[MAXBACK];
  double Mp[6];            // d M / d dG
  double theta = 0.0, dtheta = 0.0, R = 1.0;
  int iter = 0;
  const int maxIter = 50;

  while (iter < maxIter) {
    // denominators
    for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*gKin[k]*dG;

    // M(dG) and its derivative wrt dG
    for (int i = 0; i < 6; i++) {
      M[i]  = s_tr[i];
      Mp[i] = 0.0;
      for (int k = 0; k < nBack; k++) {
        double inv = 1.0/Dk[k];
        M[i]  -= alpha_n[k][i]*inv;
        Mp[i] += alpha_n[k][i]*root23*gKin[k]*inv*inv;
      }
    }
    normM = sqrt(dotT(M, M));

    // theta(dG) = 2G dG + sum_k (2/3) C_k dG / Dk ; dtheta/ddG
    theta = 2.0*G*dG;
    dtheta = 2.0*G;
    for (int k = 0; k < nBack; k++) {
      double inv = 1.0/Dk[k];
      theta  += (2.0/3.0)*Ckin[k]*dG*inv;
      dtheta += (2.0/3.0)*Ckin[k]*inv*inv;   // d(dG/Dk)/ddG = 1/Dk^2
    }

    double xiNorm = normM - theta;                      // ||xi_{n+1}||
    double pbar   = ebarP_n + root23*dG;
    R = xiNorm - root23*yieldStress(pbar);

    // dR/ddG = (M:Mp)/||M|| - dtheta - (2/3) sig_y'   (guard 1/||M||)
    double dnormM = (normM > normFloor) ? dotT(M, Mp)/normM : 0.0;
    double dR = dnormM - dtheta - (2.0/3.0)*yieldSlope(pbar);

    if (fabs(R) <= tolY) break;
    if (dR == 0.0) {
      opserr << "WARNING LadrunoJ2: singular local Jacobian (tag " << this->getTag() << ")\n";
      break;
    }
    dG -= R/dR;
    if (dG < 0.0) dG = 0.0;   // keep admissible
    iter++;
  }
  if (iter >= maxIter)
    opserr << "WARNING LadrunoJ2: local Newton did not converge (tag "
           << this->getTag() << ", |R|=" << fabs(R) << ")\n";

  // recompute final quantities at converged dG
  for (int k = 0; k < nBack; k++) Dk[k] = 1.0 + root23*gKin[k]*dG;
  for (int i = 0; i < 6; i++) {
    M[i]  = s_tr[i];
    Mp[i] = 0.0;
    for (int k = 0; k < nBack; k++) {
      double inv = 1.0/Dk[k];
      M[i]  -= alpha_n[k][i]*inv;
      Mp[i] += alpha_n[k][i]*root23*gKin[k]*inv*inv;
    }
  }
  normM = sqrt(dotT(M, M));
  if (normM <= normFloor) {
    // Degenerate relative-stress direction (||M|| -> 0): no well-defined
    // plastic flow (a converged root needs ||M|| = theta + sqrt(2/3) sig_y > 0,
    // so this means no admissible plastic solution was found). Fall back to the
    // elastic predictor at the committed state rather than dividing by ~0.
    this->setStateToCommitted();
    for (int i = 0; i < 6; i++) stress6[i] = s_tr[i];
    for (int i = 0; i < 3; i++) stress6[i] += K*tr;
    this->buildElasticTangent(Dtan);
    return;
  }
  double n[6];
  for (int i = 0; i < 6; i++) n[i] = M[i]/normM;
  double pbar = ebarP_n + root23*dG;

  // ---- state update ----
  ebarP = pbar;
  dGammaTrial = dG;
  for (int i = 0; i < 6; i++) epsP[i] = epsP_n[i] + dG*n[i];
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++)
      alpha[k][i] = (alpha_n[k][i] + (2.0/3.0)*Ckin[k]*dG*n[i]) / Dk[k];

  // ---- stress ----
  for (int i = 0; i < 6; i++) stress6[i] = s_tr[i] - 2.0*G*dG*n[i];
  for (int i = 0; i < 3; i++) stress6[i] += K*tr;

  // ---- consistent tangent (dSNPO 7.213 structure + AF dn term) ----
  // dtheta currently holds d theta/d dG at converged dG (recompute to be safe)
  dtheta = 2.0*G;
  for (int k = 0; k < nBack; k++) {
    double inv = 1.0/Dk[k];
    dtheta += (2.0/3.0)*Ckin[k]*inv*inv;
  }
  double nMp = dotT(n, Mp);
  double h   = dtheta + (2.0/3.0)*yieldSlope(pbar) - nMp;   // = -df/ddG > 0

  double Mperp[6];
  for (int i = 0; i < 6; i++) Mperp[i] = Mp[i] - nMp*n[i];

  double G2 = G*G;
  double beta1   = 1.0 - 2.0*G*dG/normM;          // coeff on 2G * IIdev
  double betaNN  = 4.0*G2*dG/normM - 4.0*G2/h;    // coeff on n (x) n
  double betaMpN = -4.0*G2*dG/(h*normM);          // coeff on Mperp (x) n

  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) {
      double oneone = (I < 3 && J < 3) ? 1.0 : 0.0;
      Dtan[I][J] = K*oneone
                 + 2.0*G*beta1*IIDEV6[I][J]
                 + betaNN*n[I]*n[J]
                 + betaMpN*Mperp[I]*n[J];
    }
}

void LadrunoJ2::buildElasticTangent(double Kt[6][6]) const
{
  const double G = shear, K = bulk;
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) {
      double oneone = (I < 3 && J < 3) ? 1.0 : 0.0;
      Kt[I][J] = K*oneone + 2.0*G*IIDEV6[I][J];
    }
}

void LadrunoJ2::setStateToCommitted(void)
{
  for (int i = 0; i < 6; i++) epsP[i] = epsP_n[i];
  ebarP = ebarP_n;
  dGammaTrial = 0.0;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha[k][i] = alpha_n[k][i];
}

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoJ2::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];  // shear = true stress
  return stressOut;
}

const Vector& LadrunoJ2::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    strainOut(a) = (full >= 3) ? 2.0*strain6[full] : strain6[full];  // tensor -> engineering
  }
  return strainOut;
}

const Matrix& LadrunoJ2::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Dtan[vmap[a]][vmap[b]];
  return tangentOut;
}

const Matrix& LadrunoJ2::getInitialTangent(void)
{
  // elastic tangent mapped to the reduced view (uncondensed, matching the
  // J2Plasticity specializations' getInitialTangent convention).
  double Kt[6][6];
  this->buildElasticTangent(Kt);
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Kt[vmap[a]][vmap[b]];
  return tangentOut;
}

const char* LadrunoJ2::getType(void) const
{
  switch (dim) {
    case DIM_PSTRAIN:    return "PlaneStrain";
    case DIM_AXISYM:     return "AxiSymmetric";
    case DIM_PLATEFIBER: return "PlateFiber";
    case DIM_PSTRESS:    return "PlaneStress";
    default:             return "ThreeDimensional";
  }
}

int LadrunoJ2::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoJ2::commitState(void)
{
  for (int i = 0; i < 6; i++) epsP_n[i] = epsP[i];
  ebarP_n  = ebarP;
  dGamma_n = dGammaTrial;
  for (int k = 0; k < nBack; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha[k][i];
  cEps22 = strain6[2];   // converged out-of-plane strain (condensed modes)
  return 0;
}

int LadrunoJ2::revertToLastCommit(void)
{
  this->setStateToCommitted();
  if (condense) strain6[2] = cEps22;
  return 0;
}

int LadrunoJ2::revertToStart(void)
{
  for (int i = 0; i < 6; i++) {
    strain6[i] = 0.0;
    stress6[i] = 0.0;
    epsP_n[i]  = 0.0;
    epsP[i]    = 0.0;
  }
  ebarP_n = 0.0; ebarP = 0.0;
  dGamma_n = 0.0; dGammaTrial = 0.0;
  cEps22 = 0.0;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) { alpha_n[k][i] = 0.0; alpha[k][i] = 0.0; }
  this->buildElasticTangent(Dtan);
  return 0;
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoJ2::getCopy(void)
{
  return new LadrunoJ2(this->getTag(), bulk, shear, sig0, Qinf, bIso, Hiso,
                       nBack, Ckin, gKin, rho, dim);
}

NDMaterial* LadrunoJ2::getCopy(const char* type)
{
  int d = -1;
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) d = DIM_3D;
  else if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0) d = DIM_PSTRAIN;
  else if (strcmp(type, "AxiSymmetric") == 0 || strcmp(type, "AxiSymmetric2D") == 0) d = DIM_AXISYM;
  else if (strcmp(type, "PlateFiber") == 0) d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) d = DIM_PSTRESS;

  if (d < 0)
    return NDMaterial::getCopy(type);   // let the base report the unsupported type

  return new LadrunoJ2(this->getTag(), bulk, shear, sig0, Qinf, bIso, Hiso,
                       nBack, Ckin, gKin, rho, d);
}

// ===========================================================================
//  parallel
// ===========================================================================
int LadrunoJ2::sendSelf(int commitTag, Channel& theChannel)
{
  // tag + 7 params + nBack + dim + 2*MAXBACK (C,gamma) + epsP(6) + ebarP(1)
  //     + alpha(MAXBACK*6) + dGamma(1) + cEps22(1)
  static Vector data(1 + 7 + 1 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1 + 1);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = bulk;
  data(c++) = shear;
  data(c++) = sig0;
  data(c++) = Qinf;
  data(c++) = bIso;
  data(c++) = Hiso;
  data(c++) = rho;
  data(c++) = nBack;
  data(c++) = dim;
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  for (int i = 0; i < 6; i++) data(c++) = epsP_n[i];
  data(c++) = ebarP_n;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) data(c++) = alpha_n[k][i];
  data(c++) = dGamma_n;
  data(c++) = cEps22;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoJ2::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static Vector data(1 + 7 + 1 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1 + 1);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  bulk  = data(c++);
  shear = data(c++);
  sig0  = data(c++);
  Qinf  = data(c++);
  bIso  = data(c++);
  Hiso  = data(c++);
  rho   = data(c++);
  nBack = (int)data(c++);
  dim   = (int)data(c++);
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  for (int i = 0; i < 6; i++) epsP_n[i] = data(c++);
  ebarP_n = data(c++);
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = data(c++);
  dGamma_n = data(c++);
  cEps22   = data(c++);

  this->setupDim();             // rebuild vmap/ncomp/condense + resize buffers
  // sync trial to committed
  this->setStateToCommitted();
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  if (condense) strain6[2] = cEps22;
  this->buildElasticTangent(Dtan);
  return 0;
}

// ===========================================================================
//  parameters / print
// ===========================================================================
int LadrunoJ2::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "K") == 0)   return param.addObject(1, this);
  if (strcmp(argv[0], "G") == 0 || strcmp(argv[0], "mu") == 0) return param.addObject(2, this);
  if (strcmp(argv[0], "rho") == 0) return param.addObject(3, this);
  if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "sig0") == 0) return param.addObject(4, this);
  if (strcmp(argv[0], "Hiso") == 0) return param.addObject(5, this);
  if (strcmp(argv[0], "Qinf") == 0) return param.addObject(6, this);
  if (strcmp(argv[0], "b") == 0)    return param.addObject(7, this);
  return -1;
}

int LadrunoJ2::updateParameter(int parameterID, Information& info)
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

int LadrunoJ2::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

void LadrunoJ2::Print(OPS_Stream& s, int flag)
{
  s << endln;
  s << "LadrunoJ2 (combined isotropic + Chaboche AF kinematic J2)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  K, G   : " << bulk << ", " << shear << endln;
  s << "  iso    : sig0=" << sig0 << " Qinf=" << Qinf << " b=" << bIso << " Hiso=" << Hiso << endln;
  s << "  nBack  : " << nBack << endln;
  for (int k = 0; k < nBack; k++)
    s << "    term " << k+1 << ": C=" << Ckin[k] << " gamma=" << gKin[k] << endln;
  s << "  rho    : " << rho << endln;
}

// ===========================================================================
//  recordable responses (stress/strain/tangent + internal state)
// ===========================================================================
Response* LadrunoJ2::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "backStress") == 0 || strcmp(a, "backstress") == 0 || strcmp(a, "alpha") == 0)
    return new MaterialResponse(this, 4, Vector(ncomp));
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 5, Vector(ncomp));
  if (strcmp(a, "equivalentPlasticStrain") == 0 || strcmp(a, "plasticStrainEq") == 0 ||
      strcmp(a, "ebarP") == 0)
    return new MaterialResponse(this, 6, Vector(1));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoJ2::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStress();
      return 0;
    case 2:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();
      return 0;
    case 3:
      if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent();
      return 0;
    case 4:                                   // total backstress, reduced view (stress convention)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) {
          double tot = 0.0;
          for (int k = 0; k < nBack; k++) tot += alpha[k][vmap[a]];
          v(a) = tot;
        }
      }
      return 0;
    case 5:                                   // plastic strain, reduced view (engineering shear)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) {
          int full = vmap[a];
          v(a) = (full >= 3) ? 2.0*epsP[full] : epsP[full];
        }
      }
      return 0;
    case 6:                                   // equivalent (accumulated) plastic strain
      if (matInfo.theVector) (*(matInfo.theVector))(0) = ebarP;
      return 0;
    default:
      return -1;
  }
}
