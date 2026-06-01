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
#include <string.h>
#include <math.h>
#include <elementAPI.h>

// ---- static return buffers -------------------------------------------------
Vector LadrunoJ2::stressV(6);
Vector LadrunoJ2::strainV(6);
Matrix LadrunoJ2::tangentM(6, 6);

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
    nBack(0), dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0),
    parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->revertToStart();
}

LadrunoJ2::LadrunoJ2(int tag, double K, double G,
                     double s0, double Qi, double b, double Hi,
                     int nb, const double* C, const double* gam, double r)
  : NDMaterial(tag, ND_TAG_LadrunoJ2),
    bulk(K), shear(G), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), rho(r),
    nBack(nb), dGamma_n(0.0), ebarP_n(0.0), ebarP(0.0), dGammaTrial(0.0),
    parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->revertToStart();
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
  strain6[0] = e(0);
  strain6[1] = e(1);
  strain6[2] = e(2);
  strain6[3] = 0.5 * e(3);   // engineering -> tensor shear
  strain6[4] = 0.5 * e(4);
  strain6[5] = 0.5 * e(5);
  this->integrate();
  return 0;
}

int LadrunoJ2::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoJ2::setTrialStrainIncr(const Vector& v)
{
  static Vector ne(6);
  ne(0) = strain6[0] + v(0);
  ne(1) = strain6[1] + v(1);
  ne(2) = strain6[2] + v(2);
  ne(3) = 2.0*strain6[3] + v(3);
  ne(4) = 2.0*strain6[4] + v(4);
  ne(5) = 2.0*strain6[5] + v(5);
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

  const double tolY = 1.0e-12 * (sig0 > 0.0 ? sig0 : 1.0);

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
  double theta = 0.0, dtheta = 0.0, normMp_dot = 0.0, R = 1.0;
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

    // dR/ddG = (M:Mp)/||M|| - dtheta - (2/3) sig_y'
    double dnormM = dotT(M, Mp)/normM;
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
  for (int i = 0; i < 6; i++) stressV(i) = stress6[i];   // shear = true stress
  return stressV;
}

const Vector& LadrunoJ2::getStrain(void)
{
  strainV(0) = strain6[0];
  strainV(1) = strain6[1];
  strainV(2) = strain6[2];
  strainV(3) = 2.0*strain6[3];   // tensor -> engineering shear
  strainV(4) = 2.0*strain6[4];
  strainV(5) = 2.0*strain6[5];
  return strainV;
}

const Matrix& LadrunoJ2::getTangent(void)
{
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) tangentM(I, J) = Dtan[I][J];
  return tangentM;
}

const Matrix& LadrunoJ2::getInitialTangent(void)
{
  double Kt[6][6];
  this->buildElasticTangent(Kt);
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) tangentM(I, J) = Kt[I][J];
  return tangentM;
}

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
  return 0;
}

int LadrunoJ2::revertToLastCommit(void)
{
  this->setStateToCommitted();
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
  LadrunoJ2* clone = new LadrunoJ2(this->getTag(), bulk, shear, sig0, Qinf, bIso, Hiso,
                                   nBack, Ckin, gKin, rho);
  return clone;
}

NDMaterial* LadrunoJ2::getCopy(const char* type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    return this->getCopy();

  opserr << "LadrunoJ2::getCopy -- type '" << type
         << "' not supported in v1 (ThreeDimensional only)\n";
  return 0;
}

// ===========================================================================
//  parallel
// ===========================================================================
int LadrunoJ2::sendSelf(int commitTag, Channel& theChannel)
{
  // tag + 7 params + nBack + 2*MAXBACK (C,gamma) + epsP(6) + ebarP(1)
  //     + alpha(MAXBACK*6) + dGamma(1)
  static Vector data(1 + 7 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1);
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
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  for (int i = 0; i < 6; i++) data(c++) = epsP_n[i];
  data(c++) = ebarP_n;
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) data(c++) = alpha_n[k][i];
  data(c++) = dGamma_n;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoJ2::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoJ2::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static Vector data(1 + 7 + 1 + 2*MAXBACK + 6 + 1 + MAXBACK*6 + 1);
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
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  for (int i = 0; i < 6; i++) epsP_n[i] = data(c++);
  ebarP_n = data(c++);
  for (int k = 0; k < MAXBACK; k++)
    for (int i = 0; i < 6; i++) alpha_n[k][i] = data(c++);
  dGamma_n = data(c++);

  // sync trial to committed
  this->setStateToCommitted();
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
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
