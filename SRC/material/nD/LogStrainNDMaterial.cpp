/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// LogStrainNDMaterial — logarithmic (Hencky) strain-space finite-strain adaptor.
// See LogStrainNDMaterial.h for the contract and the dSNPO (2008) references.
// Algorithm verified against tests/logstrain_reference.py (numpy oracle).

#include <LogStrainNDMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <Information.h>
#include <Response.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// =========================================================================== //
//  Local tensor helpers (3×3 stored row-major in double[9]; 4th-order [3][3][3][3])
// =========================================================================== //
namespace {

// Voigt order {00,11,22,01,12,20}, engineering shear. Representatives:
const int REP_I[6] = {0, 1, 2, 0, 1, 2};
const int REP_J[6] = {0, 1, 2, 1, 2, 0};

inline int VI(int i, int j) {
  // map (i,j) → Voigt index
  if (i == j) return i;                 // 0,1,2
  if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 3;
  if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 4;
  return 5;                              // (0,2)|(2,0)
}

inline double halfLog(double x)      { return 0.5 * log(x); }
inline double halfLogPrime(double x) { return 0.5 / x; }

void mat3_mul(const double A[9], const double B[9], double C[9]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += A[3*i+k] * B[3*k+j];
      C[3*i+j] = s;
    }
}

void mat3_transpose(const double A[9], double At[9]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) At[3*i+j] = A[3*j+i];
}

double mat3_det(const double A[9]) {
  return A[0]*(A[4]*A[8]-A[5]*A[7])
       - A[1]*(A[3]*A[8]-A[5]*A[6])
       + A[2]*(A[3]*A[7]-A[4]*A[6]);
}

double mat3_inv(const double A[9], double Ai[9]) {
  double det = mat3_det(A);
  double id = (det != 0.0) ? 1.0/det : 0.0;
  Ai[0] =  (A[4]*A[8]-A[5]*A[7])*id;
  Ai[1] = -(A[1]*A[8]-A[2]*A[7])*id;
  Ai[2] =  (A[1]*A[5]-A[2]*A[4])*id;
  Ai[3] = -(A[3]*A[8]-A[5]*A[6])*id;
  Ai[4] =  (A[0]*A[8]-A[2]*A[6])*id;
  Ai[5] = -(A[0]*A[5]-A[2]*A[3])*id;
  Ai[6] =  (A[3]*A[7]-A[4]*A[6])*id;
  Ai[7] = -(A[0]*A[7]-A[1]*A[6])*id;
  Ai[8] =  (A[0]*A[4]-A[1]*A[3])*id;
  return det;
}

// Cyclic Jacobi eigensolver for a symmetric 3×3 (robust at degeneracy).
// eval[a] eigenvalues; evec[3*i + a] = component i of eigenvector a.
void jacobi3(const double Ain[9], double eval[3], double evec[9]) {
  double a[3][3], v[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) { a[i][j] = Ain[3*i+j]; v[i][j] = (i==j)?1.0:0.0; }

  for (int sweep = 0; sweep < 50; sweep++) {
    double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off < 1.0e-300) break;
    const int pq[3][2] = {{0,1},{0,2},{1,2}};
    for (int r = 0; r < 3; r++) {
      int p = pq[r][0], q = pq[r][1];
      if (fabs(a[p][q]) < 1.0e-300) continue;
      double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
      double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta*theta + 1.0));
      double c = 1.0 / sqrt(t*t + 1.0);
      double s = t * c;
      double app = a[p][p], aqq = a[q][q], apq = a[p][q];
      a[p][p] = c*c*app - 2.0*s*c*apq + s*s*aqq;
      a[q][q] = s*s*app + 2.0*s*c*apq + c*c*aqq;
      a[p][q] = a[q][p] = 0.0;
      for (int k = 0; k < 3; k++) {
        if (k != p && k != q) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = a[p][k] = c*akp - s*akq;
          a[k][q] = a[q][k] = s*akp + c*akq;
        }
        double vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c*vkp - s*vkq;
        v[k][q] = s*vkp + c*vkq;
      }
    }
  }
  for (int a_ = 0; a_ < 3; a_++) {
    eval[a_] = a[a_][a_];
    for (int i = 0; i < 3; i++) evec[3*i + a_] = v[i][a_];
  }
}

// Isotropic tensor function Y(X) = Σ y(xᵢ) Eᵢ  (A.51), valid at any multiplicity.
void isoFunction(const double X[9], double (*y)(double), double Y[9]) {
  double val[3], vec[9];
  jacobi3(X, val, vec);
  for (int i = 0; i < 9; i++) Y[i] = 0.0;
  for (int a = 0; a < 3; a++) {
    double ya = y(val[a]);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Y[3*i+j] += ya * vec[3*i+a] * vec[3*j+a];
  }
}

// Derivative dY/dX of an isotropic tensor function (A.52/A.53), as D[i][j][k][l].
// Branches on eigenvalue multiplicity (the D3 degeneracy handling).
void isoFunctionDeriv(const double X[9], double (*y)(double), double (*yp)(double),
                      double D[3][3][3][3]) {
  double val[3], vec[9];
  jacobi3(X, val, vec);

  // eigenprojections E[a][i][j]
  double E[3][3][3];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        E[a][i][j] = vec[3*i+a] * vec[3*j+a];

  // building blocks: I_S and dX²/dX
  double IS[3][3][3][3], dX2[3][3][3][3];
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
      double dik=(i==k), djl=(j==l), dil=(i==l), djk=(j==k);
      IS[i][j][k][l]  = 0.5*(dik*djl + dil*djk);
      dX2[i][j][k][l] = 0.5*(dik*X[3*l+j] + dil*X[3*k+j]
                           + djl*X[3*i+k] + djk*X[3*i+l]);
    }

  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    D[i][j][k][l] = 0.0;

  // multiplicity (relative test, Remark A.2)
  double scale = fabs(val[0]);
  if (fabs(val[1]) > scale) scale = fabs(val[1]);
  if (fabs(val[2]) > scale) scale = fabs(val[2]);
  if (scale < 1.0) scale = 1.0;
  const double rtol = 1.0e-9;
  bool e01 = fabs(val[0]-val[1]) <= rtol*scale;
  bool e12 = fabs(val[1]-val[2]) <= rtol*scale;
  bool e20 = fabs(val[2]-val[0]) <= rtol*scale;

  if (e01 && e12 && e20) {                            // triple: D = y'(x) I_S
    double yp0 = yp(val[0]);
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
      D[i][j][k][l] = yp0 * IS[i][j][k][l];
    return;
  }

  if (e01 || e12 || e20) {                            // two equal (A.52 second / A.53)
    int a;                                            // singleton index
    if (e01) a = 2; else if (e12) a = 0; else a = 1;
    int c = (a+1)%3;                                  // a member of the repeated pair
    double xa = val[a], xc = val[c];
    double ya = y(xa), yc = y(xc), ypa = yp(xa), ypc = yp(xc);
    double d = xa - xc;
    double s1 = (ya - yc)/(d*d) - ypc/d;
    double s2 = 2.0*xc*(ya - yc)/(d*d) - (xa + xc)/d*ypc;
    double s3 = 2.0*(ya - yc)/(d*d*d) - (ypa + ypc)/(d*d);
    double s4 = xc*s3, s5 = xc*s3, s6 = xc*xc*s3;
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
      double dij=(i==j), dkl=(k==l);
      D[i][j][k][l] = s1*dX2[i][j][k][l] - s2*IS[i][j][k][l]
                    - s3*X[3*i+j]*X[3*k+l] + s4*X[3*i+j]*dkl
                    + s5*dij*X[3*k+l]      - s6*dij*dkl;
    }
    return;
  }

  // all distinct (A.52 first branch)
  for (int a = 0; a < 3; a++) {
    int b = (a+1)%3, c = (a+2)%3;
    double xa=val[a], xb=val[b], xc=val[c];
    double coef = y(xa) / ((xa-xb)*(xa-xc));
    double f1 = (xa-xb) + (xa-xc);
    double f2 = (xb-xc);
    double ypa = yp(xa);
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
      double EaEa = E[a][i][j]*E[a][k][l];
      double EbEb = E[b][i][j]*E[b][k][l];
      double EcEc = E[c][i][j]*E[c][k][l];
      D[i][j][k][l] += coef*( dX2[i][j][k][l] - (xb+xc)*IS[i][j][k][l]
                              - f1*EaEa - f2*(EbEb - EcEc) )
                     + ypa*EaEa;
    }
  }
}

} // anonymous namespace

// =========================================================================== //
//  Factory:  nDMaterial LogStrain $tag $innerTag                              //
// =========================================================================== //
void *OPS_LogStrainNDMaterial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial LogStrain $tag $innerTag\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial LogStrain $tag $innerTag\n";
    return 0;
  }
  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial LogStrain " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  return new LogStrainNDMaterial(iData[0], *inner);
}

// =========================================================================== //
//  Construction                                                               //
// =========================================================================== //
void LogStrainNDMaterial::setIdentity(double M[9]) {
  for (int i = 0; i < 9; i++) M[i] = 0.0;
  M[0] = M[4] = M[8] = 1.0;
}

LogStrainNDMaterial::LogStrainNDMaterial(int tag, NDMaterial &inner)
  : FiniteStrainNDMaterial(tag, ND_TAG_LogStrainNDMaterial),
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0)
{
  // an inner 3D (order-6) small-strain material is required
  if (strncmp(inner.getType(), "ThreeDimensional", 80) == 0)
    theMaterial = inner.getCopy();
  else
    theMaterial = inner.getCopy("ThreeDimensional");

  if (theMaterial == 0 || theMaterial->getOrder() != 6) {
    opserr << "LogStrainNDMaterial: inner material must be 3D (order 6)\n";
    exit(-1);
  }
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(Be_trialUpd);
  sigmaCauchy.Zero();
  henckyStrain.Zero();
  aTangent.Zero();
}

LogStrainNDMaterial::LogStrainNDMaterial()
  : FiniteStrainNDMaterial(0, ND_TAG_LogStrainNDMaterial),
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0)
{
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(Be_trialUpd);
}

LogStrainNDMaterial::~LogStrainNDMaterial()
{
  if (theMaterial) delete theMaterial;
}

// =========================================================================== //
//  The finite-strain seam: setTrialF (Box 14.3 i,ii,iv + §14.5 tangent)       //
// =========================================================================== //
int LogStrainNDMaterial::setTrialF(const Matrix &F)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Ftrial9[3*i+j] = F(i, j);

  // (i) F_Δ = F · F_n⁻¹ ;  (ii) Bᵉᵗʳ = F_Δ · Bᵉ_n · F_Δᵀ
  double Fninv[9]; mat3_inv(Fn, Fninv);
  double Fd[9];    mat3_mul(Ftrial9, Fninv, Fd);
  double tmp[9];   mat3_mul(Fd, Be_n, tmp);
  double FdT[9];   mat3_transpose(Fd, FdT);
  double Betr[9];  mat3_mul(tmp, FdT, Betr);
  // enforce symmetry (kill round-off)
  Betr[1] = Betr[3] = 0.5*(Betr[1]+Betr[3]);
  Betr[2] = Betr[6] = 0.5*(Betr[2]+Betr[6]);
  Betr[5] = Betr[7] = 0.5*(Betr[5]+Betr[7]);

  // εᵉᵗʳ = ½ ln Bᵉᵗʳ  (3×3) → engineering-shear Voigt
  double epsT[9]; isoFunction(Betr, halfLog, epsT);
  static Vector epsV(6);
  epsV(0)=epsT[0]; epsV(1)=epsT[4]; epsV(2)=epsT[8];
  epsV(3)=2.0*epsT[1]; epsV(4)=2.0*epsT[5]; epsV(5)=2.0*epsT[2];
  henckyStrain = epsV;

  // (iii) drive the UNCHANGED small-strain material with εᵉᵗʳ
  theMaterial->setTrialStrain(epsV);
  const Vector &tauV = theMaterial->getStress();   // Kirchhoff τ (6)
  const Matrix &D6   = theMaterial->getTangent();  // ∂τ/∂εᵉ (6×6)

  // (iv) Cauchy σ = τ / J
  Jdet = mat3_det(Ftrial9);
  double invJ = (Jdet != 0.0) ? 1.0/Jdet : 0.0;
  for (int a = 0; a < 6; a++) sigmaCauchy(a) = tauV(a) * invJ;

  // bᵉ to commit: elastic inner ⇒ εᵉ = εᵉᵗʳ ⇒ bᵉ = Bᵉᵗʳ. (plastic inner: the
  // seam-3 state protocol is a documented follow-up; see the header.)
  for (int i = 0; i < 9; i++) Be_trialUpd[i] = Betr[i];

  // --- spatial MATERIAL tangent  c = (1/2J)[D : L : B]  (§14.5, eq. 14.99) ---
  // (the geometric/initial-stress term −σ_il δ_jk is the element's K_geo)
  double L[3][3][3][3];
  isoFunctionDeriv(Betr, halfLog, halfLogPrime, L);   // L = ∂(½ln)/∂B
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    L[i][j][k][l] *= 2.0;                              // ∂lnB/∂B = 2·∂(½ln)/∂B

  // B_ijkl = δ_ik Be_jl + δ_jk Be_il   (14.102)
  double Bt[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) {
    double dik=(i==k), djk=(j==k);
    Bt[i][j][k][l] = dik*Betr[3*j+l] + djk*Betr[3*i+l];
  }

  // D tensor from the 6×6 (engineering Voigt) small-strain tangent
  double Dt[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) for(int l=0;l<3;l++)
    Dt[i][j][k][l] = D6(VI(i,j), VI(k,l));

  // DL_ijrs = D_ijpq L_pqrs ;  DLB_ijkl = DL_ijrs B_rskl
  double DL[3][3][3][3];
  for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int r=0;r<3;r++) for(int s=0;s<3;s++) {
    double acc = 0.0;
    for (int p=0;p<3;p++) for(int q=0;q<3;q++) acc += Dt[i][j][p][q]*L[p][q][r][s];
    DL[i][j][r][s] = acc;
  }
  // assemble into 6×6 (material spatial tangent)
  for (int I = 0; I < 6; I++) {
    int i = REP_I[I], j = REP_J[I];
    for (int J = 0; J < 6; J++) {
      int k = REP_I[J], l = REP_J[J];
      double acc = 0.0;
      for (int r=0;r<3;r++) for(int s=0;s<3;s++) acc += DL[i][j][r][s]*Bt[r][s][k][l];
      aTangent(I, J) = invJ * 0.5 * acc;
    }
  }
  return 0;
}

// =========================================================================== //
//  Query                                                                      //
// =========================================================================== //
const Vector &LogStrainNDMaterial::getStress(void)        { return sigmaCauchy; }
const Matrix &LogStrainNDMaterial::getTangent(void)       { return aTangent; }
const Vector &LogStrainNDMaterial::getStrain(void)        { return henckyStrain; }
double        LogStrainNDMaterial::getRho(void)           { return theMaterial->getRho(); }

const Matrix &LogStrainNDMaterial::getInitialTangent(void)
{
  // small-strain elastic tangent is a sound initial spatial tangent (F = I)
  return theMaterial->getInitialTangent();
}

// =========================================================================== //
//  State cycle (the adaptor owns Bᵉ_n and F_n; wraps the inner material)      //
// =========================================================================== //
int LogStrainNDMaterial::commitState(void)
{
  for (int i = 0; i < 9; i++) { Fn[i] = Ftrial9[i]; Be_n[i] = Be_trialUpd[i]; }
  return theMaterial->commitState();
}

int LogStrainNDMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int LogStrainNDMaterial::revertToStart(void)
{
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(Be_trialUpd);
  sigmaCauchy.Zero();
  henckyStrain.Zero();
  aTangent.Zero();
  Jdet = 1.0;
  return theMaterial->revertToStart();
}

// =========================================================================== //
//  Copies                                                                     //
// =========================================================================== //
NDMaterial *LogStrainNDMaterial::getCopy(void)
{
  LogStrainNDMaterial *c = new LogStrainNDMaterial(this->getTag(), *theMaterial);
  for (int i = 0; i < 9; i++) { c->Fn[i] = Fn[i]; c->Be_n[i] = Be_n[i]; }
  return c;
}

NDMaterial *LogStrainNDMaterial::getCopy(const char *type)
{
  if (strncmp(type, "ThreeDimensional", 80) == 0)
    return this->getCopy();
  opserr << "LogStrainNDMaterial::getCopy - only ThreeDimensional is supported\n";
  return 0;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int LogStrainNDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) return -1;
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) { matDbTag = theChannel.getDbTag(); theMaterial->setDbTag(matDbTag); }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) return -1;

  static Vector dataVec(18);
  for (int i = 0; i < 9; i++) { dataVec(i) = Fn[i]; dataVec(9+i) = Be_n[i]; }
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (theMaterial->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int LogStrainNDMaterial::recvSelf(int cTag, Channel &theChannel,
                                  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));

  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(dataID(1));
    if (theMaterial == 0) {
      opserr << "LogStrainNDMaterial::recvSelf - cannot create inner material classTag "
             << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(18);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  for (int i = 0; i < 9; i++) { Fn[i] = dataVec(i); Be_n[i] = dataVec(9+i); }

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) return -4;
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void LogStrainNDMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LogStrainNDMaterial\", ";
    s << "\"inner\": \"" << theMaterial->getTag() << "\"";
    s << "}";
  } else {
    s << "LogStrainNDMaterial (Hencky log-strain finite-strain adaptor), tag: "
      << this->getTag() << endln;
    s << "\tinner material tag: " << theMaterial->getTag() << endln;
  }
}

int LogStrainNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

Response *LogStrainNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  return theMaterial->setResponse(argv, argc, output);
}
