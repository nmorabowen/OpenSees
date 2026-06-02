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
// The pure numerical core lives in LogStrainKernel.h (unit-tested standalone
// against tests/logstrain_reference.py).

#include <LogStrainNDMaterial.h>
#include <LogStrainKernel.h>
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

using namespace logstrain_kernel;

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
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
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
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
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

  // Reject a non-positive Jacobian up front: a negative det F (an inverted /
  // degenerate element under a large explicit step) would otherwise give a
  // sign-flipped Cauchy σ = τ/J and a non-SPD Bᵉᵗʳ whose ½ln has -inf/NaN
  // eigenvalues — silently wrong with no diagnostic. The element (LadrunoBrick
  // -geom finite) treats a <0 return as "cut the step".
  Jdet = mat3_det(Ftrial9);
  if (Jdet <= 0.0) {
    opserr << "LogStrainNDMaterial::setTrialF - non-positive det F (" << Jdet
           << "); element inverted/degenerate\n";
    return -1;
  }

  // (i)+(ii) trial elastic left Cauchy–Green Bᵉᵗʳ and trial Hencky strain εᵉᵗʳ
  double Betr[9];
  trial_Be(Ftrial9, Fn, Be_n, Betr);
  for (int i = 0; i < 9; i++) BeTrial9[i] = Betr[i];   // kept for the §14.5 tangent
  double epsTr6[6], epsN6[6];
  hencky_voigt(Betr, epsTr6);                          // εᵉᵗʳ = ½ ln Bᵉᵗʳ
  hencky_voigt(Be_n, epsN6);                           // εᵉ_n  = ½ ln Bᵉ_n
  for (int k = 0; k < 6; k++) henckyStrain(k) = epsTr6[k];

  // (iii) PLASTIC-INNER PROTOCOL — neutralise the inner's stored εᵖ subtraction.
  // Feed ε_feed = ε_feed_n + (εᵉᵗʳ − εᵉ_n); the inner computes ε_feed − εᵖ_inner
  // = εᵉᵗʳ exactly (the invariant εᵖ_inner = ε_feed_n − εᵉ_n is maintained on
  // commit). Reduces to ε_feed = εᵉᵗʳ for an elastic inner (εᵖ_inner ≡ 0).
  static Vector epsFeedV(6);
  for (int k = 0; k < 6; k++) epsFeedV(k) = epsFeed_n[k] + (epsTr6[k] - epsN6[k]);

  theMaterial->setTrialStrain(epsFeedV);
  const Vector &tauV = theMaterial->getStress();   // Kirchhoff τ (6)
  const Matrix &D6m  = theMaterial->getTangent();  // ∂τ/∂εᵉ (6×6, elastoplastic)
  double tau6[6], D6[36];
  for (int k = 0; k < 6; k++) tau6[k] = tauV(k);
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) D6[6*I+J] = D6m(I, J);

  // (iv) Cauchy σ = τ/J and material spatial tangent c = (1/2J)[D:L:B] at Bᵉᵗʳ
  // (Jdet was computed and checked > 0 at entry)
  double sig6[6], c6[36];
  assemble_material(Betr, D6, tau6, Jdet, sig6, c6);
  for (int k = 0; k < 6; k++) sigmaCauchy(k) = sig6[k];
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) aTangent(I, J) = c6[6*I+J];

  // recover the updated elastic strain εᵉ_{n+1} = Cᵉ : τ (Cᵉ = inner elastic
  // compliance = inv of its initial tangent — exact for a linear elastic inner,
  // for both elastic and plastic steps since τ = Dᵉ:εᵉ always holds), then the
  // committed bᵉ = exp[2 εᵉ_{n+1}]. (v1 assumes a linear-elastic inner law.)
  static Matrix Ce(6, 6);
  Matrix D0(theMaterial->getInitialTangent());
  if (D0.Invert(Ce) < 0) {
    opserr << "LogStrainNDMaterial::setTrialF - inner initial tangent not invertible\n";
    return -1;
  }
  static Vector epsEnp1(6);
  epsEnp1.addMatrixVector(0.0, Ce, tauV, 1.0);     // εᵉ_{n+1} = Cᵉ τ (eng. Voigt)
  double epsEnp16[6];
  for (int k = 0; k < 6; k++) epsEnp16[k] = epsEnp1(k);
  be_from_hencky_voigt(epsEnp16, Be_trialUpd);     // bᵉ to commit

  for (int k = 0; k < 6; k++) epsFeedTrial[k] = epsFeedV(k);  // stage feed
  return 0;
}

// =========================================================================== //
//  Query                                                                      //
// =========================================================================== //
const Vector &LogStrainNDMaterial::getStress(void)  { return sigmaCauchy; }
const Matrix &LogStrainNDMaterial::getTangent(void) { return aTangent; }

// FULL 4th-order spatial constitutive modulus c_ijkl = (1/2J)[D:L:B] — the
// element's consistent-tangent channel (the 6×6 getTangent above is lossy in the
// (k,l) pair). Evaluated at the TRIAL elastic left Cauchy–Green Bᵉᵗʳ (§14.5 uses
// the trial, NOT the returned bᵉ — they differ once plastic), the inner
// small-strain tangent D, and J — all from the most recent setTrialF(). See
// FiniteStrainNDMaterial.h.
int LogStrainNDMaterial::getSpatialTangentTensor(double c[3][3][3][3])
{
  const Matrix &D6m = theMaterial->getTangent();   // ∂τ/∂εᵉ (6×6, minor-sym)
  double D6[36];
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) D6[6 * I + J] = D6m(I, J);
  spatial_tangent_full(BeTrial9, D6, Jdet, c);
  return 0;
}

const Vector &LogStrainNDMaterial::getStrain(void)  { return henckyStrain; }
double        LogStrainNDMaterial::getRho(void)     { return theMaterial->getRho(); }

const Matrix &LogStrainNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();
}

// =========================================================================== //
//  State cycle (the adaptor owns Bᵉ_n and F_n; wraps the inner material)      //
// =========================================================================== //
int LogStrainNDMaterial::commitState(void)
{
  for (int i = 0; i < 9; i++) { Fn[i] = Ftrial9[i]; Be_n[i] = Be_trialUpd[i]; }
  for (int k = 0; k < 6; k++) epsFeed_n[k] = epsFeedTrial[k];   // protocol state
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
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
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
  for (int k = 0; k < 6; k++) c->epsFeed_n[k] = epsFeed_n[k];
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

  static Vector dataVec(24);
  for (int i = 0; i < 9; i++) { dataVec(i) = Fn[i]; dataVec(9+i) = Be_n[i]; }
  for (int k = 0; k < 6; k++) dataVec(18+k) = epsFeed_n[k];
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

  static Vector dataVec(24);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  for (int i = 0; i < 9; i++) { Fn[i] = dataVec(i); Be_n[i] = dataVec(9+i); }
  for (int k = 0; k < 6; k++) epsFeed_n[k] = dataVec(18+k);

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
