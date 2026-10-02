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
#include <MaterialResponse.h>
#include <LadrunoElasticStrainProvider.h>   // Ladruno WP-144 G2: v2 elastic-strain provider mixin
#include <LadrunoMaterialStatus.h>          // Ladruno WP-144 G2 close: LADRUNO_MATERIAL_REFUSED, commit-refusal seam
#include <StagedStrainNDMaterial.h>         // Ladruno WP-144 G2 close: staged-provider construction warning
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

using namespace logstrain_kernel;

// =========================================================================== //
//  Ladruno WP-144 (G2 close): StagedStrain-over-provider construction warning  //
//                                                                             //
//  The wrapper finds the elastic-strain provider by dynamic_cast on its DIRECT //
//  inner. StagedStrain sits between LogStrain and the provider in a staged    //
//  deck (LogStrain -> StagedStrain -> LadrunoNorSand) and does not forward    //
//  the mixin, so the cast fails and the wrapper silently uses inv(D0):tau --   //
//  wrong for a pressure-dependent hyperelastic inner. InitDefGrad is not the   //
//  same case: it wraps LogStrain from OUTSIDE and never sits in between.       //
//  Called from the Tcl/Python factories only (not the ctor, which getCopy runs //
//  per Gauss point), so it prints once per command.                            //
// =========================================================================== //
void ladrunoWarnStagedProviderInner(const char *cmd, int tag, NDMaterial &inner)
{
  StagedStrainNDMaterial *staged = dynamic_cast<StagedStrainNDMaterial *>(&inner);
  if (staged == 0) return;
  if (dynamic_cast<LadrunoElasticStrainProvider *>(&inner) != 0) return;   // not the case today
  if (dynamic_cast<const LadrunoElasticStrainProvider *>(staged->getInner()) == 0) return;
  opserr << "WARNING nDMaterial " << cmd << " " << tag
         << " : the inner is a StagedStrain wrapping a LadrunoElasticStrainProvider material "
            "(e.g. LadrunoNorSand). StagedStrain does not forward the provider, so the wrapper "
            "falls back to the elastic-strain recovery inv(D0):tau, which is WRONG for a "
            "pressure-dependent hyperelastic inner and makes the committed elastic strain drift. "
            "Put the finite-strain wrapper directly on the provider (LogStrain -> LadrunoNorSand) "
            "and stage OUTSIDE it with InitDefGrad.\n";
}

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
  // Validate the inner can produce a 3D (order-6) copy BEFORE constructing, so bad
  // user input fails the command gracefully (return 0) instead of reaching the
  // constructor's hard exit(-1) and killing the interpreter/Python kernel.  // Ladruno (PLUMB-1)
  NDMaterial *probe = (strncmp(inner->getType(), "ThreeDimensional", 80) == 0)
                        ? inner->getCopy() : inner->getCopy("ThreeDimensional");
  if (probe == 0 || probe->getOrder() != 6) {
    opserr << "WARNING nDMaterial LogStrain " << iData[0]
           << " : inner nDMaterial " << iData[1]
           << " must be a 3D (order-6) material\n";
    if (probe != 0) delete probe;
    return 0;
  }
  delete probe;   // the constructor makes its own copy
  ladrunoWarnStagedProviderInner("LogStrain", iData[0], *inner);   // Ladruno WP-144 (G2 close)
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
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0), trialRefused(false)
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
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0), trialRefused(false)
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
  trialRefused = false;                                  // Ladruno WP-144 (G2 close)
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

  // Ladruno WP-144 (G2 close, owner decision 2026-10-02): PROPAGATE the inner's REFUSAL. This return
  // code used to be dropped, so a refusing inner (LadrunoNorSand past its substep cap) reached the
  // element as a SUCCESSFUL setTrialF, the points latched at commit, and the ADR-86b step cut was
  // unreachable on the finite route (-geom linear: -3 at the trial and the smaller retry returns 0;
  // -geom finite: -4 latched, retry -4). ONLY the declared sentinel, not any negative code (ADR-33/34).
  // The inner's trial is frozen at n on a refusal, so nothing below is meaningful: return BEFORE
  // touching sigmaCauchy / aTangent / Be_trialUpd / epsFeedTrial. On a 0 (or any other) return this
  // block is a no-op, so every non-refusing inner is bit-identical.
  if (theMaterial->setTrialStrain(epsFeedV) == LADRUNO_MATERIAL_REFUSED) {
    trialRefused = true;
    return LADRUNO_MATERIAL_REFUSED;     // every element tests `< 0` (LadrunoBrick::updateFinite etc.)
  }
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

  // updated elastic strain εᵉ_{n+1}, then the committed bᵉ = exp[2 εᵉ_{n+1}].
  //  v2 (Ladruno WP-144 G2, owner decision 2026-10-01): an inner that carries its own
  //  elastic strain (LadrunoElasticStrainProvider, e.g. LadrunoNorSand, whose
  //  hyperelastic K, μ are pressure-dependent so τ ≠ D0:εᵉ) PROVIDES εᵉ_{n+1}
  //  directly (engineering Voigt, its trial elastic strain).
  //  v1 (every other inner, UNCHANGED): εᵉ_{n+1} = Cᵉ : τ with Cᵉ = inner elastic
  //  compliance = inv of its initial tangent - exact for a LINEAR elastic inner,
  //  for both elastic and plastic steps since τ = Dᵉ:εᵉ always holds.
  static Vector epsEnp1(6);
  bool haveEpsE = false;
  LadrunoElasticStrainProvider *prov =              // Ladruno WP-144 G2
    dynamic_cast<LadrunoElasticStrainProvider *>(theMaterial);
  if (prov != 0) haveEpsE = prov->ladrunoGetElasticStrain(epsEnp1);
  if (!haveEpsE) {
    static Matrix Ce(6, 6);
    Matrix D0(theMaterial->getInitialTangent());
    if (D0.Invert(Ce) < 0) {
      opserr << "LogStrainNDMaterial::setTrialF - inner initial tangent not invertible\n";
      return -1;
    }
    epsEnp1.addMatrixVector(0.0, Ce, tauV, 1.0);   // εᵉ_{n+1} = Cᵉ τ (eng. Voigt)
  }
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
  // Ladruno WP-144 (G2 close): the INNER commits FIRST, and the wrapper advances bᵉ_n / F_n / the fed
  // strain only if it accepted. The two are independent updates, so on a non-refusing inner this is
  // bit-identical to the old order. A refusing inner commit (LadrunoNorSand latching a refused trial)
  // returns LADRUNO_MATERIAL_REFUSED and has already declared it to Domain::commit()
  // (ladrunoNoteCommitRefusal); advancing bᵉ here would put the wrapper one step ahead of its inner.
  int rc = theMaterial->commitState();
  if (rc == LADRUNO_MATERIAL_REFUSED) return rc;
  if (trialRefused) {
    // the trial was refused and a host committed anyway: the staged state is stale. The inner did not
    // declare it (it returned something else), so declare it here.
    ladrunoNoteCommitRefusal();
    return LADRUNO_MATERIAL_REFUSED;
  }
  for (int i = 0; i < 9; i++) { Fn[i] = Ftrial9[i]; Be_n[i] = Be_trialUpd[i]; }
  for (int k = 0; k < 6; k++) epsFeed_n[k] = epsFeedTrial[k];   // protocol state
  return rc;
}

int LogStrainNDMaterial::revertToLastCommit(void)
{
  trialRefused = false;                                          // Ladruno WP-144 (G2 close)
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
  trialRefused = false;                                          // Ladruno WP-144 (G2 close)
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

// Recorder seam (SEAM-1): the wrapper's getStress()/getStrain()/getTangent() carry
// the FINITE-STRAIN measures (Cauchy σ = τ/J, Hencky ε, spatial c). The inner
// material, queried directly, would report Kirchhoff τ, the fed strain, and the
// small-strain D — wrong by a factor det F under finite deformation. So intercept
// the generic stress/strain/tangent channels here and delegate only the
// inner-specific channels (backStress, plasticStrain, equivalentPlasticStrain,
// damage, …) to the inner.  // Ladruno (SEAM-1)
Response *LogStrainNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0) {
    if (strcmp(argv[0], "stress")  == 0 || strcmp(argv[0], "stresses") == 0 ||
        strcmp(argv[0], "Cauchy")  == 0 || strcmp(argv[0], "cauchy")   == 0)
      return new MaterialResponse(this, 1, sigmaCauchy);    // Cauchy σ (6)

    if (strcmp(argv[0], "strain")  == 0 || strcmp(argv[0], "strains")  == 0)
      return new MaterialResponse(this, 2, henckyStrain);   // Hencky εᵉᵗʳ (6)

    if (strcmp(argv[0], "tangent") == 0)
      return new MaterialResponse(this, 3, aTangent);       // spatial c (6×6)
  }
  return theMaterial->setResponse(argv, argc, output);
}

int LogStrainNDMaterial::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1: return matInfo.setVector(sigmaCauchy);
  case 2: return matInfo.setVector(henckyStrain);
  case 3: return matInfo.setMatrix(aTangent);
  default: return theMaterial->getResponse(responseID, matInfo);
  }
}
