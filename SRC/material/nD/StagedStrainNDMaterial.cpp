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
// StagedStrainNDMaterial — small-strain, dimension-general, auto-capturing staged-
// activation wrapper. See StagedStrainNDMaterial.h for the full contract, the
// composability rule (StagedStrain outermost), and the design note
// Ladruno_implementation/staged_deformation_gradiend.md.

#include <StagedStrainNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <Information.h>
#include <Response.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

// =========================================================================== //
//  Factory:  nDMaterial StagedStrain $tag $innerTag <-noInit> <-eps0 e1 e2 …>  //
// =========================================================================== //
void *OPS_StagedStrainNDMaterial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial StagedStrain $tag $innerTag "
              "<-noInit> <-eps0 e1 e2 ...>\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial StagedStrain $tag $innerTag\n";
    return 0;
  }

  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial StagedStrain " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  // GRACEFUL validation (unlike InitStrain's exit(-1)): probe a copy before
  // constructing, so a bad inner returns 0 (command fails) instead of killing the
  // interpreter. We probe the "ThreeDimensional" view because the no-arg getCopy()
  // is "subclass responsibility" (exit(-1)) for the generic material bases; a 3D
  // copy then re-derives any 2D/axisym view via the inherited getCopy(type) switch
  // (so dimension generality is preserved — see getCopy(const char*) below).  // Ladruno
  NDMaterial *probe = inner->getCopy("ThreeDimensional");
  if (probe == 0) {
    opserr << "WARNING nDMaterial StagedStrain " << iData[0]
           << " : inner nDMaterial " << iData[1] << " could not be copied\n";
    return 0;
  }
  delete probe;

  bool   active = true;
  bool   haveEps0 = false;
  Vector eps0in;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-noInit") == 0) {
      active = false;
    } else if (strcmp(opt, "-eps0") == 0) {
      int n = OPS_GetNumRemainingInputArgs();   // remaining args = the eps0 components
      if (n < 1) {
        opserr << "WARNING nDMaterial StagedStrain " << iData[0]
               << " : -eps0 needs at least one component\n";
        return 0;
      }
      double *buf = new double[n];
      if (OPS_GetDoubleInput(&n, buf) != 0) {
        opserr << "WARNING nDMaterial StagedStrain " << iData[0]
               << " : invalid -eps0 components\n";
        delete[] buf;
        return 0;
      }
      eps0in.resize(n);
      for (int i = 0; i < n; i++) eps0in(i) = buf[i];
      delete[] buf;
      haveEps0 = true;
    } else {
      opserr << "WARNING nDMaterial StagedStrain " << iData[0]
             << " : unknown option '" << opt << "'\n";
      return 0;
    }
  }

  return new StagedStrainNDMaterial(iData[0], *inner, active,
                                    haveEps0 ? &eps0in : 0);
}

// =========================================================================== //
//  Construction (buffers are sized LAZILY at the first setTrialStrain, so the   //
//  wrapper adapts to whatever order the element drives — 2D or 3D)             //
// =========================================================================== //
StagedStrainNDMaterial::StagedStrainNDMaterial(int tag, NDMaterial &inner,
                                               bool active_, const Vector *eps0Given)
  : NDMaterial(tag, ND_TAG_StagedStrainNDMaterial),
    theMaterial(0), active(active_), captured(false), eps0Explicit(false),
    eps0(), relStrain(), totalStrain()
{
  // Canonical template = the 3D view (the no-arg getCopy() is exit(-1) for the
  // generic material bases). getCopy(const char*) below re-derives the element's
  // actual view (PlaneStrain/PlaneStress/AxiSymm/3D) from this via the inner's
  // inherited getCopy(type) switch, so one tag serves 2D and 3D elements.  // Ladruno
  theMaterial = inner.getCopy("ThreeDimensional");
  if (theMaterial == 0) {
    opserr << "StagedStrainNDMaterial: failed to copy the inner material\n";
    exit(-1);   // unreachable from script — the factory probes getCopy(3D) first
  }
  if (eps0Given != 0) {            // explicit birth strain: pre-capture
    eps0 = *eps0Given;
    relStrain.resize(eps0.Size());
    totalStrain.resize(eps0.Size());
    relStrain.Zero(); totalStrain.Zero();
    captured     = true;
    eps0Explicit = true;
  }
}

StagedStrainNDMaterial::StagedStrainNDMaterial()
  : NDMaterial(0, ND_TAG_StagedStrainNDMaterial),
    theMaterial(0), active(true), captured(false), eps0Explicit(false),
    eps0(), relStrain(), totalStrain()
{
}

StagedStrainNDMaterial::~StagedStrainNDMaterial()
{
  if (theMaterial) delete theMaterial;
}

// size ε0/relStrain/totalStrain to a given order (preserve explicit ε0 if it fits)
void StagedStrainNDMaterial::sizeBuffers(int n)
{
  if (eps0.Size() != n) {
    Vector keep = eps0;
    eps0.resize(n); eps0.Zero();
    if (eps0Explicit && keep.Size() == n) eps0 = keep;   // honor matching explicit ε0
  }
  if (relStrain.Size()   != n) { relStrain.resize(n);   relStrain.Zero();   }
  if (totalStrain.Size() != n) { totalStrain.resize(n); totalStrain.Zero(); }
}

// =========================================================================== //
//  The strain seam: capture ε0 at birth, forward ε − ε0                        //
// =========================================================================== //
int StagedStrainNDMaterial::setTrialStrain(const Vector &v)
{
  const int n = v.Size();
  if (totalStrain.Size() != n) sizeBuffers(n);
  totalStrain = v;

  if (!active)                     // pass-through (-noInit)
    return theMaterial->setTrialStrain(v);

  if (!captured) {                 // auto-capture: first strain seen IS the birth strain
    eps0 = v;
    captured = true;
  } else if (eps0Explicit && eps0.Size() != n) {
    opserr << "StagedStrainNDMaterial::setTrialStrain - explicit -eps0 size "
           << eps0.Size() << " != element strain order " << n
           << "; re-capturing at birth instead\n";
    eps0 = v;                      // graceful fallback
  }

  relStrain = v;
  relStrain.addVector(1.0, eps0, -1.0);   // ε_rel = ε − ε0
  return theMaterial->setTrialStrain(relStrain);
}

int StagedStrainNDMaterial::setTrialStrain(const Vector &v, const Vector &rate)
{
  const int n = v.Size();
  if (totalStrain.Size() != n) sizeBuffers(n);
  totalStrain = v;

  if (!active)
    return theMaterial->setTrialStrain(v, rate);

  if (!captured) { eps0 = v; captured = true; }
  else if (eps0Explicit && eps0.Size() != n) { eps0 = v; }

  relStrain = v;
  relStrain.addVector(1.0, eps0, -1.0);
  return theMaterial->setTrialStrain(relStrain, rate);
}

int StagedStrainNDMaterial::setTrialStrainIncr(const Vector &v)
{
  // Reconstruct the total strain: the inner holds ε_rel = ε − ε0, so
  // ε_current = inner.getStrain() + ε0, and ε_new = ε_current + Δε.
  const int n = v.Size();
  if (totalStrain.Size() != n) sizeBuffers(n);
  Vector total(n);
  if (captured) {
    total = theMaterial->getStrain();   // ε_rel
    total.addVector(1.0, eps0, 1.0);    // + ε0  → ε_current
  } else {
    total.Zero();
  }
  total.addVector(1.0, v, 1.0);         // + Δε → ε_new
  return setTrialStrain(total);
}

int StagedStrainNDMaterial::setTrialStrainIncr(const Vector &v, const Vector &/*rate*/)
{
  return setTrialStrainIncr(v);
}

// =========================================================================== //
//  Query — delegated to the inner (it already works in ε_rel = ε − ε0)        //
// =========================================================================== //
const Vector &StagedStrainNDMaterial::getStress(void)         { return theMaterial->getStress(); }
const Matrix &StagedStrainNDMaterial::getTangent(void)        { return theMaterial->getTangent(); }
const Matrix &StagedStrainNDMaterial::getInitialTangent(void) { return theMaterial->getInitialTangent(); }
const Vector &StagedStrainNDMaterial::getStrain(void)         { return theMaterial->getStrain(); }  // ε_rel
double        StagedStrainNDMaterial::getRho(void)            { return theMaterial->getRho(); }

const char *StagedStrainNDMaterial::getType(void) const  { return theMaterial->getType(); }
int         StagedStrainNDMaterial::getOrder(void) const { return theMaterial->getOrder(); }

// =========================================================================== //
//  State cycle — ε0 is frozen at birth; the inner owns its own history         //
// =========================================================================== //
int StagedStrainNDMaterial::commitState(void)        { return theMaterial->commitState(); }
int StagedStrainNDMaterial::revertToLastCommit(void) { return theMaterial->revertToLastCommit(); }

int StagedStrainNDMaterial::revertToStart(void)
{
  // Return to the as-constructed state: auto-capture mode re-arms (a genuine
  // restart re-captures birth at the next setTrialStrain); an explicit ε0 is kept.
  if (!eps0Explicit) {
    captured = false;
    eps0.Zero();
  }
  return theMaterial->revertToStart();
}

// =========================================================================== //
//  Copies (getCopy(type) adapts the inner to the requested dimensional view)   //
// =========================================================================== //
// Both getCopy paths ADOPT an already-typed inner directly (via the empty ctor +
// theMaterial assignment) rather than routing through the main ctor, which would
// re-derive the inner to ThreeDimensional and break a 2D element's order.  // Ladruno
NDMaterial *StagedStrainNDMaterial::getCopy(void)
{
  // theMaterial is always a concrete typed subclass here (the main ctor built it
  // as ThreeDimensional), so its no-arg getCopy() is valid and preserves the type.
  NDMaterial *innerCopy = theMaterial->getCopy();
  if (innerCopy == 0) return 0;
  StagedStrainNDMaterial *c = new StagedStrainNDMaterial();
  c->setTag(this->getTag());
  c->theMaterial   = innerCopy;     // adopt — keep the inner's actual type
  c->active        = active;
  c->captured      = captured;
  c->eps0Explicit  = eps0Explicit;
  if (eps0.Size() == innerCopy->getOrder() && eps0.Size() > 0) {
    c->eps0 = eps0; c->sizeBuffers(eps0.Size());
  } else {
    c->captured = false;            // order mismatch ⇒ re-arm auto-capture
  }
  return c;
}

NDMaterial *StagedStrainNDMaterial::getCopy(const char *type)
{
  // Adapt the INNER to the element's requested view (PlaneStrain / PlaneStress /
  // AxiSymmetric / ThreeDimensional / …) via the inner's inherited getCopy(type)
  // switch, then ADOPT it directly so the wrapper carries the element's order.
  NDMaterial *innerCopy = theMaterial->getCopy(type);
  if (innerCopy == 0) {
    opserr << "StagedStrainNDMaterial::getCopy - inner material does not support "
              "type '" << type << "'\n";
    return 0;                       // graceful — no exit(-1)
  }
  StagedStrainNDMaterial *c = new StagedStrainNDMaterial();
  c->setTag(this->getTag());
  c->theMaterial = innerCopy;       // adopt the TYPED inner (do NOT re-derive to 3D)
  c->active      = active;
  // Carry an explicit ε0 only if it matches the requested order; auto-capture
  // copies start uncaptured and snapshot at their own birth.
  if (eps0Explicit && eps0.Size() == innerCopy->getOrder()) {
    c->eps0 = eps0; c->eps0Explicit = true; c->captured = true;
    c->sizeBuffers(eps0.Size());
  }
  return c;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int StagedStrainNDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) return -1;
  int dbTag = this->getDbTag();
  const int n = eps0.Size();

  static ID dataID(4);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) { matDbTag = theChannel.getDbTag(); theMaterial->setDbTag(matDbTag); }
  dataID(2) = matDbTag;
  dataID(3) = n;                    // ε0 order (may be 0 before the first capture)
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) return -1;

  Vector dataVec(3 + n);            // flags + ε0
  dataVec(0) = active       ? 1.0 : 0.0;
  dataVec(1) = captured     ? 1.0 : 0.0;
  dataVec(2) = eps0Explicit ? 1.0 : 0.0;
  for (int i = 0; i < n; i++) dataVec(3 + i) = eps0(i);
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (theMaterial->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int StagedStrainNDMaterial::recvSelf(int cTag, Channel &theChannel,
                                     FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(4);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));
  const int n = dataID(3);

  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(dataID(1));
    if (theMaterial == 0) {
      opserr << "StagedStrainNDMaterial::recvSelf - cannot create inner material classTag "
             << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  Vector dataVec(3 + n);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  active       = (dataVec(0) != 0.0);
  captured     = (dataVec(1) != 0.0);
  eps0Explicit = (dataVec(2) != 0.0);
  if (n > 0) { sizeBuffers(n); for (int i = 0; i < n; i++) eps0(i) = dataVec(3 + i); }

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) return -4;
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void StagedStrainNDMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"StagedStrainNDMaterial\", ";
    s << "\"inner\": \"" << theMaterial->getTag() << "\", ";
    s << "\"active\": " << (active ? "true" : "false") << ", ";
    s << "\"captured\": " << (captured ? "true" : "false");
    s << "}";
  } else {
    s << "StagedStrainNDMaterial (small-strain staged-activation wrapper), tag: "
      << this->getTag() << endln;
    s << "\tinner material tag: " << theMaterial->getTag() << endln;
    s << "\tactive: " << (active ? "yes" : "no")
      << ", eps0 captured: " << (captured ? "yes" : "no") << endln;
  }
}

int StagedStrainNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

int StagedStrainNDMaterial::updateParameter(int /*parameterID*/, Information &/*info*/)
{
  return 0;   // ε0 is not a runtime parameter here (set via -eps0 or auto-capture)
}

// Recorder seam: 'eps0'/'birthStrain' = the captured birth strain; 'totalStrain' =
// the strain from X_ref; everything else (stress, strain = RELATIVE ε−ε0, tangent,
// inner-specific channels) delegates to the inner.  // Ladruno
Response *StagedStrainNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0) {
    if (strcmp(argv[0], "eps0") == 0 || strcmp(argv[0], "birthStrain") == 0)
      return new MaterialResponse(this, 201, eps0);
    if (strcmp(argv[0], "totalStrain") == 0)
      return new MaterialResponse(this, 202, totalStrain);
  }
  return theMaterial->setResponse(argv, argc, output);
}

int StagedStrainNDMaterial::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
  case 201: return matInfo.setVector(eps0);
  case 202: return matInfo.setVector(totalStrain);
  default:  return theMaterial->getResponse(responseID, matInfo);
  }
}

// =========================================================================== //
//  Sensitivity (DDM) — parity with InitStrainNDMaterial                       //
// =========================================================================== //
const Vector &StagedStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int StagedStrainNDMaterial::commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
