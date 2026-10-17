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

#include <Ladruno_mutation.h>

#ifndef LADRUNO_MUTATE_CPML
#define LADRUNO_MUTATE_CPML LADRUNO_MUT_NONE
#endif

#include "LadrunoCPML.h"

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Information.h>
#include <Parameter.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <LadrunoElementShell.h>
#include <LadrunoResponseTokens.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace {
constexpr int RESP_FORCE        = 1;
constexpr int RESP_BASIC_FORCE  = 2;
constexpr int RESP_ENERGY       = 3;
constexpr int RESP_STAGE        = 4;
constexpr int RESP_STIFF        = 5;
constexpr int RESP_STIFF_INIT   = 6;
constexpr int RESP_DAMP         = 7;
constexpr int RESP_MASS         = 8;

constexpr int PARAM_STAGE       = 1;
constexpr int PARAM_E           = 2;
constexpr int PARAM_NU          = 3;
constexpr int PARAM_RHO         = 4;
constexpr int PARAM_KF          = 5;
constexpr int PARAM_PORO        = 6;
constexpr int PARAM_ALPHAB      = 7;
constexpr int PARAM_KS          = 8;
constexpr int PARAM_XPERM       = 9;
constexpr int PARAM_YPERM       = 10;
constexpr int PARAM_ZPERM       = 11;

bool resolveTopology(int ndm, int nNodes, int shapeOverride,
                     ladruno_cpml::CPMLTopology& outTopo)
{
  using ladruno_cpml::CPMLTopology;
  if (ndm == 2) {
    if (nNodes == 3) {
      outTopo = CPMLTopology::T3;
      return true;
    }
    if (nNodes == 4) {
      outTopo = CPMLTopology::Q4;
      return true;
    }
    if (nNodes == 6) {
      outTopo = (shapeOverride == 2) ? CPMLTopology::BT6 : CPMLTopology::T6;
      return true;
    }
  } else if (ndm == 3) {
    if (nNodes == 4) {
      outTopo = CPMLTopology::TET4;
      return true;
    }
    if (nNodes == 8) {
      outTopo = CPMLTopology::H8;
      return true;
    }
    if (nNodes == 10) {
      outTopo = (shapeOverride == 4) ? CPMLTopology::BTET10 : CPMLTopology::TET10;
      return true;
    }
    if (nNodes == 20) {
      outTopo = CPMLTopology::H20;
      return true;
    }
  }
  return false;
}

int shapeOverrideFromTopo(ladruno_cpml::CPMLTopology topo)
{
  using ladruno_cpml::CPMLTopology;
  switch (topo) {
    case CPMLTopology::T6:     return 1;
    case CPMLTopology::BT6:    return 2;
    case CPMLTopology::TET10:  return 3;
    case CPMLTopology::BTET10: return 4;
    default:                   return 0;
  }
}
} // namespace

void
LadrunoCPML::ensureBuffers()
{
  ndm_ = ladruno_cpml::shapes::topologyNdm(topo_);
  nNodes_ = ladruno_cpml::shapes::topologyNumNodes(topo_);
  if (biot_.enabled) {
    const int nP = (biot_.pOrder == 1)
                       ? ladruno_cpml::shapes::numVertices(topo_)
                       : nNodes_;
    nDOF_ = nNodes_ * ndm_ + nP;
    fseep_.assign(nP, 0.0);
  } else {
    nDOF_ = nNodes_ * ndm_;
    fseep_.clear();
  }

  if (nodeTags_.Size() != nNodes_) {
    nodeTags_.resize(nNodes_);
    nodeTags_.Zero();
  }
  theNodes_.assign(nNodes_, nullptr);
  coords_.assign(nNodes_ * ndm_, 0.0);
  nodalVols_.assign(nNodes_, 0.0);

  u0_.assign(nDOF_, 0.0);
  f0_.assign(nDOF_, 0.0);
  uTrial_.assign(nDOF_, 0.0);
  vTrial_.assign(nDOF_, 0.0);
  aTrial_.assign(nDOF_, 0.0);

  if (Kmat_.noRows() != nDOF_ || Kmat_.noCols() != nDOF_) {
    Kmat_.resize(nDOF_, nDOF_);
  }
  Kmat_.Zero();

  if (Cmat_.noRows() != nDOF_ || Cmat_.noCols() != nDOF_) {
    Cmat_.resize(nDOF_, nDOF_);
  }
  Cmat_.Zero();

  if (Mmat_.noRows() != nDOF_ || Mmat_.noCols() != nDOF_) {
    Mmat_.resize(nDOF_, nDOF_);
  }
  Mmat_.Zero();

  LadrunoShell::dropKi(Ki_);

  if (resid_.Size() != nDOF_) {
    resid_.resize(nDOF_);
  }
  resid_.Zero();

  if (residInc_.Size() != nDOF_) {
    residInc_.resize(nDOF_);
  }
  residInc_.Zero();

  if (basicForce_.Size() != nDOF_) {
    basicForce_.resize(nDOF_);
  }
  basicForce_.Zero();

  if (Qload_.Size() != nDOF_) {
    Qload_.resize(nDOF_);
  }
  Qload_.Zero();

  if (raScratch_.Size() != nDOF_) {
    raScratch_.resize(nDOF_);
  }
  raScratch_.Zero();

  if (rayForce_.Size() != nDOF_) {
    rayForce_.resize(nDOF_);
  }
  rayForce_.Zero();
}

void
LadrunoCPML::syncGeomParamsFromGeom()
{
  shapeOverride_ = shapeOverrideFromTopo(topo_);
  geomMode_ = (geom_.mode == ladruno_cpml::CPMLGeomMode::Box) ? 0 : 1;
  if (geomMode_ == 0) {
    if (ndm_ == 2) {
      geomParams_.resize(6);
      geomParams_(0) = geom_.box_min[0];
      geomParams_(1) = geom_.box_max[0];
      geomParams_(2) = geom_.box_min[1];
      geomParams_(3) = geom_.box_max[1];
      geomParams_(4) = geom_.L_box[0];
      geomParams_(5) = geom_.L_box[1];
    } else {
      geomParams_.resize(9);
      geomParams_(0) = geom_.box_min[0];
      geomParams_(1) = geom_.box_max[0];
      geomParams_(2) = geom_.box_min[1];
      geomParams_(3) = geom_.box_max[1];
      geomParams_(4) = geom_.box_min[2];
      geomParams_(5) = geom_.box_max[2];
      geomParams_(6) = geom_.L_box[0];
      geomParams_(7) = geom_.L_box[1];
      geomParams_(8) = geom_.L_box[2];
    }
  } else {
    if (ndm_ == 2) {
      geomParams_.resize(5);
      geomParams_(0) = geom_.x0[0];
      geomParams_(1) = geom_.x0[1];
      geomParams_(2) = geom_.n_hat[0];
      geomParams_(3) = geom_.n_hat[1];
      geomParams_(4) = geom_.L_pml;
    } else {
      geomParams_.resize(7);
      geomParams_(0) = geom_.x0[0];
      geomParams_(1) = geom_.x0[1];
      geomParams_(2) = geom_.x0[2];
      geomParams_(3) = geom_.n_hat[0];
      geomParams_(4) = geom_.n_hat[1];
      geomParams_(5) = geom_.n_hat[2];
      geomParams_(6) = geom_.L_pml;
    }
  }
}

void
LadrunoCPML::gatherKinematics(bool useCommitted,
                              double* uOut,
                              double* vOut,
                              double* aOut) const
{
  for (int a = 0; a < nNodes_; ++a) {
    if (theNodes_[a] == nullptr) {
      continue;
    }
    const Vector& disp = useCommitted ? theNodes_[a]->getDisp() : theNodes_[a]->getTrialDisp();
    const Vector& vel  = useCommitted ? theNodes_[a]->getVel()  : theNodes_[a]->getTrialVel();
    const Vector& acc  = useCommitted ? theNodes_[a]->getAccel() : theNodes_[a]->getTrialAccel();
    const int uOff = kernel_.uOffset(a);
    for (int d = 0; d < ndm_; ++d) {
      if (uOut != nullptr) {
        uOut[uOff + d] = (d < disp.Size()) ? disp(d) : 0.0;
      }
      if (vOut != nullptr) {
        vOut[uOff + d] = (d < vel.Size())  ? vel(d)  : 0.0;
      }
      if (aOut != nullptr) {
        aOut[uOff + d] = (d < acc.Size())  ? acc(d)  : 0.0;
      }
    }
  }
  if (biot_.enabled) {
    const int nP = kernel_.numPressureDOFs();
    for (int b = 0; b < nP; ++b) {
      if (theNodes_[b] == nullptr) {
        continue;
      }
      const int pIdx = kernel_.pOffset(b);
      if (pIdx < 0) {
        continue;
      }
      const Vector& disp = useCommitted ? theNodes_[b]->getDisp() : theNodes_[b]->getTrialDisp();
      const Vector& vel  = useCommitted ? theNodes_[b]->getVel()  : theNodes_[b]->getTrialVel();
      if (uOut != nullptr) {
        uOut[pIdx] = (ndm_ < disp.Size()) ? disp(ndm_) : 0.0;
      }
      if (vOut != nullptr) {
        vOut[pIdx] = (ndm_ < vel.Size())  ? vel(ndm_)  : 0.0;
      }
      if (aOut != nullptr) {
        aOut[pIdx] = 0.0;
      }
    }
  }
}

LadrunoCPML::LadrunoCPML(int tag,
                         ladruno_cpml::CPMLTopology topo,
                         const ID& nodeTags,
                         bool lumped,
                         const ladruno_cpml::CPMLMaterialParams& mat,
                         const ladruno_cpml::CPMLGradingParams& grading,
                         const ladruno_cpml::CPMLGeomConfig& geom,
                         int stage,
                         double thick,
                         const ladruno_cpml::CPMLBiotParams& biot)
    : Element(tag, ELE_TAG_LadrunoCPML),
      topo_(topo),
      lumped_(lumped),
      stage_((stage == 0) ? 0 : 1),
      initialStage_((stage == 0) ? 0 : 1),
      shapeOverride_(shapeOverrideFromTopo(topo)),
      geomMode_((geom.mode == ladruno_cpml::CPMLGeomMode::Box) ? 0 : 1),
      thick_(thick),
      charLen_(0.0),
      ndm_(ladruno_cpml::shapes::topologyNdm(topo)),
      nNodes_(ladruno_cpml::shapes::topologyNumNodes(topo)),
      nDOF_(nNodes_ * ndm_),
      valid_(false),
      lastDt_(0.0),
      skipNextUpdate_(false),
      mat_(mat),
      grading_(grading),
      geom_(geom),
      biot_(biot),
      nodeTags_(nNodes_),
      u0Latched_(false),
      Kmat_(nDOF_, nDOF_),
      Cmat_(nDOF_, nDOF_),
      Mmat_(nDOF_, nDOF_),
      Ki_(nullptr),
      resid_(nDOF_),
      residInc_(nDOF_),
      basicForce_(nDOF_),
      Qload_(nDOF_),
      raScratch_(nDOF_),
      rayForce_(nDOF_)
{
  ensureBuffers();
  syncGeomParamsFromGeom();

  if (nodeTags.Size() != nNodes_) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " expected " << nNodes_ << " node tags, got "
           << nodeTags.Size() << endln;
    return;
  }
  for (int i = 0; i < nNodes_; ++i) {
    nodeTags_(i) = nodeTags(i);
  }

  if (ndm_ == 2 && (!(thick_ > 0.0) || !std::isfinite(thick_))) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " requires positive thickness for 2D topologies (got "
           << thick_ << ")" << endln;
    return;
  }

  if (!mat_.isValid() || !grading_.isValid() || !geom_.isValid() ||
      (biot_.enabled && !biot_.isValid(ndm_, topo_))) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " received invalid C-PML material/grading/geometry/Biot parameters" << endln;
    return;
  }

  valid_ = true;
}

LadrunoCPML::LadrunoCPML(int tag,
                         int ndm,
                         const ID& nodeTags,
                         int shapeOverride,
                         double E,
                         double nu,
                         double rho,
                         double thickness,
                         int geomMode,
                         const Vector& geomParams,
                         double mDeg,
                         double R0,
                         double kappaMax,
                         double alphaMax,
                         bool lumped,
                         int initialStage,
                         const Vector& bodyForces)
    : LadrunoCPML(tag,
                  ndm,
                  nodeTags,
                  shapeOverride,
                  E,
                  nu,
                  rho,
                  thickness,
                  geomMode,
                  geomParams,
                  mDeg,
                  R0,
                  kappaMax,
                  alphaMax,
                  lumped,
                  initialStage,
                  bodyForces,
                  ladruno_cpml::CPMLBiotParams{},
                  Vector())
{
}

LadrunoCPML::LadrunoCPML(int tag,
                         int ndm,
                         const ID& nodeTags,
                         int shapeOverride,
                         double E,
                         double nu,
                         double rho,
                         double thickness,
                         int geomMode,
                         const Vector& geomParams,
                         double mDeg,
                         double R0,
                         double kappaMax,
                         double alphaMax,
                         bool lumped,
                         int initialStage,
                         const Vector& bodyForces,
                         const ladruno_cpml::CPMLBiotParams& biot,
                         const Vector& fluidBodyForces)
    : Element(tag, ELE_TAG_LadrunoCPML),
      topo_(ladruno_cpml::CPMLTopology::Q4),
      lumped_(lumped),
      stage_((initialStage == 0) ? 0 : 1),
      initialStage_((initialStage == 0) ? 0 : 1),
      shapeOverride_(shapeOverride),
      geomMode_(geomMode),
      geomParams_(geomParams),
      thick_(thickness),
      charLen_(0.0),
      ndm_(ndm),
      nNodes_(nodeTags.Size()),
      nDOF_(nodeTags.Size() * ndm),
      valid_(false),
      lastDt_(0.0),
      skipNextUpdate_(false),
      mat_(E, nu, rho),
      grading_(mDeg, R0, kappaMax, alphaMax),
      biot_(biot),
      nodeTags_(nodeTags.Size()),
      u0Latched_(false),
      Kmat_(8, 8),
      Cmat_(8, 8),
      Mmat_(8, 8),
      Ki_(nullptr),
      resid_(8),
      residInc_(8),
      basicForce_(8),
      Qload_(8),
      raScratch_(8),
      rayForce_(8)
{
  if (!resolveTopology(ndm, nodeTags.Size(), shapeOverride, topo_)) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " unsupported topology (ndm=" << ndm
           << ", nNodes=" << nodeTags.Size() << ")" << endln;
    ensureBuffers();
    return;
  }

  ensureBuffers();
  for (int i = 0; i < nNodes_; ++i) {
    nodeTags_(i) = nodeTags(i);
  }
  for (int d = 0; d < ndm_ && d < bodyForces.Size() && d < 3; ++d) {
    b_[d] = bodyForces(d);
  }
  for (int d = 0; d < ndm_ && d < fluidBodyForces.Size() && d < 3; ++d) {
    fluidBody_[d] = fluidBodyForces(d);
  }

  if (ndm_ == 2 && (!(thick_ > 0.0) || !std::isfinite(thick_))) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " requires positive thickness for 2D topologies (got "
           << thick_ << ")" << endln;
    return;
  }

  if (geomMode == 0) {
    if (ndm_ == 2 && geomParams.Size() >= 6) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makeBox2D(
          geomParams(0), geomParams(1),
          geomParams(2), geomParams(3),
          geomParams(4), geomParams(5));
    } else if (ndm_ == 3 && geomParams.Size() >= 9) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makeBox3D(
          geomParams(0), geomParams(1),
          geomParams(2), geomParams(3),
          geomParams(4), geomParams(5),
          geomParams(6), geomParams(7), geomParams(8));
    } else {
      return;
    }
  } else if (geomMode == 1) {
    if (ndm_ == 2 && geomParams.Size() >= 5) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makePlane2D(
          geomParams(0), geomParams(1),
          geomParams(2), geomParams(3),
          geomParams(4));
    } else if (ndm_ == 3 && geomParams.Size() >= 7) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makePlane3D(
          geomParams(0), geomParams(1), geomParams(2),
          geomParams(3), geomParams(4), geomParams(5),
          geomParams(6));
    } else {
      return;
    }
  } else {
    return;
  }

  if (!mat_.isValid() || !grading_.isValid() || !geom_.isValid() ||
      (biot_.enabled && !biot_.isValid(ndm_, topo_))) {
    opserr << "WARNING LadrunoCPML::LadrunoCPML - element " << tag
           << " received invalid C-PML material/grading/geometry/Biot parameters" << endln;
    return;
  }

  valid_ = true;
}

LadrunoCPML::LadrunoCPML()
    : Element(0, ELE_TAG_LadrunoCPML),
      topo_(ladruno_cpml::CPMLTopology::Q4),
      lumped_(false),
      stage_(1),
      initialStage_(1),
      shapeOverride_(0),
      geomMode_(0),
      thick_(1.0),
      charLen_(0.0),
      ndm_(2),
      nNodes_(4),
      nDOF_(8),
      valid_(false),
      lastDt_(0.0),
      skipNextUpdate_(false),
      nodeTags_(4),
      u0Latched_(false),
      Kmat_(8, 8),
      Cmat_(8, 8),
      Mmat_(8, 8),
      Ki_(nullptr),
      resid_(8),
      residInc_(8),
      basicForce_(8),
      Qload_(8),
      raScratch_(8),
      rayForce_(8)
{
  ensureBuffers();
}

LadrunoCPML::~LadrunoCPML()
{
  LadrunoShell::dropKi(Ki_);
}

int
LadrunoCPML::getNumExternalNodes() const
{
  return nNodes_;
}

const ID&
LadrunoCPML::getExternalNodes()
{
  return nodeTags_;
}

Node**
LadrunoCPML::getNodePtrs()
{
  return theNodes_.data();
}

int
LadrunoCPML::getNumDOF()
{
  return nDOF_;
}

void
LadrunoCPML::setDomain(Domain* domainPtr)
{
  if (domainPtr == nullptr) {
    std::fill(theNodes_.begin(), theNodes_.end(), nullptr);
    return;
  }

  const int nVert = ladruno_cpml::shapes::numVertices(topo_);
  for (int a = 0; a < nNodes_; ++a) {
    theNodes_[a] = domainPtr->getNode(nodeTags_(a));
    if (theNodes_[a] == nullptr) {
      opserr << "WARNING LadrunoCPML::setDomain - element " << this->getTag()
             << " node " << nodeTags_(a) << " does not exist in Domain" << endln;
      valid_ = false;
      return;
    }
    const bool isCarrier = biot_.enabled && ((biot_.pOrder == 0) || (a < nVert));
    const int expectedNDF = isCarrier ? (ndm_ + 1) : ndm_;
    if (theNodes_[a]->getNumberDOF() != expectedNDF) {
      opserr << "WARNING LadrunoCPML::setDomain - element " << this->getTag()
             << " node " << nodeTags_(a) << " has " << theNodes_[a]->getNumberDOF()
             << " DOFs, expected " << expectedNDF << endln;
      valid_ = false;
      return;
    }
    const Vector& crds = theNodes_[a]->getCrds();
    if (crds.Size() != ndm_) {
      opserr << "WARNING LadrunoCPML::setDomain - element " << this->getTag()
             << " node " << nodeTags_(a) << " has coordinate dimension " << crds.Size()
             << ", expected " << ndm_ << endln;
      valid_ = false;
      return;
    }
    for (int d = 0; d < ndm_; ++d) {
      coords_[a * ndm_ + d] = crds(d);
    }
  }

  const bool kernelOK =
      kernel_.init(topo_, lumped_, coords_.data(), mat_, grading_, geom_, biot_);
  kernel_.setStaticMode(stage_ == 0);
  if (!kernelOK || !kernel_.isValid()) {
    opserr << "WARNING LadrunoCPML::setDomain - element " << this->getTag()
           << " failed C-PML kernel initialization (singular/inverted Jacobian or invalid parameters)"
           << endln;
    valid_ = false;
    return;
  }

  nodalVols_.assign(nNodes_, 0.0);
  kernel_.computeNodalVolumes(nodalVols_.data());
  const double vol = kernel_.computeVolume();
  charLen_ = (vol > 0.0) ? std::pow(vol, 1.0 / static_cast<double>(ndm_)) : 0.0;

  valid_ = true;
  LadrunoShell::dropKi(Ki_);

  this->DomainComponent::setDomain(domainPtr);

  if (kernel_.consumeUnpackedStateFlag()) {
    lastDt_ = kernel_.getLastDt();
    skipNextUpdate_ = true;
  } else {
    this->update();
  }
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0 || Kc != nullptr) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    if (Kc != nullptr) {
      *Kc = this->getTangentStiff();
    }
  }
}

int
LadrunoCPML::commitState()
{
  int retVal = this->Element::commitState();
  if (retVal != 0) {
    opserr << "LadrunoCPML::commitState - failed in base Element::commitState()" << endln;
    return retVal;
  }
  kernel_.commit();
  return 0;
}

int
LadrunoCPML::revertToLastCommit()
{
  kernel_.revertToLastCommit();
  lastDt_ = kernel_.getLastDt();
  skipNextUpdate_ = true;
  return 0;
}

int
LadrunoCPML::revertToStart()
{
  stage_ = initialStage_;
  kernel_.revertToStart();
  kernel_.setStaticMode(stage_ == 0);
  std::fill(u0_.begin(), u0_.end(), 0.0);
  std::fill(f0_.begin(), f0_.end(), 0.0);
  u0Latched_ = false;
  lastDt_ = 0.0;
  Qload_.Zero();
  applyLoad_ = 0;
  appliedB_[0] = 0.0;
  appliedB_[1] = 0.0;
  appliedB_[2] = 0.0;
  appliedFluidBody_[0] = 0.0;
  appliedFluidBody_[1] = 0.0;
  appliedFluidBody_[2] = 0.0;
  resid_.Zero();
  residInc_.Zero();
  basicForce_.Zero();
  rayForce_.Zero();
  LadrunoShell::dropKi(Ki_);
  skipNextUpdate_ = true;
  return 0;
}

int
LadrunoCPML::update()
{
  if (skipNextUpdate_) {
    skipNextUpdate_ = false;
    return 0;
  }

  if (!valid_ || !kernel_.isValid()) {
    return 0;
  }
  for (int a = 0; a < nNodes_; ++a) {
    if (theNodes_[a] == nullptr) {
      return 0;
    }
  }

  double dt = (ops_Dt > 0.0) ? ops_Dt : 0.0;
  if (stage_ == 0) {
    lastDt_ = 0.0;
  } else if (dt > 0.0) {
    lastDt_ = dt;
  } else if (lastDt_ > 0.0) {
    dt = lastDt_;
  }

  gatherKinematics(false, uTrial_.data(), vTrial_.data(), aTrial_.data());

  if (stage_ == 0) {
    kernel_.setStaticMode(true);
    kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
  } else {
    kernel_.setStaticMode(false);
    if (u0Latched_) {
      for (int i = 0; i < nDOF_; ++i) {
        uTrial_[i] -= u0_[i];
      }
    }
    kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), dt);
  }

  return 0;
}

int
LadrunoCPML::setRayleighDampingFactors(double am, double bk, double bk0, double bkc)
{
  index = -1;
  if (Kc != nullptr) {
    delete Kc;
    Kc = nullptr;
  }
  return this->Element::setRayleighDampingFactors(am, bk, bk0, bkc);
}

const Vector&
LadrunoCPML::getRayleighDampingForces()
{
  if (!biot_.enabled) {
    return this->Element::getRayleighDampingForces();
  }
  rayForce_.Zero();
  if (!valid_ || !kernel_.isValid()) {
    return rayForce_;
  }
  if (alphaM == 0.0 && betaK == 0.0 && betaK0 == 0.0 && betaKc == 0.0) {
    return rayForce_;
  }
  const Matrix& cBase = this->Element::getDamp();
  for (int a = 0; a < nNodes_; ++a) {
    const int uBaseA = kernel_.uOffset(a);
    for (int d1 = 0; d1 < ndm_; ++d1) {
      const int r = uBaseA + d1;
      double f_r = 0.0;
      for (int b = 0; b < nNodes_; ++b) {
        const Vector& vNode = theNodes_[b]->getTrialVel();
        const int uBaseB = kernel_.uOffset(b);
        for (int d2 = 0; d2 < ndm_; ++d2) {
          const double vb = (d2 < vNode.Size()) ? vNode(d2) : 0.0;
          f_r += cBase(r, uBaseB + d2) * vb;
        }
      }
      rayForce_(r) = f_r;
    }
  }
  return rayForce_;
}

const Matrix&
LadrunoCPML::getTangentStiff()
{
  Kmat_.Zero();
  if (!valid_ || !kernel_.isValid()) {
    return Kmat_;
  }

  const double scale = scaleFactor();
  const std::vector<double>& kVec =
      (stage_ == 0 || lastDt_ <= 0.0) ? kernel_.getK0() : kernel_.getK_eff();

  for (int i = 0; i < nDOF_; ++i) {
    const size_t rowOff = static_cast<size_t>(i) * nDOF_;
    for (int j = 0; j < nDOF_; ++j) {
      Kmat_(i, j) = scale * kVec[rowOff + j];
    }
  }
  LADRUNO_MUTATE_TANGENT(CPML, Kmat_);
  return Kmat_;
}

const Matrix&
LadrunoCPML::getInitialStiff()
{
  if (Ki_ != nullptr) {
    return *Ki_;
  }

  Kmat_.Zero();
  if (valid_ && kernel_.isValid()) {
    const double scale = scaleFactor();
    const std::vector<double>& k0 = kernel_.getK0();
    for (int i = 0; i < nDOF_; ++i) {
      const size_t rowOff = static_cast<size_t>(i) * nDOF_;
      for (int j = 0; j < nDOF_; ++j) {
        Kmat_(i, j) = scale * k0[rowOff + j];
      }
    }
  }
  LADRUNO_MUTATE_TANGENT(CPML, Kmat_);
  return LadrunoShell::cacheKi(Ki_, Kmat_);
}

const Matrix&
LadrunoCPML::getDamp()
{
  Cmat_.Zero();
  if (!valid_ || !kernel_.isValid()) {
    return Cmat_;
  }

  if (stage_ == 1 && lastDt_ > 0.0) {
    const double scale = scaleFactor();
    const std::vector<double>& cVec = kernel_.getC_eff();
    for (int i = 0; i < nDOF_; ++i) {
      const size_t rowOff = static_cast<size_t>(i) * nDOF_;
      for (int j = 0; j < nDOF_; ++j) {
        Cmat_(i, j) = scale * cVec[rowOff + j];
      }
    }
    LADRUNO_MUTATE_TANGENT(CPML, Cmat_);
  } else if (biot_.enabled) {
    const double scale = scaleFactor();
    const std::vector<double>& c0 = kernel_.getC0();
    for (int i = 0; i < nDOF_; ++i) {
      const size_t rowOff = static_cast<size_t>(i) * nDOF_;
      for (int j = 0; j < nDOF_; ++j) {
        Cmat_(i, j) = scale * c0[rowOff + j];
      }
    }
    LADRUNO_MUTATE_TANGENT(CPML, Cmat_);
  }

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
    const Matrix& cBase = this->Element::getDamp();
    if (!biot_.enabled) {
      Cmat_.addMatrix(1.0, cBase, 1.0);
    } else {
      for (int a = 0; a < nNodes_; ++a) {
        const int uBaseA = kernel_.uOffset(a);
        for (int d1 = 0; d1 < ndm_; ++d1) {
          const int r = uBaseA + d1;
          for (int b = 0; b < nNodes_; ++b) {
            const int uBaseB = kernel_.uOffset(b);
            for (int d2 = 0; d2 < ndm_; ++d2) {
              const int c = uBaseB + d2;
              Cmat_(r, c) += cBase(r, c);
            }
          }
        }
      }
    }
  }
  return Cmat_;
}

const Matrix&
LadrunoCPML::getMass()
{
  Mmat_.Zero();
  if (!valid_ || !kernel_.isValid() || mat_.rho <= 0.0) {
    return Mmat_;
  }

  const double scale = scaleFactor();
  const std::vector<double>& mVec =
      (stage_ == 0) ? kernel_.getM0() : kernel_.getM_eff();

  for (int i = 0; i < nDOF_; ++i) {
    const size_t rowOff = static_cast<size_t>(i) * nDOF_;
    for (int j = 0; j < nDOF_; ++j) {
      Mmat_(i, j) = scale * mVec[rowOff + j];
    }
  }
  return Mmat_;
}

void
LadrunoCPML::zeroLoad()
{
  Qload_.Zero();
  applyLoad_ = 0;
  appliedB_[0] = 0.0;
  appliedB_[1] = 0.0;
  appliedB_[2] = 0.0;
  appliedFluidBody_[0] = 0.0;
  appliedFluidBody_[1] = 0.0;
  appliedFluidBody_[2] = 0.0;
}

int
LadrunoCPML::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  if (theLoad == nullptr) {
    return -1;
  }
  int type = 0;
  const Vector& data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
    applyLoad_ = 1;
    for (int d = 0; d < ndm_ && d < 3; ++d) {
      appliedB_[d] += loadFactor * b_[d];
      appliedFluidBody_[d] += loadFactor * fluidBody_[d];
    }
    return 0;
  } else if (type == LOAD_TAG_SelfWeight) {
    applyLoad_ = 1;
    for (int d = 0; d < ndm_ && d < data.Size() && d < 3; ++d) {
      appliedB_[d] += loadFactor * data(d) * b_[d];
      appliedFluidBody_[d] += loadFactor * data(d) * fluidBody_[d];
    }
    return 0;
  }

  opserr << "LadrunoCPML::addLoad - load type " << type
         << " unknown for element with tag: " << this->getTag() << endln;
  return -1;
}

int
LadrunoCPML::addInertiaLoadToUnbalance(const Vector& accel)
{
  if (mat_.rho <= 0.0 || !valid_ || !kernel_.isValid()) {
    return 0;
  }
  if (!biot_.enabled) {
    return LadrunoShell::addGroundInertia(Qload_,
                                          this->getMass(),
                                          theNodes_.data(),
                                          nNodes_,
                                          ndm_,
                                          lumped_,
                                          accel,
                                          raScratch_,
                                          "LadrunoCPML");
  }

  const int nVert = ladruno_cpml::shapes::numVertices(topo_);
  raScratch_.Zero();
  for (int a = 0; a < nNodes_; ++a) {
    const Vector& Raccel = theNodes_[a]->getRV(accel);
    const bool isCarrier = (biot_.pOrder == 0) || (a < nVert);
    const int need = isCarrier ? (ndm_ + 1) : ndm_;
    if (Raccel.Size() != need) {
      opserr << "LadrunoCPML::addInertiaLoadToUnbalance - element " << this->getTag()
             << ": node " << nodeTags_(a) << " RV size " << Raccel.Size()
             << " != expected " << need << endln;
      return -1;
    }
    const int uOff = kernel_.uOffset(a);
    for (int d = 0; d < ndm_; ++d) {
      raScratch_(uOff + d) = Raccel(d);
    }
  }
  Qload_.addMatrixVector(1.0, this->getMass(), raScratch_, -1.0);
  return 0;
}

const Vector&
LadrunoCPML::getResistingForce()
{
  resid_.Zero();
  if (!valid_ || !kernel_.isValid()) {
    return resid_;
  }

  const double scale = scaleFactor();
  if (stage_ == 0) {
    const std::vector<double>& rDisp =
        biot_.enabled ? kernel_.getR_disp() : kernel_.getR_int();
    for (int i = 0; i < nDOF_; ++i) {
      resid_(i) = scale * rDisp[i];
    }
  } else {
    const std::vector<double>& rDisp = kernel_.getR_disp();
    for (int i = 0; i < nDOF_; ++i) {
      resid_(i) = f0_[i] + scale * rDisp[i];
    }
  }

  LADRUNO_MUTATE_FORCE(CPML, resid_);

  const double* bUse = (applyLoad_ == 0) ? b_ : appliedB_;
  if (!biot_.enabled) {
    if (bUse[0] != 0.0 || bUse[1] != 0.0 || bUse[2] != 0.0) {
      for (int a = 0; a < nNodes_; ++a) {
        const double w = scale * nodalVols_[a];
        for (int d = 0; d < ndm_; ++d) {
          resid_(a * ndm_ + d) -= w * bUse[d];
        }
      }
    }
  } else {
    if (mat_.rho > 0.0 && (bUse[0] != 0.0 || bUse[1] != 0.0 || bUse[2] != 0.0)) {
      for (int a = 0; a < nNodes_; ++a) {
        const double w = scale * nodalVols_[a] * mat_.rho;
        const int uOff = kernel_.uOffset(a);
        for (int d = 0; d < ndm_; ++d) {
          resid_(uOff + d) -= w * bUse[d];
        }
      }
    }
    const double* fbUse = (applyLoad_ == 0) ? fluidBody_ : appliedFluidBody_;
    if (biot_.rhoF > 0.0 && (fbUse[0] != 0.0 || fbUse[1] != 0.0 || fbUse[2] != 0.0)) {
      const int nP = kernel_.numPressureDOFs();
      if (static_cast<int>(fseep_.size()) < nP) {
        fseep_.assign(nP, 0.0);
      }
      kernel_.computeSeepageLoad(fbUse, fseep_.data());
      for (int b = 0; b < nP; ++b) {
        const int pIdx = kernel_.pOffset(b);
        if (pIdx >= 0) {
          resid_(pIdx) -= scale * fseep_[b];
        }
      }
    }
  }

  resid_.addVector(1.0, Qload_, -1.0);
  return resid_;
}

const Vector&
LadrunoCPML::getResistingForceIncInertia()
{
  Vector snap = this->getResistingForce();

  if (!valid_ || !kernel_.isValid()) {
    residInc_ = snap;
    return residInc_;
  }

  const double scale = scaleFactor();
  if ((stage_ == 1 && lastDt_ > 0.0) || biot_.enabled) {
    const std::vector<double>& rVel = kernel_.getR_vel();
    Vector cVelForce(nDOF_);
    for (int i = 0; i < nDOF_; ++i) {
      cVelForce(i) = scale * rVel[i];
    }
    LADRUNO_MUTATE_FORCE(CPML, cVelForce);
    snap.addVector(1.0, cVelForce, 1.0);
  }

  if (mat_.rho > 0.0) {
    if (!biot_.enabled) {
      LadrunoShell::addNodalInertia(snap,
                                    this->getMass(),
                                    theNodes_.data(),
                                    nNodes_,
                                    ndm_,
                                    lumped_,
                                    raScratch_);
    } else {
      raScratch_.Zero();
      for (int a = 0; a < nNodes_; ++a) {
        const Vector& accel = theNodes_[a]->getTrialAccel();
        const int uOff = kernel_.uOffset(a);
        for (int d = 0; d < ndm_; ++d) {
          raScratch_(uOff + d) = (d < accel.Size()) ? accel(d) : 0.0;
        }
      }
      snap.addMatrixVector(1.0, this->getMass(), raScratch_, 1.0);
    }
  }

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
    snap.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }

  residInc_ = snap;
  return residInc_;
}

int
LadrunoCPML::sendSelf(int commitTag, Channel& theChannel)
{
  const int myDbTag = this->getDbTag();
  const int nKernelState = kernel_.numStateDoubles();

  // Single fixed-size ID (36 ints) to avoid FileDatastore (type, size, dbTag, commitTag) clobbering
  static ID idata(36);
  idata(0) = this->getTag();
  idata(1) = ndm_;
  idata(2) = nNodes_;
  idata(3) = shapeOverride_;
  idata(4) = geomMode_;
  idata(5) = lumped_ ? 1 : 0;
  idata(6) = stage_;
  idata(7) = u0Latched_ ? 1 : 0;
  idata(8) = applyLoad_;
  idata(9) = initialStage_;
  idata(10) = nKernelState;
  for (int i = 0; i < 20; ++i) {
    idata(11 + i) = (i < nNodes_) ? nodeTags_(i) : -1;
  }
  idata(31) = biot_.enabled ? 1 : 0;
  idata(32) = biot_.pOrder;
  idata(33) = biot_.formulation;
  idata(34) = biot_.stabMode;
  idata(35) = 0;

  if (theChannel.sendID(myDbTag, commitTag, idata) < 0) {
    opserr << "LadrunoCPML::sendSelf - failed to send ID data" << endln;
    return -1;
  }

  // Single Vector of size 60 + 3 * nDOF_ + nKernelState
  const int totalVecSize = 60 + 3 * nDOF_ + nKernelState;
  Vector vdata(totalVecSize);
  vdata.Zero();

  vdata(0)  = mat_.E;
  vdata(1)  = mat_.nu;
  vdata(2)  = mat_.rho;
  vdata(3)  = thick_;
  vdata(4)  = grading_.m_deg;
  vdata(5)  = grading_.R0;
  vdata(6)  = grading_.kappa_max;
  vdata(7)  = grading_.alpha_max;
  vdata(8)  = charLen_;
  vdata(9)  = lastDt_;
  vdata(10) = static_cast<double>(geomParams_.Size());
  for (int i = 0; i < 9; ++i) {
    vdata(11 + i) = (i < geomParams_.Size()) ? geomParams_(i) : 0.0;
  }
  for (int i = 0; i < 3; ++i) {
    vdata(20 + i) = b_[i];
    vdata(23 + i) = appliedB_[i];
  }
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      vdata(26 + r * 3 + c) = geom_.R[r][c];
    }
  }
  vdata(35) = geom_.x0[0];
  vdata(36) = geom_.x0[1];
  vdata(37) = geom_.x0[2];
  vdata(38) = geom_.L_pml;
  vdata(39) = static_cast<double>(nDOF_);
  vdata(40) = alphaM;
  vdata(41) = betaK;
  vdata(42) = betaK0;
  vdata(43) = betaKc;
  vdata(44) = biot_.Kf;
  vdata(45) = biot_.poro;
  vdata(46) = biot_.rhoF;
  vdata(47) = biot_.biotAlpha;
  vdata(48) = biot_.Ks;
  vdata(49) = biot_.stabValue;
  for (int i = 0; i < 3; ++i) {
    vdata(50 + i) = biot_.kbar[i];
    vdata(53 + i) = fluidBody_[i];
    vdata(56 + i) = appliedFluidBody_[i];
  }
  vdata(59) = 0.0;

  int offset = 60;
  for (int i = 0; i < nDOF_; ++i) {
    vdata(offset + i) = u0_[i];
    vdata(offset + nDOF_ + i) = f0_[i];
    vdata(offset + 2 * nDOF_ + i) = (i < Qload_.Size()) ? Qload_(i) : 0.0;
  }
  offset += 3 * nDOF_;

  if (nKernelState > 0) {
    std::vector<double> kbuf(nKernelState, 0.0);
    kernel_.packState(kbuf.data());
    for (int i = 0; i < nKernelState; ++i) {
      vdata(offset + i) = kbuf[i];
    }
  }

  if (theChannel.sendVector(myDbTag, commitTag, vdata) < 0) {
    opserr << "LadrunoCPML::sendSelf - failed to send Vector data" << endln;
    return -1;
  }

  return 0;
}

int
LadrunoCPML::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  (void)theBroker;
  const int myDbTag = this->getDbTag();

  static ID idata(36);
  if (theChannel.recvID(myDbTag, commitTag, idata) < 0) {
    opserr << "LadrunoCPML::recvSelf - failed to receive ID data" << endln;
    return -1;
  }

  this->setTag(idata(0));
  const int recvNdm = idata(1);
  const int recvNNodes = idata(2);
  shapeOverride_ = idata(3);
  geomMode_ = idata(4);
  lumped_ = (idata(5) != 0);
  stage_ = (idata(6) == 0) ? 0 : 1;
  const bool recvU0Latched = (idata(7) != 0);
  const int recvApplyLoad = idata(8);
  initialStage_ = (idata(9) == 0) ? 0 : 1;
  const int nKernelState = idata(10);

  biot_ = ladruno_cpml::CPMLBiotParams{};
  biot_.enabled = (idata(31) != 0);
  biot_.pOrder = idata(32);
  biot_.formulation = idata(33);
  biot_.stabMode = idata(34);

  if (!resolveTopology(recvNdm, recvNNodes, shapeOverride_, topo_)) {
    opserr << "LadrunoCPML::recvSelf - invalid topology in ID data" << endln;
    return -1;
  }

  ensureBuffers();
  for (int i = 0; i < nNodes_; ++i) {
    nodeTags_(i) = idata(11 + i);
  }

  const int totalVecSize = 60 + 3 * nDOF_ + nKernelState;
  Vector vdata(totalVecSize);
  if (theChannel.recvVector(myDbTag, commitTag, vdata) < 0) {
    opserr << "LadrunoCPML::recvSelf - failed to receive Vector data" << endln;
    return -1;
  }

  const double E_in   = vdata(0);
  const double nu_in  = vdata(1);
  const double rho_in = vdata(2);
  thick_              = vdata(3);
  grading_.m_deg      = vdata(4);
  grading_.R0         = vdata(5);
  grading_.kappa_max  = vdata(6);
  grading_.alpha_max  = vdata(7);
  charLen_            = vdata(8);
  const double recvLastDt = vdata(9);
  const int gpSize    = static_cast<int>(std::lround(vdata(10)));

  mat_.init(E_in, nu_in, rho_in);

  if (gpSize > 0 && gpSize <= 9) {
    geomParams_.resize(gpSize);
    for (int i = 0; i < gpSize; ++i) {
      geomParams_(i) = vdata(11 + i);
    }
  }

  for (int i = 0; i < 3; ++i) {
    b_[i]        = vdata(20 + i);
    appliedB_[i] = vdata(23 + i);
  }
  applyLoad_ = recvApplyLoad;

  if (geomMode_ == 0) {
    if (ndm_ == 2 && geomParams_.Size() >= 6) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makeBox2D(
          geomParams_(0), geomParams_(1),
          geomParams_(2), geomParams_(3),
          geomParams_(4), geomParams_(5));
    } else if (ndm_ == 3 && geomParams_.Size() >= 9) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makeBox3D(
          geomParams_(0), geomParams_(1),
          geomParams_(2), geomParams_(3),
          geomParams_(4), geomParams_(5),
          geomParams_(6), geomParams_(7), geomParams_(8));
    }
  } else {
    if (ndm_ == 2 && geomParams_.Size() >= 5) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makePlane2D(
          geomParams_(0), geomParams_(1),
          geomParams_(2), geomParams_(3),
          geomParams_(4));
    } else if (ndm_ == 3 && geomParams_.Size() >= 7) {
      geom_ = ladruno_cpml::CPMLGeomConfig::makePlane3D(
          geomParams_(0), geomParams_(1), geomParams_(2),
          geomParams_(3), geomParams_(4), geomParams_(5),
          geomParams_(6));
    }
  }

  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      geom_.R[r][c] = vdata(26 + r * 3 + c);
    }
  }
  geom_.x0[0] = vdata(35);
  geom_.x0[1] = vdata(36);
  geom_.x0[2] = vdata(37);
  geom_.L_pml = vdata(38);

  alphaM = vdata(40);
  betaK  = vdata(41);
  betaK0 = vdata(42);
  betaKc = vdata(43);

  biot_.Kf        = vdata(44);
  biot_.poro      = vdata(45);
  biot_.rhoF      = vdata(46);
  biot_.biotAlpha = vdata(47);
  biot_.Ks        = vdata(48);
  biot_.stabValue = vdata(49);
  for (int i = 0; i < 3; ++i) {
    biot_.kbar[i]          = vdata(50 + i);
    fluidBody_[i]          = vdata(53 + i);
    appliedFluidBody_[i]   = vdata(56 + i);
  }

  int offset = 60;
  for (int i = 0; i < nDOF_; ++i) {
    u0_[i]     = vdata(offset + i);
    f0_[i]     = vdata(offset + nDOF_ + i);
    Qload_(i)  = vdata(offset + 2 * nDOF_ + i);
  }
  offset += 3 * nDOF_;

  u0Latched_ = recvU0Latched;
  lastDt_    = recvLastDt;

  std::vector<double> kbuf(nKernelState, 0.0);
  for (int i = 0; i < nKernelState; ++i) {
    kbuf[i] = vdata(offset + i);
  }
  if (nKernelState > 0) {
    kernel_.setPendingState(kbuf.data(), nKernelState);
  }

  valid_ = mat_.isValid() && grading_.isValid() && geom_.isValid() &&
           (!biot_.enabled || biot_.isValid(ndm_, topo_)) &&
           (ndm_ != 2 || (thick_ > 0.0 && std::isfinite(thick_)));

  if (this->getDomain() != nullptr) {
    const Vector savedQ = Qload_;
    const std::vector<double> savedU0 = u0_;
    const std::vector<double> savedF0 = f0_;
    this->setDomain(this->getDomain());
    Qload_     = savedQ;
    u0_        = savedU0;
    f0_        = savedF0;
    u0Latched_ = recvU0Latched;
    lastDt_    = recvLastDt;
    applyLoad_ = recvApplyLoad;
    for (int i = 0; i < 3; ++i) {
      appliedB_[i]         = vdata(23 + i);
      appliedFluidBody_[i] = vdata(56 + i);
    }
  }

  skipNextUpdate_ = true;
  return 0;
}

void
LadrunoCPML::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "Element: " << this->getTag() << " Type: LadrunoCPML"
      << " Nodes: " << nodeTags_
      << " E: " << mat_.E << " nu: " << mat_.nu << " rho: " << mat_.rho
      << " stage: " << stage_ << endln;
  } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{\"name\": \"" << this->getTag()
      << "\", \"type\": \"LadrunoCPML\", \"nodes\": [";
    for (int i = 0; i < nNodes_; ++i) {
      if (i > 0) s << ", ";
      s << nodeTags_(i);
    }
    s << "], \"E\": " << mat_.E
      << ", \"nu\": " << mat_.nu
      << ", \"rho\": " << mat_.rho
      << ", \"stage\": " << stage_ << "}";
  }
}

Response*
LadrunoCPML::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  if (argc < 1 || argv == nullptr || argv[0] == nullptr) {
    return nullptr;
  }

  Response* theResponse = nullptr;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoCPML");
  output.attr("eleTag", this->getTag());
  for (int i = 0; i < nNodes_; ++i) {
    char nodeAttr[16];
    std::snprintf(nodeAttr, sizeof(nodeAttr), "node%d", i + 1);
    output.attr(nodeAttr, nodeTags_(i));
  }

  if (LadrunoResp::is(argv[0], "force")) {
    theResponse = new ElementResponse(this, RESP_FORCE, resid_);
  } else if (LadrunoResp::is(argv[0], "basicForce")) {
    theResponse = new ElementResponse(this, RESP_BASIC_FORCE, basicForce_);
  } else if (LadrunoResp::is(argv[0], "energy") ||
              std::strcmp(argv[0], "power") == 0 ||
              std::strcmp(argv[0], "externalPower") == 0) {
    theResponse = new ElementResponse(this, RESP_ENERGY, 0.0);
  } else if (std::strcmp(argv[0], "stage") == 0) {
    theResponse = new ElementResponse(this, RESP_STAGE, 0);
  } else if (LadrunoResp::is(argv[0], "stiff")) {
    theResponse = new ElementResponse(this, RESP_STIFF, Kmat_);
  } else if (LadrunoResp::is(argv[0], "stiffInitial")) {
    theResponse = new ElementResponse(this, RESP_STIFF_INIT, Kmat_);
  } else if (std::strcmp(argv[0], "damping") == 0 ||
             std::strcmp(argv[0], "damp") == 0) {
    theResponse = new ElementResponse(this, RESP_DAMP, Cmat_);
  } else if (std::strcmp(argv[0], "mass") == 0) {
    theResponse = new ElementResponse(this, RESP_MASS, Mmat_);
  }

  return LadrunoShell::finishResponse(this, theResponse, argv, argc, output);
}

int
LadrunoCPML::getResponse(int responseID, Information& eleInfo)
{
  switch (responseID) {
    case RESP_FORCE:
      return eleInfo.setVector(this->getResistingForce());
    case RESP_BASIC_FORCE: {
      basicForce_.Zero();
      if (valid_ && kernel_.isValid()) {
        const double scale = scaleFactor();
        const std::vector<double>& rInt = kernel_.getR_int();
        for (int i = 0; i < nDOF_; ++i) {
          basicForce_(i) = (stage_ == 1 ? f0_[i] : 0.0) + scale * rInt[i];
        }
      }
      return eleInfo.setVector(basicForce_);
    }
    case RESP_ENERGY:
      return eleInfo.setDouble(scaleFactor() * kernel_.getP_ext());
    case RESP_STAGE:
      return eleInfo.setInt(stage_);
    case RESP_STIFF:
      return eleInfo.setMatrix(this->getTangentStiff());
    case RESP_STIFF_INIT:
      return eleInfo.setMatrix(this->getInitialStiff());
    case RESP_DAMP:
      return eleInfo.setMatrix(this->getDamp());
    case RESP_MASS:
      return eleInfo.setMatrix(this->getMass());
    case 222222:
      return eleInfo.setVector(this->getRayleighDampingForces());
    case 444444: {
      Vector inertial = this->getResistingForceIncInertia();
      inertial -= this->getRayleighDampingForces();
      inertial -= this->getResistingForce();
      return eleInfo.setVector(inertial);
    }
    default:
      return this->Element::getResponse(responseID, eleInfo);
  }
}

int
LadrunoCPML::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1 || argv == nullptr || argv[0] == nullptr) {
    return -1;
  }

  if (argc >= 2) {
    if (argv[1] == nullptr || std::atoi(argv[1]) != this->getTag()) {
      return -1;
    }
  }

  if (std::strcmp(argv[0], "stage") == 0 ||
      std::strcmp(argv[0], "updateMaterialStage") == 0) {
    return param.addObject(PARAM_STAGE, this);
  }
  if (std::strcmp(argv[0], "E") == 0) {
    return param.addObject(PARAM_E, this);
  }
  if (std::strcmp(argv[0], "nu") == 0 || std::strcmp(argv[0], "v") == 0) {
    return param.addObject(PARAM_NU, this);
  }
  if (std::strcmp(argv[0], "rho") == 0) {
    return param.addObject(PARAM_RHO, this);
  }
  if (biot_.enabled) {
    if (std::strcmp(argv[0], "Kf") == 0) {
      return param.addObject(PARAM_KF, this);
    }
    if (std::strcmp(argv[0], "poro") == 0 || std::strcmp(argv[0], "n") == 0) {
      return param.addObject(PARAM_PORO, this);
    }
    if (std::strcmp(argv[0], "alphaB") == 0 || std::strcmp(argv[0], "biotAlpha") == 0) {
      return param.addObject(PARAM_ALPHAB, this);
    }
    if (std::strcmp(argv[0], "Ks") == 0) {
      return param.addObject(PARAM_KS, this);
    }
    if (std::strcmp(argv[0], "xPerm") == 0 || std::strcmp(argv[0], "permX") == 0) {
      return param.addObject(PARAM_XPERM, this);
    }
    if (std::strcmp(argv[0], "yPerm") == 0 || std::strcmp(argv[0], "permY") == 0) {
      return param.addObject(PARAM_YPERM, this);
    }
    if (std::strcmp(argv[0], "zPerm") == 0 || std::strcmp(argv[0], "permZ") == 0) {
      return (ndm_ == 3) ? param.addObject(PARAM_ZPERM, this) : -1;
    }
  }

  return -1;
}

int
LadrunoCPML::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case PARAM_STAGE: {
      const int newStage = (static_cast<int>(std::lround(info.theDouble)) == 0) ? 0 : 1;
      if (stage_ == 0 && newStage == 1) {
        // Latch current committed nodal displacements/pressures u_0 and elastic force f_0 = scale * K_0 * u_0
        std::fill(u0_.begin(), u0_.end(), 0.0);
        gatherKinematics(true, u0_.data(), nullptr, nullptr);
        const double scale = scaleFactor();
        const std::vector<double>& k0 = kernel_.getK0();
        for (int i = 0; i < nDOF_; ++i) {
          double fi = 0.0;
          const size_t rowOff = static_cast<size_t>(i) * nDOF_;
          for (int j = 0; j < nDOF_; ++j) {
            fi += k0[rowOff + j] * u0_[j];
          }
          f0_[i] = scale * fi;
        }
        u0Latched_ = true;
        stage_ = 1;
        lastDt_ = 0.0;
        kernel_.revertToStart();
        kernel_.setStaticMode(false);
      } else if (newStage == 0) {
        stage_ = 0;
        std::fill(u0_.begin(), u0_.end(), 0.0);
        std::fill(f0_.begin(), f0_.end(), 0.0);
        u0Latched_ = false;
        lastDt_ = 0.0;
        kernel_.revertToStart();
        kernel_.setStaticMode(true);
        gatherKinematics(true, uTrial_.data(), vTrial_.data(), aTrial_.data());
        kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
      } else {
        stage_ = 1;
        kernel_.setStaticMode(false);
      }
      return 0;
    }
    case PARAM_E: {
      if (!(info.theDouble > 0.0) || !std::isfinite(info.theDouble)) {
        return -1;
      }
      mat_.init(info.theDouble, mat_.nu, mat_.rho);
      kernel_.updateMaterialParams(mat_, lumped_);
      if (stage_ == 0 || lastDt_ <= 0.0) {
        kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
      }
      return 0;
    }
    case PARAM_NU: {
      if (!(info.theDouble > -1.0 && info.theDouble < 0.499999) || !std::isfinite(info.theDouble)) {
        return -1;
      }
      mat_.init(mat_.E, info.theDouble, mat_.rho);
      kernel_.updateMaterialParams(mat_, lumped_);
      if (stage_ == 0 || lastDt_ <= 0.0) {
        kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
      }
      return 0;
    }
    case PARAM_RHO: {
      if (!(info.theDouble > 0.0) || !std::isfinite(info.theDouble)) {
        return -1;
      }
      mat_.init(mat_.E, mat_.nu, info.theDouble);
      kernel_.updateMaterialParams(mat_, lumped_);
      if (stage_ == 0 || lastDt_ <= 0.0) {
        kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
      }
      return 0;
    }
    case PARAM_KF:
    case PARAM_PORO:
    case PARAM_ALPHAB:
    case PARAM_KS:
    case PARAM_XPERM:
    case PARAM_YPERM:
    case PARAM_ZPERM: {
      if (!biot_.enabled) {
        return -1;
      }
      ladruno_cpml::CPMLBiotParams nextBiot = biot_;
      if (parameterID == PARAM_KF)          nextBiot.Kf = info.theDouble;
      else if (parameterID == PARAM_PORO)   nextBiot.poro = info.theDouble;
      else if (parameterID == PARAM_ALPHAB) nextBiot.biotAlpha = info.theDouble;
      else if (parameterID == PARAM_KS)     nextBiot.Ks = info.theDouble;
      else if (parameterID == PARAM_XPERM)  nextBiot.kbar[0] = info.theDouble;
      else if (parameterID == PARAM_YPERM)  nextBiot.kbar[1] = info.theDouble;
      else if (parameterID == PARAM_ZPERM)  nextBiot.kbar[2] = info.theDouble;

      if (!nextBiot.isValid(ndm_, topo_)) {
        return -1;
      }
      biot_ = nextBiot;
      kernel_.updateBiotParams(biot_);
      LadrunoShell::dropKi(Ki_);
      if (stage_ == 0 || lastDt_ <= 0.0) {
        kernel_.setTrial(uTrial_.data(), vTrial_.data(), aTrial_.data(), 0.0);
      }
      return 0;
    }
    default:
      return -1;
  }
}
