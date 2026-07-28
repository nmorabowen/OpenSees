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
// Created: 07/2026
//
// LadrunoUP — unified Biot u-p saturated-porous continuum element (ADR 71 P1).
// See LadrunoUP.h for the honest-p contract table and
// Ladruno_implementation/71_ladruno_up_family_adr.md §3.1–3.5 for the theory.
//
// Family idiom donors (read/adapt/cite): LadrunoQuad.{h,cpp} (member layout,
// setDomain, update/commit cycle, setResponse, charLength), FourNodeQuadUP.cpp
// + SSPquadUP.cpp (the upstream u-p elements — what we deliberately do
// DIFFERENTLY is documented at each API below).

#include <LadrunoUP.h>
#include <LadrunoUPKernel.h>
#include <LadrunoUPShapes.h>
#include <SolidTransformation.h>        // ADR 78 — geometry-method seam
#include <SolidTransformationLinear.h>  // ADR 78 — identity method / safe fallback
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <ElementIter.h>   // stage-flip sibling broadcast (P4)
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <LadrunoResponseTokens.h>   // Ladruno — shared recorder-token aliases
#include <ElementalLoad.h>
#include <Response.h>
#include <classTags.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

// The canonical #define lands in classTags.h with WP1.C (parser/registration
// owner). This guarded fallback keeps the two work packages independently
// compilable; identical-value redefinition is a no-op once WP1.C merges.

// hard caps of the shape-provider family (BTET10: 10 nodes, H8: 8 GP)
static const int UP_MAX_NODES = 10;
static const int UP_MAX_DOF   = 40;   // 10 nodes x (3 + 1)

// ===========================================================================
//  ctors / dtor
// ===========================================================================

LadrunoUP::LadrunoUP(int tag, int ndm, const ID &nodeTags,
                     NDMaterial &theMaterial,
                     double thickness,
                     double Kf, double poro, double rhoF,
                     const Vector &perm,
                     double biotAlpha, double Ks,
                     const Vector &body, const Vector &fluidBody,
                     int formulation,
                     int pOrder,
                     bool lumped,
                     int stabMode, double stabValue,
                     bool dynSeepage,
                     int geomMethod)
 : Element(tag, ELE_TAG_LadrunoUP),
   ndm_(ndm), nNodes_(nodeTags.Size()), shapeKind_(SHAPE_NONE),
   nGP_(0), nU_(0), nP_(0), nStr_(0), nDof_(0),
   connectedExternalNodes_(nodeTags.Size()), theMaterials_(0),
   thickness_(thickness), Kf_(Kf), poro_(poro), rhoF_(rhoF),
   biotAlpha_(biotAlpha), Ks_(Ks),
   formulation_(formulation), pOrder_(pOrder), lumpedMass_(lumped),
   stabMode_(stabMode), stabValue_(stabValue), dynSeepage_(dynSeepage),
   oneOverQbar_(0.0),
   theGeom_(0), corotActive_(false),
   applyLoad_(0), Qload_(),
   kInitSolidValid_(false), kcSolidValid_(false),
   k0Valid_(false), k0Dirty_(true), stabAlpha_(0.0), stabDirty_(true),
   geomValid_(false), elemH_(0.0), charLen_(0.0)
{
  for (int n = 0; n < UP_MAX_NODES; n++)
    theNodes_[n] = 0;
  for (int i = 0; i < 3; i++) {
    kbar_[i] = 0.0;
    b_[i] = 0.0;
    fluidBody_[i] = 0.0;
    appliedB_[i] = 0.0;
    appliedFluidBody_[i] = 0.0;
  }

  for (int n = 0; n < nNodes_; n++)
    connectedExternalNodes_(n) = nodeTags(n);

  // ---- geometry method (ADR 78) --------------------------------------------
  // Parser polices -geom (corot is 3D-only); this is the belt-and-braces guard
  // for direct-construction callers. Unknown id / 2D corot degrade to linear
  // loudly rather than running a planar polar into det(H)=0.
  if (geomMethod == SolidTransformation::METHOD_COROT && ndm_ != 3) {
    opserr << "LadrunoUP::LadrunoUP - element " << tag
           << ": -geom corot is 3D-only (planar node cloud has no polar "
              "rotation); using -geom linear\n";
    geomMethod = SolidTransformation::METHOD_LINEAR;
  }
  theGeom_ = SolidTransformation::create(geomMethod);
  if (theGeom_ == 0) {
    opserr << "LadrunoUP::LadrunoUP - element " << tag
           << ": unknown geometry-method id " << geomMethod
           << "; using -geom linear\n";
    theGeom_ = new SolidTransformationLinear();
  }
  corotActive_ = (theGeom_->getMethodID() == SolidTransformation::METHOD_COROT);

  // ---- parameter intake ------------------------------------------------------
  for (int i = 0; i < ndm_; i++) {
    if (perm.Size() > i)       kbar_[i]      = perm(i);
    if (body.Size() > i)       b_[i]         = body(i);
    if (fluidBody.Size() > i)  fluidBody_[i] = fluidBody(i);
  }
  if (ndm_ == 3)
    thickness_ = 1.0;   // pinned: thickness is 2D-only

  // storage coefficient 1/Q̄ = n/K_f + (α−n)/K_s (kernel; parser polices domain,
  // element guards the one input that would produce inf and poison S silently)
  if (Kf_ > 0.0) {
    oneOverQbar_ = ladruno_up::oneOverQbar(poro_, Kf_, biotAlpha_, Ks_);
  } else {
    opserr << "LadrunoUP::LadrunoUP - element " << tag
           << ": -Kf must be > 0 (got " << Kf_ << "); storage set to 0\n";
    oneOverQbar_ = 0.0;
  }

  // ---- provider dispatch + dofMap + size-derived members + buffers ----------
  // (shared with recvSelf — see configureSizing; on a bad (ndm,nNodes) it
  // leaves shapeKind_ == SHAPE_NONE and we bail before touching materials,
  // exactly as the pre-refactor inline block did — WP1.A pin.)
  this->configureSizing();
  if (shapeKind_ == SHAPE_NONE)
    return;   // configureSizing printed the reason; setDomain refuses loudly too

  // ---- per-GP material copies (effective-stress seam) -----------------------
  const char *matType = (ndm_ == 2) ? "PlaneStrain" : "ThreeDimensional";
  theMaterials_ = new NDMaterial *[nGP_];
  for (int i = 0; i < nGP_; i++) {
    theMaterials_[i] = theMaterial.getCopy(matType);
    if (theMaterials_[i] == 0) {
      opserr << "LadrunoUP::LadrunoUP - element " << tag
             << ": failed to get a " << matType << " copy of material "
             << theMaterial.getTag() << "\n";
      exit(-1);   // family idiom (LadrunoQuad); parser pre-validates
    }
  }
}

LadrunoUP::LadrunoUP()
 : Element(0, ELE_TAG_LadrunoUP),
   ndm_(0), nNodes_(0), shapeKind_(SHAPE_NONE),
   nGP_(0), nU_(0), nP_(0), nStr_(0), nDof_(0),
   connectedExternalNodes_(0), theMaterials_(0),
   thickness_(0.0), Kf_(0.0), poro_(0.0), rhoF_(0.0),
   biotAlpha_(1.0), Ks_(0.0),
   formulation_(0), pOrder_(0), lumpedMass_(false),
   stabMode_(0), stabValue_(0.0), dynSeepage_(true),
   oneOverQbar_(0.0),
   theGeom_(new SolidTransformationLinear()), corotActive_(false),
   applyLoad_(0), Qload_(),
   kInitSolidValid_(false), kcSolidValid_(false),
   k0Valid_(false), k0Dirty_(true), stabAlpha_(0.0), stabDirty_(true),
   geomValid_(false), elemH_(0.0), charLen_(0.0)
{
  for (int n = 0; n < UP_MAX_NODES; n++)
    theNodes_[n] = 0;
  for (int i = 0; i < 3; i++) {
    kbar_[i] = 0.0;
    b_[i] = 0.0;
    fluidBody_[i] = 0.0;
    appliedB_[i] = 0.0;
    appliedFluidBody_[i] = 0.0;
  }
}

LadrunoUP::~LadrunoUP()
{
  if (theMaterials_) {
    for (int i = 0; i < nGP_; i++)
      if (theMaterials_[i])
        delete theMaterials_[i];
    delete[] theMaterials_;
  }
  if (theGeom_)
    delete theGeom_;
}

// ===========================================================================
//  construction / reconstruction sizing (ctor + recvSelf share this)
// ===========================================================================

// Provider dispatch by (ndm_, nNodes_), dofMap, every size-derived member and
// the full per-instance algebra/cache allocation — the block the ctor used to
// inline. Serialization (recvSelf) reconstructs a polymorphic element (the
// shape family varies with ndm/nNodes), so it MUST reproduce this sizing; a
// hand-copy in recvSelf would be the exact drift-bug class ADR 71 §3.5 logs
// against upstream (BrickUP/SSPbrickUP never (re)initialize massType). One
// helper, two callers, zero drift. Precondition: ndm_, nNodes_, pOrder_ and the
// physical params (Kf_/poro_/… → oneOverQbar_ is computed by the CALLER) are
// already set. On an illegal (ndm_, nNodes_): prints, leaves SHAPE_NONE, bails
// before allocating (caller returns without building materials).
void LadrunoUP::configureSizing(void)
{
  if      (ndm_ == 2 && nNodes_ == 3)  shapeKind_ = SHAPE_T3;
  else if (ndm_ == 2 && nNodes_ == 4)  shapeKind_ = SHAPE_Q4;
  else if (ndm_ == 2 && nNodes_ == 6)  shapeKind_ = SHAPE_BT6;
  else if (ndm_ == 3 && nNodes_ == 8)  shapeKind_ = SHAPE_H8;
  else if (ndm_ == 3 && nNodes_ == 10) shapeKind_ = SHAPE_BTET10;
  else {
    opserr << "LadrunoUP::configureSizing - element " << this->getTag()
           << ": no shape provider for (ndm=" << ndm_ << ", nNodes=" << nNodes_
           << "); legal pairs are (2,3) (2,4) (2,6) (3,8) (3,10)\n";
    shapeKind_ = SHAPE_NONE;
    return;   // caller (ctor/recvSelf) bails; setDomain refuses loudly too
  }

  nGP_  = this->shapeNGP();
  nU_   = nNodes_ * ndm_;
  nStr_ = (ndm_ == 2) ? 3 : 6;

  // pressure carriers: equal-order = every node; TH = the provider's vertex set
  int pCarrier[UP_MAX_NODES];
  const int *carrTH = this->shapePCarrierTH();
  for (int n = 0; n < nNodes_; n++)
    pCarrier[n] = (pOrder_ == 1) ? carrTH[n] : 1;

  uOff_.assign(nNodes_, 0);
  pOff_.assign(nNodes_, -1);
  nDof_ = ladruno_up::dofMap(nNodes_, ndm_, pCarrier, uOff_.data(), pOff_.data());
  carrierNodes_.clear();
  for (int n = 0; n < nNodes_; n++)
    if (pCarrier[n])
      carrierNodes_.push_back(n);
  nP_ = (int)carrierNodes_.size();

  // ---- size the per-instance algebra (NO class statics — ADR-40) ------------
  Qload_.resize(nDof_);            Qload_.Zero();
  K_.resize(nDof_, nDof_);         K_.Zero();
  damp_.resize(nDof_, nDof_);      damp_.Zero();
  mass_.resize(nDof_, nDof_);      mass_.Zero();
  resid_.resize(nDof_);            resid_.Zero();
  rayForce_.resize(nDof_);         rayForce_.Zero();
  K0_.resize(nDof_, nDof_);        K0_.Zero();
  KinitSolid_.resize(nU_, nU_);    KinitSolid_.Zero();
  KcSolid_.resize(nU_, nU_);       KcSolid_.Zero();
  Bscratch_.resize(nStr_, nU_);    Bscratch_.Zero();
  KsolScratch_.resize(nU_, nU_);   KsolScratch_.Zero();
  CRayScratch_.resize(nU_, nU_);   CRayScratch_.Zero();
  MsolScratch_.resize(nU_, nU_);   MsolScratch_.Zero();
  epsScratch_.resize(nStr_);       epsScratch_.Zero();
  uSolScratch_.resize(nU_);        uSolScratch_.Zero();
  fuScratch_.resize(nU_);          fuScratch_.Zero();
  pScratch_.resize(nP_);           pScratch_.Zero();

  xy_.assign(nNodes_ * ndm_, 0.0);
  NuAll_.assign(nGP_ * nNodes_, 0.0);
  dNuAll_.assign(nGP_ * nNodes_ * ndm_, 0.0);
  NpAll_.assign(nGP_ * nP_, 0.0);
  dNpAll_.assign(nGP_ * nP_ * ndm_, 0.0);
  dv_.assign(nGP_, 0.0);
  dNuBar_.assign(nNodes_ * ndm_, 0.0);
  Qblk_.assign(nU_ * nP_, 0.0);
  Hblk_.assign(nP_ * nP_, 0.0);
  Sblk_.assign(nP_ * nP_, 0.0);
  HtUnit_.assign(nP_ * nP_, 0.0);

  // corot scratch (ADR 78) — cheap; sized unconditionally so recvSelf-rebuilt
  // elements are ready regardless of which method the incoming stream carries.
  refCrdsScratch_.resize(nNodes_, 3);  refCrdsScratch_.Zero();
  curCrdsScratch_.resize(nNodes_, 3);  curCrdsScratch_.Zero();
  Rcur_.resize(3, 3);                  Rcur_.Zero();
  Rcur_(0, 0) = Rcur_(1, 1) = Rcur_(2, 2) = 1.0;
  zeroFu_.resize(nU_);                 zeroFu_.Zero();
  QrotBlk_.assign(nU_ * nP_, 0.0);
  uDcur_.resize(nU_);                  uDcur_.Zero();
  uDn_.resize(nU_);                    uDn_.Zero();
}

// ===========================================================================
//  shape-provider dispatch (stateless P0 structs, switch — WP1.A pin)
// ===========================================================================

int LadrunoUP::shapeNGP(void) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     return ladruno_up::shapes::T3::nGP;
    case SHAPE_Q4:     return ladruno_up::shapes::Q4::nGP;
    case SHAPE_H8:     return ladruno_up::shapes::H8::nGP;
    case SHAPE_BT6:    return ladruno_up::shapes::BT6::nGP;
    case SHAPE_BTET10: return ladruno_up::shapes::BTET10::nGP;
    default:           return 0;
  }
}

void LadrunoUP::shapeEval(int gp, double *Nu, double *dNu, double *detJ) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     ladruno_up::shapes::T3::eval(gp, xy_.data(), Nu, dNu, detJ); break;
    case SHAPE_Q4:     ladruno_up::shapes::Q4::eval(gp, xy_.data(), Nu, dNu, detJ); break;
    case SHAPE_H8:     ladruno_up::shapes::H8::eval(gp, xy_.data(), Nu, dNu, detJ); break;
    case SHAPE_BT6:    ladruno_up::shapes::BT6::eval(gp, xy_.data(), Nu, dNu, detJ); break;
    case SHAPE_BTET10: ladruno_up::shapes::BTET10::eval(gp, xy_.data(), Nu, dNu, detJ); break;
    default:           *detJ = 0.0; break;
  }
}

void LadrunoUP::shapeEvalP(int gp, bool taylorHood, double *Np, double *dNp) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     ladruno_up::shapes::T3::evalP(gp, xy_.data(), taylorHood, Np, dNp); break;
    case SHAPE_Q4:     ladruno_up::shapes::Q4::evalP(gp, xy_.data(), taylorHood, Np, dNp); break;
    case SHAPE_H8:     ladruno_up::shapes::H8::evalP(gp, xy_.data(), taylorHood, Np, dNp); break;
    case SHAPE_BT6:    ladruno_up::shapes::BT6::evalP(gp, xy_.data(), taylorHood, Np, dNp); break;
    case SHAPE_BTET10: ladruno_up::shapes::BTET10::evalP(gp, xy_.data(), taylorHood, Np, dNp); break;
    default: break;
  }
}

double LadrunoUP::shapeElemSizeH(void) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     return ladruno_up::shapes::T3::elemSizeH(xy_.data());
    case SHAPE_Q4:     return ladruno_up::shapes::Q4::elemSizeH(xy_.data());
    case SHAPE_H8:     return ladruno_up::shapes::H8::elemSizeH(xy_.data());
    case SHAPE_BT6:    return ladruno_up::shapes::BT6::elemSizeH(xy_.data());
    case SHAPE_BTET10: return ladruno_up::shapes::BTET10::elemSizeH(xy_.data());
    default:           return 0.0;
  }
}

const int *LadrunoUP::shapePCarrierTH(void) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     return ladruno_up::shapes::T3::pCarrierTH;
    case SHAPE_Q4:     return ladruno_up::shapes::Q4::pCarrierTH;
    case SHAPE_H8:     return ladruno_up::shapes::H8::pCarrierTH;
    case SHAPE_BT6:    return ladruno_up::shapes::BT6::pCarrierTH;
    case SHAPE_BTET10: return ladruno_up::shapes::BTET10::pCarrierTH;
    default:           return 0;
  }
}

double LadrunoUP::shapeGPWeight(int gp) const
{
  switch (shapeKind_) {
    case SHAPE_T3:     return ladruno_up::shapes::T3::gpWeight[gp];
    case SHAPE_Q4:     return ladruno_up::shapes::Q4::gpWeight[gp];
    case SHAPE_H8:     return ladruno_up::shapes::H8::gpWeight[gp];
    case SHAPE_BT6:    return ladruno_up::shapes::BT6::gpWeight[gp];
    case SHAPE_BTET10: return ladruno_up::shapes::BTET10::gpWeight[gp];
    default:           return 0.0;
  }
}

// ===========================================================================
//  topology
// ===========================================================================

int LadrunoUP::getNumExternalNodes(void) const { return nNodes_; }
const ID &LadrunoUP::getExternalNodes(void)    { return connectedExternalNodes_; }
Node **LadrunoUP::getNodePtrs(void)            { return theNodes_; }
int LadrunoUP::getNumDOF(void)                 { return nDof_; }

void LadrunoUP::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    for (int n = 0; n < UP_MAX_NODES; n++)
      theNodes_[n] = 0;
    geomValid_ = false;
    return;
  }

  if (shapeKind_ == SHAPE_NONE) {
    opserr << "LadrunoUP::setDomain - element " << this->getTag()
           << ": no shape provider selected (bad ndm/nNodes at construction); "
              "element not activated\n";
    return;
  }

  for (int n = 0; n < nNodes_; n++) {
    theNodes_[n] = theDomain->getNode(connectedExternalNodes_(n));
    if (theNodes_[n] == 0) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": node " << connectedExternalNodes_(n) << " does not exist\n";
      return;
    }
  }

  // ---- heterogeneous-ndf validation (WP3.A pin 2, STRICT) -------------------
  // Carrier nodes (dofMap gave them a p slot, pOff_>=0) need ndf = ndm+1;
  // non-carrier mid-edge nodes need ndf = ndm exactly. Equal-order lanes are
  // all-carrier so this collapses to the P1 ndf==ndm+1 check; Taylor–Hood is
  // the mixed case. A uniform-ndf TH model would leave floating mid-edge p-DOFs
  // → silent singularity — the fork rejects loudly, naming node + expectation.
  for (int n = 0; n < nNodes_; n++) {
    bool carrier = (pOff_[n] >= 0);
    int need = carrier ? (ndm_ + 1) : ndm_;
    if (theNodes_[n]->getNumberDOF() != need) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": node " << connectedExternalNodes_(n) << " has ndf = "
             << theNodes_[n]->getNumberDOF() << " but this "
             << (carrier ? "vertex/carrier" : "mid-edge")
             << " node needs ndf = " << need << " (u" << ndm_
             << (carrier ? " + p)" : ", no p)")
             << "; element not activated\n";
      return;
    }
  }

  this->DomainComponent::setDomain(theDomain);

  // ---- geometry caches (linear geometry: fixed for the analysis) ------------
  for (int a = 0; a < nNodes_; a++) {
    const Vector &crd = theNodes_[a]->getCrds();
    if (crd.Size() < ndm_) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": node " << connectedExternalNodes_(a) << " has "
             << crd.Size() << " coordinates, need " << ndm_ << "\n";
      return;
    }
    for (int i = 0; i < ndm_; i++)
      xy_[a * ndm_ + i] = crd(i);
  }

  // ---- straight-side guard (WP3.A pin 4, Bézier shapes, both pOrders) --------
  // The BT6/BTET10 providers assume an AFFINE geometry map (∇L_i exact
  // constants) — the Taylor–Hood inf-sup precondition (ADR §3.3). That holds
  // iff every mid-edge node sits at its edge midpoint. Enforce it here (the
  // providers are silent on violation). Edge→mid-edge pairing is DERIVED from
  // the Bernstein cross terms in LadrunoUPShapes.h:
  //   BT6:    N3=2ξ₁ξ₃→node3 on edge(0,1); N4=2ξ₁ξ₂→node4 on edge(1,2);
  //           N5=2ξ₂ξ₃→node5 on edge(2,0)   (node0↔ξ₃,node1↔ξ₁,node2↔ξ₂)
  //   BTET10: N4=2L1L2→node4 on (0,1); N5=2L2L3→node5 on (1,2);
  //           N6=2L1L3→node6 on (0,2); N7=2L1L4→node7 on (0,3);
  //           N8=2L3L4→node8 on (2,3); N9=2L2L4→node9 on (1,3)
  //           (node k ↔ L_{k+1}, k=0..3)
  if (shapeKind_ == SHAPE_BT6 || shapeKind_ == SHAPE_BTET10) {
    static const int bt6Edges[3][3]     = {{3,0,1},{4,1,2},{5,2,0}};
    static const int btet10Edges[6][3]  = {{4,0,1},{5,1,2},{6,0,2},
                                           {7,0,3},{8,2,3},{9,1,3}};
    const int (*edges)[3] = (shapeKind_ == SHAPE_BT6) ? bt6Edges : btet10Edges;
    int nEdge = (shapeKind_ == SHAPE_BT6) ? 3 : 6;
    for (int e = 0; e < nEdge; e++) {
      int m = edges[e][0], a = edges[e][1], b = edges[e][2];
      double d2 = 0.0, len2 = 0.0;
      for (int i = 0; i < ndm_; i++) {
        double mid = 0.5 * (xy_[a * ndm_ + i] + xy_[b * ndm_ + i]);
        double dm  = xy_[m * ndm_ + i] - mid;
        double ab  = xy_[b * ndm_ + i] - xy_[a * ndm_ + i];
        d2   += dm * dm;
        len2 += ab * ab;
      }
      double dist = sqrt(d2), edgeLen = sqrt(len2);
      if (dist > 1.0e-6 * edgeLen) {
        opserr << "LadrunoUP::setDomain - element " << this->getTag()
               << ": mid-edge node " << connectedExternalNodes_(m)
               << " is off the midpoint of edge ("
               << connectedExternalNodes_(a) << ", " << connectedExternalNodes_(b)
               << ") by " << dist << " (tol " << 1.0e-6 * edgeLen
               << " = 1e-6 * edge length " << edgeLen << "). The Bézier "
                  "Taylor–Hood providers require STRAIGHT sides (affine map); "
                  "element not activated\n";
        geomValid_ = false;
        return;
      }
    }
  }

  // ---- Jacobian sign / winding adjudication (WP3.A pin 5) --------------------
  // The providers are silent on singular Jacobians (kernel idiom) — the element
  // is the mandated guard. T3/Q4/H8/BT6: standard CCW/right-handed winding is
  // map-positive, so keep the P1 detJ>0-per-GP rejection. BTET10 is the
  // exception: the pinned node↔L map (detJ = det[v0−v3,v1−v3,v2−v3] =
  // −det[v1−v0,v2−v0,v3−v0]) makes conventionally RIGHT-handed tets evaluate
  // detJ < 0. We accept BOTH windings — the signed-J gradients (dNuAll_/dNpAll_)
  // are orientation-correct either way; only the integration MEASURE must be
  // positive, so an all-negative tet folds dv = w·|detJ|. Mixed sign or a
  // near-zero detJ (degenerate/curved) is a loud error.
  double detJall[16];
  for (int gp = 0; gp < nGP_; gp++) {
    double detJ = 0.0;
    this->shapeEval(gp, &NuAll_[gp * nNodes_], &dNuAll_[gp * nNodes_ * ndm_], &detJ);
    detJall[gp] = detJ;
  }

  double windingSign = 1.0;   // multiply detJ by this to get |detJ|
  if (shapeKind_ == SHAPE_BTET10) {
    double maxAbs = 0.0;
    for (int gp = 0; gp < nGP_; gp++)
      if (fabs(detJall[gp]) > maxAbs) maxAbs = fabs(detJall[gp]);
    if (maxAbs <= 0.0) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": degenerate BTET10 (all Gauss-point Jacobians ~ 0); "
                "element not activated\n";
      geomValid_ = false;
      return;
    }
    double eps = 1.0e-10 * maxAbs;
    int nPos = 0, nNeg = 0, nZero = 0;
    for (int gp = 0; gp < nGP_; gp++) {
      if      (detJall[gp] >  eps) nPos++;
      else if (detJall[gp] < -eps) nNeg++;
      else                         nZero++;
    }
    if (nZero > 0 || (nPos > 0 && nNeg > 0)) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": BTET10 Jacobian is mixed-sign or near-zero across Gauss "
                "points (" << nPos << " positive, " << nNeg << " negative, "
             << nZero << " near-zero) — degenerate or curved geometry; "
                "element not activated\n";
      geomValid_ = false;
      return;
    }
    if (nNeg == nGP_) {
      windingSign = -1.0;   // all-negative: fold the measure to |detJ|
      static bool windingNoteShown = false;
      if (!windingNoteShown) {
        windingNoteShown = true;
        opserr << "LadrunoUP: BTET10 element(s) supplied in the conventionally "
                  "right-handed winding evaluate detJ < 0 under the pinned "
                  "node↔barycentric map — folding the integration measure to "
                  "|detJ| (gradients are orientation-correct; printed once)\n";
      }
    }
  } else {
    for (int gp = 0; gp < nGP_; gp++) {
      if (detJall[gp] <= 0.0) {
        opserr << "LadrunoUP::setDomain - element " << this->getTag()
               << ": non-positive Jacobian detJ = " << detJall[gp] << " at GP "
               << gp << " (bad node ordering or degenerate geometry); element "
                  "not activated\n";
        geomValid_ = false;
        return;
      }
    }
  }

  double measure = 0.0;   // area (2D) | volume (3D), thickness NOT folded
  for (int gp = 0; gp < nGP_; gp++) {
    double absDetJ = windingSign * detJall[gp];   // = |detJ| (>0 by the guards)
    double w = this->shapeGPWeight(gp);
    measure += w * absDetJ;
    dv_[gp] = w * absDetJ * ((ndm_ == 2) ? thickness_ : 1.0);
    this->shapeEvalP(gp, (pOrder_ == 1), &NpAll_[gp * nP_], &dNpAll_[gp * nP_ * ndm_]);
  }

  // element-mean gradients (mean-dilatation B-bar; also the B-bar-consistent Q)
  double vol = 0.0;
  for (int k = 0; k < nNodes_ * ndm_; k++)
    dNuBar_[k] = 0.0;
  for (int gp = 0; gp < nGP_; gp++) {
    vol += dv_[gp];
    for (int k = 0; k < nNodes_ * ndm_; k++)
      dNuBar_[k] += dNuAll_[gp * nNodes_ * ndm_ + k] * dv_[gp];
  }
  if (vol > 0.0)
    for (int k = 0; k < nNodes_ * ndm_; k++)
      dNuBar_[k] /= vol;

  elemH_ = this->shapeElemSizeH();
  charLen_ = (ndm_ == 2) ? sqrt(measure) : cbrt(measure);

  geomValid_ = true;
  this->buildStaticBlocks();

  // corot (ADR 78): initialize the frame + the committed deformational state
  // from the COMMITTED nodal displacements. Fresh model: disp = 0 ⇒ R = I,
  // u_d,n = 0. After a DB restore / element migration this REBUILDS u_d,n
  // (deliberately not serialized) so the incremental storage coupling
  // differences against the true committed configuration, not zero.
  if (corotActive_) {
    for (int a = 0; a < nNodes_; a++) {
      const Vector &X = theNodes_[a]->getCrds();
      const Vector &d = theNodes_[a]->getDisp();   // COMMITTED disp
      for (int i = 0; i < 3; i++) {
        refCrdsScratch_(a, i) = X(i);
        curCrdsScratch_(a, i) = X(i) + d(i);
      }
    }
    if (theGeom_->update(nNodes_, refCrdsScratch_, curCrdsScratch_) == 0) {
      theGeom_->localizeDisp(uDn_, uDn_);
      uDcur_ = uDn_;
      theGeom_->getCurrentFrame(Rcur_);
      for (int a = 0; a < nNodes_; a++)
        for (int i = 0; i < 3; i++) {
          int row = a * 3 + i;
          for (int bb = 0; bb < nP_; bb++) {
            double s = 0.0;
            for (int j = 0; j < 3; j++)
              s += Rcur_(i, j) * Qblk_[(a * 3 + j) * nP_ + bb];
            QrotBlk_[row * nP_ + bb] = s;
          }
        }
    }
  }

  // fresh geometry ⇒ every derived cache is stale
  stabDirty_ = true;
  k0Valid_ = false;
  k0Dirty_ = true;
  kInitSolidValid_ = false;
  kcSolidValid_ = false;
}

// crack-band characteristic length — the family handshake (LadrunoQuad /
// LadrunoBrick idiom): materials that regularize softening (ASDConcrete3D
// et al.) call ops_TheActiveElement->getCharacteristicLength() on their first
// setTrialStrain; this override feeds them the geometry-true value.
double LadrunoUP::getCharacteristicLength(void)
{
  if (!geomValid_ || charLen_ <= 0.0)
    return this->Element::getCharacteristicLength();
  return charLen_;
}

// ===========================================================================
//  geometry-fixed kernel blocks
// ===========================================================================

// Q (coupling), H (Darcy), S (storage), H̃(unit-α stabilization Laplacian) are
// geometry + scalar-parameter blocks — assembled ONCE here (linear geometry)
// via the shared kernel; never duplicated in the API methods. Under bbar the
// coupling must be energy-consistent with the B-bar solid: mᵀB̄ for column
// (a,i) equals the element-MEAN gradient b̄, so Q is assembled from dNuBar_.
void LadrunoUP::buildStaticBlocks(void)
{
  std::fill(Qblk_.begin(), Qblk_.end(), 0.0);
  std::fill(Hblk_.begin(), Hblk_.end(), 0.0);
  std::fill(Sblk_.begin(), Sblk_.end(), 0.0);
  std::fill(HtUnit_.begin(), HtUnit_.end(), 0.0);

  for (int gp = 0; gp < nGP_; gp++) {
    const double *dNuQ = (formulation_ == 1) ? dNuBar_.data()
                                             : &dNuAll_[gp * nNodes_ * ndm_];
    const double *Np  = &NpAll_[gp * nP_];
    const double *dNp = &dNpAll_[gp * nP_ * ndm_];
    ladruno_up::addQ(dNuQ, Np, nNodes_, nP_, ndm_, biotAlpha_, dv_[gp],
                     Qblk_.data(), nP_);
    ladruno_up::addH(dNp, nP_, ndm_, kbar_, dv_[gp], Hblk_.data(), nP_);
    ladruno_up::addS(Np, nP_, oneOverQbar_, dv_[gp], Sblk_.data(), nP_);
    // H̃ stabilization is an EQUAL-ORDER cure only — never assembled on
    // Taylor–Hood pairs (inf-sup stable without it; WP3.A pin 3). HtUnit_
    // stays zero for TH; stabAlphaValue() also returns 0 (belt-and-braces).
    if (pOrder_ != 1)
      ladruno_up::addHtilde(dNp, nP_, ndm_, 1.0, dv_[gp], HtUnit_.data(), nP_);
  }

  // corot: seed Q_rot = Q (R = I until the first update() refreshes the frame
  // — reference configuration, so queries before any update are consistent).
  if (corotActive_)
    QrotBlk_ = Qblk_;
}

// lazy auto-α (McGann 2012 eq 73 / 2015 eq 74): α = α₀·h²/(K_s + 4G_s/3) with
// h = provider elemSizeH (largest EDGE) and the DRAINED SKELETON moduli probed
// from the material INITIAL tangent — plane strain D₀(0,0)=λ+2G, D₀(0,1)=λ,
// D₀(2,2)=G ⇒ G_s = D₀(2,2), K_s = λ + 2G/3; 3D analogous with D₀(3,3)=G.
// Cached; dirtied by updateMaterialStage / perm setParameter (WP1.A pin) so a
// staged elastic→plastic flip re-probes the stage-appropriate moduli (the SSP
// freeze-at-initial trap deliberately dodged — ADR §3.3).
double LadrunoUP::stabAlphaValue(void)
{
  if (pOrder_ != 0 || stabMode_ == 0)
    return 0.0;                       // TH pairs / -stab off: no H̃ assembled
  if (stabMode_ == 2)
    return stabValue_;                // manual raw α (upstream-style)
  if (!stabDirty_)
    return stabAlpha_;

  const Matrix &D0 = theMaterials_[0]->getInitialTangent();
  double Gs  = (ndm_ == 2) ? D0(2, 2) : D0(3, 3);
  double lam = D0(0, 1);
  double KsDrained = lam + 2.0 * Gs / 3.0;
  double den = KsDrained + 4.0 * Gs / 3.0;
  if (den <= 0.0 || elemH_ <= 0.0) {
    // kernel precondition (LadrunoUPKernel.h): a <=0 denominator would return
    // inf or a NEGATIVE α — anti-stabilization. Guarding is the element's job.
    opserr << "LadrunoUP::stabAlphaValue - element " << this->getTag()
           << ": cannot auto-compute stabilization alpha (K_s + 4G_s/3 = "
           << den << ", h = " << elemH_ << "); stabilization DISABLED for this "
              "element — check the material initial tangent\n";
    stabAlpha_ = 0.0;
  } else {
    stabAlpha_ = ladruno_up::stabAlphaAuto(elemH_, KsDrained, Gs, stabValue_);
  }
  stabDirty_ = false;
  return stabAlpha_;
}

// ===========================================================================
//  state cycle
// ===========================================================================

int LadrunoUP::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0)
    opserr << "LadrunoUP::commitState() - failed in base class\n";
  for (int i = 0; i < nGP_; i++)
    retVal += theMaterials_[i]->commitState();

  // maintain the element's OWN committed solid-K copy for βKc Rayleigh
  // (contract obligation 2 — never the base-class Kc, which snapshots
  // getTangentStiff() and would carry −Q/H into the damping).
  if (betaKc != 0.0 && geomValid_) {
    this->buildSolidK(KcSolid_, false);
    kcSolidValid_ = true;
  }

  // corot (ADR 78 §3.3): advance the committed deformational displacement —
  // the state the incremental storage coupling Qᵀ·Δu_d/Δt differences against.
  if (corotActive_)
    uDn_ = uDcur_;
  return retVal;
}

int LadrunoUP::revertToLastCommit(void)
{
  int retVal = 0;
  for (int i = 0; i < nGP_; i++)
    retVal += theMaterials_[i]->revertToLastCommit();
  return retVal;
}

int LadrunoUP::revertToStart(void)
{
  int retVal = 0;
  for (int i = 0; i < nGP_; i++)
    retVal += theMaterials_[i]->revertToStart();
  // material state (possibly a staged flip) rewinds — re-probe lazily
  stabDirty_ = true;
  k0Dirty_ = true;
  kInitSolidValid_ = false;
  kcSolidValid_ = false;
  if (corotActive_) {                 // rewind the corot committed state too
    uDn_.Zero();
    uDcur_.Zero();
  }
  return retVal;
}

// update(): trial u → strains → setTrialStrain ONLY. p NEVER enters the
// constitutive update (the effective-stress seam, ADR §1.4/§3.2) — it acts
// solely through the Q block in assembly.
//
// -geom corot (ADR 78): the per-evaluation frame refresh lives HERE (seam 0),
// then uSolScratch_ becomes the DEFORMATIONAL displacement u_d = Rᵀx⁰ − X⁰
// (seam 2) so the same reference-B strain path below feeds the material small
// strain in the corotated frame. Q_rot = R̄·Q is refreshed alongside R — the
// one per-evaluation cost the coupling blocks add.
int LadrunoUP::update(void)
{
  if (!geomValid_)
    return -1;

  if (!corotActive_) {
    for (int a = 0; a < nNodes_; a++) {
      const Vector &d = theNodes_[a]->getTrialDisp();
      for (int i = 0; i < ndm_; i++)
        uSolScratch_(a * ndm_ + i) = d(i);
    }
  } else {
    for (int a = 0; a < nNodes_; a++) {
      const Vector &X = theNodes_[a]->getCrds();
      const Vector &d = theNodes_[a]->getTrialDisp();
      for (int i = 0; i < 3; i++) {
        refCrdsScratch_(a, i) = X(i);
        curCrdsScratch_(a, i) = X(i) + d(i);
      }
    }
    if (theGeom_->update(nNodes_, refCrdsScratch_, curCrdsScratch_) != 0) {
      // collapsed/inverted cloud (det(H) <= 0): the seam degraded to a safe
      // R=I state; signal the integrator to cut the step rather than march on.
      return -1;
    }
    theGeom_->localizeDisp(uDcur_, uDcur_);   // aliasing is contract-legal
    uSolScratch_ = uDcur_;                    // strain path below reads uSolScratch_;
                                              // uDcur_ survives (getRayleighDampingForces
                                              // clobbers uSolScratch_ with velocities)
    theGeom_->getCurrentFrame(Rcur_);

    // Q_rot[(a,i), b] = Σ_j R(i,j)·Q[(a,j), b] — rows rotated to the current
    // frame; pressure columns are scalars and do not rotate (ADR 78 §3.2).
    for (int a = 0; a < nNodes_; a++)
      for (int i = 0; i < 3; i++) {
        int row = a * 3 + i;
        for (int bb = 0; bb < nP_; bb++) {
          double s = 0.0;
          for (int j = 0; j < 3; j++)
            s += Rcur_(i, j) * Qblk_[(a * 3 + j) * nP_ + bb];
          QrotBlk_[row * nP_ + bb] = s;
        }
      }
  }

  int ret = 0;
  for (int gp = 0; gp < nGP_; gp++) {
    this->formB(gp, Bscratch_);
    epsScratch_.addMatrixVector(0.0, Bscratch_, uSolScratch_, 1.0);
    ret += theMaterials_[gp]->setTrialStrain(epsScratch_);
  }
  return ret;
}

// ===========================================================================
//  assembly helpers
// ===========================================================================

// solid strain-displacement operator at GP (local solid DOF order a*ndm+i).
// formulation==1: mean-dilatation B-bar — distribute (mean − local) dilatation
// over the normal strains, 2D factor 1/2 (family LadrunoQuad), 3D factor 1/3
// (family LadrunoBrick); shear rows untouched; p-side untouched (ADR §3.4).
void LadrunoUP::formB(int gp, Matrix &B) const
{
  const double *dNu = &dNuAll_[gp * nNodes_ * ndm_];
  const bool bbar = (formulation_ == 1);
  B.Zero();

  if (ndm_ == 2) {
    for (int a = 0; a < nNodes_; a++) {
      double bx = dNu[a * 2 + 0];
      double by = dNu[a * 2 + 1];
      int cx = 2 * a, cy = 2 * a + 1;
      if (bbar) {
        double ex = 0.5 * (dNuBar_[a * 2 + 0] - bx);
        double ey = 0.5 * (dNuBar_[a * 2 + 1] - by);
        B(0, cx) = bx + ex;  B(0, cy) = ey;
        B(1, cx) = ex;       B(1, cy) = by + ey;
      } else {
        B(0, cx) = bx;
        B(1, cy) = by;
      }
      B(2, cx) = by;
      B(2, cy) = bx;
    }
  } else {
    for (int a = 0; a < nNodes_; a++) {
      double bg[3] = {dNu[a * 3 + 0], dNu[a * 3 + 1], dNu[a * 3 + 2]};
      int col[3] = {3 * a, 3 * a + 1, 3 * a + 2};
      if (bbar) {
        double e[3];
        for (int j = 0; j < 3; j++)
          e[j] = (dNuBar_[a * 3 + j] - bg[j]) / 3.0;
        for (int i = 0; i < 3; i++)          // normal rows
          for (int j = 0; j < 3; j++)
            B(i, col[j]) = ((i == j) ? bg[j] : 0.0) + e[j];
      } else {
        B(0, col[0]) = bg[0];
        B(1, col[1]) = bg[1];
        B(2, col[2]) = bg[2];
      }
      // shear rows {γxy, γyz, γzx} — OpenSees 3D Voigt (family LadrunoBrick)
      B(3, col[0]) = bg[1];  B(3, col[1]) = bg[0];
      B(4, col[1]) = bg[2];  B(4, col[2]) = bg[1];
      B(5, col[0]) = bg[2];  B(5, col[2]) = bg[0];
    }
  }
}

void LadrunoUP::buildSolidK(Matrix &Ks, bool initial)
{
  Ks.Zero();
  for (int gp = 0; gp < nGP_; gp++) {
    this->formB(gp, Bscratch_);
    const Matrix &D = initial ? theMaterials_[gp]->getInitialTangent()
                              : theMaterials_[gp]->getTangent();
    Ks.addMatrixTripleProduct(1.0, Bscratch_, D, dv_[gp]);
  }
}

// core-frame TOTAL solid force fu = ∫Bᵀσ′dV − Q·p_trial (ADR 78 §3.2): the one
// vector globalizeForce/Stiff consume under corot. Both legs use the same
// reference-B / core frame, so K/Q frame consistency is structural, and K_geo
// receives the total-stress (effective + pore) prestress — the physically
// correct load stiffness for a saturated medium. Both legs are self-
// equilibrated per component (Σ_a ∇N_a = 0), so the seam's centroid
// back-reaction stays zero. Fills pScratch_ (trial nodal p) as a side effect.
void LadrunoUP::buildCoreForce(Vector &fu)
{
  fu.Zero();
  for (int gp = 0; gp < nGP_; gp++) {
    this->formB(gp, Bscratch_);
    const Vector &sigma = theMaterials_[gp]->getStress();
    fu.addMatrixTransposeVector(1.0, Bscratch_, sigma, dv_[gp]);
  }
  for (int bb = 0; bb < nP_; bb++)
    pScratch_(bb) = theNodes_[carrierNodes_[bb]]->getTrialDisp()(ndm_);
  for (int row = 0; row < nU_; row++) {
    double s = 0.0;
    for (int bb = 0; bb < nP_; bb++)
      s += Qblk_[row * nP_ + bb] * pScratch_(bb);
    fu(row) -= s;
  }
}

// solid mass block: consistent ∫ρ NᵀN dV (ρ = material getRho() = SATURATED
// mixture density, the pinned convention — ADR §3.5), or, under -lumped, the
// HRZ row-sum diagonal (partition of unity ⇒ m_a = ∫ρ N_a dV).
void LadrunoUP::buildSolidMass(Matrix &Ms) const
{
  Ms.Zero();
  for (int gp = 0; gp < nGP_; gp++) {
    double rho = theMaterials_[gp]->getRho();
    if (rho == 0.0)
      continue;
    const double *Nu = &NuAll_[gp * nNodes_];
    if (lumpedMass_) {
      for (int a = 0; a < nNodes_; a++) {
        double m = rho * Nu[a] * dv_[gp];
        for (int i = 0; i < ndm_; i++)
          Ms(a * ndm_ + i, a * ndm_ + i) += m;
      }
    } else {
      for (int a = 0; a < nNodes_; a++)
        for (int c = 0; c < nNodes_; c++) {
          double m = rho * Nu[a] * Nu[c] * dv_[gp];
          for (int i = 0; i < ndm_; i++)
            Ms(a * ndm_ + i, c * ndm_ + i) += m;
        }
    }
  }
}

// SOLID-ONLY Rayleigh (contract obligation 2, ADR §3.2):
//   C_Ray = αM·M_s + βK·K_s(current) + βK0·K_s(initial) + βKc·K_s(committed)
// u-u block only, from the element's OWN copies. Porting reference
// SSPquadUP.cpp:437–449 — which folds betaK0/betaKc onto the CURRENT solid K;
// fixed here with true initial/committed copies.
void LadrunoUP::buildCRaySolid(Matrix &C)
{
  C.Zero();
  if (alphaM != 0.0) {
    this->buildSolidMass(MsolScratch_);
    C.addMatrix(1.0, MsolScratch_, alphaM);
  }
  if (betaK != 0.0) {
    this->buildSolidK(KsolScratch_, false);
    C.addMatrix(1.0, KsolScratch_, betaK);
  }
  if (betaK0 != 0.0) {
    if (!kInitSolidValid_) {
      this->buildSolidK(KinitSolid_, true);
      kInitSolidValid_ = true;
    }
    C.addMatrix(1.0, KinitSolid_, betaK0);
  }
  if (betaKc != 0.0) {
    if (!kcSolidValid_) {
      // first use before any commit: committed == initial trial state
      this->buildSolidK(KcSolid_, false);
      kcSolidValid_ = true;
    }
    C.addMatrix(1.0, KcSolid_, betaKc);
  }
}

void LadrunoUP::gpAccel(int gp, double *acc) const
{
  for (int i = 0; i < ndm_; i++)
    acc[i] = 0.0;
  const double *Nu = &NuAll_[gp * nNodes_];
  for (int a = 0; a < nNodes_; a++) {
    const Vector &an = theNodes_[a]->getTrialAccel();
    for (int i = 0; i < ndm_; i++)
      acc[i] += Nu[a] * an(i);
  }
}

double LadrunoUP::gpPressure(int gp) const
{
  double p = 0.0;
  const double *Np = &NpAll_[gp * nP_];
  for (int bb = 0; bb < nP_; bb++)
    p += Np[bb] * theNodes_[carrierNodes_[bb]]->getTrialDisp()(ndm_);
  return p;
}

// ===========================================================================
//  the six §3.2 contract APIs
// ===========================================================================

// [ K , −Q ; 0 , H ] — the honest-p tangent (Book eq 3.66; unsymmetric once
// paired with getDamp's Qᵀ: the deliberate divergence from upstream's
// symmetric-negative rate-DOF trick).
const Matrix &LadrunoUP::getTangentStiff(void)
{
  K_.Zero();
  if (!geomValid_)
    return K_;

  this->buildSolidK(KsolScratch_, false);

  // corot (ADR 78 §3.1/§3.2): globalize the core solid K through seam 3 with
  // the TOTAL core force (∫Bᵀσ′dV − Q·p) so K_geo carries the total-stress
  // prestress; the u-p column block becomes −R̄Q (dominant term of the exact
  // ∂[globalizeForce]/∂p, the corot v2.0 tangent policy).
  if (corotActive_) {
    this->buildCoreForce(fuScratch_);
    theGeom_->globalizeStiff(KsolScratch_, fuScratch_, KsolScratch_);
  }

  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          K_(uOff_[a] + i, uOff_[c] + j) = KsolScratch_(a * ndm_ + i, c * ndm_ + j);

  const double *Quse = corotActive_ ? QrotBlk_.data() : Qblk_.data();
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++) {
      int row = a * ndm_ + i;
      for (int bb = 0; bb < nP_; bb++)
        K_(uOff_[a] + i, pOff_[carrierNodes_[bb]]) = -Quse[row * nP_ + bb];
    }

  for (int bb = 0; bb < nP_; bb++)
    for (int cc = 0; cc < nP_; cc++)
      K_(pOff_[carrierNodes_[bb]], pOff_[carrierNodes_[cc]]) = Hblk_[bb * nP_ + cc];

  return K_;
}

// [ K₀ , −Q ; 0 , H ] with NONZERO p-rows — a family-idiom copy of upstream's
// solid-only Ki would make initial-stiffness STATICS structurally singular
// (ADR §3.2). Cached; invalidated by updateMaterialStage / perm setParameter.
// Under -geom corot this DELIBERATELY stays the reference-configuration
// operator (R = I semantics — consistent with "initial" and with LadrunoBrick
// sweep #7, which forces cur == ref before globalizing K₀; here we simply
// never let getInitialStiff touch the live frame). ADR 78 §3.6.
const Matrix &LadrunoUP::getInitialStiff(void)
{
  if (!geomValid_) {
    K_.Zero();
    return K_;
  }
  if (k0Valid_ && !k0Dirty_)
    return K0_;

  K0_.Zero();
  this->buildSolidK(KsolScratch_, true);
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          K0_(uOff_[a] + i, uOff_[c] + j) = KsolScratch_(a * ndm_ + i, c * ndm_ + j);

  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++) {
      int row = a * ndm_ + i;
      for (int bb = 0; bb < nP_; bb++)
        K0_(uOff_[a] + i, pOff_[carrierNodes_[bb]]) = -Qblk_[row * nP_ + bb];
    }

  for (int bb = 0; bb < nP_; bb++)
    for (int cc = 0; cc < nP_; cc++)
      K0_(pOff_[carrierNodes_[bb]], pOff_[carrierNodes_[cc]]) = Hblk_[bb * nP_ + cc];

  k0Valid_ = true;
  k0Dirty_ = false;
  return K0_;
}

// [ C_Ray , 0 ; Qᵀ , S + H̃ ] — the rate-multiplying blocks. NEVER the base
// Element::getDamp (it composes βK·getTangentStiff(), which contains −Q/H —
// spurious u-forces from ṗ and spurious ṗ-diffusion damping). H̃ only for
// equal-order pairs with -stab on (McGann; lands on ṗ per our contract).
const Matrix &LadrunoUP::getDamp(void)
{
  damp_.Zero();
  if (!geomValid_)
    return damp_;

  this->buildCRaySolid(CRayScratch_);
  // corot (ADR 78 §3.5): rotate the solid-only Rayleigh block as one —
  // globalizeStiff with a zero core force is the pure rotation R̄CR̄ᵀ (no
  // K_geo; damping is material, not load-stiffness). The αM·M part is exactly
  // invariant (m·I nodal blocks), so this rotates precisely the K-parts.
  if (corotActive_)
    theGeom_->globalizeStiff(CRayScratch_, zeroFu_, CRayScratch_);
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          damp_(uOff_[a] + i, uOff_[c] + j) = CRayScratch_(a * ndm_ + i, c * ndm_ + j);

  // Qᵀ — p-rows, u-cols (the transpose of the tangent's −Q partner: the two
  // Q blocks living in DIFFERENT matrices is what makes the effective
  // transient tangent unsymmetric — the honest-p price, ADR §3.2).
  // corot: (R̄Q)ᵀ — the transpose partner of the tangent's −R̄Q and the
  // dominant Newton sensitivity of the storage coupling. The RESIDUAL does
  // not use this block times velocity: getResistingForceIncInertia swaps it
  // for the incremental Qᵀ·Δu_d/Δt (ADR 78 §3.3 — the chord-defect cure).
  const double *Quse = corotActive_ ? QrotBlk_.data() : Qblk_.data();
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++) {
      int row = a * ndm_ + i;
      for (int bb = 0; bb < nP_; bb++)
        damp_(pOff_[carrierNodes_[bb]], uOff_[a] + i) = Quse[row * nP_ + bb];
    }

  double al = this->stabAlphaValue();
  for (int bb = 0; bb < nP_; bb++)
    for (int cc = 0; cc < nP_; cc++)
      damp_(pOff_[carrierNodes_[bb]], pOff_[carrierNodes_[cc]]) =
          Sblk_[bb * nP_ + cc] + al * HtUnit_[bb * nP_ + cc];

  return damp_;
}

// [ M , 0 ; 0 , 0 ] — solid mixture mass only; p-rows STAY zero in the honest-p
// contract (upstream parks −S here; ours multiplies ṗ and lives in getDamp).
const Matrix &LadrunoUP::getMass(void)
{
  mass_.Zero();
  if (!geomValid_)
    return mass_;

  this->buildSolidMass(MsolScratch_);
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          mass_(uOff_[a] + i, uOff_[c] + j) = MsolScratch_(a * ndm_ + i, c * ndm_ + j);

  return mass_;
}

// [ ∫Bᵀσ′dV − Q·p − body ; H·p − f_seep ] — NO rate terms here: static-path
// correctness (drained solid + steady seepage, the genuinely new capability)
// depends on it. Rate terms enter ONLY through getResistingForceIncInertia.
const Vector &LadrunoUP::getResistingForce(void)
{
  resid_.Zero();
  if (!geomValid_)
    return resid_;

  const double *bAct  = (applyLoad_ == 1) ? appliedB_ : b_;
  const double *fbAct = (applyLoad_ == 1) ? appliedFluidBody_ : fluidBody_;

  // ---- u-rows: internal force + coupling + body force -----------------------
  if (!corotActive_) {
    // linear (P1 path, bit-identical): internal + body accumulated per GP,
    // then −Q·p applied on the assembled u-rows.
    fuScratch_.Zero();
    for (int gp = 0; gp < nGP_; gp++) {
      this->formB(gp, Bscratch_);
      const Vector &sigma = theMaterials_[gp]->getStress();
      fuScratch_.addMatrixTransposeVector(1.0, Bscratch_, sigma, dv_[gp]);

      double rho = theMaterials_[gp]->getRho();   // saturated mixture (pinned)
      if (rho != 0.0) {
        const double *Nu = &NuAll_[gp * nNodes_];
        for (int a = 0; a < nNodes_; a++)
          for (int i = 0; i < ndm_; i++)
            fuScratch_(a * ndm_ + i) -= Nu[a] * rho * bAct[i] * dv_[gp];
      }
    }
    for (int a = 0; a < nNodes_; a++)
      for (int i = 0; i < ndm_; i++)
        resid_(uOff_[a] + i) += fuScratch_(a * ndm_ + i);

    // −Q·p on the u-rows (the pore-pressure force on the skeleton — in the
    // STATIC residual, unlike upstream where it only rides damp·vel)
    for (int bb = 0; bb < nP_; bb++)
      pScratch_(bb) = theNodes_[carrierNodes_[bb]]->getTrialDisp()(ndm_);
    for (int a = 0; a < nNodes_; a++)
      for (int i = 0; i < ndm_; i++) {
        int row = a * ndm_ + i;
        double s = 0.0;
        for (int bb = 0; bb < nP_; bb++)
          s += Qblk_[row * nP_ + bb] * pScratch_(bb);
        resid_(uOff_[a] + i) -= s;
      }
  } else {
    // corot (ADR 78 §3.2): the CORE-frame total force ∫Bᵀσ′dV − Q·p is
    // globalized as ONE vector — the K/Q frame consistency is structural.
    // (buildCoreForce also fills pScratch_, which the p-rows below reuse.)
    this->buildCoreForce(fuScratch_);
    theGeom_->globalizeForce(fuScratch_, fuScratch_);
    for (int a = 0; a < nNodes_; a++)
      for (int i = 0; i < ndm_; i++)
        resid_(uOff_[a] + i) += fuScratch_(a * ndm_ + i);

    // fixed-direction dead body load: a GLOBAL-frame quantity, kept OUT of
    // the core force and added back unrotated AFTER globalize (the
    // LadrunoBrick COROT-1 idiom).
    for (int gp = 0; gp < nGP_; gp++) {
      double rho = theMaterials_[gp]->getRho();
      if (rho == 0.0)
        continue;
      const double *Nu = &NuAll_[gp * nNodes_];
      for (int a = 0; a < nNodes_; a++)
        for (int i = 0; i < ndm_; i++)
          resid_(uOff_[a] + i) -= Nu[a] * rho * bAct[i] * dv_[gp];
    }
  }

  // ---- p-rows: H·p − f_seep --------------------------------------------------
  for (int bb = 0; bb < nP_; bb++) {
    double s = 0.0;
    for (int cc = 0; cc < nP_; cc++)
      s += Hblk_[bb * nP_ + cc] * pScratch_(cc);
    resid_(pOff_[carrierNodes_[bb]]) += s;
  }

  // f_seep: drive = fluidBody − ü_trial when -dynSeepage on (SWANDYNE residual
  // treatment, ADR §3.1 — the deliberate divergence from upstream's dead code),
  // fluidBody only when off (upstream parity; the P1 equivalence gate's leg).
  // ü is interpolated at the GP through Nu; zero in statics (trial accel = 0).
  {
    double fseep[UP_MAX_NODES];
    for (int bb = 0; bb < nP_; bb++)
      fseep[bb] = 0.0;
    double drive[3], acc[3];
    for (int gp = 0; gp < nGP_; gp++) {
      if (dynSeepage_)
        this->gpAccel(gp, acc);
      else
        acc[0] = acc[1] = acc[2] = 0.0;
      for (int i = 0; i < ndm_; i++)
        drive[i] = fbAct[i] - acc[i];
      // corot (ADR 78 §3.4): ∇Np is MATERIAL-frame, the gravity/accel drive is
      // GLOBAL-frame — pull it back: drive_mat = Rᵀ·drive_global.
      if (corotActive_) {
        double dm[3];
        for (int i = 0; i < 3; i++)
          dm[i] = Rcur_(0, i) * drive[0] + Rcur_(1, i) * drive[1]
                + Rcur_(2, i) * drive[2];
        drive[0] = dm[0];  drive[1] = dm[1];  drive[2] = dm[2];
      }
      ladruno_up::addFseep(&dNpAll_[gp * nP_ * ndm_], nP_, ndm_, kbar_,
                           rhoF_, drive, dv_[gp], fseep);
    }
    for (int bb = 0; bb < nP_; bb++)
      resid_(pOff_[carrierNodes_[bb]]) -= fseep[bb];
  }

  // applied nodal loads (inertia unbalance)
  resid_.addVector(1.0, Qload_, -1.0);
  return resid_;
}

// getResistingForce + M·ü + getDamp()·v — the element adds the rate terms
// itself (contract obligation 1): integrators feed getDamp/getMass ONLY into
// the effective tangent and never assemble damp·vel/mass·accel into the
// residual; omitting this converges Newton onto the WRONG equation (the p-row
// would degenerate to steady seepage). The full damp product carries
// [C_Ray·u̇ ; Qᵀ·u̇ + (S+H̃)·ṗ] — the p-slot "velocity" IS ṗ.
const Vector &LadrunoUP::getResistingForceIncInertia(void)
{
  this->getResistingForce();          // fills resid_
  if (!geomValid_)
    return resid_;

  double buf[UP_MAX_DOF];

  for (int a = 0; a < nNodes_; a++) {
    const Vector &an = theNodes_[a]->getTrialAccel();
    for (int i = 0; i < ndm_; i++)
      buf[uOff_[a] + i] = an(i);
    if (pOff_[a] >= 0)
      buf[pOff_[a]] = an(ndm_);       // multiplied by zero p-rows of M anyway
  }
  Vector aVec(buf, nDof_);
  resid_.addMatrixVector(1.0, this->getMass(), aVec, 1.0);

  for (int a = 0; a < nNodes_; a++) {
    const Vector &vn = theNodes_[a]->getTrialVel();
    for (int i = 0; i < ndm_; i++)
      buf[uOff_[a] + i] = vn(i);
    if (pOff_[a] >= 0)
      buf[pOff_[a]] = vn(ndm_);       // = ṗ under the honest-p contract
  }
  Vector vVec(buf, nDof_);
  resid_.addMatrixVector(1.0, this->getDamp(), vVec, 1.0);

  // corot (ADR 78 §3.3): REPLACE the p-row rate terms wholesale. Two measured
  // failure modes force this, and they pin the exact shape of the fix:
  //  (1) THE CHORD DEFECT — (R̄Q)ᵀ·u̇ contracts integrator velocities, and the
  //      chord of a finite rotation increment carries an apparent volumetric
  //      rate 2(1−cosΔθ)/Δt per step: always compressive, Q̄-amplified in
  //      undrained problems (order-1 spurious p at bearing-mechanism
  //      rotations). No velocity-linear operator can remove it (the chord
  //      contraction is representable by no spin) — only the INCREMENTAL form
  //      Qᵀ·(u_d − u_d,committed)/Δt, exactly zero under rigid motion, kills
  //      it.
  //  (2) THE MIXED-SCHEME PUMP — swapping ONLY the coupling to incremental
  //      while leaving S·ṗ on the integrator's Newmark velocity breaks the
  //      skew-symmetry of the discrete coupling pair and PUMPS the structural
  //      ringing (measured: the consolidation column's corot run grows a
  //      ±100·q p-oscillation at Δt = 0.02 where linear decays).
  // The cure is the Book's canonical u-p pairing (GN22 on u, GN11/backward
  // difference on the WHOLE p-row): under corot every p-row rate term goes
  // incremental — Qᵀ·Δu_d/Δt + (S + αH̃)·Δp/Δt — the standard
  // finite-deformation continuity discretization. Committed p comes from the
  // nodes (no extra state). getDamp keeps [(R̄Q)ᵀ, S+H̃] as the Newton
  // sensitivity surrogate (the integrator scales it by c2 = γ/(βΔt) vs the
  // true 1/Δt — a uniform γ/β factor on the p-row rate block, absorbed as an
  // inexact-tangent term under the corot v2.0 policy). Δt ≤ 0 (no transient
  // step context, e.g. post-commit reporting) falls back to the velocity
  // form, where Δu_d = Δp = 0 anyway.
  if (corotActive_) {
    Domain *dom = this->getDomain();
    double dt = (dom != 0) ? dom->getDT() : 0.0;
    if (dt > 0.0) {
      double al = this->stabAlphaValue();
      for (int bb = 0; bb < nP_; bb++) {
        double sOld = 0.0, sNew = 0.0;
        for (int a = 0; a < nNodes_; a++)
          for (int i = 0; i < 3; i++) {
            int row = a * 3 + i;
            sOld += QrotBlk_[row * nP_ + bb] * buf[uOff_[a] + i];
            sNew += Qblk_[row * nP_ + bb] * (uDcur_(row) - uDn_(row));
          }
        for (int cc = 0; cc < nP_; cc++) {
          double spc = Sblk_[bb * nP_ + cc] + al * HtUnit_[bb * nP_ + cc];
          Node *nc = theNodes_[carrierNodes_[cc]];
          sOld += spc * buf[pOff_[carrierNodes_[cc]]];          // Newmark ṗ
          sNew += spc * (nc->getTrialDisp()(ndm_) - nc->getDisp()(ndm_));
        }
        resid_(pOff_[carrierNodes_[bb]]) += sNew / dt - sOld;
      }
    }
  }

  return resid_;
}

// −M·(R·üg) with the p-slots of R·üg zeroed (belt-and-braces: M's p-rows are
// already structurally zero — ADR §3.2 table).
int LadrunoUP::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (!geomValid_)
    return 0;

  double buf[UP_MAX_DOF];
  for (int a = 0; a < nNodes_; a++) {
    const Vector &Raccel = theNodes_[a]->getRV(accel);
    // per-node ndf: carrier (vertex) nodes ndm+1, mid-edge nodes ndm (TH). A
    // uniform ndm+1 check would spuriously reject every mid-edge node.
    int need = (pOff_[a] >= 0) ? (ndm_ + 1) : ndm_;
    if (Raccel.Size() != need) {
      opserr << "LadrunoUP::addInertiaLoadToUnbalance - element " << this->getTag()
             << ": matrix and vector sizes are incompatible\n";
      return -1;
    }
    for (int i = 0; i < ndm_; i++)
      buf[uOff_[a] + i] = Raccel(i);
    if (pOff_[a] >= 0)
      buf[pOff_[a]] = 0.0;            // p-slots zeroed (pin)
  }
  Vector ra(buf, nDof_);
  Qload_.addMatrixVector(1.0, this->getMass(), ra, -1.0);
  return 0;
}

// ===========================================================================
//  solid-only Rayleigh plumbing
// ===========================================================================

// Deliberately does NOT call Element::setRayleighDampingFactors — the base
// helper snapshots Kc = getTangentStiff() (which carries −Q/H under this
// contract) and allocates the shared-static machinery that only serves the
// base getDamp/getResistingForceIncInertia paths we fully override.
int LadrunoUP::setRayleighDampingFactors(double am, double bk, double bk0, double bkc)
{
  alphaM = am;
  betaK  = bk;
  betaK0 = bk0;
  betaKc = bkc;

  if (geomValid_) {
    if (betaKc != 0.0) {
      // snapshot the CURRENT (== last-committed at this call point) solid K
      this->buildSolidK(KcSolid_, false);
      kcSolidValid_ = true;
    }
    if (betaK0 != 0.0 && !kInitSolidValid_) {
      this->buildSolidK(KinitSolid_, true);
      kInitSolidValid_ = true;
    }
  }
  return 0;
}

// shadows the non-virtual base helper (which composes βK·getTangentStiff()):
// SOLID-ONLY C_Ray times the solid trial velocities; p-slots stay zero.
const Vector &LadrunoUP::getRayleighDampingForces(void)
{
  rayForce_.Zero();
  if (!geomValid_)
    return rayForce_;

  this->buildCRaySolid(CRayScratch_);
  if (corotActive_)   // same one-block rotation as getDamp (ADR 78 §3.5)
    theGeom_->globalizeStiff(CRayScratch_, zeroFu_, CRayScratch_);
  for (int a = 0; a < nNodes_; a++) {
    const Vector &vn = theNodes_[a]->getTrialVel();
    for (int i = 0; i < ndm_; i++)
      uSolScratch_(a * ndm_ + i) = vn(i);
  }
  fuScratch_.addMatrixVector(0.0, CRayScratch_, uSolScratch_, 1.0);
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      rayForce_(uOff_[a] + i) = fuScratch_(a * ndm_ + i);
  return rayForce_;
}

// ===========================================================================
//  loads
// ===========================================================================

// zeroLoad restores the ctor -body/-fluidBody values (they are ALWAYS-ON,
// upstream behavior — applyLoad_==0 selects them in getResistingForce).
void LadrunoUP::zeroLoad(void)
{
  Qload_.Zero();
  applyLoad_ = 0;
  for (int i = 0; i < 3; i++) {
    appliedB_[i] = 0.0;
    appliedFluidBody_[i] = 0.0;
  }
}

// LOAD_TAG_SelfWeight REPLACES the active body+fluidBody with loadFactor-scaled
// ctor values — BOTH the solid rows and the seepage source (the staging knob
// the P4 gravity recipe drives, ADR §4.1). Every other elemental load type is
// rejected with a named error (surface tractions come later via standard
// patterns — NOT the SSP side-pressure route, whose index bug ADR §3.5 logs).
int LadrunoUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    if (data.Size() < ndm_) {
      opserr << "LadrunoUP::addLoad - element " << this->getTag()
             << ": SelfWeight data has " << data.Size() << " entries, need "
             << ndm_ << "\n";
      return -1;
    }
    applyLoad_ = 1;
    for (int i = 0; i < ndm_; i++) {
      appliedB_[i]         += loadFactor * data(i) * b_[i];
      appliedFluidBody_[i] += loadFactor * data(i) * fluidBody_[i];
    }
    return 0;
  }

  opserr << "LadrunoUP::addLoad - element " << this->getTag()
         << ": elemental load type " << type
         << " is not supported by LadrunoUP (only SelfWeight; use standard "
            "load patterns for surface loads)\n";
  return -1;
}

// ===========================================================================
//  parameters — the updateMaterialStage transport (contract obligation 3)
// ===========================================================================

int LadrunoUP::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // element permeabilities (live k̄ updates — upstream keep, bug-fixed wiring:
  // brickUP's zPerm→y-id defect is one of the §3.5 fix items)
  if (strcmp(argv[0], "xPerm") == 0)
    return param.addObject(3, this);
  if (strcmp(argv[0], "yPerm") == 0)
    return param.addObject(4, this);
  if (strcmp(argv[0], "zPerm") == 0) {
    if (ndm_ == 3)
      return param.addObject(5, this);
    return -1;
  }

  // per-GP material forwarding (family idiom)
  if ((strstr(argv[0], "material") != 0) && (strcmp(argv[0], "materialState") != 0)) {
    if (argc < 3)
      return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nGP_)
      return theMaterials_[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
    return -1;
  }

  // catch-all → every GP material (the upstream idiom; this is how
  // MaterialStageParameter's {"updateMaterialStage", matTag} reaches the
  // materials). When a material ACCEPTED a stage flip, co-register the element
  // itself so updateParameter(1) can dirty the stabilization-α / K₀ /
  // initial-solid-K caches (ADR §3.2 obligation 3; the SSP freeze-at-initial
  // trap). Co-registering only on acceptance keeps MaterialStageParameter's
  // find-one-element domain scan honest (it stops at the first res != -1).
  int res = -1;
  for (int i = 0; i < nGP_; i++) {
    int matRes = theMaterials_[i]->setParameter(argv, argc, param);
    if (matRes != -1)
      res = matRes;
  }
  if (res != -1 && strcmp(argv[0], "updateMaterialStage") == 0)
    param.addObject(1, this);
  return res;
}

void LadrunoUP::dirtyStageCaches(void)
{
  stabDirty_ = true;
  k0Dirty_ = true;
  kInitSolidValid_ = false;
}

int LadrunoUP::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
    case 1: {
      // updateMaterialStage flipped the materials' constitutive stage: the
      // initial tangent (→ K₀, initial solid K, auto-α skeleton moduli) is
      // stale. Recomputed lazily at next use.
      this->dirtyStageCaches();
      // SIBLING BROADCAST (P4 4.C bug fix): MaterialStageParameter registers
      // only the FIRST accepting element (its domain scan stops there —
      // MaterialStageParameter.cpp:76), but the stage flip lands in the
      // materials' shared static slot, i.e. EVERY LadrunoUP sharing the
      // matTag changed. Measured un-fixed: on a 2-element stack only the
      // scan-found element refreshed auto-α (f2 stayed 1.0000 vs the correct
      // 1.4051 — tests/test_ladruno_up_element_pdmy.py). Dirty every
      // LadrunoUP in the domain: flags only, lazily recomputed, so
      // over-dirtying elements of other matTags costs one cache rebuild.
      Domain *dom = this->getDomain();
      if (dom != 0) {
        ElementIter &els = dom->getElements();
        Element *el;
        while ((el = els()) != 0)
          if (el != this && el->getClassTag() == ELE_TAG_LadrunoUP)
            static_cast<LadrunoUP *>(el)->dirtyStageCaches();
      }
      return 0;
    }
    case 3:
    case 4:
    case 5: {
      kbar_[parameterID - 3] = info.theDouble;
      // H is geometry×k̄: rebuild it (and only it) from the cached gradients
      if (geomValid_) {
        std::fill(Hblk_.begin(), Hblk_.end(), 0.0);
        for (int gp = 0; gp < nGP_; gp++)
          ladruno_up::addH(&dNpAll_[gp * nP_ * ndm_], nP_, ndm_, kbar_,
                           dv_[gp], Hblk_.data(), nP_);
      }
      k0Dirty_ = true;    // K0_ carries H
      stabDirty_ = true;  // pinned: both caches dirty on perm updates
      return 0;
    }
    default:
      return -1;
  }
}

// ===========================================================================
//  responses — pinned IDs 1–4 (WP1.A; anything else → null response)
// ===========================================================================

Response *LadrunoUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  if (argc < 1) return 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoUP");
  output.attr("eleTag", this->getTag());

  if (LadrunoResp::is(argv[0], "stress")) {
    // 1 — EFFECTIVE stress σ′ per GP (the material's vector, nStr comps each)
    theResponse = new ElementResponse(this, 1, Vector(nGP_ * nStr_));
  } else if (LadrunoResp::is(argv[0], "stressTotal")) {
    // 2 — TOTAL stress σ = σ′ − α·m·p per GP (normal components carry the
    //     pore-pressure correction; ADR §3.1 sign convention)
    theResponse = new ElementResponse(this, 2, Vector(nGP_ * nStr_));
  } else if (strcmp(argv[0], "porePressure") == 0) {
    // 3 — per-GP pore pressure Np·p_nodal (upstream ships NOTHING like this —
    //     the §3.5 "zero element-level pore-pressure responses" fix)
    theResponse = new ElementResponse(this, 3, Vector(nGP_));
  } else if (LadrunoResp::is(argv[0], "flux")) {
    // 4 — per-GP Darcy flux −k̄·(∇p − ρ_f·drive), ndm comps each
    theResponse = new ElementResponse(this, 4, Vector(nGP_ * ndm_));
  } else if (LadrunoResp::is(argv[0], "material")) {
    // family idiom: material $gpNum <...> forwarding
    if (argc > 1) {
      int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= nGP_)
        theResponse = theMaterials_[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
    }
  // Ladruno — the rest of the family contract (force / strain / stiff /
  // stiffInitial / charLength), which this element was missing entirely.
  } else if (LadrunoResp::is(argv[0], "strain")) {
    // 5 — EFFECTIVE strain ε per GP, same layout as response 1
    theResponse = new ElementResponse(this, 5, Vector(nGP_ * nStr_));
  } else if (LadrunoResp::is(argv[0], "force")) {
    theResponse = new ElementResponse(this, 6, Vector(this->getNumDOF()));
  } else if (LadrunoResp::is(argv[0], "charLength")) {
    output.tag("ResponseType", "lch");
    theResponse = new ElementResponse(this, 7, Vector(1));
  } else if (LadrunoResp::is(argv[0], "stiff")) {
    theResponse = new ElementResponse(this, 8,
                                      Matrix(this->getNumDOF(), this->getNumDOF()));
  } else if (LadrunoResp::is(argv[0], "stiffInitial")) {
    theResponse = new ElementResponse(this, 9,
                                      Matrix(this->getNumDOF(), this->getNumDOF()));
  // Ladruno — dampingForce/dynamicForce are answered HERE, not by the base
  // chain below: Element::getRayleighDampingForces would compose
  // betaK*getTangentStiff() and inject the -Q / H coupling blocks into what it
  // calls "Rayleigh damping". Rayleigh is SOLID-ONLY on this element (class
  // header) and this element's own shadowing getRayleighDampingForces is the
  // one that honours that. inertialForce/globalForce ARE left to the base:
  // getResistingForceIncInertia is virtual, so the base already gets ours.
  } else if (LadrunoResp::is(argv[0], "dampingForce")) {
    theResponse = new ElementResponse(this, 10, Vector(this->getNumDOF()));
  } else if (LadrunoResp::is(argv[0], "dynamicForce")) {
    theResponse = new ElementResponse(this, 11, Vector(this->getNumDOF()));
  }
  // IDs 1–11 (pin: chosen clear of the plane-family's 21/stressZZ space).

  output.endTag();

  // Ladruno — base vocabulary (globalForce, dampingForce, dynamicForce,
  // inertialForce); Element::setResponse opens its own ElementOutput tag, so
  // this MUST come after endTag().
  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, output);
  return theResponse;
}

int LadrunoUP::getResponse(int responseID, Information &eleInfo)
{
  if (!geomValid_)
    return -1;

  if (responseID == 1 || responseID == 2) {
    Vector v(nGP_ * nStr_);
    int nNormal = (ndm_ == 2) ? 2 : 3;   // in-vector normal comps (plane-strain
                                         // σzz is not in the 3-vector)
    for (int gp = 0; gp < nGP_; gp++) {
      const Vector &s = theMaterials_[gp]->getStress();
      for (int k = 0; k < nStr_; k++)
        v(gp * nStr_ + k) = s(k);
      if (responseID == 2) {
        // total = effective − α·p on the normal components (σ = σ′ − α·m·p,
        // tension-positive σ / compression-positive p — ADR §3.1)
        double p = this->gpPressure(gp);
        for (int k = 0; k < nNormal; k++)
          v(gp * nStr_ + k) -= biotAlpha_ * p;
      }
    }
    return eleInfo.setVector(v);
  }

  if (responseID == 3) {
    Vector v(nGP_);
    for (int gp = 0; gp < nGP_; gp++)
      v(gp) = this->gpPressure(gp);
    return eleInfo.setVector(v);
  }

  if (responseID == 4) {
    // Darcy flux q = −k̄·(∇p − ρ_f·drive), drive = (b_f − ü) with -dynSeepage on
    // (b_f only when off) — consistent with the f_seep residual treatment.
    // corot (ADR 78 §3.4): computed in the MATERIAL frame (∇p and k̄ live
    // there; the global drive is pulled back by Rᵀ), then pushed forward by R
    // so the recorder stays global-frame. Linear: R = I, bit-identical.
    Vector v(nGP_ * ndm_);
    const double *fbAct = (applyLoad_ == 1) ? appliedFluidBody_ : fluidBody_;
    double acc[3];
    for (int gp = 0; gp < nGP_; gp++) {
      if (dynSeepage_)
        this->gpAccel(gp, acc);
      else
        acc[0] = acc[1] = acc[2] = 0.0;
      double drive[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < ndm_; i++)
        drive[i] = fbAct[i] - acc[i];
      if (corotActive_) {
        double dm[3];
        for (int i = 0; i < 3; i++)
          dm[i] = Rcur_(0, i) * drive[0] + Rcur_(1, i) * drive[1]
                + Rcur_(2, i) * drive[2];
        drive[0] = dm[0];  drive[1] = dm[1];  drive[2] = dm[2];
      }
      const double *dNp = &dNpAll_[gp * nP_ * ndm_];
      double qmat[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < ndm_; i++) {
        double gradP = 0.0;
        for (int bb = 0; bb < nP_; bb++)
          gradP += dNp[bb * ndm_ + i] *
                   theNodes_[carrierNodes_[bb]]->getTrialDisp()(ndm_);
        qmat[i] = -kbar_[i] * (gradP - rhoF_ * drive[i]);
      }
      if (corotActive_) {
        for (int i = 0; i < 3; i++)
          v(gp * ndm_ + i) = Rcur_(i, 0) * qmat[0] + Rcur_(i, 1) * qmat[1]
                           + Rcur_(i, 2) * qmat[2];
      } else {
        for (int i = 0; i < ndm_; i++)
          v(gp * ndm_ + i) = qmat[i];
      }
    }
    return eleInfo.setVector(v);
  }

  if (responseID == 5) {
    Vector v(nGP_ * nStr_);
    for (int gp = 0; gp < nGP_; gp++) {
      const Vector &e = theMaterials_[gp]->getStrain();
      for (int k = 0; k < nStr_; k++)
        v(gp * nStr_ + k) = e(k);
    }
    return eleInfo.setVector(v);
  }

  if (responseID == 6)
    return eleInfo.setVector(this->getResistingForce());

  if (responseID == 7) {
    static Vector lch(1);
    lch(0) = this->getCharacteristicLength();
    return eleInfo.setVector(lch);
  }

  if (responseID == 8)
    return eleInfo.setMatrix(this->getTangentStiff());

  if (responseID == 9)
    return eleInfo.setMatrix(this->getInitialStiff());

  if (responseID == 10)
    return eleInfo.setVector(this->getRayleighDampingForces());

  if (responseID == 11) {
    // dynamic = inertia-inclusive resisting − (solid-only) damping − static
    Vector v(this->getResistingForceIncInertia());
    v -= this->getRayleighDampingForces();
    v -= this->getResistingForce();
    return eleInfo.setVector(v);
  }

  return this->Element::getResponse(responseID, eleInfo);
}

// ===========================================================================
//  print / IO
// ===========================================================================

void LadrunoUP::Print(OPS_Stream &s, int flag)
{
  static const char *shapeNames[] = {"none", "T3", "Q4", "H8", "BT6", "BTET10"};

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoUP, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes_;
    s << "\tshape: " << shapeNames[shapeKind_]
      << "  formulation: " << (formulation_ == 1 ? "bbar" : "std")
      << "  geom: " << (theGeom_ ? theGeom_->getName() : "linear")
      << "  pOrder: " << (pOrder_ == 1 ? "linear(TH)" : "equal")
      << "  nP(carriers): " << nP_ << " / " << nNodes_ << " nodes"
      << "  nDof: " << nDof_ << "\n";
    if (ndm_ == 2)
      s << "\tthickness: " << thickness_ << "\n";
    s << "\tKf: " << Kf_ << "  poro: " << poro_ << "  rhoF: " << rhoF_
      << "  biotAlpha: " << biotAlpha_ << "  Ks(grain): " << Ks_ << "\n";
    s << "\tperm (k_hydraulic/gammaW):";
    for (int i = 0; i < ndm_; i++)
      s << " " << kbar_[i];
    s << "\n\tbody:";
    for (int i = 0; i < ndm_; i++)
      s << " " << b_[i];
    s << "  fluidBody:";
    for (int i = 0; i < ndm_; i++)
      s << " " << fluidBody_[i];
    s << "\n\tstab: " << (stabMode_ == 0 ? "off" : (stabMode_ == 1 ? "auto" : "manual"))
      << " (value " << stabValue_ << ")  lumped: " << (lumpedMass_ ? 1 : 0)
      << "  dynSeepage: " << (dynSeepage_ ? "on" : "off") << "\n";
    if (theMaterials_ != 0 && theMaterials_[0] != 0)
      theMaterials_[0]->Print(s, flag);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoUP\", ";
    s << "\"nodes\": [";
    for (int n = 0; n < nNodes_; n++) {
      s << connectedExternalNodes_(n);
      if (n < nNodes_ - 1)
        s << ", ";
    }
    s << "], ";
    if (ndm_ == 2)
      s << "\"thickness\": " << thickness_ << ", ";
    s << "\"material\": \""
      << ((theMaterials_ != 0 && theMaterials_[0] != 0) ? theMaterials_[0]->getTag() : 0)
      << "\"}";
  }
}

// ---------------------------------------------------------------------------
// Serialization (WP1.D). Ships EVERY ctor argument + the analysis-relevant
// scalar state, the topology, and the per-GP materials (class/db tag + their
// own sendSelf → committed constitutive state). Exercised by `database File`
// save/restore (Domain::sendSelf/recvSelf drive Element::send/recvSelf through
// FileDatastore-as-Channel) and by any element-migrating parallel channel.
//
// What is NOT serialized, and why (each a deliberate, documented choice — the
// opposite of the §3.5 upstream holes where massType/stab-α just went missing):
//   * geometry + kernel blocks (xy_, Nu/dNu caches, Q/H/S/H̃, dv_, dNuBar_,
//     charLen_, elemH_) — pure functions of the node coordinates; setDomain
//     rebuilds them on the receiver (recvSelf leaves geomValid_ = false).
//   * solid-K Rayleigh caches (K0_, KinitSolid_, KcSolid_) — pure functions of
//     geometry + material tangent, both restored; rebuilt LAZILY after
//     setDomain (buildCRaySolid / commitState). The Rayleigh FACTORS
//     (alphaM/betaK/betaK0/betaKc) ARE shipped so the rebuild is triggered
//     correctly; between restore and the first commit, betaKc uses the
//     "committed == current trial" fallback already coded in buildCRaySolid.
//   * transient load state (applyLoad_, appliedB_/appliedFluidBody_, Qload_) —
//     regenerated every step by zeroLoad + addLoad + addInertiaLoadToUnbalance
//     (family precedent: LadrunoQuad ships none of it). Reset to the null-ctor
//     baseline on the receiver.
// oneOverQbar_ is shipped (not recomputed) for exact identity across the wire.
// The geometry-method id (ADR 78) is shipped as data(29); the corot wrapper is
// STATELESS (R is a pure function of current nodal positions), so the id alone
// reconstructs it — same contract as LadrunoBrick.
static const int UP_NDATA = 30;   // scalar payload width (see layout below)

int LadrunoUP::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  // ---- scalar/config payload: every ctor arg + Rayleigh + storage coeff -----
  Vector data(UP_NDATA);
  data(0)  = this->getTag();
  data(1)  = ndm_;
  data(2)  = nNodes_;
  data(3)  = thickness_;
  data(4)  = Kf_;
  data(5)  = poro_;
  data(6)  = rhoF_;
  data(7)  = biotAlpha_;
  data(8)  = Ks_;
  data(9)  = kbar_[0];
  data(10) = kbar_[1];
  data(11) = kbar_[2];
  data(12) = b_[0];
  data(13) = b_[1];
  data(14) = b_[2];
  data(15) = fluidBody_[0];
  data(16) = fluidBody_[1];
  data(17) = fluidBody_[2];
  data(18) = formulation_;
  data(19) = pOrder_;
  data(20) = lumpedMass_ ? 1.0 : 0.0;
  data(21) = stabMode_;
  data(22) = stabValue_;
  data(23) = dynSeepage_ ? 1.0 : 0.0;
  data(24) = oneOverQbar_;
  data(25) = alphaM;
  data(26) = betaK;
  data(27) = betaK0;
  data(28) = betaKc;
  data(29) = theGeom_ ? theGeom_->getMethodID()
                      : SolidTransformation::METHOD_LINEAR;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoUP::sendSelf - element " << this->getTag()
           << ": failed to send data Vector\n";
    return res;
  }

  // ---- topology + per-GP material class/db tags in one ID -------------------
  // layout: [0..nNodes_-1]              node tags
  //         [nNodes_ .. +nGP_-1]        material class tags
  //         [nNodes_+nGP_ .. +nGP_-1]   material db tags
  ID idData(nNodes_ + 2 * nGP_);
  for (int n = 0; n < nNodes_; n++)
    idData(n) = connectedExternalNodes_(n);
  for (int g = 0; g < nGP_; g++) {
    idData(nNodes_ + g) = theMaterials_[g]->getClassTag();
    int matDbTag = theMaterials_[g]->getDbTag();
    if (matDbTag == 0) {                       // family idiom: assign on demand
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theMaterials_[g]->setDbTag(matDbTag);
    }
    idData(nNodes_ + nGP_ + g) = matDbTag;
  }
  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoUP::sendSelf - element " << this->getTag()
           << ": failed to send ID\n";
    return res;
  }

  // ---- per-GP materials (their own committed constitutive state) ------------
  for (int g = 0; g < nGP_; g++) {
    res = theMaterials_[g]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LadrunoUP::sendSelf - element " << this->getTag()
             << ": failed to send material at GP " << g << "\n";
      return res;
    }
  }
  return 0;
}

int LadrunoUP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  Vector data(UP_NDATA);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoUP::recvSelf - failed to receive data Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  ndm_          = (int)data(1);
  nNodes_       = (int)data(2);
  thickness_    = data(3);
  Kf_           = data(4);
  poro_         = data(5);
  rhoF_         = data(6);
  biotAlpha_    = data(7);
  Ks_           = data(8);
  kbar_[0]      = data(9);
  kbar_[1]      = data(10);
  kbar_[2]      = data(11);
  b_[0]         = data(12);
  b_[1]         = data(13);
  b_[2]         = data(14);
  fluidBody_[0] = data(15);
  fluidBody_[1] = data(16);
  fluidBody_[2] = data(17);
  formulation_  = (int)data(18);
  pOrder_       = (int)data(19);
  lumpedMass_   = (data(20) != 0.0);
  stabMode_     = (int)data(21);
  stabValue_    = data(22);
  dynSeepage_   = (data(23) != 0.0);
  oneOverQbar_  = data(24);
  alphaM        = data(25);   // solid-only Rayleigh factors (base-class members)
  betaK         = data(26);
  betaK0        = data(27);
  betaKc        = data(28);

  // geometry method (ADR 78): rebuild the stateless wrapper from its id
  if (theGeom_) {
    delete theGeom_;
    theGeom_ = 0;
  }
  theGeom_ = SolidTransformation::create((int)data(29));
  if (theGeom_ == 0) {
    opserr << "WARNING LadrunoUP::recvSelf - element " << this->getTag()
           << ": unknown geometry-method id " << (int)data(29)
           << " in the incoming stream; using linear\n";
    theGeom_ = new SolidTransformationLinear();
  }
  corotActive_ = (theGeom_->getMethodID() == SolidTransformation::METHOD_COROT);

  // rebuild every size-derived member + per-instance buffer (the SAME helper
  // the ctor uses — no drift). This sizes uOff_/pOff_/carrierNodes_/nGP_ so the
  // ID (nNodes_ + 2*nGP_) can be sized below.
  this->configureSizing();
  if (shapeKind_ == SHAPE_NONE) {
    opserr << "WARNING LadrunoUP::recvSelf - element " << this->getTag()
           << ": illegal (ndm=" << ndm_ << ", nNodes=" << nNodes_
           << ") in the incoming stream\n";
    return -1;
  }

  // transient / cache state → null-ctor baseline (setDomain + first commit
  // rebuild geometry, kernel blocks and the solid-K Rayleigh caches).
  applyLoad_ = 0;
  for (int i = 0; i < 3; i++) {
    appliedB_[i] = 0.0;
    appliedFluidBody_[i] = 0.0;
  }
  Qload_.Zero();
  geomValid_ = false;
  k0Valid_ = false;  k0Dirty_ = true;
  stabDirty_ = true;
  kInitSolidValid_ = false;  kcSolidValid_ = false;
  for (int n = 0; n < UP_MAX_NODES; n++)
    theNodes_[n] = 0;

  // ---- topology + material class/db tags ------------------------------------
  ID idData(nNodes_ + 2 * nGP_);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoUP::recvSelf - element " << this->getTag()
           << ": failed to receive ID\n";
    return res;
  }
  connectedExternalNodes_.resize(nNodes_);   // null ctor sized it to 0
  for (int n = 0; n < nNodes_; n++)
    connectedExternalNodes_(n) = idData(n);

  // ---- per-GP materials via the broker --------------------------------------
  if (theMaterials_ != 0) {                   // re-recv into a live element
    for (int g = 0; g < nGP_; g++)
      if (theMaterials_[g])
        delete theMaterials_[g];
    delete[] theMaterials_;
    theMaterials_ = 0;
  }
  theMaterials_ = new NDMaterial *[nGP_];
  for (int g = 0; g < nGP_; g++)
    theMaterials_[g] = 0;
  for (int g = 0; g < nGP_; g++) {
    int matClassTag = idData(nNodes_ + g);
    int matDbTag    = idData(nNodes_ + nGP_ + g);
    theMaterials_[g] = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterials_[g] == 0) {
      opserr << "LadrunoUP::recvSelf - element " << this->getTag()
             << ": broker could not create NDMaterial classTag " << matClassTag
             << " at GP " << g << "\n";
      return -1;
    }
    theMaterials_[g]->setDbTag(matDbTag);
    res = theMaterials_[g]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING LadrunoUP::recvSelf - element " << this->getTag()
             << ": failed to receive material at GP " << g << "\n";
      return res;
    }
  }
  return 0;
}
