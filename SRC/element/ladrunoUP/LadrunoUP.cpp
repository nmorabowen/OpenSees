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
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
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
                     bool dynSeepage)
 : Element(tag, ELE_TAG_LadrunoUP),
   ndm_(ndm), nNodes_(nodeTags.Size()), shapeKind_(SHAPE_NONE),
   nGP_(0), nU_(0), nP_(0), nStr_(0), nDof_(0),
   connectedExternalNodes_(nodeTags.Size()), theMaterials_(0),
   thickness_(thickness), Kf_(Kf), poro_(poro), rhoF_(rhoF),
   biotAlpha_(biotAlpha), Ks_(Ks),
   formulation_(formulation), pOrder_(pOrder), lumpedMass_(lumped),
   stabMode_(stabMode), stabValue_(stabValue), dynSeepage_(dynSeepage),
   oneOverQbar_(0.0),
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

  // ---- provider selection by (ndm, nNodes) — WP1.A pin -----------------------
  if      (ndm_ == 2 && nNodes_ == 3)  shapeKind_ = SHAPE_T3;
  else if (ndm_ == 2 && nNodes_ == 4)  shapeKind_ = SHAPE_Q4;
  else if (ndm_ == 2 && nNodes_ == 6)  shapeKind_ = SHAPE_BT6;
  else if (ndm_ == 3 && nNodes_ == 8)  shapeKind_ = SHAPE_H8;
  else if (ndm_ == 3 && nNodes_ == 10) shapeKind_ = SHAPE_BTET10;
  else {
    opserr << "LadrunoUP::LadrunoUP - element " << tag
           << ": no shape provider for (ndm=" << ndm_ << ", nNodes=" << nNodes_
           << "); legal pairs are (2,3) (2,4) (2,6) (3,8) (3,10)\n";
    return;   // shapeKind_ stays SHAPE_NONE; setDomain refuses loudly
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

  for (int n = 0; n < nNodes_; n++)
    connectedExternalNodes_(n) = nodeTags(n);

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

  // ---- P1 scope guard (WP1.A pin): equal-order lanes only -------------------
  // The parser already rejects 6/10-node input; guard anyway (belt-and-braces).
  if (shapeKind_ == SHAPE_BT6 || shapeKind_ == SHAPE_BTET10) {
    opserr << "LadrunoUP::setDomain - element " << this->getTag()
           << ": Bézier Taylor–Hood lanes land at P3 — heterogeneous-ndf "
              "plumbing not yet built; element not activated\n";
    return;
  }
  if (pOrder_ != 0) {
    opserr << "LadrunoUP::setDomain - element " << this->getTag()
           << ": -pOrder linear (Taylor–Hood) is reserved until P3; "
              "P1 supports equal-order only; element not activated\n";
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

  // P1: every node carries [u..., p] — ndf = ndm+1 (equal-order pin)
  for (int n = 0; n < nNodes_; n++) {
    if (theNodes_[n]->getNumberDOF() != ndm_ + 1) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": node " << connectedExternalNodes_(n) << " has ndf = "
             << theNodes_[n]->getNumberDOF() << " but the u-p element needs ndf = "
             << ndm_ + 1 << " (u" << ndm_ << " + p); element not activated\n";
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

  double measure = 0.0;   // area (2D) | volume (3D), thickness NOT folded
  for (int gp = 0; gp < nGP_; gp++) {
    double detJ = 0.0;
    this->shapeEval(gp, &NuAll_[gp * nNodes_], &dNuAll_[gp * nNodes_ * ndm_], &detJ);
    // the providers are silent on singular Jacobians (kernel idiom) — the
    // element is the mandated guard: reject detJ <= 0 BEFORE assembling.
    if (detJ <= 0.0) {
      opserr << "LadrunoUP::setDomain - element " << this->getTag()
             << ": non-positive Jacobian detJ = " << detJ << " at GP " << gp
             << " (bad node ordering or degenerate geometry); element not activated\n";
      geomValid_ = false;
      return;
    }
    double w = this->shapeGPWeight(gp);
    measure += w * detJ;
    dv_[gp] = w * detJ * ((ndm_ == 2) ? thickness_ : 1.0);
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
    ladruno_up::addHtilde(dNp, nP_, ndm_, 1.0, dv_[gp], HtUnit_.data(), nP_);
  }
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
  return retVal;
}

// update(): trial u → strains → setTrialStrain ONLY. p NEVER enters the
// constitutive update (the effective-stress seam, ADR §1.4/§3.2) — it acts
// solely through the Q block in assembly.
int LadrunoUP::update(void)
{
  if (!geomValid_)
    return -1;

  for (int a = 0; a < nNodes_; a++) {
    const Vector &d = theNodes_[a]->getTrialDisp();
    for (int i = 0; i < ndm_; i++)
      uSolScratch_(a * ndm_ + i) = d(i);
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
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          K_(uOff_[a] + i, uOff_[c] + j) = KsolScratch_(a * ndm_ + i, c * ndm_ + j);

  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++) {
      int row = a * ndm_ + i;
      for (int bb = 0; bb < nP_; bb++)
        K_(uOff_[a] + i, pOff_[carrierNodes_[bb]]) = -Qblk_[row * nP_ + bb];
    }

  for (int bb = 0; bb < nP_; bb++)
    for (int cc = 0; cc < nP_; cc++)
      K_(pOff_[carrierNodes_[bb]], pOff_[carrierNodes_[cc]]) = Hblk_[bb * nP_ + cc];

  return K_;
}

// [ K₀ , −Q ; 0 , H ] with NONZERO p-rows — a family-idiom copy of upstream's
// solid-only Ki would make initial-stiffness STATICS structurally singular
// (ADR §3.2). Cached; invalidated by updateMaterialStage / perm setParameter.
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
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++)
      for (int c = 0; c < nNodes_; c++)
        for (int j = 0; j < ndm_; j++)
          damp_(uOff_[a] + i, uOff_[c] + j) = CRayScratch_(a * ndm_ + i, c * ndm_ + j);

  // Qᵀ — p-rows, u-cols (the transpose of the tangent's −Q partner: the two
  // Q blocks living in DIFFERENT matrices is what makes the effective
  // transient tangent unsymmetric — the honest-p price, ADR §3.2)
  for (int a = 0; a < nNodes_; a++)
    for (int i = 0; i < ndm_; i++) {
      int row = a * ndm_ + i;
      for (int bb = 0; bb < nP_; bb++)
        damp_(pOff_[carrierNodes_[bb]], uOff_[a] + i) = Qblk_[row * nP_ + bb];
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

  // ---- u-rows: internal force + body force ----------------------------------
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

  // ---- −Q·p on the u-rows (the pore-pressure force on the skeleton — in the
  //      STATIC residual, unlike upstream where it only rides damp·vel) -------
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
    if (Raccel.Size() != ndm_ + 1) {
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

int LadrunoUP::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
    case 1:
      // updateMaterialStage flipped the materials' constitutive stage: the
      // initial tangent (→ K₀, initial solid K, auto-α skeleton moduli) is
      // stale. Recomputed lazily at next use.
      stabDirty_ = true;
      k0Dirty_ = true;
      kInitSolidValid_ = false;
      return 0;
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
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoUP");
  output.attr("eleTag", this->getTag());

  if (strcmp(argv[0], "stresses") == 0 || strcmp(argv[0], "stress") == 0) {
    // 1 — EFFECTIVE stress σ′ per GP (the material's vector, nStr comps each)
    theResponse = new ElementResponse(this, 1, Vector(nGP_ * nStr_));
  } else if (strcmp(argv[0], "stressesTotal") == 0 ||
             strcmp(argv[0], "stressTotal") == 0) {
    // 2 — TOTAL stress σ = σ′ − α·m·p per GP (normal components carry the
    //     pore-pressure correction; ADR §3.1 sign convention)
    theResponse = new ElementResponse(this, 2, Vector(nGP_ * nStr_));
  } else if (strcmp(argv[0], "porePressure") == 0) {
    // 3 — per-GP pore pressure Np·p_nodal (upstream ships NOTHING like this —
    //     the §3.5 "zero element-level pore-pressure responses" fix)
    theResponse = new ElementResponse(this, 3, Vector(nGP_));
  } else if (strcmp(argv[0], "flux") == 0 || strcmp(argv[0], "darcyFlux") == 0) {
    // 4 — per-GP Darcy flux −k̄·(∇p − ρ_f·drive), ndm comps each
    theResponse = new ElementResponse(this, 4, Vector(nGP_ * ndm_));
  } else if (strcmp(argv[0], "material") == 0 ||
             strcmp(argv[0], "integrPoint") == 0) {
    // family idiom: material $gpNum <...> forwarding
    if (argc > 1) {
      int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= nGP_)
        theResponse = theMaterials_[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
    }
  }
  // IDs 1–4 only (pin: chosen clear of the plane-family's 21/stressZZ space);
  // anything else falls through to the null response.

  output.endTag();
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
    Vector v(nGP_ * ndm_);
    const double *fbAct = (applyLoad_ == 1) ? appliedFluidBody_ : fluidBody_;
    double acc[3];
    for (int gp = 0; gp < nGP_; gp++) {
      if (dynSeepage_)
        this->gpAccel(gp, acc);
      else
        acc[0] = acc[1] = acc[2] = 0.0;
      const double *dNp = &dNpAll_[gp * nP_ * ndm_];
      for (int i = 0; i < ndm_; i++) {
        double gradP = 0.0;
        for (int bb = 0; bb < nP_; bb++)
          gradP += dNp[bb * ndm_ + i] *
                   theNodes_[carrierNodes_[bb]]->getTrialDisp()(ndm_);
        v(gp * ndm_ + i) = -kbar_[i] * (gradP - rhoF_ * (fbAct[i] - acc[i]));
      }
    }
    return eleInfo.setVector(v);
  }

  return -1;
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
      << "  pOrder: " << (pOrder_ == 1 ? "linear(TH)" : "equal") << "\n";
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

// P1 wave-1 STUBS — the full serialization (EVERY ctor arg incl. -lumped /
// alpha0 / dynSeepage) + MP smoke test is WP1.D (wave 2). Returning -1 keeps
// any premature MP/database use loudly broken instead of silently wrong.
int LadrunoUP::sendSelf(int commitTag, Channel &theChannel)
{
  (void)commitTag;
  (void)theChannel;
  opserr << "LadrunoUP::sendSelf - element " << this->getTag()
         << ": not implemented yet — P1 wave-2 (WP1.D)\n";
  return -1;
}

int LadrunoUP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  (void)commitTag;
  (void)theChannel;
  (void)theBroker;
  opserr << "LadrunoUP::recvSelf - element " << this->getTag()
         << ": not implemented yet — P1 wave-2 (WP1.D)\n";
  return -1;
}
