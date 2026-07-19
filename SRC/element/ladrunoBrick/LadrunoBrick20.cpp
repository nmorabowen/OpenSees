/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// LadrunoBrick20 — 20-node serendipity quadratic hexahedron (Ladruno fork).
//
// Small-strain, geometrically linear; two formulations under one class tag:
//   std (P1) — full 27-pt Gauss. Correctness anchor: reduce-to upstream
//              Twenty_Node_Brick to ~1e-12 (K, resisting force, consistent
//              mass, per-GP stress) on a distorted hex under mixed loads.
//   uri (P2) — uniform 2x2x2 reduced integration (the C3D20R analog), 8
//              material points at the Barlow points. Anchor: the P0 sympy
//              oracle (tests/hex20_reference.py) — upstream has no reduced
//              H20 to reduce to. MASS / body force / volume stay 27-pt
//              (formulation-independent BY DESIGN, ADR 72 §2.3/§3.5).
// All shape/GP/B/mass/volume math is consumed from the pure P0 kernel
// Ladruno::hex20::* (LadrunoHex20Shape.h). See
// Ladruno_implementation/72_ladruno_second_order_brick_adr.md.  // Ladruno

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <LadrunoBrick20.h>
#include <LadrunoHex20Shape.h>            // Ladruno — P0 pure H20 kernel (shape/GP/B/mass)
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Response.h>        // Ladruno — material "damage" probe (U1 advisory)
#include <Information.h>
#include <DummyStream.h>

// static data (sized 60)
Matrix  LadrunoBrick20::stiff(60, 60);
Vector  LadrunoBrick20::resid(60);
Matrix  LadrunoBrick20::mass(60, 60);

// process-once advisory / notice flags (NOT serialized; see header)  // Ladruno
bool LadrunoBrick20::advisedDamage   = false;
bool LadrunoBrick20::warnedLumped    = false;
bool LadrunoBrick20::warnedUriCoerce = false;

const char *
LadrunoBrick20::formulationName(Formulation f)
{
  switch (f) {
  case Formulation::STD:  return "std";
  case Formulation::URI:  return "uri";
  }
  return "std";
}

//null constructor (broker)
LadrunoBrick20::LadrunoBrick20()
  :Element(0, ELE_TAG_LadrunoBrick20),
   connectedExternalNodes(NEN),
   formulation(Formulation::STD),
   applyLoad(0), load(0), Ki(0), M0(0), massType(0),
   geomCached(false), badGeom(false), warnedBadUse(false), hasMass(false)
{
  for (int i = 0; i < NEN; i++) nodePointers[i] = 0;
  for (int i = 0; i < NGP; i++) {
    materialPointers[i] = 0;
    theDamping[i] = 0;
    cachedRho[i] = 0.0;
  }
  b[0] = 0.0; b[1] = 0.0; b[2] = 0.0;
  appliedB[0] = 0.0; appliedB[1] = 0.0; appliedB[2] = 0.0;
}

//full constructor
LadrunoBrick20::LadrunoBrick20(int tag,
                               int node1,  int node2,  int node3,  int node4,
                               int node5,  int node6,  int node7,  int node8,
                               int node9,  int node10, int node11, int node12,
                               int node13, int node14, int node15, int node16,
                               int node17, int node18, int node19, int node20,
                               NDMaterial &theMaterial,
                               Formulation form,
                               double b1, double b2, double b3,
                               int matype,
                               Damping *damping)
  :Element(tag, ELE_TAG_LadrunoBrick20),
   connectedExternalNodes(NEN),
   formulation(form),
   applyLoad(0), load(0), Ki(0), M0(0), massType(matype),
   geomCached(false), badGeom(false), warnedBadUse(false), hasMass(false)
{
  connectedExternalNodes(0)  = node1;   connectedExternalNodes(1)  = node2;
  connectedExternalNodes(2)  = node3;   connectedExternalNodes(3)  = node4;
  connectedExternalNodes(4)  = node5;   connectedExternalNodes(5)  = node6;
  connectedExternalNodes(6)  = node7;   connectedExternalNodes(7)  = node8;
  connectedExternalNodes(8)  = node9;   connectedExternalNodes(9)  = node10;
  connectedExternalNodes(10) = node11;  connectedExternalNodes(11) = node12;
  connectedExternalNodes(12) = node13;  connectedExternalNodes(13) = node14;
  connectedExternalNodes(14) = node15;  connectedExternalNodes(15) = node16;
  connectedExternalNodes(16) = node17;  connectedExternalNodes(17) = node18;
  connectedExternalNodes(18) = node19;  connectedExternalNodes(19) = node20;

  for (int i = 0; i < NEN; i++) nodePointers[i] = 0;

  // One material clone per formulation Gauss point (27 std / 8 uri); the
  // trailing capacity slots stay null. The factory pre-validates the 3D
  // capability (F2); this is the defensive backstop for direct-C++
  // construction — a failed clone marks the element dead (badGeom kill
  // switch) instead of exit(-1) (fork policy).  // Ladruno
  for (int i = 0; i < NGP; i++) { materialPointers[i] = 0; cachedRho[i] = 0.0; }
  for (int i = 0; i < this->nGP(); i++) {
    materialPointers[i] = theMaterial.getCopy("ThreeDimensional");
    if (materialPointers[i] == 0) {
      opserr << "LadrunoBrick20::constructor - element " << tag
             << ": failed to get a material of type ThreeDimensional (material "
             << theMaterial.getTag() << "); element marked DEAD - update() will "
                "fail and tangent/mass are zeroed.\n";
      badGeom = true;
    }
  }

  b[0] = b1; b[1] = b2; b[2] = b3;
  appliedB[0] = 0.0; appliedB[1] = 0.0; appliedB[2] = 0.0;

  for (int i = 0; i < NGP; i++) theDamping[i] = 0;
  if (damping) {
    for (int i = 0; i < this->nGP(); i++) {
      theDamping[i] = (*damping).getCopy();
      if (!theDamping[i])
        opserr << "LadrunoBrick20::LadrunoBrick20 -- failed to get copy of damping\n";
    }
  }
}

//destructor
LadrunoBrick20::~LadrunoBrick20()
{
  for (int i = 0; i < NGP; i++) {
    if (materialPointers[i]) delete materialPointers[i];
    if (theDamping[i])       delete theDamping[i];
  }
  if (load != 0) delete load;
  if (Ki != 0)   delete Ki;
  if (M0 != 0)   delete M0;
}

//set domain
void  LadrunoBrick20::setDomain(Domain *theDomain)
{
  // F9: set the base-class Domain pointer FIRST, before any early-return, so a
  // damping-init failure can never leave a stale Domain behind.  // Ladruno
  this->DomainComponent::setDomain(theDomain);

  geomCached = false;

  if (theDomain == 0) {                 // element removed from the domain
    for (int i = 0; i < NEN; i++) nodePointers[i] = 0;
    return;
  }

  for (int i = 0; i < NEN; i++)
    nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

  for (int i = 0; i < this->nGP(); i++) {
    if (theDamping[i] && theDamping[i]->setDomain(theDomain, NSTR)) {
      opserr << "LadrunoBrick20::setDomain -- Error initializing damping\n";
      return;
    }
  }

  // U1 advisory (ADR 72 §3.7): warn ONCE PER PROCESS when a softening /
  // crack-band material is attached to a quadratic element. Advisory only —
  // the run proceeds. Probe via the cached-Response pattern (the
  // LadrunoBrick::damageResponse construction), then discriminate on the
  // channel SHAPE: the crack-band concrete family reports DUAL-scalar damage
  // (ASDConcrete3D {d+,d-}, LadrunoConcrete3D {omega_t,omega_c} — size >= 2),
  // which is exactly the lch-regularized softening class the advisory targets.
  // Scalar damage channels (LadrunoJ2's Lemaitre ductile D, size 1 — exposed
  // even with the law off) and materials with no "damage" channel at all
  // (ElasticIsotropic) never fire — and never trip the process-once flag, so
  // an elastic mesh built before a concrete one cannot mask the advisory.  // Ladruno
  if (!advisedDamage && materialPointers[0] != 0) {
    DummyStream dmgStream;
    const char *dmgArgv[1] = {"damage"};
    Response *dmgResp = materialPointers[0]->setResponse(dmgArgv, 1, dmgStream);
    if (dmgResp != 0) {
      bool crackBand = false;
      if (dmgResp->getResponse() >= 0) {
        Information &info = dmgResp->getInformation();
        const Vector *d = info.theVector;
        crackBand = (d != 0 && d->Size() >= 2);
      }
      if (crackBand) {
        advisedDamage = true;           // process-once: trip ONLY when firing
        opserr << "LadrunoBrick20 (element " << this->getTag() << "): a softening / "
                  "crack-band material (dual-scalar \"damage\" channel) is attached "
                  "to a QUADRATIC element. Crack-band regularization is theory-shaky "
                  "on quadratic fields: lch = cbrt(V) may mis-size a band localizing "
                  "across a sub-cell (ADR 72 §3.7). Prefer the H8 family "
                  "(LadrunoBrick + ASDConcrete3D/LadrunoConcrete3D) for softening. "
                  "Proceeding (advisory only; printed once per process).\n";
      }
      delete dmgResp;
    }
  }

  // F11 geometry cache + F1 detJ gate (reference coords are constant, so the
  // setDomain-time validation covers every later call).  // Ladruno
  this->buildGeometryCache();

  // node set may have changed — drop the cached consistent mass
  // unconditionally, then seed the rho signature + hasMass (F14/F-1).  // Ladruno
  if (M0 != 0) { delete M0; M0 = 0; }
  this->refreshMassState();
}

//----------------------------------------------------------------------
// refreshMassState — F-1 (post-P2 adversarial gate). rho is a live NDMaterial
// parameter (e.g. ElasticIsotropic param ID 3) and a parameter update goes
// DIRECTLY to the material clones: LadrunoBrick20::setParameter forwards to
// the materials, which register THEMSELVES on the Parameter — the element's
// updateParameter never runs for it. A one-shot hasMass flag / M0 cache would
// therefore go stale (and the inertia RESIDUAL, which reads rho fresh, would
// silently disagree with the cached mass TANGENT). So every mass-path entry
// re-reads the mass-relevant densities — std: all 27 points; uri: point 0
// only, matching the mass integral's rho0 read — and drops M0 when the
// signature moved. Cost: nGP() virtual getRho() calls per entry, trivial
// against the 27-pt integration it guards.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::refreshMassState(void)
{
  const int n = (formulation == Formulation::URI) ? 1 : NGP;
  bool any = false, moved = false;
  for (int i = 0; i < n; i++) {
    const double r = (materialPointers[i] != 0) ? materialPointers[i]->getRho() : 0.0;
    if (r != 0.0) any = true;
    if (r != cachedRho[i]) { moved = true; cachedRho[i] = r; }
  }
  hasMass = any;
  if (moved && M0 != 0) { delete M0; M0 = 0; }
}

//----------------------------------------------------------------------
// buildGeometryCache — fill cachedN / cachedDNdx / cachedDV at the 27 GPs and
// validate detJ > 0 at every one (the upstream Twenty_Node_Brick::Jacobian3d
// per-GP gate). Threshold is detJ <= 0 (scale-free — never gate on an absolute
// |detJ| floor, it scales as L^3); NaN also fails the !(detJ > 0) test. On
// failure the element is marked bad: update() returns -1 and tangent/mass are
// emit-once + zeroed, so the analysis fails loudly.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::buildGeometryCache(void)
{
  geomCached = false;

  for (int i = 0; i < NEN; i++) {
    if (nodePointers[i] == 0) {
      opserr << "LadrunoBrick20::setDomain - element " << this->getTag()
             << ": node " << connectedExternalNodes(i)
             << " does not exist in the domain; element marked DEAD.\n";
      badGeom = true;
      return;
    }
  }

  double X[NEN][3];
  gatherCoords(X);

  double dN[NEN][3], J[3][3], Jinv[3][3];
  for (int L = 0; L < NGP; L++) {
    const double xi[3] = { Ladruno::hex20::GP27[L][0],
                           Ladruno::hex20::GP27[L][1],
                           Ladruno::hex20::GP27[L][2] };
    Ladruno::hex20::shape(xi, cachedN[L], dN);
    const double detJ = Ladruno::hex20::jacobian(X, dN, J);
    if (!(detJ > 0.0)) {
      opserr << "LadrunoBrick20::setDomain - element " << this->getTag()
             << ": Non-positive Jacobian: " << detJ
             << " at Gauss point " << L + 1 << " of " << NGP
             << " (xi = " << xi[0] << ", " << xi[1] << ", " << xi[2] << "). "
                "The element geometry is inverted or excessively distorted - "
                "check the node ordering and any curved / quarter-point midside "
                "nodes (a mid-edge node pushed past the quarter point makes detJ "
                "negative at the near corner). Element marked DEAD: update() "
                "will fail and tangent/mass are zeroed - the analysis will fail "
                "loudly rather than integrate a sign-flipped volume.\n";
      badGeom = true;
      return;
    }
    Ladruno::hex20::invert3(J, detJ, Jinv);
    Ladruno::hex20::cartGrad(dN, Jinv, cachedDNdx[L]);
    cachedDV[L] = Ladruno::hex20::GP27[L][3] * detJ;
  }

  // uri second cache: the 8 reduced (2x2x2 / Barlow) stiffness-residual points.
  // The 27-pt cache above is still built and gated — mass / body force /
  // volume always integrate the full rule (ADR 72 §2.3). The 8 points are
  // interior to the 27-pt hull, but gate them anyway (cheap, and NaN-safe).  // Ladruno
  if (formulation == Formulation::URI) {
    double N8[NEN];
    for (int L = 0; L < NGPU; L++) {
      const double xi[3] = { Ladruno::hex20::GP8[L][0],
                             Ladruno::hex20::GP8[L][1],
                             Ladruno::hex20::GP8[L][2] };
      Ladruno::hex20::shape(xi, N8, dN);
      const double detJ = Ladruno::hex20::jacobian(X, dN, J);
      if (!(detJ > 0.0)) {
        opserr << "LadrunoBrick20::setDomain - element " << this->getTag()
               << ": Non-positive Jacobian: " << detJ
               << " at reduced (uri) Gauss point " << L + 1 << " of " << NGPU
               << ". Element marked DEAD (see the 27-pt gate message class "
                  "above for the geometry diagnosis).\n";
        badGeom = true;
        return;
      }
      Ladruno::hex20::invert3(J, detJ, Jinv);
      Ladruno::hex20::cartGrad(dN, Jinv, cachedDNdx8[L]);
      cachedDV8[L] = Ladruno::hex20::GP8[L][3] * detJ;
    }
  }

  geomCached = true;
}

//----------------------------------------------------------------------
// cacheUsable — guard for every path that consumes the geometry cache. A bad /
// uncached element reports once (per instance) and the caller zeroes its
// output.  // Ladruno
//----------------------------------------------------------------------
bool
LadrunoBrick20::cacheUsable(const char *where)
{
  if (geomCached && !badGeom)
    return true;
  if (!warnedBadUse) {
    warnedBadUse = true;
    opserr << "LadrunoBrick20::" << where << " - element " << this->getTag()
           << (badGeom ? ": element is DEAD (non-positive Jacobian or failed "
                         "material clone - see the setDomain/constructor error "
                         "above); "
                       : ": geometry cache not built (setDomain has not run); ")
           << "returning zeros. (printed once per element)\n";
  }
  return false;
}

//----------------------------------------------------------------------
// coerceFormulationToStd — unknown-ordinal defense at the wire surface. STD
// and URI are both live (P2); an ordinal beyond the known set can only come
// from a corrupt or future-version stream — coerce to STD with a process-once
// notice instead of computing garbage.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::coerceFormulationToStd(const char *where)
{
  if (formulation == Formulation::STD || formulation == Formulation::URI)
    return;
  if (!warnedUriCoerce) {
    warnedUriCoerce = true;
    opserr << "WARNING LadrunoBrick20::" << where << " - element " << this->getTag()
           << ": unknown formulation ordinal "
           << static_cast<int>(formulation)
           << " (corrupt or future-version stream?) - coercing to std. "
              "(printed once per process)\n";
  }
  formulation = Formulation::STD;
}

int
LadrunoBrick20::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping) {
    for (int i = 0; i < this->nGP(); i++) {
      if (theDamping[i]) delete theDamping[i];
      theDamping[i] = (*damping).getCopy();
      if (!theDamping[i]) {
        opserr << "LadrunoBrick20::setDamping -- failed to get copy of damping\n";
        return -1;
      }
      if (theDamping[i]->setDomain(theDomain, NSTR)) {
        opserr << "LadrunoBrick20::setDamping -- Error initializing damping\n";
        return -2;
      }
    }
  }
  return 0;
}

int  LadrunoBrick20::getNumExternalNodes(void) const { return NEN; }
const ID &  LadrunoBrick20::getExternalNodes(void) { return connectedExternalNodes; }
Node **  LadrunoBrick20::getNodePtrs(void) { return nodePointers; }
int  LadrunoBrick20::getNumDOF(void) { return NDOF; }

//commit state (material pointers null-guarded: a dead element — failed clone —
//must not crash the domain commit)  // Ladruno
int  LadrunoBrick20::commitState(void)
{
  int success = 0;
  if ((success = this->Element::commitState()) != 0)
    opserr << "LadrunoBrick20::commitState () - failed in base class";
  for (int i = 0; i < this->nGP(); i++)
    if (materialPointers[i]) success += materialPointers[i]->commitState();
  for (int i = 0; i < this->nGP(); i++)
    if (theDamping[i]) success += theDamping[i]->commitState();
  return success;
}

int  LadrunoBrick20::revertToLastCommit(void)
{
  int success = 0;
  for (int i = 0; i < this->nGP(); i++)
    if (materialPointers[i]) success += materialPointers[i]->revertToLastCommit();
  for (int i = 0; i < this->nGP(); i++)
    if (theDamping[i]) success += theDamping[i]->revertToLastCommit();
  return success;
}

int  LadrunoBrick20::revertToStart(void)
{
  int success = 0;
  for (int i = 0; i < this->nGP(); i++)
    if (materialPointers[i]) success += materialPointers[i]->revertToStart();
  for (int i = 0; i < this->nGP(); i++)
    if (theDamping[i]) success += theDamping[i]->revertToStart();
  return success;
}

//print
void  LadrunoBrick20::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "LadrunoBrick20 (20-node serendipity hexahedron), formulation: "
      << formulationName(formulation) << "\n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Nodes: " << connectedExternalNodes;
    s << "Material Information : \n ";
    if (materialPointers[0]) materialPointers[0]->Print(s, flag);
    s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << endln;
    s << "Resisting Force (no inertia): " << this->getResistingForce();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoBrick20\", ";
    s << "\"formulation\": \"" << formulationName(formulation) << "\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
    for (int i = 1; i < NEN - 1; i++)
      s << connectedExternalNodes(i) << ", ";
    s << connectedExternalNodes(NEN - 1) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": \"" << (materialPointers[0] ? materialPointers[0]->getTag() : 0) << "\"}";
  }
}

//----------------------------------------------------------------------
// Gather the reference (undeformed) nodal coordinates into X[20][3].  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::gatherCoords(double X[NEN][3])
{
  for (int a = 0; a < NEN; a++) {
    const Vector &Xa = nodePointers[a]->getCrds();
    for (int i = 0; i < 3; i++)
      X[a][i] = Xa(i);
  }
}

//tangent stiffness
const Matrix &  LadrunoBrick20::getTangentStiff(void)
{
  formStiffness(0);
  return stiff;
}

//initial stiffness — same shared assembly, initial material tangent (F15)
const Matrix &  LadrunoBrick20::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  formStiffness(1);
  Ki = new Matrix(stiff);
  return *Ki;
}

//mass — consistent, integrated ONCE from the geometry cache and reused (F12);
//massless meshes early-out (F14); the momentum (inertia-residual) pass NEVER
//runs from here — it belongs to the residual path only.
const Matrix &  LadrunoBrick20::getMass(void)
{
  if (massType == 1 && !warnedLumped) {
    // P1 keeps the -lumped flag API-stable (parsed + serialized) but the HRZ
    // path lands P3. Fail loud-but-safe: error clearly (once per process), then
    // fall back to consistent so the run never uses an unbuilt lumped matrix.  // Ladruno
    warnedLumped = true;
    opserr << "ERROR LadrunoBrick20 (element " << this->getTag() << "): -lumped "
              "(HRZ lumped mass) is NOT yet implemented - it lands in P3 "
              "(ADR 72). Falling back to the CONSISTENT mass; results will not "
              "reflect a lumped-mass model. (printed once per process)\n";
  }

  this->refreshMassState();             // F-1: follow live rho parameter updates
  if (!hasMass || !cacheUsable("getMass")) {
    mass.Zero();
    return mass;
  }

  if (M0 == 0) {
    computeConsistentMass();
    M0 = new Matrix(mass);
  } else {
    mass = *M0;
  }
  return mass;
}

//----------------------------------------------------------------------
void  LadrunoBrick20::zeroLoad(void)
{
  if (load != 0) load->Zero();
  applyLoad = 0;
  appliedB[0] = 0.0; appliedB[1] = 0.0; appliedB[2] = 0.0;
}

int
LadrunoBrick20::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * b[0];
    appliedB[1] += loadFactor * b[1];
    appliedB[2] += loadFactor * b[2];
    return 0;
  } else if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    appliedB[2] += loadFactor * data(2) * b[2];
    return 0;
  } else {
    opserr << "LadrunoBrick20::addLoad() - ele with tag: " << this->getTag()
           << " does not deal with load type: " << type << "\n";
    return -1;
  }
}

int
LadrunoBrick20::addInertiaLoadToUnbalance(const Vector &accel)
{
  this->refreshMassState();             // F-1: follow live rho parameter updates
  if (!hasMass) return 0;               // F14 massless early-out
  if (!cacheUsable("addInertiaLoadToUnbalance")) return 0;

  // reuse the cached consistent mass (integrates once, F12)
  if (M0 == 0) {
    computeConsistentMass();
    M0 = new Matrix(mass);
  }

  int count = 0;
  for (int i = 0; i < NEN; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j = 0; j < NDF; j++)
      resid(count++) = Raccel(j);
  }

  if (load == 0) load = new Vector(NDOF);
  load->addMatrixVector(1.0, *M0, resid, -1.0);
  return 0;
}

//residual
const Vector &  LadrunoBrick20::getResistingForce(void)
{
  formResidual();
  if (load != 0) resid -= *load;
  return resid;
}

const Vector &  LadrunoBrick20::getResistingForceIncInertia(void)
{
  static Vector res(NDOF);

  formResidual();
  formInertiaResidual();

  res = resid;

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  if (load != 0) res -= *load;

  return res;
}

//----------------------------------------------------------------------
// formResidual — plain small-strain B^T sigma Gauss loop over the formulation
// rule (27-pt std / 8-pt uri), consuming the geometry cache (F11/F13: no
// kernel calls, no tangent work). The consistent body / self-weight integral
// ALWAYS runs the full 27-pt rule (a mass-class integral — formulation-
// independent like the mass, ADR 72 §2.3), in its own loop; each accumulator
// keeps the pre-P2 per-GP summation ORDER, so the std reduce-to residuals
// stay ~1e-15.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::formResidual(void)
{
  resid.Zero();

  if (!cacheUsable("getResistingForce"))
    return;

  static Vector stress(NSTR);
  static Vector dampingStress(NSTR);
  static Vector bodyForce(NDOF);

  const bool uri = (formulation == Formulation::URI);
  const int  ngp = this->nGP();

  for (int L = 0; L < ngp; L++) {
    const double dv = uri ? cachedDV8[L] : cachedDV[L];
    const double (*dNdx)[3] = uri ? cachedDNdx8[L] : cachedDNdx[L];

    stress = materialPointers[L]->getStress();

    if (theDamping[L]) {
      theDamping[L]->update(stress);
      dampingStress = theDamping[L]->getDampingForce();
      dampingStress *= dv;
    }

    stress *= dv;

    // residual: f_int += B^T sigma   (+ B^T dampingStress)
    for (int a = 0; a < NEN; a++) {
      const int c = NDF * a;
      const double bx = dNdx[a][0], by = dNdx[a][1], bz = dNdx[a][2];
      // B_a^T sigma  (sigma already scaled by dv)
      resid(c    ) += bx * stress(0) + by * stress(3) + bz * stress(5);
      resid(c + 1) += by * stress(1) + bx * stress(3) + bz * stress(4);
      resid(c + 2) += bz * stress(2) + by * stress(4) + bx * stress(5);
      if (theDamping[L]) {
        resid(c    ) += bx * dampingStress(0) + by * dampingStress(3) + bz * dampingStress(5);
        resid(c + 1) += by * dampingStress(1) + bx * dampingStress(3) + bz * dampingStress(4);
        resid(c + 2) += bz * dampingStress(2) + by * dampingStress(4) + bx * dampingStress(5);
      }
    }
  }

  // consistent body / self-weight load (subtracted from the residual) — full
  // 27-pt rule regardless of formulation (F13 hoists kept)
  const double bfx = (applyLoad == 0) ? b[0] : appliedB[0];
  const double bfy = (applyLoad == 0) ? b[1] : appliedB[1];
  const double bfz = (applyLoad == 0) ? b[2] : appliedB[2];
  if (bfx != 0.0 || bfy != 0.0 || bfz != 0.0) {
    bodyForce.Zero();
    for (int L = 0; L < NGP; L++) {
      const double dv = cachedDV[L];
      const double *N = cachedN[L];
      for (int a = 0; a < NEN; a++) {
        const int c = NDF * a;
        bodyForce(c    ) -= dv * bfx * N[a];
        bodyForce(c + 1) -= dv * bfy * N[a];
        bodyForce(c + 2) -= dv * bfz * N[a];
      }
    }
    resid += bodyForce;
  }
}

//----------------------------------------------------------------------
// formStiffness — the ONE B^T D B assembly (F15), shared by getTangentStiff
// (initialFlag=0, current material tangent) and getInitialStiff (initialFlag=1,
// initial tangent), over the formulation rule (27-pt std / 8-pt uri).
// BLOCKED per-node assembly (ADR 72 §6 P2 debt b): for each column node b,
// DB_b = D·B_b (6x3, B sparsity explicit), then each row node a accumulates
// the 3x3 block B_a^T·DB_b — no 6x60 temporaries, no zero-multiplies. No
// symmetry of D is assumed (some material tangents are unsymmetric). The
// summation reorder moves K entries only at the ~1e-16 relative level —
// re-anchored by the std reduce-to gate and the uri-vs-oracle gate.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::formStiffness(int initialFlag)
{
  stiff.Zero();

  if (!cacheUsable(initialFlag ? "getInitialStiff" : "getTangentStiff"))
    return;

  static Matrix dd(NSTR, NSTR);

  const bool uri = (formulation == Formulation::URI);
  const int  ngp = this->nGP();

  for (int L = 0; L < ngp; L++) {
    const double dv = uri ? cachedDV8[L] : cachedDV[L];
    const double (*dNdx)[3] = uri ? cachedDNdx8[L] : cachedDNdx[L];

    dd = initialFlag ? materialPointers[L]->getInitialTangent()
                     : materialPointers[L]->getTangent();
    if (theDamping[L]) dd *= theDamping[L]->getStiffnessMultiplier();
    dd *= dv;

    // stiff += B^T dd B   (dd already scaled by dv), node-pair 3x3 blocks.
    // Voigt rows {xx,yy,zz,xy,yz,zx}: B_b columns (u,v,w) carry
    //   u: {bx,0,0,by,0,bz}, v: {0,by,0,bx,bz,0}, w: {0,0,bz,0,by,bx}.
    for (int bnode = 0; bnode < NEN; bnode++) {
      const double bx = dNdx[bnode][0], by = dNdx[bnode][1], bz = dNdx[bnode][2];
      double DB[NSTR][3];
      for (int i = 0; i < NSTR; i++) {
        DB[i][0] = dd(i, 0) * bx + dd(i, 3) * by + dd(i, 5) * bz;
        DB[i][1] = dd(i, 1) * by + dd(i, 3) * bx + dd(i, 4) * bz;
        DB[i][2] = dd(i, 2) * bz + dd(i, 4) * by + dd(i, 5) * bx;
      }
      const int c = NDF * bnode;
      for (int anode = 0; anode < NEN; anode++) {
        const double ax = dNdx[anode][0], ay = dNdx[anode][1], az = dNdx[anode][2];
        const int r = NDF * anode;
        for (int j = 0; j < 3; j++) {
          stiff(r    , c + j) += ax * DB[0][j] + ay * DB[3][j] + az * DB[5][j];
          stiff(r + 1, c + j) += ay * DB[1][j] + ax * DB[3][j] + az * DB[4][j];
          stiff(r + 2, c + j) += az * DB[2][j] + ay * DB[4][j] + ax * DB[5][j];
        }
      }
    }
  }
}

//----------------------------------------------------------------------
// computeConsistentMass — consistent mass ALWAYS integrated with the full
// 27-pt rule (ADR 72 §2.3 / Abaqus AUG §28.1.1: mass uses full integration
// regardless of the stiffness formulation), from the geometry cache, into the
// static `mass` workspace. Integrated ONCE per element (getMass caches the
// result in M0, F12). massType==1 (HRZ lumped) lands P3; getMass errors and
// falls back here.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::computeConsistentMass(void)
{
  mass.Zero();

  // Under uri only 8 material clones exist but the mass rule has 27 points:
  // density is read from material point 0 (element-uniform — the 8 clones come
  // from one prototype, so this equals the std per-point read; a per-point rho
  // FIELD across the mass rule is not representable under uri).  // Ladruno
  const bool uri = (formulation == Formulation::URI);
  const double rho0 = uri ? materialPointers[0]->getRho() : 0.0;

  for (int L = 0; L < NGP; L++) {
    const double dv = cachedDV[L];
    const double *N = cachedN[L];
    const double rho = uri ? rho0 : materialPointers[L]->getRho();

    for (int a = 0; a < NEN; a++) {
      const double temp = N[a] * dv;
      const double tr = temp * rho;
      for (int bnode = 0; bnode < NEN; bnode++) {
        const double massJK = tr * N[bnode];
        for (int p = 0; p < NDF; p++)
          mass(NDF * a + p, NDF * bnode + p) += massJK;
      }
    }
  }
}

//----------------------------------------------------------------------
// formInertiaResidual — the momentum pass rho N_a N_b a_b, ACCUMULATED into
// `resid` on the residual path only (getResistingForceIncInertia); never runs
// from getMass (F12). Massless meshes early-out (F14).  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick20::formInertiaResidual(void)
{
  this->refreshMassState();             // F-1: follow live rho parameter updates
  if (!hasMass)
    return;
  if (!cacheUsable("getResistingForceIncInertia"))
    return;

  static Vector momentum(NDF);

  // 27-pt momentum pass regardless of formulation (pairs with the consistent
  // mass); uri reads rho from material point 0 — see computeConsistentMass.  // Ladruno
  const bool uri = (formulation == Formulation::URI);
  const double rho0 = uri ? materialPointers[0]->getRho() : 0.0;

  for (int L = 0; L < NGP; L++) {
    const double dv = cachedDV[L];
    const double *N = cachedN[L];
    const double rho = uri ? rho0 : materialPointers[L]->getRho();

    // momentum p = rho * sum_a N_a a_a
    momentum.Zero();
    for (int a = 0; a < NEN; a++)
      momentum.addVector(1.0, nodePointers[a]->getTrialAccel(), N[a]);
    momentum *= rho;

    for (int a = 0; a < NEN; a++) {
      const double temp = N[a] * dv;
      for (int p = 0; p < NDF; p++)
        resid(NDF * a + p) += temp * momentum(p);
    }
  }
}

//update — small-strain: set the material trial strain at every Gauss point
//(from the geometry cache). Returns -1 on a dead element so the analysis
//fails loudly (F1).
int
LadrunoBrick20::update(void)
{
  if (!cacheUsable("update"))
    return -1;

  static Vector strain(NSTR);
  int ret = 0;

  const bool uri = (formulation == Formulation::URI);
  const int  ngp = this->nGP();

  for (int L = 0; L < ngp; L++) {
    const double (*dNdx)[3] = uri ? cachedDNdx8[L] : cachedDNdx[L];

    // strain = B u  (Voigt {xx,yy,zz,xy,yz,zx}, engineering shear)
    strain.Zero();
    for (int a = 0; a < NEN; a++) {
      const Vector &ua = nodePointers[a]->getTrialDisp();
      const double bx = dNdx[a][0], by = dNdx[a][1], bz = dNdx[a][2];
      const double u0 = ua(0), u1 = ua(1), u2 = ua(2);
      strain(0) += bx * u0;
      strain(1) += by * u1;
      strain(2) += bz * u2;
      strain(3) += by * u0 + bx * u1;
      strain(4) += bz * u1 + by * u2;
      strain(5) += bz * u0 + bx * u2;
    }
    ret += materialPointers[L]->setTrialStrain(strain);
  }
  return ret;
}

//----------------------------------------------------------------------
// Reference-config element volume: V = sum of the cached w_L |J|_L (F11).
// 0 on a dead/uncached element (lch then falls back to the Element base).  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick20::computeVolume(void)
{
  if (!geomCached || badGeom)
    return 0.0;
  double V = 0.0;
  for (int L = 0; L < NGP; L++)
    V += cachedDV[L];
  return V;
}

double
LadrunoBrick20::getCharacteristicLength(void)
{
  double V = this->computeVolume();
  if (V <= 0.0)
    return Element::getCharacteristicLength();   // degenerate: fall back
  return cbrt(V);
}

//----------------------------------------------------------------------
int
LadrunoBrick20::displaySelf(Renderer &theViewer, int displayMode, float fact,
                            const char **modes, int numMode)
{
  // Corner-brick render: draw the 8 corner nodes as a hexahedron (the mid-edge
  // nodes are not part of the wireframe cube). Mirrors LadrunoBrick's drawCube.
  static Matrix coords(8, 3);
  static Vector values(8);
  static Vector vN(3);

  for (int n = 0; n < 8; n++) {
    nodePointers[n]->getDisplayCrds(vN, fact, displayMode);
    for (int i = 0; i < 3; i++)
      coords(n, i) = vN(i);
  }

  if (displayMode < 3 && displayMode > 0) {
    int index = displayMode - 1;
    for (int n = 0; n < 8; n++) {
      const Vector &stressN = materialPointers[n]->getStress();
      values(n) = stressN(index);
    }
  } else {
    for (int n = 0; n < 8; n++)
      values(n) = 0.0;
  }

  return theViewer.drawCube(coords, values, this->getTag());
}

//----------------------------------------------------------------------
Response *
LadrunoBrick20::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  char outputData[32];

  // ── LadrunoRecorder geometry self-description (contract Part A) ──────────
  // One class tag carries TWO GP layouts (std 27-pt brcshl / uri 8-pt
  // lexicographic), so the recorder cannot dispatch on the class tag alone:
  // it probes basisInfo and keys the integration rule on numGP (ADR 72 §6
  // P2 debt a). Answered BEFORE the ElementOutput tag, BezierTet10 pattern.  // Ladruno
  if (strcmp(argv[0], "basisInfo") == 0) {
    output.tag("ElementBasis");
    output.attr("topology",    "hex");
    output.attr("family",      "serendipity");
    output.attr("paramDomain", "[-1,1]");
    output.attr("rational",    0);
    output.attr("numCtrl",     NEN);          // 20 nodes
    output.attr("numGP",       this->nGP());  // 27 std / 8 uri result stations
    output.attr("orderU",      2);            // quadratic serendipity
    output.endTag();                          // ElementBasis
    return new ElementResponse(this, 101, ID(1));   // non-null sentinel
  }
  if (strcmp(argv[0], "integrationPoints") == 0)
    return new ElementResponse(this, 102, Matrix(this->nGP(), 3));
  if (strcmp(argv[0], "integrationWeights") == 0)
    return new ElementResponse(this, 103, Vector(this->nGP()));

  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoBrick20");
  output.attr("eleTag", this->getTag());
  for (int i = 1; i <= NEN; i++) {
    sprintf(outputData, "node%d", i);
    output.attr(outputData, nodePointers[i - 1]->getTag());
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
    for (int i = 1; i <= NEN; i++) {
      sprintf(outputData, "P1_%d", i); output.tag("ResponseType", outputData);
      sprintf(outputData, "P2_%d", i); output.tag("ResponseType", outputData);
      sprintf(outputData, "P3_%d", i); output.tag("ResponseType", outputData);
    }
    theResponse = new ElementResponse(this, 1, resid);

  } else if (strcmp(argv[0], "stiff") == 0 || strcmp(argv[0], "stiffness") == 0) {
    theResponse = new ElementResponse(this, 2, Matrix(NDOF, NDOF));

  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    if (argc >= 2) {                      // F5: guard before reading argv[1]
      int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= this->nGP()) {
        output.tag("GaussPoint");
        output.attr("number", pointNum);
        theResponse = materialPointers[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
        output.endTag();
      }
    }

  } else if (strcmp(argv[0], "stresses") == 0) {
    for (int i = 0; i < this->nGP(); i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma33");
      output.tag("ResponseType", "sigma12");
      output.tag("ResponseType", "sigma23");
      output.tag("ResponseType", "sigma13");
      output.endTag();
      output.endTag();
    }
    theResponse = new ElementResponse(this, 3, Vector(this->nGP() * NSTR));

  } else if (strcmp(argv[0], "strains") == 0) {
    for (int i = 0; i < this->nGP(); i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      output.tag("ResponseType", "eps11");
      output.tag("ResponseType", "eps22");
      output.tag("ResponseType", "eps33");
      output.tag("ResponseType", "eps12");
      output.tag("ResponseType", "eps23");
      output.tag("ResponseType", "eps13");
      output.endTag();
      output.endTag();
    }
    theResponse = new ElementResponse(this, 4, Vector(this->nGP() * NSTR));

  } else if (strcmp(argv[0], "stress3D6") == 0) {
    output.tag("GaussPoint");
    output.attr("number", 1);
    output.tag("NdMaterialOutput");
    output.attr("classType", materialPointers[0]->getClassTag());
    output.attr("tag", materialPointers[0]->getTag());
    output.tag("ResponseType", "sigma11");
    output.tag("ResponseType", "sigma22");
    output.tag("ResponseType", "sigma33");
    output.tag("ResponseType", "sigma12");
    output.tag("ResponseType", "sigma23");
    output.tag("ResponseType", "sigma13");
    output.endTag();
    output.endTag();
    theResponse = new ElementResponse(this, 6, Vector(NSTR));

  } else if (strcmp(argv[0], "strain3D6") == 0) {
    output.tag("GaussPoint");
    output.attr("number", 1);
    output.tag("NdMaterialOutput");
    output.attr("classType", materialPointers[0]->getClassTag());
    output.attr("tag", materialPointers[0]->getTag());
    output.tag("ResponseType", "eps11");
    output.tag("ResponseType", "eps22");
    output.tag("ResponseType", "eps33");
    output.tag("ResponseType", "eps12");
    output.tag("ResponseType", "eps23");
    output.tag("ResponseType", "eps13");
    output.endTag();
    output.endTag();
    theResponse = new ElementResponse(this, 7, Vector(NSTR));

  } else if (strcmp(argv[0], "charLength") == 0 ||
             strcmp(argv[0], "characteristicLength") == 0) {
    output.tag("ResponseType", "lch");
    theResponse = new ElementResponse(this, 9, Vector(1));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}

int
LadrunoBrick20::getResponse(int responseID, Information &eleInfo)
{
  const int ngp = this->nGP();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());

  else if (responseID == 3) {
    Vector stresses(ngp * NSTR);        // sized by formulation (27 std / 8 uri)
    int cnt = 0;
    for (int i = 0; i < ngp; i++) {
      const Vector &sigma = materialPointers[i]->getStress();
      for (int j = 0; j < NSTR; j++) stresses(cnt++) = sigma(j);
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 4) {
    Vector strains(ngp * NSTR);
    int cnt = 0;
    for (int i = 0; i < ngp; i++) {
      const Vector &eps = materialPointers[i]->getStrain();
      for (int j = 0; j < NSTR; j++) strains(cnt++) = eps(j);
    }
    return eleInfo.setVector(strains);

  } else if (responseID == 6) {
    Vector tmpStress(NSTR);
    for (int i = 0; i < ngp; i++) {
      const Vector &sigma = materialPointers[i]->getStress();
      for (int j = 0; j < NSTR; j++) tmpStress(j) += sigma(j);
    }
    tmpStress /= (double)ngp;
    return eleInfo.setVector(tmpStress);

  } else if (responseID == 7) {
    Vector tmpStrain(NSTR);
    for (int i = 0; i < ngp; i++) {
      const Vector &eps = materialPointers[i]->getStrain();
      for (int j = 0; j < NSTR; j++) tmpStrain(j) += eps(j);
    }
    tmpStrain /= (double)ngp;
    return eleInfo.setVector(tmpStrain);

  } else if (responseID == 9) {
    static Vector lch(1);
    lch(0) = this->getCharacteristicLength();
    return eleInfo.setVector(lch);

  // ── LadrunoRecorder geometry probes (contract Part A) ────────────────────
  } else if (responseID == 101) {
    // basisInfo sentinel — metadata was emitted via the stream already.
    return 0;

  } else if (responseID == 102) {
    // integrationPoints: natural (xi,eta,zeta) of the formulation GP set.
    Matrix pts(ngp, 3);
    const double (*gp)[4] = (formulation == Formulation::URI)
                              ? Ladruno::hex20::GP8 : Ladruno::hex20::GP27;
    for (int g = 0; g < ngp; g++)
      for (int j = 0; j < 3; j++)
        pts(g, j) = gp[g][j];
    return eleInfo.setMatrix(pts);

  } else if (responseID == 103) {
    // integrationWeights, same GP order (sum = 8 = reference-cube volume).
    Vector w(ngp);
    const double (*gp)[4] = (formulation == Formulation::URI)
                              ? Ladruno::hex20::GP8 : Ladruno::hex20::GP27;
    for (int g = 0; g < ngp; g++)
      w(g) = gp[g][3];
    return eleInfo.setVector(w);
  }

  return -1;
}

int
LadrunoBrick20::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  if (argc < 1) return -1;

  // damping (loop bounds from nGP() — no upstream i<4 bug class)
  if (strstr(argv[0], "damp") != 0) {
    if (argc < 2) return -1;
    for (int i = 0; i < this->nGP(); i++) {
      if (theDamping[i]) {
        int dmpRes = theDamping[i]->setParameter(argv, argc, param);
        if (dmpRes != -1) res = dmpRes;
      }
    }
    return res;
  }

  // specific material point
  if (strstr(argv[0], "material") != 0) {
    if (argc < 3) return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= this->nGP())
      return materialPointers[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  // all material points
  for (int i = 0; i < this->nGP(); i++) {
    int matRes = materialPointers[i]->setParameter(argv, argc, param);
    if (matRes != -1) res = matRes;
  }
  return res;
}

//F7: honest aggregation — any GP material failing (<0) fails the element
//update; mirrors the setParameter conventions.
int
LadrunoBrick20::updateParameter(int parameterID, Information &info)
{
  if (parameterID == -1)
    return -1;

  int res = 0;
  for (int i = 0; i < this->nGP(); i++) {
    if (materialPointers[i] == 0)
      continue;
    int matRes = materialPointers[i]->updateParameter(parameterID, info);
    if (matRes < 0)
      res = -1;
  }
  return res;
}

//----------------------------------------------------------------------
// sendSelf / recvSelf — mirror the LadrunoBrick layout, sized 20 nodes + 27
// material SLOTS (+ optional Damping). ID layout (append-only serialization):
//   [0 .. NGP-1]            27 material class-tag slots (uri fills 8, rest 0)
//   [NGP .. 2*NGP-1]        27 material db-tag slots
//   [2*NGP .. 2*NGP+NEN-1]  20 node tags
//   [2*NGP+NEN + 0]         element tag
//   [2*NGP+NEN + 1]         Rayleigh-damping flag
//   [2*NGP+NEN + 2]         Damping class tag (0 = none)
//   [2*NGP+NEN + 3]         Damping db tag
//   [2*NGP+NEN + 4]         packed: formulation ordinal + 10*massType
// Double data Vector(7): alphaM, betaK, betaK0, betaKc, b[0..2].
// The SLOT layout is formulation-independent (wire compatible); only nGP()
// materials are actually streamed, and recvSelf decodes the ordinal BEFORE
// the material loop so the count always matches the sender's.  // Ladruno
//----------------------------------------------------------------------
int  LadrunoBrick20::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();
  int matDbTag;

  const int idBase = 2 * NGP + NEN;   // = 74
  static ID idData(2 * NGP + NEN + 5);

  idData(idBase + 0) = this->getTag();
  idData(idBase + 1) = (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0) ? 1 : 0;

  // zero ALL material slots first (static ID is reused across elements — a
  // uri element must not leak the previous element's trailing slots)  // Ladruno
  for (int i = 0; i < NGP; i++) {
    idData(i) = 0;
    idData(NGP + i) = 0;
  }
  for (int i = 0; i < this->nGP(); i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(NGP + i) = matDbTag;
  }

  for (int i = 0; i < NEN; i++)
    idData(2 * NGP + i) = connectedExternalNodes(i);

  idData(idBase + 2) = 0;
  idData(idBase + 3) = 0;
  if (theDamping[0]) {
    idData(idBase + 2) = theDamping[0]->getClassTag();
    int dbTag = theDamping[0]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        for (int i = 0; i < this->nGP(); i++)
          theDamping[i]->setDbTag(dbTag);
    }
    idData(idBase + 3) = dbTag;
  }

  idData(idBase + 4) = static_cast<int>(formulation) + 10 * massType;

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoBrick20::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector dData(7);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];
  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoBrick20::sendSelf() - failed to send double data\n";
    return -1;
  }

  for (int i = 0; i < this->nGP(); i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LadrunoBrick20::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  if (theDamping[0]) {
    for (int i = 0; i < this->nGP(); i++) {
      res += theDamping[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoBrick20::sendSelf -- could not send Damping\n";
        return res;
      }
    }
  }

  return res;
}

int  LadrunoBrick20::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  const int idBase = 2 * NGP + NEN;
  static ID idData(2 * NGP + NEN + 5);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoBrick20::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(idBase + 0));

  static Vector dData(7);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoBrick20::recvSelf() - failed to recv double data\n";
    return -1;
  }
  alphaM = dData(0);
  betaK  = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  b[0]   = dData(4);
  b[1]   = dData(5);
  b[2]   = dData(6);

  for (int i = 0; i < NEN; i++)
    connectedExternalNodes(i) = idData(2 * NGP + i);

  int packed = idData(idBase + 4);
  formulation = static_cast<Formulation>(packed % 10);
  massType    = (packed / 10) % 10;

  // Unknown-ordinal wire defense (STD and URI are both live since P2; only a
  // corrupt / future-version ordinal coerces). Decoded BEFORE the material
  // loop so nGP() below matches the sender's stream.  // Ladruno
  this->coerceFormulationToStd("recvSelf");

  // geometry / formulation may have changed — the cache is rebuilt by the
  // setDomain that follows recvSelf (Domain::addElement), and the cached mass
  // and initial stiffness are stale.  // Ladruno
  geomCached = false;
  badGeom = false;
  warnedBadUse = false;
  if (M0 != 0) { delete M0; M0 = 0; }
  if (Ki != 0) { delete Ki; Ki = 0; }

  for (int i = 0; i < this->nGP(); i++) {
    int matClassTag = idData(i);
    int matDbTag = idData(NGP + i);
    if (materialPointers[i] == 0 ||
        materialPointers[i]->getClassTag() != matClassTag) {
      if (materialPointers[i]) delete materialPointers[i];
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
        opserr << "LadrunoBrick20::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
        return -1;
      }
      materialPointers[i]->setDbTag(matDbTag);
    }
    res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "LadrunoBrick20::recvSelf() - material " << i << " failed to recv itself\n";
      return res;
    }
  }
  // a reused object may hold trailing clones from a previous std life —
  // drop anything past the received formulation's count  // Ladruno
  for (int i = this->nGP(); i < NGP; i++) {
    if (materialPointers[i]) { delete materialPointers[i]; materialPointers[i] = 0; }
  }

  int dmpTag = (int)idData(idBase + 2);
  if (dmpTag) {
    for (int i = 0; i < this->nGP(); i++) {
      if (theDamping[i] == 0 || theDamping[i]->getClassTag() != dmpTag) {
        if (theDamping[i]) delete theDamping[i];
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "LadrunoBrick20::recvSelf -- could not get a Damping\n";
          return -1;
        }
      }
      theDamping[i]->setDbTag((int)idData(idBase + 3));
      res += theDamping[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoBrick20::recvSelf -- could not receive Damping\n";
        return res;
      }
    }
    for (int i = this->nGP(); i < NGP; i++) {
      if (theDamping[i]) { delete theDamping[i]; theDamping[i] = 0; }
    }
  } else {
    for (int i = 0; i < NGP; i++) {
      if (theDamping[i]) {
        delete theDamping[i];
        theDamping[i] = 0;
      }
    }
  }

  return res;
}
