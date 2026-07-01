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

// ADR-39 P1b — LadrunoContactDomain implementation. See header for design.

#include "LadrunoContactDomain.h"
#include "LadrunoContactSurface.h"
#include <OPS_Globals.h>   // opserr, endln
#include <cmath>           // sqrt
#include <algorithm>       // std::min, std::max (C2.2 Uzawa clamp)
#include <Domain.h>        // ADR-60: needsResort() looks up committed slave coords by tag
#include <Node.h>          // ADR-60
#include <Vector.h>        // ADR-60
#include <Subdomain.h>     // ADR-60 R6: serial-only refusal — detect a worker Subdomain host

LadrunoContactDomain::LadrunoContactDomain()
  : numCommits(0), numReverts(0)
{
}

LadrunoContactDomain::~LadrunoContactDomain()
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        delete theSurfaces[i];
    theSurfaces.clear();
    theContacts.clear();
}

int
LadrunoContactDomain::addSurface(LadrunoContactSurface *surf)
{
    if (surf == 0)
        return -1;
    // reject a duplicate tag
    if (getSurface(surf->getTag()) != 0) {
        opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
               << surf->getTag() << " already exists\n";
        return -1;
    }
    theSurfaces.push_back(surf);
    return 0;
}

LadrunoContactSurface *
LadrunoContactDomain::getSurface(int tag) const
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        if (theSurfaces[i]->getTag() == tag)
            return theSurfaces[i];
    return 0;
}

int
LadrunoContactDomain::addContact(int tag, int masterSurfTag, int slaveSurfTag,
                                 double kn, double kt, double mu, const double *outward,
                                 bool knAuto, double cellFrac, bool consistentTan, double muc,
                                 bool consistentNormal, double softScale,
                                 bool enableReemit, double resortFrac, int resortEvery,
                                 bool smoothNormal)
{
    if (getSurface(masterSurfTag) == 0 || getSurface(slaveSurfTag) == 0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    if (kn < 0.0) {
        // kn < 0 = attractive contact + negative-definite tangent (unstable). Reject
        // at this choke point (covers Py + Tcl), mirroring addRigidPlane's guard.
        // kn == 0 is allowed (the P1b zero-force topology path); the segment handler
        // warns + skips an inert kn<=0 SEGMENT contact. (Gate H1.)
        opserr << "WARNING LadrunoContactDomain::addContact() - penalty kn must be >= 0 (got "
               << kn << ")\n";
        return -1;
    }
    Contact c;
    c.tag = tag; c.masterSurfTag = masterSurfTag; c.slaveSurfTag = slaveSurfTag;
    c.kn = kn; c.kt = kt; c.mu = mu;
    c.knAuto = knAuto;
    c.hasOutward = (outward != 0);
    for (int d = 0; d < 3; d++) c.outward[d] = (outward != 0) ? outward[d] : 0.0;
    c.cellFrac = (cellFrac > 0.0) ? cellFrac : 1.0;   // P2.5 broad-phase cell scale
    c.consistentTan = consistentTan;                  // P3.5 friction tangent symmetry
    c.muc = (muc > 0.0) ? muc : 0.0;                  // D2 viscous stabilization (≤0 ⇒ off)
    c.consistentNormal = consistentNormal;            // B3 ∂n/∂u geometric normal tangent
    c.softScale = (softScale > 0.0) ? softScale : 0.0; // B1 SOFT=1 scale (≤0 ⇒ off, byte-identical)
    c.enableReemit = enableReemit;                     // ADR-60 finite-sliding re-emit (off ⇒ identical)
    c.resortFrac   = (resortFrac > 0.0) ? resortFrac : 0.5;
    c.resortEvery  = (resortEvery > 0) ? resortEvery : 0;
    c.smoothNormal = smoothNormal;                     // ADR-63 #4a nodal-normal smoothing (off ⇒ identical)
    theContacts.push_back(c);
    return 0;
}

// --- ADR-60: finite-sliding NTS re-emit (anchors + migration trigger). The handler registers
//     committed slave positions each handle() for opted-in contacts; Domain::commit() polls
//     needsResort() and forces a domainChange() when a slave has migrated past the band. ---
void
LadrunoContactDomain::clearReemit(void)
{
    theReemit.clear();   // handle() begin: drop the previous epoch's anchors + triggers
}

void
LadrunoContactDomain::setReemitContact(int contactTag, double band, double resortFrac, int resortEvery)
{
    ReemitContact rc;
    rc.contactTag = contactTag;
    rc.band = band;
    theReemit.push_back(rc);
    // R4/R7: ensure a PERSISTENT Trigger for this contact (one per contactTag, survives the
    // per-handle theReemit rebuild). The migration floor is fixed at the D1 default (10);
    // -resortEvery is now the FORCED cycle cadence (forceEvery = LS-DYNA BSORT), NOT the floor
    // it was previously conflated with. A first handle constructs it; later handles only refresh
    // its config (setConfig) so the accumulated cadence/arming is preserved (the R4 fix).
    double f          = (resortFrac > 0.0) ? resortFrac : 0.5;
    int    forceEvery = (resortEvery > 0) ? resortEvery : 0;
    std::map<int, LadrunoContactReemit::Trigger>::iterator it = theReemitTrig.find(contactTag);
    if (it == theReemitTrig.end())
        theReemitTrig[contactTag] = LadrunoContactReemit::Trigger(f, 0.5, 10, forceEvery);
    else
        it->second.setConfig(f, 0.5, 10, forceEvery);
}

void
LadrunoContactDomain::addReemitAnchor(int contactTag, int slaveTag, const double anchor[3])
{
    for (size_t c = 0; c < theReemit.size(); c++)
        if (theReemit[c].contactTag == contactTag) {
            ReemitAnchor a;
            a.slaveTag = slaveTag;
            a.x[0] = anchor[0]; a.x[1] = anchor[1]; a.x[2] = anchor[2];
            theReemit[c].anchors.push_back(a);
            return;
        }
}

bool
LadrunoContactDomain::reemitMembershipChanged(int contactTag, unsigned long long fingerprint)
{
    std::map<int, unsigned long long>::iterator it = theReemitFp.find(contactTag);
    if (it == theReemitFp.end()) {        // first handle for this contact — nothing to alias yet
        theReemitFp[contactTag] = fingerprint;
        return false;
    }
    if (it->second == fingerprint) return false;
    it->second = fingerprint;             // record the new ordering; report the change once
    return true;
}

void
LadrunoContactDomain::dropFrictionForContact(int contactTag)
{
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ) {
        if (it->first.c == contactTag)
            theFrictionStates.erase(it++);    // post-increment keeps the iterator valid past erase
        else
            ++it;
    }
}

// --- ADR-63 #4a: averaged nodal-normal smoothing field (per master surface). ---
void
LadrunoContactDomain::clearNormalFields(void)
{
    // Drop the per-handle scattered normals but KEEP the topological σ/fp cache (it survives across
    // handles; setNormalField recomputes σ only on a fingerprint change). A removed contact's stale
    // entry then carries no field ⇒ getSegNodalNorm returns 0 ⇒ harmless.
    for (std::map<int, NormalField>::iterator it = theNormalFields.begin();
         it != theNormalFields.end(); ++it)
        it->second.segNodalNorm.clear();
}

int
LadrunoContactDomain::setNormalField(int masterSurfTag, int nps, const int *mTags, int nSeg,
                                     const double *segCoords, const double globalSeed[3])
{
    NormalField &nf = theNormalFields[masterSurfTag];
    unsigned long long fp =
        LadrunoContactReemit::membershipFingerprint(mTags, nSeg * nps, nps);
    // (re)compute the COHERENT WINDING only when the membership changed — it is topological, so the
    // cache survives across handles and a re-mesh / element removal (new fingerprint) recomputes it
    // (and re-checks the refuse conditions: a re-mesh can change manifold-ness). R1 composes here.
    if (nf.sigma.empty() || nf.fp != fp || nf.nSeg != nSeg || nf.nps != nps) {
        nf.sigma.assign((size_t)(nSeg > 0 ? nSeg : 1), 0);
        nf.status = LadrunoContactNormalField::propagateOrientation(
            mTags, nSeg, nps, nf.sigma.empty() ? 0 : &nf.sigma[0]);
        nf.fp = fp; nf.nSeg = nSeg; nf.nps = nps;
        nf.signCaptured = false;            // membership changed ⇒ re-capture the frozen sign
    }
    if (nf.status != LadrunoContactNormalField::OK) {
        nf.segNodalNorm.clear();            // REFUSED ⇒ no field (the adapter keeps the faceted path)
        return nf.status;
    }
    // ADR-63 F1/D2 — capture the GLOBAL outward sign ONCE (first OK build) and FREEZE it. Re-deciding
    // it from the (deformed) current config every handle could flip the whole field mid-run (the R3
    // failure on the temporal axis). The nodal-normal magnitudes/directions still track deformation;
    // only the scalar sign is frozen. signConf flags a near coin-flip (seed ~⟂ field) for the handler.
    if (!nf.signCaptured) {
        nf.globalSign = LadrunoContactNormalField::voteSign(nSeg, nps, segCoords, &nf.sigma[0],
                                                            globalSeed, &nf.signConf);
        nf.signCaptured = true;
    }
    // rebuild the per-handle nodal-normal field from the current coords + the FROZEN sign
    nf.segNodalNorm.assign((size_t)nSeg * nps * 3, 0.0);
    LadrunoContactNormalField::nodalNormals(mTags, nSeg, nps, segCoords,
                                            &nf.sigma[0], nf.globalSign, &nf.segNodalNorm[0]);
    return LadrunoContactNormalField::OK;
}

double
LadrunoContactDomain::getNormalFieldSignConf(int masterSurfTag) const
{
    std::map<int, NormalField>::const_iterator it = theNormalFields.find(masterSurfTag);
    if (it == theNormalFields.end() || it->second.status != LadrunoContactNormalField::OK) return -1.0;
    return it->second.signConf;
}

const double *
LadrunoContactDomain::getSegNodalNorm(int masterSurfTag, int segIndex) const
{
    std::map<int, NormalField>::const_iterator it = theNormalFields.find(masterSurfTag);
    if (it == theNormalFields.end()) return 0;
    const NormalField &nf = it->second;
    if (nf.status != LadrunoContactNormalField::OK || nf.segNodalNorm.empty()) return 0;
    if (segIndex < 0 || segIndex >= nf.nSeg) return 0;
    return &nf.segNodalNorm[(size_t)segIndex * nf.nps * 3];
}

bool
LadrunoContactDomain::needsResort(Domain *theDomain)
{
    // O(1)-false when no contact opted in (the byte-identity contract). The Domain::commit()
    // call site also gates on !contactAugmenting (held-load Uzawa must not re-handle, gate GA-1).
    if (theReemit.empty() || theDomain == 0) return false;
    // R6 (D7 serial-only refusal): under a partitioned host this Domain is itself a worker
    // Subdomain — the migration trigger would fire uncoordinated on every rank, desyncing the
    // parallel domainChange. Refuse here too (the handler already skips anchor registration on a
    // partitioned host, so theReemit is normally empty; this is the belt-and-suspenders guard for
    // a worker that somehow registered). Sequential ⇒ the cast is null ⇒ inert (byte-identical).
    if (dynamic_cast<Subdomain *>(theDomain) != 0) return false;
    bool resort = false;
    for (size_t c = 0; c < theReemit.size(); c++) {
        ReemitContact &rc = theReemit[c];
        double dmax = 0.0;
        for (size_t i = 0; i < rc.anchors.size(); i++) {
            Node *nd = theDomain->getNode(rc.anchors[i].slaveTag);
            if (nd == 0) continue;
            const Vector &X = nd->getCrds();
            const Vector &u = nd->getDisp();   // committed disp (Domain::commit() committed it first)
            double cur[3] = { X(0) + u(0), X(1) + u(1), X(2) + u(2) };
            double d = LadrunoContactReemit::maxMigration(1, rc.anchors[i].x, cur);
            if (d > dmax) dmax = d;
        }
        // tick EVERY contact's PERSISTENT trigger each commit (cadence advances); OR the fires.
        // The trigger is keyed by contactTag in theReemitTrig (R4) — setReemitContact created it
        // this handle, so the lookup is present; guard defensively anyway.
        std::map<int, LadrunoContactReemit::Trigger>::iterator ti = theReemitTrig.find(rc.contactTag);
        if (ti != theReemitTrig.end() && ti->second.update(dmax, rc.band)) resort = true;
    }
    return resort;
}

int
LadrunoContactDomain::addMortarContact(int tag, int masterSurfTag, int slaveSurfTag,
                                       double kn, bool knAuto, double epsN, bool epsNAuto,
                                       double augTol, int maxAug, int ngp,
                                       const double *outward, double cellFrac,
                                       double mu, double epsT, bool epsTAuto,
                                       double cohesion, double tauMax, bool consistentTan,
                                       bool isTie, double muc, double softScale,
                                       bool edgeEdge, double edgeKn, bool edgeKnAuto, double edgeBand,
                                       double edgeMu, double edgeKt, double edgeCohesion,
                                       double edgeTauMax, bool edgeConsistentTan, double edgeSoftScale,
                                       bool edgeAlm, double edgeAugTol)
{
    LadrunoContactSurface *ms = getSurface(masterSurfTag);
    LadrunoContactSurface *ss = getSurface(slaveSurfTag);
    if (ms == 0 || ss == 0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    // The mortar lane integrates D over slave FACETS and M against master FACETS, so both
    // surfaces must be faceted (a SLAVE_NODES node-set carries no facet ⇒ no D). Reject the
    // wrong kinds at this choke point (covers Py + Tcl), mirroring addRigidPlane's guard.
    if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - master surface "
               << masterSurfTag << " must be -master (MASTER_SEGMENTS)\n";
        return -1;
    }
    if (ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - slave surface "
               << slaveSurfTag << " must be -slave-segments (SLAVE_SEGMENTS); the node-set "
                  "-slave kind carries no facet, so it cannot give the mortar matrix D\n";
        return -1;
    }
    if (kn < 0.0 || epsN < 0.0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - kn/epsN must be >= 0 "
                  "(got kn " << kn << ", epsN " << epsN << ")\n";
        return -1;
    }
    // C4 — a tie is an equality bond: it has no friction cone. Refuse the combination here too (the
    // command surface refuses it first, but this choke point also covers a direct API call).
    if (isTie && (mu > 0.0 || cohesion > 0.0 || tauMax > 0.0)) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - a mesh-tie (-tie) has no "
                  "friction cone; -mu/-cohesion/-tauMax are not allowed with -tie\n";
        return -1;
    }
    if (isTie && muc > 0.0) {
        // D2.2 — a tie is a permanent bond (no contact-status flips ⇒ no chatter to damp). Refuse
        // -visc with -tie at this choke point too (the command surface refuses it first).
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - -visc (viscous contact "
                  "stabilization) is not allowed with -tie (a bond has no contact-chatter regime)\n";
        return -1;
    }
    if (softScale > 0.0 && isTie) {
        // B2 (P5) — the SOFT=2 segment-based explicit penalty supports -visc (D2.2 normal damper) AND
        // -mu/-cohesion/-tauMax (Courant-stable Coulomb/Tresca over the segment overlap) under explicit;
        // implicit falls through to the regular mortar penalty/friction. Only -tie is refused at this
        // choke point too (the command surface refuses it first) — a permanent bond is not a soft penalty.
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - -soft (SOFT=2 segment-based "
                  "explicit penalty) is not allowed with -tie (a permanent bond is not a soft penalty)\n";
        return -1;
    }
    MortarContact m;
    m.tag = tag; m.masterSurfTag = masterSurfTag; m.slaveSurfTag = slaveSurfTag;
    m.kn = kn; m.knAuto = knAuto;
    m.epsN = epsN; m.epsNAuto = epsNAuto;
    m.augTol = (augTol > 0.0) ? augTol : 1e-8;
    m.maxAug = (maxAug > 0) ? maxAug : 20;
    m.ngp = (ngp >= 1) ? ngp : 2;                     // slave-facet Gauss rule order
    m.hasOutward = (outward != 0);
    for (int d = 0; d < 3; d++) m.outward[d] = (outward != 0) ? outward[d] : 0.0;
    m.cellFrac = (cellFrac > 0.0) ? cellFrac : 1.0;
    m.mu = (mu > 0.0) ? mu : 0.0;                     // C3.1 friction (≤0 ⇒ frictionless)
    m.epsT = epsT; m.epsTAuto = epsTAuto;
    m.cohesion = (cohesion > 0.0) ? cohesion : 0.0;
    m.tauMax = tauMax;                                // ≤0 ⇒ no Tresca upper cap
    m.consistentTan = consistentTan;
    m.isTie = isTie;                                  // C4 — permanent mesh-tie bond
    m.muc = (muc > 0.0) ? muc : 0.0;                  // D2.2 viscous stabilization (≤0 ⇒ off)
    m.softScale = (softScale > 0.0) ? softScale : 0.0;  // B2 (P5) SOFT=2 segment-based explicit penalty
    m.edgeEdge = edgeEdge;                             // ADR-57 E2 edge-edge fallback (off ⇒ byte-identical)
    m.edgeKn = (edgeKn > 0.0) ? edgeKn : 0.0;          // ≤0 ⇒ default to the resolved mortar penalty
    m.edgeKnAuto = edgeKnAuto;                         // size the edge penalty per master facet
    m.edgeBand = (edgeBand > 0.0) ? edgeBand : 0.0;    // ≤0 ⇒ default from the facet edge length
    m.edgeMu = (edgeMu > 0.0) ? edgeMu : 0.0;          // ADR-57 E3 edge-edge friction (≤0 ⇒ frictionless)
    m.edgeKt = (edgeKt > 0.0) ? edgeKt : 0.0;          // ≤0 ⇒ default to the edge normal penalty
    m.edgeCohesion = (edgeCohesion > 0.0) ? edgeCohesion : 0.0;
    m.edgeTauMax = edgeTauMax;                         // ≤0 ⇒ no Tresca upper cap
    m.edgeConsistentTan = edgeConsistentTan;
    m.edgeSoftScale = (edgeSoftScale > 0.0) ? edgeSoftScale : 0.0;  // ADR-57 E5 explicit SOFT (≤0 ⇒ off)
    m.edgeAlm = edgeAlm;                                // ADR-57 E6 one-scalar ALM (off ⇒ byte-identical)
    m.edgeAugTol = (edgeAugTol > 0.0) ? edgeAugTol : 1e-8;
    theMortarContacts.push_back(m);
    return 0;
}

int
LadrunoContactDomain::addRigidPlane(int tag, int slaveSurfTag,
                                    const double p0[3], const double n[3], double kn, double muc,
                                    double softScale)
{
    LadrunoContactSurface *surf = getSurface(slaveSurfTag);
    if (surf == 0) {
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - slave surface "
               << slaveSurfTag << " not defined\n";
        return -1;
    }
    if (surf->getKind() != LadrunoContactSurface::SLAVE_NODES) {
        // A MASTER_SEGMENTS surface stores flat, repeated connectivity; treating it
        // as a slave node-set would double-count contact on shared nodes. Require a
        // genuine slave node-set (the handler iterates getNodeTags() as slave nodes).
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - surface "
               << slaveSurfTag << " is not a slave node-set\n";
        return -1;
    }
    double nrm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nrm < 1e-14) {
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - zero normal\n";
        return -1;
    }
    if (kn <= 0.0) {
        // kn < 0 would make contact attractive with a negative-definite tangent
        // (unstable); kn == 0 is silently inert. Reject both at this choke point
        // (covers the Python + Tcl parsers) — mirrors the zero-normal guard above.
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - penalty kn must be "
               << "> 0 (got " << kn << ")\n";
        return -1;
    }
    RigidPlane rp;
    rp.tag = tag; rp.slaveSurfTag = slaveSurfTag; rp.kn = kn;
    rp.muc = (muc > 0.0) ? muc : 0.0;                 // D2 viscous stabilization (≤0 ⇒ off)
    rp.softScale = (softScale > 0.0) ? softScale : 0.0; // B1 SOFT=1 scale (≤0 ⇒ off, byte-identical)
    for (int d = 0; d < 3; d++) { rp.p0[d] = p0[d]; rp.n[d] = n[d] / nrm; }
    theRigidPlanes.push_back(rp);
    return 0;
}

LadrunoContactDomain::FrictionState &
LadrunoContactDomain::getOrCreateFrictionState(int contactTag, int slaveTag, int segIndex)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    // operator[] default-constructs a fresh FrictionState (engaged=false, zeroed)
    // when the key is absent — exactly the lazy-create-at-first-activation contract.
    return theFrictionStates[k];
}

// ===================================================== ADR-41 C2.2 normal ALM (λ_N)
LadrunoContactDomain::MortarNormalState &
LadrunoContactDomain::getOrCreateMortarNormalState(int contactTag, int slaveNodeTag)
{
    NodeKey k; k.c = contactTag; k.n = slaveNodeTag;
    // operator[] default-constructs a zeroed MortarNormalState (λ_N = 0) when absent.
    return theMortarNormalStates[k];
}

void
LadrunoContactDomain::accumulateMortarGap(int contactTag, int slaveNodeTag, int feTag,
                                          double gtFacet, double aFacet, double epsN)
{
    // delta update so a re-eval OVERWRITES (idempotent — the residual may be formed several
    // times per Newton iterate, and twice on the CDL first step). The node's gtGlobal/aGlobal
    // stays the running Σ over the facets this node touches (oracle T7(b): == the global g̃,a).
    PairKey fk; fk.c = contactTag; fk.s = slaveNodeTag; fk.g = feTag;
    FacetContrib &fc = theMortarFacetContribs[fk];
    NodeKey nk; nk.c = contactTag; nk.n = slaveNodeTag;
    MortarNormalState &st = theMortarNormalStates[nk];
    st.gtGlobal += gtFacet - fc.gt;
    st.aGlobal += aFacet - fc.a;
    // epsN for the commit-time Uzawa update at this node. For a fixed `-epsN val` every facet
    // writes the same value. Under `-epsN auto` the penalty is sized PER MASTER FACET, so a
    // shared slave node sees several; take the MAX (stiffest) — order-INDEPENDENT (a plain
    // overwrite would be last-writer-wins, i.e. facet-eval-order dependent). epsN is reset to 0
    // each handle() (mortarNormalGCEnd) so this max is per-analysis, never a cross-rebuild leak.
    if (epsN > st.epsN) st.epsN = epsN;
    fc.gt = gtFacet; fc.a = aFacet;
}

double
LadrunoContactDomain::getMaxMortarPenetration(void) const
{
    double pen = 0.0;
    for (std::map<NodeKey, MortarNormalState>::const_iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        const MortarNormalState &st = it->second;
        // C4: a TIE slot is not a contact node — it carries no normal gap (gtGlobal/lambdaN stay 0,
        // so it is already inert here), but skip it EXPLICITLY to mirror getMaxMortarTieResidual's
        // !isTie guard and keep this query robust if a future track ever writes those fields.
        if (st.isTie) continue;
        if (st.aGlobal <= 1e-300) continue;
        double gbar = st.gtGlobal / st.aGlobal;          // ḡ_I; < 0 ⇒ penetration
        // KKT-active only: a node held open (λ + epsN·ḡ ≥ 0) is NOT a contact violation even
        // if its raw ḡ < 0 (it would be clamped to zero pressure), so it does not count here.
        // NOTE: this uses the GLOBAL gap (the augmentation measure), whereas the per-facet
        // force/tangent active mask uses each facet's LOCAL gap; the two coincide at equilibrium
        // (where this query is read — after the converged commit), so the augTol stop criterion
        // is evaluated consistently with the converged force.
        if (st.lambdaN + st.epsN * gbar < 0.0 && -gbar > pen)
            pen = -gbar;
    }
    return pen;
}

// ===================================================== ADR-41 C4 mesh-tying (λ_tie)
void
LadrunoContactDomain::accumulateMortarTie(int contactTag, int slaveNodeTag, int feTag,
                                          const double rFacet[3], double aFacet, double epsTie)
{
    // delta update so a re-eval OVERWRITES (idempotent — the residual may be formed several times
    // per Newton iterate, and twice on the CDL first step). The node's rtGlobal/aGlobal stays the
    // running Σ over the facets this node touches. r_I is a LINEAR accumulation, so this global sum
    // is order-INDEPENDENT (the C4 shared-node resolution; unlike the friction slip MAJOR-1).
    PairKey fk; fk.c = contactTag; fk.s = slaveNodeTag; fk.g = feTag;
    FacetContrib &fc = theMortarFacetContribs[fk];
    NodeKey nk; nk.c = contactTag; nk.n = slaveNodeTag;
    MortarNormalState &st = theMortarNormalStates[nk];
    for (int d = 0; d < 3; d++) {
        st.rtGlobal[d] += rFacet[d] - fc.r[d];
        fc.r[d] = rFacet[d];
    }
    st.aGlobal += aFacet - fc.a;
    fc.a = aFacet;
    st.isTie = true;
    // epsTie for the commit-time Uzawa at this node — MAX over facets (order-independent under
    // `-epsTie auto`, which sizes per master facet), reusing the epsN slot (mutually exclusive use).
    if (epsTie > st.epsN) st.epsN = epsTie;
}

double
LadrunoContactDomain::getMaxMortarTieResidual(void) const
{
    double res = 0.0;
    for (std::map<NodeKey, MortarNormalState>::const_iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        const MortarNormalState &st = it->second;
        if (!st.isTie || st.aGlobal <= 1e-300) continue;
        for (int d = 0; d < 3; d++) {
            double rbar = std::fabs(st.rtGlobal[d] / st.aGlobal);   // |r̄_I,d|; → 0 under the tie Uzawa
            if (rbar > res) res = rbar;
        }
    }
    return res;
}

void
LadrunoContactDomain::mortarNormalGCBegin(void)
{
    liveNodeKeys.clear();
}

void
LadrunoContactDomain::mortarNormalGCMark(int contactTag, int slaveNodeTag)
{
    NodeKey k; k.c = contactTag; k.n = slaveNodeTag;
    liveNodeKeys.insert(k);
}

void
LadrunoContactDomain::mortarNormalGCEnd(void)
{
    // prune λ_N slots no live mortar pair referenced this handle() (the re-mesh leak class).
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ) {
        if (liveNodeKeys.find(it->first) == liveNodeKeys.end())
            theMortarNormalStates.erase(it++);
        else
            ++it;
    }
    liveNodeKeys.clear();
    // FE tags are reassigned each handle(), so the facet-contribution keys are NOT stable —
    // drop them all and zero the (transient) running gap sums on the survivors. The adapters
    // re-fill gtGlobal/aGlobal on the next getResidual sweep (delta from 0); only λ_N persists.
    theMortarFacetContribs.clear();
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        it->second.gtGlobal = 0.0;
        it->second.aGlobal = 0.0;
        it->second.epsN = 0.0;     // re-established (max over facets) on the next residual sweep
        // C4 — the transient tie accumulator is also re-filled each sweep (only lambdaTie persists).
        for (int d = 0; d < 3; d++) it->second.rtGlobal[d] = 0.0;
    }
}

// ===================================================== ADR-57 E2 edge-edge state
// canonical key: order each edge's two node tags so a facet edge shared by two facets (and a
// pair re-discovered with the endpoints swapped) maps to ONE slot (de-dup; never a lossy hash).
LadrunoContactDomain::EdgeEdgeState &
LadrunoContactDomain::getOrCreateEdgeEdgeState(int contactTag, int sNodeA, int sNodeB,
                                               int mNodeA, int mNodeB)
{
    EdgeKey k;
    k.c  = contactTag;
    k.sa = (sNodeA < sNodeB) ? sNodeA : sNodeB;
    k.sb = (sNodeA < sNodeB) ? sNodeB : sNodeA;
    k.ma = (mNodeA < mNodeB) ? mNodeA : mNodeB;
    k.mb = (mNodeA < mNodeB) ? mNodeB : mNodeA;
    // operator[] default-constructs a zeroed EdgeEdgeState (signN=0) when absent.
    return theEdgeEdgeStates[k];
}

void
LadrunoContactDomain::edgeGCBegin(void)
{
    liveEdgeKeys.clear();
}

void
LadrunoContactDomain::edgeGCMark(int contactTag, int sNodeA, int sNodeB, int mNodeA, int mNodeB)
{
    EdgeKey k;
    k.c  = contactTag;
    k.sa = (sNodeA < sNodeB) ? sNodeA : sNodeB;
    k.sb = (sNodeA < sNodeB) ? sNodeB : sNodeA;
    k.ma = (mNodeA < mNodeB) ? mNodeA : mNodeB;
    k.mb = (mNodeA < mNodeB) ? mNodeB : mNodeA;
    liveEdgeKeys.insert(k);
}

void
LadrunoContactDomain::edgeGCEnd(void)
{
    // prune edge-edge slots no live pair referenced this handle() (the re-mesh leak class).
    for (std::map<EdgeKey, EdgeEdgeState>::iterator it = theEdgeEdgeStates.begin();
         it != theEdgeEdgeStates.end(); ) {
        if (liveEdgeKeys.find(it->first) == liveEdgeKeys.end())
            theEdgeEdgeStates.erase(it++);
        else
            ++it;
    }
    liveEdgeKeys.clear();
    // E6 — drop the transient ALM Uzawa inputs on the survivors (gN_committed/εN are re-written each
    // getResidual when -edgeAlm is on; only signN + the friction slip persist). Mirrors
    // mortarNormalGCEnd zeroing gtGlobal/aGlobal so a re-paired slot's first commit cannot read a
    // stale gap/penalty from a prior epoch (the gate MAJOR-1 fix; λ_N persists — committed-only).
    for (std::map<EdgeKey, EdgeEdgeState>::iterator it = theEdgeEdgeStates.begin();
         it != theEdgeEdgeStates.end(); ++it) {
        it->second.gN_committed = 0.0;
        it->second.epsN = 0.0;
    }
}

// ADR-57 E6 — ‖gN‖_∞ over KKT-active edge-edge ALM pairs (the point-like getMaxMortarPenetration).
double
LadrunoContactDomain::getMaxEdgePenetration(void) const
{
    double pen = 0.0;
    for (std::map<EdgeKey, EdgeEdgeState>::const_iterator it = theEdgeEdgeStates.begin();
         it != theEdgeEdgeStates.end(); ++it) {
        const EdgeEdgeState &st = it->second;
        if (st.epsN <= 0.0) continue;                // not an ALM pair (-edgeAlm off ⇒ epsN never written)
        // KKT-active only: a pair held open (λ + εN·gN ≥ 0) is not a contact violation even if its raw
        // gN < 0 (clamped to zero pressure), mirroring getMaxMortarPenetration's active mask.
        if (st.lambdaN + st.epsN * st.gN_committed < 0.0 && -st.gN_committed > pen)
            pen = -st.gN_committed;
    }
    return pen;
}

void
LadrunoContactDomain::frictionGCBegin(void)
{
    liveKeys.clear();
    theNtsForce.clear();    // B3: drop last handle()'s NTS force snapshot (keys rebuilt this handle)
}

// B3 (P2b-2c) — NTS contact-force snapshot for the `ladrunoContactForce` query.
void
LadrunoContactDomain::setNtsForce(int contactTag, int slaveTag, int segIndex, double tn)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    theNtsForce[k] = tn;        // overwrite ⇒ holds the committed force after convergence
}

double
LadrunoContactDomain::getNtsForce(int slaveTag) const
{
    double f = 0.0;             // Σ over all (contactTag, slaveTag, segIndex) pairs of this node
    for (std::map<PairKey, double>::const_iterator it = theNtsForce.begin();
         it != theNtsForce.end(); ++it)
        if (it->first.s == slaveTag) f += it->second;
    return f;
}

// B1 (P4) — assembled nodal-mass cache (the SOFT=1 gap-mode effective mass; see the header).
void
LadrunoContactDomain::clearNodalMass(void)
{
    theNodalMass.clear();
}

void
LadrunoContactDomain::setNodalMass(int nodeTag, const double m[3])
{
    NodeMass nm; for (int d = 0; d < 3; d++) nm.m[d] = m[d];
    theNodalMass[nodeTag] = nm;
}

bool
LadrunoContactDomain::getNodalMass(int nodeTag, double m[3]) const
{
    std::map<int, NodeMass>::const_iterator it = theNodalMass.find(nodeTag);
    if (it == theNodalMass.end()) return false;
    for (int d = 0; d < 3; d++) m[d] = it->second.m[d];
    return true;
}

void
LadrunoContactDomain::frictionGCMark(int contactTag, int slaveTag, int segIndex)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    liveKeys.insert(k);
}

void
LadrunoContactDomain::frictionGCEnd(void)
{
    // erase any slot whose key the current handle() did not mark live (a re-meshed
    // or re-paired analysis leaves them behind otherwise). Marking a not-yet-created
    // key is harmless; this only prunes EXISTING slots.
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ) {
        if (liveKeys.find(it->first) == liveKeys.end())
            theFrictionStates.erase(it++);     // post-increment keeps the iterator valid
        else
            ++it;
    }
    liveKeys.clear();
}

int
LadrunoContactDomain::commit(void)
{
    // P3: promote the trial plastic slip to committed for every friction slot. The
    // counter (kept from P1b) lets a test confirm the Domain::commit() hook fires.
    numCommits++;
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ++it)
        for (int d = 0; d < 3; d++) it->second.gpT[d] = it->second.gpTtrial[d];

    // C2.2 — ONE Uzawa augmentation per commit (the EmbeddedRebar::commitState precedent):
    // λ_I ← min(0, λ_I + epsN·ḡ_I), ḡ_I = g̃_I^global / a_I^global (the just-converged trial
    // gap, accumulated across facets by the adapters this step). Across load steps this
    // augments for free on stock Newton; a held-load `analyzeAugmented` proc re-commits at
    // zero increment to drive ‖ḡ‖→augTol within a step. λ_N is mutated ONLY here, so the
    // revertToLastCommit path leaves it untouched (the committed-only invariant).
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        MortarNormalState &st = it->second;
        // C4 — mesh-tie Uzawa: λ_tie ← λ_tie + epsTie·(r_I/a_I), NO clamp (equality constraint).
        // r_I is the just-converged GLOBAL weighted relative displacement (order-independent, the
        // λ_N accumulator pattern). Drives ‖r‖ → 0 at FINITE epsTie (epsTie-INDEPENDENT bond — the
        // headline ALM win); the held-load `analyzeAugmented` proc augments it within a step. λ_tie
        // is mutated ONLY here ⇒ revertToLastCommit leaves it untouched (the committed-only invariant).
        // A tie slot has no normal KKT / friction state, so skip the rest of the loop body.
        if (st.isTie) {
            if (st.aGlobal > 1e-300)
                for (int d = 0; d < 3; d++)
                    st.lambdaTie[d] += st.epsN * (st.rtGlobal[d] / st.aGlobal);
            continue;
        }
        // C3.1 — promote the trial tangential slip (penalty friction; the EmbeddedRebar/NTS
        // FrictionState precedent). Done for every slot regardless of normal activity.
        for (int d = 0; d < 3; d++) st.gpT[d] = st.gpTtrial[d];
        // C3.3 — one tangential Uzawa step per commit: λ_T ← the returned cone-capped traction
        // (= lambdaTtrial = −tFric, written by the adapter). Drives the stick elastic creep → 0
        // at finite epsT (the EmbeddedRebar precedent, the λ_N analogue). λ_T≡0 path is inert.
        for (int d = 0; d < 3; d++) st.lambdaT[d] = st.lambdaTtrial[d];
        // C3.2 (MAJOR-2) — commit the engagement origin so a later rejected step can revert it.
        for (int d = 0; d < 3; d++) st.gT0committed[d] = st.gT0[d];
        st.engagedCommitted = st.engaged;
        if (st.aGlobal <= 1e-300) continue;          // unreferenced this step (no normal Uzawa)
        double gbar = st.gtGlobal / st.aGlobal;
        st.lambdaN = std::min(0.0, st.lambdaN + st.epsN * gbar);
    }

    // ADR-57 E2/E3/E6 — edge-edge slots. The shipped commit loop iterates only friction + mortar slots,
    // so iterating theEdgeEdgeStates here is the EXPLICIT obligation on the edge-edge PR (capstone
    // contract #1, the same trap the mortar lane had — NOT inherited). Promote the body-fixed sign
    // committed double-buffer (E2 — so a later rejected step can revert the capture) + the friction
    // trial (E3). E6 — ONE Uzawa augmentation per commit (the MortarNormalState C2.2 precedent, point-
    // like ⇒ ONE scalar gN, no shared-node accumulator): λ_N ← min(0, λ_N + εN·gN_committed). The
    // adapter writes gN_committed + εN ONLY when -edgeAlm is on; with ALM off both stay 0 ⇒ the update
    // is min(0, λ_N) and λ_N never leaves 0 ⇒ byte-identical to the E2 penalty path. λ_N is mutated
    // ONLY here ⇒ revertToLastCommit leaves it untouched (the committed-only invariant).
    for (std::map<EdgeKey, EdgeEdgeState>::iterator it = theEdgeEdgeStates.begin();
         it != theEdgeEdgeStates.end(); ++it) {
        EdgeEdgeState &st = it->second;
        st.signNcommitted = st.signN;                // body-fixed sign capture (E2)
        for (int d = 0; d < 3; d++) {                // friction (E3)
            st.gpT[d] = st.gpTtrial[d];
            st.gT0committed[d] = st.gT0[d];
        }
        st.engagedCommitted = st.engaged;
        st.lambdaN = std::min(0.0, st.lambdaN + st.epsN * st.gN_committed);   // E6 one-scalar Uzawa
    }
    return 0;
}

int
LadrunoContactDomain::revertToLastCommit(void)
{
    // P3: drop the trial back to the last committed plastic slip (a failed implicit
    // step, or a Python adaptive-step retry that calls Domain::revertToLastCommit()).
    numReverts++;
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ++it)
        for (int d = 0; d < 3; d++) it->second.gpTtrial[d] = it->second.gpT[d];
    // C3.1 — mortar friction: drop the trial slip back to committed (λ_N is committed-only ⇒
    // untouched on revert, the C2.2 invariant; only the friction trial needs reverting).
    // C3.2 (MAJOR-2) — also restore the engagement origin gT0/engaged from the committed copy, so
    // a rejected implicit step does not latch a stale origin captured at the rejected config.
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        MortarNormalState &st = it->second;
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = st.gpT[d];
            st.gT0[d] = st.gT0committed[d];
            st.lambdaTtrial[d] = st.lambdaT[d];        // C3.3 — drop the trial multiplier
        }
        st.engaged = st.engagedCommitted;
    }
    // ADR-57 E2 — edge-edge slots: restore the sign capture + the inert friction trial from the
    // committed double-buffer, so a rejected implicit step does not latch a sign/origin captured at
    // the rejected config (the C3.2 MAJOR-2 pattern). signN is mutated in getResidual (the capture),
    // so it MUST be revertible (unreachable under explicit CDL, live under implicit Newton).
    for (std::map<EdgeKey, EdgeEdgeState>::iterator it = theEdgeEdgeStates.begin();
         it != theEdgeEdgeStates.end(); ++it) {
        EdgeEdgeState &st = it->second;
        st.signN = st.signNcommitted;
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = st.gpT[d];
            st.gT0[d] = st.gT0committed[d];
        }
        st.engaged = st.engagedCommitted;
    }
    return 0;
}
