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
                                 bool knAuto, double cellFrac, bool consistentTan)
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
    theContacts.push_back(c);
    return 0;
}

int
LadrunoContactDomain::addMortarContact(int tag, int masterSurfTag, int slaveSurfTag,
                                       double kn, bool knAuto, double epsN, bool epsNAuto,
                                       double augTol, int maxAug, int ngp,
                                       const double *outward, double cellFrac,
                                       double mu, double epsT, bool epsTAuto,
                                       double cohesion, double tauMax, bool consistentTan)
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
    theMortarContacts.push_back(m);
    return 0;
}

int
LadrunoContactDomain::addRigidPlane(int tag, int slaveSurfTag,
                                    const double p0[3], const double n[3], double kn)
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
    }
}

void
LadrunoContactDomain::frictionGCBegin(void)
{
    liveKeys.clear();
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
        // C3.1 — promote the trial tangential slip (penalty friction; the EmbeddedRebar/NTS
        // FrictionState precedent). Done for every slot regardless of normal activity.
        for (int d = 0; d < 3; d++) st.gpT[d] = st.gpTtrial[d];
        if (st.aGlobal <= 1e-300) continue;          // unreferenced this step (no normal Uzawa)
        double gbar = st.gtGlobal / st.aGlobal;
        st.lambdaN = std::min(0.0, st.lambdaN + st.epsN * gbar);
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
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it)
        for (int d = 0; d < 3; d++) it->second.gpTtrial[d] = it->second.gpT[d];
    return 0;
}
