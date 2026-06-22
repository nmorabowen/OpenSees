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

// ADR-39 P1b — LadrunoContactDomain: the contact engine attached to Domain,
// parallel to it. Owns the contact surfaces + contact-pair DEFINITIONS and (P2+)
// the path-dependent pair STATE (gap0, friction). It is the single, integrator-
// agnostic home for contact state so the assembly-side adapters (LadrunoContactFE,
// created fresh every handle()) stay stateless VIEWS.
//
// Lifetime (design gate): a Domain owns at most one (nullable -> zero cost,
// byte-identical to stock). It SURVIVES domainChanged (which runs
// AnalysisModel::clearAll, not Domain::clearAll) and is deleted on Domain::clearAll
// (the ops.wipe path) + ~Domain. commit()/revertToLastCommit() are driven by the
// matching Domain hooks so path-dependent state commits/reverts integrator-
// agnostically (BLOCKER-1/2 from the design gate).
//
// P1b scope: skeleton + lifecycle plumbing, ZERO force (the narrow phase is P2).
// buildAdapterCount() reports how many (empty-connectivity, zero) adapters the
// handler should inject — one per contact definition.
//
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p1b_design.md.

#ifndef LadrunoContactDomain_h
#define LadrunoContactDomain_h

#include <vector>

class LadrunoContactSurface;
class ID;

class LadrunoContactDomain
{
  public:
    LadrunoContactDomain();
    ~LadrunoContactDomain();

    // --- model definition (from the contactSurface / contact commands) ---
    int addSurface(LadrunoContactSurface *surf);          // takes ownership
    LadrunoContactSurface *getSurface(int tag) const;
    // a contact interaction between a master and a slave surface (+ P2 params)
    int addContact(int tag, int masterSurfTag, int slaveSurfTag,
                   double kn, double kt, double mu);
    int getNumContacts(void) const { return (int)theContacts.size(); }

    // --- P2a: rigid analytical plane (point p0 + outward unit normal n) vs a
    //     slave node-set; one adapter per slave node, connectivity = {slave}. ---
    struct RigidPlane {
        int tag, slaveSurfTag;
        double p0[3], n[3];     // n stored normalized at add time
        double kn;
    };
    int addRigidPlane(int tag, int slaveSurfTag,
                      const double p0[3], const double n[3], double kn);
    int getNumRigidPlanes(void) const { return (int)theRigidPlanes.size(); }
    const RigidPlane &getRigidPlane(int i) const { return theRigidPlanes[i]; }

    // --- handler interface: how many adapters to inject this handle() ---
    int buildAdapterCount(void) const { return (int)theContacts.size(); }

    // --- lifecycle (driven by Domain::commit / revertToLastCommit) ---
    int commit(void);              // P1b: no-op on state (counter for tests)
    int revertToLastCommit(void);  // P1b: no-op on state (counter for tests)
    int getNumCommits(void) const { return numCommits; }
    int getNumReverts(void) const { return numReverts; }

  private:
    struct Contact {
        int tag, masterSurfTag, slaveSurfTag;
        double kn, kt, mu;
    };
    std::vector<LadrunoContactSurface *> theSurfaces;   // owned
    std::vector<Contact> theContacts;
    std::vector<RigidPlane> theRigidPlanes;             // P2a
    int numCommits;
    int numReverts;
};

#endif
