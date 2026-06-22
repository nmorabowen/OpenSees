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
    // a contact interaction between a master and a slave surface (+ P2 params).
    // outward (optional, null = auto): a direction toward the slave's allowed
    // half-space, used to orient the derived segment normal (design-gate BLOCKER-1);
    // null = the handler auto-derives it from the slave-vs-segment geometry.
    // knAuto: if true, kn is auto-sized per (slave,segment) pair at handle() time
    // from the owning solid element's stiffness (P2b-2b). The passed kn is then a
    // placeholder (0) and ignored; the handler computes the real value.
    int addContact(int tag, int masterSurfTag, int slaveSurfTag,
                   double kn, double kt, double mu, const double *outward = 0,
                   bool knAuto = false);
    int getNumContacts(void) const { return (int)theContacts.size(); }

    // --- P2b: faceted node-to-segment penalty contact. A Contact references a
    //     MASTER_SEGMENTS surface + a SLAVE_NODES surface; the handler builds one
    //     LadrunoContactFE per (slave node, master segment) pair. Exposed so the
    //     handler (OPS_Analysis) can read the definition. ---
    struct Contact {
        int tag, masterSurfTag, slaveSurfTag;
        double kn, kt, mu;
        bool knAuto;            // true => kn auto-sized from the master element (P2b-2b)
        bool hasOutward;        // true if an explicit orientation direction was given
        double outward[3];      // orientation direction toward the allowed half-space
    };
    const Contact &getContact(int i) const { return theContacts[i]; }

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
    std::vector<LadrunoContactSurface *> theSurfaces;   // owned
    std::vector<Contact> theContacts;
    std::vector<RigidPlane> theRigidPlanes;             // P2a
    int numCommits;
    int numReverts;
};

#endif
