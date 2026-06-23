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
#include <map>
#include <set>

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
                   bool knAuto = false, double cellFrac = 1.0,
                   bool consistentTan = false);
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
        double cellFrac;        // P2.5 bucket-sort cell = cellFrac * median seg diag
                                // (1.0 default; a huge value => 1 bucket => brute force)
        bool consistentTan;     // P3.5: true => non-symmetric consistent friction tangent
                                // (quadratic, needs FullGeneral/UmfPack); false (default)
                                // => the symmetric tangent (solver-safe on any system)
    };
    const Contact &getContact(int i) const { return theContacts[i]; }

    // --- ADR-41 C2: a MORTAR contact references a SLAVE_SEGMENTS surface + a
    //     MASTER_SEGMENTS surface; the handler (C2.1) builds a clipped-GP mortar
    //     adapter per (slave facet, master facet) pair and the frictionless
    //     commit-cycle ALM (C2.2) augments lambda_N here. C2.0 only STORES the
    //     definition (a separate vector the NTS handler loop does not read), so the
    //     command surface exists and is byte-identical/inert until C2.1 consumes it. ---
    struct MortarContact {
        int tag, masterSurfTag, slaveSurfTag;
        double kn;              // penalty / epsN seed; knAuto => sized from the solid (C2.1)
        bool   knAuto;
        double epsN;            // ALM normal penalty (C2.2); epsNAuto => = kn auto-size
        bool   epsNAuto;
        double augTol;          // Uzawa augmentation tolerance on ||g~||_inf (C2.2)
        int    maxAug;          // max augmentations per step (C2.2)
        int    ngp;             // slave-facet Gauss rule order (default 2 = 3-pt)
        bool   hasOutward;
        double outward[3];
        double cellFrac;        // broad-phase bucket-sort cell scale (shared with NTS)
        // ADR-41 C3.1 — friction (Coulomb/Tresca on the unified cone min(μN+c, τmax)). All ≤ 0 ⇒
        // the frictionless C2 path. epsT = the tangential penalty (epsTAuto ⇒ = the normal penalty).
        double mu;              // Coulomb friction coefficient
        double epsT;            // tangential penalty
        bool   epsTAuto;        // size epsT from the normal penalty at handle() time
        double cohesion;        // adhesive intercept c
        double tauMax;          // Tresca shear cap (≤0 ⇒ no upper cap)
        bool   consistentTan;   // C3.2: non-symmetric Coulomb friction tangent
    };
    int addMortarContact(int tag, int masterSurfTag, int slaveSurfTag,
                         double kn, bool knAuto, double epsN, bool epsNAuto,
                         double augTol, int maxAug, int ngp,
                         const double *outward = 0, double cellFrac = 1.0,
                         double mu = 0.0, double epsT = 0.0, bool epsTAuto = false,
                         double cohesion = 0.0, double tauMax = 0.0, bool consistentTan = false);
    int getNumMortarContacts(void) const { return (int)theMortarContacts.size(); }
    const MortarContact &getMortarContact(int i) const { return theMortarContacts[i]; }

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

    // --- P3: per-pair Coulomb friction STATE (the first path-dependent contact
    //     state). Lives on the engine (not the stateless adapter), keyed STABLY by
    //     (contactTag, slaveNodeTag, GLOBAL segment ordinal) so it survives the
    //     adapter rebuild every handle(). The adapter reads the committed plastic
    //     slip gpT + engagement origin gT0, runs the return map, and writes gpTtrial;
    //     commit()/revertToLastCommit() promote/drop the trial. ---
    struct FrictionState {
        double gpT[3];        // committed plastic slip from the engagement config
        double gpTtrial[3];   // trial plastic slip (written each getResidual)
        double gT0[3];        // ENGAGEMENT-config tangential origin (captured once)
        bool   engaged;       // has gT0 been captured at first activation
        FrictionState() : engaged(false) {
            for (int d = 0; d < 3; d++) gpT[d] = gpTtrial[d] = gT0[d] = 0.0;
        }
    };
    // lazily create + return the slot for a pair (zeroed, engaged=false if new).
    FrictionState &getOrCreateFrictionState(int contactTag, int slaveTag, int segIndex);

    // --- ADR-41 C2.2: per-GLOBAL-slave-node normal multiplier λ_N (frictionless
    //     commit-cycle ALM). Keyed (contactTag, slaveNodeTag) — NOT per facet: λ is one
    //     value per global slave node (per-facet-local λ is variationally inconsistent at
    //     shared nodes and fails the patch test). The augmented nodal pressure is
    //     p_I = min(0, λ_I + epsN·ḡ_I), ḡ_I = g̃_I^global / a_I^global where the GLOBAL
    //     weighted gap g̃_I^global = Σ_facets ∫N_I g_N dΓ and area a_I^global = Σ_facets ∫N_I dΓ
    //     are accumulated across the facets a shared slave node touches. λ_N is COMMITTED-ONLY
    //     (mutated solely in commit() ⇒ revertToLastCommit is automatically safe — the
    //     EmbeddedRebar invariant). One Uzawa step per Domain::commit(): the across-step
    //     augmentation; the within-step held-load `analyzeAugmented` proc drives ‖ḡ‖→augTol. ---
    struct MortarNormalState {
        double lambdaN;   // committed normal multiplier (≤ 0); updated ONLY in commit()
        double gtGlobal;  // Σ_facets g̃_I^facet at the current trial (incremental — see below)
        double aGlobal;   // Σ_facets a_I^facet = ∫N_I dΓ (incremental)
        double epsN;      // the penalty the adapters use at this node (for the commit update)
        // ADR-41 C3.1 — frictional mortar (folded into the SAME (contactTag,slaveNodeTag) slot
        // since it shares the key + lifecycle): per-global-slave-node elastic tangential slip.
        // gpT survives like lambdaN (committed plastic slip); gpTtrial is the transient trial
        // (a pure fn of committed state, rewritten each getResidual); gT0 is the engagement-config
        // tangential origin captured ONCE at first contact (else a late-engaging node's pre-contact
        // drift becomes a spurious stick traction — the ADR-39 P3 MAJOR-1, reused). λ_T (tangential
        // Uzawa) is C3.3 — C3.1/C3.2 ship penalty friction (λ_T≡0).
        double gpT[3];      // committed elastic tangential slip (promoted in commit())
        double gpTtrial[3]; // trial slip (written each getResidual)
        double gT0[3];      // engagement-config tangential origin (captured once)
        bool   engaged;     // has gT0 been captured at first contact activation
        MortarNormalState() : lambdaN(0.0), gtGlobal(0.0), aGlobal(0.0), epsN(0.0), engaged(false) {
            for (int d = 0; d < 3; d++) gpT[d] = gpTtrial[d] = gT0[d] = 0.0;
        }
    };
    // lazily create + return the per-node slot (zeroed if new).
    MortarNormalState &getOrCreateMortarNormalState(int contactTag, int slaveNodeTag);
    // a per-facet adapter reports its contribution g̃_I^facet / a_I^facet to a slave node,
    // keyed STABLY by the adapter's FE tag so a re-eval OVERWRITES (idempotent, never
    // double-counts). Maintains the node's gtGlobal/aGlobal as the running Σ over facets via a
    // delta update (write-then-read: call this BEFORE reading the node's global gap so the
    // calling facet's own current contribution is included). epsN is stored on the node slot.
    void accumulateMortarGap(int contactTag, int slaveNodeTag, int feTag,
                             double gtFacet, double aFacet, double epsN);
    // ‖ḡ‖_∞ restricted to PENETRATION (max over active nodes of max(0, −ḡ_I)) — the convergence
    // measure the held-load `analyzeAugmented` proc reads to stop augmenting. 0 if no contact.
    double getMaxMortarPenetration(void) const;
    int  getNumMortarNormalStates(void) const { return (int)theMortarNormalStates.size(); }

    // dead-slot GC (mirror frictionGC*): the engine survives domainChanged + across analyze()
    // calls, so a re-meshed analysis would leak old λ_N node slots. The handler rebuilds the
    // live node-set each handle(): mortarNormalGCBegin() → mortarNormalGCMark(...) per live slave
    // node → mortarNormalGCEnd() erases unmarked slots, drops ALL transient facet contributions,
    // and zeros gtGlobal/aGlobal on survivors (the adapters re-fill them next getResidual; only
    // λ_N is path-dependent and must survive). FE tags are reassigned each handle(), so the
    // facet-contribution map MUST be cleared here — its keys are not stable across rebuilds.
    void mortarNormalGCBegin(void);
    void mortarNormalGCMark(int contactTag, int slaveNodeTag);
    void mortarNormalGCEnd(void);

    // dead-slot GC (design-gate MAJOR: the engine survives domainChanged AND across
    // analyze() calls, so a re-meshed/re-paired analysis would leak old friction
    // slots — the ADR-30 theEQs leak class). The handler rebuilds the live key-set
    // each handle(): frictionGCBegin() → frictionGCMark(...) per live pair →
    // frictionGCEnd() erases any slot not marked.
    void frictionGCBegin(void);
    void frictionGCMark(int contactTag, int slaveTag, int segIndex);
    void frictionGCEnd(void);
    int  getNumFrictionStates(void) const { return (int)theFrictionStates.size(); }

    // --- lifecycle (driven by Domain::commit / revertToLastCommit) ---
    int commit(void);              // P3: gpT = gpTtrial for every slot (+ counter)
    int revertToLastCommit(void);  // P3: gpTtrial = gpT for every slot (+ counter)
    int getNumCommits(void) const { return numCommits; }
    int getNumReverts(void) const { return numReverts; }

  private:
    std::vector<LadrunoContactSurface *> theSurfaces;   // owned
    std::vector<Contact> theContacts;
    std::vector<MortarContact> theMortarContacts;       // ADR-41 C2 (separate from NTS)
    std::vector<RigidPlane> theRigidPlanes;             // P2a
    int numCommits;
    int numReverts;

    // P3 friction state store, keyed by (contactTag, slaveTag, segIndex).
    struct PairKey {
        int c, s, g;
        bool operator<(const PairKey &o) const {
            if (c != o.c) return c < o.c;
            if (s != o.s) return s < o.s;
            return g < o.g;
        }
    };
    std::map<PairKey, FrictionState> theFrictionStates;
    std::set<PairKey> liveKeys;                         // GC scratch (per handle())

    // C2.2 normal-ALM state, keyed by (contactTag, slaveNodeTag) — a 2-field key.
    struct NodeKey {
        int c, n;
        bool operator<(const NodeKey &o) const {
            if (c != o.c) return c < o.c;
            return n < o.n;
        }
    };
    std::map<NodeKey, MortarNormalState> theMortarNormalStates;
    std::set<NodeKey> liveNodeKeys;                     // mortar GC scratch (per handle())
    // transient per-facet contribution g̃_I^facet / a_I^facet keyed (contactTag, nodeTag, feTag);
    // delta-updated each getResidual, summed into the node slot's gtGlobal/aGlobal. Cleared every
    // handle() (feTag is not rebuild-stable). PairKey {c,s=nodeTag,g=feTag} reuses the 3-int key.
    struct FacetContrib { double gt, a; FacetContrib() : gt(0.0), a(0.0) {} };
    std::map<PairKey, FacetContrib> theMortarFacetContribs;
};

#endif
