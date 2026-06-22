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

// ADR-39 P1a — LadrunoContactHandler: a ConstraintHandler that injects contact
// FE_Element adapters into the assembly. handle() is the ONLY legal FE_Element
// factory in OpenSees (called inside domainChanged after AnalysisModel::clearAll),
// so the contact subsystem MUST be a ConstraintHandler to get its adapters into
// the SOE. We REPLICATE the PlainHandler DOF_Group/FE loop (not compose — a
// composed inner handler returns count3 != numFe, so the contact-FE start tag is
// unknowable, and addFE_Element silently drops duplicate tags) and then append
// the contact adapter(s) with tags continuing above the structural FE count.
//
// P1a scope: inject ONE empty-connectivity zero LadrunoContactFE to prove the
// injection is graph-neutral (bitwise-identical to no-contact) under both an
// explicit and an implicit integrator. MP constraints (equalDOF) are NOT enforced
// in P1a (Plain-style) — a warning is emitted; delegation lands in P1b.
//
// Command: constraints LadrunoContact. HANDLER_TAG_LadrunoContactHandler 33002.
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p1_design.md.

#ifndef LadrunoContactHandler_h
#define LadrunoContactHandler_h

#include <ConstraintHandler.h>

class FE_Element;
class DOF_Group;

class LadrunoContactHandler : public ConstraintHandler
{
  public:
    LadrunoContactHandler();
    ~LadrunoContactHandler();

    int handle(const ID *nodesNumberedLast = 0);
    void clearAll(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
};

#endif
