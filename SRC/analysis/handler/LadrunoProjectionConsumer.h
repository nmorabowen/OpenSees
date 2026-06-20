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

// ADR-30 P1: the abstract seam between LadrunoProjectionHandler and an explicit,
// projection-aware integrator. The handler builds the LadrunoConstraintProjector
// and pushes it through this interface (dynamic_cast on the Integrator&), so the
// handler never names a concrete integrator class and ExplicitBathe/LNVD can adopt
// the projection in P4 by implementing this one method. See
// Ladruno_implementation/30_ladruno_explicit_constraint_projection_adr.md (D3).

#ifndef LadrunoProjectionConsumer_h
#define LadrunoProjectionConsumer_h

class LadrunoConstraintProjector;

class LadrunoProjectionConsumer
{
  public:
    virtual ~LadrunoProjectionConsumer() {}

    // The handler pushes the (non-owning) projector at doneNumberingDOF(); the
    // consumer stores it and calls buildMass()/project() at the right points in
    // its leap-frog (project a0 in the starter before seeding v_{-1/2}; project
    // a_{n+1} before forming the full-step velocity).
    virtual void setConstraintProjector(LadrunoConstraintProjector *theProjector) = 0;
};

#endif
