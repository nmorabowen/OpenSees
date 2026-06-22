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

// ADR-39 P1a — LadrunoContactFE: an FE_Element adapter that injects contact
// contributions into the assembly without being backed by a Domain Element.
//
// The contact narrow phase is self-contained here: getResidual()/getTangent()
// read node trial displacements (like a real Element) and assemble F_c / K_c.
// Explicit integrators consume only getResidual() (formEleTangent is mass-only
// under CentralDifferenceLadruno); implicit consume both (tangent scaled by the
// integrator's c1). The adapter is a STATELESS VIEW — all path-dependent pair
// state lives on the Domain-owned LadrunoContactDomain (P1b), so the adapter can
// be destroyed/rebuilt every handle() with no state loss.
//
// P1a scope: EMPTY connectivity (size-0 getID/getDOFtags) + ZERO force, to prove
// the injection plumbing is graph-neutral (bitwise-identical to no-contact). With
// myEle==0 the base getResidual/getTangent exit(-1) and the base zero/add helpers
// early-return, so this subtype MUST own its own buffers and override both.
//
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p1_design.md.

#ifndef LadrunoContactFE_h
#define LadrunoContactFE_h

#include <FE_Element.h>
#include <Vector.h>
#include <Matrix.h>

class Integrator;

class LadrunoContactFE : public FE_Element
{
  public:
    // P1a: empty-connectivity adapter (numDOF_Group = 0, ndof = 0).
    LadrunoContactFE(int tag);
    ~LadrunoContactFE();

    // self-owned buffers (base buffers are unavailable when myEle == 0)
    const Vector &getResidual(Integrator *theIntegrator);
    const Matrix &getTangent(Integrator *theIntegrator);
    void addMtoTang(double fact = 1.0);   // contact pairs carry no mass -> no-op

    // With myEle==0 the base force-vector helpers print a WARNING and return the
    // shared error vector (FE_Element.cpp). They are reachable off the assembly
    // hot path (e.g. modal damping doMv -> getM_Force when mass is non-diagonal),
    // so override them to a clean size-0 return. (P2: real K/M/C force variants.)
    const Vector &getTangForce(const Vector &x, double fact = 1.0);
    const Vector &getK_Force(const Vector &x, double fact = 1.0);
    const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    const Vector &getC_Force(const Vector &x, double fact = 1.0);
    const Vector &getM_Force(const Vector &x, double fact = 1.0);

  private:
    Vector resid;   // size-0 in P1a (empty connectivity)
    Matrix tang;    // 0x0 in P1a
};

#endif
