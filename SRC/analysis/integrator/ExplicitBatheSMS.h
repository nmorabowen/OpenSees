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

#ifndef ExplicitBatheSMS_h
#define ExplicitBatheSMS_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// ExplicitBatheSMS — conventional (LS-DYNA DT2MS-style) LUMPED selective MASS
// SCALING on top of the Noh-Bathe ExplicitBathe integrator (classTag 33009; ladruno
// private band >=33000). The ExplicitBathe sibling of CentralDifferenceSMS (33007);
// see Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md (ADR-36 noted
// "ExplicitBatheSMS is a trivial follow-up once the util is proven").
//
// At domainChanged() it inflates the mass of the elements that throttle the timestep
// (per-element stable step dt_e < dtTarget) by ADDING the increment (s_e-1)*m to their
// nodes (additive nodal mass; s_e = betaK-damped closed form). Because ExplicitBathe
// assembles ONLY the mass on the RHS (system Diagonal), the injected nodal mass is
// picked up by the leap-frog automatically — no solve-path change is needed (that is
// what makes the lumped variant trivial; the CONSISTENT variant is ExplicitBatheSMS-
// Consistent, which needs the matrix-free PCG via the refineAccel hook). Re-baselined
// on every domainChanged and restored on teardown (no permanent Domain mutation).
//
// dt_e comes from the shared CriticalTimeStep kernel (self-report aware), via the
// shared Ladruno::buildMassScaling util. -cflAbort/-recompute are DOWNGRADED to
// report-only with SMS (ADR-36 MF-1, ADR-52 W1-E3a): their inherited path re-runs the
// un-augmented element eigensolve, so rather than abort we keep the run and report the
// pre-scaling dt_cr instead.
//
// Command:  integrator ExplicitBatheSMS $p $dtTarget <-maxAddedMass f>
//                      <-lump rowsum|diagonal|hrz> <-tangent> <-verbose>

#include <ExplicitBathe.h>
#include <Vector.h>
#include <map>

class Domain;

class ExplicitBatheSMS : public ExplicitBathe
{
public:
    ExplicitBatheSMS();   // broker
    ExplicitBatheSMS(double p, double dtTarget, double maxAddedMassFrac, bool verboseSMS,
                     int compute_critical_timestep, bool cflUseTangent, CTSLumping lumping);
    ~ExplicitBatheSMS();

    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

private:
    double dtTarget;
    double maxAddedMassFrac;
    bool   verboseSMS;
    CTSLumping lumpingSMS;
    bool   useTangentSMS;

    std::map<int, Vector> injected;  // per-node-tag injected diagonal ΔM (re-baseline)
    bool scaled;
    Domain *appliedDomain;
    bool warnedLimitations;

    void removeScaling(void);
};

#endif
