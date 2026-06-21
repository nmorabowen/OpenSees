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

#ifndef ExplicitBatheLNVDSMS_h
#define ExplicitBatheLNVDSMS_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// ExplicitBatheLNVDSMS — conventional (DT2MS-style) LUMPED selective MASS SCALING on
// the Noh-Bathe + FLAC local-non-viscous-damping integrator ExplicitBatheLNVD (classTag
// 33011). The LNVD sibling of ExplicitBatheSMS (33009) / CentralDifferenceSMS (33007);
// see Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.
//
// domainChanged() injects additive nodal mass (s_e-1)*m (s_e = betaK-damped closed form)
// via the shared Ladruno::buildMassScaling and restores it on teardown. ExplicitBatheLNVD
// assembles only the mass on the RHS (system Diagonal), so the injection is seen by the
// leap-frog with NO solve-path change. -cflAbort/-recompute rejected (ADR-36 MF-1).
//
// Command:  integrator ExplicitBatheLNVDSMS $p $alpha $dtTarget <-maxAddedMass f>
//                      <-lump rowsum|diagonal|hrz> <-tangent> <-verbose>

#include <ExplicitBatheLNVD.h>
#include <Vector.h>
#include <map>

class Domain;

class ExplicitBatheLNVDSMS : public ExplicitBatheLNVD
{
public:
    ExplicitBatheLNVDSMS();   // broker
    ExplicitBatheLNVDSMS(double p, double alpha_flac, double dtTarget,
                         double maxAddedMassFrac, bool verboseSMS,
                         int compute_critical_timestep, bool cflUseTangent,
                         CTSLumping lumping);
    ~ExplicitBatheLNVDSMS();

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

    std::map<int, Vector> injected;
    bool scaled;
    Domain *appliedDomain;
    bool warnedLimitations;

    void removeScaling(void);
};

#endif
