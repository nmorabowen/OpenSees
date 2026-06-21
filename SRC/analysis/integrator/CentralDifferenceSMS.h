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

#ifndef CentralDifferenceSMS_h
#define CentralDifferenceSMS_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// CentralDifferenceSMS — selective (conventional / LS-DYNA DT2MS-style) MASS
// SCALING on top of the robust explicit leap-frog CentralDifferenceLadruno
// (classTag 33007; ladruno private band >=33000). See
// Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.
//
// At domainChanged() it inflates the mass of the elements that throttle the
// timestep (per-element stable step dt_e < dtTarget) so the WHOLE model can run
// at a global step set by the bulk: for each such element, scale its mass by
// s_e=(dtTarget/dt_e)^2 (so dt_e -> dtTarget) by ADDING the increment (s_e-1)*m to
// its nodes (additive nodal mass; energy still closes because the
// EnergyBalanceRecorder KE reads nodal mass too). Re-baselined on every
// domainChanged (the injected increment is recorded and removed before recompute)
// so node masses are never permanently or doubly mutated.
//
// dt_e comes from the shared CriticalTimeStep kernel (self-report aware), so SMS
// lumps/eigensolves EXACTLY as the dt_cr query. -cflAbort and -recompute are
// REJECTED with SMS: their inherited path re-runs the element-mass eigensolve
// (which cannot see the nodal augmentation) and would spuriously hard-abort a run
// that is in fact stable (ADR-36 MF-1).
//
// Parallel: the injected lumped mass on a shared/boundary node IS reduced across
// ranks by a distributed/MPI diagonal solver (OpenSeesMP `system MPIDiagonal`,
// OpenSeesSP `system Diagonal` -> DistributedDiagonalSOE), since the explicit
// M^-1 solve reads that summed diagonal -- validated bit-identical to serial
// (ADR-36 T-MPI). The CONSISTENT (Olovsson) sibling is NOT parallel-safe: its
// matrix-free PCG (LadrunoMassScaling.h) uses rank-local inner products.

#include <CentralDifferenceLadruno.h>
#include <Vector.h>
#include <map>

class Domain;

class CentralDifferenceSMS : public CentralDifferenceLadruno
{
public:
    CentralDifferenceSMS();   // broker
    CentralDifferenceSMS(double dtTarget, double maxAddedMassFrac, bool verboseSMS,
                         int compute_critical_timestep, bool cflUseTangent,
                         CTSLumping lumping);
    ~CentralDifferenceSMS();

    // base allocates state + (optional) pre-scaling dt_cr; then the scaling pass
    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

private:
    double dtTarget;          // target global stable step (DT2MS-style)
    double maxAddedMassFrac;  // warn when added/model mass exceeds this (0 = no cap)
    bool   verboseSMS;        // report scaling summary each domainChanged
    CTSLumping lumpingSMS;    // own copy of the lump (base member is private)
    bool   useTangentSMS;     // own copy of useTangent

    std::map<int, Vector> injected;  // per-node-tag injected diagonal ΔM (re-baseline)
    bool scaled;                     // is `injected` currently committed to the Domain?
    Domain *appliedDomain;           // Domain the injection was committed to; used to
                                     //   restore in the destructor (SMS-NO-RESTORE) so
                                     //   a later stage / integrator on the same model
                                     //   sees the ORIGINAL masses.
    bool warnedLimitations;          // one-time guard for the v1-limitation warnings

    void removeScaling(void);        // subtract injected ΔM from node masses, clear map
};

#endif
