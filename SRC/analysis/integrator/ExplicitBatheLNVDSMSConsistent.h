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

#ifndef ExplicitBatheLNVDSMSConsistent_h
#define ExplicitBatheLNVDSMSConsistent_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// ExplicitBatheLNVDSMSConsistent — CONSISTENT (Olovsson) selective mass scaling on the
// Noh-Bathe + FLAC local-non-viscous-damping integrator ExplicitBatheLNVD (classTag
// 33012). The LNVD sibling of ExplicitBatheSMSConsistent (33010) /
// CentralDifferenceSMSConsistent (33008); see
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.
//
// Injects the centroidal Olovsson scaling mass M_bar_e = beta[diag(m_a) - m m^T/M_e]
// (f1 preserved) and solves the non-diagonal M_tilde a = r matrix-free
// (Ladruno::consistentPCG) at BOTH Noh-Bathe sub-step solves via the base refineAccel()
// hook. Reuses Ladruno::buildMassScalingConsistent (slave+master MP exclusion). No Node
// mutation. Requires `system Diagonal`. -cflAbort/-recompute rejected.
//
// Command:  integrator ExplicitBatheLNVDSMSConsistent $p $alpha $dtTarget
//                      <-maxAddedMass f> <-lump ...> <-tangent> <-pcgTol t>
//                      <-pcgMaxIt n> <-verbose>

#include <ExplicitBatheLNVD.h>
#include <vector>

class Domain;
namespace Ladruno { struct ConsistentBlock; }

class ExplicitBatheLNVDSMSConsistent : public ExplicitBatheLNVD
{
public:
    ExplicitBatheLNVDSMSConsistent();   // broker
    ExplicitBatheLNVDSMSConsistent(double p, double alpha_flac, double dtTarget,
                                   double maxAddedMassFrac, bool verboseSMS,
                                   int compute_critical_timestep, bool cflUseTangent,
                                   CTSLumping lumping, double pcgTol, int pcgMaxIt);
    ~ExplicitBatheLNVDSMSConsistent();

    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

protected:
    bool refinesAccel(void) const override { return true; }
    int  refineAccel(Vector &a) override;

private:
    double dtTarget;
    double maxAddedMassFrac;
    bool   verboseSMS;
    CTSLumping lumpingSMS;
    bool   useTangentSMS;
    double pcgTol;
    int    pcgMaxIt;

    std::vector<Ladruno::ConsistentBlock> *blocks;
    bool warnedLimitations;
    bool warnedSolver;
    bool warnedPCG;
};

#endif
