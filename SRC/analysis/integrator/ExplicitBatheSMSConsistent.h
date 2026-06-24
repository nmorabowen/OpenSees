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

#ifndef ExplicitBatheSMSConsistent_h
#define ExplicitBatheSMSConsistent_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// ExplicitBatheSMSConsistent — CONSISTENT (Olovsson) selective mass scaling on the
// Noh-Bathe ExplicitBathe integrator (classTag 33010; ladruno private band >=33000).
// The ExplicitBathe sibling of CentralDifferenceSMSConsistent (33008); see
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.
//
// Injects the centroidal Olovsson scaling mass M_bar_e = beta[diag(m_a) - m m^T/M_e]
// (row sums zero -> rigid translation gets no added inertia -> f1 preserved). The
// resulting M_tilde = M_lump + sum_e M_bar_e is SPD but NON-diagonal, so each of
// ExplicitBathe's TWO Noh-Bathe sub-step solves (a = M_lump^-1 r) is corrected to
// a = M_tilde^-1 r by a matrix-free Jacobi(lumped-diagonal)-preconditioned CG via the
// base refineAccel() hook (the shipped Ladruno::consistentPCG kernel). No Node mass is
// mutated. Reuses the shared sizing (dt_e / self-report / betaK / MP slave+master
// exclusion / cap) of Ladruno::buildMassScalingConsistent.
//
// Required recipe: `system Diagonal` (preconditioner source), `algorithm Linear`,
// dt <= dtTarget. -cflAbort / -recompute are DOWNGRADED to report-only with SMS (ADR-36
// MF-1, ADR-52 W1-E3a): their inherited path re-runs the un-augmented element eigensolve,
// so rather than abort we keep the run and report the pre-scaling dt_cr instead.
//
// Command:  integrator ExplicitBatheSMSConsistent $p $dtTarget <-maxAddedMass f>
//                      <-lump rowsum|diagonal|hrz> <-tangent> <-pcgTol t> <-pcgMaxIt n>
//                      <-verbose>

#include <ExplicitBathe.h>
#include <vector>

class Domain;
namespace Ladruno { struct ConsistentBlock; }

class ExplicitBatheSMSConsistent : public ExplicitBathe
{
public:
    ExplicitBatheSMSConsistent();   // broker
    ExplicitBatheSMSConsistent(double p, double dtTarget, double maxAddedMassFrac,
                               bool verboseSMS, int compute_critical_timestep,
                               bool cflUseTangent, CTSLumping lumping,
                               double pcgTol, int pcgMaxIt);
    ~ExplicitBatheSMSConsistent();

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
