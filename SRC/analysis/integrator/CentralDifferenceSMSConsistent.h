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

#ifndef CentralDifferenceSMSConsistent_h
#define CentralDifferenceSMSConsistent_h

// Written: nmb (UANDES), 06/2026  (Ladruno sibling-fork explicit integrator)
//
// CentralDifferenceSMSConsistent — CONSISTENT (Olovsson) selective mass scaling on
// top of the robust explicit leap-frog CentralDifferenceLadruno (classTag 33008;
// ladruno private integrator band >=33000). The v2 of the conventional/lumped
// CentralDifferenceSMS (33007). See
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md.
//
// Conventional (lumped) mass scaling adds an isotropic diagonal increment (s_e-1)*m to
// each throttling element's nodes — which shifts ALL global frequencies, including the
// fundamental f1 that seismic SSI cannot afford to move. This class instead injects the
// CENTROIDAL Olovsson element scaling mass
//
//     M_bar_e = beta_e * [ diag(m_a) - m m^T / M_e ]          (per spatial direction)
//
// whose row sums are zero, so M_bar_e * t_rigid = 0: rigid-body translation gets NO
// added inertia and f1 is preserved, while the relative/deformation modes that govern
// the explicit step are loaded -> dt_e -> dtTarget (the SAME betaK-damped factor
// beta_e = T^2 + 2*T*c - 1 the lumped path uses). Oracle: f1 drift < 0.3% vs 54-83%
// for lumped, at half the added mass.
//
// M_tilde = M_lump + sum_e M_bar_e is SPD but NON-diagonal, so the leap-frog can no
// longer use the trivial `system Diagonal` a = M^-1 r. The diagonal solve still runs
// (it supplies the lumped diagonal = the Jacobi preconditioner), and this class
// overrides the base refineAccel() hook to replace its result a = M_lump^-1 r by the
// matrix-free PCG solve a = M_tilde^-1 r (Ladruno::consistentPCG). No Node mass is ever
// mutated — the cross-node coupling cannot live in per-node mass, so the inject/restore
// lifecycle of the lumped sibling is not needed.
//
// Required recipe: `system Diagonal` (preconditioner source), `algorithm Linear`,
// dt <= dtTarget. -cflAbort / -recompute are DOWNGRADED to report-only with SMS (ADR-36
// MF-1, ADR-52 W1-E3a): their inherited path re-runs the un-augmented element eigensolve,
// so rather than abort we keep the run and report the pre-scaling dt_cr instead.

#include <CentralDifferenceLadruno.h>
#include <vector>

class Domain;
namespace Ladruno { struct ConsistentBlock; }

class CentralDifferenceSMSConsistent : public CentralDifferenceLadruno
{
public:
    CentralDifferenceSMSConsistent();   // broker
    CentralDifferenceSMSConsistent(double dtTarget, double maxAddedMassFrac,
                                   bool verboseSMS, int compute_critical_timestep,
                                   bool cflUseTangent, CTSLumping lumping,
                                   double pcgTol, int pcgMaxIt);
    ~CentralDifferenceSMSConsistent();

    // base allocates state + (optional) pre-scaling dt_cr; then build the consistent
    // scaling blocks (no Node mutation).
    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

protected:
    // Replace the diagonal solve a = M_lump^-1 r by a = M_tilde^-1 r (matrix-free PCG).
    bool refinesAccel(void) const override { return true; }
    int  refineAccel(Vector &a) override;

private:
    double dtTarget;          // target global stable step (DT2MS-style)
    double maxAddedMassFrac;  // warn when added/model mass exceeds this (0 = no cap)
    bool   verboseSMS;        // report scaling + PCG summary each step/domainChanged
    CTSLumping lumpingSMS;    // own copy of the lump (base member is private)
    bool   useTangentSMS;     // own copy of useTangent
    double pcgTol;            // PCG relative-residual tolerance
    int    pcgMaxIt;          // PCG iteration cap

    std::vector<Ladruno::ConsistentBlock> *blocks;  // per-element scaling blocks (owned;
                                                    //   rebuilt each domainChanged)
    bool warnedLimitations;   // one-time guard for the limitation disclosures
    bool warnedSolver;        // one-time guard for the non-Diagonal-SOE warning
    bool warnedPCG;           // one-time guard for the PCG non-convergence warning
};

#endif
