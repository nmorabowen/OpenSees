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

// Ladruno: shared refineAccel body for the THREE consistent (Olovsson) mass-scaling
// integrators (CentralDifferenceSMSConsistent 33008, ExplicitBatheSMSConsistent 33010,
// ExplicitBatheLNVDSMSConsistent 33012). Solves M_tilde a = r matrix-free, picking the
// SERIAL path (DiagonalSOE) or the DISTRIBUTED/MPI path (any LinearSOE whose
// isDistributedDiagonal() is true, i.e. MPIDiagonalSOE under OpenSeesMP `system
// MPIDiagonal`) at RUNTIME via the LinearSOE base virtuals -- the integrators live in
// the shared OpenSeesLIB and must NOT reference the MPI-only MPIDiagonalSOE type (it is
// linked only into the MP executables). See LadrunoMassScaling.h (the PCG kernels) and
// Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md (V5).
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoConsistentRefine_h
#define LadrunoConsistentRefine_h

#include <LadrunoMassScaling.h>   // ConsistentBlock + consistentPCG / consistentParPCG
#include <LinearSOE.h>
#include <DiagonalSOE.h>
#include <Vector.h>
#include <OPS_Globals.h>

#include <vector>

namespace Ladruno {

// Apply the consistent scaling mass to the incoming diagonal-solve accel `a`
// (= M_lump^-1 r), replacing it in place with M_tilde^-1 r. Returns 0 always
// (a non-fatal fallback just leaves `a` as the un-scaled diagonal solve).
// warnedSolver/warnedPCG are per-integrator latches so the warnings fire once.
inline int
applyConsistentRefine(std::vector<ConsistentBlock> &blocks, LinearSOE *theSOE, Vector &a,
                      double pcgTol, int pcgMaxIt, bool verbose,
                      bool &warnedSolver, bool &warnedPCG, const char *cls)
{
    if (theSOE == 0)
        return 0;

    // ---- DISTRIBUTED / MPI path (MPIDiagonalSOE) -------------------------------------
    // All ranks hold the same SOE type, so isDistributedDiagonal() is true on every rank
    // and the collective ops below stay in lockstep. A rank may have NO local scaled
    // elements yet still hold shared DOFs adjacent to a neighbour's M_bar, so the
    // participation decision must be GLOBAL (gActive) -- a local early-return would
    // deadlock the neighbours' assembleSharedSum / all-reduce.
    if (theSOE->isDistributedDiagonal()) {
        // Every decision to enter/skip the collective PCG MUST be global, or a rank that
        // returns early would deadlock the others' assembleSharedSum / all-reduce.
        const double *Ainv = theSOE->getScalingDiagonalA();   // factored 1/mass, GLOBAL at shared
        int neq = theSOE->getNumEqn();
        int    localActive = blocks.empty() ? 0 : 1;
        int    localBad    = (Ainv == 0 || neq != a.Size()) ? 1 : 0;
        double gActive     = theSOE->globalReduceSum((double)localActive);
        double gBad        = theSOE->globalReduceSum((double)localBad);
        if (gActive == 0.0)
            return 0;   // nothing scaled anywhere -> every rank returns together
        if (gBad > 0.0)
            return 0;   // some rank cannot apply the solve -> all bail together (no partial run)

        // weight w_i = 1/multiplicity_i (mult = #ranks sharing DOF i). mult is
        // assembleSharedSum(ones); 1 on purely-local DOFs. Rebuilt each call (one extra
        // exchange/step, cheap vs the PCG; robust to a re-partition / size change).
        std::vector<double> w((size_t)(neq > 0 ? neq : 0), 1.0);
        {
            Vector ones(neq);
            for (int i = 0; i < neq; ++i) ones(i) = 1.0;
            theSOE->assembleSharedSum(ones);
            for (int i = 0; i < neq; ++i)
                w[(size_t)i] = (ones(i) != 0.0) ? (1.0 / ones(i)) : 1.0;
        }

        Vector r(neq);
        for (int i = 0; i < neq; ++i)
            r(i) = (Ainv[i] != 0.0) ? (a(i) / Ainv[i]) : 0.0;

        LinearSOE *soe = theSOE;
        ConsistentParOps ops;
        ops.w = w.empty() ? 0 : w.data();
        ops.assembleShared = [soe](double *v, int n) { Vector vv(v, n); soe->assembleSharedSum(vv); };
        ops.globalSum      = [soe](double s) { return soe->globalReduceSum(s); };

        double relResid = 0.0;
        int iters = consistentParPCG(blocks, Ainv, neq, r, a, ops, pcgTol, pcgMaxIt, &relResid);
        if (verbose)
            opserr << "  " << cls << ": PCG " << iters << " iters (parallel), rel.resid "
                   << relResid << "\n";
        if (relResid > pcgTol && iters >= pcgMaxIt && !warnedPCG) {
            warnedPCG = true;
            opserr << "WARNING " << cls << "::refineAccel - the parallel mass-scaling PCG did NOT "
                      "converge (" << iters << " iters, rel.resid " << relResid << " > tol "
                   << pcgTol << "); consistent mass solve incomplete this step.\n";
        }
        return 0;
    }

    // ---- SERIAL / shared-memory path (DiagonalSOE) -----------------------------------
    if (blocks.empty())
        return 0;

    DiagonalSOE *theDiag = dynamic_cast<DiagonalSOE *>(theSOE);
    if (theDiag == 0) {
        if (!warnedSolver) {
            warnedSolver = true;
            opserr << "WARNING " << cls << "::refineAccel - consistent (Olovsson) mass scaling "
                      "REQUIRES `system Diagonal` (the lumped diagonal is the PCG preconditioner). "
                      "The current SOE is not a DiagonalSOE; falling back to the un-scaled diagonal "
                      "solve (the scaling mass is NOT applied -- expect the un-raised stable step).\n";
        }
        return 0;
    }

    const double *Ainv = theDiag->getDiagonalA();   // factored: 1/mass
    int neq = theDiag->getNumEqn();
    if (Ainv == 0 || neq != a.Size())
        return 0;

    Vector r(neq);
    for (int i = 0; i < neq; ++i)
        r(i) = (Ainv[i] != 0.0) ? (a(i) / Ainv[i]) : 0.0;

    double relResid = 0.0;
    int iters = consistentPCG(blocks, Ainv, neq, r, a, pcgTol, pcgMaxIt, &relResid);
    if (verbose)
        opserr << "  " << cls << ": PCG " << iters << " iters, rel.resid " << relResid << "\n";
    if (relResid > pcgTol && iters >= pcgMaxIt && !warnedPCG) {
        warnedPCG = true;
        opserr << "WARNING " << cls << "::refineAccel - the mass-scaling PCG did NOT converge ("
               << iters << " iters, rel.resid " << relResid << " > tol " << pcgTol
               << "); the consistent mass solve is incomplete this step (near-singular M_tilde?).\n";
    }
    return 0;
}

} // namespace Ladruno

#endif
