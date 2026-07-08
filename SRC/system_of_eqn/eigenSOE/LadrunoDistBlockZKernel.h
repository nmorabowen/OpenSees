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

// Ladruno: DISTRIBUTED block-real complex-shift solve kernel (ADR 43, P3c-MPI,
// L3-only per ADR R0). See 43_ladruno_feast_eigensolver_adr.md and
// feast_d2_spike/README.md.
//
// The distributed twin of LadrunoBlockZKernel: it solves the SAME symmetric
// equivalent-real 2n x 2n block system for a complex shift z = a + i*b,
//
//     [ aM-K     -bM    ] [Xr]   [ Br]
//     [ -bM   -(aM-K)   ] [Xi] = [-Bi]        (2n x 2n, symmetric indefinite)
//
// but factors/solves it with DISTRIBUTED MUMPS (dmumps, SYM=2, LDL^T) across
// the ranks of MPI_COMM_WORLD instead of serial PARDISO. This is the L3 axis of
// the ADR's three-level plan: one contour node's linear solve spread over the
// ranks; the FEAST outer loop (dfeast_srci / runFeastRci) runs REPLICATED on
// every rank and calls this kernel collectively at each RCI ijob=10/11.
//
// Distribution model (L3-only, ADR R0): every rank holds the FULL (K,M) CSR
// (the FeastEigenSOE self-assembles identically on each replicated-interpreter
// rank) and builds the full 2n block triplet list, but each rank hands MUMPS
// only its DISJOINT SLICE of the triplets (ICNTL(18)=3 distributed input), so
// the FACTORIZATION is distributed even though the operator is replicated. The
// symbolic analysis (JOB=1) runs once on the distributed structure and is
// reused across every shift (ADR R5); only the numeric factor (JOB=2) re-runs
// per contour node. The dense RHS block is assembled on the host (rank 0),
// solved (JOB=3, nrhs columns at once), and the solution is MPI_Bcast to ALL
// ranks so the replicated RCI stays in lockstep.
//
// Compiled ONLY into the parallel executable targets (OpenSeesMP/PyMP), which
// carry _MUMPS + _PARALLEL_INTERPRETERS + mpi.h + dmumps on the link line;
// NOT into OPS_SysOfEqn. FeastEigenSolver reaches it through the
// LadrunoFeastInnerSolve factory hook, which this file self-registers at load.

#ifndef LadrunoDistBlockZKernel_h
#define LadrunoDistBlockZKernel_h

#include <LadrunoFeastInnerSolve.h>
#include <mpi.h>
#include <dmumps_c.h>

class LadrunoDistBlockZKernel : public LadrunoFeastInnerSolve
{
  public:
    // CSR views borrowed from FeastEigenSOE (n equations, FULL symmetric
    // pattern, 1-based indices, ascending columns per row) — replicated on
    // every rank. The kernel copies the block STRUCTURE for its local triplet
    // slice but reads Kvals/Mvals at setShift time, so the SOE must outlive it.
    // Uses MPI_COMM_WORLD (the single L3 solve group).
    LadrunoDistBlockZKernel(int n, const int *rowPtr, const int *colInd,
                            const double *Kvals, const double *Mvals);
    ~LadrunoDistBlockZKernel();

    int setShiftBlock(double a, double b) override;
    int solveBlock(int nrhs, const double *Br, const double *Bi,
                   double *Xr, double *Xi) override;

    bool ok(void) const {return healthy;}

  private:
    int n;                  // real problem size
    int nBlk;               // 2n
    const int *rowPtr;      // borrowed (SOE-owned)
    const int *colInd;
    const double *Kvals;
    const double *Mvals;

    MPI_Comm comm;
    int rank, np;

    // local triplet slice of the 2n block LOWER triangle (1-based for MUMPS):
    int nzLoc;
    int *irnLoc;            // nzLoc row indices (1-based)
    int *jcnLoc;            // nzLoc col indices (1-based)
    double *aLoc;           // nzLoc values (recomputed per shift)
    int *srcSlot;           // nzLoc -> slot k in the n-CSR value arrays
    signed char *srcCode;   // 0: aM-K   1: -bM   2: -(aM-K)

    DMUMPS_STRUC_C id;
    bool initialized;       // MUMPS instance live (JOB=-1 done)
    bool analyzed;          // JOB=1 done (symbolic; reused across shifts)
    bool healthy;           // construction + analysis succeeded

    double *rhsHost;        // host-only dense (2n x nrhs) RHS/solution buffer
    int rhsHostCols;        // its allocated column capacity

    void terminate(void);
};

#endif
