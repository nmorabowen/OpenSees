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

// Ladruno: FEAST-RCI inner complex-shift solve interface (ADR 43, P3c). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md (S5.2, S9 R0/R1).
//
// The FEAST reverse-communication loop (dfeast_srci, runFeastRci in
// FeastEigenSolver.cpp) needs, at each contour node z = a + i*b, to factor
// (z*M - K) and solve (z*M - K) X = B for a block of complex right-hand sides.
// This pure-abstract seam lets the SAME RCI loop drive either inner solve:
//
//   - SERIAL (P3c-serial, shipped): LadrunoBlockZKernel -- the symmetric
//     equivalent-real 2n x 2n block system factored by MKL PARDISO (mtype -2),
//     entirely inside OPS_SysOfEqn.
//   - DISTRIBUTED (P3c-MPI, L3-only per ADR R0): a dmumps SYM=2 solve of the
//     same 2n block system across the ranks of a communicator. dmumps + mpi.h
//     live only in the parallel executable targets (OpenSeesMP/PyMP), NOT in
//     OPS_SysOfEqn -- so FeastEigenSolver cannot name that class directly. It
//     reaches it through the factory hook below, which the parallel target
//     registers at init and which stays null in the serial build.
//
// The interface is deliberately minimal: the RCI only ever needs to install a
// shift and solve a block. matvecZ / the -blockZGate self-test stay on the
// concrete serial kernel.

#ifndef LadrunoFeastInnerSolve_h
#define LadrunoFeastInnerSolve_h

class LadrunoFeastInnerSolve
{
  public:
    virtual ~LadrunoFeastInnerSolve() {}

    // (Re)assemble + numerically factor (z*M - K) for z = a + i*b, reusing the
    // symbolic analysis across shifts. Returns 0, or a solver error code.
    virtual int setShiftBlock(double a, double b) = 0;

    // Solve (z*M - K) X = B for nrhs columns (column-major, length-n each),
    // using the factorization the last setShiftBlock() produced. Bi == 0
    // (nullptr) means a real B. On return every participating rank holds the
    // full solution (the distributed impl broadcasts it) so a replicated RCI
    // stays in lockstep.
    virtual int solveBlock(int nrhs, const double *Br, const double *Bi,
                           double *Xr, double *Xi) = 0;

    // Force every rank to agree on an integer by taking rank 0's value
    // (collective; must be called in lockstep on all ranks). Serial default is
    // the identity. Used to pin the stochastic auto-seed m0 across ranks so the
    // FIRST distributed solve's block width (hence its MPI_Bcast length) is
    // identical everywhere — the catastrophic desync the replicated RCI must
    // avoid without relying on MKL thread-count control (which segfaults in the
    // MP MKL DLL layout). See FeastEigenSolver::solve() and the ADR R0 note.
    virtual int agreeInt(int v) { return v; }
};

// Factory hook for the DISTRIBUTED inner solve (ADR 43 P3c-MPI, L3-only).
// Defined = 0 in LadrunoFeastInnerSolve.cpp (linked into OPS_SysOfEqn, hence
// into every target). The parallel executable target (OpenSeesMP/PyMP), which
// compiles the dmumps-backed kernel with _MUMPS + _PARALLEL_INTERPRETERS +
// mpi.h, assigns this pointer at init. FeastEigenSolver::solve() consults it:
// when non-null AND the run is genuinely distributed (comm size > 1), it builds
// a distributed kernel over the shared CSR instead of the serial PARDISO one.
// Returns null if this rank cannot host a distributed solve (e.g. size == 1),
// signalling the caller to fall back to the serial kernel.
typedef LadrunoFeastInnerSolve *(*LadrunoFeastDistFactory)(
    int n, const int *rowPtr, const int *colInd,
    const double *Kvals, const double *Mvals);

extern LadrunoFeastDistFactory ladrunoFeastDistFactory;

#endif
