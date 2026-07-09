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

// Ladruno: distributed block-real complex-shift solve kernel (ADR 43, P3c-MPI,
// L3-only) — see the header for the formulation, distribution model, and role.

#include <LadrunoDistBlockZKernel.h>
#include <OPS_Globals.h>

#include <cstring>
#include <cstdio>
#include <cstdlib>

#define ICNTL(I) icntl[(I)-1]   // MUMPS 1-based control indices (as elsewhere)

// ---------------------------------------------------------------------------
// Construction: build this rank's disjoint slice of the LOWER triangle of the
// 2n x 2n block system, initialize a MUMPS instance on MPI_COMM_WORLD, and run
// the symbolic analysis (JOB=1) once — reused across every later shift (R5).
//
// MUMPS symmetric (SYM=2) expects only the LOWER triangle (global row >= col),
// matching how OpenSees's own MumpsSOE assembles symmetric matrices. Block map
// (real block r = indices 0..n-1, imag block i = n..2n-1), all lower:
//   (i,   j),   j<=i : aM-K        (code 0)   block (r,r)
//   (n+i, j),   all j: -bM         (code 1)   block (i,r), fully lower
//   (n+i, n+j), j<=i : -(aM-K)     (code 2)   block (i,i)
// ---------------------------------------------------------------------------
LadrunoDistBlockZKernel::LadrunoDistBlockZKernel(int nEqn, const int *ia,
                                                 const int *ja,
                                                 const double *kv,
                                                 const double *mv)
    : n(nEqn), nBlk(2 * nEqn), rowPtr(ia), colInd(ja), Kvals(kv), Mvals(mv),
      comm(MPI_COMM_WORLD), rank(0), np(1),
      nzLoc(0), irnLoc(0), jcnLoc(0), aLoc(0), srcSlot(0), srcCode(0),
      initialized(false), analyzed(false), healthy(false),
      rhsHost(0), rhsHostCols(0)
{
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    // ---- pass 1: total lower-triangle block nnz (global) -------------------
    long total = 0;
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col <= row) total++;      // (r,r) lower
            total++;                      // (i,r) full
            if (col <= row) total++;      // (i,i) lower
        }
    }

    // ---- this rank's contiguous slice [lo, hi) of the global triplet list --
    const long lo = (total * rank) / np;
    const long hi = (total * (rank + 1)) / np;
    nzLoc = static_cast<int>(hi - lo);
    if (nzLoc < 0) nzLoc = 0;
    irnLoc  = new int[nzLoc > 0 ? nzLoc : 1];
    jcnLoc  = new int[nzLoc > 0 ? nzLoc : 1];
    aLoc    = new double[nzLoc > 0 ? nzLoc : 1];
    srcSlot = new int[nzLoc > 0 ? nzLoc : 1];
    srcCode = new signed char[nzLoc > 0 ? nzLoc : 1];

    // ---- pass 2: emit only the local slice, in the same global order -------
    long g = 0;
    int p = 0;
    auto emit = [&](int r1, int c1, int slot, signed char code) {
        if (g >= lo && g < hi) {
            irnLoc[p]  = r1;              // 1-based for MUMPS
            jcnLoc[p]  = c1;
            srcSlot[p] = slot;
            srcCode[p] = code;
            aLoc[p]    = 0.0;
            p++;
        }
        g++;
    };
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col <= row)
                emit(row + 1, col + 1, k, 0);              // aM-K
            emit(n + row + 1, col + 1, k, 1);              // -bM
            if (col <= row)
                emit(n + row + 1, n + col + 1, k, 2);      // -(aM-K)
        }
    }

    // ---- MUMPS instance on this communicator (mirror MumpsParallelSolver) --
    memset(&id, 0, sizeof(id));
    id.par = 1;                          // host participates in the factorization
    id.sym = 2;                          // general symmetric (LDL^T)
#ifdef _OPENMPI
    id.comm_fortran = (comm == MPI_COMM_WORLD) ? 0 : MPI_Comm_c2f(comm);
#else
    id.comm_fortran = MPI_Comm_c2f(comm);
#endif
    id.job = -1;                         // initialize
    dmumps_c(&id);
    if (id.infog[0] < 0) {
        opserr << "LadrunoDistBlockZKernel -- MUMPS init failed, infog(1)="
               << id.infog[0] << "\n";
        return;
    }
    initialized = true;

    id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0;
    id.ICNTL(5) = 0;                     // assembled format
    id.ICNTL(18) = 3;                    // distributed matrix structure/values
    id.ICNTL(14) = 35;                   // workspace headroom % (retry grows it on -8/-9)
    id.ICNTL(7) = 7;                     // ordering: automatic
    // NOTE (P3c-MPI gate MINOR-1): unlike the serial PARDISO kernel (1e-8 pivot
    // perturbation, iparm[9]=8) this leaves MUMPS null-pivot detection (ICNTL24)
    // and static pivoting (CNTL) at defaults. Harmless for normal FEAST use --
    // Gauss contour nodes carry Im(z)!=0 so (zM-K) stays nonsingular -- but a
    // mode landing on a contour node could fail JOB=2 here where the serial
    // path's perturbation survives. Failure is symmetric (broadcast INFOG) so
    // it refuses cleanly, never deadlocks. Matching the cushion (ICNTL24=1 +
    // CNTL(3)) is a deferred robustness follow-up.

    id.n = nBlk;                         // GLOBAL size, set on every rank
    id.nz_loc = nzLoc;
    id.irn_loc = irnLoc;
    id.jcn_loc = jcnLoc;
    id.a_loc = aLoc;                     // values filled per shift (JOB=2)

    id.job = 1;                          // symbolic analysis (structure only)
    dmumps_c(&id);
    if (id.infog[0] < 0) {
        opserr << "LadrunoDistBlockZKernel -- MUMPS analysis failed, infog(1)="
               << id.infog[0] << " infog(2)=" << id.infog[1] << "\n";
        return;
    }
    analyzed = true;
    healthy = true;

    // proof-of-distribution hook (opt-in): lets a test confirm the DISTRIBUTED
    // kernel actually ran — a silent fallback to the serial PARDISO kernel
    // would otherwise pass a same-run differential gate unnoticed. Emits this
    // rank's disjoint triplet-slice size; the slices partition the 2n block nnz.
    if (getenv("LADRUNO_FEAST_MPI_DEBUG"))
        fprintf(stderr, "LADRUNO_FEAST_MPI rank=%d np=%d nBlk=%d nzLoc=%d\n",
                rank, np, nBlk, nzLoc);
}

LadrunoDistBlockZKernel::~LadrunoDistBlockZKernel()
{
    terminate();
    delete[] irnLoc;
    delete[] jcnLoc;
    delete[] aLoc;
    delete[] srcSlot;
    delete[] srcCode;
    delete[] rhsHost;
}

void
LadrunoDistBlockZKernel::terminate(void)
{
    if (initialized) {
        id.job = -2;
        dmumps_c(&id);
        initialized = false;
    }
}

// (Re)compute this rank's local values for z = a + i*b and numerically factor.
int
LadrunoDistBlockZKernel::setShiftBlock(double a, double b)
{
    if (!healthy) return -1;

    for (int q = 0; q < nzLoc; q++) {
        const int k = srcSlot[q];
        const double amk = a * Mvals[k] - Kvals[k];
        aLoc[q] = (srcCode[q] == 0) ? amk
                : (srcCode[q] == 1) ? -b * Mvals[k]
                                    : -amk;
    }

    id.a_loc = aLoc;
    id.job = 2;                          // numeric factorization (reuses JOB=1)
    dmumps_c(&id);

    // MUMPS INFOG(1) = -9 (real) / -8 (integer) => the internal work array is too
    // small; the documented remedy is to raise ICNTL(14) (workspace headroom %)
    // and re-run the numeric factorization. Indefinite SYM=2 3D block-real fill is
    // hard to predict (grows with rank count as the frontal tree is distributed),
    // so the fixed 20% headroom under-provisions at higher np (observed:
    // infog(1)=-9 at np=8). Retry with geometrically-growing headroom. This is
    // COLLECTIVE — every rank sees the same global INFOG(1) and retries in
    // lockstep, so it can never deadlock. ICNTL(14) is restored to the default
    // afterwards so a one-off hard shift doesn't inflate every later factor.
    const int icntl14_default = 35;
    int tries = 0;
    while ((id.infog[0] == -9 || id.infog[0] == -8) && tries < 5) {
        const int prev = id.ICNTL(14);
        id.ICNTL(14) = (prev < 60) ? 100 : prev * 2;   // 35 -> 100 -> 200 -> ...
        if (id.ICNTL(14) > 4000) break;
        opserr << "LadrunoDistBlockZKernel::setShiftBlock -- MUMPS workspace short "
               << "(infog(1)=" << id.infog[0] << ", need infog(2)=" << id.infog[1]
               << "); retrying with ICNTL(14)=" << id.ICNTL(14) << "\n";
        id.job = 2;
        dmumps_c(&id);
        tries++;
    }
    id.ICNTL(14) = icntl14_default;      // reset headroom for subsequent shifts

    if (id.infog[0] < 0) {
        opserr << "LadrunoDistBlockZKernel::setShiftBlock -- MUMPS factor "
               << "failed, infog(1)=" << id.infog[0] << " infog(2)="
               << id.infog[1] << " (shift " << a << "+" << b << "i)\n";
        return id.infog[0];
    }
    return 0;
}

// Solve (zM - K) X = B for nrhs columns. RHS assembled on the host as the
// stacked [Br; -Bi] (2n x nrhs), solved centralized, then broadcast so every
// rank returns the full solution (replicated-RCI lockstep).
int
LadrunoDistBlockZKernel::solveBlock(int nrhs, const double *Br,
                                    const double *Bi, double *Xr, double *Xi)
{
    if (!healthy) return -1;

    // size_t for the allocation (2n*nrhs can exceed INT_MAX at extreme scale);
    // the MPI_Bcast count below is int-limited, which is the real scale ceiling.
    const size_t need = static_cast<size_t>(nBlk) * nrhs;
    if (nrhs > rhsHostCols) {
        delete[] rhsHost;
        rhsHost = new double[need > 0 ? need : 1];
        rhsHostCols = nrhs;
    }

    // host assembles the stacked RHS [Br; -Bi] (column-major, ld = 2n)
    if (rank == 0) {
        for (int c = 0; c < nrhs; c++) {
            double *col = &rhsHost[static_cast<size_t>(c) * nBlk];
            for (int i = 0; i < n; i++) {
                col[i]     = Br[static_cast<size_t>(c) * n + i];
                col[n + i] = (Bi != 0) ? -Bi[static_cast<size_t>(c) * n + i]
                                       : 0.0;
            }
        }
        id.rhs = rhsHost;
        id.nrhs = nrhs;
        id.lrhs = nBlk;
    }

    id.job = 3;                          // solve (reuses the JOB=2 factor)
    dmumps_c(&id);
    if (id.infog[0] < 0) {
        opserr << "LadrunoDistBlockZKernel::solveBlock -- MUMPS solve failed, "
               << "infog(1)=" << id.infog[0] << "\n";
        return id.infog[0];
    }

    // solution is centralized on the host — broadcast so all ranks agree
    MPI_Bcast(rhsHost, static_cast<int>(need), MPI_DOUBLE, 0, comm);

    for (int c = 0; c < nrhs; c++) {
        const double *col = &rhsHost[static_cast<size_t>(c) * nBlk];
        for (int i = 0; i < n; i++) {
            Xr[static_cast<size_t>(c) * n + i] = col[i];
            Xi[static_cast<size_t>(c) * n + i] = col[n + i];
        }
    }
    return 0;
}

int
LadrunoDistBlockZKernel::agreeInt(int v)
{
    // rank 0's value wins on every rank — pins the stochastic auto-seed m0 so
    // the first solve's block width is identical everywhere (lockstep).
    MPI_Bcast(&v, 1, MPI_INT, 0, comm);
    return v;
}

// ---------------------------------------------------------------------------
// Factory + self-registration (ADR 43 P3c-MPI). Linked ONLY into the parallel
// executable target, so this hook is set exactly in the builds that can host a
// distributed solve. Returns null on a single rank so the caller falls back to
// the serial PARDISO kernel — keeping np==1 byte-identical to serial -rci.
// ---------------------------------------------------------------------------
namespace {

LadrunoFeastInnerSolve *
makeDistBlockZKernel(int n, const int *rowPtr, const int *colInd,
                     const double *Kvals, const double *Mvals)
{
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np <= 1)
        return 0;
    LadrunoDistBlockZKernel *k =
        new LadrunoDistBlockZKernel(n, rowPtr, colInd, Kvals, Mvals);
    if (!k->ok()) {
        delete k;
        return 0;
    }
    return k;
}

// Register only under the replicated-interpreter model (OpenSeesMP/PyMP), where
// every rank self-assembles the FULL CSR — the assumption the L3 distribution
// rests on. NOT under _PARALLEL_PROCESSING (OpenSeesSP), where the domain is
// partitioned and no single rank holds the whole operator. (This TU is added
// only to the _PARALLEL_INTERPRETERS targets in CMake; the guard documents and
// enforces that intent.)
#ifdef _PARALLEL_INTERPRETERS
struct DistFactoryRegistrar {
    DistFactoryRegistrar() { ladrunoFeastDistFactory = &makeDistBlockZKernel; }
};
DistFactoryRegistrar s_distFactoryRegistrar;
#endif

} // namespace
