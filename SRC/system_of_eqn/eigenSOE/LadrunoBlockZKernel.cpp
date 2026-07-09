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

// Ladruno: block-real complex-shift solve kernel (ADR 43, P3b) — see the
// header for the formulation and role.

#include <LadrunoBlockZKernel.h>
#include <OPS_Globals.h>

#include <cstring>
#include <chrono>

#if defined(_WIN32) || defined(_LADRUNO_MKL_FEAST)
extern "C" {
  // mkl_pardiso.h contract, LP64 (MKL_INT == int on this fork's build)
  void pardiso(void *pt, const int *maxfct, const int *mnum,
               const int *mtype, const int *phase, const int *n,
               const void *a, const int *ia, const int *ja,
               int *perm, const int *nrhs, int *iparm,
               const int *msglvl, void *b, void *x, int *error);
}

namespace {

double
nowSeconds(void)
{
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch())
        .count();
}

// shared PARDISO iparm setup for a solve (NOT the inertia certificate):
// matching/scaling OFF so the analysis phase is value-independent and can be
// reused verbatim across shifts (ADR R5); Bunch-Kaufman pivoting with the
// default 1e-8 perturbation handles the indefinite pivots.
void
initSolveIparm(int *iparm)
{
    std::memset(iparm, 0, 64 * sizeof(int));
    iparm[0]  = 1;   // user-supplied iparm
    iparm[1]  = 2;   // METIS fill-reduction
    iparm[9]  = 8;   // pivot perturbation 1e-8 (indefinite default)
    iparm[10] = 0;   // scaling OFF  — keeps phase-11 analysis shift-invariant
    iparm[12] = 0;   // matching OFF — ditto
}

} // namespace
#endif

LadrunoBlockZKernel::LadrunoBlockZKernel(int nEqn, const int *ia,
                                         const int *ja, const double *kv,
                                         const double *mv)
    : n(nEqn), rowPtr(ia), colInd(ja), Kvals(kv), Mvals(mv),
      nBlk(2 * nEqn), nnzBlk(0), iaBlk(0), jaBlk(0), aBlk(0),
      srcSlot(0), srcCode(0), analyzedBlk(false),
      nnzUp(0), iaZ(0), jaZ(0), aZ(0), srcSlotZ(0), analyzedZ(false),
      tFactor(0.0), tSolve(0.0)
{
    std::memset(ptBlk, 0, sizeof(ptBlk));
    std::memset(ptZ, 0, sizeof(ptZ));
    std::memset(iparmBlk, 0, sizeof(iparmBlk));
    std::memset(iparmZ, 0, sizeof(iparmZ));

    // ---- upper-triangle count of the (K,M) pattern -----------------------
    for (int row = 0; row < n; row++)
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++)
            if (colInd[k] - 1 >= row)
                nnzUp++;

    // ---- block-real pattern: 2n, upper triangle, 1-based -----------------
    // row i   (< n): [aM-K upper of row i] then [-bM ALL of row i at n+j]
    // row n+i     : [-(aM-K) upper of row i at n+j]
    const int nnzFull = rowPtr[n] - 1;
    nnzBlk = 2 * nnzUp + nnzFull;
    iaBlk = new int[nBlk + 1];
    jaBlk = new int[nnzBlk];
    aBlk = new double[nnzBlk];
    srcSlot = new int[nnzBlk];
    srcCode = new signed char[nnzBlk];

    int pos = 0;
    iaBlk[0] = 1;
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col >= row) {
                jaBlk[pos] = col + 1;          // (i, j)   <- aM-K
                srcSlot[pos] = k;
                srcCode[pos] = 0;
                pos++;
            }
        }
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            jaBlk[pos] = n + col + 1;          // (i, n+j) <- -bM
            srcSlot[pos] = k;
            srcCode[pos] = 1;
            pos++;
        }
        iaBlk[row + 1] = pos + 1;
    }
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col >= row) {
                jaBlk[pos] = n + col + 1;      // (n+i, n+j) <- -(aM-K)
                srcSlot[pos] = k;
                srcCode[pos] = 2;
                pos++;
            }
        }
        iaBlk[n + row + 1] = pos + 1;
    }

    // ---- complex-symmetric pattern: n, upper triangle, 1-based -----------
    iaZ = new int[n + 1];
    jaZ = new int[nnzUp];
    aZ = new double[2 * static_cast<size_t>(nnzUp)];
    srcSlotZ = new int[nnzUp];
    pos = 0;
    iaZ[0] = 1;
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col >= row) {
                jaZ[pos] = col + 1;
                srcSlotZ[pos] = k;
                pos++;
            }
        }
        iaZ[row + 1] = pos + 1;
    }
}

LadrunoBlockZKernel::~LadrunoBlockZKernel()
{
    releaseBlock();
    releaseComplex();
    delete[] iaBlk;
    delete[] jaBlk;
    delete[] aBlk;
    delete[] srcSlot;
    delete[] srcCode;
    delete[] iaZ;
    delete[] jaZ;
    delete[] aZ;
    delete[] srcSlotZ;
}

void
LadrunoBlockZKernel::releaseBlock(void)
{
#if defined(_WIN32) || defined(_LADRUNO_MKL_FEAST)
    if (!analyzedBlk)
        return;
    int phase = -1, error = 0, idum = 0;
    const int maxfct = 1, mnum = 1, mtype = -2, nrhs = 0, msglvl = 0;
    double ddum = 0.0;
    pardiso(ptBlk, &maxfct, &mnum, &mtype, &phase, &nBlk,
            aBlk, iaBlk, jaBlk, &idum, &nrhs, iparmBlk, &msglvl,
            &ddum, &ddum, &error);
    analyzedBlk = false;
#endif
}

void
LadrunoBlockZKernel::releaseComplex(void)
{
#if defined(_WIN32) || defined(_LADRUNO_MKL_FEAST)
    if (!analyzedZ)
        return;
    int phase = -1, error = 0, idum = 0;
    const int maxfct = 1, mnum = 1, mtype = 6, nrhs = 0, msglvl = 0;
    double ddum = 0.0;
    pardiso(ptZ, &maxfct, &mnum, &mtype, &phase, &n,
            aZ, iaZ, jaZ, &idum, &nrhs, iparmZ, &msglvl,
            &ddum, &ddum, &error);
    analyzedZ = false;
#endif
}

int
LadrunoBlockZKernel::setShiftBlock(double a, double b)
{
#if !defined(_WIN32) && !defined(_LADRUNO_MKL_FEAST)
    (void)a; (void)b;
    opserr << "LadrunoBlockZKernel -- requires MKL (Windows/oneAPI build)\n";
    return -1;
#else
    for (int p = 0; p < nnzBlk; p++) {
        const int k = srcSlot[p];
        const double amk = a * Mvals[k] - Kvals[k];
        aBlk[p] = (srcCode[p] == 0) ? amk
                : (srcCode[p] == 1) ? -b * Mvals[k]
                                    : -amk;
    }

    const int maxfct = 1, mnum = 1, mtype = -2, nrhs = 0, msglvl = 0;
    int error = 0, idum = 0;
    double ddum = 0.0;

    if (!analyzedBlk) {
        initSolveIparm(iparmBlk);
        int phase = 11;    // analysis — reused for every later shift (R5)
        pardiso(ptBlk, &maxfct, &mnum, &mtype, &phase, &nBlk,
                aBlk, iaBlk, jaBlk, &idum, &nrhs, iparmBlk, &msglvl,
                &ddum, &ddum, &error);
        if (error != 0) {
            opserr << "LadrunoBlockZKernel::setShiftBlock -- PARDISO "
                   << "analysis error " << error << "\n";
            return error;
        }
        analyzedBlk = true;
    }

    const double t0 = nowSeconds();
    int phase = 22;        // numeric factorization for this shift
    pardiso(ptBlk, &maxfct, &mnum, &mtype, &phase, &nBlk,
            aBlk, iaBlk, jaBlk, &idum, &nrhs, iparmBlk, &msglvl,
            &ddum, &ddum, &error);
    tFactor = nowSeconds() - t0;
    if (error != 0)
        opserr << "LadrunoBlockZKernel::setShiftBlock -- PARDISO "
               << "factorization error " << error << "\n";
    return error;
#endif
}

int
LadrunoBlockZKernel::setShiftComplex(double a, double b)
{
#if !defined(_WIN32) && !defined(_LADRUNO_MKL_FEAST)
    (void)a; (void)b;
    opserr << "LadrunoBlockZKernel -- requires MKL (Windows/oneAPI build)\n";
    return -1;
#else
    // zM - K = (aM - K) + i (bM), complex symmetric, upper triangle
    for (int p = 0; p < nnzUp; p++) {
        const int k = srcSlotZ[p];
        aZ[2 * p]     = a * Mvals[k] - Kvals[k];
        aZ[2 * p + 1] = b * Mvals[k];
    }

    const int maxfct = 1, mnum = 1, mtype = 6, nrhs = 0, msglvl = 0;
    int error = 0, idum = 0;
    double ddum = 0.0;

    if (!analyzedZ) {
        initSolveIparm(iparmZ);
        int phase = 11;
        pardiso(ptZ, &maxfct, &mnum, &mtype, &phase, &n,
                aZ, iaZ, jaZ, &idum, &nrhs, iparmZ, &msglvl,
                &ddum, &ddum, &error);
        if (error != 0) {
            opserr << "LadrunoBlockZKernel::setShiftComplex -- PARDISO "
                   << "analysis error " << error << "\n";
            return error;
        }
        analyzedZ = true;
    }

    const double t0 = nowSeconds();
    int phase = 22;
    pardiso(ptZ, &maxfct, &mnum, &mtype, &phase, &n,
            aZ, iaZ, jaZ, &idum, &nrhs, iparmZ, &msglvl,
            &ddum, &ddum, &error);
    tFactor = nowSeconds() - t0;
    if (error != 0)
        opserr << "LadrunoBlockZKernel::setShiftComplex -- PARDISO "
               << "factorization error " << error << "\n";
    return error;
#endif
}

int
LadrunoBlockZKernel::solveBlock(int nrhs, const double *Br, const double *Bi,
                                double *Xr, double *Xi)
{
#if !defined(_WIN32) && !defined(_LADRUNO_MKL_FEAST)
    (void)nrhs; (void)Br; (void)Bi; (void)Xr; (void)Xi;
    return -1;
#else
    if (!analyzedBlk) {
        opserr << "LadrunoBlockZKernel::solveBlock -- no factorization\n";
        return -1;
    }

    // stacked RHS [Br; -Bi] per column (see header: the symmetric form
    // negates the imaginary-part equation)
    double *rhs = new double[static_cast<size_t>(nBlk) * nrhs];
    double *sol = new double[static_cast<size_t>(nBlk) * nrhs];
    for (int c = 0; c < nrhs; c++) {
        for (int i = 0; i < n; i++) {
            rhs[static_cast<size_t>(c) * nBlk + i] =
                Br[static_cast<size_t>(c) * n + i];
            rhs[static_cast<size_t>(c) * nBlk + n + i] =
                (Bi != 0) ? -Bi[static_cast<size_t>(c) * n + i] : 0.0;
        }
    }

    const int maxfct = 1, mnum = 1, mtype = -2, msglvl = 0;
    int error = 0, idum = 0;
    int phase = 33;
    const double t0 = nowSeconds();
    pardiso(ptBlk, &maxfct, &mnum, &mtype, &phase, &nBlk,
            aBlk, iaBlk, jaBlk, &idum, &nrhs, iparmBlk, &msglvl,
            rhs, sol, &error);
    tSolve = nowSeconds() - t0;

    if (error == 0) {
        for (int c = 0; c < nrhs; c++) {
            for (int i = 0; i < n; i++) {
                Xr[static_cast<size_t>(c) * n + i] =
                    sol[static_cast<size_t>(c) * nBlk + i];
                Xi[static_cast<size_t>(c) * n + i] =
                    sol[static_cast<size_t>(c) * nBlk + n + i];
            }
        }
    } else {
        opserr << "LadrunoBlockZKernel::solveBlock -- PARDISO solve error "
               << error << "\n";
    }
    delete[] rhs;
    delete[] sol;
    return error;
#endif
}

int
LadrunoBlockZKernel::solveComplex(int nrhs, const double *Br,
                                  const double *Bi, double *Xr, double *Xi)
{
#if !defined(_WIN32) && !defined(_LADRUNO_MKL_FEAST)
    (void)nrhs; (void)Br; (void)Bi; (void)Xr; (void)Xi;
    return -1;
#else
    if (!analyzedZ) {
        opserr << "LadrunoBlockZKernel::solveComplex -- no factorization\n";
        return -1;
    }

    double *rhs = new double[2 * static_cast<size_t>(n) * nrhs];
    double *sol = new double[2 * static_cast<size_t>(n) * nrhs];
    for (int c = 0; c < nrhs; c++) {
        for (int i = 0; i < n; i++) {
            rhs[2 * (static_cast<size_t>(c) * n + i)] =
                Br[static_cast<size_t>(c) * n + i];
            rhs[2 * (static_cast<size_t>(c) * n + i) + 1] =
                (Bi != 0) ? Bi[static_cast<size_t>(c) * n + i] : 0.0;
        }
    }

    const int maxfct = 1, mnum = 1, mtype = 6, msglvl = 0;
    int error = 0, idum = 0;
    int phase = 33;
    const double t0 = nowSeconds();
    pardiso(ptZ, &maxfct, &mnum, &mtype, &phase, &n,
            aZ, iaZ, jaZ, &idum, &nrhs, iparmZ, &msglvl,
            rhs, sol, &error);
    tSolve = nowSeconds() - t0;

    if (error == 0) {
        for (int c = 0; c < nrhs; c++) {
            for (int i = 0; i < n; i++) {
                Xr[static_cast<size_t>(c) * n + i] =
                    sol[2 * (static_cast<size_t>(c) * n + i)];
                Xi[static_cast<size_t>(c) * n + i] =
                    sol[2 * (static_cast<size_t>(c) * n + i) + 1];
            }
        }
    } else {
        opserr << "LadrunoBlockZKernel::solveComplex -- PARDISO solve error "
               << error << "\n";
    }
    delete[] rhs;
    delete[] sol;
    return error;
#endif
}

void
LadrunoBlockZKernel::matvecZ(double a, double b, const double *xr,
                             const double *xi, double *yr, double *yi) const
{
    // y = (zM - K) x on the FULL pattern, complex arithmetic:
    //   yr = (aM-K) xr - (bM) xi
    //   yi = (bM) xr + (aM-K) xi
    for (int i = 0; i < n; i++) {
        double sr = 0.0, si = 0.0;
        for (int k = rowPtr[i] - 1; k < rowPtr[i + 1] - 1; k++) {
            const int j = colInd[k] - 1;
            const double amk = a * Mvals[k] - Kvals[k];
            const double bm = b * Mvals[k];
            sr += amk * xr[j] - bm * xi[j];
            si += bm * xr[j] + amk * xi[j];
        }
        yr[i] = sr;
        yi[i] = si;
    }
}
