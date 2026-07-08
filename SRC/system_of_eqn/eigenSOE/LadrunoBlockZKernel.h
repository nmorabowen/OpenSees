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

// Ladruno: block-real complex-shift solve kernel (ADR 43, P3b). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md (S5.2, S9 R1)
// and Ladruno_implementation/feast_d2_spike/README.md.
//
// FEAST's contour solves are genuinely complex: each quadrature node z = a+ib
// needs (z M - K) X = B factored and solved with K, M real symmetric. This
// kernel solves that WITHOUT a complex factorization library, via the
// symmetric equivalent-real ("K-formulation") block system
//
//     [ aM-K     -bM    ] [Xr]   [ Br]
//     [ -bM   -(aM-K)   ] [Xi] = [-Bi]        (2n x 2n, symmetric indefinite)
//
// factored LDL^T. Serial engine = MKL PARDISO mtype -2 (same library the P2
// -certify certificate uses; the distributed MUMPS SYM=2 variant rides the
// P3a sub-communicator plumbing at P3c). A native complex-symmetric path
// (PARDISO mtype 6) is included as the accuracy/cost baseline the P3b gate
// must measure against (D2 review follow-up: the equivalent-real formulation
// is known to degrade vs the native complex solve as Im(z) -> 0).
//
// Symbolic-factorization reuse across shifts (ADR R5): the block pattern
// depends only on the (K,M) sparsity, so PARDISO's analysis phase runs once
// per kernel lifetime; each new shift only re-runs the numeric phase.
//
// MKL-only (this fork: Windows/oneAPI build) — the .cpp guards on _WIN32
// exactly like FeastEigenSolver; elsewhere every method returns -1.

#ifndef LadrunoBlockZKernel_h
#define LadrunoBlockZKernel_h

class LadrunoBlockZKernel
{
  public:
    // CSR views borrowed from FeastEigenSOE: n equations, FULL symmetric
    // pattern, 1-based indices, ascending columns per row. The kernel copies
    // the pattern (building its own block/upper-triangle structures) but
    // reads Kvals/Mvals at setShift time — the SOE must outlive the kernel.
    LadrunoBlockZKernel(int n, const int *rowPtr, const int *colInd,
                        const double *Kvals, const double *Mvals);
    ~LadrunoBlockZKernel();

    // (Re)assemble + numerically factor for shift z = a + i*b.
    // First call also runs the (reused) PARDISO analysis phase.
    // Returns 0, or the PARDISO error code.
    int setShiftBlock(double a, double b);      // 2n block-real, mtype -2
    int setShiftComplex(double a, double b);    // n complex-sym, mtype  6

    // Solve (zM - K) X = B for nrhs columns (column-major, length n each).
    // Uses whichever factorization setShift*() produced last for that path.
    // Xr/Xi may alias nothing; Bi == 0 (nullptr) means a real B.
    int solveBlock(int nrhs, const double *Br, const double *Bi,
                   double *Xr, double *Xi);
    int solveComplex(int nrhs, const double *Br, const double *Bi,
                     double *Xr, double *Xi);

    // y = (zM - K) x in complex arithmetic on the full pattern — for
    // building known-solution right-hand sides and true residuals.
    void matvecZ(double a, double b, const double *xr, const double *xi,
                 double *yr, double *yi) const;

    // wall-clock seconds spent in the last numeric factorization / solve
    double lastFactorSeconds(void) const {return tFactor;}
    double lastSolveSeconds(void)  const {return tSolve;}

  private:
    int n;
    const int *rowPtr;      // borrowed (SOE-owned)
    const int *colInd;
    const double *Kvals;
    const double *Mvals;

    // --- block-real path (2n, upper triangle, 1-based CSR) ---------------
    int nBlk;               // 2n
    int nnzBlk;
    int *iaBlk;             // nBlk+1
    int *jaBlk;             // nnzBlk
    double *aBlk;           // nnzBlk (values for the current shift)
    int *srcSlot;           // nnzBlk -> source slot k in the n-CSR
    signed char *srcCode;   // 0: aM-K   1: -bM   2: -(aM-K)
    void *ptBlk[64];
    int iparmBlk[64];
    bool analyzedBlk;

    // --- complex-symmetric path (n, upper triangle, 1-based CSR) ---------
    int nnzUp;
    int *iaZ;               // n+1
    int *jaZ;               // nnzUp
    double *aZ;             // 2*nnzUp interleaved (re,im) — MKL_Complex16
    int *srcSlotZ;          // nnzUp -> source slot in the n-CSR
    void *ptZ[64];
    int iparmZ[64];
    bool analyzedZ;

    double tFactor;
    double tSolve;

    void releaseBlock(void);
    void releaseComplex(void);
};

#endif
