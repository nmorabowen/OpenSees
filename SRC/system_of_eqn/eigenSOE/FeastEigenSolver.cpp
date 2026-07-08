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

// Ladruno: FEAST band-targeted eigensolver — solver twin (ADR 43, P1). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md.

#include <FeastEigenSolver.h>
#include <FeastEigenSOE.h>
#include <LadrunoBlockZKernel.h>   // Ladruno ADR43 P3b: -blockZGate kernel
#include <LadrunoFeastInnerSolve.h> // Ladruno ADR43 P3c: RCI inner-solve seam
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

// ---------------------------------------------------------------------------
// MKL Extended Eigensolver (FEAST). MKL is on the link line only on the
// Windows/oneAPI build of this fork (the Ubuntu Zone-A CI links reference
// BLAS/LAPACK via conan, no MKL), so declare + call the kernels under _WIN32
// and provide a compiling stub elsewhere. MKL_INT == int here because the
// Windows build links the LP64 interface (-DMKL_INTERFACE_FULL=intel_lp64);
// the C-callable FEAST names carry no trailing underscore.
// ---------------------------------------------------------------------------
#ifdef _WIN32
// P3c: MKL_Complex16-layout struct for the dfeast_srci RCI (mkl_types.h is
// not included here — same minimal-declaration style as the kernels below)
typedef struct { double real; double imag; } FeastComplex16;

extern "C" {
  void feastinit(int *fpm);
  void dfeast_scsrgv(const char *uplo, const int *n,
                     const double *sa, const int *isa, const int *jsa,
                     const double *sb, const int *isb, const int *jsb,
                     int *fpm, double *epsout, int *loop,
                     const double *emin, const double *emax,
                     int *m0, double *e, double *x, int *m,
                     double *res, int *info);
  // P3c: FEAST reverse-communication interface (mkl_solvers_ee.h contract,
  // LP64: MKL_INT == int). No matrices cross this seam — the caller owns
  // the factorize/solve/matvec tasks the ijob protocol requests.
  void dfeast_srci(int *ijob, const int *n, FeastComplex16 *ze,
                   double *work, FeastComplex16 *workc,
                   double *aq, double *sq, int *fpm,
                   double *epsout, int *loop,
                   const double *emin, const double *emax,
                   int *m0, double *e, double *x, int *m,
                   double *res, int *info);
  // P2: PARDISO LDL^T for the Sturm/inertia certificate (mkl_pardiso.h;
  // _MKL_DSS_HANDLE_t is an opaque void*[64], MKL_INT = int under LP64)
  void pardiso(void *pt, const int *maxfct, const int *mnum,
               const int *mtype, const int *phase, const int *n,
               const void *a, const int *ia, const int *ja,
               int *perm, const int *nrhs, int *iparm,
               const int *msglvl, void *b, void *x, int *error);
}

namespace {

// Number of eigenvalues of the pencil (K, M) strictly below sigma =
// the number of NEGATIVE pivots of the LDL^T factorization of (K - sigma*M)
// (Sylvester inertia). Assembles the UPPER triangle of A = K - sigma*M from
// the SOE's shared full CSR pattern and factors it with MKL PARDISO
// (mtype = -2, real symmetric indefinite); the count is iparm[22]
// (1-based Fortran iparm(23): "number of negative eigenvalues").
// Returns 0 on success (negOut filled), non-zero PARDISO error otherwise.
int
inertiaBelow(int n, const int *rowPtr, const int *colInd,
             const double *Kvals, const double *Mvals, double sigma,
             int &negOut, int &perturbedOut)
{
    // upper-triangle 1-based CSR of A = K - sigma*M
    std::vector<int> ia(static_cast<size_t>(n) + 1);
    std::vector<int> ja;
    std::vector<double> a;
    ja.reserve((static_cast<size_t>(rowPtr[n]) - 1 + n) / 2 + n);
    a.reserve(ja.capacity());
    ia[0] = 1;
    for (int row = 0; row < n; row++) {
        for (int k = rowPtr[row] - 1; k < rowPtr[row + 1] - 1; k++) {
            const int col = colInd[k] - 1;
            if (col >= row) {
                ja.push_back(col + 1);
                a.push_back(Kvals[k] - sigma * Mvals[k]);
            }
        }
        ia[static_cast<size_t>(row) + 1] = static_cast<int>(ja.size()) + 1;
    }

    void *pt[64] = {0};
    int iparm[64] = {0};
    iparm[0]  = 1;    // user-supplied iparm
    iparm[1]  = 2;    // METIS fill-reduction
    iparm[9]  = 8;    // pivot perturbation 1e-8 (indefinite default)
    iparm[10] = 1;    // scaling (congruence — inertia-preserving); reduces
    iparm[12] = 1;    // matching (ditto) how often perturbation triggers
    const int maxfct = 1, mnum = 1, mtype = -2, nrhs = 0, msglvl = 0;
    int phase = 12;   // analysis + numerical factorization
    int error = 0;
    double ddum = 0.0;
    int idum = 0;

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n,
            a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    const int neg = iparm[22];        // 1-based iparm(23): negative eigenvalues
    const int perturbed = iparm[13];  // 1-based iparm(14): perturbed pivots —
                                      // nonzero means a pivot sat below the
                                      // 1e-8 perturbation floor and its SIGN
                                      // (hence the inertia) is not guaranteed

    phase = -1;       // release PARDISO memory regardless of the outcome
    int relErr = 0;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n,
            a.data(), ia.data(), ja.data(), &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &relErr);

    if (error != 0)
        return error;
    negOut = neg;
    perturbedOut = perturbed;
    return 0;
}

// P3b -blockZGate: self-test of the block-real complex-shift solve kernel
// (the FEAST-RCI inner solve of P3c) on the SOE's own CSR. For a shift a
// placed at the widest spectral gap (bounded conditioning by construction)
// it sweeps the contour flattening b = r*frac down to frac = 1e-6, solving
// (zM - K) X = B for a KNOWN X via both the symmetric 2n block-real path
// and the native complex-symmetric PARDISO baseline, and REFUSES if the
// block path's recovered-solution error exceeds 1e-8. A second sweep with
// a ON the lowest found eigenvalue is reported (not asserted) — that is the
// documented equivalent-real degradation regime (D2 review follow-up).
// Timings per factorization are printed for the cost record.
int
runBlockZGate(int n, const int *rowPtr, const int *colInd,
              const double *Kvals, const double *Mvals,
              double emin, double emax, const double *e, int m)
{
    LadrunoBlockZKernel kern(n, rowPtr, colInd, Kvals, Mvals);

    const double r = 0.5 * (emax - emin) > 0.0 ? 0.5 * (emax - emin)
                     : std::max(std::fabs(emax), 1.0);

    // shift placement: midpoint of the widest gap between adjacent found
    // eigenvalues (band edges included), so dist(a, spectrum) >= gap/2 stays
    // bounded away from zero for the asserted sweep
    std::vector<double> pts;
    pts.push_back(emin);
    for (int k = 0; k < m; k++)
        if (e[k] > emin && e[k] < emax)
            pts.push_back(e[k]);
    pts.push_back(emax);
    std::sort(pts.begin(), pts.end());
    double a = 0.5 * (emin + emax), widest = -1.0;
    for (size_t k = 0; k + 1 < pts.size(); k++) {
        const double gap = pts[k + 1] - pts[k];
        if (gap > widest) {
            widest = gap;
            a = 0.5 * (pts[k] + pts[k + 1]);
        }
    }

    // two known solutions: one real, one genuinely complex
    const int nrhs = 2;
    std::vector<double> xkr(static_cast<size_t>(n) * nrhs);
    std::vector<double> xki(static_cast<size_t>(n) * nrhs);
    for (int i = 0; i < n; i++) {
        xkr[i] = std::sin(1.0 + i);
        xki[i] = 0.0;
        xkr[static_cast<size_t>(n) + i] = std::cos(0.5 + i);
        xki[static_cast<size_t>(n) + i] = 0.7 * std::sin(2.0 + i);
    }
    std::vector<double> br(xkr.size()), bi(xkr.size());
    std::vector<double> xr(xkr.size()), xi(xkr.size());

    auto solveAndMeasure = [&](bool blockPath, double aa, double bb,
                               double &errOut, double &tFacOut) -> int {
        for (int c = 0; c < nrhs; c++)
            kern.matvecZ(aa, bb, &xkr[static_cast<size_t>(c) * n],
                         &xki[static_cast<size_t>(c) * n],
                         &br[static_cast<size_t>(c) * n],
                         &bi[static_cast<size_t>(c) * n]);
        int rc = blockPath ? kern.setShiftBlock(aa, bb)
                           : kern.setShiftComplex(aa, bb);
        if (rc != 0)
            return rc;
        tFacOut = kern.lastFactorSeconds();
        rc = blockPath
                 ? kern.solveBlock(nrhs, br.data(), bi.data(), xr.data(),
                                   xi.data())
                 : kern.solveComplex(nrhs, br.data(), bi.data(), xr.data(),
                                     xi.data());
        if (rc != 0)
            return rc;
        double num = 0.0, den = 0.0;
        for (size_t p = 0; p < xkr.size(); p++) {
            const double dr = xr[p] - xkr[p], di = xi[p] - xki[p];
            num = std::max(num, std::sqrt(dr * dr + di * di));
            const double mag = std::sqrt(xkr[p] * xkr[p] + xki[p] * xki[p]);
            den = std::max(den, mag);
        }
        errOut = num / std::max(den, 1e-30);
        return 0;
    };

    const double fracs[] = {1.0, 0.3, 0.1, 1e-2, 1e-3, 1e-4, 1e-6};
    double maxErrBlk = 0.0;
    opserr << "FeastEigenSolver -- -blockZGate sweep (a at widest gap = "
           << a << ", r = " << r << "):\n";
    opserr << "    b/r        err(block)    err(complex)  tFac(blk)  tFac(z)\n";
    for (double frac : fracs) {
        const double b = r * frac;
        double errB = 0.0, errZ = 0.0, tB = 0.0, tZ = 0.0;
        if (solveAndMeasure(true, a, b, errB, tB) != 0 ||
            solveAndMeasure(false, a, b, errZ, tZ) != 0) {
            opserr << "FeastEigenSolver -- -blockZGate: kernel failure at "
                   << "b/r = " << frac << "\n";
            return -1;
        }
        maxErrBlk = std::max(maxErrBlk, errB);
        opserr << "    " << frac << "    " << errB << "    " << errZ
               << "    " << tB << "    " << tZ << "\n";
    }

    // degradation regime (report-only): a ON a found eigenvalue (e[] is
    // unsorted at this point — any in-band eigenvalue works for this test)
    if (m > 0) {
        opserr << "  degradation regime (a on found eigenvalue " << e[0]
               << ", report-only):\n";
        for (double frac : {1e-2, 1e-4, 1e-6}) {
            const double b = r * frac;
            double errB = 0.0, errZ = 0.0, tB = 0.0, tZ = 0.0;
            const int rcB = solveAndMeasure(true, e[0], b, errB, tB);
            const int rcZ = solveAndMeasure(false, e[0], b, errZ, tZ);
            opserr << "    b/r = " << frac << ": err(block) = "
                   << (rcB == 0 ? errB : -1.0) << ", err(complex) = "
                   << (rcZ == 0 ? errZ : -1.0) << "\n";
        }
    }

    if (maxErrBlk > 1.0e-8) {
        opserr << "FeastEigenSolver -- -blockZGate FAILED: block-real "
               << "recovered-solution error " << maxErrBlk
               << " exceeds 1e-8 on the bounded-conditioning sweep\n";
        return -1;
    }
    opserr << "FeastEigenSolver -- -blockZGate PASS: max block-real error "
           << maxErrBlk << " over the b/r sweep down to 1e-6\n";
    return 0;
}

// P3c-serial: drive MKL FEAST through its reverse-communication interface
// (dfeast_srci) with LadrunoBlockZKernel as the inner complex-shift solve —
// the orchestration seam the MPI rung swaps a distributed kernel into.
//
// The ijob protocol (MKL mkl_solvers_ee.h + Intel's dfeast_sparse_rci.c
// canonical example; Fortran fpm(k) == C fpm[k-1]):
//   -1 on entry, then loop until ijob == 0:
//   10 -> factorize (ze*M - K) for the contour node ze (complex);
//         LadrunoBlockZKernel::setShiftBlock — the symmetric 2n block-real
//         LDL^T, PARDISO analysis phase reused across ALL nodes and
//         refinement loops (ADR R5).
//   11 -> solve (ze*M - K) Y = workc(n, m0), solution back INTO workc;
//         the m0 right-hand sides are genuinely complex (D2 Finding 3).
//   30 -> work(:, i:j) = K * q(:, i:j), i = fpm[23]-1 (1-based Fortran
//         fpm(24)), j-i+1 = fpm[24] columns (Fortran fpm(25)).
//   40 -> same with M.
//   -2 -> new refinement loop, nothing to do.
// The half-contour conjugate handling (skip conj(ze), conjugate the
// accumulated projector term) is INTERNAL to MKL FEAST — the caller only
// ever factorizes/solves the ze values it is handed, never a conjugate.
//
// m0 is passed BY VALUE: dfeast_srci's m0 argument is in/out and MKL
// SHRINKS it in place when the band holds fewer modes than the requested
// subspace (observed on the P3c battery: m0 15 -> 5 when the band held 5).
// The shrunken value is RCI-session bookkeeping, not a saturation verdict —
// leaking it back to the caller made the saturation-enlargement loop compare
// m against the shrunken m0 (false "saturated", plateaued retries, refusal).
// The caller's m0 stays the REQUESTED subspace size, exactly like the
// driver path's compare.
//
// Returns 0 when the RCI ran to termination (info carries FEAST's verdict,
// same codes as dfeast_scsrgv), -1 on a kernel/protocol failure (message
// already printed).
int
runFeastRci(LadrunoFeastInnerSolve &kern, int n,
            const int *rowPtr, const int *colInd,
            const double *Kvals, const double *Mvals,
            int *fpm, double *epsout, int *loop,
            const double &emin, const double &emax,
            const int m0v, double *e, double *x, int *m, double *res,
            int *info)
{
    std::vector<double> work(static_cast<size_t>(n) * m0v);
    std::vector<FeastComplex16> workc(static_cast<size_t>(n) * m0v);
    std::vector<double> aq(static_cast<size_t>(m0v) * m0v);
    std::vector<double> sq(static_cast<size_t>(m0v) * m0v);
    // split-complex staging for the kernel's (Br,Bi)->(Xr,Xi) interface
    std::vector<double> br(static_cast<size_t>(n) * m0v);
    std::vector<double> bi(static_cast<size_t>(n) * m0v);
    std::vector<double> xr(static_cast<size_t>(n) * m0v);
    std::vector<double> xi(static_cast<size_t>(n) * m0v);

    int ijob = -1;
    FeastComplex16 ze = {0.0, 0.0};
    int m0Rci = m0v;        // the RCI's own (shrinkable) in/out copy
    *info = 0;

    while (true) {
        dfeast_srci(&ijob, &n, &ze, work.data(), workc.data(),
                    aq.data(), sq.data(), fpm, epsout, loop,
                    &emin, &emax, &m0Rci, e, x, m, res, info);
        if (ijob == 0)
            break;              // finished — info carries FEAST's verdict
        if (*info != 0)
            break;              // FEAST aborted mid-protocol; report info

        switch (ijob) {
        case -2:                // new refinement loop
            break;
        case 10: {              // factorize (ze*M - K) for this contour node
            if (kern.setShiftBlock(ze.real, ze.imag) != 0) {
                opserr << "FeastEigenSolver -- -rci: block-real factorization "
                       << "failed at contour node (" << ze.real << ", "
                       << ze.imag << "i)\n";
                return -1;
            }
            break;
        }
        case 11: {              // solve (ze*M - K) Y = workc, back into workc
            // m0Rci (not m0v): after an MKL-side subspace shrink the live
            // right-hand sides are the first m0Rci columns only
            const size_t tot = static_cast<size_t>(n) * m0Rci;
            for (size_t p = 0; p < tot; p++) {
                br[p] = workc[p].real;
                bi[p] = workc[p].imag;
            }
            if (kern.solveBlock(m0Rci, br.data(), bi.data(),
                                xr.data(), xi.data()) != 0) {
                opserr << "FeastEigenSolver -- -rci: block-real solve failed "
                       << "at contour node (" << ze.real << ", " << ze.imag
                       << "i)\n";
                return -1;
            }
            for (size_t p = 0; p < tot; p++) {
                workc[p].real = xr[p];
                workc[p].imag = xi[p];
            }
            break;
        }
        case 30:                // work(:, i:j) = K * q(:, i:j)
        case 40: {              // work(:, i:j) = M * q(:, i:j)
            const double *vals = (ijob == 30) ? Kvals : Mvals;
            const int c0 = fpm[23] - 1;   // first column, 1-based fpm(24)
            const int nc = fpm[24];       // column count,          fpm(25)
            if (c0 < 0 || nc < 0 || c0 + nc > m0v) {
                opserr << "FeastEigenSolver -- -rci: matvec column window ["
                       << c0 << ", " << c0 + nc << ") outside the m0=" << m0v
                       << " subspace\n";
                return -1;
            }
            for (int c = c0; c < c0 + nc; c++) {
                const double *xc = &x[static_cast<size_t>(c) * n];
                double *yc = &work[static_cast<size_t>(c) * n];
                for (int i = 0; i < n; i++) {
                    double s = 0.0;
                    for (int k = rowPtr[i] - 1; k < rowPtr[i + 1] - 1; k++)
                        s += vals[k] * xc[colInd[k] - 1];
                    yc[i] = s;
                }
            }
            break;
        }
        default:
            opserr << "FeastEigenSolver -- -rci: unexpected ijob " << ijob
                   << " from dfeast_srci\n";
            return -1;
        }
    }

    return 0;
}

} // namespace
#endif

FeastEigenSolver::FeastEigenSolver()
:EigenSolver(EigenSOLVER_TAGS_FeastEigenSolver),
 theSOE(0), numModes(0), eigenvalues(0), eigenvectors(0),
 valueCapacity(0), vectorCapacity(0), eigenV(0)
{

}

FeastEigenSolver::~FeastEigenSolver()
{
  if (eigenvalues  != 0) delete [] eigenvalues;
  if (eigenvectors != 0) delete [] eigenvectors;
  if (eigenV       != 0) delete eigenV;
}

int
FeastEigenSolver::solve(int numModesRequested, bool generalized, bool findSmallest)
{
  if (theSOE == 0) {
    opserr << "FeastEigenSolver::solve() -- no EigenSOE has been set\n";
    return -1;
  }

  if (generalized == false) {
    opserr << "FeastEigenSolver::solve() -- only solves the generalized problem K phi = lambda M phi\n";
    return -1;
  }

#ifndef _WIN32
  opserr << "FeastEigenSolver::solve() -- FEAST requires the MKL Extended "
         << "Eigensolver, available only on the Windows/oneAPI build "
         << "(ADR 43 P1)\n";
  numModes = 0;
  return -1;
#else

  int n = theSOE->size;
  if (n <= 0) {
    numModes = 0;
    return 0;
  }

  // shared CSR pattern: K and M reuse rowPtr/colInd (isa==isb, jsa==jsb)
  const double *sa  = theSOE->Kvals;
  const double *sb  = theSOE->Mvals;
  const int    *isa = theSOE->rowPtr;
  const int    *jsa = theSOE->colInd;
  char uplo = 'F';                          // full pattern stored

  double emin = theSOE->emin;
  double emax = theSOE->emax;
  if (!(emax > emin)) {
    opserr << "FeastEigenSolver::solve() -- invalid band [" << emin << ", "
           << emax << "]; require emax > emin (lambda = (2*pi*f)^2)\n";
    return -1;
  }

  int nq        = theSOE->nq;
  int maxRefine = theSOE->maxRefine;
  int tolExp    = theSOE->tolExp;
  int fpmPrint  = theSOE->verbose ? 1 : 0;

  // ---- seed the subspace size m0 -----------------------------------------
  // 0 => auto: use FEAST's stochastic estimate (fpm[13]=2) to count the modes
  // in the band without factorizing, then over-provision (Kratos heuristic).
  int m0 = theSOE->m0;
  if (m0 <= 0) {
    int fpm[128];
    feastinit(fpm);
    fpm[0]  = fpmPrint;
    if (nq > 0)        fpm[1] = nq;
    fpm[2]  = tolExp;
    if (maxRefine > 0) fpm[3] = maxRefine;
    fpm[13] = 2;                              // stochastic estimate only

    int mGuess = (n < 40) ? n : 40;           // random-vector count for estimate
    double *eE   = new double[mGuess];
    double *xE   = new double[(size_t)n * mGuess];
    double *resE = new double[mGuess];
    double epsout = 0.0;
    int loop = 0, mEst = 0, info = 0;

    dfeast_scsrgv(&uplo, &n, sa, isa, jsa, sb, isa, jsa,
                  fpm, &epsout, &loop, &emin, &emax, &mGuess,
                  eE, xE, &mEst, resE, &info);

    delete [] eE; delete [] xE; delete [] resE;

    if (info != 0)
      opserr << "FeastEigenSolver::solve() -- stochastic mode-count estimate "
             << "failed (info=" << info << "); seeding a minimal subspace "
             << "and relying on the saturation loop\n";
    if (mEst < 0) mEst = 0;
    m0 = (int)(1.5 * mEst) + 8;               // over-provision the subspace
    if (theSOE->verbose)
      opserr << "FeastEigenSolver::solve() -- stochastic estimate ~ " << mEst
             << " mode(s) in band; seeding m0 = " << m0 << "\n";
  }
  if (m0 > n) m0 = n;
  if (m0 < 1) m0 = 1;

  // ---- saturation-enlargement loop ---------------------------------------
  const int maxEnlarge = 4;
  int enlarge = 0;
  int m = 0;
  int info = 0;
  double *e = 0, *x = 0, *res = 0;

  // P3c-serial: -rci routes the solve through the dfeast_srci RCI with the
  // block-real kernel as the inner solve (the seam the MPI rung distributes).
  // ONE kernel for the whole solve: the PARDISO analysis phase runs once and
  // is reused across every contour node, refinement loop, AND saturation
  // retry (the block pattern depends only on the (K,M) sparsity — ADR R5).
  // The m0 auto-seed estimate above deliberately stays on the driver path:
  // it is a count heuristic, not a solve, and keeping it shared makes the
  // -rci/driver eigenpair comparison seed identically.
  // ADR 43 P3c-MPI (L3-only, R0): if the parallel executable target registered
  // a distributed inner-solve factory AND it accepts this run (comm size > 1),
  // the RCI's contour solves factor/solve the 2n block system across ranks
  // (dmumps SYM=2); otherwise the serial PARDISO kernel. Either way the SAME
  // dfeast_srci outer loop runs (replicated across the ranks under MPI), so the
  // distributed spectrum is a bit-for-bit differential of the serial -rci run.
  const bool useRci = theSOE->rciFlag;
  std::unique_ptr<LadrunoFeastInnerSolve> rciKern;
  if (useRci) {
    LadrunoFeastInnerSolve *dist = 0;
    if (ladrunoFeastDistFactory != 0)
      dist = ladrunoFeastDistFactory(n, isa, jsa, sa, sb);
    if (dist != 0)
      rciKern.reset(dist);
    else
      rciKern.reset(new LadrunoBlockZKernel(n, isa, jsa, sa, sb));
  }

  while (true) {
    if (m0 > n) m0 = n;

    if (e   != 0) delete [] e;
    if (x   != 0) delete [] x;
    if (res != 0) delete [] res;
    e   = new double[m0];
    x   = new double[(size_t)n * m0];
    res = new double[m0];

    int fpm[128];
    feastinit(fpm);
    fpm[0] = fpmPrint;
    if (nq > 0)        fpm[1] = nq;
    fpm[2] = tolExp;
    if (maxRefine > 0) fpm[3] = maxRefine;

    double epsout = 0.0;
    int loop = 0;
    m = 0; info = 0;

    if (useRci) {
      if (runFeastRci(*rciKern, n, isa, jsa, sa, sb,
                      fpm, &epsout, &loop, emin, emax, m0,
                      e, x, &m, res, &info) != 0) {
        opserr << "FeastEigenSolver::solve() -- -rci: inner-solve/protocol "
               << "failure; refusing\n";
        delete [] e; delete [] x; delete [] res;
        numModes = 0;
        return -8;
      }
      if (theSOE->verbose)
        opserr << "FeastEigenSolver::solve() -- -rci: dfeast_srci finished, "
               << loop << " refinement loop(s), trace eps " << epsout
               << ", info " << info << "\n";
    } else {
      dfeast_scsrgv(&uplo, &n, sa, isa, jsa, sb, isa, jsa,
                    fpm, &epsout, &loop, &emin, &emax, &m0,
                    e, x, &m, res, &info);
    }

    if (info == 0 || info == 4) {
      // success (info==4: only a stochastic estimate was possible -> treat
      // as resolved). If the accepted count fills the subspace, the band may
      // hold more modes than m0 could capture -> enlarge and re-run.
      if (m >= m0 && m0 < n) {
        if (enlarge < maxEnlarge) {
          opserr << "FeastEigenSolver::solve() -- subspace saturated (m = m0 = "
                 << m0 << "); enlarging and retrying\n";
          m0 = (int)(1.5 * m0) + 8;
          enlarge++;
          continue;
        }
        // gate B2: a still-saturated subspace after maxEnlarge means the
        // band may hold MORE modes than we found — the completeness claim
        // would be silently false. Refuse, exactly like the info==3
        // exhausted path (same physical condition, same honesty).
        opserr << "FeastEigenSolver::solve() -- subspace STILL saturated "
               << "(m = m0 = " << m0 << ") after " << maxEnlarge
               << " enlargements; the band may hold more modes - "
               << "completeness NOT guaranteed. Re-run with an explicit "
               << "larger -m0\n";
        delete [] e; delete [] x; delete [] res;
        numModes = 0;
        return -3;
      }
      break;
    }
    else if (info == 3) {
      // subspace m0 too small
      if (m0 < n && enlarge < maxEnlarge) {
        opserr << "FeastEigenSolver::solve() -- FEAST reports m0 too small (info=3, m0="
               << m0 << "); enlarging and retrying\n";
        m0 = (int)(1.5 * m0) + 8;
        enlarge++;
        continue;
      }
      opserr << "FeastEigenSolver::solve() -- FEAST subspace still too small after "
             << maxEnlarge << " enlargements (info=3); band likely too wide for m0=n="
             << n << "\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -1;
    }
    else if (info == 1) {
      // no eigenvalue in the interval -> not an error: report zero modes
      opserr << "FeastEigenSolver::solve() -- no eigenvalues in band [" << emin
             << ", " << emax << "] (info=1)\n";
      m = 0;
      break;
    }
    else if (info == 2) {
      // gate B3: FEAST exceeded the refinement loops. Do NOT return this as
      // a clean result — gate on the residuals MKL actually computed: accept
      // only if every accepted mode meets the requested tolerance anyway
      // (slow-but-arrived); otherwise fail loudly.
      double maxRes = 0.0;
      for (int k = 0; k < m; k++)
        if (res[k] > maxRes) maxRes = res[k];
      const double tolVal = std::pow(10.0, -(double)tolExp);
      if (m > 0 && maxRes <= tolVal) {
        opserr << "FeastEigenSolver::solve() -- FEAST hit the refinement-loop "
               << "cap (info=2) but every accepted mode meets the tolerance "
               << "(max residual " << maxRes << " <= " << tolVal
               << "); accepting\n";
        break;
      }
      opserr << "FeastEigenSolver::solve() -- FEAST did NOT converge in "
             << "the refinement loop (info=2, max residual " << maxRes
             << " > " << tolVal << "); refusing the unconverged result - "
             << "raise -maxiter or loosen -tol\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -2;
    }
    else {
      // info < 0: internal error; info == 100+i: bad fpm(i); 200..202: bad n/m0/band
      opserr << "FeastEigenSolver::solve() -- MKL dfeast_scsrgv failed, info = "
             << info;
      if (info < 0)
        opserr << " (internal FEAST/inner-solver error)\n";
      else if (info >= 100 && info < 200)
        opserr << " (invalid FEAST parameter fpm[" << (info - 100) << "])\n";
      else if (info == 200)
        opserr << " (invalid search interval emin/emax)\n";
      else if (info == 201)
        opserr << " (invalid subspace size m0)\n";
      else if (info == 202)
        opserr << " (invalid problem size n)\n";
      else
        opserr << "\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return info;
    }
  }

  // ---- store accepted pairs ascending, capped at numModesRequested -------
  if (m < 0) m = 0;

  // the RCI kernel's factorization is dead weight past this point — release
  // it before -blockZGate builds ITS kernel (avoids transiently holding two
  // 2n block LDL^T factorizations; Opus gate MINOR)
  rciKern.reset();

  // ---- P2 -certify: the Sturm/inertia completeness certificate -----------
  // certified band count = neg(K - emax*M) - neg(K - emin*M), an EXACT count
  // from the factorization itself, independent of FEAST's own bookkeeping
  // (the LS-DYNA BCSLIB-EXT guarantee). Edges are nudged outward by a
  // relative delta so a mode sitting numerically ON an edge (which FEAST's
  // closed contour includes) is counted rather than left ill-posed.
  if (theSOE->certifyFlag) {
    const double eref = std::max(std::max(std::fabs(emin), std::fabs(emax)),
                                 1.0);
    const double delta = 1.0e-10 * eref;
    int negLo = 0, negHi = 0, pertLo = 0, pertHi = 0;
    const int errLo = inertiaBelow(n, theSOE->rowPtr, theSOE->colInd,
                                   theSOE->Kvals, theSOE->Mvals,
                                   emin - delta, negLo, pertLo);
    const int errHi = inertiaBelow(n, theSOE->rowPtr, theSOE->colInd,
                                   theSOE->Kvals, theSOE->Mvals,
                                   emax + delta, negHi, pertHi);
    if (errLo != 0 || errHi != 0) {
      opserr << "FeastEigenSolver::solve() -- -certify: PARDISO inertia "
             << "factorization failed (error " << errLo << "/" << errHi
             << "); cannot certify the band count\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -4;
    }
    // adversarial-gate fix: a PERTURBED pivot means (K - sigma*M) had a
    // pivot below the perturbation floor — the shift sits numerically ON an
    // eigenvalue and the pivot's sign (hence the count) is NOT guaranteed.
    // The outward nudge (1e-10 relative, lambda-units) is SMALLER than the
    // perturbation floor (1e-8 relative, stiffness-units), so it cannot
    // protect against this — detect and refuse instead of certifying a
    // possibly-off-by-one count.
    if (pertLo > 0 || pertHi > 0) {
      opserr << "FeastEigenSolver::solve() -- -certify: the factorization "
             << "used " << pertLo << "/" << pertHi << " perturbed pivot(s) "
             << "at the band edges - an eigenvalue sits numerically ON an "
             << "edge and the inertia count is not guaranteed exact. Move "
             << "the band edge away from the eigenvalue and re-run\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -6;
    }
    const int certified = negHi - negLo;
    if (certified != m) {
      opserr << "FeastEigenSolver::solve() -- -certify FAILED: the "
             << "factorization inertia counts " << certified
             << " eigenvalue(s) in the band but FEAST returned " << m
             << " - the result is NOT complete; refusing. (Raise -m0/-nq "
             << "or tighten -tol and re-run)\n";
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -5;
    }
    opserr << "FeastEigenSolver::solve() -- -certify: inertia count "
           << certified << " == FEAST count " << m
           << "; band content certified complete\n";
  }

  // ---- P3b -blockZGate: block-real complex-shift solve self-test ---------
  if (theSOE->blockZGateFlag) {
    if (runBlockZGate(n, theSOE->rowPtr, theSOE->colInd,
                      theSOE->Kvals, theSOE->Mvals, emin, emax, e, m) != 0) {
      delete [] e; delete [] x; delete [] res;
      numModes = 0;
      return -7;
    }
  }

  // sort the m accepted eigenvalues ascending, carrying a vector permutation
  std::vector<int> perm(m);
  for (int k = 0; k < m; k++) perm[k] = k;
  std::sort(perm.begin(), perm.end(),
            [&](int a, int b) { return e[a] < e[b]; });

  int report = m;
  if (numModesRequested > 0 && numModesRequested < report) {
    // gate B1: the band holds MORE modes than the analysis-side report cap
    // — retaining a subset silently would falsify the completeness claim.
    // Warn LOUDLY; the command layer's cap is max(128, -m0), so an explicit
    // larger -m0 raises it.
    opserr << "WARNING FeastEigenSolver::solve() -- the band holds " << m
           << " modes but only " << numModesRequested
           << " are being retained (analysis-side cap). The reported set is "
           << "NOT the complete band content - re-run with -m0 " << m
           << " (or larger) to retain all of them\n";
    report = numModesRequested;
  }

  if (report > valueCapacity) {
    if (eigenvalues != 0) delete [] eigenvalues;
    eigenvalues = new double[report];
    valueCapacity = report;
  }
  if ((size_t)n * report > (size_t)vectorCapacity) {
    if (eigenvectors != 0) delete [] eigenvectors;
    eigenvectors = new double[(size_t)n * report];
    vectorCapacity = n * report;
  }

  for (int k = 0; k < report; k++) {
    int src = perm[k];
    eigenvalues[k] = e[src];
    const double *xcol = &x[(size_t)src * n];
    double *dst = &eigenvectors[(size_t)k * n];
    for (int i = 0; i < n; i++)
      dst[i] = xcol[i];             // FEAST vectors are already M-orthonormal
  }

  numModes = report;

  delete [] e; delete [] x; delete [] res;

  return 0;
#endif
}

int
FeastEigenSolver::setSize(void)
{
  if (theSOE == 0) {
    opserr << "FeastEigenSolver::setSize() -- no EigenSOE has been set\n";
    return -1;
  }

  int n = theSOE->size;
  if (eigenV == 0 || eigenV->Size() != n) {
    if (eigenV != 0)
      delete eigenV;
    eigenV = new Vector(n);
    if (eigenV == 0 || eigenV->Size() != n) {
      opserr << "FeastEigenSolver::setSize() -- ran out of memory for eigenvector of size " << n << endln;
      return -2;
    }
  }

  return 0;
}

int
FeastEigenSolver::setEigenSOE(FeastEigenSOE &theFeastSOE)
{
  theSOE = &theFeastSOE;
  return 0;
}

const Vector &
FeastEigenSolver::getEigenvector(int mode)
{
  if (mode < 1 || mode > numModes) {
    // SILENT zero beyond the found count: the band defines how many modes
    // exist, but the analysis driver loops a fixed requested numEigen — the
    // command layer reconciles the count after the solve (ADR 43 P1), so an
    // over-ask here is expected, not an error.
    if (eigenV != 0)
      eigenV->Zero();
    return *eigenV;
  }

  int n = theSOE->size;
  int index = (mode - 1) * n;
  Vector &vec = *eigenV;
  if (eigenvectors != 0) {
    for (int i = 0; i < n; i++)
      vec(i) = eigenvectors[index++];
  } else {
    opserr << "FeastEigenSolver::getEigenvector() -- eigenvectors not yet computed\n";
    vec.Zero();
  }

  return *eigenV;
}

double
FeastEigenSolver::getEigenvalue(int mode)
{
  if (mode < 1 || mode > numModes) {
    // silent zero beyond the found count — see getEigenvector note
    return 0.0;
  }

  if (eigenvalues != 0)
    return eigenvalues[mode - 1];

  opserr << "FeastEigenSolver::getEigenvalue() -- eigenvalues not yet computed\n";
  return 0.0;
}

int
FeastEigenSolver::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "FeastEigenSolver::sendSelf() -- serial only (ADR 43 P1)\n";
  return -1;
}

int
FeastEigenSolver::recvSelf(int commitTag, Channel &theChannel,
                           FEM_ObjectBroker &theBroker)
{
  opserr << "FeastEigenSolver::recvSelf() -- serial only (ADR 43 P1)\n";
  return -1;
}
