/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Written: nmb (UANDES), 05/2026
//
// See CriticalTimeStep.h. Shared per-element critical-time-step eigensolve for
// the explicit Noh-Bathe integrators (hoisted out of ExplicitBathe.cpp).

#include <CriticalTimeStep.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <OPS_Globals.h>

#include <cmath>
#include <limits>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>
#endif

// LAPACK: general (DGGEV) and symmetric-definite expert (DSYGVX) generalized
// eigensolvers. The Windows OpenSees LAPACK exports uppercase, no-underscore
// names; the unix reference build exports lowercase, trailing-underscore names.
//
// NOTE (Ladruno): we use DSYGVX (the expert/range driver), NOT the simpler
// DSYGV. The bundled OTHER/LAPACK reference library shipped with this fork (the
// one the Unix/Ubuntu Zone-A build statically links) provides dsygvx.f but does
// NOT provide dsygv.f, so calling dsygv_ produces an "undefined reference"
// link error on Linux even though MKL resolves DSYGV on Windows. DSYGVX is
// already used (and proven to link) by SymmGeneralizedEigenSolver.cpp. With
// RANGE='I' and IL=IU=N it returns just the single largest eigenvalue, which is
// exactly what the critical-time-step estimate needs.
#ifdef _WIN32
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);
extern "C" int DSYGVX(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO,
                      int *N, double *A, int *LDA, double *B, int *LDB,
                      double *VL, double *VU, int *IL, int *IU,
                      double *ABSTOL, int *M, double *W, double *Z,
                      int *LDZ, double *WORK, int *LWORK, int *IWORK,
                      int *IFAIL, int *INFO);
#define LAPACK_DGGEV DGGEV
#define LAPACK_DSYGVX DSYGVX
#else
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);
extern "C" int dsygvx_(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO,
                       int *N, double *A, int *LDA, double *B, int *LDB,
                       double *VL, double *VU, int *IL, int *IU,
                       double *ABSTOL, int *M, double *W, double *Z,
                       int *LDZ, double *WORK, int *LWORK, int *IWORK,
                       int *IFAIL, int *INFO);
#define LAPACK_DGGEV dggev_
#define LAPACK_DSYGVX dsygvx_
#endif

// Largest generalized eigenvalue lambda_max of  K v = lambda M v, with M the
// diagonal (lumped) mass already packed column-major in M_data (only the
// diagonal is non-zero) and K packed column-major in K_data. Returns < 0 on
// failure. Tries DSYGVX (symmetric-definite, ~2x cheaper, no complex-eigenvalue
// fragility), requesting only the single largest eigenvalue; on failure (M not
// positive definite) falls back to DGGEV with a relative-beta threshold to
// discard near-singular eigenpairs.
static double maxGeneralizedEigenvalue(int n, double *K_data, double *M_data)
{
    // --- attempt 1: DSYGVX (needs M positive definite) --------------------
    // DSYGVX overwrites its A/B arguments, so work on copies; K_data/M_data
    // stay intact for the DGGEV fallback. With RANGE='I' and IL=IU=n we ask
    // only for the largest of the n ascending eigenvalues.
    bool mPositive = true;
    for (int i = 0; i < n; ++i)
        if (M_data[i * n + i] <= 0.0) { mPositive = false; break; }

    if (mPositive) {
        double *A = new double[n * n];
        double *B = new double[n * n];
        double *w = new double[n];        // found eigenvalues (need w[0] only)
        double *z = new double[n];        // eigenvectors not requested (jobz='N')
        int    *iwork = new int[5 * n];
        int    *ifail = new int[n];
        for (int k = 0; k < n * n; ++k) { A[k] = K_data[k]; B[k] = M_data[k]; }

        int itype = 1;      // A x = lambda B x
        char jobz = 'N';    // eigenvalues only
        char range = 'I';   // eigenvalues selected by index
        char uplo = 'U';
        int lda = n, ldb = n, ldz = n;
        int il = n, iu = n; // largest eigenvalue only (1-based, ascending)
        double vl = 0.0, vu = 0.0;  // unused when range = 'I'
        double abstol = 0.0;        // use default tolerance
        int m = 0, info = 0;

        int lwork = -1; double wkopt = 0.0;
        LAPACK_DSYGVX(&itype, &jobz, &range, &uplo, &n, A, &lda, B, &ldb,
                      &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz,
                      &wkopt, &lwork, iwork, ifail, &info);
        lwork = (info == 0) ? (int)wkopt : (8 * n);
        if (lwork < 8 * n) lwork = 8 * n;
        double *work = new double[lwork];
        LAPACK_DSYGVX(&itype, &jobz, &range, &uplo, &n, A, &lda, B, &ldb,
                      &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz,
                      work, &lwork, iwork, ifail, &info);

        double lambdaMax = -1.0;
        if (info == 0 && m >= 1) {
            // the single requested (largest) eigenvalue lands in w[0]; the
            // eigenvalues of a real symmetric-definite pencil are real.
            lambdaMax = w[0];
        }

        delete[] A; delete[] B; delete[] w; delete[] z;
        delete[] iwork; delete[] ifail; delete[] work;
        if (info == 0 && m >= 1)
            return lambdaMax;
        // else: B was not positive definite after all -> fall through to DGGEV.
    }

    // --- attempt 2: DGGEV (general, robust to indefinite/singular M) ------
    double *A = new double[n * n];
    double *B = new double[n * n];
    for (int k = 0; k < n * n; ++k) { A[k] = K_data[k]; B[k] = M_data[k]; }

    double *alphar = new double[n];
    double *alphai = new double[n];
    double *beta   = new double[n];
    double *vl = nullptr, *vr = nullptr;
    char jobvl = 'N', jobvr = 'N';
    int lda = n, ldb = n, info = 0;

    int lwork = -1; double wkopt = 0.0;
    LAPACK_DGGEV(&jobvl, &jobvr, &n, A, &lda, B, &ldb,
                 alphar, alphai, beta, vl, &lda, vr, &ldb, &wkopt, &lwork, &info);
    lwork = (info == 0) ? (int)wkopt : (8 * n);
    if (lwork < 1) lwork = 1;
    double *work = new double[lwork];
    LAPACK_DGGEV(&jobvl, &jobvr, &n, A, &lda, B, &ldb,
                 alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);

    double lambdaMax = -1.0;
    if (info == 0) {
        // Relative-beta threshold: a tiny beta means a near-infinite (spurious)
        // eigenvalue from a massless mode. Reject pairs whose |beta| is small
        // relative to the largest |beta| seen, rather than only exactly zero.
        double betaMax = 0.0;
        for (int i = 0; i < n; ++i)
            if (std::abs(beta[i]) > betaMax) betaMax = std::abs(beta[i]);
        const double betaTol = 1.0e-12 * betaMax;

        for (int i = 0; i < n; ++i) {
            if (std::abs(beta[i]) > betaTol) {
                double lambda = alphar[i] / beta[i];   // imaginary part ~ 0 for SPD pencil
                if (lambda > lambdaMax) lambdaMax = lambda;
            }
        }
    }

    delete[] A; delete[] B;
    delete[] alphar; delete[] alphai; delete[] beta; delete[] work;
    return lambdaMax;
}

CTSResult computeCriticalTimeStep(AnalysisModel *theModel,
                                  bool useTangent,
                                  CTSLumping lumping)
{
    CTSResult r;
    r.damped_dt      = std::numeric_limits<double>::infinity();
    r.undamped_dt    = std::numeric_limits<double>::infinity();
    r.damped_tag     = -1;
    r.undamped_tag   = -1;
    r.n_contributing = 0;
    r.n_scanned      = 0;

    if (theModel == 0) return r;
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) return r;

    Element *ele;
    ElementIter &elements = theDomain->getElements();

    while ((ele = elements()) != 0) {
        // --- mass first, copied to a diagonal vector ----------------------
        // Read and lump the mass BEFORE fetching the stiffness: getMass() and
        // getInitialStiff()/getTangentStiff() may both return a reference to the
        // SAME shared static matrix, so holding both at once aliases (D8).
        const Matrix &M = ele->getMass();
        int n = M.noRows();
        if (n == 0) { r.n_scanned++; continue; }

        double *mdiag = new double[n];
        if (lumping == CTSLumping::Diagonal) {
            // diagonal-of-consistent mass: M_ii (strictly positive, dimensionally
            // consistent for rotational DOFs).
            for (int i = 0; i < n; ++i)
                mdiag[i] = M(i, i);
        } else {
            // row-sum lumping: sum of each row onto the diagonal.
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int j = 0; j < n; ++j)
                    sum += M(i, j);
                mdiag[i] = sum;
            }
        }

        // --- now the stiffness (M reference no longer used) ---------------
        const Matrix &K = useTangent ? ele->getTangentStiff() : ele->getInitialStiff();
        if (K.noRows() != n) { delete[] mdiag; r.n_scanned++; continue; }

        // pack column-major for LAPACK
        double *K_data = new double[n * n];
        double *M_data = new double[n * n];
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                K_data[j * n + i] = K(i, j);
                M_data[j * n + i] = (i == j) ? mdiag[i] : 0.0;
            }
        }

        double lambdaMax = maxGeneralizedEigenvalue(n, K_data, M_data);

        delete[] mdiag;
        delete[] K_data;
        delete[] M_data;
        r.n_scanned++;

        if (lambdaMax > 0.0) {
            double w_max = std::sqrt(lambdaMax);

            Vector coefs = ele->getRayleighDampingFactors();
            double alphaM = coefs(0);
            double betaK  = coefs(1);
            double xi = 0.5 * (alphaM / w_max + betaK * w_max);

            double undamped_dt = 2.0 / w_max;
            double damped_dt   = 2.0 / w_max * (std::sqrt(1.0 + xi * xi) - xi);

            r.n_contributing++;
            if (damped_dt < r.damped_dt) {
                r.damped_dt  = damped_dt;
                r.damped_tag = ele->getTag();
            }
            if (undamped_dt < r.undamped_dt) {
                r.undamped_dt  = undamped_dt;
                r.undamped_tag = ele->getTag();
            }
        }
    }

    // --- parallel: governing step is the global minimum over all ranks ----
    // (ADR item 7). The element tags stay rank-local; only the scalar limits are
    // reduced. No-op in the sequential build.
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    {
        double local[2] = { r.damped_dt, r.undamped_dt };
        double global[2] = { local[0], local[1] };
        MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        // if our local value is no longer the governing one, the tag we hold is
        // not the global critical element -> mark it unknown to avoid lying.
        if (global[0] != r.damped_dt)   r.damped_tag   = -1;
        if (global[1] != r.undamped_dt) r.undamped_tag = -1;
        r.damped_dt   = global[0];
        r.undamped_dt = global[1];
    }
#endif

    return r;
}
