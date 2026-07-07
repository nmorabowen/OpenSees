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

// Ladruno: complex / state-space modal analysis driver (ADR 46) — P0 kernel.
// See LadrunoComplexEigen.h for the linearization and sign conventions.

#include "LadrunoComplexEigen.h"

#include <Matrix.h>
#include <OPS_Globals.h>

#include <cmath>
#include <algorithm>

#ifdef _WIN32

extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);

#endif

namespace {

// || (lambda^2 Mt + lambda Ct + Kt) z ||_2 for a complex (lambda, z).
double quadResidual(const Matrix &Mt, const Matrix &Ct, const Matrix &Kt,
                    double lr, double li,
                    const std::vector<double> &zRe,
                    const std::vector<double> &zIm)
{
    const int p = Mt.noRows();
    // lambda^2 = (lr^2 - li^2) + i (2 lr li)
    const double l2r = lr * lr - li * li;
    const double l2i = 2.0 * lr * li;
    double norm2 = 0.0;
    for (int i = 0; i < p; i++) {
        double rr = 0.0, ri = 0.0;
        for (int j = 0; j < p; j++) {
            const double m = Mt(i, j), c = Ct(i, j), k = Kt(i, j);
            const double ar = l2r * m + lr * c + k;   // Re of the (i,j) operator
            const double ai = l2i * m + li * c;       // Im of the (i,j) operator
            rr += ar * zRe[j] - ai * zIm[j];
            ri += ar * zIm[j] + ai * zRe[j];
        }
        norm2 += rr * rr + ri * ri;
    }
    return std::sqrt(norm2);
}

} // namespace

int LadrunoComplexEigen::solveReducedPencil(const Matrix &Mt, const Matrix &Ct,
                                            const Matrix &Kt,
                                            std::vector<ComplexMode> &modes,
                                            double tol)
{
    modes.clear();

    const int p = Mt.noRows();
    if (p <= 0 || Mt.noCols() != p ||
        Ct.noRows() != p || Ct.noCols() != p ||
        Kt.noRows() != p || Kt.noCols() != p) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - Mt/Ct/Kt must be "
               << "square with a common size; got " << Mt.noRows() << "x"
               << Mt.noCols() << ", " << Ct.noRows() << "x" << Ct.noCols()
               << ", " << Kt.noRows() << "x" << Kt.noCols() << "\n";
        return -1;
    }
    if (tol <= 0.0)
        tol = 1.0e-8;

    // Build the 2p x 2p pair (-B, A) column-major for LAPACK:
    //   (-B) = [-Kt 0; 0 Mt],   A = [Ct Mt; Mt 0],   (-B) w = lambda A w.
    const int n = 2 * p;
    std::vector<double> negB(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> Amat(static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            const size_t colTop = static_cast<size_t>(j) * n;       // columns 0..p-1
            const size_t colBot = static_cast<size_t>(j + p) * n;   // columns p..2p-1
            negB[colTop + i]     = -Kt(i, j);
            negB[colBot + (i+p)] =  Mt(i, j);
            Amat[colTop + i]     =  Ct(i, j);
            Amat[colBot + i]     =  Mt(i, j);
            Amat[colTop + (i+p)] =  Mt(i, j);
        }
    }

    std::vector<double> alphaR(n), alphaI(n), beta(n);
    std::vector<double> vr(static_cast<size_t>(n) * n);
    double vlDummy = 0.0;

    char jobvl = 'N';
    char jobvr = 'V';
    int  N = n, ldA = n, ldB = n, ldvl = 1, ldvr = n, info = 0;

    // workspace query, then solve
    int lwork = -1;
    double workQuery = 0.0;
#ifdef _WIN32
    DGGEV(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
          alphaR.data(), alphaI.data(), beta.data(),
          &vlDummy, &ldvl, vr.data(), &ldvr, &workQuery, &lwork, &info);
#else
    dggev_(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
           alphaR.data(), alphaI.data(), beta.data(),
           &vlDummy, &ldvl, vr.data(), &ldvr, &workQuery, &lwork, &info);
#endif
    if (info != 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - dggev workspace "
               << "query failed, INFO = " << info << "\n";
        return -2;
    }
    lwork = static_cast<int>(workQuery);
    if (lwork < 8 * n)
        lwork = 8 * n;
    std::vector<double> work(static_cast<size_t>(lwork));

#ifdef _WIN32
    DGGEV(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
          alphaR.data(), alphaI.data(), beta.data(),
          &vlDummy, &ldvl, vr.data(), &ldvr, work.data(), &lwork, &info);
#else
    dggev_(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
           alphaR.data(), alphaI.data(), beta.data(),
           &vlDummy, &ldvl, vr.data(), &ldvr, work.data(), &lwork, &info);
#endif
    if (info < 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - invalid argument "
               << -info << " passed to LAPACK dggev\n";
        return -2;
    }
    if (info > 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - LAPACK dggev "
               << "failed to converge, INFO = " << info << " (of " << n
               << " roots)\n";
        return -3;
    }

    // scale for the beta ~ 0 (infinite eigenvalue) test
    double alphaScale = 0.0, betaScale = 0.0;
    for (int k = 0; k < n; k++) {
        alphaScale = std::max(alphaScale,
                              std::hypot(alphaR[k], alphaI[k]));
        betaScale  = std::max(betaScale, std::fabs(beta[k]));
    }
    const double betaTol = tol * std::max(betaScale, alphaScale);

    // First pass: recover lambda_k, find the frequency scale for the
    // rigid-mode (omega0 ~ 0) test.
    std::vector<double> lamRe(n, 0.0), lamIm(n, 0.0);
    std::vector<bool>   infinite(n, false);
    double omegaScale = 0.0;
    for (int k = 0; k < n; k++) {
        if (std::fabs(beta[k]) <= betaTol) {
            infinite[k] = true;
            continue;
        }
        lamRe[k] = alphaR[k] / beta[k];
        lamIm[k] = alphaI[k] / beta[k];
        omegaScale = std::max(omegaScale, std::hypot(lamRe[k], lamIm[k]));
    }
    const double omegaTol = tol * std::max(omegaScale, 1.0);

    // Second pass: classify, collapse conjugate pairs, extract z.
    // LAPACK's own pair convention drives the storage layout: alphaI(k) != 0
    // marks a conjugate pair stored in consecutive columns k (real part of v,
    // eigenvalue with POSITIVE imaginary part) and k+1 (imag part of v, the
    // conjugate). The engineering tolerance omegaTol only picks the reported
    // kind, never the column layout.
    for (int k = 0; k < n; k++) {
        ComplexMode mode;
        mode.zRe.assign(static_cast<size_t>(p), 0.0);
        mode.zIm.assign(static_cast<size_t>(p), 0.0);

        if (infinite[k]) {
            // pathological (Mt singular): tag, don't divide
            mode.omega0 = 0.0;  mode.omegaD = 0.0;  mode.zeta = 0.0;
            mode.lambdaRe = 0.0; mode.lambdaIm = 0.0;
            mode.kind = RIGID;   mode.resid = 0.0;
            for (int i = 0; i < p; i++)
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
            modes.push_back(mode);
            continue;
        }

        const double lr = lamRe[k], li = lamIm[k];
        const double mag = std::hypot(lr, li);

        if (alphaI[k] != 0.0 && k + 1 < n) {
            // conjugate pair: report once (the +Im representative),
            // consume both columns
            mode.lambdaRe = lr;
            mode.lambdaIm = std::fabs(li);
            mode.omega0 = mag;
            mode.omegaD = mode.lambdaIm;
            mode.zeta   = (mag > 0.0) ? -lr / mag : 0.0;
            // near-critically-damped pairs (|Im| below the engineering
            // tolerance) are tagged OVERDAMPED but keep the pair's vector
            mode.kind = (mode.lambdaIm > omegaTol) ? UNDERDAMPED : OVERDAMPED;
            // v = VR(:,k) + i VR(:,k+1) pairs with lambda = alpha(k)/beta(k).
            // LAPACK fixes the pair layout by alphaI(k) > 0, NOT by the sign
            // of lambda: when beta(k) < 0 the column-k eigenvalue has
            // NEGATIVE imaginary part, so pairing its vector with the
            // reported +Im lambda requires the conjugate (negate zIm).
            const double imSign = (li < 0.0) ? -1.0 : 1.0;
            for (int i = 0; i < p; i++) {
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
                mode.zIm[static_cast<size_t>(i)] =
                    imSign * vr[static_cast<size_t>(k + 1) * n + i];
            }
            k++;   // skip the conjugate column
        } else {
            // exactly-real root: rigid (lambda ~ 0) or overdamped,
            // one entry per root
            mode.lambdaRe = lr;  mode.lambdaIm = 0.0;
            mode.omegaD = 0.0;
            if (mag <= omegaTol) {
                mode.kind = RIGID;
                mode.omega0 = 0.0;
                mode.zeta = 0.0;
            } else {
                mode.kind = OVERDAMPED;
                mode.omega0 = mag;
                mode.zeta = -lr / mag;   // +1 decaying, -1 growing
            }
            for (int i = 0; i < p; i++)
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
        }

        // phase normalization: largest-|z_j| entry made real-positive
        int    jmax = 0;
        double amax = -1.0;
        for (int j = 0; j < p; j++) {
            const double a = std::hypot(mode.zRe[static_cast<size_t>(j)],
                                        mode.zIm[static_cast<size_t>(j)]);
            if (a > amax) { amax = a; jmax = j; }
        }
        if (amax > 0.0) {
            const double cr = mode.zRe[static_cast<size_t>(jmax)] / amax;
            const double ci = mode.zIm[static_cast<size_t>(jmax)] / amax;
            for (int j = 0; j < p; j++) {          // z <- z * conj(c)
                const double zr = mode.zRe[static_cast<size_t>(j)];
                const double zi = mode.zIm[static_cast<size_t>(j)];
                mode.zRe[static_cast<size_t>(j)] =  zr * cr + zi * ci;
                mode.zIm[static_cast<size_t>(j)] = -zr * ci + zi * cr;
            }
        }

        mode.resid = quadResidual(Mt, Ct, Kt, mode.lambdaRe, mode.lambdaIm,
                                  mode.zRe, mode.zIm);
        modes.push_back(mode);
    }

    std::sort(modes.begin(), modes.end(),
              [](const ComplexMode &a, const ComplexMode &b) {
                  if (a.omega0 != b.omega0)
                      return a.omega0 < b.omega0;
                  return a.lambdaRe < b.lambdaRe;
              });

    return 0;
}
