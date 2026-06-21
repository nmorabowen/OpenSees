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

// ADR-30 P1 — LadrunoConstraintProjector implementation. See the header.

#include <LadrunoConstraintProjector.h>
#include <LinearSOE.h>
#include <DiagonalSOE.h>
#include <OPS_Globals.h>
#include <cmath>
#include <vector>
#include <limits>

// Ladruno (ADR-30 P6): condition-number thresholds for the per-group LtML = L^T M L.
// LtML is SPD (positive masses + full-rank L); a near-singular LtML makes project() amplify
// round-off into a garbage acceleration, which the old exact-pivot test-solve never caught.
#define LADRUNO_LTML_COND_REFUSE 1.0e12   // cond above this -> refuse (numerically singular)
#define LADRUNO_LTML_COND_WARN   1.0e8    // cond above this -> warn once (ill-conditioned)

// Ladruno (ADR-30 P6): eigenvalues of a small SYMMETRIC matrix A (n x n) by cyclic Jacobi.
// Fills eig[0..n-1] (unsorted); eigenVECTORS are not needed (used only for a condition-number
// estimate of the tiny SPD LtML, n typically <= 6). Self-contained — no LAPACK dependency, so it
// sidesteps the bundled-LAPACK symbol gaps (cf. the dsygv_ incident). Returns 0 on success.
static int
ladrunoSymEigJacobi(const Matrix &A, Vector &eig)
{
    int n = A.noRows();
    if (n <= 0 || A.noCols() != n) return -1;
    eig = Vector(n);
    if (n == 1) { eig(0) = A(0, 0); return 0; }

    std::vector<double> a((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[(size_t)i * n + j] = A(i, j);

    const int maxSweeps = 100;
    for (int sweep = 0; sweep < maxSweeps; sweep++) {
        double off = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                off += a[(size_t)i * n + j] * a[(size_t)i * n + j];
        if (off <= 1.0e-30) break;

        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = a[(size_t)p * n + q];
                if (apq == 0.0) continue;
                double app = a[(size_t)p * n + p];
                double aqq = a[(size_t)q * n + q];
                // t solves t^2 + 2*theta*t - 1 = 0 with the smaller-magnitude root; the
                // Givens rotation G (G[p][q]=s, G[q][p]=-s, diag c) zeros a'_pq in G^T A G.
                double theta = (aqq - app) / (2.0 * apq);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                           (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double s = t * c;
                // A G (columns p,q)
                for (int k = 0; k < n; k++) {
                    double g = a[(size_t)k * n + p];
                    double h = a[(size_t)k * n + q];
                    a[(size_t)k * n + p] = c * g - s * h;
                    a[(size_t)k * n + q] = s * g + c * h;
                }
                // G^T (...) (rows p,q)
                for (int k = 0; k < n; k++) {
                    double g = a[(size_t)p * n + k];
                    double h = a[(size_t)q * n + k];
                    a[(size_t)p * n + k] = c * g - s * h;
                    a[(size_t)q * n + k] = s * g + c * h;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) eig(i) = a[(size_t)i * n + i];
    return 0;
}

LadrunoConstraintProjector::LadrunoConstraintProjector()
  : groups(), massReady(false), projectICs(false), icTol(1.0e-8)
{
}

int
LadrunoConstraintProjector::enforceIC(Vector &u, Vector &v) const
{
    if (projectICs) {
        this->snapICs(u, v);
        return 0;
    }
    return this->checkIC(u, v, icTol);
}

LadrunoConstraintProjector::~LadrunoConstraintProjector()
{
}

void
LadrunoConstraintProjector::addGroup(const ID &allEqn, int nRet, const Matrix &L,
                                     const Vector &delta, const ID &fixedEqn)
{
    Group g;
    g.allEqn = allEqn;
    g.nRet   = nRet;
    g.L      = L;
    g.delta  = delta;
    g.fixedEqn = fixedEqn;
    g.m      = Vector(allEqn.Size());
    g.LtML   = Matrix(nRet, nRet);
    groups.push_back(g);
    massReady = false;
}

int
LadrunoConstraintProjector::buildMass(LinearSOE *theSOE)
{
    DiagonalSOE *theDiag = dynamic_cast<DiagonalSOE *>(theSOE);
    if (theDiag == 0) {
        opserr << "LadrunoConstraintProjector::buildMass() - the projection handler "
                  "requires `system Diagonal` (the explicit lumped-mass recipe); the "
                  "current SOE is not a DiagonalSOE.\n";
        return -1;
    }

    // The assembled diagonal IS the lumped mass M the integrator inverts. Must be
    // read BEFORE the solver factors it in place (DiagonalDirectSolver overwrites
    // A[i] with 1/A[i]); the integrator calls buildMass between formTangent and solve.
    const double *A = theDiag->getDiagonalA();
    int neq = theDiag->getNumEqn();
    if (A == 0) {
        opserr << "LadrunoConstraintProjector::buildMass() - null diagonal.\n";
        return -1;
    }

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        for (int r = 0; r < nrows; r++) {
            int e = g.allEqn(r);
            if (e < 0 || e >= neq) {
                opserr << "LadrunoConstraintProjector::buildMass() - equation " << e
                       << " out of range [0," << neq << ") in group " << gi << ".\n";
                return -1;
            }
            g.m(r) = A[e];
            if (g.m(r) <= 0.0) {
                opserr << "LadrunoConstraintProjector::buildMass() - massless DOF "
                          "(equation " << e << ", lumped mass " << g.m(r) << ") in a "
                          "constraint group. The projection handler keeps every tied DOF "
                          "in the equation set, so each needs nonzero lumped mass. For "
                          "translational DOFs add physical mass; for ROTATIONAL ties "
                          "(rigidLink -beam, rigidDiaphragm corners tie the perpendicular "
                          "rotation) add a small rotational mass (e.g. ~0.01-0.1% of the "
                          "node's translational mass). Otherwise use constraints "
                          "Penalty/Transformation (which eliminate the slave DOF).\n";
                return -1;
            }
        }

        // LtML = L^T diag(m) L   (nRet x nRet, SPD for positive masses + full-rank L)
        int nRet = g.nRet;
        g.LtML.Zero();
        for (int p = 0; p < nRet; p++) {
            for (int q = 0; q < nRet; q++) {
                double s = 0.0;
                for (int r = 0; r < nrows; r++)
                    s += g.L(r, p) * g.m(r) * g.L(r, q);
                g.LtML(p, q) = s;
            }
        }

        // Ladruno (ADR-30 P6 / O1): condition-number gate. LtML is SPD; estimate
        // cond = lambda_max/lambda_min via a self-contained symmetric Jacobi eigensolve.
        // A near-singular LtML (a barely-dependent retained direction, a tiny relative mass,
        // a near-redundant constraint) passes an exact-pivot test-solve but makes project()
        // amplify round-off into a garbage acceleration -> REFUSE above COND_REFUSE (this also
        // catches the exact-singular case lambda_min -> 0), WARN once above COND_WARN.
        if (nRet > 0) {
            Vector eig;
            if (ladrunoSymEigJacobi(g.LtML, eig) < 0) {
                opserr << "LadrunoConstraintProjector::buildMass() - eigenvalue estimate "
                          "failed for Lt M L in group " << gi << ".\n";
                return -1;
            }
            double lmax = eig(0), lmin = eig(0);
            for (int k = 1; k < eig.Size(); k++) {
                if (eig(k) > lmax) lmax = eig(k);
                if (eig(k) < lmin) lmin = eig(k);
            }
            double cond = (lmin > 0.0) ? (lmax / lmin)
                                       : std::numeric_limits<double>::infinity();
            if (lmin <= 0.0 || cond > LADRUNO_LTML_COND_REFUSE) {
                opserr << "LadrunoConstraintProjector::buildMass() - near-singular Lt M L in "
                          "group " << gi << " (condition number " << cond << ", min eigenvalue "
                       << lmin << "). A retained direction has (almost) no mass, or the "
                          "constraints are redundant/cyclic. Fix the constraint set, add mass, "
                          "or use constraints Penalty/Transformation.\n";
                return -1;
            }
            if (cond > LADRUNO_LTML_COND_WARN)
                opserr << "LadrunoConstraintProjector::buildMass() - NOTE Lt M L in group " << gi
                       << " is ill-conditioned (condition number " << cond << "); the projected "
                          "acceleration may lose precision. Check for near-redundant constraints "
                          "or widely disparate tied masses.\n";
        }
    }

    massReady = true;
    return 0;
}

void
LadrunoConstraintProjector::project(Vector &a)
{
    if (!massReady) {
        opserr << "LadrunoConstraintProjector::project() - mass not built; buildMass() "
                  "must run after the first M-only formTangent and before project().\n";
        return;
    }

    // P3: (re)size + zero the tie-force cache to the current system size.
    if (tieForce.Size() != a.Size())
        tieForce = Vector(a.Size());
    tieForce.Zero();

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        Group &g = groups[gi];

        // orphaned slaves (every retained master SP-fixed) -> a_c = C*0 = 0
        for (int i = 0; i < g.fixedEqn.Size(); i++)
            a(g.fixedEqn(i)) = 0.0;   // tie force on these rows not cached (mass not held)

        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;
        if (nRet == 0 || nrows == 0)
            continue;

        // capture a_raw BEFORE overwriting (needed for the tie force f = M(a_raw-a_proj))
        Vector araw(nrows);
        for (int r = 0; r < nrows; r++)
            araw(r) = a(g.allEqn(r));

        // rhs = L^T (M .* a_raw)
        Vector rhs(nRet);
        rhs.Zero();
        for (int p = 0; p < nRet; p++) {
            double s = 0.0;
            for (int r = 0; r < nrows; r++)
                s += g.L(r, p) * g.m(r) * araw(r);
            rhs(p) = s;
        }

        // ar = (L^T M L)^-1 rhs   (solve on a copy; Solve may modify the matrix)
        Vector ar(nRet);
        Matrix tmp(g.LtML);
        if (tmp.Solve(rhs, ar) < 0) {
            opserr << "LadrunoConstraintProjector::project() - solve failed in group "
                   << gi << "\n";
            continue;
        }

        // a_full = L ar, scatter back (overwrites BOTH retained and constrained),
        // and cache the tie force f = m (a_raw - a_proj) per row.
        for (int r = 0; r < nrows; r++) {
            double aproj = 0.0;
            for (int p = 0; p < nRet; p++)
                aproj += g.L(r, p) * ar(p);
            int e = g.allEqn(r);
            tieForce(e) = g.m(r) * (araw(r) - aproj);
            a(e) = aproj;
        }
    }
}

double
LadrunoConstraintProjector::tieForceAtEqn(int eqn) const
{
    if (eqn < 0 || eqn >= tieForce.Size())
        return 0.0;
    return tieForce(eqn);
}

int
LadrunoConstraintProjector::checkIC(const Vector &u, const Vector &v, double icTol) const
{
    int nViol = 0;

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        const Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;

        // scale tolerance by the group's displacement/velocity magnitude
        double umax = 0.0, vmax = 0.0;
        for (int r = 0; r < nrows; r++) {
            double au = fabs(u(g.allEqn(r)));  if (au > umax) umax = au;
            double av = fabs(v(g.allEqn(r)));  if (av > vmax) vmax = av;
        }
        double tolU = icTol * (1.0 + umax);
        double tolV = icTol * (1.0 + vmax);

        // constrained rows: u_c - C u_r ?= delta ;  v_c - C v_r ?= 0
        for (int j = 0; j < nrows - nRet; j++) {
            int row = nRet + j;
            int ec  = g.allEqn(row);
            double Cur = 0.0, Cvr = 0.0;
            for (int p = 0; p < nRet; p++) {
                Cur += g.L(row, p) * u(g.allEqn(p));
                Cvr += g.L(row, p) * v(g.allEqn(p));
            }
            double ru = u(ec) - Cur - g.delta(j);
            double rv = v(ec) - Cvr;
            if (fabs(ru) > tolU || fabs(rv) > tolV) {
                opserr << "LadrunoConstraintProjector::checkIC() - initial state "
                          "violates a projected constraint at equation " << ec
                       << " (disp residual " << ru << ", vel residual " << rv
                       << "). The constrained DOF and its retained master(s) must "
                          "start on the constraint manifold; fix the ICs or pass "
                          "-projectICs.\n";
                nViol++;
            }
        }

        // orphaned (SP-fixed-master) slaves must start at zero
        for (int i = 0; i < g.fixedEqn.Size(); i++) {
            int ec = g.fixedEqn(i);
            if (fabs(u(ec)) > tolU || fabs(v(ec)) > tolV) {
                opserr << "LadrunoConstraintProjector::checkIC() - constrained DOF "
                          "(equation " << ec << ") tied to an SP-fixed master must "
                          "start at zero (disp " << u(ec) << ", vel " << v(ec)
                       << ").\n";
                nViol++;
            }
        }
    }

    return nViol;
}

void
LadrunoConstraintProjector::snapICs(Vector &u, Vector &v) const
{
    for (int gi = 0; gi < (int)groups.size(); gi++) {
        const Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;

        for (int j = 0; j < nrows - nRet; j++) {
            int row = nRet + j;
            int ec  = g.allEqn(row);
            double Cur = 0.0, Cvr = 0.0;
            for (int p = 0; p < nRet; p++) {
                Cur += g.L(row, p) * u(g.allEqn(p));
                Cvr += g.L(row, p) * v(g.allEqn(p));
            }
            u(ec) = Cur + g.delta(j);
            v(ec) = Cvr;
        }
        for (int i = 0; i < g.fixedEqn.Size(); i++) {
            u(g.fixedEqn(i)) = 0.0;
            v(g.fixedEqn(i)) = 0.0;
        }
    }
}
