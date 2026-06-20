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

        // singularity guard: a rank-deficient LtML = a massless/redundant retained
        // direction. Probe with a test solve on a copy (Solve may factor in place).
        if (nRet > 0) {
            Matrix tmp(g.LtML);
            Vector e0(nRet), x0(nRet);
            e0(0) = 1.0;
            if (tmp.Solve(e0, x0) < 0) {
                opserr << "LadrunoConstraintProjector::buildMass() - singular Lt M L in "
                          "group " << gi << " (a retained direction has no mass, or a "
                          "redundant/cyclic constraint). Fix the constraint set or add "
                          "mass.\n";
                return -1;
            }
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
