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

// ADR-62 — LadrunoTie kinematic mesh-tie generator. See LadrunoTie.h for the design.
// This file holds the geometry generator (LadrunoTie::generate) and the command
// front-end (OPS_LadrunoTie, registered in both the Python and Tcl interpreters).

#include "LadrunoTie.h"

#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <ElementIter.h>
#include <EQ_Constraint.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

// header-only, OpenSees-free geometry reused verbatim from the contact engine
#include <LadrunoContactProjection.h>   // closest-point projection + shape weights (ADR-39/41)
#include <LadrunoContactBucketSort.h>   // bucket-sort broad phase (ADR-39 P2.5)

#include <cmath>
#include <set>
#include <vector>
#include <string.h>

// --------------------------------------------------------------------------- //
// helpers
// --------------------------------------------------------------------------- //

// Padded 3D reference coords of a node (z=0 for a 2D node). Returns false if the
// node is not in the domain.
static bool ltNodeCoords3(Domain *dom, int tag, double x[3])
{
    Node *nd = dom->getNode(tag);
    if (nd == 0) return false;
    const Vector &X = nd->getCrds();
    x[0] = x[1] = x[2] = 0.0;
    for (int d = 0; d < X.Size() && d < 3; d++) x[d] = X(d);
    return true;
}

// --------------------------------------------------------------------------- //
// generator
// --------------------------------------------------------------------------- //

int
LadrunoTie::generate(Domain *dom,
                     const ID &slaves,
                     const ID &mfn, int nps,
                     const ID &dofsIn, double tolFrac)
{
    if (dom == 0) {
        opserr << "WARNING LadrunoTie - domain is not defined\n";
        return -1;
    }
    if (nps != 3 && nps != 4) {
        opserr << "WARNING LadrunoTie - nps must be 3 (tri-3) or 4 (quad-4), got " << nps << "\n";
        return -1;
    }
    if (slaves.Size() == 0) {
        opserr << "WARNING LadrunoTie - no slave nodes given\n";
        return -1;
    }
    if (mfn.Size() == 0 || (mfn.Size() % nps) != 0) {
        opserr << "WARNING LadrunoTie - master facet node count " << mfn.Size()
               << " is not a positive multiple of nps " << nps << "\n";
        return -1;
    }
    const int nf = mfn.Size() / nps;

    // tied DOFs: default to the spatial dimension of the first slave node.
    ID dofs = dofsIn;
    if (dofs.Size() == 0) {
        Node *s0 = dom->getNode(slaves(0));
        if (s0 == 0) {
            opserr << "WARNING LadrunoTie - slave node " << slaves(0) << " not in domain\n";
            return -1;
        }
        int ndm = s0->getCrds().Size();
        dofs.resize(ndm);
        for (int d = 0; d < ndm; d++) dofs(d) = d + 1;   // 1-based
    }

    // --- master facets: flat reference coords (nf*nps*3), node-disjoint set, sizes ---
    std::vector<double> seg((size_t)nf * nps * 3, 0.0);
    std::set<int> masterSet;
    std::vector<double> facetSize(nf, 0.0);
    for (int f = 0; f < nf; f++) {
        double Xf[4][3];
        for (int i = 0; i < nps; i++) {
            int tag = mfn(f * nps + i);
            double x[3];
            if (!ltNodeCoords3(dom, tag, x)) {
                opserr << "WARNING LadrunoTie - master facet node " << tag << " not in domain\n";
                return -1;
            }
            masterSet.insert(tag);
            for (int d = 0; d < 3; d++) { seg[((size_t)f * nps + i) * 3 + d] = x[d]; Xf[i][d] = x[d]; }
        }
        // characteristic facet size = max pairwise node distance (the diagonal)
        double dmax = 0.0;
        for (int i = 0; i < nps; i++)
            for (int j = i + 1; j < nps; j++) {
                double dx = Xf[i][0]-Xf[j][0], dy = Xf[i][1]-Xf[j][1], dz = Xf[i][2]-Xf[j][2];
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d > dmax) dmax = d;
            }
        facetSize[f] = dmax;
    }

    // --- broad phase over the master facets (reference config; a tie is static) ---
    LadrunoContactBucketSort::Grid grid(nf, nps, seg.data());
    std::vector<int> cand(nf);

    // --- BLOCKER-2 pre-scan: nodes carrying lumped mass (nodal mass OR an incident
    //     element whose mass matrix has a nonzero diagonal). At model-build the
    //     element mass is already formable from rho+geometry (same call the
    //     projection handler's consistentMassGuard makes at handle() time). ---
    std::set<int> massed;
    ElementIter &eIter = dom->getElements();
    Element *elePtr;
    while ((elePtr = eIter()) != 0) {
        const Matrix &M = elePtr->getMass();
        int n = M.noRows();
        if (n == 0) continue;
        double mx = 0.0;
        for (int p = 0; p < n; p++) if (std::fabs(M(p,p)) > mx) mx = std::fabs(M(p,p));
        if (mx <= 0.0) continue;
        const ID &en = elePtr->getExternalNodes();
        for (int k = 0; k < en.Size(); k++) massed.insert(en(k));
    }

    std::set<int> doneSlaves;
    int emitted = 0;

    for (int si = 0; si < slaves.Size(); si++) {
        int s = slaves(si);

        // BLOCKER-1a: a slave listed twice would be doubly constrained.
        if (doneSlaves.count(s)) {
            opserr << "WARNING LadrunoTie - slave node " << s << " is listed more than once. v1 ties "
                      "each slave node to exactly one facet (no double constraints). Remove the duplicate.\n";
            return -1;
        }
        doneSlaves.insert(s);

        // BLOCKER-1b: slave and master surfaces must be node-disjoint (else MP-chain).
        if (masterSet.count(s)) {
            opserr << "WARNING LadrunoTie - node " << s << " is BOTH a slave and a master facet node. "
                      "The projection handler refuses MP-chains; the slave and master surfaces must be "
                      "node-disjoint.\n";
            return -1;
        }

        Node *sNode = dom->getNode(s);
        double xs[3];
        if (sNode == 0 || !ltNodeCoords3(dom, s, xs)) {
            opserr << "WARNING LadrunoTie - slave node " << s << " not in domain\n";
            return -1;
        }

        // BLOCKER-2: every tied slave DOF must carry lumped mass (the projection keeps
        // slave DOFs in the equation set, so a massless tied DOF is singular).
        bool hasMass = massed.count(s) != 0;
        if (!hasMass) {
            const Matrix &Mn = sNode->getMass();
            for (int p = 0; p < Mn.noRows(); p++) if (std::fabs(Mn(p,p)) > 0.0) { hasMass = true; break; }
        }
        if (!hasMass) {
            opserr << "WARNING LadrunoTie - tied slave node " << s << " carries no mass (no nodal mass and "
                      "no incident element with mass). Define LadrunoTie AFTER the elements/masses; the "
                      "projection handler requires lumped mass on every tied DOF.\n";
            return -1;
        }

        // --- nearest in-bounds facet (broad phase first, brute-force fallback) ---
        int    bestF = -1;
        double bestDist = 0.0, bestXi = 0.0, bestEta = 0.0;
        for (int pass = 0; pass < 2 && bestF < 0; pass++) {
            int nc = (pass == 0) ? grid.candidates(xs, cand.data()) : nf;
            for (int ci = 0; ci < nc; ci++) {
                int f = (pass == 0) ? cand[ci] : ci;
                double X[4][3];
                for (int i = 0; i < nps; i++)
                    for (int d = 0; d < 3; d++) X[i][d] = seg[((size_t)f * nps + i) * 3 + d];
                double xi, eta;
                if (LadrunoContactProjection::project(nps, X, xs, xi, eta, 1e-12, 30) != 0)
                    continue;   // out-of-bounds or degenerate facet
                double N[4], dNx[4], dNe[4], xb[3];
                LadrunoContactProjection::shape(nps, xi, eta, N, dNx, dNe);
                LadrunoContactProjection::interp(nps, N, X, xb);
                double dd = std::sqrt((xs[0]-xb[0])*(xs[0]-xb[0]) +
                                      (xs[1]-xb[1])*(xs[1]-xb[1]) +
                                      (xs[2]-xb[2])*(xs[2]-xb[2]));
                if (bestF < 0 || dd < bestDist) { bestF = f; bestDist = dd; bestXi = xi; bestEta = eta; }
            }
        }
        if (bestF < 0) {
            opserr << "WARNING LadrunoTie - slave node " << s << " does not project in-bounds onto any "
                      "master facet. Check the master surface covers the slave surface (nps/facets correct).\n";
            return -1;
        }

        // OQ-3: require conforming-at-interface geometry (no IC snapping).
        double tol = tolFrac * (facetSize[bestF] > 0.0 ? facetSize[bestF] : 1.0);
        if (bestDist > tol) {
            opserr << "WARNING LadrunoTie - slave node " << s << " projects " << bestDist
                   << " off the master surface (tol " << tol << " = " << tolFrac
                   << " x facet size). v1 requires conforming-at-interface geometry: the slave must lie on "
                      "the master surface (off-manifold ICs are not snapped). Move the node onto the "
                      "interface, or relax -tol if the offset is intentional.\n";
            return -1;
        }

        // --- shape weights at the projection; emit one EQ_Constraint per tied DOF ---
        double X[4][3];
        for (int i = 0; i < nps; i++)
            for (int d = 0; d < 3; d++) X[i][d] = seg[((size_t)bestF * nps + i) * 3 + d];
        double N[4], dNx[4], dNe[4];
        LadrunoContactProjection::shape(nps, bestXi, bestEta, N, dNx, dNe);

        // drop near-zero weights so a corner/edge collocation does not spuriously couple
        // (and group-connect) master nodes that carry no share of this slave.
        int cnt = 0;
        for (int i = 0; i < nps; i++) if (std::fabs(N[i]) > 1e-12) cnt++;

        for (int di = 0; di < dofs.Size(); di++) {
            int dof1 = dofs(di);                 // 1-based
            if (dof1 < 1 || dof1 > sNode->getNumberDOF()) {
                opserr << "WARNING LadrunoTie - slave node " << s << " has no DOF " << dof1
                       << " (numberDOF " << sNode->getNumberDOF() << ")\n";
                return -1;
            }
            ID rNode(cnt), rDOF(cnt);
            Vector Ccr(cnt);
            int p = 0;
            for (int i = 0; i < nps; i++)
                if (std::fabs(N[i]) > 1e-12) {
                    rNode(p) = mfn(bestF * nps + i);
                    rDOF(p)  = dof1 - 1;          // 0-based, stored
                    Ccr(p)   = N[i];              // u_s = sum N_i u_{m,i}  (Ccr = N, no negation)
                    p++;
                }
            // ModelBuilder EQ_Constraint ctor stores `constraint` = Ccr directly.
            EQ_Constraint *eq = new EQ_Constraint(s, dof1 - 1, Ccr, rNode, rDOF);
            if (eq == 0 || dom->addEQ_Constraint(eq) == false) {
                opserr << "WARNING LadrunoTie - failed to add EQ_Constraint for slave node " << s
                       << " dof " << dof1 << "\n";
                if (eq != 0) delete eq;
                return -1;
            }
            emitted++;
        }
    }

    opserr << "LadrunoTie: emitted " << emitted << " EQ_Constraint(s) tying " << (int)doneSlaves.size()
           << " slave node(s) to " << nf << " master facet(s). Enforce with `constraints LadrunoProjection`.\n";
    return emitted;
}

// --------------------------------------------------------------------------- //
// command front-end
//
//   LadrunoTie -slaveNodes <ns> s1 ... sns
//              -masterFacets <nps> <nf> (nf groups of nps node tags)
//              [-dof <nd> d1 ... dnd]      (default: 1..ndm of the slave nodes)
//              [-tol <frac>]              (default 1e-6 * facet size)
//
// Flat token streams (matching equationConstraint/contactSurface); both Tcl and
// Python pass node tags as separate arguments, not a single braced/list arg.
// --------------------------------------------------------------------------- //

int OPS_LadrunoTie()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) {
        opserr << "WARNING LadrunoTie - domain is not defined\n";
        return -1;
    }
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING want - LadrunoTie -slaveNodes ns s1.. -masterFacets nps nf m.. "
                  "<-dof nd d1..> <-tol frac>\n";
        return -1;
    }

    ID slaveNodes(0), masterFacetNodes(0), dofs(0);
    int nps = 0;
    double tolFrac = 1.0e-6;
    bool haveSlaves = false, haveMaster = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-slaveNodes") == 0 || strcmp(opt, "-slave") == 0) {
            int one = 1, ns = 0;
            if (OPS_GetIntInput(&one, &ns) < 0 || ns < 1) {
                opserr << "WARNING LadrunoTie -slaveNodes - need a positive count\n";
                return -1;
            }
            slaveNodes.resize(ns);
            for (int i = 0; i < ns; i++) {
                int v;
                if (OPS_GetIntInput(&one, &v) < 0) {
                    opserr << "WARNING LadrunoTie -slaveNodes - could not read node tag " << i + 1 << "\n";
                    return -1;
                }
                slaveNodes(i) = v;
            }
            haveSlaves = true;
        }
        else if (strcmp(opt, "-masterFacets") == 0 || strcmp(opt, "-master") == 0) {
            int one = 1, nf = 0;
            if (OPS_GetIntInput(&one, &nps) < 0 || (nps != 3 && nps != 4)) {
                opserr << "WARNING LadrunoTie -masterFacets - need nps (3 or 4)\n";
                return -1;
            }
            if (OPS_GetIntInput(&one, &nf) < 0 || nf < 1) {
                opserr << "WARNING LadrunoTie -masterFacets - need a positive facet count\n";
                return -1;
            }
            masterFacetNodes.resize(nf * nps);
            for (int i = 0; i < nf * nps; i++) {
                int v;
                if (OPS_GetIntInput(&one, &v) < 0) {
                    opserr << "WARNING LadrunoTie -masterFacets - could not read facet node " << i + 1
                           << " of " << nf * nps << "\n";
                    return -1;
                }
                masterFacetNodes(i) = v;
            }
            haveMaster = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            int one = 1, nd = 0;
            if (OPS_GetIntInput(&one, &nd) < 0 || nd < 1) {
                opserr << "WARNING LadrunoTie -dof - need a positive count\n";
                return -1;
            }
            dofs.resize(nd);
            for (int i = 0; i < nd; i++) {
                int v;
                if (OPS_GetIntInput(&one, &v) < 0 || v < 1) {
                    opserr << "WARNING LadrunoTie -dof - could not read DOF " << i + 1 << "\n";
                    return -1;
                }
                dofs(i) = v;
            }
        }
        else if (strcmp(opt, "-tol") == 0) {
            int one = 1;
            if (OPS_GetDoubleInput(&one, &tolFrac) < 0 || tolFrac < 0.0) {
                opserr << "WARNING LadrunoTie -tol - need a non-negative number\n";
                return -1;
            }
        }
        else {
            opserr << "WARNING LadrunoTie - unknown option " << opt << "\n";
            return -1;
        }
    }

    if (!haveSlaves || !haveMaster) {
        opserr << "WARNING LadrunoTie - both -slaveNodes and -masterFacets are required\n";
        return -1;
    }

    int rc = LadrunoTie::generate(theDomain, slaveNodes, masterFacetNodes, nps, dofs, tolFrac);
    return (rc < 0) ? -1 : 0;
}
