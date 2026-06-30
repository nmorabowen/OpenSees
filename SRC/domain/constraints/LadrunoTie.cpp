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
#include <LadrunoMortarKernel.h>        // clip + Gauss -> mortar D, M, weighted gap (ADR-41 C1)

#include <cmath>
#include <map>
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
// P2 — integral-mortar generator
// --------------------------------------------------------------------------- //

int
LadrunoTie::generateMortar(Domain *dom,
                           const ID &sfn, int npsS,
                           const ID &mfn, int npsM,
                           const ID &dofsIn, double tolFrac, const double *outward)
{
    if (dom == 0) { opserr << "WARNING LadrunoTie - domain is not defined\n"; return -1; }
    if (npsS != 3 && npsS != 4) {
        opserr << "WARNING LadrunoTie -mortar - slave nps must be 3 or 4, got " << npsS << "\n"; return -1; }
    if (npsM != 3 && npsM != 4) {
        opserr << "WARNING LadrunoTie -mortar - master nps must be 3 or 4, got " << npsM << "\n"; return -1; }
    if (sfn.Size() == 0 || (sfn.Size() % npsS) != 0) {
        opserr << "WARNING LadrunoTie -mortar - slave facet node count " << sfn.Size()
               << " not a positive multiple of npsS " << npsS << "\n"; return -1; }
    if (mfn.Size() == 0 || (mfn.Size() % npsM) != 0) {
        opserr << "WARNING LadrunoTie -mortar - master facet node count " << mfn.Size()
               << " not a positive multiple of npsM " << npsM << "\n"; return -1; }
    const int nfS = sfn.Size() / npsS;
    const int nfM = mfn.Size() / npsM;

    // --- unique slave / master node lists (tag <-> local index) ---
    std::map<int,int> sIdx, mIdx;
    std::vector<int> sTag, mTag;
    for (int i = 0; i < sfn.Size(); i++)
        if (sIdx.find(sfn(i)) == sIdx.end()) { sIdx[sfn(i)] = (int)sTag.size(); sTag.push_back(sfn(i)); }
    for (int i = 0; i < mfn.Size(); i++)
        if (mIdx.find(mfn(i)) == mIdx.end()) { mIdx[mfn(i)] = (int)mTag.size(); mTag.push_back(mfn(i)); }
    const int Ns = (int)sTag.size(), Nm = (int)mTag.size();

    // BLOCKER-1: node-disjoint slave / master surfaces (else MP-chain).
    for (int i = 0; i < Ns; i++)
        if (mIdx.find(sTag[i]) != mIdx.end()) {
            opserr << "WARNING LadrunoTie -mortar - node " << sTag[i] << " is on BOTH the slave and "
                      "master surface. They must be node-disjoint (the handler refuses MP-chains).\n";
            return -1;
        }

    // tied DOFs default = spatial dim of the first slave node.
    ID dofs = dofsIn;
    if (dofs.Size() == 0) {
        Node *s0 = dom->getNode(sTag[0]);
        if (s0 == 0) { opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[0] << " not in domain\n"; return -1; }
        int ndm = s0->getCrds().Size();
        dofs.resize(ndm);
        for (int d = 0; d < ndm; d++) dofs(d) = d + 1;
    }

    // --- master facet reference coords (flat nfM*npsM*3) + average outward normal ---
    std::vector<double> segM((size_t)nfM * npsM * 3, 0.0);
    double refDir[3] = {0.0, 0.0, 0.0};
    for (int f = 0; f < nfM; f++) {
        double X[4][3];
        for (int i = 0; i < npsM; i++) {
            double x[3];
            if (!ltNodeCoords3(dom, mfn(f * npsM + i), x)) {
                opserr << "WARNING LadrunoTie -mortar - master facet node " << mfn(f*npsM+i) << " not in domain\n";
                return -1;
            }
            for (int d = 0; d < 3; d++) { segM[((size_t)f*npsM+i)*3+d] = x[d]; X[i][d] = x[d]; }
        }
        // accumulate the (unoriented) facet normal for a default refDir
        double a1[3], a2[3], nrm[3];
        for (int d = 0; d < 3; d++) { a1[d] = X[1][d]-X[0][d]; a2[d] = X[npsM-1][d]-X[0][d]; }
        nrm[0]=a1[1]*a2[2]-a1[2]*a2[1]; nrm[1]=a1[2]*a2[0]-a1[0]*a2[2]; nrm[2]=a1[0]*a2[1]-a1[1]*a2[0];
        double L = std::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
        if (L > 1e-300) for (int d = 0; d < 3; d++) refDir[d] += nrm[d]/L;
    }
    if (outward != 0) { for (int d = 0; d < 3; d++) refDir[d] = outward[d]; }
    double Lr = std::sqrt(refDir[0]*refDir[0]+refDir[1]*refDir[1]+refDir[2]*refDir[2]);
    if (Lr < 1e-12) {
        opserr << "WARNING LadrunoTie -mortar - could not infer an outward direction (master normals "
                  "cancel; the surface may be closed/curved). Pass -outward ox oy oz.\n";
        return -1;
    }
    for (int d = 0; d < 3; d++) refDir[d] /= Lr;

    // --- assemble global D (Ns x Ns), M (Ns x Nm), coverage a_I, weighted gap g_I ---
    // Brute-force slave x master facet pairs: the clip returns instantly for a non-
    // overlapping pair, and an exhaustive sweep guarantees the mortar integral captures
    // EVERY overlap (no silent coverage gap). The ADR-39 bucket-sort is the large-
    // interface optimization (the coverage guard below would flag any missed overlap).
    Matrix D(Ns, Ns), M(Ns, Nm);
    D.Zero(); M.Zero();
    // cover[I]  = INT N_I over the slave/MASTER overlap (the bonded measure);
    // fullCov[I]= INT N_I over the whole slave surface (master-independent, via a
    //             self-clip of each slave facet) = the tributary area node I *should*
    //             have if fully covered. cover/fullCov < 1 => the slave protrudes past
    //             the master there (a partially-bonded node) => refuse.
    std::vector<double> cover(Ns, 0.0), fullCov(Ns, 0.0), gap(Ns, 0.0);
    for (int fs = 0; fs < nfS; fs++) {
        double Xs[4][3];
        for (int i = 0; i < npsS; i++) {
            double x[3];
            if (!ltNodeCoords3(dom, sfn(fs*npsS+i), x)) {
                opserr << "WARNING LadrunoTie -mortar - slave facet node " << sfn(fs*npsS+i) << " not in domain\n";
                return -1;
            }
            for (int d = 0; d < 3; d++) Xs[i][d] = x[d];
        }
        // slave facet self-mass (clip the facet against itself = the full facet) ->
        // fullCov[I] += Sum_J INT N_I N_J = INT N_I over this facet.
        LadrunoMortarKernel::PairResult prSelf;
        if (LadrunoMortarKernel::integratePair(npsS, Xs, npsS, Xs, refDir, prSelf, 4) == 0)
            for (int a = 0; a < npsS; a++) {
                int I = sIdx[sfn(fs*npsS+a)];
                for (int b = 0; b < npsS; b++) fullCov[I] += prSelf.D[a][b];
            }
        for (int fm = 0; fm < nfM; fm++) {
            double Xm[4][3];
            for (int i = 0; i < npsM; i++)
                for (int d = 0; d < 3; d++) Xm[i][d] = segM[((size_t)fm*npsM+i)*3+d];
            LadrunoMortarKernel::PairResult pr;
            if (LadrunoMortarKernel::integratePair(npsS, Xs, npsM, Xm, refDir, pr, 4) != 0)
                continue;                              // empty/degenerate overlap
            for (int a = 0; a < npsS; a++) {
                int I = sIdx[sfn(fs*npsS+a)];
                gap[I] += pr.g[a];
                for (int b = 0; b < npsS; b++) {
                    int J = sIdx[sfn(fs*npsS+b)];
                    D(I,J) += pr.D[a][b];
                    cover[I] += pr.D[a][b];            // Sum_J D_IJ = INT N_I dGamma
                }
                for (int k = 0; k < npsM; k++)
                    M(I, mIdx[mfn(fm*npsM+k)]) += pr.M[a][k];
            }
        }
    }

    // --- coverage + conforming-gap guards (mortar OQ-3 / uncovered-node refusal) ---
    double maxCover = 0.0;
    for (int I = 0; I < Ns; I++) if (cover[I] > maxCover) maxCover = cover[I];
    if (maxCover <= 1e-300) {
        opserr << "WARNING LadrunoTie -mortar - the slave and master surfaces do not overlap "
                  "(zero integrated area). Check the facets are coincident and the nps/winding correct.\n";
        return -1;
    }
    double facetScale = std::sqrt(maxCover);          // ~ a facet edge length
    for (int I = 0; I < Ns; I++) {
        // a slave node must be covered over essentially its WHOLE tributary area, else
        // the slave surface protrudes beyond the master there (a partial / extrapolated
        // bond) and/or D's row is near-singular (DGESV only flags an EXACT zero pivot).
        // Require cover/fullCov >= 1 - 1e-3 (0.1% slack for clip/quadrature noise).
        if (fullCov[I] > 1e-300 && cover[I] < (1.0 - 1.0e-3) * fullCov[I]) {
            opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[I] << " is only "
                   << (cover[I] / fullCov[I]) * 100.0 << "% covered by the master surface (the slave "
                      "surface protrudes beyond the master). A mesh-tie needs the master surface to cover "
                      "the whole slave surface. Extend the master surface or trim the slave.\n";
            return -1;
        }
        if (cover[I] < 1e-6 * maxCover) {             // backstop: a (near-)zero D row
            opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[I] << " is not covered by the "
                      "master surface (integrated measure ~0 => D singular). Extend the master surface.\n";
            return -1;
        }
        double gbar = std::fabs(gap[I]) / cover[I];
        if (gbar > tolFrac * facetScale) {
            opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[I] << " sits " << gbar
                   << " off the master surface (tol " << tolFrac * facetScale << "). v1 requires "
                      "conforming-at-interface geometry (coincident surfaces). Move the slave onto the "
                      "interface, or relax -tol.\n";
            return -1;
        }
    }

    // --- BLOCKER-2: every tied slave node carries lumped mass ---
    std::set<int> massed;
    ElementIter &eIter = dom->getElements();
    Element *elePtr;
    while ((elePtr = eIter()) != 0) {
        const Matrix &Me = elePtr->getMass();
        int n = Me.noRows();
        if (n == 0) continue;
        double mx = 0.0;
        for (int p = 0; p < n; p++) if (std::fabs(Me(p,p)) > mx) mx = std::fabs(Me(p,p));
        if (mx <= 0.0) continue;
        const ID &en = elePtr->getExternalNodes();
        for (int k = 0; k < en.Size(); k++) massed.insert(en(k));
    }
    for (int I = 0; I < Ns; I++) {
        bool hasMass = massed.count(sTag[I]) != 0;
        if (!hasMass) {
            Node *nd = dom->getNode(sTag[I]);
            const Matrix &Mn = nd->getMass();
            for (int p = 0; p < Mn.noRows(); p++) if (std::fabs(Mn(p,p)) > 0.0) { hasMass = true; break; }
        }
        if (!hasMass) {
            opserr << "WARNING LadrunoTie -mortar - tied slave node " << sTag[I] << " carries no mass "
                      "(no nodal mass, no incident element with mass). Define LadrunoTie AFTER the "
                      "elements/masses; the projection handler requires lumped mass on tied DOFs.\n";
            return -1;
        }
    }

    // --- condense P = D^{-1} M (one SPD solve at setup) ---
    Matrix P(Ns, Nm);
    if (D.Solve(M, P) < 0) {
        opserr << "WARNING LadrunoTie -mortar - failed to solve D P = M (the slave-interface mass D is "
                  "singular; an uncovered slave node likely slipped the coverage guard).\n";
        return -1;
    }
    // Post-solve sanity: every transfer row must satisfy partition of unity (Sum_k P_Ik = 1),
    // which is algebraically exact for P = D^{-1}M (D 1 = M 1). A row that drifts from 1 means
    // the solve was ill-conditioned (a near-singular D that DGESV accepted) => refuse rather
    // than emit a constant-field-violating tie.
    for (int I = 0; I < Ns; I++) {
        double rs = 0.0;
        for (int k = 0; k < Nm; k++) rs += P(I, k);
        if (std::fabs(rs - 1.0) > 1e-6) {
            opserr << "WARNING LadrunoTie -mortar - the mortar transfer for slave node " << sTag[I]
                   << " violates partition of unity (Sum P = " << rs << " != 1): the interface matrix D "
                      "is ill-conditioned. Check for sliver overlaps / near-degenerate slave facets.\n";
            return -1;
        }
    }

    // --- emit one EQ_Constraint per slave node per tied DOF: u_s = sum_k P_sk u_{m,k} ---
    int emitted = 0;
    for (int I = 0; I < Ns; I++) {
        int s = sTag[I];
        Node *sNode = dom->getNode(s);
        if (sNode == 0) { opserr << "WARNING LadrunoTie -mortar - slave node " << s << " not in domain\n"; return -1; }

        // master columns with a non-negligible weight (drop true zeros only; P is dense)
        int cnt = 0;
        for (int k = 0; k < Nm; k++) if (std::fabs(P(I,k)) > 1e-12) cnt++;
        if (cnt == 0) {
            opserr << "WARNING LadrunoTie -mortar - slave node " << s << " produced an empty transfer row.\n";
            return -1;
        }

        for (int di = 0; di < dofs.Size(); di++) {
            int dof1 = dofs(di);
            if (dof1 < 1 || dof1 > sNode->getNumberDOF()) {
                opserr << "WARNING LadrunoTie -mortar - slave node " << s << " has no DOF " << dof1 << "\n";
                return -1;
            }
            ID rNode(cnt), rDOF(cnt);
            Vector Ccr(cnt);
            int p = 0;
            for (int k = 0; k < Nm; k++)
                if (std::fabs(P(I,k)) > 1e-12) {
                    rNode(p) = mTag[k];
                    rDOF(p)  = dof1 - 1;
                    Ccr(p)   = P(I,k);                 // u_s = sum_k P_sk u_{m,k}  (Ccr = P row)
                    p++;
                }
            EQ_Constraint *eq = new EQ_Constraint(s, dof1 - 1, Ccr, rNode, rDOF);
            if (eq == 0 || dom->addEQ_Constraint(eq) == false) {
                opserr << "WARNING LadrunoTie -mortar - failed to add EQ_Constraint for slave node " << s
                       << " dof " << dof1 << "\n";
                if (eq != 0) delete eq;
                return -1;
            }
            emitted++;
        }
    }

    opserr << "LadrunoTie (mortar): emitted " << emitted << " EQ_Constraint(s) tying " << Ns
           << " slave node(s) to " << Nm << " master node(s) over " << nfS << " slave facet(s). "
              "Enforce with `constraints LadrunoProjection`.\n";
    return emitted;
}

// --------------------------------------------------------------------------- //
// command front-end
//
//  P1 (collocation, default):
//   LadrunoTie -slaveNodes <ns> s1.. -masterFacets <npsM> <nfM> m..
//              [-dof <nd> d1..] [-tol <frac>]
//
//  P2 (integral mortar):
//   LadrunoTie -mortar -slaveFacets <npsS> <nfS> s.. -masterFacets <npsM> <nfM> m..
//              [-dof <nd> d1..] [-tol <frac>] [-outward ox oy oz]
//
// Flat token streams (matching equationConstraint/contactSurface); both Tcl and
// Python pass node tags as separate arguments, not a single braced/list arg.
// --------------------------------------------------------------------------- //

// read "<n> v1..vn" into id; returns false on a parse/count error.
static bool ltReadCountedID(ID &id, const char *what)
{
    int one = 1, n = 0;
    if (OPS_GetIntInput(&one, &n) < 0 || n < 1) {
        opserr << "WARNING LadrunoTie " << what << " - need a positive count\n";
        return false;
    }
    id.resize(n);
    for (int i = 0; i < n; i++) {
        int v;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING LadrunoTie " << what << " - could not read entry " << i + 1 << "\n";
            return false;
        }
        id(i) = v;
    }
    return true;
}

// read "<nps> <nf> (nps*nf tags)" into id; sets nps; returns false on error.
static bool ltReadFacets(ID &id, int &nps, const char *what)
{
    int one = 1, nf = 0;
    if (OPS_GetIntInput(&one, &nps) < 0 || (nps != 3 && nps != 4)) {
        opserr << "WARNING LadrunoTie " << what << " - need nps (3 or 4)\n";
        return false;
    }
    if (OPS_GetIntInput(&one, &nf) < 0 || nf < 1) {
        opserr << "WARNING LadrunoTie " << what << " - need a positive facet count\n";
        return false;
    }
    id.resize(nf * nps);
    for (int i = 0; i < nf * nps; i++) {
        int v;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING LadrunoTie " << what << " - could not read facet node " << i + 1
                   << " of " << nf * nps << "\n";
            return false;
        }
        id(i) = v;
    }
    return true;
}

int OPS_LadrunoTie()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) {
        opserr << "WARNING LadrunoTie - domain is not defined\n";
        return -1;
    }
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING want - LadrunoTie -slaveNodes ns s1.. -masterFacets npsM nf m.. "
                  "<-dof nd d1..> <-tol frac>   |   LadrunoTie -mortar -slaveFacets npsS nf s.. "
                  "-masterFacets npsM nf m.. <-dof..> <-tol..> <-outward ox oy oz>\n";
        return -1;
    }

    ID slaveNodes(0), slaveFacetNodes(0), masterFacetNodes(0), dofs(0);
    int npsM = 0, npsS = 0;
    double tolFrac = 1.0e-6;
    double outward[3] = {0.0, 0.0, 0.0};
    bool mortar = false, haveSlaveNodes = false, haveSlaveFacets = false;
    bool haveMaster = false, haveOutward = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-mortar") == 0) {
            mortar = true;
        }
        else if (strcmp(opt, "-slaveNodes") == 0 || strcmp(opt, "-slave") == 0) {
            if (!ltReadCountedID(slaveNodes, "-slaveNodes")) return -1;
            haveSlaveNodes = true;
        }
        else if (strcmp(opt, "-slaveFacets") == 0) {
            if (!ltReadFacets(slaveFacetNodes, npsS, "-slaveFacets")) return -1;
            haveSlaveFacets = true;
        }
        else if (strcmp(opt, "-masterFacets") == 0 || strcmp(opt, "-master") == 0) {
            if (!ltReadFacets(masterFacetNodes, npsM, "-masterFacets")) return -1;
            haveMaster = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            if (!ltReadCountedID(dofs, "-dof")) return -1;
        }
        else if (strcmp(opt, "-tol") == 0) {
            int one = 1;
            if (OPS_GetDoubleInput(&one, &tolFrac) < 0 || tolFrac < 0.0) {
                opserr << "WARNING LadrunoTie -tol - need a non-negative number\n";
                return -1;
            }
        }
        else if (strcmp(opt, "-outward") == 0) {
            int three = 3;
            if (OPS_GetDoubleInput(&three, outward) < 0) {
                opserr << "WARNING LadrunoTie -outward - need ox oy oz\n";
                return -1;
            }
            haveOutward = true;
        }
        else {
            opserr << "WARNING LadrunoTie - unknown option " << opt << "\n";
            return -1;
        }
    }

    int rc;
    if (mortar) {
        if (!haveSlaveFacets || !haveMaster) {
            opserr << "WARNING LadrunoTie -mortar - both -slaveFacets and -masterFacets are required\n";
            return -1;
        }
        rc = LadrunoTie::generateMortar(theDomain, slaveFacetNodes, npsS, masterFacetNodes, npsM,
                                        dofs, tolFrac, haveOutward ? outward : 0);
    } else {
        if (!haveSlaveNodes || !haveMaster) {
            opserr << "WARNING LadrunoTie - both -slaveNodes and -masterFacets are required "
                      "(did you mean -mortar with -slaveFacets?)\n";
            return -1;
        }
        rc = LadrunoTie::generate(theDomain, slaveNodes, masterFacetNodes, npsM, dofs, tolFrac);
    }
    return (rc < 0) ? -1 : 0;
}
