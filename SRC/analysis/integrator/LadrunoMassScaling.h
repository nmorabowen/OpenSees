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

// Ladruno: selective (conventional / DT2MS-style) mass-scaling machinery for the
// explicit integrators. See Ladruno_implementation/36_ladruno_selective_mass_scaling_adr.md.
//
// For every element whose per-element stable step dt_e < dt_target, scale the
// element's mass by s_e = (dt_target/dt_e)^2 so dt_e -> dt_target, implemented by
// ADDING the fictitious increment (s_e-1)*m to the element's NODES (additive
// nodal mass, M2). The per-element dt_e comes from the shared CriticalTimeStep
// kernel (lumpElementMass + elementCriticalDt -- self-report aware, ADR-36 SMS-1),
// so SMS lumps/eigensolves EXACTLY as the dt_cr query does.
//
// Nodal injection (not SOE injection): nodal mass is the canonical M every
// consumer reads (the leap-frog M^-1 AND the EnergyBalanceRecorder KE), so energy
// still closes with the scaled mass. The injected increment is recorded per node
// so it can be re-baselined on a re-entrant domainChanged and restored on
// integrator destruction (M5) -- never a permanent mutation of shared Domain state.
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoMassScaling_h
#define LadrunoMassScaling_h

#include <CriticalTimeStep.h>   // CTSLumping, lumpElementMass, elementCriticalDt
#include <AnalysisModel.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <OPS_Globals.h>

#include <map>
#include <cmath>

namespace Ladruno {

struct MassScalingReport {
    double addedMass;     // sum of injected translational increments
    double modelMass;     // total ELEMENT translational lumped mass (the denominator;
                          //   nodal getMass() is zero on element-rho models, so the
                          //   cap MUST use the element mass — SMS-CAP-DEAD fix)
    int    nScaled;       // elements that received scaling (dt_e < dt_target)
    int    nElems;        // elements scanned
    double minDtScaled;   // smallest un-scaled dt_e that was below target (diagnostic)
    int    nSelfReport;   // throttling elements SKIPPED because their bound is mass-
                          //   independent (self-reported) — scaling can't help them
    int    nMismatch;     // elements skipped due to a non-node-major / DOF mismatch
};

// Build the per-node fictitious-mass increment (additive diagonal) into `injected`
// (node tag -> Vector(ndf)) for the given target step. Does NOT touch the Domain;
// call applyMassScaling to commit. `injected` must be empty (freshly re-baselined).
inline MassScalingReport
buildMassScaling(AnalysisModel *theModel, double dtTarget, CTSLumping lumping,
                 bool useTangent, std::map<int, Vector> &injected)
{
    MassScalingReport rep; rep.addedMass = 0.0; rep.modelMass = 0.0;
    rep.nScaled = 0; rep.nElems = 0; rep.minDtScaled = dtTarget;
    rep.nSelfReport = 0; rep.nMismatch = 0;
    if (theModel == 0 || dtTarget <= 0.0) return rep;
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) return rep;

    Element *ele;
    ElementIter &elements = theDomain->getElements();
    while ((ele = elements()) != 0) {
        rep.nElems++;

        // mass first (D8), lumped via the SAME shared kernel the dt_cr query uses
        const Matrix &M = ele->getMass();
        int n = M.noRows();
        if (n == 0) continue;
        double *mdiag = new double[n];
        lumpElementMass(ele, M, lumping, mdiag);

        // --- model-mass denominator (SMS-CAP-DEAD): sum each element's translational
        //     lumped mass. Element masses are distinct physical contributions, so the
        //     sum over elements is the structural translational mass (no double count).
        //     This is the cap denominator; nodal getMass() is zero on element models.
        {
            Node **nds = ele->getNodePtrs();
            int nn = ele->getNumExternalNodes();
            int p = 0;
            for (int a = 0; a < nn && p < n; ++a) {
                Node *ndp = nds ? nds[a] : 0;
                if (ndp == 0) break;
                int ndmA = ndp->getCrds().Size();
                int ndf = ndp->getNumberDOF();
                for (int k = 0; k < ndf && p < n; ++k, ++p)
                    if (k < ndmA) rep.modelMass += mdiag[p];
            }
        }

        // --- SMS-SELFREPORT: an element with a positive self-reported bound carries a
        //     MASS-INDEPENDENT stability limit (e.g. a bipenalty/kinematic coupling:
        //     2*sqrt(m_internal/k)). Adding nodal mass cannot move it, so do NOT
        //     pretend to scale it — skip and count it (the integrator warns).
        if (ele->getExplicitCriticalTimeStep() > 0.0) {
            double selfDt = ele->getExplicitCriticalTimeStep();
            if (selfDt < dtTarget) rep.nSelfReport++;
            delete[] mdiag; continue;
        }

        // self-report-aware per-element stable step
        double dt_e = elementCriticalDt(ele, useTangent, mdiag, n);
        if (dt_e <= 0.0 || dt_e >= dtTarget) { delete[] mdiag; continue; }

        // s_e = (dt_target/dt_e)^2 ; add (s_e-1)*m to the element's nodes by lumped share
        double ratio = dtTarget / dt_e;
        double factor = ratio * ratio - 1.0;   // > 0

        // --- SMS-PARTIAL-INJECT: stage this element's increments LOCALLY, then commit
        //     to the shared `injected` map ONLY if the full DOF walk maps node-by-node
        //     (pos == n). A null node or non-node-major layout aborts the element with
        //     NOTHING committed (vs. the old half-injected, un-rolled-back state).
        std::map<int, Vector> staged;
        Node **nodes = ele->getNodePtrs();
        int numNodes = ele->getNumExternalNodes();
        int pos = 0;
        bool ok = (nodes != 0);
        double addedTrans = 0.0;
        for (int a = 0; a < numNodes && ok && pos < n; ++a) {
            Node *nd = nodes[a];
            if (nd == 0) { ok = false; break; }
            int tag = nd->getTag();
            int ndfN = nd->getNumberDOF();
            int ndmA = nd->getCrds().Size();
            Vector lv(ndfN);
            for (int k = 0; k < ndfN && pos < n; ++k, ++pos) {
                double dInc = factor * mdiag[pos];
                lv(k) = dInc;
                if (k < ndmA) addedTrans += dInc;
            }
            std::map<int, Vector>::iterator s = staged.find(tag);
            if (s == staged.end()) staged.insert(std::make_pair(tag, lv));
            else s->second += lv;   // node repeated within this element (degenerate)
        }
        delete[] mdiag;

        if (!ok || pos != n) { rep.nMismatch++; continue; }   // discard staged

        // commit: merge staged into the shared injected map (accumulate on shared nodes)
        for (std::map<int, Vector>::iterator s = staged.begin(); s != staged.end(); ++s) {
            std::map<int, Vector>::iterator it = injected.find(s->first);
            if (it == injected.end())
                injected.insert(*s);
            else
                it->second += s->second;
        }
        rep.addedMass += addedTrans;
        rep.nScaled++;
        if (dt_e < rep.minDtScaled) rep.minDtScaled = dt_e;
    }
    return rep;
}

// Commit (+1) or restore (-1) the recorded increments to the Domain's node masses.
// sign = +1 adds the increment (apply), sign = -1 subtracts it (re-baseline /
// teardown). Reads each node's current ndf x ndf mass, adjusts the diagonal, and
// writes it back (Node::setMass overwrites).
inline void
applyMassScaling(Domain *theDomain, const std::map<int, Vector> &injected, double sign)
{
    if (theDomain == 0) return;
    for (std::map<int, Vector>::const_iterator it = injected.begin();
         it != injected.end(); ++it) {
        Node *nd = theDomain->getNode(it->first);
        if (nd == 0) continue;
        const Vector &inj = it->second;
        Matrix m = nd->getMass();          // copy of the node's ndf x ndf mass
        int ndf = m.noRows();
        int lim = (inj.Size() < ndf) ? inj.Size() : ndf;
        for (int k = 0; k < lim; ++k)
            m(k, k) += sign * inj(k);
        nd->setMass(m);
    }
}

} // namespace Ladruno

#endif
