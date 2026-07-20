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
#include <MP_Constraint.h>      // Ladruno: constraint-exclusion guard (slave-node hazard)
#include <MP_ConstraintIter.h>
#include <LoadPattern.h>        // Ladruno (ADR-73 P3b): overlay collection (Kadd seam)
#include <LoadPatternIter.h>
#include <classTags.h>          // PATTERN_TAG_LadrunoPorousOverlay
#include <LadrunoPorousOverlay.h>  // Ladruno (ADR-73 P3b): getUndrainedAugmentation
                                   // (no cycle: the overlay header includes only
                                   //  LoadPattern/Vector/ID; CriticalTimeStep.cpp
                                   //  already includes it the same way)
#include <FE_Element.h>         // Ladruno: consistent (Olovsson) path — element eqn map
#include <FE_EleIter.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <OPS_Globals.h>

#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <functional>   // Ladruno: parallel consistent PCG (ConsistentParOps callbacks)

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
    double minDtSelfReport; // smallest self-reported bound below dtTarget (it still
                          //   GOVERNS after scaling — mass cannot move it; <=0 none)
    int    nMismatch;     // elements skipped due to a non-node-major / DOF mismatch
    int    nConstrained;  // sub-target elements SKIPPED because they touch an MP-constrained
                          //   (slave) node — injected mass would transfer through T^T M T and
                          //   the dt boost would not land; they remain GOVERNING (see below)
    double minDtConstrained; // smallest dt_e among the excluded constrained elements (the
                          //   step that still governs because they were not scaled; <=0 none)
    int    nOverlayAugScaled; // Ladruno (ADR-73 P3b): overlay-owned (undrained-augmented)
                          //   elements that actually received scaling. The CONSISTENT
                          //   integrator qualifies its "stable after scaling" report line
                          //   when > 0 (Olovsson under-delivery, ADR-73 SS12 P3b item 5).
};

// --- Ladruno (ADR-73 P3b §3b.4): SMS + porous-overlay composability.
//     LUMPED builder = the composability DELIVERABLE: it prices the UNDRAINED
//     pencil AND delivers dtTarget (real nodal mass injection scales the full
//     coupled mode; certified march measured stable). CONSISTENT builder =
//     PRICED-BUT-WARNED: it prices the same undrained pencil (strictly better
//     reporting) but the Olovsson centroid-preserving M_bar leaves the
//     undrained coupling mode's rigid-translation component unscaled, so it
//     under-delivers dtTarget on overlay cells — a one-time LOUD warning fires
//     when any overlay-owned element is scaled (ADR-73 §12 P3b, measured).
//     BOTH builders collect the LadrunoPorousOverlay
//     patterns once per build and, per element, sum the getUndrainedAugmentation
//     blocks dK_e = Q_e S_e^-1 Q_e^T into a dense Matrix passed to
//     elementCriticalDt(..., &Kaug) — SMS then sizes against the UNDRAINED
//     per-element pencil. The closed-form scale s = T^2 + 2*T*c is UNCHANGED and
//     remains exact: dK_e is mass- and state-independent, so lambda_max still
//     scales as 1/s under mass injection — the augmented dt_e is the correct
//     sizing input with zero formula changes. Residual not-augmented cases stay
//     loud: size mismatch advises once + sizes DRAINED (mirroring
//     computeCriticalTimeStep), self-reported bounds win UNCORRECTED (the CTS
//     scan already advises for that combo), and a not-ready snapshot warns
//     loudly below (unreachable in practice — the snapshot fires in setDomain at
//     pattern creation, long before analysis setup; a FAILED snapshot aborts
//     every commit anyway). History: pre-P3b no caller passed Kadd and SMS +
//     overlay was UNSUPPORTED (blanket warnIfOverlayPresentSMS, RETIRED here) —
//     the discrete undrained factor measured ~26x on the e72 soft soil, ABOVE
//     the 13-21x material range (mode-shape excess) — see ADR-73 §12 P3/P3b.
inline void collectPorousOverlaysSMS(Domain *theDomain,
                                     std::vector<LadrunoPorousOverlay*> &overlays)
{
    overlays.clear();
    if (theDomain == 0) return;
    LoadPattern *lp;
    LoadPatternIter &lps = theDomain->getLoadPatterns();
    while ((lp = lps()) != 0)
        if (lp->getClassTag() == PATTERN_TAG_LadrunoPorousOverlay) {
            LadrunoPorousOverlay *ov = static_cast<LadrunoPorousOverlay*>(lp);
            if (!ov->snapshotReady())
                opserr << "WARNING LadrunoMassScaling -- porous overlay "
                       << ov->getTag()
                       << " has NO completed geometry snapshot at SMS sizing "
                          "time: its elements are sized against the DRAINED "
                          "pencil (optimistic). This should be unreachable -- "
                          "check the overlay snapshot errors above.\n";
            overlays.push_back(ov);
        }
}

// Per-element augmentation fold: sum dK_e over the overlays that own `ele`
// (overlapping overlays sum, ⟨A-13⟩). Returns true if Kaug holds >= 1 block.
// A size mismatch (exotic ndf layout) advises once and skips that block — the
// element is then sized against its DRAINED pencil, never a silently wrong
// augmented one (the computeCriticalTimeStep behavior, mirrored).
inline bool overlayAugmentationSMS(const std::vector<LadrunoPorousOverlay*> &overlays,
                                   Element *ele, int n, Matrix &Kaug)
{
    static bool sizeMismatchWarned = false;   // one-time, process-wide
    bool augmented = false;
    for (size_t ov = 0; ov < overlays.size(); ++ov) {
        Matrix Kadd;
        if (!overlays[ov]->getUndrainedAugmentation(ele->getTag(), Kadd))
            continue;
        if (Kadd.noRows() != n) {
            if (!sizeMismatchWarned) {
                sizeMismatchWarned = true;
                opserr << "LadrunoMassScaling - NOTE overlay "
                       << overlays[ov]->getTag() << " augmentation for element "
                       << ele->getTag() << " is " << Kadd.noRows() << "x"
                       << Kadd.noRows() << " but the element pencil is "
                       << n << "x" << n << " (non-first-ndm DOF layout?): "
                          "augmentation SKIPPED for such elements - they are "
                          "sized against the DRAINED pencil (optimistic). "
                          "(printed once)\n";
            }
            continue;
        }
        if (!augmented) { Kaug = Kadd; augmented = true; }
        else              Kaug += Kadd;
    }
    return augmented;
}

// One-time INFO line when sizing augmented >= 1 element (battery-greppable
// honesty, pins §3b.4; static latch like the retired warning — one
// process-wide print).
inline void noteOverlayAugmentedSMS(int nAugmented)
{
    static bool infoShown = false;
    if (infoShown || nAugmented <= 0) return;
    infoShown = true;
    opserr << "LadrunoMassScaling: SMS sizing priced the UNDRAINED pencil for "
           << nAugmented << " overlay-owned elements (ADR-73 P3b)\n";
}

// Build the per-node fictitious-mass increment (additive diagonal) into `injected`
// (node tag -> Vector(ndf)) for the given target step. Does NOT touch the Domain;
// call applyMassScaling to commit. `injected` must be empty (freshly re-baselined).
inline MassScalingReport
buildMassScaling(AnalysisModel *theModel, double dtTarget, CTSLumping lumping,
                 bool useTangent, std::map<int, Vector> &injected)
{
    MassScalingReport rep; rep.addedMass = 0.0; rep.modelMass = 0.0;
    rep.nScaled = 0; rep.nElems = 0; rep.minDtScaled = dtTarget;
    rep.nSelfReport = 0; rep.minDtSelfReport = -1.0; rep.nMismatch = 0;
    rep.nConstrained = 0; rep.minDtConstrained = -1.0; rep.nOverlayAugScaled = 0;
    if (theModel == 0 || dtTarget <= 0.0) return rep;
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) return rep;

    // --- SMS-CONSTRAINTS (ADR-36 v1.1): collect the MP-constrained (SLAVE) node tags. A
    //     slave DOF is eliminated by the handler and any injected mass is redistributed to
    //     the retained node through T^T M T, so the per-element dt boost would NOT land —
    //     scaling such an element silently mis-distributes mass and under-delivers dt.
    //     EXCLUDE those elements (skip + count + report they still govern). SP/fix nodes are
    //     NOT a hazard (mass on a removed DOF is inert; refinement at a fixed support — the
    //     motivating case — must still scale), so only MP slaves are collected here.
    std::set<int> constrainedNodes;
    {
        MP_Constraint *theMP;
        MP_ConstraintIter &mps = theDomain->getMPs();
        while ((theMP = mps()) != 0)
            constrainedNodes.insert(theMP->getNodeConstrained());
    }

    // Ladruno (ADR-73 P3b §3b.4): collect porous overlays once per build (the
    // computeCriticalTimeStep idiom) — sizing prices the UNDRAINED pencil below.
    std::vector<LadrunoPorousOverlay*> overlays;
    collectPorousOverlaysSMS(theDomain, overlays);
    int nAugmented = 0;

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
            if (selfDt < dtTarget) {
                rep.nSelfReport++;
                if (rep.minDtSelfReport < 0.0 || selfDt < rep.minDtSelfReport)
                    rep.minDtSelfReport = selfDt;
            }
            delete[] mdiag; continue;
        }

        // --- Ladruno (ADR-73 P3b §3b.4): fold the overlay UNDRAINED augmentation
        //     into the sizing pencil (Kadd seam). Sizing math below is unchanged —
        //     the closed form s = T^2 + 2*T*c is scale-exact under the state-
        //     independent augmentation.
        Matrix Kaug;
        bool augmented = overlayAugmentationSMS(overlays, ele, n, Kaug);
        if (augmented) nAugmented++;

        // self-report-aware per-element UNDAMPED stable step (= 2/w0, w0 = sqrt(lambdaMax))
        double dt_e = elementCriticalDt(ele, useTangent, mdiag, n,
                                        augmented ? &Kaug : 0);
        if (dt_e <= 0.0) { delete[] mdiag; continue; }

        // --- SMS-BETAK: size against the DAMPED step. Stiffness-proportional Rayleigh
        //     damping shrinks the explicit stable step (xi = betaK*w/2 GROWS with
        //     w), so sizing against the undamped 2/w0 UNDER-scales -> still unstable at
        //     dtTarget. With c = betaK/dt_e (= 0.5*betaK*w0) the betaK-damped step at
        //     mass-scale s is dt_d(s) = (2/w0)(sqrt(s + c^2) - c) -- a CLOSED-FORM that
        //     inverts to s = T^2 + 2*T*c, T = dtTarget/dt_e. (Mass-proportional alphaM is
        //     intentionally excluded: it does NOT reduce the high-frequency step that
        //     governs explicit stability, and folding it in across scales is non-monotonic.
        //     Mirrors the betaK term of computeCriticalTimeStep's damped estimate.)
        //     Reduces EXACTLY to the undamped s=(dtTarget/dt_e)^2 when betaK=0 -> no-damping
        //     models are byte-identical. betaK is the TOTAL stiffness-proportional
        //     coefficient betaK + betaKinit + betaKcomm (stiffnessRayleighBeta): all three
        //     slots damp identically at high frequency, and reading coefs(1) alone silently
        //     UNDER-scaled `rayleigh 0 0 betaKinit 0` models.
        double betaK = stiffnessRayleighBeta(ele);
        double cDamp = betaK / dt_e;                                   // = 0.5*betaK*w0
        double dtDamped = dt_e * (std::sqrt(1.0 + cDamp * cDamp) - cDamp);   // current (s=1) damped step
        if (dtDamped >= dtTarget) { delete[] mdiag; continue; }        // already stable incl. betaK damping

        // --- SMS-CONSTRAINTS: this sub-target element WOULD be scaled, but if any of its
        //     nodes is an MP-constrained slave, skip it (mass would not land through the
        //     constraint). Count it and remember its (damped) dt — it remains GOVERNING, so
        //     the integrator can honestly tell the user dtTarget is not delivered for it.
        if (!constrainedNodes.empty()) {
            Node **cnds = ele->getNodePtrs();
            int cnn = ele->getNumExternalNodes();
            bool touchesConstrained = false;
            for (int a = 0; a < cnn && cnds; ++a) {
                if (cnds[a] != 0 && constrainedNodes.count(cnds[a]->getTag())) {
                    touchesConstrained = true; break;
                }
            }
            if (touchesConstrained) {
                rep.nConstrained++;
                if (rep.minDtConstrained < 0.0 || dtDamped < rep.minDtConstrained)
                    rep.minDtConstrained = dtDamped;
                delete[] mdiag; continue;
            }
        }

        // betaK-damped scale s = T^2 + 2*T*c (closed-form inverse of dt_d(s)=dtTarget); add
        // (s-1)*m to the element's nodes by lumped share.
        double T = dtTarget / dt_e;
        double factor = T * T + 2.0 * T * cDamp - 1.0;   // = s - 1 > 0

        // --- SMS-NDF-PREWALK (review 2026-07-02 hardening): reject a non-node-major
        //     layout UP FRONT. The staging loop below is bounded by pos < n, so an
        //     element whose nodes' TOTAL ndf exceeds n would stop mid-walk with
        //     pos == n exactly and COMMIT misaligned mass (the pos != n check after
        //     the loop cannot see it). Pre-walk the full sum, like the consistent
        //     builder's first pass, and skip + count unless it maps exactly.
        Node **nodes = ele->getNodePtrs();
        int numNodes = ele->getNumExternalNodes();
        {
            int sumNdf = 0;
            bool walkOK = (nodes != 0);
            for (int a = 0; a < numNodes && walkOK; ++a) {
                if (nodes[a] == 0) { walkOK = false; break; }
                sumNdf += nodes[a]->getNumberDOF();
            }
            if (!walkOK || sumNdf != n) {
                rep.nMismatch++; delete[] mdiag; continue;
            }
        }

        // --- SMS-PARTIAL-INJECT: stage this element's increments LOCALLY, then commit
        //     to the shared `injected` map ONLY if the full DOF walk maps node-by-node
        //     (pos == n). A null node or non-node-major layout aborts the element with
        //     NOTHING committed (vs. the old half-injected, un-rolled-back state).
        std::map<int, Vector> staged;
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
        if (augmented) rep.nOverlayAugScaled++;   // Ladruno (ADR-73 P3b)
        if (dtDamped < rep.minDtScaled) rep.minDtScaled = dtDamped;   // governing (damped) step
    }
    noteOverlayAugmentedSMS(nAugmented);   // Ladruno (ADR-73 P3b): one-time INFO
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

// ==========================================================================
//  CONSISTENT (Olovsson) selective mass scaling — ADR 38 (v2 of ADR 36)
// ==========================================================================
//
// Instead of the isotropic diagonal increment (s-1)*m_a that lumped scaling adds
// to each node (which shifts ALL global frequencies, incl. f1), the consistent path
// injects the CENTROIDAL Olovsson element scaling mass
//
//     M_bar_e = beta_e * [ diag(m_a) - m m^T / M_e ]      (per spatial direction)
//
// whose row sums are ZERO -> M_bar_e * t_rigid = 0: rigid-body translation receives
// no added inertia, so f1 is preserved while the relative/deformation modes (which
// govern the explicit step) are loaded -> dt_e -> dtTarget exactly as in v1. The same
// beta_e = T^2 + 2*T*c - 1 (betaK-damped, T = dtTarget/dt_e, c = betaK/dt_e) is reused.
//
// The resulting M_tilde = M_lump + sum_e M_bar_e is SPD but NON-diagonal (the
// -m m^T/M_e term couples the element's nodes), so it CANNOT live in per-node
// Node::setMass and the trivial `system Diagonal` a = M^-1 r no longer applies. The
// integrator solves M_tilde a = r matrix-free by Jacobi(lumped-diagonal)-preconditioned
// CG (consistentPCG below) — the lumped diagonal is the SOE diagonal the explicit run
// already assembles. Proven by the numpy oracle (Ladruno_implementation/
// mass_scaling_consistent/oracle_olovsson_sms.py): 12-21 CG iters, f1 drift < 0.3% vs
// 54-83% for lumped, at the same dtTarget. See 38_ladruno_consistent_mass_scaling_adr.md.

// One scaled element's consistent scaling mass + its global equation-number map.
//   eleTag  : the element's tag (Ladruno V4 — keys the EnergyBalanceRecorder energy
//             conduit so the recorder can add 1/2 v^T M_bar v per element; the matvec
//             itself never needs it).
//   eqn(i)  : SOE equation number of local DOF i (node-major); < 0 if the DOF is
//             constrained/fixed (eliminated) — skipped in the mat-vec.
//   Mbar    : n x n centroidal scaling mass (beta already applied; nonzero only on
//             same-direction translational DOF pairs). Node-major local-DOF layout.
struct ConsistentBlock {
    int    eleTag;
    ID     eqn;
    Matrix Mbar;
    ConsistentBlock(int tag, const ID &e, const Matrix &m) : eleTag(tag), eqn(e), Mbar(m) {}
};

// Build the per-element consistent scaling blocks for the given target step. Iterates
// the AnalysisModel's FE_Elements (for the equation-number map getID()), reusing the
// SAME per-element sizing the lumped path uses (dt_e via the shared CriticalTimeStep
// kernel, self-report skip, betaK-damped beta, MP-slave constraint exclusion). Does NOT
// touch the Domain or any Node mass — the blocks are integrator-owned scratch consumed
// only by consistentMatVec / consistentPCG. `blocks` must be empty.
inline MassScalingReport
buildMassScalingConsistent(AnalysisModel *theModel, double dtTarget, CTSLumping lumping,
                           bool useTangent, std::vector<ConsistentBlock> &blocks)
{
    MassScalingReport rep; rep.addedMass = 0.0; rep.modelMass = 0.0;
    rep.nScaled = 0; rep.nElems = 0; rep.minDtScaled = dtTarget;
    rep.nSelfReport = 0; rep.minDtSelfReport = -1.0; rep.nMismatch = 0;
    rep.nConstrained = 0; rep.minDtConstrained = -1.0; rep.nOverlayAugScaled = 0;
    if (theModel == 0 || dtTarget <= 0.0) return rep;
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) return rep;

    // Ladruno (ADR-73 P3b §3b.4): collect porous overlays once per build — the
    // consistent path prices the SAME undrained pencil as the lumped one (the
    // pre-P3b verifier caught this path silently drained-priced too). NOTE the
    // consistent path prices-but-WARNS: delivery is not honored for overlay
    // cells (see the under-delivery warning at the end of this builder).
    std::vector<LadrunoPorousOverlay*> overlays;
    collectPorousOverlaysSMS(theDomain, overlays);
    int nAugmented = 0;
    int nAugScaled = 0;   // augmented AND sub-target (received an M_bar block)

    // MP-constrained node tags — STRICTER than the lumped path: exclude any sub-target
    // element touching EITHER a slave (getNodeConstrained) OR a master/retained
    // (getNodeRetained) MP node. Rationale: the consistent path is matrix-free in the
    // EQUATION-number space (FE_Element::getID()), so it bypasses the constraint handler's
    // assembly. Under the default TransformationConstraintHandler an element attached to a
    // constrained node has a TRANSFORMED FE id — a different basis AND a different size
    // (a retained node's TransformationDOF_Group absorbs the slave's DOFs, so getID().Size()
    // > the element's own n) — so scattering the ORIGINAL-basis n x n M_bar through that id
    // is both wrong-basis and an out-of-bounds read. (The lumped path can stay slave-only
    // because it injects into Node mass, which the handler then transforms correctly via
    // T^T M T.) Excluded elements are counted (nConstrained) and reported as still governing.
    // This also keeps M_bar off every projector (ADR-30) equation, so the lumped-diagonal
    // projector mass stays exact.
    std::set<int> constrainedNodes;
    {
        MP_Constraint *theMP;
        MP_ConstraintIter &mps = theDomain->getMPs();
        while ((theMP = mps()) != 0) {
            constrainedNodes.insert(theMP->getNodeConstrained());
            constrainedNodes.insert(theMP->getNodeRetained());
        }
    }

    FE_Element *fePtr;
    FE_EleIter &fes = theModel->getFEs();
    while ((fePtr = fes()) != 0) {
        Element *ele = fePtr->getElement();
        if (ele == 0) continue;
        rep.nElems++;

        const Matrix &M = ele->getMass();
        int n = M.noRows();
        if (n == 0) continue;
        double *mdiag = new double[n];
        lumpElementMass(ele, M, lumping, mdiag);

        // model-mass denominator (sum of element translational lumped mass), same as v1.
        Node **nds = ele->getNodePtrs();
        int nn = ele->getNumExternalNodes();
        {
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

        // self-report: mass-independent bound -> scaling can't help; skip + count.
        if (ele->getExplicitCriticalTimeStep() > 0.0) {
            double selfDt = ele->getExplicitCriticalTimeStep();
            if (selfDt < dtTarget) {
                rep.nSelfReport++;
                if (rep.minDtSelfReport < 0.0 || selfDt < rep.minDtSelfReport)
                    rep.minDtSelfReport = selfDt;
            }
            delete[] mdiag; continue;
        }

        // Ladruno (ADR-73 P3b §3b.4): fold the overlay UNDRAINED augmentation
        // (Kadd seam) — identical to the lumped builder.
        Matrix Kaug;
        bool augmented = overlayAugmentationSMS(overlays, ele, n, Kaug);
        if (augmented) nAugmented++;

        double dt_e = elementCriticalDt(ele, useTangent, mdiag, n,
                                        augmented ? &Kaug : 0);
        if (dt_e <= 0.0) { delete[] mdiag; continue; }

        // betaK-damped sizing (identical closed form to the lumped path; betaK is the
        // summed betaK + betaKinit + betaKcomm, see stiffnessRayleighBeta).
        double betaK = stiffnessRayleighBeta(ele);
        double cDamp = betaK / dt_e;
        double dtDamped = dt_e * (std::sqrt(1.0 + cDamp * cDamp) - cDamp);
        if (dtDamped >= dtTarget) { delete[] mdiag; continue; }   // already stable

        // constraint exclusion (MP slave): skip + count + remember governing dt.
        if (!constrainedNodes.empty()) {
            bool touches = false;
            for (int a = 0; a < nn && nds; ++a)
                if (nds[a] != 0 && constrainedNodes.count(nds[a]->getTag())) { touches = true; break; }
            if (touches) {
                rep.nConstrained++;
                if (rep.minDtConstrained < 0.0 || dtDamped < rep.minDtConstrained)
                    rep.minDtConstrained = dtDamped;
                delete[] mdiag; continue;
            }
        }

        double T = dtTarget / dt_e;
        double beta = T * T + 2.0 * T * cDamp - 1.0;   // = s - 1 > 0

        // ---- build the centroidal scaling mass M_bar (n x n) -------------------
        // Walk the element node-by-node to recover, for each translational direction d,
        // the local DOF positions and their lumped masses; then for each node pair (a,b)
        //   M_bar[loc_a_d, loc_b_d] += beta*(m_a_d * delta_ab - m_a_d*m_b_d / M_e_d).
        // Non-translational (rotational / u-p) DOFs get no scaling mass (M_bar row = 0),
        // mirroring the lumped path which only adds to k < ndmA.
        Matrix Mbar(n, n);
        Mbar.Zero();

        // first pass: per-node DOF base offset + translational dim (node-major layout).
        // pos must reach exactly n (node-major mass), else this element's layout is not
        // node-major -> skip (count nMismatch), consistent with the lumped partial-inject.
        // ndmOf is clamped per node by the node's OWN ndf (review 2026-07-02 hardening):
        // a node with ndf < ndm would otherwise let the direction loop below index
        // base[a]+d past its DOF block — into the next node's block, and off the end of
        // Mbar/mdiag on the last node. No behavior change for ndf >= ndm (all solids).
        int *base = new int[nn];
        int *ndmOf = new int[nn];
        int pos = 0; bool ok = (nds != 0);
        for (int a = 0; a < nn && ok; ++a) {
            Node *ndp = nds[a];
            if (ndp == 0) { ok = false; break; }
            base[a] = pos;
            int ndmA = ndp->getCrds().Size();
            int ndfA = ndp->getNumberDOF();
            ndmOf[a] = (ndfA < ndmA) ? ndfA : ndmA;
            pos += ndfA;
        }
        if (!ok || pos != n) {
            delete[] base; delete[] ndmOf; delete[] mdiag;
            rep.nMismatch++; continue;
        }

        // common translational dimension across the element's nodes (min ndm seen).
        int ndmE = ndmOf[0];
        for (int a = 1; a < nn; ++a) if (ndmOf[a] < ndmE) ndmE = ndmOf[a];

        double addedTrans = 0.0;
        for (int d = 0; d < ndmE; ++d) {
            double Med = 0.0;
            for (int a = 0; a < nn; ++a) Med += mdiag[base[a] + d];
            if (Med <= 0.0) continue;   // massless direction -> no scaling mass
            for (int a = 0; a < nn; ++a) {
                double ma = mdiag[base[a] + d];
                int ia = base[a] + d;
                for (int b = 0; b < nn; ++b) {
                    double mb = mdiag[base[b] + d];
                    int ib = base[b] + d;
                    double val = beta * (((a == b) ? ma : 0.0) - ma * mb / Med);
                    Mbar(ia, ib) += val;
                }
                addedTrans += beta * (ma - ma * ma / Med);   // diagonal trace contribution
            }
        }
        delete[] base; delete[] ndmOf; delete[] mdiag;

        // DEFENSIVE: the matvec scatters the n x n M_bar through this id, so the id MUST
        // be the plain node-major map of size n. A surviving element touches no MP node, so
        // under the standard handlers getID().Size()==n; but guard anyway against any FE-layer
        // DOF expansion (transformed/exotic handler) -> skip (count nMismatch) rather than
        // read M_bar out of bounds.
        const ID &feID = fePtr->getID();
        if (feID.Size() != n) { rep.nMismatch++; continue; }

        // store the block with this element's tag (Ladruno V4 energy conduit) and
        // equation-number map.
        blocks.push_back(ConsistentBlock(ele->getTag(), feID, Mbar));
        if (augmented) { nAugScaled++; rep.nOverlayAugScaled++; }   // overlay-owned element actually scaled (ADR-73 P3b)
        rep.addedMass += addedTrans;
        rep.nScaled++;
        if (dtDamped < rep.minDtScaled) rep.minDtScaled = dtDamped;
    }
    // Ladruno (ADR-73 §12 P3b, measured): consistent-SMS UNDER-DELIVERY on
    // overlay cells. The Olovsson centroid-preserving M_bar scales only the
    // non-rigid element modes, but the overlay's undrained coupling mode
    // carries a large rigid-translation component that stays UNSCALED — the
    // coupled frequency scales by LESS than sqrt(s), so consistent SMS cannot
    // deliver dtTarget against the undrained pencil (the certified march
    // diverges geometrically from step 1). The lumped builder injects real
    // nodal mass and CAN deliver. Pricing stays undrained above (strictly
    // better reporting); delivery is warned loudly here.
    if (nAugScaled > 0) {
        static bool consistentUnderDeliveryWarned = false;   // one-time latch
        if (!consistentUnderDeliveryWarned) {
            consistentUnderDeliveryWarned = true;
            opserr << "WARNING LadrunoMassScaling (consistent) -- " << nAugScaled
                   << " overlay-owned element(s) scaled: the Olovsson "
                      "centroid-preserving blocks under-scale the undrained "
                      "COUPLING mode (measured under-delivery, ADR-73 §12 P3b) "
                      "-- the certified dtTarget is NOT honored for them; use "
                      "lumped SMS (CentralDifferenceSMS) or size dt from the "
                      "overlay-aware criticalTimeStep() report. (printed once)\n";
        }
    }
    noteOverlayAugmentedSMS(nAugmented);   // Ladruno (ADR-73 P3b): one-time INFO
    return rep;
}

// Matrix-free apply  y = M_tilde x = (diag) .* x + sum_e M_bar_e x_e.
//   diagMinv : the SOE diagonal AFTER the in-place factor (DiagonalDirectSolver stores
//              1/mass), length neq. The lumped diagonal is 1/diagMinv[i]; this is also
//              the Jacobi preconditioner. Pass it directly to avoid recomputing.
// y and x are full SOE-length vectors. Constrained DOFs (eqn < 0) are skipped.
inline void
consistentMatVec(const std::vector<ConsistentBlock> &blocks,
                 const double *diagMinv, int neq, const Vector &x, Vector &y)
{
    for (int i = 0; i < neq; ++i)
        y(i) = (diagMinv[i] != 0.0 ? (x(i) / diagMinv[i]) : 0.0);   // diag mass = 1/diagMinv
    for (size_t e = 0; e < blocks.size(); ++e) {
        const ConsistentBlock &blk = blocks[e];
        const ID &id = blk.eqn;
        const Matrix &Mb = blk.Mbar;
        int nl = id.Size();
        for (int a = 0; a < nl; ++a) {
            int ea = id(a);
            if (ea < 0 || ea >= neq) continue;
            double acc = 0.0;
            for (int b = 0; b < nl; ++b) {
                int eb = id(b);
                if (eb < 0 || eb >= neq) continue;
                acc += Mb(a, b) * x(eb);
            }
            y(ea) += acc;
        }
    }
}

// Solve M_tilde a = r matrix-free by Jacobi(diag)-preconditioned CG, in place: `a`
// enters as the diagonal solve M_lump^-1 r (the starting guess AND the way r is
// recovered by the caller) and exits as M_tilde^-1 r. diagMinv = 1/lumped-mass
// (the SOE factored diagonal = the preconditioner). Returns the iteration count;
// writes the achieved relative residual to *relResid if non-null.
inline int
consistentPCG(const std::vector<ConsistentBlock> &blocks,
              const double *diagMinv, int neq, const Vector &r, Vector &a,
              double tol = 1.0e-10, int maxit = 200, double *relResid = 0)
{
    Vector Ax(neq), z(neq), p(neq), Ap(neq);
    // start from a = M_lump^-1 r (already in `a` on entry); residual = r - M_tilde a
    consistentMatVec(blocks, diagMinv, neq, a, Ax);
    Vector res(neq);
    for (int i = 0; i < neq; ++i) res(i) = r(i) - Ax(i);
    for (int i = 0; i < neq; ++i) z(i) = diagMinv[i] * res(i);
    p = z;
    double rz = res ^ z;
    double rnorm0 = r.Norm();
    if (rnorm0 <= 0.0) rnorm0 = 1.0;
    int k = 0;
    for (k = 1; k <= maxit; ++k) {
        consistentMatVec(blocks, diagMinv, neq, p, Ap);
        double pAp = p ^ Ap;
        if (pAp <= 0.0) break;                  // breakdown guard (M_tilde SPD => pAp>0)
        double alpha = rz / pAp;
        a.addVector(1.0, p, alpha);
        res.addVector(1.0, Ap, -alpha);
        if (res.Norm() <= tol * rnorm0) break;
        for (int i = 0; i < neq; ++i) z(i) = diagMinv[i] * res(i);
        double rzNew = res ^ z;
        double bcg = rzNew / rz;
        for (int i = 0; i < neq; ++i) p(i) = z(i) + bcg * p(i);
        rz = rzNew;
    }
    if (relResid) *relResid = res.Norm() / rnorm0;
    return (k > maxit) ? maxit : k;   // loop exits at maxit+1 if never converged; clamp
}

// ===========================================================================
// Ladruno: DISTRIBUTED (MPI) consistent PCG — the parallel sibling of the two
// functions above. Solves M_tilde a = r matrix-free across MPI ranks where the
// model is manually partitioned (OpenSeesMP, `system MPIDiagonal`) and a
// partition-boundary node DOF is SHARED by several ranks.
//
// The whole correctness trick is one weight vector  w_i = 1 / multiplicity_i
// (multiplicity = number of ranks that hold DOF i; 1 for purely-local DOFs):
//
//   * Matrix-free matvec  y = M_tilde x:  the lumped diagonal (the SOE's
//     getDiagonalA() is the GLOBAL, already cross-rank-reduced 1/mass at shared
//     DOFs and is therefore replicated on every sharing rank) is applied
//     WEIGHTED by w_i, the off-diagonal M_bar_e blocks (each element lives on
//     exactly one rank) are applied in FULL, and then the shared entries of y
//     are summed across ranks (assembleShared). The sum turns w_i * (full diag)
//     over multiplicity ranks back into the full diagonal exactly once, while
//     the off-diagonal contributions of all element-owning ranks accumulate.
//   * Global inner products  <x,y> = globalSum( sum_i w_i x_i y_i ):  the weight
//     makes each shared DOF count once across the ranks that hold it.
//
// All CG control scalars (rz, pAp, residual norm) are GLOBAL reductions, so the
// iteration count is identical on every rank -> the per-iteration assembleShared
// / globalSum collectives stay in lockstep (no deadlock). With w == null,
// assembleShared == null and globalSum == null this reduces EXACTLY to the
// serial consistentMatVec / consistentPCG above (the serial callers keep using
// those; this path is only taken when a distributed SOE is detected).
// ===========================================================================
struct ConsistentParOps {
    const double *w;                          // length neq; w[i]=1/mult (1 if non-shared). null => all ones
    std::function<void(double *, int)> assembleShared;  // sum v[0..n) shared entries across ranks, in place. null => no-op
    std::function<double(double)>      globalSum;        // sum a scalar across ranks. null => identity
    ConsistentParOps() : w(0) {}
};

inline double
consistentParDot(const Vector &x, const Vector &y, int neq, const ConsistentParOps &ops)
{
    const double *w = ops.w;
    double s = 0.0;
    for (int i = 0; i < neq; ++i)
        s += (w ? w[i] : 1.0) * x(i) * y(i);
    return ops.globalSum ? ops.globalSum(s) : s;
}

// y = M_tilde x  (distributed). x must be consistent (identical on every rank) at
// shared DOFs on entry; y is consistent at shared DOFs on exit.
inline void
consistentParMatVec(const std::vector<ConsistentBlock> &blocks, const double *diagMinv,
                    int neq, const Vector &x, Vector &y, const ConsistentParOps &ops)
{
    const double *w = ops.w;
    for (int i = 0; i < neq; ++i) {
        double diag = (diagMinv[i] != 0.0) ? (x(i) / diagMinv[i]) : 0.0;
        y(i) = (w ? w[i] : 1.0) * diag;       // diagonal weighted -> assembled to full once
    }
    for (size_t e = 0; e < blocks.size(); ++e) {
        const ConsistentBlock &blk = blocks[e];
        const ID &id = blk.eqn;
        const Matrix &Mb = blk.Mbar;
        int nl = id.Size();
        for (int a = 0; a < nl; ++a) {
            int ea = id(a);
            if (ea < 0 || ea >= neq) continue;
            double acc = 0.0;
            for (int b = 0; b < nl; ++b) {
                int eb = id(b);
                if (eb < 0 || eb >= neq) continue;
                acc += Mb(a, b) * x(eb);
            }
            y(ea) += acc;                     // off-diagonal applied in FULL (element is rank-local)
        }
    }
    if (ops.assembleShared)
        ops.assembleShared(&y(0), neq);       // sum shared entries across ranks
}

inline int
consistentParPCG(const std::vector<ConsistentBlock> &blocks, const double *diagMinv,
                 int neq, const Vector &r, Vector &a, const ConsistentParOps &ops,
                 double tol = 1.0e-10, int maxit = 200, double *relResid = 0)
{
    Vector Ax(neq), z(neq), p(neq), Ap(neq), res(neq);
    consistentParMatVec(blocks, diagMinv, neq, a, Ax, ops);   // a enters as M_lump^-1 r
    for (int i = 0; i < neq; ++i) res(i) = r(i) - Ax(i);
    for (int i = 0; i < neq; ++i) z(i) = diagMinv[i] * res(i);
    p = z;
    double rz = consistentParDot(res, z, neq, ops);
    double rnorm0 = std::sqrt(consistentParDot(r, r, neq, ops));
    if (rnorm0 <= 0.0) rnorm0 = 1.0;
    int k = 0;
    for (k = 1; k <= maxit; ++k) {
        consistentParMatVec(blocks, diagMinv, neq, p, Ap, ops);
        double pAp = consistentParDot(p, Ap, neq, ops);
        if (pAp <= 0.0) break;                // breakdown guard (M_tilde SPD => pAp>0 globally)
        double alpha = rz / pAp;
        a.addVector(1.0, p, alpha);
        res.addVector(1.0, Ap, -alpha);
        double resnorm = std::sqrt(consistentParDot(res, res, neq, ops));
        if (resnorm <= tol * rnorm0) break;
        for (int i = 0; i < neq; ++i) z(i) = diagMinv[i] * res(i);
        double rzNew = consistentParDot(res, z, neq, ops);
        double bcg = rzNew / rz;
        for (int i = 0; i < neq; ++i) p(i) = z(i) + bcg * p(i);
        rz = rzNew;
    }
    if (relResid) *relResid = std::sqrt(consistentParDot(res, res, neq, ops)) / rnorm0;
    return (k > maxit) ? maxit : k;
}

} // namespace Ladruno

#endif
