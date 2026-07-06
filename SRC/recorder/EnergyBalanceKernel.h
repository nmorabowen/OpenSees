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

#ifndef EnergyBalanceKernel_h
#define EnergyBalanceKernel_h

// EnergyBalanceKernel.h — the ONE definition of the structural-dynamics
// energy-balance math (ADR D8). Header-only / inline so both the global
// EnergyBalanceRecorder and ladruno::EnergyBalanceSource link against a single
// copy. Lives in namespace ebkernel (NOT ladruno) so it is usable from the
// global-namespace recorder AND from the ladruno recorder code.
//
// 6 components, in this order: KE, IE, DW, ULW, RES, ERR where
//   KE  = 1/2 v^T M v              (instantaneous; element + nodal mass)
//   IE  = int F_resisting^T v dt   (trapezoidal)
//   DW  = int v^T C v dt           (trapezoidal; element + nodal Rayleigh damping)
//   ULW = int v^T P_ext dt         (trapezoidal; Node::getUnbalancedLoad)
//   RES = ULW - (KE + IE + DW)     (closure residual)
//   ERR = 100 * RES / E_ref        (E_ref = running max of |KE|+|IE|+|DW|+|ULW|,
//                                    per scope)
//
// The two per-entity helpers (addElementEnergy / addNodeEnergy) are lifted
// VERBATIM from the EnergyBalanceRecorder anonymous namespace, made inline.
// Every documented quirk is preserved:
//   * element mass is COPIED before getDamp() (elements share a static scratch
//     Matrix between getMass/getDamp);
//   * velocities gathered DOF-by-DOF with zero-fill;
//   * empty (0x0) mass/damp matrices are skipped;
//   * F = element resisting force dotted with v.

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <MeshRegion.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>

#include <LadrunoMassScalingEnergy.h>   // Ladruno V4: consistent (Olovsson) scaling-mass KE
#include <LadrunoEnergyChannels.h>      // Ladruno ADR-69: named energy channels (v2)

#include <cmath>

namespace ebkernel {

// Single source of truth for the component count. The recorder still defines
// EBR_NUM_ENERGY_COMPONENTS (== 6) for its own response sizing; this constant
// must agree with it.
static const int NUM_ENERGY_COMPONENTS = 6;

// Accumulate the energy contributions of a single element. vel is a caller-
// supplied scratch vector (>= numDOF) so the hot loop performs no allocation.
//   kinetic       += 1/2 v^T M v   (instantaneous kinetic energy)
//   internal_rate += F^T v         (internal power; caller integrates over dt)
//   damping_rate  += v^T C v       (damping power;  caller integrates over dt)
inline void addElementEnergy(Element *ele, Vector &vel,
                             double &kinetic,
                             double &internal_rate,
                             double &damping_rate)
{
    // COPY the mass: many elements return a reference to a shared static matrix
    // from BOTH getMass() and getDamp() (and getDamp() internally calls getMass/
    // getTangentStiff), so a reference here would be clobbered by getDamp() below
    // and the kinetic energy would be computed from the damping matrix.
    Matrix M = ele->getMass();
    const Matrix &C = ele->getDamp();
    const Vector &F = ele->getResistingForce();
    const int numExternalNodes = ele->getNumExternalNodes();
    const int numDOF = ele->getNumDOF();
    Node **elenodes = ele->getNodePtrs();

    int cnt = 0;
    for (int i = 0; i < numExternalNodes && cnt < numDOF; ++i) {
        const Vector &node_vel = elenodes[i]->getVel();
        for (int j = 0; j < node_vel.Size() && cnt < numDOF; ++j)
            vel(cnt++) = node_vel(j);
    }
    for (; cnt < numDOF; ++cnt)
        vel(cnt) = 0.0;

    // Some elements return an empty (0x0) mass or damping matrix when they
    // carry none; only accumulate when the matrix is the expected size.
    const bool haveMass = (M.noRows() == numDOF && M.noCols() == numDOF);
    const bool haveDamp = (C.noRows() == numDOF && C.noCols() == numDOF);
    if (haveMass)
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);

    // Ladruno V4 — CONSISTENT (Olovsson) selective mass scaling injects a cross-node
    // element scaling mass M_bar_e that lives only inside the integrator (never in
    // Node/Element getMass), so getMass() above sees only the lumped diagonal. Add the
    // missing 1/2 v^T M_bar_e v here (vel is already this element's node-major velocity)
    // so the consistent-path energy balance closes. The registry is empty for the lumped
    // path and every base integrator -> active() is false -> no extra work, no double
    // counting. See SRC/analysis/integrator/LadrunoMassScalingEnergy.h.
    {
        const Ladruno::MassScalingEnergyRegistry &reg =
            Ladruno::MassScalingEnergyRegistry::instance();
        if (reg.active())
            kinetic += reg.elementScalingKE(ele->getTag(), vel, numDOF);
    }

    if (haveDamp)
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                damping_rate += C(i, j) * vel(i) * vel(j);

    if (F.Size() == numDOF) {
        double dot = 0.0;
        for (int i = 0; i < numDOF; ++i)
            dot += F(i) * vel(i);
        internal_rate += dot;
    }
}

// Accumulate the nodal contributions: lumped-mass kinetic energy, nodal
// Rayleigh damping (alphaM * nodalMass), and the work rate of the external
// applied nodal load. M is fully consumed before getDamp() is called because
// Node::getMass()/getDamp() can share the same scratch matrix when massless.
inline void addNodeEnergy(Node *node,
                          double &kinetic,
                          double &damping_rate,
                          double &unbalanced_rate)
{
    const Vector &vel = node->getVel();
    const int n = vel.Size();

    const Matrix &M = node->getMass();
    if (M.noRows() == n && M.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);

    const Matrix &C = node->getDamp();
    if (C.noRows() == n && C.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                damping_rate += C(i, j) * vel(i) * vel(j);

    const Vector &Pext = node->getUnbalancedLoad();
    if (Pext.Size() == n)
        unbalanced_rate += vel ^ Pext;
}

// Per-scope accumulator: encapsulates the trapezoidal work integrals + the
// closure/ERR arithmetic, reproducing EnergyBalanceRecorder::record() exactly.
struct EnergyAccumulator {
    double internal_energy, damping_work, unbalanced_load_work;
    double prev_internal_rate, prev_damping_rate, prev_unbalanced_rate;
    double eref;          // running max of |KE|+|IE|+|DW|+|ULW| for THIS scope

    EnergyAccumulator() { reset(); }

    void reset() {
        internal_energy = 0.0;
        damping_work = 0.0;
        unbalanced_load_work = 0.0;
        prev_internal_rate = 0.0;
        prev_damping_rate = 0.0;
        prev_unbalanced_rate = 0.0;
        eref = 0.0;
    }

    // Integrate one step. `firstRecord` => only seed the prev rates (the first
    // recorded step uses the rectangle rule so no increment is dropped, exactly
    // matching the recorder), thereafter trapezoidal.
    // Fills out[6] = {KE, IE, DW, ULW, RES, ERR}.
    void step(double dT, bool firstRecord,
              double ke, double internal_rate, double damping_rate, double unbalanced_rate,
              double out[NUM_ENERGY_COMPONENTS])
    {
        if (firstRecord) {
            // first recorded step: rectangle rule (there is no previous-step
            // rate yet). This reproduces EnergyBalanceRecorder's original
            // record() EXACTLY (its firstRecord branch is rectangle, not
            // trapezoidal/seed-only — verified against the pre-refactor source).
            // Skipping it would drop one increment of work and leave a constant
            // residual offset. NB the caller passes dT = thisStepTime - prevTime;
            // EnergyBalanceSource starts prevTime at 0, so its first-step dT
            // matches EnergyBalanceRecorder's (time_last seeded from the start
            // time) for the usual t0==0 analysis. A dynamic stage that starts at
            // t0!=0 (e.g. after a gravity stage without a time reset) would take
            // a larger first increment in the source than in the text sidecar;
            // that is a documented v2 caveat, not a closure error (RES tracks it).
            internal_energy      += internal_rate   * dT;
            damping_work         += damping_rate    * dT;
            unbalanced_load_work += unbalanced_rate * dT;
        } else {
            // trapezoidal rule
            internal_energy      += 0.5 * (prev_internal_rate   + internal_rate)   * dT;
            damping_work         += 0.5 * (prev_damping_rate    + damping_rate)    * dT;
            unbalanced_load_work += 0.5 * (prev_unbalanced_rate + unbalanced_rate) * dT;
        }

        prev_internal_rate   = internal_rate;
        prev_damping_rate    = damping_rate;
        prev_unbalanced_rate = unbalanced_rate;

        const double sum = ke + internal_energy + damping_work;
        const double res = unbalanced_load_work - sum;
        // energy scale = running max of summed component magnitudes (does not
        // collapse to ~0 when KE and incremental IE cancel, e.g. free vibration)
        const double mag = std::fabs(ke) + std::fabs(internal_energy)
                         + std::fabs(damping_work) + std::fabs(unbalanced_load_work);
        if (mag > eref) eref = mag;
        const double err = (eref > 1.0e-16) ? 100.0 * res / eref : 0.0;

        out[0] = ke;
        out[1] = internal_energy;
        out[2] = damping_work;
        out[3] = unbalanced_load_work;
        out[4] = res;
        out[5] = err;
    }
};

// ---------------------------------------------------------------------------
// Ladruno ADR-69 (v2): channel-aware accumulator with element/nodal split.
//
// Differences from the V1 EnergyAccumulator above (which is frozen — the
// legacy 6-column path must stay byte-identical):
//   * KE and DW are carried as ELEMENT and NODAL parts separately. Under MPI
//     the element parts reduce correctly by naive cross-rank summation
//     (elements are uniquely partitioned); the nodal parts are multiply-
//     counted on shared partition-boundary nodes and are quarantined in
//     their own columns so an offline reducer can treat them deliberately.
//   * Reads the ADR-69 energy channels: ABSORB_LEAK L (the resisting-sign
//     injection work that pollutes the raw IE integral) and LNVD_WORK D.
//     Both are cumulative registry totals; this struct snapshots a baseline
//     at the first record and uses deltas, so a fresh recorder starts at
//     zero regardless of registry history.
//   * Balance:  IE_display = IE_raw - L        (rebucket: truthful IE)
//               E_inject   = -L                (source-side, positive = in)
//               RES = (ULW + E_inject)
//                     - (KE + IE_display + DW + E_lnvd)
//     The L terms cancel in RES (rebucketing never fakes closure); the LNVD
//     term genuinely closes what v1 leaked into RES.
//
// Output slot order (fixed; the recorder skips undeclared channel slots when
// writing): KE_ele KE_nod IE DW_ele DW_nod ULW E_inject E_lnvd RES ERR.
// ---------------------------------------------------------------------------
static const int NUM_V2_SLOTS = 10;
enum V2Slot {
    V2_KE_ELE = 0, V2_KE_NOD, V2_IE, V2_DW_ELE, V2_DW_NOD, V2_ULW,
    V2_E_INJECT, V2_E_LNVD, V2_RES, V2_ERR
};

struct EnergyAccumulatorV2 {
    double internal_energy;                       // raw IE (may contain leak)
    double damping_work_ele, damping_work_nod;
    double unbalanced_load_work;
    double prev_internal_rate;
    double prev_damping_rate_ele, prev_damping_rate_nod;
    double prev_unbalanced_rate;
    double eref;                                  // running max, this scope
    bool   baselined;
    double base_leak, base_lnvd;                  // registry snapshot at first record

    EnergyAccumulatorV2() { reset(); }

    void reset() {
        internal_energy = 0.0;
        damping_work_ele = 0.0;
        damping_work_nod = 0.0;
        unbalanced_load_work = 0.0;
        prev_internal_rate = 0.0;
        prev_damping_rate_ele = 0.0;
        prev_damping_rate_nod = 0.0;
        prev_unbalanced_rate = 0.0;
        eref = 0.0;
        baselined = false;
        base_leak = 0.0;
        base_lnvd = 0.0;
    }

    // Same integration convention as V1 step(): first recorded step uses the
    // rectangle rule, thereafter trapezoidal. Channel totals are read from
    // the registry by the caller and passed in (kernel stays registry-free at
    // this level so the arithmetic is testable in isolation).
    void step(double dT, bool firstRecord,
              double ke_ele, double ke_nod,
              double internal_rate,
              double damping_rate_ele, double damping_rate_nod,
              double unbalanced_rate,
              double leak_total, double lnvd_total,
              double out[NUM_V2_SLOTS])
    {
        if (firstRecord) {
            internal_energy      += internal_rate    * dT;
            damping_work_ele     += damping_rate_ele * dT;
            damping_work_nod     += damping_rate_nod * dT;
            unbalanced_load_work += unbalanced_rate  * dT;
        } else {
            internal_energy      += 0.5 * (prev_internal_rate    + internal_rate)    * dT;
            damping_work_ele     += 0.5 * (prev_damping_rate_ele + damping_rate_ele) * dT;
            damping_work_nod     += 0.5 * (prev_damping_rate_nod + damping_rate_nod) * dT;
            unbalanced_load_work += 0.5 * (prev_unbalanced_rate  + unbalanced_rate)  * dT;
        }

        prev_internal_rate    = internal_rate;
        prev_damping_rate_ele = damping_rate_ele;
        prev_damping_rate_nod = damping_rate_nod;
        prev_unbalanced_rate  = unbalanced_rate;

        if (!baselined) {
            base_leak = leak_total;
            base_lnvd = lnvd_total;
            baselined = true;
        }
        const double L  = leak_total - base_leak;   // resisting-sign leak
        const double Dl = lnvd_total - base_lnvd;   // LNVD dissipation (>= 0)

        const double ke      = ke_ele + ke_nod;
        const double dw      = damping_work_ele + damping_work_nod;
        const double ie_disp = internal_energy - L;
        const double e_inj   = -L;

        const double res = (unbalanced_load_work + e_inj)
                         - (ke + ie_disp + dw + Dl);

        const double mag = std::fabs(ke) + std::fabs(ie_disp) + std::fabs(dw)
                         + std::fabs(unbalanced_load_work)
                         + std::fabs(e_inj) + std::fabs(Dl);
        if (mag > eref) eref = mag;
        const double err = (eref > 1.0e-16) ? 100.0 * res / eref : 0.0;

        out[V2_KE_ELE]   = ke_ele;
        out[V2_KE_NOD]   = ke_nod;
        out[V2_IE]       = ie_disp;
        out[V2_DW_ELE]   = damping_work_ele;
        out[V2_DW_NOD]   = damping_work_nod;
        out[V2_ULW]      = unbalanced_load_work;
        out[V2_E_INJECT] = e_inj;
        out[V2_E_LNVD]   = Dl;
        out[V2_RES]      = res;
        out[V2_ERR]      = err;
    }
};

// Sum instantaneous KE + the three work rates over the whole domain. Mirrors
// EnergyBalanceRecorder::record()'s whole-model element+node loops: for
// elements accumulate KE+IE+DW, for nodes accumulate KE+DW+ULW. `velScratch`
// is a caller-owned Vector sized to the max element DOF (avoids per-call alloc).
inline void sweepDomain(Domain* dom, Vector& velScratch,
                        double& ke, double& internal_rate,
                        double& damping_rate, double& unbalanced_rate)
{
    ke = 0.0; internal_rate = 0.0; damping_rate = 0.0; unbalanced_rate = 0.0;
    if (dom == 0)
        return;

    {
        Element *ele;
        ElementIter &elements = dom->getElements();
        while ((ele = elements()) != 0) {
            double eke = 0., eie = 0., edw = 0.;
            addElementEnergy(ele, velScratch, eke, eie, edw);
            ke += eke; internal_rate += eie; damping_rate += edw;
        }
    }
    {
        Node *node;
        NodeIter &nodes = dom->getNodes();
        while ((node = nodes()) != 0) {
            double nke = 0., ndw = 0., nulw = 0.;
            addNodeEnergy(node, nke, ndw, nulw);
            ke += nke; damping_rate += ndw; unbalanced_rate += nulw;
        }
    }
}

// Same, restricted to a MeshRegion's elements + nodes (element/node pointers
// resolved via the domain). Mirrors the recorder's per-region accumulation
// (elements -> KE+IE+DW, nodes -> KE+DW+ULW).
inline void sweepRegion(Domain* dom, MeshRegion* reg, Vector& velScratch,
                        double& ke, double& internal_rate,
                        double& damping_rate, double& unbalanced_rate)
{
    ke = 0.0; internal_rate = 0.0; damping_rate = 0.0; unbalanced_rate = 0.0;
    if (dom == 0 || reg == 0)
        return;

    const ID &elems = reg->getElements();
    for (int i = 0; i < elems.Size(); ++i) {
        Element *ele = dom->getElement(elems(i));
        if (ele != 0) {
            double eke = 0., eie = 0., edw = 0.;
            addElementEnergy(ele, velScratch, eke, eie, edw);
            ke += eke; internal_rate += eie; damping_rate += edw;
        }
    }
    const ID &rnodes = reg->getNodes();
    for (int i = 0; i < rnodes.Size(); ++i) {
        Node *node = dom->getNode(rnodes(i));
        if (node != 0) {
            double nke = 0., ndw = 0., nulw = 0.;
            addNodeEnergy(node, nke, ndw, nulw);
            ke += nke; damping_rate += ndw; unbalanced_rate += nulw;
        }
    }
}

} // namespace ebkernel

#endif // EnergyBalanceKernel_h
