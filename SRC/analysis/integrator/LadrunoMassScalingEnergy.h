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

#ifndef LadrunoMassScalingEnergy_h
#define LadrunoMassScalingEnergy_h

// Ladruno V4 — consistent (Olovsson) selective mass-scaling energy conduit.
// See Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md and the
// mass-scaling developer handoff (the "V4 — consistent-mass energy closure" item).
//
// THE PROBLEM. CONSISTENT (Olovsson) mass scaling injects a centroidal element
// scaling mass  M_bar_e = beta_e [diag(m_a) - m m^T / M_e]  that lives ONLY inside
// the integrator (a matrix-free PCG operand; it is NEVER written to Node/Element
// mass, because its cross-node coupling cannot live in a per-node diagonal). The
// EnergyBalanceRecorder builds its kinetic energy from Node::getMass()/Element::
// getMass(), so it sees the lumped diagonal M_lump but NOT the cross-node M_bar_e.
// The recorded KE is therefore 1/2 v^T M_lump v, short of the true run kinetic
// energy 1/2 v^T M_tilde v = 1/2 v^T (M_lump + sum_e M_bar_e) v, and the consistent
// run's energy balance does not close (the lumped sibling's does — it injects into
// Node mass, which the recorder reads).
//
// THE SEAM. The recorder holds only a Domain* — it has no handle on the active
// integrator. This process-global registry is that missing recorder<->integrator
// hook, kept minimal: the active CONSISTENT-SMS integrator PUBLISHES its per-element
// node-major M_bar blocks here at the end of domainChanged and CLEARS them on
// teardown; the shared energy kernel (EnergyBalanceKernel.h) QUERIES it per element
// and adds the missing 1/2 v^T M_bar_e v. M_bar_e is stored in the element's
// node-major local-DOF layout — exactly the order ebkernel::addElementEnergy already
// gathers nodal velocities — so the kernel can form the quadratic with the velocity
// vector it already builds, with no equation-number bookkeeping.
//
// SCOPE / LIFETIME. There is exactly one active transient integrator per process, so
// a global registry is the minimal correct seam. The lumped SMS path and every base
// integrator leave the registry empty (active() == false), so their recorded KE is
// untouched and there is no double counting. An `owner` token (the publishing
// integrator's `this`) guards clear(): a stale teardown of a superseded integrator
// cannot wipe a newer one's blocks. NOT parallel-shared-node aware — that is the
// separate T-MPI work item; the registry mirrors the sequential per-rank blocks.
//
// ADR-69 P3 (KE_ms advisory). The registry ALSO carries the LUMPED (DT2MS-style)
// SMS injection: the per-node additive diagonal ΔM the lumped integrators commit
// to Node mass (CentralDifferenceSMS / ExplicitBathe -sms). Unlike the consistent
// M_bar (which the recorder's getMass() CANNOT see, so its KE is ADDED to the
// balance), the lumped ΔM is already inside Node::getMass() — the recorder's KE
// is correct and nothing is added. The nodal side exists so the -v2 recorder can
// MEASURE the fictitious added-mass kinetic-energy share as the KE_ms advisory
// column (LS-DYNA GLSTAT-style fidelity metric): KE_ms = 1/2 v^T ΔM v (nodal,
// lumped) + sum_e 1/2 v^T M_bar_e v (element, consistent). Same owner-guarded
// publish/clear lifecycle, kept on its own owner token because the lumped and
// consistent publishers are different integrator objects.
//
// Header pulls only <map>, <Matrix.h>, <Vector.h> (no Domain/Element) so it is cheap
// to include from the header-only energy kernel without a dependency cycle.

#include <map>
#include <Matrix.h>
#include <Vector.h>

namespace Ladruno {

class MassScalingEnergyRegistry {
public:
    static MassScalingEnergyRegistry &instance();

    // Publish/replace the active consistent-SMS integrator's per-element node-major
    // scaling-mass blocks (element tag -> M_bar). `owner` identifies the publisher so
    // a later clear() from a superseded integrator cannot wipe a newer publisher.
    void publish(const void *owner, const std::map<int, Matrix> &blocks);

    // Clear ONLY if `owner` is the current publisher; a no-op otherwise. Call from the
    // integrator's destructor / recvSelf so the registry never outlives its blocks.
    void clear(const void *owner);

    // Fast guard for the recorder's hot loop: true iff some integrator has published
    // a non-empty block set. When false the energy kernel does zero map lookups.
    bool active() const { return !theBlocks.empty(); }

    // 1/2 v^T M_bar_e v for element `eleTag`; `vel` is the element's node-major
    // velocity vector (length >= numDOF, exactly as ebkernel::addElementEnergy
    // builds it). Returns 0.0 if the element carries no scaling mass.
    double elementScalingKE(int eleTag, const Vector &vel, int numDOF) const;

    // ---- ADR-69 P3: lumped (DT2MS-style) SMS nodal injection ----
    // Publish/replace the active lumped-SMS integrator's per-node additive diagonal
    // ΔM (node tag -> Vector(ndf)), exactly the `injected` map it committed to the
    // Domain via applyMassScaling(+1). Same wholesale-replace + owner-guard
    // semantics as publish()/clear() above, on an independent owner token.
    void publishNodal(const void *owner, const std::map<int, Vector> &diag);
    void clearNodal(const void *owner);

    // Fast guard for the recorder's nodal loop (mirrors active()).
    bool activeNodal() const { return !theNodalDiag.empty(); }

    // True iff EITHER side holds blocks — the -v2 recorder's KE_ms column gate.
    bool anyActive() const { return active() || activeNodal(); }

    // 1/2 v^T ΔM v for node `nodeTag` (diagonal quadratic in the node's own
    // velocity). Returns 0.0 if the node carries no injected mass.
    double nodalScalingKE(int nodeTag, const Vector &vel) const;

private:
    MassScalingEnergyRegistry() : owner(0), nodalOwner(0) {}
    MassScalingEnergyRegistry(const MassScalingEnergyRegistry &);
    MassScalingEnergyRegistry &operator=(const MassScalingEnergyRegistry &);

    const void *owner;                 // current publisher (clear() is owner-guarded)
    std::map<int, Matrix> theBlocks;   // element tag -> node-major M_bar

    const void *nodalOwner;            // current lumped publisher (ADR-69 P3)
    std::map<int, Vector> theNodalDiag; // node tag -> injected diagonal ΔM
};

} // namespace Ladruno

#endif // LadrunoMassScalingEnergy_h
