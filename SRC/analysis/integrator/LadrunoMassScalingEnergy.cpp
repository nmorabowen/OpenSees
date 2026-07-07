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

// Ladruno V4 — consistent (Olovsson) mass-scaling energy conduit. The single
// translation unit that owns the process-global block registry; the header-only
// energy kernel and the consistent-SMS integrators link against these symbols.
// See LadrunoMassScalingEnergy.h.

#include <LadrunoMassScalingEnergy.h>

namespace Ladruno {

MassScalingEnergyRegistry &
MassScalingEnergyRegistry::instance()
{
    static MassScalingEnergyRegistry theRegistry;   // thread-safe init (C++11+)
    return theRegistry;
}

void
MassScalingEnergyRegistry::publish(const void *o, const std::map<int, Matrix> &blocks)
{
    // Wholesale replace: domainChanged rebuilds the block set from scratch, so the
    // newest publish is authoritative for its owner.
    owner = o;
    theBlocks = blocks;
}

void
MassScalingEnergyRegistry::clear(const void *o)
{
    // Owner-guarded: a stale teardown of a superseded integrator must not wipe the
    // blocks a newer integrator has since published.
    if (o == owner) {
        theBlocks.clear();
        owner = 0;
    }
}

// ---- ADR-69 P3: lumped (DT2MS-style) SMS nodal injection ----

void
MassScalingEnergyRegistry::publishNodal(const void *o, const std::map<int, Vector> &diag)
{
    // Wholesale replace, mirroring publish(): every domainChanged re-baselines and
    // rebuilds `injected` from scratch, so the newest publish is authoritative.
    nodalOwner = o;
    theNodalDiag = diag;
}

void
MassScalingEnergyRegistry::clearNodal(const void *o)
{
    if (o == nodalOwner) {
        theNodalDiag.clear();
        nodalOwner = 0;
    }
}

double
MassScalingEnergyRegistry::nodalScalingKE(int nodeTag, const Vector &vel) const
{
    std::map<int, Vector>::const_iterator it = theNodalDiag.find(nodeTag);
    if (it == theNodalDiag.end())
        return 0.0;

    const Vector &dm = it->second;
    const int n = (dm.Size() < vel.Size()) ? dm.Size() : vel.Size();

    double ke = 0.0;
    for (int i = 0; i < n; ++i)
        ke += 0.5 * dm(i) * vel(i) * vel(i);
    return ke;
}

double
MassScalingEnergyRegistry::elementScalingKE(int eleTag, const Vector &vel, int numDOF) const
{
    std::map<int, Matrix>::const_iterator it = theBlocks.find(eleTag);
    if (it == theBlocks.end())
        return 0.0;

    const Matrix &Mbar = it->second;
    int n = Mbar.noRows();
    if (n > numDOF)     n = numDOF;       // defensive: never read past the element DOFs
    if (n > vel.Size()) n = vel.Size();   //            nor past the supplied velocity

    double ke = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            ke += 0.5 * Mbar(i, j) * vel(i) * vel(j);
    return ke;
}

} // namespace Ladruno
