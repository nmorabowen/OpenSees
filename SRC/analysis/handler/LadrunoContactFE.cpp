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

// ADR-39 P1a — LadrunoContactFE implementation (empty-connectivity zero adapter).
// See LadrunoContactFE.h for the design rationale.

#include "LadrunoContactFE.h"

LadrunoContactFE::LadrunoContactFE(int tag)
  : FE_Element(tag, /*numDOF_Group=*/0, /*ndof=*/0),
    resid(0), tang(0, 0)
{
    // myDOF_Groups and myID are size 0 (empty connectivity): the adapter adds NO
    // edges to the DOF graph, so the numberer permutation is untouched and the
    // result is bitwise-identical to no-contact. P2 declares real per-segment
    // connectivity (and switches the gate to 1e-12).
}

LadrunoContactFE::~LadrunoContactFE()
{
}

const Vector &
LadrunoContactFE::getResidual(Integrator *)
{
    // P1a: zero force. (P2: project pairs from this adapter's connected nodes'
    // trial displacements and assemble F_c here.)
    resid.Zero();
    return resid;
}

const Matrix &
LadrunoContactFE::getTangent(Integrator *)
{
    // P1a: zero. Implicit assembly does addKtToTang(c1), so P2's real tangent
    // must be returned already scaled by the integrator's c1 factor.
    tang.Zero();
    return tang;
}

void
LadrunoContactFE::addMtoTang(double)
{
    // contact pairs contribute no mass (mass lives on the real structural
    // elements at the same nodes); no-op keeps the explicit LHS = lumped M only.
}

// Size-0 returns for the off-hot-path force variants (base would warn on myEle==0).
const Vector &LadrunoContactFE::getTangForce(const Vector &, double) { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getK_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getKi_Force(const Vector &, double)  { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getC_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getM_Force(const Vector &, double)   { resid.Zero(); return resid; }
