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

// ADR-39 P2b — LadrunoContactKernel: header-only, OpenSees-free pure-function NTS
// penalty NORMAL-LAW + COULOMB FRICTION math (traction, return map, consistent tangent).
// The closest-point PROJECTION geometry was EXTRACTED to LadrunoContactProjection.h
// (ADR-41 A2, capstone Track A2) — ONE projection for both the NTS path and the mortar
// Gauss-point loop. Mirrors LadrunoJ2Kernel.h so both OPS_Analysis (the FE adapter) and
// OPS_Domain (the surface) can include it with no link inversion. NO OpenSees types —
// only raw double[] + <cmath>.
//
// Validated oracle-first in Ladruno_implementation/contact_prototypes/proto_p2b_nts.py
// (7/7: projection, winding-flip-immune normal, penalty, B-operator, self-equilibrium,
// symmetric kn*BᵀB + ∂n/∂u tangent vs FD). 3D, master segment = tri-3 or quad-4.
//
// Conventions (match shipped P2a, verified): penetration ⇔ gap < 0; traction
// tn = kn*<−gap>₊ (Macaulay); getResidual returns +tn over the gap operator B so a
// penetrating slave is pushed back along +n (restoring, never attract). The assembled
// tangent is K = −∂r/∂q = +kn*BᵀB (main term; ∂n/∂u block is P2b-2). Outward normal is
// DERIVED + oriented from a reference point, NEVER trusted from node winding
// (design-gate BLOCKER-1).
//
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p2b_design.md.

#ifndef LadrunoContactKernel_h
#define LadrunoContactKernel_h

#include <cmath>

namespace LadrunoContactKernel {

// Geometry → LadrunoContactProjection.h (ADR-41 A2); Coulomb/Tresca friction →
// LadrunoFrictionKernel.h (ADR-41 A1). This header now owns ONLY the normal penalty law.

// penalty traction tn = kn*<−gap>₊ . Internal Macaulay clamp (gate KMF-3): returns
// 0 for gap>=0 so a stray call can never produce a non-physical (adhesive) traction.
inline double traction(double kn, double gap) { return (gap < 0.0) ? kn * (-gap) : 0.0; }

} // namespace LadrunoContactKernel

#endif
