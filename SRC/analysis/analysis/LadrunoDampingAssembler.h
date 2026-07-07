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

// Ladruno: reduced-operator projection for complex modal analysis (ADR 46 P2).
// See Ladruno_implementation/46_ladruno_complex_modal_adr.md §4.6 Route B.
//
//   Projects the model's damping and mass operators onto the real undamped
//   modal basis Phi from the last `eigen` run WITHOUT ever assembling a
//   global matrix: accumulate element-by-element
//
//       Ct~ += Phi_e^T C_e Phi_e     (C_e = Element::getDamp() — Rayleigh
//                                     per-element factors + material/damper
//                                     overrides, exactly the transient's C)
//       Mt~ += Phi_e^T M_e Phi_e     (M_e = Element::getMass())
//
//   plus the nodal contributions (lumped nodal mass and its alphaM Rayleigh
//   term via Node::getDamp()/getMass()). Phi_e is gathered from the nodal
//   eigenvector storage (Node::getEigenvectors) in the element's node/DOF
//   ordering — no DOF numbering, no SOE, no integrator.
//
//   The stiffness projection is NOT assembled: K phi_a = w_a^2 M phi_a gives
//   Kt~ = Mt~ * diag(w^2) EXACTLY (symmetrized to kill roundoff), for any
//   basis normalization and even for repeated eigenvalues (vectors of
//   distinct eigenvalues are automatically M-orthogonal; within a repeated
//   eigenspace Lambda is a multiple of I, so Mt~ and Lambda commute).
//   Passing the FULL projected Mt~ (not an assumed identity) to the QZ
//   kernel is what makes the repeated-eigenvalue case exact without any
//   re-orthonormalization — see LEDGER_quirks (ADR 46 P2 pre-warning).
//
//   Route B strictly contains Route A: for a pure global-Rayleigh model the
//   element/node getDamp() sums reproduce the closed form (the P2 test gate
//   pins B == A at 1e-10), and scoped (region/element) Rayleigh plus
//   betaKinit/betaKcomm terms are captured exactly because getDamp() already
//   encodes them per element.
//
//   CONTRACT / STALENESS (read this before trusting zeta):
//   - C is exactly what getDamp() returns: elements whose damping does NOT
//     flow through getDamp() are invisible — modalDamping (integrator-side;
//     the driver warns), element 'damping' framework objects (folded into
//     getResistingForceIncInertia), HHT-alpha numerical damping, and any
//     element with its -doRayleigh flag left at its default 0 (Truss and the
//     zeroLength family DEFAULT OFF — see LEDGER_quirks).
//   - M~ and C~ are read from the LIVE model while Lambda is the SNAPSHOT
//     spectrum of the last eigen run. Editing the model after eigen (mass,
//     stiffness, state) makes them inconsistent — and because Kt = M~ Lambda
//     ECHOES the snapshot, the reported omega0 still equals the stale
//     frequencies, masking the staleness in exactly the field a user checks.
//     Re-run eigen after any model change (same convention as
//     modalProperties).
//   - Eigen must have run under a constraint handler that distributes
//     eigenvectors to constrained/slave DOFs (the default Transformation
//     handler does; 'constraints Plain' + MP constraints writes ZERO slave
//     rows and would misproject).

#ifndef LadrunoDampingAssembler_h
#define LadrunoDampingAssembler_h

class Domain;
class Matrix;

class LadrunoDampingAssembler
{
  public:
    // Project the model's M and C onto the first p modes of the last eigen
    // run and synthesize Kt from the eigenvalues (see header note). Mt, Ct,
    // Kt must be p x p (they are zeroed here). Returns 0 on success, < 0 on
    // error (no/short spectrum, a node with fewer stored modes than p, or an
    // element whose damping/mass dimension does not match its nodes' DOFs).
    static int projectReduced(Domain &theDomain, int p,
                              Matrix &Mt, Matrix &Ct, Matrix &Kt);
};

#endif
