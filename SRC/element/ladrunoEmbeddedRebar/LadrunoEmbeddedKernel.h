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

// Ladruno: shared coupling-kernel for the embedded-element family (ADR 23 §D1).
//
//   A COMPILED HELPER translation unit (LadrunoEmbeddedKernel.{h,cpp}) holding the
//   constitutive-agnostic numerics shared by `LadrunoEmbeddedRebar` (ELE 33005) and
//   the planned general `LadrunoEmbeddedNode` (ELE 33006). It is NOT a header-only
//   OpenSees-free leaf (cf. LadrunoJ2Kernel): its functions operate on framework
//   types (Matrix/Vector) and take framework handles (Node**) as explicit
//   parameters — the host-Element / Domain orchestration and the warnings stay in
//   each element (each owns its own messaging + ndf/DOF bookkeeping).
//
//   What lives here (the genuinely reusable, bit-reproducible core):
//     - maxAbsDiagonal            : the host-stiffness/mass scale reader (auto-kt, ω_host).
//     - effectiveCouplingStiffness: k_eff = max(k_axial0, k_t)  (the bounding spring).
//     - massPenaltyDtcr/Wcap      : the bipenalty m_p budgets (ADR 20 §10.6 / ADR 23 D5).
//     - criticalTimeStep          : the self-reported explicit step 2√(m_p/k_eff).
//     - computeGap                : g = u_c − Σ_i N_i u_host,i  (current-config rel. disp).
//     - assembleResistingForce    : P = Bᵀt  with B = [ I | −N_1 I … −N_M I ].
//     - assembleTangent           : K += Bᵀ D B  (c = [1, −N_1, …, −N_M], block c_p c_q D).
//
//   assembleResistingForce/assembleTangent are the TRANSLATIONAL blocks (the U tie);
//   the node element's UR/UP rows are layered on top in its own assembly (ADR 23 §3).
//
//   See Ladruno_implementation/23_ladruno_embedded_node_adr.md (§D1, §13 M3/KX-1) and
//   Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoEmbeddedKernel_h
#define LadrunoEmbeddedKernel_h

class Matrix;
class Vector;
class Node;

namespace LadrunoEmbedded {

// --- pure data math (framework-data-aware, handle-free) ---

// max_i |M(i,i)| — the host-element stiffness/mass scale (~E·lch / lumped mass).
// 0.0 for an empty matrix.
double maxAbsDiagonal(const Matrix& M);

// k_eff = max(kAxial0, kt): the stiffest coupled spring sets the bounded bipenalty
// mode. Inputs are assumed non-negative (callers pass |·|).
double effectiveCouplingStiffness(double kAxial0, double kt);

// bipenalty mass penalty (ADR 20 §10.6 / ADR 23 D5):
//   -dtcr dt : m_p = k_eff·(dt/2)²   (explicit budget; no host query).
double massPenaltyDtcr(double kEff, double dt);
//   -wcap β  : m_p = k_eff/(β²·ω_host²),  ω_host² = ‖K_host‖_∞/‖M_host‖_∞.
double massPenaltyWcap(double kEff, double beta, double omega2Host);

// self-reported explicit critical step 2√(m_p/k_eff); −1 ("no opinion") if either
// input is non-positive.
double criticalTimeStep(double mPenalty, double kEff);

// --- framework-handle-taking helpers ---

// gap g = u_constrained − Σ_i N_i u_host,i, read from the nodes' TRIAL displacements.
// `nodes` is [constrained, host_1 … host_M]; `N` holds the M host weights. Resizes
// and fills `g` (size ndm).
void computeGap(Node** nodes, const Vector& N, int nHost, int ndm, Vector& g);

// P = Bᵀt with B = [ I | −N_1 I … −N_M I ]: block 0 (constrained) = +t, block (1+i)
// (host i) = −N_i t. ASSIGNS the (1+M)·ndm touched entries (caller need not pre-zero).
// `t` is the ndm-vector traction; nHost = N.Size().
void assembleResistingForce(Vector& P, const Vector& N, const Vector& t, int ndm);

// K += Bᵀ D B with c = [1, −N_1, …, −N_M]: block (p,q) = c_p c_q D (ndm×ndm). ADDS
// into K (caller pre-zeroes); skips c_p c_q == 0 blocks. nHost = N.Size().
void assembleTangent(Matrix& K, const Vector& N, const Matrix& D, int ndm);

} // namespace LadrunoEmbedded

#endif // LadrunoEmbeddedKernel_h
