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

// ADR-62 — LadrunoTie: a setup-time KINEMATIC mesh-tie GENERATOR (no enforcement
// code).  A non-conforming explicit mesh-tie is a permanent, non-sliding bond, so
// the pairing is frozen once at model-build: each slave node is projected onto the
// nearest master facet (closest-point projection), the master facet shape weights
// N_i(xi_bar) are evaluated at the projection, and one ordinary OpenSees
// EQ_Constraint per slave node per translational DOF is emitted:
//
//        u_s = sum_i N_i u_{m,i}            (1*u_s + sum_i (-N_i) u_{m,i} = 0)
//
// The user then runs `constraints LadrunoProjection` (ADR-30, the shipped explicit
// handler) — or any EQ-aware handler (Lagrange/Penalty for implicit/static) — to
// enforce the tie.  Enforcement is EXACT, dt_cr-neutral, momentum-conserving, and
// adds no fictitious mass: the constructive successor to the shelved penalty route
// (ADR-61).  See Ladruno_implementation/62_ladruno_kinematic_mesh_tie_adr.md.
//
// Reused, unchanged: LadrunoContactProjection (ADR-39/41 closest-point projection +
// shape weights) and LadrunoContactBucketSort (ADR-39 broad phase).  Scope (P1):
// solid–solid, node-to-segment collocation only (nps = 3 tri / 4 quad); shells and
// integral-mortar are deferred.

#ifndef LadrunoTie_h
#define LadrunoTie_h

class Domain;
class ID;

namespace LadrunoTie {

// Pair each slave node to its nearest master facet, evaluate the master shape
// weights N_i at the projection, and emit one EQ_Constraint per slave node per
// tied DOF (u_s = sum_i N_i u_{m,i}).
//
//   slaveNodes        : the slave-surface node tags (each tied to ONE facet).
//   masterFacetNodes  : flat node tags, nps per facet, row-major [facet][node].
//   nps               : nodes per master facet (3 tri-3 / 4 quad-4).
//   dofs              : 1-based translational DOFs to tie (e.g. {1,2,3} for a solid).
//   tolFrac           : conforming tolerance — refuse if a slave projects farther
//                       than tolFrac * facet-size off the master surface (OQ-3:
//                       require conforming-at-interface geometry, no IC snapping).
//
// Returns the number of EQ_Constraints emitted (>= 0) on success, or -1 on a named
// refusal: BLOCKER-1 (slave/master not node-disjoint, or a slave listed twice =>
// MP-chain / double-constraint the projection handler refuses), BLOCKER-2 (a tied
// slave DOF carries no lumped mass), OQ-3 (slave off the master manifold), or a
// geometry/lookup error.
int generate(Domain *theDomain,
             const ID &slaveNodes,
             const ID &masterFacetNodes, int nps,
             const ID &dofs, double tolFrac);

// P2 — integral-mortar variant. Instead of point collocation, tie the slave
// SURFACE to the master surface in the weak (integral) sense: assemble the global
// mortar operators over the clipped overlap
//        D_IJ = INT N_I^s N_J^s dGamma   (slave-slave, interface consistent mass)
//        M_IK = INT N_I^s phi_K^m dGamma (slave-master)
// (reusing LadrunoMortarKernel::integratePair verbatim), CONDENSE once at setup to
//        u_s = P u_m ,  P = D^{-1} M
// and emit one EQ_Constraint per slave node per tied DOF (u_s = sum_k P_sk u_{m,k}).
// Pre-inverting D globally makes every P row tie to MASTER nodes only => no MP-chains
// => the shipped projection handler enforces it unchanged (variationally consistent,
// optimal for non-matching interfaces).  Both surfaces are FACETED (tri-3/quad-4).
//
//   slaveFacetNodes  : flat tags, npsS per slave facet, row-major [facet][node].
//   masterFacetNodes : flat tags, npsM per master facet.
//   outward          : interface outward direction (orients the mortar normal); if 0,
//                      the average master-facet normal is used.
//
// Returns #EQ_Constraints emitted (>= 0), or -1 on a named refusal: surfaces not
// node-disjoint, an uncovered slave node (=> singular D), a non-conforming gap
// (surfaces not coincident, OQ-3), or a massless tied DOF.
int generateMortar(Domain *theDomain,
                   const ID &slaveFacetNodes, int npsS,
                   const ID &masterFacetNodes, int npsM,
                   const ID &dofs, double tolFrac, const double *outward);

} // namespace LadrunoTie

#endif
