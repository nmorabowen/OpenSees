/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoUPKernel — the PURE numerical core of the shared Biot u-p (solid
// displacement + pore pressure) mixed-continuum assembly (ADR 71). Plain doubles
// + <math.h>, NO OpenSees dependency, so it is unit-tested standalone against the
// numpy oracle (tests/ladruno_up_reference.py) WITHOUT building OpenSees, and
// reused by every geometry (T3 / Q4 / H8 / Bézier-T6 / Bézier-Tet10 …) — each
// supplies only its shape values + gradients + GP rule (LadrunoUPShapes.h); the
// coupled Biot block mechanics live here, written once (the ADR-70 idiom, mirror
// of SRC/element/ladrunoPlane/LadrunoFiniteStrain2DKernel.h).
//
// ----------------------------------------------------------------------------
// LAYOUTS (row-major, caller-owned, caller-ZEROED before the GP loop):
//
//   Nu[a]              a = 0 … nNu-1   solid (displacement) shape values
//   dNu[a*ndm + i]     ∂N_a/∂x_i       solid cartesian gradient (node-major)
//   Np[b]              b = 0 … nNp-1   pressure shape values (COMPACTED over
//                                      carriers for Taylor–Hood: nNp = #carriers)
//   dNp[b*ndm + i]     ∂N^p_b/∂x_i     pressure cartesian gradient (node-major)
//   dv                 w · detJ · thick(2D)  |  w · detJ (3D)   folded by CALLER
//   kbar[i]            i = 0 … ndm-1   diagonal Darcy permeability k_hydraulic/γ_w
//   drive[i]           i = 0 … ndm-1   seepage drive: b (static/body) OR
//                                      b − ü_trial (dynamic seepage), built by CALLER
//
// Every add*() ACCUMULATES into a caller-owned, caller-zeroed buffer; ld* is the
// leading dimension (row stride) of the destination matrix. The kernel neither
// allocates nor validates: geometry guards (detJ > 0, non-degenerate Jacobian)
// are the PROVIDER's / ELEMENT's job — this header stays OpenSees-free and silent
// on failure, exactly like LadrunoFiniteStrain2DKernel.h. A caller that feeds a
// singular Jacobian's gradients inherits a silent-garbage trap; guard upstream.
//
// Block shapes (all in the effective-stress / honest-p contract, ADR §3.2):
//   Q   (coupling)      : (nNu·ndm) rows × nNp cols   Q[(a·ndm+i), b]
//   H   (Darcy)         : nNp × nNp                   H[b, c]
//   S   (storage)       : nNp × nNp                   S[b, c]
//   Htilde(stabilizn.)  : nNp × nNp                   Ht[b, c]  (equal-order only)
//   fseep(seepage src.) : nNp                         fseep[b]
//
// References: Zienkiewicz & Shiomi (1984); Zienkiewicz et al., "Computational
// Geomechanics" (1999) §2–3, eqs 2.17 / 3.66; McGann–Arduino–Mackenzie-Helnwein
// (2012 eq 73, 2015 eq 74) for the direct-α stabilization. Block definitions and
// the honest-p arrangement are ADR-71 §3.1–3.3.

#ifndef LadrunoUPKernel_h
#define LadrunoUPKernel_h

#include <math.h>

namespace ladruno_up {

// --------------------------------------------------------------------------- //
//  Q — volumetric coupling block  ∫ Bᵀ α m N_p dV                             //
//                                                                             //
//  Q[(a·ndm+i), b] += dNu[a·ndm+i] · biotAlpha · Np[b] · dv                   //
//                                                                             //
//  Rows are per-node-interleaved solid DOFs (node a, component i); cols are   //
//  the (compacted) pressure carriers. ldQ = nNp.                             //
// --------------------------------------------------------------------------- //
inline void addQ(const double* dNu, const double* Np, int nNu, int nNp, int ndm,
                 double biotAlpha, double dv, double* Q, int ldQ) {
  for (int a = 0; a < nNu; a++)
    for (int i = 0; i < ndm; i++) {
      int row = a * ndm + i;
      double g = dNu[row] * biotAlpha * dv;
      for (int b = 0; b < nNp; b++)
        Q[row * ldQ + b] += g * Np[b];
    }
}

// --------------------------------------------------------------------------- //
//  H — Darcy permeability (Laplacian) block  ∫ (∇N_p)ᵀ k̄ ∇N_p dV  (diag k̄)   //
//                                                                             //
//  H[b, c] += Σ_i dNp[b·ndm+i] · kbar[i] · dNp[c·ndm+i] · dv                  //
//  Symmetric; ldH = nNp.                                                      //
// --------------------------------------------------------------------------- //
inline void addH(const double* dNp, int nNp, int ndm, const double* kbar,
                 double dv, double* H, int ldH) {
  for (int b = 0; b < nNp; b++)
    for (int c = 0; c < nNp; c++) {
      double s = 0.0;
      for (int i = 0; i < ndm; i++)
        s += dNp[b * ndm + i] * kbar[i] * dNp[c * ndm + i];
      H[b * ldH + c] += s * dv;
    }
}

// --------------------------------------------------------------------------- //
//  S — storage / compressibility block  ∫ N_pᵀ (1/Q̄) N_p dV                   //
//                                                                             //
//  S[b, c] += Np[b] · oneOverQbar · Np[c] · dv   (symmetric; ldS = nNp)       //
// --------------------------------------------------------------------------- //
inline void addS(const double* Np, int nNp, double oneOverQbar, double dv,
                 double* S, int ldS) {
  for (int b = 0; b < nNp; b++) {
    double g = Np[b] * oneOverQbar * dv;
    for (int c = 0; c < nNp; c++)
      S[b * ldS + c] += g * Np[c];
  }
}

// --------------------------------------------------------------------------- //
//  Htilde — pressure-Laplacian stabilization block (equal-order pairs)        //
//           ∫ (∇N_p)ᵀ α_stab ∇N_p dV   (isotropic scalar α_stab)              //
//                                                                             //
//  Same form as addH with a scalar in place of the diagonal k̄. Lands in the  //
//  damp p-p slot (S* = S + H̃, ADR §3.3). Symmetric; ldHt = nNp.              //
// --------------------------------------------------------------------------- //
inline void addHtilde(const double* dNp, int nNp, int ndm, double stabAlpha,
                      double dv, double* Ht, int ldHt) {
  for (int b = 0; b < nNp; b++)
    for (int c = 0; c < nNp; c++) {
      double s = 0.0;
      for (int i = 0; i < ndm; i++)
        s += dNp[b * ndm + i] * dNp[c * ndm + i];
      Ht[b * ldHt + c] += stabAlpha * s * dv;
    }
}

// --------------------------------------------------------------------------- //
//  fseep — seepage source vector  ∫ (∇N_p)ᵀ k̄ ρ_f drive dV                    //
//                                                                             //
//  fseep[b] += Σ_i dNp[b·ndm+i] · kbar[i] · rhoF · drive[i] · dv              //
//  drive = b (static/body) or b − ü_trial (dynamic seepage), built by CALLER. //
// --------------------------------------------------------------------------- //
inline void addFseep(const double* dNp, int nNp, int ndm, const double* kbar,
                     double rhoF, const double* drive, double dv, double* fseep) {
  for (int b = 0; b < nNp; b++) {
    double s = 0.0;
    for (int i = 0; i < ndm; i++)
      s += dNp[b * ndm + i] * kbar[i] * rhoF * drive[i];
    fseep[b] += s * dv;
  }
}

// --------------------------------------------------------------------------- //
//  Auto stabilization coefficient (McGann et al. 2012 eq 73 / 2015 eq 74):    //
//    α = α₀ · h² / (K_s + 4·G_s/3)                                            //
//  h = largest element dimension = longest EDGE (Provider::elemSizeH, 0.E-F1).//
//  KsDrained, GsDrained = DRAINED SKELETON bulk & shear moduli from the       //
//  material's initial isotropic tangent — NOT the grain bulk K_s that         //
//  oneOverQbar() takes below (name-collision trap: same letter in the papers, //
//  two different physical moduli).                                            //
//  PRECONDITION: KsDrained + 4·GsDrained/3 > 0. A ≤0 denominator returns inf  //
//  or a NEGATIVE α — anti-stabilization that silently DEstabilizes addHtilde. //
//  Guarding is the ELEMENT's job (P1 parser/setup); the kernel stays silent.  //
// --------------------------------------------------------------------------- //
inline double stabAlphaAuto(double h, double KsDrained, double GsDrained,
                            double alpha0) {
  return alpha0 * h * h / (KsDrained + 4.0 * GsDrained / 3.0);
}

// --------------------------------------------------------------------------- //
//  Storage coefficient 1/Q̄ = n/K_f + (α − n)/K_s   (ADR §3.1, Book eq 2.17).  //
//  Ks here is the GRAIN bulk modulus (cf. the drained-skeleton KsDrained in   //
//  stabAlphaAuto above). Ks <= 0 encodes incompressible grains (K_s → ∞):     //
//  the second term vanishes.                                                  //
//  Expected domain (policed by the P1 parser, not here): 0 < n <= biotAlpha   //
//  <= 1 and Kf > 0 — outside it the formula returns a physically meaningless  //
//  (possibly negative) storage coefficient.                                   //
// --------------------------------------------------------------------------- //
inline double oneOverQbar(double n, double Kf, double biotAlpha, double Ks) {
  double val = n / Kf;
  if (Ks > 0.0)
    val += (biotAlpha - n) / Ks;
  return val;
}

// --------------------------------------------------------------------------- //
//  Per-node-interleaved DOF map for one element [u…, p?] with p last on each   //
//  carrier node (ADR §3.5 "keep"). pCarrier[n] ∈ {0,1}. Fills:                //
//    uOff[n] = first u-DOF slot of node n   (always ndm slots)                //
//    pOff[n] = the node's p slot, or −1 if it carries no pressure             //
//  Returns the total DOF count. Caller sizes uOff/pOff to nNodes.             //
// --------------------------------------------------------------------------- //
inline int dofMap(int nNodes, int ndm, const int* pCarrier, int* uOff, int* pOff) {
  int slot = 0;
  for (int n = 0; n < nNodes; n++) {
    uOff[n] = slot;
    slot += ndm;
    if (pCarrier[n]) {
      pOff[n] = slot;
      slot += 1;
    } else {
      pOff[n] = -1;
    }
  }
  return slot;
}

} // namespace ladruno_up

#endif
