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

// Ladruno: 1D bond-slip tau-s UniaxialMaterial (CEB-FIP Model Code 2010 backbone).
//
//   Relates local bond stress tau to bar-concrete slip s for the bar-axial slot
//   of the embedded-reinforcement element (LadrunoEmbeddedRebar, ADR 20 / D4).
//   v1 = MONOTONIC backbone + elastic unload/reload (cyclic degradation deferred
//   to v2). Sign-symmetric: tau(-s) = -tau(s).
//
//   Backbone (MC2010), with the ADR-D4 must-fix regularizations baked in:
//     0    <= |s| < s0 : LINEAR  tau = k0*|s|         (D4.1 initial-slip reg.:
//                        kills the power-law dtau/ds -> inf singularity at s->0;
//                        s0 defaults to 0.1*s1, k0 = tau_max*(s0/s1)^alpha / s0)
//     s0   <= |s| < s1 : tau = tau_max*(|s|/s1)^alpha   (ascending power law)
//     s1   <= |s| < s2 : tau = tau_max                  (plateau)
//     s2   <= |s| < s3 : tau = tau_max - (tau_max-tau_f)*(|s|-s2)/(s3-s2)
//                        (softening, NEGATIVE tangent -> caller must use
//                        displacement/arc-length/IMPLEX control, D4.2)
//     |s|  >= s3       : tau = tau_f                    (residual friction)
//
//   D4.3 fracture-energy regularization: if Gf > 0 it OVERRIDES s3 so the
//   softening branch dissipates the bond fracture energy per unit interface area
//   Gf = 0.5*(tau_max - tau_f)*(s3 - s2)  =>  s3 = s2 + 2*Gf/(tau_max - tau_f).
//   (The element scales tau by perimeter*L_trib for nodal force; feeding Gf here
//   fixes the per-area dissipation so the embedded response is mesh-objective.)
//
//   Unload/reload: elastic with the initial stiffness k0, clamped to the
//   envelope at the committed peak |slip|. Tangent is finite everywhere (k0 off
//   the envelope; the branch slope on it). getInitialTangent = k0.
//
//   classTag MAT_TAG_LadrunoBondSlip = 33002 (Ladruno private uniaxial band,
//   after LadrunoUniaxialJ2=33000, LadrunoRebarBuckling=33001).
//   See Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md (D4).
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoBondSlip_h
#define LadrunoBondSlip_h

#include <UniaxialMaterial.h>

class LadrunoBondSlip : public UniaxialMaterial
{
 public:
  LadrunoBondSlip(int tag, double tau_max, double s1, double s2, double s3,
                  double tau_f, double alpha, double Gf, double s0);
  LadrunoBondSlip();
  ~LadrunoBondSlip();

  const char* getClassType(void) const { return "LadrunoBondSlip"; }

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void)         { return Tslip; }
  double getStress(void)         { return Tstress; }
  double getTangent(void)        { return Ttangent; }
  double getInitialTangent(void) { return k0; }

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial* getCopy(void);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  // --- backbone parameters (s3 may be overridden by Gf) ---
  double tau_max, s1, s2, s3, tau_f, alpha, Gf, s0;
  double k0;                 // initial (linear-segment) stiffness = envelope'(0+)

  void setDerived(void);     // compute s0 default, k0, and Gf-override of s3
  double envelope(double a) const;      // backbone for a = |slip| >= 0
  double envelopeTangent(double a) const;

  // --- committed state ---
  double Cslip;              // committed slip
  double Cstress;            // committed bond stress
  double CslipMaxAbs;        // max |slip| reached on the envelope (damage memory)

  // --- trial state ---
  double Tslip;
  double Tstress;
  double Ttangent;
  double TslipMaxAbs;
};

#endif
