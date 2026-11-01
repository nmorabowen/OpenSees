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

// Ladruno: discrete cohesive moment–rotation law for the LadrunoDispBeamColumn
// Tier-2 embedded strong-discontinuity hinge (ADR 32 §Tier-2). A UniaxialMaterial
// whose "strain" is the rotation jump [[theta]] (kappa) and whose "stress" is the
// cohesive moment M([[theta]]). The hinge carries a *fracture energy* Gf per hinge
// directly (units of energy = moment * rotation) — there is NO characteristic
// length to calibrate, which is the whole point of the discrete-hinge formulation
// vs. the Tier-1 lch crack-band smearing.
//
// Law (rigid–softening cohesive with irreversible secant-to-origin damage):
//
//   pre-peak  (kappa <= kappa0):  M = Kpen * kappa            (near-rigid penalty:
//                                  the hinge stays CLOSED until M reaches Mc)
//   softening (kappa >  kappa0):  exponential  M = Mc * exp(-a (kappa-kappa0))
//                                  linear       M = Mc (1 - (kappa-kappa0)/kappaf)
//
// The softening branch is calibrated so the TOTAL area under the monotonic
// envelope equals Gf EXACTLY:
//   kappa0 = Mc/Kpen                       (jump at activation)
//   pre-peak (recoverable) area = Mc^2/(2 Kpen)
//   exp :  a      = Mc / (Gf - Mc^2/(2 Kpen))   ->  total = Mc^2/(2Kpen) + Mc/a = Gf
//   lin :  kappaf = 2 Gf/Mc - Mc/Kpen           ->  total = Mc^2/(2Kpen) + Mc kappaf/2 = Gf
//
// GUARDED PENALTY FLOOR: both calibrations require Gf > Mc^2/(2 Kpen), i.e.
//   Kpen > Kpen_floor = Mc^2 / (2 Gf).
// Below the floor the pre-peak elastic energy alone would exceed Gf and the
// softening exponent/length goes non-physical. The default penalty is
//   Kpen = penaltyRatio * Kpen_floor   (penaltyRatio default 1000)
// so the recoverable pre-peak energy is 1/penaltyRatio of Gf (≈0.1%): near-rigid,
// floor-safe by construction. A user -penalty K below the floor is clamped with a
// loud opserr warning.
//
// Irreversibility: the max |kappa| (kappaMax) is a monotone history variable; on
// unload/reload the response is the secant Ksec = Menv(kappaMax)/kappaMax to the
// origin (isotropic damage), so the hinge never recovers strength and never stores
// residual moment at zero jump.
//
// classTag MAT_TAG_LadrunoCohesiveHinge = 33003 (Ladruno private uniaxial band).
// See Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoCohesiveHinge_h
#define LadrunoCohesiveHinge_h

#include <UniaxialMaterial.h>

class LadrunoCohesiveHinge : public UniaxialMaterial
{
 public:
  enum SofteningShape { EXP = 0, LINEAR = 1 };

  LadrunoCohesiveHinge(int tag, double Mc, double Gf,
                       double Kpen = -1.0,            // <0 -> default penaltyRatio*floor
                       int shape = EXP,
                       double penaltyRatio = 1000.0,
                       double P0 = 0.0,
                       double alphaPM = 0.5,
                       bool implex = false,
                       double P0t = 0.0,
                       double Pbal = 0.0,
                       double muBal = 1.0);
  LadrunoCohesiveHinge();
  ~LadrunoCohesiveHinge();

  const char* getClassType(void) const { return "LadrunoCohesiveHinge"; }

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void)         { return Tstrain; }
  double getStress(void)         { return Tstress; }
  double getTangent(void)        { return Ttangent; }
  double getInitialTangent(void) { return (P0 > 0.0) ? Kpen * getPMMultiplier() : Kpen; }   // near-rigid: hinge starts closed
  double getEnergy(void)         { return Twork; }   // path work int M d[[theta]]

  void setAxialForce(double P)   { P_axial = P; }
  double getAxialForce(void) const { return P_axial; }
  double getPMMultiplier(void) const;
  double getPMDerivative(void) const;
  void setPM(double p0, double apm) { P0 = p0; alphaPM = apm; P0t = 0.0; Pbal = 0.0; }
  void setPMBal(double p0c, double p0t, double pbal, double mubal) {
    P0 = p0c; P0t = p0t; Pbal = pbal; muBal = mubal;
  }
  bool hasPM(void) const         { return P0 > 0.0; }
  bool isImplex(void) const      { return implex; }
  void setImplex(bool flag)      { implex = flag; }
  double getImplexError(void) const;

  bool ladrunoThreadSafeUpdate(void) const { return true; }
  bool ladrunoThreadSafeForm(void) const   { return true; }

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial* getCopy(void);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);

 private:
  // calibrate kappa0 / aExp / kappaf from (Mc, Gf, Kpen, shape); enforce the floor.
  void computeDerived(void);
  // envelope moment & consistent slope at a non-negative max-jump km
  void envelope(double km, double& Menv, double& dMenv) const;

  // material parameters
  double Mc;            // cohesive moment capacity (peak)
  double Gf;            // fracture energy per hinge (energy units)
  double Kpen;          // near-rigid penalty stiffness
  bool   userKpen;      // true if user explicitly specified -penalty / Kpen > 0
  int    shape;         // EXP or LINEAR
  double penaltyRatio;  // default-penalty multiple of the floor (record for getCopy)
  double P0;            // axial compressive capacity for P-M interaction (<=0 disables)
  double alphaPM;       // parabolic P-M interaction coefficient
  double P0t;           // tensile axial capacity for -pmBal (>0 enables 4-parameter RC balance law)
  double Pbal;          // compressive balance axial load for -pmBal (0 < Pbal < P0)
  double muBal;         // balance moment ratio M_bal / M_0 for -pmBal
  bool   implex;        // Oliver et al. (2008) IMPL-EX integration flag

  // transient axial force supplied by parent beam element
  double P_axial;

  // derived (rebuilt by computeDerived, never serialized)
  double kappa0;        // Mc/Kpen   -- jump at activation
  double aExp;          // exponential softening rate
  double kappaf;        // linear softening jump span (Mc -> 0)

  // committed history
  double CkappaMax;     // max |kappa| ever seen (irreversible damage driver)
  double CkappaMax_old; // previous step committed max |kappa| (for IMPL-EX extrapolation)
  double Cstrain;       // committed jump   (work-integral anchor)
  double Cstress;       // committed moment (work-integral anchor)
  double Cwork;         // committed path work int M d[[theta]]
  double dtime_n_commit;
  int    implexCommittedSteps;

  // trial state
  double Tstrain;
  double Tstress;
  double Ttangent;
  double TkappaMax;
  double Twork;
};

#endif
