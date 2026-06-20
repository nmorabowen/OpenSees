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

// Ladruno: COUPLED biaxial cohesive moment-rotation law for the
// LadrunoDispBeamColumn3d Tier-2 embedded hinge (ADR 34). An NDMaterial of
// ORDER 2: strain = the rotation-jump vector [[theta]] = [alpha_z, alpha_y],
// stress = the cohesive moment vector [M_z, M_y], tangent = the 2x2 dM/dalpha.
// A discrete cohesive law with an Mz-My INTERACTION SURFACE (ellipse) carrying a
// fracture energy per bending axis (Gf_z, Gf_y) directly — no characteristic
// length. The coupled generalization of the scalar LadrunoCohesiveHinge (33003):
// it reduces EXACTLY to that scalar law on each pure bending axis.
//
// Formulation (mixed-mode rigid-softening + isotropic secant damage):
//   per-axis penalty Kpen_i = penaltyRatio*Mc_i^2/(2 Gf_i) (floored at Mc_i^2/2Gf_i),
//   activation jump alpha0_i = Mc_i/Kpen_i, normalized jump e_i = alpha_i/alpha0_i.
//   ELLIPTICAL norm r = sqrt(e_z^2 + e_y^2): r<=1 -> CLOSED (M_i = Kpen_i alpha_i,
//   onset ellipse sqrt((Mz/Mcz)^2+(My/Mcy)^2)=1); r>1 -> damage. History r_max.
//   Mode mix w_z = e_z^2/r^2, w_y = e_y^2/r^2; S = w_z Mc_z^2/Kpen_z + w_y Mc_y^2/Kpen_y
//   (peak-scale, LINEAR -> exact pure-axis reduction). Fracture energy by Benzeggagh-Kenane:
//   Gf_mix = Gf_z + (Gf_y - Gf_z) w_y^eta (eta=1 -> linear; eta!=1 makes D mix-dependent and
//   the 2x2 tangent's off-radial mode-mix term live). Effective 1D law T_eff(r) (T_eff*dr = M*dalpha):
//   pre-peak T_eff = S r; softening calibrated to int T_eff dr = Gf_mix
//     EXP    T_eff = S exp(-A(r-1)),  A = S/(Gf_mix - S/2)
//     LINEAR T_eff = S(1-(r-1)/r_f),  r_f = 2(Gf_mix - S/2)/S
//   Isotropic secant damage D = 1 - T_eff(r_max)/(S r_max); M_i = (1-D) Kpen_i alpha_i.
//   Reduces to LadrunoCohesiveHinge(Mc_i, Gf_i) on each pure axis (peak, shape, Gf).
//
// classTag ND_TAG_LadrunoCohesiveHingeBiaxial = 33004 (Ladruno private ND band).
// See Ladruno_implementation/34_ladruno_cohesive_hinge_biaxial_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoCohesiveHingeBiaxial_h
#define LadrunoCohesiveHingeBiaxial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class LadrunoCohesiveHingeBiaxial : public NDMaterial
{
 public:
  enum SofteningShape { EXP = 0, LINEAR = 1 };

  LadrunoCohesiveHingeBiaxial(int tag, double Mcz, double Gfz, double Mcy, double Gfy,
                              int shape = EXP, double penaltyRatio = 1000.0,
                              double bkEta = 1.0);
  LadrunoCohesiveHingeBiaxial();
  ~LadrunoCohesiveHingeBiaxial();

  const char* getClassType(void) const { return "LadrunoCohesiveHingeBiaxial"; }

  int setTrialStrain(const Vector& v);
  int setTrialStrain(const Vector& v, const Vector& r);
  const Vector& getStress(void)        { return Tstress; }
  const Vector& getStrain(void)        { return Tstrain; }
  const Matrix& getTangent(void)       { return Ttangent; }
  const Matrix& getInitialTangent(void){ return Kinit; }

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial* getCopy(void);
  NDMaterial* getCopy(const char* type);
  const char* getType(void) const  { return "BeamHingeBiaxial"; }
  int getOrder(void) const         { return 2; }

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

  double getEnergy(void) { return Twork; }   // path work int M . dalpha

 private:
  void computeDerived(void);   // per-axis Kpen/alpha0 from (Mc_i, Gf_i), enforce the floor
  // effective 1D law at the current mix (w_z,w_y) and history r_max:
  //   returns S, Gf_mix, T_eff(r_max), dT_eff/dr(r_max)
  void effectiveLaw(double wz, double wy, double rMax,
                    double& S, double& Teff, double& dTeff) const;

  // material parameters
  double Mcz, Gfz, Mcy, Gfy;
  int    shape;
  double penaltyRatio;
  double bkEta;            // Benzeggagh-Kenane mode-mix exponent (1 = linear Gf interpolation)

  // derived (rebuilt by computeDerived, never serialized)
  double Kpenz, Kpeny;     // per-axis near-rigid penalties
  double a0z, a0y;         // per-axis activation jumps Mc_i/Kpen_i
  Matrix Kinit;            // diag(Kpenz, Kpeny) — closed-hinge initial tangent

  // committed history
  double CrMax;            // max elliptical norm r ever seen (irreversible damage driver)
  Vector Cstrain;          // committed jump (work anchor)
  Vector Cstress;          // committed moment (work anchor)
  double Cwork;            // committed path work int M . dalpha

  // trial state
  Vector Tstrain;          // [alpha_z, alpha_y]
  Vector Tstress;          // [M_z, M_y]
  Matrix Ttangent;         // 2x2 dM/dalpha
  double TrMax;
  double Twork;
};

#endif
