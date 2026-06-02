/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: rebar-buckling WRAPPER UniaxialMaterial.
//
// Applies Dhakal-Maekawa (2002) reinforcing-bar buckling on top of ANY bare-bar
// UniaxialMaterial (e.g. LadrunoUniaxialJ2, Steel02, Steel4). DM is a post-hoc
// modification of the bare-bar stress at the same strain:
//     sigma_buckled = r(eps, lsr) * sigma_bare        (compression branch)
// so this is a stress-modifying wrapper that leaves the wrapped material a pure
// oracle (its plastic state evolves on the TRUE strain underneath). Tension and
// pre-onset compression pass through unchanged.
//
// The buckled average compressive stress departs the bare curve at a
// slenderness-dependent onset strain eStar(lsr) and softens linearly toward a
// residual ~ -0.2 fy. Mirrors ReinforcingSteel::Buckled_stress_Dhakal (the port
// oracle), but as a composable layer rather than baked into one monolithic
// material; consistent tangent is ANALYTIC (upstream finite-differences it).
//
// v1: monotonic compression buckling (from the tensile-strain anchor e_cross).
// Full cyclic re-straightening transitions are a documented follow-up (B4).
//
// classTag MAT_TAG_LadrunoRebarBuckling = 33001 (Ladruno uniaxial band).
// See Ladruno_implementation/14_ladruno_rebar_buckling_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoRebarBuckling_h
#define LadrunoRebarBuckling_h

#include <UniaxialMaterial.h>

class LadrunoRebarBuckling : public UniaxialMaterial
{
 public:
  enum { MODEL_DM = 0, MODEL_GA = 1 };
  static const double LSR_MIN;     // slenderness below which buckling is disabled

  LadrunoRebarBuckling(int tag, UniaxialMaterial& bar, double lsr,
                       double fy, double E, int model = MODEL_DM,
                       double alpha = 1.0);
  LadrunoRebarBuckling();
  ~LadrunoRebarBuckling();

  const char* getClassType(void) const { return "LadrunoRebarBuckling"; }

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void)          { return theBar->getStrain(); }
  double getStress(void)          { return Tstress; }
  double getTangent(void)         { return Ttangent; }
  double getInitialTangent(void)  { return theBar->getInitialTangent(); }

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
  UniaxialMaterial* theBar;     // wrapped bare-bar material (owned)

  // parameters
  double lsr;                   // slenderness L/D (= s/d); <= LSR_MIN => pass-through
  double fy;                    // yield stress (for eStar / eyp)
  double E;                     // elastic modulus (for eyp = fy/E)
  int    model;                 // MODEL_DM | MODEL_GA
  double alpha;                 // DM amplification (beta in ReinforcingSteel)

  // committed buckling state (the ONLY state this wrapper owns)
  double CmaxTenStrain;         // max tensile strain reached (e_cross anchor)
  double CfStarL;               // bare stress anchored at the buckling onset eStar
  int    Conset;                // 0 = not yet buckled this excursion, 1 = past onset

  // trial mirrors
  double TmaxTenStrain;
  double TfStarL;
  int    Tonset;

  // trial state
  double Tstress;
  double Ttangent;

  // helpers
  double eStarStrain(void) const;                  // onset strain eStar(lsr,fy,E) (compressive, <0)
  // DM buckled average stress at relative strain e (= eps - e_cross); fills the
  // consistent tangent dSigma/dEps in tanOut (uses the anchored CfStarL/TfStarL).
  double dmStress(double e, double sBare, double kBare, double fStarL,
                  double& tanOut) const;
};

#endif
