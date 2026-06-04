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

// Ladruno: rebar-buckling wrapper UniaxialMaterial.
//
//   A stress-modifying wrapper over any tension-compression UniaxialMaterial
//   ("the bare bar"). In compression past a slenderness-dependent onset strain
//   it applies the Dhakal-Maekawa (2002) buckled-average degradation:
//
//       sigma_buckled = r(e, lambda) * sigma_bare        (r in [resid, 1])
//
//   so the wrapped core stays byte-untouched (its V7 oracle property is
//   preserved) and the geometric buckling is a clean opt-in overlay. Tension
//   and pre-onset compression pass through unchanged; lsr <= 0 is the exact
//   identity gate (pure pass-through). The buckled average is referenced to the
//   tensile-excursion anchor (max tensile strain ε_max), so a re-loaded bar
//   straightens back toward the bare curve.
//
//   This is the "core = bare bar (a pure oracle), wrapper = the geometric
//   buckling average on top" decomposition. It is designed for the Ladruno
//   uniaxial steel (LadrunoUniaxialJ2) but composes over Steel02/Steel4 too,
//   and is stackable: Fatigue ∘ RebarBuckling ∘ J2.
//
//   Model:    Dhakal-Maekawa (-model dm, default) or Gomes-Appleton (-model ga).
//             Oracle: the verbatim formulae ported from
//             ReinforcingSteel::Buckled_stress_Dhakal / _Gomes(ess, fss).
//   Tangent:  analytic consistent tangent  dσ/dε = r·k_bare + σ_bare·(∂r/∂ε)
//             (improves on ReinforcingSteel, which finite-differences it).
//
//   classTag MAT_TAG_LadrunoRebarBuckling = 33001 (Ladruno private uniaxial
//   band; sibling of LadrunoUniaxialJ2 = 33000).
//   See Ladruno_implementation/14_ladruno_rebar_buckling_adr.md.
//   Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoRebarBuckling_h
#define LadrunoRebarBuckling_h

#include <UniaxialMaterial.h>

class LadrunoRebarBuckling : public UniaxialMaterial
{
 public:
  // buckling-law selector
  enum { MODEL_DM = 1, MODEL_GA = 2 };
  // re-straightening span mode (DECISION 3): L_rs = c*|e_rev-e_cross_rs| (floor eY)
  //                                     vs   L_rs = f(lambda)*eY
  enum { RS_C = 0, RS_LAMBDA = 1 };
  // cyclic branch (v2): the regime of the converged step
  enum { BR_PASS = 0, BR_BUCKLING = 1, BR_RESTRAIGHTEN = 2 };

  LadrunoRebarBuckling(int tag, UniaxialMaterial& bar,
                       double lsr, double alpha,
                       double gaReduction, double gaFsuFrac,
                       double fy, double E, int model,
                       int restraightenMode, double restraightenC);
  LadrunoRebarBuckling();
  ~LadrunoRebarBuckling();

  const char* getClassType(void) const { return "LadrunoRebarBuckling"; }

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void)          { return Tstrain; }
  double getStrainRate(void);
  double getStress(void)          { return Tstress; }
  double getTangent(void)         { return Ttangent; }
  double getDampTangent(void);
  double getInitialTangent(void)  { return (theBar != 0) ? theBar->getInitialTangent() : E; }

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial* getCopy(void);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  UniaxialMaterial* theBar;     // wrapped bare-bar (owns a getCopy)

  // model parameters
  int    model;                 // MODEL_DM | MODEL_GA
  double lsr;                   // slenderness s/d (LDratio); lsr<=0 => pass-through
  double alpha;                 // DM residual-shape factor (ReinforcingSteel beta, 0.75..1.0)
  double gaReduction;           // GA blend r in [0,1] (0=full GA buckling, 1=none); == RS GABuck r
  double gaFsuFrac;             // GA fsu_fraction (== RS GABuck gama); blends fsup anchor stress
  double fy;                    // yield stress  (for eY = fy/E and the -0.2 fy floor)
  double E;                     // Young's modulus (for eY and the deep-branch slope)
  int    restraightenMode;      // RS_C | RS_LAMBDA  (v2 cyclic span form, DECISION 3)
  double restraightenC;         // RS_C: L_rs = max(c*|e_rev-e_cross_rs|, eY)

  // ---- committed state (v1 anchor + onset cache, KEPT) ----
  double CmaxTenStrain;         // ε_max anchor (max engineering strain seen)
  double CmaxTenStress;         // bare stress at the ε_max anchor (for e_cross)
  double CfStarL;               // cached backbone stress at the onset strain
  double CfStarLcross;          // the e_cross for which CfStarL was computed

  // ---- committed cyclic state (v2: σ = σ_bare - g, g in [0, g_law]) ----
  int    Cbranch;               // BR_PASS | BR_BUCKLING | BR_RESTRAIGHTEN
  double Cstrain;               // committed engineering strain (direction detection)
  double Cstress;               // committed wrapper output stress (s_rev source)
  double Cg;                    // committed residual deficit g >= 0 below bare
  double Cg_law;                // committed v1 buckling-law gap g_law = sBare - sigma_v1
  // RESTRAIGHTEN episode latches (set at the BUCKLING->RESTRAIGHTEN reversal):
  double Ce_rev;                // reversal strain
  double Cs_rev;                // reversal (buckled) stress = continuity anchor
  double Ce_cross_rs;           // frozen anchor crossing = CmaxTenStrain - CmaxTenStress/E
  double Cf_bare_cross;         // σ_bare(Ce_cross_rs) via backboneStressAt, latched at reversal
  double CL_rs;                 // re-straightening span (§9.4)
  // re-buckle carry latches (set at RESTRAIGHTEN->BUCKLING; Ce_reentry=LRB_INVALID => no carry):
  double Ce_reentry;            // strain at re-compression
  double Cg_law_reentry;        // g_law at re-entry (the offset reference)
  double Cg_reentry;            // deficit g at re-entry (the carried value)
  // deepest-buckle envelope (monotone; bounds the carry re-merge):
  double Ce_bow;                // deepest committed buckle strain (LRB_INVALID => none)
  double Cg_bow;               // g_law at Ce_bow (the value g re-merges to)

  // ---- trial twins of all of the above ----
  double TmaxTenStrain, TmaxTenStress, TfStarL, TfStarLcross;
  int    Tbranch;
  double Tg, Tg_law;
  double Te_rev, Ts_rev, Te_cross_rs, Tf_bare_cross, TL_rs;
  double Te_reentry, Tg_law_reentry, Tg_reentry;
  double Te_bow, Tg_bow;

  // ---- trial output ----
  double Tstrain;
  double Tstress;
  double Ttangent;
  double Tr;                    // effective reduction factor σ/σ_bare (for the "reduction" response)

  // ---- helpers ----
  // v1 monotonic buckling-law value (DM or GA) -> (sig, tan); updates the fStarL cache.
  void   lawV1(double eps, double sBare, double kBare, double& sig, double& tan);
  void   applyBucklingDM(double eps, double sBare, double kBare, double& sig, double& tan);
  void   applyBucklingGA(double eps, double sBare, double kBare, double& sig, double& tan);
  double gaFactor(double eps, double e_cross) const;              // closed-form GA factor(eps)
  double backboneStressAt(double engStrain);                     // committed-clone probe
  double computeLrs(double e_rev, double e_cross_rs) const;       // §9.4 span (mode-selected)
  bool   onsetReached(double eps, double e_cross_live) const;     // law-specific PASS->BUCKLING gate
  void   copyCommittedToTrial(void);                             // C->T for the cyclic latches
  // BUCKLING output: sigma = sBare - g (residual-gap reconciliation, with re-buckle carry);
  // updates Tstress/Ttangent/Tr/Tg and deepens the Te_bow/Tg_bow envelope.
  void   bucklingBranch(double eps, double sBare, double kBare, double sigV1, double kV1);
  // RESTRAIGHTEN output: Phase-1 elastic unload then Phase-2 smoothstep gap-close.
  void   restraightenBranch(double eps, double sBare, double kBare);

  int parameterID;
};

#endif
