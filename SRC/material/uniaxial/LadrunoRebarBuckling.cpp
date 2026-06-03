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

// Ladruno: rebar-buckling wrapper UniaxialMaterial. See LadrunoRebarBuckling.h
// and Ladruno_implementation/14_ladruno_rebar_buckling_adr.md for the model.
//
// The Dhakal-Maekawa buckled-average degradation is ported verbatim from
// ReinforcingSteel::Buckled_stress_Dhakal(ess, fss) but works in engineering
// strain on a generic wrapped bar (ReinforcingSteel works in log strain on its
// own monolithic backbone), and re-implements the cyclic bookkeeping with its
// own minimal state (an ε_max anchor) instead of the ReinforcingSteel
// branch-state machine. The consistent tangent is analytic (product rule)
// rather than finite-differenced.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoRebarBuckling.h>
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

static const double LRB_INVALID = 1.0e30;   // "fStarL not yet computed" sentinel

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoRebarBuckling tag matTag -lsr s/d
//                    <-model dm|ga> <-alpha a> <-reduction r> <-fsufrac g>
//                    <-fy fy> <-E E>
// ===========================================================================
void* OPS_LadrunoRebarBuckling(void)
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoRebarBuckling tag? matTag? -lsr s/d? "
           << "<-model dm|ga> <-alpha a?> <-reduction r?> <-fsufrac g?> "
           << "<-fy fy?> <-E E?>\n";
    return 0;
  }

  int idata[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, idata) < 0) {
    opserr << "WARNING invalid LadrunoRebarBuckling tag/matTag\n";
    return 0;
  }

  double lsr   = 0.0;          // off by default (identity gate)
  double alpha = 1.0;          // DM residual-shape factor (ReinforcingSteel beta)
  double gaReduction = 0.0;    // GA blend r: 0 = full GA buckling (wrapper default), 1 = none
  double gaFsuFrac   = 0.5;    // GA fsu_fraction (ReinforcingSteel GABuck "gama" default)
  double fy    = 0.0;          // queried/required when lsr>0
  double E     = 0.0;          // <=0 => take theBar->getInitialTangent()
  int    model = LadrunoRebarBuckling::MODEL_DM;
  int    restraightenMode = LadrunoRebarBuckling::RS_C;   // v2 cyclic span form (DECISION 3)
  double restraightenC    = 1.0;                          // RS_C fraction of the compressive excursion

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-lsr") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &lsr) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -lsr wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-alpha") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &alpha) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -alpha wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-reduction") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &gaReduction) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -reduction wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-fsufrac") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &gaFsuFrac) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -fsufrac wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-fy") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &fy) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -fy wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-E") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &E) < 0) {
        opserr << "WARNING LadrunoRebarBuckling: -E wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-model") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "dm") == 0 || strcmp(law, "DM") == 0) {
        model = LadrunoRebarBuckling::MODEL_DM;
      } else if (strcmp(law, "ga") == 0 || strcmp(law, "GA") == 0) {
        model = LadrunoRebarBuckling::MODEL_GA;
      } else {
        opserr << "WARNING LadrunoRebarBuckling: unknown -model '" << law
               << "' (dm or ga)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-restraighten") == 0) {
      // v2 cyclic re-straightening span (DECISION 3):
      //   -restraighten lambda            -> L_rs = f(lambda)*eY
      //   -restraighten c <value>         -> L_rs = max(value*|e_rev-e_cross_rs|, eY)
      const char* mode = OPS_GetString();
      if (strcmp(mode, "lambda") == 0 || strcmp(mode, "lsr") == 0) {
        restraightenMode = LadrunoRebarBuckling::RS_LAMBDA;
      } else if (strcmp(mode, "c") == 0) {
        restraightenMode = LadrunoRebarBuckling::RS_C;
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &restraightenC) < 0) {
          opserr << "WARNING LadrunoRebarBuckling: -restraighten c wants a value\n";
          return 0;
        }
      } else {
        opserr << "WARNING LadrunoRebarBuckling: unknown -restraighten '" << mode
               << "' (lambda or c <value>)\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoRebarBuckling: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* bar = OPS_getUniaxialMaterial(idata[1]);
  if (bar == 0) {
    opserr << "WARNING LadrunoRebarBuckling: wrapped material " << idata[1]
           << " does not exist (for material " << idata[0] << ")\n";
    return 0;
  }

  if (E <= 0.0)
    E = bar->getInitialTangent();

  if (gaReduction < 0.0) gaReduction = 0.0;
  if (gaReduction > 1.0) gaReduction = 1.0;

  // DM needs fy (onset eStar/eY + the -0.2 fy floor); GA does not (it uses the
  // fsup anchor stress + lsr + reduction/fsu_fraction). Both need E>0 (e_cross).
  if (lsr > 0.0 && E <= 0.0) {
    opserr << "WARNING LadrunoRebarBuckling: with -lsr>0 a positive E is required "
           << "(via -E or the wrapped material's initial tangent); material "
           << idata[0] << "\n";
    return 0;
  }
  if (lsr > 0.0 && model == LadrunoRebarBuckling::MODEL_DM && fy <= 0.0) {
    opserr << "WARNING LadrunoRebarBuckling: -model dm with -lsr>0 needs a "
           << "positive -fy; material " << idata[0] << "\n";
    return 0;
  }

  UniaxialMaterial* theMat =
    new LadrunoRebarBuckling(idata[0], *bar, lsr, alpha,
                             gaReduction, gaFsuFrac, fy, E, model,
                             restraightenMode, restraightenC);
  if (theMat == 0) {
    opserr << "WARNING LadrunoRebarBuckling: failed to allocate material\n";
    return 0;
  }
  return theMat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
// Initialize ALL cyclic state to the virgin PASS configuration. Inline (NOT via
// revertToStart, which forwards to theBar and would wipe a getCopy-preserved
// state -- the classic wrapper footgun; FatigueMaterial/MinMaxMaterial too).
#define LRB_INIT_STATE()                                                        \
  do {                                                                          \
    CmaxTenStrain = CmaxTenStress = CfStarL = 0.0; CfStarLcross = LRB_INVALID;  \
    TmaxTenStrain = TmaxTenStress = TfStarL = 0.0; TfStarLcross = LRB_INVALID;  \
    Cbranch = Tbranch = BR_PASS;                                                \
    Cstrain = Cstress = Cg = Cg_law = 0.0;                                      \
    Ce_rev = Cs_rev = Ce_cross_rs = Cf_bare_cross = CL_rs = 0.0;                \
    Ce_reentry = LRB_INVALID; Cg_law_reentry = Cg_reentry = 0.0;               \
    Ce_bow = LRB_INVALID; Cg_bow = 0.0;                                         \
    Tg = Tg_law = 0.0;                                                          \
    Te_rev = Ts_rev = Te_cross_rs = Tf_bare_cross = TL_rs = 0.0;                \
    Te_reentry = LRB_INVALID; Tg_law_reentry = Tg_reentry = 0.0;               \
    Te_bow = LRB_INVALID; Tg_bow = 0.0;                                         \
    Tstrain = Tstress = 0.0; Tr = 1.0;                                          \
  } while (0)

LadrunoRebarBuckling::LadrunoRebarBuckling(int tag, UniaxialMaterial& bar,
                                           double lsr_, double alpha_,
                                           double gaReduction_, double gaFsuFrac_,
                                           double fy_, double E_, int model_,
                                           int restraightenMode_, double restraightenC_)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), model(model_), lsr(lsr_), alpha(alpha_),
    gaReduction(gaReduction_), gaFsuFrac(gaFsuFrac_), fy(fy_), E(E_),
    restraightenMode(restraightenMode_), restraightenC(restraightenC_),
    parameterID(0)
{
  theBar = bar.getCopy();
  if (theBar == 0) {
    opserr << "LadrunoRebarBuckling::LadrunoRebarBuckling -- failed to get a "
           << "copy of the wrapped material\n";
  }
  if (E <= 0.0 && theBar != 0)
    E = theBar->getInitialTangent();
  LRB_INIT_STATE();
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : E;
}

LadrunoRebarBuckling::LadrunoRebarBuckling()
  : UniaxialMaterial(0, MAT_TAG_LadrunoRebarBuckling),
    theBar(0), model(MODEL_DM), lsr(0.0), alpha(1.0),
    gaReduction(0.0), gaFsuFrac(0.5), fy(0.0), E(0.0),
    restraightenMode(RS_C), restraightenC(1.0),
    parameterID(0)
{
  // theBar is built by recvSelf via the broker.
  LRB_INIT_STATE();
  Ttangent = 0.0;
}

LadrunoRebarBuckling::~LadrunoRebarBuckling()
{
  if (theBar)
    delete theBar;
}

// ===========================================================================
//  committed-clone backbone probe -- the bare stress at an arbitrary strain,
//  reached monotonically from the last committed state. Material-agnostic.
// ===========================================================================
double LadrunoRebarBuckling::backboneStressAt(double engStrain)
{
  if (theBar == 0)
    return 0.0;
  UniaxialMaterial* probe = theBar->getCopy();   // copies theBar's trial state
  if (probe == 0)
    return 0.0;
  probe->revertToLastCommit();                   // back to the committed anchor
  probe->setTrialStrain(engStrain);              // push to the onset strain
  double s = probe->getStress();
  delete probe;
  return s;
}

// ===========================================================================
//  Dhakal-Maekawa buckled-average overlay (sets Tstress / Ttangent / Tr)
//  Ported from ReinforcingSteel::Buckled_stress_Dhakal(ess, fss); the clean
//  r*sigma_bare branch (the monotonic / TBranchNum%4>1 form).
// ===========================================================================
void LadrunoRebarBuckling::applyBucklingDM(double eps, double sBare, double kBare,
                                          double& sig, double& tan)
{
  const double eY = fy / E;                        // yield strain
  const double e_cross = TmaxTenStrain - TmaxTenStress / E;   // tensile-unload anchor
  const double e = eps - e_cross;                  // strain from the anchor

  // pre-onset (incl. tension): exact pass-through
  if (e >= -eY) {
    sig = sBare; tan = kBare;
    return;
  }

  const double kfac = sqrt(eY * 2000.0);           // = sqrt(fy/E*2000)
  double eStarMag = 55.0 - 2.3 * kfac * lsr;
  if (eStarMag < 7.0) eStarMag = 7.0;
  const double eStar = -eStarMag * eY;             // onset strain (negative)

  // backbone stress at the onset strain (cached per anchor). Mirrors upstream
  // ReinforcingSteel::Buckled_stress_Dhakal, which evaluates Backbone_f(eStar)
  // with eStar treated as the absolute backbone strain (it equals the onset's
  // absolute strain when e_cross==0, the monotonic-from-virgin case).
  if (TfStarLcross != e_cross) {
    TfStarL = backboneStressAt(eStar);
    TfStarLcross = e_cross;
  }
  const double fStarL = TfStarL;
  if (fabs(fStarL) < 1.0e-12) {                    // degenerate bare backbone => no overlay
    sig = sBare; tan = kBare;
    return;
  }

  double fStar = fStarL * alpha * (1.1 - 0.016 * kfac * lsr);
  if (fStar > -0.2 * fy) fStar = -0.2 * fy;        // residual floor

  if (e >= eStar) {
    // intermediate branch:  -eY > e >= eStar
    const double ratio = 1.0 - (1.0 - fStar / fStarL) * (e + eY) / (eStar + eY);
    const double dratio = -(1.0 - fStar / fStarL) / (eStar + eY);   // d(ratio)/d(eps)
    sig = sBare * ratio;
    tan = kBare * ratio + sBare * dratio;
  } else {
    // deep branch:  e < eStar  (ramp toward the -0.2 fy floor)
    const double num = fStar - 0.02 * E * (e - eStar);
    double ave = sBare * num / fStarL;
    if (ave > -0.2 * fy) {
      // floored: flat residual => (near) zero tangent
      sig = -0.2 * fy;
      tan = 0.0;
    } else {
      sig = ave;
      tan = (kBare * num + sBare * (-0.02 * E)) / fStarL;
    }
  }
}

// ===========================================================================
//  Gomes-Appleton buckled-average overlay
//  Ported from ReinforcingSteel::Buckled_stress_Gomes(ess, fss). The upstream
//  stress function hardcodes its own beta=1, gama=0.1, Dft=0.25 (shadowing the
//  user GABuck "beta"), so only lsr, reduction (GABuck r) and fsu_fraction
//  (GABuck gama) plus the fsup anchor stress actually drive the result.
// ===========================================================================
double LadrunoRebarBuckling::gaFactor(double eps, double e_cross) const
{
  const double C = 9.42477796076938;               // 3*pi (upstream constant)
  double denom = e_cross - eps;                     // > 0 in the buckled region
  if (denom < 1.0e-300) denom = 1.0e-300;
  const double fs_buck = sqrt(32.0 / denom) / (C * lsr);
  double betaLoc = 1.0;
  const double gamaLoc = 0.1, Dft = 0.25;
  const double stress_diff = fabs(fs_buck - 1.0);
  if (stress_diff <= Dft) betaLoc = 1.0 - gamaLoc * (Dft - stress_diff) / Dft;
  const double m = (fs_buck < 1.0) ? fs_buck : 1.0;
  return m * betaLoc * (1.0 - gaReduction) + gaReduction;
}

void LadrunoRebarBuckling::applyBucklingGA(double eps, double sBare, double kBare,
                                          double& sig, double& tan)
{
  const double e_cross = TmaxTenStrain - TmaxTenStress / E;   // tensile-unload anchor

  // GA gate: pass-through unless compressed past the anchor crossing
  if (eps >= e_cross) {
    sig = sBare; tan = kBare;
    return;
  }

  const double fsup = TmaxTenStress;               // peak tensile (anchor) stress
  const double ff   = gaFsuFrac;                    // fsu_fraction
  const double factor = this->gaFactor(eps, e_cross);

  //   sigma = fsup*ff - (factor+ff)*(fsup*ff - sBare)/(1+ff)
  sig = fsup * ff - (factor + ff) * (fsup * ff - sBare) / (1.0 + ff);

  // Consistent tangent dsigma/deps = rho*k_bare + (dsigma/dfactor)*dfactor/deps,
  //   rho = dsigma/dsBare = (factor+ff)/(1+ff),
  //   dsigma/dfactor = -(fsup*ff - sBare)/(1+ff),
  //   dfactor/deps from a central FD of the closed-form scalar factor (cheap, no
  // material probe). More accurate than ReinforcingSteel's Buckled_mod_Gomes,
  // which FDs the whole stress with sBare frozen then adds k_bare (i.e. uses
  // k_bare in place of rho*k_bare).
  const double rho = (factor + ff) / (1.0 + ff);
  const double dSig_dFactor = -(fsup * ff - sBare) / (1.0 + ff);
  const double h = 1.0e-8;
  double epsm = eps - h, epsp = eps + h;
  if (epsp >= e_cross) epsp = eps;                 // stay inside the buckled region
  const double dfac = (this->gaFactor(epsp, e_cross) - this->gaFactor(epsm, e_cross))
                      / (epsp - epsm);
  tan = rho * kBare + dSig_dFactor * dfac;
}

// ===========================================================================
//  lawV1 -- the v1 monotonic buckling-law value (dispatch DM/GA), into (sig,tan).
//  This is the "g_law" source: g_law = sBare - sig. For the identity/pre-onset
//  region it returns (sBare,kBare) so g_law = 0. Mutates the fStarL cache only.
// ===========================================================================
void LadrunoRebarBuckling::lawV1(double eps, double sBare, double kBare,
                                 double& sig, double& tan)
{
  if (model == MODEL_GA)
    this->applyBucklingGA(eps, sBare, kBare, sig, tan);
  else
    this->applyBucklingDM(eps, sBare, kBare, sig, tan);
}

// ===========================================================================
//  computeLrs -- the re-straightening span (DECISION 3, §9.4), mode-selected.
// ===========================================================================
double LadrunoRebarBuckling::computeLrs(double e_rev, double e_cross_rs) const
{
  const double eY = fy / E;
  if (restraightenMode == RS_LAMBDA) {
    // f(lambda) = 0.5 + 0.10*max(0, lambda-5), clamped to [0.5, 6.0]
    double f = 0.5 + 0.10 * (lsr > 5.0 ? (lsr - 5.0) : 0.0);
    if (f < 0.5) f = 0.5;
    if (f > 6.0) f = 6.0;
    return f * eY;
  }
  // RS_C: fraction of the compressive excursion, floored at eY
  double span = restraightenC * fabs(e_rev - e_cross_rs);
  return (span > eY) ? span : eY;
}

// ===========================================================================
//  onsetReached -- the law-specific PASS->BUCKLING gate (live anchor crossing).
// ===========================================================================
bool LadrunoRebarBuckling::onsetReached(double eps, double e_cross_live) const
{
  const double eY = fy / E;
  if (model == MODEL_GA)
    return (eps < e_cross_live);                   // any compression past the anchor
  return ((eps - e_cross_live) < -eY);             // DM: past the yield offset
}

// ===========================================================================
//  state setting
// ===========================================================================
// ---------------------------------------------------------------------------
//  copy the committed cyclic latches into their trial twins (steady-state
//  default; the state machine overrides on a transition).
// ---------------------------------------------------------------------------
void LadrunoRebarBuckling::copyCommittedToTrial(void)
{
  Tbranch        = Cbranch;
  Te_rev         = Ce_rev;        Ts_rev         = Cs_rev;
  Te_cross_rs    = Ce_cross_rs;   Tf_bare_cross  = Cf_bare_cross;
  TL_rs          = CL_rs;
  Te_reentry     = Ce_reentry;    Tg_law_reentry = Cg_law_reentry;
  Tg_reentry     = Cg_reentry;
  Te_bow         = Ce_bow;        Tg_bow         = Cg_bow;
}

// ---------------------------------------------------------------------------
//  BUCKLING branch:  sigma = sBare - g,  g in [0, g_law].
//   no carry  : g = g_law  (=> sigma = sigma_v1, BIT-IDENTICAL v1, tangent k_v1)
//   carry     : g = g_law - Delta*(1 - S(phi)),  Delta = Cg_law_reentry - Cg_reentry,
//               phi = (g_law - Cg_law_reentry)/(Cg_bow - Cg_law_reentry), S = smoothstep,
//               re-merging onto g_law at the deepest committed buckle e_bow (phi=1).
// ---------------------------------------------------------------------------
void LadrunoRebarBuckling::bucklingBranch(double eps, double sBare, double kBare,
                                          double sigV1, double kV1)
{
  double g, tang;
  const bool carry = (Te_reentry != LRB_INVALID);

  if (!carry) {
    g    = Tg_law;                                 // exact v1 (g == g_law)
    tang = kV1;
  } else if (Te_bow == LRB_INVALID || eps <= Te_bow) {
    // re-compressed to / beyond the deepest prior bow: re-merged onto the v1
    // envelope (and the envelope extends below). Clear the carry.
    g = Tg_law; tang = kV1; Te_reentry = LRB_INVALID;
  } else {
    // Carry: smoothstep the deficit in STRAIN from the carried value Tg_reentry
    // (at Te_reentry) to the deepest-bow gap Tg_bow (at Te_bow), then follow v1
    // beyond e_bow. Strain-based (NOT g_law-based) so it is C0 even when
    // re-entry happens on the tension side where g_law == 0 but the bow deficit
    // Tg_reentry > 0.  xi = 0 at re-entry, 1 at e_bow.
    const double span = Te_reentry - Te_bow;       // > 0 (e_bow is more compressive)
    if (span <= 1.0e-12) {
      g = Tg_law; tang = kV1; Te_reentry = LRB_INVALID;
    } else {
      double xi = (Te_reentry - eps) / span;
      if (xi < 0.0) xi = 0.0;
      if (xi > 1.0) xi = 1.0;
      const double S  = xi * xi * (3.0 - 2.0 * xi);    // smoothstep 3xi^2 - 2xi^3
      const double Sp = 6.0 * xi * (1.0 - xi);         // dS/dxi
      g = Tg_reentry + (Tg_bow - Tg_reentry) * S;
      const double dg_deps = (Tg_bow - Tg_reentry) * Sp * (-1.0 / span);
      tang = kBare - dg_deps;
    }
  }

  if (g < 0.0) g = 0.0;                             // raise above bare is non-negative
  Tg       = g;
  Tstress  = sBare + g;                             // sigma = sigma_bare + raise
  Ttangent = tang;
  Tr       = (sBare != 0.0) ? Tstress / sBare : 1.0;

  // deepen the monotone bow envelope (the carry re-merge target)
  if (Te_bow == LRB_INVALID || eps < Te_bow) {
    Te_bow = eps;
    Tg_bow = Tg_law;
  }
}

// ---------------------------------------------------------------------------
//  RESTRAIGHTEN branch:  Phase 1 (elastic unload at k_rev) then Phase 2
//  (smoothstep gap-close).  Uses the TRIAL episode latches (Te_rev/Ts_rev/
//  Te_cross_rs/Tf_bare_cross/TL_rs), which are committed-derived.
// ---------------------------------------------------------------------------
void LadrunoRebarBuckling::restraightenBranch(double eps, double sBare, double kBare)
{
  // Track the LIVE bare curve plus a buckling raise g that is HELD on the
  // compression side (no straightening yet) and DECAYED to 0 over TL_rs on the
  // tension side (eps >= Te_cross_rs). g_rev = Ts_rev is the raise at reversal.
  //   sigma = sigma_bare(eps) + g,   g in [0, g_rev].
  const double g_rev = Ts_rev;
  double g, tang;

  if (eps < Te_cross_rs) {
    // Phase 1 -- hold the raise, recover parallel to bare (yields with the bar)
    g    = g_rev;
    tang = kBare;
  } else {
    // Phase 2 -- tension-side straightening (smoothstep decay of the raise)
    double q = (TL_rs > 0.0) ? (eps - Te_cross_rs) / TL_rs : 1.0;
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    const double B  = 1.0 - q * q * (3.0 - 2.0 * q);   // 1 - smoothstep, B(0)=1, B(1)=0
    const double Bp = -6.0 * q * (1.0 - q);            // dB/dq
    g    = g_rev * B;
    tang = kBare + g_rev * Bp / ((TL_rs > 0.0) ? TL_rs : 1.0);
  }

  if (g < 0.0) g = 0.0;
  Tg       = g;
  Tstress  = sBare + g;
  Ttangent = tang;
  Tr = (sBare != 0.0) ? Tstress / sBare : 1.0;
}

// ===========================================================================
//  state setting -- the v2 cyclic state machine (committed-latch / idempotent).
//  Output is recomputed from scratch every call as a pure function of the
//  COMMITTED state + the trial strain; transitions persist only in commitState.
// ===========================================================================
int LadrunoRebarBuckling::setTrialStrain(double strain, double strainRate)
{
  if (theBar == 0)
    return -1;

  int res = theBar->setTrialStrain(strain, strainRate);
  double sBare = theBar->getStress();
  double kBare = theBar->getTangent();
  Tstrain = strain;

  // v1 tensile-strain anchor (live, path-dependent BY MANDATE -- keeps the
  // PASS->BUCKLING onset + v1 BUCKLING bit-identical to v1). Invalidate the
  // fStarL cache if the anchor moves.
  if (strain > TmaxTenStrain) {
    TmaxTenStrain = strain;
    TmaxTenStress = sBare;
    TfStarLcross  = LRB_INVALID;
  }

  // seed the trial latches from committed (steady state); transitions override
  this->copyCommittedToTrial();

  // identity gate (B0): exact pass-through. DM also needs fy>0; GA does not.
  if (lsr <= 0.0 || E <= 0.0 || (model == MODEL_DM && fy <= 0.0)) {
    Tstress = sBare; Ttangent = kBare; Tr = 1.0;
    Tg = 0.0; Tg_law = 0.0; Tbranch = BR_PASS;
    return res;
  }

  const double eY   = fy / E;
  const double tol  = 1.0e-12 + 1.0e-8 * eY;
  const double dEps = strain - Cstrain;
  const double e_cross_live = TmaxTenStrain - TmaxTenStress / E;

  // v1 law value (g_law source), always evaluated. g is the RAISE of the
  // buckled stress ABOVE bare: g = sigma - sigma_bare >= 0 (buckling reduces the
  // compressive magnitude, i.e. raises sigma toward 0). So g_law = sigma_v1 - sBare.
  double sigV1, kV1;
  this->lawV1(strain, sBare, kBare, sigV1, kV1);
  Tg_law = sigV1 - sBare;
  if (Tg_law < 0.0) Tg_law = 0.0;

  if (Cbranch == BR_PASS) {
    if (this->onsetReached(strain, e_cross_live)) {
      Tbranch = BR_BUCKLING;
      Te_reentry = LRB_INVALID;                    // virgin buckle: no carry
      this->bucklingBranch(strain, sBare, kBare, sigV1, kV1);
    } else {
      Tbranch = BR_PASS;
      Tstress = sBare; Ttangent = kBare; Tr = 1.0; Tg = 0.0;
    }
  }
  else if (Cbranch == BR_BUCKLING) {
    if (dEps > tol) {
      // reversal -> RESTRAIGHTEN: latch the episode from committed state.
      // Ts_rev holds h_rev = the RAISE (sigma - sigma_bare) at the reversal,
      // which is exactly the committed g (Cg). The straighten-start Te_cross_rs
      // is the frozen tensile anchor crossing; the raise is held until then,
      // then decayed to 0 over TL_rs. (No elastic-line probe is needed.)
      Tbranch        = BR_RESTRAIGHTEN;
      Te_rev         = Cstrain;
      Ts_rev         = Cg;                         // h_rev = committed raise at reversal
      Te_cross_rs    = CmaxTenStrain - CmaxTenStress / E;
      TL_rs          = this->computeLrs(Te_rev, Te_cross_rs);
      Tf_bare_cross  = 0.0;                        // reserved (elastic-line probe removed)
      Te_reentry     = LRB_INVALID;                // leaving any carry
      this->restraightenBranch(strain, sBare, kBare);
    } else {
      Tbranch = BR_BUCKLING;
      this->bucklingBranch(strain, sBare, kBare, sigV1, kV1);
    }
  }
  else { // Cbranch == BR_RESTRAIGHTEN
    const double q = (CL_rs > 0.0) ? (strain - Ce_cross_rs) / CL_rs : 1.0;
    if (q >= 1.0) {
      // fully straightened -> PASS (bow gone: reset the episode + envelope)
      Tbranch = BR_PASS;
      Tstress = sBare; Ttangent = kBare; Tr = 1.0; Tg = 0.0;
      Te_reentry = LRB_INVALID; Te_bow = LRB_INVALID; Tg_bow = 0.0;
    } else if (dEps < -tol) {
      // re-compress before completion -> BUCKLING with residual-gap carry
      Tbranch        = BR_BUCKLING;
      Te_reentry     = Cstrain;
      Tg_law_reentry = Cg_law;
      Tg_reentry     = Cg;
      this->bucklingBranch(strain, sBare, kBare, sigV1, kV1);
    } else {
      Tbranch = BR_RESTRAIGHTEN;
      this->restraightenBranch(strain, sBare, kBare);
    }
  }
  return res;
}

double LadrunoRebarBuckling::getStrainRate(void)
{
  return (theBar != 0) ? theBar->getStrainRate() : 0.0;
}

double LadrunoRebarBuckling::getDampTangent(void)
{
  return (theBar != 0) ? theBar->getDampTangent() : 0.0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRebarBuckling::commitState(void)
{
  // v1 anchor + onset cache
  CmaxTenStrain = TmaxTenStrain;
  CmaxTenStress = TmaxTenStress;
  CfStarL       = TfStarL;
  CfStarLcross  = TfStarLcross;
  // v2 cyclic state (the authoritative branch + latch transition persists here)
  Cbranch       = Tbranch;
  Cstrain       = Tstrain;
  Cstress       = Tstress;
  Cg            = Tg;
  Cg_law        = Tg_law;
  Ce_rev        = Te_rev;        Cs_rev         = Ts_rev;
  Ce_cross_rs   = Te_cross_rs;   Cf_bare_cross  = Tf_bare_cross;
  CL_rs         = TL_rs;
  Ce_reentry    = Te_reentry;    Cg_law_reentry = Tg_law_reentry;
  Cg_reentry    = Tg_reentry;
  Ce_bow        = Te_bow;        Cg_bow         = Tg_bow;
  return (theBar != 0) ? theBar->commitState() : 0;
}

int LadrunoRebarBuckling::revertToLastCommit(void)
{
  TmaxTenStrain = CmaxTenStrain;
  TmaxTenStress = CmaxTenStress;
  TfStarL       = CfStarL;
  TfStarLcross  = CfStarLcross;
  Tbranch       = Cbranch;
  Tstrain       = Cstrain;
  Tstress       = Cstress;
  Tg            = Cg;
  Tg_law        = Cg_law;
  Te_rev        = Ce_rev;        Ts_rev         = Cs_rev;
  Te_cross_rs   = Ce_cross_rs;   Tf_bare_cross  = Cf_bare_cross;
  TL_rs         = CL_rs;
  Te_reentry    = Ce_reentry;    Tg_law_reentry = Cg_law_reentry;
  Tg_reentry    = Cg_reentry;
  Te_bow        = Ce_bow;        Tg_bow         = Cg_bow;
  // Ttangent/Tr are recomputed on the next setTrialStrain.
  return (theBar != 0) ? theBar->revertToLastCommit() : 0;
}

int LadrunoRebarBuckling::revertToStart(void)
{
  // Here a full reset (incl. the bar) IS intended -- this is NOT the getCopy
  // footgun path, so forwarding to theBar->revertToStart() is correct.
  LRB_INIT_STATE();
  Ttangent = (theBar != 0) ? theBar->getInitialTangent() : E;
  Tr = 1.0;
  return (theBar != 0) ? theBar->revertToStart() : 0;
}

UniaxialMaterial* LadrunoRebarBuckling::getCopy(void)
{
  LadrunoRebarBuckling* theCopy =
    new LadrunoRebarBuckling(this->getTag(), *theBar, lsr, alpha,
                             gaReduction, gaFsuFrac, fy, E, model,
                             restraightenMode, restraightenC);

  // committed (v1 anchor/cache + v2 cyclic)
  theCopy->CmaxTenStrain = CmaxTenStrain;  theCopy->CmaxTenStress = CmaxTenStress;
  theCopy->CfStarL       = CfStarL;        theCopy->CfStarLcross  = CfStarLcross;
  theCopy->Cbranch       = Cbranch;
  theCopy->Cstrain       = Cstrain;        theCopy->Cstress       = Cstress;
  theCopy->Cg            = Cg;             theCopy->Cg_law        = Cg_law;
  theCopy->Ce_rev        = Ce_rev;         theCopy->Cs_rev        = Cs_rev;
  theCopy->Ce_cross_rs   = Ce_cross_rs;    theCopy->Cf_bare_cross = Cf_bare_cross;
  theCopy->CL_rs         = CL_rs;
  theCopy->Ce_reentry    = Ce_reentry;     theCopy->Cg_law_reentry= Cg_law_reentry;
  theCopy->Cg_reentry    = Cg_reentry;
  theCopy->Ce_bow        = Ce_bow;         theCopy->Cg_bow        = Cg_bow;
  // trial twins
  theCopy->TmaxTenStrain = TmaxTenStrain;  theCopy->TmaxTenStress = TmaxTenStress;
  theCopy->TfStarL       = TfStarL;        theCopy->TfStarLcross  = TfStarLcross;
  theCopy->Tbranch       = Tbranch;
  theCopy->Tg            = Tg;             theCopy->Tg_law        = Tg_law;
  theCopy->Te_rev        = Te_rev;         theCopy->Ts_rev        = Ts_rev;
  theCopy->Te_cross_rs   = Te_cross_rs;    theCopy->Tf_bare_cross = Tf_bare_cross;
  theCopy->TL_rs         = TL_rs;
  theCopy->Te_reentry    = Te_reentry;     theCopy->Tg_law_reentry= Tg_law_reentry;
  theCopy->Tg_reentry    = Tg_reentry;
  theCopy->Te_bow        = Te_bow;         theCopy->Tg_bow        = Tg_bow;
  theCopy->Tstrain       = Tstrain;        theCopy->Tstress       = Tstress;
  theCopy->Ttangent      = Ttangent;       theCopy->Tr            = Tr;

  return theCopy;
}

// ===========================================================================
//  parallel
// ===========================================================================
// Serialization schema version (v1 was implicitly 1, 11-wide Vector + 3-wide ID;
// v2 writes 2, with the version in dataID(3) -- NOT data(0), which holds model).
static const int LRB_SCHEMA_VERSION = 2;
static const int LRB_DATA_LEN       = 28;   // v2 Vector width

int LadrunoRebarBuckling::sendSelf(int cTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(4);
  dataID(0) = this->getTag();
  dataID(1) = theBar->getClassTag();
  int matDbTag = theBar->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theBar->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  dataID(3) = LRB_SCHEMA_VERSION;                 // schema version (read before the Vector)
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send ID\n";
    return -1;
  }

  static Vector data(LRB_DATA_LEN);
  // 0..10 : v1 layout, UNCHANGED
  data(0)  = model;        data(1)  = lsr;          data(2)  = alpha;
  data(3)  = fy;           data(4)  = E;
  data(5)  = CmaxTenStrain; data(6) = CmaxTenStress;
  data(7)  = CfStarL;      data(8)  = CfStarLcross;
  data(9)  = gaReduction;  data(10) = gaFsuFrac;
  // 11.. : v2 additions
  data(11) = restraightenMode;  data(12) = restraightenC;
  data(13) = Cbranch;      data(14) = Cstrain;      data(15) = Cstress;
  data(16) = Cg;           data(17) = Cg_law;
  data(18) = Ce_rev;       data(19) = Cs_rev;       data(20) = Ce_cross_rs;
  data(21) = Cf_bare_cross; data(22) = CL_rs;
  data(23) = Ce_reentry;   data(24) = Cg_law_reentry; data(25) = Cg_reentry;
  data(26) = Ce_bow;       data(27) = Cg_bow;
  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send Vector\n";
    return -2;
  }

  if (theBar->sendSelf(cTag, theChannel) < 0) {
    opserr << "LadrunoRebarBuckling::sendSelf - failed to send the bar\n";
    return -3;
  }
  return 0;
}

int LadrunoRebarBuckling::recvSelf(int cTag, Channel& theChannel,
                                   FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(4);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv ID "
           << "(a v1 3-wide stream is not readable by this v2 build)\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  const int version = int(dataID(3));
  if (version != LRB_SCHEMA_VERSION) {
    opserr << "LadrunoRebarBuckling::recvSelf - unsupported serialization version "
           << version << " (need " << LRB_SCHEMA_VERSION << ")\n";
    return -3;
  }

  if (theBar == 0) {
    int matClassTag = int(dataID(1));
    theBar = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theBar == 0) {
      opserr << "LadrunoRebarBuckling::recvSelf - failed to create bar with "
             << "classTag " << matClassTag << endln;
      return -2;
    }
  }
  theBar->setDbTag(dataID(2));

  static Vector data(LRB_DATA_LEN);
  if (theChannel.recvVector(dbTag, cTag, data) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv Vector\n";
    return -4;
  }
  model = int(data(0));    lsr   = data(1);          alpha = data(2);
  fy    = data(3);         E     = data(4);
  CmaxTenStrain = data(5); CmaxTenStress = data(6);
  CfStarL       = data(7); CfStarLcross  = data(8);
  gaReduction   = data(9); gaFsuFrac     = data(10);
  restraightenMode = int(data(11)); restraightenC = data(12);
  Cbranch = int(data(13)); Cstrain = data(14);       Cstress = data(15);
  Cg      = data(16);      Cg_law  = data(17);
  Ce_rev  = data(18);      Cs_rev  = data(19);        Ce_cross_rs = data(20);
  Cf_bare_cross = data(21); CL_rs  = data(22);
  Ce_reentry = data(23);   Cg_law_reentry = data(24); Cg_reentry = data(25);
  Ce_bow  = data(26);      Cg_bow  = data(27);

  // trial = committed (full state recomputed on next setTrialStrain)
  TmaxTenStrain = CmaxTenStrain; TmaxTenStress = CmaxTenStress;
  TfStarL = CfStarL;             TfStarLcross  = CfStarLcross;
  Tbranch = Cbranch;
  Tg = Cg;                       Tg_law = Cg_law;
  Te_rev = Ce_rev;               Ts_rev = Cs_rev;
  Te_cross_rs = Ce_cross_rs;     Tf_bare_cross = Cf_bare_cross;  TL_rs = CL_rs;
  Te_reentry = Ce_reentry;       Tg_law_reentry = Cg_law_reentry; Tg_reentry = Cg_reentry;
  Te_bow = Ce_bow;               Tg_bow = Cg_bow;
  Tstrain = Cstrain; Tstress = Cstress; Ttangent = E; Tr = 1.0;

  if (theBar->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoRebarBuckling::recvSelf - failed to recv the bar\n";
    return -5;
  }
  return 0;
}

// ===========================================================================
//  Print
// ===========================================================================
void LadrunoRebarBuckling::Print(OPS_Stream& s, int flag)
{
  const char* brName = (Cbranch == BR_BUCKLING) ? "BUCKLING"
                     : (Cbranch == BR_RESTRAIGHTEN) ? "RESTRAIGHTEN" : "PASS";
  const char* rsName = (restraightenMode == RS_LAMBDA) ? "lambda" : "c";

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoRebarBuckling\", ";
    s << "\"material\": \"" << theBar->getTag() << "\", ";
    s << "\"model\": \"" << (model == MODEL_DM ? "dm" : "ga") << "\", ";
    s << "\"lsr\": " << lsr << ", ";
    s << "\"alpha\": " << alpha << ", ";
    s << "\"reduction\": " << gaReduction << ", ";
    s << "\"fsufrac\": " << gaFsuFrac << ", ";
    s << "\"fy\": " << fy << ", ";
    s << "\"E\": " << E << ", ";
    s << "\"restraighten\": \"" << rsName << "\", ";
    s << "\"restraightenC\": " << restraightenC << ", ";
    s << "\"branch\": \"" << brName << "\"}";
    return;
  }

  s << "LadrunoRebarBuckling, tag: " << this->getTag() << endln;
  s << "  wrapped material tag: " << theBar->getTag() << endln;
  s << "  model: " << (model == MODEL_DM ? "Dhakal-Maekawa" : "Gomes-Appleton") << endln;
  if (model == MODEL_DM)
    s << "  lsr (s/d): " << lsr << "  alpha: " << alpha << "  fy: " << fy
      << "  E: " << E << endln;
  else
    s << "  lsr (s/d): " << lsr << "  reduction: " << gaReduction
      << "  fsufrac: " << gaFsuFrac << "  E: " << E << endln;
  s << "  restraighten: " << rsName;
  if (restraightenMode == RS_C) s << " (c=" << restraightenC << ")";
  s << "  branch: " << brName << endln;
  if (Cbranch != BR_PASS)
    s << "    g=" << Cg << "  g_law=" << Cg_law
      << "  e_rev=" << Ce_rev << "  s_rev=" << Cs_rev
      << "  e_cross_rs=" << Ce_cross_rs << "  L_rs=" << CL_rs
      << "  e_bow=" << Ce_bow << endln;
  s << "  wrapped material:" << endln;
  theBar->Print(s, flag);
}

// ===========================================================================
//  parameters
// ===========================================================================
int LadrunoRebarBuckling::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "lsr") == 0) {
    param.setValue(lsr);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "alpha") == 0) {
    param.setValue(alpha);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0], "restraightenC") == 0) {
    param.setValue(restraightenC);
    return param.addObject(3, this);
  }
  // forward everything else to the wrapped bar (so element setParameter
  // can still reach the core material's parameters)
  return (theBar != 0) ? theBar->setParameter(argv, argc, param) : -1;
}

int LadrunoRebarBuckling::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case 1: lsr           = info.theDouble; return 0;
  case 2: alpha         = info.theDouble; return 0;
  case 3: restraightenC = info.theDouble; return 0;
  default: return -1;
  }
}

// ===========================================================================
//  responses
// ===========================================================================
Response* LadrunoRebarBuckling::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc == 0)
    return 0;

  // Responses we don't recognize are forwarded to the wrapped bar so core
  // quantities (plasticStrain, backStress, ...) stay reachable. Decide FIRST,
  // so forwarding does not emit a spurious empty UniaxialMaterialOutput block.
  const bool mine =
    (strcmp(argv[0], "stress") == 0) ||
    (strcmp(argv[0], "tangent") == 0) ||
    (strcmp(argv[0], "strain") == 0) ||
    (strcmp(argv[0], "stressStrain") == 0) ||
    (strcmp(argv[0], "stressANDstrain") == 0) ||
    (strcmp(argv[0], "reduction") == 0) ||
    (strcmp(argv[0], "buckling") == 0);

  if (!mine)
    return (theBar != 0) ? theBar->setResponse(argv, argc, s) : 0;

  Response* theResponse = 0;

  s.tag("UniaxialMaterialOutput");
  s.attr("matType", this->getClassType());
  s.attr("matTag", this->getTag());

  if (strcmp(argv[0], "stress") == 0) {
    s.tag("ResponseType", "sigma11");
    theResponse = new MaterialResponse(this, 1, this->getStress());
  }
  else if (strcmp(argv[0], "tangent") == 0) {
    s.tag("ResponseType", "C11");
    theResponse = new MaterialResponse(this, 2, this->getTangent());
  }
  else if (strcmp(argv[0], "strain") == 0) {
    s.tag("ResponseType", "eps11");
    theResponse = new MaterialResponse(this, 3, this->getStrain());
  }
  else if ((strcmp(argv[0], "stressStrain") == 0) ||
           (strcmp(argv[0], "stressANDstrain") == 0)) {
    s.tag("ResponseType", "sig11");
    s.tag("ResponseType", "eps11");
    theResponse = new MaterialResponse(this, 4, Vector(2));
  }
  else if ((strcmp(argv[0], "reduction") == 0) ||
           (strcmp(argv[0], "buckling") == 0)) {
    s.tag("ResponseType", "r");
    theResponse = new MaterialResponse(this, 5, Tr);
  }

  s.endTag();
  return theResponse;
}

int LadrunoRebarBuckling::getResponse(int responseID, Information& matInfo)
{
  static Vector stressStrain(2);
  switch (responseID) {
  case 1: matInfo.setDouble(this->getStress());  return 0;
  case 2: matInfo.setDouble(this->getTangent()); return 0;
  case 3: matInfo.setDouble(this->getStrain());  return 0;
  case 4:
    stressStrain(0) = this->getStress();
    stressStrain(1) = this->getStrain();
    matInfo.setVector(stressStrain);
    return 0;
  case 5: matInfo.setDouble(Tr); return 0;
  default: return -1;
  }
}
