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

// Ladruno: LadrunoHS -- Hardening Soil (Schanz, Vermeer & Bonnier 1999) and
// Hardening Soil with Small-Strain Stiffness (HSsmall, Benz 2007) nDMaterial
// with selectable Lode-angle cross-section (-lodeShape mc vs -lodeShape mn,
// Panteghini & Lagioia 2014).
//
// See LadrunoHS.h, LadrunoHSKernel.h, and Ladruno_implementation/210_ladruno_hs_plan.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoHS.h>
#include <LadrunoMaterialStatus.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

using namespace LadrunoHSKernel;

namespace {

// Ladruno ADR-92 P2-7c pattern: upstream MaterialStageParameter::setDomain()
// exits its element loop after the FIRST matching element, so only eid=1
// receives updateParameter(1, info). Storing the stage keyed by material tag
// (with reference-counted cleanup in ~LadrunoHS() when ops.wipe() destroys all
// instances of a tag) allows every other element's Gauss point to lazily run
// applyStageTransition() at the top of integrate3D() without getCopy() clobbering
// an active stage transition.
std::unordered_map<int, int>& ladrunoHSTagStageMap() {
  static auto* s_map = new std::unordered_map<int, int>();
  return *s_map;
}

std::unordered_map<int, int>& ladrunoHSTagRefCountMap() {
  static auto* s_refs = new std::unordered_map<int, int>();
  return *s_refs;
}

static void* parseLadrunoHSCommon(bool presetSmallStrain) {
  const int minArgs = presetSmallStrain ? 10 : 8;
  if (OPS_GetNumRemainingInputArgs() < minArgs) {
    opserr << "WARNING LadrunoHS: insufficient arguments. Want: "
           << (presetSmallStrain
                   ? "nDMaterial LadrunoHSSmall tag E50ref Eoedref Eurref m c phi psi G0ref gamma07 <options...>\n"
                   : "nDMaterial LadrunoHS tag E50ref Eoedref Eurref m c phi psi <options...>\n");
    return nullptr;
  }

  int tag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING LadrunoHS: invalid material tag\n";
    return nullptr;
  }

  double core[7];
  numData = 7;
  if (OPS_GetDoubleInput(&numData, core) < 0) {
    opserr << "WARNING LadrunoHS (tag " << tag
           << "): expected 7 core parameters: E50ref Eoedref Eurref m c phi psi\n";
    return nullptr;
  }

  HSParams p;
  p.E50_ref  = core[0];
  p.Eoed_ref = core[1];
  p.Eur_ref  = core[2];
  p.m        = core[3];
  p.c        = core[4];
  p.phi_deg  = core[5];
  p.psi_deg  = core[6];

  if (presetSmallStrain) {
    double ss[2];
    numData = 2;
    if (OPS_GetDoubleInput(&numData, ss) < 0) {
      opserr << "WARNING LadrunoHSSmall (tag " << tag
             << "): expected G0ref gamma07 after psi\n";
      return nullptr;
    }
    p.use_small_strain = true;
    p.G0_ref           = ss[0];
    p.gamma_07         = ss[1];
  }

  double rho = 0.0;
  int stage = 1;
  bool implex = false;
  bool hasSig0 = false;
  double sig0Comp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (flag == nullptr) {
      opserr << "WARNING LadrunoHS (tag " << tag << "): null option flag\n";
      return nullptr;
    }

    if (std::strcmp(flag, "-nu") == 0 || std::strcmp(flag, "-nu_ur") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.nu_ur) < 0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -nu wants nu_ur\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-pref") == 0 || std::strcmp(flag, "-p_ref") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.p_ref) < 0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -pref wants p_ref\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-Rf") == 0 || std::strcmp(flag, "-rf") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.Rf) < 0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -Rf wants Rf\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-K0nc") == 0 || std::strcmp(flag, "-k0nc") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.K0nc) < 0 || p.K0nc <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -K0nc wants K0nc > 0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-sigmat") == 0 || std::strcmp(flag, "-sigma_t") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.sigma_t) < 0 || p.sigma_t < 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -sigmat wants sigma_t >= 0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-OCR") == 0 || std::strcmp(flag, "-ocr") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.OCR) < 0 || p.OCR < 1.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -OCR wants OCR >= 1.0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-pop") == 0 || std::strcmp(flag, "-POP") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.pop) < 0 || p.pop < 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -pop wants pop >= 0.0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-e0") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.e0) < 0 || p.e0 <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -e0 wants e0 > 0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-emax") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.emax) < 0 || p.emax <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -emax wants emax > 0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-smallStrain") == 0 || std::strcmp(flag, "-hssmall") == 0) {
      double ss[2];
      numData = 2;
      if (OPS_GetDoubleInput(&numData, ss) < 0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -smallStrain wants G0ref gamma07\n";
        return nullptr;
      }
      p.use_small_strain = true;
      p.G0_ref           = ss[0];
      p.gamma_07         = ss[1];
    } else if (std::strcmp(flag, "-lodeShape") == 0 || std::strcmp(flag, "-lode") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -lodeShape wants mc or mn\n";
        return nullptr;
      }
      const char* s = OPS_GetString();
      if (std::strcmp(s, "mc") == 0 || std::strcmp(s, "MC") == 0 ||
          std::strcmp(s, "mohr-coulomb") == 0) {
        p.lode_shape = LODE_MC;
      } else if (std::strcmp(s, "mn") == 0 || std::strcmp(s, "MN") == 0 ||
                 std::strcmp(s, "matsuoka-nakai") == 0) {
        p.lode_shape = LODE_MN;
      } else {
        opserr << "WARNING LadrunoHS (tag " << tag << "): unknown -lodeShape '" << s
               << "' (expected 'mc' or 'mn')\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-stiffStress") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -stiffStress wants sigma3 or p\n";
        return nullptr;
      }
      const char* s = OPS_GetString();
      if (std::strcmp(s, "sigma3") == 0 || std::strcmp(s, "sig3") == 0) {
        p.stiff_stress = STIFF_SIGMA3;
      } else if (std::strcmp(s, "p") == 0 || std::strcmp(s, "mean") == 0) {
        p.stiff_stress = STIFF_MEAN_P;
      } else {
        opserr << "WARNING LadrunoHS (tag " << tag << "): unknown -stiffStress '" << s
               << "' (expected 'sigma3' or 'p')\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-alphaCap") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.alpha_cap) < 0 || p.alpha_cap <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -alphaCap wants alpha > 0\n";
        return nullptr;
      }
      p.user_alpha_cap = p.alpha_cap;
    } else if (std::strcmp(flag, "-betaCap") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.beta_cap) < 0 || p.beta_cap <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -betaCap wants beta > 0\n";
        return nullptr;
      }
      p.user_beta_cap = p.beta_cap;
    } else if (std::strcmp(flag, "-pMinRatio") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p.p_min_ratio) < 0 || p.p_min_ratio <= 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -pMinRatio wants ratio > 0\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-sig0") == 0) {
      double p0 = 0.0;
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &p0) < 0 || p0 < 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -sig0 wants compressive pressure p0 >= 0\n";
        return nullptr;
      }
      hasSig0 = true;
      sig0Comp[0] = sig0Comp[1] = sig0Comp[2] = p0;
      sig0Comp[3] = sig0Comp[4] = sig0Comp[5] = 0.0;
    } else if (std::strcmp(flag, "-sig03") == 0) {
      double s3[3];
      numData = 3;
      if (OPS_GetDoubleInput(&numData, s3) < 0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -sig03 wants s11 s22 s33 (compressive > 0)\n";
        return nullptr;
      }
      hasSig0 = true;
      sig0Comp[0] = s3[0];
      sig0Comp[1] = s3[1];
      sig0Comp[2] = s3[2];
      sig0Comp[3] = sig0Comp[4] = sig0Comp[5] = 0.0;
    } else if (std::strcmp(flag, "-stage") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &stage) < 0 || (stage != 0 && stage != 1)) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -stage wants 0 (elastic) or 1 (plastic)\n";
        return nullptr;
      }
    } else if (std::strcmp(flag, "-implex") == 0) {
      implex = true;
    } else if (std::strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0 || rho < 0.0) {
        opserr << "WARNING LadrunoHS (tag " << tag << "): -rho wants rho >= 0\n";
        return nullptr;
      }
    } else {
      opserr << "WARNING LadrunoHS (tag " << tag << "): unknown option '" << flag << "'\n";
      return nullptr;
    }
  }

  if (!p.finalize()) {
    opserr << "WARNING LadrunoHS (tag " << tag
           << "): invalid material parameter combination (check E50ref, Eoedref, Eurref, m, phi, nu, Rf, G0ref, gamma07)\n";
    return nullptr;
  }

  ladrunoHSTagStageMap()[tag] = stage;
  return new LadrunoHS(tag, p, rho, LadrunoHS::DIM_3D, stage, implex, hasSig0, sig0Comp);
}

}  // namespace

void* OPS_LadrunoHS(void) {
  return parseLadrunoHSCommon(false);
}

void* OPS_LadrunoHSSmall(void) {
  return parseLadrunoHSCommon(true);
}

// ===========================================================================
//  Constructors & Destructor
// ===========================================================================
LadrunoHS::LadrunoHS()
  : NDMaterial(0, ND_TAG_LadrunoHS),
    params_(),
    rho_(0.0),
    dim_(DIM_3D),
    ncomp_(6),
    condense_(false),
    stage_(1),
    implex_(false),
    has_sig0_(false),
    regime_n_(REGIME_ELASTIC),
    regime_trial_(REGIME_ELASTIC),
    cEps22_comp_(0.0),
    refused_(false),
    registered_in_tag_map_(false) {
  params_.finalize();
  for (int i = 0; i < 6; ++i) {
    sig0_comp_[i]      = 0.0;
    dsig_p_n_[i]       = 0.0;
    dsig_p_trial_[i]   = 0.0;
    sig_impl_trial_[i] = 0.0;
  }
  this->setupDim();
  this->revertToStart();
}

LadrunoHS::LadrunoHS(int tag,
                     const HSParams& params,
                     double rho,
                     int dimMode,
                     int stage,
                     bool implex,
                     bool hasSig0,
                     const double sig0Comp[6])
  : NDMaterial(tag, ND_TAG_LadrunoHS),
    params_(params),
    rho_(rho),
    dim_(dimMode),
    ncomp_(6),
    condense_(false),
    stage_(stage),
    implex_(implex),
    has_sig0_(hasSig0),
    regime_n_(REGIME_ELASTIC),
    regime_trial_(REGIME_ELASTIC),
    cEps22_comp_(0.0),
    refused_(false),
    registered_in_tag_map_(true) {
  params_.finalize();
  for (int i = 0; i < 6; ++i) {
    sig0_comp_[i]      = (hasSig0 && sig0Comp != nullptr) ? sig0Comp[i] : 0.0;
    dsig_p_n_[i]       = 0.0;
    dsig_p_trial_[i]   = 0.0;
    sig_impl_trial_[i] = 0.0;
  }
  this->setupDim();
  if (condense_) {
    sig0_comp_[2] = 0.0;
  }
  auto& refs = ladrunoHSTagRefCountMap();
  if (refs[this->getTag()] == 0) {
    ladrunoHSTagStageMap()[this->getTag()] = stage_;
  }
  refs[this->getTag()] += 1;
  this->revertToStart();
}

LadrunoHS::~LadrunoHS() {
  if (registered_in_tag_map_) {
    int t = this->getTag();
    auto& refs = ladrunoHSTagRefCountMap();
    auto it = refs.find(t);
    if (it != refs.end()) {
      it->second -= 1;
      if (it->second <= 0) {
        refs.erase(it);
        ladrunoHSTagStageMap().erase(t);
      }
    }
  }
}

void LadrunoHS::setupDim(void) {
  // Full 6-component ordering: 0:11, 1:22, 2:33, 3:12, 4:23, 5:31
  switch (dim_) {
    case DIM_PSTRAIN:
      ncomp_ = 3;
      vmap_[0] = 0; vmap_[1] = 1; vmap_[2] = 3;
      condense_ = false;
      break;
    case DIM_AXISYM:
      ncomp_ = 4;
      vmap_[0] = 0; vmap_[1] = 1; vmap_[2] = 2; vmap_[3] = 3;
      condense_ = false;
      break;
    case DIM_PLATEFIBER:
      ncomp_ = 5;
      vmap_[0] = 0; vmap_[1] = 1; vmap_[2] = 3; vmap_[3] = 4; vmap_[4] = 5;
      condense_ = true;
      break;
    case DIM_PSTRESS:
      ncomp_ = 3;
      vmap_[0] = 0; vmap_[1] = 1; vmap_[2] = 3;
      condense_ = true;
      break;
    case DIM_3D:
    default:
      ncomp_ = 6;
      for (int a = 0; a < 6; ++a) vmap_[a] = a;
      condense_ = false;
      break;
  }
  stressOut_.resize(ncomp_);
  strainOut_.resize(ncomp_);
  tangentOut_.resize(ncomp_, ncomp_);
  initTangentOut_.resize(ncomp_, ncomp_);
}

void LadrunoHS::condenseTangent(double D[6][6]) const {
  double d22 = D[2][2];
  if (std::abs(d22) < 1.0e-30) return;
  double col[6], row[6];
  for (int i = 0; i < 6; ++i) {
    col[i] = D[i][2];
    row[i] = D[2][i];
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      D[i][j] -= col[i] * row[j] / d22;
    }
  }
}

void LadrunoHS::buildInitialTangent6x6(double Kt[6][6]) const {
  double s0 = 1.0;
  if (has_sig0_) {
    double p_n = 0.0, q_n = 0.0, x_n = 1.0, Th_n = 1.0;
    double dq[6], dTh[6];
    invariants_and_lode_grad(params_, sig0_comp_, p_n, q_n, x_n, Th_n, dq, dTh);
    double ds_dp = 0.0, ds_dq = 0.0;
    stiffness_scale_and_derivs(params_, p_n, Th_n * q_n, s0, ds_dp, ds_dq);
  }
  double K0 = params_.Kur_ref * s0;
  double G0 = (params_.use_small_strain && params_.G0_ref > params_.Gur_ref)
                  ? (params_.G0_ref * s0)
                  : (params_.Gur_ref * s0);
  build_isotropic_elastic_tangent(K0, G0, Kt);
  if (condense_) {
    this->condenseTangent(Kt);
  }
}

// ===========================================================================
//  Constitutive Integration (Stage 0 Elastic / Stage 1 HS / IMPL-EX)
// ===========================================================================
void LadrunoHS::applyStageTransition(int newStage) {
  if (stage_ == 0 && newStage == 1) {
    // Transition from elastic gravity stage (0) to elastoplastic HS (1):
    // Initialize p_p and gamma_p from the committed gravity stress st_n_.sig
    // while preserving the committed strain st_n_.eps, deviatoric strain e_dev, and void_ratio!
    double eps_saved[6];
    for (int i = 0; i < 6; ++i) eps_saved[i] = st_n_.eps[i];
    double e_saved = st_n_.void_ratio;
    initialize_state_from_stress(params_, st_n_.sig, st_n_);
    for (int i = 0; i < 6; ++i) {
      st_n_.eps[i]       = eps_saved[i];
      sig_impl_trial_[i] = st_n_.sig[i];
      dsig_p_n_[i]       = 0.0;
      dsig_p_trial_[i]   = 0.0;
    }
    double eps_v_mean = (eps_saved[0] + eps_saved[1] + eps_saved[2]) / 3.0;
    st_n_.e_dev[0][0] = eps_saved[0] - eps_v_mean;
    st_n_.e_dev[1][1] = eps_saved[1] - eps_v_mean;
    st_n_.e_dev[2][2] = eps_saved[2] - eps_v_mean;
    st_n_.e_dev[0][1] = st_n_.e_dev[1][0] = 0.5 * eps_saved[3];
    st_n_.e_dev[1][2] = st_n_.e_dev[2][1] = 0.5 * eps_saved[4];
    st_n_.e_dev[0][2] = st_n_.e_dev[2][0] = 0.5 * eps_saved[5];
    st_n_.void_ratio = e_saved;
    st_trial_ = st_n_;
  }
  stage_ = newStage;
}

int LadrunoHS::integrate3D(const double epsComp[6]) {
  auto itStage = ladrunoHSTagStageMap().find(this->getTag());
  if (itStage != ladrunoHSTagStageMap().end() && itStage->second != stage_) {
    this->applyStageTransition(itStage->second);
  }
  if (stage_ == 0) {
    // Stage 0: pressure-independent linear elastic gravity initialization
    for (int i = 0; i < 6; ++i) {
      if (!std::isfinite(epsComp[i])) return -1;
    }
    st_trial_ = st_n_;
    double deps[6];
    for (int i = 0; i < 6; ++i) {
      st_trial_.eps[i] = epsComp[i];
      deps[i]          = epsComp[i] - st_n_.eps[i];
    }
    double deps_v = deps[0] + deps[1] + deps[2];
    st_trial_.void_ratio = st_n_.void_ratio - (1.0 + st_n_.void_ratio) * deps_v;

    double eps_v_mean = (epsComp[0] + epsComp[1] + epsComp[2]) / 3.0;
    st_trial_.e_dev[0][0] = epsComp[0] - eps_v_mean;
    st_trial_.e_dev[1][1] = epsComp[1] - eps_v_mean;
    st_trial_.e_dev[2][2] = epsComp[2] - eps_v_mean;
    st_trial_.e_dev[0][1] = st_trial_.e_dev[1][0] = 0.5 * epsComp[3];
    st_trial_.e_dev[1][2] = st_trial_.e_dev[2][1] = 0.5 * epsComp[4];
    st_trial_.e_dev[0][2] = st_trial_.e_dev[2][0] = 0.5 * epsComp[5];

    build_isotropic_elastic_tangent(params_.Kur_ref, params_.Gur_ref, Dtan_);
    for (int i = 0; i < 6; ++i) {
      double dsig_i = 0.0;
      for (int j = 0; j < 6; ++j) dsig_i += Dtan_[i][j] * deps[j];
      st_trial_.sig[i]   = st_n_.sig[i] + dsig_i;
      sig_impl_trial_[i] = st_trial_.sig[i];
      dsig_p_trial_[i]   = 0.0;
    }
    regime_trial_ = REGIME_ELASTIC;
    return 0;
  }

  // Stage 1: Elastoplastic Hardening Soil (HS / HSsmall)
  double deps_norm_sq = 0.0;
  for (int i = 0; i < 6; ++i) {
    double di = epsComp[i] - st_n_.eps[i];
    deps_norm_sq += di * di;
  }
  ActiveRegime reg = return_map_3d_comp(params_, st_n_, epsComp, st_trial_, Dtan_);
  if (reg == REGIME_REFUSED) {
    return -1;
  }
  regime_trial_ = (reg == REGIME_ELASTIC && deps_norm_sq <= 1.0e-30) ? regime_n_ : reg;

  for (int i = 0; i < 6; ++i) {
    sig_impl_trial_[i] = st_trial_.sig[i];
    dsig_p_trial_[i]   = 0.0;
  }

  if (implex_) {
    // Oliver et al. (2008) IMPL-EX:
    //   1. Evaluate SPD paraelastic tangent C_e(sig_n) at the committed stress st_n_.sig.
    //   2. Compute elastic predictor stress sig_tr = sig_n + C_e(sig_n) : deps.
    //   3. Record current step's plastic stress correction dsig_p_trial_ = sig_tr - sig_impl_trial_.
    //   4. Expose extrapolated stress sig_tilde = sig_tr - dsig_p_n_ as st_trial_.sig for global
    //      equilibrium while committing sig_impl_trial_ into st_n_.sig at commitState().
    double norm_sig_n = 0.0;
    for (int i = 0; i < 6; ++i) norm_sig_n += st_n_.sig[i] * st_n_.sig[i];
    double s_n = 1.0;
    if (norm_sig_n > 1.0e-24) {
      double p_n = 0.0, q_n = 0.0, x_n = 1.0, Th_n = 1.0;
      double dq[6], dTh[6];
      invariants_and_lode_grad(params_, st_n_.sig, p_n, q_n, x_n, Th_n, dq, dTh);
      double ds1 = 0.0, ds2 = 0.0;
      stiffness_scale_and_derivs(params_, p_n, Th_n * q_n, s_n, ds1, ds2);
    }
    build_isotropic_elastic_tangent(params_.Kur_ref * s_n, params_.Gur_ref * s_n, Dtan_);
    double deps[6];
    for (int i = 0; i < 6; ++i) deps[i] = epsComp[i] - st_n_.eps[i];
    for (int i = 0; i < 6; ++i) {
      double dsig_e = 0.0;
      for (int j = 0; j < 6; ++j) dsig_e += Dtan_[i][j] * deps[j];
      double sig_tr_i  = st_n_.sig[i] + dsig_e;
      dsig_p_trial_[i] = sig_tr_i - sig_impl_trial_[i];
      st_trial_.sig[i] = sig_tr_i - dsig_p_n_[i];
    }
  }

  for (int i = 0; i < 6; ++i) {
    if (!std::isfinite(st_trial_.sig[i])) return -1;
    for (int j = 0; j < 6; ++j) {
      if (!std::isfinite(Dtan_[i][j])) return -1;
    }
  }
  return 0;
}

int LadrunoHS::setTrialStrain(const Vector& v) {
  if (v.Size() != ncomp_) {
    refused_ = true;
    return LADRUNO_MATERIAL_REFUSED;
  }

  double epsComp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (condense_) {
    epsComp[2] = st_trial_.eps[2];
  }
  // OpenSees continuum strain is tension-positive; kernel is compression-positive:
  for (int a = 0; a < ncomp_; ++a) {
    double val = v(a);
    if (!std::isfinite(val)) {
      refused_ = true;
      return LADRUNO_MATERIAL_REFUSED;
    }
    epsComp[vmap_[a]] = -val;
  }

  if (!condense_) {
    if (this->integrate3D(epsComp) < 0) {
      refused_ = true;
      return LADRUNO_MATERIAL_REFUSED;
    }
    refused_ = false;
    return 0;
  }

  // Condensed modes (DIM_PSTRESS / DIM_PLATEFIBER): enforce sig_22 == 0 via Newton on epsComp[2]:
  const int maxIt = 25;
  bool conv = false;
  for (int it = 0; it < maxIt; ++it) {
    if (this->integrate3D(epsComp) < 0) {
      refused_ = true;
      return LADRUNO_MATERIAL_REFUSED;
    }
    double d22 = Dtan_[2][2];
    double smag = 0.0;
    for (int i = 0; i < 6; ++i) smag += st_trial_.sig[i] * st_trial_.sig[i];
    smag = std::sqrt(smag);
    double tol22 = 1.0e-10 * std::max(1.0, smag);
    if (std::abs(st_trial_.sig[2]) <= tol22) {
      conv = true;
      break;
    }
    if (std::abs(d22) < 1.0e-30) break;
    epsComp[2] -= st_trial_.sig[2] / d22;
  }
  if (!conv) {
    refused_ = true;
    return LADRUNO_MATERIAL_REFUSED;
  }
  this->condenseTangent(Dtan_);
  refused_ = false;
  return 0;
}

int LadrunoHS::setTrialStrain(const Vector& v, const Vector&) {
  return this->setTrialStrain(v);
}

int LadrunoHS::setTrialStrainIncr(const Vector& v) {
  if (v.Size() != ncomp_) {
    refused_ = true;
    return LADRUNO_MATERIAL_REFUSED;
  }
  Vector ne(ncomp_);
  for (int a = 0; a < ncomp_; ++a) {
    ne(a) = -st_trial_.eps[vmap_[a]] + v(a);
  }
  return this->setTrialStrain(ne);
}

int LadrunoHS::setTrialStrainIncr(const Vector& v, const Vector&) {
  return this->setTrialStrainIncr(v);
}

// ===========================================================================
//  Element-Facing Stress, Strain & Tangent Getters (OpenSees Sign Convention)
// ===========================================================================
const Vector& LadrunoHS::getStress(void) {
  for (int a = 0; a < ncomp_; ++a) {
    stressOut_(a) = -st_trial_.sig[vmap_[a]];
  }
  return stressOut_;
}

double LadrunoHS::getStressZZ(void) {
  if (dim_ == DIM_PSTRAIN) {
    return -st_trial_.sig[2];
  }
  return NDMaterial::getStressZZ();
}

const Vector& LadrunoHS::getStrain(void) {
  for (int a = 0; a < ncomp_; ++a) {
    strainOut_(a) = -st_trial_.eps[vmap_[a]];
  }
  return strainOut_;
}

const Matrix& LadrunoHS::getTangent(void) {
  for (int a = 0; a < ncomp_; ++a) {
    for (int b = 0; b < ncomp_; ++b) {
      tangentOut_(a, b) = Dtan_[vmap_[a]][vmap_[b]];
    }
  }
  return tangentOut_;
}

const Matrix& LadrunoHS::getInitialTangent(void) {
  double Kt[6][6];
  this->buildInitialTangent6x6(Kt);
  for (int a = 0; a < ncomp_; ++a) {
    for (int b = 0; b < ncomp_; ++b) {
      initTangentOut_(a, b) = Kt[vmap_[a]][vmap_[b]];
    }
  }
  return initTangentOut_;
}

const char* LadrunoHS::getType(void) const {
  switch (dim_) {
    case DIM_PSTRAIN:    return "PlaneStrain";
    case DIM_AXISYM:     return "AxiSymmetric";
    case DIM_PLATEFIBER: return "PlateFiber";
    case DIM_PSTRESS:    return "PlaneStress";
    default:             return "ThreeDimensional";
  }
}

// ===========================================================================
//  State Commit / Revert Cycle
// ===========================================================================
int LadrunoHS::commitState(void) {
  if (refused_) {
    ladrunoNoteCommitRefusal();
    return LADRUNO_MATERIAL_REFUSED;
  }
  st_n_        = st_trial_;
  if (implex_) {
    for (int i = 0; i < 6; ++i) {
      st_n_.sig[i] = sig_impl_trial_[i];
      dsig_p_n_[i] = dsig_p_trial_[i];
    }
  } else {
    for (int i = 0; i < 6; ++i) {
      dsig_p_n_[i] = 0.0;
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      Dtan_n_[i][j] = Dtan_[i][j];
    }
  }
  regime_n_    = regime_trial_;
  cEps22_comp_ = st_trial_.eps[2];
  return 0;
}

int LadrunoHS::revertToLastCommit(void) {
  st_trial_     = st_n_;
  regime_trial_ = regime_n_;
  if (condense_) {
    st_trial_.eps[2] = cEps22_comp_;
  }
  for (int i = 0; i < 6; ++i) {
    sig_impl_trial_[i] = st_n_.sig[i];
    dsig_p_trial_[i]   = dsig_p_n_[i];
    for (int j = 0; j < 6; ++j) {
      Dtan_[i][j] = Dtan_n_[i][j];
    }
  }
  refused_ = false;
  return 0;
}

int LadrunoHS::revertToStart(void) {
  if (has_sig0_) {
    initialize_state_from_stress(params_, sig0_comp_, st_n_);
  } else {
    // Initialize preconsolidation pressure p_p at reference K0nc stress while keeping sig = 0:
    const double sigRef[6] = {
        params_.p_ref,
        params_.p_ref * params_.K0nc,
        params_.p_ref * params_.K0nc,
        0.0, 0.0, 0.0
    };
    initialize_state_from_stress(params_, sigRef, st_n_);
    for (int i = 0; i < 6; ++i) st_n_.sig[i] = 0.0;
    st_n_.gamma_p = 0.0;
  }
  st_trial_     = st_n_;
  regime_n_     = REGIME_ELASTIC;
  regime_trial_ = REGIME_ELASTIC;
  cEps22_comp_  = 0.0;
  refused_      = false;
  for (int i = 0; i < 6; ++i) {
    dsig_p_n_[i]       = 0.0;
    dsig_p_trial_[i]   = 0.0;
    sig_impl_trial_[i] = st_n_.sig[i];
  }
  this->buildInitialTangent6x6(Dtan_);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      Dtan_n_[i][j] = Dtan_[i][j];
    }
  }
  return 0;
}

// ===========================================================================
//  Copies (Propagating State & Refusal Latch per Quirk Rule)
// ===========================================================================
NDMaterial* LadrunoHS::getCopy(void) {
  auto itStage = ladrunoHSTagStageMap().find(this->getTag());
  int activeStage = (itStage != ladrunoHSTagStageMap().end()) ? itStage->second : this->stage_;
  LadrunoHS* copy = new LadrunoHS(
      this->getTag(), params_, rho_, dim_, activeStage, implex_, has_sig0_, sig0_comp_);
  copy->st_n_         = this->st_n_;
  copy->st_trial_     = this->st_trial_;
  copy->regime_n_     = this->regime_n_;
  copy->regime_trial_ = this->regime_trial_;
  copy->cEps22_comp_  = this->cEps22_comp_;
  copy->refused_      = this->refused_;
  for (int i = 0; i < 6; ++i) {
    copy->dsig_p_n_[i]       = this->dsig_p_n_[i];
    copy->dsig_p_trial_[i]   = this->dsig_p_trial_[i];
    copy->sig_impl_trial_[i] = this->sig_impl_trial_[i];
    for (int j = 0; j < 6; ++j) {
      copy->Dtan_[i][j]   = this->Dtan_[i][j];
      copy->Dtan_n_[i][j] = this->Dtan_n_[i][j];
    }
  }
  return copy;
}

NDMaterial* LadrunoHS::getCopy(const char* type) {
  int d = -1;
  if (std::strcmp(type, "ThreeDimensional") == 0 || std::strcmp(type, "3D") == 0) {
    d = DIM_3D;
  } else if (std::strcmp(type, "PlaneStrain") == 0 || std::strcmp(type, "PlaneStrain2D") == 0) {
    d = DIM_PSTRAIN;
  } else if (std::strcmp(type, "AxiSymmetric") == 0 || std::strcmp(type, "AxiSymmetric2D") == 0) {
    d = DIM_AXISYM;
  } else if (std::strcmp(type, "PlateFiber") == 0) {
    d = DIM_PLATEFIBER;
  } else if (std::strcmp(type, "PlaneStress") == 0 || std::strcmp(type, "PlaneStress2D") == 0) {
    d = DIM_PSTRESS;
  }
  if (d < 0) {
    return NDMaterial::getCopy(type);
  }

  auto itStage = ladrunoHSTagStageMap().find(this->getTag());
  int activeStage = (itStage != ladrunoHSTagStageMap().end()) ? itStage->second : this->stage_;
  LadrunoHS* copy = new LadrunoHS(
      this->getTag(), params_, rho_, d, activeStage, implex_, has_sig0_, sig0_comp_);
  if (copy->condense_ && !this->condense_) {
    // Condensed copy from a 3D prototype: keep the plane-stress / plate-fiber
    // initialization from revertToStart() (where sig0_comp_[2] == 0.0 and st_n_.sig[2] == 0.0)
    // and propagate only the refusal latch.
    copy->refused_ = this->refused_;
    return copy;
  }
  copy->st_n_         = this->st_n_;
  copy->st_trial_     = this->st_trial_;
  copy->regime_n_     = this->regime_n_;
  copy->regime_trial_ = this->regime_trial_;
  copy->cEps22_comp_  = this->cEps22_comp_;
  copy->refused_      = this->refused_;
  for (int i = 0; i < 6; ++i) {
    copy->dsig_p_n_[i]       = this->dsig_p_n_[i];
    copy->dsig_p_trial_[i]   = this->dsig_p_trial_[i];
    copy->sig_impl_trial_[i] = this->sig_impl_trial_[i];
  }
  copy->buildInitialTangent6x6(copy->Dtan_);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      copy->Dtan_n_[i][j] = copy->Dtan_[i][j];
    }
  }
  return copy;
}

// ===========================================================================
//  Wire Serialization (sendSelf / recvSelf)
// ===========================================================================
static constexpr int kLadrunoHSWireSize = 75;

int LadrunoHS::sendSelf(int commitTag, Channel& theChannel) {
  static Vector data(kLadrunoHSWireSize);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = params_.E50_ref;
  data(c++) = params_.Eoed_ref;
  data(c++) = params_.Eur_ref;
  data(c++) = params_.m;
  data(c++) = params_.c;
  data(c++) = params_.phi_deg;
  data(c++) = params_.psi_deg;
  data(c++) = params_.nu_ur;
  data(c++) = params_.p_ref;
  data(c++) = params_.Rf;
  data(c++) = params_.K0nc;
  data(c++) = params_.sigma_t;
  data(c++) = params_.OCR;
  data(c++) = params_.pop;
  data(c++) = params_.e0;
  data(c++) = params_.emax;
  data(c++) = params_.use_small_strain ? 1.0 : 0.0;
  data(c++) = params_.G0_ref;
  data(c++) = params_.gamma_07;
  data(c++) = static_cast<double>(params_.lode_shape);
  data(c++) = static_cast<double>(params_.stiff_stress);
  data(c++) = params_.alpha_cap;
  data(c++) = params_.beta_cap;
  data(c++) = params_.user_alpha_cap;
  data(c++) = params_.user_beta_cap;
  data(c++) = params_.p_min_ratio;
  data(c++) = rho_;
  data(c++) = static_cast<double>(dim_);
  data(c++) = static_cast<double>(stage_);
  data(c++) = implex_ ? 1.0 : 0.0;
  data(c++) = has_sig0_ ? 1.0 : 0.0;
  for (int i = 0; i < 6; ++i) data(c++) = sig0_comp_[i];

  // Committed state st_n_ + IMPL-EX correction + flags (37 doubles):
  for (int i = 0; i < 6; ++i) data(c++) = st_n_.sig[i];
  for (int i = 0; i < 6; ++i) data(c++) = st_n_.eps[i];
  data(c++) = st_n_.gamma_p;
  data(c++) = st_n_.p_p;
  data(c++) = st_n_.void_ratio;
  data(c++) = st_n_.e_dev[0][0];
  data(c++) = st_n_.e_dev[1][1];
  data(c++) = st_n_.e_dev[2][2];
  data(c++) = st_n_.e_dev[0][1];
  data(c++) = st_n_.e_dev[1][2];
  data(c++) = st_n_.e_dev[0][2];
  data(c++) = st_n_.H_hist[0][0];
  data(c++) = st_n_.H_hist[1][1];
  data(c++) = st_n_.H_hist[2][2];
  data(c++) = st_n_.H_hist[0][1];
  data(c++) = st_n_.H_hist[1][2];
  data(c++) = st_n_.H_hist[0][2];
  data(c++) = st_n_.is_virgin ? 1.0 : 0.0;
  data(c++) = static_cast<double>(regime_n_);
  data(c++) = cEps22_comp_;
  data(c++) = refused_ ? 1.0 : 0.0;
  for (int i = 0; i < 6; ++i) data(c++) = dsig_p_n_[i];

  if (c != kLadrunoHSWireSize) {
    opserr << "LadrunoHS::sendSelf - wire size mismatch (" << c << " vs " << kLadrunoHSWireSize << ")\n";
    return -1;
  }
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoHS::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoHS::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&) {
  static Vector data(kLadrunoHSWireSize);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoHS::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  if (registered_in_tag_map_) {
    int oldTag = this->getTag();
    auto& refs = ladrunoHSTagRefCountMap();
    auto it = refs.find(oldTag);
    if (it != refs.end()) {
      it->second -= 1;
      if (it->second <= 0) {
        refs.erase(it);
        ladrunoHSTagStageMap().erase(oldTag);
      }
    }
    registered_in_tag_map_ = false;
  }
  this->setTag(static_cast<int>(data(c++)));
  params_.E50_ref          = data(c++);
  params_.Eoed_ref         = data(c++);
  params_.Eur_ref          = data(c++);
  params_.m                = data(c++);
  params_.c                = data(c++);
  params_.phi_deg          = data(c++);
  params_.psi_deg          = data(c++);
  params_.nu_ur            = data(c++);
  params_.p_ref            = data(c++);
  params_.Rf               = data(c++);
  params_.K0nc             = data(c++);
  params_.sigma_t          = data(c++);
  params_.OCR              = data(c++);
  params_.pop              = data(c++);
  params_.e0               = data(c++);
  params_.emax             = data(c++);
  params_.use_small_strain = (data(c++) != 0.0);
  params_.G0_ref           = data(c++);
  params_.gamma_07         = data(c++);
  params_.lode_shape       = static_cast<int>(data(c++));
  params_.stiff_stress     = static_cast<int>(data(c++));
  params_.alpha_cap        = data(c++);
  params_.beta_cap         = data(c++);
  params_.user_alpha_cap   = data(c++);
  params_.user_beta_cap    = data(c++);
  params_.p_min_ratio      = data(c++);
  rho_                     = data(c++);
  dim_                     = static_cast<int>(data(c++));
  stage_                   = static_cast<int>(data(c++));
  implex_                  = (data(c++) != 0.0);
  has_sig0_                = (data(c++) != 0.0);
  for (int i = 0; i < 6; ++i) sig0_comp_[i] = data(c++);

  params_.finalize();
  this->setupDim();

  for (int i = 0; i < 6; ++i) st_n_.sig[i] = data(c++);
  for (int i = 0; i < 6; ++i) st_n_.eps[i] = data(c++);
  st_n_.gamma_p    = data(c++);
  st_n_.p_p        = data(c++);
  st_n_.void_ratio = data(c++);
  st_n_.e_dev[0][0] = data(c++);
  st_n_.e_dev[1][1] = data(c++);
  st_n_.e_dev[2][2] = data(c++);
  st_n_.e_dev[0][1] = st_n_.e_dev[1][0] = data(c++);
  st_n_.e_dev[1][2] = st_n_.e_dev[2][1] = data(c++);
  st_n_.e_dev[0][2] = st_n_.e_dev[2][0] = data(c++);
  st_n_.H_hist[0][0] = data(c++);
  st_n_.H_hist[1][1] = data(c++);
  st_n_.H_hist[2][2] = data(c++);
  st_n_.H_hist[0][1] = st_n_.H_hist[1][0] = data(c++);
  st_n_.H_hist[1][2] = st_n_.H_hist[2][1] = data(c++);
  st_n_.H_hist[0][2] = st_n_.H_hist[2][0] = data(c++);
  st_n_.is_virgin  = (data(c++) != 0.0);
  regime_n_        = static_cast<ActiveRegime>(static_cast<int>(data(c++)));
  cEps22_comp_     = data(c++);
  refused_         = (data(c++) != 0.0);
  for (int i = 0; i < 6; ++i) {
    dsig_p_n_[i]       = data(c++);
    dsig_p_trial_[i]   = dsig_p_n_[i];
    sig_impl_trial_[i] = st_n_.sig[i];
  }
  if (c != kLadrunoHSWireSize) {
    opserr << "LadrunoHS::recvSelf - wire size mismatch (" << c << " vs " << kLadrunoHSWireSize << ")\n";
    return -1;
  }

  st_trial_     = st_n_;
  regime_trial_ = regime_n_;
  ladrunoHSTagRefCountMap()[this->getTag()] += 1;
  registered_in_tag_map_ = true;
  ladrunoHSTagStageMap()[this->getTag()] = stage_;
  this->buildInitialTangent6x6(Dtan_);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      Dtan_n_[i][j] = Dtan_[i][j];
    }
  }
  return 0;
}

// ===========================================================================
//  Parameters (updateMaterialStage / stage / material properties)
// ===========================================================================
int LadrunoHS::setParameter(const char** argv, int argc, Parameter& param) {
  if (argc < 1) return -1;
  const char* a = argv[0];
  if (argc >= 3 && (std::strcmp(argv[0], "material") == 0 || std::strcmp(argv[0], "-material") == 0)) {
    int matTag = std::atoi(argv[1]);
    if (matTag != this->getTag()) return -1;
    a = argv[2];
  } else if (argc >= 2) {
    char* endptr = nullptr;
    long matTag = std::strtol(argv[1], &endptr, 10);
    if (endptr != argv[1] && endptr != nullptr && *endptr == '\0') {
      if (static_cast<int>(matTag) != this->getTag()) return -1;
    }
  }
  if (std::strcmp(a, "updateMaterialStage") == 0) {
    return param.addObject(1, this);
  }
  if (std::strcmp(a, "materialState") == 0 ||
      std::strcmp(a, "stage") == 0 ||
      std::strcmp(a, "firstStage") == 0) {
    return param.addObject(13, this);
  }
  if (std::strcmp(a, "E50") == 0 || std::strcmp(a, "E50ref") == 0) return param.addObject(2, this);
  if (std::strcmp(a, "Eoed") == 0 || std::strcmp(a, "Eoedref") == 0) return param.addObject(3, this);
  if (std::strcmp(a, "Eur") == 0 || std::strcmp(a, "Eurref") == 0) return param.addObject(4, this);
  if (std::strcmp(a, "c") == 0 || std::strcmp(a, "cohesion") == 0) return param.addObject(5, this);
  if (std::strcmp(a, "phi") == 0) return param.addObject(6, this);
  if (std::strcmp(a, "psi") == 0) return param.addObject(7, this);
  if (std::strcmp(a, "G0") == 0 || std::strcmp(a, "G0ref") == 0) return param.addObject(8, this);
  if (std::strcmp(a, "gamma07") == 0) return param.addObject(9, this);
  if (std::strcmp(a, "OCR") == 0) return param.addObject(10, this);
  if (std::strcmp(a, "pop") == 0) return param.addObject(11, this);
  if (std::strcmp(a, "rho") == 0) return param.addObject(12, this);
  return -1;
}

int LadrunoHS::updateParameter(int parameterID, Information& info) {
  auto recalibrateCapAndFinalize = [this]() {
    params_.alpha_cap = params_.user_alpha_cap;
    params_.beta_cap  = params_.user_beta_cap;
    params_.finalize();
  };
  switch (parameterID) {
    case 1:
    case 13: {
      int newStage = (parameterID == 1) ? info.theInt : static_cast<int>(info.theDouble);
      ladrunoHSTagStageMap()[this->getTag()] = newStage;
      this->applyStageTransition(newStage);
      return 0;
    }
    case 2: params_.E50_ref  = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 3: params_.Eoed_ref = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 4: params_.Eur_ref  = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 5: params_.c        = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 6: params_.phi_deg  = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 7: params_.psi_deg  = info.theDouble; recalibrateCapAndFinalize(); return 0;
    case 8: params_.G0_ref   = info.theDouble; params_.finalize(); return 0;
    case 9: params_.gamma_07 = info.theDouble; params_.finalize(); return 0;
    case 10: params_.OCR     = info.theDouble; return 0;
    case 11: params_.pop     = info.theDouble; return 0;
    case 12: rho_            = info.theDouble; return 0;
    default: return -1;
  }
}

void LadrunoHS::Print(OPS_Stream& s, int) {
  s << "LadrunoHS Material (tag: " << this->getTag() << ", type: " << this->getType() << ")" << endln;
  s << "  E50_ref=" << params_.E50_ref << " Eoed_ref=" << params_.Eoed_ref
    << " Eur_ref=" << params_.Eur_ref << " m=" << params_.m << endln;
  s << "  c=" << params_.c << " phi=" << params_.phi_deg << " psi=" << params_.psi_deg
    << " Rf=" << params_.Rf << " K0nc=" << params_.K0nc << endln;
  s << "  lodeShape=" << (params_.lode_shape == LODE_MN ? "mn" : "mc")
    << " stiffStress=" << (params_.stiff_stress == STIFF_MEAN_P ? "p" : "sigma3")
    << " alpha_cap=" << params_.alpha_cap << " beta_cap=" << params_.beta_cap << endln;
  if (params_.use_small_strain) {
    s << "  HSsmall: G0_ref=" << params_.G0_ref << " gamma_07=" << params_.gamma_07 << endln;
  }
}

// ===========================================================================
//  Recordable Responses
// ===========================================================================
Response* LadrunoHS::setResponse(const char** argv, int argc, OPS_Stream& s) {
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (std::strcmp(a, "stress") == 0 || std::strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (std::strcmp(a, "strain") == 0 || std::strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (std::strcmp(a, "tangent") == 0 || std::strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (std::strcmp(a, "state") == 0 || std::strcmp(a, "internal") == 0)
    return new MaterialResponse(this, 4, Vector(7));
  if (std::strcmp(a, "plasticShearStrain") == 0 ||
      std::strcmp(a, "gamma_p") == 0 ||
      std::strcmp(a, "gammap") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (std::strcmp(a, "preconsolidationPressure") == 0 ||
      std::strcmp(a, "p_p") == 0 ||
      std::strcmp(a, "pp") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (std::strcmp(a, "voidRatio") == 0 || std::strcmp(a, "e") == 0)
    return new MaterialResponse(this, 7, Vector(1));
  if (std::strcmp(a, "regime") == 0 || std::strcmp(a, "activeSurface") == 0)
    return new MaterialResponse(this, 8, Vector(1));
  if (std::strcmp(a, "capParams") == 0 || std::strcmp(a, "alphaCap") == 0)
    return new MaterialResponse(this, 9, Vector(2));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoHS::getResponse(int responseID, Information& matInfo) {
  switch (responseID) {
    case 1:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStress();
      return 0;
    case 2:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();
      return 0;
    case 3:
      if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent();
      return 0;
    case 4:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        double p_mean = 0.0, q = 0.0, x_sin3 = 1.0, Theta = 1.0;
        double dq[6], dTh[6];
        invariants_and_lode_grad(params_, st_trial_.sig, p_mean, q, x_sin3, Theta, dq, dTh);
        v(0) = p_mean;
        v(1) = q;
        v(2) = Theta * q;
        v(3) = st_trial_.gamma_p;
        v(4) = st_trial_.p_p;
        v(5) = st_trial_.void_ratio;
        v(6) = static_cast<double>(regime_trial_);
      }
      return 0;
    case 5:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = st_trial_.gamma_p;
      return 0;
    case 6:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = st_trial_.p_p;
      return 0;
    case 7:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = st_trial_.void_ratio;
      return 0;
    case 8:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = static_cast<double>(regime_trial_);
      return 0;
    case 9:
      if (matInfo.theVector) {
        (*(matInfo.theVector))(0) = params_.alpha_cap;
        (*(matInfo.theVector))(1) = params_.beta_cap;
      }
      return 0;
    default:
      return -1;
  }
}
