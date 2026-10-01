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

// Implementation of the LadrunoNorSand nDMaterial shell. See LadrunoNorSand.h for the design, the
// Voigt <-> tensor mapping and the refusal contract; the mathematics is LadrunoNorSandKernel.h.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "LadrunoNorSand.h"
#include "LadrunoNorSand3D.h"
#include "LadrunoNorSandPlaneStrain.h"
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <LadrunoMaterialStatus.h>   // LADRUNO_MATERIAL_REFUSED, WP-99 commit-refusal seam
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits>
#include <cmath>

using ladruno_norsand::Params;
using ladruno_norsand::State;
using ladruno_norsand::StepInfo;

// ===========================================================================
//  The 20 double-valued parameters, in ONE table shared by the parser and the wire format.
//  The position in this table is the wire offset (nsw::PARAMS + i); do not reorder.
// ===========================================================================
namespace {

struct ParamEntry { const char* flag; double Params::* mem; };

const int NPD = 20;
const ParamEntry kParam[NPD] = {
  {"-p0",           &Params::p0},          //  0  reference pressure (< 0)
  {"-kappa_hat",    &Params::kappa_hat},   //  1  elastic compressibility
  {"-eps_v0",       &Params::eps_v0},      //  2  reference elastic volumetric strain (at p0)
  {"-mu0",          &Params::mu0},         //  3  shear modulus
  {"-alpha0",       &Params::alpha0},      //  4  pressure/shear coupling
  {"-M",            &Params::M},           //  5  critical stress ratio, compression
  {"-N",            &Params::N},           //  6  F curvature
  {"-N_bar",        &Params::N_bar},       //  7  Q curvature
  {"-rho",          &Params::rho},         //  8  F ellipticity (NOT the mass density: that is -density)
  {"-rho_bar",      &Params::rho_bar},     //  9  Q ellipticity
  {"-chi",          &Params::chi},         // 10  maximum-dilatancy coefficient (< 0)
  {"-h",            &Params::h},           // 11  hardening constant
  {"-lambda_tilde", &Params::lambda_tilde},// 12  paper CSL slope
  {"-v_c0",         &Params::v_c0},        // 13  paper CSL intercept
  {"-e0",           &Params::e0},          // 14  fork CSL
  {"-lambda_c",     &Params::lambda_c},    // 15  fork CSL
  {"-xi",           &Params::xi},          // 16  fork CSL
  {"-p_a",          &Params::p_a},         // 17  fork CSL reference pressure
  {"-c1",           &Params::c1},          // 18  cap blend start
  {"-c2",           &Params::c2},          // 19  cap blend end
};
enum { I_P0 = 0, I_KAPPA, I_EPSV0, I_MU0, I_ALPHA0, I_M, I_N, I_NBAR, I_RHO, I_RHOBAR, I_CHI, I_H,
       I_LTILDE, I_VC0, I_E0, I_LC, I_XI, I_PA, I_C1, I_C2 };

// sendSelf / recvSelf wire layout (offsets into the one Vector of LadrunoNorSand::WIRE_LEN doubles)
namespace nsw {
enum {
  TAG = 0,          // material tag
  DIM = 1,          // dimension mode (informational: the class fixes it)
  DENSITY = 2,      // mass density
  PARAMS = 3,       // 20 doubles, kParam order            [3 .. 22]
  CSL = 23,         // csl_mode 0 paper | 1 fork
  ZETA = 24,        // zeta 0 WW | 1 GA
  CAP = 25,         // cap 0 none | 1 planar | 2 smooth
  SIGMA0 = 26,      // initial stress, kernel order        [26 .. 31]
  V0 = 32,          // initial specific volume
  PI0 = 33,         // initial image pressure (meaningful if PI0GIVEN)
  PI0GIVEN = 34,    // 1 = -pi0 supplied, 0 = on the yield surface
  SC = 35,          // committed State (12)                 [35 .. 46]
  ST = 47,          // trial State (12)                     [47 .. 58]
  EPSC = 59,        // committed total strain (tensor)     [59 .. 64]
  EPST = 65,        // trial total strain (tensor)         [65 .. 70]
  LATCHED = 71,
  TRIALREF = 72,    // the trial was refused
  LASTREF = 73,     // last refusal code (-1: non-finite)
  NREF = 74,        // refused trials since revertToStart
  NSUB = 75,        // accepted steps that needed substepping
  INFO = 76,        // last StepInfo (7)                    [76 .. 82]
  END = 83
};
// State packing (12): eps_e[6], pi_i, v, v0, eps_p_v, eps_p_s, D_last
const int STATE_LEN = 12;
}  // namespace nsw

static_assert(nsw::END == LadrunoNorSand::WIRE_LEN, "LadrunoNorSand wire layout out of sync with WIRE_LEN");
static_assert(nsw::PARAMS + NPD == nsw::CSL, "parameter block size");
static_assert(nsw::SC + nsw::STATE_LEN == nsw::ST && nsw::ST + nsw::STATE_LEN == nsw::EPSC, "state block size");

void packState(double* a, const State& s)
{
  for (int i = 0; i < 6; i++) a[i] = s.eps_e[i];
  a[6] = s.pi_i; a[7] = s.v; a[8] = s.v0; a[9] = s.eps_p_v; a[10] = s.eps_p_s; a[11] = s.D_last;
}

void unpackState(const double* a, State& s)
{
  for (int i = 0; i < 6; i++) s.eps_e[i] = a[i];
  s.pi_i = a[6]; s.v = a[7]; s.v0 = a[8]; s.eps_p_v = a[9]; s.eps_p_s = a[10]; s.D_last = a[11];
}

// ---------------------------------------------------------------------------
// psi = v - v_c(pressure), sheet 144a section 6 (S.22), both CSL modes. NaN if pressure >= 0.
// (pressure = pi_i gives psi_i; pressure = p gives the state parameter psi.)
// ---------------------------------------------------------------------------
double cslPsi(const Params& p, double v, double pressure)
{
  if (!(pressure < 0.0)) return std::numeric_limits<double>::quiet_NaN();
  if (p.csl_mode == 0)
    return v - p.v_c0 + p.lambda_tilde * std::log(-pressure);
  return (v - 1.0) - p.e0 + p.lambda_c * std::pow(-pressure / p.p_a, p.xi);
}

const char* refusalName(int code)
{
  switch (code) {
    case ladruno_norsand::LOCAL_NOCONV:       return "local Newton did not converge";
    case ladruno_norsand::LOCAL_LINESEARCH:   return "local line search failed";
    case ladruno_norsand::PI_NOBRACKET:       return "pi_i root not bracketed";
    case ladruno_norsand::PI_NOCONV:          return "pi_i solve did not converge";
    case ladruno_norsand::B_NONPOS:           return "pi_i* guard B <= 0";
    case ladruno_norsand::P_OR_PI_NONNEG:     return "p or pi_i not negative";
    case ladruno_norsand::NEGATIVE_DLAMBDA:   return "negative plastic multiplier";
    case ladruno_norsand::SUBSTEPS_EXHAUSTED: return "substeps exhausted (2^8)";
    case -1:                                  return "non-finite trial strain or kernel output";
    default:                                  return "unknown";
  }
}

bool allFinite(const double* a, int n)
{
  for (int i = 0; i < n; i++) if (!std::isfinite(a[i])) return false;
  return true;
}

}  // namespace

// ===========================================================================
//  constructors / destructor
// ===========================================================================
static void zeroParams(Params& p)
{
  p = Params();      // value-initialise (zero unless the kernel declares defaults)
  p.csl_mode = 0; p.zeta = 0; p.cap = 0;
}

// null constructor (broker / recvSelf)
LadrunoNorSand::LadrunoNorSand()
  : NDMaterial(0, ND_TAG_LadrunoNorSand),
    density(0.0), v0init(0.0), pi0init(std::numeric_limits<double>::quiet_NaN()), pi0given(false),
    initOk(false), dim(DIM_3D), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), nRefusals(0), nSubstepped(0)
{
  zeroParams(kp);
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  for (int i = 0; i < 6; i++) { sigma0[i] = 0.0; epsC[i] = 0.0; epsT[i] = 0.0; sigT[i] = 0.0; }
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = 0.0;
  this->setupDim();
}

LadrunoNorSand::LadrunoNorSand(int tag, const Params& p, const double sig0[6],
                               double v0, double pi0, double dens)
  : NDMaterial(tag, ND_TAG_LadrunoNorSand),
    kp(p), density(dens), v0init(v0), pi0init(pi0), pi0given(std::isfinite(pi0)),
    initOk(false), dim(DIM_3D), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), nRefusals(0), nSubstepped(0)
{
  for (int i = 0; i < 6; i++) sigma0[i] = sig0[i];
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  this->setupDim();
  this->buildInitialState();
}

// derived-class constructors
LadrunoNorSand::LadrunoNorSand(int classTag, int dimMode)
  : NDMaterial(0, classTag),
    density(0.0), v0init(0.0), pi0init(std::numeric_limits<double>::quiet_NaN()), pi0given(false),
    initOk(false), dim(dimMode), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), nRefusals(0), nSubstepped(0)
{
  zeroParams(kp);
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  for (int i = 0; i < 6; i++) { sigma0[i] = 0.0; epsC[i] = 0.0; epsT[i] = 0.0; sigT[i] = 0.0; }
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = 0.0;
  this->setupDim();
}

LadrunoNorSand::LadrunoNorSand(int tag, int classTag, const Params& p, const double sig0[6],
                               double v0, double pi0, double dens, int dimMode)
  : NDMaterial(tag, classTag),
    kp(p), density(dens), v0init(v0), pi0init(pi0), pi0given(std::isfinite(pi0)),
    initOk(false), dim(dimMode), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), nRefusals(0), nSubstepped(0)
{
  for (int i = 0; i < 6; i++) sigma0[i] = sig0[i];
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  this->setupDim();
  this->buildInitialState();
}

LadrunoNorSand::~LadrunoNorSand() {}

// ---------------------------------------------------------------------------
//  dimensional view: reduced element vector order -> kernel (== full Voigt) index
//    kernel / full order: 0:00 1:11 2:22 3:01 4:12 5:02
// ---------------------------------------------------------------------------
void LadrunoNorSand::setupDim(void)
{
  if (dim == DIM_PSTRAIN) {
    ncomp = 3; vmap[0] = 0; vmap[1] = 1; vmap[2] = 3;
    vmap[3] = vmap[4] = vmap[5] = -1;
  } else {
    ncomp = 6;
    for (int a = 0; a < 6; a++) vmap[a] = a;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
  stressOut.Zero(); strainOut.Zero(); tangentOut.Zero();
}

// Build s0 from (kp, sigma0, v0init, pi0init) and put committed = trial = s0, strain = 0.
void LadrunoNorSand::buildInitialState(void)
{
  std::string msg;
  double pi0 = pi0given ? pi0init : std::numeric_limits<double>::quiet_NaN();
  int rc = ladruno_norsand::initialState(kp, sigma0, v0init, pi0, s0, msg);
  initOk = (rc == 0);
  initMsg = msg;
  if (!initOk) s0 = State();
  sC = s0; sT = s0;
  for (int i = 0; i < 6; i++) { epsC[i] = 0.0; epsT[i] = 0.0; }
  if (initOk) {
    ladruno_norsand::stress(kp, sC, sigT);
    ladruno_norsand::elasticTangent(kp, sC, CT);
  } else {
    for (int i = 0; i < 6; i++) sigT[i] = 0.0;
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = 0.0;
  }
  lastInfo = StepInfo();
}

// ===========================================================================
//  strain interface
// ===========================================================================
int LadrunoNorSand::setTrialStrain(const Vector& e)
{
  // OpenSees reduced vector (engineering shear) -> kernel tensor components.
  for (int i = 0; i < 6; i++) epsT[i] = 0.0;
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double val = e(a);
    if (full >= 3) val *= 0.5;     // engineering gamma -> tensor eps
    epsT[full] = val;
  }
  return this->integrate();
}

int LadrunoNorSand::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoNorSand::setTrialStrainIncr(const Vector& v)
{
  // current trial strain (reduced, engineering) + increment
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double cur = (full >= 3) ? 2.0 * epsT[full] : epsT[full];
    ne(a) = cur + v(a);
  }
  return this->setTrialStrain(ne);
}

int LadrunoNorSand::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ---------------------------------------------------------------------------
//  one kernel step from the COMMITTED state to the trial strain.
//  Returns 0, or LADRUNO_MATERIAL_REFUSED (committed state untouched, trial frozen at n).
// ---------------------------------------------------------------------------
int LadrunoNorSand::integrate(void)
{
  int refusal = 0;                 // 0 = accepted

  if (latched) {
    refusal = (lastRefusal != 0) ? lastRefusal : -1;       // stays refused until revertToStart()
  } else if (!allFinite(epsT, 6)) {
    refusal = -1;                  // a diverged Newton iterate: never hand NaN to the kernel
  } else {
    double deps[6];
    for (int i = 0; i < 6; i++) deps[i] = epsT[i] - epsC[i];

    State np1 = sC;
    double sig[6];
    double C[6][6];
    StepInfo info = StepInfo();
    int rc = ladruno_norsand::step(kp, sC, deps, np1, sig, C, info);
    lastInfo = info;

    if (rc != 0 || info.refusal != 0) {
      refusal = (rc != 0) ? rc : info.refusal;
    } else {
      // finite-ness of everything we are about to adopt (NaN-blind checks are a known trap)
      bool fin = allFinite(sig, 6) && allFinite(&C[0][0], 36) &&
                 allFinite(np1.eps_e, 6) && std::isfinite(np1.pi_i) && std::isfinite(np1.v) &&
                 std::isfinite(np1.v0) && std::isfinite(np1.eps_p_v) && std::isfinite(np1.eps_p_s) &&
                 std::isfinite(np1.D_last);
      if (!fin) refusal = -1;
      else {
        sT = np1;
        for (int i = 0; i < 6; i++) sigT[i] = sig[i];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = C[i][j];
        trialRefused = false;
        if (info.substeps > 1) nSubstepped++;
        return 0;
      }
    }
  }

  // ---- refused: committed state untouched, trial frozen at the committed state ----
  trialRefused = true;
  lastRefusal = refusal;
  nRefusals++;
  sT = sC;
  ladruno_norsand::stress(kp, sC, sigT);
  ladruno_norsand::elasticTangent(kp, sC, CT);
  if (!warnedTrial) {
    warnedTrial = true;
    opserr << "WARNING LadrunoNorSand tag " << this->getTag() << ": the return map REFUSED this trial strain ("
           << refusalName(refusal) << "; local iters " << lastInfo.local_iters << ", pi iters "
           << lastInfo.pi_iters << ", substeps " << lastInfo.substeps << "). The committed state is untouched"
           << " and the trial stress/tangent are the committed ones; the step must be cut (code "
           << LADRUNO_MATERIAL_REFUSED << "). Further refusals of this point are counted silently"
           << " (`refusal` response)." << endln;
  }
  return LADRUNO_MATERIAL_REFUSED;
}

// ===========================================================================
//  element-facing accessors (Voigt, engineering shear)
// ===========================================================================
const Vector& LadrunoNorSand::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = sigT[vmap[a]];   // shear stress: same number
  return stressOut;
}

double LadrunoNorSand::getStressZZ(void)
{
  if (dim == DIM_PSTRAIN) return sigT[2];
  return NDMaterial::getStressZZ();   // NaN = not applicable
}

const Vector& LadrunoNorSand::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    strainOut(a) = (full >= 3) ? 2.0 * epsT[full] : epsT[full];   // tensor -> engineering
  }
  return strainOut;
}

// T[a][b] = d sigma_voigt_a / d eps_voigt_b = C[a][b] * w_b ; rows unchanged, shear columns halved.
const Matrix& LadrunoNorSand::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) {
      int fb = vmap[b];
      tangentOut(a, b) = CT[vmap[a]][fb] * ((fb >= 3) ? 0.5 : 1.0);
    }
  return tangentOut;
}

const Matrix& LadrunoNorSand::getInitialTangent(void)
{
  // the hyperelastic tangent at the INITIAL state (not the base default getTangent(): -initial must
  // stay a genuine initial-stiffness iteration)
  double C[6][6];
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) C[i][j] = 0.0;
  if (initOk) ladruno_norsand::elasticTangent(kp, s0, C);
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) {
      int fb = vmap[b];
      tangentOut(a, b) = C[vmap[a]][fb] * ((fb >= 3) ? 0.5 : 1.0);
    }
  return tangentOut;
}

const char* LadrunoNorSand::getType(void) const
{
  return (dim == DIM_PSTRAIN) ? "PlaneStrain" : "ThreeDimensional";
}

int LadrunoNorSand::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
void LadrunoNorSand::restoreTrialFromCommitted(void)
{
  sT = sC;
  for (int i = 0; i < 6; i++) epsT[i] = epsC[i];
  if (initOk) {
    ladruno_norsand::stress(kp, sC, sigT);
    ladruno_norsand::elasticTangent(kp, sC, CT);
  }
  trialRefused = false;
}

int LadrunoNorSand::commitState(void)
{
  // WP-99 (F7) commit-refusal seam. A refused trial that reaches commitState was DISCARDED by the
  // host element (Domain::commit() drops commitState's return code), so the refusal is declared
  // out of band and this point latches until revertToStart(). It must not commit a frozen state.
  if (latched || trialRefused) {
    if (!latched) {
      latched = true;
      if (!warnedLatch) {
        warnedLatch = true;
        opserr << "WARNING LadrunoNorSand tag " << this->getTag()
               << ": a REFUSED update (" << refusalName(lastRefusal) << ") reached commitState -- the host"
               << " element discarded the material's return code. The commit is ABORTED ("
               << LADRUNO_MATERIAL_REFUSED << ") and this point LATCHES: every further trial and commit is"
               << " refused until revertToStart(). Use an element that forwards the code (so the step is"
               << " cut instead of committed)." << endln;
      }
    }
    ladrunoNoteCommitRefusal();                 // Ladruno WP-99 (F7)
    this->restoreTrialFromCommitted();
    return LADRUNO_MATERIAL_REFUSED;
  }

  sC = sT;
  for (int i = 0; i < 6; i++) epsC[i] = epsT[i];
  return 0;
}

int LadrunoNorSand::revertToLastCommit(void)
{
  // the latch survives: only revertToStart() clears it
  this->restoreTrialFromCommitted();
  return 0;
}

int LadrunoNorSand::revertToStart(void)
{
  // back to the INITIAL state (sigma0, v0, pi0); v0 is part of the state and is restored with it
  if (initOk) {
    sC = s0; sT = s0;
    ladruno_norsand::stress(kp, sC, sigT);
    ladruno_norsand::elasticTangent(kp, sC, CT);
  } else {
    for (int i = 0; i < 6; i++) sigT[i] = 0.0;
  }
  for (int i = 0; i < 6; i++) { epsC[i] = 0.0; epsT[i] = 0.0; }
  trialRefused = false;
  latched = false;
  warnedTrial = false;
  warnedLatch = false;
  lastRefusal = 0;
  nRefusals = 0;
  nSubstepped = 0;
  lastInfo = StepInfo();
  return 0;
}

// ===========================================================================
//  copies (history isolation: every instance owns all of its state; nothing is static or shared)
// ===========================================================================
void LadrunoNorSand::copyFrom(const LadrunoNorSand& o)
{
  this->setTag(o.getTag());
  kp = o.kp;
  density = o.density;
  for (int i = 0; i < 6; i++) {
    sigma0[i] = o.sigma0[i]; epsC[i] = o.epsC[i]; epsT[i] = o.epsT[i]; sigT[i] = o.sigT[i];
  }
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = o.CT[i][j];
  v0init = o.v0init; pi0init = o.pi0init; pi0given = o.pi0given;
  initOk = o.initOk; initMsg = o.initMsg;
  s0 = o.s0; sC = o.sC; sT = o.sT;
  trialRefused = o.trialRefused;
  latched = o.latched;               // a latched point must hand out latched copies
  warnedTrial = o.warnedTrial; warnedLatch = o.warnedLatch;
  lastRefusal = o.lastRefusal; nRefusals = o.nRefusals; nSubstepped = o.nSubstepped;
  lastInfo = o.lastInfo;
  // dim / ncomp / vmap / output buffers belong to the clone's own class: not copied
}

NDMaterial* LadrunoNorSand::getCopy(void)
{
  LadrunoNorSand* c = new LadrunoNorSand();
  c->copyFrom(*this);
  return c;
}

NDMaterial* LadrunoNorSand::getCopy(const char* type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
    LadrunoNorSand3D* c = new LadrunoNorSand3D();
    c->copyFrom(*this);
    return c;
  }
  if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0) {
    LadrunoNorSandPlaneStrain* c = new LadrunoNorSandPlaneStrain();
    c->copyFrom(*this);
    return c;
  }
  return NDMaterial::getCopy(type);   // let the base report the unsupported type
}

// ===========================================================================
//  parallel / database
// ===========================================================================
int LadrunoNorSand::sendSelf(int commitTag, Channel& theChannel)
{
  Vector data(WIRE_LEN);
  data.Zero();
  data(nsw::TAG) = this->getTag();
  data(nsw::DIM) = dim;
  data(nsw::DENSITY) = density;
  for (int i = 0; i < NPD; i++) data(nsw::PARAMS + i) = kp.*(kParam[i].mem);
  data(nsw::CSL) = kp.csl_mode;
  data(nsw::ZETA) = kp.zeta;
  data(nsw::CAP) = kp.cap;
  for (int i = 0; i < 6; i++) data(nsw::SIGMA0 + i) = sigma0[i];
  data(nsw::V0) = v0init;
  data(nsw::PI0) = pi0given ? pi0init : 0.0;
  data(nsw::PI0GIVEN) = pi0given ? 1.0 : 0.0;
  {
    double a[nsw::STATE_LEN];
    packState(a, sC); for (int i = 0; i < nsw::STATE_LEN; i++) data(nsw::SC + i) = a[i];
    packState(a, sT); for (int i = 0; i < nsw::STATE_LEN; i++) data(nsw::ST + i) = a[i];
  }
  for (int i = 0; i < 6; i++) { data(nsw::EPSC + i) = epsC[i]; data(nsw::EPST + i) = epsT[i]; }
  data(nsw::LATCHED) = latched ? 1.0 : 0.0;
  data(nsw::TRIALREF) = trialRefused ? 1.0 : 0.0;
  data(nsw::LASTREF) = lastRefusal;
  data(nsw::NREF) = nRefusals;
  data(nsw::NSUB) = nSubstepped;
  data(nsw::INFO + 0) = lastInfo.refusal;
  data(nsw::INFO + 1) = lastInfo.plastic;
  data(nsw::INFO + 2) = lastInfo.vertex;
  data(nsw::INFO + 3) = lastInfo.cap_active;
  data(nsw::INFO + 4) = lastInfo.local_iters;
  data(nsw::INFO + 5) = lastInfo.pi_iters;
  data(nsw::INFO + 6) = lastInfo.substeps;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoNorSand::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoNorSand::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  Vector data(WIRE_LEN);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoNorSand::recvSelf - failed to recv vector\n";
    return -1;
  }
  this->setTag((int)data(nsw::TAG));
  density = data(nsw::DENSITY);
  for (int i = 0; i < NPD; i++) kp.*(kParam[i].mem) = data(nsw::PARAMS + i);
  kp.csl_mode = (int)data(nsw::CSL);
  kp.zeta = (int)data(nsw::ZETA);
  kp.cap = (int)data(nsw::CAP);
  for (int i = 0; i < 6; i++) sigma0[i] = data(nsw::SIGMA0 + i);
  v0init = data(nsw::V0);
  pi0given = (data(nsw::PI0GIVEN) != 0.0);
  pi0init = pi0given ? data(nsw::PI0) : std::numeric_limits<double>::quiet_NaN();

  // rebuild s0 (the initial state is a function of the inputs), then overwrite committed/trial
  this->setupDim();
  this->buildInitialState();

  double a[nsw::STATE_LEN];
  for (int i = 0; i < nsw::STATE_LEN; i++) a[i] = data(nsw::SC + i);
  unpackState(a, sC);
  for (int i = 0; i < nsw::STATE_LEN; i++) a[i] = data(nsw::ST + i);
  unpackState(a, sT);
  for (int i = 0; i < 6; i++) { epsC[i] = data(nsw::EPSC + i); epsT[i] = data(nsw::EPST + i); }
  latched = (data(nsw::LATCHED) != 0.0);
  trialRefused = (data(nsw::TRIALREF) != 0.0);
  lastRefusal = (int)data(nsw::LASTREF);
  nRefusals = (int)data(nsw::NREF);
  nSubstepped = (int)data(nsw::NSUB);
  lastInfo.refusal = (int)data(nsw::INFO + 0);
  lastInfo.plastic = (int)data(nsw::INFO + 1);
  lastInfo.vertex = (int)data(nsw::INFO + 2);
  lastInfo.cap_active = (int)data(nsw::INFO + 3);
  lastInfo.local_iters = (int)data(nsw::INFO + 4);
  lastInfo.pi_iters = (int)data(nsw::INFO + 5);
  lastInfo.substeps = (int)data(nsw::INFO + 6);

  // The trial STRESS is a function of the trial state; the trial TANGENT is not sent, so until the
  // next setTrialStrain the tangent is the hyperelastic one at the trial state.
  if (initOk) {
    ladruno_norsand::stress(kp, sT, sigT);
    ladruno_norsand::elasticTangent(kp, sT, CT);
  }
  return 0;
}

// ===========================================================================
//  print / parameter echo
// ===========================================================================
static const char* cslName(int m) { return m == 0 ? "paper (v_c = v_c0 - lambda_tilde ln(-p))"
                                                   : "fork (e_c = e0 - lambda_c (-p/p_a)^xi)"; }
static const char* zetaName(int z) { return z == 0 ? "WW (Willam-Warnke)" : "GA (Gudehus-Argyris)"; }
static const char* capName(int c) { return c == 0 ? "none" : (c == 1 ? "planar" : "smooth"); }

void LadrunoNorSand::echoParameters(OPS_Stream& s) const
{
  s << "LadrunoNorSand tag " << this->getTag() << " (WP-144; NorSand, Andrade & Borja 2006 form; classTag "
    << this->getClassTag() << ")" << endln;
  s << "  elastic (BA06 energy): p0=" << kp.p0 << " kappa_hat=" << kp.kappa_hat << " eps_v0=" << kp.eps_v0
    << " mu0=" << kp.mu0 << " alpha0=" << kp.alpha0 << endln;
  s << "  surface: M=" << kp.M << " N=" << kp.N << " N_bar=" << kp.N_bar << " rho=" << kp.rho
    << " rho_bar=" << kp.rho_bar << "  (rho = ellipticity; the mass density is -density)" << endln;
  s << "  dilatancy/hardening: chi=" << kp.chi << " h=" << kp.h << endln;
  s << "  CSL: " << cslName(kp.csl_mode) << ": ";
  if (kp.csl_mode == 0) s << "lambda_tilde=" << kp.lambda_tilde << " v_c0=" << kp.v_c0 << endln;
  else                  s << "e0=" << kp.e0 << " lambda_c=" << kp.lambda_c << " xi=" << kp.xi << " p_a=" << kp.p_a << endln;
  s << "  Lode shape zeta: " << zetaName(kp.zeta) << "; Q-cap: " << capName(kp.cap);
  if (kp.cap != 0) s << " (c1=" << kp.c1 << " c2=" << kp.c2 << ")";
  s << endln;
  s << "  initial state: sigma0 = [" << sigma0[0] << " " << sigma0[1] << " " << sigma0[2] << " "
    << sigma0[3] << " " << sigma0[4] << " " << sigma0[5] << "] (11 22 33 12 23 13), v0=" << v0init
    << ", pi_i0=" << s0.pi_i << (pi0given ? " (given)" : " (on the yield surface)")
    << ", psi_i0=" << cslPsi(kp, s0.v, s0.pi_i) << endln;
  s << "  density=" << density << endln;
  s << "  NON-symmetric consistent tangent: use an unsymmetric solver. A refused return map reaches the"
    << " element as " << LADRUNO_MATERIAL_REFUSED << " (trial) or aborts the commit and latches the point"
    << " (a host that discards the code); work per step is bounded by the kernel caps." << endln;
}

void LadrunoNorSand::Print(OPS_Stream& s, int /*flag*/)
{
  s << endln;
  this->echoParameters(s);
  s << "  committed: pi_i=" << sC.pi_i << " v=" << sC.v << " v0=" << sC.v0 << " eps_p_v=" << sC.eps_p_v
    << " eps_p_s=" << sC.eps_p_s << endln;
  s << "  refused trials=" << nRefusals << " substepped steps=" << nSubstepped
    << " latched=" << (latched ? 1 : 0) << endln;
}

// ===========================================================================
//  recordable responses
//    stress, strain, tangent         element-facing (Voigt, engineering shear)
//    state          [pi_i, psi_i, v, v0, eps_p_v, eps_p_s]
//    D              dissipation of the last step (>= 0)
//    refusal        [last refusal code (-1 non-finite), refused trials, latched]
//    substeps       [substeps of the last step, steps that needed substepping]
//    stepInfo       [refusal, plastic, vertex, cap_active, local_iters, pi_iters, substeps]
//    psi            state parameter psi = v - v_c(p)
//    elasticStrain  elastic strain, Voigt engineering shear
// ===========================================================================
Response* LadrunoNorSand::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "state") == 0)
    return new MaterialResponse(this, 4, Vector(6));
  if (strcmp(a, "D") == 0 || strcmp(a, "dissipation") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "refusal") == 0)
    return new MaterialResponse(this, 6, Vector(3));
  if (strcmp(a, "substeps") == 0)
    return new MaterialResponse(this, 7, Vector(2));
  if (strcmp(a, "stepInfo") == 0)
    return new MaterialResponse(this, 8, Vector(7));
  if (strcmp(a, "psi") == 0)
    return new MaterialResponse(this, 9, Vector(1));
  if (strcmp(a, "elasticStrain") == 0 || strcmp(a, "elasticStrains") == 0)
    return new MaterialResponse(this, 10, Vector(ncomp));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoNorSand::getResponse(int responseID, Information& matInfo)
{
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
        v(0) = sT.pi_i;
        v(1) = cslPsi(kp, sT.v, sT.pi_i);
        v(2) = sT.v;
        v(3) = sT.v0;
        v(4) = sT.eps_p_v;
        v(5) = sT.eps_p_s;
      }
      return 0;
    case 5:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = sT.D_last;
      return 0;
    case 6:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = lastRefusal; v(1) = nRefusals; v(2) = latched ? 1.0 : 0.0;
      }
      return 0;
    case 7:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = lastInfo.substeps; v(1) = nSubstepped;
      }
      return 0;
    case 8:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = lastInfo.refusal; v(1) = lastInfo.plastic; v(2) = lastInfo.vertex; v(3) = lastInfo.cap_active;
        v(4) = lastInfo.local_iters; v(5) = lastInfo.pi_iters; v(6) = lastInfo.substeps;
      }
      return 0;
    case 9:
      if (matInfo.theVector) {
        double pm = (sigT[0] + sigT[1] + sigT[2]) / 3.0;
        (*(matInfo.theVector))(0) = cslPsi(kp, sT.v, pm);
      }
      return 0;
    case 10:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) {
          int full = vmap[a];
          v(a) = (full >= 3) ? 2.0 * sT.eps_e[full] : sT.eps_e[full];
        }
      }
      return 0;
    default:
      return -1;
  }
}

// ===========================================================================
//  OPS parser (Tcl and Python)
//
//   nDMaterial LadrunoNorSand tag
//       -p0 p0 -kappa_hat kh -mu0 mu0 [-eps_v0 e] [-alpha0 a]
//       -M M -N N [-N_bar Nb] -rho rho [-rho_bar rb] -chi chi -h h
//       [-csl paper|fork]   paper: -lambda_tilde lt -v_c0 vc0     fork: -e0 e0 -lambda_c lc -xi xi [-p_a pa]
//       [-zeta WW|GA] [-cap none|planar|smooth [-c1 c1] [-c2 c2]]
//       -v0 v0 [-sigma0 s11 s22 s33 s12 s23 s13] [-pi0 pi_i0] [-density d]
//
//   defaults: N_bar = N, rho_bar = rho, eps_v0 = alpha0 = 0, -csl paper, -zeta WW, -cap none,
//   -sigma0 = isotropic p0, -pi0 = on the yield surface, -density 0.
//   -rho is the ELLIPTICITY of F (the model parameter); the mass density is -density.
//   The kernel's validate() is called: a refused parameter set is a hard error, rho > rho_bar a warning.
// ===========================================================================
static void nsUsage()
{
  opserr << "Want: nDMaterial LadrunoNorSand tag? -p0 p0? -kappa_hat kh? -mu0 mu0? <-eps_v0 e?> <-alpha0 a?>"
         << " -M M? -N N? <-N_bar Nb?> -rho rho? <-rho_bar rb?> -chi chi? -h h?"
         << " <-csl paper|fork> (paper: -lambda_tilde lt? -v_c0 vc0?  fork: -e0 e0? -lambda_c lc? -xi xi? <-p_a pa?>)"
         << " <-zeta WW|GA> <-cap none|planar|smooth> <-c1 c1?> <-c2 c2?>"
         << " -v0 v0? <-sigma0 s11? s22? s33? s12? s23? s13?> <-pi0 pi_i0?> <-density d?>" << endln;
}

static bool ieq(const char* a, const char* b)
{
  for (; *a && *b; a++, b++)
    if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return false;
  return *a == 0 && *b == 0;
}

// Parse everything after the command name; returns the prototype (a 3D-capable LadrunoNorSand).
static void* parseLadrunoNorSand(void)
{
  if (OPS_GetNumRemainingInputArgs() < 1) { nsUsage(); return 0; }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING LadrunoNorSand: invalid tag\n";
    nsUsage();
    return 0;
  }

  Params p;
  zeroParams(p);
  // defaults for the parameters of the INACTIVE modes (never echoed); the active ones are required
  p.eps_v0 = 0.0; p.alpha0 = 0.0;
  p.lambda_tilde = 0.0135; p.v_c0 = 1.81;
  p.e0 = 0.83; p.lambda_c = 0.027; p.xi = 0.45; p.p_a = 101.325;
  p.c1 = 0.05; p.c2 = 0.15;

  bool seen[NPD];
  for (int i = 0; i < NPD; i++) seen[i] = false;
  double sigma0[6] = {0, 0, 0, 0, 0, 0};
  bool haveSigma0 = false, haveV0 = false, havePi0 = false;
  double v0 = 0.0, pi0 = std::numeric_limits<double>::quiet_NaN(), density = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (flag == 0) break;

    int idx = -1;
    for (int i = 0; i < NPD; i++)
      if (strcmp(flag, kParam[i].flag) == 0) { idx = i; break; }

    if (idx >= 0) {
      double d;
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &d) < 0) {
        opserr << "WARNING LadrunoNorSand: " << flag << " wants a number\n";
        return 0;
      }
      p.*(kParam[idx].mem) = d;
      seen[idx] = true;
    }
    else if (strcmp(flag, "-csl") == 0 || strcmp(flag, "-zeta") == 0 || strcmp(flag, "-cap") == 0) {
      const char* val = OPS_GetString();
      if (val == 0) { opserr << "WARNING LadrunoNorSand: " << flag << " wants a mode\n"; return 0; }
      if (strcmp(flag, "-csl") == 0) {
        if (ieq(val, "paper")) p.csl_mode = 0;
        else if (ieq(val, "fork")) p.csl_mode = 1;
        else { opserr << "WARNING LadrunoNorSand: -csl wants paper|fork, got '" << val << "'\n"; return 0; }
      } else if (strcmp(flag, "-zeta") == 0) {
        if (ieq(val, "WW")) p.zeta = 0;
        else if (ieq(val, "GA")) p.zeta = 1;
        else { opserr << "WARNING LadrunoNorSand: -zeta wants WW|GA, got '" << val << "'\n"; return 0; }
      } else {
        if (ieq(val, "none")) p.cap = 0;
        else if (ieq(val, "planar")) p.cap = 1;
        else if (ieq(val, "smooth")) p.cap = 2;
        else { opserr << "WARNING LadrunoNorSand: -cap wants none|planar|smooth, got '" << val << "'\n"; return 0; }
      }
    }
    else if (strcmp(flag, "-sigma0") == 0) {
      numData = 6;
      if (OPS_GetDoubleInput(&numData, sigma0) < 0) {
        opserr << "WARNING LadrunoNorSand: -sigma0 wants s11 s22 s33 s12 s23 s13\n";
        return 0;
      }
      haveSigma0 = true;
    }
    else if (strcmp(flag, "-v0") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &v0) < 0) { opserr << "WARNING LadrunoNorSand: -v0 wants a number\n"; return 0; }
      haveV0 = true;
    }
    else if (strcmp(flag, "-pi0") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &pi0) < 0) { opserr << "WARNING LadrunoNorSand: -pi0 wants a number\n"; return 0; }
      havePi0 = true;
    }
    else if (strcmp(flag, "-density") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &density) < 0) { opserr << "WARNING LadrunoNorSand: -density wants a number\n"; return 0; }
    }
    else {
      opserr << "WARNING LadrunoNorSand: unknown flag '" << flag << "'\n";
      nsUsage();
      return 0;
    }
  }

  // ---- defaults that depend on other parameters ----
  if (!seen[I_NBAR])   p.N_bar = p.N;
  if (!seen[I_RHOBAR]) p.rho_bar = p.rho;
  if (p.cap == 1 && !seen[I_C2] && seen[I_C1]) p.c2 = p.c1;      // planar: c1 = c2

  // ---- required parameters ----
  {
    static const int reqAll[] = {I_P0, I_KAPPA, I_MU0, I_M, I_N, I_RHO, I_CHI, I_H};
    static const int reqPaper[] = {I_LTILDE, I_VC0};
    static const int reqFork[] = {I_E0, I_LC, I_XI};
    bool missing = false;
    for (size_t i = 0; i < sizeof(reqAll) / sizeof(int); i++)
      if (!seen[reqAll[i]]) { opserr << "WARNING LadrunoNorSand: missing required " << kParam[reqAll[i]].flag << "\n"; missing = true; }
    if (p.csl_mode == 0) {
      for (size_t i = 0; i < sizeof(reqPaper) / sizeof(int); i++)
        if (!seen[reqPaper[i]]) { opserr << "WARNING LadrunoNorSand: -csl paper needs " << kParam[reqPaper[i]].flag << "\n"; missing = true; }
    } else {
      for (size_t i = 0; i < sizeof(reqFork) / sizeof(int); i++)
        if (!seen[reqFork[i]]) { opserr << "WARNING LadrunoNorSand: -csl fork needs " << kParam[reqFork[i]].flag << "\n"; missing = true; }
    }
    if (p.cap != 0 && !seen[I_C1]) { opserr << "WARNING LadrunoNorSand: -cap planar|smooth needs -c1\n"; missing = true; }
    if (p.cap == 2 && !seen[I_C2]) { opserr << "WARNING LadrunoNorSand: -cap smooth needs -c2\n"; missing = true; }
    if (!haveV0) { opserr << "WARNING LadrunoNorSand: missing required -v0 (initial specific volume)\n"; missing = true; }
    if (missing) { nsUsage(); return 0; }
  }

  // ---- the owner-approved parameter refusals ----
  {
    std::string msg;
    bool warnRho = false;
    int rc = ladruno_norsand::validate(p, msg, warnRho);
    if (rc != 0) {
      opserr << "WARNING LadrunoNorSand tag " << tag << ": parameter set REFUSED (code " << rc << "): "
             << msg.c_str() << endln;
      return 0;
    }
    // on success msg carries the kernel's warnings (rho > rho_bar, chi > 0), if any
    if (warnRho || !msg.empty())
      opserr << "WARNING LadrunoNorSand tag " << tag << ": " << msg.c_str()
             << (warnRho ? " [rho > rho_bar: only a warning; the dissipation guarantee still holds]" : "") << endln;
  }

  if (!haveSigma0)
    for (int i = 0; i < 3; i++) sigma0[i] = p.p0;     // isotropic at the reference pressure

  LadrunoNorSand* mat = new LadrunoNorSand(tag, p, sigma0, v0, havePi0 ? pi0 : std::numeric_limits<double>::quiet_NaN(), density);
  if (!mat->initOK()) {
    opserr << "WARNING LadrunoNorSand tag " << tag << ": the initial state was REFUSED: " << mat->initMessage() << endln;
    delete mat;
    return 0;
  }
  mat->echoParameters(opserr);
  return mat;
}

void* OPS_LadrunoNorSand(void)
{
  return parseLadrunoNorSand();
}
