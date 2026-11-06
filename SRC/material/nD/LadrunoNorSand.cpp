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
#include <cstdio>
#include <string>

using ladruno_norsand::Params;
using ladruno_norsand::State;
using ladruno_norsand::StepInfo;

// ===========================================================================
//  The 24 double-valued parameters, in ONE table used by the parser. The first NPD_BASE = 20 are the
//  round-0 block of the wire (offset nsw::PARAMS + i; do not reorder); the last four (HAR k, g, n_e and the
//  floor p_min, WP-144 round 3) are appended on the wire at nsw::XPARAMS + (i - NPD_BASE), so the older offsets
//  did not move.
// ===========================================================================
namespace {

struct ParamEntry { const char* flag; double Params::* mem; };

const int NPD_BASE = 20;
const int NPD = 24;
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
  {"-k",            &Params::k},           // 20  HAR bulk stiffness factor (energy HAR only)
  {"-g",            &Params::g},           // 21  HAR shear stiffness factor (energy HAR only)
  {"-n",            &Params::n_e},         // 22  HAR pressure exponent (energy HAR only; default 1/2)
  {"-pmin",         &Params::p_min},       // 23  p' floor (default 5e-3 p_ref; 0 = off)
};
enum { I_P0 = 0, I_KAPPA, I_EPSV0, I_MU0, I_ALPHA0, I_M, I_N, I_NBAR, I_RHO, I_RHOBAR, I_CHI, I_H,
       I_LTILDE, I_VC0, I_E0, I_LC, I_XI, I_PA, I_C1, I_C2, I_K, I_G, I_NE, I_PMIN };

// sendSelf / recvSelf wire layout (offsets into the one Vector of LadrunoNorSand::WIRE_LEN doubles)
//
// What is NOT on the wire, and why (WP-144 G2 finding, LEDGER_quirks "Domain::recvSelf calls update()"):
// the TRIAL stress (6) and the trial consistent tangent (36) are not sent. Domain::recvSelf calls
// element->update() on every element after its recvSelf, which re-integrates each point from its
// restored COMMITTED state to the restored nodal strain and so recomputes both: no public route can
// observe the 42 doubles, and carrying them would only suggest a guarantee the host does not give.
// recvSelf rebuilds them from the restored trial State instead (the hyperelastic stress and tangent
// of sT), so a restored point is self-consistent even before that update(). Everything the update()
// does NOT recompute (committed State, trial State, strains, the latch and the refusal record)
// stays on the wire.
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
  PI0 = 33,         // initial image pressure
  SC = 34,          // committed State (12)                 [34 .. 45]
  ST = 46,          // trial State (12)                     [46 .. 57]
  EPSC = 58,        // committed total strain (tensor)     [58 .. 63]
  EPST = 64,        // trial total strain (tensor)         [64 .. 69]
  LATCHED = 70,
  TRIALREF = 71,    // the trial was refused
  LASTREF = 72,     // last refusal code (-1: non-finite)
  LASTFIN = 73,     // finest-level cause of that refusal (StepInfo.finest; 0 none)
  LASTFINSUB = 74,  // ... and its sub-reason (StepInfo.finest_sub)
  NREF = 75,        // refused trials since revertToStart
  NSUB = 76,        // COMMITTED steps that needed substepping
  INFO = 77,        // last StepInfo (9): refusal plastic vertex cap_active local_iters pi_iters substeps finest finest_sub
  // ---- round 3 (HAR energy, p' floor): appended, so every offset above is unchanged ----
  ENERGY = 86,      // 0 BA06 | 1 HAR
  XPARAMS = 87,     // k, g, n_e, p_min (kParam[NPD_BASE ..])  [87 .. 90]
  SCF = 91,         // committed floor counters (6)           [91 .. 96]
  STF = 97,         // trial floor counters (6)               [97 .. 102]
  INFOF = 103,      // last StepInfo floor block (6): floor_tr floor_post at_floor deps_f_v W_f E_f  [103 .. 108]
  EFC = 109,        // cumulative floor energy E_f over the committed history
  EFT = 110,        // ... including the trial step
  PI0AUTO = 111,    // pi_i0 came from the unified rule (S.53)
  END = 112
};
// State packing (12): eps_e[6], pi_i, v, v0, eps_p_v, eps_p_s, D_last
const int STATE_LEN = 12;
const int INFO_LEN = 9;
// Floor counters packing (6): eps_f_v, W_f, n_f_tr, n_f_post, n_f_init, at_floor
const int FLOOR_LEN = 6;
const int INFOF_LEN = 6;
}  // namespace nsw

static_assert(nsw::END == LadrunoNorSand::WIRE_LEN, "LadrunoNorSand wire layout out of sync with WIRE_LEN");
static_assert(nsw::PARAMS + NPD_BASE == nsw::CSL, "parameter block size");
static_assert(nsw::XPARAMS + (NPD - NPD_BASE) == nsw::SCF, "round-3 parameter block size");
static_assert(nsw::SCF + nsw::FLOOR_LEN == nsw::STF && nsw::STF + nsw::FLOOR_LEN == nsw::INFOF, "floor counter block size");
static_assert(nsw::INFOF + nsw::INFOF_LEN == nsw::EFC && nsw::EFC + 1 == nsw::EFT && nsw::EFT + 1 == nsw::PI0AUTO &&
              nsw::PI0AUTO + 1 == nsw::END, "round-3 tail block size");
static_assert(nsw::INFO + nsw::INFO_LEN == nsw::ENERGY, "round-3 blocks follow the StepInfo block");
static_assert(nsw::PI0 + 1 == nsw::SC, "initial-state block size");
static_assert(nsw::SC + nsw::STATE_LEN == nsw::ST && nsw::ST + nsw::STATE_LEN == nsw::EPSC, "state block size");
static_assert(nsw::EPSC + 6 == nsw::EPST && nsw::EPST + 6 == nsw::LATCHED, "strain block size");
static_assert(nsw::NSUB + 1 == nsw::INFO, "StepInfo block position");
// FE_Datastore keys a sent Vector by its SIZE (LEDGER_quirks "FE_Datastore keys a sent Vector by its SIZE"). This class
// sends exactly ONE Vector under its dbTag and commitTag and the NDMaterial base sends none, so there is nothing to
// collide with. If a base block, a subclass block or a second vector is ever added, its length must differ from
// WIRE_LEN: extend this list with that length (SANISAND's LWIRE_SIZE != 97 is the pattern).
static_assert(LadrunoNorSand::WIRE_LEN > 0 && nsw::END == 112, "LadrunoNorSand wire length changed: re-check the FE_Datastore size-uniqueness note above");

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

// the cumulative floor counters of a State (sheet 9.7), a separate wire block so the State block did not move
void packFloor(double* a, const State& s)
{
  a[0] = s.eps_f_v; a[1] = s.W_f; a[2] = s.n_f_tr; a[3] = s.n_f_post; a[4] = s.n_f_init; a[5] = s.at_floor;
}

void unpackFloor(const double* a, State& s)
{
  s.eps_f_v = a[0]; s.W_f = a[1];
  s.n_f_tr = (int)a[2]; s.n_f_post = (int)a[3]; s.n_f_init = (int)a[4]; s.at_floor = (int)a[5];
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

// The sub-reason of the finest-level refusal (kernel EvalErr name), "" if none.
const char* subName(int sub)
{
  return ladruno_norsand::detail::evalErrName(sub);
}

// "finest level: <cause> (<sub-reason>)", empty when the kernel recorded none (shell-level refusals).
std::string finestText(int fin, int finSub)
{
  if (fin == 0) return std::string();
  std::string t = std::string("finest level: ") + refusalName(fin);
  if (finSub != 0) { t += " ("; t += subName(finSub); t += ")"; }
  return t;
}

// Yield function F(sigma, pi_i) of the sheet (S.12), via the kernel's own flow evaluation. NaN if it
// cannot be evaluated (p >= 0, pi_i >= 0, eigen-decomposition failure).
double yieldF(const Params& p, const double sig[6], double pi)
{
  namespace kd = ladruno_norsand::detail;
  double S[3][3], w[3], V[3][3];
  kd::t6_to_m(sig, S);
  if (!kd::eig_sym3(S, w, V)) return std::numeric_limits<double>::quiet_NaN();
  kd::Invariants inv;
  kd::invariants(w, inv);
  kd::Flow fl;
  if (kd::flow(p, inv, pi, fl) != kd::EE_NONE) return std::numeric_limits<double>::quiet_NaN();
  return fl.F;
}

// Start-state classification, in units of p_ref (|p0| under BA06, p_a under HAR: ladruno_norsand::pRef), ONE band: |F0| <= 1e-6 |p0| is ON the surface (accepted with a
// WARNING); F0 > 1e-6 |p0| is OUTSIDE and REFUSED. The two constants are the same number by construction, so
// no start state is both "ON" and refused.
const double F0_ON_SURFACE_REL = 1.0e-6;
const double F0_OUTSIDE_REL = F0_ON_SURFACE_REL;

// shell-level initial-state refusal codes (the kernel uses 100 + the validate code, and 101..106)
enum { INIT_PI0_MISSING = 107, INIT_F0_NAN = 108, INIT_OUTSIDE = 109 };

// parser-level refusal codes of the energy option (sheet 2.4; the kernel's validate() owns 1..26, initialState 101..108)
enum { PARSE_BA06_WITH_HAR = 201,     // a BA06 constant (-p0 -kappa_hat -eps_v0 -mu0 -alpha0) given with -energy HAR
       PARSE_HAR_WITH_BA06 = 202,     // a HAR constant (-k -g -n -G0 -nu -e_ref) given with the BA06 energy
       PARSE_HAR_AMBIGUOUS = 203,     // both the direct (-k -g) and the DM04 (-G0 -nu) form given
       PARSE_HAR_INCOMPLETE = 204,    // -k without -g (or the reverse), -G0 without -nu (or the reverse), -e_ref alone
       PARSE_DM04_RANGE = 206,        // -G0 <= 0, -nu outside (-1, 1/2) or -e_ref <= -1
       PARSE_PI0_BOTH = 207 };        // -pi0 and -pi0_auto both given

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
    density(0.0), v0init(0.0), pi0init(std::numeric_limits<double>::quiet_NaN()),
    F0init(std::numeric_limits<double>::quiet_NaN()), pi0Auto(false),
    initOk(false), dim(DIM_3D), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), lastFinest(0), lastFinestSub(0), nRefusals(0), nSubstepped(0), efC(0.0), efT(0.0)
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
    kp(p), density(dens), v0init(v0), pi0init(pi0), F0init(std::numeric_limits<double>::quiet_NaN()), pi0Auto(false),
    initOk(false), dim(DIM_3D), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), lastFinest(0), lastFinestSub(0), nRefusals(0), nSubstepped(0), efC(0.0), efT(0.0)
{
  for (int i = 0; i < 6; i++) sigma0[i] = sig0[i];
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  this->setupDim();
  this->buildInitialState();
}

// derived-class constructors
LadrunoNorSand::LadrunoNorSand(int clsTag, int dimMode)
  : NDMaterial(0, clsTag),
    density(0.0), v0init(0.0), pi0init(std::numeric_limits<double>::quiet_NaN()),
    F0init(std::numeric_limits<double>::quiet_NaN()), pi0Auto(false),
    initOk(false), dim(dimMode), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), lastFinest(0), lastFinestSub(0), nRefusals(0), nSubstepped(0), efC(0.0), efT(0.0)
{
  zeroParams(kp);
  s0 = State(); sC = State(); sT = State(); lastInfo = StepInfo();
  for (int i = 0; i < 6; i++) { sigma0[i] = 0.0; epsC[i] = 0.0; epsT[i] = 0.0; sigT[i] = 0.0; }
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = 0.0;
  this->setupDim();
}

LadrunoNorSand::LadrunoNorSand(int tag, int clsTag, const Params& p, const double sig0[6],
                               double v0, double pi0, double dens, int dimMode)
  : NDMaterial(tag, clsTag),
    kp(p), density(dens), v0init(v0), pi0init(pi0), F0init(std::numeric_limits<double>::quiet_NaN()), pi0Auto(false),
    initOk(false), dim(dimMode), ncomp(6),
    trialRefused(false), latched(false), warnedTrial(false), warnedLatch(false),
    lastRefusal(0), lastFinest(0), lastFinestSub(0), nRefusals(0), nSubstepped(0), efC(0.0), efT(0.0)
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
// pi0init = NaN selects the UNIFIED pi_i0 rule (sheet 5.4 (S.53), owner decision (d)); the resolved value then
// replaces the NaN (so recvSelf, which rebuilds s0 from the stored inputs, sees a finite pi_i0). A finite pi_i0
// is the deck's and overrides the rule. +-Inf is not a pi_i0: refused like a missing one.
// The p' floor (sheet 9.7) may project the initial state (State::n_f_init = 1): the state then sits at
// sigma(eps^e_f), not at the deck's sigma0, and F0 is evaluated at THAT stress.
void LadrunoNorSand::buildInitialState(void)
{
  std::string msg;
  F0init = std::numeric_limits<double>::quiet_NaN();
  int rc;
  const bool wantAuto = std::isnan(pi0init);
  if (wantAuto) pi0Auto = true;
  if (!wantAuto && !std::isfinite(pi0init)) {
    // -pi0 is REQUIRED: the kernel's old "NaN = on the yield surface" default is not offered by the shell
    // (a silent on-surface start makes the first loading step plastic by construction). -pi0_auto asks for
    // the unified rule of sheet 5.4 explicitly.
    rc = INIT_PI0_MISSING;
    msg = "-pi0 (the initial image pressure pi_i0) is required and must be finite and < 0";
  } else {
    rc = ladruno_norsand::initialState(kp, sigma0, v0init, pi0init, s0, msg);   // NaN: PI0_UNIFIED (default rule)
  }
  if (rc == 0) {
    if (wantAuto) pi0init = s0.pi_i;
    // the stress the initial state really has: the deck's sigma0, or sigma(eps^e_f) if the floor projected it
    double sigEff[6];
    if (s0.n_f_init != 0) ladruno_norsand::stress(kp, s0, sigEff);
    else for (int i = 0; i < 6; i++) sigEff[i] = sigma0[i];
    // F(sigma, pi_i0): an initial stress OUTSIDE the surface is inadmissible
    F0init = yieldF(kp, sigEff, pi0init);
    const double tol = F0_OUTSIDE_REL * ladruno_norsand::pRef(kp);   // F0 within +-tol: ON the surface
    if (!std::isfinite(F0init)) {
      rc = INIT_F0_NAN;
      msg = "the yield function F(sigma0, pi_i0) could not be evaluated";
    } else if (F0init > tol) {
      rc = INIT_OUTSIDE;
      char buf[512];
      std::snprintf(buf, sizeof(buf),
                    "the initial stress lies OUTSIDE the yield surface of -pi0 = %.17g: F(sigma0, pi_i0) = %.6g > %.3g"
                    " (= %.1g p_ref)", pi0init, F0init, tol, F0_OUTSIDE_REL);
      msg = buf;
      // hint: the pi_i0 of the surface THROUGH sigma0 (a more negative -pi0 starts inside it): the pre-round-3
      // rule, eta* = eta_init, NOT the unified one (whose surface passes through eta* = max(eta_init, c2 M))
      State sOn; std::string m2;
      if (ladruno_norsand::initialState(kp, sigma0, v0init, std::numeric_limits<double>::quiet_NaN(), sOn, m2,
                                        ladruno_norsand::PI0_LEGACY) == 0) {
        std::snprintf(buf, sizeof(buf), "; the surface through sigma0 is at pi_i = %.17g: give -pi0 <= that"
                      " (more negative starts inside)", sOn.pi_i);
        msg += buf;
      }
    }
  }
  initOk = (rc == 0);
  initMsg = msg;
  if (!initOk) s0 = State();
  sC = s0; sT = s0;
  efC = 0.0; efT = 0.0;
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
  int fin = 0, finSub = 0;         // finest-level cause of a kernel refusal (StepInfo)
  bool newCause = true;            // false: a latched point repeats the recorded cause

  if (latched) {
    refusal = (lastRefusal != 0) ? lastRefusal : -1;       // stays refused until revertToStart()
    newCause = false;
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
      fin = info.finest; finSub = info.finest_sub;
    } else {
      // finite-ness of everything we are about to adopt (NaN-blind checks are a known trap)
      bool allOk = allFinite(sig, 6) && allFinite(&C[0][0], 36) &&
                 allFinite(np1.eps_e, 6) && std::isfinite(np1.pi_i) && std::isfinite(np1.v) &&
                 std::isfinite(np1.v0) && std::isfinite(np1.eps_p_v) && std::isfinite(np1.eps_p_s) &&
                 std::isfinite(np1.D_last) && std::isfinite(np1.eps_f_v) && std::isfinite(np1.W_f) &&
                 std::isfinite(info.E_f);
      if (!allOk) refusal = -1;
      else {
        sT = np1;
        for (int i = 0; i < 6; i++) sigT[i] = sig[i];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) CT[i][j] = C[i][j];
        efT = efC + info.E_f;        // the floor energy of this increment (a diagnostic; never fed back, sheet 9.7)
        trialRefused = false;
        return 0;
      }
    }
  }

  // ---- refused: committed state untouched, trial frozen at the committed state ----
  trialRefused = true;
  lastRefusal = refusal;
  if (newCause) { lastFinest = fin; lastFinestSub = finSub; }
  nRefusals++;
  sT = sC;
  efT = efC;                 // a refused increment counts nothing (sheet 9.7)
  ladruno_norsand::stress(kp, sC, sigT);
  ladruno_norsand::elasticTangent(kp, sC, CT);
  if (!warnedTrial) {
    warnedTrial = true;
    const std::string fx = finestText(lastFinest, lastFinestSub);
    opserr << "WARNING LadrunoNorSand tag " << this->getTag() << ": the return map REFUSED this trial strain ("
           << refusalName(refusal) << (fx.empty() ? "" : "; ") << fx.c_str() << "; local iters "
           << lastInfo.local_iters << ", pi iters " << lastInfo.pi_iters << ", substeps " << lastInfo.substeps
           << "). The committed state is untouched"
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

// LadrunoElasticStrainProvider: the TRIAL elastic strain (the one belonging to sigT), tensor -> engineering Voigt.
// A refused trial leaves sT = sC (frozen at the committed state), so this stays consistent with the stress
// getStress() returns. // Ladruno WP-144 G2
bool LadrunoNorSand::ladrunoGetElasticStrain(Vector& epsE) const
{
  if (!initOk || epsE.Size() != 6 || !allFinite(sT.eps_e, 6)) return false;
  // p' floor (sheet 9.7): sT.eps_e is the POST-FLOORED elastic strain, the one that reproduces sigT exactly, so
  // the wrapper rebuilds b^e from the floored strain (Delta eps^f joins the inelastic part of the split). Under HAR
  // it lies in dom Psi by construction (a trial outside it is a floor event or a refusal, never committed); the
  // check below only keeps a wrapper from ever receiving an out-of-domain strain.
  if (kp.energy == 1 &&
      !(ladruno_norsand::detail::har_estar(kp, (sT.eps_e[0] + sT.eps_e[1]) + sT.eps_e[2]) > 0.0)) return false;
  for (int i = 0; i < 6; i++) epsE(i) = (i >= 3) ? 2.0 * sT.eps_e[i] : sT.eps_e[i];
  return true;
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
  // the hyperelastic tangent at the INITIAL state s0 (not the base default getTangent(): -initial must
  // stay a genuine initial-stiffness iteration). s0 is fixed at construction (buildInitialState) and is
  // touched by nothing afterwards (not by a step, a commit, a revert or a latch), so this is the elastic
  // tangent at (sigma0, v0, pi_i0) however much plastic history the point has. A G2 test kills the
  // "return getTangent()" mutant (tests/test_ladruno_norsand.py).
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
  efT = efC;
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
               << ": a REFUSED update (" << refusalName(lastRefusal)
               << (lastFinest != 0 ? "; " : "") << finestText(lastFinest, lastFinestSub).c_str()
               << ") reached commitState -- the host"
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

  // nSubstepped counts COMMITTED steps that needed substepping (lastInfo is the converged trial's);
  // counting per trial would count every Newton re-integration of the same step.
  if (lastInfo.substeps > 1) nSubstepped++;
  sC = sT;
  efC = efT;
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
  efC = 0.0; efT = 0.0;
  trialRefused = false;
  latched = false;
  warnedTrial = false;
  warnedLatch = false;
  lastRefusal = 0;
  lastFinest = 0; lastFinestSub = 0;
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
  v0init = o.v0init; pi0init = o.pi0init; F0init = o.F0init; pi0Auto = o.pi0Auto;
  efC = o.efC; efT = o.efT;          // the floor-energy census of the committed path: carried by value (like nSubstepped)
  initOk = o.initOk; initMsg = o.initMsg;
  s0 = o.s0; sC = o.sC; sT = o.sT;
  nSubstepped = o.nSubstepped;       // history census of the committed path: carried by value
  // REFUSAL STATE IS NOT INHERITED (WP-144 G2 decision, LEDGER_quirks "getCopy of a latched LadrunoNorSand").
  // A clone starts with no latch, no refused trial and an empty refusal record, whatever the source holds.
  //   * The latch exists to stop THIS integration point from continuing after a commit the host discarded.
  //     A refused trial never reaches the committed state (the commit is aborted and the trial restored from
  //     the committed state), so the source's committed state, which is all the clone copies, is a valid
  //     state: there is nothing for the clone to be protected from.
  //   * A clone that inherited the latch would make a FRESH element refuse every step from its first update()
  //     and fail with no refusal of its own to explain it (the cause belongs to another point).
  //   * It is the opposite of a "work already done" latch (ASDConcrete3D regularization): not copying that one
  //     redoes the work, not copying this one redoes nothing.
  // sendSelf/recvSelf is a different case: it moves THE SAME point (a restart), so it keeps the latch.
  trialRefused = false;
  latched = false;
  warnedTrial = false; warnedLatch = false;
  lastRefusal = 0; lastFinest = 0; lastFinestSub = 0;
  nRefusals = 0;
  lastInfo = StepInfo();             // the last attempt's info belongs to the source's (possibly refused) trial
  // the trial of a refused source is frozen at its committed state: make that explicit for the clone
  if (o.trialRefused || o.latched) this->restoreTrialFromCommitted();
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
  for (int i = 0; i < NPD_BASE; i++) data(nsw::PARAMS + i) = kp.*(kParam[i].mem);
  for (int i = NPD_BASE; i < NPD; i++) data(nsw::XPARAMS + (i - NPD_BASE)) = kp.*(kParam[i].mem);
  data(nsw::ENERGY) = kp.energy;
  data(nsw::CSL) = kp.csl_mode;
  data(nsw::ZETA) = kp.zeta;
  data(nsw::CAP) = kp.cap;
  for (int i = 0; i < 6; i++) data(nsw::SIGMA0 + i) = sigma0[i];
  data(nsw::V0) = v0init;
  data(nsw::PI0) = pi0init;
  {
    double a[nsw::STATE_LEN];
    packState(a, sC); for (int i = 0; i < nsw::STATE_LEN; i++) data(nsw::SC + i) = a[i];
    packState(a, sT); for (int i = 0; i < nsw::STATE_LEN; i++) data(nsw::ST + i) = a[i];
    double f[nsw::FLOOR_LEN];
    packFloor(f, sC); for (int i = 0; i < nsw::FLOOR_LEN; i++) data(nsw::SCF + i) = f[i];
    packFloor(f, sT); for (int i = 0; i < nsw::FLOOR_LEN; i++) data(nsw::STF + i) = f[i];
  }
  data(nsw::EFC) = efC;
  data(nsw::EFT) = efT;
  data(nsw::PI0AUTO) = pi0Auto ? 1.0 : 0.0;
  for (int i = 0; i < 6; i++) { data(nsw::EPSC + i) = epsC[i]; data(nsw::EPST + i) = epsT[i]; }
  // the trial stress and tangent are NOT sent (see the wire-layout note at nsw): Domain::recvSelf's update() recomputes them
  data(nsw::LATCHED) = latched ? 1.0 : 0.0;
  data(nsw::TRIALREF) = trialRefused ? 1.0 : 0.0;
  data(nsw::LASTREF) = lastRefusal;
  data(nsw::LASTFIN) = lastFinest;
  data(nsw::LASTFINSUB) = lastFinestSub;
  data(nsw::NREF) = nRefusals;
  data(nsw::NSUB) = nSubstepped;
  data(nsw::INFO + 0) = lastInfo.refusal;
  data(nsw::INFO + 1) = lastInfo.plastic;
  data(nsw::INFO + 2) = lastInfo.vertex;
  data(nsw::INFO + 3) = lastInfo.cap_active;
  data(nsw::INFO + 4) = lastInfo.local_iters;
  data(nsw::INFO + 5) = lastInfo.pi_iters;
  data(nsw::INFO + 6) = lastInfo.substeps;
  data(nsw::INFO + 7) = lastInfo.finest;
  data(nsw::INFO + 8) = lastInfo.finest_sub;
  data(nsw::INFOF + 0) = lastInfo.floor_tr;
  data(nsw::INFOF + 1) = lastInfo.floor_post;
  data(nsw::INFOF + 2) = lastInfo.at_floor;
  data(nsw::INFOF + 3) = lastInfo.deps_f_v;
  data(nsw::INFOF + 4) = lastInfo.W_f;
  data(nsw::INFOF + 5) = lastInfo.E_f;

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
  for (int i = 0; i < NPD_BASE; i++) kp.*(kParam[i].mem) = data(nsw::PARAMS + i);
  for (int i = NPD_BASE; i < NPD; i++) kp.*(kParam[i].mem) = data(nsw::XPARAMS + (i - NPD_BASE));
  kp.energy = (int)data(nsw::ENERGY);
  kp.csl_mode = (int)data(nsw::CSL);
  kp.zeta = (int)data(nsw::ZETA);
  kp.cap = (int)data(nsw::CAP);
  for (int i = 0; i < 6; i++) sigma0[i] = data(nsw::SIGMA0 + i);
  v0init = data(nsw::V0);
  pi0init = data(nsw::PI0);

  // rebuild s0 (the initial state is a function of the inputs), then overwrite committed/trial
  this->setupDim();
  this->buildInitialState();

  double a[nsw::STATE_LEN];
  for (int i = 0; i < nsw::STATE_LEN; i++) a[i] = data(nsw::SC + i);
  unpackState(a, sC);
  for (int i = 0; i < nsw::STATE_LEN; i++) a[i] = data(nsw::ST + i);
  unpackState(a, sT);
  {
    double f[nsw::FLOOR_LEN];
    for (int i = 0; i < nsw::FLOOR_LEN; i++) f[i] = data(nsw::SCF + i);
    unpackFloor(f, sC);
    for (int i = 0; i < nsw::FLOOR_LEN; i++) f[i] = data(nsw::STF + i);
    unpackFloor(f, sT);
  }
  efC = data(nsw::EFC);
  efT = data(nsw::EFT);
  pi0Auto = (data(nsw::PI0AUTO) != 0.0);
  for (int i = 0; i < 6; i++) { epsC[i] = data(nsw::EPSC + i); epsT[i] = data(nsw::EPST + i); }
  latched = (data(nsw::LATCHED) != 0.0);
  trialRefused = (data(nsw::TRIALREF) != 0.0);
  lastRefusal = (int)data(nsw::LASTREF);
  lastFinest = (int)data(nsw::LASTFIN);
  lastFinestSub = (int)data(nsw::LASTFINSUB);
  nRefusals = (int)data(nsw::NREF);
  nSubstepped = (int)data(nsw::NSUB);
  lastInfo.refusal = (int)data(nsw::INFO + 0);
  lastInfo.plastic = (int)data(nsw::INFO + 1);
  lastInfo.vertex = (int)data(nsw::INFO + 2);
  lastInfo.cap_active = (int)data(nsw::INFO + 3);
  lastInfo.local_iters = (int)data(nsw::INFO + 4);
  lastInfo.pi_iters = (int)data(nsw::INFO + 5);
  lastInfo.substeps = (int)data(nsw::INFO + 6);
  lastInfo.finest = (int)data(nsw::INFO + 7);
  lastInfo.finest_sub = (int)data(nsw::INFO + 8);
  lastInfo.floor_tr = (int)data(nsw::INFOF + 0);
  lastInfo.floor_post = (int)data(nsw::INFOF + 1);
  lastInfo.at_floor = (int)data(nsw::INFOF + 2);
  lastInfo.deps_f_v = data(nsw::INFOF + 3);
  lastInfo.W_f = data(nsw::INFOF + 4);
  lastInfo.E_f = data(nsw::INFOF + 5);
  // sigT / CT are not on the wire: rebuild them from the restored trial State (hyperelastic stress and tangent
  // of sT; the consistent plastic tangent of the source's last step is NOT reproduced, and the host's own
  // update() recomputes both anyway).
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
  if (kp.energy == 1) {
    const double nuAxis = (3.0 * kp.k - 2.0 * kp.g) / (6.0 * kp.k + 2.0 * kp.g);
    s << "  elastic (HAR energy, Houlsby-Amorosi-Rojas 2005): k=" << kp.k << " g=" << kp.g << " n=" << kp.n_e
      << " p_a=" << kp.p_a << "  (isotropic axis: K = k p_a (|p|/p_a)^n, G = g p_a (|p|/p_a)^n, nu = " << nuAxis
      << "; p = -p_a at zero elastic strain; the BA06 constants are not read)" << endln;
  } else {
    s << "  elastic (BA06 energy): p0=" << kp.p0 << " kappa_hat=" << kp.kappa_hat << " eps_v0=" << kp.eps_v0
      << " mu0=" << kp.mu0 << " alpha0=" << kp.alpha0 << endln;
  }
  s << "  energy option: " << (kp.energy == 1 ? "HAR" : "BA06 (default; paper mode)") << ", p_ref=" << ladruno_norsand::pRef(kp)
    << (kp.energy == 1 ? " (= p_a)" : " (= |p0|)") << ": F_tol, r4 and the floor default scale with it" << endln;
  if (kp.p_min > 0.0)
    s << "  p' floor: ON, p_min=" << kp.p_min << " (" << kp.p_min / ladruno_norsand::pRef(kp) << " p_ref): strain-space"
      << " projection (sheet 9.7 S.48) at the trial and after convergence, never inside the local Newton; it never"
      << " refuses and is COUNTED (responses `floor`, `floorEnergy`); the consistent tangent is the EXACT one (zero bulk"
      << " stiffness at a floored state, no regularisation)" << endln;
  else
    s << "  p' floor: OFF (-pmin 0): a trial with p >= 0 (HAR: outside the elastic domain) is a REFUSAL, as before round 3"
      << endln;
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
    << ", pi_i0=" << s0.pi_i << ", psi_i0=" << cslPsi(kp, s0.v, s0.pi_i)
    << ", F0=" << F0init
    << (F0init > F0_ON_SURFACE_REL * ladruno_norsand::pRef(kp) ? " (OUTSIDE the surface)"
        : (std::fabs(F0init) <= F0_ON_SURFACE_REL * ladruno_norsand::pRef(kp) ? " (ON the surface)" : " (inside the surface)")) << endln;
  if (pi0Auto)
    s << "  pi_i0 from the unified rule (sheet 5.4 (S.53)): the surface through (p_init, max(eta_init, c2 M)), c2 = 0 without a cap"
      << endln;
  if (s0.n_f_init != 0)
    s << "  the initial state was PROJECTED by the p' floor (n_f_init = 1, eps^f_v = " << s0.eps_f_v << "): sigma0 is replaced by"
      << " sigma(eps^e_f) (p = -p_min); the first equilibrium iteration absorbs the O(p_min) imbalance" << endln;
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
  s << "  refused trials=" << nRefusals << " committed substepped steps=" << nSubstepped
    << " latched=" << (latched ? 1 : 0) << endln;
  s << "  p' floor: p_min=" << kp.p_min << " at_floor=" << sC.at_floor << " n_f_tr=" << sC.n_f_tr << " n_f_post=" << sC.n_f_post
    << " n_f_init=" << sC.n_f_init << " eps_f_v=" << sC.eps_f_v << " W_f=" << sC.W_f << " E_f(cumulative)=" << efC << endln;
}

// ===========================================================================
//  recordable responses
//    stress, strain, tangent         element-facing (Voigt, engineering shear)
//    state          [pi_i, psi_i, v, v0, eps_p_v, eps_p_s]
//    D              dissipation of the last step (>= 0)
//    refusal        [last refusal code (-1 non-finite), refused trials, latched, finest, finest_sub]
//    substeps       [substeps of the last step, committed steps that needed substepping]
//    stepInfo       [refusal, plastic, vertex, cap_active, local_iters, pi_iters, substeps, finest, finest_sub,
//                    floor_tr, floor_post]   (the last two: sheet 9.7, floor events of the last step, summed over its sub-increments)
//    psi            state parameter psi = v - v_c(p)
//    elasticStrain  elastic strain, Voigt engineering shear (the POST-FLOORED one: it reproduces the stress)
//    floor          (alias floored) [at_floor, n_f_tr, n_f_post, eps_f_v, W_f] of the TRIAL state (sheet 9.7): at_floor = 1 if
//                   p > -p_min (1 + 1e-10); n_f_* the CUMULATIVE numbers of trial / post floor events; eps_f_v = sum of the
//                   floor's volumetric strain; W_f = p_min eps_f_v, the bound on the energy the floor injected. Always counted.
//    floorEnergy    [E_f of the last increment, E_f cumulative over the history incl. the trial increment, W_f]: the floor
//                   energy (S.52) from the closed-form Psi (BA06 (S.4) / HAR (S.4h)), a DIAGNOSTIC (never fed back);
//                   0 <= E_f <= W_f. The initial-state event (n_f_init) is in W_f but has no E_f.
//    floorInit      [n_f_init] 1 if the p' floor projected the initial state
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
    return new MaterialResponse(this, 6, Vector(5));
  if (strcmp(a, "substeps") == 0)
    return new MaterialResponse(this, 7, Vector(2));
  if (strcmp(a, "stepInfo") == 0)
    return new MaterialResponse(this, 8, Vector(11));
  if (strcmp(a, "psi") == 0)
    return new MaterialResponse(this, 9, Vector(1));
  if (strcmp(a, "elasticStrain") == 0 || strcmp(a, "elasticStrains") == 0)
    return new MaterialResponse(this, 10, Vector(ncomp));
  if (strcmp(a, "floor") == 0 || strcmp(a, "floored") == 0)
    return new MaterialResponse(this, 11, Vector(5));
  if (strcmp(a, "floorEnergy") == 0)
    return new MaterialResponse(this, 12, Vector(3));
  if (strcmp(a, "floorInit") == 0)
    return new MaterialResponse(this, 13, Vector(1));

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
        v(3) = lastFinest; v(4) = lastFinestSub;
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
        v(7) = lastInfo.finest; v(8) = lastInfo.finest_sub;
        v(9) = lastInfo.floor_tr; v(10) = lastInfo.floor_post;
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
    case 11:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = sT.at_floor; v(1) = sT.n_f_tr; v(2) = sT.n_f_post; v(3) = sT.eps_f_v; v(4) = sT.W_f;
      }
      return 0;
    case 12:
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = lastInfo.E_f; v(1) = efT; v(2) = sT.W_f;
      }
      return 0;
    case 13:
      if (matInfo.theVector) (*(matInfo.theVector))(0) = sT.n_f_init;
      return 0;
    default:
      return -1;
  }
}

// ===========================================================================
//  OPS parser (Tcl and Python)
//
//   nDMaterial LadrunoNorSand tag
//       <-energy BA06|HAR>                                                        (default BA06)
//         BA06: -p0 p0 -kappa_hat kh -mu0 mu0 [-eps_v0 e] [-alpha0 a]
//         HAR : -k k -g g [-n n] -p_a pa                (the stiffness factors, directly), or
//               -G0 G0 -nu nu [-e_ref e] -p_a pa        (the DM04 mapping of sheet 2.3 / 15, exact on the axis)
//       -M M -N N [-N_bar Nb] -rho rho [-rho_bar rb] -chi chi -h h
//       [-csl paper|fork]   paper: -lambda_tilde lt -v_c0 vc0     fork: -e0 e0 -lambda_c lc -xi xi [-p_a pa]
//       [-zeta WW|GA] [-cap none|planar|smooth [-c1 c1] [-c2 c2]]
//       [-pmin p]                                                                  (default 5e-3 p_ref; 0 = off)
//       -v0 v0 (-pi0 pi_i0 | -pi0_auto) [-sigma0 s11 s22 s33 s12 s23 s13] [-density d]
//
//   defaults: N_bar = N, rho_bar = rho, eps_v0 = alpha0 = 0, HAR n = 1/2, e_ref = v0 - 1, -csl paper, -zeta WW,
//   -cap none, -sigma0 = isotropic p0 (HAR: p0 := -p_a), -density 0, -pmin = 5e-3 p_ref with p_ref = |p0| (BA06) or
//   p_a (HAR) (sheet 9.7: 0.5 / 0.505 kPa on the K2 / TIMs sets).
//   REQUIRED: -p0 -kappa_hat -mu0 (BA06) or -k -g / -G0 -nu, -p_a (HAR), -M -N -rho -chi -h, the active CSL's constants
//   (fork: -e0 -lambda_c -xi AND -p_a: it multiplies a stress, so a unit-blind default would be silently wrong), -v0, and
//   -pi0 (no default onto the yield surface: F(sigma0, pi_i0) is computed; a start OUTSIDE the surface is refused, an
//   on-surface start warned) OR -pi0_auto (the unified rule of sheet 5.4 (S.53); giving both is refused).
//   p_a is ONE flag: the HAR energy's reference pressure and the fork CSL's are the same number (TIMs: 101 kPa). Under
//   BA06 with the paper CSL it is unused (a note is printed).
//   Energy refusals (sheet 2.4; codes 201..207, then the kernel's validate()): a BA06 constant (-p0 -kappa_hat
//   -eps_v0 -mu0 -alpha0) with -energy HAR; a HAR constant (-k -g -n -G0 -nu -e_ref) with BA06; both stiffness forms
//   (-k/-g and -G0/-nu); one of a pair; an out-of-range DM04 input; then k <= 0, g <= 0, n outside [0, 1), p_a <= 0,
//   p_min < 0 (kernel codes 21..25). HAR replaces alpha0: none of the five BA06 values is read, so none is accepted.
//   Units: every stress-like input (-p0, -mu0, -p_a, -sigma0, -pi0, -pmin) and the paper CSL's -v_c0 (the intercept
//   of v_c = v_c0 - lambda_tilde ln(-p)) are in ONE consistent stress unit chosen by the model; k, g, n, G0, nu, e_ref
//   are dimensionless.
//   -rho is the ELLIPTICITY of F (the model parameter); the mass density is -density.
//   The kernel's validate() is called: a refused parameter set is a hard error, rho > rho_bar a warning.
// ===========================================================================
static void nsUsage()
{
  opserr << "Want: nDMaterial LadrunoNorSand tag? <-energy BA06|HAR>"
         << " (BA06: -p0 p0? -kappa_hat kh? -mu0 mu0? <-eps_v0 e?> <-alpha0 a?>;"
         << " HAR: -k k? -g g? <-n n?> -p_a pa?  or  -G0 G0? -nu nu? <-e_ref e?> -p_a pa?)"
         << " -M M? -N N? <-N_bar Nb?> -rho rho? <-rho_bar rb?> -chi chi? -h h?"
         << " <-csl paper|fork> (paper: -lambda_tilde lt? -v_c0 vc0?  fork: -e0 e0? -lambda_c lc? -xi xi? -p_a pa?)"
         << " <-zeta WW|GA> <-cap none|planar|smooth> <-c1 c1?> <-c2 c2?> <-pmin p?>"
         << " -v0 v0? (-pi0 pi_i0? | -pi0_auto) <-sigma0 s11? s22? s33? s12? s23? s13?> <-density d?>"
         << " (-pi0 or -pi0_auto is REQUIRED; fork CSL and HAR also require -p_a; -pmin 0 switches the floor off)" << endln;
}

static bool ieq(const char* a, const char* b)
{
  for (; *a && *b; a++, b++)
    if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return false;
  return *a == 0 && *b == 0;
}

// Copy a token returned by OPS_GetString into a local buffer AT ONCE (the LadrunoSANISAND pattern): the
// openseespy backend returns a pointer into a temporary it has already released, and returns the literal
// "Invalid String Input!" (still consuming the token) when the argument is a number. Nothing may be
// fetched between OPS_GetString and this copy.
static void copyTok(char* dst, int cap, const char* src)
{
  int i = 0;
  if (src != 0) while (i < cap - 1 && src[i] != '\0') { dst[i] = src[i]; i++; }
  dst[i] = '\0';
}

// A parser-level refusal (sheet 2.4): printed with its code, the command returns 0.
static void* nsRefuse(int tag, int code, const std::string& why)
{
  opserr << "WARNING LadrunoNorSand tag " << tag << ": parameter set REFUSED (code " << code << "): " << why.c_str() << endln;
  return 0;
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
  double initSigma[6] = {0, 0, 0, 0, 0, 0};
  bool haveSigma0 = false, haveV0 = false, havePi0 = false, wantPi0Auto = false;
  double v0 = 0.0, pi0 = std::numeric_limits<double>::quiet_NaN(), massDensity = 0.0;
  // the DM04 form of the HAR stiffness (sheet 2.3): -G0 -nu [-e_ref]
  double dmG0 = 0.0, dmNu = 0.0, dmEref = 0.0;
  bool haveG0 = false, haveNu = false, haveEref = false;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* rawFlag = OPS_GetString();      // consumes exactly one token in every backend
    if (rawFlag == 0) break;                    // classic-Tcl backend returns 0 when exhausted
    char flagTok[64];
    copyTok(flagTok, (int)sizeof(flagTok), rawFlag);
    const char* flag = flagTok;

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
    else if (strcmp(flag, "-csl") == 0 || strcmp(flag, "-zeta") == 0 || strcmp(flag, "-cap") == 0 ||
             strcmp(flag, "-energy") == 0) {
      const char* rawVal = OPS_GetString();
      if (rawVal == 0) { opserr << "WARNING LadrunoNorSand: " << flag << " wants a mode\n"; return 0; }
      char valTok[32];
      copyTok(valTok, (int)sizeof(valTok), rawVal);
      const char* val = valTok;
      if (strcmp(flag, "-energy") == 0) {
        if (ieq(val, "BA06")) p.energy = 0;
        else if (ieq(val, "HAR")) p.energy = 1;
        else { opserr << "WARNING LadrunoNorSand: -energy wants BA06|HAR, got '" << val << "'\n"; return 0; }
      } else if (strcmp(flag, "-csl") == 0) {
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
    else if (strcmp(flag, "-G0") == 0 || strcmp(flag, "-nu") == 0 || strcmp(flag, "-e_ref") == 0) {
      double d;
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &d) < 0) { opserr << "WARNING LadrunoNorSand: " << flag << " wants a number\n"; return 0; }
      if (strcmp(flag, "-G0") == 0) { dmG0 = d; haveG0 = true; }
      else if (strcmp(flag, "-nu") == 0) { dmNu = d; haveNu = true; }
      else { dmEref = d; haveEref = true; }
    }
    else if (strcmp(flag, "-sigma0") == 0) {
      numData = 6;
      if (OPS_GetDoubleInput(&numData, initSigma) < 0) {
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
    else if (strcmp(flag, "-pi0_auto") == 0) {
      wantPi0Auto = true;                       // no value: the unified rule of sheet 5.4 (S.53)
    }
    else if (strcmp(flag, "-density") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &massDensity) < 0) { opserr << "WARNING LadrunoNorSand: -density wants a number\n"; return 0; }
    }
    else {
      opserr << "WARNING LadrunoNorSand: unknown flag '" << flag << "'\n";
      nsUsage();
      return 0;
    }
  }

  // ---- energy option (sheet 2.4): a constant of the other energy is REFUSED, never ignored ----
  const bool har = (p.energy == 1);
  if (har) {
    static const int baIdx[] = {I_P0, I_KAPPA, I_EPSV0, I_MU0, I_ALPHA0};
    std::string given;
    for (size_t i = 0; i < sizeof(baIdx) / sizeof(int); i++)
      if (seen[baIdx[i]]) { given += " "; given += kParam[baIdx[i]].flag; }
    if (!given.empty())
      return nsRefuse(tag, PARSE_BA06_WITH_HAR,
                      "-energy HAR REPLACES the BA06 constants (alpha0 included) and does not read them; given:" + given +
                      ". Remove them: they are refused, never ignored.");
  } else {
    std::string given;
    if (seen[I_K]) given += " -k";
    if (seen[I_G]) given += " -g";
    if (seen[I_NE]) given += " -n";
    if (haveG0) given += " -G0";
    if (haveNu) given += " -nu";
    if (haveEref) given += " -e_ref";
    if (!given.empty())
      return nsRefuse(tag, PARSE_HAR_WITH_BA06,
                      "these are HAR constants and the BA06 energy does not read them; given:" + given +
                      ". Add -energy HAR, or remove them.");
  }
  if (wantPi0Auto && havePi0)
    return nsRefuse(tag, PARSE_PI0_BOTH, "-pi0 and -pi0_auto both given: -pi0 is the deck's value, -pi0_auto asks for the"
                    " unified rule of sheet 5.4; give one");

  bool dm04Pending = false;      // -G0 -nu given, e_ref needs -v0 (default v0 - 1), which may not be given yet
  if (har) {
    const bool haveKG = seen[I_K] || seen[I_G];
    const bool haveDM = haveG0 || haveNu || haveEref;
    if (haveKG && haveDM)
      return nsRefuse(tag, PARSE_HAR_AMBIGUOUS, "give the HAR stiffness either directly (-k, -g) or by the DM04 mapping"
                      " (-G0, -nu, -e_ref), not both");
    if (haveKG && !(seen[I_K] && seen[I_G]))
      return nsRefuse(tag, PARSE_HAR_INCOMPLETE, "-k and -g go together (HAR needs both stiffness factors)");
    if (haveDM) {
      if (!(haveG0 && haveNu))
        return nsRefuse(tag, PARSE_HAR_INCOMPLETE, "-G0 and -nu go together (the DM04 mapping needs both; -e_ref is optional)");
      if (haveV0 || haveEref) {
        const double eref = haveEref ? dmEref : v0 - 1.0;
        if (!(std::isfinite(dmG0) && dmG0 > 0.0) || !(std::isfinite(dmNu) && dmNu > -1.0 && dmNu < 0.5) ||
            !(std::isfinite(eref) && eref > -1.0)) {
          char buf[256];
          std::snprintf(buf, sizeof(buf), "DM04 mapping out of range: need -G0 > 0, -1 < -nu < 1/2, e_ref > -1 (got G0 = %.6g,"
                        " nu = %.6g, e_ref = %.6g)", dmG0, dmNu, eref);
          return nsRefuse(tag, PARSE_DM04_RANGE, buf);
        }
        const double f = (2.97 - eref) * (2.97 - eref) / (1.0 + eref);          // DM04 f(e) = (2.97 - e)^2 / (1 + e)
        p.g = dmG0 * f;
        p.k = p.g * 2.0 * (1.0 + dmNu) / (3.0 * (1.0 - 2.0 * dmNu));              // exact on the isotropic axis
        seen[I_K] = true; seen[I_G] = true;
        opserr << "LadrunoNorSand tag " << tag << ": HAR stiffness from the DM04 mapping (sheet 2.3/15): G0=" << dmG0
               << " nu=" << dmNu << " e_ref=" << eref << (haveEref ? "" : " (= v0 - 1)") << " -> g=" << p.g << " k=" << p.k
               << " (exact on the isotropic axis; off the axis the moduli differ from DM04's)" << endln;
      } else {
        dm04Pending = true;
      }
    }
    if (!seen[I_NE]) p.n_e = 0.5;       // TIMs / Toyoura
    p.p0 = -p.p_a;                      // sheet 2.4: p0 := -p_a (the reference pressure of the 9.1 scalings)
  } else if (seen[I_PA] && p.csl_mode == 0) {
    opserr << "LadrunoNorSand tag " << tag << ": NOTE -p_a is not used (BA06 energy with the paper CSL)" << endln;
  }

  // ---- defaults that depend on other parameters ----
  if (!seen[I_NBAR])   p.N_bar = p.N;
  if (!seen[I_RHOBAR]) p.rho_bar = p.rho;
  if (p.cap == 1 && !seen[I_C2] && seen[I_C1]) p.c2 = p.c1;      // planar: c1 = c2
  if (!seen[I_PMIN]) p.p_min = ladruno_norsand::defaultPmin(p);    // 5e-3 p_ref (sheet 9.7); -pmin 0 switches it off

  // ---- required parameters ----
  {
    static const int reqAll[] = {I_M, I_N, I_RHO, I_CHI, I_H};
    static const int reqBA06[] = {I_P0, I_KAPPA, I_MU0};
    static const int reqPaper[] = {I_LTILDE, I_VC0};
    static const int reqFork[] = {I_E0, I_LC, I_XI, I_PA};
    bool missing = false;
    for (size_t i = 0; i < sizeof(reqAll) / sizeof(int); i++)
      if (!seen[reqAll[i]]) { opserr << "WARNING LadrunoNorSand: missing required " << kParam[reqAll[i]].flag << "\n"; missing = true; }
    if (har) {
      if (!seen[I_K] && !seen[I_G] && !dm04Pending)
        { opserr << "WARNING LadrunoNorSand: -energy HAR needs -k and -g (or the DM04 form -G0 and -nu)\n"; missing = true; }
      if (dm04Pending)
        { opserr << "WARNING LadrunoNorSand: the DM04 mapping needs -e_ref or -v0 (e_ref defaults to v0 - 1)\n"; missing = true; }
      if (!seen[I_PA])
        { opserr << "WARNING LadrunoNorSand: -energy HAR needs -p_a (the reference pressure; one flag shared with the fork CSL)\n"; missing = true; }
    } else {
      for (size_t i = 0; i < sizeof(reqBA06) / sizeof(int); i++)
        if (!seen[reqBA06[i]]) { opserr << "WARNING LadrunoNorSand: missing required " << kParam[reqBA06[i]].flag << "\n"; missing = true; }
    }
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
    if (!havePi0 && !wantPi0Auto) {
      opserr << "WARNING LadrunoNorSand: missing required -pi0 (initial image pressure pi_i0 < 0; there is no default"
             << " onto the yield surface: F(sigma0, pi_i0) must be <= 0, a more negative -pi0 starts further inside;"
             << " or give -pi0_auto for the unified rule of sheet 5.4)\n";
      missing = true;
    }
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
    for (int i = 0; i < 3; i++) initSigma[i] = p.p0;     // isotropic at the reference pressure (HAR: -p_a)

  // -pi0_auto: pi0 stays NaN, which the constructor reads as "the unified rule" (sheet 5.4)
  LadrunoNorSand* mat = new LadrunoNorSand(tag, p, initSigma, v0, wantPi0Auto ? std::numeric_limits<double>::quiet_NaN() : pi0,
                                           massDensity);
  if (!mat->initOK()) {
    opserr << "WARNING LadrunoNorSand tag " << tag << ": the initial state was REFUSED: " << mat->initMessage() << endln;
    delete mat;
    return 0;
  }
  if (std::fabs(mat->initF0()) <= F0_ON_SURFACE_REL * ladruno_norsand::pRef(p))
    opserr << "WARNING LadrunoNorSand tag " << tag << ": the initial state is ON the yield surface (F0 = " << mat->initF0()
           << ", |F0| <= " << F0_ON_SURFACE_REL << " p_ref): the first loading increment is plastic from the start."
           << " A more negative -pi0 starts inside the surface." << endln;
  mat->echoParameters(opserr);
  return mat;
}

void* OPS_LadrunoNorSand(void)
{
  return parseLadrunoNorSand();
}
