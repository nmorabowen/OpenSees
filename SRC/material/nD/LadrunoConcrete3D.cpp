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

// Ladruno: CDPM2-grade solid-concrete plastic-damage nDMaterial (3D). See LadrunoConcrete3D.h and
// Ladruno_implementation/31_ladruno_concrete3d_adr.md. classTag 33017.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoConcrete3D.h>
#include <LadrunoConcrete3DKernel.h>   // the header-only OpenSees-free CDPM2 kernel
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <Element.h>                   // ops_TheActiveElement->getCharacteristicLength()
#include <OPS_Globals.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

using Ladruno::Concrete3D::Params;
using Ladruno::Concrete3D::State;

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoConcrete3D tag? E? nu? fc? ft? Gf? Gc?
//       <-e e? | -kupfer ratio?> <-Df Df?> <-As As?> <-rho rho?>
//       <-hardening qh0? Hp?> <-ductility Ah? Bh? Ch? Dh?>
//       <-lch lch?> <-autoRegularization>
//  (fc, ft are POSITIVE magnitudes; compression is negative in the model.)
// ===========================================================================
void* OPS_LadrunoConcrete3D(void)
{
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoConcrete3D tag? E? nu? fc? ft? Gf? Gc? "
           << "<-e e? | -kupfer ratio?> <-Df Df?> <-As As?> <-rho rho?> "
           << "<-hardening qh0? Hp?> <-ductility Ah? Bh? Ch? Dh?> <-lch lch?> <-autoRegularization> "
           << "<-implex> <-eta eta?> <-ctTemper none|alphat|proj> <-hoop K? <-hoopFy fy?>>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoConcrete3D tag\n";
    return 0;
  }

  double d6[6];
  numData = 6;
  if (OPS_GetDoubleInput(&numData, d6) < 0) {
    opserr << "WARNING invalid LadrunoConcrete3D E nu fc ft Gf Gc\n";
    return 0;
  }
  double E = d6[0], nu = d6[1], fc = d6[2], ft = d6[3], Gf = d6[4], Gc = d6[5];

  // defaults (oracle make_material / ADR 6): e from Kupfer 1.16; CDPM2 published ductility table
  double ecc = -1.0;            // <0 => derive from Kupfer ratio below
  double kupfer = 1.16;
  bool haveE = false;
  double Df = 1.0, As = 2.0, rho = 0.0;
  double qh0 = 0.3, Hp = 0.5, Ah = 0.08, Bh = 0.003, Ch = 2.0, Dh = 1.0e-6;
  double lch = 1.0;
  bool autoReg = false;
  bool implex = false;
  double eta = 0.0;
  int ctTemper = 0;            // P2h compression->tension damage temper: 0=none 1=alphat 2=proj
  double hoopK = 0.0;          // P5 confined-fiber (BeamFiber view): passive hoop confining stiffness
  double hoopFy = 1.0e30;      //                                     hoop yield (caps the confining pressure)

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-e") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &ecc) < 0) { opserr << "WARNING LadrunoConcrete3D: -e wants e\n"; return 0; }
      haveE = true;
    }
    else if (strcmp(flag, "-kupfer") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &kupfer) < 0) { opserr << "WARNING LadrunoConcrete3D: -kupfer wants fcc/fc\n"; return 0; }
    }
    else if (strcmp(flag, "-Df") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &Df) < 0) { opserr << "WARNING LadrunoConcrete3D: -Df wants Df\n"; return 0; }
    }
    else if (strcmp(flag, "-As") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &As) < 0) { opserr << "WARNING LadrunoConcrete3D: -As wants As\n"; return 0; }
    }
    else if (strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) { opserr << "WARNING LadrunoConcrete3D: -rho wants rho\n"; return 0; }
    }
    else if (strcmp(flag, "-hardening") == 0) {
      double h[2]; numData = 2;
      if (OPS_GetDoubleInput(&numData, h) < 0) { opserr << "WARNING LadrunoConcrete3D: -hardening wants qh0 Hp\n"; return 0; }
      qh0 = h[0]; Hp = h[1];
    }
    else if (strcmp(flag, "-ductility") == 0) {
      double h[4]; numData = 4;
      if (OPS_GetDoubleInput(&numData, h) < 0) { opserr << "WARNING LadrunoConcrete3D: -ductility wants Ah Bh Ch Dh\n"; return 0; }
      Ah = h[0]; Bh = h[1]; Ch = h[2]; Dh = h[3];
    }
    else if (strcmp(flag, "-lch") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &lch) < 0 || lch <= 0.0) { opserr << "WARNING LadrunoConcrete3D: -lch wants lch > 0\n"; return 0; }
    }
    else if (strcmp(flag, "-autoRegularization") == 0) {
      autoReg = true;
    }
    else if (strcmp(flag, "-implex") == 0) {
      implex = true;
    }
    else if (strcmp(flag, "-eta") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &eta) < 0 || eta < 0.0) { opserr << "WARNING LadrunoConcrete3D: -eta wants eta >= 0\n"; return 0; }
    }
    else if (strcmp(flag, "-ctTemper") == 0) {
      const char* mode = OPS_GetString();
      if      (strcmp(mode, "none") == 0)   ctTemper = 0;
      else if (strcmp(mode, "alphat") == 0) ctTemper = 1;
      else if (strcmp(mode, "proj") == 0)   ctTemper = 2;
      else { opserr << "WARNING LadrunoConcrete3D: -ctTemper wants {none|alphat|proj}\n"; return 0; }
    }
    else if (strcmp(flag, "-hoop") == 0) {
      // P5 confined-fiber view (§4.6): the passive transverse-steel hoop stiffness. Used ONLY by the
      // BeamFiber view (getCopy("BeamFiber"), e.g. fibers of an NDFiberSection3d); inert for the solid
      // 3D / plane views. Optional yield via a separate -hoopFy (default = no yield).
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &hoopK) < 0 || hoopK < 0.0) { opserr << "WARNING LadrunoConcrete3D: -hoop wants K >= 0\n"; return 0; }
    }
    else if (strcmp(flag, "-hoopFy") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &hoopFy) < 0 || hoopFy <= 0.0) { opserr << "WARNING LadrunoConcrete3D: -hoopFy wants fy > 0\n"; return 0; }
    }
    else {
      opserr << "WARNING LadrunoConcrete3D: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  // ---- foot-gun guards (the kernel ASSUMES valid params; reject at the boundary) ----
  if (E <= 0.0)                      { opserr << "WARNING LadrunoConcrete3D: E must be > 0\n"; return 0; }
  if (nu < 0.0 || nu >= 0.5)         { opserr << "WARNING LadrunoConcrete3D: nu must be in [0,0.5)\n"; return 0; }
  if (fc <= 0.0 || ft <= 0.0)        { opserr << "WARNING LadrunoConcrete3D: fc, ft must be > 0 (POSITIVE magnitudes)\n"; return 0; }
  if (ft >= fc)                      { opserr << "WARNING LadrunoConcrete3D: need ft < fc (m0 = 3(fc^2-ft^2)/(fc ft) e/(e+1) > 0)\n"; return 0; }
  if (Gf <= 0.0 || Gc <= 0.0)        { opserr << "WARNING LadrunoConcrete3D: Gf, Gc must be > 0\n"; return 0; }
  if (As < 1.0)                      { opserr << "WARNING LadrunoConcrete3D: As must be >= 1 (ductility amplitude)\n"; return 0; }
  if (kupfer <= 1.0)                 { opserr << "WARNING LadrunoConcrete3D: -kupfer fcc/fc must be > 1\n"; return 0; }

  if (!haveE)
    ecc = Ladruno::Concrete3D::eccentricityFromKupfer(fc, ft, kupfer);   // (0.5,1] by construction
  if (!(ecc > 0.5 && ecc <= 1.0)) {
    opserr << "WARNING LadrunoConcrete3D: eccentricity e=" << ecc << " must be in (0.5, 1]\n";
    return 0;
  }

  // The CDPM2 consistent tangent is NON-SYMMETRIC (non-associated flow, and ~2% even associated
  // from the semi-implicit theta-freeze) — a symmetric solver will under-converge or diverge.
  // The Tier-2 IMPL-EX secant (-implex) is the degraded-elastic D_dam(w~):C0 — symmetric-part SPD on
  // SINGLE-SIGN principal states (so a symmetric solver works there) but it can lose definiteness on
  // mixed-sign high-omega states (a tension crack carrying lateral compression); an unsymmetric solver
  // is the safe default in both tiers.
  if (implex)
    opserr << "LadrunoConcrete3D (tag " << tag << "): -implex (Tier-2) — the secant is the degraded-elastic "
           << "D_dam(w~):C0, symmetric-part SPD on single-sign states; an unsymmetric solver is still "
           << "the safe default (mixed-sign high-omega can lose definiteness).\n";
  else
    opserr << "LadrunoConcrete3D (tag " << tag << "): the consistent tangent is NON-SYMMETRIC; "
           << "use an unsymmetric solver (e.g. `system FullGeneral` or `system UmfPack`).\n";

  // Duvaut-Lions viscoplastic relaxation (-eta): relaxes the inviscid plastic return toward the elastic
  // trial by beta = dt/(eta+dt) (eta in TIME units; dt = ops_Dt). eta=0 (or no time increment) => inviscid,
  // byte-identical to Tier-1. v1 applies -eta in the Tier-1 path only; under -implex the implicit solve
  // stays inviscid (the -eta + -implex composition is deferred — warn).
  if (eta > 0.0 && implex)
    opserr << "LadrunoConcrete3D (tag " << tag << "): -eta with -implex — the Duvaut-Lions relaxation is "
           << "applied to the Tier-1 implicit path only; the IMPL-EX extrapolation runs INVISCID (the "
           << "combined mode is not yet validated).\n";
  else if (eta > 0.0)
    opserr << "LadrunoConcrete3D (tag " << tag << "): -eta = " << eta << " (Duvaut-Lions, Tier-1) — relaxes "
           << "with beta = dt/(eta+dt); needs a positive time increment (transient or pseudo-time), else inviscid.\n";

  // P2h compression->tension damage temper (-ctTemper): shields the tensile damage history during
  // compression so a subsequent tension reload is not pre-damaged to ~0 (literal CDPM2, the default
  // 'none', does pre-damage it). 'alphat' (w_t=1-alpha_c) keeps both monotonic backbones exact;
  // 'proj' (tensile-stress-projected plastic-strain fraction) lightly softens the tension backbone.
  if (ctTemper == 1)
    opserr << "LadrunoConcrete3D (tag " << tag << "): -ctTemper alphat — compression->tension damage "
           << "coupling tempered by w_t = 1 - alpha_c (tension restored after compression).\n";
  else if (ctTemper == 2)
    opserr << "LadrunoConcrete3D (tag " << tag << "): -ctTemper proj — compression->tension damage "
           << "coupling tempered by the tensile-stress-projected plastic-strain fraction.\n";

  if (hoopK > 0.0) {
    opserr << "LadrunoConcrete3D (tag " << tag << "): -hoop K=" << hoopK << (hoopFy < 1.0e29 ? "" : " (no yield)")
           << " — confined-fiber view (§4.6); active only through getCopy(\"BeamFiber\") (e.g. NDFiberSection3d).\n";
    if (implex || eta > 0.0)
      opserr << "LadrunoConcrete3D (tag " << tag << "): the confined-fiber (BeamFiber) view runs Tier-1 implicit "
             << "ONLY — -implex / -eta are INERT in this view (the 3D / plane views still use them).\n";
  }

  NDMaterial* mat = new LadrunoConcrete3D(tag, E, nu, fc, ft, Gf, Gc, ecc, Df, As,
                                          qh0, Hp, Ah, Bh, Ch, Dh, rho, lch, autoReg, implex, eta, ctTemper,
                                          hoopK, hoopFy);
  if (mat == 0) {
    opserr << "WARNING LadrunoConcrete3D: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoConcrete3D::LadrunoConcrete3D()
  : NDMaterial(0, ND_TAG_LadrunoConcrete3D),
    E(0.0), nu(0.0), fc(0.0), ft(0.0), Gf(0.0), Gc(0.0), ecc(0.0), m0(0.0),
    Df(0.0), As(2.0), qh0(0.3), Hp(0.5), Ah(0.08), Bh(0.003), Ch(2.0), Dh(1.0e-6),
    rho(0.0), lchFixed(1.0), autoReg(false), implex(false), eta(0.0), ctTemper(0),
    hoopK(0.0), hoopFy(1.0e30),
    dim(DIM_3D), ncomp(6), condense(false), confined(false), cEps22(0.0),
    kp_n(0.0), etmax_n(0.0), kdt1_n(0.0), kdt2_n(0.0), kdc_n(0.0), kdc1_n(0.0), kdc2_n(0.0),
    sigtmax_n(0.0), sigcmax_n(0.0),
    wt_n(0.0), wc_n(0.0), dwt_n(0.0), dwc_n(0.0), dtn_n(0.0),
    kp_t(0.0), etmax_t(0.0), kdt1_t(0.0), kdt2_t(0.0), kdc_t(0.0), kdc1_t(0.0), kdc2_t(0.0),
    sigtmax_t(0.0), sigcmax_t(0.0),
    wt_t(0.0), wc_t(0.0), dwt_t(0.0), dwc_t(0.0), dtn_t(0.0),
    omegaT(0.0), omegaC(0.0), lastStatus(0),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->setupDim();
  this->revertToStart();
}

LadrunoConcrete3D::LadrunoConcrete3D(int tag, double E_, double nu_, double fc_, double ft_,
                                     double Gf_, double Gc_, double e_, double Df_, double As_,
                                     double qh0_, double Hp_, double Ah_, double Bh_, double Ch_, double Dh_,
                                     double rho_, double lch_, bool autoReg_, bool implex_, double eta_,
                                     int ctTemper_, double hoopK_, double hoopFy_, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoConcrete3D),
    E(E_), nu(nu_), fc(fc_), ft(ft_), Gf(Gf_), Gc(Gc_), ecc(e_),
    m0(Ladruno::Concrete3D::m0Of(fc_, ft_, e_)),
    Df(Df_), As(As_), qh0(qh0_), Hp(Hp_), Ah(Ah_), Bh(Bh_), Ch(Ch_), Dh(Dh_),
    rho(rho_), lchFixed(lch_), autoReg(autoReg_), implex(implex_), eta(eta_), ctTemper(ctTemper_),
    hoopK(hoopK_), hoopFy(hoopFy_),
    dim(dimMode), ncomp(6), condense(false), confined(false), cEps22(0.0),
    kp_n(0.0), etmax_n(0.0), kdt1_n(0.0), kdt2_n(0.0), kdc_n(0.0), kdc1_n(0.0), kdc2_n(0.0),
    sigtmax_n(0.0), sigcmax_n(0.0),
    wt_n(0.0), wc_n(0.0), dwt_n(0.0), dwc_n(0.0), dtn_n(0.0),
    kp_t(0.0), etmax_t(0.0), kdt1_t(0.0), kdt2_t(0.0), kdc_t(0.0), kdc1_t(0.0), kdc2_t(0.0),
    sigtmax_t(0.0), sigcmax_t(0.0),
    wt_t(0.0), wc_t(0.0), dwt_t(0.0), dwc_t(0.0), dtn_t(0.0),
    omegaT(0.0), omegaC(0.0), lastStatus(0),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->setupDim();
  this->revertToStart();
}

LadrunoConcrete3D::~LadrunoConcrete3D() {}

// ===========================================================================
//  dimensional view setup: reduced-vector order -> full 6-comp tensor index
//  (full ordering 0:00 1:11 2:22 3:01 4:12 5:02; the out-of-plane normal is
//  always index 2, so condensation acts on row/col 2.)
// ===========================================================================
void LadrunoConcrete3D::setupDim(void)
{
  confined = false;
  switch (dim) {
    case DIM_PSTRAIN:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_AXISYM:     ncomp = 4; { int m[4]={0,1,2,3};     for(int a=0;a<4;a++) vmap[a]=m[a]; } condense=false; break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5};   for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    // P5 confined-fiber (BeamFiber): retained {00,01,02} = axial + 2 transverse shears; the lateral
    // block {11,22,12} is condensed by driveConfinedFiber vs the hoop (NOT the sigma_22=0 single-index
    // path => `condense` stays false; `confined` selects the dedicated kernel call in integrate()).
    case DIM_BEAMFIBER:  ncomp = 3; { int m[3]={0,3,5};       for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=false; confined=true; break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

// Static condensation of the out-of-plane (33, index 2) dof so sigma_22 = 0:
//   Dtan6[I][J] -= Dtan6[I][2] Dtan6[2][J] / Dtan6[2][2].
// Done in the kernel TENSOR convention; since index 2 is a NORMAL component its
// column is unscaled, and column-halving (engineering shear) commutes with this
// rank-1 update, so getTangent halves the shear columns of the condensed matrix.
void LadrunoConcrete3D::condenseTangent(void)
{
  double d22 = Dtan6[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan6[i][2]; row[i] = Dtan6[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan6[I][J] -= col[I]*row[J]/d22;
}

// ===========================================================================
//  parameter pack + the kernel call
// ===========================================================================
void LadrunoConcrete3D::integrate(bool doTangent)
{
  Params p;
  p.E = E; p.nu = nu; p.rho = rho;
  p.fc = fc; p.ft = ft;
  p.e = ecc; p.m0 = m0;
  p.Gf = Gf; p.Gc = Gc;
  p.Df = Df; p.As = As;
  p.qh0 = qh0; p.Hp = Hp; p.Ah = Ah; p.Bh = Bh; p.Ch = Ch; p.Dh = Dh;
  p.eta = eta; p.implex = implex;                 // Tier-1 (default, +Duvaut-Lions -eta) / Tier-2 IMPL-EX (-implex)
  p.ctTemper = ctTemper;                          // P2h compression->tension damage temper (none/alphat/proj)

  // crack-band: lch from the active element (mesh-objective) when -autoRegularization is on,
  // else the fixed -lch the input was calibrated for.
  double lch = lchFixed;
  if (autoReg && ops_TheActiveElement != 0) {
    double l = ops_TheActiveElement->getCharacteristicLength();
    if (l > 0.0) lch = l;
  }
  p.lch = lch; p.lch_ref = lch;

  State in;
  for (int i = 0; i < 6; i++) { in.eps[i] = eps_n[i]; in.sig[i] = sig_n[i]; in.sigEff[i] = sigEff_n[i]; }
  in.kp = kp_n; in.et_max = etmax_n;
  in.kdt1 = kdt1_n; in.kdt2 = kdt2_n;
  in.kdc = kdc_n; in.kdc1 = kdc1_n; in.kdc2 = kdc2_n;
  in.sigtMax = sigtmax_n; in.sigcMax = sigcmax_n;     // P2g monotone drive history (no-heal cyclic damage)
  // IMPL-EX committed bookkeeping (the extrapolation source); harmless when !implex
  in.wt = wt_n; in.wc = wc_n; in.dwt = dwt_n; in.dwc = dwc_n; in.dt_n = dtn_n;
  for (int i = 0; i < 6; i++) in.depl[i] = depl_n[i];

  // current time increment: the IMPL-EX extrapolation ratio r = dt/dt_n (-implex) AND the Duvaut-Lions
  // relaxation factor beta = dt/(eta+dt) (-eta) both read it (ASDConcrete3D uses the same global).
  // Ignored when neither -implex nor -eta is active (inviscid Tier-1).
  double dt = ops_Dt;

  State out;
  double sigEffImpl[6];
  if (confined) {
    // P5 confined-fiber view: the axial (retained {00,01,02}) strain is imposed in strain6; the lateral
    // block {11,22,12} is condensed by a nested Newton vs the passive hoop spring (§4.6). driveConfinedFiber
    // overwrites strain6[{1,2,4}] with the converged lateral strains + writes the NOMINAL stress + the
    // condensed tangent on the retained comps. omegaT/omegaC (recorders) = the implicit damage in out.
    lastStatus = Ladruno::Concrete3D::driveConfinedFiber(p, strain6, in, out, stress6, sigEffImpl, Dtan6,
                                                         doTangent, hoopK, hoopFy, dt);
    omegaT = out.wt; omegaC = out.wc;
  } else {
    lastStatus = Ladruno::Concrete3D::returnMap(p, strain6, in, out, stress6, sigEffImpl, Dtan6,
                                                doTangent, dt, /*hardening=*/true, &omegaT, &omegaC);
  }

  for (int i = 0; i < 6; i++) sigEff6[i] = out.sigEff[i];
  kp_t = out.kp; etmax_t = out.et_max;
  kdt1_t = out.kdt1; kdt2_t = out.kdt2;
  kdc_t = out.kdc; kdc1_t = out.kdc1; kdc2_t = out.kdc2;
  sigtmax_t = out.sigtMax; sigcmax_t = out.sigcMax;   // P2g monotone drive history
  // trial IMPL-EX bookkeeping (committed on commitState)
  wt_t = out.wt; wc_t = out.wc; dwt_t = out.dwt; dwc_t = out.dwc; dtn_t = out.dt_n;
  for (int i = 0; i < 6; i++) depl_t[i] = out.depl[i];

  if (lastStatus != 0)
    opserr << "WARNING LadrunoConcrete3D: return map did not converge (tag "
           << this->getTag() << ") — step-cut\n";
}

// ===========================================================================
//  strain interface (element passes ENGINEERING shear; kernel wants TENSOR shear)
//
//  The element vector is the reduced view; map it into the full 6-comp tensor
//  strain (the kernel is always 3D). Components NOT in vmap are zero (e.g. eps_22
//  in plane strain, the transverse shears in 2D); for condensed modes (PlaneStress
//  / PlateFiber) eps_22 is an internal unknown solved so sigma_22 = 0 — preserve
//  its running value as the Newton starting guess.
// ===========================================================================
int LadrunoConcrete3D::setTrialStrain(const Vector& e)
{
  if (confined) {
    // P5 confined-fiber (BeamFiber): retained {00,01,02} from the element; seed the condensed lateral
    // block {11,22,12} from the committed strain (warm start for the inner Newton, which lives inside
    // driveConfinedFiber). One integrate() does the whole nested lateral solve + the condensed tangent.
    strain6[1] = eps_n[1]; strain6[2] = eps_n[2]; strain6[4] = eps_n[4];
    for (int a = 0; a < ncomp; a++) {
      int full = vmap[a];
      double val = e(a);
      if (full >= 3) val *= 0.5;        // engineering -> tensor shear
      strain6[full] = val;
    }
    this->integrate(true);
    return 0;
  }

  double eps22 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps22;

  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double val = e(a);
    if (full >= 3) val *= 0.5;          // engineering -> tensor shear
    strain6[full] = val;
  }

  if (!condense) {
    this->integrate(true);
    return 0;
  }

  // enforce sigma_22 = 0: Newton on eps_22 (= strain6[2]); dSNPO sec 9.2.3. The
  // CDPM2 tangent can soften/lose definiteness post-peak (snap-back) — guard the
  // pivot and warn rather than diverge silently.
  const int maxIt = 30;
  for (int it = 0; it < maxIt; it++) {
    this->integrate(true);
    double d22 = Dtan6[2][2];
    double smag = 0.0;
    for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-9 * (smag > 1.0 ? smag : 1.0);
    if (fabs(stress6[2]) <= tol22) break;
    if (fabs(d22) < 1.0e-300) break;
    strain6[2] -= stress6[2] / d22;
    if (it == maxIt - 1)
      opserr << "WARNING LadrunoConcrete3D: sigma_22 condensation did not converge (tag "
             << this->getTag() << ", |s22|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return 0;
}

int LadrunoConcrete3D::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoConcrete3D::setTrialStrainIncr(const Vector& v)
{
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    double cur = (full >= 3) ? 2.0 * strain6[full] : strain6[full];   // tensor -> engineering
    ne(a) = cur + v(a);
  }
  return this->setTrialStrain(ne);
}

int LadrunoConcrete3D::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoConcrete3D::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];   // true tensor stress == engineering
  return stressOut;
}

double LadrunoConcrete3D::getStressZZ(void)
{
  // plane strain drops sigma_zz from getStress(); it lives in stress6[2]
  if (dim == DIM_PSTRAIN)
    return stress6[2];
  return NDMaterial::getStressZZ();   // NaN = not applicable
}

const Vector& LadrunoConcrete3D::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) {
    int full = vmap[a];
    strainOut(a) = (full >= 3) ? 2.0 * strain6[full] : strain6[full];  // tensor -> engineering shear
  }
  return strainOut;
}

const Matrix& LadrunoConcrete3D::getTangent(void)
{
  // kernel tangent is in the TENSOR convention (dsig_ij = 2G deps_ij); the element wants
  // d(sigma)/d(engineering strain), so the shear COLUMNS (deps_eng = 2 deps_tensor) are halved.
  // For condensed modes Dtan6 is already statically condensed (sigma_22 = 0) by setTrialStrain.
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) {
      int fb = vmap[b];
      tangentOut(a, b) = (fb < 3) ? Dtan6[vmap[a]][fb] : 0.5 * Dtan6[vmap[a]][fb];
    }
  return tangentOut;
}

const Matrix& LadrunoConcrete3D::getInitialTangent(void)
{
  Params p;
  p.E = E; p.nu = nu;
  double Cel[6][6];
  Ladruno::Concrete3D::elasticC(p, Cel);
  // condensed views (PlaneStress / PlateFiber) get the proper plane-stress elastic
  // modulus as K0 (the rank-1 33-condensation of the isotropic C) — a correct initial
  // stiffness matters for the softening-fragile Newton path.
  if (condense) {
    double d22 = Cel[2][2];
    if (fabs(d22) > 1.0e-300) {
      double col[6], rw[6];
      for (int i = 0; i < 6; i++) { col[i] = Cel[i][2]; rw[i] = Cel[2][i]; }
      for (int I = 0; I < 6; I++)
        for (int J = 0; J < 6; J++) Cel[I][J] -= col[I]*rw[J]/d22;
    }
  } else if (confined) {
    // P5 BeamFiber: condense the elastic C over the lateral block {11,22,12}. At eps_lat=0 the hoop is
    // slack (hoopStiffness=0), so K0 is the free-reduction BeamFiber elastic stiffness on {00,01,02}.
    const int L[3] = {1, 2, 4};
    double Kll[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Kll[i][j] = Cel[L[i]][L[j]];
    double Tmat[3][6];
    for (int k = 0; k < 6; k++) {
      double rhs[3] = { Cel[L[0]][k], Cel[L[1]][k], Cel[L[2]][k] }, tcol[3];
      if (!Ladruno::Concrete3D::solve3(Kll, rhs, tcol)) { tcol[0]=tcol[1]=tcol[2]=0.0; }
      for (int i = 0; i < 3; i++) Tmat[i][k] = tcol[i];
    }
    double Ccond[6][6];
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) {
      double v = Cel[i][j];
      for (int a = 0; a < 3; a++) v -= Cel[i][L[a]]*Tmat[a][j];
      Ccond[i][j] = v;
    }
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) Cel[i][j] = Ccond[i][j];
  }
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) {
      int fb = vmap[b];
      tangentOut(a, b) = (fb < 3) ? Cel[vmap[a]][fb] : 0.5 * Cel[vmap[a]][fb];
    }
  return tangentOut;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoConcrete3D::commitState(void)
{
  for (int i = 0; i < 6; i++) { eps_n[i] = strain6[i]; sig_n[i] = stress6[i]; sigEff_n[i] = sigEff6[i]; }
  kp_n = kp_t; etmax_n = etmax_t;
  kdt1_n = kdt1_t; kdt2_n = kdt2_t;
  kdc_n = kdc_t; kdc1_n = kdc1_t; kdc2_n = kdc2_t;
  sigtmax_n = sigtmax_t; sigcmax_n = sigcmax_t;   // P2g monotone drive history
  wt_n = wt_t; wc_n = wc_t; dwt_n = dwt_t; dwc_n = dwc_t; dtn_n = dtn_t;
  for (int i = 0; i < 6; i++) depl_n[i] = depl_t[i];
  cEps22 = strain6[2];               // converged out-of-plane strain (condensed modes)
  return 0;
}

int LadrunoConcrete3D::revertToLastCommit(void)
{
  // restore the trial buffers to the committed state (a re-issued setTrialStrain overwrites them)
  for (int i = 0; i < 6; i++) { strain6[i] = eps_n[i]; stress6[i] = sig_n[i]; sigEff6[i] = sigEff_n[i]; }
  if (condense) strain6[2] = cEps22;
  kp_t = kp_n; etmax_t = etmax_n;
  kdt1_t = kdt1_n; kdt2_t = kdt2_n;
  kdc_t = kdc_n; kdc1_t = kdc1_n; kdc2_t = kdc2_n;
  sigtmax_t = sigtmax_n; sigcmax_t = sigcmax_n;   // P2g monotone drive history
  wt_t = wt_n; wc_t = wc_n; dwt_t = dwt_n; dwc_t = dwc_n; dtn_t = dtn_n;
  for (int i = 0; i < 6; i++) depl_t[i] = depl_n[i];
  omegaT = 0.0; omegaC = 0.0;
  Params p; p.E = E; p.nu = nu;
  Ladruno::Concrete3D::elasticC(p, Dtan6);
  if (condense) this->condenseTangent();
  return 0;
}

int LadrunoConcrete3D::revertToStart(void)
{
  for (int i = 0; i < 6; i++) {
    eps_n[i] = 0.0; sig_n[i] = 0.0; sigEff_n[i] = 0.0;
    strain6[i] = 0.0; stress6[i] = 0.0; sigEff6[i] = 0.0;
    depl_n[i] = 0.0; depl_t[i] = 0.0;
  }
  kp_n = 0.0; etmax_n = 0.0; kdt1_n = 0.0; kdt2_n = 0.0; kdc_n = 0.0; kdc1_n = 0.0; kdc2_n = 0.0;
  kp_t = 0.0; etmax_t = 0.0; kdt1_t = 0.0; kdt2_t = 0.0; kdc_t = 0.0; kdc1_t = 0.0; kdc2_t = 0.0;
  sigtmax_n = sigcmax_n = sigtmax_t = sigcmax_t = 0.0;   // P2g monotone drive history
  wt_n = wc_n = dwt_n = dwc_n = dtn_n = 0.0;
  wt_t = wc_t = dwt_t = dwc_t = dtn_t = 0.0;
  omegaT = 0.0; omegaC = 0.0; lastStatus = 0;
  cEps22 = 0.0;
  Params p; p.E = E; p.nu = nu;
  Ladruno::Concrete3D::elasticC(p, Dtan6);
  if (condense) this->condenseTangent();
  return 0;
}

const char* LadrunoConcrete3D::getType(void) const
{
  switch (dim) {
    case DIM_PSTRAIN:    return "PlaneStrain";
    case DIM_AXISYM:     return "AxiSymmetric";
    case DIM_PLATEFIBER: return "PlateFiber";
    case DIM_PSTRESS:    return "PlaneStress";
    case DIM_BEAMFIBER:  return "BeamFiber";
    default:             return "ThreeDimensional";
  }
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoConcrete3D::getCopy(void)
{
  return new LadrunoConcrete3D(this->getTag(), E, nu, fc, ft, Gf, Gc, ecc, Df, As,
                               qh0, Hp, Ah, Bh, Ch, Dh, rho, lchFixed, autoReg, implex, eta, ctTemper,
                               hoopK, hoopFy, dim);
}

NDMaterial* LadrunoConcrete3D::getCopy(const char* type)
{
  int d = -1;
  if      (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)            d = DIM_3D;
  else if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)      d = DIM_PSTRAIN;
  else if (strcmp(type, "AxiSymmetric") == 0 || strcmp(type, "AxiSymmetric2D") == 0)    d = DIM_AXISYM;
  else if (strcmp(type, "PlateFiber") == 0)                                             d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0)      d = DIM_PSTRESS;
  else if (strcmp(type, "BeamFiber") == 0)                                              d = DIM_BEAMFIBER;

  if (d < 0)
    return NDMaterial::getCopy(type);   // let the base report the unsupported type

  return new LadrunoConcrete3D(this->getTag(), E, nu, fc, ft, Gf, Gc, ecc, Df, As,
                               qh0, Hp, Ah, Bh, Ch, Dh, rho, lchFixed, autoReg, implex, eta, ctTemper,
                               hoopK, hoopFy, d);
}

// ===========================================================================
//  parallel / serialization (flat Vector — the kernel state is all fixed-size scalars)
// ===========================================================================
static const int LC3D_NDATA = 1 + 18 + 1 + 1 + 25 + 2 + 1 + 1 + 1 + 1 + 11 + 2;
// tag +18 params +autoReg +dim +25 state +2 P2g(sigtmax,sigcmax) +cEps22 +implex +eta
// +1 ctTemper(P2h) +IMPL-EX committed(wt,wc,dwt,dwc,dtn + depl[6]) +2 P5 hoop(hoopK,hoopFy)

int LadrunoConcrete3D::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(LC3D_NDATA);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = E; data(c++) = nu; data(c++) = fc; data(c++) = ft; data(c++) = Gf; data(c++) = Gc;
  data(c++) = ecc; data(c++) = m0; data(c++) = Df; data(c++) = As;
  data(c++) = qh0; data(c++) = Hp; data(c++) = Ah; data(c++) = Bh; data(c++) = Ch; data(c++) = Dh;
  data(c++) = rho; data(c++) = lchFixed;
  data(c++) = autoReg ? 1.0 : 0.0;
  data(c++) = (double)dim;
  for (int i = 0; i < 6; i++) data(c++) = eps_n[i];
  for (int i = 0; i < 6; i++) data(c++) = sig_n[i];
  for (int i = 0; i < 6; i++) data(c++) = sigEff_n[i];
  data(c++) = kp_n; data(c++) = etmax_n;
  data(c++) = kdt1_n; data(c++) = kdt2_n;
  data(c++) = kdc_n; data(c++) = kdc1_n; data(c++) = kdc2_n;
  data(c++) = sigtmax_n; data(c++) = sigcmax_n;   // P2g monotone drive history
  data(c++) = cEps22;
  data(c++) = implex ? 1.0 : 0.0;
  data(c++) = eta;
  data(c++) = (double)ctTemper;                    // P2h compression->tension damage temper mode
  data(c++) = wt_n; data(c++) = wc_n; data(c++) = dwt_n; data(c++) = dwc_n; data(c++) = dtn_n;
  for (int i = 0; i < 6; i++) data(c++) = depl_n[i];
  data(c++) = hoopK; data(c++) = hoopFy;           // P5 confined-fiber hoop spring (DIM_BEAMFIBER)

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoConcrete3D::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoConcrete3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(LC3D_NDATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoConcrete3D::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  E = data(c++); nu = data(c++); fc = data(c++); ft = data(c++); Gf = data(c++); Gc = data(c++);
  ecc = data(c++); m0 = data(c++); Df = data(c++); As = data(c++);
  qh0 = data(c++); Hp = data(c++); Ah = data(c++); Bh = data(c++); Ch = data(c++); Dh = data(c++);
  rho = data(c++); lchFixed = data(c++);
  autoReg = (data(c++) != 0.0);
  dim = (int)data(c++);
  for (int i = 0; i < 6; i++) eps_n[i] = data(c++);
  for (int i = 0; i < 6; i++) sig_n[i] = data(c++);
  for (int i = 0; i < 6; i++) sigEff_n[i] = data(c++);
  kp_n = data(c++); etmax_n = data(c++);
  kdt1_n = data(c++); kdt2_n = data(c++);
  kdc_n = data(c++); kdc1_n = data(c++); kdc2_n = data(c++);
  sigtmax_n = data(c++); sigcmax_n = data(c++);   // P2g monotone drive history
  cEps22 = data(c++);
  implex = (data(c++) != 0.0);
  eta = data(c++);
  ctTemper = (int)data(c++);                        // P2h compression->tension damage temper mode
  wt_n = data(c++); wc_n = data(c++); dwt_n = data(c++); dwc_n = data(c++); dtn_n = data(c++);
  for (int i = 0; i < 6; i++) depl_n[i] = data(c++);
  hoopK = data(c++); hoopFy = data(c++);            // P5 confined-fiber hoop spring (DIM_BEAMFIBER)

  this->setupDim();             // rebuild vmap/ncomp/condense + resize the return buffers
  this->revertToLastCommit();   // sync trial buffers + tangent to the received committed state
  return 0;
}

// ===========================================================================
//  print
// ===========================================================================
void LadrunoConcrete3D::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoConcrete3D (CDPM2-grade solid concrete, 3D)" << endln;
  s << "  tag      : " << this->getTag() << endln;
  s << "  E, nu    : " << E << ", " << nu << endln;
  s << "  fc, ft   : " << fc << ", " << ft << "   (e=" << ecc << ", m0=" << m0 << ")" << endln;
  s << "  Gf, Gc   : " << Gf << ", " << Gc << "   (lch=" << lchFixed
    << (autoReg ? ", auto-regularized" : "") << ")" << endln;
  s << "  Df, As   : " << Df << ", " << As << endln;
  s << "  harden   : qh0=" << qh0 << " Hp=" << Hp
    << "   ductility Ah=" << Ah << " Bh=" << Bh << " Ch=" << Ch << " Dh=" << Dh << endln;
  s << "  damage   : omega_t=" << omegaT << " omega_c=" << omegaC << "   kappa_p=" << kp_t << endln;
  s << "  rho      : " << rho << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoConcrete3D::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "effectiveStress") == 0 || strcmp(a, "sigEff") == 0)
    return new MaterialResponse(this, 4, Vector(ncomp));
  if (strcmp(a, "damage") == 0 || strcmp(a, "omega") == 0)
    return new MaterialResponse(this, 5, Vector(2));
  if (strcmp(a, "kappaP") == 0 || strcmp(a, "kappa_p") == 0 || strcmp(a, "kappaPlastic") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 7, Vector(ncomp));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoConcrete3D::getResponse(int responseID, Information& matInfo)
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
    case 4:                                   // effective (undamaged) stress, reduced view
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int a = 0; a < ncomp; a++) v(a) = sigEff6[vmap[a]];   // stress unscaled
      }
      return 0;
    case 5:                                   // dual damage [omega_t, omega_c]
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = omegaT; v(1) = omegaC;
      }
      return 0;
    case 6:                                   // hardening variable kappa_p
      if (matInfo.theVector) (*(matInfo.theVector))(0) = kp_t;
      return 0;
    case 7:                                   // plastic strain eps - C^-1 sigEff (engineering shear, reduced view)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        double epl[6];
        Params p; p.E = E; p.nu = nu;
        Ladruno::Concrete3D::plasticStrain6(sigEff6, strain6, p, epl);
        for (int a = 0; a < ncomp; a++) {
          int full = vmap[a];
          v(a) = (full >= 3) ? 2.0 * epl[full] : epl[full];
        }
      }
      return 0;
    default:
      return -1;
  }
}
