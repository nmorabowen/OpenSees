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

// Ladruno RC plastic-damage nDMaterial — see LadrunoRCConcrete.h. classTag 33015.

#include "LadrunoRCConcrete.h"
#include <elementAPI.h>
#include <LadrunoOptSpec.h>   // WP-167: fail-closed option table
#include <Channel.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <OPS_Globals.h>   // ops_TheActiveElement (Phase 3b lch latch)
#include <LadrunoMaterialStatus.h>   // Ladruno (ADR-86b/94) - LADRUNO_MATERIAL_REFUSED
#include <Element.h>       // Element::getCharacteristicLength()
#include <string.h>
#include <math.h>
#include <vector>

using namespace ladruno_rc_kernel;

// ===========================================================================
//  parser:  nDMaterial LadrunoRCConcrete $tag $E $nu
//             -Ce {..} -Cs {..} <-Cd {..}>  (compression backbone: strain/nom-stress/damage)
//             -Te {..} -Ts {..} <-Td {..}>  (tension backbone)
//             <-Kc $Kc> <-beta> <-betaFloor $f> <-lublinerReduced>
//             <-tensStiff {vc|cm}> <-tensStiffC $c> <-tensStiffAlpha $a>  (Phase 3; vc c default 200, C2)
//             <-betaC $C> <-crackedNu $nu>                                 (Ladruno C2)
//             <-addonPlane {membrane|spatial}> <-crackGate {strain|damage}>  (Ladruno WP-241)
//             <-ties $rhoY $rhoZ $fyt $Es> <-tieHardening $b>  (BeamFiber view, Ladruno WP-243)
//             <-crackCheck sum $smx $smt>                          (BeamFiber view, Ladruno WP-243)
//             <-lateralNoRoot {refuse|modelC}>                     (BeamFiber view, Ladruno WP-247)
//             <-rho $rho> <-secant | -numericalTangent>
// ===========================================================================
void* OPS_LadrunoRCConcrete(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "nDMaterial LadrunoRCConcrete error: too few args.\n"
           << "  want: tag E nu -Ce {..} -Cs {..} <-Cd {..}> -Te {..} -Ts {..} <-Td {..}>"
           << " <-Kc Kc> <-beta> <-betaFloor f> <-lublinerReduced> <-rho rho> <-secant|-numericalTangent>\n";
    return 0;
  }

  int numData = 1;
  int tag = 0;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "nDMaterial LadrunoRCConcrete error: invalid tag.\n";
    return 0;
  }
  double E = 0.0, nu = 0.0;
  if (OPS_GetDoubleInput(&numData, &E) < 0)  { opserr << "LadrunoRCConcrete: invalid E.\n";  return 0; }
  if (OPS_GetDoubleInput(&numData, &nu) < 0) { opserr << "LadrunoRCConcrete: invalid nu.\n"; return 0; }

  std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
  double Kc = 2.0 / 3.0, betaFloor = 0.1, rho = 0.0;
  bool betaOn = false, lubRed = false;
  int  tanMode = 0;
  // Phase 2a/2b: aggregate-interlock shear retention (default off)
  bool   interlockOn = false, interlockCyclic = false;
  double aggSize = 16.0, crackStrain = 0.0, crackSpacing = 0.0, lch = 0.0, betaSrMin = 0.01;
  int    shearRetMode = 0;       // 0 mcft (default) | 1 const | 2 dsfm | 3 rots
  double shearRetFactor = 0.4;   // const-mode retention mu in (0,1]
  // Phase 2b.2b: second orthogonal crack (X-cracking) + slip-driven interlock wear
  bool   xcrackOn = false;
  double degKappa = 0.5, degSlipRef = 0.01, degMin = 0.1;
  // Phase 4: IMPL-EX (implicit-explicit) robustness for cyclic softening
  bool   implexOn = false, implexCtrl = false;
  double implexAlpha = 1.0, implexErrTol = 0.05, implexTimeRedLim = 0.01;
  // Phase 3: tension stiffening (default off => baseline-identical)
  int    tensStiffMode = 0;       // 0 off | 1 vc (Bentz) | 2 cm (Collins-Mitchell)
  // Ladruno (C2): the vc default is the Vecchio-Collins 1986 MCFT constant 200
  // (f1 = ft/(1+sqrt(200 e1))); it was 500 (Collins-Mitchell 1991, = the cm law) before C2.
  // Pass -tensStiffC 500 for the pre-C2 vc curve.
  double tensStiffC = 200.0, tensStiffAlpha = 1.0;
  bool   tensStiffCGiven = false;
  // Ladruno (C2): MCFT softening coefficient + cracked-state Poisson ratio (defaults = pre-C2)
  double betaC = 170.0, crackedNu = -1.0;
  // Ladruno WP-241: add-on plane (0 membrane | 1 spatial) and crack gate (-1 = view default | 0 strain | 1 damage)
  int addonPlane = 0, crackGate = -1;
  // Ladruno WP-243: BeamFiber-view ties (default none: a free-lateral fibre)
  double tieRhoY = 0.0, tieRhoZ = 0.0, tieFy = 0.0, tieEs = 0.0, tieB = 0.0;
  bool   tiesGiven = false;
  int    crackCheck = 0; double ccSmx = 0.0, ccSmt = 0.0;   // WP-243 Model C crack check
  int    lateralNoRoot = 0;                                // WP-247: 0 refuse | 1 modelC
  // Phase 3b: crack-band (Bazant-Oh) regularization (default off => baseline-identical)
  bool   autoReg = false;
  double lchRef = 1.0;

  auto readList = [](std::vector<double>& v) {
    v.clear();
    int nd = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      double item;
      int oldn = OPS_GetNumRemainingInputArgs();
      if (OPS_GetDoubleInput(&nd, &item) < 0) {
        int newn = OPS_GetNumRemainingInputArgs();
        if (newn < oldn) OPS_ResetCurrentInputArg(-1);   // un-consume the flag/string
        break;
      }
      v.push_back(item);
    }
  };

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    // WP-167: fail closed. A token the family table does not declare is REFUSED (the
    // command fails and names it), never ignored; a declared option must have its values.
    const ladruno_opt::OptSpec* spec = ladruno_opt::find(ladruno_opt::kLadrunoRCOptions, opt);
    if (spec == 0) { ladruno_opt::reportUnknown("LadrunoRCConcrete", tag, opt, ladruno_opt::kLadrunoRCOptions); return 0; }
    if (!ladruno_opt::haveValues("LadrunoRCConcrete", tag, *spec)) return 0;
    if      (strcmp(opt, "-Ce") == 0) readList(Ce);
    else if (strcmp(opt, "-Cs") == 0) readList(Cs);
    else if (strcmp(opt, "-Cd") == 0) readList(Cd);
    else if (strcmp(opt, "-Te") == 0) readList(Te);
    else if (strcmp(opt, "-Ts") == 0) readList(Ts);
    else if (strcmp(opt, "-Td") == 0) readList(Td);
    else if (strcmp(opt, "-Kc") == 0)        { int nd = 1; if (OPS_GetDoubleInput(&nd, &Kc) < 0)        { opserr << "LadrunoRCConcrete: -Kc needs a value.\n";        return 0; } }
    else if (strcmp(opt, "-betaFloor") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &betaFloor) < 0) { opserr << "LadrunoRCConcrete: -betaFloor needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-rho") == 0)       { int nd = 1; if (OPS_GetDoubleInput(&nd, &rho) < 0)       { opserr << "LadrunoRCConcrete: -rho needs a value.\n";       return 0; } }
    else if (strcmp(opt, "-beta") == 0)              betaOn = true;
    else if (strcmp(opt, "-lublinerReduced") == 0)   lubRed = true;
    else if (strcmp(opt, "-secant") == 0)            tanMode = 1;
    else if (strcmp(opt, "-numericalTangent") == 0)  tanMode = 2;
    else if (strcmp(opt, "-interlock") == 0)         interlockOn = true;
    else if (strcmp(opt, "-cyclic") == 0)            interlockCyclic = true;
    else if (strcmp(opt, "-agg") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &aggSize) < 0)      { opserr << "LadrunoRCConcrete: -agg needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-crackStrain") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackStrain) < 0)  { opserr << "LadrunoRCConcrete: -crackStrain needs a value.\n";  return 0; } }
    else if (strcmp(opt, "-crackSpacing") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackSpacing) < 0) { opserr << "LadrunoRCConcrete: -crackSpacing needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-lch") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &lch) < 0)          { opserr << "LadrunoRCConcrete: -lch needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-betaSrMin") == 0)    { int nd = 1; if (OPS_GetDoubleInput(&nd, &betaSrMin) < 0)    { opserr << "LadrunoRCConcrete: -betaSrMin needs a value.\n";    return 0; } }
    else if (strcmp(opt, "-xcrack") == 0)            xcrackOn = true;
    else if (strcmp(opt, "-degKappa") == 0)     { int nd = 1; if (OPS_GetDoubleInput(&nd, &degKappa) < 0)     { opserr << "LadrunoRCConcrete: -degKappa needs a value.\n";     return 0; } }
    else if (strcmp(opt, "-degSlipRef") == 0)   { int nd = 1; if (OPS_GetDoubleInput(&nd, &degSlipRef) < 0)   { opserr << "LadrunoRCConcrete: -degSlipRef needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-degMin") == 0)       { int nd = 1; if (OPS_GetDoubleInput(&nd, &degMin) < 0)       { opserr << "LadrunoRCConcrete: -degMin needs a value.\n";       return 0; } }
    else if (strcmp(opt, "-implex") == 0)            implexOn = true;
    else if (strcmp(opt, "-implexAlpha") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &implexAlpha) < 0)   { opserr << "LadrunoRCConcrete: -implexAlpha needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-implexControl") == 0) {
      implexCtrl = true; int nd = 1;
      if (OPS_GetDoubleInput(&nd, &implexErrTol) < 0 || OPS_GetDoubleInput(&nd, &implexTimeRedLim) < 0) {
        opserr << "LadrunoRCConcrete: -implexControl needs $errTol $timeReductionLimit.\n"; return 0; }
    }
    // Phase 2b.2c: -shearRetention {mcft|const|dsfm|rots} crack-shear retention curve.
    else if (strcmp(opt, "-shearRetention") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "mcft")  == 0) shearRetMode = 0;
      else if (m && strcmp(m, "const") == 0) shearRetMode = 1;
      else if (m && strcmp(m, "dsfm")  == 0) shearRetMode = 2;
      else if (m && strcmp(m, "rots")  == 0) shearRetMode = 3;
      else { opserr << "LadrunoRCConcrete: -shearRetention needs {mcft|const|dsfm|rots}.\n"; return 0; }
    }
    else if (strcmp(opt, "-shearRetFactor") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &shearRetFactor) < 0) { opserr << "LadrunoRCConcrete: -shearRetFactor needs a value.\n"; return 0; } }
    // Phase 3: -tensStiff {vc|cm} tension-stiffening mode (+ -tensStiffC c, -tensStiffAlpha a).
    else if (strcmp(opt, "-tensStiff") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "vc") == 0) tensStiffMode = 1;
      else if (m && strcmp(m, "cm") == 0) tensStiffMode = 2;
      else { opserr << "LadrunoRCConcrete: -tensStiff needs {vc|cm}.\n"; return 0; }
    }
    else if (strcmp(opt, "-tensStiffC") == 0)     { int nd = 1; tensStiffCGiven = true; if (OPS_GetDoubleInput(&nd, &tensStiffC) < 0)     { opserr << "LadrunoRCConcrete: -tensStiffC needs a value.\n";     return 0; } }
    // Ladruno (C2): -betaC C  (beta = 1/(0.8 + C eps1); default 170; V&C 1986: 0.34/|eps'c|)
    else if (strcmp(opt, "-betaC") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &betaC) < 0)          { opserr << "LadrunoRCConcrete: -betaC needs a value.\n";          return 0; } }
    // Ladruno (C2): -crackedNu nu  (Poisson ratio once eps1 >= eps_cr; default off = elastic nu kept)
    else if (strcmp(opt, "-crackedNu") == 0)      { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackedNu) < 0)      { opserr << "LadrunoRCConcrete: -crackedNu needs a value.\n";      return 0; } }
    // Ladruno WP-241: -addonPlane {membrane|spatial}  (where beta / TS / cracked-nu read eps1)
    else if (strcmp(opt, "-addonPlane") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "membrane") == 0) addonPlane = 0;
      else if (m && strcmp(m, "spatial")  == 0) addonPlane = 1;
      else { opserr << "LadrunoRCConcrete: -addonPlane needs {membrane|spatial}.\n"; return 0; }
    }
    // Ladruno WP-241: -crackGate {strain|damage}  (what opens the cracked-nu latch and the TS floor)
    else if (strcmp(opt, "-crackGate") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "strain") == 0) crackGate = 0;
      else if (m && strcmp(m, "damage") == 0) crackGate = 1;
      else { opserr << "LadrunoRCConcrete: -crackGate needs {strain|damage}.\n"; return 0; }
    }
    else if (strcmp(opt, "-tensStiffAlpha") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &tensStiffAlpha) < 0) { opserr << "LadrunoRCConcrete: -tensStiffAlpha needs a value.\n"; return 0; } }
    // Phase 3b: -autoRegularization $lch_ref  (crack-band Bazant-Oh regularization).
    else if (strcmp(opt, "-autoRegularization") == 0) {
      autoReg = true; int nd = 1;
      if (OPS_GetDoubleInput(&nd, &lchRef) < 0) { opserr << "LadrunoRCConcrete: -autoRegularization needs $lch_ref.\n"; return 0; }
    }
    // Ladruno WP-243: -ties rhoY rhoZ fyt Es  (smeared ties of the BeamFiber view, per lateral direction)
    else if (strcmp(opt, "-ties") == 0) {
      double v[4]; int nd = 4;
      if (OPS_GetDoubleInput(&nd, v) < 0) { opserr << "LadrunoRCConcrete: -ties needs rhoY rhoZ fyt Es.\n"; return 0; }
      tieRhoY = v[0]; tieRhoZ = v[1]; tieFy = v[2]; tieEs = v[3]; tiesGiven = true;
    }
    // Ladruno WP-243: -tieHardening b  (bilinear kinematic hardening ratio of the ties)
    // Ladruno WP-243: -crackCheck sum smx smt  (Model C crack check, BeamFiber view; fixed arity 3,
    // so absent = off)
    else if (strcmp(opt, "-crackCheck") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if (m && strcmp(m, "sum") == 0) {
        int nd = 2; double v[2];
        if (OPS_GetDoubleInput(&nd, v) < 0) { opserr << "LadrunoRCConcrete: -crackCheck sum needs $smx $smt.\n"; return 0; }
        crackCheck = 1; ccSmx = v[0]; ccSmt = v[1];
      }
      else { opserr << "LadrunoRCConcrete: -crackCheck needs sum $smx $smt.\n"; return 0; }
    }
    // Ladruno WP-247: -lateralNoRoot {refuse|modelC}  (BeamFiber view: when the lateral balance has no
    // root, refuse the trial, or keep Model C's bisection end point)
    else if (strcmp(opt, "-lateralNoRoot") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if (m && strcmp(m, "refuse") == 0) lateralNoRoot = 0;
      else if (m && strcmp(m, "modelC") == 0) lateralNoRoot = 1;
      else { opserr << "LadrunoRCConcrete: -lateralNoRoot needs refuse or modelC.\n"; return 0; }
    }
    else if (strcmp(opt, "-tieHardening") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &tieB) < 0) { opserr << "LadrunoRCConcrete: -tieHardening needs a value.\n"; return 0; } }
    // declared in kLadrunoRCOptions but no branch above: the table and the ladder drifted
    else { ladruno_opt::reportUnhandled("LadrunoRCConcrete", tag, opt); return 0; }
  }

  // -cyclic implies -interlock; -xcrack implies -cyclic (the X-crack/wear law lives inside
  // the cyclic friction-slip block, which lives inside the interlock block).
  if (xcrackOn) interlockCyclic = true;
  // a non-default shear-retention curve is an interlock feature; const/dsfm parametrize the
  // CYCLIC slip stiffness, so they imply -interlock -cyclic. rots only needs -interlock
  // (it disables the fixed-crack shear). mcft (0) implies nothing (it is the default).
  if (shearRetMode == 1 || shearRetMode == 2) interlockCyclic = true;
  if (shearRetMode != 0) interlockOn = true;
  if (interlockCyclic && !interlockOn) interlockOn = true;

  // tension stiffening: c sits under a sqrt (sigma_ts = ft/(1+sqrt(c*e1))), so c must be
  // positive or it produces NaN once the crack gate opens (cm mode uses a fixed 500).
  if (tensStiffMode == 1 && tensStiffC <= 0.0) {
    opserr << "nDMaterial LadrunoRCConcrete error: -tensStiffC must be > 0.\n";
    return 0;
  }
  // cm mode hard-codes the Collins-Mitchell 500 coefficient, so -tensStiffC is inert there;
  // warn rather than silently ignore a user-supplied value (footgun: they likely wanted vc).
  if (tensStiffMode == 2 && tensStiffCGiven && tensStiffC != 500.0)   // Ladruno (C2): only a USER value
    opserr << "WARNING nDMaterial LadrunoRCConcrete: -tensStiffC is ignored in cm mode "
              "(cm uses the fixed Collins-Mitchell 500); use -tensStiff vc for a tunable c.\n";
  if (betaC <= 0.0) {   // Ladruno (C2)
    opserr << "nDMaterial LadrunoRCConcrete error: -betaC must be > 0.\n";
    return 0;
  }
  if (crackedNu >= 0.5) {   // Ladruno (C2); a negative value means off
    opserr << "nDMaterial LadrunoRCConcrete error: -crackedNu must be in [0, 0.5).\n";
    return 0;
  }
  // Ladruno WP-241: the fixed-crack interlock is an x-y object (2D normal, scalar slip); it has
  // no meaning on the spatial add-on plane, so the combination fails closed (implications above
  // already turned -cyclic/-xcrack/-shearRetention into -interlock).
  if (addonPlane == 1 && interlockOn) {
    opserr << "nDMaterial LadrunoRCConcrete error: -addonPlane spatial cannot be combined with "
              "-interlock/-cyclic/-xcrack/-shearRetention (the fixed crack is x-y only).\n";
    return 0;
  }
  // spatial defaults to the damage gate. An unset gate on the membrane plane stays -1: the kernel
  // reads it as strain (it tests crackGate == 1), and the BeamFiber view resolves it to damage
  // (Ladruno WP-243), so an explicit -crackGate strain is still honoured there.
  if (crackGate < 0 && addonPlane == 1) crackGate = 1;
  // Ladruno WP-243: ties
  if (tiesGiven) {
    if (tieRhoY < 0.0 || tieRhoZ < 0.0 || tieRhoY >= 1.0 || tieRhoZ >= 1.0) {
      opserr << "nDMaterial LadrunoRCConcrete error: -ties rhoY rhoZ must be in [0, 1).\n"; return 0;
    }
    if (!(tieFy > 0.0) || !(tieEs > 0.0)) {
      opserr << "nDMaterial LadrunoRCConcrete error: -ties fyt and Es must be > 0.\n"; return 0;
    }
  }
  if (crackCheck == 1 && !(ccSmx > 0.0 && ccSmt > 0.0)) {
    opserr << "nDMaterial LadrunoRCConcrete error: -crackCheck sum $smx $smt must be > 0.\n"; return 0;
  }
  if (tieB < 0.0 || tieB >= 1.0) {
    opserr << "nDMaterial LadrunoRCConcrete error: -tieHardening b must be in [0, 1).\n"; return 0;
  }
  if (autoReg && lchRef <= 0.0) {
    opserr << "nDMaterial LadrunoRCConcrete error: -autoRegularization $lch_ref must be > 0.\n";
    return 0;
  }

  if (Ce.size() < 2 || Cs.size() != Ce.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Ce and -Cs of equal length (>=2).\n";
    return 0;
  }
  if (Te.size() < 2 || Ts.size() != Te.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Te and -Ts of equal length (>=2).\n";
    return 0;
  }

  Params P;
  P.E = E; P.nu = nu; P.Kc = Kc; P.betaFloor = betaFloor;
  P.cdf = 0.0; P.eta = 0.0;
  P.betaOn = betaOn; P.lublinerTCReduced = lubRed; P.tangentMode = tanMode;
  P.interlockOn = interlockOn; P.interlockCyclic = interlockCyclic;
  P.shearRetMode = shearRetMode; P.shearRetFactor = shearRetFactor;
  P.aggSize = aggSize; P.crackStrain = crackStrain; P.crackSpacing = crackSpacing;
  P.lch = lch; P.betaSrMin = betaSrMin; P.sqrtFc = 0.0;
  P.xcrackOn = xcrackOn; P.degKappa = degKappa; P.degSlipRef = degSlipRef; P.degMin = degMin;
  P.implex = implexOn; P.implexAlpha = implexAlpha; P.implexControl = implexCtrl;
  P.implexErrTol = implexErrTol; P.implexTimeRedLim = implexTimeRedLim;
  P.tensStiffMode = tensStiffMode; P.tensStiffC = tensStiffC;
  P.tensStiffAlpha = tensStiffAlpha; P.ftPeak = 0.0;   // ftPeak set by setupParams
  P.autoReg = autoReg; P.lchRef = lchRef;
  P.betaC = betaC; P.crackedNu = crackedNu;   // Ladruno (C2)
  P.addonPlane = addonPlane; P.crackGate = crackGate;   // Ladruno WP-241
  P.tieRhoY = tieRhoY; P.tieRhoZ = tieRhoZ; P.tieFy = tieFy; P.tieEs = tieEs; P.tieB = tieB;   // WP-243
  P.crackCheck = crackCheck; P.ccSmx = ccSmx; P.ccSmt = ccSmt;                                 // WP-243
  P.lateralNoRoot = lateralNoRoot;                                                            // WP-247

  // build backbones via the faithful ASDConcrete3D HardeningLaw c-tor + adjust()
  // (elastic-consistent q). -Cd/-Td are optional -> pad with zeros to match length.
  std::vector<double> Cdf(Ce.size(), 0.0), Tdf(Te.size(), 0.0);
  for (size_t i = 0; i < Cd.size() && i < Cdf.size(); ++i) Cdf[i] = Cd[i];
  for (size_t i = 0; i < Td.size() && i < Tdf.size(); ++i) Tdf[i] = Td[i];
  buildBackbone(P.hc, E, Ce.data(), Cs.data(), Cdf.data(), (int)Ce.size());
  buildBackbone(P.ht, E, Te.data(), Ts.data(), Tdf.data(), (int)Te.size());
  setupParams(P);

  return new LadrunoRCConcrete(tag, P, rho, LadrunoRCConcrete::DIM_3D);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoRCConcrete::LadrunoRCConcrete()
  : NDMaterial(0, ND_TAG_LadrunoRCConcrete),
    rho(0.0), dim(DIM_3D), ncomp(6), condense(false), cEps33(0.0), confined(false), bfIters(0), bfWarned(0), ccBinds(0), lateralRoot(0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    regularizationDone(false), regLch(0.0), regWarned(false),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  // safe defaults until recvSelf populates P
  P.E = 1.0; P.nu = 0.0; P.Kc = 2.0/3.0; P.fcft_ratio = 5.0; P.betaFloor = 0.1;
  P.cdf = 0.0; P.eta = 0.0; P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
  P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4; P.aggSize = 16.0;
  P.crackStrain = 0.0; P.crackSpacing = 0.0; P.lch = 0.0; P.betaSrMin = 0.01; P.sqrtFc = 0.0;
  P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
  P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
  P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
  P.tensStiffMode = 0; P.tensStiffC = 500.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
  P.autoReg = false; P.lchRef = 1.0;
  P.ht.n = 0; P.hc.n = 0;
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::LadrunoRCConcrete(int tag, const Params& P_, double rho_, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoRCConcrete),
    P(P_), rho(rho_), dim(dimMode), ncomp(6), condense(false), cEps33(0.0), confined(false), bfIters(0), bfWarned(0), ccBinds(0), lateralRoot(0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    regularizationDone(false), regLch(0.0), regWarned(false),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::~LadrunoRCConcrete() {}

// ===========================================================================
//  dimensional view setup
// ===========================================================================
void LadrunoRCConcrete::setupDim(void)
{
  confined = false;
  switch (dim) {
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};     for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    // Ladruno WP-243: BeamFiber. Retained (eps_xx, gamma_xy, gamma_xz); the lateral block {1,2,4}
    // is solved against the ties by driveBeamFiber (not the single-index sigma_33 = 0 path).
    case DIM_BEAMFIBER:  ncomp = 3; { int m[3]={0,3,5};     for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=false; confined=true; break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5}; for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

void LadrunoRCConcrete::condenseTangent(void)
{
  double d22 = Dtan[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan[i][2]; row[i] = Dtan[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] -= col[I]*row[J]/d22;
}

// ===========================================================================
//  strain interface  (ENGINEERING shear — no tensor conversion)
// ===========================================================================
extern double ops_Dt;   // current (pseudo-)time increment, from OPS_Globals

void LadrunoRCConcrete::integrate(bool do_implex, double tfac)
{
  status = returnMap3D(P, strain6, histN, stress6, Dtan, histTr, true, 0, do_implex, tfac);
}

double LadrunoRCConcrete::implexTimeFactor(void) const
{
  // Extrapolation factor tf = (dt_n / dt_{n-1}) * alpha. In DYNAMIC analysis with a
  // smooth time step this is ~alpha. In STATIC analysis the "time" is the load factor,
  // whose increment is erratic (and resets across loadConst), so the raw ratio can be
  // garbage (huge/negative) and would blow up the threshold extrapolation. Guard it:
  // fall back to alpha whenever the ratio is non-finite/non-positive, and clamp the
  // acceleration to a safe band so a one-step time spike cannot detonate the damage.
  if (!P.implex || !commitDone) return P.implexAlpha;
  if (dtime_n_commit <= 0.0 || dtime_n <= 0.0) return P.implexAlpha;
  double tf = dtime_n / dtime_n_commit * P.implexAlpha;
  if (!(tf > 0.0)) return P.implexAlpha;        // NaN / non-positive guard
  const double tfMax = 2.0 * (P.implexAlpha > 0.0 ? P.implexAlpha : 1.0);
  if (tf > tfMax) tf = tfMax;
  return tf;
}

// Phase 3b: latch the characteristic length once and regularize the softening backbones so the
// dissipated energy is mesh-objective (clone of ASDConcrete3D's first-setTrialStrain latch). The
// lch is the element's getCharacteristicLength() (already EAS-aware on ASDShellQ4) unless -lch was
// given. LOUD FAILURE (no silent fallback, ADR D5) if -autoRegularization is on but no lch can be
// resolved. Latched once (regularizationDone); getCopy runs before analysis so copies regularize
// per-element. Default off => never called.
int LadrunoRCConcrete::regularizeIfNeeded(void)
{
  if (!P.autoReg || regularizationDone) return 0;
  double lch = (P.lch > 0.0) ? P.lch
             : (ops_TheActiveElement ? ops_TheActiveElement->getCharacteristicLength() : 0.0);
  if (!(lch > 0.0)) {
    // LOUD FAILURE: do NOT latch (regularizationDone stays false) — otherwise an algorithm
    // step-retry would find it "done" and silently proceed on the un-regularized backbone. We
    // re-fail every step; the message is gated to once (regWarned) to avoid spam.
    if (!regWarned) {
      opserr << "FATAL LadrunoRCConcrete (tag " << this->getTag()
             << "): -autoRegularization requested but no characteristic length is available "
                "(no active element and no -lch). Refusing to run with a non-objective softening "
                "law — supply -lch $l or use the material inside an element.\n";
      regWarned = true;
    }
    regLch = 0.0;
    return -1;
  }
  regLch = lch;
  ladruno_rc_kernel::regularize(P.ht, lch, P.lchRef, P.E);
  ladruno_rc_kernel::regularize(P.hc, lch, P.lchRef, P.E);
  regularizationDone = true;                  // latch only on SUCCESS
  return 0;
}

int LadrunoRCConcrete::setTrialStrain(const Vector& e)
{
  // ADR-94: a refused trial returns the sentinel, never a bare -1 (hosts check the sentinel only). Nothing below
  // touches the COMMITTED state (histN, cEps33): only strain6/stress6/Dtan/histTr/status, which revertToLastCommit() resets.
  if (this->regularizeIfNeeded() < 0) { status = STATUS_NO_CONVERGE; return LADRUNO_MATERIAL_REFUSED; }
  if (P.implex) {
    dtime_n = ops_Dt;
    if (!commitDone) { dtime_0 = dtime_n; dtime_n_commit = dtime_n; }
  }
  const bool   dox = P.implex;
  const double tf  = this->implexTimeFactor();

  if (confined) {
    // Ladruno WP-243: the lateral block starts from the committed lateral strain inside the driver
    for (int i = 0; i < 6; i++) strain6[i] = 0.0;
    for (int a = 0; a < ncomp; a++) strain6[vmap[a]] = e(a);
    bool ccb = false, rsc = false, nra = false;
    status = driveBeamFiber(P, strain6, histN, tieEpN, stress6, Dtan, histTr, tieEpTr, tieSig, &bfIters, dox, tf, &ccb,
                            &rsc, &nra);
    ccBinds = ccb ? 1 : 0;
    lateralRoot = nra ? 2 : (rsc ? 1 : 0);   // WP-247
    if (status == STATUS_NO_CONVERGE) {
      if (bfWarned < 5) {
        opserr << "WARNING LadrunoRCConcrete (tag " << this->getTag() << "): BeamFiber lateral balance "
                  "against the ties did not converge; the trial is refused (eps_xx=" << e(0)
               << ", gamma_xy=" << e(1) << ", gamma_xz=" << e(2) << ")\n";
        if (++bfWarned == 5) opserr << "  (further BeamFiber non-convergence messages for this point suppressed)\n";
      }
      return LADRUNO_MATERIAL_REFUSED;
    }
    return 0;
  }

  double eps33 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps33;
  for (int a = 0; a < ncomp; a++) strain6[vmap[a]] = e(a);   // engineering shear, no factor

  if (!condense) { this->integrate(dox, tf); return (status == STATUS_NO_CONVERGE) ? LADRUNO_MATERIAL_REFUSED : 0; }

  // enforce sigma_33 = 0: guarded Newton on eps_33 (= strain6[2]); dSNPO sec 9.4
  const int maxIt = 25;
  double prevAbs = 1.0e300;
  bool converged = false;
  for (int it = 0; it < maxIt; it++) {
    this->integrate(dox, tf);
    double smag = 0.0; for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-10 * (smag > 1.0 ? smag : 1.0);
    double abs33 = fabs(stress6[2]);
    if (abs33 <= tol22) { converged = true; break; }
    double d22 = Dtan[2][2];
    double step;
    if (fabs(d22) < 1.0e-14 * P.E) step = (stress6[2] > 0.0 ? 1.0 : -1.0) * 1.0e-6;  // bisection nudge
    else                          step = stress6[2] / d22;
    double lambda = (abs33 > prevAbs) ? 0.5 : 1.0;   // damp if residual grew
    strain6[2] -= lambda * step;
    prevAbs = abs33;
  }
  if (!converged) {
    status = STATUS_NO_CONVERGE;
    opserr << "WARNING LadrunoRCConcrete: sigma_33 condensation did not converge (tag "
           << this->getTag() << ", |s33|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return (status == STATUS_NO_CONVERGE) ? LADRUNO_MATERIAL_REFUSED : 0;
}

int LadrunoRCConcrete::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v)
{
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) ne(a) = strain6[vmap[a]] + v(a);
  return this->setTrialStrain(ne);
}
int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoRCConcrete::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];   // true stress
  return stressOut;
}
const Vector& LadrunoRCConcrete::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) strainOut(a) = strain6[vmap[a]];   // engineering shear
  return strainOut;
}
const Matrix& LadrunoRCConcrete::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Dtan[vmap[a]][vmap[b]];
  return tangentOut;
}
const Matrix& LadrunoRCConcrete::getInitialTangent(void)
{
  double C0[6][6];
  elasticTangent(P.E, P.nu, C0);
  if (confined) beamFiberInitialTangent(P, C0);   // Ladruno WP-243: condensed against rho Es
  // Ladruno WP-241: the condensed views (PlaneStress, PlateFiber: sigma_33 = 0) return the
  // Schur complement on index 2, as getTangent does, so -initial Newton iterates against the
  // plane-stress E/(1-nu^2) and not the plane-strain-like lambda + 2 mu of the raw C0 block.
  if (condense) {
    double d22 = C0[2][2], col[6], row[6];
    for (int i = 0; i < 6; i++) { col[i] = C0[i][2]; row[i] = C0[2][i]; }
    for (int I = 0; I < 6; I++)
      for (int J = 0; J < 6; J++) C0[I][J] -= col[I]*row[J]/d22;
  }
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = C0[vmap[a]][vmap[b]];
  return tangentOut;
}

const char* LadrunoRCConcrete::getType(void) const
{
  switch (dim) {
    case DIM_PSTRESS:    return "PlaneStress";
    case DIM_PLATEFIBER: return "PlateFiber";
    case DIM_BEAMFIBER:  return "BeamFiber";
    default:             return "ThreeDimensional";
  }
}
int LadrunoRCConcrete::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRCConcrete::commitState(void)
{
  if (P.implex) {
    // IMPL-EX commit: the trial currently holds the EXPLICIT response (extrapolated,
    // frozen damage). Re-integrate IMPLICITLY at the converged strain to advance the
    // TRUE thresholds for next step's extrapolation and measure the implex error.
    double dt_ex = histTr.dt_bar, dc_ex = histTr.dc_bar;
    this->integrate(false, 1.0);                       // implicit pass -> histTr = true state
    implexError = fmax(fabs(dt_ex - histTr.dt_bar), fabs(dc_ex - histTr.dc_bar));
    // roll n -> n-1: the previously committed thresholds become the 'old' generation
    histTr.xt_old = histN.xt; histTr.xc_old = histN.xc; histTr.eps1_old = histN.eps1;
    dtime_n_commit = dtime_n;
    // advisory error control: warn if the implex extrapolation error is large (the user
    // should reduce the time/step size). Automatic dt-reduction is deferred.
    if (P.implexControl && implexError > P.implexErrTol)
      opserr << "WARNING LadrunoRCConcrete (tag " << this->getTag() << "): IMPL-EX error "
             << implexError << " > tol " << P.implexErrTol << " — reduce the step size.\n";
  }
  histN = histTr;
  cEps33 = strain6[2];
  tieEpN[0] = tieEpTr[0]; tieEpN[1] = tieEpTr[1];   // WP-243
  commitDone = true;
  return 0;
}
int LadrunoRCConcrete::revertToLastCommit(void)
{
  histTr = histN;
  if (condense) strain6[2] = cEps33;
  tieEpTr[0] = tieEpN[0]; tieEpTr[1] = tieEpN[1];   // WP-243
  return 0;
}
int LadrunoRCConcrete::revertToStart(void)
{
  histZero(histN);
  histZero(histTr);
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  cEps33 = 0.0;
  for (int d = 0; d < 2; d++) { tieEpN[d] = tieEpTr[d] = tieSig[d] = 0.0; }   // WP-243
  bfIters = 0; ccBinds = 0; lateralRoot = 0;
  implexError = 0.0; dtime_n = dtime_n_commit = dtime_0 = 0.0; commitDone = false;
  this->resetViewTangent();
  return 0;
}

// elastic tangent of the current view: raw C0, its sigma_33 = 0 Schur complement (WP-241), or the
// BeamFiber condensation against the elastic ties (WP-243)
void LadrunoRCConcrete::resetViewTangent(void)
{
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  if (confined) beamFiberInitialTangent(P, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  if (condense) this->condenseTangent();   // Ladruno WP-241: a condensed view starts condensed
}

// ===========================================================================
//  copies
// ===========================================================================
// NOTE (Phase 3b): regularize() mutates P.ht/P.hc in place, so we must propagate the
// regularization LATCH to the copy (mirror ASDConcrete3D's copy-ctor, which copies
// regularization_done). Otherwise a copy made AFTER this instance regularized would carry the
// already-stretched backbone but a reset flag, and re-regularize (double-scale). In the normal
// host flow getCopy runs before any setTrialStrain so the latch is false here, but propagating
// it makes copy-of-a-used-instance safe too.
NDMaterial* LadrunoRCConcrete::getCopy(void)
{
  LadrunoRCConcrete* c = new LadrunoRCConcrete(this->getTag(), P, rho, dim);
  c->regularizationDone = regularizationDone; c->regLch = regLch;
  return c;
}
NDMaterial* LadrunoRCConcrete::getCopy(const char* type)
{
  int d = -1;
  if      (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) d = DIM_3D;
  else if (strcmp(type, "PlateFiber") == 0)                                  d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) d = DIM_PSTRESS;
  else if (strcmp(type, "BeamFiber") == 0)                                   d = DIM_BEAMFIBER;   // WP-243
  if (d < 0) return NDMaterial::getCopy(type);
  Params Pv = P;
  if (d == DIM_BEAMFIBER) {
    // Ladruno WP-243: the BeamFiber view selects its own add-on mode. A fibre's shear sits in
    // gamma_xy AND gamma_xz, so beta / TS / cracked-nu read eps1 from the full strain (spatial),
    // and an unset crack gate becomes damage (the strain gate opens in plain compression, which a
    // beam fibre always sees). The fixed-crack interlock is x-y only: it is switched off here.
    // (A null copy would exit() inside the stock fibre section, so the view adapts and says so.)
    Pv.addonPlane = 1;
    if (Pv.crackGate < 0) Pv.crackGate = 1;
    if (Pv.interlockOn) {
      static bool warned = false;
      if (!warned) {
        opserr << "WARNING LadrunoRCConcrete (tag " << this->getTag() << "): -interlock/-cyclic/-xcrack/"
                  "-shearRetention are x-y only and are switched off in the BeamFiber view.\n";
        warned = true;
      }
      Pv.interlockOn = false; Pv.interlockCyclic = false; Pv.xcrackOn = false;
    }
  }
  LadrunoRCConcrete* c = new LadrunoRCConcrete(this->getTag(), Pv, rho, d);
  c->regularizationDone = regularizationDone; c->regLch = regLch;
  return c;
}

// ===========================================================================
//  parallel  (serialize params + backbones + committed history)
// ===========================================================================
static const int RC_SCHEMA_VERSION = 9;    // bump when the wire layout changes (hard-checked in recvSelf); v2 = +IMPL-EX; v3 = +shearRetFactor; v4 = +tension stiffening; v5 = +crack-band regularization; v6 = +betaC,crackedNu,nuCracked (C2); v7 = +addonPlane,crackGate,tsScale (WP-241); v8 = +ties, tie plastic strains, crackCheck (WP-243); v9 = +lateralNoRoot (WP-247)
static const int RC_NSCALAR = 1 /*schemaVersion*/ + 3 /*tag,dim,rho*/
                            + 10 /*E,nu,Kc,fcft,betaFloor,cdf,eta,betaOn,lubRed,tanMode*/
                            + 9 /*interlockOn,shearRetMode,shearRetFactor,aggSize,crackStrain,crackSpacing,lch,betaSrMin,sqrtFc*/
                            + 1 /*interlockCyclic*/
                            + 4 /*xcrackOn,degKappa,degSlipRef,degMin*/
                            + 5 /*implex,implexAlpha,implexControl,implexErrTol,implexTimeRedLim*/
                            + 4 /*tensStiffMode,tensStiffC,tensStiffAlpha,ftPeak*/
                            + 4 /*autoReg,lchRef,regularizationDone,regLch*/
                            + 5 /*dtime_n,dtime_n_commit,dtime_0,commitDone,implexError*/
                            + 2 /*betaC,crackedNu (C2)*/
                            + 2 /*addonPlane,crackGate (WP-241)*/
                            + 5 /*tieRhoY,tieRhoZ,tieFy,tieEs,tieB (WP-243)*/
                            + 3 /*crackCheck,ccSmx,ccSmt (WP-243)*/
                            + 1 /*lateralNoRoot (WP-247)*/;
static const int RC_BACK = 1 + 3*MAXPTS;   // n + x[]+y[]+q[]
static const int RC_HIST = 6 + 6 + 6        // stress_eff, strain, (xt,xc,dt_bar,dc_bar,beta,eps1)
                         + 5                // + (cracked,crackC,crackS,wmax,betaSr)
                         + 2                // + Phase-2b (tauCr,gammaCr)
                         + 2                // + Phase-2b.2b (cracked2,slipCum)
                         + 3                // + Phase-4 IMPL-EX (xt_old,xc_old,eps1_old)
                         + 1                // + C2 cracked-nu latch (nuCracked)
                         + 1;               // + WP-241 damage-gated TS onset scale (tsScale)
static const int RC_DATA = RC_NSCALAR + 2*RC_BACK + RC_HIST + 1 /*cEps33*/ + 2 /*tieEpN (WP-243)*/;

int LadrunoRCConcrete::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(RC_DATA);
  int c = 0;
  data(c++) = RC_SCHEMA_VERSION;
  data(c++) = this->getTag();
  data(c++) = dim;
  data(c++) = rho;
  data(c++) = P.E; data(c++) = P.nu; data(c++) = P.Kc; data(c++) = P.fcft_ratio;
  data(c++) = P.betaFloor; data(c++) = P.cdf; data(c++) = P.eta;
  data(c++) = P.betaOn ? 1.0 : 0.0;
  data(c++) = P.lublinerTCReduced ? 1.0 : 0.0;
  data(c++) = P.tangentMode;
  data(c++) = P.interlockOn ? 1.0 : 0.0;
  data(c++) = P.interlockCyclic ? 1.0 : 0.0;
  data(c++) = P.shearRetMode; data(c++) = P.shearRetFactor;
  data(c++) = P.aggSize; data(c++) = P.crackStrain; data(c++) = P.crackSpacing;
  data(c++) = P.lch; data(c++) = P.betaSrMin; data(c++) = P.sqrtFc;
  data(c++) = P.xcrackOn ? 1.0 : 0.0;
  data(c++) = P.degKappa; data(c++) = P.degSlipRef; data(c++) = P.degMin;
  data(c++) = P.implex ? 1.0 : 0.0;
  data(c++) = P.implexAlpha;
  data(c++) = P.implexControl ? 1.0 : 0.0;
  data(c++) = P.implexErrTol; data(c++) = P.implexTimeRedLim;
  data(c++) = P.tensStiffMode; data(c++) = P.tensStiffC;
  data(c++) = P.tensStiffAlpha; data(c++) = P.ftPeak;
  data(c++) = P.autoReg ? 1.0 : 0.0; data(c++) = P.lchRef;
  data(c++) = regularizationDone ? 1.0 : 0.0; data(c++) = regLch;
  data(c++) = dtime_n; data(c++) = dtime_n_commit; data(c++) = dtime_0;
  data(c++) = commitDone ? 1.0 : 0.0; data(c++) = implexError;
  data(c++) = P.betaC; data(c++) = P.crackedNu;   // v6 (C2)
  data(c++) = P.addonPlane; data(c++) = P.crackGate;   // v7 (WP-241)
  data(c++) = P.tieRhoY; data(c++) = P.tieRhoZ; data(c++) = P.tieFy; data(c++) = P.tieEs; data(c++) = P.tieB;   // v8 (WP-243)
  data(c++) = P.crackCheck; data(c++) = P.ccSmx; data(c++) = P.ccSmt;   // v8 (WP-243)
  data(c++) = P.lateralNoRoot;                                          // v9 (WP-247)
  data(c++) = P.ht.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.q[i];
  data(c++) = P.hc.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.q[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.stress_eff[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.strain[i];
  data(c++) = histN.xt; data(c++) = histN.xc;
  data(c++) = histN.dt_bar; data(c++) = histN.dc_bar;
  data(c++) = histN.beta; data(c++) = histN.eps1;
  data(c++) = histN.cracked; data(c++) = histN.crackC; data(c++) = histN.crackS;
  data(c++) = histN.wmax; data(c++) = histN.betaSr;
  data(c++) = histN.tauCr; data(c++) = histN.gammaCr;
  data(c++) = histN.cracked2; data(c++) = histN.slipCum;
  data(c++) = histN.xt_old; data(c++) = histN.xc_old; data(c++) = histN.eps1_old;
  data(c++) = histN.nuCracked;   // v6 (C2)
  data(c++) = histN.tsScale;     // v7 (WP-241)
  data(c++) = cEps33;
  data(c++) = tieEpN[0]; data(c++) = tieEpN[1];   // v8 (WP-243)

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoRCConcrete::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(RC_DATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  int ver = (int)data(c++);
  if (ver != RC_SCHEMA_VERSION) {
    opserr << "LadrunoRCConcrete::recvSelf - unsupported schema version " << ver
           << " (expected " << RC_SCHEMA_VERSION << ")\n";
    return -1;
  }
  this->setTag((int)data(c++));
  dim = (int)data(c++);
  rho = data(c++);
  P.E = data(c++); P.nu = data(c++); P.Kc = data(c++); P.fcft_ratio = data(c++);
  P.betaFloor = data(c++); P.cdf = data(c++); P.eta = data(c++);
  P.betaOn = (data(c++) != 0.0);
  P.lublinerTCReduced = (data(c++) != 0.0);
  P.tangentMode = (int)data(c++);
  P.interlockOn = (data(c++) != 0.0);
  P.interlockCyclic = (data(c++) != 0.0);
  P.shearRetMode = (int)data(c++); P.shearRetFactor = data(c++);
  P.aggSize = data(c++); P.crackStrain = data(c++); P.crackSpacing = data(c++);
  P.lch = data(c++); P.betaSrMin = data(c++); P.sqrtFc = data(c++);
  P.xcrackOn = (data(c++) != 0.0);
  P.degKappa = data(c++); P.degSlipRef = data(c++); P.degMin = data(c++);
  P.implex = (data(c++) != 0.0);
  P.implexAlpha = data(c++);
  P.implexControl = (data(c++) != 0.0);
  P.implexErrTol = data(c++); P.implexTimeRedLim = data(c++);
  P.tensStiffMode = (int)data(c++); P.tensStiffC = data(c++);
  P.tensStiffAlpha = data(c++); P.ftPeak = data(c++);
  P.autoReg = (data(c++) != 0.0); P.lchRef = data(c++);
  regularizationDone = (data(c++) != 0.0); regLch = data(c++);
  dtime_n = data(c++); dtime_n_commit = data(c++); dtime_0 = data(c++);
  commitDone = (data(c++) != 0.0); implexError = data(c++);
  P.betaC = data(c++); P.crackedNu = data(c++);   // v6 (C2)
  P.addonPlane = (int)data(c++); P.crackGate = (int)data(c++);   // v7 (WP-241)
  P.tieRhoY = data(c++); P.tieRhoZ = data(c++); P.tieFy = data(c++); P.tieEs = data(c++); P.tieB = data(c++);   // v8 (WP-243)
  P.crackCheck = (int)data(c++); P.ccSmx = data(c++); P.ccSmt = data(c++);   // v8 (WP-243)
  P.lateralNoRoot = (int)data(c++);                                          // v9 (WP-247)
  P.ht.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.q[i] = data(c++);
  P.hc.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.q[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.stress_eff[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.strain[i] = data(c++);
  histN.xt = data(c++); histN.xc = data(c++);
  histN.dt_bar = data(c++); histN.dc_bar = data(c++);
  histN.beta = data(c++); histN.eps1 = data(c++);
  histN.cracked = data(c++); histN.crackC = data(c++); histN.crackS = data(c++);
  histN.wmax = data(c++); histN.betaSr = data(c++);
  histN.tauCr = data(c++); histN.gammaCr = data(c++);
  histN.cracked2 = data(c++); histN.slipCum = data(c++);
  histN.xt_old = data(c++); histN.xc_old = data(c++); histN.eps1_old = data(c++);
  histN.nuCracked = data(c++);   // v6 (C2)
  histN.tsScale = data(c++);     // v7 (WP-241)
  cEps33 = data(c++);
  tieEpN[0] = data(c++); tieEpN[1] = data(c++);   // v8 (WP-243)

  this->setupDim();
  histTr = histN;
  tieEpTr[0] = tieEpN[0]; tieEpTr[1] = tieEpN[1]; tieSig[0] = tieSig[1] = 0.0;
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  if (condense) strain6[2] = cEps33;
  this->resetViewTangent();
  return 0;
}

void LadrunoRCConcrete::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoRCConcrete (RC plastic-damage + MCFT compression softening)" << endln;
  s << "  tag   : " << this->getTag() << endln;
  s << "  E, nu : " << P.E << ", " << P.nu;
  if (P.crackedNu >= 0.0) s << "   crackedNu=" << P.crackedNu << (histN.nuCracked >= 0.5 ? " (cracked)" : "");
  s << endln;
  s << "  Kc    : " << P.Kc << "   fc/ft ratio: " << P.fcft_ratio << endln;
  s << "  beta  : " << (P.betaOn ? "ON" : "off") << "  floor=" << P.betaFloor << "  C=" << P.betaC
    << "  lublinerReduced=" << (P.lublinerTCReduced ? 1 : 0) << endln;
  s << "  interlock: " << (P.interlockOn ? "ON" : "off")
    << (P.interlockCyclic ? " (cyclic)" : "") << (P.xcrackOn ? " (xcrack)" : "")
    << "  a_g=" << P.aggSize
    << "  eps_cr=" << P.crackStrain << "  s_theta=" << P.crackSpacing
    << "  betaSrMin=" << P.betaSrMin;
  if (P.xcrackOn) s << "  degKappa=" << P.degKappa << " degSlipRef=" << P.degSlipRef << " degMin=" << P.degMin;
  if (P.interlockOn) {
    const char* srn = (P.shearRetMode == 1) ? "const" : (P.shearRetMode == 2) ? "dsfm"
                    : (P.shearRetMode == 3) ? "rots" : "mcft";
    s << "  shearRetention=" << srn;
    if (P.shearRetMode == 1) s << "(mu=" << P.shearRetFactor << ")";
  }
  s << endln;
  s << "  implex: " << (P.implex ? "ON" : "off");
  if (P.implex) s << "  alpha=" << P.implexAlpha
                  << (P.implexControl ? " (control)" : "") << "  lastError=" << implexError;
  s << endln;
  s << "  tensStiff: " << (P.tensStiffMode == 1 ? "vc" : P.tensStiffMode == 2 ? "cm" : "off");
  if (P.tensStiffMode == 1) s << "  c=" << P.tensStiffC;
  if (P.tensStiffMode == 2) s << "  alpha=" << P.tensStiffAlpha;
  if (P.tensStiffMode != 0) s << "  ft=" << P.ftPeak;
  s << endln;
  s << "  autoRegularization: " << (P.autoReg ? "ON" : "off");
  if (P.autoReg) s << "  lch_ref=" << P.lchRef << (regularizationDone ? "  lch=" : "  (pending)");
  if (P.autoReg && regularizationDone) s << regLch;
  s << endln;
  s << "  view  : " << this->getType() << " (order " << ncomp << ")" << endln;
  if (confined || P.tieRhoY > 0.0 || P.tieRhoZ > 0.0)   // WP-243
    s << "  ties  : rhoY=" << P.tieRhoY << "  rhoZ=" << P.tieRhoZ << "  fyt=" << P.tieFy
      << "  Es=" << P.tieEs << "  b=" << P.tieB << (confined ? "" : "  (BeamFiber view only)") << endln;
  if (P.crackCheck == 1)                                 // WP-243
    s << "  crackCheck: sum  smx=" << P.ccSmx << "  smt=" << P.ccSmt << "  a_g=" << P.aggSize
      << (confined ? "" : "  (BeamFiber view only)") << endln;
  if (P.lateralNoRoot == 1)                              // WP-247
    s << "  lateralNoRoot: modelC" << (confined ? "" : "  (BeamFiber view only)") << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoRCConcrete::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];
  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "damage") == 0 || strcmp(a, "Damage") == 0)
    return new MaterialResponse(this, 4, Vector(2));   // (dt_bar, dc_bar)
  if (strcmp(a, "beta") == 0 || strcmp(a, "compressionSoftening") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "eps1") == 0 || strcmp(a, "transverseTensileStrain") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "betaShear") == 0 || strcmp(a, "shearRetention") == 0)
    return new MaterialResponse(this, 7, Vector(1));
  if (strcmp(a, "crackState") == 0 || strcmp(a, "crackNormal") == 0)
    return new MaterialResponse(this, 8, Vector(4));   // (cracked, crackC, crackS, wmax)
  if (strcmp(a, "crackShear") == 0)
    return new MaterialResponse(this, 9, Vector(2));   // (tauCr, gammaCr) Phase-2b
  if (strcmp(a, "xcrackState") == 0)
    return new MaterialResponse(this, 10, Vector(2));  // (cracked2, slipCum) Phase-2b.2b
  if (strcmp(a, "implexError") == 0 || strcmp(a, "ImplexError") == 0)
    return new MaterialResponse(this, 11, Vector(1));  // IMPL-EX |dt_ex - dt_im| Phase-4
  if (strcmp(a, "nuCracked") == 0 || strcmp(a, "crackedNu") == 0)
    return new MaterialResponse(this, 12, Vector(1));  // C2 cracked-nu latch (0/1)
  // Ladruno WP-243 (BeamFiber view; zeros in the other views)
  if (strcmp(a, "lateralStrain") == 0)
    return new MaterialResponse(this, 13, Vector(3));  // (eps_yy, eps_zz, gamma_yz)
  if (strcmp(a, "tieStress") == 0)
    return new MaterialResponse(this, 14, Vector(2));  // bar stress (y, z)
  if (strcmp(a, "lateralIterations") == 0)
    return new MaterialResponse(this, 15, Vector(1));
  if (strcmp(a, "crackCheck") == 0)
    return new MaterialResponse(this, 16, Vector(1));  // 1 = the Model C cap binds
  if (strcmp(a, "lateralRoot") == 0)                 // WP-247
    return new MaterialResponse(this, 17, Vector(1));  // 0 Newton | 1 root scan | 2 -lateralNoRoot modelC
  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoRCConcrete::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1: if (matInfo.theVector) *(matInfo.theVector) = this->getStress();  return 0;
    case 2: if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();  return 0;
    case 3: if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent(); return 0;
    case 4: if (matInfo.theVector) { Vector& v = *(matInfo.theVector); v(0) = histTr.dt_bar; v(1) = histTr.dc_bar; } return 0;
    case 5: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.beta;   return 0;
    case 6: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.eps1;   return 0;
    case 7: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.betaSr; return 0;
    case 8: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.cracked; v(1) = histTr.crackC; v(2) = histTr.crackS; v(3) = histTr.wmax; } return 0;
    case 9: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.tauCr; v(1) = histTr.gammaCr; } return 0;
    case 10: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.cracked2; v(1) = histTr.slipCum; } return 0;
    case 11: if (matInfo.theVector) (*(matInfo.theVector))(0) = implexError; return 0;
    case 12: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.nuCracked; return 0;
    case 13: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = confined ? strain6[1] : 0.0; v(1) = confined ? strain6[2] : 0.0; v(2) = confined ? strain6[4] : 0.0; } return 0;
    case 14: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = confined ? tieSig[0] : 0.0; v(1) = confined ? tieSig[1] : 0.0; } return 0;
    case 15: if (matInfo.theVector) (*(matInfo.theVector))(0) = bfIters; return 0;
    case 16: if (matInfo.theVector) (*(matInfo.theVector))(0) = ccBinds; return 0;
    case 17: if (matInfo.theVector) (*(matInfo.theVector))(0) = lateralRoot; return 0;
    default: return -1;
  }
}
