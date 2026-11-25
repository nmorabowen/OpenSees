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

// Ladruno WP-254 (program #1131): see LadrunoMCFTStrip.h.

#include "LadrunoMCFTStrip.h"

#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <LadrunoOptSpec.h>

#include <string.h>
#include <math.h>
#include <cmath>

using namespace ladruno_mcft;

void* OPS_LadrunoMCFTStrip()
{
  static const char* OWNER = "LadrunoMCFTStrip";
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "ERROR nDMaterial LadrunoMCFTStrip: want  nDMaterial LadrunoMCFTStrip $tag $fc $fp $ep "
              "$e20 -crackSpacing $smx $smy [-ties $rho $fyt $fut] [-tieSteel $Es $esh $esu] "
              "[-agg $a] [-betaC $C] [-cvci $c] [-Gxz $G] [-tangentStep $h]\n";
    return 0;
  }
  int tag, one = 1;
  if (OPS_GetIntInput(&one, &tag) < 0) {
    opserr << "ERROR nDMaterial LadrunoMCFTStrip: invalid $tag\n";
    return 0;
  }
  double d4[4]; int four = 4;
  if (OPS_GetDoubleInput(&four, d4) < 0) {
    opserr << "ERROR nDMaterial LadrunoMCFTStrip " << tag << ": want four numbers $fc $fp $ep $e20\n";
    return 0;
  }
  Strip S;
  S.conc.fc = d4[0]; S.conc.fp = d4[1]; S.conc.ep = d4[2]; S.conc.e20 = d4[3];
  double Gxz = -1.0, hTan = 1e-7;
  bool haveSpacing = false, haveTies = false, haveTieSteel = false, haveGxz = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    const ladruno_opt::OptSpec* o = ladruno_opt::find(ladruno_opt::kLadrunoMCFTStripOptions, opt);
    if (o == 0) { ladruno_opt::reportUnknown("LadrunoMCFTStrip", tag, opt, ladruno_opt::kLadrunoMCFTStripOptions); return 0; }
    if (!ladruno_opt::haveValues(OWNER, tag, *o)) return 0;
    double v[3];
    int n = o->nargs;
    if (n > 0 && OPS_GetDoubleInput(&n, v) < 0) {
      opserr << "ERROR nDMaterial LadrunoMCFTStrip " << tag << ": " << opt << " wants " << o->nargs
             << " number(s)\n";
      return 0;
    }
    if (strcmp(opt, "-crackSpacing") == 0) {
      S.conc.smx = v[0]; S.conc.smy = v[1]; haveSpacing = true;
    }
    else if (strcmp(opt, "-ties") == 0) {
      S.rho = v[0]; S.tie.fy = v[1]; S.tie.fu = v[2]; haveTies = true;
    }
    else if (strcmp(opt, "-tieSteel") == 0) {
      S.tie.Es = v[0]; S.tie.esh = v[1]; S.tie.esu = v[2]; haveTieSteel = true;
    }
    else if (strcmp(opt, "-agg") == 0)         S.conc.agg = v[0];
    else if (strcmp(opt, "-betaC") == 0)       S.conc.kbeta = v[0];
    else if (strcmp(opt, "-cvci") == 0)        S.conc.cvci = v[0];
    else if (strcmp(opt, "-Gxz") == 0)         { Gxz = v[0]; haveGxz = true; }
    else if (strcmp(opt, "-tangentStep") == 0) hTan = v[0];
    else { ladruno_opt::reportUnhandled(OWNER, tag, opt); return 0; }
  }

  const Concrete& C = S.conc;
  const char* bad = 0;
  if (!(C.fc > 0.0))                            bad = "$fc must be > 0";
  else if (!(C.fp > 0.0) || !(C.ep > 0.0))      bad = "$fp and $ep are compression magnitudes, > 0";
  else if (!(C.e20 > C.ep))                     bad = "$e20 must exceed $ep (linear softening from the peak)";
  else if (!haveSpacing)                        bad = "-crackSpacing $smx $smy is required";
  else if (!(C.smx > 0.0) || !(C.smy > 0.0))    bad = "-crackSpacing wants two positive lengths";
  else if (!(C.agg >= 0.0))                     bad = "-agg must be >= 0";
  else if (!(C.kbeta >= 0.0))                   bad = "-betaC must be >= 0";
  else if (!(C.cvci > 0.0))                     bad = "-cvci must be > 0";
  else if (haveTies && !(S.rho >= 0.0))         bad = "-ties: $rho must be >= 0";
  else if (haveTies && S.rho > 0.0 && (!(S.tie.fy > 0.0) || !(S.tie.fu >= S.tie.fy)))
                                                bad = "-ties: want $fyt > 0 and $fut >= $fyt";
  else if (haveTieSteel && !haveTies)           bad = "-tieSteel belongs to -ties";
  else if (!(S.tie.Es > 0.0) || !(S.tie.esh > S.tie.fy / S.tie.Es) || !(S.tie.esu > S.tie.esh))
                                                bad = "-tieSteel: want $Es > 0 and fyt/Es < $esh < $esu";
  else if (haveGxz && !(Gxz >= 0.0))           bad = "-Gxz must be >= 0";
  else if (!(hTan > 0.0))                       bad = "-tangentStep must be > 0";
  if (bad != 0) {
    opserr << "ERROR nDMaterial LadrunoMCFTStrip " << tag << ": " << bad << "\n";
    return 0;
  }
  S.conc.derive();
  if (!(S.conc.fp / S.conc.ep < S.conc.Ec)) {
    opserr << "ERROR nDMaterial LadrunoMCFTStrip " << tag << ": the secant to the peak, fp/ep = "
           << S.conc.fp / S.conc.ep << ", must be below Ec = 4700 sqrt(fc) = " << S.conc.Ec
           << " (Popovics)\n";
    return 0;
  }
  if (!haveGxz) Gxz = 0.5 * S.conc.Ec;
  return new LadrunoMCFTStrip(tag, S, Gxz, hTan, 3);
}

LadrunoMCFTStrip::LadrunoMCFTStrip(int tag, const Strip& s, double G, double h, int ord)
  : NDMaterial(tag, ND_TAG_LadrunoMCFTStrip), strip(s), Gxz(G), hTan(h), order(ord),
    hC(), eyT(0.0), oT(), foundT(0), tanValid(false),
    eps(ord), epsC(ord), sig(ord), D(ord, ord), Dinit(ord, ord)
{
  evaluateTrial();
}

LadrunoMCFTStrip::LadrunoMCFTStrip()
  : NDMaterial(0, ND_TAG_LadrunoMCFTStrip), strip(), Gxz(0.0), hTan(1e-7), order(3),
    hC(), eyT(0.0), oT(), foundT(0), tanValid(false),
    eps(3), epsC(3), sig(3), D(3, 3), Dinit(3, 3)
{
}

LadrunoMCFTStrip::~LadrunoMCFTStrip()
{
}

void
LadrunoMCFTStrip::evaluateTrial(void)
{
  if (!(strip.conc.Ec > 0.0)) {   // an empty (broker) object
    sig.Zero(); oT = Out(); eyT = 0.0; foundT = 0; tanValid = false;
    return;
  }
  foundT = evaluate(strip, hC, eps(0), eps(1), eyT, oT);
  sig(0) = oT.sx;
  sig(1) = oT.txy;
  if (order == 3) sig(2) = Gxz * eps(2);
  tanValid = false;
}

int
LadrunoMCFTStrip::setTrialStrain(const Vector& v)
{
  if (v.Size() != order) {
    opserr << "LadrunoMCFTStrip " << this->getTag() << ": strain of size " << v.Size()
           << " for an order-" << order << " view\n";
    return -1;
  }
  eps = v;
  evaluateTrial();
  for (int i = 0; i < order; i++)
    if (!std::isfinite(sig(i))) {
      opserr << "LadrunoMCFTStrip " << this->getTag() << ": non-finite stress at strain " << v;
      return -1;
    }
  return 0;
}

int LadrunoMCFTStrip::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int
LadrunoMCFTStrip::setTrialStrainIncr(const Vector& dv)
{
  Vector v(eps);
  v.addVector(1.0, dv, 1.0);
  return this->setTrialStrain(v);
}

int LadrunoMCFTStrip::setTrialStrainIncr(const Vector& dv, const Vector&) { return this->setTrialStrainIncr(dv); }

const Matrix&
LadrunoMCFTStrip::getTangent(void)
{
  if (tanValid) return D;
  D.Zero();
  if (strip.conc.Ec > 0.0) {
    double K[2][2];
    tangentFD(strip, hC, eps(0), eps(1), oT, hTan, hTan, K);
    D(0, 0) = K[0][0]; D(0, 1) = K[0][1];
    D(1, 0) = K[1][0]; D(1, 1) = K[1][1];
    if (order == 3) D(2, 2) = Gxz;
  }
  tanValid = true;
  return D;
}

const Matrix&
LadrunoMCFTStrip::getInitialTangent(void)
{
  Dinit.Zero();
  Dinit(0, 0) = strip.conc.Ec;
  Dinit(1, 1) = 0.5 * strip.conc.Ec;
  if (order == 3) Dinit(2, 2) = Gxz;
  return Dinit;
}

const Vector& LadrunoMCFTStrip::getStress(void) { return sig; }
const Vector& LadrunoMCFTStrip::getStrain(void) { return eps; }

int
LadrunoMCFTStrip::commitState(void)
{
  commit(strip, hC, eyT, oT);
  epsC = eps;
  // the committed state is the trial state under the new memory: re-evaluate so a
  // tangent asked for before the next trial strain uses the new memory
  evaluateTrial();
  return 0;
}

int
LadrunoMCFTStrip::revertToLastCommit(void)
{
  eps = epsC;
  evaluateTrial();
  return 0;
}

int
LadrunoMCFTStrip::revertToStart(void)
{
  hC = Hist();
  eps.Zero(); epsC.Zero();
  evaluateTrial();
  return 0;
}

NDMaterial*
LadrunoMCFTStrip::getCopy(void)
{
  LadrunoMCFTStrip* c = new LadrunoMCFTStrip(this->getTag(), strip, Gxz, hTan, order);
  c->hC = hC;
  c->eps = eps; c->epsC = epsC;
  c->evaluateTrial();
  return c;
}

NDMaterial*
LadrunoMCFTStrip::getCopy(const char* type)
{
  if (strcmp(type, "BeamFiber") == 0 || strcmp(type, "TimoshenkoFiber") == 0)
    return new LadrunoMCFTStrip(this->getTag(), strip, Gxz, hTan, 3);
  if (strcmp(type, "BeamFiber2d") == 0 || strcmp(type, "TimoshenkoFiber2d") == 0)
    return new LadrunoMCFTStrip(this->getTag(), strip, Gxz, hTan, 2);
  opserr << "ERROR nDMaterial LadrunoMCFTStrip " << this->getTag() << ": no '" << type
         << "' view -- the MCFT strip is a beam fibre (BeamFiber, BeamFiber2d) only (refused)\n";
  return 0;
}

const char*
LadrunoMCFTStrip::getType(void) const
{
  return (order == 3) ? "BeamFiber" : "BeamFiber2d";
}

int LadrunoMCFTStrip::getOrder(void) const { return order; }

// responses: stress, strain, tangent; "mcft" = (eps_y, e1, e2, w, f_st, f1, f2, beta,
// theta, bracketed); "eps_y" / "transverseStrain"; "crackWidth"; "tieStress"
Response*
LadrunoMCFTStrip::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];
  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "mcft") == 0 || strcmp(a, "stripState") == 0)
    return new MaterialResponse(this, 4, Vector(10));
  if (strcmp(a, "eps_y") == 0 || strcmp(a, "transverseStrain") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "crackWidth") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "tieStress") == 0)
    return new MaterialResponse(this, 7, Vector(1));
  return NDMaterial::setResponse(argv, argc, s);
}

int
LadrunoMCFTStrip::getResponse(int responseID, Information& info)
{
  switch (responseID) {
    case 1: if (info.theVector) *(info.theVector) = this->getStress();  return 0;
    case 2: if (info.theVector) *(info.theVector) = this->getStrain();  return 0;
    case 3: if (info.theMatrix) *(info.theMatrix) = this->getTangent(); return 0;
    case 4: if (info.theVector) {
              Vector& v = *(info.theVector);
              v(0) = eyT; v(1) = oT.e1; v(2) = oT.e2; v(3) = oT.w; v(4) = oT.fst;
              v(5) = oT.f1; v(6) = oT.f2; v(7) = oT.beta; v(8) = oT.theta; v(9) = foundT;
            }
            return 0;
    case 5: if (info.theVector) (*(info.theVector))(0) = eyT;    return 0;
    case 6: if (info.theVector) (*(info.theVector))(0) = oT.w;   return 0;
    case 7: if (info.theVector) (*(info.theVector))(0) = oT.fst; return 0;
    default: return NDMaterial::getResponse(responseID, info);
  }
}

// wire: ID(2) = (tag, order); Vector(kNData) = parameters, committed memory, strains
static const int kNData = 28;

int
LadrunoMCFTStrip::sendSelf(int commitTag, Channel& theChannel)
{
  ID idata(2);
  idata(0) = this->getTag();
  idata(1) = order;
  int res = theChannel.sendID(this->getDbTag(), commitTag, idata);
  if (res < 0) { opserr << "LadrunoMCFTStrip::sendSelf - failed to send ID\n"; return res; }
  Vector d(kNData);
  d.Zero();
  const Concrete& C = strip.conc;
  d(0) = C.fc;  d(1) = C.fp;  d(2) = C.ep;  d(3) = C.e20; d(4) = C.smx; d(5) = C.smy;
  d(6) = C.agg; d(7) = C.kbeta; d(8) = C.cvci;
  d(9) = strip.rho; d(10) = strip.tie.fy; d(11) = strip.tie.fu; d(12) = strip.tie.Es;
  d(13) = strip.tie.esh; d(14) = strip.tie.esu;
  d(15) = Gxz; d(16) = hTan;
  d(17) = hC.ey; d(18) = hC.em1; d(19) = hC.em2; d(20) = hC.ept; d(21) = hC.used ? 1.0 : 0.0;
  for (int i = 0; i < order; i++) { d(22 + i) = epsC(i); d(25 + i) = eps(i); }
  res = theChannel.sendVector(this->getDbTag(), commitTag, d);
  if (res < 0) { opserr << "LadrunoMCFTStrip::sendSelf - failed to send data\n"; return res; }
  return 0;
}

int
LadrunoMCFTStrip::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  ID idata(2);
  int res = theChannel.recvID(this->getDbTag(), commitTag, idata);
  if (res < 0) { opserr << "LadrunoMCFTStrip::recvSelf - failed to recv ID\n"; return res; }
  this->setTag(idata(0));
  if (idata(1) != 2 && idata(1) != 3) {
    opserr << "LadrunoMCFTStrip::recvSelf - bad order " << idata(1) << "\n";
    return -1;
  }
  if (idata(1) != order) {
    order = idata(1);
    eps = Vector(order); epsC = Vector(order); sig = Vector(order);
    D = Matrix(order, order); Dinit = Matrix(order, order);
  }
  Vector d(kNData);
  res = theChannel.recvVector(this->getDbTag(), commitTag, d);
  if (res < 0) { opserr << "LadrunoMCFTStrip::recvSelf - failed to recv data\n"; return res; }
  Concrete& C = strip.conc;
  C.fc = d(0); C.fp = d(1); C.ep = d(2); C.e20 = d(3); C.smx = d(4); C.smy = d(5);
  C.agg = d(6); C.kbeta = d(7); C.cvci = d(8);
  C.derive();
  strip.rho = d(9); strip.tie.fy = d(10); strip.tie.fu = d(11); strip.tie.Es = d(12);
  strip.tie.esh = d(13); strip.tie.esu = d(14);
  Gxz = d(15); hTan = d(16);
  hC.ey = d(17); hC.em1 = d(18); hC.em2 = d(19); hC.ept = d(20); hC.used = d(21) != 0.0;
  for (int i = 0; i < order; i++) { epsC(i) = d(22 + i); eps(i) = d(25 + i); }
  evaluateTrial();
  return 0;
}

void
LadrunoMCFTStrip::Print(OPS_Stream& s, int flag)
{
  const Concrete& C = strip.conc;
  s << "LadrunoMCFTStrip, tag: " << this->getTag() << " (" << this->getType()
    << " MCFT strip, one transverse unknown)" << endln;
  s << "  fc " << C.fc << "  Ec " << C.Ec << "  fp " << C.fp << "  ep " << C.ep << "  e20 " << C.e20
    << "  smx " << C.smx << "  smy " << C.smy << "  agg " << C.agg << "  betaC " << C.kbeta
    << "  cvci " << C.cvci << endln;
  s << "  ties rho " << strip.rho << "  fyt " << strip.tie.fy << "  fut " << strip.tie.fu
    << "  Es " << strip.tie.Es << "  esh " << strip.tie.esh << "  esu " << strip.tie.esu
    << "  Gxz " << Gxz << "  tangentStep " << hTan << endln;
  s << "  committed: eps_y " << hC.ey << "  em1 " << hC.em1 << "  em2 " << hC.em2
    << "  ept " << hC.ept << "  used " << (hC.used ? 1 : 0) << endln;
}
