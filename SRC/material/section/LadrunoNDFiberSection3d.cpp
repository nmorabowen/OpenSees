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

// Ladruno WP-245 (program #1131, rcbf-3): see LadrunoNDFiberSection3d.h.

#include "LadrunoNDFiberSection3d.h"

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <LadrunoOptSpec.h>

#include <string.h>
#include <math.h>

using ladruno_shear::Spec;

void* OPS_LadrunoNDFiberSection3d()
{
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "ERROR section LadrunoNDFiber: want  section LadrunoNDFiber $tag -shearProfile "
              "uniform|parabolic ...\n";
    return 0;
  }
  int tag, one = 1;
  if (OPS_GetIntInput(&one, &tag) < 0) {
    opserr << "ERROR section LadrunoNDFiber: invalid tag\n";
    return 0;
  }
  Spec spec;
  bool compCentroid = true;
  if (LadrunoNDFiberSection3d::parseOptions(tag, spec, compCentroid) < 0) return 0;
  return new LadrunoNDFiberSection3d(tag, 30, spec, compCentroid);
}

int
LadrunoNDFiberSection3d::parseOptions(int tag, Spec& S, bool& compCentroid)
{
  static const char* OWNER = "LadrunoNDFiber";
  bool haveProfile = false, haveDepth = false, haveResultant = false, haveAlpha = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    const ladruno_opt::OptSpec* o = ladruno_opt::find(ladruno_opt::kLadrunoNDFiberOptions, opt);
    if (o == 0) { ladruno_opt::reportUnknown("LadrunoNDFiber", tag, opt, ladruno_opt::kLadrunoNDFiberOptions); return -1; }
    if (!ladruno_opt::haveValues(OWNER, tag, *o)) return -1;
    if (strcmp(opt, "-shearProfile") == 0) {
      const char* v = OPS_GetString();
      if      (strcmp(v, "uniform") == 0)   S.profile = ladruno_shear::UNIFORM;
      else if (strcmp(v, "parabolic") == 0) S.profile = ladruno_shear::PARABOLIC;
      else {
        opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearProfile '" << v
               << "' -- want uniform or parabolic\n";
        return -1;
      }
      haveProfile = true;
    }
    else if (strcmp(opt, "-shearDepth") == 0) {
      double d[2]; int two = 2;
      if (OPS_GetDoubleInput(&two, d) < 0 || !(d[0] > 0.0) || !(d[1] > 0.0)) {
        opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearDepth wants two positive "
                  "extents $hy $hz\n";
        return -1;
      }
      S.hy = d[0]; S.hz = d[1];
      haveDepth = true;
    }
    else if (strcmp(opt, "-shearResultant") == 0) {
      const char* v = OPS_GetString();
      if      (strcmp(v, "modelC") == 0)    S.resultant = ladruno_shear::MODELC;
      else if (strcmp(v, "conjugate") == 0) S.resultant = ladruno_shear::CONJUGATE;
      else {
        opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearResultant '" << v
               << "' -- want modelC or conjugate\n";
        return -1;
      }
      haveResultant = true;
    }
    else if (strcmp(opt, "-alpha") == 0) {
      int n1 = 1;
      if (OPS_GetDoubleInput(&n1, &S.alpha) < 0 || !(S.alpha > 0.0)) {
        opserr << "ERROR section LadrunoNDFiber " << tag << ": -alpha wants a positive value\n";
        return -1;
      }
      haveAlpha = true;
    }
    else if (strcmp(opt, "-noCentroid") == 0) {
      compCentroid = false;
    }
    else { ladruno_opt::reportUnhandled(OWNER, tag, opt); return -1; }
  }

  if (!haveProfile) {
    opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearProfile uniform|parabolic is "
              "required (no silent default)\n";
    return -1;
  }
  if (S.profile == ladruno_shear::PARABOLIC) {
    if (!haveDepth) {
      opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearProfile parabolic needs "
                "-shearDepth $hy $hz (the section extent along y and z)\n";
      return -1;
    }
    if (haveAlpha) {
      opserr << "ERROR section LadrunoNDFiber " << tag << ": -alpha is the uniform profile's "
                "shape factor; the parabolic profile takes none\n";
      return -1;
    }
  } else {
    if (haveDepth || haveResultant) {
      opserr << "ERROR section LadrunoNDFiber " << tag << ": -shearDepth and -shearResultant "
                "belong to -shearProfile parabolic\n";
      return -1;
    }
  }
  return 0;
}

LadrunoNDFiberSection3d::LadrunoNDFiberSection3d()
  : NDFiberSection3d(0, SEC_TAG_LadrunoNDFiberSection3d, 0, 1.0, true), spec(), geomState(0)
{
}

LadrunoNDFiberSection3d::LadrunoNDFiberSection3d(int tag, int sizeHint, const Spec& S, bool compCentroid)
  : NDFiberSection3d(tag, SEC_TAG_LadrunoNDFiberSection3d, sizeHint > 0 ? sizeHint : 30,
                     S.profile == ladruno_shear::UNIFORM ? S.alpha : 1.0, compCentroid),
    spec(S), geomState(0)
{
}

bool
LadrunoNDFiberSection3d::checkGeometry(void)
{
  if (geomState != 0) return geomState > 0;
  geomState = 1;
  if (spec.profile != ladruno_shear::PARABOLIC) return true;
  for (int i = 0; i < numFibers; i++) {
    const double y = matData[3*i] - yBar, z = matData[3*i+1] - zBar;
    if (!ladruno_shear::inside(spec, y, z)) {
      opserr << "ERROR section LadrunoNDFiber " << this->getTag() << ": fibre " << i
             << " at (y, z) = (" << y << ", " << z << ") from the centroid lies outside the "
                "parabolic band |y| <= " << 0.5*spec.hy << ", |z| <= " << 0.5*spec.hz
             << " (-shearDepth); refused\n";
      geomState = -1;
      return false;
    }
  }
  return true;
}

int
LadrunoNDFiberSection3d::setTrialSectionDeformation(const Vector &deforms)
{
  e = deforms;
  for (int i = 0; i < 6; i++) sData[i] = 0.0;
  for (int i = 0; i < 36; i++) kData[i] = 0.0;
  if (!checkGeometry()) return -1;

  double ev[6], a[3][6], b[3][6], epsv[3], sig[3], D[3][3], S6[6] = {0,0,0,0,0,0}, K[6][6];
  for (int i = 0; i < 6; i++) ev[i] = deforms(i);
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) K[i][j] = 0.0;

  static Vector eps(3);
  int res = 0;
  for (int f = 0; f < numFibers; f++) {
    const double y = matData[3*f] - yBar, z = matData[3*f+1] - zBar, A = matData[3*f+2];
    ladruno_shear::maps(spec, y, z, a, b);
    ladruno_shear::strain(a, ev, epsv);
    eps(0) = epsv[0]; eps(1) = epsv[1]; eps(2) = epsv[2];
    res += theMaterials[f]->setTrialStrain(eps);
    const Vector& st = theMaterials[f]->getStress();
    const Matrix& Dt = theMaterials[f]->getTangent();
    for (int r = 0; r < 3; r++) {
      sig[r] = st(r);
      for (int c = 0; c < 3; c++) D[r][c] = Dt(r, c);
    }
    ladruno_shear::accumulate(A, a, b, sig, D, S6, K);
  }
  // kData is column-major (Matrix storage)
  for (int i = 0; i < 6; i++) {
    sData[i] = S6[i];
    for (int j = 0; j < 6; j++) kData[j*6 + i] = K[i][j];
  }
  return res;
}

void
LadrunoNDFiberSection3d::assembleFromMaterials(void)
{
  double a[3][6], b[3][6], sig[3], D[3][3], S6[6] = {0,0,0,0,0,0}, K[6][6];
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) K[i][j] = 0.0;
  for (int f = 0; f < numFibers; f++) {
    const double y = matData[3*f] - yBar, z = matData[3*f+1] - zBar, A = matData[3*f+2];
    ladruno_shear::maps(spec, y, z, a, b);
    const Vector& st = theMaterials[f]->getStress();
    const Matrix& Dt = theMaterials[f]->getTangent();
    for (int r = 0; r < 3; r++) {
      sig[r] = st(r);
      for (int c = 0; c < 3; c++) D[r][c] = Dt(r, c);
    }
    ladruno_shear::accumulate(A, a, b, sig, D, S6, K);
  }
  for (int i = 0; i < 6; i++) {
    sData[i] = S6[i];
    for (int j = 0; j < 6; j++) kData[j*6 + i] = K[i][j];
  }
}

const Matrix&
LadrunoNDFiberSection3d::getInitialTangent(void)
{
  static double kInit[36];
  static Matrix ki(kInit, 6, 6);
  ki.Zero();
  double a[3][6], b[3][6], D[3][3], K[6][6];
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) K[i][j] = 0.0;
  for (int f = 0; f < numFibers; f++) {
    const double y = matData[3*f] - yBar, z = matData[3*f+1] - zBar, A = matData[3*f+2];
    ladruno_shear::maps(spec, y, z, a, b);
    const Matrix& Dt = theMaterials[f]->getInitialTangent();
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++) D[r][c] = Dt(r, c);
    ladruno_shear::accumulate(A, a, b, 0, D, 0, K);
  }
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) ki(i, j) = K[i][j];
  return ki;
}

int
LadrunoNDFiberSection3d::revertToLastCommit(void)
{
  int err = 0;
  for (int f = 0; f < numFibers; f++) err += theMaterials[f]->revertToLastCommit();
  assembleFromMaterials();
  return err;
}

int
LadrunoNDFiberSection3d::revertToStart(void)
{
  int err = 0;
  for (int f = 0; f < numFibers; f++) err += theMaterials[f]->revertToStart();
  assembleFromMaterials();
  return err;
}

SectionForceDeformation*
LadrunoNDFiberSection3d::getCopy(void)
{
  LadrunoNDFiberSection3d* c = new LadrunoNDFiberSection3d();
  c->setTag(this->getTag());
  c->spec = spec;
  c->geomState = geomState;
  c->numFibers = numFibers;
  c->sizeFibers = numFibers;
  if (numFibers != 0) {
    c->theMaterials = new NDMaterial *[numFibers];
    c->matData = new double [numFibers*3];
    for (int i = 0; i < numFibers; i++) {
      c->matData[i*3]   = matData[i*3];
      c->matData[i*3+1] = matData[i*3+1];
      c->matData[i*3+2] = matData[i*3+2];
      c->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber");
      if (c->theMaterials[i] == 0) {
        opserr << "LadrunoNDFiberSection3d::getCopy -- failed to get a BeamFiber copy of fibre "
               << i << "'s material\n";
        delete c;
        return 0;
      }
    }
  }
  c->e = e;
  c->QzBar = QzBar; c->QyBar = QyBar; c->Abar = Abar;
  c->yBar = yBar;   c->zBar = zBar;
  c->computeCentroid = computeCentroid;
  c->alpha = alpha;
  c->parameterID = parameterID;
  for (int i = 0; i < 6; i++)  c->sData[i] = sData[i];
  for (int i = 0; i < 36; i++) c->kData[i] = kData[i];
  c->sectionIntegr = 0;
  return c;
}

// Wire: the stock payload, then one Vector(7) with the profile. 7 so it never
// shares a size with the stock 3*numFibers fibre vector on a database channel.
int
LadrunoNDFiberSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = NDFiberSection3d::sendSelf(commitTag, theChannel);
  if (res < 0) return res;
  static Vector v(7);
  v(0) = 1.0;                 // schema
  v(1) = spec.profile;
  v(2) = spec.resultant;
  v(3) = spec.hy;
  v(4) = spec.hz;
  v(5) = spec.alpha;
  v(6) = computeCentroid ? 1.0 : 0.0;
  res = theChannel.sendVector(this->getDbTag(), commitTag, v);
  if (res < 0) opserr << "LadrunoNDFiberSection3d::sendSelf - failed to send the shear profile\n";
  return res;
}

int
LadrunoNDFiberSection3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = NDFiberSection3d::recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) return res;
  static Vector v(7);
  res = theChannel.recvVector(this->getDbTag(), commitTag, v);
  if (res < 0) {
    opserr << "LadrunoNDFiberSection3d::recvSelf - failed to recv the shear profile\n";
    return res;
  }
  if ((int)v(0) != 1) {
    opserr << "LadrunoNDFiberSection3d::recvSelf - unknown wire schema " << v(0) << "\n";
    return -1;
  }
  spec.profile   = (int)v(1);
  spec.resultant = (int)v(2);
  spec.hy = v(3); spec.hz = v(4); spec.alpha = v(5);
  alpha = (spec.profile == ladruno_shear::UNIFORM) ? spec.alpha : 1.0;
  geomState = 0;
  return 0;
}

void
LadrunoNDFiberSection3d::Print(OPS_Stream &s, int flag)
{
  s << "\nLadrunoNDFiberSection3d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid (y,z): " << yBar << ' ' << zBar << endln;
  if (spec.profile == ladruno_shear::PARABOLIC)
    s << "\tShear profile: parabolic, hy = " << spec.hy << ", hz = " << spec.hz
      << ", resultant " << (spec.resultant == ladruno_shear::MODELC ? "modelC (V = sum tau A, unsymmetric)"
                                                                     : "conjugate (V = sum tau phi A)")
      << endln;
  else
    s << "\tShear profile: uniform, alpha = " << spec.alpha << endln;
  if (flag == 1) {
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y,z) = " << matData[3*i] << ' ' << matData[3*i+1];
      s << "\nArea = " << matData[3*i+2] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

int
LadrunoNDFiberSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc >= 1 && strstr(argv[0], "alpha") != 0) {
    opserr << "ERROR section LadrunoNDFiber " << this->getTag()
           << ": 'alpha' is not a parameter of this section (refused)\n";
    return -1;
  }
  return NDFiberSection3d::setParameter(argv, argc, param);
}

const Vector&
LadrunoNDFiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  if (spec.profile == ladruno_shear::PARABOLIC) {
    static Vector zero(6);
    static bool warned = false;
    if (!warned) {
      opserr << "ERROR section LadrunoNDFiber: DDM sensitivity is not implemented for the "
                "parabolic profile (returns zero)\n";
      warned = true;
    }
    zero.Zero();
    return zero;
  }
  return NDFiberSection3d::getStressResultantSensitivity(gradIndex, conditional);
}

const Matrix&
LadrunoNDFiberSection3d::getInitialTangentSensitivity(int gradIndex)
{
  if (spec.profile == ladruno_shear::PARABOLIC) {
    static Matrix zero(6, 6);
    zero.Zero();
    return zero;
  }
  return NDFiberSection3d::getInitialTangentSensitivity(gradIndex);
}

int
LadrunoNDFiberSection3d::commitSensitivity(const Vector& dedh, int gradIndex, int numGrads)
{
  if (spec.profile == ladruno_shear::PARABOLIC) return -1;
  return NDFiberSection3d::commitSensitivity(dedh, gradIndex, numGrads);
}
