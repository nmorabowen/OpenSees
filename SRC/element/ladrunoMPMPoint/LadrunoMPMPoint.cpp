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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 11/2026
//
// LadrunoMPMPoint — see the header (WP-232f, mode A; interface
// tests/mpm_implicit_interface.md §1).

#include "LadrunoMPMPoint.h"
#include "../../domain/pattern/ladrunoMPM/LadrunoMPMCore.h"
#include "../ladrunoUP/LadrunoUPKernel.h"      // WP-232l: the u-p blocks at a particle (P3)

#include <Domain.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementalLoad.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <Ladruno_mutation.h>
#include <LadrunoElementShell.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <string.h>

namespace {

// Shared scratch, keyed by the number of DOFs (at most 64 nodes x 3: WP-232i's
// cubic B-spline; 27 x 3 before). The pool is static and allocates a size only
// when one is first used, so no element stores more (232i F2). The element is
// serial (no threading opt-in), so one scratch per size is safe.
const int MAX_DOF = 256;            // WP-232l: 64 nodes x (3 + p)

Matrix &scratchK(int n)
{
  static Matrix *pool[MAX_DOF + 1] = {0};
  if (pool[n] == 0) pool[n] = new Matrix(n, n);
  return *pool[n];
}

Matrix &scratchM(int n)
{
  static Matrix *pool[MAX_DOF + 1] = {0};
  if (pool[n] == 0) pool[n] = new Matrix(n, n);
  return *pool[n];
}

// WP-232l: the two-phase layout's full-size returns (node-major [u, p]), apart from the
// compact u-u scratch they are assembled from
Matrix &scratchKF(int n)
{
  static Matrix *pool[MAX_DOF + 1] = {0};
  if (pool[n] == 0) pool[n] = new Matrix(n, n);
  return *pool[n];
}

Matrix &scratchMF(int n)
{
  static Matrix *pool[MAX_DOF + 1] = {0};
  if (pool[n] == 0) pool[n] = new Matrix(n, n);
  return *pool[n];
}

Matrix &scratchC(int n)
{
  static Matrix *pool[MAX_DOF + 1] = {0};
  if (pool[n] == 0) pool[n] = new Matrix(n, n);
  return *pool[n];
}

// σ (Voigt, tensor shear) -> 3x3 row-major tensor
void voigtToTensor(int nd, const double *s, double t[9])
{
  for (int k = 0; k < 9; k++) t[k] = 0.0;
  if (nd == 3) {
    t[0] = s[0]; t[4] = s[1]; t[8] = s[2];
    t[1] = t[3] = s[3]; t[5] = t[7] = s[4]; t[2] = t[6] = s[5];
  } else {
    t[0] = s[0]; t[4] = s[1]; t[1] = t[3] = s[2];
  }
}

}  // namespace

bool LadrunoMPMPoint::diagnose = false;

LadrunoMPMPoint::LadrunoMPMPoint(int tag, int particle, int hostTag, int nd, NDMaterial *ownedMat)
  : Element(tag, ELE_TAG_LadrunoMPMPoint), p(particle), hostPatternTag(hostTag), ndm(nd),
    nv(nd == 3 ? 6 : 3), nen(0), mat(ownedMat), st(ARMED), needUpdate(true), trialValid(false),
    alMode(0), nCorner(0), alSuppressed(false), cellH(0.0), Jc(1.0), alRel(0.0),
    connect(0), J1(1.0), vol(0.0), mp(0.0), payloadDbTag(0), P(0), Pinc(0), Ki(0),
    ndfN(nd), tpOn(false), smallK(false), stabA(0.0), tpAlpha(1.0), tpRhoF(0.0), affOn(false),
    piStart_(0.0)
{
  for (int k = 0; k < 12; k++) affC[k] = 0.0;
  tpK[0] = tpK[1] = tpK[2] = 0.0;
  for (int k = 0; k < 8; k++) corner[k] = -1;
  for (int k = 0; k < 9; k++) Fci[k] = (k % 4 == 0) ? 1.0 : 0.0;
  for (int k = 0; k < 6; k++) depsT[k] = 0.0;
  for (int k = 0; k < 24; k++) gHat[k] = 0.0;
  diagTan = diagJr = 1.0;
  for (int k = 0; k < 9; k++) Frel[k] = F1[k] = R1[k] = (k % 4 == 0) ? 1.0 : 0.0;
  for (int k = 0; k < 6; k++) epsT[k] = sigHatT[k] = sigT[k] = 0.0;
  bf[0] = bf[1] = bf[2] = 0.0;
}

LadrunoMPMPoint::LadrunoMPMPoint()
  : Element(0, ELE_TAG_LadrunoMPMPoint), p(-1), hostPatternTag(0), ndm(2), nv(3), nen(0),
    mat(0), st(RESTORED), needUpdate(true), trialValid(false), alMode(0), nCorner(0),
    alSuppressed(false), cellH(0.0), Jc(1.0), alRel(0.0), connect(0), J1(1.0), vol(0.0),
    mp(0.0), payloadDbTag(0), P(0), Pinc(0), Ki(0), ndfN(2), tpOn(false), smallK(false), stabA(0.0),
    tpAlpha(1.0), tpRhoF(0.0), affOn(false), piStart_(0.0)
{
  for (int k = 0; k < 12; k++) affC[k] = 0.0;
  tpK[0] = tpK[1] = tpK[2] = 0.0;
  for (int k = 0; k < 8; k++) corner[k] = -1;
  for (int k = 0; k < 9; k++) Fci[k] = (k % 4 == 0) ? 1.0 : 0.0;
  for (int k = 0; k < 6; k++) depsT[k] = 0.0;
  for (int k = 0; k < 24; k++) gHat[k] = 0.0;
  diagTan = diagJr = 1.0;
  for (int k = 0; k < 9; k++) Frel[k] = F1[k] = R1[k] = (k % 4 == 0) ? 1.0 : 0.0;
  for (int k = 0; k < 6; k++) epsT[k] = sigHatT[k] = sigT[k] = 0.0;
  bf[0] = bf[1] = bf[2] = 0.0;
}

LadrunoMPMPoint::~LadrunoMPMPoint()
{
  LadrunoMPMPointHost *h = host();
  if (h != 0) h->mpmElementGone(p);
  delete mat;
  delete P;
  delete Pinc;
  delete Ki;
}

NDMaterial *
LadrunoMPMPoint::releaseMaterial(void)
{
  NDMaterial *m = mat;
  mat = 0;
  return m;
}

void
LadrunoMPMPoint::setRayleigh(double am, double bk)
{
  this->Element::setRayleighDampingFactors(am, bk, 0.0, 0.0);   // zero clears them (item 9)
}

void
LadrunoMPMPoint::resize(void)
{
  const int n = nen * ndm, nf = nen * ndfN;     // WP-232l: the residuals carry every node DOF
  if (P == 0 || P->Size() != nf) { delete P; P = new Vector(nf); }
  if (Pinc == 0 || Pinc->Size() != nf) { delete Pinc; Pinc = new Vector(nf); }
  g.assign(n, 0.0);
  fInt.assign(n, 0.0);
  uT.assign(n, 0.0);
  if ((int)ghostK.size() != nen) ghostK.assign(nen, 0.0);
  if ((int)S.size() != nen) S.assign(nen, 0.0);
  if ((int)dS.size() != n) dS.assign(n, 0.0);
  if ((int)vHat.size() != n) vHat.assign(n, 0.0);
  if ((int)aHat.size() != n) aHat.assign(n, 0.0);
}

void
LadrunoMPMPoint::refetchNodes(void)
{
  Domain *dom = this->getDomain();
  nodes.assign(nen, (Node *)0);
  if (dom == 0) return;
  for (int k = 0; k < nen; k++) {
    nodes[k] = dom->getNode(connect(k));
    if (nodes[k] == 0)
      opserr << "LadrunoMPMPoint " << this->getTag() << " -- grid node " << connect(k)
             << " is not in the Domain (MPM host bug)\n";
  }
}

void
LadrunoMPMPoint::setDomain(Domain *theDomain)
{
  this->DomainComponent::setDomain(theDomain);
  refetchNodes();
}

bool
LadrunoMPMPoint::connectivityEquals(const ID &tags) const
{
  if (tags.Size() != nen) return false;
  for (int k = 0; k < nen; k++)
    if (tags(k) != connect(k)) return false;
  return true;
}

int
LadrunoMPMPoint::repoint(const ID &tags, const std::vector<double> &Sk,
                         const std::vector<double> &dSk, const std::vector<char> &inactive,
                         const std::vector<double> &vH, const std::vector<double> &aH,
                         const std::vector<double> &Sck, const std::vector<double> &dSck,
                         const int *corners, int nCorners, double cellh)
{
  nen = tags.Size();
  if (nen < 1 || nen * ndm > MAX_DOF) {
    opserr << "LadrunoMPMPoint " << this->getTag() << " -- " << nen
           << " support nodes is out of range (MPM host bug)\n";
    return -1;
  }
  connect = tags;
  S = Sk;
  dS = dSk;
  if (!Sck.empty() && ((int)Sck.size() != nen * ndm || (int)dSck.size() != nen * ndm * ndm)) {
    opserr << "LadrunoMPMPoint " << this->getTag() << " -- per-component support of the wrong"
           << " size (MPM host bug)\n";
    return -1;
  }
  Sc = Sck;              // WP-232i B8: empty = the scalar path
  dSc = dSck;
  if (Sc.empty()) { std::vector<double>().swap(gc); }
  else gc.assign(nen * ndm * ndm, 0.0);
  vHat = vH;
  aHat = aH;
  // WP-232j (J11): the corner slots of the particle's cell, or none (off). A slot
  // < 0 anywhere = the fallback (J5), kept as corner[0] = -1.
  alMode = 0;
  nCorner = 0;
  for (int k = 0; k < 8; k++) corner[k] = -1;
  if (corners != 0 && nCorners == (1 << ndm) && cellh > 0.0) {
    alMode = 1;
    nCorner = nCorners;
    cellH = cellh;
    bool all = true;
    for (int k = 0; k < nCorners; k++) {
      if (corners[k] < 0 || corners[k] >= nen) all = false;
      corner[k] = corners[k];
    }
    if (!all) corner[0] = -1;
  }
  ghostK.assign(nen, 0.0);
  // WP-232l: a repoint is single-phase until setTwoPhase says otherwise
  ndfN = ndm;
  tpOn = false;
  stabA = 0.0;
  affOn = false;
  piStart_ = 0.0;
  std::vector<double>().swap(Sp); std::vector<double>().swap(dSp); std::vector<double>().swap(pRef);
  std::vector<double>().swap(Qb); std::vector<double>().swap(Hb); std::vector<double>().swap(Sb);
  std::vector<double>().swap(Htb); std::vector<double>().swap(Gs);
  LadrunoMPMPointHost *h = host();
  smallK = false;
  if (h != 0) {
    for (int k = 0; k < nen; k++)
      if (k < (int)inactive.size() && inactive[k]) ghostK[k] = h->mpmGhostSpring(p);
    mp = h->mpmParticles().m[p];
    smallK = h->mpmSmallKinematics();
  }
  resize();
  if (this->getDomain() != 0) refetchNodes();
  delete Ki;              // the initial stiffness belongs to the old support (dropKi)
  Ki = 0;
  index = -1;             // Element's Rayleigh scratch is keyed to numDOF (interface D1.4)
  st = ARMED;
  needUpdate = true;
  trialValid = false;
  return 0;
}

// ---------------------------------------------------------------------------
//  The trial state (interface D1.5)
// ---------------------------------------------------------------------------
int
LadrunoMPMPoint::computeTrial(void)
{
  LadrunoMPMPointHost *h = host();
  ladruno_mpm::Particles &PP = h->mpmParticles();
  const int nd = ndm, nn = nd * nd;
  for (int k = 0; k < nen; k++) {
    if (nodes[k] == 0) return -5;
    const Vector &u = nodes[k]->getTrialDisp();
    for (int a = 0; a < nd; a++) uT[nd * k + a] = u(a);
  }
  // F_rel = I + Σ u_a ⊗ ∇S_a ;  F_1 = F_rel · F_n
  for (int a = 0; a < nd; a++)
    for (int c = 0; c < nd; c++) {
      double s = (a == c) ? 1.0 : 0.0;
      if (Sc.empty())
        for (int k = 0; k < nen; k++) s += uT[nd * k + a] * dS[nd * k + c];
      else                                   // B8: component a's own gradient
        for (int k = 0; k < nen; k++) s += uT[nd * k + a] * dSc[(nd * k + a) * nd + c];
      Frel[nd * a + c] = s;
    }
  const double *Fn = &PP.F[nn * p];
  double F1loc[9], R1loc[9], deps[6], J;
  // WP-232j (J1, J9): the kernel takes dFbar F_n; Frel stays the actual increment
  // (the spatial gradients g_a, the trial centre, G2P's position all use it).
  double Fuse[9];
  for (int k = 0; k < nn; k++) Fuse[k] = Frel[k];
  alRel = 0.0;
  if (alActive()) {
    double Fc[9], gN[8][3];
    cornerGradients(gN);
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++) {
        double s = (a == c) ? 1.0 : 0.0;
        for (int k = 0; k < nCorner; k++) s += uT[nd * corner[k] + a] * gN[k][c];
        Fc[nd * a + c] = s;
      }
    if (ladruno_mpm::fbarIncrement(nd, Frel, Fc, Fuse) != 0) return -2;
    Jc = ladruno_mpm::detN(nd, Fc);
    if (nd == 3) ladruno_hypo::inv3(Fc, Fci);
    else {
      Fci[0] = Fc[3] / Jc;  Fci[1] = -Fc[1] / Jc;
      Fci[2] = -Fc[2] / Jc; Fci[3] = Fc[0] / Jc;
    }
    for (int k = 0; k < nCorner; k++)           // C1: ghat_k = gradN_k . dF_c^-1
      for (int j = 0; j < nd; j++) {
        double s = 0.0;
        for (int L = 0; L < nd; L++) s += gN[k][L] * Fci[nd * L + j];
        gHat[nd * k + j] = s;
      }
    diagJr = Jc / ladruno_mpm::detN(nd, Frel);
    double num = 0.0, den = 0.0;
    for (int k = 0; k < nn; k++) {
      num = std::max(num, std::fabs(Fuse[k] - Frel[k]));
      den = std::max(den, std::fabs(Frel[k] - ((k % (nd + 1) == 0) ? 1.0 : 0.0)));
    }
    alRel = (den > 0.0) ? num / den : num;
  }
  for (int a = 0; a < nd; a++)
    for (int c = 0; c < nd; c++) {
      double s = 0.0;
      for (int e = 0; e < nd; e++) s += Fuse[nd * a + e] * Fn[nd * e + c];
      F1loc[nd * a + c] = s;
    }
  int kst;
  if (smallK) {                          // WP-232l (P6): sym(dF - I), R = I, V = V0
    kst = ladruno_mpm::smallIncrement(nd, Fuse, deps, R1loc) ? 0 : -1;
    J = 1.0;
  } else {
    kst = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1loc, deps, R1loc, &J)
                    : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1loc, deps, R1loc, &J);
  }
  if (kst != 0) return -2;
  for (int k = 0; k < nv; k++) depsT[k] = deps[k];
  // the trial centre stays in the real grid (interface A7)
  double x[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < nd; a++) {
    double s = PP.x[nd * p + a];
    if (Sc.empty())
      for (int k = 0; k < nen; k++) s += S[k] * uT[nd * k + a];
    else
      for (int k = 0; k < nen; k++) s += Sc[nd * k + a] * uT[nd * k + a];
    x[a] = s;
  }
  if (!h->mpmTrialAdmissible(p, x)) return -3;
  // the constitutive update on the accumulated unrotated strain
  static Vector e3(3), e6(6);
  Vector &e = (nv == 6) ? e6 : e3;
  for (int k = 0; k < nv; k++) {
    epsT[k] = PP.eps[nv * p + k] + deps[k];
    e(k) = epsT[k];
  }
  if (mat->setTrialStrain(e) < 0) {
    LadrunoMPMPointHost *hh = host();             // WP-232k (D14): the host names the material
    if (hh != 0) hh->mpmTrialRefused(p);
    return -4;
  }
  if (diagnose) {                         // WP-232j: `mpm debug failDump`
    // COPY the current tangent before asking for the initial one: a material may return
    // both through ONE member buffer (DruckerPragerPlaneStrain.cpp:129-164 fills `tangent`
    // in both getters), and a held reference then reads the initial tangent twice -- the
    // 4fc274d51 dumps' constant ratio 1.0 (quirk WP-232j-ndmaterial-tangent-getters-share-a-buffer)
    const Matrix Dt(mat->getTangent());
    const Matrix &D0 = mat->getInitialTangent();
    double a2 = 0.0, b2 = 0.0;
    for (int r = 0; r < Dt.noRows(); r++)
      for (int c = 0; c < Dt.noCols(); c++) { a2 += Dt(r, c) * Dt(r, c); b2 += D0(r, c) * D0(r, c); }
    diagTan = (b2 > 0.0) ? std::sqrt(a2 / b2) : 1.0;
  }
  const Vector &sm = mat->getStress();
  for (int k = 0; k < nv; k++) sigHatT[k] = PP.sig0[nv * p + k] + sm(k);
  if (smallK) for (int k = 0; k < nv; k++) sigT[k] = sigHatT[k];   // WP-232l (P6)
  else if (nd == 3) ladruno_hypo::pushStress6(R1loc, sigHatT, sigT);
  else ladruno_hypo::pushStressPlaneStrain(R1loc, sigHatT, sigT);
  for (int k = 0; k < nn; k++) { F1[k] = F1loc[k]; R1[k] = R1loc[k]; }
  J1 = J;
  vol = J * PP.V0[p];
  // spatial gradients g_a = ∇S_a · F_rel⁻¹ (WP-232l, small: g_a = ∇S_a, F_rel⁻¹ taken as I)
  double Fi[9];
  if (smallK) {
    for (int k = 0; k < 9; k++) Fi[k] = 0.0;
    for (int a = 0; a < nd; a++) Fi[nd * a + a] = 1.0;
  } else if (nd == 3) {
    if (!(ladruno_hypo::inv3(Frel, Fi) > 0.0)) return -2;
  } else {
    const double det = Frel[0] * Frel[3] - Frel[1] * Frel[2];
    if (!(det > 0.0)) return -2;
    Fi[0] = Frel[3] / det;  Fi[1] = -Frel[1] / det;
    Fi[2] = -Frel[2] / det; Fi[3] = Frel[0] / det;
  }
  double sg[9];
  voigtToTensor(nd, sigT, sg);
  if (!Sc.empty()) {                     // B8: g'_{k,i} = grad S'_{k,i} . F_rel^-1
    for (int k = 0; k < nen; k++)
      for (int i = 0; i < nd; i++) {
        const double *dsi = &dSc[(nd * k + i) * nd];
        double *gi = &gc[(nd * k + i) * nd];
        for (int j = 0; j < nd; j++) {
          double s = 0.0;
          for (int K = 0; K < nd; K++) s += dsi[K] * Fi[nd * K + j];
          gi[j] = s;
        }
        double s = 0.0;
        for (int j = 0; j < nd; j++) s += sg[3 * i + j] * gi[j];
        fInt[nd * k + i] = (LADRUNO_MUTATE_MPM == LADRUNO_MUT_ZERO) ? 0.0 : vol * s;
      }
    addConjugateForce(sg);                 // WP-232j (C1)
    trialValid = true;
    needUpdate = false;
    return 0;
  }
  for (int k = 0; k < nen; k++) {
    for (int j = 0; j < nd; j++) {
      double s = 0.0;
      for (int K = 0; K < nd; K++) s += dS[nd * k + K] * Fi[nd * K + j];
      g[nd * k + j] = s;
    }
    for (int i = 0; i < nd; i++) {
      double s = 0.0;
      for (int j = 0; j < nd; j++) s += sg[3 * i + j] * g[nd * k + j];
      // MPM mutation family, ZERO: the internal force is dropped (ADR 232 §8)
      fInt[nd * k + i] = (LADRUNO_MUTATE_MPM == LADRUNO_MUT_ZERO) ? 0.0 : vol * s;
    }
  }
  addConjugateForce(sg);                   // WP-232j (C1)
  trialValid = true;
  needUpdate = false;
  return 0;
}

// WP-232j (C1, interface §11): the internal force work-conjugate to the F-bar kinematics,
//   f_a = w [sigma g_a + pbar (ghat_a - g_a)],  pbar = tr_d(sigma) / d,
// ghat on the cell's corners only: Sum_a f_a . v_a = w sigma : lbar, lbar = l + beta I_d.
// The per-component (B8) path subtracts pbar g'_{a,i,i}. Nothing off or for a fallback.
void
LadrunoMPMPoint::addConjugateForce(const double *sg)
{
  if (!alActive() || LADRUNO_MUTATE_MPM == LADRUNO_MUT_ZERO) return;
  const int nd = ndm;
  double pbar = 0.0;
  for (int i = 0; i < nd; i++) pbar += sg[3 * i + i];
  pbar /= nd;
  for (int k = 0; k < nen; k++)
    for (int i = 0; i < nd; i++)
      fInt[nd * k + i] -= vol * pbar * (Sc.empty() ? g[nd * k + i] : gc[(nd * k + i) * nd + i]);
  for (int k = 0; k < nCorner; k++)
    for (int i = 0; i < nd; i++) fInt[nd * corner[k] + i] += vol * pbar * gHat[nd * k + i];
}

int
LadrunoMPMPoint::update(void)
{
  if (st == RESTORED) {        // Domain::addElement's update() after a restore
    st = COMMITTED;
    return 0;
  }
  LadrunoMPMPointHost *h = host();
  // review item 1: armed AND inside the host's `mpm analyze` / tangentCheck window; a plain
  // `analyze` after a failed `mpm analyze` finds the elements still ARMED.
  if (h == 0 || st != ARMED || !h->mpmInStep()) {
    static int nWarn = 0;
    if (nWarn++ < 3)
      opserr << "ERROR LadrunoMPMPoint " << this->getTag() << " -- "
             << (h == 0 ? "its MPM host is gone (remove loadPattern)"
                        : "a step not prepared by `mpm analyze`")
             << ": mode-A particle elements move only through `mpm analyze` (ADR 232 §5.3)\n";
    return -1;
  }
  const int failMode = h->mpmFailThisAttempt();     // `mpm debug failAttempt`
  if (failMode == 1) {
    trialValid = false;
    return -1;
  }
  const int rc = computeTrial();
  if (failMode == 2) {                 // after the trial: the material carries a trial state
    trialValid = false;
    return -1;
  }
  if (rc < 0) {
    trialValid = false;
    const char *why = (rc == -2) ? "det(F_mid) <= 0 or det(F) <= 0 (inverted or collapsed)"
                    : (rc == -3) ? "its trial centre leaves the grid (interface A7)"
                    : (rc == -4) ? "its material failed to integrate the trial strain"
                                 : "a support node is missing";
    opserr << "LadrunoMPMPoint " << this->getTag() << " (particle index " << p
           << ") -- trial refused: " << why << "\n";
    return -1;
  }
  return 0;
}

int
LadrunoMPMPoint::commitState(void)
{
  int rc = this->Element::commitState();
  if (mat != 0) mat->commitState();
  LadrunoMPMPointHost *h = host();
  if (h != 0 && st == ARMED) h->mpmParticleCommitted(p);   // the batched G2P (D1.9)
  return rc;
}

int
LadrunoMPMPoint::revertToLastCommit(void)
{
  // Freeze decision 7: restore the material's trial state and re-arm, so the
  // driver's retry starts from the committed state.
  int rc = (mat != 0) ? mat->revertToLastCommit() : 0;
  if (st == ARMED) {
    needUpdate = true;
    trialValid = false;
  }
  return rc;
}

int
LadrunoMPMPoint::revertToStart(void)
{
  // review item 7: a plain `reset` (Domain::revertToStart) would revert the material while the
  // host keeps the particle's deformed x / F / strain -- plastic history lost silently. While a
  // host owns the particle, refuse it loudly; `mpm reset` reverts particles and materials together.
  if (host() != 0) {
    static int nWarn = 0;
    if (nWarn++ < 1)
      opserr << "WARNING LadrunoMPMPoint -- `reset` does not reach MPM particles: their state is"
                " the host's; use `mpm reset` (particles and materials together). Nothing reverted.\n";
    return 0;
  }
  int rc = (mat != 0) ? mat->revertToStart() : 0;
  needUpdate = true;
  trialValid = false;
  return rc;
}

void
LadrunoMPMPoint::g2p(ladruno_mpm::Particles &PP, double flip, bool quasiStatic)
{
  const int nd = ndm, nn = nd * nd;
  for (int a = 0; a < nd; a++) {
    double dx = 0.0;
    if (Sc.empty())
      for (int k = 0; k < nen; k++) dx += S[k] * nodes[k]->getTrialDisp()(a);
    else
      for (int k = 0; k < nen; k++) dx += Sc[nd * k + a] * nodes[k]->getTrialDisp()(a);
    // MPM mutation family, SCALE: the G2P position increment x1.5 (ADR 232 §8)
    if (LADRUNO_MUTATE_MPM == LADRUNO_MUT_SCALE) dx *= ladruno_mpm::MUT_SCALE_FACTOR;
    PP.x[nd * p + a] += dx;
  }
  for (int k = 0; k < nn; k++) { PP.F[nn * p + k] = F1[k]; PP.R[nn * p + k] = R1[k]; }
  PP.V[p] = vol;
  for (int k = 0; k < nv; k++) {
    PP.eps[nv * p + k] = epsT[k];
    PP.sigHat[nv * p + k] = sigHatT[k];
    PP.sig[nv * p + k] = sigT[k];
  }
  if (!quasiStatic) {
    for (int a = 0; a < nd; a++) {
      double dv = 0.0, vpic = 0.0, da = 0.0, apic = 0.0;
      for (int k = 0; k < nen; k++) {
        const double v = nodes[k]->getTrialVel()(a), ac = nodes[k]->getTrialAccel()(a);
        const double Sk = Sc.empty() ? S[k] : Sc[nd * k + a];
        dv += Sk * (v - vHat[nd * k + a]);
        vpic += Sk * v;
        da += Sk * (ac - aHat[nd * k + a]);
        apic += Sk * ac;
      }
      double &vp = PP.v[nd * p + a];
      double &ap = PP.a[nd * p + a];
      vp = flip * (vp + dv) + (1.0 - flip) * vpic;
      ap = flip * (ap + da) + (1.0 - flip) * apic;
    }
  }
  // WP-232l (P5): p, pd, pdd <- alpha_p (x + sum Sp (x_a - xref_a)) + (1 - alpha_p) sum Sp x_a,
  // xref the unconditioned projection (pRef)
  if (tpOn && (int)PP.pw.size() > p) {
    LadrunoMPMPointHost *h = host();
    const double ap = (h != 0) ? h->mpmPressureBlend() : 1.0;
    double inc[3] = {0.0, 0.0, 0.0}, pic[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < nen; k++) {
      const double val[3] = {nodes[k]->getTrialDisp()(ndm), nodes[k]->getTrialVel()(ndm),
                             nodes[k]->getTrialAccel()(ndm)};
      for (int q = 0; q < 3; q++) {
        inc[q] += Sp[k] * (val[q] - pRef[3 * k + q]);
        pic[q] += Sp[k] * val[q];
      }
    }
    // WP-232l fix 2: the interpolated pressure carries the known affine term pi; the FLIP
    // increment adds pi_end - pi_start, the PIC part pi_end (pi-dot neglected: interface 7b)
    if (affOn) {
      double g3[3];
      const double piEnd = piNow(g3);
      inc[0] += piEnd - piStart_;
      pic[0] += piEnd;
    }
    double *tgt[3] = {&PP.pw[p], &PP.pwd[p], &PP.pwdd[p]};
    for (int q = 0; q < 3; q++) {
      if (quasiStatic && q > 0) break;
      *tgt[q] = ap * (*tgt[q] + inc[q]) + (1.0 - ap) * pic[q];
    }
  }
  st = COMMITTED;
}

int
LadrunoMPMPoint::tangentCheck(double scale, double eps, double &rel, double &normK)
{
  if (st != ARMED || mat == 0) return -1;
  if (ndfN != ndm) return -1;        // WP-232l: the u-only check; a two-phase epoch skips it
  const int n = nen * ndm;
  std::vector<Vector> saved;
  for (int k = 0; k < nen; k++) saved.push_back(nodes[k]->getTrialDisp());
  std::vector<double> U(n);
  for (int j = 0; j < n; j++) U[j] = scale * std::sin(1.7 * j + 0.3);
  Vector u(ndm);
  auto setU = [&](const std::vector<double> &W) {
    for (int k = 0; k < nen; k++) {
      for (int a = 0; a < ndm; a++) u(a) = W[ndm * k + a];
      nodes[k]->setTrialDisp(u);
    }
  };
  auto restore = [&]() {
    for (int k = 0; k < nen; k++) nodes[k]->setTrialDisp(saved[k]);
    mat->revertToLastCommit();
    needUpdate = true;
    trialValid = false;
  };
  setU(U);
  if (this->update() < 0) { restore(); return -1; }
  Matrix K(this->getTangentStiff());
  Matrix KFD(n, n);
  int bad = 0;
  for (int j = 0; j < n; j++) {
    const double u0 = U[j];
    U[j] = u0 + eps;
    setU(U);
    bad += this->update() < 0;
    Vector Pp(this->getResistingForce());
    U[j] = u0 - eps;
    setU(U);
    bad += this->update() < 0;
    Vector Pm(this->getResistingForce());
    U[j] = u0;
    for (int i = 0; i < n; i++) KFD(i, j) = (Pp(i) - Pm(i)) / (2.0 * eps);
  }
  double d2 = 0.0, k2 = 0.0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      d2 += (K(i, j) - KFD(i, j)) * (K(i, j) - KFD(i, j));
      k2 += K(i, j) * K(i, j);
    }
  restore();
  if (bad) return -1;
  normK = std::sqrt(k2);
  rel = (k2 > 0.0) ? std::sqrt(d2 / k2) : 0.0;
  return 0;
}

// ---------------------------------------------------------------------------
//  Force, tangent, mass (interface D1.6, D1.7)
// ---------------------------------------------------------------------------
void
LadrunoMPMPoint::assembleForce(void)
{
  if (ndfN != ndm) { assembleForceUP(); return; }     // WP-232l
  LadrunoMPMPointHost *h = host();
  if (h != 0) {
    const double *b = h->mpmBodyForce(p);
    for (int a = 0; a < ndm; a++) bf[a] = b[a];
  }
  for (int k = 0; k < nen; k++)
    for (int a = 0; a < ndm; a++) {
      const int r = ndm * k + a;
      const double Sk = Sc.empty() ? S[k] : Sc[ndm * k + a];   // B8
      (*P)(r) = fInt[r] - Sk * bf[a] + ghostK[k] * uT[r];
    }
}

const Vector &
LadrunoMPMPoint::getResistingForce(void)
{
  if (P == 0) resize();
  if (st == ARMED) {
    if (needUpdate || !trialValid) this->update();
    if (trialValid) assembleForce();
  }
  return *P;
}

const Vector &
LadrunoMPMPoint::getResistingForceIncInertia(void)
{
  if (P == 0) resize();
  if (st != ARMED) return *Pinc;     // COMMITTED: the converged cache
  if (ndfN != ndm) {
    // WP-232l (interface §1, LadrunoUP.cpp:1682-1708): R + M a + Damp v, every node DOF; the
    // p-slot "velocity" is pdot. Damp carries the solid-only Rayleigh (obligation 2).
    Vector R(this->getResistingForce());
    const int nf = nen * ndfN;
    Vector acc(nf), vel(nf);
    for (int k = 0; k < nen; k++) {
      const Vector &ak = nodes[k]->getTrialAccel();
      const Vector &vk = nodes[k]->getTrialVel();
      for (int c = 0; c < ndfN; c++) { acc(ndfN * k + c) = ak(c); vel(ndfN * k + c) = vk(c); }
    }
    R.addMatrixVector(1.0, this->getMass(), acc, 1.0);
    const Matrix C(this->getDamp());
    R.addMatrixVector(1.0, C, vel, 1.0);
    *Pinc = R;
    return *Pinc;
  }
  // snapshot the static residual into a LOCAL before the Rayleigh forces
  // (betaK's getTangentStiff may re-enter; element guide, lint L1)
  Vector R(this->getResistingForce());
  const Matrix &M = this->getMass();
  Vector acc(nen * ndm);
  LadrunoShell::addNodalInertia(R, M, nodes.data(), nen, ndm, true, acc);
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    R += this->getRayleighDampingForces();
  *Pinc = R;
  return *Pinc;
}

void
LadrunoMPMPoint::tangent(Matrix &K, bool initial)
{
  K.Zero();
  LadrunoMPMPointHost *h = host();
  const int nd = ndm;
  double Rm[9], volm, sg[9];
  const double *grad;
  if (initial) {
    if (h == 0) return;
    ladruno_mpm::Particles &PP = h->mpmParticles();
    for (int k = 0; k < nd * nd; k++) Rm[k] = PP.R[nd * nd * p + k];
    volm = PP.V[p];
    grad = dS.data();
    for (int k = 0; k < 9; k++) sg[k] = 0.0;
  } else {
    for (int k = 0; k < nd * nd; k++) Rm[k] = R1[k];
    volm = vol;
    grad = g.data();
    voigtToTensor(nd, sigT, sg);
  }
  const Matrix &Dm = initial ? mat->getInitialTangent() : mat->getTangent();
  double Dc[36];
  if (smallK) {                          // WP-232l (P6): R = I, the material tangent as is
    for (int r = 0; r < nv; r++)
      for (int c = 0; c < nv; c++) Dc[nv * r + c] = Dm(r, c);
  } else if (nd == 3) {
    double Dmat[36];
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++) Dmat[6 * r + c] = Dm(r, c);
    ladruno_hypo::pushTangent6(Rm, Dmat, Dc);
  } else {
    double Dmat[9];
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++) Dmat[3 * r + c] = Dm(r, c);
    ladruno_hypo::pushTangentPlaneStrain(Rm, Dmat, Dc);
  }
  if (!Sc.empty()) {                     // WP-232i B8
    tangentPerComponent(K, initial, Rm, volm, sg, Dc);
    addFbarTangent(K, initial, volm, sg, Dc);      // WP-232j (J10)
    return;
  }
  // material part: B_a rows {xx,yy,zz,xy,yz,zx} (3D) or {xx,yy,xy} (2D)
  double Ba[6][3], Bb[6][3], DB[6][3];
  for (int a = 0; a < nen; a++) {
    const double *ga = &grad[nd * a];
    for (int r = 0; r < nv; r++) for (int c = 0; c < 3; c++) Ba[r][c] = 0.0;
    if (nd == 3) {
      Ba[0][0] = ga[0]; Ba[1][1] = ga[1]; Ba[2][2] = ga[2];
      Ba[3][0] = ga[1]; Ba[3][1] = ga[0];
      Ba[4][1] = ga[2]; Ba[4][2] = ga[1];
      Ba[5][0] = ga[2]; Ba[5][2] = ga[0];
    } else {
      Ba[0][0] = ga[0]; Ba[1][1] = ga[1];
      Ba[2][0] = ga[1]; Ba[2][1] = ga[0];
    }
    for (int b = 0; b < nen; b++) {
      const double *gb = &grad[nd * b];
      for (int r = 0; r < nv; r++) for (int c = 0; c < 3; c++) Bb[r][c] = 0.0;
      if (nd == 3) {
        Bb[0][0] = gb[0]; Bb[1][1] = gb[1]; Bb[2][2] = gb[2];
        Bb[3][0] = gb[1]; Bb[3][1] = gb[0];
        Bb[4][1] = gb[2]; Bb[4][2] = gb[1];
        Bb[5][0] = gb[2]; Bb[5][2] = gb[0];
      } else {
        Bb[0][0] = gb[0]; Bb[1][1] = gb[1];
        Bb[2][0] = gb[1]; Bb[2][1] = gb[0];
      }
      for (int r = 0; r < nv; r++)
        for (int c = 0; c < nd; c++) {
          double s = 0.0;
          for (int q = 0; q < nv; q++) s += Dc[nv * r + q] * Bb[q][c];
          DB[r][c] = s;
        }
      for (int i = 0; i < nd; i++)
        for (int k = 0; k < nd; k++) {
          double s = 0.0;
          for (int r = 0; r < nv; r++) s += Ba[r][i] * DB[r][k];
          K(nd * a + i, nd * b + k) += s * volm;
        }
      // geometric (initial-stress) term: δ_ik g_a·σ·g_b V (none under `small`, WP-232l P6)
      if (!initial && !smallK) {
        double s = 0.0;
        for (int j = 0; j < nd; j++)
          for (int l = 0; l < nd; l++) s += ga[j] * sg[3 * j + l] * gb[l];
        for (int i = 0; i < nd; i++) K(nd * a + i, nd * b + i) += s * volm;
      }
    }
    for (int i = 0; i < nd; i++) K(nd * a + i, nd * a + i) += ghostK[a];   // ghost spring (D2.3)
  }
  addFbarTangent(K, initial, volm, sg, Dc);        // WP-232j (J10)
}

// WP-232j (interface §11, J10 re-derived for C1): the anti-locking part of the consistent
// tangent of f_ai = w [sigma_ij g_aj + pbar D_ai], D_ai = ghat_ai - g_ai, under the model that
// makes the shipped assembly exact (delta(J sigma) = J (c:d + l sigma + sigma l^T), applied to
// the F-bar kinematics lbar = l + beta I_d, beta = (1/d) Sum_bk D_bk du_bk, w the carried volume):
//   K += h D^T                                   (sigma's response to beta)
//      + D_ai w [ (q^T B_b)_k + (2/d)(sigma g_b)_k + (s + 2 pbar / d) D_bk ]   (d(w pbar))
//      + w pbar (g_ak g_bi - ghat_ak ghat_bi)    (d D_ai)
//   h_ai = w Sum_j Q_ij g_aj,  Q = (1/d)(c:I_d + 2 sigma),  q_col = (1/d) Sum_{r<d} c[r][col],
//   s = (1/d) Sum_{col<d} q_col,  pbar = tr_d(sigma) / d.
// The per-component (B8) path uses g'_{a,i,.} for component i: g_aj -> g'_{a,i,j},
// D_ai = ghat_ai - g'_{a,i,i}, (sigma g_b)_k -> Sum_m g'_{b,k,m} sigma_mk, g_ak g_bi ->
// g'_{a,i,k} g'_{b,k,i}. Initial: dF_c = I, sigma = 0, the initial tangent, the step-start
// gradients. Nothing for a fallback particle or with it off. Plane strain: d = 2, in plane.
void
LadrunoMPMPoint::fillB(int a, bool initial, double B[6][3]) const
{
  const int nd = ndm;
  const bool pc = !Sc.empty();
  const double *grad = pc ? (initial ? dSc.data() : gc.data()) : (initial ? dS.data() : g.data());
  auto G = [&](int i) { return pc ? &grad[(nd * a + i) * nd] : &grad[nd * a]; };
  for (int r = 0; r < 6; r++) for (int c = 0; c < 3; c++) B[r][c] = 0.0;
  const double *g0 = G(0), *g1 = G(1), *g2 = (nd == 3) ? G(2) : 0;
  if (nd == 3) {
    B[0][0] = g0[0]; B[1][1] = g1[1]; B[2][2] = g2[2];
    B[3][0] = g0[1]; B[3][1] = g1[0];
    B[4][1] = g1[2]; B[4][2] = g2[1];
    B[5][0] = g0[2]; B[5][2] = g2[0];
  } else {
    B[0][0] = g0[0]; B[1][1] = g1[1];
    B[2][0] = g0[1]; B[2][1] = g1[0];
  }
}

void
LadrunoMPMPoint::addFbarTangent(Matrix &K, bool initial, double volm, const double *sg,
                                const double *Dc)
{
  if (!alActive()) return;
  const int nd = ndm, n = nen * nd;
  const bool pc = !Sc.empty();
  const double *grad = pc ? (initial ? dSc.data() : gc.data()) : (initial ? dS.data() : g.data());
  // gradient of component i of node a, entry j
  auto GR = [&](int a, int i, int j) { return pc ? grad[(nd * a + i) * nd + j] : grad[nd * a + j]; };
  double Qv[6], Q[9], q[6], s = 0.0, pbar = 0.0;
  for (int r = 0; r < nv; r++) {
    double t = 0.0;
    for (int c = 0; c < nd; c++) t += Dc[nv * r + c];
    Qv[r] = t;
  }
  for (int c = 0; c < nv; c++) {
    double t = 0.0;
    for (int r = 0; r < nd; r++) t += Dc[nv * r + c];
    q[c] = t / nd;
  }
  for (int c = 0; c < nd; c++) s += q[c];
  s /= nd;
  voigtToTensor(nd, Qv, Q);
  for (int k = 0; k < 9; k++) Q[k] = (Q[k] + 2.0 * sg[k]) / nd;
  for (int i = 0; i < nd; i++) pbar += sg[3 * i + i];
  pbar /= nd;
  // ghat per node slot (0 off the corners): initial dF_c = I
  static std::vector<double> gh, h, D, row;
  gh.assign(n, 0.0);
  h.assign(n, 0.0);
  D.assign(n, 0.0);
  row.assign(n, 0.0);
  double gN[8][3];
  cornerGradients(gN);
  for (int k = 0; k < nCorner; k++)
    for (int j = 0; j < nd; j++) gh[nd * corner[k] + j] = initial ? gN[k][j] : gHat[nd * k + j];
  for (int a = 0; a < nen; a++)
    for (int i = 0; i < nd; i++) {
      double t = 0.0;
      for (int j = 0; j < nd; j++) t += Q[3 * i + j] * GR(a, i, j);
      h[nd * a + i] = volm * t;
      D[nd * a + i] = gh[nd * a + i] - GR(a, i, i);
    }
  // row_bk = w [ (q^T B_b)_k + (2/d)(sigma g_b)_k + (s + 2 pbar/d) D_bk ]
  double B[6][3];
  for (int b = 0; b < nen; b++) {
    fillB(b, initial, B);
    for (int k = 0; k < nd; k++) {
      double t = 0.0;
      for (int c = 0; c < nv; c++) t += q[c] * B[c][k];
      double sgb = 0.0;
      for (int m = 0; m < nd; m++) sgb += GR(b, k, m) * sg[3 * m + k];
      row[nd * b + k] = volm * (t + 2.0 / nd * sgb + (s + 2.0 * pbar / nd) * D[nd * b + k]);
    }
  }
  for (int a = 0; a < nen; a++)
    for (int i = 0; i < nd; i++) {
      const int r = nd * a + i;
      for (int b = 0; b < nen; b++)
        for (int k = 0; k < nd; k++) {
          const int c = nd * b + k;
          double v = h[r] * D[c] + D[r] * row[c];
          if (pbar != 0.0)
            v += volm * pbar * (GR(a, i, k) * GR(b, k, i) - gh[nd * a + k] * gh[nd * b + i]);
          K(r, c) += v;
        }
    }
}

// The Q1 gradients at the cell centre in the host's corner order
// (ladruno_mpm::cellCorners: bit d of k = the +h side on axis d).
void
LadrunoMPMPoint::cornerGradients(double gN[8][3]) const
{
  const double a = (ndm == 3) ? 0.25 / cellH : 0.5 / cellH;
  for (int k = 0; k < 8; k++)
    for (int d = 0; d < 3; d++) gN[k][d] = (d < ndm) ? (((k >> d) & 1) ? a : -a) : 0.0;
}

void
LadrunoMPMPoint::antiLockingTrial(double &rel, double &depsMax) const
{
  rel = alRel;
  depsMax = 0.0;
  for (int k = 0; k < nv; k++) depsMax = std::max(depsMax, std::fabs(depsT[k]));
}

double
LadrunoMPMPoint::trialStressChange(void) const
{
  LadrunoMPMPointHost *h = host();
  if (h == 0) return 0.0;
  const ladruno_mpm::Particles &PP = h->mpmParticles();
  double m = 0.0;
  for (int k = 0; k < nv; k++) m = std::max(m, std::fabs(sigHatT[k] - PP.sigHat[nv * p + k]));
  return m;
}

// WP-232i B8: the tangent with each displacement component's own gradient,
// g'_{a,i} (= T^T K T for the images' reduction). Same algebra as tangent().
void
LadrunoMPMPoint::tangentPerComponent(Matrix &K, bool initial, const double *, double volm,
                                     const double *sg, const double *Dc)
{
  const int nd = ndm;
  const double *grad = initial ? dSc.data() : gc.data();   // [(nd*a + i)*nd + j]
  auto G = [&](int a, int i) { return &grad[(nd * a + i) * nd]; };
  auto fillB = [&](int a, double B[6][3]) {
    for (int r = 0; r < nv; r++) for (int c = 0; c < 3; c++) B[r][c] = 0.0;
    const double *g0 = G(a, 0), *g1 = G(a, 1), *g2 = (nd == 3) ? G(a, 2) : 0;
    if (nd == 3) {
      B[0][0] = g0[0]; B[1][1] = g1[1]; B[2][2] = g2[2];
      B[3][0] = g0[1]; B[3][1] = g1[0];
      B[4][1] = g1[2]; B[4][2] = g2[1];
      B[5][0] = g0[2]; B[5][2] = g2[0];
    } else {
      B[0][0] = g0[0]; B[1][1] = g1[1];
      B[2][0] = g0[1]; B[2][1] = g1[0];
    }
  };
  double Ba[6][3], Bb[6][3], DB[6][3];
  for (int a = 0; a < nen; a++) {
    fillB(a, Ba);
    for (int b = 0; b < nen; b++) {
      fillB(b, Bb);
      for (int r = 0; r < nv; r++)
        for (int c = 0; c < nd; c++) {
          double s = 0.0;
          for (int q = 0; q < nv; q++) s += Dc[nv * r + q] * Bb[q][c];
          DB[r][c] = s;
        }
      for (int i = 0; i < nd; i++)
        for (int k = 0; k < nd; k++) {
          double s = 0.0;
          for (int r = 0; r < nv; r++) s += Ba[r][i] * DB[r][k];
          K(nd * a + i, nd * b + k) += s * volm;
        }
      if (!initial && !smallK)         // geometric term: g'_{a,i} . sigma . g'_{b,i} V
        for (int i = 0; i < nd; i++) {
          const double *ga = G(a, i), *gb = G(b, i);
          double s = 0.0;
          for (int j = 0; j < nd; j++)
            for (int l = 0; l < nd; l++) s += ga[j] * sg[3 * j + l] * gb[l];
          K(nd * a + i, nd * b + i) += s * volm;
        }
    }
    for (int i = 0; i < nd; i++) K(nd * a + i, nd * a + i) += ghostK[a];   // ghost spring (D2.3)
  }
}

const Matrix &
LadrunoMPMPoint::getTangentStiff(void)
{
  if (ndfN != ndm) {
    // WP-232l (P3): [K_u, -Q; 0, H] on the node-major [u, p] layout (LadrunoUP.cpp:1360-1408)
    const int nf = nen * ndfN;
    Matrix &K = scratchKF(nf > 0 ? nf : 1);
    K.Zero();
    if (nf == 0) return K;
    if (st == ARMED && (needUpdate || !trialValid)) this->update();
    Matrix &Kc = scratchK(nen * ndm);
    compactTangent(Kc, false);
    for (int k = 0; k < nen; k++)
      for (int a = 0; a < ndm; a++)
        for (int l = 0; l < nen; l++)
          for (int b = 0; b < ndm; b++) K(dofU(k, a), dofU(l, b)) = Kc(ndm * k + a, ndm * l + b);
    if (tpOn)
      for (int k = 0; k < nen; k++)
        for (int b = 0; b < nen; b++) {
          for (int i = 0; i < ndm; i++) K(dofU(k, i), dofP(b)) = -Qb[(ndm * k + i) * nen + b];
          K(dofP(k), dofP(b)) = Hb[k * nen + b];
        }
    return K;
  }
  const int n = nen * ndm;
  Matrix &K = scratchK(n > 0 ? n : 1);
  if (n == 0) return K;
  if (st == ARMED && (needUpdate || !trialValid)) this->update();
  tangent(K, false);
  return K;
}

void
LadrunoMPMPoint::compactTangent(Matrix &Kc, bool initial)
{
  tangent(Kc, initial);
}

const Matrix &
LadrunoMPMPoint::getInitialStiff(void)
{
  if (ndfN != ndm) {
    // WP-232l: [K0_u, -Q; 0, H] -- p rows NOT zero (ADR-71 §3.2, LadrunoUP.cpp:1417-1448)
    const int nf = nen * ndfN;
    if (Ki == 0) {
      Ki = new Matrix(nf > 0 ? nf : 1, nf > 0 ? nf : 1);
      if (nf > 0 && mat != 0) {
        Matrix Kc(nen * ndm, nen * ndm);
        compactTangent(Kc, true);
        for (int k = 0; k < nen; k++)
          for (int a = 0; a < ndm; a++)
            for (int l = 0; l < nen; l++)
              for (int b = 0; b < ndm; b++) (*Ki)(dofU(k, a), dofU(l, b)) = Kc(ndm * k + a, ndm * l + b);
        if (tpOn)
          for (int k = 0; k < nen; k++)
            for (int b = 0; b < nen; b++) {
              for (int i = 0; i < ndm; i++) (*Ki)(dofU(k, i), dofP(b)) = -Qb[(ndm * k + i) * nen + b];
              (*Ki)(dofP(k), dofP(b)) = Hb[k * nen + b];
            }
      }
    }
    return *Ki;
  }
  const int n = nen * ndm;
  if (Ki == 0) {
    Ki = new Matrix(n > 0 ? n : 1, n > 0 ? n : 1);
    if (n > 0 && mat != 0) tangent(*Ki, true);
  }
  return *Ki;
}

const Matrix &
LadrunoMPMPoint::getMass(void)
{
  // one mass model: lumped, M_aa = S_a m_p I (the row sum of the consistent
  // N m_p Nᵀ, Σ S_a = 1). The inertia residual uses this same matrix (D1.7).
  if (ndfN != ndm) {                    // WP-232l: [M_u, 0; 0, 0]
    const int nf = nen * ndfN;
    Matrix &M = scratchMF(nf > 0 ? nf : 1);
    M.Zero();
    for (int k = 0; k < nen; k++)
      for (int a = 0; a < ndm; a++) M(dofU(k, a), dofU(k, a)) = S[k] * mp;
    return M;
  }
  const int n = nen * ndm;
  Matrix &M = scratchM(n > 0 ? n : 1);
  M.Zero();
  for (int k = 0; k < nen; k++)
    for (int a = 0; a < ndm; a++) M(ndm * k + a, ndm * k + a) = S[k] * mp;
  return M;
}

// WP-232l (P3): [C_Ray, 0; Q^T, S + alpha_st H~] (LadrunoUP.cpp:1450-1503). C_Ray =
// alphaM M_u + betaK K_u(current), u-u only: never Element::getDamp in a two-phase epoch
// (it would compose betaK getTangentStiff, which carries -Q and H; ADR-71 §3.2, obl. 2).
const Matrix &
LadrunoMPMPoint::getDamp(void)
{
  if (ndfN == ndm) return this->Element::getDamp();
  const int nf = nen * ndfN;
  Matrix &C = scratchC(nf > 0 ? nf : 1);
  C.Zero();
  if (nf == 0) return C;
  if (alphaM != 0.0)
    for (int k = 0; k < nen; k++)
      for (int a = 0; a < ndm; a++) C(dofU(k, a), dofU(k, a)) += alphaM * S[k] * mp;
  if (betaK != 0.0) {
    if (st == ARMED && (needUpdate || !trialValid)) this->update();
    Matrix &Kc = scratchK(nen * ndm);
    compactTangent(Kc, false);
    for (int k = 0; k < nen; k++)
      for (int a = 0; a < ndm; a++)
        for (int l = 0; l < nen; l++)
          for (int b = 0; b < ndm; b++) C(dofU(k, a), dofU(l, b)) += betaK * Kc(ndm * k + a, ndm * l + b);
  }
  if (tpOn)
    for (int k = 0; k < nen; k++)
      for (int b = 0; b < nen; b++) {
        for (int i = 0; i < ndm; i++) C(dofP(b), dofU(k, i)) = Qb[(ndm * k + i) * nen + b];
        C(dofP(k), dofP(b)) = Sb[k * nen + b] + stabA * Htb[k * nen + b];
      }
  return C;
}

// WP-232l (interface §1): R_static = [f' - S bf + ghost - Q p ; H p - f_seep], every node DOF
void
LadrunoMPMPoint::assembleForceUP(void)
{
  LadrunoMPMPointHost *h = host();
  if (h != 0) {
    const double *b = h->mpmBodyForce(p);
    for (int a = 0; a < ndm; a++) bf[a] = b[a];
  }
  P->Zero();
  for (int k = 0; k < nen; k++)
    for (int a = 0; a < ndm; a++) {
      const int r = ndm * k + a;
      const double Sk = Sc.empty() ? S[k] : Sc[ndm * k + a];
      (*P)(dofU(k, a)) = fInt[r] - Sk * bf[a] + ghostK[k] * uT[r];
    }
  if (!tpOn) return;
  pT.assign(nen, 0.0);
  for (int k = 0; k < nen; k++) pT[k] = nodes[k]->getTrialDisp()(ndm);
  for (int k = 0; k < nen; k++)
    for (int i = 0; i < ndm; i++) {
      double s = 0.0;
      for (int b = 0; b < nen; b++) s += Qb[(ndm * k + i) * nen + b] * pT[b];
      (*P)(dofU(k, i)) -= s;
    }
  const double *drv = (h != 0) ? h->mpmFluidDrive(p) : 0;
  for (int b = 0; b < nen; b++) {
    double s = 0.0;
    for (int c = 0; c < nen; c++) s += Hb[b * nen + c] * pT[c];
    if (drv != 0)
      for (int i = 0; i < ndm; i++) s -= Gs[b * ndm + i] * drv[i];
    (*P)(dofP(b)) += s;
  }
  // WP-232l fix 2 (interface 7b): the affine image's known pressure pi in the coupling force
  // (-alpha V dS'_{a,i} pi) and in the Darcy flux (V dS^p_b . kbar grad pi)
  if (affOn && h != 0) {
    double gpi[3] = {0, 0, 0};
    const double pi = piNow(gpi);
    const double Vp = h->mpmParticles().V[p];
    for (int k = 0; k < nen; k++)
      for (int i = 0; i < ndm; i++) {
        const double dn = Sc.empty() ? dS[ndm * k + i] : dSc[(ndm * k + i) * ndm + i];
        (*P)(dofU(k, i)) -= tpAlpha * Vp * dn * pi;
      }
    for (int b = 0; b < nen; b++) {
      double s = 0.0;
      for (int i = 0; i < ndm; i++) s += dSp[ndm * b + i] * tpK[i] * gpi[i];
      (*P)(dofP(b)) += Vp * s;
    }
  }
}

// WP-232l fix 2: pi = rhoF sum_d drive_d a_d (drive = lambda b_f, the host's), its gradient
double
LadrunoMPMPoint::piNow(double *grad) const
{
  for (int i = 0; i < ndm; i++) grad[i] = 0.0;
  LadrunoMPMPointHost *h = host();
  const double *drv = (h != 0) ? h->mpmFluidDrive(p) : 0;
  if (!affOn || drv == 0) return 0.0;
  double pi = 0.0;
  for (int d = 0; d < ndm; d++) {
    pi += tpRhoF * drv[d] * affC[d * (ndm + 1)];
    for (int i = 0; i < ndm; i++) grad[i] += tpRhoF * drv[d] * affC[d * (ndm + 1) + 1 + i];
  }
  return pi;
}

// `mpm debug particleBlocks` (tests/mpm_twophase_interface.md §7a Q1, Q2)
int
LadrunoMPMPoint::writeBlocks(const double sc[6], std::ostream &os, double lambda, long long id)
{
  LadrunoMPMPointHost *h = host();
  if (st != ARMED || mat == 0 || h == 0 || (int)nodes.size() != nen) return -1;
  for (int k = 0; k < nen; k++) if (nodes[k] == 0) return -1;
  const int nf = nen * ndfN;
  std::vector<Vector> sd, sv, sa;
  for (int k = 0; k < nen; k++) {
    sd.push_back(nodes[k]->getTrialDisp());
    sv.push_back(nodes[k]->getTrialVel());
    sa.push_back(nodes[k]->getTrialAccel());
  }
  std::vector<double> tv[6];             // u, p, v, pd, a, pdd as written
  for (int k = 0; k < nen; k++) {
    Vector d(ndfN), v(ndfN), a(ndfN);
    for (int c = 0; c < ndfN; c++) {
      const int j = ndfN * k + c;
      const double s = std::sin(1.7 * j + 0.3);
      const int kind = (c < ndm) ? 0 : 1;
      d(c) = sc[kind] * s;
      v(c) = sc[2 + kind] * s;
      a(c) = sc[4 + kind] * s;
      tv[kind].push_back(d(c));
      tv[2 + kind].push_back(v(c));
      tv[4 + kind].push_back(a(c));
    }
    nodes[k]->setTrialDisp(d);
    nodes[k]->setTrialVel(v);
    nodes[k]->setTrialAccel(a);
  }
  auto restore = [&]() {
    for (int k = 0; k < nen; k++) {
      nodes[k]->setTrialDisp(sd[k]);
      nodes[k]->setTrialVel(sv[k]);
      nodes[k]->setTrialAccel(sa[k]);
    }
    mat->revertToLastCommit();
    needUpdate = true;
    trialValid = false;
  };
  needUpdate = true;
  trialValid = false;
  if (this->update() < 0) { restore(); return -1; }
  const Matrix K(this->getTangentStiff());
  const Matrix C(this->getDamp());
  const Matrix M(this->getMass());
  const Vector Rs(this->getResistingForce());
  const Vector R(this->getResistingForceIncInertia());
  const ladruno_mpm::Particles &PP = h->mpmParticles();
  os << std::setprecision(17);
  auto vec = [&](const std::vector<double> &x) {
    os << "[";
    for (size_t i = 0; i < x.size(); i++) os << (i ? ", " : "") << x[i];
    os << "]";
  };
  auto rows = [&](const std::vector<double> &x, int w) {   // x as rows of width w
    os << "[";
    for (int r = 0; w > 0 && r < (int)x.size() / w; r++) {
      os << (r ? ", " : "") << "[";
      for (int c = 0; c < w; c++) os << (c ? ", " : "") << x[w * r + c];
      os << "]";
    }
    os << "]";
  };
  auto mat2 = [&](const Matrix &A) {
    os << "[";
    for (int r = 0; r < nf; r++) {
      os << (r ? ", " : "") << "[";
      for (int c = 0; c < nf; c++) os << (c ? ", " : "") << A(r, c);
      os << "]";
    }
    os << "]";
  };
  auto vec2 = [&](const Vector &x) {
    os << "[";
    for (int i = 0; i < nf; i++) os << (i ? ", " : "") << x(i);
    os << "]";
  };
  std::vector<double> crd;
  for (int k = 0; k < nen; k++) {
    const Vector &x = nodes[k]->getCrds();
    for (int a = 0; a < ndm; a++) crd.push_back(a < x.Size() ? x(a) : 0.0);
  }
  const std::vector<double> &SpO = tpOn ? Sp : S, &dSpO = tpOn ? dSp : dS;
  os << "{\"particle\": " << p << ", \"id\": " << id << ", \"ndm\": " << ndm << ", \"ndf\": "
     << ndfN << ", \"nen\": " << nen << ", \"lambda\": " << lambda << ", \"V\": " << PP.V[p]
     << ", \"m\": " << mp << ", \"alpha_st\": " << stabA << ", \"two_phase\": " << (tpOn ? 1 : 0);
  {
    double gpi[3];                       // fix 2: the known affine term at the particle (0 = none)
    os << ", \"pi\": " << piNow(gpi);
  }
  os << ",\n \"nodes\": ";
  os << "[";
  for (int k = 0; k < nen; k++) os << (k ? ", " : "") << connect(k);
  os << "],\n \"coords\": ";
  rows(crd, ndm);
  os << ",\n \"S\": "; vec(S);
  os << ",\n \"dS\": "; rows(dS, ndm);
  os << ",\n \"Sp\": "; vec(SpO);
  os << ",\n \"dSp\": "; rows(dSpO, ndm);
  os << ",\n \"trial\": {\"u\": "; rows(tv[0], ndm);
  os << ", \"p\": "; vec(tv[1]);
  os << ", \"v\": "; rows(tv[2], ndm);
  os << ", \"pd\": "; vec(tv[3]);
  os << ", \"a\": "; rows(tv[4], ndm);
  os << ", \"pdd\": "; vec(tv[5]);
  os << "},\n \"K\": "; mat2(K);
  os << ",\n \"C\": "; mat2(C);
  os << ",\n \"M\": "; mat2(M);
  os << ",\n \"R\": "; vec2(R);
  os << ",\n \"R_static\": "; vec2(Rs);
  os << "}\n";
  restore();
  return 0;
}

int
LadrunoMPMPoint::setTwoPhase(int ndfNode, const ladruno_mpm::TwoPhaseProps *tp, double h,
                             const std::vector<double> &Spk, const std::vector<double> &dSpk,
                             const std::vector<double> &pRefk, const double *aff, double piStart)
{
  if (ndfNode != ndm && ndfNode != ndm + 1) return -1;
  ndfN = ndfNode;
  tpOn = (tp != 0) && ndfN == ndm + 1;
  if (nen * ndfN > MAX_DOF) {
    opserr << "LadrunoMPMPoint " << this->getTag() << " -- " << nen
           << " support nodes is out of range for a two-phase element (MPM host bug)\n";
    return -1;
  }
  resize();
  delete Ki;
  Ki = 0;
  index = -1;
  if (!tpOn) return 0;
  if ((int)Spk.size() != nen || (int)dSpk.size() != nen * ndm || (int)pRefk.size() != 3 * nen) {
    opserr << "LadrunoMPMPoint " << this->getTag() << " -- two-phase weights of the wrong size"
           << " (MPM host bug)\n";
    return -1;
  }
  Sp = Spk;
  dSp = dSpk;
  pRef = pRefk;
  tpAlpha = tp->alpha;
  tpRhoF = tp->rhoF;
  for (int i = 0; i < 3; i++) tpK[i] = tp->kbar[i];
  affOn = false;
  piStart_ = piStart;
  if (aff != 0)
    for (int k = 0; k < ndm * (ndm + 1); k++) {
      affC[k] = aff[k];
      if (aff[k] != 0.0) affOn = true;
    }
  LadrunoMPMPointHost *hh = host();
  // V_p: the step-start volume (V0 under `small`, J_n V0 under hypo; interface §1)
  const double Vp = (hh != 0) ? hh->mpmParticles().V[p] : 0.0;
  // the u-row gradient of Q is component i's own (B-spline folds, P9)
  std::vector<double> dNu(nen * ndm);
  for (int k = 0; k < nen; k++)
    for (int i = 0; i < ndm; i++)
      dNu[ndm * k + i] = Sc.empty() ? dS[ndm * k + i] : dSc[(ndm * k + i) * ndm + i];
  Qb.assign(nen * ndm * nen, 0.0);
  Hb.assign(nen * nen, 0.0);
  Sb.assign(nen * nen, 0.0);
  Htb.assign(nen * nen, 0.0);
  Gs.assign(nen * ndm, 0.0);
  ladruno_up::addQ(dNu.data(), Sp.data(), nen, nen, ndm, tp->alpha, Vp, Qb.data(), nen);
  ladruno_up::addH(dSp.data(), nen, ndm, tp->kbar, Vp, Hb.data(), nen);
  ladruno_up::addS(Sp.data(), nen, ladruno_up::oneOverQbar(tp->n, tp->Kf, tp->alpha, tp->Ks), Vp,
                   Sb.data(), nen);
  ladruno_up::addHtilde(dSp.data(), nen, ndm, 1.0, Vp, Htb.data(), nen);
  // G[b, i] = dSp_{b,i} kbar_i rhoF Vp: f_seep[b] = sum_i G[b, i] drive_i (addFseep, unit drives)
  for (int i = 0; i < ndm; i++) {
    double unit[3] = {0.0, 0.0, 0.0};
    unit[i] = 1.0;
    std::vector<double> f(nen, 0.0);
    ladruno_up::addFseep(dSp.data(), nen, ndm, tp->kbar, tp->rhoF, unit, Vp, f.data());
    for (int b = 0; b < nen; b++) Gs[b * ndm + i] = f[b];
  }
  // alpha_st (P1 / §1): off, manual, or auto = alpha0 h^2 / (K_d + 4 G_d / 3) from the
  // material's initial tangent, as LadrunoUP::stabAlphaValue reads it (LadrunoUP.cpp:749-772)
  stabA = 0.0;
  if (tp->stabMode == 2) {
    stabA = tp->stabValue;
  } else if (tp->stabMode == 1 && mat != 0) {
    const Matrix &D0 = mat->getInitialTangent();
    const double Gd = (ndm == 2) ? D0(2, 2) : D0(3, 3);
    const double lam = D0(0, 1);
    const double Kd = lam + 2.0 * Gd / 3.0;
    if (Kd + 4.0 * Gd / 3.0 > 0.0 && h > 0.0) stabA = ladruno_up::stabAlphaAuto(h, Kd, Gd, tp->stabValue);
  }
  return 0;
}

void
LadrunoMPMPoint::zeroLoad(void)
{
}

int
LadrunoMPMPoint::addLoad(ElementalLoad *, double)
{
  opserr << "WARNING LadrunoMPMPoint " << this->getTag()
         << " -- element loads are refused on MPM particles; use `mpm pforce` or `mpm gravity`\n";
  return -1;
}

// ---------------------------------------------------------------------------
//  Restart (interface §6)
// ---------------------------------------------------------------------------
int
LadrunoMPMPoint::sendSelf(int commitTag, Channel &theChannel)
{
  const int dataTag = this->getDbTag();
  if (mat == 0) {
    opserr << "LadrunoMPMPoint::sendSelf -- element " << this->getTag() << " has no material\n";
    return -1;
  }
  int matDbTag = mat->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0) mat->setDbTag(matDbTag);
  }
  // FileDatastore files by (type, size, dbTag, commitTag): the connectivity and
  // the Rayleigh factors go under a second dbTag so no size can coincide with
  // the 10-entry header (LEGACY quirk "FileDatastore silently CLOBBERS").
  if (payloadDbTag == 0) payloadDbTag = theChannel.getDbTag();
  static ID hdr(10);
  hdr(0) = this->getTag();
  hdr(1) = p;
  hdr(2) = hostPatternTag;
  hdr(3) = ndm;
  hdr(4) = nen;
  hdr(5) = mat->getClassTag();
  hdr(6) = matDbTag;
  hdr(7) = payloadDbTag;
  hdr(8) = ndfN;                       // WP-232l: DOFs per node (ndm, or ndm + 1)
  hdr(9) = 0;
  if (theChannel.sendID(dataTag, commitTag, hdr) < 0) {
    opserr << "LadrunoMPMPoint::sendSelf -- failed to send the header\n";
    return -1;
  }
  if (nen > 0 && theChannel.sendID(payloadDbTag, commitTag, connect) < 0) {
    opserr << "LadrunoMPMPoint::sendSelf -- failed to send the connectivity\n";
    return -1;
  }
  static Vector rv(5);
  rv(0) = alphaM; rv(1) = betaK; rv(2) = betaK0; rv(3) = betaKc; rv(4) = mp;
  if (theChannel.sendVector(payloadDbTag, commitTag, rv) < 0) {
    opserr << "LadrunoMPMPoint::sendSelf -- failed to send the Rayleigh factors\n";
    return -1;
  }
  if (mat->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoMPMPoint::sendSelf -- the material failed to send itself\n";
    return -1;
  }
  return 0;
}

int
LadrunoMPMPoint::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  const int dataTag = this->getDbTag();
  static ID hdr(10);
  if (theChannel.recvID(dataTag, commitTag, hdr) < 0) {
    opserr << "LadrunoMPMPoint::recvSelf -- failed to receive the header\n";
    return -1;
  }
  this->setTag(hdr(0));
  p = hdr(1);
  hostPatternTag = hdr(2);
  ndm = hdr(3);
  nv = (ndm == 3) ? 6 : 3;
  nen = hdr(4);
  const int matClass = hdr(5), matDbTag = hdr(6);
  payloadDbTag = hdr(7);
  ndfN = (hdr(8) == ndm + 1) ? ndm + 1 : ndm;   // WP-232l (0 in an older database: ndm)
  tpOn = false;                                  // the next `mpm analyze` re-points and rebuilds
  connect = ID(nen > 0 ? nen : 1);
  if (nen > 0 && theChannel.recvID(payloadDbTag, commitTag, connect) < 0) {
    opserr << "LadrunoMPMPoint::recvSelf -- failed to receive the connectivity\n";
    return -1;
  }
  if (nen == 0) connect = ID(0);
  static Vector rv(5);
  if (theChannel.recvVector(payloadDbTag, commitTag, rv) < 0) {
    opserr << "LadrunoMPMPoint::recvSelf -- failed to receive the Rayleigh factors\n";
    return -1;
  }
  alphaM = rv(0); betaK = rv(1); betaK0 = rv(2); betaKc = rv(3); mp = rv(4);
  if (mat == 0 || mat->getClassTag() != matClass) {
    delete mat;
    mat = theBroker.getNewNDMaterial(matClass);
    if (mat == 0) {
      opserr << "LadrunoMPMPoint::recvSelf -- the broker has no NDMaterial " << matClass << "\n";
      return -1;
    }
  }
  mat->setDbTag(matDbTag);
  if (mat->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "LadrunoMPMPoint::recvSelf -- the material failed to receive itself\n";
    return -1;
  }
  // every cache is dropped; the next `mpm analyze` re-points the element
  S.assign(nen, 0.0);
  ghostK.assign(nen, 0.0);
  resize();
  delete Ki;
  Ki = 0;
  index = -1;
  trialValid = false;
  needUpdate = true;
  if (this->getDomain() != 0) {            // into a live element
    refetchNodes();
    st = COMMITTED;
  } else {
    st = RESTORED;                         // through the broker, before addElement
  }
  return 0;
}

void
LadrunoMPMPoint::Print(OPS_Stream &s, int)
{
  s << "LadrunoMPMPoint " << this->getTag() << " (MPM particle index " << p << ", host pattern "
    << hostPatternTag << ", ndm " << ndm << ", " << nen << " support nodes, state "
    << (st == ARMED ? "armed" : (st == COMMITTED ? "committed" : "restored")) << ")\n";
}

Response *
LadrunoMPMPoint::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  return this->Element::setResponse(argv, argc, output);
}

int
LadrunoMPMPoint::getResponse(int responseID, Information &eleInfo)
{
  return this->Element::getResponse(responseID, eleInfo);
}

int
LadrunoMPMPoint::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1 || mat == 0) return -1;
  // WP-232k (tests/mpm_soil_interface.md D9 (c)): MaterialStageParameter stops at the first
  // element that answers (MaterialStageParameter.cpp:79), so a particle that answered alone
  // would leave every sibling to SANISAND's lazy flip, which skips Elastic2Plastic. A soil
  // stage request is handed to the host, which registers itself and walks every particle.
  if (argc >= 2 && strcmp(argv[0], "updateMaterialStage") == 0 &&
      param.getClassTag() == PARAMETER_TAG_MaterialStageParameter) {
    LadrunoMPMPointHost *h = host();
    if (h != 0) {
      const int r = h->mpmStageParameter(argv, argc, param);
      if (r != -2) return r;
    }
  }
  // `parameter ... element <tag> material 1 <args>` or a forall broadcast
  if (LadrunoShell::isMaterialPointToken(argv[0]))
    return LadrunoShell::forwardToMaterialPoint(&mat, 1, true, argv, argc, param);
  return LadrunoShell::forwardToMaterials(&mat, 1, argv, argc, param);
}
