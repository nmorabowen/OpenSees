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

#include <cmath>
#include <string.h>

namespace {

// Shared scratch, keyed by the number of DOFs (at most 27 nodes x 3). The
// element is serial (no threading opt-in), so one scratch per size is safe.
const int MAX_DOF = 81;

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

LadrunoMPMPoint::LadrunoMPMPoint(int tag, int particle, int hostTag, int nd, NDMaterial *ownedMat)
  : Element(tag, ELE_TAG_LadrunoMPMPoint), p(particle), hostPatternTag(hostTag), ndm(nd),
    nv(nd == 3 ? 6 : 3), nen(0), mat(ownedMat), st(ARMED), needUpdate(true), trialValid(false),
    connect(0), J1(1.0), vol(0.0), mp(0.0), payloadDbTag(0), P(0), Pinc(0), Ki(0)
{
  for (int k = 0; k < 9; k++) Frel[k] = F1[k] = R1[k] = (k % 4 == 0) ? 1.0 : 0.0;
  for (int k = 0; k < 6; k++) epsT[k] = sigHatT[k] = sigT[k] = 0.0;
  bf[0] = bf[1] = bf[2] = 0.0;
}

LadrunoMPMPoint::LadrunoMPMPoint()
  : Element(0, ELE_TAG_LadrunoMPMPoint), p(-1), hostPatternTag(0), ndm(2), nv(3), nen(0),
    mat(0), st(RESTORED), needUpdate(true), trialValid(false), connect(0), J1(1.0), vol(0.0),
    mp(0.0), payloadDbTag(0), P(0), Pinc(0), Ki(0)
{
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
  if (am != 0.0 || bk != 0.0)
    this->Element::setRayleighDampingFactors(am, bk, 0.0, 0.0);
}

void
LadrunoMPMPoint::resize(void)
{
  const int n = nen * ndm;
  if (P == 0 || P->Size() != n) { delete P; P = new Vector(n); }
  if (Pinc == 0 || Pinc->Size() != n) { delete Pinc; Pinc = new Vector(n); }
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
                         const std::vector<double> &vH, const std::vector<double> &aH)
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
  vHat = vH;
  aHat = aH;
  ghostK.assign(nen, 0.0);
  LadrunoMPMPointHost *h = host();
  if (h != 0) {
    for (int k = 0; k < nen; k++)
      if (k < (int)inactive.size() && inactive[k]) ghostK[k] = h->mpmGhostSpring(p);
    mp = h->mpmParticles().m[p];
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
      for (int k = 0; k < nen; k++) s += uT[nd * k + a] * dS[nd * k + c];
      Frel[nd * a + c] = s;
    }
  const double *Fn = &PP.F[nn * p];
  double F1loc[9], R1loc[9], deps[6], J;
  for (int a = 0; a < nd; a++)
    for (int c = 0; c < nd; c++) {
      double s = 0.0;
      for (int e = 0; e < nd; e++) s += Frel[nd * a + e] * Fn[nd * e + c];
      F1loc[nd * a + c] = s;
    }
  const int kst = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1loc, deps, R1loc, &J)
                            : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1loc, deps, R1loc, &J);
  if (kst != 0) return -2;
  // the trial centre stays in the real grid (interface A7)
  double x[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < nd; a++) {
    double s = PP.x[nd * p + a];
    for (int k = 0; k < nen; k++) s += S[k] * uT[nd * k + a];
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
  if (mat->setTrialStrain(e) < 0) return -4;
  const Vector &sm = mat->getStress();
  for (int k = 0; k < nv; k++) sigHatT[k] = PP.sig0[nv * p + k] + sm(k);
  if (nd == 3) ladruno_hypo::pushStress6(R1loc, sigHatT, sigT);
  else ladruno_hypo::pushStressPlaneStrain(R1loc, sigHatT, sigT);
  for (int k = 0; k < nn; k++) { F1[k] = F1loc[k]; R1[k] = R1loc[k]; }
  J1 = J;
  vol = J * PP.V0[p];
  // spatial gradients g_a = ∇S_a · F_rel⁻¹
  double Fi[9];
  if (nd == 3) {
    if (!(ladruno_hypo::inv3(Frel, Fi) > 0.0)) return -2;
  } else {
    const double det = Frel[0] * Frel[3] - Frel[1] * Frel[2];
    if (!(det > 0.0)) return -2;
    Fi[0] = Frel[3] / det;  Fi[1] = -Frel[1] / det;
    Fi[2] = -Frel[2] / det; Fi[3] = Frel[0] / det;
  }
  double sg[9];
  voigtToTensor(nd, sigT, sg);
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
  trialValid = true;
  needUpdate = false;
  return 0;
}

int
LadrunoMPMPoint::update(void)
{
  if (st == RESTORED) {        // Domain::addElement's update() after a restore
    st = COMMITTED;
    return 0;
  }
  LadrunoMPMPointHost *h = host();
  if (h == 0 || st != ARMED) {
    static int nWarn = 0;
    if (nWarn++ < 3)
      opserr << "ERROR LadrunoMPMPoint " << this->getTag() << " -- "
             << (h == 0 ? "its MPM host is gone (remove loadPattern)"
                        : "a step not prepared by `mpm analyze`")
             << ": mode-A particle elements move only through `mpm analyze` (ADR 232 §5.3)\n";
    return -1;
  }
  if (h->mpmFailThisAttempt()) {       // `mpm debug failAttempt`
    trialValid = false;
    return -1;
  }
  const int rc = computeTrial();
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
    for (int k = 0; k < nen; k++) dx += S[k] * nodes[k]->getTrialDisp()(a);
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
        dv += S[k] * (v - vHat[nd * k + a]);
        vpic += S[k] * v;
        da += S[k] * (ac - aHat[nd * k + a]);
        apic += S[k] * ac;
      }
      double &vp = PP.v[nd * p + a];
      double &ap = PP.a[nd * p + a];
      vp = flip * (vp + dv) + (1.0 - flip) * vpic;
      ap = flip * (ap + da) + (1.0 - flip) * apic;
    }
  }
  st = COMMITTED;
}

// ---------------------------------------------------------------------------
//  Force, tangent, mass (interface D1.6, D1.7)
// ---------------------------------------------------------------------------
void
LadrunoMPMPoint::assembleForce(void)
{
  LadrunoMPMPointHost *h = host();
  if (h != 0) {
    const double *b = h->mpmBodyForce(p);
    for (int a = 0; a < ndm; a++) bf[a] = b[a];
  }
  for (int k = 0; k < nen; k++)
    for (int a = 0; a < ndm; a++) {
      const int r = ndm * k + a;
      (*P)(r) = fInt[r] - S[k] * bf[a] + ghostK[k] * uT[r];
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
  if (nd == 3) {
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
      // geometric (initial-stress) term: δ_ik g_a·σ·g_b V
      if (!initial) {
        double s = 0.0;
        for (int j = 0; j < nd; j++)
          for (int l = 0; l < nd; l++) s += ga[j] * sg[3 * j + l] * gb[l];
        for (int i = 0; i < nd; i++) K(nd * a + i, nd * b + i) += s * volm;
      }
    }
    for (int i = 0; i < nd; i++) K(nd * a + i, nd * a + i) += ghostK[a];   // ghost spring (D2.3)
  }
}

const Matrix &
LadrunoMPMPoint::getTangentStiff(void)
{
  const int n = nen * ndm;
  Matrix &K = scratchK(n > 0 ? n : 1);
  if (n == 0) return K;
  if (st == ARMED && (needUpdate || !trialValid)) this->update();
  tangent(K, false);
  return K;
}

const Matrix &
LadrunoMPMPoint::getInitialStiff(void)
{
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
  const int n = nen * ndm;
  Matrix &M = scratchM(n > 0 ? n : 1);
  M.Zero();
  for (int k = 0; k < nen; k++)
    for (int a = 0; a < ndm; a++) M(ndm * k + a, ndm * k + a) = S[k] * mp;
  return M;
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
  hdr(8) = 0;
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
  // `parameter ... element <tag> material 1 <args>` or a forall broadcast
  if (LadrunoShell::isMaterialPointToken(argv[0]))
    return LadrunoShell::forwardToMaterialPoint(&mat, 1, true, argv, argc, param);
  return LadrunoShell::forwardToMaterials(&mat, 1, argv, argc, param);
}
