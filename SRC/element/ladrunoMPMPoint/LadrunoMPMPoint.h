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
// LadrunoMPMPoint — the mode-A MPM particle element (WP-232f, ADR 232 §5.3;
// frozen interface tests/mpm_implicit_interface.md §1). ELE_TAG 33026.
//
// One element per particle: one integration point, one owned NDMaterial copy.
// Its nodes are the active GIMP support of the particle at the STEP-START
// position (2 or 3 per axis: 4-9 nodes in 2D, 8-27 in 3D; linear hat 4 / 8;
// WP-232i: bspline2 <= 3 per axis, bspline3 <= 4, so <= 64 nodes in 3D);
// every grid node has ndf = ndm. Not user-creatable: the MPM host
// (LadrunoMPMPattern) creates it, re-points it in place every step and owns the
// particle's kinematic state (ladruno_mpm::Particles); the element writes that
// state at G2P.
//
// STATES. ARMED: re-pointed for this step, update() computes the trial state.
// COMMITTED: the host's batched G2P ran; update() refuses loudly (a step not
// prepared by `mpm analyze`), the resisting forces return the converged cache.
// RESTORED: built by the broker; the one update() Domain::addElement makes
// answers 0, then COMMITTED. A host that is gone leaves the element ORPHAN
// (the shared link's host pointer is null): update() refuses loudly.
//
// KINEMATICS (interface D1.5). Step reference = the reset grid (u = 0 at every
// node after P2G). F_rel = I + Σ u_a ⊗ ∇S_a, F_1 = F_rel·F_n, then the
// incremental hypo-kernel entry (Hughes–Winget, unrotated frame); the force and
// tangent are LadrunoBrick's hypo assembly at one point with g_a = ∇S_a·F_rel⁻¹
// and V = J_1 V_0 (D1.6). Mass: lumped M_aa = S_a m_p I (D1.7). Rayleigh:
// alphaM and betaK; betaK0 / betaKc are refused by the host (D1.8).
//
// G2P (D1.9, ADR §5.3 amendment 1): commitState() commits the material and
// counts in with the host; the host's N-th count runs g2p() on every particle
// unless a commit refusal is pending — then nothing moves.
//
// ANTI-LOCKING (WP-232j, tests/mpm_antilocking_interface.md J7-J12): with the
// host's `mpm antiLocking fbar`, repoint() also receives the slots of the 2^d
// corners of the particle's cell (always inside its own support, J3). The trial
// increment becomes dFbar = (det dF_c / det dF)^(1/d) dF (J1), the residual
// weight the carried volume det(dFbar F_n) V_0, the internal force the one
// work-conjugate to that kinematics, f_a = w [sigma g_a + pbar (ghat_a - g_a)]
// (interface §11, C1: Sum f.v = w sigma : lbar), and the tangent its consistent
// linearisation. No corners = the shipped arithmetic, bit for bit.
//
// The element overrides no threading opt-in: it stays off every WP-107 / WP-192
// path (its scratch is static, shared by size).

#ifndef LadrunoMPMPoint_h
#define LadrunoMPMPoint_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

#include <iosfwd>
#include <memory>
#include <vector>

class NDMaterial;
class Node;
class Parameter;
class Channel;
class FEM_ObjectBroker;

namespace ladruno_mpm { struct Particles; }

// The host's side of the particle elements (implemented by LadrunoMPMPattern).
class LadrunoMPMPointHost
{
 public:
  virtual ~LadrunoMPMPointHost() {}
  virtual ladruno_mpm::Particles &mpmParticles(void) = 0;
  // body force on particle p: λ (m_p b + f_p), ndm entries (live, the host's λ)
  virtual const double *mpmBodyForce(int p) = 0;
  // the trial centre x stays admissible (interface A7)
  virtual bool   mpmTrialAdmissible(int p, const double *x) = 0;
  // the debug hook `mpm debug failAttempt`: 0 no, 1 fail before the trial, 2 after it
  virtual int    mpmFailThisAttempt(void) = 0;
  // inside the host's `mpm analyze` / tangentCheck: the only window update() may run in
  virtual bool   mpmInStep(void) = 0;
  // the batched G2P: every ARMED particle element counts in at commit
  virtual void   mpmParticleCommitted(int p) = 0;
  // the element is being destroyed (its material with it)
  virtual void   mpmElementGone(int p) = 0;
  // ghost regularization spring of particle p (interface D2.3); 0 = none
  virtual double mpmGhostSpring(int p) = 0;
  // WP-232k (tests/mpm_soil_interface.md D9 (c)): `updateMaterialStage` reached this
  // particle element through MaterialStageParameter. The host registers ITSELF and walks
  // every particle (0), or answers -2: not a soil stage request, the element forwards as before.
  virtual int    mpmStageParameter(const char **, int, Parameter &) { return -2; }
  // WP-232k (D14): particle p's material refused the trial strain (the host names it)
  virtual void   mpmTrialRefused(int) {}
  // WP-232l (tests/mpm_twophase_interface.md): `mpm kinematics small` (P6); the pore-pressure
  // G2P blend alpha_p (P5); particle p's seepage drive lambda b_f (ndm values, P1)
  virtual bool   mpmSmallKinematics(void) { return false; }
  virtual double mpmPressureBlend(void) { return 1.0; }
  virtual const double *mpmFluidDrive(int) { return 0; }
};

namespace ladruno_mpm { struct TwoPhaseProps; }

// Shared liveness token: the host nulls `host` in its destructor.
struct LadrunoMPMHostLink {
  LadrunoMPMPointHost *host = nullptr;
};

class LadrunoMPMPoint : public Element
{
 public:
  enum State { ARMED = 0, COMMITTED = 1, RESTORED = 2 };

  // particle index p, the host pattern's tag, ndm, an OWNED material copy
  LadrunoMPMPoint(int tag, int particle, int hostTag, int ndm, NDMaterial *ownedMat);
  LadrunoMPMPoint();
  ~LadrunoMPMPoint();

  const char *getClassType(void) const { return "LadrunoMPMPoint"; }

  // ---- the host's interface ---------------------------------------------------
  void setLink(const std::shared_ptr<LadrunoMPMHostLink> &l) { link = l; }
  int  particle(void) const { return p; }
  int  hostTag(void) const { return hostPatternTag; }
  NDMaterial *material(void) { return mat; }
  // Hand the material back to the host and forget it (the element is about to
  // be removed and deleted by the host, `mpm reset`).
  NDMaterial *releaseMaterial(void);
  // Re-point in place (interface D1.4): node tags, S_a, ∇S_a (ndm per node,
  // spatial at the step-start position), inactive flags (ghost mode), the
  // projected step-start nodal velocity / acceleration (ndm per node).
  // Refetches Node* when the element is in a Domain. Leaves the element ARMED.
  // WP-232i (tests/mpm_bspline_interface.md B8): Sc, dSc are the per-component
  // S'_{k,a} (ndm per node) and grad S'_{k,a} (ndm x ndm per node, component-major)
  // of a support folded by images with a sign -1; empty for every other particle,
  // which keeps the scalar path bit for bit. S then holds the unsigned weight W
  // (the lumped mass), dS its gradient.
  int  repoint(const ID &tags, const std::vector<double> &S, const std::vector<double> &dS,
              const std::vector<char> &inactive, const std::vector<double> &vHat,
              const std::vector<double> &aHat, const std::vector<double> &Sc,
              const std::vector<double> &dSc, const int *corners = 0, int nCorners = 0,
              double cellH = 0.0);
  // WP-232j: 0 not averaged (no corners), 1 averaged, 2 fallback (a corner absent)
  int  antiLockingPath(void) const { return alMode == 0 ? 0 : (corner[0] < 0 ? 2 : 1); }
  // `mpm debug tangentCheck -antiLocking off|fbar`: 1 suppresses the F-bar for the check
  void setAntiLockingSuppressed(bool s) { alSuppressed = s; needUpdate = true; trialValid = false; }
  // `mpm debug affineCheck` (L1): the last trial's dF (actual) and dFbar, and its deps
  void antiLockingTrial(double &rel, double &depsMax) const;
  // max |sigHat_trial - sigHat_n| of the last trial (affineCheck: the spurious stress, L1b)
  double trialStressChange(void) const;
  // WP-232j failure diagnosis (`mpm debug failDump`): when set, every trial records the
  // material's tangent norm against its initial one (an apex / perfectly plastic point
  // returns a regularised, near-zero tangent).
  static bool diagnose;
  double diagTangentRatio(void) const { return diagTan; }
  double diagVolumeRatio(void) const { return diagJr; }       // det dF_c / det dF (1 off)
  const double *trialStress(void) const { return sigT; }
  const Vector *lastResidual(void) const { return P; }
  bool connectivityEquals(const ID &tags) const;
  State state(void) const { return st; }
  // The host's batched G2P (interface D1.9): writes particle p of P from the
  // converged trial state; SCALE mutant on the position increment. flip = α_F;
  // quasiStatic leaves v and a alone. The element is COMMITTED afterwards.
  void g2p(ladruno_mpm::Particles &P, double flip, bool quasiStatic);
  // `mpm debug tangentCheck`: ||K - K_FD||_F / ||K||_F at the trial displacement
  // scale * sin(1.7 j + 0.3) (central differences of step eps); state restored.
  int  tangentCheck(double scale, double eps, double &rel, double &normK);
  // Internal force V σ·g_a at node slot k (ndm entries) of the last trial /
  // converged state — the cases' cond_force (interface §5).
  const double *internalForceAt(int k) const { return &fInt[ndm * k]; }
  bool refusedRayleigh(void) const { return betaK0 != 0.0 || betaKc != 0.0; }
  void setRayleigh(double am, double bk);
  // WP-232l (tests/mpm_twophase_interface.md P3, P4): after repoint, in a two-phase epoch.
  // ndfNode = ndm + 1 (every grid node carries p last); tp = the body's fluid constants, or 0
  // for a single-phase particle (zero p rows); h = the grid spacing (the stabilization
  // length); Sp, dSp = the p weights (nen, nen x ndm); pRef = the unconditioned projected
  // p, pd, pdd per support node (3 x nen), the G2P references. The u-p blocks are built here
  // from the step-start weights and volume.
  // WP-232l fix 2 (interface 7b): aff = the affine-image coefficients a_d, grad a_d per axis
  // (ndm x (ndm + 1), 0 = none); piStart = the known pressure term pi at the projection
  int  setTwoPhase(int ndfNode, const ladruno_mpm::TwoPhaseProps *tp, double h,
                   const std::vector<double> &Sp, const std::vector<double> &dSp,
                   const std::vector<double> &pRef, const double *aff = 0, double piStart = 0.0);
  bool twoPhase(void) const { return tpOn; }
  // `mpm debug particleBlocks` (§7a Q1/Q2): set the pattern trial state, form K, C, M, R and
  // R_static, write the JSON object, restore. 0, or -1 (message printed).
  int  writeBlocks(const double scales[6], std::ostream &os, double lambda, long long id);

  // ---- Element ----------------------------------------------------------------
  int getNumExternalNodes(void) const { return nen; }
  const ID &getExternalNodes(void) { return connect; }
  Node **getNodePtrs(void) { return nodes.empty() ? 0 : nodes.data(); }
  int getNumDOF(void) { return nen * ndfN; }      // WP-232l: ndfN = ndm (+ 1 in a two-phase epoch)
  void setDomain(Domain *theDomain);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);
  // WP-232l: [C_Ray(solid), 0; Q^T, S*] in a two-phase epoch (ADR-71 §3.2 obligation 2);
  // otherwise Element::getDamp, bit for bit
  const Matrix &getDamp(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &eleInfo);
  int setParameter(const char **argv, int argc, Parameter &param);

 private:
  int  p, hostPatternTag, ndm, nv, nen;
  NDMaterial *mat;
  std::shared_ptr<LadrunoMPMHostLink> link;
  State st;
  bool  needUpdate;      // re-pointed or reverted, no update() since
  bool  trialValid;      // the trial cache holds a successful update()
  // WP-232j (J8): alMode 0 off; 1 corners given (corner[0] < 0 = fallback, J5)
  int    alMode, nCorner, corner[8];
  bool   alSuppressed;
  double cellH, Jc, alRel;   // alRel: |dFbar - dF|_max / |dF - I|_max of the last trial
  double Fci[9];          // dF_c^-1 of the last trial (the tangent's g-hat)
  double depsT[6];        // the last trial's unrotated strain increment (affineCheck)
  double gHat[24];        // C1: ghat_k = gradN_k(x_c) . dF_c^-1 per corner (the conjugate force)
  double diagTan, diagJr; // failure diagnosis (diagnose)
  bool   alActive(void) const { return alMode != 0 && !alSuppressed && corner[0] >= 0; }

  ID connect;
  std::vector<Node *> nodes;
  std::vector<double> S, dS, vHat, aHat, ghostK;   // per support node (ghostK: spring, 0 = none)
  // WP-232i B8: per-component S', grad S' and spatial g' (empty = the scalar path)
  std::vector<double> Sc, dSc, gc;

  // trial cache (interface D1.5 / D1.6)
  double Frel[9], F1[9], R1[9], J1;
  double epsT[6], sigHatT[6], sigT[6];
  std::vector<double> g;          // spatial gradients g_a, ndm per node
  std::vector<double> fInt;       // V σ g_a
  std::vector<double> uT;         // trial nodal displacement (ghost springs)
  double bf[3];                   // the body force of the last evaluation
  double vol;                     // V = J_1 V_0
  double mp;                      // particle mass (lumped mass, the restart)
  int    payloadDbTag;            // second dbTag: connectivity + factors (no size clash)

  Vector *P;                      // resisting force, sized numDOF
  Vector *Pinc;                   // incl. inertia and Rayleigh
  Matrix *Ki;                     // initial stiffness of this support

  // WP-232l (two-phase, interface §1 P3-P6)
  int    ndfN;                    // DOFs per node: ndm, or ndm + 1 in a two-phase epoch
  bool   tpOn;                    // this particle carries the u-p blocks
  bool   smallK;                  // `mpm kinematics small` (set at repoint)
  double stabA;                   // alpha_st
  std::vector<double> Sp, dSp, pRef;      // p weights, gradients; p, pd, pdd references
  std::vector<double> Qb, Hb, Sb, Htb, Gs; // Q[(ndm k + i) nen + b], H, S, H~ [b nen + c], G[b ndm + i]
  std::vector<double> pT;                  // trial nodal p of the last assembly
  double tpAlpha, tpRhoF, tpK[3];          // the body's alpha, rhoF, kbar (fix 2's known term)
  double affC[12];                         // fix 2: a_d and grad a_d, ndm x (ndm + 1)
  bool   affOn;
  double piStart_;
  double piNow(double *grad) const;        // fix 2: pi and grad pi at the current drive
  // WP-232l (review of record): the axes whose normal component a plate carries; their lumped
  // mass is the signed, translation-consistent weight S'_a (B-spline affine plate image)
  unsigned plateMask_ = 0u;
 public:
  void   setPlateMask(unsigned m) { plateMask_ = m; }
 private:
  double massW(int k, int a) const
  { return ((plateMask_ & (1u << a)) && !Sc.empty()) ? Sc[ndm * k + a] : S[k]; }
  int  dofU(int k, int a) const { return ndfN * k + a; }
  int  dofP(int k) const { return ndfN * k + ndm; }
  void compactTangent(Matrix &Kc, bool initial);
  void assembleForceUP(void);

  LadrunoMPMPointHost *host(void) const { return (link && link->host) ? link->host : 0; }
  void resize(void);
  int  computeTrial(void);
  void assembleForce(void);
  void tangent(Matrix &K, bool initial);
  void tangentPerComponent(Matrix &K, bool initial, const double *Rm, double volm,
                           const double *sg, const double *Dc);
  void refetchNodes(void);
  void addFbarTangent(Matrix &K, bool initial, double volm, const double *sg, const double *Dc);
  void fillB(int a, bool initial, double B[6][3]) const;
  void addConjugateForce(const double *sg);
  void cornerGradients(double gN[8][3]) const;
};

#endif
