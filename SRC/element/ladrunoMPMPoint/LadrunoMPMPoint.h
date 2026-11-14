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
// position (2 or 3 per axis: 4-9 nodes in 2D, 8-27 in 3D; linear hat 4 / 8);
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
// The element overrides no threading opt-in: it stays off every WP-107 / WP-192
// path (its scratch is static, shared by size).

#ifndef LadrunoMPMPoint_h
#define LadrunoMPMPoint_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

#include <memory>
#include <vector>

class NDMaterial;
class Node;
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
  // the debug hook `mpm debug failAttempt`: this attempt fails before commit
  virtual bool   mpmFailThisAttempt(void) = 0;
  // the batched G2P: every ARMED particle element counts in at commit
  virtual void   mpmParticleCommitted(int p) = 0;
  // the element is being destroyed (its material with it)
  virtual void   mpmElementGone(int p) = 0;
  // ghost regularization spring of particle p (interface D2.3); 0 = none
  virtual double mpmGhostSpring(int p) = 0;
};

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
  int  repoint(const ID &tags, const std::vector<double> &S, const std::vector<double> &dS,
              const std::vector<char> &inactive, const std::vector<double> &vHat,
              const std::vector<double> &aHat);
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

  // ---- Element ----------------------------------------------------------------
  int getNumExternalNodes(void) const { return nen; }
  const ID &getExternalNodes(void) { return connect; }
  Node **getNodePtrs(void) { return nodes.empty() ? 0 : nodes.data(); }
  int getNumDOF(void) { return nen * ndm; }
  void setDomain(Domain *theDomain);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

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

  ID connect;
  std::vector<Node *> nodes;
  std::vector<double> S, dS, vHat, aHat, ghostK;   // per support node (ghostK: spring, 0 = none)

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

  LadrunoMPMPointHost *host(void) const { return (link && link->host) ? link->host : 0; }
  void resize(void);
  int  computeTrial(void);
  void assembleForce(void);
  void tangent(Matrix &K, bool initial);
  void refetchNodes(void);
};

#endif
