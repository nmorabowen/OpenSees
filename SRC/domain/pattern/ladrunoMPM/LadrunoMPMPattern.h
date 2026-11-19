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
// LadrunoMPMPattern — the host of the explicit MPM engine (WP-232 mode F,
// ADR 232 §5.2). PATTERN_TAG_LadrunoMPM = 33027.
//
// WHY A LoadPattern (ADR 232 §5.2, after ADR-73's LadrunoPorousOverlay): with no
// Domain edit it gets deletion by `wipe` (Domain::clearAll) and by
// `remove loadPattern`; the `parameter <t> loadPattern <tag> ...` target, which
// it forwards to its particles' materials (W12: the stock stage walks reach
// Domain elements only); and its own TimeSeries (gravity ramps). It is NOT an
// element: an element with no DOFs makes FE_Element exit(-1) at the first
// analysis set-up (W15).
//
// WHAT IT OWNS: the pure engine core (LadrunoMPMCore.h — grid, particles, the
// MUSL step) and one NDMaterial copy per particle. Outside `mpm analyze` its
// applyLoad does nothing, so a later FE analysis in the same Domain is unaffected.
//
// ONE STEP (ADR 232 §5.2), driven by analyze():
//   1. refuse an inadmissible particle or Δt > s·Δt_c
//   2. Domain::applyLoad(t + Δt): sets the domain time and ops_Dt; this
//      pattern's applyLoad reads its series factor
//   3–7. the core: P2G, grid update and conditions, MUSL with the kernel entry
//      and ops_TheActiveElement = 0 (W9), commit the particle materials (a
//      declared refusal aborts the run, W17), G2P
//   8. Domain::commit(): committed time, commit tag, recorders — after G2P (W11)
//   9. the profiler's "step" scope and OPS_PROFILE_STEP row
// and Domain::flushRecorders() at the end.
//
// v1 REFUSES: free FE DOFs in the Domain; any material outside the admitted
// table; LadrunoJ2 -autoRegularization (a length consumer); initial stress on a
// non-elastic material (Q2: no yield-consistent seeding in v1); database
// save / restore (sendSelf / recvSelf); OpenSeesSP / OpenSeesMP (the command).
//
// MODE A (WP-232f, ADR 232 §5.3, tests/mpm_implicit_interface.md): `mpm mode
// implicit | quasiStatic` makes every particle a LadrunoMPMPoint element and the
// active support nodes Domain nodes, re-pointed every step; the user's
// Transient / Static analysis solves the step (LadrunoMPMImplicit.cpp). The
// host implements LadrunoMPMPointHost: it owns the particle state, runs the
// batched G2P at the last particle commit, and owns the grid conditions as its
// own SP_Constraints. Database save / restore is admitted in mode A only.

#ifndef LadrunoMPMPattern_h
#define LadrunoMPMPattern_h

#include <LoadPattern.h>
#include "LadrunoMPMCore.h"
#include "LadrunoMPMIO.h"
#include "../../../element/ladrunoMPMPoint/LadrunoMPMPoint.h"

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class NDMaterial;
class TimeSeries;
class Channel;
class FEM_ObjectBroker;
class OPS_Stream;
class Parameter;
class Information;
class StaticAnalysis;
class DirectIntegrationAnalysis;
class StaticIntegrator;
class TransientIntegrator;
class ConstraintHandler;
class LinearSOE;
class SP_Constraint;
class Node;

// The user's analysis, as each engine's bridge sees it (interface §4): classic
// Tcl from its globals (commands.cpp), Python from the `cmds` singleton
// (OpenSeesCommands.cpp). Mode F needs none.
struct LadrunoMPMAnalysisRefs {
  StaticAnalysis            *sa = 0;
  DirectIntegrationAnalysis *ta = 0;
  StaticIntegrator          *si = 0;
  TransientIntegrator       *ti = 0;
  ConstraintHandler         *handler = 0;
  bool                       bridged = false;   // an engine bridge filled it
  // WP-232j (K3 refusal): the user's system, and PARDISO's matType when it is the
  // PARDISO SOE (-1 otherwise; only the bridges compile with _PARDISO)
  LinearSOE                 *soe = 0;
  int                        pardisoMatType = -1;
};

class LadrunoMPMPattern : public LoadPattern, public LadrunoMPMPointHost
{
 public:
  explicit LadrunoMPMPattern(int tag);
  LadrunoMPMPattern();
  ~LadrunoMPMPattern();

  // ---- LoadPattern overrides -----------------------------------------------
  void         applyLoad(double pseudoTime = 0.0);
  void         setTimeSeries(TimeSeries *theSeries);
  double       getLoadFactor(void);
  int          sendSelf(int commitTag, Channel &theChannel);
  int          recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void         Print(OPS_Stream &s, int flag = 0);
  LoadPattern *getCopy(void);
  int          setParameter(const char **argv, int argc, Parameter &param);
  int          updateParameter(int parameterID, Information &info);
  int          activateParameter(int parameterID);

  // ---- the engine, driven by the `mpm` command ------------------------------
  ladruno_mpm::Engine &engine() { return eng; }
  bool hasStepped() const { return stepped; }
  // Incremented by every `mpm reset`: the .ladruno particle group's
  // PARTICLE_EPOCH stage attribute (232a contract 2, phase c2).
  int  particleEpoch() const { return epoch; }

  // One particle with its own material copy (admission checked). Returns 0, or
  // -1 with the message in err.
  int addParticle(int64_t id, int body, int matTag, double rho, double V0,
                  const double *x, const double *lp, const double *v,
                  const double *sig0, std::string &err);
  // Lattice seeding of a body (interface §5). shape 0 box, 1 disc, 2 ball.
  // Returns the count, or -1.
  int seedBody(int shape, const double *geom, int ppc, int body, int matTag,
               double rho, std::string &err);
  // Initial stress on every particle of a body (Voigt, tensor shear). Refused on
  // a non-elastic material (Q2). Returns the count, or -1.
  int setBodyStress(int body, const std::vector<double> &sig, std::string &err);
  // Initial velocity of every particle of a body: uniform v, or (sine) v_comp =
  // v0 sin(beta (X_axis - xref)), other components 0. Returns the count, or -1.
  int setBodyVelocity(int body, const double *v, std::string &err);
  // WP-232i (tests/mpm_bspline_interface.md B11, B12): one flat list, ndm per
  // particle in storage order. The body force (per unit mass, unscaled) may be set
  // between any two `mpm analyze` calls; the velocity only before the first step.
  int setParticleBodyForce(const std::vector<double> &b, std::string &err);
  int setParticleVelocity(const std::vector<double> &v, std::string &err);
  // WP-232i: a B-spline grid face with both `fix` and `wall` is refused (B5); a
  // gridvel box reaching the ghost layer behind a `fix` face too (B10). "" = fine.
  std::string bsplineFaceConflict() const;
  int setBodyVelocitySine(int body, int comp, int axis, double v0, double beta,
                          double xref, std::string &err);
  int setBodyStressGeostatic(int body, double yTop, double rho, double g, double K0,
                             std::string &err);

  // ---- phase c2: the two contracts (LadrunoMPMIO.h) ----------------------------
  // `mpm particles -file <path>`: load a Contract-1 set, all or nothing. The path
  // is opened as given. Returns the count, or -1 (every problem printed).
  int loadParticleSet(const std::string &path);
  // `mpm output ...`: the .ladruno particle-group writer. 0, or -1 with err.
  int setOutput(const ladruno_mpm_io::WriterConfig &cfg, std::string &err);
  void closeOutput(void);
  // `mpm units <L> <F> <T> [-system <id>]`: the units Contract 2 carries when no
  // Contract-1 file declared them.
  int setUnits(const std::string &L, const std::string &F, const std::string &T,
               const std::string &system, std::string &err);
  const std::vector<std::string> &loadedFiles(void) const { return sourcePaths; }
  // Test hooks for the W17 path (`mpm debug injectRefusal | staleRefusal`):
  // a refusal declared at the next particle commit, or a stale count left in
  // the WP-99 counter before it (which must NOT abort).
  void injectCommitRefusal(void);
  void injectStaleRefusal(void);
  bool isAborted(void) const { return aborted; }

  double criticalTimeStep(void);
  // Run nSteps explicit steps of dt. 0 on success; -1 on a refusal or failure
  // (message printed).
  int    analyze(int nSteps, double dt);

  // ---- mode A (WP-232f) -----------------------------------------------------------
  enum Mode { MODE_EXPLICIT = 0, MODE_IMPLICIT = 1, MODE_QUASISTATIC = 2 };
  int    setMode(Mode m, std::string &err);
  Mode   mode(void) const { return theMode; }
  // `mpm analyze` in mode A (interface §3): nSteps of dt through the user's analysis
  int    analyzeImplicit(int nSteps, double dt, const LadrunoMPMAnalysisRefs &A);
  // `mpm gridvel <name> box <lo> <hi> v <value|free ..> [start t0] [track axis]
  //  [translate axis] [profile t0 Tr t1]` (WP-232k D17: `translate` moves both bounds,
  //  `profile` is a trapezoidal velocity factor, present from t = 0)
  int    addGridVel(const std::string &name, const double *lo, const double *hi,
                    const double *v, unsigned freeMask, double start, int track, std::string &err,
                    int translate = -1, const double *profile = 0);
  // `mpm condForce <name>`: ndm doubles at the last accepted step
  int    condForce(const std::string &name, std::vector<double> &out, std::string &err);
  // `mpm pforce <body> box <lo> <hi> <f>`: returns the count of particles
  int    addPointForce(int body, const double *lo, const double *hi, const double *f, std::string &err);
  // WP-232k (S1): `mpm gridvel <name> box <lo> <hi> affine <L: ndm x ndm> [ref <x: ndm>]` -- every
  // node of the set prescribed v = L (x_node - ref), all components; `mpm gridvel <name>
  // setAffine <L>` changes L between `mpm analyze` calls (the harness's per-step dF path)
  int    addGridVelAffine(const std::string &name, const double *lo, const double *hi,
                          const double *L, const double *ref, std::string &err);
  int    setGridVelAffine(const std::string &name, const double *L, std::string &err);
  void   setMassCutoffUser(double c) { eng.par.masscut = c; masscutUser = true; }
  // emptyNodes: -1 exclude (held at zero increment), EMPTY_RENORMALIZE exclude with
  // the kept shape functions renormalised, >= 0 ghost spring factor gamma
  static constexpr double EMPTY_EXCLUDE = -1.0;
  static constexpr double EMPTY_RENORMALIZE = -2.0;
  void   setEmptyNodes(double gamma) { ghostGamma = gamma; }
  void   setCut(double factor, int maxCuts) { cutFactor = factor; maxCutCount = maxCuts; }
  int    setTagBases(int nodeBase, int eleBase, std::string &err);
  void   setRayleighMPM(double am, double bk);
  void   setDebugFail(int step, int attempts, bool afterTrial = false)
    { failStep = step; failAttempts = attempts; failAfterTrial = afterTrial; }
  void   setDebugReseed(int m) { reseedMode = m; }               // 0 default, 1 always, 2 skip
  // alOverride: -1 the deck's setting, 0 anti-locking suppressed for the check (WP-232j L2)
  int    debugTangentCheck(double scale, double eps, std::vector<double> &out, std::string &err,
                           int alOverride = -1);
  // WP-232j (L1): `mpm debug affineCheck <A> <b>`, both modes
  int    debugAffineCheck(const double *A, const double *b, std::vector<double> &out, std::string &err);
  // WP-232j (K3): `mpm antiLocking off|fbar`, fixed for an epoch like the mode
  int    setAntiLocking(int on, std::string &err);
  // WP-232j (J13, L8): `mpm debug packVersion 2|3` -- the restart pack's layout for the next
  // database save (2 = the 232i layout, no anti-locking field: the v2-restores-as-off test)
  int    setPackVersion(int v);
  // WP-232j (L4 r2 diagnosis): `mpm debug failDump <path> [<k>] | off` -- on every mode-A
  // attempt that fails before commit, append the k particles with the largest residual at
  // their support nodes (free DOFs), their trial stress, tangent ratio and apex flag.
  int    setFailDump(const std::string &path, int k);
  // `mpm info` counters (interface §5)
  long   infoCounter(const std::string &what, bool &ok) const;
  const std::vector<double> &dtHistory(void) const { return dtLog; }
  void   gridNodeTags(std::vector<double> &out) const;

  // ---- WP-232k: soil particles (tests/mpm_soil_interface.md) ----------------------
  // D7: the one stage walk, both modes. Every particle material with the tag (and every
  // LadrunoSANISAND particle, whose stage flag is class-wide) registers through its own
  // setParameter(`materialState`, tag) on a host-local Parameter, then update(stage).
  // Returns the number of particle materials that answered, or -1 with err.
  int    stageWalk(int matTag, int stage, std::string &err);
  // D10: "" or why `mpm analyze` is refused (SANISAND's class-wide stage moved by a
  // route that did not reach the particles)
  std::string stageGuard(void) const;
  // D15: "" or the first particle material whose tangent may be unsymmetric
  // (LadrunoSANISAND or LadrunoHS without -implex)
  std::string unsymmetricSoil(void) const;
  // AMEND A5: "" or why a B-spline gridvel box, at time t (it may move: track / translate),
  // reaches the ghost layer behind a `fix` face. staticAllTimes (the up-front check, review
  // item 4): a box that does not move is checked whether or not it is present yet (a future
  // `start`); moving boxes are checked where they are at t.
  std::string gridvelGhostConflict(double t, bool staticAllTimes = false) const;
  // AMEND A4: a SANISAND or HS particle material exists
  bool   anySoil(void) const;
  // D16 queries: sigma_zz per particle (plane strain: the material's; 3D: sig[2]);
  // a material response per particle, flat (sizeOnly: its length per particle);
  // `mpm info materials`: [tag, classTag, scheme, implex, stage, count] per tag
  int    particleSigZZ(std::vector<double> &out, std::string &err);
  int    particleResponse(const std::vector<std::string> &args, bool sizeOnly,
                          std::vector<double> &out, std::string &err);
  void   materialsInfo(std::vector<double> &out);
  // D13: "particle id <id>, nDMaterial <class> tag <t>" of the first refusing particle
  std::string refuserText(int p) const;

  // ---- WP-232l: two-phase particles (tests/mpm_twophase_interface.md) ---------------
  // P1: `mpm twoPhase <body> ...` -- the body's fluid constants; before the first step
  int    setTwoPhaseBody(int body, const ladruno_mpm::TwoPhaseProps &tp, std::string &err);
  // P8 `mpm pfix <face> <p0>`; P10 `mpm plate <face> <F>`; P10b `mpm traction ...`
  int    setPfix(int face, double p0, std::string &err);
  int    setPlate(int face, double F, std::string &err);
  int    addTraction(const std::string &name, int face, const double *lo, const double *hi,
                     const double *t, std::string &err);
  // P6 `mpm kinematics small|hypo`; P5 `mpm pressureBlend <alpha_p>`
  int    setKinematics(int small, std::string &err);
  int    setPressureBlend(double a, std::string &err);
  // §7a Q5: initial particle state (flat, storage order); `mpm porePressure <body> uniform <p>`
  int    setParticlePressure(const std::vector<double> &v, std::string &err);
  int    setParticleAcceleration(const std::vector<double> &a, std::string &err);
  int    setBodyPorePressure(int body, double p, std::string &err);
  // P12: `mpm info twoPhase <body>`
  int    twoPhaseInfo(int body, std::vector<double> &out, std::string &err);
  // §7a: `mpm debug particleBlocks <p> <path> <6 scales>`; `mpm debug assemble <path> <dt>`
  int    debugParticleBlocks(int p, const std::string &path, const double *scales, std::string &err);
  int    debugAssemble(const std::string &path, double dt, const LadrunoMPMAnalysisRefs &A,
                       std::string &err);
  bool   anyTwoPhase(void) const;
  bool   anyPlate(void) const;
  bool   isTwoPhaseParticle(int p) const;

  // ---- LadrunoMPMPointHost --------------------------------------------------------
  bool   mpmSmallKinematics(void) { return eng.par.smallKin != 0; }
  double mpmPressureBlend(void) { return eng.par.pBlend; }
  const double *mpmFluidDrive(int p);
  int    mpmStageParameter(const char **argv, int argc, Parameter &param);   // D9 (c)
  void   mpmTrialRefused(int p);                                             // D14
  ladruno_mpm::Particles &mpmParticles(void) { return eng.P; }
  const double *mpmBodyForce(int p);
  bool   mpmTrialAdmissible(int p, const double *x);
  int    mpmFailThisAttempt(void) { return failNow; }
  bool   mpmInStep(void) { return armedWindow; }
  void   mpmParticleCommitted(int p);
  void   mpmElementGone(int p);
  double mpmGhostSpring(int p);
  // Particles back to the state at the first `mpm analyze` (materials
  // revertToStart). The Domain time is the user's (`reset` / `setTime`).
  int    reset(void);

 private:
  class Mats;                         // the MaterialPoints adapter
  ladruno_mpm::Engine     eng;
  std::vector<NDMaterial *> mats;     // one owned copy per particle
  std::vector<char>       elastic;    // per particle: the initial-stress offset is exact
  Mats                   *adapter;
  TimeSeries             *series;     // owned by LoadPattern; non-owning view
  double                  lambda;
  bool                    inStep;
  bool                    stepped;
  int64_t                 nextLatticeId;
  // the initial state, captured at the first step, restored by reset()
  bool                    haveSnapshot;
  int                     epoch;
  ladruno_mpm::Particles  snapshot;

  // phase c2 state
  ladruno_mpm_io::ParticleWriter *writer;
  double                  thickness;       // 2D; 1 for lattice-seeded particles
  bool                    thicknessFixed;
  std::string             unitL, unitF, unitT, unitSystem;
  bool                    unitsFromFile;
  std::vector<std::string> sourceSetIds, sourcePaths;
  std::map<int, std::string> bodyNames;
  int                     stepsSinceOutput;
  double                  nextOutputTime;
  std::vector<int>        idOrder;         // storage index of each particle, by id

  bool                    aborted;         // a commit refusal aborted the run (W17); `mpm reset` clears it

  // ---- mode A state (WP-232f) -----------------------------------------------------
  struct GridVel {
    std::string name;
    double lo[3], hi[3], v[3];
    unsigned freeMask;          // bit d: component d free
    double start;               // absent before start
    int    track;               // -1, or the axis whose lower bound moves with v
    int    translate = -1;      // WP-232k D17: -1, or the axis along which BOTH bounds move
    bool   profiled = false;    // WP-232k D17: present from t = 0, velocity v * phi(t)
    double pT0 = 0.0, pTr = 0.0, pT1 = 0.0;   // the trapezoid: ramp up from pT0 over pTr,
                                               // ramp down ending at pT1
    bool   affine = false;      // WP-232k (S1): v = L (x - ref) per node, every component
    double L[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, ref[3] = {0, 0, 0};
    // the prescribed component a at a node with coordinates c, at step factor fac
    double valueAt(int a, const double *c, int nd, double fac) const;
    bool   presentAt(double t) const;
    double displacementFactor(double t) const;      // Phi(t) (profile) or max(t - start, 0)
    double stepFactor(double t, double h) const;    // (Phi(t+h) - Phi(t)) / h, or 1
    void   boxAt(double t, double lo_[3], double hi_[3]) const;
    std::vector<int> nodes;     // the node set of the current / last accepted step
    std::vector<double> force;  // cond_force at the last accepted step
  };
  Mode                    theMode;
  bool                    masscutUser;     // `mpm massCutoff` given: wins in both modes
  double                  ghostGamma;      // < 0 exclude (phase 1), else ghost spring factor γ
  double                  cutFactor;
  int                     maxCutCount;
  double                  rayAlphaM, rayBetaK;
  std::vector<GridVel>    gridVels;
  std::vector<double>     pforce;          // ndm per particle, scaled by λ
  std::shared_ptr<LadrunoMPMHostLink> link;
  std::vector<LadrunoMPMPoint *> eles;     // per particle; null before the first mode-A step
  std::vector<char>       matOwnedByEle;   // per particle: the element owns mats[p]
  std::unordered_map<long, int> liveTag;   // engine lattice index -> live grid node tag
  bool                    basesSet;
  int                     nodeBase, nextNodeTag, eleBase;
  std::vector<int>        hostSPTags;      // this step's grid-condition SPs, by tag (owned by the
                                           // pattern; `wipe` may delete them first, so never by pointer)
  struct SpecRow { long ext; int dof; double value; bool constant;
                   bool operator==(const SpecRow &o) const
                   { return ext == o.ext && dof == o.dof && value == o.value && constant == o.constant; } };
  std::vector<SpecRow>    lastSpec;        // the SP set those tags realise, in lattice terms
  int                     hostPayloadDbTag;
  bool                    elementsRestored; // restored with elements in the Domain: relink first
  double                  stepT0;          // Domain time at the start of the attempt
  int                     failNow;         // 0, 1 fail before the trial, 2 after it (review item 13)
  bool                    failAfterTrial;
  bool                    armedWindow;     // inside `mpm analyze` / tangentCheck (review item 1)
  std::unordered_set<int> liveTags;        // the live grid node tags (review item 10)
  int                     failStep, failAttempts, reseedMode;
  int                     armedCount, committedCount;
  bool                    g2pRan, g2pSkipped;
  long                    nG2P, nAccepted, nRetries, nDomainChanges, nConnChanges, epochSteps;
  std::vector<double>     dtLog;
  std::vector<double>     bfBuf;
  std::vector<unsigned char> alPend;       // WP-232j: this attempt's anti-locking path per particle
  int                     packVersionOut = 0;   // WP-232j: 0 = the current layout (one save only)
  std::string             failDumpPath;          // WP-232j: empty = off
  int                     failDumpK = 0;
  void                    writeFailDump(double t, int attempt);
  // mode F's affine check with the particles' materials (defined where Mats is complete)
  void                    engineAffineCheck(const double *A, const double *b, std::vector<double> &o);

  // ---- WP-232l state ------------------------------------------------------------
  std::map<int, ladruno_mpm::TwoPhaseProps> twoPhase;   // body id -> fluid constants
  int                     plateTag[6] = {0, 0, 0, 0, 0, 0};   // the plate unknowns' nodes (0 none)
  struct Traction { std::string name; int face; double lo[3], hi[3], t[3]; };
  std::vector<Traction>   tractions;
  std::vector<std::pair<int, std::vector<double>>> tractionLoads;   // this step: node tag, force
  std::vector<double>     fdBuf;                   // mpmFluidDrive scratch
  int    prepareStepUP(double h, std::string &err);   // a two-phase epoch or a plate (§2, §3)
  int    tagOfExt(long ext) const;                     // a lattice id or plateId(f) -> live tag
  void   retirePlates(void);
  std::string twoPhaseAdmission(const LadrunoMPMAnalysisRefs &A) const;

  double modeAMassCutoff(void) const;
  int    prepareStep(double h, std::string &err);
  int    createElements(std::string &err);
  int    relinkElements(std::string &err);
  void   removeHostSPs(void);
  void   retireNode(long ext);
  void   teardownModeA(void);
  void   updateCondForces(void);
  int    admitAnalysis(const LadrunoMPMAnalysisRefs &A, std::string &err);
  int    writeAfterStep(std::string &err);
  int    packState(Vector &v);
  int    unpackState(const Vector &v);
  struct ParamRef { Domain *domain; int tag; Parameter *param; };
  std::vector<ParamRef>   paramRefs;       // Parameters that hold our materials (F2)

  void detachParameters(void);
  // WP-232k (D9, D10, D13)
  struct StageReg { int matTag; bool intPath; };
  std::vector<StageReg>   stageRegs;       // Parameter ids STAGE_PARAM_BASE + index
  static constexpr int    STAGE_PARAM_BASE = 1000000;
  // the Parameter id of the (matTag, intPath) registration, reusing an existing one (review
  // item 5: the list does not grow with every `updateMaterialStage`)
  int                     stageRegId(int matTag, bool intPath);
  int                     sanisandStage = 0;   // the stage the last walk gave SANISAND particles
  std::map<int, int>      hsWalked;            // AMEND A1: HS tag -> the stage its last walk gave
  int                     restoredPackVersion = 0;   // review item 3: the last unpack's layout
  // WP-232k (review item 7): the soil class tags, in one place
  static bool             isSanisandTag(int classTag);
  static bool             isHSTag(int classTag);
  static bool             isSoilTag(int classTag);
  int                     refusalBase = 0;     // mode A: the WP-99 count at the attempt start
  int                     firstRefuserA = -1;  // mode A: the first particle whose commit refused
  bool soilStageRequest(const char **argv, int argc, int &matTag) const;
  // WP-232k (D18, #1109): mode F's prescribed grid velocities for the step [t, t + h]
  int  buildPrescribed(double t, double h, std::string &err);
  // WP-232k (D14): the first particle whose material refused the trial strain in this attempt
  int  firstTrialRefuserA = -1;
  int  admitMaterial(int matTag, NDMaterial *&copy, bool &isElastic, std::string &err);
  bool domainHasFreeDOF(std::string &which);
  bool isGridNode(int tag) const;
  void nameBody(int body);
  int  openOutputStage(std::string &err);
  int  writeOutputNow(std::string &err);
};

// The host in the Domain, or 0 (v1 allows one per Domain).
LadrunoMPMPattern *ladrunoMPMFindHost(Domain *theDomain);

// `mpm criticalTimeStep | particleData <f> | info time|summary` as doubles
// (OPS_LadrunoMPM.cpp): 1 handled, 0 not such a query, -1 error (printed).
// The classic-Tcl bridge uses it to publish full-precision Tcl doubles.
int ladrunoMPMDoubleQuery(Domain *dom, const std::vector<std::string> &args,
                          std::vector<double> &out);

#endif
