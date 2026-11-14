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
  // `mpm gridvel <name> box <lo> <hi> v <value|free ..> [start t0] [track axis]`
  int    addGridVel(const std::string &name, const double *lo, const double *hi,
                    const double *v, unsigned freeMask, double start, int track, std::string &err);
  // `mpm condForce <name>`: ndm doubles at the last accepted step
  int    condForce(const std::string &name, std::vector<double> &out, std::string &err);
  // `mpm pforce <body> box <lo> <hi> <f>`: returns the count of particles
  int    addPointForce(int body, const double *lo, const double *hi, const double *f, std::string &err);
  void   setMassCutoffUser(double c) { eng.par.masscut = c; masscutUser = true; }
  void   setEmptyNodes(double gamma) { ghostGamma = gamma; }      // < 0: exclude
  void   setCut(double factor, int maxCuts) { cutFactor = factor; maxCutCount = maxCuts; }
  int    setTagBases(int nodeBase, int eleBase, std::string &err);
  void   setRayleighMPM(double am, double bk) { rayAlphaM = am; rayBetaK = bk; }
  void   setDebugFail(int step, int attempts) { failStep = step; failAttempts = attempts; }
  void   setDebugReseed(int m) { reseedMode = m; }               // 0 default, 1 always, 2 skip
  // `mpm info` counters (interface §5)
  long   infoCounter(const std::string &what, bool &ok) const;
  const std::vector<double> &dtHistory(void) const { return dtLog; }
  void   gridNodeTags(std::vector<double> &out) const;

  // ---- LadrunoMPMPointHost --------------------------------------------------------
  ladruno_mpm::Particles &mpmParticles(void) { return eng.P; }
  const double *mpmBodyForce(int p);
  bool   mpmTrialAdmissible(int p, const double *x);
  bool   mpmFailThisAttempt(void) { return failNow; }
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
  bool                    failNow;
  int                     failStep, failAttempts, reseedMode;
  int                     armedCount, committedCount;
  bool                    g2pRan, g2pSkipped;
  long                    nG2P, nAccepted, nRetries, nDomainChanges, nConnChanges, epochSteps;
  std::vector<double>     dtLog;
  std::vector<double>     bfBuf;

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
