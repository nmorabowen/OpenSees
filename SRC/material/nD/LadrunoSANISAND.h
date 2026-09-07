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

// Ladruno: LadrunoSANISAND -- Manzari-Dafalias (2004) with the two low-stress
// constants made settable, carried on the wire, and echoed at construction.
//
// `ManzariDafalias::initialize()` assigns, unconditionally and with no way to
// override:
//
//     m_Pmin      = 1.0e-4 * m_P_atm;   //  0.0101 kPa at P_atm = 101
//     m_Presidual = 1.0e-2 * m_P_atm;   //  1.01   kPa at P_atm = 101
//
// `m_Presidual` is threaded through ~30 mean-stress evaluations as
// `p = one3*GetTrace(stress) + m_Presidual`, including the one that feeds
// `GetPSI()`, so it displaces the state parameter psi and with it M^b, M^d and
// the dilatancy. It acts as an apparent cohesion c = p_r*tan(phi) ~ 0.95 kPa:
// negligible at 100 kPa, +20 % on the model's own bounding-surface identity at
// 1 kPa. It is also NOT on the vanilla wire (`Vector(97)` carries `m_Pmin` at
// data(96) and nothing for `m_Presidual`), so an OpenSeesMP worker or a
// database-restored material runs at p_r = 0 while the serial process beside it
// runs at 1.01 kPa.
//
// This class changes exactly two constants and nothing else. Defaults follow
// NTUASand02 (Gorini): p_residual = 0 (a cohesionless sand has no cohesion) and
// p_min = 1.0e-3 * P_atm.
//
// DESIGN NOTE (do not "clean up"):
//   `m_Presidual` and `m_Pmin` are PROTECTED DATA read at run time by every base
//   integrator, so this subclass never has to override an integrator -- it only
//   has to win the LAST WRITE. Every method it must control for that is already
//   virtual through NDMaterial/MovableObject: revertToStart, sendSelf, recvSelf,
//   both getCopy forms, Print.
//
//   `initialize()` is NOT virtual and MUST NOT BECOME VIRTUAL. The declaration
//   below is deliberately a *shadow*: the base constructor's `initialize()` call
//   still binds to the base version (correct -- the base must run its own init),
//   and the derived constructor re-applies the two constants afterwards.
//   `ManzariDafaliasRO` (UWmaterials/ManzariDafaliasRO.h:93-97) shadows
//   `initialize()` and two `GetElasticModuli` overloads with identical
//   signatures; adding `virtual` to any of those in the base would silently turn
//   those shadows into overrides and make every existing RO deck run
//   Ramberg-Osgood elasticity from inside the base integrators.
//
// Vanilla `ManzariDafalias` is not edited.
// See Ladruno_implementation/86_ladruno_sanisand_adr.md.
// classTags 33019 (base) / 33020 (3D) / 33021 (PlaneStrain).
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoSANISAND_h
#define LadrunoSANISAND_h

// Quote-include so it resolves against THIS file's own directory first --
// SRC/material/nD/, where fork-authored ND materials live flat (the LadrunoJ2
// convention, ADR 86 sec.4.2) -- which is why the path carries the
// `UWmaterials/` segment: the base lives one level down, in
// SRC/material/nD/UWmaterials/. Same spelling FEM_ObjectBrokerAllClasses.cpp uses.
// (An earlier version of this comment said THIS file lives in UWmaterials/. It
// does not; the include resolved correctly anyway, so nothing depended on the
// error, but on this class the comments have been wrong more often than the code.)
#include "UWmaterials/ManzariDafalias.h"

// ===========================================================================
//  Ladruno (ADR-92 P1): the IMPL-EX option block.
//
//  Carried as ONE aggregate rather than seven more trailing constructor
//  arguments. The four constructors already take up to 28; adding seven more
//  would have to be repeated in both full constructors, both wrapper
//  constructors and getCopy(const char*) -- exactly the duplication
//  sanitiseLadrunoInputs() exists to undo. The parser (and getCopy(const char*))
//  construct first and then call setLadrunoImplexOptions(); the wrappers'
//  getCopy(void) is a memberwise `*clone = *this` and carries this struct for
//  free.
//
//  DEFAULTS ARE "OFF" AND THAT IS LOAD-BEARING: with `enabled == false` every
//  code path added by ADR 92 is skipped by a single branch in
//  ladrunoTrialUpdate() / commitState() / revertToLastCommit(), so a deck that
//  does not say -implex is byte-identical to the pre-ADR-92 build.
// ===========================================================================
struct LadrunoImplexOptions            // Ladruno (ADR-92 P1)
{
    // dt source (ADR-92 D2). `pseudo` = ops_Dt, the ASDConcrete3D behaviour and
    // the default; `strain` = the norm of the strain increment; `user` = a value
    // the deck supplies with `-implexDt user $dt` (or at run time through
    // setParameter "implexDt").
    enum DtSource { DT_PSEUDO = 0, DT_STRAIN = 1, DT_USER = 2 };

    // Ladruno ADR-92 P2-1: what -implexControl does at the reduction floor, where
    // there is nothing left to cut. `implicit` (the DEFAULT) delivers the
    // COMPANION's stress for that step at that Gauss point, under the frozen Ce
    // -- the operator stays SPD, the committed state cannot carry an O(1) error,
    // and the self-sustaining gap-closing loop ADR 93's 2026-09-07 census
    // measured has nothing left to feed on. `accept` is the pre-P2 behaviour
    // (commit the extrapolation whatever its error). `refuse` returns the
    // sentinel at the floor too -- an honest wall instead of a creeping curve.
    enum FloorMode { FLOOR_IMPLICIT = 0, FLOOR_ACCEPT = 1, FLOOR_REFUSE = 2 };

    bool   enabled;          // -implex
    bool   control;          // -implexControl
    double errorTol;         // -implexControl $tol
    double reductionLimit;   // -implexControl ... $reductionLimit
    double alpha;            // -implexAlpha
    int    dtSource;         // -implexDt {pseudo|strain|user}
    double dtUser;           // -implexDt user $dt
    int    floorMode;        // -implexFloor {implicit|accept|refuse}  Ladruno ADR-92 P2-1
    bool   guard;            // -implexGuard {on|off}                  Ladruno ADR-92 P2-2

    // Ladruno ADR-92 P2-6: the TRIAL-time graded fallback. Esmeralda 146569
    // (dense q10) measured the leg crawling to its subdivision budget because
    // W7 refuses a trial extrapolated with the full f the instant its error
    // crosses tol, even when the committed predecessor gave no advance warning
    // (the P2-2 guard only catches a predecessor that ALREADY showed Kp <= 0 or
    // a reversal). ON (the DEFAULT) tries f = 0 for THIS trial, before
    // refusing, whenever the primary extrapolation's error is past tol and the
    // reduction floor has not been reached; OFF reproduces the pre-P2-6
    // behaviour (refuse immediately, exactly as P2-1..P2-5 shipped it).
    bool   trialGuard;       // -implexTrialGuard {on|off}             Ladruno ADR-92 P2-6

    // errorTol default: measured 2026-09-06 (_adr92_p1_bvp_gate_rerun.md sweep)
    // -- 0.05 fails on reach, 0.1 is the tightest tolerance that beats the
    // implicit control's depth under 5% mean deviation. (WP-92d)
    LadrunoImplexOptions()
      : enabled(false), control(false), errorTol(0.1), reductionLimit(0.01),
        alpha(1.0), dtSource(DT_PSEUDO), dtUser(0.0),
        floorMode(FLOOR_IMPLICIT), guard(true),
        trialGuard(true) {}   // Ladruno ADR-92 P2 / P2-6
};

class LadrunoSANISAND : public ManzariDafalias
{
  public:

    // full constructor, explicit classTag -- used by the parser and by the
    // LadrunoSANISAND3D / LadrunoSANISANDPlaneStrain wrappers.
    // Defaults of the five optional integration args match the base's
    // equivalent (int tag, int classTag, ...) form: 1, 0, 1, 1e-7, 1e-7.
    LadrunoSANISAND(int tag, int classTag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
                    double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0,
                    double nd, double z_max, double cz, double mDen,
                    int integrationScheme = 1, int tangentType = 0, int JacoType = 1,
                    double TolF = 1.0e-7, double TolR = 1.0e-7,
                    double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0,
                    int maxSubsteps = 0,
                    double reversalTol = 1.0e-10,    // Ladruno ADR-92 P2-5
                    double reversalRel = 0.05);      // Ladruno ADR-92 P2-5b

    // full constructor, classTag defaults to ND_TAG_LadrunoSANISAND.
    // Defaults of the five optional integration args match the base's
    // equivalent (int tag, ...) form: 2, 2, 1, 1e-7, 1e-7.
    LadrunoSANISAND(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
                    double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0,
                    double nd, double z_max, double cz, double mDen,
                    int integrationScheme = 2, int tangentType = 2, int JacoType = 1,
                    double TolF = 1.0e-7, double TolR = 1.0e-7,
                    double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0,
                    int maxSubsteps = 0,
                    double reversalTol = 1.0e-10,    // Ladruno ADR-92 P2-5
                    double reversalRel = 0.05);      // Ladruno ADR-92 P2-5b

    // specific-type null constructor (used by the wrappers' null constructors)
    LadrunoSANISAND(int classTag);

    // null constructor (FEM_ObjectBroker)
    LadrunoSANISAND();

    // destructor
    ~LadrunoSANISAND();

    // --- the six overrides. See ADR 86 section 4.4; none is optional. -------

    // Replicates the base ops_InitialStateAnalysis guard but re-initialises
    // through THIS class's initialize(), so p_residual is not restored to
    // 1.0e-2*P_atm mid-analysis.
    int revertToStart(void);

    // The base getCopy(const char*) hardcodes `new ManzariDafalias3D(...)` /
    // `new ManzariDafaliasPlaneStrain(...)`, which would strip these settings at
    // every Gauss point.
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);

    // Base Vector(97) wire format unchanged; one extra Vector(26) follows it
    // (Vector(4) at ADR-86, 5 at ADR-86b, 22 at ADR-92 P1, 25 at ADR-92 P2,
    // 26 at ADR-92 P2-5).
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    // --- ADR-92 P1: IMPL-EX --------------------------------------------------

    // The companion return and the extrapolation history live here, so the
    // commit hook is a seventh override. `ManzariDafalias::commitState()` is
    // virtual through NDMaterial; when -implex is off this delegates to it
    // unconditionally and the deck is byte-identical.
    int commitState(void);              // Ladruno (ADR-92 P1)

    // The base's is an empty `return 0`. Under -implex the trial members hold an
    // EXTRAPOLATED state that has no companion behind it, so a rejected step
    // must put them back on the committed ones and re-arm the step. Guarded on
    // the flag -- with -implex off this is the base's empty body verbatim.
    int revertToLastCommit(void);       // Ladruno (ADR-92 P1)

    // Set after construction (see the struct's note). Validates the option set
    // against this deck's IntScheme and -maxSubsteps and refuses what ADR 92 D3
    // refuses; returns 0 on accept, -1 if the options were rejected (in which
    // case IMPL-EX is left OFF and the material runs implicit).
    // `verbose` is false on the getCopy / recvSelf paths: those run once per
    // Gauss point, and the echo rule for this class (ADR 86 sec.4.4) is
    // once per DECK command, not once per integration point.
    int  setLadrunoImplexOptions(const LadrunoImplexOptions &opt,
                                 bool verbose = true);               // Ladruno (ADR-92 P1)
    const LadrunoImplexOptions &getLadrunoImplexOptions(void) const { return mImplexOpt; }

    // `implexError` / `avgImplexError`, on the ASDConcrete3DMaterial.cpp
    // :2073-2077 template, plus this material's own per-point detail response.
    int setParameter(const char **argv, int argc, Parameter &param);  // Ladruno (ADR-92 P1)
    int updateParameter(int parameterID, Information &info);          // Ladruno (ADR-92 P1)

    // --- ADR-86b: the substep cap's two public seams -------------------------

    // Ladruno (ADR-86b): "substeps" / "substepsME" -- the ModifiedEuler substep
    // count of the LAST material update at this integration point, and whether
    // the cap fired. Response IDs are in the 3308x band so they cannot collide
    // with the base's 1..8. Added on the BASE Ladruno class so both wrappers
    // inherit them (neither overrides setResponse).
    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int       getResponse(int responseID, Information &matInformation);

    // --- ADR-92 P1: the one entry point both wrappers' setTrialStrain use -----
    //
    // Replaces the bare `this->integrate(); return this->ladrunoUpdateStatus();`
    // pair. With -implex OFF it IS that pair, in that order, so both wrappers
    // stay byte-identical; with -implex ON it runs the extrapolated update
    // instead. Public because the wrappers call it from their own
    // setTrialStrain, which is public, and because a test may want to drive it.
    int ladrunoTrialUpdate(void);   // Ladruno (ADR-92 P1)

  protected:

    // Ladruno (ADR-86b): the status the wrappers' setTrialStrain returns.
    // 0 = the update integrated; -1 = the substep cap fired, the strain
    // increment was NOT integrated, and the committed state is untouched.
    // Both LadrunoSANISAND3D and LadrunoSANISANDPlaneStrain call this rather
    // than reading mSubstepCapHitInME directly, so the rule lives in one place.
    int ladrunoUpdateStatus(void) const;   // Ladruno

    double mPresidualInput;   // Ladruno: p_residual as given by the user (>= 0)
    double mPminInput;        // Ladruno: p_min as given; < 0 is the SENTINEL for
                              //          "resolve to 1.0e-3 * P_atm" (P_atm is not
                              //          known until the base ctor has run)
    int    mHonorTolR;        // Ladruno: the DECK-LEVEL request, 0 or 1. Wired in PR-3
                              //          to the base's `bool mHonorTolRInME` (opened by
                              //          PR-2), which ModifiedEuler() reads as
                              //          `TolE = mHonorTolRInME ? mTolR : 1e-4`.
                              //          0 = vanilla's hardcoded 1e-4 substep error
                              //          tolerance; 1 = honour the deck's TolR.
                              //          Two names on purpose: this is the request,
                              //          mHonorTolRInME is the base-side seam it acts
                              //          on, and keeping them distinct keeps the two
                              //          greppable apart.
    int    mMaxSubsteps;      // Ladruno (ADR-86b): the DECK-LEVEL request, >= 0.
                              //          0 = UNCAPPED = vanilla, the default, so every
                              //          existing deck stays byte-identical. A positive
                              //          value is wired to the base's `int
                              //          mMaxSubstepsInME` seam, which
                              //          ManzariDafalias::ModifiedEuler() reads as a
                              //          substep-COUNT cap and on which it FAILS the
                              //          update rather than force-accepting. Same
                              //          two-name convention as mHonorTolR above: this
                              //          is the request, mMaxSubstepsInME is the
                              //          base-side seam it acts on.

    // Ladruno ADR-92 P2-5: absolute strain-increment threshold below which
    // ManzariDafalias::integrate()'s unconditional loading-reversal reset
    // (mAlpha_in := mAlpha_n on a sign test with no magnitude guard,
    // ManzariDafalias.cpp:1002-1013) is undone by
    // ladrunoGuardReversalNoise(). Measured on Esmeralda 146458: on a
    // zero-increment step (a `LoadControl 0.0` hold) the Gauss-point strain
    // increment is round-off (1e-12..1e-16), its sign is noise, and the
    // "reversal" it triggers commits at 28-54% of 34,560 points, sending the
    // plastic modulus h ~ 1/|(alpha - alpha_in):n| toward infinity.
    // Default 1.0e-10 (absolute strain); 0 disables the guard. Applies
    // unconditionally -- with -implex OFF exactly as with it ON -- because
    // the defect lives in the base's integrate(), not in ADR-92's own code.
    double mReversalTol;      // Ladruno ADR-92 P2-5

    // Ladruno ADR-92 P2-5b: P2-5's ABSOLUTE floor cannot work by itself -- a
    // LoadControl(0.0) hold's per-Gauss-point strain increment tracks the
    // SOLVER's tolerance, not round-off. Measured on the R3 footing (1600
    // GPs, LoadControl 0.0): median ~4e-9, max 6.4e-8 (IMPL-EX) / 1.4e-6
    // (implicit), so 99-100% of points exceed 1e-10 and mAlpha_in still
    // resets at 42% (IMPL-EX) / 9.5% (implicit) of them; at 1e-7 the
    // implicit arm still resets 2.2%. -reversalRel scales a SECOND,
    // RELATIVE threshold off mDEpsNormCommit (below): the guard predicate in
    // ladrunoGuardReversalNoise() becomes
    //
    //     ||d_eps|| < max(mReversalTol, mReversalRel * mDEpsNormCommit)
    //
    // Default 0.05: a hold's increment is <= 1e-2 of the previous real
    // step's, a genuine reversal step is ~1x it, and a halved retry is
    // 0.5x -- 0.05 separates a hold from either with margin. 0 disables the
    // relative part and recovers P2-5's pure absolute test exactly.
    double mReversalRel;      // Ladruno ADR-92 P2-5b

    // Ladruno ADR-92 P2-5b: the COMMITTED reference the relative threshold
    // above reads -- ||eps_n - eps_{n-1}||, GetNorm_Cov (the same convention
    // ladrunoGuardReversalNoise() uses for d_eps: mEpsilon is strain, i.e.
    // covariant/engineering-shear-component, not stress; see the base's own
    // comments at ManzariDafalias.cpp:5279-5299). Updated at EVERY commit
    // EXCEPT a zero-increment one (composes with P2-3: on the implicit path
    // a hold is `ops_Dt == 0.0`; on the IMPL-EX path it reuses P2-3's own
    // `mImplexDt == 0.0` predicate) -- a hold must not overwrite the
    // reference the NEXT step's threshold reads, on the exact same
    // reasoning P2-3 already applies to mImplexDtCommit / mImplexDEpsP.
    // Zeroed on revertToStart (ladrunoImplexInitState); carried through
    // getCopy(const char*) and the wire (sendSelf/recvSelf) like any other
    // per-Gauss-point committed history variable.
    double mDEpsNormCommit;   // Ladruno ADR-92 P2-5b

    // =======================================================================
    //  Ladruno (ADR-92 P1): IMPL-EX state.
    //
    //  D1 (measured at P0, 18-22x over the dGamma form on the deck-default
    //  integrator) extrapolates the PLASTIC STRAIN, and that is why this block
    //  needs NO hook into vanilla: eps_p(n) = mEpsilon_n - mEpsilonE_n is
    //  already committed and exact on the base, so the only new committed
    //  quantity is its increment. ADR 92 D5: vanilla footprint ZERO.
    // =======================================================================

    LadrunoImplexOptions mImplexOpt;   // Ladruno (ADR-92 P1): the deck's request

    // d_eps_p(n) = eps_p(n) - eps_p(n-1), COMMITTED. The one new history
    // variable. Zero until the first commit after the elastic-plastic stage
    // flip, which is what makes the first plastic step a pure elastic
    // prediction (and what the P0 oracle's `Implex.__init__` does).
    Vector mImplexDEpsP;               // Ladruno (ADR-92 P1)

    // Per-step bookkeeping.
    double mImplexDt;         // dt_{n+1}, frozen at the FIRST trial call of the step
    double mImplexDtCommit;   // dt_n, the frozen value of the last committed step
    double mImplexDt0;        // Ladruno ADR-92 fix (red/blue B2): the MAGNITUDE |dt|
                              //          of the first non-zero dt seen -- the
                              //          -implexControl reduction floor, which the
                              //          parser documents as a fraction of it. It armed
                              //          only from a POSITIVE dt, so on a
                              //          settlement-driven deck it stayed 0.0 and the
                              //          floor was dead.
    double mImplexFactor;     // f = (dt_{n+1}/dt_n)*alpha, frozen with mImplexDt.
                              //          Ladruno ADR-92 fix (red/blue B1): the ratio is
                              //          formed whenever dt_n != 0, sign-consistently,
                              //          so a monotonically negative clock gets the
                              //          positive ratio it is entitled to instead of
                              //          collapsing to alpha.
    bool   mImplexStepArmed;  // true until the first trial call after a commit/revert
    bool   mImplexTrialDone;  // the last trial pass was an EXTRAPOLATED one, so
                              // commitState() owes a companion return. False on
                              // stage 0 and with -implex off, which is what keeps
                              // gravity and the LoadControl 0.0 hold bit-identical.

    // Error accounting, all measured at commit against the companion return.
    // The deviatoric and volumetric splits exist because ADR-92 P1 section 4
    // makes the clamp-vs-bound decision on them: the p_min clamp repairs
    // tr(sigma~) only, so if the deviatoric error is the same order the
    // pressure error was, the clamp is the wrong fix and a bound on f replaces
    // it. Measured, reported, not asserted.
    // Ladruno ADR-92 P2-2: the softening / reversal guard, ARMED AT COMMIT and
    // consumed by the NEXT step's arm. `reversal` = the committed mAlpha_in_n
    // moved in this commit (the base's loading-reversal detector fired, so the
    // previous plastic increment belongs to a branch the material has just left);
    // `softening` = the plastic modulus Kp at the NEW committed state is <= 0.
    // Either makes d_eps_p(n) the wrong thing to extrapolate -- the seat replay
    // (_adr93_seat_replay.md) measured err 0.4625 there, and f = 0 on that step
    // brings it to 0.029, under tol, with no refusal. The two components are kept
    // separately only so Print() can say WHICH fired.
    bool   mImplexGuardArmed;      // Ladruno ADR-92 P2-2: reversal || softening
    bool   mImplexGuardReversal;   // Ladruno ADR-92 P2-2
    bool   mImplexGuardSoftening;  // Ladruno ADR-92 P2-2

    double mImplexError;      // ||sigma~ - sigma_impl|| / (||sigma_impl|| + P_atm*||eps_n||)
    double mImplexErrorDev;   // the deviatoric part of the same quotient
    double mImplexErrorVol;   // the volumetric part (sqrt(3)*|dp|) of the same quotient
    bool   mImplexClampFired; // the p_min clamp acted on the LAST extrapolation
    long   mImplexClampCount; // how often it has acted at this integration point

    // Ladruno ADR-92 P2-7: has the elastic (0) -> plastic (nonzero) stage flip
    // ALREADY run its zero-increment drift-absorbing companion return? Guards
    // updateParameter()'s W5 hook against re-absorbing on a redundant repeat
    // call (`updateMaterialStage 1` issued again while already at stage 1),
    // which would otherwise re-run integrate()+commitState() on an unchanged
    // committed state every time. NOT touched by ladrunoImplexInitState() --
    // its lifetime is the plastic stage itself, not a single step -- so it is
    // set false in the constructors/initialize() below and reset false only
    // when the stage goes back to 0 (see updateParameter()).
    bool   mImplexStageFlipAbsorbed;   // Ladruno ADR-92 P2-7

    // SHADOW of the non-virtual ManzariDafalias::initialize(). Same signature on
    // purpose -- see the DESIGN NOTE above. DO NOT add `virtual` here or in the
    // base.
    void initialize();

    // --- ADR-92 P1 internals -------------------------------------------------

    // Sizes mImplexDEpsP and zeroes every IMPL-EX scalar. Called from all four
    // constructors and from revertToStart via initialize(); it does NOT touch
    // mImplexOpt, because the option set is the deck's request and survives a
    // revertToStart exactly as p_residual does.
    void ladrunoImplexInitState(void);         // Ladruno (ADR-92 P1)

    // true only when an extrapolated update is actually owed: the flag is on AND
    // the material is past the elastic stage. mElastFlag is a STATIC shared by
    // every instance and is flipped by `updateMaterialStage`, so this is the
    // whole of W5's stage-0 inertness.
    bool ladrunoImplexActive(void) const;      // Ladruno (ADR-92 P1)

    // The extrapolated update: freezes Ce at the committed mean stress, advances
    // sigma_n incrementally, clamps at p_min, and (with -implexControl) runs the
    // in-step implicit and refuses past tolerance.
    int  ladrunoImplexTrial(void);             // Ladruno (ADR-92 P1)

    // The companion return, at commitState only.
    int  ladrunoImplexCommit(void);            // Ladruno (ADR-92 P1)

    // dt_{n+1} for this step from the configured source, and f from it. Called
    // ONCE per step (mImplexStepArmed), so f is a constant within the step and
    // d(sigma~)/d(eps) is exactly Ce -- acceptance 1. Under `-implexDt strain`
    // that freeze is not a convenience, it is the reason the tangent identity
    // survives at all.
    void ladrunoImplexArmStep(void);           // Ladruno (ADR-92 P1)

    // Put every trial member back on its committed twin. Used by
    // revertToLastCommit and by the -implexControl in-step probe, which has to
    // undo the implicit return it just ran before handing the element the
    // extrapolated state.
    void ladrunoRestoreTrialFromCommitted(void);   // Ladruno (ADR-92 P1)

    // Ladruno ADR-92 P2-5 / P2-5b: undoes ManzariDafalias::integrate()'s
    // loading-reversal reset when the strain increment that triggered it was
    // round-off/solver noise, not a real reversal. Call AFTER integrate() (or
    // the -implexControl probe's own call to it) returns, while mEpsilon /
    // mEpsilon_n still hold the increment integrate() just read. See the
    // mReversalTol / mReversalRel member comments for the defect this repairs
    // and the relative threshold P2-5b added.
    //
    // Returns true when this call classified mEpsilon - mEpsilon_n as noise
    // (below the combined threshold) and therefore undid the reset -- read by
    // ladrunoImplexCommit() to gate the P2-2 reversal/softening detection at
    // commit with the SAME criterion; the trial-time call sites discard it.
    //
    // Ladruno ADR-92 P2-5c: `ops_Dt == 0.0` is checked FIRST and is
    // sufficient on its own, ahead of (and regardless of) mReversalTol /
    // mReversalRel -- a hold is a fact the DOMAIN reports, not something the
    // material has to infer from a strain norm, so it is not vulnerable to
    // the same undershoot P2-5b fixed for far-field / early-push points
    // whose own last increment was itself tiny (R3 footing: 136/1600
    // IMPL-EX, 88/1600 implicit still reset at the relative floor). Under
    // `-implexDt user` a user-supplied dt of exactly 0 is ALSO a hold by
    // this rule (mImplexDt tracks ops_Dt only for the DT_PSEUDO default; a
    // user or strain dt source reads a different `ops_Dt == 0.0` fact --
    // still true when the domain itself is not advancing time). A
    // genuinely load-controlled hold under DisplacementControl has d(lambda)
    // = 0 rather than ops_Dt = 0 and is NOT caught by this rule -- it is a
    // different kind of hold, unaddressed here.
    bool ladrunoGuardReversalNoise(void);          // Ladruno ADR-92 P2-5 / P2-5b / P2-5c

    // implexError and its deviatoric / volumetric split, on ADR 92 section 2's
    // definition. `epsRef` is the strain the denominator is scaled by: the P0
    // oracle uses the NEW committed strain, so the commit path passes
    // mEpsilon_n and the in-step -implexControl probe passes mEpsilon (the same
    // vector, one commit earlier).
    void ladrunoImplexMeasureError(const Vector &sigTilde,
                                   const Vector &sigImplicit,
                                   const Vector &epsRef);      // Ladruno (ADR-92 P1)

    // Write Ce(p_n) into all three tangent slots, so getTangent() returns the
    // operator the extrapolated stress was actually built with whatever TanType
    // the deck asked for. Returns K and G through the reference arguments.
    void ladrunoImplexFreezeTangent(double &K, double &G);   // Ladruno (ADR-92 P1)

    // Ladruno ADR-92 P2-2: the plastic modulus at the CURRENT COMMITTED state,
    // Kp = (2/3) p h (b:n), computed from ONE ManzariDafalias::GetStateDependent()
    // call on (mSigma_n, mAlpha_n, mFabric_n, mVoidRatio, mAlpha_in_n). This is
    // the base's own expression, verbatim from ManzariDafalias.cpp:1373 /
    // GetElastoPlasticTangent(:4954) -- the source-true seam, not a proxy. Called
    // ONCE per commit, only with -implex on and -implexGuard on. Not const:
    // GetStateDependent() is not.
    double ladrunoImplexCommittedKp(void);   // Ladruno ADR-92 P2-2

    // The single "win the last write" helper. Called from every constructor
    // (after the base ctor has returned), from revertToStart via initialize(),
    // and from recvSelf. Sets m_Presidual, m_Pmin AND (PR-3) the base's
    // mHonorTolRInME seam -- all three are base-side state that something else
    // may have written first, and all three are re-asserted in exactly one place.
    void applyLadrunoConstants(void);   // Ladruno

    // Per-construction echo (NOT latched -- see ADR 86 section 4.4).
    void echoLadrunoConstants(void);    // Ladruno

    // Ladruno (ADR-86 PR-3): the three input checks, in one place rather than
    // copied into both full constructors. Called before applyLadrunoConstants().
    void sanitiseLadrunoInputs(int tag);   // Ladruno

    // Ladruno (ADR-86 PR-3): true when this deck's IntScheme actually routes to
    // ManzariDafalias::ModifiedEuler(), the ONE function that reads the
    // mHonorTolRInME seam. On every other scheme `-honorTolR 1` is stored, echoed
    // and wired -- and inert, which is the defect class this ADR exists to fix, so
    // the echo says so. See the dispatch table at the definition.
    bool schemeReachesModifiedEuler(void) const;   // Ladruno
};

#endif
