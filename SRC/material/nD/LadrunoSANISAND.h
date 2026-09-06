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
                    int maxSubsteps = 0);   // Ladruno

    // full constructor, classTag defaults to ND_TAG_LadrunoSANISAND.
    // Defaults of the five optional integration args match the base's
    // equivalent (int tag, ...) form: 2, 2, 1, 1e-7, 1e-7.
    LadrunoSANISAND(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
                    double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0,
                    double nd, double z_max, double cz, double mDen,
                    int integrationScheme = 2, int tangentType = 2, int JacoType = 1,
                    double TolF = 1.0e-7, double TolR = 1.0e-7,
                    double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0,
                    int maxSubsteps = 0);   // Ladruno

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

    // Base Vector(97) wire format unchanged; one extra Vector(4) follows it.
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    // --- ADR-86b: the substep cap's two public seams -------------------------

    // Ladruno (ADR-86b): "substeps" / "substepsME" -- the ModifiedEuler substep
    // count of the LAST material update at this integration point, and whether
    // the cap fired. Response IDs are in the 3308x band so they cannot collide
    // with the base's 1..8. Added on the BASE Ladruno class so both wrappers
    // inherit them (neither overrides setResponse).
    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int       getResponse(int responseID, Information &matInformation);

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

    // SHADOW of the non-virtual ManzariDafalias::initialize(). Same signature on
    // purpose -- see the DESIGN NOTE above. DO NOT add `virtual` here or in the
    // base.
    void initialize();

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
