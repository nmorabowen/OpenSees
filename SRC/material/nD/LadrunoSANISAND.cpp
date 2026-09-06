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

// Ladruno: LadrunoSANISAND -- Manzari-Dafalias with settable, wired and echoed
// low-stress constants (p_residual, p_min). See LadrunoSANISAND.h for the design
// note and Ladruno_implementation/86_ladruno_sanisand_adr.md for the audit that
// motivated it.
//
// The underlying constitutive model is that of Alborz Ghofrani and Pedro Arduino
// (U. Washington), after Dafalias & Manzari (2004). Nothing in that model is
// changed here.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "LadrunoSANISAND.h"
#include "LadrunoSANISAND3D.h"            // Ladruno: written by stage S2
#include "LadrunoSANISANDPlaneStrain.h"   // Ladruno: written by stage S2

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>                  // ops_InitialStateAnalysis (declared at OPS_Globals.h:75)
#include <classTags.h>                    // ND_TAG_LadrunoSANISAND{,3D,PlaneStrain} -- added by stage S3
#include <elementAPI.h>
#include <MaterialResponse.h>             // Ladruno (ADR-86b): the -maxSubsteps diagnostic response
#include <Information.h>
#include <LadrunoMaterialStatus.h>        // Ladruno (ADR-86b): LADRUNO_MATERIAL_REFUSED

#include <string.h>
#include <math.h>

// ===========================================================================
//  OPS parser
//
//    nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi \
//        $P_atm $m $h0 $ch $nb $A0 $nd $z_max $cz $Rho                         \
//        <$IntScheme $TanType $JacoType $TolF $TolR>                           \
//        <-Presidual $pr> <-Pmin $pmin> <-honorTolR 0|1> <-maxSubsteps $n>
//
//  The first 18 positional doubles and the 5 positional optionals occupy the
//  SAME SLOTS as `nDMaterial ManzariDafalias`, so a deck migrates by renaming
//  the command. ONE default differs, deliberately, as of ADR-86b:
//
//      TanType  vanilla parser 0 (ELASTIC tangent) -> here 2 (CONSISTENT)
//
//  Rationale (ADR-90 GATE U, LEDGER_quirks): `ManzariDafalias3D::getTangent()`
//  returns `mCe` for TanType 0, so a deck that emits only the 18 positional
//  parameters hands `algorithm Newton` an elastic tangent -- modified Newton
//  with linear convergence, dressed as full Newton. It is invisible on a
//  zero-free-DOF material-point deck (there are no equations to iterate) and
//  cost ~7x of wall time on the first boundary-value problem this fork put the
//  material into. Vanilla `OPS_ManzariDafaliasMaterial` KEEPS its own default
//  of 0 (`ManzariDafalias.cpp:93`) -- it is not changed, because every existing
//  vanilla deck and golden file depends on it. So the two parsers now disagree
//  on this one slot, ON PURPOSE, and both say so.
//  `mCep_Consistent` is UNSYMMETRIC under non-associated flow: pair TanType 2
//  with an unsymmetric solver (`system Pardiso -matrixType 0`, UmfPack, ...).
//
//  ARGUMENT ARITHMETIC -- this parser deliberately does NOT compute
//  `numArgs - k` to decide how many optionals to read. `OPS_SAniSandMSMaterial`
//  (UANDESmaterials/SAniSandMS.cpp ~:134, ~:143) does, and its arithmetic is off
//  by one in two places, which SILENTLY DROPS TolF and TolR (ADR 86 section 7.1).
//  Instead every remaining token is consumed one at a time and classified, so the
//  count never has to be right. Hand-checked for 19..25 total args:
//    19 -> tag + 18 doubles,      0 positional optionals, all 5 defaulted
//    20 -> + IntScheme
//    21 -> + TanType
//    22 -> + JacoType
//    23 -> + TolF      (the leg SAniSandMS drops)
//    24 -> + TolR      (the leg SAniSandMS drops)
//    25 -> hard error: more than 5 positional optionals
//  and with flags mixed in, e.g. 22 = tag + 18 + IntScheme + "-Pmin" + $pmin.
// ===========================================================================

static int numLadrunoSANISANDMaterials = 0;

// Ladruno (ADR-86 PR-3): IntScheme numbers, re-declared here because
// ManzariDafalias's own INT_* macros are defined in ManzariDafalias.cpp (:38-49)
// and are therefore not visible outside that TU. Prefixed so they can never
// collide with the base's if that ever changes. Only the three this file needs.
// KEEP IN SYNC with ManzariDafalias.cpp:38-49 -- there is no compile-time link
// between the two, which is why schemeReachesModifiedEuler() spells out the
// dispatch it is asserting rather than trusting the names.
#define INT_LSANISAND_MAXENE_MFE     0    // ManzariDafalias INT_MAXENE_MFE
#define INT_LSANISAND_ModifiedEuler  1    // ManzariDafalias INT_ModifiedEuler
#define INT_LSANISAND_RungeKutta45  45    // ManzariDafalias INT_RungeKutta45

void *
OPS_LadrunoSANISAND(void)
{
    // One-time authorship credit (this is NOT the parameter echo -- that one is
    // per construction and lives in the constructors; see ADR 86 section 4.4).
    if (numLadrunoSANISANDMaterials == 0)
        opserr << "LadrunoSANISAND nDMaterial - Ladruno subclass of ManzariDafalias "
               << "(model: A.Ghofrani, P.Arduino, U.Washington); settable p_residual/p_min. ADR 86.\n";
    numLadrunoSANISANDMaterials++;

    // Counted BEFORE the tag is consumed, and used ONLY for the minimum check.
    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 19) {
        opserr << "Want: nDMaterial LadrunoSANISAND tag? G0? nu? e_init? Mc? c? lambda_c? e0? ksi?"
               << " P_atm? m? h0? Ch? nb? A0? nd? z_max? cz? Rho?"
               << " <IntScheme? TanType? JacoType? TolF? TolR?>"
               << " <-Presidual pr?> <-Pmin pmin?> <-honorTolR 0|1>"
               << " <-maxSubsteps n?>" << endln;
        return 0;
    }

    int    tag;
    double dData[18];
    double oData[5];

    oData[0] = 1;          // IntScheme   ) slots identical to
    // Ladruno (ADR-86b): TanType default 0 -> 2. This is the ONE positional
    // default that differs from OPS_ManzariDafaliasMaterial, and it differs on
    // purpose -- see the long note at the head of this parser. 0 = mCe
    // (ELASTIC), 1 = mCep, 2 = mCep_Consistent.
    oData[1] = 2;          // TanType     ) OPS_ManzariDafaliasMaterial's,
    oData[2] = 1;          // JacoType    ) and all but TanType carry the
    oData[3] = 1.0e-7;     // TolF        ) same DEFAULT too
    oData[4] = 1.0e-7;     // TolR        )

    double presidual = 0.0;    // Ladruno: default -- a cohesionless sand has no cohesion
    double pmin      = -1.0;   // Ladruno: sentinel -- resolve to 1.0e-3 * P_atm in the ctor
    int    honorTolR = 0;      // Ladruno: default = vanilla's hardcoded ModifiedEuler
                               //          substep tolerance 1e-4 (ADR-86 PR-3)
    int    maxSubsteps = 0;    // Ladruno (ADR-86b): 0 = UNCAPPED = vanilla, so every
                               //          existing deck is byte-identical

    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "WARNING invalid nDMaterial LadrunoSANISAND material tag" << endln;
        return 0;
    }

    numData = 18;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid material data for nDMaterial LadrunoSANISAND material with tag: "
               << tag << endln;
        return 0;
    }

    // --- remaining tokens: positional optionals first, then `-flag` keywords ---
    int  nPos     = 0;
    bool seenFlag = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char *rawArg = OPS_GetString();   // consumes exactly one token in every backend
        if (rawArg == 0)                        // classic-Tcl backend returns 0 when exhausted
            break;

        // Copy immediately: the openseespy backend returns a pointer into a
        // temporary that it has already Py_DECREF'd, and it returns the literal
        // "Invalid String Input!" (having still consumed the token) whenever the
        // argument is a number rather than a string.
        char argTok[64];
        int  ic = 0;
        while (ic < 63 && rawArg[ic] != '\0') { argTok[ic] = rawArg[ic]; ic++; }
        argTok[ic] = '\0';

        if (strcmp(argTok, "-Presidual") == 0 || strcmp(argTok, "-presidual") == 0) {
            seenFlag = true;
            numData  = 1;
            if (OPS_GetDoubleInput(&numData, &presidual) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -Presidual wants one value" << endln;
                return 0;
            }
            if (presidual < 0.0) {
                // A negative would collide with the -Pmin sentinel convention and
                // is physically meaningless (it is an apparent cohesion).
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -Presidual must be >= 0 (got " << presidual
                       << "). p_residual is an apparent cohesion c = p_r*tan(phi)." << endln;
                return 0;
            }
        }
        else if (strcmp(argTok, "-Pmin") == 0 || strcmp(argTok, "-pmin") == 0) {
            seenFlag = true;
            numData  = 1;
            if (OPS_GetDoubleInput(&numData, &pmin) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -Pmin wants one value" << endln;
                return 0;
            }
            if (pmin <= 0.0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -Pmin must be > 0 (got " << pmin
                       << "). Omit the flag to get the default 1.0e-3*P_atm." << endln;
                return 0;
            }
        }
        else if (strcmp(argTok, "-honorTolR") == 0 || strcmp(argTok, "-honortolr") == 0) {
            seenFlag = true;
            numData  = 1;
            if (OPS_GetIntInput(&numData, &honorTolR) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -honorTolR wants 0 or 1" << endln;
                return 0;
            }
            if (honorTolR != 0 && honorTolR != 1) {
                // Ladruno (ADR-86 PR-3): the seam is wired now, but it is a BOOLEAN
                // seam -- the base member it drives is a `bool`, so 2 and -1 would
                // both silently mean 1. Refuse them rather than quietly widening.
                // (PR-1 refused every value but 0, because the seam did not exist
                // yet and a flag that claims to have done something it did not do is
                // the exact defect class this ADR exists to fix. PR-2 opened the seam
                // in vanilla; PR-3 connects it, so 1 is now a real request.)
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -honorTolR wants 0 or 1 (got " << honorTolR
                       << "). 0 = vanilla's hardcoded ModifiedEuler substep tolerance"
                       << " 1e-4; 1 = honour this deck's TolR." << endln;
                return 0;
            }
        }
        else if (strcmp(argTok, "-maxSubsteps") == 0 || strcmp(argTok, "-maxsubsteps") == 0) {
            // Ladruno (ADR-86b): the substep-COUNT cap on ModifiedEuler.
            seenFlag = true;
            numData  = 1;
            if (OPS_GetIntInput(&numData, &maxSubsteps) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -maxSubsteps wants one non-negative integer" << endln;
                return 0;
            }
            if (maxSubsteps < 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -maxSubsteps must be >= 0 (got " << maxSubsteps
                       << "). 0 means UNCAPPED, which is vanilla's behaviour." << endln;
                return 0;
            }
        }
        else {
            // Not one of our flags, so it must be a positional optional.
            if (seenFlag) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": positional optional '" << argTok << "' may not follow a -flag."
                       << " Order is: <IntScheme TanType JacoType TolF TolR> then -flags." << endln;
                return 0;
            }
            if (nPos >= 5) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": too many positional optional arguments (max 5:"
                       << " IntScheme TanType JacoType TolF TolR), at '" << argTok << "'" << endln;
                return 0;
            }
            OPS_ResetCurrentInputArg(-1);   // rewind the token we just consumed
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &oData[nPos]) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": unrecognized option '" << argTok << "'."
                       << " Expected a numeric positional optional or one of"
                       << " -Presidual / -Pmin / -honorTolR / -maxSubsteps" << endln;
                return 0;
            }
            nPos++;
        }
    }

    NDMaterial *theMaterial =
        new LadrunoSANISAND(tag, ND_TAG_LadrunoSANISAND,
                            dData[0],  dData[1],  dData[2],  dData[3],  dData[4],  dData[5],
                            dData[6],  dData[7],  dData[8],  dData[9],  dData[10], dData[11],
                            dData[12], dData[13], dData[14], dData[15], dData[16], dData[17],
                            (int)oData[0], (int)oData[1], (int)oData[2], oData[3], oData[4],
                            presidual, pmin, honorTolR, maxSubsteps);                 // Ladruno

    if (theMaterial == 0)
        opserr << "WARNING ran out of memory for nDMaterial LadrunoSANISAND material with tag: "
               << tag << endln;

    return theMaterial;
}

// ===========================================================================
//  constructors / destructor
//
//  In every one of them the base constructor runs first and calls the BASE
//  `initialize()` (correct -- the base must run its own init), which sets
//  m_Pmin = 1.0e-4*P_atm and m_Presidual = 1.0e-2*P_atm. applyLadrunoConstants()
//  then wins the last write.
// ===========================================================================

LadrunoSANISAND::LadrunoSANISAND(int tag, int classTag, double G0, double nu, double e_init, double Mc,
    double c, double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch,
    double nb, double A0, double nd, double z_max, double cz, double mDen,
    int integrationScheme, int tangentType, int JacoType, double TolF, double TolR,
    double Presidual, double Pmin, int honorTolR, int maxSubsteps)
  : ManzariDafalias(tag, classTag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch,
                    nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType, JacoType, TolF, TolR),
    mPresidualInput(Presidual),
    mPminInput(Pmin),
    mHonorTolR(honorTolR),
    mMaxSubsteps(maxSubsteps)                                                         // Ladruno
{
    // Defensive input sanitising -- the parser already rejects these, but the
    // wrappers and getCopy() also reach this constructor.
    this->sanitiseLadrunoInputs(tag);   // Ladruno (ADR-86 PR-3)

    this->applyLadrunoConstants();
    this->echoLadrunoConstants();
}

LadrunoSANISAND::LadrunoSANISAND(int tag, double G0, double nu, double e_init, double Mc,
    double c, double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch,
    double nb, double A0, double nd, double z_max, double cz, double mDen,
    int integrationScheme, int tangentType, int JacoType, double TolF, double TolR,
    double Presidual, double Pmin, int honorTolR, int maxSubsteps)
  : ManzariDafalias(tag, ND_TAG_LadrunoSANISAND, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm,
                    m, h0, ch, nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType,
                    JacoType, TolF, TolR),
    mPresidualInput(Presidual),
    mPminInput(Pmin),
    mHonorTolR(honorTolR),
    mMaxSubsteps(maxSubsteps)                                                         // Ladruno
{
    this->sanitiseLadrunoInputs(tag);   // Ladruno (ADR-86 PR-3)

    this->applyLadrunoConstants();
    this->echoLadrunoConstants();
}

// specific-type null constructor -- used by the wrappers' null constructors and,
// through them, by FEM_ObjectBrokerAllClasses. P_atm is 0 here, so there is
// nothing meaningful to echo; recvSelf restores the real values.
LadrunoSANISAND::LadrunoSANISAND(int classTag)
  : ManzariDafalias(classTag),
    mPresidualInput(0.0),
    mPminInput(-1.0),
    mHonorTolR(0),
    mMaxSubsteps(0)                                                                   // Ladruno
{
    this->applyLadrunoConstants();
}

// null constructor (FEM_ObjectBroker)
LadrunoSANISAND::LadrunoSANISAND()
  : ManzariDafalias(ND_TAG_LadrunoSANISAND),
    mPresidualInput(0.0),
    mPminInput(-1.0),
    mHonorTolR(0),
    mMaxSubsteps(0)                                                                   // Ladruno
{
    this->applyLadrunoConstants();
}

LadrunoSANISAND::~LadrunoSANISAND()
{
}

// Ladruno (ADR-86 PR-3): input sanitising, in ONE place.
//
// Both full constructors used to carry byte-identical copies of these three
// checks, and PR-3 made that worse before it made it better -- adding the
// honorTolR check took the count from four duplicated blocks to six. The header's
// own design note praises "one place" for the base-side WRITES
// (applyLadrunoConstants); this is the same rule applied to the input side, which
// is where a future rule change would otherwise have to be made twice.
//
// The parser already rejects all three cases with a hard error, so nothing a DECK
// can write reaches these. They exist because the wrappers and getCopy(const
// char*) also reach this constructor, and a clone must not be able to smuggle in a
// value the parser would have refused.
void
LadrunoSANISAND::sanitiseLadrunoInputs(int tag)
{
    if (mPresidualInput < 0.0) {
        opserr << "WARNING LadrunoSANISAND tag " << tag << ": p_residual = " << mPresidualInput
               << " < 0 is meaningless; using 0." << endln;
        mPresidualInput = 0.0;
    }
    if (mPminInput == 0.0) {
        opserr << "WARNING LadrunoSANISAND tag " << tag
               << ": p_min = 0 disables the low-stress clamp; using the default 1.0e-3*P_atm." << endln;
        mPminInput = -1.0;   // back to the sentinel
    }
    if (mHonorTolR != 0 && mHonorTolR != 1) {
        // The base member this drives is a `bool`, so any nonzero would collapse to
        // 1 silently. Refuse rather than widen.
        opserr << "WARNING LadrunoSANISAND tag " << tag << ": honorTolR = " << mHonorTolR
               << " is not 0 or 1; using 0 (vanilla's hardcoded ModifiedEuler"
               << " substep tolerance 1e-4)." << endln;
        mHonorTolR = 0;
    }
    // Ladruno (ADR-86b): the base seam is an `int` substep COUNT, so a negative
    // value is not "a different cap", it is a nonsense one -- and because the
    // guard tests `> 0`, a negative would silently mean "uncapped" instead of
    // being refused. Refuse it here for the same reason the parser does.
    if (mMaxSubsteps < 0) {
        opserr << "WARNING LadrunoSANISAND tag " << tag << ": maxSubsteps = " << mMaxSubsteps
               << " < 0 is meaningless; using 0 (UNCAPPED, vanilla's behaviour)." << endln;
        mMaxSubsteps = 0;
    }
}

// ===========================================================================
//  the "win the last write" helper + the echo
// ===========================================================================

// Ladruno: the whole class, in four lines. Every base integrator reads
// m_Presidual / m_Pmin as protected data at run time, so re-asserting them after
// anything that may have recomputed them is sufficient -- no integrator override
// is needed, and none is wanted.
//
// PR-3 adds the third line. `mHonorTolRInME` is the flag seam PR-2 opened in
// vanilla (ManzariDafalias.h, read once in ModifiedEuler() as
// `TolE = mHonorTolRInME ? mTolR : 1e-4`); it is a protected base member set
// `false` in all four ManzariDafalias constructors and written nowhere else in
// that class, so vanilla stays bit-identical and THIS is the only writer. It
// belongs here rather than in the constructor bodies for the same reason the
// other two do: applyLadrunoConstants() is also reached from initialize() (hence
// revertToStart) and from recvSelf, and a seam re-asserted in only some of those
// places is a seam that silently reverts -- which is the exact failure mode the
// `revertToStart` override exists to prevent for m_Presidual.
//
// NOT set here, deliberately: the sibling seam `mUseCurrentVoidRatioInG`
// (ADR 86 sec.7.3). Wiring it moves a CALIBRATED quantity -- Gorini's
// G0 = 264.32 was fitted against the frozen m_e_init form -- so it stays open
// pending the decision recorded in 86_ladruno_sanisand_pr3_tripwire_memo.md.
void
LadrunoSANISAND::applyLadrunoConstants(void)
{
    m_Presidual     = mPresidualInput;
    m_Pmin          = (mPminInput < 0.0) ? 1.0e-3 * m_P_atm : mPminInput;
    mHonorTolRInME  = (mHonorTolR != 0);   // Ladruno (ADR-86 PR-3): the seam, wired
    mMaxSubstepsInME = mMaxSubsteps;       // Ladruno (ADR-86b): the substep-count cap
}

// Ladruno (ADR-86 PR-3): does this deck's IntScheme actually reach the code that
// reads the seam?
//
// `mHonorTolRInME` is read at EXACTLY ONE site: ManzariDafalias::ModifiedEuler().
// So on a scheme that never calls ModifiedEuler, `-honorTolR 1` is accepted,
// stored, echoed, wired -- and inert. That is the "a flag claims to have done
// something it did not do" defect this ADR was written about, so it is warned
// about rather than left for the user to discover.
//
// Traced through the dispatch in ManzariDafalias.cpp, and the answer is NOT the
// obvious one:
//   scheme  1  INT_ModifiedEuler  explicit_integrator -> ModifiedEuler      REACHES
//   scheme  0  INT_MAXENE_MFE     MaxEnergyInc        -> ModifiedEuler      REACHES
//   anything not in {0..9, 45}    explicit_integrator -> default:           REACHES
//   scheme  2  INT_BackwardEuler  integrate() branches to BackwardEuler_CPPM  no
//   scheme  3  INT_RungeKutta                         -> RungeKutta4          no
//   scheme  5  INT_ForwardEuler                       -> ForwardEuler         no
//   scheme 45  INT_RungeKutta45                       -> RungeKutta45         no  (it
//              already uses `TolE = mTolR` unconditionally -- the seam exists
//              because ModifiedEuler was the outlier that did NOT)
//   schemes 4, 6                  MaxEnergyInc -> ForwardEuler / RungeKutta4  no
//   schemes 7, 8, 9               MaxStrainInc -> ForwardEuler in BOTH switch
//              branches, including `case INT_MAXSTR_MFE` -- so 7 does NOT reach
//              ModifiedEuler despite its name. Read the switch, not the name.
bool
LadrunoSANISAND::schemeReachesModifiedEuler(void) const
{
    // mScheme is `char unsigned` in the base, so s is in [0, 255] -- there is no
    // negative case to test for, and writing one would be dead code that reads as
    // a guard.
    const int s = (int)mScheme;
    if (s == INT_LSANISAND_ModifiedEuler || s == INT_LSANISAND_MAXENE_MFE)
        return true;
    // Scheme numbers the base's switch does not enumerate fall through
    // explicit_integrator's `default:`, which is ModifiedEuler.
    if (s > 9 && s != INT_LSANISAND_RungeKutta45)
        return true;
    return false;
}

// Ladruno: PER-DECK-MATERIAL echo (ADR 86 section 4.4).
//
// Deliberately NOT latched behind a `static bool printedOnce` -- a latched
// message is observable only by whichever test in a session runs first, which is
// exactly how the p_residual defect stayed invisible. Every `nDMaterial
// LadrunoSANISAND ...` command echoes, every time, in every process.
//
// But NOT once per Gauss point either. `getCopy(const char*)` runs a full
// constructor for EVERY integration point, so echoing unconditionally put one
// 207-byte line per Gauss point on stderr: measured at 10 lines for a single
// stdBrick and 514 for 64 of them, i.e. ~400k lines (~83 MB) before the first
// analysis step on a 50k-element model. That is not what section 4.4 asks for --
// its requirement is that the message not be LATCHED, and echoing once per deck
// command satisfies that in full.
//
// The discriminator is the class tag, which needs no signature change: the
// object the parser builds carries ND_TAG_LadrunoSANISAND, while every
// Gauss-point copy is a wrapper (ND_TAG_LadrunoSANISAND{3D,PlaneStrain}). So the
// deck-level material announces itself and its copies stay quiet. `Print()` is
// unguarded and still reports the constants for ANY instance on demand, so a
// per-Gauss-point record is still obtainable when someone actually wants one.
void
LadrunoSANISAND::echoLadrunoConstants(void)
{
    if (this->getClassTag() != ND_TAG_LadrunoSANISAND)   // Ladruno: copies stay quiet
        return;

    opserr << "LadrunoSANISAND tag " << this->getTag()
           << ": p_residual = " << m_Presidual;
    if (mPresidualInput == 0.0)
        opserr << " (default, cohesionless)";
    else
        opserr << " (user)";

    opserr << ", p_min = " << m_Pmin;
    if (mPminInput < 0.0)
        opserr << " (default = 1e-3*P_atm, P_atm = " << m_P_atm << ")";
    else
        opserr << " (user)";

    // Ladruno (ADR-86 PR-3): the honoured ModifiedEuler substep error tolerance.
    // Named as a NUMBER, not as a flag state -- "honorTolR = 1" tells the reader
    // what was asked for, "TolE = 1e-06" tells them what the integrator ran.
    opserr << ", honorTolR = " << mHonorTolR
           << " (ModifiedEuler substep TolE = "
           << (mHonorTolR ? mTolR : 1.0e-4)
           << (mHonorTolR ? ", this deck's TolR" : ", vanilla's hardcoded value")
           << ")";

    // Ladruno (ADR-86b): the TANGENT the deck will actually run. This is echoed
    // because the parser default MOVED (0 -> 2) and a silent default change is
    // exactly the defect class this class exists to make impossible. Named as a
    // tangent, not as a number, so the reader does not have to remember the map.
    opserr << ", TanType = " << (int)mTangType << " ("
           << ((int)mTangType == 0 ? "ELASTIC mCe -- Newton degrades to modified Newton"
              : ((int)mTangType == 1 ? "continuum mCep" : "consistent mCep_Consistent (unsymmetric)"))
           << ")";

    // Ladruno (ADR-86b): the substep-count cap. Say the NUMBER and say what
    // happens when it fires -- "0" alone reads like "no substeps".
    opserr << ", maxSubsteps = " << mMaxSubsteps;
    if (mMaxSubsteps == 0)
        opserr << " (UNCAPPED, vanilla: dT_min = 1e-6 is the only bound)";
    else
        opserr << " (ModifiedEuler FAILS the update past this many substeps,"
                  " so the step can be cut)";

    // ADR 86 D5a (still open) / D5b (repaired in vanilla by PR-2).
    // The PR-1 text here read "D_factor ... ships UNCHANGED and is still
    // kPa-dimensional in this PR". PR-2 non-dimensionalised the sigmoid at all
    // four sites, so that sentence became FALSE the moment PR-2 merged and stayed
    // in the echo -- a stale claim in exactly the channel this class exists to
    // make trustworthy. Corrected in PR-3.
    // NB the wording avoids the literal token `p_residual = <value>`: the battery
    // counts echo LINES by that signature to prove the echo is not latched, and
    // an earlier draft of this sentence containing "once p_residual = 0" doubled
    // that count. Prose in a machine-read stream has to stay out of the machine's
    // way.
    opserr << " [D5b: D_factor sigmoid non-dimensionalised in vanilla (PR-2),"
              " an exact no-op at P_atm = 101; D5a: whether it should exist at all"
              " once the residual pressure is zero is still OPEN]";

    opserr << endln;

    // Ladruno (ADR-86 PR-3): -honorTolR 1 on a scheme that never calls
    // ModifiedEuler() is accepted, stored and wired -- and does nothing. Say so
    // at construction rather than letting the deck's author infer it from a
    // runtime that did not change.
    if (mHonorTolR != 0 && !this->schemeReachesModifiedEuler()) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -honorTolR 1 has NO EFFECT with IntScheme " << (int)mScheme
               << ". The seam it sets (ManzariDafalias mHonorTolRInME) is read at"
               << " exactly one site, inside ModifiedEuler(), and this scheme does"
               << " not route there."
               << (((int)mScheme == INT_LSANISAND_RungeKutta45)
                     ? " IntScheme 45 (RungeKutta45) already honours TolR"
                       " unconditionally -- ModifiedEuler was the outlier that did"
                       " not, which is why the seam exists."
                     : " Use IntScheme 1 (ModifiedEuler) if you want it.")
               << endln;
    }

    // Ladruno (ADR-86b): -maxSubsteps reads the SAME single site as -honorTolR
    // (inside ModifiedEuler), so it has exactly the same inertness hazard and
    // gets exactly the same warning. Kept as a separate `if` rather than folded
    // into the one above so that a deck setting only one of the two flags is
    // told about the one it set.
    if (mMaxSubsteps != 0 && !this->schemeReachesModifiedEuler()) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -maxSubsteps " << mMaxSubsteps << " has NO EFFECT with IntScheme "
               << (int)mScheme << ". The seam it sets (ManzariDafalias mMaxSubstepsInME)"
               << " is read at exactly one site, inside ModifiedEuler(), and this scheme"
               << " does not route there. Use IntScheme 1 (ModifiedEuler) if you want it."
               << endln;
    }
}

// SHADOW of ManzariDafalias::initialize() -- non-virtual, same signature, on
// purpose. See the DESIGN NOTE in LadrunoSANISAND.h. This is what stops the base
// recomputing p_r from P_atm anywhere the derived type is statically known.
void
LadrunoSANISAND::initialize()
{
    ManzariDafalias::initialize();   // the base must still run its own init
    this->applyLadrunoConstants();   // Ladruno: and we take the last write
}

// ===========================================================================
//  revertToStart -- replicates ManzariDafalias::revertToStart (:465-476) but
//  routes through THIS class's initialize().
//
//  Without this override the base version runs, and its `this->initialize()`
//  binds statically to ManzariDafalias::initialize() inside the base TU, which
//  restores m_Presidual = 1.0e-2*P_atm MID-ANALYSIS. It looks like the most
//  skippable of the six overrides. It is not.
// ===========================================================================
int
LadrunoSANISAND::revertToStart(void)
{
    // added: C.McGann, U.Washington for InitialStateAnalysis
    if (ops_InitialStateAnalysis) {
        // do nothing, keep state variables from last step
    } else {
        // normal call for revertToStart (not initialStateAnalysis)
        this->initialize();          // Ladruno: LadrunoSANISAND::initialize(), not the base's
    }

    return 0;
}

// ===========================================================================
//  getCopy
// ===========================================================================

NDMaterial *
LadrunoSANISAND::getCopy(void)
{
    // Mirrors ManzariDafalias::getCopy(void): the dimension-specific subclass is
    // responsible. LadrunoSANISAND3D / LadrunoSANISANDPlaneStrain override this.
    opserr << "LadrunoSANISAND::getCopy -- subclass responsibility\n";
    exit(-1);
    return 0;
}

NDMaterial *
LadrunoSANISAND::getCopy(const char *type)
{
    // Ladruno: the base version hardcodes `new ManzariDafaliasPlaneStrain(...)`
    // and `new ManzariDafalias3D(...)` (ManzariDafalias.cpp:392-410), so
    // inheriting it would silently strip p_residual/p_min at EVERY Gauss point --
    // the deck would echo the settings once and then run vanilla everywhere.
    if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
        LadrunoSANISANDPlaneStrain *clone;
        clone = new LadrunoSANISANDPlaneStrain(this->getTag(), m_G0, m_nu, m_e_init, m_Mc,
                       m_c, m_lambda_c, m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0,
                       m_nd, m_z_max, m_cz, massDen, mScheme, mTangType, mJacoType, mTolF, mTolR,
                       mPresidualInput, mPminInput, mHonorTolR, mMaxSubsteps);        // Ladruno
        return clone;
    } else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
        LadrunoSANISAND3D *clone;
        clone = new LadrunoSANISAND3D(this->getTag(), m_G0, m_nu, m_e_init, m_Mc, m_c, m_lambda_c,
                       m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz,
                       massDen, mScheme, mTangType, mJacoType, mTolF, mTolR,
                       mPresidualInput, mPminInput, mHonorTolR, mMaxSubsteps);        // Ladruno
        return clone;
    } else {
        opserr << "LadrunoSANISAND::getCopy failed to get copy: " << type << endln;
        return 0;
    }
}

// ===========================================================================
//  sendSelf / recvSelf
//
//  The base Vector(97) wire format is UNCHANGED -- mixed-build parallel and
//  database restore against vanilla ManzariDafalias keep working, and the class
//  tag already separates the two materials. One extra Vector(4) follows it on the
//  same channel:
//
//      data(0) = mPresidualInput   (as given by the user)
//      data(1) = mPminInput        (as given; < 0 = "1e-3*P_atm" sentinel)
//      data(2) = m_Presidual       (the RESOLVED value actually in force)
//      data(3) = (double)mHonorTolR
//      data(4) = (double)mMaxSubsteps   (Ladruno, ADR-86b)
//
//  ADR-86b widened this Vector(4) to a Vector(5). Both ends of any channel are
//  the SAME BUILD in every supported workflow (an MP job runs one binary; a
//  FileDatastore is read back by the engine that wrote it), so this is a
//  fork-internal format and no compatibility shim is owed. What it DOES break is
//  reading a datastore written by a pre-ADR-86b build -- that is a size mismatch
//  on a different FileDatastore file, i.e. a loud failure, not a silent one.
//
//  data(2) is redundant with data(0) today. It is sent anyway so a restored
//  material can be compared against what the sender was actually running without
//  re-deriving it -- that comparison is the whole point of ADR 86 section 3, in
//  which serial and MP were found to be running different soils with nothing
//  warning.
//
//  Vectors of different sizes go to different FileDatastore files
//  (FileDatastore.cpp:823-856 keys the file by Vector::Size()), so a Vector(4)
//  cannot collide with the base's Vector(97) under the same dbTag/commitTag.
//  On a stream channel (MPI) the two sends are simply ordered, and recvSelf
//  reads them back in the same order.
// ===========================================================================

int
LadrunoSANISAND::sendSelf(int commitTag, Channel &theChannel)
{
    int res = ManzariDafalias::sendSelf(commitTag, theChannel);
    if (res < 0) {
        opserr << "WARNING: LadrunoSANISAND::sendSelf - base class failed to send" << endln;
        return -1;
    }

    static Vector ladrunoData(5);                                                     // Ladruno

    ladrunoData(0) = mPresidualInput;
    ladrunoData(1) = mPminInput;
    ladrunoData(2) = m_Presidual;
    ladrunoData(3) = (double)mHonorTolR;
    ladrunoData(4) = (double)mMaxSubsteps;   // Ladruno (ADR-86b)

    res = theChannel.sendVector(this->getDbTag(), commitTag, ladrunoData);
    if (res < 0) {
        opserr << "WARNING: LadrunoSANISAND::sendSelf - failed to send Ladruno constants"
               << " to channel" << endln;
        return -1;
    }

    return 0;
}

int
LadrunoSANISAND::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = ManzariDafalias::recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "WARNING: LadrunoSANISAND::recvSelf - base class failed to receive" << endln;
        return -1;
    }

    static Vector ladrunoData(5);                                                     // Ladruno

    res = theChannel.recvVector(this->getDbTag(), commitTag, ladrunoData);
    if (res < 0) {
        opserr << "WARNING: LadrunoSANISAND::recvSelf - failed to receive Ladruno constants"
               << " from channel" << endln;
        return -1;
    }

    mPresidualInput = ladrunoData(0);
    mPminInput      = ladrunoData(1);
    m_Presidual     = ladrunoData(2);   // overwritten by applyLadrunoConstants below;
                                        // restored first so a future divergence is visible
    mHonorTolR      = (int)ladrunoData(3);
    mMaxSubsteps    = (int)ladrunoData(4);   // Ladruno (ADR-86b)

    // The base recvSelf restored m_Pmin from its own data(96) and never re-runs
    // initialize(); we take the last write here.
    this->applyLadrunoConstants();                                                    // Ladruno

    return 0;
}

// ===========================================================================
//  Ladruno (ADR-86b): the substep-cap status, and its diagnostic response
// ===========================================================================

// The rule the wrappers' setTrialStrain returns, in ONE place.
//
// `mSubstepCapHitInME` is reset at the top of every ManzariDafalias::integrate()
// and set only by ModifiedEuler()'s cap, so it describes THIS material update
// and no earlier one. When it is set the strain increment was not integrated:
// the trial state is partially updated (harmless -- nothing reads it after a
// failed update) and the COMMITTED state is untouched, because integrate()
// writes only trial members.
//
// LADRUNO_MATERIAL_REFUSED, not a bare -1, and the distinction is load-bearing:
// see the long note in SRC/material/LadrunoMaterialStatus.h. In one line --
// `ASDConcrete3DMaterial` already returns a negative code to mean "an inner
// iteration missed, here is my best state", and the fork's own rule (ADR-33/34,
// LEDGER_quirks) is that such a code must NOT fail the step, because doing so
// makes softening analyses fragile. Measured: propagating any `< 0` in
// LadrunoBrick killed two long-green ASDConcrete mesh-objectivity gates. This
// code means something else -- "the increment was not integrated at all" -- so
// it gets its own value and the element propagates only that.
//
// NOT applied to schemes that never call ModifiedEuler: on those the flag can
// never be set, so this is 0 by construction rather than by a branch.
int
LadrunoSANISAND::ladrunoUpdateStatus(void) const
{
    return mSubstepCapHitInME ? LADRUNO_MATERIAL_REFUSED : 0;
}

// Ladruno (ADR-86b): "substeps" / "substepsME" -- what the last update cost.
//
// Cheap by construction: two ints already maintained by the integrator, read
// through the ordinary MaterialResponse path, so nothing is computed for a deck
// that does not ask. Defined on the BASE Ladruno class, so LadrunoSANISAND3D and
// LadrunoSANISANDPlaneStrain (neither of which overrides setResponse) inherit
// it, and anything else falls through to ManzariDafalias::setResponse unchanged.
//
// The id only has to miss the base's 1..8; it is not a class tag and nothing may
// derive one from it. Named for what it reports rather than for its digits, and
// `constexpr` rather than a macro for the same reasons as the refusal code.
constexpr int LadrunoSanisandSubstepResponseID = 33086;   // Ladruno (ADR-86b)

Response *
LadrunoSANISAND::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc > 0 && (strcmp(argv[0], "substeps") == 0 ||
                     strcmp(argv[0], "substepsME") == 0 ||
                     strcmp(argv[0], "ladrunoSubsteps") == 0)) {
        static Vector probe(2);
        return new MaterialResponse(this, LadrunoSanisandSubstepResponseID, probe);
    }
    return ManzariDafalias::setResponse(argv, argc, output);
}

int
LadrunoSANISAND::getResponse(int responseID, Information &matInformation)
{
    if (responseID == LadrunoSanisandSubstepResponseID) {
        static Vector out(2);
        out(0) = (double)mSubstepsTakenInME;   // substeps spent in the LAST update
        out(1) = mSubstepCapHitInME ? 1.0 : 0.0;
        return matInformation.setVector(out);
    }
    return ManzariDafalias::getResponse(responseID, matInformation);
}

// ===========================================================================
//  Print -- a record that cannot state what it ran is the thing this override
//  exists to prevent (ADR 86 section 4.4, and D5a/D5b in section 7.2 / 7.2.1).
// ===========================================================================
void
LadrunoSANISAND::Print(OPS_Stream &s, int flag)
{
    s << "LadrunoSANISAND Material, tag: " << this->getTag() << endln;

    // ManzariDafalias::getType() is "subclass responsibility" and calls exit(-1),
    // so only ask for it on a dimension-specific instance. The object the material
    // registry holds carries the base tag and has no type.
    if (this->getClassTag() == ND_TAG_LadrunoSANISAND)
        s << "  Type: (dimensionless prototype; getCopy(\"3D\"|\"PlaneStrain\") makes the working copy)" << endln;
    else
        s << "  Type: " << this->getType() << endln;

    s << "  Base model: ManzariDafalias (Dafalias & Manzari 2004; A.Ghofrani, P.Arduino, UW)" << endln;
    s << "  p_residual = " << m_Presidual
      << "   (input " << mPresidualInput
      << (mPresidualInput == 0.0 ? ", default: cohesionless)" : ", user)") << endln;
    s << "  p_min      = " << m_Pmin;
    if (mPminInput < 0.0)
        s << "   (default = 1e-3*P_atm, P_atm = " << m_P_atm << ")" << endln;
    else
        s << "   (input " << mPminInput << ", user)" << endln;
    s << "  P_atm      = " << m_P_atm << ",  IntScheme = " << (int)mScheme
      << ",  TanType = " << (int)mTangType << ",  JacoType = " << (int)mJacoType << endln;
    s << "  TolF       = " << mTolF << ",  TolR = " << mTolR
      << ",  honorTolR = " << mHonorTolR << endln;
    // Ladruno (ADR-86b): the substep cap and the state it left behind.
    s << "  maxSubsteps = " << mMaxSubsteps
      << (mMaxSubsteps == 0 ? "  (UNCAPPED = vanilla; ModifiedEuler is bounded only by"
                              " dT_min = 1e-6, i.e. up to 1e6 substeps per update)"
                            : "  (ModifiedEuler returns FAILURE past this count instead of"
                              " force-accepting; the committed state is left untouched)") << endln;
    s << "             last update: " << mSubstepsTakenInME << " ModifiedEuler substep(s)"
      << (mSubstepCapHitInME ? ", CAP HIT (that update did not integrate)" : "") << endln;
    // Ladruno (ADR-86b): the same inertness note the -honorTolR block below carries.
    // Both flags drive seams read at EXACTLY ONE site, inside ModifiedEuler(), so on
    // a scheme that never routes there the cap is stored, echoed, wired -- and does
    // nothing. NB IntScheme 7 is called INT_MAXSTR_MFE and does NOT reach
    // ModifiedEuler: MaxStrainInc has no case for it and falls through to
    // ForwardEuler (ManzariDafalias.cpp:1199-1207). Read the switch, not the name;
    // schemeReachesModifiedEuler() encodes the switch and is correct as written.
    if (mMaxSubsteps != 0 && !this->schemeReachesModifiedEuler())
        s << "             NOTE: IntScheme " << (int)mScheme << " does not route to"
             " ModifiedEuler(), so -maxSubsteps is INERT on this deck." << endln;
    // A record that also states what the cap would COST if it fired: the element
    // must propagate the refusal. LadrunoBrick does; Brick / BrickUP / QuadUP /
    // stdBrick discard it, and under those a capped run is invalid, not merely
    // un-cut. The material cannot see its element, so this is a statement, not a
    // check.
    if (mMaxSubsteps != 0)
        s << "             NOTE: a cap is only safe under an element that PROPAGATES"
             " a material refusal (today: LadrunoBrick)." << endln;
    // Ladruno (ADR-86 PR-3): the seam is wired now (it was "inactive in PR-1").
    // Report the tolerance the integrator ACTUALLY ran with, and whether this
    // deck's scheme even reaches the site that reads it -- a record that says
    // "honorTolR = 1" while the scheme never calls ModifiedEuler is a record that
    // overstates what happened.
    s << "             ModifiedEuler substep TolE = " << (mHonorTolR ? mTolR : 1.0e-4)
      << (mHonorTolR ? "  (this deck's TolR, via the ManzariDafalias mHonorTolRInME seam)"
                     : "  (vanilla's hardcoded 1e-4; -honorTolR 1 selects TolR instead)") << endln;
    if (mHonorTolR != 0 && !this->schemeReachesModifiedEuler())
        s << "             NOTE: IntScheme " << (int)mScheme << " does not route to"
             " ModifiedEuler(), so -honorTolR 1 is INERT on this deck." << endln;

    // ADR 86 section 7.2 / 7.2.1.  D5b is DONE (vanilla, PR-2); D5a is open.
    // The PR-1 text this replaces asserted the sigmoid "ships UNCHANGED" and "is
    // still kPa-DIMENSIONAL in this PR", then contradicted itself two lines later
    // by reporting the PR-2 repair.  Both halves cannot be true after PR-2 merged.
    s << "  ADR 86 D5b (DONE, vanilla, PR-2): the D_factor dilatancy sigmoid is" << endln;
    s << "             non-dimensionalised at all four sites. For p < 0.05*P_atm the model" << endln;
    s << "             applies D_factor = 1/(1 + exp(7.6349 - 7.2713*101.0/P_atm * p))." << endln;
    s << "             The shipped form multiplied a RAW stress by the bare literal 7.2713," << endln;
    s << "             which therefore carried units of 1/stress: at a true confinement of" << endln;
    s << "             1 kPa the factor was 0.410 in kPa, 1.000 in Pa (never fires) and" << endln;
    s << "             0.0005 in MPa (total suppression). The repair is an EXACT no-op" << endln;
    s << "             wherever P_atm = 101 (kPa); it changes only Pa/MPa decks, and a deck" << endln;
    s << "             declaring P_atm = 101.3 shifts ~1.3% in D_factor at 1 kPa." << endln;
    s << "  ADR 86 D5a (OPEN): whether the sigmoid should exist AT ALL once p_residual = 0" << endln;
    s << "             is a separate modelling question and is NOT settled by the units" << endln;
    s << "             repair. This class does not change its shape. Half-suppression sits" << endln;
    s << "             at p = 7.6349/7.2713 = 1.050 kPa, within 4% of vanilla's p_residual" << endln;
    s << "             of 1.01 kPa -- see ADR 86 section 7.2.1 and the PR-3 tripwire memo." << endln;
}
