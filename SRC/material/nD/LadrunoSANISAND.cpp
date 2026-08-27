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

#include <string.h>
#include <math.h>

// ===========================================================================
//  OPS parser
//
//    nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi \
//        $P_atm $m $h0 $ch $nb $A0 $nd $z_max $cz $Rho                         \
//        <$IntScheme $TanType $JacoType $TolF $TolR>                           \
//        <-Presidual $pr> <-Pmin $pmin> <-honorTolR 0|1>
//
//  The first 18 positional doubles and the 5 positional optionals are IDENTICAL
//  to `nDMaterial ManzariDafalias`, so a deck migrates by renaming the command.
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
               << " <-Presidual pr?> <-Pmin pmin?> <-honorTolR 0|1>" << endln;
        return 0;
    }

    int    tag;
    double dData[18];
    double oData[5];

    oData[0] = 1;          // IntScheme   ) defaults identical to
    oData[1] = 0;          // TanType     ) OPS_ManzariDafaliasMaterial,
    oData[2] = 1;          // JacoType    ) so a renamed deck behaves
    oData[3] = 1.0e-7;     // TolF        ) identically
    oData[4] = 1.0e-7;     // TolR        )

    double presidual = 0.0;    // Ladruno: default -- a cohesionless sand has no cohesion
    double pmin      = -1.0;   // Ladruno: sentinel -- resolve to 1.0e-3 * P_atm in the ctor
    int    honorTolR = 0;      // Ladruno: PR-2 seam

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
            if (honorTolR != 0) {
                // Ladruno: PR-2 seam. The flag would be one line at
                // ManzariDafalias.cpp:1320 (`TolE = mHonorTolRInME ? mTolR : 1e-4`),
                // and PR-1 does not edit vanilla ManzariDafalias (ADR 86 D3/D7).
                // Failing loudly beats accepting a flag that does nothing -- a flag
                // that claims to have done something it did not do is exactly the
                // class of defect this ADR exists to fix.
                opserr << "ERROR nDMaterial LadrunoSANISAND tag " << tag
                       << ": -honorTolR " << honorTolR << " is NOT WIRED YET."
                       << " ADR 86 PR-2 opened the base flag seam (ManzariDafalias mHonorTolRInME,"
                       << " read in ModifiedEuler), but nothing sets it from this class or this"
                       << " parser -- that wiring is PR-3. Accepting the flag now would claim an"
                       << " error tolerance the integrator is not honouring, which is the exact"
                       << " defect class ADR 86 exists to fix. Use -honorTolR 0, or omit the flag."
                       << endln;
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
                       << " -Presidual / -Pmin / -honorTolR" << endln;
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
                            presidual, pmin, honorTolR);                              // Ladruno

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
    double Presidual, double Pmin, int honorTolR)
  : ManzariDafalias(tag, classTag, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch,
                    nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType, JacoType, TolF, TolR),
    mPresidualInput(Presidual),
    mPminInput(Pmin),
    mHonorTolR(honorTolR)                                                             // Ladruno
{
    // Defensive input sanitising -- the parser already rejects these, but the
    // wrappers and getCopy() also reach this constructor.
    if (mPresidualInput < 0.0) {                                                      // Ladruno
        opserr << "WARNING LadrunoSANISAND tag " << tag << ": p_residual = " << mPresidualInput
               << " < 0 is meaningless; using 0." << endln;
        mPresidualInput = 0.0;
    }
    if (mPminInput == 0.0) {                                                          // Ladruno
        opserr << "WARNING LadrunoSANISAND tag " << tag
               << ": p_min = 0 disables the low-stress clamp; using the default 1.0e-3*P_atm." << endln;
        mPminInput = -1.0;   // back to the sentinel
    }
    if (mHonorTolR != 0) {                                                            // Ladruno: PR-2 seam
        opserr << "ERROR LadrunoSANISAND tag " << tag << ": honorTolR = " << mHonorTolR
               << " is not implemented in PR-1 (ADR 86 D7 puts the ModifiedEuler mTolR"
               << " seam in PR-2); forcing 0." << endln;
        mHonorTolR = 0;
    }

    this->applyLadrunoConstants();
    this->echoLadrunoConstants();
}

LadrunoSANISAND::LadrunoSANISAND(int tag, double G0, double nu, double e_init, double Mc,
    double c, double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch,
    double nb, double A0, double nd, double z_max, double cz, double mDen,
    int integrationScheme, int tangentType, int JacoType, double TolF, double TolR,
    double Presidual, double Pmin, int honorTolR)
  : ManzariDafalias(tag, ND_TAG_LadrunoSANISAND, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm,
                    m, h0, ch, nb, A0, nd, z_max, cz, mDen, integrationScheme, tangentType,
                    JacoType, TolF, TolR),
    mPresidualInput(Presidual),
    mPminInput(Pmin),
    mHonorTolR(honorTolR)                                                             // Ladruno
{
    if (mPresidualInput < 0.0) {                                                      // Ladruno
        opserr << "WARNING LadrunoSANISAND tag " << tag << ": p_residual = " << mPresidualInput
               << " < 0 is meaningless; using 0." << endln;
        mPresidualInput = 0.0;
    }
    if (mPminInput == 0.0) {                                                          // Ladruno
        opserr << "WARNING LadrunoSANISAND tag " << tag
               << ": p_min = 0 disables the low-stress clamp; using the default 1.0e-3*P_atm." << endln;
        mPminInput = -1.0;
    }
    if (mHonorTolR != 0) {                                                            // Ladruno: PR-2 seam
        opserr << "ERROR LadrunoSANISAND tag " << tag << ": honorTolR = " << mHonorTolR
               << " is not implemented in PR-1 (ADR 86 D7 puts the ModifiedEuler mTolR"
               << " seam in PR-2); forcing 0." << endln;
        mHonorTolR = 0;
    }

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
    mHonorTolR(0)                                                                     // Ladruno
{
    this->applyLadrunoConstants();
}

// null constructor (FEM_ObjectBroker)
LadrunoSANISAND::LadrunoSANISAND()
  : ManzariDafalias(ND_TAG_LadrunoSANISAND),
    mPresidualInput(0.0),
    mPminInput(-1.0),
    mHonorTolR(0)                                                                     // Ladruno
{
    this->applyLadrunoConstants();
}

LadrunoSANISAND::~LadrunoSANISAND()
{
}

// ===========================================================================
//  the "win the last write" helper + the echo
// ===========================================================================

// Ladruno: the whole class, in three lines. Every base integrator reads
// m_Presidual / m_Pmin as protected data at run time, so re-asserting them after
// anything that may have recomputed them is sufficient -- no integrator override
// is needed, and none is wanted.
void
LadrunoSANISAND::applyLadrunoConstants(void)
{
    m_Presidual = mPresidualInput;
    m_Pmin      = (mPminInput < 0.0) ? 1.0e-3 * m_P_atm : mPminInput;
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

    // ADR 86 D5a/D5b (amended 2026-08-26): declared, not fixed, in PR-1.
    opserr << " [D5a/D5b: D_factor dilatancy sigmoid ships UNCHANGED and is still"
              " kPa-dimensional in this PR]"
           << endln;
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
                       mPresidualInput, mPminInput, mHonorTolR);                      // Ladruno
        return clone;
    } else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
        LadrunoSANISAND3D *clone;
        clone = new LadrunoSANISAND3D(this->getTag(), m_G0, m_nu, m_e_init, m_Mc, m_c, m_lambda_c,
                       m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz,
                       massDen, mScheme, mTangType, mJacoType, mTolF, mTolR,
                       mPresidualInput, mPminInput, mHonorTolR);                      // Ladruno
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

    static Vector ladrunoData(4);                                                     // Ladruno

    ladrunoData(0) = mPresidualInput;
    ladrunoData(1) = mPminInput;
    ladrunoData(2) = m_Presidual;
    ladrunoData(3) = (double)mHonorTolR;

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

    static Vector ladrunoData(4);                                                     // Ladruno

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

    // The base recvSelf restored m_Pmin from its own data(96) and never re-runs
    // initialize(); we take the last write here.
    this->applyLadrunoConstants();                                                    // Ladruno

    return 0;
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
      << ",  honorTolR = " << mHonorTolR << " (ADR 86 PR-2 seam, inactive in PR-1)" << endln;

    // ADR 86 D5a/D5b / section 7.2, 7.2.1 (amended 2026-08-26) -- declared here,
    // deliberately not fixed in PR-1.
    s << "  ADR 86 D5a/D5b: the D_factor dilatancy sigmoid ships UNCHANGED from ManzariDafalias." << endln;
    s << "             It is still kPa-DIMENSIONAL in this PR: for p < 0.05*P_atm the model" << endln;
    s << "             applies D_factor = 1/(1 + exp(7.6349 - 7.2713*p)) in which 7.2713" << endln;
    s << "             multiplies a raw stress and is a bare literal, so it carries units" << endln;
    s << "             of 1/stress. At a true confinement of 1 kPa the factor is 0.410 in kPa," << endln;
    s << "             1.000 in Pa (never fires) and 0.0005 in MPa (total suppression). These" << endln;
    s << "             results are as intended ONLY if this deck's stress unit is kPa." << endln;
    s << "             D5b: the UNIT dependence is repaired in vanilla ManzariDafalias in PR-2" << endln;
    s << "             (b = 7.2713*101.0 on p/P_atm, all four sites; a no-op at P_atm = 101)." << endln;
    s << "             D5a: whether the sigmoid should exist at all once p_residual = 0 is a" << endln;
    s << "             separate, still-open modelling question -- see ADR 86 section 7.2.1." << endln;
}
