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
#include <Parameter.h>                   // Ladruno (ADR-92 P1): setParameter/updateParameter
#include <LadrunoMaterialStatus.h>        // Ladruno (ADR-86b): LADRUNO_MATERIAL_REFUSED

#include <string.h>
#include <math.h>

// ===========================================================================
//  OPS parser
//
//    nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi \
//        $P_atm $m $h0 $ch $nb $A0 $nd $z_max $cz $Rho                         \
//        <$IntScheme $TanType $JacoType $TolF $TolR>                           \
//        <-Presidual $pr> <-Pmin $pmin> <-honorTolR 0|1> <-maxSubsteps $n>     \
//        <-implex> <-implexControl $tol $reductionLimit> <-implexAlpha $a>      \
//        <-implexDt pseudo|strain|user <$dt>>
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
               << " <-maxSubsteps n?>"
               << " <-implex> <-implexControl tol? reductionLimit?>"
               << " <-implexAlpha a?> <-implexDt pseudo|strain|user <dt?>>"     // Ladruno (ADR-92)
               << endln;
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

    // Ladruno (ADR-92 P1): every default here is "IMPL-EX off", which is what
    // makes an existing SANISAND deck byte-identical.
    LadrunoImplexOptions implexOpt;                                             // Ladruno (ADR-92)
    // Tracked explicitly rather than inferred from "is any field non-default":
    // `-implexAlpha 1.0` and `-implexDt pseudo` written WITHOUT -implex are
    // default-valued, and inferring would let exactly those two slip past the
    // check that refuses an option nothing will read.
    bool sawImplexToken = false;                                                // Ladruno (ADR-92)

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
        // ------------------------------------------------------------------
        //  Ladruno (ADR-92 P1): the IMPL-EX flags.
        //
        //  Parsed here, VALIDATED in setLadrunoImplexOptions() -- the checks
        //  need mScheme and mMaxSubsteps, which only exist once the material is
        //  built, and a clone must not be able to smuggle past a rule the deck
        //  was refused (the sanitiseLadrunoInputs argument, one level up).
        // ------------------------------------------------------------------
        else if (strcmp(argTok, "-implex") == 0 || strcmp(argTok, "-implEx") == 0 ||
                 strcmp(argTok, "-implEX") == 0) {
            seenFlag = true;
            sawImplexToken    = true;
            implexOpt.enabled = true;
        }
        else if (strcmp(argTok, "-implexControl") == 0 || strcmp(argTok, "-implexcontrol") == 0) {
            seenFlag = true;
            sawImplexToken    = true;
            implexOpt.control = true;
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &implexOpt.errorTol) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexControl wants two values: $tol $reductionLimit" << endln;
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &implexOpt.reductionLimit) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexControl wants two values: $tol $reductionLimit" << endln;
                return 0;
            }
            if (implexOpt.errorTol <= 0.0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexControl tolerance must be > 0 (got " << implexOpt.errorTol
                       << "). P0 measured the extrapolation error at 5e-4 (nominal increment)"
                       << " to 1.26 (p0 = 5 kPa at 5e-4); 0.05 is the ADR-92 working value."
                       << endln;
                return 0;
            }
            if (implexOpt.reductionLimit <= 0.0 || implexOpt.reductionLimit > 1.0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexControl reductionLimit must be in (0, 1] (got "
                       << implexOpt.reductionLimit
                       << "). It is the floor below which the material stops refusing,"
                          " as a fraction of the FIRST dt seen." << endln;
                return 0;
            }
        }
        else if (strcmp(argTok, "-implexAlpha") == 0 || strcmp(argTok, "-implexalpha") == 0) {
            seenFlag = true;
            sawImplexToken = true;
            numData  = 1;
            if (OPS_GetDoubleInput(&numData, &implexOpt.alpha) != 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexAlpha wants one value" << endln;
                return 0;
            }
            if (implexOpt.alpha < 0.0) {
                // A negative alpha extrapolates BACKWARDS along the plastic flow.
                // There is no reading of IMPL-EX in which that is what the deck
                // meant, and it would silently manufacture an unloading history.
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexAlpha must be >= 0 (got " << implexOpt.alpha
                       << "). It scales the extrapolation of the committed plastic-strain"
                          " increment; 1.0 is the standard IMPL-EX, 0.0 disables the"
                          " extrapolation and leaves a purely elastic predictor." << endln;
                return 0;
            }
        }
        else if (strcmp(argTok, "-implexDt") == 0 || strcmp(argTok, "-implexdt") == 0) {
            seenFlag = true;
            sawImplexToken = true;
            const char *rawMode = OPS_GetString();
            if (rawMode == 0) {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexDt wants one of pseudo|strain|user" << endln;
                return 0;
            }
            char modeTok[32];
            int  mc = 0;
            while (mc < 31 && rawMode[mc] != '\0') { modeTok[mc] = rawMode[mc]; mc++; }
            modeTok[mc] = '\0';

            if (strcmp(modeTok, "pseudo") == 0) {
                implexOpt.dtSource = LadrunoImplexOptions::DT_PSEUDO;
            }
            else if (strcmp(modeTok, "strain") == 0) {
                implexOpt.dtSource = LadrunoImplexOptions::DT_STRAIN;
            }
            else if (strcmp(modeTok, "user") == 0) {
                implexOpt.dtSource = LadrunoImplexOptions::DT_USER;
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &implexOpt.dtUser) != 0) {
                    opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                           << ": -implexDt user wants the value after it: -implexDt user $dt"
                           << endln;
                    return 0;
                }
                if (implexOpt.dtUser < 0.0) {
                    opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                           << ": -implexDt user $dt must be >= 0 (got " << implexOpt.dtUser
                           << ")." << endln;
                    return 0;
                }
            }
            else {
                opserr << "WARNING nDMaterial LadrunoSANISAND tag " << tag
                       << ": -implexDt wants pseudo|strain|user, got '" << modeTok << "'."
                       << " pseudo = ops_Dt (the default); strain = the norm of the strain"
                          " increment; user = a fixed value given after the keyword." << endln;
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
                       << " -Presidual / -Pmin / -honorTolR / -maxSubsteps /"
                       << " -implex / -implexControl / -implexAlpha / -implexDt" << endln;
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

    if (theMaterial == 0) {
        opserr << "WARNING ran out of memory for nDMaterial LadrunoSANISAND material with tag: "
               << tag << endln;
        return 0;
    }

    // Ladruno (ADR-92 P1): the IMPL-EX request, applied after construction (see
    // the LadrunoImplexOptions note in the header for why it is not a
    // constructor argument). setLadrunoImplexOptions() is where D3's scheme
    // refusals live; a refusal there means the deck asked for something ADR 92
    // does not qualify, and the material must NOT quietly run implicit under a
    // flag that says otherwise -- so the parser fails the command.
    if (sawImplexToken) {
        if (((LadrunoSANISAND *)theMaterial)->setLadrunoImplexOptions(implexOpt) != 0) {
            delete theMaterial;
            return 0;
        }
    }

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

    this->ladrunoImplexInitState();     // Ladruno (ADR-92 P1)
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

    this->ladrunoImplexInitState();     // Ladruno (ADR-92 P1)
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
    this->ladrunoImplexInitState();     // Ladruno (ADR-92 P1)
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
    this->ladrunoImplexInitState();     // Ladruno (ADR-92 P1)
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
    // Ladruno (ADR-92 P1): a revertToStart puts the material back at step 0, so
    // the extrapolation history has to go with it -- d_eps_p(n) from a discarded
    // load path is the single most misleading thing this class could carry
    // forward. The OPTION SET survives, exactly as p_residual does.
    this->ladrunoImplexInitState();
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
        // Ladruno (ADR-92 P1): the IMPL-EX request is not a constructor argument
        // (header note), so it is transferred here. Quietly -- this runs once per
        // Gauss point. The HISTORY is deliberately NOT transferred: a fresh
        // integration point starts with d_eps_p = 0, which makes its first
        // plastic step a pure elastic prediction, exactly as at the stage flip.
        clone->setLadrunoImplexOptions(mImplexOpt, false);                            // Ladruno
        return clone;
    } else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
        LadrunoSANISAND3D *clone;
        clone = new LadrunoSANISAND3D(this->getTag(), m_G0, m_nu, m_e_init, m_Mc, m_c, m_lambda_c,
                       m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz,
                       massDen, mScheme, mTangType, mJacoType, mTolF, mTolR,
                       mPresidualInput, mPminInput, mHonorTolR, mMaxSubsteps);        // Ladruno
        clone->setLadrunoImplexOptions(mImplexOpt, false);                            // Ladruno (ADR-92 P1)
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
//  Ladruno (ADR-92 P1) widened it again, 5 -> 22, for the IMPL-EX option set and
//  its ONE history variable:
//
//      data(5)      = (double)mImplexOpt.enabled
//      data(6)      = (double)mImplexOpt.control
//      data(7)      = mImplexOpt.errorTol
//      data(8)      = mImplexOpt.reductionLimit
//      data(9)      = mImplexOpt.alpha
//      data(10)     = (double)mImplexOpt.dtSource
//      data(11)     = mImplexOpt.dtUser
//      data(12)     = mImplexDt          )
//      data(13)     = mImplexDtCommit    ) the step state -- WITHOUT dt_n the
//      data(14)     = mImplexDt0         ) restored point would extrapolate by
//      data(15)     = mImplexFactor      ) the wrong factor on its first step
//      data(16..21) = mImplexDEpsP(0..5) = d_eps_p(n), the D1 history
//
//  d_eps_p is COMMITTED state and has to cross the wire for the same reason
//  mAlpha_n does: an MP worker that receives a material mid-analysis and
//  restarts its extrapolation from zero silently runs a different constitutive
//  law from the rank beside it -- which is the ADR-86 section 3 defect, in a new
//  variable. mImplexError / the clamp counters are diagnostics and are NOT sent.
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

    static Vector ladrunoData(22);                                                    // Ladruno (ADR-92 P1)

    ladrunoData(0) = mPresidualInput;
    ladrunoData(1) = mPminInput;
    ladrunoData(2) = m_Presidual;
    ladrunoData(3) = (double)mHonorTolR;
    ladrunoData(4) = (double)mMaxSubsteps;   // Ladruno (ADR-86b)

    // Ladruno (ADR-92 P1)
    ladrunoData(5)  = mImplexOpt.enabled ? 1.0 : 0.0;
    ladrunoData(6)  = mImplexOpt.control ? 1.0 : 0.0;
    ladrunoData(7)  = mImplexOpt.errorTol;
    ladrunoData(8)  = mImplexOpt.reductionLimit;
    ladrunoData(9)  = mImplexOpt.alpha;
    ladrunoData(10) = (double)mImplexOpt.dtSource;
    ladrunoData(11) = mImplexOpt.dtUser;
    ladrunoData(12) = mImplexDt;
    ladrunoData(13) = mImplexDtCommit;
    ladrunoData(14) = mImplexDt0;
    ladrunoData(15) = mImplexFactor;
    for (int i = 0; i < 6; i++)
        ladrunoData(16 + i) = mImplexDEpsP(i);

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

    static Vector ladrunoData(22);                                                    // Ladruno (ADR-92 P1)

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

    // Ladruno (ADR-92 P1). setLadrunoImplexOptions() is used rather than a raw
    // assignment so a restored material passes the SAME D3 scheme checks the
    // deck did -- if the sender's IntScheme did not survive the base recvSelf,
    // this is where that shows up, loudly, instead of running an unqualified
    // companion. Quiet, because recvSelf runs once per received material.
    {
        LadrunoImplexOptions opt;
        opt.enabled        = (ladrunoData(5) != 0.0);
        opt.control        = (ladrunoData(6) != 0.0);
        opt.errorTol       = ladrunoData(7);
        opt.reductionLimit = ladrunoData(8);
        opt.alpha          = ladrunoData(9);
        opt.dtSource       = (int)ladrunoData(10);
        opt.dtUser         = ladrunoData(11);
        this->ladrunoImplexInitState();
        if (opt.enabled && this->setLadrunoImplexOptions(opt, false) != 0) {
            opserr << "WARNING: LadrunoSANISAND::recvSelf - the received -implex option set"
                   << " was REFUSED by this build's ADR-92 D3 checks; this material will run"
                   << " the implicit companion only. The sending process was running IMPL-EX."
                   << endln;
        }
        mImplexDt       = ladrunoData(12);
        mImplexDtCommit = ladrunoData(13);
        mImplexDt0      = ladrunoData(14);
        mImplexFactor   = ladrunoData(15);
        for (int i = 0; i < 6; i++)
            mImplexDEpsP(i) = ladrunoData(16 + i);
    }

    // The base recvSelf restored m_Pmin from its own data(96) and never re-runs
    // initialize(); we take the last write here.
    this->applyLadrunoConstants();                                                    // Ladruno

    return 0;
}

// ===========================================================================
//  Ladruno (ADR-92 P1): IMPL-EX for LadrunoSANISAND
//
//  THE OPERATOR, in the material's own (compression-positive) convention, with
//  `n` the last committed step:
//
//      f        = (dt_{n+1} / dt_n) * implexAlpha
//      sigma~   = sigma_n + Ce(p_n) : ( (eps_{n+1} - eps_n) - f * d_eps_p(n) )
//      dsigma~/deps = Ce(p_n)                       (constant within the step)
//
//  and, at commitState ONLY, the true return from state `n` with the actual
//  strain increment, which produces sigma_implicit, the whole internal state,
//  and the next d_eps_p.
//
//  THREE THINGS THAT LOOK LIKE DETAILS AND ARE NOT:
//
//  1. INCREMENTAL, NOT TOTAL. SANISAND's elasticity is hypoelastic -- the base
//     integrates dsigma = Ce(p):deps_e with the moduli read at the COMMITTED
//     stress -- so sigma~ must be advanced FROM sigma_n and never rebuilt from a
//     total elastic strain. P0 measured the total form returning p = 1.11 kPa
//     against a committed 5.00 kPa (78 % error) at a ZERO strain increment. That
//     form is the negative control for this file, not an alternative to it.
//
//  2. Ce IS FROZEN AT THE COMMITTED MEAN STRESS. If the moduli are re-read on
//     the extrapolated stress the delivered tangent is no longer the operator the
//     stress was built with, and acceptance 1 (tangent identity to machine
//     precision; P0 measured 3.5e-11 on the oracle) fails. All three tangent
//     slots are written, so the identity holds whatever TanType the deck asked
//     for -- under -implex TanType is inert and the echo says so.
//
//  3. THE EXTRAPOLATED STRESS IS CLAMPED AT p_min. The floor bounds the
//     COMMITTED stress; it does not bound f*d_eps_p(n), and P0 measured sigma~
//     crossing into tension (min p = -1.37 kPa on a committed +0.0101 kPa state,
//     first order in the step but O(1)-O(10) relative). The code's own device --
//     `sigma~ = dev(sigma~) + p_min*I1` -- is applied here, and the DEVIATORIC
//     error is measured beside the pressure error at every commit, because ADR-92
//     P1 section 4 decides clamp-vs-bound-on-f on those two numbers.
//
//  WHY THE PLASTIC STRAIN AND NOT dGamma (D1, measured at P0: 18-22x better on
//  the deck-default integrator): mDGamma is the step total only under
//  BackwardEuler_CPPM. Under every substepped explicit scheme it is the LAST
//  SUBSTEP's multiplier, and at the corner the substep count swings by orders of
//  magnitude between steps, so extrapolating it extrapolates noise.
//  eps_p(n) = mEpsilon_n - mEpsilonE_n is committed and exact on the base, which
//  is why this whole file needs NO hook into vanilla ManzariDafalias (D5:
//  vanilla footprint ZERO).
// ===========================================================================

// Process-wide IMPL-EX error accounting, on the
// `ASDConcrete3DMaterial::GlobalParameters` template (:307-342). Anonymous
// namespace: this is one process's diagnostic accumulator, not an interface, and
// nothing outside this translation unit may reach it.
namespace {

class LadrunoImplexGlobals                                    // Ladruno (ADR-92 P1)
{
  public:
    static LadrunoImplexGlobals &instance(void)
    {
        static LadrunoImplexGlobals theInstance;
        return theInstance;
    }
    void accumulate(double e)
    {
        if (e > maxError)
            maxError = e;
        sumError += e;
        count++;
    }
    double getMaxError(void) const { return maxError; }
    // The average SINCE THE LAST QUERY, and it resets the accumulator -- the
    // same contract ASDConcrete3D's getAverageError() has, spelled out here
    // because "average" alone would be read as "since the analysis began".
    double takeAverageError(void)
    {
        double a = (count > 0) ? (sumError / (double)count) : 0.0;
        sumError = 0.0;
        count    = 0;
        return a;
    }

  private:
    LadrunoImplexGlobals() : maxError(0.0), sumError(0.0), count(0) {}
    LadrunoImplexGlobals(const LadrunoImplexGlobals &);
    LadrunoImplexGlobals &operator=(const LadrunoImplexGlobals &);

    double maxError;
    double sumError;
    long   count;
};

} // anonymous namespace

// Zeroes the IMPL-EX state and sizes the one history vector. Called from all
// four constructors, from initialize() (hence revertToStart), from recvSelf, and
// at the elastoplastic stage flip -- W5's "the history initialises at the stage
// switch, not before it". It does NOT touch mImplexOpt: the option set is the
// deck's request and survives a revertToStart exactly as p_residual does.
void
LadrunoSANISAND::ladrunoImplexInitState(void)
{
    if (mImplexDEpsP.Size() != 6)
        mImplexDEpsP.resize(6);
    mImplexDEpsP.Zero();

    mImplexDt         = 0.0;
    mImplexDtCommit   = 0.0;
    mImplexDt0        = 0.0;
    mImplexFactor     = 0.0;
    mImplexStepArmed  = true;
    mImplexTrialDone  = false;
    mImplexError      = 0.0;
    mImplexErrorDev   = 0.0;
    mImplexErrorVol   = 0.0;
    mImplexClampFired = false;
    mImplexClampCount = 0;
}

// ---------------------------------------------------------------------------
//  The option set, and D3's refusals.
//
//  Validated HERE and not in the parser because the checks need mScheme and
//  mMaxSubsteps, which do not exist until the material is built -- and because
//  getCopy(const char*) and recvSelf also reach this function, so a clone or a
//  restored material cannot smuggle past a rule the deck was refused. Same
//  argument as sanitiseLadrunoInputs(), one level up.
// ---------------------------------------------------------------------------
int
LadrunoSANISAND::setLadrunoImplexOptions(const LadrunoImplexOptions &opt, bool verbose)
{
    if (!opt.enabled) {
        // A control tolerance, an alpha or a dt source with no -implex to read
        // them is the "a flag claims to have done something it did not do"
        // defect this class exists to make impossible. Refuse, do not ignore.
        if (opt.control || opt.alpha != 1.0 ||
            opt.dtSource != LadrunoImplexOptions::DT_PSEUDO) {
            opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
                   << ": -implexControl / -implexAlpha / -implexDt were given without"
                      " -implex. Nothing would read them. Add -implex, or remove them."
                   << endln;
            return -1;
        }
        mImplexOpt = opt;
        this->ladrunoImplexInitState();
        return 0;
    }

    const int s = (int)mScheme;

    // D3, as REVERSED by P0. The companion must be an integrator whose failure
    // mode is bounded and observable.
    if (s != 1 && s != 2) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -implex REFUSES IntScheme " << s << ".";
        if (s == 3 || s == 5 || s == 7 || s == 8 || s == 9) {
            opserr << " ADR 92 D3: schemes 3/5/7/8/9 carry no error control on the"
                      " substep, so the IMPL-EX companion could not report a failed"
                      " return and -implexControl would have nothing to refuse.";
            if (s == 5)
                opserr << " Scheme 5 additionally runs ForwardEuler, whose `r` is"
                          " identically zero (a block-scoped redeclaration,"
                          " ManzariDafalias.cpp ForwardEuler): it silently drops both"
                          " volumetric coupling terms and cost 26 % in mobilised"
                          " strength when P0 measured it.";
        }
        else {
            opserr << " ADR 92 qualified only IntScheme 1 (ModifiedEuler, the deck"
                      " default) and IntScheme 2 (BackwardEuler_CPPM) as IMPL-EX"
                      " companions; nothing else was measured at P0, and an"
                      " unqualified companion is a wrong answer that passes every"
                      " gate.";
        }
        opserr << " Use IntScheme 1 with -maxSubsteps." << endln;
        return -1;
    }

    // The companion runs once per committed step and there is no global Newton
    // left to walk it through a seizure, so it MUST be able to fail. On scheme 1
    // that is precisely what -maxSubsteps buys (ADR-86b / #792 T1); without it
    // ModifiedEuler is bounded only by dT_min = 1e-6, i.e. up to 1e6 return maps
    // in one commit, and ADR-90 GATE U measured single updates of 34 minutes.
    if (s == 1 && mMaxSubsteps <= 0) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -implex on IntScheme 1 REQUIRES -maxSubsteps > 0 (got "
               << mMaxSubsteps << "). The IMPL-EX companion runs at commitState,"
                  " where no global Newton is left to react, so it must be able to"
                  " FAIL rather than force-accept at dT_min. See ADR 92 D3 and"
                  " ADR-90 GATE U." << endln;
        return -1;
    }

    if (s == 2 && verbose) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -implex with IntScheme 2 is PERMITTED but is not the ADR-92"
                  " default. P0 measured 58-74 % of BackwardEuler_CPPM's calls on the"
                  " low-confinement corner path taking its low-p branch, whose Newton"
                  " is disabled by a literal errFlag = 0 and which therefore integrates"
                  " by explicit_integrator (ModifiedEuler) -- so scheme 2 is not an"
                  " implicit return where this campaign's problem lives, and it costs a"
                  " 19-unknown Newton everywhere else." << endln;
        if (mMaxSubsteps <= 0)
            opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
                   << ": -implex with IntScheme 2 and -maxSubsteps 0. Scheme 2's own"
                      " retry ladder ends in explicit_integrator, so the companion can"
                      " still spend an unbounded number of ModifiedEuler substeps in one"
                      " commit. -maxSubsteps is strongly recommended." << endln;
    }

    mImplexOpt = opt;
    this->ladrunoImplexInitState();

    if (verbose) {
        opserr << "LadrunoSANISAND tag " << this->getTag()
               << ": IMPL-EX ON. companion IntScheme " << s
               << (s == 1 ? " (ModifiedEuler, capped at " : " (BackwardEuler_CPPM, cap ")
               << mMaxSubsteps << " substeps), alpha = " << mImplexOpt.alpha
               << ", dt = "
               << (mImplexOpt.dtSource == LadrunoImplexOptions::DT_PSEUDO ? "pseudo (ops_Dt)"
                  : (mImplexOpt.dtSource == LadrunoImplexOptions::DT_STRAIN ? "strain increment norm"
                                                                           : "user"));
        if (mImplexOpt.dtSource == LadrunoImplexOptions::DT_USER)
            opserr << " " << mImplexOpt.dtUser;
        if (mImplexOpt.control)
            opserr << ", -implexControl tol = " << mImplexOpt.errorTol
                   << " reductionLimit = " << mImplexOpt.reductionLimit
                   << " (refuses with LADRUNO_MATERIAL_REFUSED = " << LADRUNO_MATERIAL_REFUSED
                   << ", which only an element that PROPAGATES it can act on -- today"
                      " LadrunoBrick)";
        else
            opserr << ", -implexControl OFF (the extrapolation error is measured and"
                      " reported at every commit but NOTHING refuses a step; P0 measured"
                      " IMPL-EX unusable from d_eps = 5e-4 at p0 = 5 kPa, so at a"
                      " low-confinement corner the control is a requirement, not an option)";
        opserr << ". TanType " << (int)mTangType
               << " is INERT under -implex: the delivered tangent is Ce(p_n), which is"
                  " symmetric -- the non-associated consistent tangent disappears."
               << endln;
        opserr << "LadrunoSANISAND tag " << this->getTag()
               << ": READING HAZARD (ADR 92 section 8) -- every equilibrium on this"
                  " material is an equilibrium of the EXTRAPOLATED stress. Any limit"
                  " point read off an -implex leg must be confirmed on the implicit"
                  " material, and implexError must be printed beside the verdict."
               << endln;
    }

    return 0;
}

// ---------------------------------------------------------------------------
//  Stage-0 inertness (W5), in one expression.
//
//  mElastFlag is a STATIC shared by every ManzariDafalias instance and is
//  flipped for all of them at once by `updateMaterialStage`. While it is 0 the
//  base takes its elastic_integrator branch and -implex is inert, so gravity and
//  the `LoadControl 0.0` re-equilibration are bit-identical with the flag on and
//  off. There is no branch to get wrong: the extrapolated path is simply not
//  reachable.
// ---------------------------------------------------------------------------
bool
LadrunoSANISAND::ladrunoImplexActive(void) const
{
    return mImplexOpt.enabled && (mElastFlag != 0);
}

// ---------------------------------------------------------------------------
//  The one entry point both wrappers' setTrialStrain uses.
//
//  With -implex OFF this is `this->integrate(); return this->ladrunoUpdateStatus();`
//  in that order -- the two statements ADR-86b left in the wrappers -- so every
//  existing SANISAND deck is byte-identical.
// ---------------------------------------------------------------------------
int
LadrunoSANISAND::ladrunoTrialUpdate(void)
{
    if (!this->ladrunoImplexActive()) {
        mImplexTrialDone = false;   // so commitState() takes the base path
        this->integrate();
        return this->ladrunoUpdateStatus();
    }
    return this->ladrunoImplexTrial();
}

// ---------------------------------------------------------------------------
//  dt_{n+1} and the extrapolation factor, computed ONCE per step.
//
//  The freeze is not a micro-optimisation. Under `-implexDt strain` the source
//  is the norm of the strain increment, so an f recomputed on every trial call
//  would be a FUNCTION OF eps and d(sigma~)/d(eps) would no longer be Ce -- the
//  tangent identity (acceptance 1) would fail by exactly the term nobody would
//  think to look for. Freezing f at the first trial call after a commit makes it
//  a constant within the step for every dt source, which is what the identity
//  needs and what the P0 oracle does (one extrapolate() call per step).
//
//  D2's guards:
//    dt_{n+1} == 0   a hold. f = 0: no time advanced, so no plastic flow is
//                    predicted, and sigma~ = sigma_n + Ce:d_eps exactly.
//    dt_n == 0       the ratio is undefined (first step after a hold, or the
//                    first step of all). Fall back to f = alpha, which is what
//                    the P0 oracle's default does; on the first plastic step
//                    d_eps_p is zero anyway, so the fallback is unobservable there.
//    dt_{n+1} < 0    REFUSED in ladrunoImplexTrial(), not here.
// ---------------------------------------------------------------------------
void
LadrunoSANISAND::ladrunoImplexArmStep(void)
{
    double dt = 0.0;

    if (mImplexOpt.dtSource == LadrunoImplexOptions::DT_STRAIN) {
        Vector dEps(6);
        dEps = mEpsilon;
        dEps.addVector(1.0, mEpsilon_n, -1.0);
        dt = this->GetNorm_Cov(dEps);
    }
    else if (mImplexOpt.dtSource == LadrunoImplexOptions::DT_USER) {
        dt = mImplexOpt.dtUser;
    }
    else {
        // ops_Dt is the domain's `currentTime - committedTime`, i.e. the pseudo
        // time increment of the step. Under LoadControl on a prescribed-settlement
        // pattern it is proportional to the settlement increment, which is D2's
        // whole argument for the default.
        dt = ops_Dt;
    }

    mImplexDt = dt;
    if (mImplexDt0 <= 0.0 && dt > 0.0)
        mImplexDt0 = dt;

    if (mImplexDtCommit > 0.0)
        mImplexFactor = (dt / mImplexDtCommit) * mImplexOpt.alpha;
    else
        mImplexFactor = mImplexOpt.alpha;

    if (dt == 0.0)
        mImplexFactor = 0.0;

    mImplexStepArmed = false;
}

// Put every TRIAL member back on its committed twin. mEpsilon is deliberately
// NOT touched: it is the element's own input and the caller (the -implexControl
// probe) still needs it. revertToLastCommit() restores it separately.
void
LadrunoSANISAND::ladrunoRestoreTrialFromCommitted(void)
{
    mSigma     = mSigma_n;
    mEpsilonE  = mEpsilonE_n;
    mAlpha     = mAlpha_n;
    mFabric    = mFabric_n;
    mAlpha_in  = mAlpha_in_n;
    mDGamma    = mDGamma_n;
    // mVoidRatio has no committed twin on the base -- it is a single copy that
    // every integrator overwrites in place and commitState() re-derives. Re-derive
    // it here from the committed strain for the same reason, because
    // getVoidRatio() is public and a probe run must not leave it on a trial value.
    mVoidRatio = m_e_init - (1 + m_e_init) * this->GetTrace(mEpsilon_n);
}

// Ce(p_n) into all three tangent slots.
//
// It does NOT write mK / mG, and that omission is deliberate. Those two are the
// moduli integrate() hands its integrators, and vanilla's own contract for them
// is "whatever commitState() last derived from the committed stress" -- which
// this function's expression reproduces at every step EXCEPT the one right after
// an elastoplastic stage flip, where commitState() ran with mElastFlag == 0 (no
// sqrt(p/P_atm) factor) and this call would run with it == 1. Writing them would
// therefore silently change the companion's first plastic step relative to the
// implicit run it is supposed to be compared against. The extrapolation operator
// itself IS built with the stage-1 moduli, which is correct for sigma~ and keeps
// the tangent identity exact; the -implexControl probe saves and restores mK/mG
// around its integrate() so nothing it touches escapes either.
//
// mCe IS written, and that is not optional: integrate() reads it for the
// loading-reversal test (`trialDirection = mCe*(mEpsilon - mEpsilon_n)`), and the
// P0 oracle does that test with the COMMITTED stiffness, so writing it is what
// makes the companion reproduce the oracle rather than nearly reproduce it.
void
LadrunoSANISAND::ladrunoImplexFreezeTangent(double &K, double &G)
{
    this->GetElasticModuli(mSigma_n, mVoidRatio, K, G);
    mCe             = this->GetStiffness(K, G);
    mCep            = mCe;
    mCep_Consistent = mCe;
}

// The IMPL-EX error and its two parts, on ADR 92 section 2's definition:
//
//     implexError = ||sigma~ - sigma_impl|| / ( ||sigma_impl|| + P_atm*||eps|| )
//
//  with the contravariant norm on BOTH terms -- including on the strain, which
//  is what the P0 oracle does (`norm_contr(m.eps_n)`), so parity is exact rather
//  than nearly so.
//
//  The split: sigma = dev + p*I1 and dev : I1 = 0, so
//  ||d_sigma||^2 = ||d_dev||^2 + 3*(dp)^2 exactly, and the two reported numbers
//  are the two legs of that identity over the same denominator. ADR-92 P1
//  section 4 decides clamp-vs-bound on whether the deviatoric leg is the same
//  order the pressure leg was before the clamp -- which is why this is measured
//  at every commit and not only when someone asks.
void
LadrunoSANISAND::ladrunoImplexMeasureError(const Vector &sigTilde, const Vector &sigImplicit,
                                           const Vector &epsRef)
{
    Vector d(6);
    d = sigTilde;
    d.addVector(1.0, sigImplicit, -1.0);

    double den = this->GetNorm_Contr(sigImplicit) + m_P_atm * this->GetNorm_Contr(epsRef);
    if (den < small)
        den = small;

    Vector dDev(6);
    dDev = this->GetDevPart(d);

    mImplexError    = this->GetNorm_Contr(d) / den;
    mImplexErrorDev = this->GetNorm_Contr(dDev) / den;
    mImplexErrorVol = sqrt(3.0) * fabs(one3 * this->GetTrace(d)) / den;
}

// ---------------------------------------------------------------------------
//  The extrapolated update -- the response the element sees, every global
//  iteration.
//
//  WHERE -implexControl COMPUTES ITS IMPLICIT (ADR-92 P1 section 4, second open
//  question), and why this file computes it HERE rather than a-posteriori at
//  commit:
//
//    * only the in-step site can refuse the step it measured. The a-posteriori
//      site can only shrink the NEXT step, so the step whose error it measured
//      is already committed -- and P0 measured errors of 1.26 (p0 = 5 kPa,
//      d_eps = 5e-4) where the answer is not merely inaccurate but wrong
//      (q/p 0.09 against 2.07);
//    * only the in-step site can discover a COMPANION FAILURE before commit. If
//      the substep cap fires inside the commit-time return there is no good
//      move left: committing that state is garbage and refusing at commit
//      leaves the analysis believing a step advanced that did not;
//    * the cost objection -- ASD pays one implicit return per setTrialStrain --
//      is bounded by IMPL-EX's own structural claim: on a frozen operator the
//      global step is linear, so it is 1-2 returns per step against the 125
//      state-determination passes ADR-90 GATE U measured. If the BVP gate shows
//      the ladder still firing at CP1's rate, that claim is refuted and THIS is
//      the first thing to move.
//
//  The decision is recorded in the plan's section 4 and the BVP gate is its
//  falsifier. The a-posteriori number is produced for free either way: the
//  companion runs at every commit, so implexError is reported with -implexControl
//  off as well.
// ---------------------------------------------------------------------------
int
LadrunoSANISAND::ladrunoImplexTrial(void)
{
    // commitState() owes a companion return after this.
    mImplexTrialDone = true;

    // W6: freeze f for the step -- but ONLY from a trial call that actually
    // carries a strain increment.
    //
    // The naive `if (mImplexStepArmed) armStep();` is WRONG, and it is wrong on
    // the DEFAULT dt source. `Domain::revertToLastCommit()` does not stop at
    // reverting: it sets `dT = 0.0`, re-applies the committed load and ends
    // `return this->update();` (Domain.cpp), i.e. it pushes ONE zero-increment
    // state determination through every material at `ops_Dt == 0`. That call
    // would arm the step at `dt = 0`, hence `f = 0`, and CONSUME the arm -- so
    // the retried step that follows would run its whole ladder rung with the
    // plastic extrapolation switched off, silently, and `implexError` would be
    // the error of an operator the deck never asked for. A failed step followed
    // by a revert is an ordinary event in every adaptive run, and this is the
    // same trap `LEDGER_quirks` already records for rate-dependent materials
    // ("Domain::revertToLastCommit() sets dT = 0 and re-applies the load").
    //
    // A zero strain increment therefore takes f = 0 for THAT evaluation --
    // sigma~ = sigma_n exactly, which is what W6 requires of a `LoadControl 0.0`
    // hold -- and leaves the step ARMED for the first call that moves. On a
    // genuine hold every call is a zero increment, so f stays 0 throughout and
    // the hold is unchanged. Under `-implexDt pseudo` this cannot perturb any
    // step that does move: ops_Dt is the same on every call of a step, so which
    // call arms it is unobservable.
    if (mImplexStepArmed) {
        Vector dArm(6);
        dArm = mEpsilon;
        dArm.addVector(1.0, mEpsilon_n, -1.0);
        if (this->GetNorm_Cov(dArm) > 0.0) {
            this->ladrunoImplexArmStep();       // clears mImplexStepArmed
        } else {
            mImplexDt     = 0.0;
            mImplexFactor = 0.0;                // no strain advanced, no plastic flow predicted
        }
    }

    // D2's refusal, by symptom. A material has no handle on the integrator, so
    // "refuse DisplacementControl and arc length" is enforced where it actually
    // bites: a pseudo-time increment that has gone NEGATIVE is a load factor that
    // has turned round past a limit point, and dt_{n+1}/dt_n is then not the
    // extrapolation the operator assumes. NOTE THE GAP, it is real: a
    // DisplacementControl run BEFORE its limit point has d(lambda) > 0 and is NOT
    // caught here. Such a deck must pass -implexDt user (or strain).
    if (mImplexDt < 0.0) {
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": -implex REFUSES this step -- the pseudo-time increment is negative ("
               << mImplexDt << "). dt_{n+1}/dt_n is the extrapolation factor itself, so a"
                  " load factor that has turned round (DisplacementControl or arc length"
                  " past a limit point) makes it a wrong answer that would pass every"
                  " gate. Use -implexDt user or -implexDt strain on such a deck."
               << endln;
        this->ladrunoRestoreTrialFromCommitted();
        double Kr = 0.0, Gr = 0.0;
        this->ladrunoImplexFreezeTangent(Kr, Gr);
        return LADRUNO_MATERIAL_REFUSED;
    }

    // ---- -implexControl: the in-step implicit -----------------------------
    // A plain local, NOT a function-local static: a static would be shared by
    // every Gauss point in the process, and this one is live across a call to
    // integrate(). The allocation is noise beside the implicit return it holds.
    Vector sigImplicit(6);
    if (mImplexOpt.control) {
        // The probe must be invisible: mK / mG are vanilla's committed moduli and
        // integrate() overwrites them with the last substep's, which the
        // commit-time companion would then inherit.
        const double Kprobe = mK;
        const double Gprobe = mG;

        this->integrate();          // clobbers every trial member; restored below

        mK = Kprobe;
        mG = Gprobe;

        if (mSubstepCapHitInME) {
            // The companion could not integrate this increment. That is the one
            // thing -implexControl exists to catch before the step commits.
            this->ladrunoRestoreTrialFromCommitted();
            double Kf = 0.0, Gf = 0.0;
            this->ladrunoImplexFreezeTangent(Kf, Gf);
            return LADRUNO_MATERIAL_REFUSED;
        }

        sigImplicit = mSigma;
        this->ladrunoRestoreTrialFromCommitted();
    }

    // ---- W1: the extrapolated stress, INCREMENTAL --------------------------
    double K = 0.0, G = 0.0;
    this->ladrunoImplexFreezeTangent(K, G);

    Vector dEps(6);
    dEps = mEpsilon;
    dEps.addVector(1.0, mEpsilon_n,  -1.0);          // d_eps
    dEps.addVector(1.0, mImplexDEpsP, -mImplexFactor);  // d_eps - f*d_eps_p(n)

    mSigma = mSigma_n;
    mSigma.addVector(1.0, mCe * dEps, 1.0);

    // ---- W2: the p_min floor clamp ----------------------------------------
    // The same device the base applies to the committed stress (the pre-loop
    // guard in ModifiedEuler and its four siblings). The test omits m_Presidual
    // on both sides exactly as those sites do -- it cancels there and it cancels
    // here.
    //
    // AND A CONSEQUENCE THE ADR DID NOT CARRY, recorded because acceptance 1 has
    // to know about it: the clamp is a PROJECTION, so wherever it fires
    // d(sigma~)/d(eps) is the deviatoric projection of Ce, not Ce -- the tangent
    // identity holds exactly everywhere the clamp is idle and nowhere it is not.
    // The delivered operator stays Ce, which keeps it symmetric and positive
    // definite (the reason IMPL-EX was asked for), so this is a knowingly
    // inexact tangent at the floor, not a broken one. It is also a NEW argument
    // for ADR-92 P1 section 4's alternative: bounding f keeps sigma~ admissible
    // WITHOUT a projection and therefore keeps the identity exact everywhere.
    // See the section-4 record in _adr92_p1_execution_plan.md.
    mImplexClampFired = false;
    if (one3 * this->GetTrace(mSigma) < m_Pmin) {
        Vector sTilde(6);
        sTilde = this->GetDevPart(mSigma);
        mSigma = sTilde;
        mSigma.addVector(1.0, mI1, m_Pmin);
        mImplexClampFired = true;
        mImplexClampCount++;
    }

    // The trial elastic strain that goes with the extrapolated stress, so
    // getEStrain() / getPStrain() report the extrapolated split rather than a
    // stale committed one. Nothing in the base READS the trial mEpsilonE (the
    // integrators take mEpsilonE_n), so this is a recorder-facing consistency
    // fix and cannot feed back into the companion.
    mEpsilonE = mEpsilonE_n;
    mEpsilonE.addVector(1.0, dEps, 1.0);

    // mAlpha, mFabric, mAlpha_in and mDGamma are DELIBERATELY untouched: no
    // fabric is accumulated and no loading reversal is detected on an
    // extrapolated state (ADR 92 section 2). They advance only at commit.

    // ---- W7: the refusal ---------------------------------------------------
    if (mImplexOpt.control) {
        this->ladrunoImplexMeasureError(mSigma, sigImplicit, mEpsilon);
        if (mImplexError > mImplexOpt.errorTol) {
            // Below the reduction limit there is nothing left to cut, so refusing
            // again would only turn a bounded inaccuracy into a dead analysis.
            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
                return LADRUNO_MATERIAL_REFUSED;
        }
    }

    return 0;
}

// ---------------------------------------------------------------------------
//  The companion return, at commitState only.
// ---------------------------------------------------------------------------
int
LadrunoSANISAND::ladrunoImplexCommit(void)
{
    // The extrapolated stress the element actually equilibrated with.
    Vector sigTilde(6);
    sigTilde = mSigma;

    // eps_p(n), before the companion advances the committed state.
    Vector epsPOld(6);
    epsPOld = mEpsilon_n;
    epsPOld.addVector(1.0, mEpsilonE_n, -1.0);

    // The true return from state n with the actual strain increment.
    this->integrate();

    if (mSubstepCapHitInME) {
        // Reachable only with -implexControl OFF (with it on, the same failure was
        // caught in setTrialStrain and the step was refused). There is no good move
        // here and the message says so rather than picking one quietly: the state
        // ModifiedEuler left is NOT an answer, committing it would manufacture
        // history, and committing sigma~ instead would manufacture a different one.
        opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
               << ": the -implex COMPANION hit the -maxSubsteps cap at commitState."
                  " The committed state is NOT advanced and this commit is refused ("
               << LADRUNO_MATERIAL_REFUSED << "). This is the ADR-92 section 8 risk"
                  " measured: under -implex the global step is solved on the elastic"
                  " operator, so the increment handed to the commit-time return can be"
                  " larger and differently directed than implicit Newton would have"
                  " found. Re-run with -implexControl, which refuses such a step BEFORE"
                  " it converges." << endln;
        this->ladrunoRestoreTrialFromCommitted();
        mSigma = sigTilde;          // leave the element holding what it equilibrated with
        double Kf = 0.0, Gf = 0.0;
        this->ladrunoImplexFreezeTangent(Kf, Gf);
        return LADRUNO_MATERIAL_REFUSED;
    }

    Vector sigImplicit(6);
    sigImplicit = mSigma;

    // The committed state is the IMPLICIT one -- standard IMPL-EX, and the same
    // choice ASDConcrete3DMaterial::commitState() makes. Only the EQUILIBRIUM was
    // found on the extrapolated stress.
    int res = ManzariDafalias::commitState();

    // W3: the one new history variable, d_eps_p(n+1) = eps_p(n+1) - eps_p(n).
    // After the base commit, mEpsilon_n / mEpsilonE_n are the NEW committed pair.
    mImplexDEpsP = mEpsilon_n;
    mImplexDEpsP.addVector(1.0, mEpsilonE_n, -1.0);
    mImplexDEpsP.addVector(1.0, epsPOld,     -1.0);

    mImplexDtCommit  = mImplexDt;
    mImplexStepArmed = true;
    mImplexTrialDone = false;

    // The denominator uses the NEW committed strain, matching the P0 oracle's
    // `Implex.commit()` (which measures after `m.commit()`), so parity is exact.
    this->ladrunoImplexMeasureError(sigTilde, sigImplicit, mEpsilon_n);
    LadrunoImplexGlobals::instance().accumulate(mImplexError);

    // Hand the next step -- and anything that asks for a tangent between now and
    // the next setTrialStrain -- the frozen operator at the new committed p.
    double K = 0.0, G = 0.0;
    this->ladrunoImplexFreezeTangent(K, G);

    return res;
}

// ---------------------------------------------------------------------------
//  commitState / revertToLastCommit
// ---------------------------------------------------------------------------
int
LadrunoSANISAND::commitState(void)
{
    // mImplexTrialDone is false with -implex off AND on stage 0, so gravity and
    // the LoadControl 0.0 hold take the base path verbatim.
    if (!mImplexTrialDone)
        return ManzariDafalias::commitState();

    return this->ladrunoImplexCommit();
}

int
LadrunoSANISAND::revertToLastCommit(void)
{
    int res = ManzariDafalias::revertToLastCommit();   // an empty `return 0` today

    if (!mImplexOpt.enabled)
        return res;

    // Under -implex the trial members hold an EXTRAPOLATED state with no
    // companion behind it. A rejected step must put them back, and must re-arm
    // the step so the next attempt recomputes f from the new dt -- which is the
    // whole point of a subdivision.
    mEpsilon = mEpsilon_n;
    this->ladrunoRestoreTrialFromCommitted();

    mImplexTrialDone  = false;
    mImplexStepArmed  = true;
    mImplexDt         = mImplexDtCommit;
    mImplexClampFired = false;

    // The -implexControl probe's integrate() writes these two, and nothing on
    // the IMPL-EX path calls ladrunoUpdateStatus() to clear them. Left sticky
    // they would make the `substeps` response and Print() report a cap hit
    // belonging to a step that was thrown away. Diagnostic-only -- integrate()
    // resets both at its top -- but a diagnostic that survives a revert is a
    // diagnostic that lies.
    mSubstepsTakenInME = 0;
    mSubstepCapHitInME = false;

    double K = 0.0, G = 0.0;
    this->ladrunoImplexFreezeTangent(K, G);

    return res;
}

// ---------------------------------------------------------------------------
//  setParameter / updateParameter -- `implexError` and `avgImplexError` on the
//  ASDConcrete3DMaterial.cpp:2073-2077 template, plus the one writable knob.
//
//  IDs sit in the same 3308x band as the ADR-86b response id, for the same
//  reason: they only have to miss the base's 1..9, they are NOT class tags, and
//  nothing may derive one from them.
// ---------------------------------------------------------------------------
static const int LadrunoSanisandImplexErrorParamID    = 33087;   // Ladruno (ADR-92 P1)
static const int LadrunoSanisandAvgImplexErrorParamID = 33088;   // Ladruno (ADR-92 P1)
static const int LadrunoSanisandImplexDtParamID       = 33089;   // Ladruno (ADR-92 P1)

int
LadrunoSANISAND::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc > 0) {
        // Deliberately NOT guarded by argv[1] == this->getTag(), matching the ASD
        // template: these two are PROCESS-wide accumulators, not per-material
        // state, and a deck asks for them once.
        if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
            param.setValue(LadrunoImplexGlobals::instance().getMaxError());
            return param.addObject(LadrunoSanisandImplexErrorParamID, this);
        }
        if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
            param.setValue(LadrunoImplexGlobals::instance().takeAverageError());
            return param.addObject(LadrunoSanisandAvgImplexErrorParamID, this);
        }
        if (strcmp(argv[0], "implexDt") == 0 || strcmp(argv[0], "implexDT") == 0) {
            param.setValue(mImplexOpt.dtUser);
            return param.addObject(LadrunoSanisandImplexDtParamID, this);
        }
    }
    return ManzariDafalias::setParameter(argv, argc, param);
}

int
LadrunoSANISAND::updateParameter(int parameterID, Information &info)
{
    // Read-only diagnostics: a `parameter` object bound to them exists so a deck
    // can RECORD them; writing one would falsify the accumulator it reports.
    if (parameterID == LadrunoSanisandImplexErrorParamID ||
        parameterID == LadrunoSanisandAvgImplexErrorParamID)
        return 0;

    if (parameterID == LadrunoSanisandImplexDtParamID) {
        if (info.theDouble < 0.0) {
            opserr << "WARNING LadrunoSANISAND tag " << this->getTag()
                   << ": implexDt must be >= 0 (got " << info.theDouble << "); ignored."
                   << endln;
            return 0;
        }
        mImplexOpt.dtUser  = info.theDouble;
        mImplexOpt.dtSource = LadrunoImplexOptions::DT_USER;
        return 0;
    }

    int res = ManzariDafalias::updateParameter(parameterID, info);

    // W5: the history initialises AT THE STAGE FLIP, not before it. ids 1 and 5
    // are the base's `updateMaterialStage` and `materialState`; both write the
    // static mElastFlag and, on 1, run Elastic2Plastic(). Anything IMPL-EX
    // accumulated during the elastic stage is the wrong history to extrapolate
    // from, and there is no reason a stage-0 dt should set the ratio for the
    // first plastic step either.
    if (res == 0 && mImplexOpt.enabled && (parameterID == 1 || parameterID == 5)) {
        if (mElastFlag != 0)
            this->ladrunoImplexInitState();
    }

    return res;
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

// Ladruno (ADR-92 P1): the IMPL-EX diagnostics. Same band, same rule -- these
// are response ids, not class tags, and nothing may derive one from them.
//
// `implexError` is the PER-INTEGRATION-POINT value at the last commit, which is
// what a recorder wants; `avgImplexError` is the process accumulator on the
// ASDConcrete3D template. `implexDetail` carries the split ADR-92 P1 section 4
// needs to decide clamp-vs-bound, plus the clamp's own counters -- gate 4's
// evidence, readable from a deck instead of a debugger.
constexpr int LadrunoSanisandImplexErrorResponseID    = 33090;   // Ladruno (ADR-92 P1)
constexpr int LadrunoSanisandAvgImplexErrorResponseID = 33091;   // Ladruno (ADR-92 P1)
constexpr int LadrunoSanisandImplexDetailResponseID   = 33092;   // Ladruno (ADR-92 P1)

Response *
LadrunoSANISAND::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc > 0 && (strcmp(argv[0], "substeps") == 0 ||
                     strcmp(argv[0], "substepsME") == 0 ||
                     strcmp(argv[0], "ladrunoSubsteps") == 0)) {
        static Vector probe(2);
        return new MaterialResponse(this, LadrunoSanisandSubstepResponseID, probe);
    }
    // Ladruno (ADR-92 P1)
    if (argc > 0 && (strcmp(argv[0], "implexError") == 0 ||
                     strcmp(argv[0], "ImplexError") == 0)) {
        static Vector probe1(1);
        return new MaterialResponse(this, LadrunoSanisandImplexErrorResponseID, probe1);
    }
    if (argc > 0 && (strcmp(argv[0], "avgImplexError") == 0 ||
                     strcmp(argv[0], "AvgImplexError") == 0)) {
        static Vector probe1(1);
        return new MaterialResponse(this, LadrunoSanisandAvgImplexErrorResponseID, probe1);
    }
    if (argc > 0 && (strcmp(argv[0], "implexDetail") == 0 ||
                     strcmp(argv[0], "ImplexDetail") == 0)) {
        static Vector probe6(6);
        return new MaterialResponse(this, LadrunoSanisandImplexDetailResponseID, probe6);
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
    // Ladruno (ADR-92 P1)
    if (responseID == LadrunoSanisandImplexErrorResponseID) {
        static Vector out1(1);
        out1(0) = mImplexError;
        return matInformation.setVector(out1);
    }
    if (responseID == LadrunoSanisandAvgImplexErrorResponseID) {
        static Vector out1(1);
        out1(0) = LadrunoImplexGlobals::instance().takeAverageError();
        return matInformation.setVector(out1);
    }
    if (responseID == LadrunoSanisandImplexDetailResponseID) {
        static Vector out6(6);
        out6(0) = mImplexError;                        // total
        out6(1) = mImplexErrorDev;                     // deviatoric leg
        out6(2) = mImplexErrorVol;                     // volumetric leg (sqrt(3)|dp|)
        out6(3) = mImplexClampFired ? 1.0 : 0.0;       // p_min clamp on the LAST pass
        out6(4) = (double)mImplexClampCount;           // ... and how often, ever
        out6(5) = mImplexFactor;                       // f, frozen for this step
        return matInformation.setVector(out6);
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

    // Ladruno (ADR-92 P1): the IMPL-EX record. A curve produced on this material
    // cannot be read without these numbers -- see ADR 92 section 8's reading
    // hazard and section 10's disclosure text.
    if (!mImplexOpt.enabled) {
        s << "  IMPL-EX     = OFF (ADR 92; -implex not given, so this material is"
             " byte-identical to the pre-ADR-92 build)" << endln;
    } else {
        s << "  IMPL-EX     = ON (ADR 92 P1).  alpha = " << mImplexOpt.alpha
          << ",  dt source = "
          << (mImplexOpt.dtSource == LadrunoImplexOptions::DT_PSEUDO ? "pseudo (ops_Dt)"
             : (mImplexOpt.dtSource == LadrunoImplexOptions::DT_STRAIN ? "strain-increment norm"
                                                                       : "user"))
          << ",  f (this step) = " << mImplexFactor << endln;
        s << "              dt_n+1 = " << mImplexDt << ", dt_n = " << mImplexDtCommit
          << ", dt_0 = " << mImplexDt0
          << (mElastFlag == 0 ? "   [stage 0: INERT, the elastic integrator is running]" : "")
          << endln;
        s << "              -implexControl = " << (mImplexOpt.control ? "ON" : "OFF");
        if (mImplexOpt.control)
            s << " (tol " << mImplexOpt.errorTol << ", reductionLimit "
              << mImplexOpt.reductionLimit << "; refuses with "
              << LADRUNO_MATERIAL_REFUSED << ")";
        s << endln;
        s << "              implexError (last commit) = " << mImplexError
          << "   [deviatoric " << mImplexErrorDev
          << ", volumetric " << mImplexErrorVol << "]" << endln;
        s << "              p_min clamp on sigma~: "
          << (mImplexClampFired ? "FIRED on the last pass" : "not on the last pass")
          << ", " << mImplexClampCount << " time(s) at this point" << endln;
        s << "              READING HAZARD: every equilibrium here is an equilibrium of"
             " the EXTRAPOLATED" << endln;
        s << "              stress, not of the converged one. A limit point read off this"
             " leg must be" << endln;
        s << "              confirmed on the implicit material, and implexError printed"
             " beside it." << endln;
        s << "              The step size is the regularization parameter and it has no"
             " length in it:" << endln;
        s << "              every band and post-peak branch read off an -implex leg is"
             " regularized by dt." << endln;
    }

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
