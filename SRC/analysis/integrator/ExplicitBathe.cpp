/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024
// Ladruno ADR-52 W1-E2 (nmb, 06/2026): collapsed the six ExplicitBathe* class tags
// (ExplicitBathe / *LNVD / *SMS / *SMSConsistent / *LNVDSMS / *LNVDSMSConsistent)
// into this one class with orthogonal -lnvd / -sms / -consistent flags. The six
// historical command names + class tags remain as deprecated aliases for one release.

#include <ExplicitBathe.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>         // Ladruno (W1-E2, LNVD): base formUnbalance()
#include <Vector.h>
#include <Matrix.h>                        // Ladruno (W1-E2, consistent): M_bar publish
#include <ID.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>                          // Ladruno (ADR-30 P5): tie-force scatter
#include <NodeIter.h>
#include <LadrunoConstraintProjector.h>    // Ladruno (ADR-30 P5): constraint projection
#include <LadrunoMassScaling.h>            // Ladruno (W1-E2, -sms): build/apply scaling
#include <LadrunoConsistentRefine.h>       // Ladruno (W1-E2, -consistent): PCG refine
#include <LadrunoMassScalingEnergy.h>      // Ladruno (W1-E2, -consistent): ALLVD blocks
#include <DiagonalSOE.h>                   // Ladruno (W1-E2, -consistent): SOE-type check

#include <classTags.h>

#define OPS_Export

// LAPACK eigenvalue solver declarations (kept for ABI parity; the eigensolve itself
// now lives in the shared CriticalTimeStep.{h,cpp}).
#ifdef _WIN32
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);
#else
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);
#endif

// =====================================================================================
//  Ladruno (W1-E2): flag <-> legacy-classTag bijection.
//  The three capability axes (lnvd, sms, consistent) are orthogonal, but "consistent"
//  is only meaningful with "sms". The six valid combos map 1:1 onto the six historical
//  integrator tags, so a serialized object reports the exact legacy tag and any saved
//  DB / parallel recvSelf reconstructs through the existing broker switch.
// =====================================================================================
int ExplicitBathe::tagForFlags(bool lnvd, bool sms, bool consistent)
{
    if (!sms)        consistent = false;   // consistent requires sms (defensive)
    if (!lnvd && !sms)                return INTEGRATOR_TAGS_ExplicitBathe;              // 33000
    if ( lnvd && !sms)                return INTEGRATOR_TAGS_ExplicitBatheLNVD;          // 33002
    if (!lnvd &&  sms && !consistent) return INTEGRATOR_TAGS_ExplicitBatheSMS;           // 33009
    if (!lnvd &&  sms &&  consistent) return INTEGRATOR_TAGS_ExplicitBatheSMSConsistent; // 33010
    if ( lnvd &&  sms && !consistent) return INTEGRATOR_TAGS_ExplicitBatheLNVDSMS;       // 33011
    return INTEGRATOR_TAGS_ExplicitBatheLNVDSMSConsistent;                               // 33012
}

void ExplicitBathe::flagsForTag(int classTag, bool &lnvd, bool &sms, bool &consistent)
{
    switch (classTag) {
    case INTEGRATOR_TAGS_ExplicitBathe:
        lnvd = false; sms = false; consistent = false; break;
    case INTEGRATOR_TAGS_ExplicitBatheLNVD:
        lnvd = true;  sms = false; consistent = false; break;
    case INTEGRATOR_TAGS_ExplicitBatheSMS:
        lnvd = false; sms = true;  consistent = false; break;
    case INTEGRATOR_TAGS_ExplicitBatheSMSConsistent:
        lnvd = false; sms = true;  consistent = true;  break;
    case INTEGRATOR_TAGS_ExplicitBatheLNVDSMS:
        lnvd = true;  sms = true;  consistent = false; break;
    case INTEGRATOR_TAGS_ExplicitBatheLNVDSMSConsistent:
        lnvd = true;  sms = true;  consistent = true;  break;
    default:
        lnvd = false; sms = false; consistent = false; break;
    }
}

// Ladruno (W1-E2): broker factory — build a default object with the flag set decoded
// from a (possibly retired) integrator classTag. recvSelf() then fills the parameters.
ExplicitBathe *ExplicitBathe::makeForBroker(int classTag)
{
    bool lnvd = false, sms = false, consistent = false;
    flagsForTag(classTag, lnvd, sms, consistent);
    // dtTarget/maxAddedMass/pcg defaults match the historical broker defaults; the real
    // values arrive via recvSelf. compute_critical_timestep starts at 0 (recvSelf is the
    // sole post-broker initializer, matching the old per-class default ctors).
    return new ExplicitBathe(0.0, /*cct*/0, /*verbose*/false, /*cflAbort*/false,
                             /*div*/0.0, /*cflUseTangent*/false, /*cflRecompute*/0,
                             CTSLumping::RowSum, lnvd, /*alpha*/0.0,
                             sms, consistent, /*dtTarget*/0.0, /*maxAddedMass*/0.05,
                             /*pcgTol*/1.0e-10, /*pcgMaxIt*/200);
}

// =====================================================================================
//  OPS factories. The plain factory is the primary command; the other five are
//  DEPRECATED aliases retained one release (each forces a fixed flag set). Their
//  historical positional grammars + validation messages are preserved verbatim so
//  existing fork scripts keep parsing.
// =====================================================================================

// Usage: integrator ExplicitBathe <p> <computeCriticalTimestep>
//                                  <-verbose> <-cfl> <-cflAbort> <-divergence f>
//                                  <-tangent> <-recompute N> <-lump rowsum|diagonal|hrz>
//                                  <-lnvd <alpha>> <-sms <dtTarget>> <-consistent>
//                                  <-maxAddedMass f> <-pcgTol t> <-pcgMaxIt n>
void *OPS_ExplicitBathe(void) {
    double p = 0.54;            // default: good high-frequency dissipation
    int compute_critical_timestep = 0;
    bool verbose = false;
    bool cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    CTSLumping lumping = CTSLumping::RowSum;
    // Ladruno (W1-E2): orthogonal capability flags on the unified command.
    bool useLNVD = false;       double alpha_flac = 0.8;     // FLAC default if -lnvd given
    bool useSMS = false, useConsistent = false;
    double dtTarget = 0.0, maxAddedMassFrac = 0.05;
    double pcgTol = 1.0e-10;    int pcgMaxIt = 200;

    // p is the leading numeric positional. Read it with the typed getter so it
    // works under BOTH Tcl and OpenSeesPy (OPS_GetString mis-reads a numeric arg).
    if (OPS_GetNumRemainingInputArgs() > 0) {
        int nd = 1;
        if (OPS_GetDoubleInput(&nd, &p) < 0) {
            opserr << "WARNING ExplicitBathe - could not read p (give p first, "
                      "e.g. integrator ExplicitBathe 0.54)\n";
            return 0;
        }
    }
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verbose = true;
        } else if (strcmp(arg, "-cflAbort") == 0) {
            cflAbort = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-cfl") == 0 ||
                   strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-divergence") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetDoubleInput(&nd, &divergenceFactor);
            }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-recompute") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetIntInput(&nd, &cflRecomputeEvery);
            }
            if (cflRecomputeEvery <= 0)
                opserr << "WARNING ExplicitBathe -recompute expects a positive integer N "
                          "(steps between dt_cr refreshes); dt_cr will be computed once\n";
            cflUseTangent = true;            // recomputing K0 every N steps is pointless
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING ExplicitBathe - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping rowsum)\n";
            }
        } else if (strcmp(arg, "-lnvd") == 0) {
            // Ladruno (W1-E2): FLAC local non-viscous damping; optional alpha follows
            // (default stays 0.8). Use the PROVEN contact `-soft` optional-value idiom:
            // peek as a STRING (OPS_GetString advances on BOTH Tcl + Py), classify by the
            // leading '-' (value vs next flag) — NOT by strtod, because under openseespy
            // OPS_GetString returns "Invalid String Input!" for a numeric arg, so a strtod
            // test would mis-classify it as non-numeric and silently DROP the alpha (the
            // documented "getString mis-reads a numeric Py arg" quirk). Then un-read and, if
            // it was a value, read it Python-safely with OPS_GetDoubleInput.
            useLNVD = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *peek = OPS_GetString();
                bool isFlag = (peek != 0 && peek[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!isFlag) {
                    int nd = 1; double a;
                    if (OPS_GetDoubleInput(&nd, &a) < 0) {
                        opserr << "WARNING ExplicitBathe -lnvd - need an alpha value or omit "
                                  "it (default 0.8)\n";
                        return 0;
                    }
                    alpha_flac = a;
                }
            }
        } else if (strcmp(arg, "-sms") == 0) {
            // Ladruno (W1-E2): selective mass scaling; REQUIRES a dtTarget.
            useSMS = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                if (OPS_GetDoubleInput(&nd, &dtTarget) < 0)
                    opserr << "WARNING ExplicitBathe -sms expects a positive dtTarget\n";
            }
        } else if (strcmp(arg, "-consistent") == 0) {
            useConsistent = true;
        } else if (strcmp(arg, "-maxAddedMass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int nd = 1; OPS_GetDoubleInput(&nd, &maxAddedMassFrac); }
        } else if (strcmp(arg, "-pcgTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int nd = 1; OPS_GetDoubleInput(&nd, &pcgTol); }
        } else if (strcmp(arg, "-pcgMaxIt") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int nd = 1; OPS_GetIntInput(&nd, &pcgMaxIt); }
        } else {
            opserr << "WARNING ExplicitBathe - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBathe - p must be in (0,1); typical 0.50-0.95, "
                  "recommended 0.54 (got " << p << ")\n";
        return 0;
    }
    if (useLNVD && (alpha_flac < 0.0 || alpha_flac >= 1.0)) {
        opserr << "WARNING ExplicitBathe -lnvd alpha must be in [0,1); FLAC default 0.8 "
                  "(got " << alpha_flac << ")\n";
        return 0;
    }
    if (useConsistent && !useSMS) {
        opserr << "WARNING ExplicitBathe -consistent requires -sms <dtTarget>; ignoring "
                  "-consistent\n";
        useConsistent = false;
    }
    if (useSMS && dtTarget <= 0.0) {
        opserr << "WARNING ExplicitBathe -sms expects a positive dtTarget\n";
        return 0;
    }
    if (useSMS) {
        compute_critical_timestep = 1;   // report the pre-scaling dt_cr (as the SMS aliases do)
        // W1-E3a (MF-1): under mass scaling the per-element eigensolve cannot see the nodal
        // augmentation, so -cflAbort/-recompute on the un-augmented pencil would (wrongly)
        // abort a run that -sms made stable at dt<=dtTarget. DOWNGRADE to report-only —
        // mirrors OPS_ExplicitBatheSMS_impl and the class docstring — instead of aborting.
        if (cflAbort || cflRecomputeEvery > 0) {
            opserr << "NOTE ExplicitBathe - -cflAbort/-recompute are downgraded to REPORT-ONLY "
                      "under -sms mass scaling (the un-augmented element pencil cannot see the "
                      "scaling mass, MF-1); the pre-scaling dt_cr is still reported.\n";
            cflAbort = false;
            cflRecomputeEvery = 0;
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBathe(p, compute_critical_timestep, verbose, cflAbort,
                          divergenceFactor, cflUseTangent, cflRecomputeEvery, lumping,
                          useLNVD, alpha_flac, useSMS, useConsistent,
                          dtTarget, maxAddedMassFrac, pcgTol, pcgMaxIt);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    return theIntegrator;
}

// DEPRECATED alias (was ExplicitBatheLNVD). Usage: integrator ExplicitBatheLNVD <p>
//   <alpha> <computeCriticalTimestep> <-verbose> <-cfl> <-cflAbort> <-divergence f>
void *OPS_ExplicitBatheLNVD(void) {
    double p = 0.54;
    double alpha_flac = 0.8;          // classic FLAC default
    int compute_critical_timestep = 0;
    bool verbose = false, cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    CTSLumping lumping = CTSLumping::RowSum;

    int navail = OPS_GetNumRemainingInputArgs();
    if (navail >= 2) {
        int nd = 2; double vals[2];
        if (OPS_GetDoubleInput(&nd, vals) < 0) {
            opserr << "WARNING ExplicitBatheLNVD - could not read p and alpha "
                      "(e.g. integrator ExplicitBatheLNVD 0.54 0.8)\n";
            return 0;
        }
        p = vals[0]; alpha_flac = vals[1];
    } else if (navail == 1) {
        int nd = 1;
        if (OPS_GetDoubleInput(&nd, &p) < 0) {
            opserr << "WARNING ExplicitBatheLNVD - could not read p\n";
            return 0;
        }
    }
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verbose = true;
        } else if (strcmp(arg, "-cflAbort") == 0) {
            cflAbort = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-cfl") == 0 ||
                   strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-divergence") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetDoubleInput(&nd, &divergenceFactor);
            }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-recompute") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetIntInput(&nd, &cflRecomputeEvery);
            }
            if (cflRecomputeEvery <= 0)
                opserr << "WARNING ExplicitBatheLNVD -recompute expects a positive integer N "
                          "(steps between dt_cr refreshes); dt_cr will be computed once\n";
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING ExplicitBatheLNVD - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping rowsum)\n";
            }
        } else {
            opserr << "WARNING ExplicitBatheLNVD - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVD - p must be in (0,1), recommended 0.54 "
                  "(got " << p << ")\n";
        return 0;
    }
    if (alpha_flac < 0.0 || alpha_flac >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVD - alpha must be in [0,1); FLAC default "
                  "0.8 (got " << alpha_flac << ")\n";
        return 0;
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBathe(p, compute_critical_timestep, verbose, cflAbort,
                          divergenceFactor, cflUseTangent, cflRecomputeEvery, lumping,
                          /*useLNVD*/true, alpha_flac);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheLNVD integrator\n";
    return theIntegrator;
}

// Shared parser body for the four DEPRECATED selective-mass-scaling aliases. `isLNVD`
// adds FLAC damping (an alpha positional follows dtTarget); `isConsistent` selects the
// Olovsson matrix-free PCG path (extra -pcgTol/-pcgMaxIt options).
static void *OPS_ExplicitBatheSMS_impl(const char *name, bool isLNVD, bool isConsistent)
{
    int nNeed = isLNVD ? 3 : 2;   // (p, alpha, dtTarget) or (p, dtTarget)
    if (OPS_GetNumRemainingInputArgs() < nNeed) {
        opserr << "WARNING integrator " << name
               << (isLNVD ? " $p $alpha $dtTarget <options>" : " $p $dtTarget <options>")
               << " - needs the Noh-Bathe p" << (isLNVD ? ", FLAC alpha," : "")
               << " and a target time step\n";
        return 0;
    }

    double pin[3] = {0.0, 0.0, 0.0};
    int nLead = nNeed;
    if (OPS_GetDoubleInput(&nLead, pin) < 0) {
        opserr << "WARNING " << name << " - could not read the leading positionals\n";
        return 0;
    }
    double p          = pin[0];
    double alpha_flac = isLNVD ? pin[1] : 0.0;
    double dtTarget   = isLNVD ? pin[2] : pin[1];
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING " << name << " - $p must be in (0,1) (Noh-Bathe sub-step)\n";
        return 0;
    }
    if (isLNVD && (alpha_flac < 0.0 || alpha_flac >= 1.0)) {
        opserr << "WARNING " << name << " - $alpha must be in [0,1) (FLAC local damping)\n";
        return 0;
    }
    if (dtTarget <= 0.0) {
        opserr << "WARNING " << name << " - $dtTarget must be a positive double\n";
        return 0;
    }

    double maxAddedMassFrac = 0.05;
    bool   verbose = false;
    bool   cflUseTangent = false;
    CTSLumping lumping = CTSLumping::Diagonal;   // matches the system Diagonal run
    double pcgTol = 1.0e-10;
    int    pcgMaxIt = 200;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verbose = true;
        } else if (strcmp(arg, "-maxAddedMass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetDoubleInput(&n1, &maxAddedMassFrac); }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
        } else if (strcmp(arg, "-cfl") == 0 || strcmp(arg, "-criticalTimestep") == 0) {
            // dt_cr is always reported under SMS; accept the flag for compatibility.
        } else if (isConsistent && strcmp(arg, "-pcgTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetDoubleInput(&n1, &pcgTol); }
        } else if (isConsistent && strcmp(arg, "-pcgMaxIt") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) { int n1 = 1; OPS_GetIntInput(&n1, &pcgMaxIt); }
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING " << name << " - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else if (strcmp(arg, "-cflAbort") == 0 || strcmp(arg, "-recompute") == 0) {
            // W1-E3a (ADR-52): do NOT refuse the run. Under SMS these flags cannot do
            // their job -- the element-mass eigensolve can't see the nodal augmentation,
            // so an abort/recompute on the un-augmented pencil would be wrong (MF-1).
            // DOWNGRADE to report-only: keep the integrator and report the pre-scaling dt_cr.
            opserr << "NOTE " << name << " - " << arg << " is downgraded to REPORT-ONLY "
                      "under mass scaling (no abort/recompute on the un-augmented element "
                      "pencil, MF-1); the pre-scaling dt_cr will be reported each "
                      "domainChanged.\n";
            verbose = true;
        } else {
            opserr << "WARNING " << name << " - unknown option " << arg << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBathe(p, /*cct*/1, verbose, /*cflAbort*/false, /*div*/0.0,
                          cflUseTangent, /*cflRecompute*/0, lumping,
                          isLNVD, alpha_flac, /*useSMS*/true, isConsistent,
                          dtTarget, maxAddedMassFrac, pcgTol, pcgMaxIt);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating " << name << " integrator\n";
    return theIntegrator;
}

void *OPS_ExplicitBatheSMS(void) {
    return OPS_ExplicitBatheSMS_impl("ExplicitBatheSMS", /*isLNVD*/false, /*isConsistent*/false);
}
void *OPS_ExplicitBatheSMSConsistent(void) {
    return OPS_ExplicitBatheSMS_impl("ExplicitBatheSMSConsistent", false, true);
}
void *OPS_ExplicitBatheLNVDSMS(void) {
    return OPS_ExplicitBatheSMS_impl("ExplicitBatheLNVDSMS", true, false);
}
void *OPS_ExplicitBatheLNVDSMSConsistent(void) {
    return OPS_ExplicitBatheSMS_impl("ExplicitBatheLNVDSMSConsistent", true, true);
}

// =====================================================================================
//  Constructors
// =====================================================================================

// Default constructor (broker for classTag 33000)
ExplicitBathe::ExplicitBathe()
    : ExplicitBathe(INTEGRATOR_TAGS_ExplicitBathe, 0.0, 0, false, false, 0.0,
                    false, 0, CTSLumping::RowSum,
                    false, 0.0, false, false, 0.0, 0.05, 1.0e-10, 200)
{}

// Public unified constructor — classTag is derived from the flag combo.
ExplicitBathe::ExplicitBathe(double _p, int compute_critical_timestep_,
                             bool verbose_, bool cflAbort_, double divergenceFactor_,
                             bool cflUseTangent_, int cflRecomputeEvery_,
                             CTSLumping lumping_,
                             bool useLNVD_, double alpha_flac_,
                             bool useSMS_, bool useConsistent_,
                             double dtTarget_, double maxAddedMassFrac_,
                             double pcgTol_, int pcgMaxIt_)
    : ExplicitBathe(tagForFlags(useLNVD_, useSMS_, useConsistent_),
                    _p, compute_critical_timestep_, verbose_, cflAbort_, divergenceFactor_,
                    cflUseTangent_, cflRecomputeEvery_, lumping_,
                    useLNVD_, alpha_flac_, useSMS_, useConsistent_,
                    dtTarget_, maxAddedMassFrac_, pcgTol_, pcgMaxIt_)
{}

// Master constructor (the real init). q0/q1/q2 from p:
//   q1 = (1 - 2p) / (2p(1-p)); q2 = 0.5 - p*q1; q0 = -q1 - q2 + 0.5
ExplicitBathe::ExplicitBathe(int classTag, double _p, int compute_critical_timestep_,
                             bool verbose_, bool cflAbort_, double divergenceFactor_,
                             bool cflUseTangent_, int cflRecomputeEvery_,
                             CTSLumping lumping_,
                             bool useLNVD_, double alpha_flac_,
                             bool useSMS_, bool useConsistent_,
                             double dtTarget_, double maxAddedMassFrac_,
                             double pcgTol_, int pcgMaxIt_)
    : TransientIntegrator(classTag),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), R_tdt(0),
      updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(verbose_), cflAbort(cflAbort_), divergenceFactor(divergenceFactor_),
      prevKE(0.0), firstStep(true),
      cflUseTangent(cflUseTangent_), cflRecomputeEvery(cflRecomputeEvery_),
      cflStepCount(0), cflFirstComputation(true), lumping(lumping_),
      // Ladruno (W1-E2): LNVD
      useLNVD(useLNVD_), alpha_flac(alpha_flac_), lastUnbalanceNorm(-1.0),
      // Ladruno (W1-E2): SMS (consistent requires sms)
      useSMS(useSMS_), useConsistent(useSMS_ && useConsistent_),
      dtTarget(dtTarget_), maxAddedMassFrac(maxAddedMassFrac_),
      lumpingSMS(lumping_), useTangentSMS(cflUseTangent_),
      pcgTol(pcgTol_), pcgMaxIt(pcgMaxIt_),
      injected(), scaled(false), appliedDomain(0), blocks(0),
      warnedLimitations(false), warnedSolver(false), warnedPCG(false),
      theProjector(0), massBuilt(false)   // Ladruno (ADR-30 P5)
{
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    if (verbose) {
        opserr << "ExplicitBathe: p = " << p
               << ", compute_critical_timestep = " << compute_critical_timestep;
        if (useLNVD) opserr << ", -lnvd alpha = " << alpha_flac;
        if (useSMS)  opserr << ", -sms dtTarget = " << dtTarget
                            << (useConsistent ? " (consistent)" : " (lumped)");
        opserr << endln;
    }
}

// Destructor - clean up allocated memory + any applied mass scaling
ExplicitBathe::~ExplicitBathe() {
    // Ladruno (W1-E2, -sms lumped): undo the additive nodal injection on teardown.
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    // Ladruno (W1-E2, -consistent): release the energy-recorder blocks (owner-guarded).
    if (useConsistent)
        Ladruno::MassScalingEnergyRegistry::instance().clear(this);
    if (blocks != 0) delete blocks;

    if (U_t) delete U_t;
    if (V_t) delete V_t;
    if (A_t) delete A_t;
    if (U_tpdt) delete U_tpdt;
    if (V_tpdt) delete V_tpdt;
    if (V_fake) delete V_fake;
    if (A_tpdt) delete A_tpdt;
    if (U_tdt) delete U_tdt;
    if (V_tdt) delete V_tdt;
    if (A_tdt) delete A_tdt;
    if (R_tdt) delete R_tdt;
}

// Ladruno (ADR-30 P5): the projection handler pushes its (non-owning) projector.
// Re-read diag(M) next starter (massBuilt=false) since a new projector implies a
// new constraint set / numbering. theProjector is owned by LadrunoProjectionHandler.
void ExplicitBathe::setConstraintProjector(LadrunoConstraintProjector *pr)
{
    theProjector = pr;
    massBuilt = false;
}

// =====================================================================================
//  Ladruno (W1-E2): selective mass scaling helpers (inert unless useSMS)
// =====================================================================================
void ExplicitBathe::removeScaling(void)
{
    if (scaled && appliedDomain != 0 && !injected.empty())
        Ladruno::applyMassScaling(appliedDomain, injected, -1.0);
    injected.clear();
    scaled = false;
    appliedDomain = 0;
}

int ExplicitBathe::applyMassScalingSMS(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
    if (theModel == 0 || theDomain == 0) {
        opserr << "ExplicitBathe(-sms)::domainChanged - missing model/domain\n";
        return -1;
    }

    if (!useConsistent) {
        // ---- lumped (DT2MS-style) additive nodal injection ----
        if (!warnedLimitations) {
            warnedLimitations = true;
            opserr << "ExplicitBathe(-sms): v1 limitations -- (1) elements touching an "
                      "MP-constraint SLAVE node are EXCLUDED from scaling (they remain "
                      "governing). (2) sizing accounts for betaK Rayleigh damping "
                      "(closed-form s=T^2+2*T*betaK/dt_e); alphaM is not folded in. (3) in "
                      "PARALLEL the injected lumped mass on a shared node IS summed across "
                      "ranks by a distributed/MPI diagonal solver (`system MPIDiagonal` / "
                      "OpenSeesSP `system Diagonal`).\n";
        }
        Ladruno::MassScalingReport rep =
            Ladruno::buildMassScaling(theModel, dtTarget, lumpingSMS, useTangentSMS, injected);
        Ladruno::applyMassScaling(theDomain, injected, +1.0);
        scaled = true;
        appliedDomain = theDomain;

        double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
        bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
        if (verbose || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
            opserr << "ExplicitBathe(-sms): dtTarget=" << dtTarget
                   << " scaled " << rep.nScaled << "/" << rep.nElems
                   << " elements; added mass " << (100.0 * frac) << "% of model mass"
                   << "; PRE-SCALING dt_cr estimate=" << rep.minDtScaled
                   << " (governing un-scaled element step; AFTER scaling the run is stable "
                      "at dt <= dtTarget=" << dtTarget << ")\n";
        }
        if (rep.nSelfReport > 0)
            opserr << "WARNING ExplicitBathe(-sms): " << rep.nSelfReport
                   << " throttling element(s) have a MASS-INDEPENDENT (self-reported) bound "
                      "and were NOT scaled -- lower dtTarget to satisfy them.\n";
        if (rep.nMismatch > 0)
            opserr << "WARNING ExplicitBathe(-sms): " << rep.nMismatch
                   << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
        if (rep.nConstrained > 0)
            opserr << "WARNING ExplicitBathe(-sms): " << rep.nConstrained
                   << " sub-target element(s) touch an MP-constrained (slave) node and were "
                      "EXCLUDED -- they still GOVERN at dt_e=" << rep.minDtConstrained
                   << " < dtTarget=" << dtTarget << ".\n";
        if (over)
            opserr << "WARNING ExplicitBathe(-sms): added mass " << (100.0 * frac)
                   << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac)
                   << "% (proceeding; the scaled inertia shifts global frequencies)\n";
        return 0;
    }

    // ---- consistent (Olovsson) matrix-free PCG ----
    if (blocks == 0) blocks = new std::vector<Ladruno::ConsistentBlock>();
    if (!warnedLimitations) {
        warnedLimitations = true;
        opserr << "ExplicitBathe(-sms -consistent): consistent (Olovsson) scaling on "
                  "Noh-Bathe -- (1) elements touching an MP-constraint slave OR master node "
                  "are EXCLUDED. (2) sizing accounts for betaK damping; alphaM is not folded "
                  "in. (3) both Noh-Bathe sub-steps solve M_tilde a = r by matrix-free PCG -- "
                  "REQUIRES `system Diagonal` (serial) or `system MPIDiagonal` (OpenSeesMP). "
                  "(4) in PARALLEL the distributed PCG reduces shared/boundary-node coupling "
                  "across ranks -- use `system MPIDiagonal`.\n";
    }

    Ladruno::MassScalingReport rep =
        Ladruno::buildMassScalingConsistent(theModel, dtTarget, lumpingSMS,
                                            useTangentSMS, *blocks);

    double frac = (rep.modelMass > 0.0) ? rep.addedMass / rep.modelMass : 0.0;
    bool over = (maxAddedMassFrac > 0.0 && frac > maxAddedMassFrac);
    if (verbose || over || rep.nSelfReport > 0 || rep.nMismatch > 0 || rep.nConstrained > 0) {
        opserr << "ExplicitBathe(-sms -consistent): dtTarget=" << dtTarget
               << " scaled " << rep.nScaled << "/" << rep.nElems
               << " elements (Olovsson); added mass " << (100.0 * frac) << "% of model mass"
               << "; PRE-SCALING dt_cr estimate=" << rep.minDtScaled
               << " (governing un-scaled element step; AFTER scaling the run is stable "
                  "at dt <= dtTarget=" << dtTarget << ")\n";
    }
    if (rep.nSelfReport > 0)
        opserr << "WARNING ExplicitBathe(-sms -consistent): " << rep.nSelfReport
               << " throttling element(s) have a MASS-INDEPENDENT bound and were NOT scaled.\n";
    if (rep.nMismatch > 0)
        opserr << "WARNING ExplicitBathe(-sms -consistent): " << rep.nMismatch
               << " element(s) had a non-node-major / DOF-mismatch mass and were SKIPPED.\n";
    if (rep.nConstrained > 0)
        opserr << "WARNING ExplicitBathe(-sms -consistent): " << rep.nConstrained
               << " sub-target element(s) touch an MP-constrained node and were EXCLUDED -- "
                  "they still GOVERN at dt_e=" << rep.minDtConstrained << ".\n";
    if (over)
        opserr << "WARNING ExplicitBathe(-sms -consistent): added mass " << (100.0 * frac)
               << "% exceeds -maxAddedMass cap " << (100.0 * maxAddedMassFrac) << "%\n";

    {
        LinearSOE *soe = this->getLinearSOE();
        bool okSOE = (dynamic_cast<DiagonalSOE *>(soe) != 0) ||
                     (soe != 0 && soe->isDistributedDiagonal());
        if (rep.nScaled > 0 && !okSOE)
            opserr << "WARNING ExplicitBathe(-sms -consistent): " << rep.nScaled
                   << " element(s) need scaling but the SOE is NOT `system Diagonal`/"
                      "`system MPIDiagonal` -- the consistent mass CANNOT be applied and dt="
                   << dtTarget << " will run UNSCALED -> expect INSTABILITY.\n";
    }

    // Ladruno V4: publish the node-major M_bar blocks (keyed by element tag) for the
    // EnergyBalanceRecorder. Empty when nothing is scaled -> registry stays inactive.
    {
        std::map<int, Matrix> ebBlocks;
        for (size_t i = 0; i < blocks->size(); ++i)
            ebBlocks[(*blocks)[i].eleTag] = (*blocks)[i].Mbar;
        Ladruno::MassScalingEnergyRegistry::instance().publish(this, ebBlocks);
    }
    return 0;
}

bool ExplicitBathe::refinesAccel(void) const
{
    return useSMS && useConsistent;
}

int ExplicitBathe::refineAccel(Vector &a)
{
    if (!useSMS || !useConsistent || blocks == 0)
        return 0;
    return Ladruno::applyConsistentRefine(*blocks, this->getLinearSOE(), a,
                                          pcgTol, pcgMaxIt, verbose,
                                          warnedSolver, warnedPCG,
                                          "ExplicitBathe(-sms -consistent)");
}

// =====================================================================================
//  Ladruno (W1-E2, -lnvd): unified formUnbalance() + FLAC local non-viscous damping
// =====================================================================================
int ExplicitBathe::formUnbalance(void) {
    int res = this->IncrementalIntegrator::formUnbalance();
    if (res < 0) return res;

    if (useLNVD) {
        LinearSOE *theSOE = this->getLinearSOE();
        if (theSOE != 0)
            lastUnbalanceNorm = theSOE->getB().pNorm(0);
        if (alpha_flac > 0.0)
            this->addLocalDamping();
    }
    return res;
}

void ExplicitBathe::addLocalDamping(void) {
    LinearSOE *theSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theSOE == 0 || theModel == 0) return;

    const Vector &B = theSOE->getB();
    int n = B.Size();
    if (n == 0) return;

    // gather the global trial velocity at the instant being solved
    Vector vglob(n);
    vglob.Zero();
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        const Vector &v = dofPtr->getTrialVel();
        for (int i = 0; i < id.Size(); ++i) {
            int loc = id(i);
            if (loc >= 0 && loc < n)
                vglob(loc) = v(i);
        }
    }

    // r_i <- r_i - alpha*|r_i|*sign(v_i)  (opposes motion, vanishes at v->0)
    Vector damped(B);
    for (int i = 0; i < n; ++i) {
        double v = vglob(i);
        double sgn = (v > 0.0) ? 1.0 : ((v < 0.0) ? -1.0 : 0.0);
        damped(i) = B(i) - alpha_flac * std::abs(B(i)) * sgn;
    }
    theSOE->setB(damped);
}

// Handle domain changes (mesh modifications, initial conditions, etc.)
int ExplicitBathe::domainChanged() {
    // Ladruno (W1-E2, -sms lumped): re-baseline onto the ORIGINAL masses before re-sizing.
    if (useSMS && !useConsistent)
        removeScaling();
    // Ladruno (W1-E2, -sms consistent): drop the stale blocks before rebuilding.
    if (useSMS && useConsistent && blocks != 0)
        blocks->clear();

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBathe::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    if (size == 0) {
        opserr << "ExplicitBathe::domainChanged - invalid size\n";
        return -1;
    }

    // Allocate or reallocate state vectors if size has changed
    if (!U_t || U_t->Size() != size) {
        // Delete old vectors
        if (U_t) delete U_t;
        if (V_t) delete V_t;
        if (A_t) delete A_t;
        if (U_tpdt) delete U_tpdt;
        if (V_tpdt) delete V_tpdt;
        if (V_fake) delete V_fake;
        if (A_tpdt) delete A_tpdt;
        if (U_tdt) delete U_tdt;
        if (V_tdt) delete V_tdt;
        if (A_tdt) delete A_tdt;

        // Create new vectors
        U_t = new Vector(size);
        V_t = new Vector(size);
        A_t = new Vector(size);
        U_tpdt = new Vector(size);
        V_tpdt = new Vector(size);
        V_fake = new Vector(size);
        A_tpdt = new Vector(size);
        U_tdt = new Vector(size);
        V_tdt = new Vector(size);
        A_tdt = new Vector(size);

        // Verify allocation was successful
        if (!U_t || !V_t || !A_t || !U_tpdt || !V_tpdt || !V_fake ||
            !A_tpdt || !U_tdt || !V_tdt || !A_tdt) {
            opserr << "ExplicitBathe::domainChanged - out of memory\n";
            return -1;
        }
    }

    // Initialize state vectors from committed DOF values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;

    while ((dofPtr = theDOFs()) != nullptr) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        // Initialize displacements
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*U_t)(loc) = disp(i);
                (*U_tpdt)(loc) = disp(i);
                (*U_tdt)(loc) = disp(i);
            }
        }

        // Initialize velocities
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*V_t)(loc) = vel(i);
            }
        }

        // Initialize accelerations
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*A_t)(loc) = accel(i);
            }
        }
    }

    // Ladruno (ADR-30 P5): the projected-mass cache is invalid after a domain change;
    // re-read diag(M) at the next starter. Then check (or snap) IC compliance of the
    // committed (u0, v0) against every constraint group — U_t/V_t hold them now.
    // enforceIC() needs no mass (it uses L, delta only), so it is safe here.
    massBuilt = false;
    if (theProjector != 0) {
        int nViol = theProjector->enforceIC(*U_t, *V_t);
        if (nViol > 0) {
            opserr << "ExplicitBathe::domainChanged - " << nViol
                   << " initial-condition constraint violation(s); aborting. Fix the "
                      "ICs or add -projectICs to the constraints command.\n";
            return -1;
        }
    }

    // Initialize critical timestep values
    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    // Reset computation flag if it was already computed; also reset the
    // periodic-recompute counter so the cadence restarts cleanly after a
    // mesh change / staged construction.
    if (compute_critical_timestep == 2) {
        compute_critical_timestep = 1;
    }
    cflStepCount = 0;

    // Ladruno (W1-E2, -sms): (re)build + apply the selective mass scaling for this mesh.
    if (useSMS) {
        if (applyMassScalingSMS() < 0)
            return -1;
    }

    return 0;
}

// The per-element critical-time-step eigensolve lives in the shared
// SRC/analysis/integrator/CriticalTimeStep.{h,cpp} (pulled in via ExplicitBathe.h).
// See computeCriticalTimeStep().

// Advance to a new time step (first Noh-Bathe sub-step: t -> t + p*dt)
int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;

    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
        return -1;
    }

    // Each step performs exactly two solves; reset here so a failed/retried
    // step cannot poison the next step's update() counter (matches CentralDifference).
    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();

    // Cold-start diagnostic: an explicit scheme started from rest with zero
    // committed initial acceleration produces no motion on the first step if a
    // load is already applied. Warn once so the user runs an initial
    // equilibrium step or supplies a consistent a0.
    if (firstStep) {
        firstStep = false;
        if (A_t->pNorm(0) == 0.0 && V_t->pNorm(0) == 0.0)
            opserr << "ExplicitBathe - NOTE starting from rest with zero initial "
                      "acceleration; if a load is already applied, run one "
                      "equilibrium step (or set the initial acceleration) so a0 is "
                      "consistent.\n";
    }

    // Ladruno (ADR-30 P5): build the projector's diag(M) once per stage (re-armed by
    // massBuilt=false at domainChanged / setConstraintProjector) and project the
    // committed initial acceleration onto the constraint manifold. We must read M from
    // the assembled Diagonal SOE BEFORE the algorithm's solve() factors it in place
    // (DiagonalDirectSolver overwrites A[i] with 1/A[i]); so we form the (mass) tangent
    // here ourselves. The algorithm re-forms M before its own solve, so the extra
    // assembly is harmless (cost: one mass assembly on the first step of a stage).
    if (theProjector != 0 && !massBuilt) {
        LinearSOE *theLinSOE = this->getLinearSOE();
        if (this->formTangent(CURRENT_TANGENT) < 0) {
            opserr << "ExplicitBathe::newStep() - starter: formTangent failed (projector)\n";
            return -3;
        }
        if (theProjector->buildMass(theLinSOE) < 0) {
            opserr << "ExplicitBathe::newStep() - starter: projector buildMass failed "
                      "(massless DOF / non-Diagonal SOE)\n";
            return -3;
        }
        massBuilt = true;
        theProjector->project(*A_t);   // project committed a0 onto the manifold
    }

    // Critical-time-step estimate (per-element generalized eigenproblem).
    // damped_minimum_critical_timestep is the CONSERVATIVE central-difference
    // limit; the Noh-Bathe scheme is stable up to ~EB_NB_STABILITY_FACTOR times it.
    // With -recompute N it is refreshed every N steps from the tangent stiffness
    // so it tracks softening/stiffening in nonlinear runs.
    if (cflRecomputeEvery > 0 && compute_critical_timestep == 2) {
        if (++cflStepCount >= cflRecomputeEvery) {
            compute_critical_timestep = 1;
            cflStepCount = 0;
        }
    }
    if (compute_critical_timestep == 1) {
        // shared eigensolve (DSYGV/DGGEV, relative-beta, chosen lumping, MPI_MIN
        // reduced across ranks). Fresh CTSResult, so no manual reset needed.
        CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
        damped_minimum_critical_timestep   = r.damped_dt;
        undamped_minimum_critical_timestep = r.undamped_dt;
        damped_critical_element_tag        = r.damped_tag;
        undamped_critical_element_tag      = r.undamped_tag;
        compute_critical_timestep = 2;  // mark as computed

        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (cflFirstComputation || verbose) {
            opserr << "ExplicitBathe: critical time step estimate"
                   << (cflUseTangent ? " (tangent)" : "") << "\n"
                   << "  central-difference limit (conservative): "
                   << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n"
                   << "  Noh-Bathe limit (~" << EB_NB_STABILITY_FACTOR << "x): " << dt_nb << "\n";
            if (deltaT > dt_nb)
                opserr << "  WARNING dt = " << deltaT
                       << " exceeds the Noh-Bathe stability limit - expect INSTABILITY.\n";
            else if (deltaT > damped_minimum_critical_timestep)
                opserr << "  note: dt exceeds the conservative limit but is within the "
                          "Noh-Bathe stable range.\n";
        }
        cflFirstComputation = false;
    }

    // Optional hard CFL guard.
    if (cflAbort && damped_minimum_critical_timestep > 0.0) {
        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (deltaT > dt_nb) {
            opserr << "ExplicitBathe::newStep() - ABORT: dt = " << deltaT
                   << " > Noh-Bathe stability limit " << dt_nb << "\n";
            return -2;
        }
    }

    // Optional per-step reporting (off by default).
    if (verbose && damped_minimum_critical_timestep > 0.0) {
        const double f = deltaT / damped_minimum_critical_timestep;
        opserr << "ExplicitBathe::newStep() dt = " << deltaT
               << " dt_cd = " << damped_minimum_critical_timestep
               << " factor = " << f << (f < 1.0 ? " [OK]" : " [> CD limit]") << endln;
    }

    // Compute integration coefficients for this time step
    a0 = p * deltaT;
    a1 = std::pow(p * deltaT, 2) / 2.0;
    a2 = a0 / 2.0;
    a3 = (1.0 - p) * deltaT;
    a4 = std::pow((1.0 - p) * deltaT, 2) / 2.0;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // Predict displacement and velocity at t + p*dt
    *U_tpdt = *U_t;
    U_tpdt->addVector(1.0, *V_t, a0);      // += v_t * p*dt
    U_tpdt->addVector(1.0, *A_t, a1);      // += a_t * (p*dt)^2/2

    *V_fake = *V_t;
    V_fake->addVector(1.0, *A_t, a0);      // += a_t * p*dt

    A_tpdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + p * deltaT;

    if (theModel->updateDomain(newtime, p * deltaT) < 0) {
        opserr << "ExplicitBathe::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

// Update the response quantities. Called twice per time step:
//   1st call: after solving at t + p*dt, prepare for t + dt
//   2nd call: after solving at t + dt, finalize the time step
int ExplicitBathe::update(const Vector &U) {
    updateCount++;
    if (updateCount > 2) {
        opserr << "WARNING ExplicitBathe::update() - called more than twice in a step.\n";
        opserr << "  ExplicitBathe is an explicit scheme and requires 'algorithm Linear' "
                  "(exactly 2 solves per step).\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitBathe::update() - no AnalysisModel set\n";
        return -2;
    }

    if (U_t == 0) {
        opserr << "WARNING ExplicitBathe::update() - domainChange() failed or not called\n";
        return -3;
    }

    if (U.Size() != A_t->Size()) {
        opserr << "WARNING ExplicitBathe::update() - Vector size mismatch\n";
        return -4;
    }

    LinearSOE *theLinSOE = this->getLinearSOE();

    // Store acceleration at t + p*dt
    *A_tpdt = U;
    // Ladruno (W1-E2, -consistent): consistent mass scaling refines this sub-step-1
    // accel (a = M_tilde^-1 r) from the just-factored diagonal SOE. No-op otherwise.
    if (this->refinesAccel())
        this->refineAccel(*A_tpdt);
    // Ladruno (ADR-30 P5): project the sub-step-1 accel onto the constraint manifold
    // (after any consistent-mass refine) before it feeds V_tpdt / U_tdt.
    if (theProjector != 0)
        theProjector->project(*A_tpdt);

    // Update velocity at t + p*dt (corrected)
    // v_{t+p*dt} = v_t + (a_t + a_{t+p*dt}) * p*dt/2
    *V_tpdt = *V_t;
    V_tpdt->addVector(1.0, *A_t, a2);      // += a_t * p*dt/2
    V_tpdt->addVector(1.0, *A_tpdt, a2);   // += a_{t+p*dt} * p*dt/2

    // Prepare fake velocity for response setting
    *V_fake = *V_tpdt;
    V_fake->addVector(1.0, *A_tpdt, a3);   // += a_{t+p*dt} * (1-p)*dt

    // Predict displacement at t + dt
    // u_{t+dt} = u_{t+p*dt} + v_{t+p*dt} * (1-p)*dt + a_{t+p*dt} * ((1-p)*dt)^2/2
    *U_tdt = *U_tpdt;
    U_tdt->addVector(1.0, *V_tpdt, a3);    // += v_{t+p*dt} * (1-p)*dt
    U_tdt->addVector(1.0, *A_tpdt, a4);    // += a_{t+p*dt} * ((1-p)*dt)^2/2

    A_tdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + (1.0 - p) * deltaT;

    if (theModel->updateDomain(newtime, (1.0 - p) * deltaT) < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -3;
    }

    // Solve for acceleration at t + dt (formUnbalance() adds FLAC damping when -lnvd)
    this->formUnbalance();
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();
    // Ladruno (W1-E2, -consistent): refine this sub-step-2 accel. The SOE diagonal is
    // still the factored 1/mass (DiagonalDirectSolver factors once; "just solve" reuse),
    // so refineAccel's r-recovery is valid here too.
    if (this->refinesAccel())
        this->refineAccel(*A_tdt);
    // Ladruno (ADR-30 P5): project the sub-step-2 accel (the load-bearing one — it
    // becomes A_t at commit and drives next step's predictor → manifold invariance).
    // The NaN/Inf circuit breaker below deliberately runs on the PROJECTED accel.
    if (theProjector != 0)
        theProjector->project(*A_tdt);

    // Circuit breaker 1: catch a blown-up (NaN/Inf) acceleration before it
    // silently propagates through the rest of the run.
    const double A_max = A_tdt->pNorm(0);
    if (A_max != A_max || A_max == std::numeric_limits<double>::infinity()) {
        opserr << "ExplicitBathe::update() - ABORT: non-finite acceleration "
                  "(NaN/Inf) - the integration has diverged (dt likely too large).\n";
        return -5;
    }

    // Circuit breaker 2 (opt-in): abort on runaway kinetic energy growth.
    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*V_tpdt) ^ (*V_tpdt));   // velocity-based KE proxy
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "ExplicitBathe::update() - ABORT: kinetic-energy proxy grew by "
                   << (ke / prevKE) << "x in one step (> " << divergenceFactor
                   << ") - spurious energy growth / instability.\n";
            return -6;
        }
        if (ke > 0.0) prevKE = ke;
    }

    if (verbose) {
        opserr << "ExplicitBathe::update() max|a| = " << A_max;
        if (useLNVD) opserr << "  |unbalance| = " << lastUnbalanceNorm;
        opserr << endln;
    }

    // Final velocity update at t + dt using all three accelerations
    // v_{t+dt} = v_{t+p*dt} + q0*a_t*(1-p)*dt + (0.5+q1)*a_{t+p*dt}*(1-p)*dt + q2*a_{t+dt}*(1-p)*dt
    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);       // += q0 * a_t * (1-p)*dt
    V_tdt->addVector(1.0, *A_tpdt, a6);    // += (0.5+q1) * a_{t+p*dt} * (1-p)*dt
    V_tdt->addVector(1.0, *A_tdt, a7);     // += q2 * a_{t+dt} * (1-p)*dt

    // Set final response
    theModel->setResponse(*U_tdt, *V_tdt, *A_tdt);

    if (theModel->updateDomain() < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

// Form element tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

// Form nodal tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// Commit the state for this time step
int ExplicitBathe::commit() {
    updateCount = 0;  // Reset update counter for next step

    // Update state vectors: t+dt becomes new t
    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBathe::commit() - no AnalysisModel set\n";
        return -1;
    }

    // Ladruno (ADR-30 P4/P5): scatter the projector's per-equation tie force
    // f = M(a_raw - a_proj) (cached during the step's last project() = sub-step 2)
    // onto the nodes, so the node-based LadrunoRecorder can record a constraintTieForce
    // field under the Noh-Bathe integrators too. Done before commitDomain() (recorders
    // run after commit). Untied / SP-fixed DOFs -> 0.
    if (theProjector != 0 && theProjector->isMassReady()) {
        Domain *theDomain = theModel->getDomainPtr();
        if (theDomain != 0) {
            NodeIter &theNodes = theDomain->getNodes();
            Node *nodePtr;
            while ((nodePtr = theNodes()) != 0) {
                DOF_Group *dg = nodePtr->getDOF_GroupPtr();
                if (dg == 0) continue;
                const ID &id = dg->getID();
                int ndf = nodePtr->getNumberDOF();
                Vector tf(ndf);
                for (int d = 0; d < ndf && d < id.Size(); d++) {
                    int eqn = id(d);
                    tf(d) = (eqn >= 0) ? theProjector->tieForceAtEqn(eqn) : 0.0;
                }
                nodePtr->setProjectionTieForce(tf);
            }
        }
    }

    return theModel->commitDomain();
}

// Get current velocity (for modal damping interface)
const Vector &ExplicitBathe::getVel() {
    return *V_t;
}

// Conservative (central-difference) critical time step; valid once a step has
// run with critical-time-step computation enabled. <=0 if not yet computed.
double ExplicitBathe::getCriticalTimeStep(void) const {
    if (compute_critical_timestep == 0)
        return -1.0;
    if (damped_minimum_critical_timestep == std::numeric_limits<double>::infinity())
        return -1.0;
    return damped_minimum_critical_timestep;
}

// Ladruno (W1-E2, was ExplicitBatheLNVD): |r|_inf at the most recent solve (<0 until
// the first solve, or if -lnvd is off — only the LNVD path tracks it).
double ExplicitBathe::getUnbalanceNorm(void) const {
    return lastUnbalanceNorm;
}

// =====================================================================================
//  Serialization. The full parameter superset is sent every time (fixed size); the
//  {lnvd, sms, consistent} flags travel implicitly through the integrator classTag,
//  which the broker uses (makeForBroker) to construct the right object before recvSelf.
// =====================================================================================
int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(11);
    data(0)  = p;
    data(1)  = (double)(int)lumping;          // 0=RowSum, 1=Diagonal, 2=HRZ
    data(2)  = alpha_flac;
    data(3)  = dtTarget;
    data(4)  = maxAddedMassFrac;
    data(5)  = (double)(int)lumpingSMS;
    data(6)  = useTangentSMS ? 1.0 : 0.0;
    data(7)  = pcgTol;
    data(8)  = (double)pcgMaxIt;
    data(9)  = verbose ? 1.0 : 0.0;
    data(10) = (double)compute_critical_timestep;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(11);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);
    {   // decode lumping; default to RowSum (this integrator's default)
        int lc = (int)data(1);
        lumping = (lc == 2) ? CTSLumping::HRZ
                : (lc == 1) ? CTSLumping::Diagonal
                            : CTSLumping::RowSum;
    }
    alpha_flac       = data(2);
    dtTarget         = data(3);
    maxAddedMassFrac = data(4);
    {
        int lc = (int)data(5);
        lumpingSMS = (lc == 2) ? CTSLumping::HRZ
                   : (lc == 1) ? CTSLumping::Diagonal
                               : CTSLumping::RowSum;
    }
    useTangentSMS = (data(6) != 0.0);
    pcgTol        = data(7);
    pcgMaxIt      = (int)data(8);
    verbose       = (data(9) != 0.0);
    compute_critical_timestep = (int)data(10);

    // The {useLNVD, useSMS, useConsistent} flags are set at construction from the
    // classTag (makeForBroker); a fresh-receive resets the transient SMS bookkeeping.
    injected.clear();
    scaled = false;
    appliedDomain = 0;
    if (blocks != 0) blocks->clear();
    if (useConsistent)
        Ladruno::MassScalingEnergyRegistry::instance().clear(this);
    warnedLimitations = false;
    warnedSolver = false;
    warnedPCG = false;
    lastUnbalanceNorm = -1.0;

    // Recalculate integration coefficients from received p
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    return 0;
}

// Print integrator information
void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method";
    if (useLNVD) stream << " + FLAC Local Non-Viscous Damping";
    if (useSMS)  stream << (useConsistent ? " + consistent (Olovsson) mass scaling"
                                          : " + lumped selective mass scaling");
    stream << " (classTag " << this->getClassTag() << ")\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  Damping parameter p: " << p << "\n";
    stream << "  Integration coefficients: q0 = " << q0
           << ", q1 = " << q1 << ", q2 = " << q2 << "\n";
    if (useLNVD)
        stream << "  alpha (FLAC local damping): " << alpha_flac << "\n";
    if (useSMS) {
        stream << "  -sms dtTarget = " << dtTarget
               << ", -maxAddedMass = " << maxAddedMassFrac << "\n";
        if (useConsistent)
            stream << "  M_tilde = M_lump + sum_e beta_e[diag(m)-mm^T/M_e], matrix-free PCG "
                      "(tol=" << pcgTol << ", maxIt=" << pcgMaxIt << ") at both sub-steps\n";
    }
    if (compute_critical_timestep > 0) {
        stream << "  Critical timestep (damped): " << damped_minimum_critical_timestep << "\n";
        stream << "  Critical timestep (undamped): " << undamped_minimum_critical_timestep << "\n";
    }
}
