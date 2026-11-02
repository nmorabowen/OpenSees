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

// Ladruno WP-221: `ladrunoSymmetryGuard` -- read the half-storage asymmetry
// guard of the active linear system (LadrunoSymmetryGuard.h) from a script.
//
// Header-only for the same link reason as LadrunoSolverQuery.h (#729): the
// classic Tcl engine (SRC/tcl/commands.cpp) and the DL engine
// (OpenSeesCommands.cpp: openseespy, DL Tcl) both include this one copy, and each
// binds the elementAPI calls to its own backend. The SOE is a PARAMETER because
// each engine keeps its own `theSOE`.
//
//   ladrunoSymmetryGuard            -> the full record, in this order:
//       latched worstDeviation scale assembly checkedAssemblies checks
//       strict matrixType maxRelative refused
//   ladrunoSymmetryGuard <field>    -> that one value
//
//   latched            1 once an element matrix exceeded tolerance (1e-8 relative)
//   worstDeviation     max |m(i,j) - m(j,i)| of the matrix that latched it
//   scale              max |m(i,j)| of that matrix
//   assembly           tangent assembly of the current pattern where it latched;
//                      -1 while unlatched
//   checkedAssemblies  tangent assemblies in which the guard actually compared
//                      at least one matrix. 0 means "never looked", NOT "clean":
//                      without -strictSymmetry the guard SAMPLES (3 of every 64
//                      assemblies), so an assembly can go unchecked
//   checks             element matrices compared
//   strict             1 under -strictSymmetry
//   matrixType         what is half-stored: Pardiso's -matrixType (0 = full CSR,
//                      the guard is off); ProfileSPD reports 1 (always upper)
//   maxRelative        worst worstDev/scale over every check, latched or not
//   refused            solves refused under -strictSymmetry
//
// It fails (non-zero, so openseespy raises) when there is no system yet or the
// system is not one that half-stores with a guard -- a script asking "was my
// tangent symmetric?" must never get a silent "0" from a solver that never
// looked.

#ifndef LadrunoSymmetryGuardQuery_h
#define LadrunoSymmetryGuardQuery_h

#include <LinearSOE.h>
#include <ProfileSPDLinSOE.h>
#ifdef _PARDISO
#include <PARDISOGenLinSOE.h>
#endif
#include <LadrunoSymmetryGuard.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <string.h>

namespace Ladruno {

inline int symmetryGuardCommand(LinearSOE *soe)
{
    if (soe == 0) {
        opserr << "WARNING ladrunoSymmetryGuard - no system has been defined yet\n";
        return -1;
    }

    const LadrunoSymmetryGuard *g = 0;
    int matType = 0;
    if (soe->getClassTag() == LinSOE_TAGS_ProfileSPDLinSOE) {
        g = &((ProfileSPDLinSOE *)soe)->getSymmetryGuard();
        matType = 1;
    }
#ifdef _PARDISO
    else if (soe->getClassTag() == LinSOE_TAGS_PARDISOGenLinSOE) {
        PARDISOGenLinSOE *p = (PARDISOGenLinSOE *)soe;
        g = &p->getSymmetryGuard();
        matType = p->getMatType();
    }
#endif
    if (g == 0) {
        opserr << "WARNING ladrunoSymmetryGuard - the active system has no "
                  "half-storage asymmetry guard (only `system Pardiso` and "
                  "`system ProfileSPD` carry one)\n";
        return -1;
    }

    const int nField = 10;
    double rec[nField] = {
        (double)g->latched, g->worstDev, g->scale, (double)g->pass,
        (double)g->checkedPasses, (double)g->checks, (double)g->strict,
        (double)matType, g->maxRel, (double)g->refused};
    static const char *names[nField] = {
        "latched", "worstDeviation", "scale", "assembly", "checkedAssemblies",
        "checks", "strict", "matrixType", "maxRelative", "refused"};

    if (OPS_GetNumRemainingInputArgs() < 1) {
        int n = nField;
        if (OPS_SetDoubleOutput(&n, rec, false) < 0) return -1;
        return 0;
    }

    const char *sub = OPS_GetString();
    for (int i = 0; i < nField; i++)
        if (strcmp(sub, names[i]) == 0) {
            int n = 1;
            if (OPS_SetDoubleOutput(&n, &rec[i], true) < 0) return -1;
            return 0;
        }
    if (strcmp(sub, "tolerance") == 0) {
        double tol = LadrunoSymmetryGuard::tolerance();
        int n = 1;
        if (OPS_SetDoubleOutput(&n, &tol, true) < 0) return -1;
        return 0;
    }
    opserr << "WARNING ladrunoSymmetryGuard - unknown field '" << sub
           << "' (use latched|worstDeviation|scale|assembly|checkedAssemblies|"
              "checks|strict|matrixType|maxRelative|refused|tolerance)\n";
    return -1;
}

}  // namespace Ladruno

#endif
