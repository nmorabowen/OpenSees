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

// Ladruno (#729): the ONE copy of the `ladrunoArcLength` / `ladrunoDR` subcommand
// dispatch, shared by the two command engines.
//
// It lives in a header, not in a .cpp, for a link reason that is not obvious and
// cost a failed build to find. The natural home was OpenSeesCommands.cpp, and the
// CMake notes even say `OPS_InterpPyCmds` is linked by the Tcl OpenSees/SP/MP
// exes -- but being ON THE LINK LINE is not the same as being LINKED IN. Nothing
// in the classic engine referenced any symbol in that TU, so the linker never
// pulled the object at all. The moment SRC/tcl/commands.cpp externs one function
// from it, the object IS pulled, and it brings the whole file: ~40 LNK2005
// duplicate-symbol errors against elementAPI_TCL.cpp (both define the ops_get*_ /
// ops_set*_ elementAPI backend) plus unresolved OPS_SparsePython* externals.
// Measured 2026-08-11.
//
// A header-only inline core has no object of its own, so it adds no link surface,
// and each including TU binds the elementAPI calls to ITS OWN backend --
// elementAPI_TCL.cpp for classic Tcl, elementAPI.cpp for the DL engine -- which is
// exactly the behaviour wanted: the same dispatch reading and publishing through
// whichever interpreter is asking. Same extraction precedent as
// LadrunoFictitiousMass.h / LadrunoJ2Kernel.h.
//
// The integrator is a PARAMETER, and that is load-bearing. The DL engine's
// OPS_Ladruno*Cmd() wrappers read the `cmds` singleton, which only that engine
// constructs; classic Tcl owns its integrators in commands.cpp's own
// theStaticIntegrator / theTransientIntegrator globals. A classic bridge that
// called the singleton form would compile, link, and silently answer NOTHING
// (`cmds == 0` -> return 0 -> no output written).
//
// Header-only. Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoSolverQuery_h
#define LadrunoSolverQuery_h

#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <LadrunoArcLength.h>
#include <LadrunoDynamicRelaxation.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <string.h>

namespace Ladruno {

// ladrunoArcLength [sub [value]] -- runtime control of the active LadrunoArcLength
// static integrator (Layer-B cut-and-retry, ADR-20 §4.3). Exposes the scalar radius
// mutators and read-only state to a script so a failed step can be re-attempted
// without reconstructing the integrator:
//   ok = analyze(1)
//   while ok != 0:
//       ladrunoArcLength('reduceStep', 0.5);  ok = analyze(1)
// Subcommands (value-taking): reduceStep f | increaseStep f | setArcLength v |
//   scaleCVisc f. Queries (return a scalar): arcLength | deltaLambdaStep |
//   currentLambda | sign | deltaUstepNorm | dissipationRatio | dissipatedEnergy |
//   referenceEnergy. With no subcommand, returns arcLength.
inline int arcLengthCommand(StaticIntegrator *si)
{
    if (si == 0 || si->getClassTag() != INTEGRATOR_TAGS_LadrunoArcLength) {
        opserr << "WARNING ladrunoArcLength - the active static integrator is not "
                  "a LadrunoArcLength (set `integrator LadrunoArcLength ...` first)\n";
        return -1;
    }
    LadrunoArcLength *la = (LadrunoArcLength *)si;

    // no subcommand => return the current arc length
    if (OPS_GetNumRemainingInputArgs() < 1) {
        double al = la->getArcLength();
        int nd = 1;
        OPS_SetDoubleOutput(&nd, &al, true);
        return 0;
    }

    const char *sub = OPS_GetString();

    // value-taking mutators
    if (strcmp(sub, "reduceStep") == 0 || strcmp(sub, "increaseStep") == 0 ||
        strcmp(sub, "setArcLength") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING ladrunoArcLength " << sub << " - expects one value\n";
            return -1;
        }
        int nd = 1;
        double v = 0.0;
        if (OPS_GetDoubleInput(&nd, &v) < 0) {
            opserr << "WARNING ladrunoArcLength " << sub << " - failed to read value\n";
            return -1;
        }
        int res = 0;
        if (strcmp(sub, "reduceStep") == 0)        res = la->reduceStep(v);
        else if (strcmp(sub, "increaseStep") == 0) res = la->increaseStep(v);
        else                                       res = la->setArcLength(v);
        if (res < 0) return -1;
        double al = la->getArcLength();            // echo the resulting radius
        int n = 1;
        OPS_SetDoubleOutput(&n, &al, true);
        return 0;
    }

    // value-taking actuator (ADR-31 R-RAMPDOWN): scale the viscous coefficient;
    // echoes the resulting dissipation ratio so the driver sees the effect.
    if (strcmp(sub, "scaleCVisc") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING ladrunoArcLength scaleCVisc - expects one value\n";
            return -1;
        }
        int nd = 1;
        double v = 0.0;
        if (OPS_GetDoubleInput(&nd, &v) < 0) {
            opserr << "WARNING ladrunoArcLength scaleCVisc - failed to read value\n";
            return -1;
        }
        if (la->scaleCVisc(v) < 0) return -1;
        double r = la->getStabilizationDissipationRatio();
        int n = 1;
        OPS_SetDoubleOutput(&n, &r, true);
        return 0;
    }

    // read-only queries
    double out = 0.0;
    if      (strcmp(sub, "arcLength") == 0)        out = la->getArcLength();
    else if (strcmp(sub, "deltaLambdaStep") == 0)  out = la->getDeltaLambdaStep();
    else if (strcmp(sub, "currentLambda") == 0)    out = la->getCurrentLambda();
    else if (strcmp(sub, "sign") == 0)             out = (double)la->getSignLastDeltaLambdaStep();
    else if (strcmp(sub, "deltaUstepNorm") == 0)   out = la->getDeltaUstepNorm();
    // ADR-31 Layer-1.5 stabilization-energy gate
    else if (strcmp(sub, "dissipationRatio") == 0) out = la->getStabilizationDissipationRatio();
    else if (strcmp(sub, "dissipatedEnergy") == 0) out = la->getStabilizationDissipatedEnergy();
    else if (strcmp(sub, "referenceEnergy") == 0)  out = la->getReferenceStrainEnergy();
    else {
        opserr << "WARNING ladrunoArcLength - unknown subcommand '" << sub
               << "' (use reduceStep|increaseStep|setArcLength|scaleCVisc|arcLength|"
                  "deltaLambdaStep|currentLambda|sign|deltaUstepNorm|"
                  "dissipationRatio|dissipatedEnergy|referenceEnergy)\n";
        return -1;
    }
    int n = 1;
    OPS_SetDoubleOutput(&n, &out, true);
    return 0;
}

// ladrunoDR <sub> -- runtime query of the active LadrunoDynamicRelaxation (33005);
// the rung-5 settling / micro-burst signals (ADR-31 R-DR-ENERGY). Read-only.
//   ladrunoDR residualNorm    -> ||f_ext - f_int||_inf  (mass-free settling gate)
//   ladrunoDR kineticEnergy   -> 1/2 v^T M* v           (micro-burst signal)
//   ladrunoDR stabilityMargin -> (omega_max*dt/2)^2 of the mass in use, measured
//        against the LIVE tangent (note 83 §3, #728). <= 1 stable; == massSafety^2
//        on an unchanged tangent; > 1 = marching at/over the explicit boundary,
//        where DR amplifies round-off and can relax to a silently WRONG state.
//        NEGATIVE = NOT MEASURED, and must never be read as "safe": -1 = no
//        gershgorin mass (lumped/unity have no such bound), -2 = probe disabled
//        (-marginEvery 0).
// The robust-solve driver reads residualNorm/residualNorm0 each DR step to decide
// when the dynamics excursion has relaxed to a quasi-static rest state -- the
// physical-mass EnergyBalance KE is ~0 on DR's pseudo-mass models, so the gate
// MUST be this force residual, not a KE ratio.
inline int drCommand(TransientIntegrator *ti)
{
    if (ti == 0 || ti->getClassTag() != INTEGRATOR_TAGS_LadrunoDynamicRelaxation) {
        opserr << "WARNING ladrunoDR - the active transient integrator is not a "
                  "LadrunoDynamicRelaxation (set `integrator LadrunoDynamicRelaxation "
                  "...` first)\n";
        return -1;
    }
    LadrunoDynamicRelaxation *dr = (LadrunoDynamicRelaxation *)ti;

    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING ladrunoDR - expects a subcommand "
                  "(residualNorm|kineticEnergy|stabilityMargin)\n";
        return -1;
    }
    const char *sub = OPS_GetString();
    double out = 0.0;
    if      (strcmp(sub, "residualNorm") == 0)  out = dr->getResidualNorm();
    else if (strcmp(sub, "kineticEnergy") == 0) out = dr->getKineticEnergy();
    // Ladruno (#728, merged into the shared dispatch here): the -massSafety
    // silent-divergence detector. It landed in OpenSeesCommands.cpp while this
    // branch was extracting that function into this header, so the two PRs
    // conflicted on exactly this line -- and a merge that took the extraction
    // wholesale would have DROPPED the subcommand while compiling cleanly.
    else if (strcmp(sub, "stabilityMargin") == 0) out = dr->getStabilityMargin();
    else {
        opserr << "WARNING ladrunoDR - unknown subcommand '" << sub
               << "' (use residualNorm|kineticEnergy|stabilityMargin)\n";
        return -1;
    }
    int n = 1;
    OPS_SetDoubleOutput(&n, &out, true);
    return 0;
}

} // namespace Ladruno

#endif
