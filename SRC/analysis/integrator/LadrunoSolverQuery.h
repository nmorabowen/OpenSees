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
#include <DirectIntegrationAnalysis.h>   // Ladruno WP-240: ladrunoDR relax
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <string.h>
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>     // Ladruno WP-240: drInstallMPIReductions
#endif

namespace Ladruno {

// Ladruno WP-240 (from jaabell/OpenSees 3c2f585b5, J. Abell): point the
// LadrunoDynamicRelaxation reductions at MPI_Allreduce. Called where each engine
// creates the integrator. It compiles to nothing outside the
// _PARALLEL_INTERPRETERS targets (OpenSeesMP.exe, openseesmp.pyd), which are the
// only ones that include this header with MPI. One copy for both engines, for
// the reason at the top of this file. Upstream sets the pointers inline in
// SRC/tcl/commands.cpp only; this fork also has the Python engine of
// openseesmp.pyd, which never runs that file.
inline void drInstallMPIReductions(void)
{
#ifdef _PARALLEL_INTERPRETERS
    int up = 0, down = 0;
    MPI_Initialized(&up);
    MPI_Finalized(&down);
    if (!up || down)
        return;                 // no MPI to reduce over: keep the identity
    drAllreduceSum = [](double x) {
        double y = 0.0;
        MPI_Allreduce(&x, &y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return y;
    };
    drAllreduceMax = [](double x) {
        double y = 0.0;
        MPI_Allreduce(&x, &y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return y;
    };
#endif
}

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
//
// Ladruno WP-240 (from jaabell/OpenSees 9a429237e, J. Abell):
//   ladrunoDR maxDisplacement -> max |u| over the equation DOFs
//   ladrunoDR relax -tolResidual r <-tolKE k> -maxSteps n <-runaway f> <-uRef u>
//                   <-checkEvery m> <-dt dt>
//        march until rest, runaway or n steps in ONE call, instead of one
//        interpreter round trip per step:
//          rest     residualNorm < r (and kineticEnergy < k when given),
//                   checked every step;
//          runaway  max|u| > f * uRef, checked every m steps (default 100);
//                   uRef defaults to max|u| after the first step of this call.
//        Returns {status steps residualNorm kineticEnergy maxDisplacement} with
//        status 0 rest, 1 maxSteps, 2 runaway, 3 failed step. dt defaults to 1.
//        The three values are global in a parallel run, so every process takes
//        the same decision at the same step.
//        It loops the stock DirectIntegrationAnalysis::analyze(1, dt): no new
//        analysis type and no copy of the driver (ADR 21 section 7.3 stands).
//        The analysis is a PARAMETER for the reason the integrator is one.
inline int drRelax(LadrunoDynamicRelaxation *dr, DirectIntegrationAnalysis *ta)
{
    if (ta == 0) {
        opserr << "WARNING ladrunoDR relax - no transient analysis (use `analysis Transient`)\n";
        return -1;
    }
    double tolR = -1.0, tolKE = -1.0, runaway = 0.0, uRef = -1.0, dt = 1.0;
    int maxSteps = -1, checkEvery = 100;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        int nd = 1;
        int bad = 0;
        if (OPS_GetNumRemainingInputArgs() < 1 && opt != 0 && opt[0] == '-') {
            opserr << "WARNING ladrunoDR relax - " << opt << " expects a value\n";
            return -1;
        }
        if      (strcmp(opt, "-tolResidual") == 0) bad = OPS_GetDoubleInput(&nd, &tolR);
        else if (strcmp(opt, "-tolKE") == 0)       bad = OPS_GetDoubleInput(&nd, &tolKE);
        else if (strcmp(opt, "-runaway") == 0)     bad = OPS_GetDoubleInput(&nd, &runaway);
        else if (strcmp(opt, "-uRef") == 0)        bad = OPS_GetDoubleInput(&nd, &uRef);
        else if (strcmp(opt, "-dt") == 0)          bad = OPS_GetDoubleInput(&nd, &dt);
        else if (strcmp(opt, "-maxSteps") == 0)    bad = OPS_GetIntInput(&nd, &maxSteps);
        else if (strcmp(opt, "-checkEvery") == 0)  bad = OPS_GetIntInput(&nd, &checkEvery);
        else {
            opserr << "WARNING ladrunoDR relax - unknown option '" << opt << "'\n";
            return -1;
        }
        if (bad < 0) {
            opserr << "WARNING ladrunoDR relax - failed to read the value of " << opt << "\n";
            return -1;
        }
    }
    if (!(tolR > 0.0) || maxSteps < 1) {
        opserr << "WARNING ladrunoDR relax - needs -tolResidual r > 0 and -maxSteps n >= 1\n";
        return -1;
    }
    if (!(dt > 0.0)) {
        opserr << "WARNING ladrunoDR relax - -dt must be positive\n";
        return -1;
    }
    if (checkEvery < 1) checkEvery = 1;
    // uRef by default: max|u| after the first step, when the integrator has been
    // (re)seeded from the domain
    const bool uRefFromFirstStep = runaway > 0.0 && !(uRef > 0.0);

    int status = 1, steps = 0;
    while (steps < maxSteps) {
        if (ta->analyze(1, dt) < 0) { status = 3; ++steps; break; }
        ++steps;
        if (steps == 1 && uRefFromFirstStep) {
            uRef = dr->getMaxDisplacement();
            if (!(uRef > 0.0)) uRef = 1.0e-12;  // from rest at u = 0
        }
        if (dr->getResidualNorm() < tolR && (tolKE < 0.0 || dr->getKineticEnergy() < tolKE)) {
            status = 0; break;
        }
        if (runaway > 0.0 && steps % checkEvery == 0
            && dr->getMaxDisplacement() > runaway * uRef) {
            status = 2; break;
        }
    }
    double out[5] = {(double)status, (double)steps, dr->getResidualNorm(),
                     dr->getKineticEnergy(), dr->getMaxDisplacement()};
    int n = 5;
    OPS_SetDoubleOutput(&n, out, false);
    return 0;
}

inline int drCommand(TransientIntegrator *ti, DirectIntegrationAnalysis *ta = 0)
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
                  "(residualNorm|kineticEnergy|stabilityMargin|maxDisplacement|relax)\n";
        return -1;
    }
    const char *sub = OPS_GetString();
    if (strcmp(sub, "relax") == 0)               // Ladruno WP-240
        return drRelax(dr, ta);
    double out = 0.0;
    if      (strcmp(sub, "residualNorm") == 0)  out = dr->getResidualNorm();
    else if (strcmp(sub, "kineticEnergy") == 0) out = dr->getKineticEnergy();
    else if (strcmp(sub, "maxDisplacement") == 0) out = dr->getMaxDisplacement();   // Ladruno WP-240
    // Ladruno (#728, merged into the shared dispatch here): the -massSafety
    // silent-divergence detector. It landed in OpenSeesCommands.cpp while this
    // branch was extracting that function into this header, so the two PRs
    // conflicted on exactly this line -- and a merge that took the extraction
    // wholesale would have DROPPED the subcommand while compiling cleanly.
    else if (strcmp(sub, "stabilityMargin") == 0) out = dr->getStabilityMargin();
    else {
        opserr << "WARNING ladrunoDR - unknown subcommand '" << sub
               << "' (use residualNorm|kineticEnergy|stabilityMargin|maxDisplacement|relax)\n";
        return -1;
    }
    int n = 1;
    OPS_SetDoubleOutput(&n, &out, true);
    return 0;
}

} // namespace Ladruno

#endif
