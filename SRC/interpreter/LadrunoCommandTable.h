// Ladruno WP-168 (apeGmsh#1490, panel P9 R10): THE table of fork-only
// interpreter commands. X-macro: deliberately no include guard. Every includer
// defines LADRUNO_COMMAND(name, dl, classic), includes this file, and undefines
// it again. Only the three Ladruno_registerCommands hooks include it:
//
//   SRC/tcl/LadrunoCommandsClassicTcl.h          classic Tcl (commands.cpp:
//                                                 OpenSees.exe, OpenSeesSP, OpenSeesMP)
//   SRC/interpreter/LadrunoCommandsTclWrapper.h  DL Tcl engine (TclWrapper.cpp)
//   SRC/interpreter/LadrunoCommandsPython.h      Python (PythonWrapper.cpp:
//                                                 opensees.pyd, openseesmp.pyd)
//
// One row per command:
//
//   LADRUNO_COMMAND("name", dl, classic)
//
//   dl       the no-arg OPS_* entry point the two DL engines (TclWrapper,
//            Python) call after plumbing argv into the elementAPI cursor. Every
//            fork command in those engines has exactly that shape, so the hook
//            generates their bridge. LADRUNO_NONE = not registered there.
//   classic  the classic-Tcl Tcl_CmdProc defined in SRC/tcl/commands.cpp (or
//            declared in SRC/tcl/commands.h). Classic Tcl keeps its own bridges
//            on purpose: some bind ENGINE STATE the shared OPS_* form cannot see
//            (ladrunoArcLength/ladrunoDR pass theStaticIntegrator /
//            theTransientIntegrator; the no-arg OPS_Ladruno*Cmd() forms read the
//            DL-only `cmds` singleton and would answer nothing, #729), and the DL
//            entry points live in OpenSeesCommands.cpp, which must never be pulled
//            into the classic link (LadrunoSolverQuery.h, trap 1).
//            LADRUNO_NONE = not registered in classic Tcl.
//
// A column is expanded ONLY in its own engine's translation unit, so a classic
// proc name never reaches the Python build and a DL entry point never reaches
// the classic link.
//
// Adding a command: one row here, nothing in the three upstream files.
// tests/test_ladruno_command_registry.py proves that OpenSees.exe and
// opensees.pyd answer every row exactly as its columns say (present AND
// absent), and fails on any Ladruno command registered outside these hooks.
//
// The LADRUNO_NONE gaps are recorded, not accidental. Closing one means writing
// the missing bridge and flipping the column; the parity test then demands it.

// ---- ADR-39 / ADR-41 / ADR-57: contact family ------------------------------
LADRUNO_COMMAND("contactSurface",            OPS_LadrunoContactSurface,        ladrunoContactSurface)
LADRUNO_COMMAND("contact",                   OPS_LadrunoContact,               ladrunoContact)
LADRUNO_COMMAND("contactPlane",              OPS_LadrunoContactPlane,          ladrunoContactPlane)
LADRUNO_COMMAND("ladrunoContactInfo",        OPS_LadrunoContactInfo,           ladrunoContactInfo)
LADRUNO_COMMAND("ladrunoContactForce",       OPS_LadrunoContactForce,          ladrunoContactForce)
LADRUNO_COMMAND("ladrunoMortarPenetration",  OPS_LadrunoMortarPenetration,     ladrunoMortarPenetration)
LADRUNO_COMMAND("ladrunoMortarTieResidual",  OPS_LadrunoMortarTieResidual,     ladrunoMortarTieResidual)
LADRUNO_COMMAND("ladrunoEdgePenetration",    OPS_LadrunoEdgePenetration,       ladrunoEdgePenetration)
LADRUNO_COMMAND("ladrunoBeginAugment",       OPS_LadrunoBeginAugment,          ladrunoBeginAugment)
LADRUNO_COMMAND("ladrunoEndAugment",         OPS_LadrunoEndAugment,            ladrunoEndAugment)

// ---- ADR-30 / ADR-62: mesh ties ---------------------------------------------
LADRUNO_COMMAND("ladrunoProjectionTieForce", OPS_LadrunoProjectionTieForce,    LADRUNO_NONE)
LADRUNO_COMMAND("LadrunoTie",                OPS_LadrunoTie,                   LADRUNO_NONE)

// ---- provenance, runtime knobs, diagnostics ---------------------------------
LADRUNO_COMMAND("ladrunoBuild",              OPS_LadrunoBuild,                 ladrunoBuild)
LADRUNO_COMMAND("ladrunoThreads",            OPS_LadrunoThreads,               ladrunoThreads)                // WP-107
LADRUNO_COMMAND("ladrunoMutation",           OPS_LadrunoMutation,              ladrunoMutation)               // ADR-87 D2
LADRUNO_COMMAND("ladrunoSANISANDReplay",     OPS_LadrunoSANISANDReplay,        ladrunoSANISANDReplay)         // WP-127
LADRUNO_COMMAND("profiler",                  OPS_profiler,                     TclCommand_profiler)
LADRUNO_COMMAND("ladrunoNumbering",          LADRUNO_NONE,                     TclCommand_ladrunoNumbering)   // ADR-74 N0

// ---- ADR-52 W1-I1b: trial-state probes --------------------------------------
LADRUNO_COMMAND("ladrunoTrialResidualNorm",  OPS_LadrunoTrialResidualNorm,     LADRUNO_NONE)
LADRUNO_COMMAND("ladrunoSetNodeTrial",       OPS_LadrunoSetNodeTrial,          LADRUNO_NONE)

// ---- solver-state queries (ADR-20 / ADR-31 / ADR-80) ------------------------
LADRUNO_COMMAND("ladrunoArcLength",          OPS_LadrunoArcLengthCmd,          ladrunoArcLength)              // #729
LADRUNO_COMMAND("ladrunoDR",                 OPS_LadrunoDRCmd,                 ladrunoDR)                     // #729
LADRUNO_COMMAND("ladrunoLoadControl",        OPS_LadrunoLoadControlCmd,        LADRUNO_NONE)                  // ADR-80 S1

// ---- analysis drivers --------------------------------------------------------
LADRUNO_COMMAND("LadrunoStaggeredAnalyze",   OPS_LadrunoStaggeredAnalyze,      TclCommand_ladrunoStaggeredAnalyze)   // ADR-73 P2
LADRUNO_COMMAND("criticalTimeStep",          OPS_criticalTimeStep,             LADRUNO_NONE)
LADRUNO_COMMAND("complexEigen",              OPS_complexEigen,                 LADRUNO_NONE)                  // ADR-46
LADRUNO_COMMAND("modalResponseHistory",      OPS_LadrunoModalResponseHistory,  modalResponseHistory)          // ADR-44 P1a
LADRUNO_COMMAND("frequencyResponse",         OPS_LadrunoFrequencyResponse,     frequencyResponse)             // ADR-44 P2
LADRUNO_COMMAND("steadyStateDynamics",       OPS_LadrunoSteadyStateDynamics,   steadyStateDynamics)           // ADR-44 P2
LADRUNO_COMMAND("randomResponse",            OPS_LadrunoRandomResponse,        randomResponse)                // ADR-44 P3
