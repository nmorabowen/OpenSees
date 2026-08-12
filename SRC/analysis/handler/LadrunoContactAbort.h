// LADRUNO-HEADER-START
// ==========================================================================
//  Ladruno — a research fork of OpenSees
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-78 P1 — the single exit for a contact failure inside
// LadrunoContactHandler::handle().
//
// This lives in its OWN translation unit for one reason: it is the only
// parallel-sensitive code in the contact subsystem, and it MUST be compiled
// with the consuming target's `_PARALLEL_PROCESSING` / `_PARALLEL_INTERPRETERS`
// define. LadrunoContactHandler.cpp is compiled into the OPS_Analysis OBJECT
// library, which is built ONCE with no parallel define and folded into every
// target -- so an `#ifdef _PARALLEL_*` block written there compiles to nothing
// in every build, including OpenSeesMP. That was measured, not assumed: the
// guarded MPI_Abort silently vanished and the mutation deck hung exactly as it
// had before the "fix".
//
// This file is therefore listed in OPS_MPI_PER_TARGET_SOURCES and added to
// each executable/module directly, mirroring Patch 8
// (OPS_INTERP_TCL_PER_TARGET_SOURCES) and Patch 9 (OPS_InterpPyCmds*). Keeping
// it to ~one function keeps the duplicated compilation trivial.

#ifndef LadrunoContactAbort_h
#define LadrunoContactAbort_h

// Returns -1 always, so callers read as `return ladrunoContactFatal();`.
//
// Serial builds: that is all it does -- StaticAnalysis /
// DirectIntegrationAnalysis::domainChanged() report the failure and the run
// stops, exactly as before.
//
// MPI builds with np > 1: returning is NOT enough. The failing rank leaves
// handle() while its peers walk into the next collective (the numberer's
// gather) and block on a rank that never arrives -- the job hangs until an
// external timeout kills it, stranding orphaned processes. On a 240-rank run
// that burns the whole allocation and buries the FATAL line in one log. So the
// job is torn down instead.
int ladrunoContactFatal();

#endif
