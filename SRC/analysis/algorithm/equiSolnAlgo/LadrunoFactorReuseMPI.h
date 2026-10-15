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

// Ladruno WP-188 -- the rank agreement for factorization reuse.
//
// This lives in its OWN translation unit, listed in OPS_MPI_PER_TARGET_SOURCES
// (root CMakeLists.txt), for the reason recorded at LadrunoContactAbort.h: the
// equiSolnAlgo sources are compiled ONCE into the OPS_Analysis OBJECT library with
// neither _PARALLEL_PROCESSING nor _PARALLEL_INTERPRETERS, so an `#ifdef _PARALLEL_*`
// written in LadrunoFactorReuse.cpp would compile to nothing in every binary,
// OpenSeesMP included. Compiled per target, the MP binaries get the MPI_Allreduce
// and the serial ones the identity.
//
// Why a reduction at all: whether a rank forms the tangent decides whether its
// distributed SOE asks MUMPS for job 5 (factor + solve) or job 3 (solve only). The
// inputs are rank-identical in practice (the test norms are computed from the
// broadcast solution), but a rank that disagreed would deadlock the collective
// factorization, so the decision is reduced, not assumed.
//
// Scope: MPI_COMM_WORLD. A `system Mumps -commSplit` run (several independent
// models on sub-communicators) is not supported with reuse on (ADR-188 §5).

#ifndef LadrunoFactorReuseMPI_h
#define LadrunoFactorReuseMPI_h

// In place: flags[i] = max over all ranks. Serial builds, MPI not initialised or
// already finalised, or a single rank: no-op. Returns 0.
int ladrunoFactorReuseAgreeMax(int *flags, int n);

// WP-193: the same for doubles (the measured costs of the -reuseCost model).
int ladrunoFactorReuseAgreeMaxDouble(double *vals, int n);

// This process's rank in MPI_COMM_WORLD (0 in serial builds or outside MPI).
int ladrunoFactorReuseRank(void);

#endif
