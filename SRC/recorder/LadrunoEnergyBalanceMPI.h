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

// Ladruno -- EnergyBalance -v2: optional columns agreed across ranks.
//
// This lives in its OWN translation unit, listed in OPS_MPI_PER_TARGET_SOURCES
// (root CMakeLists.txt), for the reason recorded at LadrunoFactorReuseMPI.h: the
// recorder sources are compiled ONCE into the OPS_Recorder OBJECT library with
// neither _PARALLEL_PROCESSING nor _PARALLEL_INTERPRETERS, so an `#ifdef _PARALLEL_*`
// written in EnergyBalanceRecorder.cpp would compile to nothing in every binary.
//
// Why a reduction at all: each -v2 optional column (E_inject, E_lnvd, E_modal,
// E_hg, KE_ms, DW_rel/W_ref) was switched on by a RANK-LOCAL test (the energy
// registry's declared set, the hourglass probe, the SMS registry). Under
// interpreter-per-rank MPI a rank that owns no DRM/absorbing element wrote and
// echoed one model column fewer than its peers, and per-rank readers that
// require identical layouts failed. OR-reducing the flags gives every rank the
// same layout; a rank with no local contribution writes 0 in that column.
//
// Scope: _PARALLEL_INTERPRETERS only (openseesmp, OpenSeesMP Tcl), where every
// rank builds its own recorder and reaches its first record() inside the same
// analyze call. PartitionedDomain targets (_PARALLEL_PROCESSING) get the
// identity: there the recorder lives on P0 only (its sendSelf is empty), so a
// collective would hang. COLLECTIVE: under openseesmp every rank must define
// the -v2 recorder, or the first record() blocks (quirks ledger).

#ifndef LadrunoEnergyBalanceMPI_h
#define LadrunoEnergyBalanceMPI_h

// In place: flags[i] = max over all ranks of MPI_COMM_WORLD. Serial builds,
// _PARALLEL_PROCESSING builds, MPI not initialised or already finalised, or a
// single rank: no-op. Returns 0.
int ladrunoEnergyBalanceAgreeMax(int *flags, int n);

#endif
