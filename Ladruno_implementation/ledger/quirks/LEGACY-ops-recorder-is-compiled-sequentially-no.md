---
wp: LEGACY
title: "OPS_Recorder is compiled sequentially (no _PARALLEL_PROCESSING) for ALL targets"
legacy_seq: 37
---
## OPS_Recorder is compiled sequentially (no `_PARALLEL_PROCESSING`) for ALL targets

`SRC/recorder/LadrunoRecorder.cpp` (and the frozen `MPCORecorder.cpp`) live in the shared
`OPS_Recorder` static lib, which CMake compiles **once with the sequential define
set** (no `-D_PARALLEL_PROCESSING`, no `-D_MUMPS`) and links into OpenSeesPy,
OpenSeesMP **and** openseesmp alike. Consequence: any `#ifdef _PARALLEL_PROCESSING`
code in a recorder is **dead** — it compiles out in every artifact, including the
MP ones. The frozen `MPCORecorder`'s only `_PARALLEL_PROCESSING` block (the
per-step stage-stamp `MPI_Allreduce`) is therefore never active in this build; its
per-partition output works purely via `sendSelf`/`recvSelf` + `send_self_count`
(always compiled, no MPI). A recorder that needs its MPI rank cannot call
`MPI_Comm_rank` — read the launcher's `PMI_RANK`/`PMI_SIZE` env instead (Intel MPI
and MS-MPI both export them per rank; verified `mpiexec -n N` sets them).

Also: the openseespy parallel model here is `_PARALLEL_INTERPRETERS` (one
independent interpreter per rank, manual `getPID`/`getNP` partition — see
`example_mpi_paralleltruss_explicit.py`), NOT `_PARALLEL_PROCESSING`
(PartitionedDomain + `sendSelf` broadcast). The `partition` command (auto domain
decomposition, the path that *would* exercise `sendSelf`) is **METIS-4-blocked**
(`OPS_partition` returns an error; needs METIS 5 / `OPS_HAVE_METIS5`) — but ONLY
the openseespy `partition` command under `_PARALLEL_INTERPRETERS`. **Correction
(WP-163 review, 2026-10-03):** the Tcl `OpenSeesSP.exe` partitions automatically at
the first `analyze`/`eigen` (`SRC/tcl/commands.cpp` `partitionModel()` → `new Metis`,
the METIS-4 legacy API in `MetisWrapper.cpp`), so the PartitionedDomain + `sendSelf`
recorder path DOES ship in this build — untested, not unreachable (e.g. `-G energy`
null-derefs on a ShadowSubdomain there; WP-163 M1).
