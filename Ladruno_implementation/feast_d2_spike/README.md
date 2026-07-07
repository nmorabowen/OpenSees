---
title: "D2 spike record — MKL FEAST RCI over an MPI sub-communicator via MumpsParallelSOE"
project: Ladruno
type: spike-record
related:
  - "[[43_ladruno_feast_eigensolver_adr]]"
  - "[[45_ladruno_modal_family_roadmap_adr]]"
  - "[[LEDGER_quirks]]"
date: 2026-07-07
---

# D2 spike — decision gate before ADR 43 P-C

**Question:** can MKL FEAST RCI (`dfeast_srci`) route its per-contour-point complex
shifted linear solves through the fork's existing `MumpsParallelSOE` on an
`MPI_Comm_split` sub-communicator, or does the complex-shift problem force vendoring
PFEAST 4.0?

This was a **spike, not a PR** — no C++ was built or run. The MPI toy-model run
originally scoped (~1k-DOF frame, 2–4 ranks) turned out to be unnecessary: reading the
actual `MumpsParallelSolver`/`MumpsParallelSOE` source and the local MUMPS build state
surfaced two independent, concrete blockers that settle the question without needing a
runtime experiment. A numpy check filled the one genuinely open numerical question
(does the ADR's "real solver per conjugate pair" claim hold).

## Finding 1 — `MumpsParallelSolver` hardcodes `MPI_COMM_WORLD`

`MumpsParallelSolver::MumpsParallelSolver(int mpi_comm, int ICNTL7, int ICNTL14)`
(`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp:54-64`) accepts `mpi_comm`
but never stores it — no member is set. `initializeMumps()` (`:93-105`) hardcodes
`id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD)` on the Intel-MPI path this fork's
Windows/oneAPI build uses (`_OPENMPI` path is `0`, MUMPS's own WORLD); the
`MPI_Comm_rank/size` probe two lines later also reads `MPI_COMM_WORLD` directly. No
other path in the file reaches a non-WORLD communicator. **Verified by direct source
read**, cross-checked by an Opus adversarial-verification pass.

A bounded plumbing bug, not an architectural wall — recorded as a load-bearing gotcha
in [[LEDGER_quirks]] and scoped as ADR 43 P3a.

## Finding 2 — no complex MUMPS build; `dmumps` only

`mumps-build/CMakeCache.txt:417` → `arith:UNINITIALIZED=d` (double real only).
`mumps-install/lib/` contains `dmumps.lib`, `mumps_common.lib`, `pord.lib` — no
`zmumps.lib`. `MumpsParallelSolver.{h,cpp}` `#include <dmumps_c.h>` exclusively, use
`DMUMPS_STRUC_C`, call only `dmumps_c()`. There is no complex code path in the SRC
mumps wrapper today.

## Finding 3 — ADR 43 §5.2's "real solver per conjugate pair" claim was wrong

Refuted by an Opus literature cross-check against the FEAST v3/v4 User Guides and the
MKL `?feast_srci` reference: at `ijob=10/11` FEAST hands the caller a **genuinely
complex** workspace and asks for a complex factorization/solve of $A_z=z_jM-K$ — every
individual contour solve is complex. The Hermitian conjugate-pairing (`fpm(2)`) only
halves the *count* of solves (skip $\bar z_j$, conjugate the result); it does not make
any one solve real. Even MKL's convenience driver `dfeast_scsrgv` doesn't dodge this —
it calls complex PARDISO (`mtype=6`) internally.

`block_real_conjugate_check.py` (this folder) numerically confirms the *correct*
version of the trick: for real symmetric $K,M$ and complex shift $z=a+bi$, solving
$(zM-K)X=B$ (real $B$) is exact (1e-15 relative error against a direct complex solve)
to forming the symmetric real $2n\times2n$ block-augmented system

$$
\begin{bmatrix}aM-K & -bM\\-bM & -(aM-K)\end{bmatrix}
\begin{bmatrix}X_r\\X_i\end{bmatrix}=\begin{bmatrix}B\\0\end{bmatrix}
$$

and recombining $X=X_r+iX_i$. This is a real solve **twice the size** of the original
problem (roughly 2× the nonzero count), not the original $n\times n$ `MumpsParallelSOE`
call with a swapped-in complex shift.

Run it: `python block_real_conjugate_check.py` (needs numpy; no OpenSees build
required).

## Decision — Outcome A-prime (see ADR 43 §9 R1, umbrella §5 D2)

Stay on MKL FEAST as the P-C driver. **Do not vendor PFEAST 4.0** — neither finding
argues for it; PFEAST buys nothing these two fixes don't, and adds a Fortran build
dependency. Instead, P3 splits into three sub-phases:

- **P3a** — fix the sub-communicator plumbing (Finding 1): store the passed comm,
  `MPI_Comm_c2f` it (not WORLD), thread it through `MumpsParallelSOE` construction.
  Gate: 2+ independent sub-comms each solve a distinct real system concurrently with no
  cross-talk (this *is* the originally-scoped D2 toy test — now a real P3a unit test
  instead of a throwaway spike).
- **P3b** — new inner SOE for the complex contour solve. Recommended: the symmetric
  real $2n\times2n$ block-augmentation (Finding 3), LDLᵀ via the existing real `dmumps`
  (`SYM=2`) — **zero new external dependency**. Record a cost comparison against a
  `zmumps`-based alternative in the P3 validation bundle before committing (block-real
  costs ~2× factorization work per contour point; don't assume it wins without
  measuring).
- **P3c** — `FeastEigenSolver` orchestration: `MPI_Comm_split` into per-quadrature
  sub-comms, each hosting a P3a/P3b inner SOE; `Allreduce`-sum the projector. Gate:
  strong/weak scaling, bit-comparable eigenpairs to the P1 serial run.

P4 (SP/MP build-flag unification) is unchanged by this spike.

## Verification

An Opus subagent independently re-read the source (not just my line numbers) and
cross-checked Findings 1–3 against the FEAST v3/v4 User Guides and the MKL
`?feast_srci`/`?feast_scsrgv` references. All three findings: **CONFIRMED**. Report
archived in this session's transcript; citations: FEAST v4 User Guide
(arxiv.org/pdf/2002.04807), MKL `?feast_srci`/`?feast_hcsrgv` reference, FEAST v3 User
Guide (arxiv.org/pdf/1203.4031).
