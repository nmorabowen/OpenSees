# OpenSees modal / eigenvalue analysis — ground-truth baseline (current state)

Audit of the upstream + Ladruno-fork source under `SRC/`, focused on how
modal/eigen analysis works today and what the PARALLEL (OpenSeesSP / OpenSeesMP /
MPI) story actually is. All citations are `file:line` against the worktree at
`C:\Users\nmora\Github\OpenSees_Compile\OpenSees\SRC`.

> TL;DR of the parallel reality (read section C for the full trace): `eigen` is
> wired end-to-end into the partitioned/subdomain machinery and the shift-invert
> `(K−σM)⁻¹` solve delegates to whatever `LinearSOE` is active — which *can* be a
> distributed `MumpsParallelSOE`. **BUT** the ARPACK Lanczos driver itself is the
> **serial** `dsaupd_`/`dseupd_` (no PARPACK), so the Krylov iteration runs on one
> rank with distributed linear-solves/matvecs underneath. There is **no
> distributed eigensolver** (no PARPACK / SLEPc / ScaLAPACK / distributed Lanczos).
> Practically: the parallel-eigen path exists but is lightly used, only reachable
> from the Tcl `_PARALLEL_PROCESSING` build (`MumpsParallelSOE` is gated to the
> `_PARALLEL_INTERPRETERS`/MP build), and the post-processing
> (`modalProperties`/`responseSpectrum`) is serial-only.

---

## (A) Serial eigen pipeline + flags

### A.1 The `eigen` command and how a solver is selected

Two parsers exist (modern "xara" runtime and legacy/parallel runtime); they agree:

- **xara runtime**: `SRC/runtime/commands/analysis/analysis.cpp:250` (`eigenAnalysis`).
  Default `typeSolver = EigenSOE_TAGS_ArpackSOE` (`analysis.cpp:269`). Flag map:
  - `frequency` / `-frequency` / `generalized` / `-generalized` → generalized GEVP, default (`analysis.cpp:277-281`)
  - `standard` / `-standard` → `generalizedAlgo=false` + `EigenSOE_TAGS_SymBandEigenSOE` (`analysis.cpp:283-287`)
  - `-findLargest` → `findSmallest=false` (`analysis.cpp:289-290`)
  - `genBandArpack` / `-genBandArpack` / `genBandArpackEigen` → `EigenSOE_TAGS_ArpackSOE` (`analysis.cpp:292-296`)
  - `symmBandLapack` / `-symmBandLapack` → `EigenSOE_TAGS_SymBandEigenSOE` (`analysis.cpp:298-302`)
  - `fullGenLapack` / `-fullGenLapack` → `EigenSOE_TAGS_FullGenEigenSOE` (`analysis.cpp:304-308`)
  - Then `builder->newEigenAnalysis(typeSolver, shift)` (`analysis.cpp:326`) and
    `builder->eigen(...)` (`analysis.cpp:328`).
- **legacy / openseespy / parallel runtime**: `SRC/interpreter/OpenSeesCommands.cpp:270`
  (`OpenSeesCommands::eigen`). Same tags, plus an extra case
  `EigenSOE_TAGS_SymmGeneralizedEigenSOE` that is **disabled on Windows**
  (`OpenSeesCommands.cpp:340-346`: "SymmGeneralizedEigenSolver not currently
  compiled for Windows").

SOE construction: `SRC/runtime/runtime/BasicAnalysisBuilder.cpp:871`
(`newEigenAnalysis`) builds `SymBandEigenSOE` / `FullGenEigenSOE` / `ArpackSOE`
(`BasicAnalysisBuilder.cpp:894-903`) and immediately wires the linear SOE into it
(`theEigenSOE->setLinearSOE(*theSOE)`, `BasicAnalysisBuilder.cpp:910`).

`shift` is currently hard-zero in the xara parser (`analysis.cpp:271`, never
reassigned) — there is no exposed `-shift`/`-sigma` flag in that path; the shift
plumbing exists at the SOE level (`ArpackSOE(double shift)`) but the command does
not surface it. (Gap, see E.)

### A.2 What each EigenSOE/EigenSolver actually does

| SOE (classTag) | Solver | Algorithm | Matrices | Selected by |
|---|---|---|---|---|
| `ArpackSOE` (default) | `ArpackSolver` | **ARPACK** symmetric implicitly-restarted Lanczos (`dsaupd_`/`dseupd_`), shift-invert mode 3, `bmat='G'` | wraps a **`LinearSOE`** for the `(K−σM)⁻¹` solve; M held separately (diagonal fast path or element matvec) | default / `genBandArpack` |
| `BandArpackSOE` | `BandArpackSolver` | **ARPACK** Lanczos + **LAPACK banded** factor/solve (`dgbtrf_`/`dgbtrs_`/`dgbsv_`) for shift-invert | banded K, separate M | (legacy; not in xara flag map) |
| `SymArpackSOE` | `SymArpackSolver` | ARPACK Lanczos + sparse-symmetric factor | sparse-sym | (legacy) |
| `SymBandEigenSOE` | `SymBandEigenSolver` | **LAPACK** dense-banded symmetric eigensolver (`dsbgvx`-family) — full spectrum, no Krylov | symmetric banded K, M | `standard` / `symmBandLapack` |
| `FullGenEigenSOE` | `FullGenEigenSolver` | **LAPACK** full general GEVP (`dggev`) | full dense K, M | `fullGenLapack` |
| `SymmGeneralizedEigenSOE` | `SymmGeneralizedEigenSolver` | LAPACK symmetric-definite GEVP | symmetric | only via `OpenSeesCommands.cpp` (non-Windows) |

ARPACK confirmation: `SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp:120-130`
(`dsaupd_`/`dseupd_` externs), driver loop `ArpackSolver.cpp:221-291`,
`mode=3` (shift-invert) at `ArpackSolver.cpp:207`, `bmat='G'` at
`ArpackSolver.cpp:200`. BandArpack LAPACK banded externs:
`BandArpackSolver.cpp:98-117`.

`ArpackSolver` **only solves the generalized problem** — it errors out if
`generalized==false` (`ArpackSolver.cpp:138-141`). The `standard`/`-standard`
flag therefore *must* route to `SymBandEigenSOE` (which it does), not ARPACK.

### A.3 Assembly path (serial)

`StaticAnalysis::eigen` (`SRC/analysis/analysis/StaticAnalysis.cpp:244`) and the
twin in `BasicAnalysisBuilder::eigen` (`BasicAnalysisBuilder.cpp:916`):
- `theEigenSOE->zeroA(); zeroM()` (`StaticAnalysis.cpp:275-276`)
- form K: per element `zeroTangent(); addKtToTang(1.0); addA(getTangent(0), id)`
  (`StaticAnalysis.cpp:285-293`)
- form M: per element/DOF `addMtoTang(1.0); addM(...)`
  (`StaticAnalysis.cpp:300-323`)
- `theEigenSOE->solve(numMode, generalized, findSmallest)` (`StaticAnalysis.cpp:330`)
- store eigenvalues + push eigenvectors into nodes via
  `theAnalysisModel->setEigenvector` (`StaticAnalysis.cpp:339-345`).

**KEY:** the assembled "A" is always `Kt` (`addKtToTang`) and "M" is always mass
(`addMtoTang`). There is **no Kg (geometric stiffness) assembly** in this path —
see section D.

`EigenIntegrator::formEleTangK` likewise calls only `addKtToTang(1.0)`
(`SRC/analysis/integrator/EigenIntegrator.cpp:191-196`). Algorithms:
`FrequencyAlgo` (generalized, `formK`+`formM`) and `StandardEigenAlgo` (standard,
`formK` only) under `SRC/analysis/algorithm/eigenAlgo/`.

---

## (B) modalProperties / responseSpectrum — capabilities & limits

Two implementations of `DomainModalProperties` exist:
`SRC/runtime/commands/analysis/modal/DomainModalProperties.{cpp,h}` (active, ASDEA
Petracca/Camata) and an older sibling `SRC/domain/domain/DomainModalProperties.*`.
Command registration is via `G3_AddTclDomainCommands`
(`SRC/runtime/commands/domain/commands.cpp:121`, `modalProperties`) and in the
legacy Tcl table (`SRC/tcl/commands.cpp:1033-1035`), and in Python. The driver is
`SRC/runtime/commands/analysis/modal/modal.cpp` (`modalProperties`, line 32;
`responseSpectrumAnalysis`, line 102).

### B.1 What `DomainModalProperties::compute()` produces
- **Total mass** per global DOF (translational + rotational), both full and
  free-DOF variants (`DomainModalProperties.cpp:555-567`).
- **Center of mass** from free-DOF masses (`.cpp:485-523`).
- **Generalized masses** `diag(VᵀMV)` (`.cpp:602`).
- **Modal participation factors** `(VᵀMR)/diag(VᵀMV)` (`.cpp:646`), R = rigid-body
  influence vector.
- **Effective modal masses** `(VᵀMR)²/diag(VᵀMV)` (`.cpp:649`).
- **Modal participation mass ratios** (cumulative % of free mass) (`.cpp:654-668`).
- **Parallel-axis rotational mass** augmentation about the COM (`.cpp:531-553`).

Requires `eigen` to have been run first — reads `domain->getEigenvalues()`
(`.cpp:298`, error check `.cpp:294-295`) and nodal eigenvectors
(`node->getEigenvectors()`, `.cpp:452`). Mass is read from `element->getMass()`
(`.cpp:439`) and `node->getMass()` (`.cpp:453`); a consistent mass is diagonalized
to a lumped form internally (`.cpp:354-415`). Restricted to ndm 2 or 3
(`.cpp:84-85`); excludes pressure DOFs (`.cpp:133-135, 267`).

### B.2 ResponseSpectrumAnalysis
- **Per-mode static superposition only** — NO built-in CQC/SRSS. Comment is
  explicit (`modal.cpp` header / ASDEA design). Per mode it sets nodal modal
  displacements `u = V·Γ·Sa(T)/λ` via `node->setTrialDisp(...)` and commits the
  domain once per mode (`ResponseSpectrumAnalysis.cpp:~242-277`); the **user must
  combine** the per-mode results externally (SRSS/CQC done in script/recorder).
- Acts along ONE global direction (`$dir`, 1..ndf); no arbitrary direction vector.
- Spectrum given either as a `timeSeries` or `-Tn/-fn` + `-Sa` lists with linear
  interpolation (`modal.cpp:241-322`).

### B.3 Limitations
- **Constraints/diaphragms**: mass on eliminated (slave, `id==-1`) DOFs is simply
  skipped (`DomainModalProperties.cpp:137-138, 359-369`); there is no transfer of
  constrained mass to masters. Rigid-diaphragm mass accounting is therefore only
  as good as the user's lumping; off-diagonal mass coupling is dropped in the
  rotational augmentation.
- **No damping** enters modalProperties (undamped eigenvectors only).
- **Serial-only**: no `processID`, channels, or MPI anywhere in either
  `DomainModalProperties.*` or `ResponseSpectrumAnalysis.*`. In a partitioned
  OpenSeesSP run these post-processors see only the rank-0 master domain, not the
  subdomain elements → wrong/partial mass and participation.

---

## (C) PARALLEL eigen reality — the definitive answer

### C.1 Is there a distributed EigenSOE / EigenSolver?
**No.** There is no PARPACK, SLEPc, ScaLAPACK, or distributed Lanczos in the tree.
The only Krylov driver is the **serial** ARPACK `dsaupd_`/`dseupd_`
(`ArpackSolver.cpp:120-130`, `BandArpackSolver.cpp:108-117`). grep for
`parpack`/`pdsaupd` returns nothing. The "parallel" eigen story is entirely:
*serial Lanczos on one rank + a distributed linear-solve / matvec underneath.*

### C.2 What actually happens when `eigen` runs in a partitioned domain

The machinery genuinely exists and is wired end-to-end:

1. `PartitionedDomain::eigenAnalysis` (`SRC/domain/domain/partitioned/PartitionedDomain.cpp:1179`)
   runs the base `Domain::eigenAnalysis` on the master (`:1207`) and then loops
   subdomains calling `theSub->eigenAnalysis(...)` (`:1216`).
2. `ShadowSubdomain::eigenAnalysis` (`SRC/domain/subdomain/ShadowSubdomain.cpp:1425`)
   and `ShadowSubdomain::setAnalysisEigenSOE` (`:1170`) ship the request/SOE to the
   remote actor.
3. `ActorSubdomain` handles `ShadowActorSubdomain_setAnalysisEigenSOE`
   (`SRC/domain/subdomain/ActorSubdomain.cpp:749`): it builds a **fresh** EigenSOE
   on the remote via `theBroker->getNewEigenSOE(theType)` (`:751`), `recvObject`s
   it (`:754`), and `setAnalysisEigenSOE`s it into the local analysis (`:755`);
   the eigen run is dispatched at `ActorSubdomain.cpp:139-150`.
4. **Note on serialization**: `ArpackSOE::sendSelf/recvSelf`
   (`ArpackSOE.cpp:269-372`) send only a `processID`/channel-handshake ID — **not**
   the matrices. Each subdomain assembles its own K/M locally into its own
   `ArpackSOE` + `LinearSOE`. The matrices stay distributed.

### C.3 Where the `(K−σM)⁻¹` solve happens
`ArpackSOE` wraps a `LinearSOE* theSOE` (`ArpackSOE.h:77`,
`ArpackSOE::setLinearSOE`, `ArpackSOE.cpp:382-386`). Inside the Lanczos loop the
shift-invert solve is `theSOE->solve()` (`ArpackSolver.cpp:258` for `ido=-1`,
`:276` for `ido=1`), and the `M·v` product + cross-rank reduction is done over the
ArpackSOE channels (`ArpackSolver::myMv`, `ArpackSolver.cpp:543-562`).

Crucially the active `LinearSOE` is propagated into the EigenSOE automatically:
`StaticAnalysis::setLinearSOE` and `setEigenSOE` both call
`theEigenSOE->setLinearSOE(*theSOE)` (`StaticAnalysis.cpp:549, 577`). So **if the
active `system` is a distributed `MumpsParallelSOE`, the shift-invert solve is a
genuine parallel MUMPS solve** (matrix A kept distributed, X/B on all ranks —
`MumpsParallelSOE.h:30-31`). The Lanczos coordinate iteration (`dsaupd_`) still
runs serially on the master; `ArpackSOE::checkSameInt` (`ArpackSOE.cpp:388-426`)
keeps the `ido` reverse-communication flag synchronized across ranks so all ranks
enter the distributed solve together.

### C.4 How the parallel build is wired
- **OpenSeesSP** = `_PARALLEL_PROCESSING` (PartitionedDomain / master-slave
  subdomains), Tcl main `SRC/tcl/mpiMain.cpp` (e.g. `mpiMain.cpp:184-291`,
  `OPS_PARALLEL_PROCESSING`/`OPS_NUM_SUBDOMAINS` set ~`:218,:229`).
- **OpenSeesMP** = `_PARALLEL_INTERPRETERS` (replicated interpreter, every rank
  runs the same script), main `SRC/tcl/mpiParameterMain.cpp:222-349`;
  `PythonMPIModule.cpp` is the MP-Python module.
- CMake gates: `OpenSeesCommands.cpp` / `OpenSeesMiscCommands.cpp` are the sources
  carrying the `_PARALLEL_*` guards (per root `CMakeLists.txt`).
- **Gotcha / mismatch**: `MumpsParallelSOE` is created only under a
  `_PARALLEL_INTERPRETERS` guard in `OPS_MumpsSolver`
  (`OpenSeesCommands.cpp:~4079-4094`, with `-ICNTL7`/`-ICNTL14`/`-matrixType`
  flags) — i.e. it's reachable from the **MP** (replicated-interpreter) Python/Tcl
  build, where it fetches MPI rank/channels and calls
  `setProcessID`/`setChannels`. The **SP** subdomain path (`getNewEigenSOE` +
  per-subdomain `ArpackSOE`) instead relies on each subdomain's own serial
  `LinearSOE`. So the two "parallel" stories don't compose cleanly: SP gives you
  distributed assembly but per-subdomain serial solves under a serial Lanczos; MP
  gives you a distributed MUMPS solve but every rank redundantly runs the whole
  (serial) Lanczos. Neither is a true distributed eigensolver.

### C.5 Does modalProperties / responseSpectrum work in parallel?
**Effectively no.** Both are serial (no MPI awareness, C.3 of section B). Under SP
they only see the master domain; under MP each rank would recompute on its own
replicated model (only meaningful if the whole model is replicated, which defeats
the purpose of partitioning). There is no parallel reduction of participation
factors / effective modal mass.

---

## (D) Buckling / geometric-stiffness eigen status

**Not implemented as an eigenproblem.** There is no `K·φ = λ·Kg·φ` path.

- Eigen assembly forms only `Kt` and mass `M` (`StaticAnalysis.cpp:285-323`,
  `EigenIntegrator.cpp:191-196`); `addM` always means mass (`addMtoTang`).
- The `eigen` parser comment lists "2 - buckling" for `generalizedAlgo`
  (`analysis.cpp:265-268`) but `generalizedAlgo` is a **bool** — buckling is
  **dead code / never wired**. Same in `OpenSeesCommands.cpp`.
- The element API *has* geometric-stiffness hooks but they are not used by eigen:
  `FE_Element::addKgToTang` (`SRC/analysis/fe_ele/FE_Element.cpp:423-436`) →
  `Element::getGeometricTangentStiff` which returns a **zero** matrix in the base
  class (only `PFEMElement2DCompressible` overrides). `addKgToTang` is used in
  PFEM sensitivity, never in any eigen algorithm.
- P-Δ / corotational transforms fold geometric effects directly into `Kt` at the
  current state (e.g. `PDeltaCrdTransf2d.cpp`), so the only route to "buckling"
  today is a nonlinear incremental analysis watching the tangent go singular — not
  a linearized buckling eigenvalue about a prestressed state.
- **No** `LadrunoBuckle` or any buckling command/class. `LadrunoRebarBuckling`
  (`SRC/material/uniaxial/LadrunoRebarBuckling.h`) is a *uniaxial material*
  (Dhakal-Maekawa σ-degradation), unrelated to structural buckling eigen.

---

## (E) Concrete gaps a design must close — ranked

1. **No distributed eigensolver (the headline gap).** ARPACK is serial
   (`dsaupd_`); for large partitioned models the Lanczos iteration is a serial
   bottleneck and the SP path doesn't even use a parallel linear solve per
   subdomain. *Close with:* PARPACK (drop-in `pdsaupd_`/`pdseupd_` mirroring the
   reverse-communication loop already in `ArpackSolver.cpp`), or SLEPc/ScaLAPACK.
   The reverse-communication structure + `checkSameInt` sync already in place make
   PARPACK the lowest-friction option.

2. **SP vs MP parallel-eigen incoherence.** `MumpsParallelSOE` is `_PARALLEL_INTERPRETERS`-gated
   (MP) while the subdomain `ArpackSOE` distribution is `_PARALLEL_PROCESSING` (SP).
   A design must pick one consistent model: distributed-assembly + distributed
   solve + distributed Krylov, all in the same build. Today you cannot get
   "partitioned domain + parallel `(K−σM)⁻¹` + parallel Lanczos" together.

3. **modalProperties / responseSpectrum are serial-only.** They read only the
   master domain and have zero MPI awareness — wrong results under SP. Needs a
   parallel reduction of `VᵀMV`, `VᵀMR`, total mass across ranks/subdomains.

4. **No linearized buckling eigenproblem (`K φ = λ Kg φ`).** Hooks exist
   (`addKgToTang`, `getGeometricTangentStiff`) but return zero and are never wired
   to eigen. A buckling design needs (a) elements to return a real Kg, (b) an
   eigen assembly that forms Kg into the "M" slot, (c) a command flag (revive the
   dead "buckling" option as a real mode), (d) a prestress state capture.

5. **No CQC/SRSS modal combination built in.** `responseSpectrumAnalysis` is
   per-mode static superposition only; combination is left to the user. A complete
   modal-analysis offering should provide CQC (needs modal damping ratios) and
   SRSS natively.

6. **Shift/sigma not exposed in the xara `eigen` command.** `shift` is hard-zero
   (`analysis.cpp:271`); the SOE supports it but there's no `-shift`/`-sigma`
   flag, hurting convergence for clustered/high modes and for stiff-but-soft modes.

7. **Constraint / rigid-diaphragm mass handling in modalProperties is naive.**
   Eliminated-DOF mass is dropped, off-diagonal mass coupling ignored. For
   diaphragm-heavy building models the effective modal mass can be wrong.

8. **`SymmGeneralizedEigenSOE` disabled on Windows** (`OpenSeesCommands.cpp:340-346`)
   — a portability gap if that solver is wanted as a serial dense fallback.

---

### Appendix — file map

- Serial eigen SOEs/solvers: `SRC/system_of_eqn/eigenSOE/` (Arpack*, BandArpack*,
  SymArpack*, SymBand*, FullGen*, SymmGeneralized*).
- Eigen analysis/integrator/algorithms: `SRC/analysis/analysis/EigenAnalysis.*`,
  `SRC/analysis/integrator/EigenIntegrator.*`,
  `SRC/analysis/algorithm/eigenAlgo/{FrequencyAlgo,StandardEigenAlgo}.*`.
- `eigen` parsers: `SRC/runtime/commands/analysis/analysis.cpp:250`,
  `SRC/interpreter/OpenSeesCommands.cpp:270`.
- Modal post-processing: `SRC/runtime/commands/analysis/modal/`
  (`DomainModalProperties.*`, `ResponseSpectrumAnalysis.*`, `modal.cpp`); legacy
  sibling `SRC/domain/domain/DomainModalProperties.*`.
- Parallel: `SRC/domain/domain/partitioned/PartitionedDomain.cpp:1179`,
  `SRC/domain/subdomain/{ShadowSubdomain,ActorSubdomain,Subdomain}.cpp`,
  `SRC/system_of_eqn/linearSOE/mumps/MumpsParallel{SOE,Solver}.*`,
  `SRC/tcl/mpiMain.cpp` (SP), `SRC/tcl/mpiParameterMain.cpp` (MP).
- Ladruno modal-relevant: `SRC/analysis/integrator/CentralDifferenceLadruno.*`
  (explicit; `getVel()` modal-damping hook ~`:539`),
  `SRC/analysis/integrator/CriticalTimeStep.*` (per-element generalized eigensolve
  `dsygvx_`/`dggev` for dt_cr — stability only, not modal output). The fictitious-
  mass family (`LadrunoDynamicRelaxation`, `LadrunoArcLength`,
  `LadrunoFictitiousMass.h`) is quasi-static, not modal.
