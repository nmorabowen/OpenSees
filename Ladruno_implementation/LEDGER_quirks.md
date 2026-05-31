---
title: Ledger — OpenSees quirks & gotchas
project: Ladruno
tags:
  - ledger
  - quirks
  - gotchas
---

# Ledger — OpenSees quirks & gotchas we learned

Surprising, undocumented, or bug-prone behaviours of upstream OpenSees that
cost us time. Recording them here so we (and future us) stop re-discovering
them. This is observation-only — fixes we actually applied are tracked in
[[LEDGER_vanilla_files]] / [[LEDGER_implementations]].

## Conventions

- **One section per quirk.** Title = the symptom you'd search for.
- State: *what bites*, *why*, *workaround/status*, and the *date* learned.
- If a quirk drove a code change, cross-link the ledger row / PR.
- Deep build/toolchain quirks may live in
  [[../Ladruno_internal/01_compilation_journal]]; link rather than duplicate.

## Quirks

### MPCORecorder `exit(-1)` kills the kernel
- **Bites:** ~25 raw libc `exit(-1)` calls in upstream `MPCORecorder.cpp` hard-kill
  the Jupyter kernel on any recorder error and leave the `.mpco` unflushed.
- **Why:** upstream uses `exit()` for error handling; no exception path. The file
  has been byte-identical to `master` and frozen ~4.5 years.
- **Status:** diagnose-only as of 2026-05-18. A local `ladruno` patch is the only
  fix; deferred. The `MPCO_Ladruno` rewrite (`mpcol`) avoids the pattern.

### `patch_banner.py` path went stale after workspace consolidation
- **Bites:** `patch_banner.py` computed `TARGET = ROOT/OpenSees/SRC/...`, assuming an
  inner `OpenSees/` dir. After the 2026-05-30 consolidation `SRC/` sits directly
  under the repo root, so the script silently targeted a non-existent path.
- **Why:** the script predated moving the workspace into the repo root.
- **Status:** **fixed** — paths recomputed to `ROOT/SRC/...`; the script now also
  patches the feature list into both `tclMain.cpp` and `PythonModule.cpp`.

### MPCO is nodal-blind without an element-result parity gate
- **Bites:** a recorder can pass all *nodal* parity checks yet silently drop or
  mis-tag *element/material* results (e.g. Lagrange quad/tri `setResponse` tags).
- **Why:** nodal and element result paths are independent; testing one doesn't
  cover the other.
- **Status:** addressed by the element-result parity gate ([#14](https://github.com/nmorabowen/OpenSees/pull/14)) and the upstream `setResponse` fix ([#7](https://github.com/nmorabowen/OpenSees/pull/7)).

### Damping shrinks the explicit critical time step — and βK is a trap
- **Bites:** adding damping *reduces* `dt_cr` in both coupled and explicit modes;
  stiffness-proportional damping (`βK`) is especially punishing under explicit
  integration.
- **Why:** βK inflates the effective high-frequency content that sets `dt_cr`.
- **Status:** use mass-proportional (`αM`) damping in explicit; captured in the
  CentralDifferenceLadruno plan (`project_robust_central_difference`).

### `partition` needs the METIS 5 API (Patch 9)
- **Bites:** the openseespy `partition` command fails / mislinks without the
  METIS 5 API path.
- **Why:** upstream assumes an older METIS; the MP build links METIS 5.
- **Status:** guarded with a clear error message in `OpenSeesMiscCommands.cpp`
  ([#1](https://github.com/nmorabowen/OpenSees/pull/1)); `OPS_HAVE_METIS5` follow-up noted.

### HDF5 won't create a group under a dataset — the "writeSections" red herring
- **Bites:** `mpcoLadruno` emitted ~3 non-fatal `HDF5-DIAG` errors to stderr for
  every custom-rule element (Lobatto beam, MVLEM, BezierTri6). The error
  signature (`invalid location identifier` / `dset_id is not a dataset ID`) was
  misattributed to `writeSections` for a while.
- **Why:** `writeModelElements` made `MODEL/ELEMENTS/<name>` the CONNECTIVITY
  **dataset**, then tried to `H5Gcreate` a `QUADRATURE` child **group** under it.
  HDF5 datasets cannot parent groups, so the handle was invalid and the
  `GP_PARAM`/`GP_WEIGHT` writes under it failed — silently losing the quadrature
  group (data was otherwise intact; GP_X attr carried the coords). The opserr
  buffering made the errors *look* like they came from `writeSections`, which is
  clean. Reliably-ordered `fprintf(stderr)` pinned them to `writeModelElements`.
- **Status:** **fixed.** [#16](https://github.com/nmorabowen/OpenSees/pull/16)
  removed the broken block as a stopgap (GP_X-only); [#18](https://github.com/nmorabowen/OpenSees/pull/18)
  did the real schema-v1 fix — `<name>` is now a GROUP holding `CONNECTIVITY` +
  `QUADRATURE`/{`GP_PARAM`,`GP_WEIGHT`}. See [[LEDGER_implementations]] MPCO_Ladruno row.

### `Response`/`Information` stores a Matrix and a Vector separately — `getData()` only returns the Vector
- **Bites:** an element `setResponse` that returns a `Matrix` (via `setMatrix`,
  e.g. BezierTri6 `integrationPoints` → `Matrix(nGP,2)`) is invisible to code
  reading `response->getInformation().getData()` (that returns the **Vector**
  slot, which is empty/unrelated). The recorder's legacy `getCustomGaussPointLocations`
  read the Vector and silently missed the matrix-typed barycentric GP coords.
- **Why:** `Information` is a tagged union (`theType` ∈ {…, VectorType, MatrixType}
  with separate `theVector`/`theMatrix` pointers); `getData()` hardcodes the Vector.
- **Status:** **handled** — check `info.theType == MatrixType` and read `*info.theMatrix`
  for multi-dim parametric rules ([#18](https://github.com/nmorabowen/OpenSees/pull/18)).

### Gauss-point ordering is per-element, NOT a standard tensor order
- **Bites:** there is no global "Gauss point #k → natural coords" convention across
  OpenSees solid elements. Deriving GP coordinates from a `(geometry, rule)` table and
  assuming the usual lexicographic tensor order silently mis-maps results (stress at
  GP `k` paired with the wrong location) — invisible to nodal parity and to any
  value-only check (the frozen MPCO recorder never writes GP coords; STKO holds a
  hardcoded per-rule table reader-side).
- **Why:** each element hardcodes its own integration-point loop. Verified from sources:
  `FourNodeQuad` (4-pt 2×2) walks **counter-clockwise** `(−,−),(+,−),(+,+),(−,+)`
  (FourNodeQuad.cpp:298-305) — *not* the lexicographic `(−,−),(+,−),(−,+),(+,+)`;
  `Brick`/stdBrick (8-pt) walks nested `for i{for j{for k}}` ⇒ x-outer..z-inner
  **lexicographic** (Brick.cpp:536-542); `NineNodeQuad` (9-pt) walks 4 CCW corners →
  4 CCW edge-mids → center (NineNodeQuad.cpp:127-144), a serendipity-style order, not
  a 3×3 tensor sweep. So even two "quad GL" elements need not share GP order.
- **Status:** for `MPCO_Ladruno`, the standard-quadrature `GP_PARAM[k]` MUST equal the
  element's own k-th GP natural coords (so it pairs with result `gauss_id k`); the
  table is verified against each canonical element's source, and the belt-and-suspenders
  `GLOBAL_GP_COORDS` round-trip oracle (`x(GP_PARAM[k])` vs the C++-computed global GP)
  catches any ordering/basis mismatch at write time. Learned 2026-05-30 building the
  standard-rule QUADRATURE table.

### `FourNodeTetrahedron` has 1 Gauss point, not 4 (the recorder comment lies)
- **Bites:** the recorder's `getGeometryAndIntRuleByClassTag` comment reads
  "4-node tetrahedron with **4 gp**" and maps it to `Tetrahedron_GaussLegendre_1`; an
  implementer trusting that would write 4 GP coordinates for a 1-GP element.
- **Why:** `FourNodeTetrahedron.cpp` defines `sg[]={0.25}` (a single abscissa) and its
  Gauss loop is collapsed: `i = j = k = 0; // Just one Gauss point in a tet`
  (FourNodeTetrahedron.cpp:226,574). A linear tet integrates exactly with one centroid
  point `(¼,¼,¼)` in barycentric coords; the "4 gp" comment is stale/wrong.
- **Status:** table uses **1 GP** for the tet rule. Authoritative GP count comes from
  the element source / the OutputDescriptor response walk, never the recorder comment.
  Learned 2026-05-30.

### `openseessp` is an unbuilt subsystem (Python has no SP)
- **Bites:** expecting an `import openseessp` analogous to `openseesmp`.
- **Why:** the SP parallel engine exists only on the Tcl side; the Python engine
  never had an SP build. Its build trace was fully removed from the fork.
- **Status:** by design — use `openseesmp`. Rationale in `02_openseespymp.md`.

## Quirk: bundled `OTHER/LAPACK` ships `dsygvx` but NOT `dsygv` (Linux link break)

The fork's in-tree reference LAPACK (`OTHER/LAPACK/`, statically linked by the
Unix/Ubuntu Zone-A build) is a *curated subset* of LAPACK, not the full library.
It provides the **expert/range** symmetric-definite generalized eigensolver
`dsygvx.f` (line ~147 of `OTHER/LAPACK/CMakeLists.txt`) but does **not** provide
the simple driver `dsygv.f`. It also ships `dsyevx` but not `dsyev`, `dggev` but
the `x`-suffixed expert drivers are the ones actually bundled for the symmetric
path.

Consequence: fork code that calls `dsygv_` compiles per-TU fine and links fine
on **Windows** (full MKL resolves `DSYGV`), but fails the **Linux** link with
`undefined reference to 'dsygv_'`. This silently blocked the entire ladruno
Zone-A CI Build step (and thus every downstream runtime test battery) after
`CriticalTimeStep.cpp` (CentralDifferenceLadruno, PR #22) introduced a `dsygv_`
call — even though the recorder/MPCO work was unrelated.

**Rule for fork code:** for a symmetric-definite generalized eigenproblem, call
the bundled **`dsygvx_`** (expert driver), not `dsygv_`. To get just the largest
eigenvalue use `RANGE='I'` with `IL=IU=N`; the value lands in `W[0]`. Mirror the
existing, proven usage in
`SRC/system_of_eqn/eigenSOE/SymmGeneralizedEigenSolver.cpp`. Adding `dsygv.f` to
the bundle instead is the wrong fix: it cascades (dsygv calls `dsyev`, also not
bundled) and edits the upstream LAPACK source list.
*(Fixed in CriticalTimeStep.cpp via `fix/criticaltimestep-dsygvx-link`.)*

## OPS_Recorder is compiled sequentially (no `_PARALLEL_PROCESSING`) for ALL targets

`SRC/recorder/MPCORecorderLadruno.cpp` (and `MPCORecorder.cpp`) live in the shared
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
(`OPS_partition` returns an error; needs METIS 5 / `OPS_HAVE_METIS5`), so the
broadcast path is not runtime-testable in this build.
