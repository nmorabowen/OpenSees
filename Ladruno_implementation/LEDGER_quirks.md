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

### `Truss` (and most elements) default to a LUMPED mass — `-lump diagonal` ≡ `-lump rowsum` on them
- **Bites:** the ExplicitBathe / CriticalTimeStep `-lump <rowsum|diagonal>` option
  only changes `dt_cr` when the element's `getMass()` returns a **consistent**
  (non-diagonal) matrix. The `Truss` element defaults to `cMass=0` (lumped:
  `0.5*rho*L` on the diagonal, zero off-diagonals). For a diagonal matrix,
  "diagonal-of-consistent" (take `M_ii`) and "row-sum" (sum each row onto the
  diagonal) are identical, so both lumpings return the same `dt_cr` (= `le/c` for
  a 2-node bar). A test that asserts `dt_diag < dt_rowsum` on a default lumped bar
  is **vacuous and fails** (the two tie).
- **Why:** `Truss::getMass()` returns lumped mass unless built with `-cMass 1`;
  the consistent branch is `rho*A*Le/6 * [[2,1],[1,2]]`. Only then do the lumpings
  diverge: diagonal-of-consistent = `rho*A*Le/3`, row-sum = `rho*A*Le/2`, giving
  `dt_diag/dt_rowsum = sqrt((Le/3)/(Le/2)) = sqrt(2/3) ≈ 0.816` (per-element pencil).
- **How it surfaced:** `tests/test_explicitBathe_integrator.py::test_eb_lump_diagonal`
  failed on Zone-A (Ubuntu) CI with `diagonal=0.02000 vs rowsum=0.02000` (a tie at
  `le/c`). The Linux result is **correct**; the source `_verify_explicit.py` test 11
  had the same faulty premise (built a lumped bar yet asserted the inequality).
- **Workaround/status:** build the bar with `-cMass 1` (consistent mass) for any
  test that means to exercise the diagonal-vs-rowsum *difference*. Fixed in the
  Zone-A port (PR #40) and mirrored back into `Ladruno_scripts/_verify_explicit.py`.
  *(Learned 2026-05-31, Zone-A ExplicitBathe battery.)*

### MPCORecorder `exit(-1)` kills the kernel
- **Bites:** ~25 raw libc `exit(-1)` calls in upstream `MPCORecorder.cpp` hard-kill
  the Jupyter kernel on any recorder error and leave the `.mpco` unflushed.
- **Why:** upstream uses `exit()` for error handling; no exception path. The file
  has been byte-identical to `master` and frozen ~4.5 years.
- **Status:** diagnose-only as of 2026-05-18. A local `ladruno` patch is the only
  fix; deferred. The Ladruno recorder rewrite (`ladruno`) avoids the pattern.

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
- **Bites:** `ladruno` emitted ~3 non-fatal `HDF5-DIAG` errors to stderr for
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
  `QUADRATURE`/{`GP_PARAM`,`GP_WEIGHT`}. See [[LEDGER_implementations]] Ladruno row.

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
- **Status:** for `Ladruno`, the standard-quadrature `GP_PARAM[k]` MUST equal the
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
(`OPS_partition` returns an error; needs METIS 5 / `OPS_HAVE_METIS5`), so the
broadcast path is not runtime-testable in this build.

### Ladruno recorder nested result-group names need the `<display>` parent pre-created
- **Bites:** writing an element result whose `schema.name` is nested
  (`"<display>/<bucket>"`, e.g. `stress/204-FourNodeQuad[201:0:0]`) under a parent
  that does NOT already contain the intermediate `<display>` group → the
  `H5Gcreate` in `createResultGroup` returns an invalid handle and the datasets
  underneath silently fail (HDF5-DIAG noise / missing data).
- **Why:** the recorder's group-creation property list (`h_group_proplist`) is a
  GCPL with `CRT_ORDER_TRACKED|INDEXED` but carries **no create-intermediate-group
  LINK property**, and `createResultGroup` passes `lcpl=H5P_DEFAULT`. `H5Gcreate`
  only creates the *final* path component; any intermediate must already exist.
- **Status/workaround:** the time-series `StreamingSink` path is safe because the
  recorder pre-creates `ON_ELEMENTS/<display>` at init. The `EnvelopeSink` is
  self-contained, so it pre-creates each prefix segment of `m_name` inside
  `writeEnvelope` (the per-flush `H5Ldelete` only deletes the leaf link, so the
  intermediate persists). This is why element envelopes were latently broken until
  PR #45 — only flat node/domain names had ever been written. Learned 2026-05-31.

### `decode_columns` in `ladruno_format.py` must flatten COLUMN_MAP arrays — the recorder writes them 2-D `[k×1]`
- **Bites:** `int(mult[i])` (and the other per-block scalars) in `decode_columns`
  throws `TypeError: only 0-dimensional arrays can be converted to Python scalars`
  under numpy 2.x when reading a **real recorder** `.ladruno`.
- **Why:** the C++ writer stores each per-block COLUMN_MAP array via
  `createAndWrite(vec, k, 1)` → **2-D `[k×1]`**, so `arr[i]` is a `(1,)`-array, not
  a scalar. `make_synthetic.py` writes them 1-D `[k]`, which masked it; and the
  element-**parity** gate keys results by flat column index and never calls
  `decode_columns`, so no test exercised it on real 2-D output until PR #45's
  element-envelope checker.
- **Status:** fixed (PR #45) — `decode_columns` `.reshape(-1)`s GAUSS_ID/SECTION_TAG/
  FIBER_ID/NUM_COMP/MULTIPLICITY (LEVELS is consumed via `np.atleast_1d`, left as-is).
  Learned 2026-05-31.

### h5py reads of a freshly-written `.mpco`/`.ladruno` HANG on HDF5 file locking
- **Bites:** a venv-python checker calling `h5py.File(path, "r")` on a file the
  build-python just wrote (in the same gate run) **hangs indefinitely** — no error,
  just blocks (e.g. `parity_check.py` stuck forever).
- **Why:** HDF5's default file locking; the writer's lock/superblock state isn't
  cleared promptly on a synced/Temp FS, so the reader blocks acquiring the lock.
- **Workaround:** set `HDF5_USE_FILE_LOCKING=FALSE` in the checker's environment
  (`$env:HDF5_USE_FILE_LOCKING="FALSE"`) → opens instantly (parity 80/80). Apply to
  ALL `.ladruno`/`.mpco` read steps after a recorder run. Learned 2026-05-31.

### Parallel build OOM (`cl.exe` C1060 "out of heap") on the giant template TUs under RAM pressure
- **Bites:** with low free RAM (~1–2 GB of 28), `cmake --build ... -j8`/`-j16` dies
  with `fatal error C1060: compiler is out of heap space` on the huge template TUs
  (`OPS_AllASDPlasticMaterial3Ds.cpp`, `MPCORecorder.cpp`, `LadrunoRecorder.cpp`),
  and even ordinary TUs get OS-killed (ninja `FAILED: [code=2]` with no compiler
  diagnostic = the process was killed, not a code error).
- **Workaround:** compile the monsters **serially first** —
  `ninja -j1 CMakeFiles/OPS_Material.dir/SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp.obj`
  (and the two MPCO recorder objs) — then `cmake --build build\build\Release
  --target OpenSeesPy -j2` for the rest (ninja resumes cached objs). Don't assume a
  `code=2` with no error text is a code bug; check free RAM first. Learned 2026-05-31.

### `getCommitTag()` is a GLOBAL monotonic counter — it does NOT reset on `wipe()`
- **Bites:** any per-step recorder/series that uses `Domain::getCommitTag()` as its
  step axis (the analysis monitor's `STEP`, the profiler per-step series). Across
  several `analyze()` runs in one interpreter session — even with `wipe()` between
  them — the commitTag keeps climbing (run 1 → 0..199, run 2 → 200..399, ...). It is
  NOT a within-analysis 0-based step index.
- **Implication:** don't assert absolute step values or compare step arrays across
  runs by value; compare the *stride* (`np.diff(step) == every`) instead. For a live
  viewer, treat `STEP` as a monotonic id, not "step N of this analysis."
- Learned 2026-05-31 building the analysis monitor (`08_analysis_monitor.md`).

### `recorder ladruno` is NOT wired into the classic Tcl `OpenSees.exe` (was; now fixed)
- **Bites:** `recorder ladruno ...` works from OpenSeesPy/openseesmp and the
  interpreter-based Tcl (`TclWrapper`→`OPS_Recorder`, the shared map in
  `OpenSeesOutputCommands.cpp`), but the **classic** Tcl `OpenSees.exe`
  (`commands.cpp`→`addRecorder`→`TclAddRecorder` in `TclRecorderCommands.cpp`)
  hardcodes its recorder dispatch in a *separate* file that the rename PR never
  touched — it had `mpco`/`vtkhdf`/`gmsh`/`EnergyBalance` but no `ladruno`. So
  `recorder ladruno` raised "recorder type ladruno is unknown" only in `OpenSees.exe`.
- **Fix:** added the `else if (strcmp(argv[1],"ladruno")==0)` branch (+ extern
  `OPS_LadrunoRecorder`) mirroring the `mpco` block in `TclRecorderCommands.cpp`.
  Lesson: there are TWO recorder-dispatch tables (shared `OPS_Recorder` map for
  Py/interpreter-Tcl; hardcoded `TclAddRecorder` for classic Tcl) — wire new
  recorders into BOTH. Learned 2026-05-31.

### Ladruno recorder `modesOfVibration` (eigen) output writes no data — modal `DATA` group/dataset collision
- **Bites:** `recorder ladruno -N modesOfVibration` after `ops.eigen(n)` creates the
  `MODEL_STAGE[*]/RESULTS/ON_NODES/MODES_OF_VIBRATION(U)` group with a valid schema
  (ID, COMPONENTS) but `DATA` stays empty `(0, nNodes, nComp)` and no `MODE_k`
  datasets appear; HDF5-DIAG "can't synchronously write data / Write failed" errors
  fire. Happens for ANY model (reproduced with a bare elasticBeamColumn portal —
  NOT the known fiber-section `writeSections` noise), under both `ops.record()` and
  an `analyze()` step.
- **Why:** `LadrunoRecorder::recordModeChannel` (LadrunoRecorder.cpp ~1692) calls
  `ch.sink->begin()`, which (StreamingSink) creates `.../MODES_OF_VIBRATION(U)/DATA`
  as a **chunked dataset** (normal time-series layout), then tries to
  `h5::group::create` a **group** at `DATA/STEP_<step>` with `MODE_k` datasets
  under it (the MPCO modal layout). A group cannot be created beneath an existing
  dataset → the HDF5 calls fail, no modal data is written.
- **Status:** **FIXED 2026-05-31.** `recordModeChannel` no longer calls the StreamingSink
  `begin()`. It now owns the modal init, mirroring frozen `ResultRecorderModesOfVibration::record`:
  once per stage (idempotent via `H5Lexists`) it creates the result group
  (`h5::group::createResultGroup`) + `ID` dataset + `DATA` **group**, then per step writes
  `DATA/STEP_<step>/MODE_<k>` datasets with MODE/LAMBDA/OMEGA/FREQUENCY/PERIOD attrs. The
  validator `ladruno_format.py::_check_data_shape` was taught the modal layout (DATA = group
  of STEP_<step> groups of MODE_<k> datasets). **Modal eigenvectors now match frozen mpco to
  1e-12**; the EIGEN gate is promoted into the counted regression battery and `eigen_check.py`
  does a real modal value-parity diff vs the ref `.mpco`. Files: `SRC/recorder/LadrunoRecorder.cpp`,
  `ladruno_format.py`, `eigen_check.py`, `run_regression.bat`. Found + fixed 2026-05-31 by the
  new eigen coverage gate (the test scheme catching, then confirming the fix of, a real bug).

### Energy-balance check script `energy_check.py` was stale vs the chunked `DATA` layout (fixed)
- **Bites:** `energy_check.py` failed with `TypeError: Only 1D arrays allowed for
  fancy indexing` — its `_read_result` iterated `grp["DATA"]` as a per-step *group*,
  but the recorder writes `ON_DOMAIN/ON_REGIONS energyBalance DATA` as a chunked
  `[T×nrows×ncomp]` **dataset** (the standardized streaming layout; recorder output
  is correct).
- **Fix:** read via `lf.iter_step_slices(grp)` (the canonical slicer the parity
  checks already use; tolerates both chunked and legacy `DATA/STEP_k` layouts).
  After the fix the energy kernel matches the EnergyBalance text sidecar to ~5e-9.
  Learned 2026-05-31.
