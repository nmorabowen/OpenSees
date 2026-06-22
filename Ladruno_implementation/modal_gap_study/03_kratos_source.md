# Kratos Multiphysics — Eigen / Buckling / Harmonic / Parallel-Eigen source dossier

Portable open-source template study for a new OpenSees modal/buckling/FRF feature.

Scope: how Kratos `StructuralMechanicsApplication` + `LinearSolversApplication`
implement (1) modal/eigen analysis, (2) linear buckling (prebuckling), (3)
harmonic / frequency-response, (4) complex / damped modal, (5) parallel
eigensolving. Each section ends with an **OpenSees porting note**.

## Sources read

Local skill refs (full):
- `C:\Users\nmora\.claude\skills\kratos\SKILL.md`
- `C:\Users\nmora\.claude\skills\kratos\references\structural-mechanics.md`
- `...\references\core-framework.md`
- `...\references\developer-guide.md`

GitHub source / docs (fetched from `raw.githubusercontent.com/KratosMultiphysics/Kratos/master/...`,
and `gh api search/code` for path discovery):
- `applications/StructuralMechanicsApplication/custom_strategies/custom_strategies/eigensolver_strategy.hpp`
- `.../custom_strategies/custom_strategies/prebuckling_strategy.hpp`
- `.../custom_strategies/custom_strategies/direct_harmonic_analysis_strategy.hpp`
- `.../custom_strategies/custom_strategies/modal_harmonic_analysis_strategy.hpp`
- `applications/StructuralMechanicsApplication/python_scripts/structural_mechanics_eigensolver.py`
- `.../python_scripts/structural_mechanics_prebuckling_solver.py`
- `applications/LinearSolversApplication/README.md`
- `applications/LinearSolversApplication/custom_solvers/eigensystem_solver.h`
- `applications/LinearSolversApplication/custom_solvers/feast_eigensystem_solver.h`
- `applications/LinearSolversApplication/external_libraries/FEAST/4.0/include/pfeast*.h` (existence/parallel confirmation via `gh api`)
- `docs/pages/Kratos/For_Developers/Solvers/FEAST_Solver.md`

Key negative finding (verified by `gh api search/code q='Anasazi repo:KratosMultiphysics/Kratos'`
→ **zero hits**, and `q='eigen path:applications/TrilinosApplication'` → **zero hits**):
**Kratos has NO Trilinos/Anasazi/SLEPc eigensolver.** Its only distributed-eigen
path is the bundled **PFEAST** (MPI FEAST). Parallelism elsewhere is shared-memory
OpenMP. This is itself an important data point for the OpenSees decision.

---

## 0. The architectural pattern (shared by all four analyses)

Kratos cleanly separates the modal/buckling/FRF analyses into **three reusable
layers**, all driven from JSON with no recompile:

1. **Scheme** (`Scheme<...>`) — decides *what each element contributes* by setting
   a `BUILD_LEVEL` flag in `ProcessInfo`. `EigensolverDynamicScheme` makes elements
   emit either their mass matrix (`BUILD_LEVEL==1`) or stiffness matrix
   (`BUILD_LEVEL==2`) through the *same* `CalculateLocalSystem` entry point.
2. **BuilderAndSolver** — assembles the global sparse matrix for whichever
   `BUILD_LEVEL` is active, applies constraints (`ApplyConstraints`) and Dirichlet
   BCs (`ApplyDirichletConditions`) **without changing the DOF/sparsity pattern**
   (critical: K and M must share an identical pattern for the generalized solve).
3. **Strategy** (`SolvingStrategy`) — orchestrates: build M, build K, hand the two
   sparse matrices to a *generalized eigensolver* (which is just the strategy's
   `LinearSolver` plug-in), then write results back to `EIGENVALUE_VECTOR`
   (ProcessInfo) and `EIGENVECTOR_MATRIX` (nodal, dense per-node block).

The generalized eigensolver is **not** a bespoke strategy method — it is a
`LinearSolver` subclass whose `Solve(K, M, values, vectors)` overload takes *two*
matrices. This is the single cleanest idea to lift: the eigensolver is a swappable
backend behind one interface.

**OpenSees porting note.** OpenSees already mirrors layers 1–3:
`EigenAnalysis`/`EigenIntegrator`/`EigenSOE`/`EigenSolver` (`SRC/analysis/...`,
`SRC/system_of_eqn/eigenSOE/`). The Kratos `BUILD_LEVEL` trick maps onto
`Integrator::formK()` vs `formM()` already present in `EigenIntegrator`. The lift
is mostly conceptual: keep M and K on an identical sparsity, and treat the
eigensolver as a pluggable backend (OpenSees has `ArpackSOE`/`ArpackSolver`,
`SymBandEigenSolver`, `FullGenEigenSolver` — the slot exists).

---

## 1. EIGENSOLVER STRATEGY (`eigensolver_strategy.hpp`)

### Assembly of M and K
In `SolveSolutionStep()`:

```cpp
// MASS
rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
TSparseSpace::SetToZero(rMassMatrix);
pBuilderAndSolver->Build(pScheme, rModelPart, rMassMatrix, b);   // b discarded

// STIFFNESS
rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
TSparseSpace::SetToZero(rStiffnessMatrix);
pBuilderAndSolver->Build(pScheme, rModelPart, rStiffnessMatrix, b);
```

Both go through the *same* `Build()`; the scheme (`EigensolverDynamicScheme`)
reads `BUILD_LEVEL` and asks each element for `CalculateMassMatrix` vs
`CalculateLeftHandSide`/local stiffness. Constraints and Dirichlet are applied to
*both* identically so the patterns match.

### Generalized eigen solve (K φ = λ M φ)
```cpp
pBuilderAndSolver->GetLinearSystemSolver()->Solve(
    rStiffnessMatrix,   // K
    rMassMatrix,        // M
    Eigenvalues,        // λ  (= ω² for undamped structural)
    Eigenvectors);      // Φ
```
The actual algorithm lives entirely in the plugged-in `LinearSolver` (§5). The
strategy is backend-agnostic.

### Result storage
- `rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = Eigenvalues;`
- per-node dense block `node.GetValue(EIGENVECTOR_MATRIX)` = (n_modes × n_node_dofs)
  slab of Φ. (This is how the harmonic solvers later re-read the modes — §3.)

### Number of modes / shift
**Not** in the strategy. `number_of_eigenvalues` and `shift` are parameters of the
eigensolver backend (`eigensolver_settings`, §5). Defaults (Python layer): solver
`"spectra_sym_g_eigs_shift"`, `number_of_eigenvalues: 5`, `scheme_type: "dynamic"`.

### Mass normalization
Optional `MassNormalizeEigenvectors()` (flag
`normalize_eigenvectors_with_mass_matrix`): for each mode φⱼ compute
`mⱼ = φⱼᵀ M φⱼ`, scale `φⱼ ← φⱼ / √mⱼ` so `φᵀ M φ = I`. Also an optional
`compute_modal_decomposition` flag to emit modal participation data.

### Python wiring (`structural_mechanics_eigensolver.py`)
- `_CreateScheme()` → `StructuralMechanicsApplication.EigensolverDynamicScheme()`
- `_CreateLinearSolver()` → `eigen_solver_factory.ConstructSolver(eigensolver_settings)`
- `_CreateSolutionStrategy()` → `EigensolverStrategy(...,
  mass_matrix_diagonal_value, stiffness_matrix_diagonal_value,
  compute_modal_decomposition, normalize_eigenvectors_with_mass_matrix)`
- `mass_matrix_diagonal_value` / `stiffness_matrix_diagonal_value`: tiny diagonal
  values injected on *constrained* DOFs so the rows aren't structurally singular —
  a practical robustness trick.

**OpenSees porting note.** OpenSees `eigen N` already does K φ = λ M φ via Arpack
(`ArpackSOE`/`ArpackSolver`) and stores modes on nodes (`Node::getEigenvectors`).
The transferable pieces: (a) the explicit `normalize_eigenvectors_with_mass_matrix`
post-pass (OpenSees normalizes differently per solver — making it an explicit,
documented option is cleaner); (b) the `*_matrix_diagonal_value` regularization on
constrained DOFs to avoid singular eigen rows; (c) emitting modes in a stable
per-node block that downstream FRF can re-read — OpenSees lacks a built-in modal-FRF
consumer, so this is the hook to add.

---

## 2. PREBUCKLING STRATEGY (`prebuckling_strategy.hpp`) — the clean buckling reference

This is the most valuable single file for the OpenSees buckling gap, because
OpenSees has no first-class linear-buckling procedure. Kratos implements **Jia &
Mang's** consistent-linearization buckling, *not* the textbook "two-snapshot
secant" form most codes use.

### Two stiffness snapshots, differenced
The strategy holds:
- `mpStiffnessMatrix` — current load step's tangent stiffness
- `mpStiffnessMatrixPrevious` — previous step's tangent stiffness

Each load increment is driven by a genuine **nonlinear static solve**:
```cpp
pBuilderAndSolver->BuildAndSolve(pScheme, rModelPart, rStiffnessMatrix, rDx, rRHS);
```
so geometric/stress stiffening accumulates through the element's own geometric
nonlinearity (no separate "Kg element routine" needed — the tangent already
contains K0+Kg at the current state).

After convergence of a step, it stores the previous tangent and forms the
**stiffness difference** (the discrete derivative of K w.r.t. load):
```cpp
if (mLoadStepIteration % 2 == 0)
    rStiffnessMatrixPrevious = rStiffnessMatrix;        // snapshot reference
...
rStiffnessMatrix = rStiffnessMatrixPrevious - rStiffnessMatrix;   // ΔK ≈ -dK/dλ · Δλ
```

### Eigenproblem actually solved
```cpp
pGetEigenSolver()->GetLinearSystemSolver()->Solve(
    rStiffnessMatrixPrevious,   // reference tangent  K_ref  (≈ K0 at the step)
    rStiffnessMatrix,           // ΔK  (stiffness change between steps)
    Eigenvalues, Eigenvectors);
```
i.e. **K_ref φ = λ ΔK φ.** The eigenvalue λ is a multiplier on the *incremental*
stiffness change, which is the consistent linearization of the buckling condition
det(K(λ)) = 0 near the current load. This avoids the inaccuracy of the naive
`(K0 + λ Kg) φ = 0` form when Kg is taken from a single far preload.

### Critical-load extraction & incremental driver
```cpp
mLambda = Eigenvalues(0) * delta_load_multiplier;          // smallest eig → critical
for (i ...) Eigenvalues[i] = mLambdaPrev + Eigenvalues[i]*delta_load_multiplier;
AssignVariables(Eigenvalues, Eigenvectors);                // → EIGENVALUE_VECTOR / EIGENVECTOR_MATRIX
```
The outer loop is an **eigenvalue-tracking continuation**:
- iteration 1: apply `mInitialLoadIncrement·(λ+λprev)`
- odd iterations: tiny probe increment `mSmallLoadIncrement·λprev`, re-solve eig,
  test convergence `|mLambda/mLambdaPrev| < mConvergenceRatio`
- even iterations: path-following advance `λprev += mPathFollowingStep·λ`
- loop continues until the smallest eigenvalue stabilizes → that is the critical
  buckling load factor. Load applied via `UpdateLoadConditions()` using `TIME` as
  the load multiplier.

### Python wiring (`structural_mechanics_prebuckling_solver.py`)
Needs **TWO** linear solvers:
- static path: `ResidualBasedEliminationBuilderAndSolver` + `_GetLinearSolver()`
  (forces *elimination* B&S; warns if you ask for block — buckling needs DOFs
  eliminated, not penalized)
- eigen: `_CreateLinearSolverEigen()` → `eigen_solver_factory.ConstructSolver(eigensolver_settings)`
- scheme: `ResidualBasedIncrementalUpdateStaticScheme`
- convergence criterion forced to `displacement_criterion`
- `buckling_settings` defaults: `initial_load_increment: 1.0`,
  `small_load_increment: 0.0005`, `path_following_step: 0.5`,
  `convergence_ratio: 0.05`; eigensolver default `"eigen_eigensystem"`,
  `number_of_eigenvalues: 5`.

**OpenSees porting note.** This is a near-drop-in recipe for an OpenSees
`bucklingAnalysis`:
1. OpenSees elements with geometric nonlinearity (corotational/PDelta transforms,
   `forceBeamColumn`, shells) already return a tangent that is K0+Kg at the current
   state — so, exactly like Kratos, you do **not** need a separate Kg assembler;
   run a `staticIntegrator`/`integrator LoadControl` step, grab the tangent from the
   `LinearSOE`/`AnalysisModel`, advance the load a hair, grab the second tangent,
   and difference them.
2. Feed (`K_ref`, `ΔK`) to the existing `eigen` machinery (`ArpackSOE` does the
   shift-invert generalized solve already). The only new C++ is a small "buckling
   integrator/analysis" that (a) captures two tangents around a load level and
   (b) wraps them as the generalized pair, plus the eigenvalue-tracking outer loop.
3. The smallest λ × current load factor = critical load; eigenvector = buckling
   mode. Reuse `Node::setEigenvectors` for the shape so existing recorders work.
4. Simpler v1: skip continuation, do a single two-snapshot difference at a chosen
   reference load (classic linear buckling). The Kratos continuation loop is the
   accuracy upgrade for v2.

---

## 3. HARMONIC / FREQUENCY-RESPONSE

Kratos ships **two** FRF strategies. Both treat the excitation angular frequency Ω
as `ProcessInfo[TIME]` and rely on an **external frequency-sweep loop** (the
`AnalysisStage` time loop advances "time" = frequency; each step is one Ω). Neither
strategy contains the sweep itself — the sweep is the outer stage loop.

### 3a. DIRECT harmonic (`direct_harmonic_analysis_strategy.hpp`)
Assembles a true **complex** dynamic stiffness per frequency in
`BuildComplexDynamicSystem(omega)`:
```cpp
AddRealMatrixToComplex(rA, rK, ComplexType( 1.0,        0.0));     //  K
AddRealMatrixToComplex(rA, rM, ComplexType(-omega*omega,0.0));     // -Ω²M
if (mAssembleDampingMatrix)
    AddRealMatrixToComplex(rA, rC, ComplexType(0.0, omega));       // +iΩC
// Rayleigh added directly without forming C:
if (mRayleighAlpha) AddRealMatrixToComplex(rA, rM, ComplexType(0.0, omega*mRayleighAlpha));
if (mRayleighBeta ) AddRealMatrixToComplex(rA, rK, ComplexType(0.0, omega*mRayleighBeta ));
```
→ **D(Ω) = K − Ω²M + iΩ(αM+βK+C)**. Then a single **complex** linear solve per Ω:
```cpp
mpComplexLinearSolver->Solve(*mpDynamicMatrix, *mpComplexSolution, *mpComplexRHS);
```
RHS supports separate real and imaginary load sub-model-parts
(`rb[i] = ComplexType(F_real[i], F_imag[i])`). Output: real part → displacement
DOF, imaginary part → `DISPLACEMENT_IMAGINARY`. `RAYLEIGH_ALPHA`/`RAYLEIGH_BETA`
read from `Properties`.

### 3b. MODAL harmonic (`modal_harmonic_analysis_strategy.hpp`)
Reuses the prior eigen run. In `Initialize()` it reads
`EIGENVALUE_VECTOR` and assembles a system-size×n_modes modal matrix `mpModalMatrix`
from each node's `EIGENVECTOR_MATRIX`. Per frequency:
```cpp
ComplexType factor( eigenvalues[i] - pow(Omega,2),
                    2*modal_damping*sqrt(eigenvalues[i])*Omega );   // ωⱼ²-Ω² + i·2ζⱼωⱼΩ
mode_weight = inner_prod(modal_vector, f) / factor;                 // qⱼ = φⱼᵀF / factor
...
modal_displacement[eq] += mode_weight * eigenvector_component;      // u = Σ qⱼ φⱼ
```
Per-mode damping blends sources:
```cpp
modal_damping = mSystemDamping
              + mRayleighAlpha/(2*sqrt(eigenvalues[i]))
              + mRayleighBeta*sqrt(eigenvalues[i])/2
              + mMaterialDampingRatios[i];      // strain-energy-weighted, per submodelpart
```
(`SYSTEM_DAMPING_RATIO` = uniform ζ; eigenvalue = ωⱼ²; material damping computed
per-mode from sub-model-part strain energy.) Output stores **magnitude** → DOF and
**phase** (arg) → reaction slot.

### Random / PSD
**Not present.** No PSD/random-vibration response in either strategy (no
Lyapunov/PSD integration). FRF is deterministic per-frequency only.

**OpenSees porting note.** OpenSees has *no* native FRF. The **modal** route (3b)
is by far the cheaper add and dovetails with the existing `eigen` command: after
`eigen N`, read `Node::getEigenvectors`, compute ωⱼ from the returned eigenvalues,
and evaluate `qⱼ(Ω) = φⱼᵀF / (ωⱼ²−Ω²+i·2ζⱼωⱼΩ)` over a Python/Tcl frequency loop —
this can be a **pure post-processing layer with zero C++** (a Python helper in
`Ladruno_implementation/`), mirroring 3b exactly. The **direct** route (3a) needs a
complex SOE/solver, which OpenSees lacks (its solvers are real); that is a larger
C++ effort (a `ComplexLinearSOE`) and should be deferred. Damping blend formula and
the magnitude/phase output convention port verbatim. Random/PSD is a genuine gap in
both codes — a differentiator if added.

---

## 4. COMPLEX / DAMPED MODAL

- No quadratic/state-space damped eigenproblem ([K, C, M] → 2N complex modes)
  anywhere in `StructuralMechanicsApplication`. Damped modal response is handled
  *approximately* via the real undamped modes + modal damping ratios (§3b), the
  standard proportional-damping assumption.
- True complex eigenvalues are available only through the **complex FEAST**
  backend `ComplexFEASTGeneralEigensystemSolver` (`solver_type: "feast_complex"`),
  which solves a *general* (non-symmetric, complex) A x = λ B x inside a complex
  contour — usable for a damped/non-symmetric pencil if you assemble it, but Kratos
  ships no strategy that forms the quadratic-eigenproblem companion linearization
  to feed it. So: the *solver* can do complex; no *structural workflow* uses it.

**OpenSees porting note.** Same situation as OpenSees, which also has only
real-symmetric modal + proportional modal damping. If complex/damped modes are
wanted, the portable pattern is: linearize the quadratic pencil
(λ²M + λC + K)φ=0 into a 2N×2N first-order generalized problem and hand it to a
*general* complex eigensolver. Kratos shows the backend half (complex FEAST) but
not the linearization — so this remains greenfield in both; design the companion
linearization once, target whatever general eigensolver each code has.

---

## 5. THE EIGENSOLVER BACKENDS (`LinearSolversApplication`)

All are `LinearSolver` subclasses with a `Solve(K, M, values, vectors)` overload.
Registered `solver_type` strings and problem classes:

| solver_type | class | problem | domain | backend |
|---|---|---|---|---|
| `eigen_eigensystem` | `EigensystemSolver` | K φ = λ M φ, K sym / M SPD | real | Eigen |
| `spectra_sym_g_eigs_shift` | `SpectraSymGEigsShiftSolver` | K φ = λ M φ, sym/SPD, shift-invert | real | Spectra |
| `feast` | `FEASTGeneralEigensystemSolver` | A x = λ B x (general) | real | Intel MKL / bundled FEAST |
| `feast_complex` | `ComplexFEASTGeneralEigensystemSolver` | A x = λ B x (general) | complex | FEAST |

### 5a. `EigensystemSolver` (default for buckling) — Bathe subspace iteration
Implements **subspace iteration** straight out of Bathe, *Finite Element
Procedures* (~p.954). For K φ = λ M φ:
1. pick nc-dim subspace (nc > n_roots)
2. iterate: solve `K r = M r⁽ʲ⁾` (one LU of **K only**, *not* K−σM — no shift),
   form reduced `ar = rᵀK r`, `br = rᵀM r`
3. solve the small dense projected problem with
   `Eigen::GeneralizedSelfAdjointEigenSolver`
4. converge on relative eigenvalue change vs `tolerance`.
Params (defaults): `number_of_eigenvalues:1`, `max_iteration:1000`,
`tolerance:1e-6`, `normalize_eigenvectors:false`, `echo_level:1`,
`use_mkl_if_available:true` (PARDISO else Eigen SparseLU).

### 5b. `SpectraSymGEigsShiftSolver` (default for modal) — Spectra shift-invert
Wraps the **Spectra** library's symmetric generalized shift-invert eigensolver
(implicitly-restarted Lanczos / Krylov-Schur family). Params:
`number_of_eigenvalues`, `shift` (compute smallest eigenvalues > shift —
the shift-invert target), `max_iteration`, `echo_level`. Best for "lowest few
modes of a large sparse SPD pencil".

### 5c. `FEASTGeneralEigensystemSolver` — contour integration
FEAST = **density-matrix / contour-integration** eigensolver, explicitly *not*
Krylov/Lanczos (Kratos docs quote this). Finds all eigenpairs whose eigenvalue lies
inside a search window:
- real: interval `[e_min, e_max]` (for structural, ω²=λ so this is a frequency band)
- complex (`feast_complex`): circular contour `e_mid_re/e_mid_im` + radius `e_r`

Subspace size `M0` auto-set: `2·n_eig` if extremal search, else `1.5·n_eig`.
The wrapper calls the FEAST Fortran kernel once (no user RCI loop — RCI is internal
to FEAST):
```cpp
feast(&UPLO,&N, A,IA,JA, B,IB,JB, fpm,&epsout,&loop, Emin,Emax,&M0, E,X,&M,res,&info);
// dfeast_scsrgv / zfeast_scsrgv (sym), dfeast_gcsrgv / zfeast_gcsrgv (general)
// fpm[0]=echo, fpm[2]=-log10(tol), fpm[3]=max_iter, fpm[14]=right-vecs,
// fpm[39]=±1 lowest/highest
```
Params: `e_min,e_max` (or contour), `number_of_eigenvalues`, `subspace_size`,
`symmetric`, `search_lowest_eigenvalues`/`search_highest_eigenvalues`, `tolerance`
(doc default 1e-12), `max_iteration` (doc default 20). Requires
`-DUSE_EIGEN_FEAST=ON` + MKL (or the bundled FEAST 4.0 in
`external_libraries/FEAST/`).

Backend selection in Python: `eigen_solver_factory.ConstructSolver(settings)`
dispatches on `solver_type`.

**OpenSees porting note.** OpenSees's generalized-eigen backends are Arpack
(`ArpackSOE`/`ArpackSolver`, shift-invert Lanczos — the analog of 5b),
`SymBandEigenSolver`, and `FullGenEigenSolver`. So OpenSees already *has* the
shift-invert Lanczos path Kratos calls "default modal"; no new backend is needed
for basic modal. The genuinely new capability Kratos exposes is **FEAST contour
integration**: "give me every mode whose frequency falls in band [f1, f2]" — which
is exactly the query seismic/HVAC FRF users want and Arpack does awkwardly (you
shift near a target and hope). FEAST 4.0 is BSD-ish and *bundled inside Kratos as
portable Fortran* — it could be vendored into OpenSees `OTHER/` the same way LAPACK
is, then wrapped as a new `EigenSolver`. This is the cleanest single porting target
for a "modes-in-a-band" feature.

---

## 6. PARALLEL EIGENSOLVING — the most valuable pattern

**Finding: Kratos has no Trilinos/Anasazi/SLEPc/PETSc-SLEPc eigensolver.**
(`gh api search/code` for `Anasazi` and for `eigen` under `TrilinosApplication`
both returned zero hits; the TrilinosApplication README documents only linear
solvers — Epetra/AztecOO/Amesos/ML — no eigen.) So the Kratos "parallel eigen"
story is *not* the distributed-Krylov route one might expect. It is two-pronged:

### 6a. Shared-memory (the common case)
The Eigen/Spectra/serial-FEAST backends parallelize via **OpenMP** inside the dense
projected solves and the sparse factorization (MKL PARDISO threaded). The Kratos
`block_for_each` machinery threads the *assembly* of K and M. This covers
single-node multicore, which is how the vast majority of Kratos eigen runs go.

### 6b. Distributed-memory: **PFEAST** (MPI FEAST) — the real distributed pattern
The bundled FEAST 4.0 includes its **parallel** layer (confirmed present:
`external_libraries/FEAST/4.0/include/pfeast.h`, `pfeast_sparse.h`,
`pfeast_tools.h`, and Fortran `src/sparse/pdzfeast_sparse.f90`,
`pdzifeast_sparse.f90`, `pdzfeast_pev_sparse.f90`). PFEAST's parallelism is
**three-level MPI**, which is the pattern worth emulating:

1. **L1 — search-interval / contour parallelism**: split the spectrum into several
   sub-intervals (or contours), one MPI sub-communicator per interval. Embarrassingly
   parallel across intervals.
2. **L2 — quadrature-point parallelism**: FEAST evaluates a contour integral
   `(zⱼ B − A)⁻¹` at ~8–16 Gauss points on the complex contour; each quadrature
   point's *independent linear solve* is mapped to its own MPI sub-communicator.
   This is the natural, near-perfect parallelism FEAST is famous for — the modes
   come from a sum over independent solves, no global Krylov synchronization.
3. **L3 — within each linear solve**: the distributed sparse solve at one
   quadrature point runs across the ranks of that sub-communicator (MKL
   CPardiso / MUMPS).
Then a *small* dense Rayleigh-Ritz on the M0-dim subspace (reduced, cheap, can be
replicated). FEAST iterates the whole thing (the "loop" counter) until residuals
drop — a handful of outer iterations, *not* hundreds of Lanczos restarts.

Why this beats distributed Lanczos/Anasazi for OpenSees's purposes: the dominant
cost (the linear solves) is **already** what OpenSees does well in parallel
(`OpenSeesMP` + MUMPS via `MumpsParallelSolver`/`MumpsSOE`). FEAST turns
"distributed eigen" into "many independent distributed *linear* solves + a tiny
dense problem" — reusing exactly the distributed-solver infrastructure OpenSees
already has, with MPI-subgroup orchestration on top.

**OpenSees porting note (highest value).** OpenSees has no parallel eigensolver;
`eigen` under `OpenSeesMP` is effectively serial/replicated. The PFEAST pattern is
the cleanest open-source template to close that gap *without* adopting Trilinos:
- OpenSees already ships **MUMPS** distributed solves (`SRC/system_of_eqn/linearSOE/mumps/`,
  `MumpsParallelSOE`/`MumpsParallelSolver`) and runs under MPI (`OpenSeesMP`,
  `PartitionedDomain`).
- A `ParallelFEASTEigenSolver` would: (a) form the complex contour around the target
  frequency band, (b) split `MPI_COMM_WORLD` into per-quadrature-point
  sub-communicators, (c) at each contour point solve `(zⱼM − K) X = M Y` with the
  *existing* MumpsParallel solver on that sub-comm, (d) sum into the subspace,
  (e) reduced Rayleigh-Ritz, (f) iterate. The bundled FEAST/PFEAST Fortran is
  vendorable (like LAPACK in `OTHER/`); its RCI form lets you supply OpenSees's own
  distributed solver as the inner solve.
- This reuses the proven `Mumps*` + `PartitionedDomain` stack, adds only the
  MPI-subgroup splitting + contour bookkeeping, and gives "all modes in a frequency
  band, in parallel" — directly relevant to large soil-structure / DRM models.
- Contrast: a distributed block-Lanczos (Anasazi-style) would require a parallel
  Krylov basis with global orthogonalization each step — much more invasive and
  exactly the route Kratos itself *declined* to take. The Kratos evidence is a
  strong vote for the FEAST/contour pattern over distributed-Krylov for a code that
  already has good distributed *linear* solvers.

---

## 7. Condensed crosswalk (Kratos → OpenSees)

| Capability | Kratos mechanism | OpenSees today | Port verdict |
|---|---|---|---|
| Modal K φ=λMφ | `EigensolverStrategy` + Spectra shift-invert | `eigen` + Arpack | parity; lift `normalize`+diag-regularization options |
| Mass-normalize modes | `MassNormalizeEigenvectors` flag | implicit/per-solver | make explicit option |
| Linear buckling | `PrebucklingStrategy` (Jia-Mang, two-tangent difference + continuation) | **none** | **add**: difference two static tangents → existing eigen pair; reuse corot/PDelta Kg |
| Direct FRF | complex D(Ω)=K−Ω²M+iΩC, complex solve/Ω | **none** | needs complex SOE → defer |
| Modal FRF | `qⱼ=φⱼᵀF/(ωⱼ²−Ω²+i2ζⱼωⱼΩ)` superposition | **none** | **add as pure post-proc** on `eigen` output (zero C++) |
| Modal damping blend | system+Rayleigh+per-mode material ζ | n/a | port formula |
| Complex/damped modal | only via `feast_complex`; no QEP linearization | none | greenfield both; design companion linearization |
| Modes-in-a-band | FEAST contour integration (bundled Fortran) | none | **vendor FEAST → new EigenSolver** |
| Parallel eigen | **PFEAST** 3-level MPI (interval/quadrature/solve); **no Anasazi** | none | **vendor PFEAST; inner solve = existing MumpsParallel; split MPI sub-comms** |

---

## 8. Bottom line for the OpenSees feature

1. **Buckling** is the cheapest high-value add and Kratos hands you a tested recipe:
   difference two geometric-nonlinear static tangents and feed the pair to the eigen
   machinery OpenSees already has. No new Kg element code.
2. **Modal FRF** can ship as a zero-C++ Python/Tcl post-processor over `eigen`.
3. **Parallel eigen** should follow **FEAST/PFEAST**, *not* Trilinos/Anasazi —
   Kratos itself made that choice. FEAST recasts distributed eigen as independent
   distributed *linear* solves over contour quadrature points, which maps perfectly
   onto OpenSees's existing MUMPS-parallel + PartitionedDomain stack. Vendor the
   bundled FEAST Fortran (LAPACK-style) and wrap it as an `EigenSolver` /
   `ParallelFEASTEigenSolver`.
