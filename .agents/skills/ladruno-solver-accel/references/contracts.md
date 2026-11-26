# WP-260 Cross-Pillar Contracts, Precision Policy & Tier-3 Traps

Summary of ADR §4–§6. The ADR text is authoritative; contracts freeze when the owner accepts the ADR.

## 1. Six Cross-Pillar C++ Interface Contracts

| # | Seam & Files | Signature / Contract | Rules |
|---|---|---|---|
| **1a** | **Stress-only evaluation** (`SRC/material/nD/NDMaterial.h`, `FiniteStrainNDMaterial.h`) | `virtual int setTrialStrainStressOnly(const Vector&)`, `setTrialFStressOnly(const Matrix&)`, `setTrialFincrStressOnly(const Matrix&)`, `virtual bool ladrunoSupportsStressOnly() const` (default `false`) | Defaults forward to `setTrialStrain` / `setTrialF` / `setTrialFincr`. Covers **every** tangent channel: `getTangent()`, `getSpatialTangentTensor()`, and wrappers reading an inner tangent. The stress-only trial caches what the tangent needs; producing the tangent later is a pure function of that cache (never re-runs the return map from committed state, never returns the previous iteration's tangent). Writes only the instance's own members; only for `ladrunoThreadSafeForm()` materials with one thread per instance during form. `tangentValid_` resets on commit/revert. |
| **1b** | **IMPL-EX control on the existing API** | Reuse `getImplexError()` + `LadrunoImplexClock`; add `virtual double getImplexError() const` only where missing (ASD materials keep their `const Vector&` accessor plus a differently-named `double` adapter); add `supportsImplex()`, `setImplexMode(bool enable, double alpha = 1.0)` | A material declared with `-implex`/`-implexControl` keeps its own clock; `setImplexMode` on it is a warn-once no-op. Never two extrapolations. `getImplexElasticDrift()` dropped. |
| **2** | **Effective tangent key** (`SRC/analysis/integrator/IncrementalIntegrator.h`) | `virtual bool getEffectiveTangentKey(double &c1, double &c2, double &c3) const { return false; }` | `LadrunoFactorReuse` keys on `(SOE, integrator, neq, c1, c2, c3)` when true, today's `(SOE, integrator, neq, dt)` when false. `TRBDF2`/`TRBDF3` should implement it. |
| **3** | **Substep commit without recorders** (`SRC/domain/domain/Domain.h/.cpp`) — gated by D6 | `int suppressRecorderDepth = 0` + RAII `Domain::ScopedRecorderSuppression`; skips only the recorder loop and `commitTag++` | **Separate from `contactAugmenting`** (that flag also drives contact augment-mode commit and gates the ADR-60 re-emit). Per side effect under a substep commit: element/node commit runs; WP-99 refusal latch runs; contact commit runs with substep `dT`, never augment mode; `needsResort` re-emit deferred to the step-final commit; `LadrunoPorousOverlay::onDomainCommit` **skipped** (fluid advances once per full step); recorders/`commitTag` skipped. Rollback: don't commit substep 1 (hold as trial), or keep a step-start snapshot for **all** materials. Two Newton solves per step need a fork-owned substep driver, not the integrator calling the algorithm. |
| **4** | **Fast CSR scatter & inexact Newton** (`SRC/system_of_eqn/linearSOE/LinearSOE.h`, `LinearSOESolver.h`) | Non-virtual `LinearSOE::setSize(Graph&)` → subclass `doSetSize(Graph&)` → `++patternGeneration_` (NVI); `locateEntry(row,col) const`, `addAMapped(const double*, int nr, const IndexT *map, double fact=1.0)`, `supportsMappedScatter() const`; `LinearSOESolver::setRelativeTolerance(double eta)` | Plans keyed on `(const LinearSOE*, patternGeneration_)`. A subclass still overriding `setSize` directly reports `supportsMappedScatter()==false`. `IndexT` = the SOE's own index type (`int64_t` for ILP64/`pardiso_64`). Default `addAMapped` falls back to `addA`; default `setRelativeTolerance` is ignored. |
| **5** | **64-byte SoA `LadrunoGPStore`** | 64-byte-aligned `eps[6][N]`, `sig[6][N]`, `Dtan[36][N]` (null when stress-only), optional lazy `statePrevCommit`, `syncToDomain()`, `evalBatchStressOnly(W, eps, state, sig)`, two-pass compacted `plasticIdx[]` | Keep the virtual `NDMaterial*` fallback. Zero allocations inside `evalBatch`. `*Kernel.h` stay plain C++: no allocation, no `opserr`, no statics, no virtual calls. `LADRUNO_HD` (`__host__ __device__`) only if D3 puts GPU in scope. |
| **6** | **Windows toolchain / OpenMP** (`CMakeLists.txt`) | Single OpenMP runtime (`libiomp5md`, zero `vcomp140.dll`); `/fp:precise` / `-ffp-contract=off` for `G-BYTE`; `mimalloc` + `/GL /LTCG` per D4 | **No global `/arch:AVX2`** (ADR-77 V0). Verify with `dumpbin /dependents` on every Windows DLL/EXE in `dist\`. |

---

## 2. Unified Mixed-Precision Policy

- **Allowed (opt-in `G-OPT` / `G-EQUIV`):**
  - Implicit direct solvers: FP32 MUMPS/PARDISO factorization (**for PARDISO `iparm[27]=1`, keep `float*` shadows `A_f`, `b_f`, `x_f`**) + **FP64 outer GMRES-IR** (Carson–Higham, up to $\kappa \lesssim 10^{11}$ with automatic FP64 fallback).
  - Explicit dynamics: FP32 element kinematics / elastic / CPML / J2 **force evaluation** + **FP64 nodal state accumulation** (`u += dt * v_half`).
- **Strictly FP64:**
  - Any state accumulated over steps: nodal state, **CPML auxiliary memory variables** (recursive convolution), plastic/damage internal variables.
  - `LadrunoSANISAND` (`p' -> p_min`), `LadrunoConcrete3D` softening, `LadrunoHS` corners, and un-accumulated long-duration wave time-stepping.

---

## 3. Tier-3 Traps (Do Not Build)

1. Parareal / MGRIT / PFASST for hyperbolic wave propagation or hysteretic plasticity.
2. Dynamic h-AMR / p-adaptivity inside OpenSees `Domain` / `LinearSOE`.
3. High-order implicit SDIRK / ESDIRK (order 3–5) for seismic/contact dynamics.
4. Unpreconditioned JFNK on stiff heterogeneous SSI/RC/contact.
5. Rosenbrock / W-methods over Global IMPL-EX.
6. Implicit GPU offload before SoA/assembly fixes.
7. Stock `HHTGeneralizedExplicit` as a diagonal-mass explicit wave integrator (`C` added to LHS).
8. Global `/arch:AVX2` (ADR-77 V0: no speed, "do not add the flag") and standalone `LadrunoJ2` `tang_flag` or closed-form plane-stress J2 (ADR-68 T1/T2).

## 4. Closed Unless the Owner Re-Opens

- **Implicit Bathe** (ADR-77 G3, closed with the benchmark on record). Re-opens only via D6's pre-registered shared-factor benchmark clearing ADR-77's ≥1.3× bar.
