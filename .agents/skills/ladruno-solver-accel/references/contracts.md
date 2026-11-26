# WP-260 Locked Cross-Pillar Contracts, Precision Policy & Tier-3 Traps

## 1. Six Locked Cross-Pillar C++ Interface Contracts

| # | Seam & Files | Exact Signature / Contract | Rules |
|---|---|---|---|
| **1** | **Stress-Only & Unified IMPL-EX Hooks** (`SRC/material/nD/NDMaterial.h`, `FiniteStrainNDMaterial.h`) | `virtual int setTrialStrainStressOnly(const Vector &v);` `virtual int setTrialFStressOnly(const Matrix &F);` `virtual bool ladrunoSupportsStressOnly() const;` Plus unified IMPL-EX: `supportsImplex()`, `setImplexMode(bool enable, double alpha = 1.0)`, `getImplexStressGap(double &gapNorm2, double &refNorm2)`, `getImplexElasticDrift()`. | Default `setTrialStrainStressOnly` forwards to `setTrialStrain`. Ladruno implementations compute stress/state with `doTangent = false` and set `tangentValid_ = false` (reset on commit/revert). If `getTangent()` is called when `!tangentValid_`, materialize `Dtan` lazily on demand (never return stale tangent). |
| **2** | **Effective Tangent Key** (`SRC/analysis/integrator/IncrementalIntegrator.h`) | `virtual bool getEffectiveTangentKey(double &c1, double &c2, double &c3) const { return false; }` | `LadrunoFactorReuse` compares `(c1, c2, c3)` when true (falling back to `dt` comparison when false). Allows `LadrunoImplicitBathe` ($\gamma_0(\rho_\infty)$) to share 1 LU/AMG factor across both substeps. |
| **3** | **Commit-vs-Recorder Decoupling** (`SRC/domain/domain/Domain.h/.cpp`) | Ref-counted `int suppressRecorderDepth = 0` + RAII `Domain::ScopedRecorderSuppression` (generalizes `contactAugmenting`). | Commits material/contact state cleanly during composite substeps (`LadrunoImplicitBathe`) and subcycled explicit/IMEX zones without firing recorders or bumping `commitTag`. |
| **4** | **Fast CSR Scatter & Inexact-Newton Hooks** (`SRC/system_of_eqn/linearSOE/LinearSOE.h`) | `uint64_t patternGeneration_ = 0;` incremented in `setSize()`; `locateEntry(int row, int col) const`, `addAMapped(const double *mColMajor, int nr, const int32_t *map, double fact = 1.0)`, `supportsMappedScatter() const`; `setRelativeTolerance(double eta)`. | Default fallbacks preserve 100% compatibility with un-upgraded SOEs. Map rebuilds automatically when `plan.patternGen != soe.patternGeneration_`. |
| **5** | **64-Byte SoA `LadrunoGPStore` & `LADRUNO_HD` Kernels** | Pure POD `*Kernel.h` headers compiled `__host__ __device__`; `LadrunoGPStore` holds 64-byte-aligned `eps[6][N]`, `sig[6][N]`, `Dtan[36][N]` (null when `needTangent=false`), optional lazy `statePrevCommit` (for $O(1)$ Substep-2 rollback), `syncToDomain()`, `evalBatchStressOnly(W, eps, state, sig)`, and Pass-1 SIMD predictor + Pass-2 compacted `plasticIdx[]`. | Always retain virtual `NDMaterial*` fallback for non-Ladruno materials. Zero dynamic allocations inside `evalBatch`. |
| **6** | **Windows Toolchain / OpenMP / ISA Contract** (`CMakeLists.txt`) | Single OpenMP runtime (`libiomp5md`, zero `vcomp140.dll`); `/fp:precise` / `-ffp-contract=off` for `G-BYTE` gates; AVX2 + `mimalloc` + `/GL /LTCG` governed by ADR §8 D2/D4. | Verify with `dumpbin /dependents` on Windows DLL/EXE targets. |

---

## 2. Unified Mixed-Precision Policy

- **Allowed & High-ROI (Opt-In `G-OPT` / `G-EQUIV`):**
  - Implicit direct solvers: FP32 MUMPS/PARDISO factorization (**for PARDISO `iparm[27]=1`, MUST maintain `float*` shadow `A_f`, `b_f`, `x_f`!**) + **FP64 outer GMRES-IR** (Carson–Higham, up to $\kappa \lesssim 10^{11}$ with automatic FP64 fallback).
  - Explicit dynamics: FP32 element kinematics / elastic / CPML / J2 force kernels + **FP64 nodal state accumulation** (`u += dt * v_half`).
- **Strictly Forbidden (FP64 Only):**
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
8. Standalone `LadrunoJ2` `tang_flag` or closed-form plane-stress J2 (demoted by ADR-68 T1/T2).
