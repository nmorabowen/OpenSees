---
title: "WP-260 orchestration plan — work packages, waves, parallelism, and dual-agent protocol"
project: Ladruno
status: decisions D1-D7 recorded (2026-11-26); ready for owner merge
program_issue: "#1156"
companion: 260_ladruno_solver_acceleration_adr
updated: 2026-11-26
---

# WP-260 Orchestration Plan (`Solver & Wave Acceleration`)

## 0. Dual-Agent Mechanics (Gemini Author + Claude Reviewer)

- **Gemini / Antigravity** discovers workspace skills in `.agents/skills/ladruno-solver-accel/` and runs in-workspace worktrees under `.agents/worktrees/wp260*` (gitignored via `.gitignore`: `.agents/*` except `.agents/skills/`), keeping all file edits, `build.bat`, static gates and `pytest` inside its worktree.
- **Claude** does not author WP-260 code. It runs the adversarial bug-finding review on each Draft PR before ready-flip. That review is **not** the ADR-87 D4 out-of-family review: D4 is the owner-run handoff pack + committed verdict (ADR §7, decision D7).
- **GitHub holds all live state** (`Ladruno_internal/ISSUES.md`): one `program` issue (#1156), board `#937`, and one Draft PR per WP (`wp/260<x>-<slug>`) with `### Handoff` comments.
- **Local Windows builds are serialized:** never run two `Ladruno_scripts\build.bat` builds at once on the same workstation; run static gates and coding in parallel, then queue local 5-target builds one worktree at a time before firing parallel Zone-A Ubuntu CI runs.

---

## 0.5 Selection Rule — Measure First, Then Pick 2–3 (review recommendation, adopted by the owner 2026-11-25)

The work-package map below is a **menu, not a schedule**. Every speedup figure in it is a hypothesis (literature or code inspection) until `260-w0` measures where the time actually goes on the D1 decks. Earlier fork ADRs found that steps which looked expensive took little of the runtime once measured (ADR-68: the J2 tangent flag 3.9 %; ADR-77: `/arch:AVX2` nothing, implicit Bathe below the bar), so the program commits work only behind measurements.

1. **`260-w0` first.** Nothing beyond `260-w0` starts until its phase-share table (§`260-w0` brief) is posted on the `program` issue.
2. **Selection gate.** A work package is scheduled only if, on at least one D1 deck:
   - the phase it targets takes **p ≥ 10 %** of the run (from the 1-thread `-deep` profile), and
   - its projected whole-run gain **G = 1 / ((1 − p) + p / s) ≥ 1.10×**, where `s` is the WP's own *conservative* phase speedup (the low end of its hypothesis, halved if unmeasured anywhere).
   The `260-w0` post ranks every Wave-1/2 package by `G` per deck.
3. **Pick 2–3 at a time.** The owner schedules the top 2–3 packages by `G` (ties broken by risk: `G-BYTE` before `G-EQUIV` before `G-OPT`). The rest stay on the menu.
4. **Ride-alongs.** A cheap `G-BYTE` hygiene item (≤ 2 days, e.g. the dead loop in `KrylovAccelerator.cpp`) may ride along with a selected package that already owns the same files, without passing the gate.
5. **Re-rank after each merge.** Each merged package posts its measured whole-run gain on the D1 decks; `260-w0`'s table is updated and the menu re-ranked before the next pick. A package whose measured `G` lands below half its projection triggers a short "why" note on the `program` issue before the next pick.
6. **Wave 3 needs a written case.** A Wave-3 package starts only after the Wave-1/2 packages feeding it have merged, its phase passes the gate on the re-ranked table, and the owner has answered its decision (D3 for GPU).

**First-pick expectation (to be confirmed by `260-w0`, not assumed):** factorization (`soe.factor`/`soe.symbolic`) usually dominates large 3D implicit decks → `260a` (PARDISO levers) / `260g`; element + material time (`elem_by_type`, `commit.elements`) dominates finite-strain and soil decks → `260c1` / `260c2` / `260f`; `update`/`formUnbalance` dominate explicit decks → `260b`; `dc.setSize` / `soe.addA` show up on Mode-A MPM → `260a`'s MUMPS items.

---

## 1. Work Package Map & File Ownership

Ownership is disjoint **between WP-260 packages**. It is **not** disjoint from the open PRs #1147/#1153/#1154 (ADR §2.2): a package that owns a file those PRs edit is marked *blocked-by* and is cut from `origin/ladruno` only after that PR merges.

### Wave 0 — Bootstrap & Phase-0 Profiler Baseline
| WP | Branch | Content | Dependencies | Gate |
|---|---|---|---|---|
| **260** | `wp/260-solver-accel-program` | Program ADR, this plan, `ladruno-solver-accel` skill, `.gitignore`/`.agents/` setup, `AGENTS.md` row, quirk-pointer lint widened to `.agents/skills/`, 3 quirk fragments | None (merges first) | Static gates green; Claude review; owner D1–D7; ADR-87 D4 verdict if D7 says so |
| **260-w0** | `wp/260w0-phase0-profiles` | Profiler sweep on the D1 decks (`origin/ladruno` + `2CLD12`): phase-share table, phase→WP map and the §0.5 ranking by projected whole-run gain, posted on the `program` issue (no `SRC/` edits) | `260` + D1 | Table + ranking posted; owner picks the first 2–3 WPs |

### Wave 1 — Surgical Wins (Mostly `G-BYTE`) — candidates for the first pick (§0.5)
| WP | Branch | Scope & File Ownership | Blocked by | Speedup hypothesis (unmeasured; replaced after `260-w0`) | Gate |
|---|---|---|---|---|---|
| **260a** | `wp/260a-soe-scatter-hygiene` | `MumpsSOE::addA` binary search + skip `rowA++/--` on `JOB=3` + `std::sort` `setSize()` (incl. local-vertex `MumpsParallelSOE`); `PARDISOGenLinSOE::setSize` `std::sort`; `KrylovAccelerator.cpp:200-216` dead loop; fix the `PARDISOGenLinSolver.cpp:366-368` comment (re-prove `-deterministic`) & expose opt-in `iparm[23]/iparm[1]/iparm[24]=2/iparm[10/12]` (refused or proved under CNR; oneMKL 2024.2.2 floor via `check_mkl_compat.py`); `AcceleratedNewton::sendSelf` reuse serialization | — | 3–20× `Mumps`/`setSize`, 1.3–2× PARDISO factor (to be measured), faster `KrylovNewton` | `G-BYTE` (default) + `G-OPT` (flags) + MKL compat |
| **260b** | `wp/260b-explicit-wave-hygiene` | `ExplicitBathe` CDL mass-cache + `-commitSolveState` + hoisted `-lnvd`; ADR-67 P1.1/P1.2′ consistent-SMS & `MPIDiagonalSOE` alloc/halo hoist + rank-1 $O(n)$ centroidal matvec + `CriticalTimeStep` massless check; T7 inertia skips (`LadrunoElementShell.h`, `Brick20`, `CPML`); `LadrunoBrick::formUri`/BV $c_d$ cache; `ConstraintProjector` Cholesky $L$; global $O(N)$ $\alpha_M$ (`TransientIntegrator.cpp`); `LadrunoCPML` OpenMP opt-in | — | 1.4–2.2× explicit (`CDL`/`ExplicitBathe`), 4–8× consistent-SMS PCG | `G-BYTE` / `G-EQUIV` |
| **260c1** | `wp/260c1-logstrain-stressonly` | `LogStrainNDMaterial::setTrialF()` 5-bug fix; Contract #1a on `NDMaterial.h`/`FiniteStrainNDMaterial.h` + `LogStrainNDMaterial`/`LogStrainKernel.h` (all tangent channels incl. `getSpatialTangentTensor`); `LadrunoUP` OpenMP opt-in after guarding the `:599` static | — | 3–6× finite-strain GP | `G-BYTE` / `G-EQUIV` |
| **260c2** | `wp/260c2-concrete-hs-implex` | Contract #1a on `LadrunoConcrete3D*`/`LadrunoHS*`/`LadrunoRCConcrete*`; `doTangent && !mp.implex` (`LadrunoConcrete3DKernel.h:3253`); commit-deferred `-implex` in `Concrete3D` & `HS` (`Ge_implex` freeze) on the existing `LadrunoImplexClock` (Contract #1b) | **#1153, #1154** | 1.5–4× unused-tangent skips; 4–15× `-implex` | `G-BYTE` / `G-EQUIV` |
| **260d** | `wp/260d-adaptive-transient` | `LadrunoAdaptiveTransient` (HAFTOL half-step LTE + quantized $2^k/\sqrt{2}^k$ rungs) + Contract #2 `getEffectiveTangentKey` (incl. `TRBDF2`/`TRBDF3`) + vectorized `addModalDampingForce`. **Owned:** `SRC/analysis/integrator/IncrementalIntegrator.{h,cpp}`, `TRBDF2.{h,cpp}`, `TRBDF3.{h,cpp}`, new `LadrunoAdaptiveTransient*`, `SRC/analysis/algorithm/equiSolnAlgo/LadrunoFactorReuse*`. `LadrunoImplicitBathe` + Contract #3 are **not** in 260d: they wait for D6 (260d-B below) | — | 2–5× adaptive runs (to be measured) | `G-EQUIV` / `G-OPT` |
| **260e** | `wp/260e-allocator-lto-omp` | Unify Windows OpenMP on `libiomp5md` (remove the `vcomp` dual pool, `CMakeLists.txt:774-781`, verified via `dumpbin`); evaluate/wire `mimalloc` + `/GL /LTCG` (per D4). **No `/arch:AVX2`** (ADR-77 V0) | D4 | to be measured (#933's heap-contention finding is the hypothesis) | `G-BYTE` + `dumpbin` check |

**Gated (not scheduled until its decision):**
- **260d-B** `wp/260dB-implicit-bathe` — D6 benchmark first (shared factor across substeps under `LadrunoFactorReuse`, on the ADR-77 T3 deck + a contact/wave deck, ≥1.3× bar). Only if it clears: `LadrunoImplicitBathe` + Contract #3 (`ScopedRecorderSuppression`, commit-side-effect table, substep driver, no-commit substep 1). Owns `SRC/domain/domain/Domain.{h,cpp}` and the new integrator/driver files.

### Wave 2 — Core Kernels, Solvers & SSI Upgrades (3–8 Weeks) — menu, scheduled by §0.5
- **`260f` (3D/MPM/MCFT Threading + Persistent Assembly):** De-static & opt-in `LadrunoBrick*`, `Bezier*`, `LadrunoSolidShell`, Mode-A Implicit MPM, Model C `LadrunoNDFiberSection3d`/`MCFTStrip`; `LadrunoFastScatter` (Contract #4, NVI `setSize`/`doSetSize`) + persistent `LadrunoThreadedAssembly` + colored scatter + sparse `TransformationFE` $T^T K T$ + `Domain` snapshot cache. *Blocked by #1147 (MPM, `LadrunoThreads`) and #1153 (fiber sections).*
- **`260g` (Mixed-Precision GMRES-IR & Sym-Preconditioned Krylov):** FP32 PARDISO (`iparm[27]=1` `float*` shadow) + MKL RCI `dfgmres` FP64 GMRES-IR; `LadrunoMumpsMixed` GMRES-IR upgrade; `LadrunoSymPrecondKrylovSOE` ($K_{\text{sym}}$ factor + FGMRES). oneMKL 2024.2.2 floor via `check_mkl_compat.py`.
- **`260h` (Analytic Return Maps & Outer Loop):** `Concrete3D` analytic `damagedTangent` & Voigt lift + Cardano/Smith 3×3 eigensolver; `LadrunoHS` analytic R4/R5 Jacobian & IFT Newton–Brent; `SANISAND` `SAS-RK23/45`; `ASDPlastic3D` closed-form VM/DP; Eisenstat–Walker Inexact Newton (Contract #4 `LinearSOESolver::setRelativeTolerance`) + `LadrunoAndersonNewton` + two-loop Powell-damped `LadrunoLBFGS`; optional modal Woodbury & static PTC. *Blocked by #1153/#1154 (SasME, ASDConcrete3D, Concrete3D, HS).*
- **`260i` (Global IMPL-EX + `LadrunoSuperelement` + `H5DRM` v2):** `LadrunoImplexAlgo`/`LadrunoImplexLinear` on Contract #1b; `LadrunoSuperelement` wrapping `OPS_LadrunoCMS`; velocity-consistent & Biot $u\text{–}p$ `H5DRMLoadPattern` v2.
- **`260j` (`LadrunoContact` Acceleration):** Slave-aware corner-span bucket sort / BVH for 3D Mortar (`LadrunoContactHandler.cpp:2345-2352`), MPI `-reemit`, threaded `LadrunoContactFE`, ALM active-set audit.

### Wave 3 — Transformative Engines & HPC Frontier (2–4+ Months) — menu, needs a written case (§0.5 rule 6)
- **`260k` (`LadrunoGPStore` SoA + 7×7 Schur SANISAND CPPM):** 64-byte SoA + two-pass active-set compaction + Morton/RCM internal cache ordering.
- **`260l` (`LadrunoExplicitEngine` + Clustered LTS):** WP-232 Mode-F `LoadPattern` flat explicit engine + flat `CPML`/`Contact`/`EmbeddedNode`/`PorousOverlay` + power-of-two LTS subcycling + FP32-force/FP64-state (CPML memory variables FP64).
- **`260m` (`LadrunoAMGSOE` + Monolithic Biot Schur + Cluster Direct):** Refreshed Boost-free `OTHER/AMGCL` (6 rigid-body near-nullspace modes + 3×3 block SA-AMG + monolithic `LadrunoUP` fixed-stress Schur); `LadrunoClusterPardisoSOE`; weighted METIS `vwgts`; optional `cuDSS` (per D3).
- **`260n` (`LadrunoIMEX` Spatial Partitioning):** Phase A Hughes–Liu pruned-LU IMEX → Phase B Prakash–Hjelmstad dual-Schur multi-rate subcycling atop `260l`.
- **`260o` (`LadrunoSpectralHex` + Pure-$M$ Integrators + GPU Explicit/MPM):** $p=3\text{–}5$ GLL SEM + pure-$M$ Hulbert–Chung / 2-stage Soares + `LADRUNO_HD` CUDA/HIP backend (only if D3 puts GPU in scope).

---

## 2. Landing & Review Cadence

1. Merge **`WP-260`** first after Claude review + owner D1–D7 answers (+ the D4 verdict file if D7 requires it).
2. Run **`260-w0`** and post the phase-share table + ranking on the `program` issue. The owner picks the top **2–3** packages by projected gain (§0.5); only those launch, in parallel worktrees under `.agents/worktrees/`. **`260c2`** can be picked only after #1153/#1154 merge; **`260d-B`** only after D6; **`260e`** only after D4.
3. Serialize local `build.bat` per worktree → run the `ladruno.yml` static gates (+ the MKL 2024.2 compat job for MKL edits) + local pytest + before/after profiles on the D1 decks → push Draft PRs → run Zone-A Ubuntu CI in parallel → post `### Review request` → Claude review (§3) → address findings → owner merge commit.
4. After each merge: post the measured whole-run gain, update the `260-w0` table, re-rank the menu (§0.5 rule 5), then pick the next package.

---

## 3. Review Protocol (Claude Reviews Every WP-260 PR)

Gemini authors; **Claude reviews**. The full checklist is `.agents/skills/ladruno-solver-accel/references/review-protocol.md`. In short:

1. **Request.** When the warrant package is green (static gates, local 5-target build, Zone-A, before/after profiles), the author posts a `### Review request — 260<x> — <head sha>` comment on the draft PR: what changed, contracts touched, the gate type claimed per change, and the measured table.
2. **Review.** Claude reviews the PR at that SHA against the checklist (scope and ownership, correctness gate proof, contract conformance, performance evidence, fork rules, threading) and posts **one** PR review comment with numbered findings and a verdict: `ready`, `changes needed`, or `stop — ADR amendment needed`.
3. **Response.** The author answers every finding by number: `fixed in <sha>` or `declined: <rationale>`. Claude re-checks the fixes and posts an updated verdict. No ready-flip while any finding is unanswered.
4. **Program-level review.** After each merge Claude also reviews the re-ranked menu (§0.5 rule 5) before the next pick, and flags any package whose measured gain fell far below its projection.

