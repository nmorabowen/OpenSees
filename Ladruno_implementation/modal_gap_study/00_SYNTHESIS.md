# OpenSees modal-analysis gap — cross-code synthesis

Deep-theory dive comparing OpenSees (Ladruno fork) modal/eigen capabilities against
**Abaqus**, **Kratos**, and **LS-DYNA**, including the **parallel (MP/MPI) eigensolver**
dimension. Source dossiers:

- [`01_opensees_current_state.md`](01_opensees_current_state.md) — ground truth (file:line audit of the fork)
- [`02_abaqus_theory.md`](02_abaqus_theory.md) — Abaqus Theory Guide §2.5 / §2.3 / §2.14 formulations
- [`03_kratos_source.md`](03_kratos_source.md) — Kratos C++ strategies + PFEAST (open-source template)
- [`04_lsdyna_theory.md`](04_lsdyna_theory.md) — LS-DYNA implicit/complex eigen + MPP

---

## 1. Ground truth — what OpenSees actually does today

- **Serial eigen**: `eigen` drives **serial ARPACK** (`dsaupd_`/`dseupd_`, reverse-communication
  Lanczos) through `ArpackSolver`. Backends: `-genBandArpack` (default), `-fullGenLapack`,
  `-symmBandLapack`. The shift-invert `(K−σM)⁻¹` solve delegates to whatever `LinearSOE` is active.
- **Parallel eigen — the nuance**: the RC plumbing *is* distributed. `PartitionedDomain::eigenAnalysis`
  fans to subdomains; each keeps its matrices local; the `(K−σM)⁻¹` solve becomes a real **parallel
  MUMPS** solve **if** the active SOE is `MumpsParallelSOE`. **But the Krylov driver is still serial
  ARPACK** — there is **no PARPACK / SLEPc / ScaLAPACK / distributed Lanczos**. And the two parallel
  models do not compose: `MumpsParallelSOE` is `_PARALLEL_INTERPRETERS` (OpenSeesMP) while subdomain
  `ArpackSOE` distribution is `_PARALLEL_PROCESSING` (OpenSeesSP). You cannot get
  partitioned-domain + parallel solve + parallel Lanczos in one build.
- **Post-processing**: `modalProperties` (participation, effective mass, total-mass check) and
  `responseSpectrumAnalysis` exist (Petracca) but are **serial-only** (see only the master domain
  under SP). No native CQC/SRSS combination command.
- **Buckling**: hooks exist but **`Kφ=λKgφ` returns zero — the "buckling" path is dead code**. No
  `LadrunoBuckle`.
- **Complex / damped modal**: **none**. Real symmetric undamped only.
- **Frequency domain (FRF/SSD/random)**: **none**.

---

## 2. The two convergent algorithms (every code agrees)

### 2a. Complex / state-space modal — DON'T solve the full quadratic problem
Both mature codes avoid the full $N$-dim quadratic eigenproblem $(\lambda^2M+\lambda C+K)\psi=0$:

- **Abaqus `*COMPLEX FREQUENCY`**: reuse the real undamped basis $\Phi$ (from a prior FREQUENCY
  step), **project** $\tilde M=\Phi^TM\Phi,\ \tilde C=\Phi^TC\Phi,\ \tilde K=\Phi^TK\Phi$ (with
  $\tilde C$ *full* because damping is non-classical), and solve the small $2p\times2p$
  **state-space generalized nonsymmetric (QZ)** problem. Then
  $\omega_0=|\lambda|,\ \omega_d=\mathrm{Im}\,\lambda,\ \zeta=-\mathrm{Re}\,\lambda/|\lambda|$,
  complex modes $\psi=\Phi z$.
- **LS-DYNA `LCPACK=3`**: makes the eigenproblem nonsymmetric and **switches to ARPACK**, folding
  damping (+ gyroscopic terms) into first-order state-space form → complex eigenpairs.
  **Eigenvalue `LCPACK=3` is SMP-only.**

→ **For OpenSees this is a low-core-risk quick win.** `eigen` + `modalProperties` already give
$\Phi$ and $\phi^TM\phi$. New pieces: an *assembled-C accessor*, a small **LAPACK `dggev`/`zggev`
QZ wrapper** on the $2p\times2p$ reduced problem, and complex post-processing. No assembly-core or
big-eigensolver change.

### 2b. Linear buckling / prestressed modal — Kratos gives the cleanest template
- **Abaqus**: rides the *linear-perturbation-about-base-state* architecture; stress stiffening is
  free because state carries across steps; solves $(K_0+\lambda K_\Delta)\phi=0$.
- **Kratos `prebuckling_strategy`** (Jia–Mang consistent linearization): run genuine **nonlinear
  static steps** (so each element's own geometric nonlinearity already builds $K_0+K_g$ in its
  tangent — **no separate $K_g$ assembler needed**), snapshot two tangents, **difference them**
  $\Delta K=K_\mathrm{ref}-K$, and solve the pair **$K_\mathrm{ref}\,\phi=\lambda\,\Delta K\,\phi$**
  with the standard eigensolver + an eigenvalue-tracking continuation loop.

→ **For OpenSees**: difference two static tangents from corot/PDelta elements, feed the pair to the
existing ARPACK `eigen`. Only new C++ is a small buckling analysis/integrator (`LadrunoBuckle`).
Reuses everything.

---

## 3. The parallel / MP answer (the user's key question)

**Both Kratos and LS-DYNA deliberately avoid a distributed Krylov eigensolver. They reduce
distributed eigen to *many independent distributed LINEAR solves*** — which is exactly the
capability OpenSees already has (`MumpsParallelSOE`).

- **LS-DYNA MPP**: no special distributed Lanczos. Parallelism lives **entirely in the distributed
  sparse factorization** of the shift-invert solves (`LSOLVR=7` multifrontal, or `LSOLVR=30` =
  **MUMPS**), with ParMETIS orderings. (Fast Lanczos `EIGMTH=103` is the MPP path for many
  approximate NVH modes.)
- **Kratos**: **no Trilinos/Anasazi/SLEPc eigensolver at all** (confirmed). Its distributed-eigen
  path is bundled **PFEAST (MPI FEAST)** — contour-integration eigensolver with **three-level
  parallelism**: search-interval × contour-quadrature-point × inner distributed linear solve. Each
  quadrature point is an *independent* distributed linear solve $(z_jM-K)^{-1}$, then a tiny dense
  Rayleigh–Ritz.

**Strategic implication — FEAST is a triple-win for OpenSees.** A FEAST/PFEAST `EigenSolver`
backend would:
1. **Parallelize naturally** over MPI subgroups (per contour point), each inner solve = an existing
   `MumpsParallelSOE` solve → **resolves the SP/MP incoherence** by routing both through the same
   linear-solve stack.
2. **Band-target** ("give me all modes in $[f_1,f_2]$") — the contour *is* the band. Closes the
   no-band-targeting + no-Sturm-guarantee robustness gaps in one move (FEAST's subspace count is
   self-certifying).
3. **Handle complex/damped** via complex FEAST contours → also a path to the §2a gap at scale.

The lower-friction alternative is **PARPACK** (reuse the existing RC loop), but both reference codes
*chose not to* go distributed-Krylov — a strong vote for the FEAST pattern.

---

## 4. Frequency domain (FRF / SSD / random)

- **Modal route** (cheap, portable): once `eigen` gives $\Phi,\omega_\alpha$, FRF / steady-state /
  random (input CSD → output PSD → RMS) and modal-superposition transient (exact piecewise-linear
  SDOF integration, Abaqus TG §2.5.3) are **post-processors with ~zero C++** — scriptable on top of
  `eigen` + `modalProperties`, then optionally a C++ `LadrunoModalResponse`.
- **Direct route** (deferred): direct SSD needs a complex linear solver $(K-\Omega^2M+i\Omega C)$.

---

## 5. Proposed ADR family + build order

| ADR | Scope | Reuses | Effort | Teacher |
|---|---|---|---|---|
| **A — Complex/state-space modal** | Modal-projection QZ on $\Phi^T\{M,C,K\}\Phi$; complex $\omega_d,\zeta$, complex modes | `eigen`, `modalProperties`, LAPACK `zggev` | **S** (quick win) | Abaqus §2.5.1 / LS-DYNA LCPACK=3 |
| **B — Prestressed modal + buckling** | `LadrunoBuckle`: two-tangent difference $K_\mathrm{ref}\phi=\lambda\Delta K\phi$ + tracking | corot/PDelta $K_g$, ARPACK `eigen` | **S–M** | Kratos `prebuckling_strategy` / Abaqus `*BUCKLE` |
| **C — Robust + parallel eigensolver** | **FEAST/PFEAST `EigenSolver`**: band-targeted, self-certifying, MPI via per-contour `MumpsParallelSOE`; unifies SP/MP; also complex contours | `MumpsParallelSOE`, PartitionedDomain | **L** (strategic) | Kratos PFEAST / LS-DYNA MPP factorization |
| **D — Frequency domain** | Modal FRF/SSD/random + modal-superposition transient (post-processors); direct SSD deferred | `eigen`, `modalProperties` | **M** | Abaqus §2.5.6–8 / LS-DYNA `*FREQUENCY_DOMAIN_*` |

**Recommended order: A → B → C → D.**
A validates the complex-modal direction cheaply and serially; B is high-demand plumbing; C is the
big strategic unifier that simultaneously fixes **parallel + band-targeting + Sturm robustness**
(and re-hosts A's complex case at scale); D layers cheaply on top.

**Cross-cutting robustness item** (fold into C, usable earlier): expose the **Sturm-sequence
inertia count** (negative pivots of the factorized $(K-\sigma M)$) from the Mumps/Band solvers so
any shift-invert path can *certify no modes were missed* — the guarantee Abaqus/LS-DYNA have and
ARPACK-OpenSees lacks.

---

## 6. Strategic importance — is this load-bearing?

Before committing PR effort, we asked whether this family is *load-bearing* (infrastructure others
build on) or a set of *leaf features* (end deliverables nothing depends on). The honest answer:
**the value is concentrated, not even.** Two ADRs are genuinely load-bearing; two are deliverables
that ride them.

### 6.1 Infrastructure vs leaf classification

| ADR | Role | Load-bearing? | What it unlocks *beyond its own ADR* |
|---|---|---|---|
| **43** FEAST + SP/MP fix | **Substrate** | **Yes — high** | The eigensolver every other modal capability rides; **and** the SP/MP parallel-composition fix is *general* parallel infrastructure (helps any large partitioned analysis, not just modal). |
| **40** Complex modal | **Domain-enabling** | **Yes — for our portfolio** | Correct modal damping for SSI/DRM/isolation/dampers (all non-classically damped); prerequisite for reduced-order modeling of those systems. |
| **42** Buckling / Kg | Standalone analysis | Modest | One real cross-link: the geometric-stiffness eigenpath could feed limit-point / bifurcation detection into the arc-length solvers (ADRs 20/22). |
| **44** Frequency domain | **Deliverable** | No | Consumes 46/43; nothing downstream builds on it. Pure end-feature (NVH, response spectrum, random). |

### 6.2 The dependency map

```
        ADR 43  (robust + parallel eigensolver, SP/MP fix)   ← foundation
            │     everything modal rides the eigensolver
   ┌────────┼─────────────────────┐
   ▼        ▼                     ▼
 ADR 46   ADR 42               ADR 44
 complex  buckling/Kg          FRF/SSD/random
   │        │                    ▲ (consumes 46/43; leaf)
   ▼        ▼
 unlocks   feeds limit-point detection
 beyond    into arc-length (ADR 20/22)
 modal ↓
 ROM / Craig-Bampton substructuring (candidate ADR 47)
```

### 6.3 The strongest "everything depends on it" argument

**Modal eigen sits upstream of every damped time-history run.** You cannot set Rayleigh damping
without modal frequencies, and today you cannot *trust* the frequencies of a large or partitioned
model (serial-only ARPACK, no completeness guarantee). So this is not a side feature — it sits under
NLTHA, SSI/DRM, and the apeSees/apeGmsh validation pipeline. That is the real load-bearing claim,
and it holds.

### 6.4 Honest counterpoints (do not oversell)

- **44 is a deliverable, not infrastructure** — build it last, or only when a project needs it.
- **42 is largely standalone** — valuable (stability is a missing analysis type) but a weak unlock-multiplier.
- **Modal analysis is more self-contained** than, say, a material kernel ten elements consume. The
  unlock multiplier comes from (a) the **parallel-composition fix in 43** and (b) **enabling ROM /
  substructuring + trustworthy damping calibration** — not from a sprawling dependency tree.

### 6.5 Revised center of gravity for sequencing

The earlier "A→B→C→D" order optimizes for *cheapest-first proof*. If the bar is *unlock the most*,
the weight shifts toward **43 (infrastructure)** and **46 (domain)**:

- **ADR 46 first** — cheap serial proof the complex-modal direction is correct (reuses existing
  `eigen`; validates against isolated/damped models). Low risk, directly serves the research.
- **ADR 43 is the strategic investment** — the eigensolver + parallel fix is what makes the whole
  family (and large-model NLTHA) usable; the SP/MP fix has value independent of modal analysis.
- **42 and 44 are opportunistic** — build when a specific project asks.

### 6.6 Forward note — the biggest genuine unlock (candidate ADR 47)

The single capability that would most strengthen the load-bearing case is **reduced-order modeling /
Craig–Bampton component-mode substructuring** (flagged as future in §2b of the Abaqus dossier). It
rides directly on this family (needs a trustworthy modal basis + parallel eigen) and is the enabler
for fast large SSI / real-time hybrid simulation. **Proposed as a candidate ADR 47 — not yet
written; pending decision.**

---

## 7. Execution plan

The phased rollout (P-A…P-G), cross-cutting decisions (assembled-`C`, MKL-FEAST vs PFEAST, `-shift`
exposure, vanilla-footprint policy), program gates, and risk register live in the umbrella program
ADR **[[45_ladruno_modal_family_roadmap_adr]]**. Headline sequence: **P-A** ADR 46 serial (cheap
proof) → **P-B** ADR 43 serial MKL-FEAST (substrate, zero new dep) → **P-C** ADR 43 parallel +
SP/MP unification → **P-D** ADR 42 buckling (opportunistic) → **P-E** ADR 43 complex contours
(re-host 46 at scale) → **P-F** ADR 44 frequency domain → *(P-G)* ADR 47 ROM.
