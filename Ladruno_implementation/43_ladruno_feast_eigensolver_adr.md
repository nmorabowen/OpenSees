---
title: "ADR 43 — FEAST/PFEAST band-targeted parallel eigensolver: design spec"
project: Ladruno
type: ADR / design spec
status: draft
priority: high
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"     # §3 — the parallel/MP answer; this ADR is its heart (track "C")
  - "[[modal_gap_study/01_opensees_current_state]]" # ground-truth file:line audit of the current parallel-eigen plumbing
  - "[[modal_gap_study/03_kratos_source]]" # PRIMARY template — bundled PFEAST, 3-level MPI, no Anasazi/Trilinos eigen
  - "[[modal_gap_study/04_lsdyna_theory]]" # MPP eigen = distributed factorization, not distributed Lanczos; Sturm/inertia certification
  - "[[46_ladruno_complex_modal_adr]]"     # complex/damped state-space modal — re-host its complex case at scale via complex contours
  - "[[42_ladruno_buckling_adr]]"          # linear buckling — benefits from band-targeting + Sturm certification
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
  - "[[LEDGER_quirks]]"
tags:
  - adr
  - solver
  - parallel
  - eigen
  - modal
  - feast
  - mpi
  - mumps
updated: 2026-06-22
---

# ADR 43 — `FeastEigenSOE` / `FeastEigenSolver` (band-targeted parallel eigensolver)

> **Strategic role (load-bearing assessment — see [[modal_gap_study/00_SYNTHESIS]] §6).**
> **The substrate — highest load-bearing of the family.** This is infrastructure, not a feature:
> (1) it is the eigensolver every other modal capability (46/42/44) rides; (2) the SP/MP
> parallel-composition fix is *general* parallel infrastructure that helps any large partitioned
> analysis, not just modal. Combined with the fact that modal eigen sits upstream of every damped
> time-history run (Rayleigh-damping calibration), this ADR is what makes the whole family — and
> large-model NLTHA — actually trustworthy at scale. **The strategic investment** (large build);
> sequence after the cheap ADR-46 proof.

**Status:** in progress — **P1 (serial MKL-FEAST, `eigen -feast`) MERGED
([#515](https://github.com/nmorabowen/OpenSees/pull/515)); P2 (`-certify` Sturm/inertia via
PARDISO on the SOE's own CSR + banner) MERGED
([#517](https://github.com/nmorabowen/OpenSees/pull/517)); gates G-A/G-B CLOSED. D2 spike
RUN 2026-07-07 — see §5.2/§9 R1: MKL FEAST confirmed as the P3 driver, no PFEAST vendoring;
P3 scoped into P3a/P3b/P3c below. P3a (sub-communicator plumbing: `setCommunicator` on
`MumpsParallelSolver` + `system Mumps -commSplit` hook) BUILT + gated same day —
4-rank/2-group concurrent-solve gate green (`feast_d2_spike/p3a_commsplit_gate.py`),
Opus adversarial gate clean (no blockers; WORLD/tag-0 channel-envelope subtlety
documented in [[LEDGER_quirks]], true envelope isolation deferred to P3c-MPI).
P3b (block-real complex-shift kernel `LadrunoBlockZKernel`, symmetric 2n×2n
PARDISO mtype −2 with reused analysis phase + native mtype 6 baseline, gated via
`eigen -feast -blockZGate`) MERGED ([#527](https://github.com/nmorabowen/OpenSees/pull/527)).
P3c-serial (`eigen -feast … -rci`: the `dfeast_srci` RCI orchestration with
`LadrunoBlockZKernel` as the inner contour solve — the seam the MPI rung
distributes; driver-parity gate λ ≤1e-8 rel + MAC ≥0.999 + `-certify`
composition; MKL RCI QUIRK: `dfeast_srci` SHRINKS its in/out `m0` in place
when the band holds fewer modes than the subspace — keep the caller's m0
by value or the saturation-enlargement compare plateaus and refuses)
MERGED ([#530](https://github.com/nmorabowen/OpenSees/pull/530)). P3c-MPI (L3-only
per R0): `eigen -feast … -rci` under openseesmp routes the RCI inner solve
through `LadrunoDistBlockZKernel` — distributed `dmumps` SYM=2 of the 2n block
system across the ranks (via a new `LadrunoFeastInnerSolve` seam + self-registering
factory) with the `dfeast_srci` outer loop replicated + solution broadcast for
lockstep; gate `feast_d2_spike/p3c_mpi_gate.py` green at mpiexec −n 2/4 (dist
spectrum == serial oracle to 3–5e-13, all ranks agree, ΦᵀMΦ==I); Opus gate +
PR pending.** classTags **33022**
(`FeastEigenSOE`) + **33023** (`FeastEigenSolver`) **ACTIVE in `SRC/classTags.h`**. P1 deviations
from this draft (deliberate, ledgered): the packed CSR driver `dfeast_scsrgv` instead of the §5.2
RCI seam — the RCI's complex shifted solves cannot route through a *real* inner `LinearSOE`, so
that seam moves to P3 alongside the D2 MPI spike; and `solve()` is `_WIN32`-guarded (Ubuntu Zone-A
CI links reference LAPACK, no MKL — root CMake's MKL branch defines no macro to gate on).
This is the strategic unifier of the modal-analysis effort: it simultaneously closes
the *distributed-eigen-at-scale*, the *SP/MP-don't-compose*, the *no-band-targeting*,
and the *no-completeness-guarantee* gaps identified in
[[modal_gap_study/00_SYNTHESIS|the synthesis]] (§3) — in one architecture.

> [!info] What this feature is, in one line
> A **contour-integration** eigensolver (FEAST) wrapped as an OpenSees
> `EigenSOE`/`EigenSolver` pair that answers *"give me **all** modes in the frequency
> band [f₁, f₂]"* by turning one distributed eigenproblem into **many independent
> distributed `(zM−K)⁻¹` linear solves** — each routed through the **existing**
> `MumpsParallelSOE` over an MPI sub-communicator — plus a tiny dense Rayleigh–Ritz.
> Band-targeted, self-certifying (subspace count vs Sturm/inertia), embarrassingly
> parallel, and a natural path to complex/damped modes at scale.

---

## 1. Driver & goal — OpenSees cannot do distributed modal at scale

The ground-truth audit ([[modal_gap_study/01_opensees_current_state|dossier 01]],
§C, file:line-pinned) establishes the parallel reality precisely:

1. **There is no distributed eigensolver.** The only Krylov driver is **serial
   ARPACK** (`dsaupd_`/`dseupd_`, reverse-communication Lanczos) in
   `SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp` — the RC loop runs on **one rank**
   (driver while-loop `ArpackSolver.cpp:221-291`, shift-invert `mode=3` at `:207`,
   `bmat='G'` at `:200`). `grep parpack|pdsaupd` returns nothing. No PARPACK, no
   SLEPc, no ScaLAPACK, no distributed Lanczos.

2. **The "parallel-eigen" path is serial-Lanczos-over-distributed-solves.** The
   RC plumbing *is* wired end-to-end into the partitioned machinery
   (`PartitionedDomain::eigenAnalysis`,
   `SRC/domain/domain/partitioned/PartitionedDomain.cpp:1179`; subdomains via
   `ShadowSubdomain::eigenAnalysis` / `ActorSubdomain::getNewEigenSOE`,
   `ActorSubdomain.cpp:749-755`), and the shift-invert `(K−σM)⁻¹` solve delegates to
   whatever `LinearSOE` is active (`ArpackSolver.cpp:258` for `ido=-1`, `:276` for
   `ido=1`). So if that `LinearSOE` is a distributed `MumpsParallelSOE`, the inner
   solve *is* a genuine parallel MUMPS solve — **but the Lanczos coordinate iteration
   `dsaupd_` still runs serially on the master**, with `ArpackSOE::checkSameInt`
   (`ArpackSOE.cpp:388-426`) merely keeping the `ido` flag synchronized across ranks.

3. **The two parallel models do not compose.** `MumpsParallelSOE` is created **only
   under `_PARALLEL_INTERPRETERS`** (`OpenSeesCommands.cpp:4079-4094` — the OpenSeesMP
   replicated-interpreter build), while the subdomain `ArpackSOE` distribution path is
   `_PARALLEL_PROCESSING` (OpenSeesSP, `SRC/tcl/mpiMain.cpp`). You **cannot** get
   *partitioned domain + parallel `(K−σM)⁻¹` + parallel coordinate iteration* in one
   build. SP gives distributed assembly but per-subdomain serial solves under a serial
   Lanczos; MP gives a distributed MUMPS solve but every rank redundantly re-runs the
   whole serial Lanczos. Neither is a true distributed eigensolver.

4. **No band targeting.** The `eigen` command takes `numEigen` near a single `shift`
   (and the xara parser hard-zeros `shift` at `analysis.cpp:271`). There is no
   *"all modes in [f₁, f₂]"* query — the seismic/NVH question — and no way to know you
   got them all.

5. **No completeness guarantee.** ARPACK-OpenSees has **no Sturm/inertia interval
   certification**. Both reference codes do: LS-DYNA's BCSLIB-EXT inertia count
   ([[modal_gap_study/04_lsdyna_theory|dossier 04]] §1) and FEAST's self-certifying
   subspace count. The user must request `numEigen` and trust the result.

6. **Serial-only post-processing.** `DomainModalProperties` /
   `ResponseSpectrumAnalysis` have zero MPI awareness; under SP they see only the
   master domain → wrong participation/effective-mass.

**Goal.** A distributed, band-targeted, self-certifying generalized eigensolver that
**reuses the distributed linear-solve infrastructure OpenSees already has** rather
than adding a distributed Krylov basis — exactly the architectural bet both Kratos and
LS-DYNA made.

---

## 2. Decision summary

**Adopt FEAST (contour integration), not PARPACK (distributed Krylov).**

FEAST recasts the generalized eigenproblem $K\phi=\lambda M\phi$ as the action of a
**spectral projector built from a contour integral of the resolvent**
$(zM-K)^{-1}$. Numerically, the contour integral becomes a small sum over
quadrature points, and **each quadrature point is an independent linear solve with a
different complex shift $z_j$.** This is the decisive property:

- **the dominant cost is the linear solves** — precisely what OpenSees already does
  well in parallel via MUMPS (`MumpsParallelSOE`/`MumpsParallelSolver`,
  `SRC/system_of_eqn/linearSOE/mumps/`);
- **there is no global Krylov synchronization** — no distributed orthogonalization,
  no parallel Lanczos restart, none of the invasiveness that a PARPACK/Anasazi block
  basis demands.

**Both reference codes deliberately avoided distributed Krylov.**
- **Kratos** has *no* Trilinos/Anasazi/SLEPc eigensolver at all (verified by
  `gh api search/code`, [[modal_gap_study/03_kratos_source|dossier 03]] §6); its only
  distributed-eigen path is **bundled PFEAST** with three-level MPI parallelism.
- **LS-DYNA MPP** has *no* special distributed Lanczos; parallelism lives **entirely
  in the distributed sparse factorization** of the shift-invert solves
  (`LSOLVR=7` multifrontal / `LSOLVR=30` = **MUMPS**, ParMETIS ordering;
  [[modal_gap_study/04_lsdyna_theory|dossier 04]] §4).

Two codes, two independent teams, the same verdict: **for a code that already has a
good distributed *linear* solver, reduce distributed eigen to many independent
distributed linear solves.** PARPACK is the lower-friction *port* (it reuses the
existing RC loop nearly verbatim) but it is the architecture both teachers rejected.
This ADR follows the teachers.

### The build decision: MKL FEAST vs vendored PFEAST

**The OpenSees build already links Intel MKL** (oneAPI; per the build journal and
`project_build_env_*` notes). **MKL's Extended Eigensolver IS FEAST** — the same
Polizzi contour-integration kernel, exposed as `?feast_*` (e.g. `dfeast_scsrgv`,
`zfeast_gcsrgv`) and the matrix-free RCI forms (`?feast_srci`/`?feast_grci`).
This is the central integration question of the ADR:

| Option | What | Pros | Cons |
|---|---|---|---|
| **A — MKL FEAST** (serial + RCI) | call the FEAST kernels bundled in MKL we already link | **zero new external dep**; threaded; the RCI form lets us supply OpenSees's *own* distributed solve as the inner solve | MKL's *distributed* FEAST is via **MKL Cluster** (CPardiso), **not** PFEAST's clean MPI-subgroup model; binding the inner solve to `MumpsParallelSOE` still needs the **RCI** path; licensing/redistribution of MKL parallel layer |
| **B — vendored PFEAST** | drop the BSD-ish FEAST/PFEAST 4.0 Fortran into `OTHER/FEAST/` (exactly like `OTHER/ARPACK`, `OTHER/LAPACK`, `OTHER/MUMPS`) | the **three-level MPI** model is first-class and documented; the RCI lets us route the inner solve to `MumpsParallelSOE`; no MKL-parallel-layer licensing entanglement; matches the Kratos template 1:1 | a new vendored Fortran lib + CMake target; PFEAST's MPI-comm bookkeeping to wire; build-matrix surface |

**Recommended split (de-risks the decision):**
- **P1–P2 (serial): Option A.** Use MKL FEAST directly — it is already on the link
  line, so the serial band-targeted solver and the Sturm/inertia certification cost
  **zero new dependency**. Validate the whole contour/Rayleigh–Ritz machine serially
  against ARPACK first.
- **P3+ (distributed): re-evaluate A-RCI vs B at the gate.** The deciding factor is
  whether MKL's RCI FEAST (`?feast_srci`) cleanly lets each quadrature point's solve
  be an *OpenSees* `MumpsParallelSOE` solve on an MPI sub-communicator. If yes → stay
  on A (no vendoring). If the MKL parallel layer fights the sub-communicator model →
  vendor PFEAST (Option B), which is purpose-built for it. **This is the key open
  decision (see §9, R1).**

**Net decision:** *Build the FEAST EigenSOE/Solver pair; start on MKL FEAST (no new
dep); keep the inner solve behind the RCI seam so the distributed backend
(MKL-cluster vs vendored PFEAST) is swappable at P3 without touching the OpenSees-side
architecture.*

---

## 3. Scope-fence & classTags

`#define EigenSOE_TAGS_FeastEigenSOE 33022`  *(EigenSOE tag namespace; current max is
`SparsePythonCOOEigenSOE`=9 — 33022 sits in the Ladruno private ≥33000 band, mirroring
the per-namespace convention used for `INTEGRATOR_TAGS_ExplicitDifferenceStatic`=33001,
`HANDLER_TAG`/`RECORDER_TAGS` 33001…)*

`#define EigenSOLVER_TAGS_FeastEigenSolver 33023`  *(EigenSOLVER tag namespace; current
max is `GeneralArpackSolver`=6)*

Both **RESERVED, not yet built** — record in `LEDGER_implementations.md` at design time
so a sibling ADR cannot collide.

**In scope (this ADR):**
- **Real generalized band-targeted eigen:** $K\phi=\lambda M\phi$, $K$ symmetric,
  $M$ SPD, *all* modes with $\lambda\in[\lambda_1,\lambda_2]=[(2\pi f_1)^2,(2\pi
  f_2)^2]$, via a real contour (interval) FEAST.
- **Serial first** (MKL FEAST, P1), **then MPI** via per-contour-point
  `MumpsParallelSOE` inner solves (P3).
- **Sturm/inertia certification** (P2): expose the negative-pivot count of the
  factorized $(K-\sigma M)$ from the Mumps/Band solvers so any shift-invert path can
  certify no modes were missed — and use it to cross-check FEAST's subspace count.
- **Unify the SP/MP build gating** (P4): a single coherent
  partitioned-domain + distributed-solve + contour-parallel model in one build.

**NOT in scope (handed to siblings / later):**
- **Full complex-contour damped eigen** — re-host of [[46_ladruno_complex_modal_adr|ADR
  46]]'s state-space complex modal at scale. FEAST's complex contours
  (`zfeast_*`) are the *mechanism*; the quadratic-pencil linearization is ADR 46's.
  P5 coordinates with ADR 46; not built here.
- **AMLS / component-mode synthesis** (the EIGMTH=101-class "thousands of NVH modes"
  route) — a separate future ADR.
- **Parallel `modalProperties` / CQC-SRSS** — needed for a complete parallel modal
  offering but tracked separately ([[modal_gap_study/00_SYNTHESIS|synthesis]] gap E.3).

---

## 4. Formulation

### 4.1 The contour-integral spectral projector

For the generalized pencil $(K, M)$ with $M$ SPD, the eigenvalues $\lambda_i$ are real.
Let $\Gamma$ be a positively-oriented closed contour in the complex plane enclosing
exactly the target interval $[\lambda_1, \lambda_2]$ on the real axis. The **spectral
projector** onto the invariant subspace of the enclosed eigenvalues is the Cauchy
contour integral of the resolvent:

$$
P \;=\; \frac{1}{2\pi i}\oint_{\Gamma} (zM - K)^{-1}\,M\;dz .
$$

$P$ is the $M$-orthogonal projector onto $\mathrm{span}\{\phi_i : \lambda_i \in
[\lambda_1,\lambda_2]\}$: $P\phi_i=\phi_i$ for enclosed modes and $P\phi_j=0$
otherwise. **All band physics is in the contour** — the band *is* the contour.

### 4.2 Numerical quadrature on an ellipse

Apply $P$ to a tall random matrix $Y\in\mathbb{R}^{n\times m_0}$ (with $m_0 \ge$ the
number of modes expected in the band) to obtain a basis $Q = PY$ for the subspace.
$P$ is never formed; the integral is evaluated by $N_q$-point quadrature on an ellipse
over $[\lambda_1,\lambda_2]$ (FEAST default: Gauss–Legendre, $N_q\!\approx\!8\!-\!16$;
trapezoidal on a circle is the alternative). With nodes $z_j$ and weights $\omega_j$
on the contour (the ellipse has center $(\lambda_1{+}\lambda_2)/2$, real semi-axis
$(\lambda_2{-}\lambda_1)/2$, and an aspect ratio $r$ flattening it onto the real axis):

$$
Q \;\approx\; \sum_{j=1}^{N_q} \omega_j\,(z_j M - K)^{-1} M\,Y .
$$

**Each term is one linear solve with a complex shift $z_j$:**

$$
(z_j M - K)\,X_j \;=\; M\,Y, \qquad Q \;=\; \sum_{j=1}^{N_q}\omega_j X_j .
$$

The solves for different $j$ share *nothing* but the right-hand side $MY$ — different
factorizations, no coupling. **→ embarrassingly parallel over the quadrature points.**

### 4.3 Rayleigh–Ritz reduction

Project the full pencil onto the $m_0$-dimensional subspace spanned by $Q$:

$$
A_Q = Q^{\!\top} K\,Q \in \mathbb{R}^{m_0\times m_0}, \qquad
B_Q = Q^{\!\top} M\,Q \in \mathbb{R}^{m_0\times m_0},
$$

and solve the **tiny dense** generalized eigenproblem
$A_Q\,W = B_Q\,W\,\Lambda$ with LAPACK (`dsygv`/`dggev`) on the reduced $m_0\times m_0$
matrices. Ritz pairs whose $\lambda$ falls inside $[\lambda_1,\lambda_2]$ are accepted;
the Ritz vectors lift back as $\Phi = Q W$. Residuals
$\lVert K\phi_i - \lambda_i M\phi_i\rVert$ drive a **FEAST outer refinement loop**
(re-apply $P$ to the current $\Phi$, typically 2–4 iterations) until convergence — *not*
hundreds of Lanczos restarts.

### 4.4 Subspace-size self-certification

FEAST is self-certifying: it returns $m$, the count of accepted eigenvalues actually
found in the band, and a stochastic estimate of the true count. The contract:

- **$m_0$ must over-estimate** the modes in the band. If $m < m_0$ comfortably →
  the band is fully resolved and **provably complete** (the subspace was large enough
  to capture every enclosed mode).
- **If $m$ saturates at $m_0$** → the band held more modes than the subspace could
  hold; **enlarge $m_0$ and re-run**. (Practical estimate, mirroring Kratos:
  $m_0 = \lceil 1.5\,m_{\text{guess}}\rceil$ for an interior band,
  $2\,m_{\text{guess}}$ for an extremal band.)

This closes **both** robustness gaps at once: band-targeting *and* the
"no-missed-modes" guarantee. **Cross-check with Sturm/inertia (P2):** factor
$(K-\sigma M)$ at $\sigma=\lambda_1$ and $\sigma=\lambda_2$; the difference in
negative-pivot counts (matrix inertia) is the *exact* number of eigenvalues in the
band. FEAST's $m$ must equal it — an independent certificate from the factorization
itself (the LS-DYNA BCSLIB-EXT guarantee, [[modal_gap_study/04_lsdyna_theory|dossier
04]] §1).

### 4.5 Three-level parallelism map (the heart of the ADR)

Mirroring PFEAST ([[modal_gap_study/03_kratos_source|dossier 03]] §6b), the MPI ranks
split into a three-level hierarchy:

| Level | What is parallelized | OpenSees mapping |
|---|---|---|
| **L1 — interval / contour** | split a wide $[f_1,f_2]$ into sub-bands, one contour each; embarrassingly parallel across sub-bands | one MPI sub-communicator group per sub-band (optional; v1 = one band) |
| **L2 — quadrature point** | the $N_q$ independent solves $(z_jM-K)X_j=MY$; **the natural FEAST parallelism** — sum over independent solves, no global Krylov sync | one MPI sub-communicator per quadrature point $z_j$ |
| **L3 — inner distributed solve** | the single distributed sparse solve $(z_jM-K)X_j=MY$ across the ranks of *its* sub-comm | **the existing `MumpsParallelSOE`** on that sub-comm (the entire reuse) |

Then a **cheap dense Rayleigh–Ritz** (§4.3) on the $m_0$-dim reduced pencil — small
enough to replicate or do on one rank, then broadcast. Why each quadrature solve is
independent → **embarrassingly parallel**: the projector is a *sum* (§4.2); term $j$
needs only $z_j$ and the shared $MY$; there is no data dependence between terms and no
iterative coupling — the polar opposite of a distributed Lanczos basis that must
re-orthogonalize against all prior vectors every step.

---

## 5. Architecture & API

### 5.1 Where the classes slot in

OpenSees already mirrors the Kratos Scheme/B&S/Strategy layering
([[modal_gap_study/03_kratos_source|dossier 03]] §0): `EigenAnalysis` /
`EigenIntegrator` / `EigenSOE` / `EigenSolver`
(`SRC/analysis/...`, `SRC/system_of_eqn/eigenSOE/`). The new pair drops into that
hierarchy beside `ArpackSOE`/`ArpackSolver`:

```
EigenSOE  (base, SRC/system_of_eqn/eigenSOE/EigenSOE.h)
  ├── ArpackSOE        (33022-sibling slot; wraps a LinearSOE for (K-σM)⁻¹)
  └── FeastEigenSOE    (33022)  ── holds K, M handles + contour params + m0
EigenSolver (base, SRC/system_of_eqn/eigenSOE/EigenSolver.h)
  ├── ArpackSolver     (serial RC Lanczos)
  └── FeastEigenSolver (33023) ── drives the contour quadrature + Rayleigh-Ritz
```

`FeastEigenSOE`, like `ArpackSOE`, **wraps a `LinearSOE* theInnerSOE`** for the
resolvent solves (`ArpackSOE.h:77`, `ArpackSOE::setLinearSOE` at
`ArpackSOE.cpp:382-386` is the exact precedent). `StaticAnalysis::setLinearSOE` /
`setEigenSOE` already auto-propagate the active `LinearSOE` into the EigenSOE
(`StaticAnalysis.cpp:549, 577`) — so wiring a distributed inner solve is *free*: set
`system Mumps` and the contour solves become parallel MUMPS solves, identical to how
ArpackSOE inherits its distributed inner solve today.

### 5.2 Inner-solve reuse — the RCI seam

`FeastEigenSolver::solve()` runs the **FEAST RCI** form (`?feast_srci`): FEAST hands
back a complex shift $z_j$ and a right-hand side and asks the caller to solve
$(z_jM-K)X=B$. The caller assembles the shifted matrix into the inner SOE and calls
`theInnerSOE->solve()` — **the same delegation ArpackSolver does** at
`ArpackSolver.cpp:258,276`. The difference from the ARPACK RC loop is only the *flag
protocol*: FEAST's RCI `ijob` replaces ARPACK's `ido`; the body is the same
"assemble → `theInnerSOE->solve()` → hand X back" pattern.

> [!warning] **CORRECTED 2026-07-07 (D2 spike, [[_modal_family_handoff]]).** This
> paragraph previously claimed "for the real-symmetric band case FEAST's `srci` keeps
> the real-arithmetic structure (it pairs $z_j$ and $\bar z_j$), so a real distributed
> solver (MUMPS) can be used per conjugate pair" — **refuted by an adversarial
> literature check (Opus, cross-referencing the FEAST v3/v4 User Guides and the MKL
> `?feast_srci` reference)**. At `ijob=10/11` FEAST hands the caller a **genuinely
> complex** workspace (`zwork`) and asks for a complex factorization/solve of
> $A_z=z_jM-K$ — every individual contour solve is complex, full stop. The Hermitian
> conjugate-pairing (`fpm(2)`) only **halves the count of solves** (skip $\bar z_j$,
> conjugate the result) — it does not make any *one* solve real. `dfeast_scsrgv`'s
> convenience form doesn't dodge this either: the observable behavior is consistent
> with MKL performing a **complex PARDISO** factorization (`mtype=6`,
> complex-symmetric) per contour point internally — an inference about undocumented
> MKL internals, but the documented RCI contract (complex `zwork`, complex $A_z$
> factorization requested from the caller) is sufficient on its own. There is no
> real-only-$n\times n$ path in FEAST, reference or vendor. See §9 R1 for the
> corrected decision.

### 5.3 MPI communicator splitting

At P3, `FeastEigenSolver` splits the run's communicator: `MPI_Comm_split` carves
`MPI_COMM_WORLD` into $N_q$ sub-communicators (one per quadrature point, L2); each
sub-comm hosts one `MumpsParallelSOE` (L3) that factorizes and solves $(z_jM-K)$ across
its ranks. Results are `MPI_Allreduce`-summed into $Q$ (the projector sum, §4.2). The
reduced $m_0\times m_0$ Rayleigh–Ritz is replicated/broadcast. This is the only genuinely
new MPI bookkeeping; everything below it (the distributed factor/solve) is the proven
`Mumps*` stack.

### 5.4 Command

```
eigen -feast $fmin $fmax  [-m0 $m0] [-nq $Nq] [-tol $tol] [-maxiter $n] \
                          [-contour ellipse|circle] [-aspect $r] [-certify]
```

- `$fmin $fmax` — frequency band (Hz); internally $\lambda_{1,2}=(2\pi f)^2$.
- `-m0` — subspace size (default heuristic §4.4; auto-enlarge on saturation).
- `-nq` — quadrature points (default 8).
- `-certify` — also run the Sturm/inertia count and assert it equals FEAST's $m$.
- Returns the eigenvalues (as today, `domain->getEigenvalues()`); eigenvectors pushed
  to nodes via `theAnalysisModel->setEigenvector` (`StaticAnalysis.cpp:339-345`,
  unchanged) so existing recorders/`modalProperties` consume them verbatim.

Parsers: the band path is **added** beside the existing flag map in **both**
`SRC/runtime/commands/analysis/analysis.cpp:250` (xara) and
`SRC/interpreter/OpenSeesCommands.cpp:270` (legacy/parallel) — selecting
`EigenSOE_TAGS_FeastEigenSOE` and stashing `fmin/fmax/m0` on the SOE.

### 5.5 How this unifies the SP/MP split (the concrete fix)

The incoherence ([[modal_gap_study/01_opensees_current_state|dossier 01]] §C.4) is:
`MumpsParallelSOE` is `_PARALLEL_INTERPRETERS`-only (`OpenSeesCommands.cpp:4079`) while
the subdomain `ArpackSOE` distribution is `_PARALLEL_PROCESSING`. FEAST dissolves it
because **its only distributed need is the inner linear solve** — there is no
distributed *eigen* coordinate iteration to gate. Concretely:

- The contour-parallel orchestration (L1/L2) is **eigensolver-internal MPI**
  (`MPI_Comm_split`), independent of which OpenSees parallel-build macro is set.
- The inner solve (L3) is `MumpsParallelSOE`, which already works under MP.
- So a **single build** that exposes `MumpsParallelSOE` + the FEAST orchestrator gives
  *partitioned-or-replicated assembly + distributed solve + contour-parallel eigen* —
  the combination that is impossible today. P4 is the build-flag surgery that makes
  `MumpsParallelSOE` reachable alongside the contour orchestrator in one coherent
  configuration (candidate: a `_PARALLEL_FEAST`/unified guard that compiles the Mumps
  parallel SOE + the comm-splitting without requiring the full replicated-interpreter
  model). Because the eigen coordinate work is *not* distributed, the choice of
  SP-style partition vs MP-style replication becomes orthogonal to the eigensolve.

---

## 6. OpenSees integration points (file:line)

| Concern | File:line | Action |
|---|---|---|
| **The RC loop to mirror** | `SRC/system_of_eqn/eigenSOE/ArpackSolver.cpp:221-291` (`ido` while-loop; `ido=-1`→solve `:247-262`, `ido=1`→solve `:264-282`) | `FeastEigenSolver::solve()` replicates this as a FEAST-RCI `ijob` loop: same "assemble shifted system → `theInnerSOE->solve()` → return X" body |
| **The `(K−σM)⁻¹` delegation** | `ArpackSolver.cpp:258` (`ido=-1`), `:276` (`ido=1`); SOE wrap at `ArpackSOE.h:77`, `ArpackSOE::setLinearSOE` `ArpackSOE.cpp:382-386` | `FeastEigenSOE` wraps the same `LinearSOE*`; contour solve = assemble $(z_jM-K)$ + `theInnerSOE->solve()` |
| **MumpsParallelSOE reuse** | `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.{h,cpp}`, `MumpsParallelSolver.{h,cpp}`; constructed at `OpenSeesCommands.cpp:4079-4094` | each L2 sub-comm hosts one of these as the inner SOE (no change to Mumps itself) |
| **`_PARALLEL_PROCESSING` vs `_PARALLEL_INTERPRETERS` gating** | `OpenSeesCommands.cpp:4079` (`#ifdef _PARALLEL_INTERPRETERS` around Mumps**Parallel**SOE); SP main `SRC/tcl/mpiMain.cpp:184-291`; MP main `SRC/tcl/mpiParameterMain.cpp:222-349` | P4: introduce a coherent guard so `MumpsParallelSOE` + contour comm-split compile together (LEDGER_vanilla_files row) |
| **Command registration** | xara parser `SRC/runtime/commands/analysis/analysis.cpp:250` (flag map `:277-308`, build at `:326-328`); legacy/parallel `SRC/interpreter/OpenSeesCommands.cpp:270`; SOE construction `SRC/runtime/runtime/BasicAnalysisBuilder.cpp:871-910` | add `-feast` branch in both parsers + `newEigenAnalysis` case building `FeastEigenSOE` and wiring inner SOE (`setLinearSOE`) |
| **EigenSOE registration for parallel send** | `FEM_ObjectBrokerAllClasses::getNewEigenSOE` `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp:3378-3396` (switch currently only `ArpackSOE`); also `TclPackageClassBroker.cpp:2027`, `FEM_ObjectBroker.cpp:366` | add `case EigenSOE_TAGS_FeastEigenSOE: return new FeastEigenSOE();` so `ActorSubdomain::getNewEigenSOE` (`ActorSubdomain.cpp:751`) can reconstruct it remotely |
| **sendSelf/recvSelf** | precedent `ArpackSOE::sendSelf/recvSelf` `ArpackSOE.cpp:269-372` (sends only a `processID`/channel-handshake ID — **not** the matrices; `checkSameInt` `:388-426` syncs the RC flag) | `FeastEigenSOE::sendSelf/recvSelf` follows the same contract **plus** marshals the scalar contour params (`fmin,fmax,m0,Nq,tol`); matrices stay distributed, assembled locally per rank |
| **Assembly (K into A, M into M)** | `StaticAnalysis::eigen` `SRC/analysis/analysis/StaticAnalysis.cpp:244` (`zeroA/zeroM` `:275-276`, K via `addKtToTang`+`addA` `:285-293`, M via `addMtoTang`+`addM` `:300-323`, `solve(numMode,gen,findSmallest)` `:330`); `EigenIntegrator::formEleTangK` `EigenIntegrator.cpp:191-196` | **unchanged** — FEAST consumes the same assembled $K,M$; only `FeastEigenSOE::solve` differs. (Note: K is always `Kt` here — band-of-frequencies, not buckling; the buckling $K_g$-into-M variant is ADR 42) |
| **Vendoring precedent (if Option B)** | `OTHER/ARPACK/` (Fortran, own `CMakeLists.txt`), `OTHER/LAPACK/`, `OTHER/MUMPS/` | `OTHER/FEAST/` added the same way + a `USE_FEAST`/`USE_PFEAST` CMake gate |

---

## 7. Phased roadmap + gates

| Phase | Deliverable | Gate |
|---|---|---|
| **P1 — serial MKL-FEAST EigenSolver** | `FeastEigenSOE`/`FeastEigenSolver` (33022/33023) using MKL `dfeast_scsrgv` (or `dfeast_srci` RCI through the existing inner `LinearSOE`); `eigen -feast fmin fmax`; band → ellipse contour → quadrature → Rayleigh–Ritz | **same eigenpairs as ARPACK** on a medium serial model within $\le10^{-6}$ on $\lambda$ and MAC $\ge0.999$ on $\phi$ for every mode in the band (§8.1) |
| **P2 — Sturm/inertia certification** | expose negative-pivot/inertia count from the Band/Mumps factorization; `-certify` asserts FEAST $m$ == inertia($\lambda_2$)−inertia($\lambda_1$) | hand-built model with **known closely-spaced modes straddling the band edge**: FEAST $m$ equals the inertia count; deliberately under-sized $m_0$ is **detected** (saturation flagged), not silently wrong (§8.2) |
| **P3 — MPI per-contour parallel** | Three sub-phases (D2-driven, §9 R1): **P3a** fix `MumpsParallelSolver`/`SOE` to honor a passed sub-communicator instead of hardcoded `MPI_COMM_WORLD`; **P3b** new symmetric $2n\times2n$ block-real inner SOE for the complex contour solve (LDLᵀ via existing `dmumps`, `SYM=2`); **P3c** `FeastEigenSolver` orchestration, two rungs — **P3c-serial** (`-rci`): the `dfeast_srci` RCI ijob loop (10=factorize z_jM−K, 11=solve, 30/40 matvecs) with ONE `LadrunoBlockZKernel` per solve (setShiftBlock per contour node, PARDISO analysis phase reused across nodes/refinement-loops/saturation-retries); **P3c-MPI (L3-only, per R0 panel decision)**: a distributed 2n×2n SYM=2 block-real `dmumps` solve on ONE sub-communicator (P3a plumbing), wired as the `dfeast_srci` inner solve; the validated `runFeastRci` outer loop runs REPLICATED across the sub-comm; solution broadcast to all ranks (lockstep RCI); symbolic analysis reused across shifts; + the scoped envelope-isolation fix. **L2 deferred** (§R0): `MPI_Comm_split` per-quadrature sub-comms + `Allreduce`-projector is a later multiplier, not this rung | P3a: 2+ sub-comms solve independent systems concurrently, no cross-talk; P3b: block-real solve matches a direct complex solve (numpy oracle) to solver precision **including a contour aspect-ratio sweep asserting residuals as Im(z)→0** (equivalent-real conditioning degrades near the real axis), plus a measured cost comparison vs zmumps; P3c-serial: `-rci` reproduces the `dfeast_scsrgv` driver-path eigenpairs (λ ≤1e-8 rel, MAC ≥0.999) on the frame battery + `-certify` composition; P3c-MPI: `mpiexec -n 2/4` distributed run == serial `-rci` spectrum bit-comparable (differential test vs the shipped oracle) + MUMPS strong-scaling instrumented on a representative mesh (§8.3) |
| **P4 — unify SP/MP build gating** | one coherent build where `MumpsParallelSOE` + the contour orchestrator compile together; comm-split orthogonal to partition/replication | the *impossible-today* combination — partitioned/replicated assembly + distributed solve + contour-parallel eigen — runs green in **one** binary (`OpenSeesSP` and/or `OpenSeesMP`); LEDGER_vanilla_files documents the guard change |
| **P5 — complex contours** *(coordinate with [[46_ladruno_complex_modal_adr|ADR 46]])* | `zfeast_*` complex contour for the damped/non-symmetric pencil; needs a complex inner SOE | re-host ADR 46's complex/state-space modal at scale; gated on ADR 46's quadratic-pencil linearization landing first |

Adversarial-gate policy (per [[feedback_adversarial_gate_when]]): **P1–P2 warrant the
full multi-agent gate** (novel-to-the-fork contour math + a completeness *claim*);
**P3–P4 are core/parallel-build edits → also gated** (touching the parallel mains is
exactly the "core-or-vanilla edit" trigger).

---

## 8. Validation / oracle battery

### 8.1 Equivalence to ARPACK (P1)
On a medium model (e.g. a meshed 3D frame, $n\sim10^4$–$10^5$ DOF), run
`eigen N` (ARPACK, lowest $N$) and `eigen -feast 0 fmax` with $f_\text{max}$ chosen to
enclose exactly those $N$ modes. **Assert:** same set of $\lambda_i$ ($\le10^{-6}$
rel.), MAC $\ge0.999$ mode-by-mode, $\Phi^\top M\Phi=I$ (mass-normalization).

### 8.2 Band-targeting correctness (P1–P2)
Compute the **full spectrum** with `fullGenLapack` (`FullGenEigenSOE`, dense
`dggev`) on a small model. Pick several bands $[f_1,f_2]$ — interior, edge-straddling,
empty, and one with a degenerate/closely-spaced pair on the boundary. **Assert:**
FEAST returns *exactly* the modes whose $f\in[f_1,f_2]$, no more, no fewer; the
$m_0$-saturation case is *flagged and auto-enlarged*, never silently truncated.

### 8.3 Parallel scaling (P3)
Partitioned model under MPI. **Strong scaling:** fixed model, increasing ranks
mapped to quadrature sub-comms — wall-time drop tracks the embarrassingly-parallel
quadrature (near-linear up to $N_q$ groups, then L3 takes over). **Weak scaling:**
grow model + ranks together. **Assert:** eigenpairs match the P1 serial run; report
the strong/weak curves in the bundle.

### 8.4 Completeness via inertia (P2)
For each band in 8.2, independently count eigenvalues via the Sturm/inertia difference
of the factorized $(K-\sigma M)$ at the two band edges. **Assert:** equals FEAST's $m$.
This is the *guarantee* ARPACK-OpenSees lacks and both teachers have.

**Oracle bundle:** `Ladruno_implementation/feast_eigensolver_validation/` (README +
figures/ + code/), mirroring the `lemaitre_validation/` and concrete validation
precedents.

---

## 9. Risk register

> [!important] **R0 — DECIDED 2026-07-08 (design review panel): the P3c-MPI rung is
> L3-only (distributed inner solve on the validated `dfeast_srci` outer loop); true
> L2 quadrature parallelism is a deferred follow-on, not this rung.**
> After P3c-serial shipped (#530), building the MPI rung surfaced a hard fact: MKL's
> `dfeast_srci` visits the $N_q$ contour nodes **sequentially inside its RCI loop and
> accumulates the projector internally** (confirmed against Intel's canonical
> `dfeast_sparse_rci.c`). There is **no seam** to hand different quadrature nodes to
> different communicators — so the §5.3 "`MPI_Comm_split` per quadrature point,
> `Allreduce` the projector" (L2) is achievable **only by abandoning `dfeast_srci`
> and re-implementing FEAST's entire contour outer loop** (Gauss–Legendre ellipse
> nodes/weights + Jacobian, projector, Rayleigh–Ritz, refinement, trace convergence,
> $m_0$ management, spurious-mode filtering). A three-seat expert panel (FEAST/PFEAST
> authority "PET-CE", parallel-HPC/MUMPS-scaling, OpenSees-fork architecture) was
> **unanimous for L3-only** as this rung. Load-bearing reasons:
> 1. **L3 is the unavoidable foundation for BOTH.** Pure-L2 caps at $N_q$ ranks
>    ($\le 8$, and $\sim\times4$ after half-contour conjugate symmetry); to spend
>    more ranks L2 must nest L3 *inside* each sub-comm — which is exactly the
>    distributed-MUMPS work L3-only builds. So L3 is on the critical path regardless;
>    L2 is a later *multiplier* on top (the PFEAST $L2\times L3$ nest), and P3a's
>    `setCommunicator(MPI_Comm)` already parameterizes L3 on an arbitrary sub-comm so
>    the composition is not foreclosed.
> 2. **L3 owns the target regime.** For $n\sim10^5$–$10^6$ structural models the
>    dominant cost is one big sparse symmetric factorization (the PFEAST paper
>    foregrounds the *distributed linear solver*, not L2); the 2n block augmentation
>    (P3b) doubles the factored dimension, pushing harder toward L3. L2-without-L3 is
>    useless at $n\sim10^6$ — a single node's LDLᵀ factor won't fit in one rank's RAM.
>    L3's scaling ceiling *grows with n* and covers the 16–64-rank budget; L2's is
>    pinned at $\sim\times4$.
> 3. **Memory:** L3-only holds ONE distributed factorization at a time (sequential
>    nodes, reused buffers); L2 holds $N_q$ simultaneous 2n factorizations spread
>    across the machine — the difference between "runs" and "OOM" for a $2n=2\times10^6$
>    3D factor.
> 4. **Gate & reuse:** L3-only reuses the shipped `runFeastRci` **verbatim** and its
>    correctness gate is a **bit-for-bit differential test against the serial `-rci`
>    oracle** (same MKL outer loop, same arithmetic; only the inner solve swaps
>    PARDISO→distributed MUMPS, whose residual FEAST's refinement absorbs). L2 would
>    discard `runFeastRci` and demand a full from-scratch numerical qualification of
>    a novel contour engine — the fork's largest novel-math surface. B's one genuine
>    edge (dropping MKL FEAST would let the distributed path run on the Linux MP
>    build too) buys only a *secondary* target: the blessed MP build is Windows/oneAPI
>    with MKL already linked (verified — `dist/openseesmp/` ships `mkl_core` +
>    `impi.dll`, so `dfeast_srci` is callable there).
> **P3c-MPI (L3) scope:** a distributed block-real complex-shift solve (2n×2n SYM=2
> `dmumps` on a sub-communicator via the P3a plumbing, distributed triplet input,
> solution broadcast to all ranks so the replicated RCI stays in lockstep, symbolic
> analysis reused across shifts) wired as the `dfeast_srci` inner solve; the
> replicated RCI runs on one sub-comm; plus the scoped envelope-isolation fix and an
> `mpiexec -n 2/4` differential test vs the serial `-rci` spectrum. **L2 (deferred):**
> `MPI_Comm_split` the world into $N_q$ node-groups, run the L3 kernel inside each,
> `Allreduce` the projector — pursued only if a profiled workload shows a $>\times4$
> quadrature win L3 can't deliver.

> [!question]- **R1 — RESOLVED 2026-07-07 (D2 spike): stay on MKL FEAST + build two
> OpenSees-side inner-solve fixes; do not vendor PFEAST.**
> Original question: does MKL's **RCI** FEAST (`?feast_srci`) cleanly let each
> quadrature solve be an *OpenSees* `MumpsParallelSOE` solve on an arbitrary MPI
> sub-communicator? **D2 spike findings** (code read + numerical check + Opus
> literature cross-check, full record in [[_modal_family_handoff]] and the D2 bundle):
>
> 1. **`MumpsParallelSolver` hardcodes `MPI_COMM_WORLD`, not an arbitrary sub-comm.**
>    The `mpi_comm`-taking constructor (`MumpsParallelSolver.cpp:54-64`) accepts the
>    parameter and **discards it** — no member stores it. `initializeMumps()`
>    (`:93-105`) sets `id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD)` unconditionally
>    on the Intel-MPI path this fork builds (`_OPENMPI` path is `0`, MUMPS's own
>    WORLD), and the rank/size probe two lines later reads `MPI_COMM_WORLD` directly.
>    **This is a real, bounded plumbing bug**, not an architectural dead end: store the
>    passed communicator, `MPI_Comm_c2f` it, thread it through `MumpsParallelSOE`
>    construction (which today never passes one either).
> 2. **Every contour solve is genuinely complex — the §5.2 "real solver per conjugate
>    pair" claim was wrong**, refuted against the FEAST v3/v4 User Guides and the MKL
>    `?feast_srci` reference (§5.2 above has the corrected paragraph and citations).
>    `dfeast_scsrgv`'s convenience driver doesn't dodge this either — it calls complex
>    PARDISO internally. The fork's local MUMPS build is real-only
>    (`mumps-build/CMakeCache.txt:417` `arith=d`; `mumps-install/lib/` has `dmumps.lib`
>    only, no `zmumps.lib`) and `MumpsParallelSolver.{h,cpp}` `#include <dmumps_c.h>`
>    exclusively — there is no complex code path in the SRC mumps wrapper today.
>
> **Decision — "Outcome A-prime"**: neither the clean "yes" (Outcome A, ship as-is) nor
> vendoring PFEAST (Outcome B) is correct. MKL FEAST's RCI is architecturally fine and
> stays the driver (**no PFEAST vendoring** — it would buy nothing these two fixes
> don't, and adds a Fortran build dependency); but P3 must build:
>   (a) the sub-communicator plumbing fix in `MumpsParallelSolver`/`MumpsParallelSOE`
>       (item 1), verified with 2+ independent sub-comms each solving a distinct real
>       system concurrently — the actual "D2 spike toy" test, now a P3 unit test
>       instead of a throwaway experiment; and
>   (b) a complex inner solve. **Recommended: the symmetric real $2n\times 2n$
>       block-augmentation** (numerically verified exact to 1e-15 relative error, see
>       [[_modal_family_handoff]]): for $z=a+bi$, $(zM-K)X=B$ (real $B$) is equivalent to
>       $\begin{bmatrix}aM-K & -bM\\-bM & -(aM-K)\end{bmatrix}\begin{bmatrix}X_r\\X_i\end{bmatrix}=\begin{bmatrix}B\\0\end{bmatrix}$,
>       symmetric indefinite (MUMPS `SYM=2`/LDLᵀ), reusing the existing real `dmumps` —
>       **zero new external dependency**, at the cost of a new inner-SOE class and
>       ~2× the factorization work of a hypothetical *complex* $n\times n$ solve
>       (~8× the current *real* $n\times n$ solve in dense-flop terms; sparse fill is
>       a third number — measure it). The zmumps alternative (build complex MUMPS + a
>       `ComplexMumpsParallelSOE`) is the fallback if block-real proves too slow or
>       too ill-conditioned in practice — the P3 validation bundle must record BOTH a
>       cost comparison AND a conditioning sweep: the equivalent-real formulation is
>       known to degrade vs the native complex system as $\mathrm{Im}(z_j)\to0$
>       (flattened-ellipse contour points near the real axis), so assert residuals
>       across the contour aspect-ratio range, not just mid-contour.
> This is de-risked exactly as planned: P1–P2 (serial) are unaffected either way.

> [!question] **R2 — MPI integer size / MKL threading model (ABI).**
> FEAST/MKL are Fortran; MUMPS and MKL each have LP64 vs ILP64 builds. A mismatch
> between OpenSees's MPI int width, MKL's `?feast` int width, and MUMPS's int width is
> the classic Fortran/MKL ABI trap (cf. `BUILD_GOTCHAS.md`, dSAUPD/dsygvx history).
> Also: MKL's internal threading (`MKL_NUM_THREADS`) nested under MPI sub-comms can
> oversubscribe. Pin int width + threading in the compilation journal.

> [!question] **R3 — Contour choice & $m_0$ guess.**
> Too-tight an ellipse aspect ratio $r$ slows convergence near band edges; too-loose
> captures spurious modes outside $[f_1,f_2]$ (rejected by Ritz filtering but wasteful).
> $m_0$ under-guess → saturation (handled by §4.4 auto-enlarge); over-guess → wasted
> solves. Defaults from Kratos ($1.5\times$/$2\times$) + auto-enlarge mitigate.

> [!question] **R4 — Subspace not converged (FEAST refinement loop).**
> The outer FEAST loop (§4.3) may need several iterations for ill-conditioned pencils
> (near-degenerate modes on the contour, soft modes). Expose `-maxiter`/`-tol`; warn
> on non-convergence (mirror ArpackSolver's `info<0` diagnostics
> `ArpackSolver.cpp:293+`). Degenerate eigenvalue *on* the contour is the pathological
> case — nudge the contour or widen the band.

> [!question] **R5 — MUMPS factorization reuse across contour points.**
> Each $z_j$ needs its *own* factorization of $(z_jM-K)$ (different shift) → $N_q$
> factorizations per FEAST outer iteration. They cannot share a factor, but the
> **symbolic** factorization (sparsity pattern of $K,M$) is identical across all $z_j$
> and all outer iterations → factor it once, reuse the analysis phase (MUMPS `JOB=1`)
> for every numeric factorization (`JOB=2`). This is a real, large speedup and a
> design requirement on the inner-SOE reuse, not an afterthought.

> [!question] **R6 — SP/MP build-flag surgery risk (P4).**
> Touching `OpenSeesCommands.cpp:4079`'s `_PARALLEL_INTERPRETERS` guard and the
> parallel mains is core/vanilla territory — high blast radius (it gates *all* of
> OpenSeesMP). Must be additive (a new coherent configuration), not a rewrite of the
> existing SP/MP behavior; full adversarial gate + LEDGER_vanilla_files row mandatory.

- **Backwards compatibility:** purely additive — a new `-feast` branch + two new
  classTags. Existing `eigen`, ARPACK, `modalProperties` untouched; default behavior
  byte-identical.
- **Windows note:** `SymmGeneralizedEigenSOE` is already disabled on Windows
  (`OpenSeesCommands.cpp:340-346`); MKL FEAST *is* available on the Windows oneAPI
  build (it is part of MKL), so P1 is reachable on the primary dev box — a portability
  advantage of Option A over a vendored-Fortran Option B that would need the Windows
  Fortran toolchain.

---

## 10. Ledger / header / PR obligations

- **`LEDGER_implementations.md`** — at design-time add two RESERVED rows:
  `FeastEigenSOE` (EigenSOE, classTag **33022**, files
  `SRC/system_of_eqn/eigenSOE/FeastEigenSOE.{h,cpp}`, status *design/ADR 43*) and
  `FeastEigenSolver` (EigenSolver, classTag **33023**, files
  `FeastEigenSolver.{h,cpp}`). Flip to *active* per phase.
- **`LEDGER_vanilla_files.md`** — rows for every upstream file touched:
  `FEM_ObjectBrokerAllClasses.cpp` (+ `TclPackageClassBroker.cpp`,
  `FEM_ObjectBroker.cpp` — the `getNewEigenSOE` switch), the two `eigen` parsers
  (`analysis.cpp`, `OpenSeesCommands.cpp`), `BasicAnalysisBuilder.cpp`
  (`newEigenAnalysis`), and at P4 the `_PARALLEL_*` guard in `OpenSeesCommands.cpp` +
  parallel mains. Mark each edit in-source with a `// Ladruno ADR43 …` comment so the
  table is reconstructable via `grep -rn "Ladruno" SRC/`.
- **`LEDGER_quirks.md`** — record: ArpackSOE `sendSelf` ships only a handshake ID not
  the matrices (the distributed-assembly contract FEAST must follow); MUMPS symbolic-
  factor-reuse across contour shifts (R5); MKL/MUMPS int-width ABI (R2).
- **Banner** — add a `Ladruno_scripts/banner_features.txt` line (e.g.
  *"FEAST band-targeted parallel eigensolver (eigen -feast)"*) once a phase ships,
  then `python Ladruno_scripts/patch_banner.py` + rebuild. (Every `shipped` row gets a
  banner line.)
- **Header stamp** — run `Ladruno_scripts/stamp_headers.py` on the four new files and
  add them to its GLOBS (the four-author LADRUNO header is non-optional for
  fork-authored files, per [[feedback_always_stamp_header]]).
- **`Ladruno_internal/01_compilation_journal.md`** — record the FEAST/MKL (and, if
  Option B, vendored PFEAST) build dependency: MKL link already present; document the
  `?feast_*` symbols, int-width pin, and (P3) the MPI-comm-split / MUMPS-per-subcomm
  recipe.
- **PR base** — base on **`ladruno`** (default branch), one logical PR per phase;
  verify each prior PR `state==MERGED` before stacking (per
  [[feedback_stranded_commits_after_automerge]] and the global stacked-PR `--base`
  pitfall).

---

## 11. Cross-references

- [[modal_gap_study/00_SYNTHESIS]] **§3** — "the parallel / MP answer"; this ADR *is*
  track **C** of the proposed ADR family (the strategic unifier). §5 build order:
  A→B→**C**→D.
- [[modal_gap_study/01_opensees_current_state]] — the file:line ground truth this ADR
  is grounded in (§C parallel reality, §E ranked gaps 1/2/3 closed here).
- [[modal_gap_study/03_kratos_source]] §2/§5c/§6 — the PFEAST template (3-level MPI,
  bundled Fortran, no Anasazi); §5c FEAST RCI/`fpm` knobs.
- [[modal_gap_study/04_lsdyna_theory]] §1/§4 — MPP eigen = distributed factorization +
  Sturm/inertia certification (the completeness guarantee we replicate in P2).
- [[46_ladruno_complex_modal_adr]] (**ADR 46**) — complex/state-space damped modal;
  P5 here re-hosts it at scale via `zfeast_*` complex contours.
- [[42_ladruno_buckling_adr]] (**ADR 42**) — linear buckling $K\phi=\lambda K_g\phi$;
  benefits directly from FEAST band-targeting (buckling-load *band*) + Sturm
  certification of the critical cluster. (ADR 44: parallel `modalProperties`/CQC — the
  post-processing companion — tracked as synthesis gap E.3/E.5.)
