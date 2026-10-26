---
title: "ADR 213 — Modal basis completeness: shared Krylov-Ritz kernel + residualModes (WP-213)"
project: Ladruno
type: ADR
status: "PR 2 implemented (draft PR) — kernel C + residualModes D; ritzVectors (E) is PR 3"
priority: medium
owner: nmora
related:
  - "[[44_ladruno_frequency_domain_adr]] (§9 deferred the missing-mass correction)"
  - "[[188_ladruno_factor_reuse_adr]] (factor reuse on the analysis LinearSOE)"
  - "WP-211 / PR #1020 (modalProperties keeps MP-retained mass under Transformation)"
  - "WP-212 / PR #1021 (RSA mass-shortfall guard, -requireMass)"
  - "apeWorkbench opensees-dynamic-analysis, branch work/claude/rsa-gap-study: studies/rsa-gaps/{gap-study.md, fork-implementation-guide.md, v1/}"
tags: [adr, modal, response-spectrum, missing-mass, ritz, asce7, wp-213]
updated: 2026-10-26
---

# ADR 213 — Modal basis completeness for RSA

## 1. Context

- **Codes.** ASCE 7-22 §12.9.1.1 asks for modes carrying 100 % of the mass;
  modes under 0.05 s may be represented by one rigid mode at 0.05 s, and 90 %
  per direction is the exception. NEC-SE-DS §6.2.2e asks for ≥ 90 %.
- **Models.** apeGmsh shell and solid models spread mass over hundreds of
  local modes: 90 % can take hundreds of eigenpairs, 100 % is out of reach
  without a residual term.
- **Practice.** CSI appends static-correction (missing-mass) modes to the
  eigenvectors, or uses load-dependent Ritz vectors (CARM §19.2.5, §19.3;
  Wilson ch. 14). NASTRAN calls them residual vectors.
- **The fork.** ADR 44 §9 deferred the correction. PR 1 (WP-212, #1021)
  reports the shortfall; this ADR adds the correction (PR 2) and prepares
  Ritz vectors (PR 3).

## 2. Decisions

### D1 — one kernel (`LadrunoRitzKernel`)

`SRC/analysis/analysis/LadrunoRitzKernel.{h,cpp}`: block Krylov–Ritz on
(K, M), shared by `residualModes` (this PR) and `ritzVectors` (PR 3).

- **K** is what `eigen` assembles (every FE's `addKtToTang(1.0)`, constraint
  FEs included) at the current state, so gravity + `loadConst` + PDelta gives
  K − Kg. It is assembled into the **analysis' LinearSOE**; the first solve
  factors it and every later solve reuses the factor (the `factored` flag ADR
  188 relies on). The first batch of solves is checked against an
  element-level K·x (relative residual < 1e-4), which refuses an SOE that
  cannot hold K (diagonal / explicit systems, a stale size).
- **M·x** goes through `FE_Element::getM_Force` / `DOF_Group::getM_Force` in
  equation numbering (the inertia-load path): Transformation gives TᵀMT,
  Penalty / Lagrange use the handler's own objects, multiplier rows carry no
  mass, FEs without an Element (constraint, contact) are skipped (WP-183 rule).
- **Orthogonalization**: deflate against an optional basis Φ (copied and
  M-orthonormalized internally, so any eigen normalization works) and the
  earlier Krylov columns, twice ("twice is enough"), then normalize; a column
  whose M-norm fell below `dropTol ×` its norm right after the solve is
  dropped and counted per starting vector.
- **Rayleigh–Ritz**: K̃ = XᵀKX (element sweep), LAPACK `dsyev`; a
  non-positive or non-finite Ritz value is refused.
- **API for PR 3**: `Options{dropTol, maxBlocks, maxVectors, cyclesPerStart,
  solveStart}` and `Result{psi, lambda, keptPerStart, droppedPerStart,
  numDropped, maxSolveResidual}`; helpers `applyM`, `applyK`,
  `rigidInfluence(dir)` (zero at SP-fixed DOFs), `domainModes`, `install`.
  `solveStart = true` is LDR (block 0 = K⁻¹F, loads); `false` takes F as shapes.
- **No classTag**: the kernel and the command are not domain objects
  (nothing `TaggedObject` / `MovableObject`), so the manifest gate has nothing
  to reserve.

### D2 — install into the domain's mode set

Pseudo-modes are appended through `AnalysisModel::setNumEigenvectors` /
`setEigenvector` / `setEigenvalues`, the calls `eigen` makes.
`TransformationDOF_Group::setEigenvector` rebuilds the constrained nodes.
`modalProperties`, `responseSpectrumAnalysis` (per mode, `-combine`,
`-increment`), `modalResponseHistory`, `frequencyResponse`, recorders and
apeGmsh's `spectrum_after_gravity` use them unchanged. `residualModes`
unsets the domain's modal properties, so a stale `modalProperties` cannot be
paired with the new set.

**Mode-kind registry** (`LadrunoModeKinds`, decided in the guide's §11 q1):
Ladruno side, keyed by `Domain*`, no `Domain` member. A record is valid only
while the domain's eigenvalue vector is bit-identical to the snapshot taken
when it was written, so a new `eigen` (or `wipe`, a FEAST trim) clears the
tags without touching the eigen path. `residualModes -kinds` returns them
(0 eigen, 1 residual, 2 Ritz). Tagging the rows of `modalProperties -print`
is deferred: `DomainModalProperties.cpp` is WP-211's file.

### D3 — the residual space is span{r_d, K⁻¹M r_d} (deviation from the guide)

The guide (§6) used one block from F = M ι_d: the classical static
correction x_d = K⁻¹M ι_d deflated against Φ. Two facts measured on the V3
podium model (numpy, before any C++) changed the design:

1. **One static-correction vector cannot reach 100 % effective mass.** Its
   captured mass is (xᵀMr)² / (xᵀMx · rᵀMr) of the residual mass rᵀMr
   (r = ι − ΦΦᵀMι), which by Cauchy–Schwarz is 1 only when x ∥ r, i.e. when
   all missing modes share one frequency. On V3 it is 93.8 % of the
   residual at every stiffness scale (total 94.9 % with 3 modes). What
   reaches 100 % is the *static* participation (CSI's "static load
   participation ratio"), not the `modalProperties` ratio.
2. **Starting the Krylov sequence from the shape ι_d fixes it.** Block 0 =
   r_d (the mass the modes miss), block 1 = K⁻¹M r_d = x_d (because
   K⁻¹MΦ ⊂ span Φ, this is exactly the static correction). With r_d in the
   span the effective mass is 100 % by construction; with x_d in the span
   the static response to M ι_d is reproduced exactly. On V3 the base shear
   error with 3 eigenmodes drops from +4.0 % (one vector, `kbase` 1000) /
   +1.6 % (`kbase` 300) to +0.55 %.

So `residualModes` runs the kernel with `solveStart = false`, two blocks, up
to 2 pseudo-modes per direction (a column already spanned is dropped — e.g.
a direction a single mode carries completely).

Default drop tolerance for `residualModes` is `1e-4` on the M-norm ratio (a
dropped column carries < 1e-8 of the direction's mass). It has to sit above
the eigenvector error, not at round-off: on V1 with `-dir 1 2` the X residual
is pure eigenvector noise (3.6e-6 of ‖ι_x‖_M, 1e-11 of the mass, from
`eigen -fullGenLapack` under Lagrange and Penalty alike). Kept at `1e-6`, that
noise vector matched the genuine torsional content of the Y static correction
on the massive DOFs, which then dropped; but the noise carries arbitrary
values on the massless DOFs, which the M-norm cannot see, so the Ritz pair
came out at T = 0.033 s instead of 0.313 s. The kernel's own default stays
`1e-10` (LDR load starts, PR 3).

### D3b — periods: `rayleigh` (default) or `-period $T`

- `rayleigh` keeps the Ritz values.
- `-period $T` sets λ = (2π/T)² on every pseudo-mode (the guide's ASCE
  §12.9.1.1 "rigid mode at 0.05 s"). **Limitation:** the RSA displacement of
  a mode is φ Γ Sa(T)/λ, so a fixed λ changes the amplitude by λ_R/λ_T, not
  only the Sa lookup. On V3 (pseudo-mode Rayleigh periods 0.31 s and
  0.16 s), `-period 0.05` gives 0.085 × the reference base shear;
  `rayleigh` gives 1.006 ×. The two bracket the reference, as the guide
  asked, but `-period 0.05` is only conservative when the Rayleigh periods
  are *shorter* than 0.05 s (the case it is meant for: genuinely rigid
  missing mass). The modal-mass ratio is unaffected (the shape is). The
  physically right ASCE rigid mode — the static response times Sa(0.05) —
  needs RSA to look Sa up at an overridden period while dividing by the true
  λ: a `ResponseSpectrumAnalysis` change, left for after WP-212 merges
  (§6 follow-up 1).

### D4 — Ritz vectors as an alternative to `eigen` (PR 3)

`ritzVectors $nvec -dir ... -pattern ... -cycles ...` will call the same
kernel with Φ empty, `solveStart = true`, `maxVectors = nvec`,
`cyclesPerStart`, tag the result `RITZ` and install it with `install()`.
Nothing in this PR implements it.

### D5 — no rigid/periodic split in v1

No Gupta split (`-combine GMC`, CARM §20.5.1, Gupta & Chen 1984); the
pseudo-modes carry the rigid mass and combine by CQC like any mode (their
periods are far from the dynamic ones, so ρ ≈ 0 against them).

### Rejected — an `EigenSOE` subclass

An `EigenSOE` sees K and M but not node directions (ι_d), load patterns
(PR 3) or the domain's existing modes, and `residualModes` must **append** to
a set an earlier `eigen` produced (FEAST included). A command over the
analysis' LinearSOE keeps both eigen paths and every downstream consumer
untouched.

## 3. Command

```tcl
residualModes -dir 1 <2 <3>> <-period rayleigh|$T> <-replace> <-tol $t>
residualModes -kinds
```

- After `eigen` / `eigen -feast`; then `modalProperties` again before RSA.
- Returns the appended eigenvalues (empty when everything was already
  spanned: the set is then left as it was).
- `-replace` drops the earlier residual modes first; without it a second call
  appends only what the current set (residual modes included) misses.
- Directions are global translations 1..ndm with free mass; ι_d is zero on
  SP-fixed DOFs under every handler.
- Refused: no modes, no LinearSOE, a partitioned (OpenSeesSP) domain
  (`Domain::ladrunoThreadedUpdateAllowed()` is false exactly there).
- Prints, per direction, the free mass the modes carry and what the
  pseudo-modes add, and warns if a pseudo-mode's Rayleigh period exceeds the
  shortest retained eigen-period (the missing mass is then not rigid at this
  truncation).
- Wired once in `SRC/interpreter/LadrunoCommandTable.h` (Python, DL Tcl,
  classic Tcl via the `residualModes` bridge in `SRC/tcl/commands.cpp`).

## 4. Verification (`tests/test_wp213_residual_modes.py`, numpy oracles)

| Check | Result |
|---|---|
| V3 (a) 3 eigenmodes | modal mass 18 %, base shear 0.058 × reference |
| V3 (b) + `residualModes -dir 1` | modal mass 100 % (1e-10), base shear 1.0055 × reference; Ritz values = numpy kernel (1e-9), shapes up to sign; ΨᵀMΨ = I, ΨᵀMΦ = 0 (1e-10), ΨᵀKΨ = Λ |
| V3 (c) `-period 0.05` | λ = (2π/0.05)²; base shear = numpy oracle of the same rule (1e-8), 0.085 × reference; brackets with (b) |
| full-size basis | 3-storey chain: 1 mode + residual = `eigen -fullGenLapack 3` (1e-10) |
| `-replace`, re-run, `eigen` clears tags | kinds `[0,0,0,1,1]` → unchanged → `[0,0,0,0]` |
| rigid diaphragm, mass on constrained nodes | Penalty and Lagrange: modal mass 100 % X and Y, eigenvalues and CQC base shears = 3-DOF reference (1e-6), handlers agree (1e-7); Transformation: strict xfail until WP-211; constrained-node recovery and M-orthonormality pass now |
| preloaded PDelta column | Ritz values = preloaded eigenvalues 2, 3 (1e-8), differ from unloaded |
| `modalResponseHistory` | 3 modes + residual within 10 % of all modes (3 modes alone > 50 % off) |
| classic Tcl | `OpenSees.exe` returns the same eigenvalues and kinds as Python |

## 5. Consequences

- The analysis' LinearSOE holds K and its factor after `residualModes`; the
  command calls the algorithm's `domainChanged()` so WP-188 factor reuse
  cannot start a step on it. Every other algorithm forms its tangent at the
  start of a step.
- Lagrange needs a solver that factors an indefinite KKT matrix (UmfPack,
  Mumps, Pardiso); ProfileSPD will fail the solve.
- Under Penalty / Lagrange a node with mass on an SP-fixed DOF: ι is zero
  there, while `modalProperties` counts that mass as free (WP-211 note), so
  the ratio stays below 100 % by exactly that mass.

## 6. Follow-ups

1. RSA: Sa at an overridden period for RESIDUAL modes, true λ in the
   amplitude (makes `-period 0.05` the ASCE rigid mode). Touches
   `ResponseSpectrumAnalysis.cpp`; after WP-212.
2. `modalProperties -print` tags R/Z (after WP-211).
3. PR 3: `ritzVectors` on the same kernel.
4. apeGmsh: `apeSees.residual_modes(dirs, period=...)`, and read 100 % from
   `modalProperties` only with residual modes present.
