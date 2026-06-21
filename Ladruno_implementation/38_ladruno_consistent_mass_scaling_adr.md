---
title: Consistent (Olovsson) selective mass scaling — the v2 of CentralDifferenceSMS
project: Ladruno
status: IMPLEMENTED — Zone-A 6/6 (oracle + C++ integrator built & passing; pending PR)
priority: high
owner: nmora
tags:
  - implementation
  - adr
  - integrator
  - explicit
  - mass-scaling
  - olovsson
---

# Consistent (Olovsson) selective mass scaling

> **The v2 of the mass-scaling stack.** Lumped (conventional, `DT2MS`-style) scaling
> shipped as ADR 36 (`CentralDifferenceSMS`, classTag 33007). This ADR is the
> **consistent / Olovsson** increment parked there (decision M4, open-question
> "lumped vs consistent"): inject a *consistent* element scaling-mass that preserves
> the global frequencies at the same %added-mass, instead of an isotropic diagonal.
> **Phase 0 (numpy oracle) is complete and proves the formulation** — see
> [[mass_scaling_consistent/README]] and `oracle_olovsson_sms.py`.

## What

For every sub-target element (`dt_e < dtTarget`), instead of adding the isotropic
diagonal increment `(s−1)·m_a` to each node (lumped v1), add the **centroidal Olovsson
scaling mass**

```
M_bar_e = beta_e · [ diag(m_a) − m m^T / M_e ]            (per spatial direction)
beta_e  = (dtTarget/dt_e)^2 − 1     (the SAME factor as lumped v1; betaK-damped variant
                                     s = T^2 + 2*T*c carries over verbatim)
```

with `m = (m_1..m_n)` the element's lumped nodal masses and `M_e = Σ m_a`. The row sums
of `M_bar_e` are zero, so `M_bar_e · t_rigid = 0`: **rigid-body translation receives no
added inertia**. The added inertia loads only the relative/deformation modes that govern
the explicit step, so the global low modes (f1 — non-negotiable for seismic SSI) are
preserved while `dt_e → dtTarget` exactly as in v1.

The scaled global mass `M̃ = M_lump + Σ_e M_bar_e` is SPD but **non-diagonal** (the
`−m m^T/M_e` term couples the element's nodes), so the explicit leap-frog can no longer
use the trivial `system Diagonal` `a = M⁻¹r`. The solve is done **matrix-free** by
preconditioned CG (Phase 2), with the lumped diagonal as the Jacobi preconditioner.

**In scope:** consistent scaling on the explicit central-difference integrator; reuse of
the v1 sizing (`dt_e`, self-report skip, betaK-damped factor, MP-slave constraint
exclusion, `-maxAddedMass` cap); matrix-free PCG solve; `T-CONSISTENT` validation.
**Out of scope:** `c·K` stiffness-proportional scaling (the centroidal form is element-
local, PSD, and needs no K — preferred); consistent scaling under MPI shared nodes
(rides on the same T-MPI gap as v1); the `ExplicitBatheSMS` clone.

## Why

ADR-36 M4 / open-question: *"Lumped shifts global frequencies; for seismic SSI the
fundamental period is non-negotiable. … If f1 drift exceeds 1% at the scaling needed for
real 3D SSI meshes, elevate consistent scaling to a v2 ADR."* The Phase-0 oracle confirms
the failure mode quantitatively: at an aggressive target, **lumped scaling shifts f1 by
−54% to −83%** while **Olovsson holds f1 to −0.1% to −0.2%** — a >5× (here ~400×)
improvement, at ~half the total added mass. This is exactly the "10× more aggressive
scaling" headroom M4 anticipated.

## How — formulation (Phase 0, PROVEN)

Centroidal Olovsson, per the oracle. For a 2-node equal-mass bar this is
`beta·(m/4)·[[1,−1],[−1,1]]`: the `[1,1]` (rigid) mode is untouched, the `[1,−1]`
(relative) mode mass scales by `(1+beta)` ⇒ `omega_e → omega_e/√(1+beta)` ⇒
`dt_e → dt_e·√(1+beta)`, identical to the lumped factor but without moving the low modes.

The solve, matrix-free PCG (`oracle_olovsson_sms.py::pcg`):

```
M̃ x = b ,   apply M̃ x = M_lump·x + Σ_e M_bar_e · x_e   (gather/scatter element DOFs)
precond     = diag(M_lump)^{-1}
```

Oracle: 12–21 CG iters to ≤1e-11 on the test bars (M̄ is a small perturbation of the
dominant lumped diagonal, so the preconditioned operator is well-conditioned ~`1+beta_max`).

## Where (C++, Phases 1–2 — pending)

- `SRC/analysis/integrator/LadrunoMassScaling.h` — **extend**: a `buildMassScalingConsistent`
  that emits per-element `M_bar_e` blocks (node-tag-keyed, with the element's DOF map) in
  addition to / instead of the diagonal `injected` map. Reuse the existing dt_e / guards.
- the SMS integrator — **the solve override** (Phase 2): a matrix-free PCG in the leap-frog
  step. The lumped diagonal still comes from the assembled `system Diagonal` (preconditioner
  + the `M_lump` part of the operator); the `M_bar_e` mat-vec maps element node DOFs to SOE
  equation numbers via `FE_Element::getID()` (the same map `CriticalTimeStep` uses).
- `classTags.h` / brokers — IF a new class (see Decision V-PKG): `INTEGRATOR_TAGS_… 33008`
  (next free). 33008 is also `ELE_TAG_LadrunoTri6`-adjacent — verify against `classTags.h`
  + the G2 within-family gate at allocation time.
- ledger rows, banner line, header stamp — per CLAUDE.md, in the merging PR.

## Decisions

| # | Decision | Rationale |
|---|----------|-----------|
| V1 | **Centroidal `M_bar_e = beta[diag(m_a) − m m^T/M_e]`**, NOT `c·K` | Element-local, PSD, no stiffness fetch, exact rigid-mode annihilation; the textbook Olovsson form. `c·K` couples more widely and risks indefiniteness. |
| V2 | **Matrix-free PCG inside the integrator**, lumped-diagonal preconditioner | Keeps the user recipe (`system Diagonal` still supplies the preconditioner + M_lump); no new SOE/solver class; localizes the consistent complexity to the integrator. Oracle: few iters. |
| V3 | Reuse v1 sizing verbatim (`dt_e`, betaK `s=T²+2Tc`, self-report, MP-slave exclusion, cap) | The scale factor `beta` is identical; only the *placement* of the mass changes. One sizing path for both modes. |
| V-PKG | **DECIDED (2026-06-20): new subclass `CentralDifferenceSMSConsistent`, classTag 33008**, sibling of `CentralDifferenceSMS` off `CentralDifferenceLadruno`; shares the sizing util. | The per-step accel solve is driven by the `Linear` algorithm (`theSOE->solve()` → `update()`), not inside the integrator, so the consistent solve must be injected at the two sites that consume the diagonal result (the `newStep()` starter + `update()`). A tiny protected `virtual refineAccel(Vector&)` no-op hook in the base gives v1 a **byte-identical** path for free (vs. a `-consistent` runtime branch we'd have to prove identical). The subclass also never mutates `Node` mass (the coupling can't live there), so it skips v1's inject/restore-on-destruct lifecycle entirely. |
| V4 | Energy: `EnergyBalanceRecorder` KE must include the `M_bar` term | The recorder sums node/element `getMass()` (diagonal) — it will NOT see the cross-node `M_bar`. Either expose `½ vᵀ M̃ v` from the integrator or document the KE gap. Decide in Phase 3. |

## Open questions / risks

- **V-PKG packaging** (above) — the one decision blocking Phase 1.
- **Energy closure (V4)** — the consistent KE is not visible to the diagonal-summing
  recorder; needs an integrator-side KE or a documented gap. The Phase-0 oracle does not
  exercise the recorder.
- **Projector interaction (ADR-30)** — `LadrunoProjectionHandler` reads `diag(M)` from the
  SOE for `(LᵀML)⁻¹LᵀM`. Under consistent mass the projector's mass is no longer the true
  M̃; v1 already excludes MP-slave elements from scaling (so scaled elements are projector-
  free), but document the interaction.
- **First-step starter** — the starter `a_0 = M⁻¹(P_0 − Cv_0 − F_int)` must also go through
  the PCG.

## Implementation log

- **2026-06-20 — Phase 0 complete.** numpy oracle `oracle_olovsson_sms.py` proves C1
  (f1 preservation, >5× vs lumped), C2 (less mass, same step), C3 (PCG 12–21 iters).
  Formulation V1 + solve V2 locked.
- **2026-06-20 — Phases 1–3 complete (C++ built & passing).** V-PKG decided: new subclass
  `CentralDifferenceSMSConsistent` (33008). Kernel: `buildMassScalingConsistent` /
  `consistentMatVec` / `consistentPCG` in `LadrunoMassScaling.h`. Base gained a no-op
  `refineAccel(Vector&)` hook (lumped path **byte-identical** — v1 regression 24/24 still
  green); the subclass overrides it to recover `r` from the factored DiagonalSOE
  (`A[i]=1/mass`) and replace the diagonal solve by the matrix-free PCG. Registered at all
  6 sites + banner + ledgers + manifest row. **Zone-A `T-CONSISTENT` 6/6:** supra-stable
  necessity (CD diverges @dtTarget); **f1 −0.17% vs lumped −53.4%** via transient FFT on the
  oracle Case-A bar (matches the numpy oracle −0.135% / −54.18%); NO nodal-mass mutation (vs
  lumped control); reduce-to-base bit-identical to CDL; PCG 3 iters on the small chain.
  **Still open:** V4 (energy-recorder KE does not see the `M_bar` term — diagonal-summing
  recorder, documented gap); projector(ADR-30)+consistent interaction documented not tested;
  MPI shared-node reduction (shared with the lumped T-MPI gap). NEXT: an adversarial review
  pass, then the energy-KE reconciliation (V4).
- **2026-06-20 — adversarial review folded (3-lens, refute-by-default).** Lens 1 (kernel
  math): CONFIRMED CORRECT — Mbar matches the oracle incl. row-sum-zero, betaK closed form
  inverts exactly + reduces to the oracle at betaK=0, matVec gather/scatter index-correct,
  PCG solves M̃a=r to machine precision on SPD M̃, r-recovery exact in all solve modes
  (incl. `-factorOnce`). One **CRITICAL** found + FIXED: under the default
  `TransformationConstraintHandler`, an element touching an MP *retained/master* node has a
  transformed `FE_Element::getID()` whose size EXCEEDS the element's `n` → pairing it with the
  `n×n` M_bar was an out-of-bounds read (fires on the exact SSI/diaphragm models this targets).
  Fix: exclude elements touching EITHER a slave OR a master MP node (`getNodeConstrained` +
  `getNodeRetained`) + a defensive `getID().Size()==n` skip. This ALSO resolves the MEDIUM
  projector(ADR-30)×consistent mass-mismatch (no M_bar lands on any projector equation now).
  Plus: the non-Diagonal-SOE fallback is now flagged up-front in `domainChanged` (it would run
  UNSCALED→unstable, not merely slow); PCG non-convergence is surfaced (one-time warning);
  PCG iter-count off-by-one on non-convergence clamped. Regression
  `test_consistent_excludes_constraint_touching_elements` (TCONS-7). Zone-A now **7/7**;
  v1 lumped regression still byte-identical (24/24).
- **2026-06-20 — V4 decision: documented limitation for v1 (recorder coupling deferred).**
  The `EnergyBalanceRecorder` computes KE = ½vᵀMv from element+nodal `getMass()`, which does
  NOT see the matrix-free `M_bar`; a correct consistent KE needs a recorder↔integrator hook
  (coupling = v2). The leap-frog still conserves the true M̃-energy; the recorder's RES just
  doesn't close for the consistent integrator (it does for lumped). Documented in the manifest,
  handoff, and [[LEDGER_quirks]]; the gap is small in the design regime (low-participation
  scaled elements). NEXT (v2): the recorder KE hook, then T-MPI.
- **2026-06-20 — ExplicitBathe SMS family added (33009 lumped, 33010 consistent).** The
  same kernel (`buildMassScaling*` / `consistentPCG`) and `refineAccel()` hook pattern,
  ported to the Noh-Bathe `ExplicitBathe`. `ExplicitBathe` gained a protected classTag ctor
  + the no-op `refineAccel()` hook called at BOTH sub-step solves (`A_tpdt` from the external
  algorithm solve, `A_tdt` from the inline second solve); default no-op keeps `ExplicitBathe`
  byte-identical (base regression green). `ExplicitBatheSMS` (lumped) needs NO solve hook —
  ExplicitBathe assembles only the mass on the RHS, so nodal injection is seen directly (the
  ADR-36 "trivial follow-up"). `ExplicitBatheSMSConsistent` overrides `refineAccel` with the
  matrix-free PCG; the inline second solve reuses the factored DiagonalSOE (`isAfactored` ⇒
  1/mass at both sites, so r-recovery is valid). Measured identical to the CDL family: f1
  **−0.17% (consistent) vs −53.41% (lumped)** via transient FFT. Command takes the Noh-Bathe
  `$p` first. Zone-A `tests/test_explicitBatheSMS_integrator.py`. LNVD (separate base) is a
  clean follow-up, not done here.
- **2026-06-20 — LNVD SMS family added (33011 lumped, 33012 consistent).** The same pattern
  ported to `ExplicitBatheLNVD` (Noh-Bathe + FLAC local non-viscous damping; a separate
  `TransientIntegrator` base, not an `ExplicitBathe` subclass). Same two-sub-step solve
  structure ⇒ same hook placement (`A_tpdt`, `A_tdt`); protected classTag ctor + no-op
  `refineAccel` (base LNVD + lumped LNVDSMS byte-identical, regression green). The FLAC `alpha`
  is orthogonal to the mass-scaling sizing. Command takes `$p $alpha` first. Measured identical:
  f1 **−0.17% vs −53.41%**. Zone-A `tests/test_explicitBatheLNVDSMS_integrator.py`. The
  explicit-integrator SMS axis (CentralDifference / ExplicitBathe / ExplicitBatheLNVD) is now
  complete for both lumped and consistent.
- **2026-06-21 — V5: distributed (MPI) consistent PCG — DONE.** The serial `consistentPCG`
  is rank-local (local inner products, no shared-DOF `M̄` exchange), so the consistent path
  was wrong at partition boundaries under OpenSeesMP. V5 adds a distributed PCG
  (`consistentParPCG`/`consistentParMatVec` + `ConsistentParOps` in `LadrunoMassScaling.h`)
  built on ONE ownership-free weight `wᵢ=1/multiplicityᵢ` (mult = #ranks sharing DOF i):
  the matvec applies the GLOBAL (replicated) lumped diagonal WEIGHTED + off-diagonal `M̄ₑ`
  in FULL then sum-assembles shared DOFs across ranks (diagonal → full-once, off-diagonal
  accumulates); global inner products use the same `w` + all-reduce; every CG control scalar
  is global ⇒ identical iteration count on all ranks ⇒ collectives stay in lockstep (no
  deadlock). Collapses byte-identically to the serial PCG at `np=1` (no neighbours, `w≡1`,
  identity reduce) ⇒ serial path untouched. **Architecture:** the 3 consistent integrators
  live in the shared `OpenSeesLIB` and cannot reference the MP-only `MPIDiagonalSOE`
  (`_PARALLEL_INTERPRETERS` + that SOE exist only in the MP executables), so dispatch goes
  through 4 new `LinearSOE` base virtuals (serial no-op defaults; `MPIDiagonalSOE` overrides
  `getScalingDiagonalA`/`assembleSharedSum`/`globalReduceSum`/`isDistributedDiagonal`). The
  3 `refineAccel` bodies are factored into the shared `LadrunoConsistentRefine.h` (collective-
  safe early-returns). **Validation** (`Ladruno_implementation/mass_scaling_mpi/`, local
  2-rank, NOT CI-gated — single-process CI): MPI `np=2` == serial `DiagonalSOE`+`consistentPCG`
  gold reference (max abs diff 0 to recorder precision, final 1.746e-3) AND `np=1`==`np=2` AND
  consistent (1.746e-3) ≠ lumped (2.35e-3) ⇒ genuinely Olovsson, partition-invariant. Serial
  Zone-A 36/36 green (no regression from the `refineAccel` refactor). Two bugs caught & fixed
  en route: `build.bat` honored only `%1` (built only the first target → stale MP binary), and
  the A/B comparator passed two diverged runs (now rejects non-finite). With the lumped path
  already parallel-correct (T-MPI), the parallel SMS axis is COMPLETE. Lumped/consistent under
  OpenSeesSP (`system Diagonal`→DistributedDiagonalSOE) is the remaining (deferred) increment —
  it would need the analogous Channel-based assemble, behind the same virtuals.
