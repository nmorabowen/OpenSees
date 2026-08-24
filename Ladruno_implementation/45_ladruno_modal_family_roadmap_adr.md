---
title: "ADR 45 — Modal-analysis family: implementation roadmap & sequencing plan"
project: Ladruno
type: ADR / program plan (umbrella over ADRs 46/42/43/44)
status: draft — planning, NO code
priority: high
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"            # the cross-code theory + §6 load-bearing assessment this plan executes
  - "[[46_ladruno_complex_modal_adr]]"            # member: complex/state-space modal (33019)
  - "[[42_ladruno_buckling_adr]]"                 # member: prestressed modal + linear buckling (33021)
  - "[[43_ladruno_feast_eigensolver_adr]]"        # member: band-targeted parallel eigensolver (33022/33023) — the substrate
  - "[[44_ladruno_frequency_domain_adr]]"         # member: FRF/SSD/random + modal transient (33024)
  - "[[modal_gap_study/01_opensees_current_state]]" # ground-truth file:line audit
  - "[[66_ladruno_solidshell_adr]]"               # classTag boundary: 33020 taken by LadrunoSolidShell (shipped #482; we skip it)
  - "[[Ladruno_explicit_roadmap]]"                # sibling precedent: a program roadmap spanning several ADRs
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags: [adr, program-plan, roadmap, modal, eigen, complex-modes, buckling, feast, parallel, frequency-domain, sequencing]
updated: 2026-07-08
---

> [!warning] The `status:` above is STALE — this ADR has shipped
> Its frontmatter still carries the pre-implementation value. Trust
> [[LEDGER_implementations]] for *does it work / which PR*, and this ADR for *why*.
> Flagged 2026-08-23 by a ledger audit; see [[README]] §Conventions. Remove this
> banner when `status:` is corrected.

# ADR 45 — Modal-analysis family: implementation roadmap & sequencing plan

**Status:** draft. **Planning only — no code.** This is the *umbrella ADR* governing the rollout of
the modal-analysis family (ADRs **46 / 42 / 43 / 44**, candidate **47**). The per-feature ADRs
specify *what each feature is and how it works*; this ADR decides *what we build, in what order, and
which cross-cutting decisions gate the program*. Theory lives in
[[modal_gap_study/00_SYNTHESIS]] (+ the four code dossiers).

> [!info] One-line thesis
> Bring OpenSees modal analysis to commercial-code parity (complex/damped modes, prestressed +
> buckling, robust band-targeted **parallel** eigensolving, frequency-domain response) by building
> **one shared eigensolver substrate** (ADR 43) and a **cheap serial complex-modal proof** (ADR 46)
> first, then layering the opportunistic deliverables (42, 44) and a future ROM extension (47).

---

## 1. Why a program ADR (not just the four feature ADRs)

The four feature ADRs were scoped independently, but they are **not independent to build**:

- They **share a substrate** — every one rides the eigensolver, and three of four ride an
  *assembled-operator* contract (`M`, `C`, `K`) that OpenSees only partially exposes today.
- They **share cross-cutting decisions** — the assembled-`C` accessor (46, reused by 44 and FEAST
  complex), the `-shift` exposure on the `eigen` command (42, touches shared code), the
  **MKL-FEAST-vs-vendored-PFEAST** build-dependency call (43), and the vanilla-footprint policy
  (CQC edit in 44, SP/MP build-flag surgery in 43).
- They have a **non-obvious optimal order** — the cheapest ADR (46) and the most load-bearing ADR
  (43) are *different* ADRs, so "build cheapest first" and "build foundation first" disagree. That
  tension needs a decision, recorded once, here.

Precedent: [[Ladruno_explicit_roadmap]] plays the same umbrella role for the explicit-dynamics ADRs.

---

## 2. The family at a glance

| ADR | Feature | classTag | Effort | Load-bearing (synthesis §6) |
|---|---|---|---|---|
| **46** | `LadrunoComplexEigen` — complex/state-space modal (non-classical damping) | 33019 | **S** | Domain-enabling (SSI/DRM/isolation/dampers) |
| **42** | `LadrunoBuckle` — prestressed modal + linear buckling | 33021 | **S–M** | Standalone analysis (modest) |
| **43** | `FeastEigenSOE`/`FeastEigenSolver` — band-targeted **parallel** eigensolver | 33022/33023 | **L** | **Substrate (highest)** + general SP/MP fix |
| **44** | `LadrunoModalResponse` — FRF/SSD/random + modal transient | 33024 | **M**¹ | Deliverable (none downstream) |
| *(ROM)* | *ROM / Craig–Bampton substructuring (candidate, future; ADR number assigned when drafted)²* | *TBD* | *L* | *Rides the family; biggest forward unlock* |

¹ ADR 44's "M" reflects *breadth* (three sub-deliverables: FRF, SSD/random, modal transient), not
core risk — per-core-risk it is the cheapest member ("almost zero new core", ADR 44 §2).
² Originally "candidate ADR 47"; **47 has since been taken** by
[[47_ladruno_contact_deferrals_adr]] (the family's third numbering collision — see §9). We no
longer pre-assign numbers to unwritten ADRs; the ROM candidate gets the next free number when its
ADR is actually drafted.

`33020` is **deliberately skipped** — taken by `LadrunoSolidShell`, which has since **shipped**
([[66_ladruno_solidshell_adr]], [#482](https://github.com/nmorabowen/OpenSees/pull/482);
`classTags.h:939`). The family's own tags 33019/33021–33024 remain undefined in `SRC/classTags.h`
(re-verified 2026-07-06).
ADR **41** is the unrelated mortar/ALM contact ADR already on `ladruno`.

---

## 3. Dependency graph (what blocks what)

```
                ADR 43  FEAST eigensolver (MP-parallel)        ← SUBSTRATE
                  │   (band-target · Sturm · distributed dmumps under MP; SP unsupported)
   ┌──────────────┼───────────────────────────┐
   │ (serial eigen already enough)             │ (parallel + complex contours)
   ▼              ▼                             ▼
 ADR 46        ADR 42                        ADR 46 @ scale
 complex       buckling/Kg                   (= P-E, ADR 43 P5 complex contours —
 (serial OK)   (serial OK)                    NOT delivered by the real-parallel P-C)
   │              │                               │
   └──────┬───────┘                               │
          ▼                                       ▼
        ADR 44  frequency domain (FRF/SSD/random) — consumes 46 + the eigensolver
          │
          ▼
        ROM / Craig–Bampton (future, ADR number TBD) — needs trustworthy basis + parallel eigen
```

Key reading: **46, 42, 44 can each ship on the *existing* serial ARPACK `eigen`.** ADR 43 is not a
*blocker* for a first version of any of them — it is the **scale + robustness + parallel** upgrade
that makes them trustworthy on large/partitioned models and re-hosts 46's complex case via complex
contours.

---

## 4. Sequencing decision

Two defensible orders (synthesis §5 vs §6.5):

- **Cheapest-first proof:** 46 → 42 → 43 → 44. Validates direction with the smallest spend.
- **Unlock-the-most:** 43 → (46, 42) → 44. Builds the load-bearing substrate first.

**Decision — a hybrid that front-loads the cheap proof, then the substrate:**

| Phase | Work | Why here | Depends on |
|---|---|---|---|
| **P-A** | **ADR 46 P0–P3** — complex modal, *serial*, on existing `eigen` | Cheapest, directly serves the research portfolio (isolation/dampers/SSI); de-risks the projection+QZ approach before any big build | existing `eigen`, `modalProperties` |
| **P-B** | **ADR 43 P1–P2** — serial **MKL-FEAST** eigensolver (band-target + Sturm) | The substrate; **zero new build dep** (MKL already linked); validated vs ARPACK | MKL |
| **P-C** | **ADR 43 P3** — **distributed inner solve under MP** (per-contour `(z_jM−K)` via distributed `dmumps` on the MP communicator) — **SHIPPED #532**; MP is the single blessed parallel config, **P4 SP/MP-unification DE-SCOPED** | The strategic payoff: large-model modal in the MP build; the reusable distributed-inner-solve seam | P-B, `dmumps`, `OpenSeesMP`/PyMP |
| **P-D** | **ADR 42** — prestressed modal + linear buckling | Opportunistic; rides serial eigen, gains band/Sturm from P-B; can jump ahead of P-C if a project needs it | corot/PDelta Kg, `eigen` |
| **P-E** | **ADR 43 P5** — complex contours (re-host ADR 46 at scale) | Unifies the complex case onto the parallel substrate | P-A, P-C |
| **P-F** | **ADR 44** — frequency domain (FRF/SSD/random, modal transient) — **SHIPPED**: P1a `modalResponseHistory` #537, P1b RSA `-combine` #539, P2 `frequencyResponse`/`steadyStateDynamics` #544 (+ classic-Tcl #546), P3 `randomResponse` PSD→RMS [#552](https://github.com/nmorabowen/OpenSees/pull/552); **gate G-F MET** (all five items, incl. random-vs-Monte-Carlo). Remaining follow-ups: nodal-force `-load`, cross-PSD input | Deliverable layer; ~~build when a project asks~~ built | P-A, eigen |
| *(P-G)* | *ROM / Craig–Bampton (candidate; ADR number TBD — see §2 note ²)* | *Future; the biggest forward unlock* | *whole family* |

**Rationale.** P-A buys confidence + an immediately useful research deliverable for ~S effort. P-B
starts the substrate at *zero dependency cost* (the most-load-bearing work that carries no build
risk). P-C is the one large, build-risky step (MPI distributed inner solve) and is deliberately
sequenced after the substrate is proven serially. **De-scope note (2026-07-08):** P-C originally
bundled a P4 "SP/MP build-flag unification"; that was **dropped**. MP (`OpenSeesMP`/PyMP) is the
single blessed parallel configuration, the distributed FEAST solve already runs there (#532), and
there is no partitioned-SP + distributed-solve combination anyone wants — so the high-blast-radius
build-flag surgery is *not* built. SP is explicitly unsupported for FEAST. P-D/P-F are demand-driven.

---

## 5. Cross-cutting decisions (program-level — decide once, here)

These appear in multiple ADRs; resolving them at the program level prevents divergent choices.

### D1 — Assembled-`C` accessor (from ADR 46; reused by 44, 43-complex)
OpenSees has **no accessor that returns an isolated assembled global `C`**. Precise audit
(2026-07-06): `FE_Element::addCtoTang()` *does* exist (`FE_Element.cpp:383` — folds
`Element::getDamp()` into the tangent), but only ever into the **combined** transient effective
tangent `aM+bC+cK`; there is no path that exports a standalone `C` the way `EigenSOE` gets `K`
and `M`, and no `formEleTangC` anywhere in `SRC/`. (Don't state this as "no `addC`" — a reviewer
grepping `addC` will find `addCtoTang` and reject the premise.) **Decision:** build a single
`LadrunoDampingAssembler` path in P-A and **reuse it verbatim** in ADR 44 and ADR 43's complex
contours. v1 contract = *exactly* `getDamp()` + Rayleigh; **warn** (don't silently absorb) when
`modalDamping`/HHT numerical damping is active. Note the closed form `C̃=αM̃+βK̃` applies to the
**Rayleigh part only**, and only because the projection basis is the real undamped modes
(`M̃=I`, `K̃=diag(ω²)`); element/material dampers require the full `ΦᵀCΦ` projection (ADR 46
§4.6 Route B), never a closed form. Owner: P-A.

### D2 — MKL-FEAST vs vendored PFEAST (from ADR 43 — the big build call) — **RESOLVED 2026-07-07**
Intel MKL's Extended Eigensolver **is** FEAST and the build already links MKL → **serial P-B costs
zero new dependency.** The open question was whether MKL's RCI form lets each contour solve run as an
OpenSees `MumpsParallelSOE` solve on an arbitrary MPI sub-communicator. **D2 spike run** (source read
+ numerical check + Opus literature cross-check, full record [[_modal_family_handoff]] and ADR 43
§5.2/§9 R1): found two independent, concrete facts, neither of which favors vendoring PFEAST.
(1) `MumpsParallelSolver` **hardcodes `MPI_COMM_WORLD`** — the `mpi_comm` constructor argument is
silently discarded (`MumpsParallelSolver.cpp:54-64,97,104-105`) — a bounded plumbing bug, not an
architectural wall. (2) **Every FEAST contour solve is genuinely complex**, refuting ADR 43's old
§5.2 claim that a real solver suffices "per conjugate pair" (checked against the FEAST v3/v4 User
Guides and MKL's `?feast_srci` reference — even MKL's own convenience driver calls complex PARDISO
internally); the fork's local MUMPS build is real-only (`arith=d`, no `zmumps.lib`). **Decision:
stay on MKL FEAST — no PFEAST vendoring.** P-C/P3 must instead ship two OpenSees-side fixes: (a) make
`MumpsParallelSolver`/`SOE` honor a passed sub-communicator, and (b) a complex inner solve —
recommended as a symmetric real $2n\times2n$ block-augmented system (LDLᵀ via the existing real
`dmumps`, `SYM=2`, numerically verified exact to 1e-15) rather than building `zmumps`. Zero new
external dependency either way. Owner: P-C (now unblocked — see ADR 43 phased roadmap P3a/b/c).

### D3 — `eigen` `-shift` exposure (from ADR 42)
Buckling needs a non-zero shift (`ΔK` indefinite) but `eigen` currently hard-zeros the shift
(`shift = 0.0` at `OpenSeesCommands.cpp:2187`, no `-shift` token in the option loop). The plumbing
below is complete and unreachable: `ArpackSOE(shift)` stores and applies a non-zero shift
correctly (`ArpackSOE.cpp:47,215,264`) — so exposure is parser-only work, no solver change.
**Decision:** expose `-shift` on the new `buckling` command in P-D; **only** un-hard-zero the shared
`eigen` path if P-D shows a concrete need (minimize vanilla touch). Owner: P-D.

### D4 — Vanilla-footprint policy (from ADR 44 + ADR 43)
ADR 44's CQC/SRSS touches Petracca's upstream `ResponseSpectrumAnalysis`; ADR 43 touches the
`EigenSOE` base + parallel mains (the SP/MP gate). Sharpening (source audit 2026-07-06): the
upstream `ResponseSpectrumAnalysis` deliberately performs **no combination at all** — it computes
per-mode modal displacements and explicitly defers SRSS/CQC to the user
(`SRC/runtime/commands/analysis/modal/ResponseSpectrumAnalysis.cpp:96-100`, mirrored in
`SRC/analysis/analysis/ResponseSpectrumAnalysis.cpp:362`). ADR 44 would therefore be **adding**
combination logic, not editing existing math — which makes the fork-authored-sibling route even
cheaper than the ADR assumed. **Decision:** prefer **fork-authored siblings** that *read* committed
state where feasible (ADR 44 path); where an upstream edit is unavoidable (ADR 43 SP/MP
build-flag), mark with `// Ladruno ADR43` and record in [[LEDGER_vanilla_files]] **in the same
PR**. Owner: each phase.

---

## 6. Milestones & exit gates (program view; details in each ADR)

| Gate (phase) | Proven by |
|---|---|
| **G-A** (P-A) complex modal correct (serial) | 2-DOF non-classical closed-form complex modes; base-isolated stick model; vs `scipy.linalg.eig` on projected matrices; vs log-dec from a decay history |
| **G-B** (P-B) FEAST serial = ARPACK | Same eigenpairs on a medium model; band-targeting counts *all* modes in `[f₁,f₂]`; Sturm/inertia completeness |
| **G-C** (P-C) distributed MP solve = serial | **MET by #532** (L3-only MP distributed `dmumps` inner solve): `mpiexec -n 2/4` distributed spectrum == serial `-rci` oracle (1.5–3.3e-13, MAC≥0.999), lockstep + ΦᵀMΦ=I. **Scope narrowed** to MP-vs-serial (SP de-scoped); "SP vs MP" identity no longer a gate criterion. Optional L2 quadrature scaling deferred (demand-driven) |
| **G-D** (P-D) buckling/prestressed | Euler `P_cr=π²EI/(KL)²` (multiple BCs), plate buckling, string-tension frequency shift |
| **G-E** (P-E) complex modes at scale | Complex-FEAST eigenpairs (λ, ψ up to scaling/phase) match ADR 46's serial projection result on a non-classically-damped model within the projection's own convergence envelope; conjugate-pair completeness (every λ paired with λ̄, no orphans in the contour); serial complex-contour == parallel complex-contour spectrum |
| **G-F** (P-F) frequency domain | SDOF/2-DOF FRF vs analytic; modal transient == direct Newmark (linear); random RMS vs Monte-Carlo; CQC reproduces a published closely-spaced-modes example (ADR 44 P1c); `-combine`-absent path **byte-identical** to current `responseSpectrumAnalysis` (ADR 44 P1d — protects the D4 vanilla surface) |

Each gate is the ship gate for its phase (P-G's gate is defined when its ADR is drafted). No phase
merges without its gate green.

---

## 7. Program risk register

| # | Risk | Mitigation |
|---|---|---|
| R1 | **FEAST build dependency** — **RESOLVED**: no vendoring; MKL FEAST stays, gated on two OpenSees-side P3 fixes (D2 spike, 2026-07-07) | P3a comm-split plumbing fix in `MumpsParallelSolver`; P3b symmetric 2n×2n block-real inner SOE (existing `dmumps`, zero new dep); both unit-tested before P3c orchestration |
| R2 | **SP/MP build-flag surgery** (`_PARALLEL_PROCESSING` vs `_PARALLEL_INTERPRETERS`) — **RETIRED 2026-07-08**: P4 de-scoped, so the surgery is never performed | MP is the single blessed parallel config; distributed FEAST inner solve compiles into the MP targets only behind the `LadrunoFeastInnerSolve` seam (#532), upstream guards untouched; SP unsupported for FEAST — no gate change to isolate |
| R3 | **MKL ABI** (MPI integer size, threading model) | Pin in [[Ladruno_internal]] compilation journal at P-C; mirror existing MKL usage |
| R4 | **Assembled-`C` scope creep** (D1) | v1 = exactly `getDamp()`+Rayleigh; warn on others; widen only with sign-off |
| R5 | **classTag OR ADR-number collision** during the open window (it has happened **three times**: 41→42 shift, 40→46 shift, and "candidate 47" taken by [[47_ladruno_contact_deferrals_adr]]) | Re-verify 33019/33021/33022/33023/33024 vs fresh `ladruno` HEAD before each merge ([[feedback_stale_pr_ledger_ci]]); **never pre-assign ADR numbers to unwritten ADRs** — take the next free number at drafting time |
| R6 | **Adversarial-gate need** | Run the full multi-agent gate for the *novel math* phases (P-A complex projection, P-C parallel); skip for mechanical phases per [[feedback_adversarial_gate_when]] |

---

## 8. Build / CI / ledger obligations across the family

- **classTag reservations** (this PR, docs-only): 33019 (46), 33021 (42), 33022/33023 (43),
  33024 (44) — RESERVED in [[LEDGER_implementations]]; **enter `SRC/classTags.h` only when each
  implementation merges**. 33020 skipped — now occupied by the shipped LadrunoSolidShell
  ([[66_ladruno_solidshell_adr]], `classTags.h:939`).
- **Per-shipping-phase:** add the `SRC/classTags.h` define; flip the ledger row RESERVED→active;
  stamp the LADRUNO header on new source ([[feedback_always_stamp_header]]); add a
  `Ladruno_scripts/banner_features.txt` line **only when the feature actually ships**; record any
  upstream edit in [[LEDGER_vanilla_files]] (for ADR 43: `getNewEigenSOE` broker switch + the two
  `eigen` parsers + `BasicAnalysisBuilder` — **not** the parallel mains, since P4 is de-scoped and
  the distributed inner solve sits behind the `LadrunoFeastInnerSolve` seam; `ResponseSpectrumAnalysis`
  for ADR 44 if D4 takes the edit path).
- **Build deps:** ADR 43 links MKL FEAST (already on the link line) — no PFEAST vendoring (D2 resolved);
  note the MP distributed-solve recipe in [[Ladruno_internal]] compilation journal.
- **PRs:** one logical phase per PR; base on `ladruno`; verify branch is current before each push
  ([[feedback_stranded_commits_after_automerge]]).

---

## 9. Decision log / provenance

- **2026-06-21/22.** Family scoped from a cross-code deep-theory dive (Abaqus / Kratos / LS-DYNA
  skills + LS-DYNA manuals), four parallel research dossiers, and a load-bearing analysis
  ([[modal_gap_study/00_SYNTHESIS]] §6). Convergent findings: complex modal = project-to-real-modes
  QZ (Abaqus + LS-DYNA both avoid the full quadratic); parallel eigen = many independent distributed
  *linear* solves (Kratos PFEAST + LS-DYNA distributed factorization both avoid distributed Krylov).
- **Numbering rebased** onto `ladruno` after ADR 41 (mortar/ALM contact) landed: family shifted
  41→42, 42→43, 43→44; this roadmap took 45. **Then a concurrent `40_ladruno_performance_adr`
  landed on `ladruno`, colliding on 40 → complex-modal moved 40→46 and the ROM candidate 46→47.**
- This ADR + the four feature ADRs + the theory study ship together as a **docs-only PR**
  ([#351](https://github.com/nmorabowen/OpenSees/pull/351)); no `SRC/` change, no banner line yet.
- **2026-07-08 — P-C shipped + FEAST runway closed.** ADR 43 P3 landed across #524 (P3a
  comm-split) / #527 (P3b block-real kernel) / #530 (P3c-serial RCI) / #532 (P3c-MPI-L3 distributed
  `dmumps` inner solve). **Decision: MP (`OpenSeesMP`/PyMP) is the single blessed parallel
  configuration; SP (`_PARALLEL_PROCESSING`) is explicitly unsupported for FEAST; the P4 SP/MP
  build-flag unification is DE-SCOPED** (not deferred) — the distributed solve already runs in MP,
  so the "impossible-today" partitioned-SP + distributed-solve combination P4 targeted is not
  wanted. Gate **G-C is MET by #532** rather than pending an SP-unification PR; program risk R2
  retired. The only open ADR-43 item is the optional, demand-driven L2 quadrature-parallel
  multiplier (ADR 43 §9 R0). The FEAST substrate line is functionally complete.
- **2026-07-06 — adversarial review applied** (source-verified against `SRC/` + all four member
  ADRs; all 8 code-state premises held, 2 needed rewording):
  - **Added the missing G-E gate** — P-E (complex contours at scale, the riskiest math phase) had
    no ship criterion; the gate table also skipped from G-D to G-F.
  - **G-F extended** with ADR 44's P1c (CQC closely-spaced) + P1d (byte-identical `-combine`-absent
    regression) so the program gate protects the one sensitive vanilla surface.
  - **D1 reworded**: "no `addC`" was refutable — `FE_Element::addCtoTang()` exists
    (`FE_Element.cpp:383`) but only folds C into the combined effective tangent; the true gap is
    the absence of an *isolated* assembled-C export. Also pinned that the closed form is
    Rayleigh-only (undamped-mode basis); dampers need full `ΦᵀCΦ`.
  - **D3/D4 strengthened with file:line ground truth**: `eigen` shift hard-zeroed at
    `OpenSeesCommands.cpp:2187` with complete unreachable plumbing in `ArpackSOE`; upstream
    `ResponseSpectrumAnalysis` performs **no** SRSS/CQC (defers to user) — ADR 44 *adds*
    combination, it doesn't edit it.
  - **Staleness swept**: LadrunoSolidShell shipped and occupies 33020 under its own
    [[66_ladruno_solidshell_adr]] (was cited via ADR 19); "candidate ADR 47" collided with the
    contact-deferrals ADR → ROM candidate is now number-TBD-at-drafting (R5 widened to cover
    ADR-number collisions); fixed ADR 43's renumbering leftovers ("46/41/43" banner, "ADR 40"
    re-host reference, `// Ladruno ADR42` marker → ADR43).

---

## 10. Cross-references

[[modal_gap_study/00_SYNTHESIS]] (theory + load-bearing) · [[46_ladruno_complex_modal_adr]] ·
[[42_ladruno_buckling_adr]] · [[43_ladruno_feast_eigensolver_adr]] ·
[[44_ladruno_frequency_domain_adr]] · ROM/Craig–Bampton candidate (ADR number TBD, §2 note ²) ·
[[Ladruno_explicit_roadmap]] (umbrella-ADR precedent).
