---
title: "ADR 45 — Modal-analysis family: implementation roadmap & sequencing plan"
project: Ladruno
type: ADR / program plan (umbrella over ADRs 40/42/43/44)
status: draft — planning, NO code
priority: high
owner: nmora
related:
  - "[[modal_gap_study/00_SYNTHESIS]]"            # the cross-code theory + §6 load-bearing assessment this plan executes
  - "[[40_ladruno_complex_modal_adr]]"            # member: complex/state-space modal (33019)
  - "[[42_ladruno_buckling_adr]]"                 # member: prestressed modal + linear buckling (33021)
  - "[[43_ladruno_feast_eigensolver_adr]]"        # member: band-targeted parallel eigensolver (33022/33023) — the substrate
  - "[[44_ladruno_frequency_domain_adr]]"         # member: FRF/SSD/random + modal transient (33024)
  - "[[modal_gap_study/01_opensees_current_state]]" # ground-truth file:line audit
  - "[[19_ladruno_rc_shell_adr]]"                 # classTag boundary: 33020 reserved for LadrunoSolidShell (we skip it)
  - "[[Ladruno_explicit_roadmap]]"                # sibling precedent: a program roadmap spanning several ADRs
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags: [adr, program-plan, roadmap, modal, eigen, complex-modes, buckling, feast, parallel, frequency-domain, sequencing]
updated: 2026-06-22
---

# ADR 45 — Modal-analysis family: implementation roadmap & sequencing plan

**Status:** draft. **Planning only — no code.** This is the *umbrella ADR* governing the rollout of
the modal-analysis family (ADRs **40 / 42 / 43 / 44**, candidate **46**). The per-feature ADRs
specify *what each feature is and how it works*; this ADR decides *what we build, in what order, and
which cross-cutting decisions gate the program*. Theory lives in
[[modal_gap_study/00_SYNTHESIS]] (+ the four code dossiers).

> [!info] One-line thesis
> Bring OpenSees modal analysis to commercial-code parity (complex/damped modes, prestressed +
> buckling, robust band-targeted **parallel** eigensolving, frequency-domain response) by building
> **one shared eigensolver substrate** (ADR 43) and a **cheap serial complex-modal proof** (ADR 40)
> first, then layering the opportunistic deliverables (42, 44) and a future ROM extension (46).

---

## 1. Why a program ADR (not just the four feature ADRs)

The four feature ADRs were scoped independently, but they are **not independent to build**:

- They **share a substrate** — every one rides the eigensolver, and three of four ride an
  *assembled-operator* contract (`M`, `C`, `K`) that OpenSees only partially exposes today.
- They **share cross-cutting decisions** — the assembled-`C` accessor (40, reused by 44 and FEAST
  complex), the `-shift` exposure on the `eigen` command (42, touches shared code), the
  **MKL-FEAST-vs-vendored-PFEAST** build-dependency call (43), and the vanilla-footprint policy
  (CQC edit in 44, SP/MP build-flag surgery in 43).
- They have a **non-obvious optimal order** — the cheapest ADR (40) and the most load-bearing ADR
  (43) are *different* ADRs, so "build cheapest first" and "build foundation first" disagree. That
  tension needs a decision, recorded once, here.

Precedent: [[Ladruno_explicit_roadmap]] plays the same umbrella role for the explicit-dynamics ADRs.

---

## 2. The family at a glance

| ADR | Feature | classTag | Effort | Load-bearing (synthesis §6) |
|---|---|---|---|---|
| **40** | `LadrunoComplexEigen` — complex/state-space modal (non-classical damping) | 33019 | **S** | Domain-enabling (SSI/DRM/isolation/dampers) |
| **42** | `LadrunoBuckle` — prestressed modal + linear buckling | 33021 | **S–M** | Standalone analysis (modest) |
| **43** | `FeastEigenSOE`/`FeastEigenSolver` — band-targeted **parallel** eigensolver | 33022/33023 | **L** | **Substrate (highest)** + general SP/MP fix |
| **44** | `LadrunoModalResponse` — FRF/SSD/random + modal transient | 33024 | **M** | Deliverable (none downstream) |
| *(46)* | *ROM / Craig–Bampton substructuring (candidate, future)* | *TBD* | *L* | *Rides the family; biggest forward unlock* |

`33020` is **deliberately skipped** — reserved for `LadrunoSolidShell` ([[19_ladruno_rc_shell_adr]]).
ADR **41** is the unrelated mortar/ALM contact ADR already on `ladruno`.

---

## 3. Dependency graph (what blocks what)

```
                ADR 43  FEAST eigensolver + SP/MP fix          ← SUBSTRATE
                  │   (band-target · Sturm · MPI via MumpsParallelSOE)
   ┌──────────────┼───────────────────────────┐
   │ (serial eigen already enough)             │ (parallel + complex contours)
   ▼              ▼                             ▼
 ADR 40        ADR 42                        ADR 40 @ scale
 complex       buckling/Kg                   (re-host via complex FEAST)
 (serial OK)   (serial OK)                        │
   │              │                               │
   └──────┬───────┘                               │
          ▼                                       ▼
        ADR 44  frequency domain (FRF/SSD/random) — consumes 40 + the eigensolver
          │
          ▼
        ADR 46  ROM / Craig–Bampton (future) — needs trustworthy basis + parallel eigen
```

Key reading: **40, 42, 44 can each ship on the *existing* serial ARPACK `eigen`.** ADR 43 is not a
*blocker* for a first version of any of them — it is the **scale + robustness + parallel** upgrade
that makes them trustworthy on large/partitioned models and re-hosts 40's complex case via complex
contours.

---

## 4. Sequencing decision

Two defensible orders (synthesis §5 vs §6.5):

- **Cheapest-first proof:** 40 → 42 → 43 → 44. Validates direction with the smallest spend.
- **Unlock-the-most:** 43 → (40, 42) → 44. Builds the load-bearing substrate first.

**Decision — a hybrid that front-loads the cheap proof, then the substrate:**

| Phase | Work | Why here | Depends on |
|---|---|---|---|
| **P-A** | **ADR 40 P0–P3** — complex modal, *serial*, on existing `eigen` | Cheapest, directly serves the research portfolio (isolation/dampers/SSI); de-risks the projection+QZ approach before any big build | existing `eigen`, `modalProperties` |
| **P-B** | **ADR 43 P1–P2** — serial **MKL-FEAST** eigensolver (band-target + Sturm) | The substrate; **zero new build dep** (MKL already linked); validated vs ARPACK | MKL |
| **P-C** | **ADR 43 P3–P4** — **parallel** (MPI per-contour via `MumpsParallelSOE`) + **SP/MP unification** | The strategic payoff: large-model modal + the *general* parallel-composition fix | P-B, `MumpsParallelSOE`, PartitionedDomain |
| **P-D** | **ADR 42** — prestressed modal + linear buckling | Opportunistic; rides serial eigen, gains band/Sturm from P-B; can jump ahead of P-C if a project needs it | corot/PDelta Kg, `eigen` |
| **P-E** | **ADR 43 P5** — complex contours (re-host ADR 40 at scale) | Unifies the complex case onto the parallel substrate | P-A, P-C |
| **P-F** | **ADR 44** — frequency domain (FRF/SSD/random, modal transient) | Deliverable layer; build when a project asks | P-A, eigen |
| *(P-G)* | *ADR 46 — ROM / Craig–Bampton* | *Future; the biggest forward unlock* | *whole family* |

**Rationale.** P-A buys confidence + an immediately useful research deliverable for ~S effort. P-B
starts the substrate at *zero dependency cost* (the most-load-bearing work that carries no build
risk). P-C is the one large, build-risky step (MPI + build-flag surgery) and is deliberately
sequenced after the substrate is proven serially. P-D/P-F are demand-driven.

---

## 5. Cross-cutting decisions (program-level — decide once, here)

These appear in multiple ADRs; resolving them at the program level prevents divergent choices.

### D1 — Assembled-`C` accessor (from ADR 40; reused by 44, 43-complex)
OpenSees has **no assembled damping-matrix accessor** (audit: no `addC`/`formEleTangC`; only
`Element::getDamp()`). **Decision:** build a single `LadrunoDampingAssembler` path in P-A (Rayleigh
projected in closed form as `αM̃+βK̃`; element/material dampers via `getDamp()`), and **reuse it
verbatim** in ADR 44 and ADR 43's complex contours. v1 contract = *exactly* `getDamp()` + Rayleigh;
**warn** (don't silently absorb) when `modalDamping`/HHT numerical damping is active. Owner: P-A.

### D2 — MKL-FEAST vs vendored PFEAST (from ADR 43 — the big build call)
Intel MKL's Extended Eigensolver **is** FEAST and the build already links MKL → **serial P-B costs
zero new dependency.** Open question is whether MKL's RCI form lets each contour solve run as an
OpenSees `MumpsParallelSOE` solve on an arbitrary MPI sub-communicator. **Decision:** stay on
**MKL-FEAST for P-B**; **defer the MKL-vs-vendored-PFEAST call to the P-C gate** (P-B/P-A are
identical either way). If MKL's distributed layer resists user sub-communicators, vendor PFEAST 4.0
into `OTHER/FEAST/` (ARPACK precedent in `OTHER/ARPACK/`). Owner: P-C gate.

### D3 — `eigen` `-shift` exposure (from ADR 42)
Buckling needs a non-zero shift (`ΔK` indefinite) but `eigen` currently hard-zeros the shift.
**Decision:** expose `-shift` on the new `buckling` command in P-D; **only** un-hard-zero the shared
`eigen` path if P-D shows a concrete need (minimize vanilla touch). Owner: P-D.

### D4 — Vanilla-footprint policy (from ADR 44 + ADR 43)
ADR 44's CQC/SRSS edits Petracca's upstream `ResponseSpectrumAnalysis`; ADR 43 touches the
`EigenSOE` base + parallel mains (the SP/MP gate). **Decision:** prefer **fork-authored siblings**
that *read* committed state where feasible (ADR 44 path); where an upstream edit is unavoidable
(ADR 43 SP/MP build-flag), mark with `// Ladruno` and record in [[LEDGER_vanilla_files]] **in the
same PR**. Owner: each phase.

---

## 6. Milestones & exit gates (program view; details in each ADR)

| Gate | Proven by |
|---|---|
| **G-A** complex modal correct (serial) | 2-DOF non-classical closed-form complex modes; base-isolated stick model; vs `scipy.linalg.eig` on projected matrices; vs log-dec from a decay history |
| **G-B** FEAST serial = ARPACK | Same eigenpairs on a medium model; band-targeting counts *all* modes in `[f₁,f₂]`; Sturm/inertia completeness |
| **G-C** parallel scaling + SP/MP unified | Strong/weak scaling on a partitioned model; identical spectrum SP vs MP vs serial |
| **G-D** buckling/prestressed | Euler `P_cr=π²EI/(KL)²` (multiple BCs), plate buckling, string-tension frequency shift |
| **G-F** frequency domain | SDOF/2-DOF FRF vs analytic; modal transient == direct Newmark (linear); random RMS vs Monte-Carlo |

Each gate is the ship gate for its phase. No phase merges without its gate green.

---

## 7. Program risk register

| # | Risk | Mitigation |
|---|---|---|
| R1 | **FEAST build dependency** (MKL RCI sub-communicator limits → vendoring PFEAST) | Defer to P-C gate; P-A/P-B unaffected; ARPACK vendoring precedent |
| R2 | **SP/MP build-flag surgery** (`_PARALLEL_PROCESSING` vs `_PARALLEL_INTERPRETERS`) is general-infra risk | Land serially first (P-B); isolate the gate change; test SP==MP==serial spectrum |
| R3 | **MKL ABI** (MPI integer size, threading model) | Pin in [[Ladruno_internal]] compilation journal at P-C; mirror existing MKL usage |
| R4 | **Assembled-`C` scope creep** (D1) | v1 = exactly `getDamp()`+Rayleigh; warn on others; widen only with sign-off |
| R5 | **classTag collision** during the open window | Re-verify 33019/33021/33022/33023/33024 vs fresh `ladruno` HEAD before each merge ([[feedback_stale_pr_ledger_ci]]) |
| R6 | **Adversarial-gate need** | Run the full multi-agent gate for the *novel math* phases (P-A complex projection, P-C parallel); skip for mechanical phases per [[feedback_adversarial_gate_when]] |

---

## 8. Build / CI / ledger obligations across the family

- **classTag reservations** (this PR, docs-only): 33019 (40), 33021 (42), 33022/33023 (43),
  33024 (44) — RESERVED in [[LEDGER_implementations]]; **enter `SRC/classTags.h` only when each
  implementation merges**. 33020 skipped (LadrunoSolidShell).
- **Per-shipping-phase:** add the `SRC/classTags.h` define; flip the ledger row RESERVED→active;
  stamp the LADRUNO header on new source ([[feedback_always_stamp_header]]); add a
  `Ladruno_scripts/banner_features.txt` line **only when the feature actually ships**; record any
  upstream edit in [[LEDGER_vanilla_files]] (expect: `EigenSOE` base + parallel mains for ADR 43;
  `ResponseSpectrumAnalysis` for ADR 44 if D4 takes the edit path).
- **Build deps:** ADR 43 P-C may add FEAST/PFEAST — note in [[Ladruno_internal]] compilation journal.
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
  41→42, 42→43, 43→44; ROM candidate 45→**46** (this ADR took 45).
- This ADR + the four feature ADRs + the theory study ship together as a **docs-only PR**
  ([#351](https://github.com/nmorabowen/OpenSees/pull/351)); no `SRC/` change, no banner line yet.

---

## 10. Cross-references

[[modal_gap_study/00_SYNTHESIS]] (theory + load-bearing) · [[40_ladruno_complex_modal_adr]] ·
[[42_ladruno_buckling_adr]] · [[43_ladruno_feast_eigensolver_adr]] ·
[[44_ladruno_frequency_domain_adr]] · candidate ADR 46 (ROM/Craig–Bampton) ·
[[Ladruno_explicit_roadmap]] (umbrella-ADR precedent).
