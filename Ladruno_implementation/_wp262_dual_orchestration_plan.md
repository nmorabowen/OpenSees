---
title: "Dual-number program — orchestration plan (work packages, lanes, models, protocol)"
project: dual-numbers-opensees
status: approved by the owner 2026-11-28 ("Sí, apruebo") with the §8 defaults — see §8
orchestrator: "the authoring session (Opus 5.5), the §8 default the owner approved on 2026-11-28"
companion: "research bench gitAPE/dual-numbers-opensees: research/00_overview.md (vision, levels L0-L3), 04_porous_biot.md, 05_downsides.md; operating layer .claude/skills/ladruno-dual"
updated: 2026-11-28
---

# Dual-number program: orchestration plan

The research bench (this folder) settled **what** dual numbers are for in the Ladruno fork
(`research/00_overview.md`). WP-262 proved **how**:
- the `LadrunoDual.h` header;
- the J2 kernels templated, bit-identical;
- tangent oracles at round-off;
- the first DDM sensitivities, `LadrunoUniaxialJ2` in a fiber beam, matching FD at 3e-10.

This plan turns the roadmap into work packages that several agents can run in parallel without
colliding, and names where the live state lives. It follows the fork's program precedent
(MPM, `Ladruno_implementation/_wp232_orchestration_plan.md`, `ladruno-mpm` skill, program #1075).

## 0. Where things live

| What | Home |
|---|---|
| Program status, WP claims, decisions | a fork **`program` issue**, "Program: dual numbers (WP-262…)" — **to open (§8 item 2)** |
| One WP's state and handoffs | its **draft PR** (`### Handoff — WP-<n> — <date> — <model> @ <effort>`) |
| The fork-wide view | board row on **#937** (already there for WP-262) |
| The operating layer for sessions | a fork skill **`ladruno-dual`** (`.claude/skills/ladruno-dual/`), written in wave 0 |
| Research, measurements, prototypes, this plan | this habitat (`research/`, `prototypes/`, `plan/`) |
| apeGmsh slices | apeGmsh's own program protocol (T coordinator, its board); this plan only recommends |

**WP numbers are claimed on the program issue, never guessed.** Below, the WPs carry
**slugs** (`D-RC`, `D-C3D`, …); each gets its number when the orchestrator claims it at kickoff.

## 1. The map

```
WAVE 0  WP-262 (in flight: header + J2 + UniaxialJ2 DDM) ── merged first
        + program issue + ladruno-dual skill + M0 (MSVC bench)
          │
WAVE 1  (parallel; disjoint files)
          ├─ D-RC      RC returnMap3D Dual<6> tangent + Newton-iteration measurement   ──► DECISION D2
          ├─ D-C3D     Concrete3D FD driver gradients + RebarBuckling dfactor → dual
          ├─ D-IMPLEX  generic IMPL-EX tangent-identity gate + non-zero mutation IDENT floor
          ├─ D-UP1     LadrunoUP kernel templated + Dual<32> hypo Jacobian oracle (dropped blocks measured)
          └─ D-SENS2   DDM for the other fork uniaxials (CohesiveHinge, BondSlip, RebarBuckling)
          │
WAVE 2  (gated)
          ├─ D-HS      HS regime residuals templated; Dual IFT oracle for sens[2][4]      (after M1)
          ├─ D-SAN     SANISAND kernel extraction + substep-chain tangent  (H)            (after D2 = go)
          ├─ D-COROT   SolidTransformation corot seam consistent tangent (Dual<ndof>)     (after D-UP1 lesson)
          ├─ D-MORTAR  LadrunoMortarKernel linearizePair via Dual<24>
          ├─ D-UPDDM   LadrunoUP getResistingForceSensitivity / commitSensitivity
          └─ G1        apeGmsh DDM plumbing + refusal catalogue (apeGmsh ADR)            (after 262 is released)
          │
WAVE 3  (decided at wave-2 close)
          ├─ D-FSL     adaptive fixed-stress L from a dual K_dr                           (needs D-HS or D-SAN)
          ├─ D-UNSAT   unsaturated Biot physics mode (own ADR, H)                         (DECISION D4)
          ├─ D-LINT    quirk-lint rule fd-in-production-tangent                           (after D-C3D)
          ├─ D-FS2D    ladruno_fs2d / hypo element oracles
          ├─ G2        material-point harness from Python                                (needs WP-248 ndTest in Python)
          └─ G3        calibration with exact Jacobians (habitats, not apeGmsh core)
```

Cross-cutting measurement work packages:
- **M0:** the `Dual<N>` benchmarks under MSVC `cl.exe /O2`. All the numbers so far are g++.
- **M1:** a Newton-iteration and wall-time harness on fork models, reused by D-RC, D-COROT and
  D-UP1. apeGmsh's Pardiso factorisation counts (its ADR 0106) are the reference metric.

## 2. Work packages

| Slug | Level | Content | Files (fork) | Agent | `/code-review` | Done when |
|---|---|---|---|---|---|---|
| **WP-262** | L0+L3 | `LadrunoDual.h`; J2 + Lemaitre kernels templated; dual tangent oracles; `LadrunoUniaxialJ2Kernel.h` + DDM | in PR #1172 | this session | `high` | Zone-A green; out-of-family critic read; owner merges |
| **M0** | — | `Dual<N>` benches (N = 1, 6, 24, 32) under MSVC `/O2` and g++, inside the fork's toolchain | habitat `prototypes/` | Sonnet · medium | — | a table of ns per call by compiler, committed with the build flags |
| **M1** | — | iteration / wall-time harness: runs a deck twice (tangent A vs B), reports Newton iterations per step, factorisations, total time, with the build hash | `Ladruno_scripts/` or `tests/_testbed/` (new file) | Sonnet · medium | `low` | reproduces a known count on an existing deck |
| **D-RC** | L1 | `LadrunoRCKernel.h` `returnMap3D` templated; `-tangent dual` (new `tanMode`) replacing `-numericalTangent`, valid under IMPL-EX; oracle vs the secant | `LadrunoRCKernel.h`, `LadrunoRCConcrete.cpp`, `LadrunoRCFiniteStrain.cpp`, `LadrunoOptSpec.h` | Opus · high | `high` | bit-identity of the double path; dual = FD at round-off; **M1 measurement on an RC panel and wall deck, before vs after** |
| **D-C3D** | L1 | `dscalarDsig` / `dscalarDsigAxis` driver gradients → `Dual<6>`; the axis case keeps its branch; `LadrunoRebarBuckling` `dfactor/dε` FD → `Dual<1>` | `LadrunoConcrete3DKernel.h`, `LadrunoRebarBuckling.cpp` | Opus · high | `high` | replay lanes bit-exact where the FD is gone; damaged tangent agrees in MPa/kPa/Pa with no step constant; per-call time vs PR-877 |
| **D-IMPLEX** | L0 | one harness: per IMPL-EX material, dual derivative of the extrapolated map (clock latched) = shipped secant; `mutation_gate.py` IDENT floor > 0 | `tests/` (new), `Ladruno_scripts/mutation_gate.py`, `LadrunoImplexClock.h` (read-only) | Opus · high | `medium` | a planted live `df/dε` fails the gate; every family's IDENT floor > 0 |
| **D-UP1** | L0 | `LadrunoUPKernel.h` + the hypo GP residual templated; `Dual<32>` H8 element Jacobian oracle | `SRC/element/ladrunoUP/*` | Opus · high | `high` | size of each dropped block reported; iteration count full vs dominant-term on a large-strain consolidation (M1) |
| **D-SENS2** | L3 | DDM for `LadrunoCohesiveHinge`, `LadrunoBondSlip`, `LadrunoRebarBuckling` by the WP-262 recipe | those three classes (+ kernels extracted as needed) | Opus · medium (pattern follows 262) | `high` | per material: kernel bit-identity, DDM-through-history = FD, NOHIST fails; a beam or zeroLength end-to-end vs FD |
| **D-HS** | L0→L2 | HS regime residuals templated; `Dual<2+4>` IFT oracle replacing the hand `sens[2][4]` check | `LadrunoHSKernel.h` | Opus · high | `high` | oracle = analytic over the 5 regimes + tension corner + WP-210 states |
| **D-SAN** | L1 | extract `LadrunoSANISANDKernel.h` from vanilla-derived MD code; dual through ModifiedEuler / SAS-ME; `Dual<19>` CPPM Jacobian | `LadrunoSANISAND*` (fork), MD untouched | Opus · high, one agent across rounds | `high` + out-of-family verdict | algorithmic tangent of the chain; iteration count vs TanType 2 (M1) |
| **D-COROT** | L1 | `SolidTransformationCorot` consistent tangent (the deferred v2.1 terms) via `Dual<ndof>` | `SolidTransformationCorot.cpp` (fork seam) | Opus · high | `high` | quadratic convergence on the corot battery (M1) |
| **D-MORTAR** | L1 | `linearizePair` (stub returns −1) via `Dual<24>` | `LadrunoMortarKernel.h` | Opus · high | `high` | FD = dual on the contact patch tests; Newton count on a mortar lane |
| **D-UPDDM** | L3 | `LadrunoUP` element sensitivities by parameter-seeded duals (k̄, α, n₀, K_f) | `SRC/element/ladrunoUP/*` | Opus · high | `high` | Terzaghi `∂p/∂k̄` vs closed form and FD |
| **G1** | L3 | apeGmsh: `ops.sensitivity.*`, `-NS` on the Ladruno recorder primitive, reader entries, archive zone, `SENSITIVITY_CAPABLE` refusal table | apeGmsh | Opus · high (ADR), Sonnet · medium (plumbing) | apeGmsh rules | frame model `dU/dθ` = FD through the bridge; a soil model is refused, not zeroed |

Waves 3 rows (D-FSL, D-UNSAT, D-LINT, D-FS2D, G2, G3) are scoped at the wave-2 close from what
waves 1–2 measured; D-UNSAT starts with its own ADR.

## 3. Parallelism rules

- **One writer per file set.** Wave 1 is chosen so the five WPs touch disjoint files (RC, Concrete3D +
  RebarBuckling, tests + mutation gate, LadrunoUP, three uniaxials). **D-SENS2 and D-C3D both touch
  `LadrunoRebarBuckling.cpp`**: D-C3D owns it, and D-SENS2 does RebarBuckling last, after D-C3D merges.
- **Shared header.** `LadrunoDual.h` is owned by WP-262 until it merges. After that, a change to it is
  its own small PR, never part of a feature WP, so two workers never edit it at once. Workers who need
  a new primitive ask the orchestrator.
- **One worktree per WP**, branch `wp/<n>-<slug>` from fresh `origin/ladruno`; never the shared checkout.
- **Builds are launched by the orchestrator** in the worker's worktree. A subagent's own background build
  never wakes it. Up to two concurrent `build.bat` runs: this machine already slows when a build runs
  beside pytest.

## 4. Roles

| Role | Model · effort | Why |
|---|---|---|
| Orchestrator: sequencing, claims, builds, verification, handoffs | **this session (Opus 5.5)** or a fresh **Fable · high**; **max** at D2 / D3 and when adjudicating a review | fork precedent (ADR-90, ADR-94, WP-232) |
| Kernel templating + oracle authors (D-RC, D-C3D, D-UP1, D-HS, D-SAN, D-COROT, D-MORTAR, D-UPDDM) | **Opus · high**, one agent per WP continued across review rounds | bit-identity and tangent algebra are correctness-critical |
| Recipe ports (D-SENS2) | **Opus · medium** | follows the WP-262 recipe, but each material has its own branch structure |
| Harnesses (M0, M1) | **Sonnet · medium** | written once, run many times; the orchestrator launches the runs |
| Critic at each ready-flip | the **other family**: Fable · high for Opus authors | critics must not share the author's context |
| Read-only sweeps | `Explore`, "very thorough" | keeps the orchestrator's context lean |
| Out-of-family verdict (ADR-87 D4) | run by the owner at D-SAN and at the first L1 production change (D-RC) | the citable review |

## 5. Orchestrator protocol

1. **Boot:**
   - read #937, the program issue, and the WP's draft PR (its last `### Handoff` is the resume point);
   - check `ListAgents` and the `in-flight` label for a live duplicate;
   - `git fetch origin`, and read code from `origin/ladruno`.
2. **Claim:** post the next WP number on the program issue, open the branch and worktree, and open the draft PR on day one.
3. **Dispatch:**
   - brief from the skill's `references/wp-briefs.md`;
   - `subagent_type: general-purpose`, the roster's `model` and `effort`, `isolation: "worktree"`;
   - independent workers go in one message;
   - add `in-flight`.
4. **Verify — never accept a worker's word:**
   - re-run the cheapest gates yourself: bit-identity, the oracle, the static gates;
   - open every binary probe with `ladrunoBuild()` = the worker's head;
   - a reported number names the build hash, the deck, the output path and the date.
5. **Review:** one out-of-family critic per WP at ready-flip, briefed with the diff and the plan row, not with your conclusions. Then `/code-review high` for C++.
6. **Land:**
   - flip to ready only when the row's done-when holds;
   - `gh workflow run ladruno.yml --ref <branch>` for Zone-A;
   - read CI with the `ccd_pr` tools and Auto-fix; never poll;
   - merges follow §8 item 5.
7. **Record:** a ledger fragment per WP; a quirk fragment for every trap; the board row and program issue refreshed with each handoff; this habitat's `research/` updated with each measurement.

## 6. Standing rules for every worker (verbatim in each brief)

1. **Template, don't fork.** Write the kernel once as `template <class T>` (plus `P` for constants
   when parameters are seeded). Keep the `double` entry point under its old name and signature.
2. **Bit-identity is the first gate.** Before any derivative work, show the `double` path unchanged:
   - either diff the existing C++ checks' `%.17g` output against the pre-WP headers (`git show origin/ladruno:`);
   - or compare bit for bit against a verbatim transcription of the old code.
3. **Branches read `value(x)`.** Tolerances are `double`. No implicit `Dual → double` anywhere.
4. **Local Newton loops get the IFT step:** one value-zeroed Newton step at the converged point
   (`zero_value`). Vector local solves use the converged Jacobian solve.
5. **Every oracle has a break-the-gate leg:** a mutant the oracle must flag (symmetrized tangent,
   dropped history, dropped term).
6. **Kinks keep their dedicated handling:**
   - apex, hydrostatic axis, repeated eigenvalues, `sqrt(0)`;
   - never push a dual through an eigen-solver at a degenerate state.
7. **Production replaces only an FD, a secant or a stub** (research/00 §10). An analytic tangent
   that passes the oracle stays.
8. **Zero vanilla files** unless the WP's row says otherwise; a vanilla touch needs a `// Ladruno`
   marker and a vanilla fragment.
9. **No `-march=native`** in g++ checks on Windows (MinGW AVX spill crash, prototypes README).
10. **Measure, don't promise:** a convergence or speed claim carries M1/M0 numbers.
11. Ledger fragments per WP; never edit `LEDGER_*.md`; quote the code you claim about.

## 7. Decision points

| ID | When | Question | Input |
|---|---|---|---|
| **D1** | now | Program scope: research/00 §10 (material point / element only; in-house `LadrunoDual.h`; production replaces only FD / secant / stub) | the research bench |
| **D2** | D-RC done | Does an exact tangent cut Newton iterations at model scale (RC panel and wall, hardening branch vs softening)? If **yes**: L1 proceeds (D-SAN, D-COROT, D-MORTAR as production). If **no**: L1 shrinks to FD and stub replacement; the program stands on L0 + L3 | M1 numbers |
| **D3** | after D2 | SANISAND kernel extraction (high effort, vanilla-derived code) — go / no-go | D2 + D-HS lessons |
| **D4** | wave-2 close | Unsaturated Biot as a new physics mode of `LadrunoUP` (own ADR) | D-UP1 measurement, demand |
| **D5** | wave-2 close | apeGmsh G1 goes ahead (needs a fork release with WP-262 and its `FORK_PIN` bump) | WP-262 merged + released |

## 8. Go / no-go items for the owner

Owner's "Sí, apruebo" on 2026-11-28 approved the plan with the defaults below.

1. ~~Approve §2 and §4.~~ Approved as written.
2. ~~Open the `program` issue claiming WP-262 onward.~~ Opened by the orchestrator; linked from #937.
3. ~~Choose the orchestrator session.~~ The authoring session (Opus 5.5).
4. ~~Approve the wave-0 bootstrap.~~ This plan and the `ladruno-dual` skill ride WP-262's branch (#1172).
5. Merges: **owner-only** (not delegated).
6. ~~Decide D1.~~ Approved: research/00 §10 scope. D2–D5 come back at their gates.

**Wave-1 claims (2026-11-28):** WP-263 D-RC, WP-264 D-C3D, WP-265 D-IMPLEX, WP-266 D-UP1,
WP-267 D-SENS2. M0 / M1: M1 is folded into WP-263 (its `rc_newton_count.py` harness) and WP-266
(its consolidation deck); M0 runs at the wave-1 close on the orchestrator's builds.
