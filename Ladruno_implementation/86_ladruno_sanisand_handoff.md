---
title: LadrunoSANISAND — session handoff (after PR-1 and PR-2)
project: Ladruno
status: active
owner: nmora
tags:
  - handoff
  - adr
  - material
  - sanisand
  - manzari-dafalias
---

# LadrunoSANISAND handoff

State after PR-1 ([#767](https://github.com/nmorabowen/OpenSees/pull/767), merged) and PR-2
([#768](https://github.com/nmorabowen/OpenSees/pull/768)). Spec is
[[86_ladruno_sanisand_adr]]; read it first, **including its correction boxes** — several of its
original numbers were refuted by measurement and corrected in place.

---

## 1. Read this before you touch anything

Four traps cost real time in the first two PRs. Each is cheap to avoid and expensive to rediscover.

**A plain `pytest` tests the WRONG BINARY.** `opensees_venv` carries `ladruno_opensees.pth`, which at
interpreter startup puts `C:\Program Files\Ladruno\OpenSees\bin` on `sys.path[0]` and eagerly
imports `opensees`. Measured: `test_manzari_safety_pack.py` is *2 failed, 1 passed* against the
installed engine and *3 passed* against the worktree one. Use a rebinding plugin and **confirm the
bound path on stderr every run**. Never modify anything under `C:\Program Files\Ladruno\`.

**Your binary can be stale, and green tests will not tell you.** A rebuild that finishes in under a
second did nothing. A header change to `ManzariDafalias.h` triggers a **~1500-TU, >10-minute**
rebuild. At the end of PR-2 the working binary was ahead of `SRC/` and nobody noticed until the
final rebuild recompiled 1492 TUs. **Rebuild from committed HEAD before recording any number.**

**CI does not run what you think it runs.** No workflow sets `--runslow` or `LADRUNO_RUN_SLOW`, so
anything marked `@pytest.mark.slow` **never runs in CI**. Before PR-1's fix, `_opensees_exe()` looked
only in `dist/bin`, which CI never creates — so the classic-Tcl gate skipped on every push and had
*never once executed*. Marking a test slow is equivalent to deleting it from CI.

**A skip and a pass look identical in a summary line.** `WORKFLOW_GOTCHAS` §7 records the G1 Tcl
suite reporting `OK (0 checks passed)` for months. When you arm a gate, expect it to go red — that
is evidence it exists.

---

## 2. The architecture, in one page

`LadrunoSANISAND : public ManzariDafalias`. Not a copy, not a wrapper — **inheritance**. When a deck
uses it, `integrate()` / `ModifiedEuler()` are `ManzariDafalias`'s functions, the same compiled code
vanilla runs. Our ~600 lines contain **no physics**.

What it changes: two `protected` members (`ManzariDafalias.h:162-163`) that ~30 sites in the physics
read at run time. Vanilla's `initialize()` hardcodes them; we let it run and then overwrite —
***win the last write***. That is also why G1 is bit-identical: given vanilla's constants, it is
vanilla's code on vanilla's data.

The six overrides exist only to stop the base taking those values back (re-`initialize()`,
`revertToStart` mid-analysis, a wire format with no slot for `p_r`, `getCopy` hardcoding vanilla
wrapper types at every Gauss point, and `Print` so the record can state what it ran).

**The ceiling.** Inheritance reaches `protected` data and `virtual` methods. `ModifiedEuler`,
`integrate` and `elastic_integrator` are **non-virtual** and dispatched through a pointer-to-member
explicitly typed to `ManzariDafalias`, so a subclass cannot substitute them — and you **cannot** fix
that by adding `virtual`, because `ManzariDafaliasRO` (`ManzariDafaliasRO.h:93-97`) *shadows*
`initialize()` and two `GetElasticModuli` overloads with identical signatures. Promoting them would
silently convert those shadows into overrides and change every existing Ramberg-Osgood deck.

### The three lanes (D9)

| lane | for | vanilla decks |
|---|---|---|
| **1. subclass** | anything reachable via protected data | untouched |
| **2. vanilla bugfix** | genuine **errors** | move — disclosed |
| **3. flag seam in vanilla** | **opinions**, and calibration-moving errors | bit-identical |

**The test is: error, or opinion?** A flag seam is a protected member on the base, default = vanilla,
read at the defect site, set only from a derived constructor. Two exist: `mHonorTolRInME` and
`mUseCurrentVoidRatioInG`. **Neither is wired to the subclass or the parser yet.**

**The fork tripwire.** The subclass stops being viable the day we change the *formulation* rather
than fix an error — the `alpha` re-set, `Elastic2Plastic`'s `M_c` repair, D5a's sigmoid shape, §7.3's
`m_e_init`. All are non-virtual and entangled, and changing them in vanilla imposes a modelling
opinion on every existing deck. **D5a is queued for PR-3, so this decision is one PR away.** Open it
consciously, not by drift.

---

## 3. What shipped

**PR-1** — the class, both wrappers, all six overrides, five registration sites, ND tags
33019/33020/33021, the battery, ledgers, manifest, banner, header stamps. Vanilla untouched.

**PR-2** — five vanilla commits (`mTolR` seam, clamp made observable, `D_factor` non-dimensionalised
at four sites, void-ratio interpolant fixed at three, clamp `p` self-consistency), one declined as
ambiguous, one flag seam, and a confine-first test deck.

### Constants (re-measured on an engine rebuilt from PR-2 HEAD)

| | ramp deck | confined deck |
|---|---|---|
| G1 equivalence | **0.0 bit-identical** | 0.0 bit-identical |
| non-vacuity (1.01 vs 0) | 5.579874e-02 | **9.325164e-02** |
| isolate-the-constant (1.01 vs 0.5) | 2.692933e-02 | **4.241113e-02** |
| PlaneStrain sensitivity | 3.177732e-02 | **1.091772e-01** |
| PlaneStrain mirror margin | 0.1866 | **19.157** |
| `\|eps_pl\|` at `p_r = 0` | **0.0 exactly** | 5.675e-03 |
| Outside Bounding / leg | **8** | **0** |

**The interpolant fix (PR-2 commit 4) moves `IntScheme 45` decks by ~0.6 %.** RK4 does not move (its loop reads a separate `CurVoidRatio`); **RK45 does** (11 in-loop reads of `NextVoidRatio`, 0 of `CurVoidRatio`). Measured by A/B rebuild on a NON-isochoric triaxial: reldiff **6.33e-3**. PR-2 initially reported zero because the confined test deck is exactly isochoric and the two anchors differ by exactly trace(d_eps) — a measurement that was true on its deck and falsely generalised. Caught by adversarial review before merge. **Lesson: an inertness measurement is only as general as the deck's strain path.**

Low-`p` gate: vanilla **+18.1258 %**, `p_r = 0` **+0.0754 %**, invariant `err·p_end ~= p_r` on both.
`-Presidual` plastic threshold: between **0.0230 and 0.0240 kPa**.

---

## 4. Two mechanisms you must know

**The radial freeze.** Every *elastic* step re-sets `NextAlpha = GetDevPart(NextStress)/p` with `p`
including `m_Presidual` (`ManzariDafalias.cpp:974-977`). On a perfectly proportional path with
`p_r = 0`, `s/p` is constant along the ray, so `alpha` is frozen at exactly the current stress ratio
and `f = ||s − p·alpha|| − √(2/3)·m·p` collapses to `−√(2/3)·m·p < 0` **forever**. The material never
yields, on any integration scheme. Not a defect of the `p_r = 0` default — `NTUASand02` ships
`p_r = 0` with identical code and is used for staged foundation work where the stress direction
changes constantly. It is a property of **proportional decks**, and a user can hit it with a
single-element calibration deck. Pinned by `test_radial_ramp_with_pr0_never_yields`.

**`M_c` inflation.** On the 3D ramp deck every leg — including vanilla's — starts the plastic stage
outside the bounding surface, and `Elastic2Plastic` silently raises `M_c` from 1.3309 to 1.99878, 8×
per leg. Any η/M^b measured there is against a surface the material stopped using. Zero on the
confined deck, and now asserted.

---

## 5. Owed work

**PR-3 candidates**
- Wire `-honorTolR 1` through to `mHonorTolRInME`; the parser currently rejects `1` by design.
- A **behavioural** `-Pmin` gate. Our default is `1e-3·P_atm`, **10× vanilla's**, is demonstrably not
  inert at low stress (via `Elastic2Plastic`'s `3*m_Pmin` reset), and every test pins `-Pmin` to
  vanilla's value. PR-2's clamp diagnostic is the observable that finally makes this writable.
- `LadrunoSANISAND3D::getCopy(void)` is uncovered (proven by mutation).
- **D5a** — should the `D_factor` sigmoid exist at all once `p_r = 0`? Modelling question. **This is
  the fork tripwire.**
- **§7.3** `m_e_init` — seam is open, unwired. Wiring it moves a calibrated quantity: `G0 = 264.32`
  was fitted against the frozen form.

**Subclass strengthening (proposed, not done)**
- A **silent-elasticity diagnostic**: override `commitState()`, and if the plastic stage is active
  while `mDGamma` has been identically zero for N consecutive commits, warn (throttled, process
  budget — *not* per instance; `getCopy` makes every Gauss point an instance). Catches the radial
  freeze for any user on any deck. Entirely subclass-side.
- Extend the echo/`Print`: `p_r = 0` costs more than cohesion — it also removes a **numerical
  regulariser** (the low-`p` leg fails at 400 steps where vanilla survives; needs 1200) and permits
  the radial degeneracy.

**Separate PR — `SAniSandMS`**
- §7.1: `numData = numArgs - 19` at `:134` and `numData -= 5` at `:143` silently drop `TolF`/`TolR`
  for 4–5 optionals. The ADR's line numbers are **exact**; a proposed "correction" to `:132`/`:141`
  was itself wrong.
- Its `D_factor` trigger is `0.001*P_atm`, **not** `0.05` — over that window the factor runs
  4.8e-4…1.0e-3, a ~1000× kill of dilatancy. A *shape* defect on top of the units one.
- Its `RungeKutta4:1550` carries the same interpolant bug.

**Do NOT re-propose (D8):** upstreaming. Considered and declined at fork level. Prof. Gorini is to be
informed; the upstream-PR call is **his**.

---

## 6. Method that worked

**Adversarial review by a different model found something real every single time** — a false
"self-arming" coverage claim, a CI gate that had never executed, three undefended overrides, a false
"8 events" claim, and two wrong premises of the orchestrator's. ADR §9's staffing note is not
ceremony.

**Mutation testing is what separates a gate from a decoration.** Three overrides had no defence at
all and read as covered. And the naive PlaneStrain assertion — "the stress should be compressive" —
*passes the mutant*, because the sign-flipped material still returns compression (−0.73 kPa). Only
the **deviator sign** discriminates. Reading the code would not have told you that.

**Name the deck, the leg and the constants when claiming a diagnostic fires.** Two agents reported
"8 events" about two unrelated diagnostics in one PR. Number collisions are how false claims survive
review.

**G1 is a pair, not a test.** Under the `getCopy(const char*)` mutation — every Gauss point silently
vanilla — G1 *passed*, comparing vanilla with vanilla. The non-vacuity companion caught it. The ADR
asked only for G1; the companion was added on suspicion and is the half that cannot be fooled.

---

## Log

- 2026-08-27 — Written at the end of the PR-2 session.
