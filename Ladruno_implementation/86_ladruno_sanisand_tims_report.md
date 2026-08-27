---
title: The Hidden Cohesion — evidence report for the TIMs project team
project: Ladruno
status: issued
date: 2026-08-27
audience: TIMs project team
tags:
  - report
  - evidence
  - material
  - sanisand
  - manzari-dafalias
---

# The Hidden Cohesion

**Evidence report for the TIMs project team — 27 August 2026**

OpenSees' `ManzariDafalias` adds a hardcoded residual pressure to every mean-stress evaluation.
In a cohesionless sand it behaves as roughly **0.95 kPa of apparent cohesion** — negligible at
100 kPa, dominant below 10, and invisible because every run converges and reports success.

| | |
|---|---|
| Fork | Ladruno / OpenSees, branch `ladruno` |
| Shipped in | PR #767 (`4870f802c`), PR #768 (`9922dfb1f`) |
| Specification | `Ladruno_implementation/86_ladruno_sanisand_adr.md` (ADR 86) |
| Status | Merged |

---

## 1. What this means for your results

> **Above about 50 kPa of confinement, this does not change your numbers meaningfully.**
> The effect scales as `p_r / p'`, so it is under 2 % at 50–200 kPa.
>
> **Below about 10 kPa, your results carry an apparent cohesion nobody declared**, rising to
> roughly **+20 % at 1 kPa**. And if you compared an `OpenSeesMP` run against a serial one, or
> restored a model from a database, **those two runs were not the same material.**

None of this announced itself. Every leg converged. This is a silent-wrong-answer defect, not a
robustness one — which is exactly why it survived this long.

---

## 2. The constant, and what it does

`ManzariDafalias::initialize()` assigns, unconditionally and with no way to override it from a deck:

```cpp
m_Pmin      = 1.0e-4 * m_P_atm;   //  0.0101 kPa at P_atm = 101
m_Presidual = 1.0e-2 * m_P_atm;   //  1.01   kPa at P_atm = 101
```

`m_Presidual` is described in the header as a *"small residual pressure (due to cohesion)"* and is
threaded through **thirty** mean-stress evaluations as
`p = one3*GetTrace(stress) + m_Presidual`. It was not a constructor argument, not settable from the
deck, not echoed anywhere, and — see §3 — not on the wire.

Crucially it reaches `GetPSI()` as well. So the residual pressure perturbs **ψ**, the state
parameter the whole critical-state formulation is organised around, and through ψ it displaces
`M^b`, `M^d` and the dilatancy together. Dafalias & Manzari (2004) put no residual pressure there.

### Measured against the model's own identity

The cleanest instrument needs no experiment at all: at the end of a drained shear the formulation
satisfies `η = M^b(ψ)` by construction. Any departure is the instrument reading itself.

| p₀ (kPa) | 100 | 10 | 3 | **1** | 0.5 | 0.2 |
|---|---|---|---|---|---|---|
| departure from its own bounding surface | +1.0 % | +3.7 % | ~+8 % | **+20.1 %** | ~+30 % | **~+52 %** |

The invariant `err × p'_end ≈ p_r` reproduces across the whole sweep, which identifies the constant
without needing any modelling context: at p₀ = 1 kPa, `err × p'_end = 1.0136` against
`p_r = 1.0100`.

**Why "+52 %" understates it.** η saturates at a geometric ceiling in a triaxial test. In physical
terms the model returns mobilised friction angles up to **83°** where its own `M^b` says 47°, and
the deviator barely changes while the confinement that should produce it falls fiftyfold.
**A strength independent of confinement is a cohesion.**

---

## 3. Findings

### 3.1 Serial and parallel were running different soils — *fixed in the new class*

`sendSelf` serialises into `Vector data(97)`; `m_Pmin` occupies the last slot.
**`m_Presidual` has no slot at all.** On the broker path the null constructor sets `m_P_atm = 0` and
then calls `initialize()`, so `m_Presidual = 1.0e-2 × 0.0 = 0`; `recvSelf` restores `m_P_atm` but
never re-runs `initialize()`.

So an `OpenSeesMP` worker, or any database-restored material, ran at `p_r = 0` while the serial
process beside it ran at 1.01 kPa. It is also path-dependent: `revertToStart()` calls `initialize()`
mid-life, so the constant reappears.

> **Measured** — a restored vanilla material changes its stress **5.5 %** the instant it is
> restored, and 4.2 % by the end of the path. Any MP-vs-serial comparison on this material at low
> confinement was invalid.

### 3.2 A substep interpolant anchored to the wrong end — *fixed in vanilla*

Inside the substepping loops the pseudo-time `T` sweeps 0→1 across the strain increment, so the void
ratio at `T` must be evaluated at `CurStrain + T·dStrain`. All three substepping integrators
evaluated at `NextStrain + T·dStrain` — a full increment ahead at `T = 0`. Void ratio feeds ψ, which
feeds `M^b`.

**Independently corroborated.** `NTUASand02` (Gorini) contains the same bug, found and fixed on
`21/09/2021`, dated in the source with the buggy line left commented directly above the fix — and
written character-for-character as we wrote it. `PM4Sand` and `PM4Silt`, same folder and lineage,
always had it right.

> **Blast radius, measured** — `IntScheme 1` (ModifiedEuler, the default) moves **2.8e-6** relative.
> `IntScheme 3` does not move. **`IntScheme 45` moves 6.33e-3** — about 0.6 %. If you ran RK45
> decks, your results shift by that order.

### 3.3 A dilatancy sigmoid that changes meaning with your unit system — *fixed in vanilla*

Below `p < 0.05·P_atm` the model scales dilatancy by `1/(1 + exp(7.6349 − 7.2713·p))`. The trigger
scales with the unit system; **the constant does not** — `7.2713` multiplies a raw stress and is a
bare literal, so it silently carries units of 1/stress.

> **At 1 kPa true confinement** the factor is **0.410 in kPa**, **1.000 in Pa** (never fires),
> **0.0005 in MPa** (total suppression). Three different materials from one input file, and OpenSees
> has no unit system to catch it. Now non-dimensionalised at all four sites; an exact no-op at
> `P_atm = 101`.

### 3.4 A stress clamp that had never printed for anyone — *now observable*

When mean pressure falls under the floor, the integrators silently *rebuild* the stress tensor. The
notice was gated behind `debugFlag`, a `static const bool = false` with **61 dead guards** in the
file. The clamp changes the answer and said nothing.

> **Now** — throttled to 11 lines per site, independent of mesh size. It does not fire on ordinary
> decks; it fires when `p` genuinely reaches the floor, which is what you want to know about.

### 3.5 The bounding surface gets silently inflated 50 % — *found, not changed*

If a deck's stress ratio is already above `M_c` when it switches to the plastic stage,
`Elastic2Plastic` raises `M_c` to accommodate it — from the calibrated **1.3309** to **1.99878** on
our own test deck, once per Gauss point.

> **What it means for you** — if you apply deviatoric load during the elastic stage and then flip to
> plastic, you may be running a material whose bounding surface is not the one you calibrated.
> **Confine hydrostatically first, then switch, then shear.** That ordering produces zero inflation;
> we rebuilt our own test decks around it.

### 3.6 On a perfectly proportional path with `p_r = 0`, the sand never yields — *characterised, pinned by a test*

Every *elastic* step re-sets the back-stress ratio to `s/(p + p_r)`. On a single proportional strain
ramp the stress path is exactly radial, so with `p_r = 0` that ratio is constant along the ray — the
yield surface is welded to the stress and `f < 0` forever. Zero plastic strain, on every integration
scheme.

> **Where you could hit this** — a single-element calibration deck running staged gravity and then a
> proportional ramp. It is a property of *proportional* decks, not of `p_r = 0`: `NTUASand02` ships
> `p_r = 0` with identical code and is used for staged foundation work, where the stress direction
> changes constantly. Threshold measured between `-Presidual` 0.0230 and 0.0240 kPa.


### 3.7 The elastic moduli never evolve with density — *switch opened, deliberately not wired*

`commitState` does everything right: it commits the strain, recomputes the current void ratio as
`mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon)`, and hands it to `GetElasticModuli`
as the parameter `en`.

**`GetElasticModuli` then ignores `en` and uses `m_e_init`.** So it is not that the void ratio is
initial rather than committed — the committed value exists, is correct, and is passed in. It is
discarded at the point of use. Directly above that line, commented out, sits the correct form,
including a more careful secant branch using both `en` and `en1` for a finite volumetric increment.

Dafalias & Manzari (2004) p.623 make `e` in `(2.97-e)^2/(1+e)` the **current** void ratio.
`SAniSandMS` and `NTUASand02` both use the current one. `b0`, in this same file, uses the current
one.

Because `K = (2/3)(1+nu)/(1-2nu) * G`, **both** elastic moduli are frozen at the initial-state value
for the whole analysis: a densifying sand keeps its loose stiffness, a dilating one keeps its dense
stiffness.

**Magnitude.** At your `e_init = 0.6944` the sensitivity is `dG/G = -1.47 * de`:

| change in void ratio | error in G *and* K |
|---|---|
| 0.002 (our low-`p` verification deck) | −0.3 % |
| 0.02 | −2.9 % |
| 0.05 | −7.3 % |
| 0.10 | **−14.7 %** |

> **Negligible on a monotonic, low-strain path; material wherever the void ratio actually moves** —
> cyclic loading, densification, liquefaction. It compounds with §2 at low confinement, where both
> distort stiffness and strength at once.

**Why it was not simply fixed.** `G0 = 264.32` was fitted **against the frozen form**. Correcting
`G` without revisiting `G0` would preserve one error by means of another — a differently-wrong
stiffness under a calibration that no longer matches it. PR-2 therefore opened a switch,
`mUseCurrentVoidRatioInG`, default off, vanilla bit-identical, with all three `GetElasticModuli`
overloads seamed together so they cannot drift apart — and **deliberately left it unwired**.

> **This one is yours to decide.** Flipping it is a one-line change plus a refit of `G0`. Nobody
> outside the calibration can make that call.

---

## 4. What was implemented

`nDMaterial LadrunoSANISAND` — a thin subclass of `ManzariDafalias`. It contains **no physics**. The
5,229 lines of constitutive code run unchanged; the subclass only makes the two low-stress constants
settable, carries them across the wire, and echoes them at construction.

```tcl
nDMaterial LadrunoSANISAND $tag <18 params> <IntScheme TanType JacoType TolF TolR> \
    <-Presidual $pr> <-Pmin $pmin>
```

The first 18 positional arguments and the five optional ones are **identical to `ManzariDafalias`**,
so a deck migrates by renaming the command. Defaults are `p_residual = 0.0` and
`p_min = 1.0e-3·P_atm` — exactly the values `NTUASand02` chose, where the zero is documented as
deliberate: *a cohesionless sand has no cohesion.*

### Why a new class rather than changing the default

Changing `ManzariDafalias`'s own default would silently alter every archived result, and exposing
the constant there would need a wire-format change that breaks mixed-build parallel runs. A separate
name keeps every past result reproducible under the name it was measured with, and makes *"which
soil did this number describe"* a question about a class name rather than a build hash.

> **The claim this licenses.** Given vanilla's two constants, `LadrunoSANISAND` returns
> **bit-identical** committed stresses to `ManzariDafalias` — relative difference exactly `0.0`.
> That is what permits the sentence *"these two runs differ in exactly one constant and are
> otherwise the same code."* It is the property an A/B experiment needs, and a forked copy could
> never provide it.

### Verification

| Check | Result |
|---|---|
| Equivalence to vanilla, given vanilla's constants | **0.0 — bit-identical** |
| Sensitivity to `p_residual` alone (1.01 → 0.5, both plastic) | 4.24e-02 |
| Departure at p₀ = 1.01 kPa — vanilla constants | **+18.1258 %** |
| Departure at p₀ = 1.01 kPa — `p_r = 0` | **+0.0754 %** |
| Invariant `err × p_end` vs `p_r` | 0.9903 / 1.0100 |
| Database round-trip carries `p_residual` | yes |
| Bounding-surface inflation on the confined deck | 0 events |

Setting `p_r = 0` removes the departure essentially completely — **+0.075 %, zero to measurement**.
That result corrected the specification, which had predicted a residual departure on the grounds
that the dilatancy sigmoid would still suppress ~59 % of dilatancy there. The sigmoid scales how the
state *travels* to the bounding surface, not the identity that holds once it is there.

---

## 5. What to do

1. **Triage by confinement, not by deck.** Above ~50 kPa, no action. Below ~10 kPa, treat past
   results as carrying an undeclared apparent cohesion and re-run. The rule of thumb is the
   invariant: the error in stress ratio is approximately `p_r / p'_end`.

2. **Treat any past MP-vs-serial comparison at low confinement as invalid.** Not approximate —
   different materials. The same applies to anything restored from a `FileDatastore`. Re-run those
   under `LadrunoSANISAND`, where the constant is on the wire.

3. **Migrate by renaming the command, and pin both constants explicitly.** Write `-Presidual` and
   `-Pmin` in the deck even when you want the defaults. The material echoes what it is running at
   construction, so the record states it. Note the new default `-Pmin` is `1e-3·P_atm`, **ten times
   vanilla's** — pin it in any A/B comparison so only one variable moves.

4. **Confine hydrostatically before switching to the plastic stage.** Applying deviatoric load
   during the elastic stage and then flipping can inflate `M_c` by 50 % (§3.5). Confine first,
   switch, then shear.

5. **Expect small shifts in existing vanilla results.** The interpolant repair moves `IntScheme 1`
   decks by 2.8e-6 and **`IntScheme 45` decks by ~0.6 %**. If you bisect a change of that order to
   this work, that is the cause — and it is a correction, not a regression.

---

## 6. Still open

Recorded rather than fixed, so nothing here is hidden behind a green build:

- The frozen elastic moduli of §3.7 — the switch is open and unwired, awaiting a calibration
  decision that is yours to make.
- Whether the dilatancy sigmoid should exist at all once `p_residual = 0` is a modelling question
  for the model's authors. The two look like overlapping patches for the same problem applied at
  different times.
- `SAniSandMS` carries the same interpolant bug, an argument-parsing bug that silently drops
  `TolF`/`TolR`, and a sigmoid trigger of `0.001·P_atm` rather than `0.05` — over that window it
  kills dilatancy by a factor of about 1000. Scheduled for its own change with its own regression
  tests.

**On upstreaming.** These repairs are genuinely upstreamable, and every one we hold is a maintenance
liability for this fork. The decision recorded is that we do not open an upstream pull request
ourselves — **Prof. Gorini is to be informed, and whether it is worth making is his call.** He found
the interpolant bug independently in 2021; because that fix stayed private, it was rediscovered from
scratch five years later.

---

## Provenance

Specification and full audit trail in `Ladruno_implementation/86_ladruno_sanisand_adr.md`, with
corrections recorded in place. Defect provenance in `LEDGER_quirks.md`; every upstream file touched
is listed in `LEDGER_vanilla_files.md`. Engineering handoff for the next session in
`Ladruno_implementation/86_ladruno_sanisand_handoff.md`. Shipped in PR #767 and PR #768.

Figures in §2 marked as a sweep are from the originating audit. All figures in §4 were re-measured
on an engine rebuilt from the merged source.
