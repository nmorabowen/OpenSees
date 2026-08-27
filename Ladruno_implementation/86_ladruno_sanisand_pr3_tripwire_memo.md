---
title: LadrunoSANISAND — the two fork-tripwire decisions (D5a and §7.3), measured but NOT taken
project: Ladruno
status: decision pending — owner's call
priority: high
owner: nmora
tags:
  - adr
  - decision
  - material
  - sanisand
  - manzari-dafalias
---

# D5a and §7.3 — what they would cost, and what would settle them

> **This memo changes no code and takes no decision.** It was written during ADR-86 PR-3
> under an explicit instruction to investigate and defer. Spec is
> [[86_ladruno_sanisand_adr]]; the architecture note it depends on is
> [[86_ladruno_sanisand_handoff]] §2.
>
> Both items change the **formulation** rather than fix an error, which is where handoff §2
> says the subclass architecture stops being viable. Neither is reachable from a subclass:
> `GetStateDependent` and all three `GetElasticModuli` overloads are **non-virtual**
> (`ManzariDafalias.h:335-342`) and — per D4 and the `ManzariDafaliasRO` shadow hazard — must
> not be made virtual. So each is a vanilla edit or a flag seam, never an override.

---

## 1. The headline: the two low-`p` devices are QUANTITATIVELY COUPLED, and our default changes which one dominates

This was the ADR's suspicion (§7.2.1: *"two overlapping patches for the same problem"*). It is
now a measurement.

`GetStateDependent` computes **one** `p` and hands it to **both** the state parameter and the
dilatancy sigmoid (`ManzariDafalias.cpp:4914`):

```cpp
double p = one3 * GetTrace(stress) + m_Presidual;   // <-- includes p_residual
p = (p < small) ? small : p;
psi = GetPSI(e, p);                                 // ... and then, at :4957
if (p < 0.05 * m_P_atm)
    D_factor = 1.0 / (1.0 + exp(7.6349 - 7.2713 * 101.0 / m_P_atm * p));
```

So `p_residual` **shifts the sigmoid's argument to the right by 1.01 kPa**, and the sigmoid's
own half-suppression point is

```
p_half = 7.6349 / 7.2713 = 1.0500 kPa      against      p_residual = 1.0100 kPa
```

— the same number to within 4 %. Whether that is design or coincidence is a provenance
question (§4 below), but the consequence is arithmetic:

| | `D_factor` floor as true `p` → 0 | at true `p` = 1 kPa |
|---|---|---|
| vanilla, `p_r` = 1.01 (argument `p`+1.01) | `1/(1+e^(7.6349−7.3440))` = **0.4278** | 0.9991 |
| Ladruno default, `p_r` = 0 (argument `p`) | `1/(1+e^7.6349)` = **4.830e-4** | 0.4101 |
| ratio | **886×** | 2.44× |

**`p_residual` was bounding `D_factor` from below at 0.43.** It is not possible, in vanilla,
for this sigmoid to suppress dilatancy by more than a factor of ~2.3 — no matter how low the
true confinement. Setting `p_residual = 0`, which is exactly what `LadrunoSANISAND` does by
default, removes that bound and lets the same unmodified sigmoid suppress dilatancy by up to
~2070×.

Measured on the battery's confine-first deck (10 confinement + 40 isochoric deviatoric steps,
Gorini's constants, `-Pmin` pinned to vanilla's, IntScheme 1), reading committed `p` per step
and evaluating the shipped sigmoid on the `p` the code actually passes it:

| leg | steps with `p` < 5.05 kPa (the trigger) | min `D_factor` on the path |
|---|---|---|
| `p_r` = 1.01 (vanilla) | 36 / 50 | **0.7227** |
| `p_r` = 0 (Ladruno default) | 34 / 50 | **0.0016821** |

Same strain path, same 40 steps, **430× more dilatancy suppression** — from a constant we did
not touch, unmasked by one we did.

### Why this matters for the decision

D5a was filed as a general, family-wide modelling question, deferrable because it affects
everyone equally. It does not affect everyone equally. **On a `p_r = 0` deck the sigmoid is
the dominant remaining low-`p` device, and `p_r = 0` is our default.** The ADR's own
correction box (§5 test 5) established that `D_factor` does not perturb the `η = M^b(ψ)`
identity — it shapes the *path*, not the endpoint — and independently measured the two legs
reaching `p_end` 3.29 vs 5.46 kPa. That path difference is this factor.

This does **not** say the sigmoid is wrong. It says the sigmoid is now load-bearing for our
class in a way it was not for vanilla, and that "defer, it is family-wide" is no longer a
complete description of the risk.

---

## 2. D5a — what would settle it

The code's only stated rationale is the comment above the block:
`// Apply a factor to D so it doesn't go very big when p is small`. That rationale is
**not obviously satisfied by the code**, and this is the first thing to resolve:

- `D = A·(d:n)` with `A = A0(1 + ⟨z:n⟩)` bounded by `z_max`, and
  `α_θ^d = g·M_c·exp(n_d ψ) − m` bounded because `e_c(p) → e_0` as `p → 0`. So **`D` itself is
  bounded as `p → 0`** — it does not "go very big".
- What *is* singular at low `p` is the plastic modulus: `K_p = (2/3)·p·h·(b:n)` with
  `h = b_0/⟨…⟩` and `b_0 ∝ p^(−1/2)`, so `K_p ~ p^(1/2) → 0`.

So the sigmoid reads like a **stability regulariser on the plastic multiplier**, mislabelled as
a bound on `D`. If that is right, `p_residual` and `D_factor` are two patches on the same
singularity applied at different times — and removing one may make the other necessary rather
than redundant. **This is a hypothesis, not a finding.**

### The experiments, in the order that makes each one cheap

Each is an A/B on a **throwaway build**, never a shipped change.

1. **Ablation.** Build with `D_factor ≡ 1` at all four live sites. Re-run, at `p_r = 0` and at
   `p_r = 1.01`:
   - does every step converge, and at what step count (the ADR already records that the
     `p_r = 0` low-`p` leg fails at 400 steps and needs 1200 — that is the regulariser
     question in miniature);
   - the `η/M^b` identity error at the end (expected: unchanged — `D_factor` does not touch
     the endpoint);
   - `p_end`, and the number of ModifiedEuler `dT == dT_min` force-accepts (the
     `mUseElasticTan = true` branch at `ManzariDafalias.cpp:1626`), which is the direct
     measure of "did the integrator struggle".
   **Read:** if convergence and the identity are unchanged and only `p_end` moves, the sigmoid
   is a path-shaping modelling choice and D5a is a calibration question. If the ablated build
   stops converging, it is a numerical regulariser and removing it is not on the table.

2. **The 2×2 factorial.** `p_r ∈ {0, 1.01}` × `D_factor ∈ {on, off}`, one deck, four legs.
   This is the measurement that actually answers *"are they overlapping patches?"* — if the
   sigmoid's effect size is large at `p_r = 0` and negligible at `p_r = 1.01`, they overlap and
   §1's arithmetic is the whole story. If both matter independently, they do not.

3. **Where the band actually bites.** §1 measured 34-36 of 50 steps inside `p < 0.05·P_atm` on
   a deck whose `p_end` is ~10 kPa. Repeat at `p₀` = 100 kPa, the regime real decks run in,
   and report the fraction. If it is zero, D5a is a low-confinement-only question and can be
   scoped to a flag rather than a formulation change.

4. **The shape, not just the switch.** Even granting the sigmoid, is `p/P_atm` the right
   non-dimensionalisation? PR-2 made it `p/P_atm` because that is what the bare constants
   implied at `P_atm = 101`, which is a *units* repair and deliberately not a *shape* claim.
   A state-based reference (e.g. `p` against the current `p_c` or against `m_Pmin`) is a
   different model. Not settleable by A/B on one deck.

### What is NOT settleable by measurement

Where `7.6349` and `7.2713` came from. There is no derivation in the code, no citation in
Dafalias & Manzari (2004), and both sibling implementations carry the identical bare literals —
`SAniSandMS.cpp:2627` and `NTUASand02`. That is a **provenance question for the model's
authors**, and per D8 the route is: inform Prof. Gorini, and let the call be his. Do not
reverse-engineer intent from constants.

### The cost if we act

`D_factor` appears at **five** places in `ManzariDafalias.cpp`. Four are live and were
non-dimensionalised by PR-2 (`:3147`, `:3976`, `:4455` in the analytical Jacobian, each with its
own `temp1` derivative, and `:4958` in `GetStateDependent`). Changing the *shape* means changing
all four **plus the three analytical derivatives**, which is where a shape change stops being a
one-line edit: get a derivative wrong and the Newton iteration degrades silently.

> **A fifth site, found during PR-3 and NOT part of PR-2's "all four".**
> `ManzariDafalias.cpp:3556-3563`, inside `NewtonSol2`, carries a *different* sigmoid:
> trigger `p < 0.001·m_P_atm`, and
> `be = 207232.6584 * 2.0 * m_Pmin; temp1 = exp(20.72326584 - be*p)`.
> Two things are wrong with it and one thing saves it.
> * `be·p` carries **stress²**, so it is dimensionally worse than the four PR-2 repaired.
> * Its steepness is **proportional to `m_Pmin`**. Working backwards, `207232.6584 = 20.72326584/1e-4`
>   and `20.72326584 = −ln(1e-9)`, so `be` reduces to `2·20.723·P_atm` **only when `m_Pmin` has
>   its vanilla value `1e-4·P_atm`**. The constant silently assumes vanilla's floor;
>   `LadrunoSANISAND`'s `-Pmin` default of `1e-3·P_atm` would make it 10× steeper.
> * **It is unreachable.** `NewtonSol2` is called only from `NewtonIter3`, and `NewtonIter3` has
>   no caller anywhere in `SRC/` (verified by grep across the tree; `SAniSandMS` carries its own
>   unrelated pair). So PR-2's claim is **correct about behaviour and incomplete about source**.
> Recorded, not fixed. If `NewtonIter3` is ever wired up, this becomes live and it is a
> `-Pmin`-coupled shape defect, not a units one.

---

## 3. §7.3 — `m_e_init` in the elastic modulus: what it costs the calibration

Dafalias & Manzari (2004) p. 623 give

```
G = G0 · p_at · (2.97 − e)² / (1 + e) · (p/p_at)^½        with e the CURRENT void ratio
```

All three `GetElasticModuli` overloads use `m_e_init` (the INITIAL void ratio) at six lines,
while being handed the current void ratio as parameter `en` and not using it; the correct line
sits commented out three lines above the first; `b0` in the same file uses the current `e`; and
both siblings (`SAniSandMS:2295`, `NTUASand02`) use `en`. **The seam is already open** —
`mUseCurrentVoidRatioInG`, opened by PR-2, `false` in every `ManzariDafalias` constructor,
read as `const double& eG = mUseCurrentVoidRatioInG ? en : m_e_init` at
`ManzariDafalias.cpp:4727`, `:4747`, `:4767`. PR-3 deliberately does **not** write it (see the
note in `LadrunoSANISAND::applyLadrunoConstants`).

### The cost, in one number

Write `f(e) = (2.97 − e)²/(1 + e)`. Then

```
d ln f / de  =  −2/(2.97 − e)  −  1/(1 + e)  =  −1.4691   at e = 0.6944
e = e_init − (1 + e_init)·tr(ε)                       so   de = −1.6944 · d tr(ε)
=>  dG/G  =  −1.4691 · de  =  +2.489 · d tr(ε)
```

So **G moves by about 2.5 % per 1 % of volumetric strain**, in the direction that compaction
softens the material.

| volumetric strain `tr(ε)` | `e` | `f(e)` | `G` vs the frozen form |
|---|---|---|---|
| 0 (the frozen form) | 0.694400 | 3.05616 | — |
| −0.5 % (compaction) | 0.702872 | 3.01835 | **−1.237 %** |
| −1 % | 0.711344 | 2.98101 | **−2.459 %** |
| +1 % (dilation) | 0.677456 | 3.13317 | **+2.520 %** |

*(Recomputed 2026-08-27 after an adversarial review found two rows of an earlier draft
drifting from the formula directly above them — −1.22 / −2.46 / +2.47, against the −1.237
/ −2.459 / +2.520 the expression actually gives. The linear coefficient and the
derivative were correct in that draft; only the worked rows were off. Note the third row
is **not** the mirror of the second: `f` is nonlinear in `e`, which is the whole reason a
single `G0` rescale cannot absorb this change.)*

### Why this cannot be absorbed into `G0`

`G0 = 264.32` was fitted against the frozen form. The obvious repair — rescale `G0` to
compensate — **works at exactly one void ratio and nowhere else**, because the correction is
state-dependent and the state moves during the analysis. Concretely:

- On an **isochoric** path `tr(ε) = 0` throughout, so `e ≡ e_init` and the two forms are
  identical. The battery's confine-first deck is isochoric by construction and measures
  `e` = 0.694384750 on all three legs to nine digits — so **the current battery cannot see
  this change at all**, and a green battery would not be evidence of anything.
- On a **drained** path reaching a few percent volumetric strain, no single `G0` reproduces the
  old `G` over the history. Correcting `G` without refitting `G0` against the same experimental
  data preserves one error by means of another — which is exactly why D6 declined it and why
  D9 routes calibration-moving errors to a seam rather than to a bugfix.

### What would settle it

1. **Is Gorini's `G0` refittable?** This is not a code question. Wiring the seam is a
   one-line change (`applyLadrunoConstants` sets `mUseCurrentVoidRatioInG`); the work is
   re-running the calibration that produced 264.32 against the corrected form. **Whoever owns
   that fit has to be asked before the seam is wired**, or every archived result under
   `G0 = 264.32` becomes ambiguous.
2. **A non-isochoric gate must exist first.** The battery cannot currently detect the change
   (see above). A drained-triaxial deck with a few percent volumetric strain is a prerequisite
   for the seam, not a follow-up to it.
3. **Effect size on a real deck.** Run the A/B on that drained deck with the seam on vs off and
   report the change in `p_end`, `η/M^b` and the initial stiffness. The table above bounds it at
   a few percent in `G`; what that does to a stress path is not derivable from it.

### The asymmetry worth noting

Unlike D5a, §7.3 is **not** a modelling opinion — the paper is unambiguous and both siblings
follow it. It is a genuine error whose only obstacle is a calibration. That makes it the
*easier* of the two to justify and the *harder* to schedule: the code change is one line and
the prerequisite is somebody else's fit.

---

## 4. Recommendation (advisory — the decision is the owner's)

1. **Neither in PR-3.** Done; PR-3 shipped the three low-risk items and this memo.
2. **Raise §1 with Prof. Gorini, and put the 886× number in front of him.** It is the single
   most decision-relevant fact either item produced, it concerns a class we ship with
   `p_r = 0` by default, and under D8 the modelling call is his. This is the one item that
   should not wait for a spare PR.
3. **D5a: run experiment 1 (ablation) and experiment 2 (the 2×2) before deciding anything.**
   They are throwaway builds and a few hours. Deciding the shape question without them is
   guessing; the ADR has already had one such prediction refuted by measurement (§5 test 5).
4. **§7.3: build the non-isochoric gate first, then ask about the fit.** Wiring the seam
   before there is a test that can see it would be a change nothing could verify.
5. **Interim disclosure, which PR-3 ships:** `Print` now states that `D5b` is done and `D5a`
   open, and names the 1.050 kPa half-suppression point against vanilla's 1.01 kPa
   `p_residual`, so any archived run says on its face which regime it was in.

---

## Log

- 2026-08-27 — Written during ADR-86 PR-3, under instruction to investigate both tripwire items
  and change nothing. All numbers measured on an engine rebuilt from the PR-3 worktree
  (`dist/bin/opensees.pyd`, hash `647a60ad…`), except the §3 table, which is arithmetic on the
  published `G` expression and is reproducible with a calculator.
