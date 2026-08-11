# Note 83 — The tangent-free path probe: does a quadratic element finish the Prandtl collapse when no tangent is formed?

**PR:** [#727](https://github.com/nmorabowen/OpenSees/pull/727)
**Kind:** MEASUREMENT. Ships no feature, changes no element, changes no
integrator.
**Predecessors:** [[81_quadratic_hex_limit_load_measurement]] and
[[82_quadratic_path_loss_diagnosis]]. Read note 82 §0 and §8 first; **their
retractions are load-bearing and nothing they withdrew is re-derived here.**
This note answers note 82 §8 item 3 — the experiment that note recommended be
run before anything else.
**Harness:** `Ladruno_files/testbed/hypo_bearing/quad_dr_probe.py`
(+ `quad_dr_summary.py`, `dr_elastic_check.py`), built on `h20_prandtl.py`
(mesh, consistent Q8 surcharge, material, cone, oracle, surcharge step and all
of its controls — **imported, not copied**).
**Engine:** `ladrunoBuild() = bbb2ad732a2bd1ae714fa6d710196b112c67258d`
**Problem:** Prandtl–Reissner, `q_u = q0·N_q = 138.907 kPa` exact,
φ_txc = 20°, q₀ = 10 kPa, γ = 0, ν = 0.45, Chen & Han plane-strain cone match,
non-associated (ρ̄ = 0). Unchanged from notes 81 and 82 and from the merged R3
gate.

> **ATTRIBUTION AND HOLD.** Where this note leans on results credited to the
> **TIMs element-technology campaign** it does so through notes 81/82, which
> cite them with credit. TIMs material is internal-only and under the same hold
> as the unfiled UW report: do not carry their numbers or phrasing into anything
> public or upstream-facing.

---

## 0. The answer, and the two things that almost made it the wrong answer

**THE BINARY QUESTION.** Note 82 established that the linear relieved hex
reaches the exact Prandtl answer on a genuine plateau while *every* quadratic
leg — hex and tet, Lagrange and Bernstein, relieved and locked — stops early on
a path-controller limit at s/B ≈ 0.001–0.011, and that at the point of failure
the quadratic element's tangent is **flat** and its mechanism is **as developed**
as the element that goes on to succeed. Four explanations are dead. The one
remaining fork: is the wall in the **Newton path**, or in the **discretization /
material return map**? Dynamic relaxation settles it because it removes the
tangent from the question entirely — `LadrunoDynamicRelaxation` (33005) never
forms, factorizes or trusts one; `formTangent` pokes a diagonal fictitious mass
onto the SOE and the solve degenerates to an `M*⁻¹` apply.

**THE ANSWER, in one line: the wall is in the NEWTON PATH.** The linear control
reproduces the established answer under DR — **0.9968 of exact, on a plateau,
`TARGET`, a capacity on all three clauses**, with every one of its ten hold
checks landing within 0.45 % of the static Newton curve (§5.1). And the
quadratic element, on the same mesh under the same tangent-free path, **walks
straight through s/B = 0.01009 — the settlement where the Newton path died —
and reaches settled, gate-passing static equilibria at 1.11× and 1.24× beyond
it, with nothing singular happening at the wall** (§5.2). Two DR configurations
disagreeing in three controller settings agree there to **0.21 %**.

**What this note does NOT establish**, and the distinction is the whole of §5.3:
the quadratic leg **did not finish the problem**. It stopped at s/B = 0.0138 on
the **wall clock**, against the control's 0.0500 and its plateau — a
**CONTROLLER ALLOWANCE**, so no quadratic number here is an element ceiling. The
binary question is answered in its first half only: *the wall is gone*; whether
a quadratic element plateaus at the exact answer is still untested. A flattening
in the last three holds is consistent both with a real ceiling near 0.70 of
exact and with the relaxation degrading, and three points cannot separate them.

Before all of it, the two negatives — because both produced a *fully converged,
entirely wrong* number first, and either one alone would have been published as
a result.

1. **The Gershgorin fictitious mass puts the march exactly ON the
   central-difference stability boundary (§3).** `-mass gershgorin` sets
   `m_i = (dt²/4)·Σ_j|K_ij|`, which bounds `ω_max·dt ≤ 2` — and for an FE
   stiffness that bound is very nearly *attained*, so the march amplifies
   round-off instead of damping it. Measured: at the nominal step the zero-push
   hold departs from an exact equilibrium and reaches a residual of **87 kN on a
   300 kN problem**. A safety factor is mandatory, it is **not** exposed as an
   option, and its value is element-dependent.
2. **A displacement increment that is STEPPED rather than RAMPED silently
   rewrites the load history (§4).** Stepping the increment and relaxing gives a
   genuine static equilibrium — the residual really does reach 1e-13 kN — but it
   is the equilibrium of a *different path*: the discontinuity launches a
   transient through a mesh whose mass is **fictitious**, so it has no physical
   wave speed, and on an elastic-plastic material the overshoot leaves plastic
   strain behind. Measured: **−24 % to −29 % against the static Newton answer at
   matched settlement, with a perfectly converged residual at every point.**

Both were caught only because the harness runs a control that compares DR
against a state we already know (control DR-0 in §2, and §4.1). Neither would have been caught by
looking at residuals, and the second would have produced a self-consistent
quadratic-vs-linear comparison in which *both* legs were wrong by ~25 %.

---

## 1. Why DR is the decisive probe — and a correction to the acceptance rule

DR reaches static equilibrium by integrating a fictitious damped transient to
rest: `M*ü + C*u̇ + f_int(u) = f_ext`, with `M*` and the pseudo-step chosen for
fastest decay, not physical accuracy. When the transient stops, what survives is
`f_int(u) = f_ext`. No tangent is involved anywhere. So:

- if a quadratic element **finishes** under DR, the wall was in the Newton path,
  not in the discretization, and quadratic elements are usable for limit loads
  via a tangent-free path;
- if it **walls** under DR too, the entire Newton/conditioning family of
  explanations is eliminated and the Drucker–Prager apex / return-map probe
  (note 82 §8 item 1) becomes the next move.

### 1.1 "It never settles at collapse" is the WRONG acceptance rule here

The natural DR acceptance is: below collapse the structure settles at each
imposed increment; at collapse it never settles, it keeps accelerating under the
unbalanced load. **That rule is correct under LOAD control and inapplicable
here**, and the reason is structural rather than a detail:

> this problem is **displacement-controlled** — a rigid footing, `sp()` on dof 3
> of every footing node. The collapse mechanism carries footing settlement, so
> under displacement control the mechanism is **driven, not free**: it is not in
> the span of the free DOFs at all. The constrained operator therefore stays
> non-singular through collapse — which is exactly what note 82 §7.3 *measured*
> for the linear element, σ_min flat at 1.2e-3 all the way to a fully formed
> mechanism on a dead-flat plateau. Every increment settles, at and past
> collapse.

Changing the loading to load control to recover the slogan would have changed
the problem away from the established legs, and comparability with notes 81/82
is worth more than a tidier acceptance test. So the acceptance is built on what
displacement control actually offers, and it is stated before any leg was run:

- each increment is **held** and relaxed to static rest, so every recorded
  `(s, q)` point is a genuine **static equilibrium** point, directly comparable
  with the static campaign's converged points;
- **"settled" = the true static unbalance** `‖f_ext − f_int‖_∞ < 1e-6 × 300 kN
  = 3.0e-4 kN`, i.e. ten times tighter in absolute force than the static
  ladder's own `NormUnbalance` acceptance;
- the **collapse signature is still the q(s) PLATEAU**, unchanged from the
  static campaign, so the two are judged on the same footing;
- an increment that runs out of DR iterations is the DR analogue of the static
  controller's **step floor**: an ALLOWANCE, never a capacity.

**A design fault caught in passing, and worth recording.** The first version of
the settle gate also required "the footing reaction stopped moving over a
chunk". That gate is **chunk-length dependent** — shorten the chunk and less
happens inside it, so the same state passes. A gate whose verdict moves with a
reporting parameter is measuring the reporting parameter. The residual is
chunk-free and is now the whole gate; `dR` is written to the CSV and decides
nothing.

### 1.2 Capacity rule

The merged R3 gate's three clauses, transposed:

1. **plateau** — tail dq/ds < 2 % of the initial tangent;
2. **termination mode `TARGET`**. `ITERCAP`, `DIVERGED`, `WALL` and `TRUNCATED`
   are seizure modes and are never a capacity;
3. **free advance** — over the final tenth of the run no increment spent more
   than 50 % of its iteration allowance.

Any leg failing a clause is a **CONTROLLER ALLOWANCE and the allowance is
named**.

## 2. Controls

Controls 0–4 are `h20_prandtl.py`'s, imported and unchanged, so the expensive
load algebra cannot fork. Control DR-0 is new and is the one that carries this
note.

| control | what it tests | h8bbar h0 = 0.5 | h20uri h0 = 0.5 |
|---|---|---|---|
| **0** provenance | `ladrunoBuild()` recorded before any number | `bbb2ad73…` | `bbb2ad73…` |
| **2** consistent surcharge | Q4 / Q8 weights vs analytic, sum = top area | 26 faces, 30.000000000 m², err 3.1e-17·A | 26 faces, **54 nodes carry a NEGATIVE weight**, err 1.2e-16·A |
| **3** resultant identity | Σ R_z = q₀·Σw | 300.000000 vs 300.000000 (+0.0000000 %) | 300.000000 vs 300.000000 (+0.0000000 %) |
| **1** 1-D elastic stress patch | the **distribution**, which the identity is structurally blind to | σ_zz 4.5e-14, σ_xx 3.1e-14, τ 1.4e-14 | σ_zz 2.4e-13, σ_xx 2.9e-13, τ 1.5e-13 |
| **DR-0** zero-push hold | DR reproduces the static state it was handed | **R = 12.500000 vs 12.500000, +0.000000 %**, \|du\| = 0, KE 6.0e-31, res 2.9e-13 kN | **R = 10.833333 vs 10.833333, +0.000000 %**, \|du\| = 0, KE 6.6e-30, res 2.1e-12 kN |

**Control DR-0, and why it is run as a long hold rather than a quick look.**
After the static surcharge step, the footing's *current* displacement is
re-imposed as an SP and DR is run. The state is already a static equilibrium, so
DR must not move it. Run for a few hundred steps this passes trivially. Run for
**3000–6000** steps it catches §3 — which is the entire reason the hold length
is a named parameter (`--dr0-steps`, default 3000) rather than "until it looks
quiet".

*What would have to be true for DR-0 to pass while DR is broken?* DR-0 tests a
state that is already in equilibrium, so it cannot detect an error that only
appears once the material yields — and §4 is exactly such an error. DR-0 is
necessary and **not** sufficient, and §4.1 is the control that covers the gap.

## 3. NEGATIVE ONE — the Gershgorin mass sits on the stability boundary

`buildGershgorinDiagonal` (`SRC/analysis/integrator/LadrunoFictitiousMass.h`)
builds `m_i = (dt²/4)·Σ_j|K^e_ij|` accumulated over elements. By Gershgorin
`λ_max(M*⁻¹K) ≤ 4/dt²`, hence `ω_max·dt ≤ 2` — the central-difference stability
limit, **with equality**, not with margin. At the limit the amplification matrix
is defective and round-off is amplified rather than damped.

The integrator exposes **no safety factor**. It can be applied from the script
side, because the mass is sized for the integrator's `-dt` while the march uses
the `analyze` dt: passing `-dt D` and marching at `dt = f·D` gives
`ω_max·dt ≤ 2f` and scales the per-step relaxation by `f²`. `--dtfac` is that
`f`.

**Measured** — zero-push hold, `h8bbar` h0 = 1.0, from an exact static
equilibrium (the only thing varied is `--dtfac`):

| dtfac | steps held | KE at end | res at end (kN) | verdict |
|---|---|---|---|---|
| 1.00 | 6000 | 7.5e+01 | **8.7e+01** | **UNSTABLE** |
| 0.90 | 6000 | 2.7e-02 | **6.9e+01** | **UNSTABLE** |
| 0.85 | 3000 | 8.8e-02 | **1.4e+01** | **UNSTABLE** |
| 0.75 | 3000 | 1.2e-31 | 8.3e-14 | stable |
| 0.70 | 3000 | 1.1e-31 | 2.0e-13 | stable |
| 0.60 | 3000 | 1.7e-31 | 2.5e-14 | stable |
| 0.50 | 6000 | 3.1e-31 | 4.1e-14 | stable |

The stable rows are at **round-off** — res 1e-13 kN on a 300 kN problem, and the
footing reaction printing as `12.500000` for six thousand steps. The unstable
rows grow from that same round-off: at `dtfac = 1.0` the residual rose from
2.6e-13 kN to 8.7e+01 kN over 5000 steps, an amplification of ~1.006 per step,
and the reaction drifted from 15.000000 to 10–13 kN as the spurious oscillation
drove Gauss points past yield.

**The threshold is element-dependent.** `h20uri` h0 = 0.5 fails the same hold at
`dtfac = 0.75` (KE 2.9, res 8.0e+01 kN after 1000 steps) and passes at 0.50
(KE 6.6e-30, res 2.1e-12 kN over 3000 steps). So the linear hex tolerates
~0.8 and the quadratic hex ~0.5–0.75. That is a legitimate per-leg setting — an
explicit core has an element-dependent stable step — and both legs below are run
at **dtfac = 0.5**, inside both thresholds, with the value recorded.

**This is a defect in `LadrunoDynamicRelaxation`, not in this measurement**, and
it is written up as such in `LEDGER_quirks.md`. Its practical reach is wider
than this note: `robust_drive.py` rung 5 and `torture_dynamics.py` both call
`integrator("LadrunoDynamicRelaxation")` with the default `-dt`, i.e. at
`dtfac = 1`. Short excursions will not show it; long ones will.

## 4. NEGATIVE TWO — a STEPPED displacement increment rewrites the load history

Holding the imposed settlement and relaxing produces an exact static
equilibrium. It does **not** produce the equilibrium of the *quasi-static path*.
A displacement discontinuity at the footing launches a transient through a mesh
whose mass is **fictitious** — it has no physical wave speed, and the near-field
elements absorb the whole jump in a thin layer. On an elastic material that
rings and relaxes away harmlessly. On an **elastic-plastic** material the
overshoot leaves plastic strain behind, and the settled state is permanently
softer.

**Measured, against the static Newton leg on the identical mesh at matched
settlement** (`h8bbar` h0 = 0.5, ds = 1 mm stepped, `dtfac` 0.5): the DR curve
runs **−28.9 %** at s/B = 0.0005, decaying monotonically to **−8.1 %** at
s/B = 0.0455 — while reporting a settled residual of 1e-6 to 1e-13 kN at every
single point.

### 4.1 The control that separated "wrong plumbing" from "wrong path"

A factor of ~2 in the elastic range cannot be plasticity, so the mechanics were
isolated from the material: the same mesh, the same consistent surcharge and the
same imposed 1 mm footing displacement, driven once by static Newton and once by
DR, on an **ElasticIsotropic** material where the answer is unique and path-free
(`dr_elastic_check.py`). It also prints the footing displacement actually
reached, because "the imposed displacement is only partly applied" and "the
structure is half as stiff" produce the same reaction and are otherwise
indistinguishable.

| elastic, h0 = 0.5, ds = 1 mm | static Newton | DR | DR/static | \|u − u_target\| |
|---|---|---|---|---|
| `h8bbar` R (kN) | 12.500000 → **28.620369** | 12.500000 → **28.620369** | **1.000000** | 0.00e+00 |
| `h20uri` R (kN) | 10.833333 → **26.708338** | 10.833333 → **26.708446** | **1.000007** | 0.00e+00 |

**The plumbing is exact.** The SP is applied in full, the reaction is measured
identically, and the two paths find the same elastic equilibrium to seven
figures. So the softening is real spurious plastic straining — a quasi-staticness
violation — and not a harness error.

### 4.2 The remedy, and its convergence

Each held increment is applied as `--ramp N` equal sub-jumps separated by
`--rampevery K` DR steps, then relaxed. The control variable is the resulting
**displacement rate per DR step**.

**Measured** (`h8bbar` h0 = 0.5, `dtfac` 0.5, each leg compared against the
static Newton curve at its own end settlement):

| ramp | rate (m/step) | s/B at end | DR q | static q | **diff** | DR steps/inc |
|---|---|---|---|---|---|---|
| **1 (stepped)** | — (discontinuity) | 0.0050 | 61.876 | 81.947 | **−24.49 %** | 1 500 |
| 10 × 40, ds 2e-4 | 5.0e-7 | 0.0012 | 46.414 | 47.453 | −2.19 % | 1 860 |
| 50 × 40 | 5.0e-7 | 0.0030 | 65.839 | 66.989 | −1.72 % | 3 460 |
| 200 × 40 | 1.25e-7 | 0.0010 | 43.156 | 43.508 | −0.81 % | 9 210 |
| **50 × 100** | **2.0e-7** | 0.0015 | 51.674 | 51.985 | **−0.60 %** | 6 067 |

Two things this says, and one it does not.

- The error is controlled by the **rate**, not by the increment size: `ds = 2e-4`
  with 10 sub-jumps and `ds = 1e-3` with 50 sub-jumps have the same rate and land
  within 0.5 % of each other (−2.19 % vs −1.72 %) despite a 5× difference in
  `ds`. Shrinking `ds` alone is the expensive way to buy quasi-staticness.
- It converges toward the Newton path from **below**, monotonically in rate, but
  **it does not reach zero**, and the sub-jumps are themselves mini-shocks.
- The **quadratic element is markedly more rate-sensitive than the linear one**:
  at rate 5.0e-7 the linear leg tracks static to −0.3 %/−3.7 % over its first two
  increments while the quadratic reads −5.9 %/−11.2 %/−7.7 %. Any DR comparison
  across element orders at a *fixed* rate is therefore comparing two different
  amounts of contamination, which is exactly the kind of shared-mechanism
  agreement notes 81/82 kept getting caught by.

### 4.3 The remedy that was adopted: a genuinely CONTINUOUS ramp, bought with the dt null-control

The sub-jump ramp still steps, just more finely, and it pays a `domainChanged`
per sub-jump — which is where most of its wall clock goes. The **dt null-control
of §6.1 makes a better construction available**: with `-mass gershgorin` the
march depends only on the *ratio* of the analyze `dt` to the integrator's `-dt`,
never on their absolute size. The pseudo-time step is therefore **free to be
repurposed as the loading clock**. Set `-dt = rate/dtfac`, march at `dt = rate`,
and give the push pattern a `Linear` series carrying `sp(node, 3, −1.0)`:

    uz = −1.0 · t,   t₀ = −uz₀   ⇒   uz = uz₀ − n·rate  after n DR steps.

The dynamics are bit-identical to the same run at `dt = 1.0`; what is bought is
that the imposed displacement is now **continuous in the DR step**, with no
sub-jumps and no `domainChanged` anywhere in the march.

A continuously ramped point is quasi-static, not static — it carries whatever
inertial lag the rate leaves behind — so the lag is **measured rather than
assumed away**: at intervals the ramp is frozen, a `Constant`-series pattern is
installed at the current displacement, and the state is relaxed to the same
residual gate every held leg uses (`hold_check`). The headline always comes from
a **final hold**, i.e. from a genuine static equilibrium point.

**Measured** (h0 = 0.5, `dtfac` 0.5, against the static Newton curve; the
continuous ramp errs on the *high* side because the lag is a lag):

| element | rate (m/step) | tracking vs static, first → last sample |
|---|---|---|
| `h8bbar` | 2.0e-6 | −2.63 → +0.05 → +0.21 → −1.77 % |
| `h8bbar` | **5.0e-7** | +3.06 → +1.55 → +0.98 → +0.54 → +0.01 → −0.04 → +0.33 → **−0.42 %** |
| `h20uri` | 2.0e-6 | +16.55 → +6.35 → +6.85 → +8.63 % |
| `h20uri` | **5.0e-7** | +9.00 → +4.52 → +1.51 → +2.20 → −0.87 → +0.43 → +1.15 → **+0.54 %** |

At **rate 5.0e-7 both elements settle onto the static Newton curve to within
about ±1 %** after the start-up transient — a fivefold better result than the
sub-jump ramp achieved at the same rate, and cheaper. That is the operating
point for §5, and `--holdevery` reports the inertial lag along the way.

*What would have to be true for these checks to pass while the result is wrong?*
They compare against the static Newton leg, so they inherit that leg's
correctness — but that leg is the one the merged R3 gate certifies against an
exact oracle, the strongest anchor available. They are also blind to any error
that affects DR and Newton identically; §2 control 1 and the exact oracle cover
that. And the tracking check only constrains DR **where the static leg exists** —
for the quadratic element that is s/B ≤ 0.0101, i.e. exactly up to the wall and
no further. Beyond it, the quadratic DR leg is unanchored, and §5 says so.

## 5. THE LEG TABLE — the deliverable

All legs: `h0 = 0.5`, non-associated, continuous ramp, `-mass gershgorin`,
`-damping kinetic`. `q_exact = 138.907 kPa`. **Headline q is taken from a HOLD**
— a relaxed static equilibrium — never from a ramped sample.

| leg | DOF | mode | rate (m/step) | dtfac | end s/B | **q_hold** | **q/q_exact** | tail % | plateau | **CAPACITY** |
|---|---|---|---|---|---|---|---|---|---|---|
| **`h8bbar` — THE CONTROL** | 2592 | **TARGET** | 5.0e-7 | 0.50 | **0.0500** | **138.462** | **0.9968** | −0.896 | **yes** | **YES** |
| **`h20uri` — the quadratic** | 8814 | **WALL** (wall-clock) | 2.5e-7 | 0.25 | 0.0138 | 97.552 @ 0.0125 | 0.7023 | (see §5.3) | NO | **no — ALLOWANCE: wall-clock** |
| `h20uri` (withdrawn, §5.4) | 8814 | ITERCAP | 5.0e-7 | 0.50 | 0.0160 | 93.274 @ 0.0100 | 0.6715 | — | NO | no — ALLOWANCE: DR instability |
| `h8bbar`, STEPPED (§4) | 2592 | TARGET | — | 0.50 | 0.0600 | — | 0.9288 | 0.658 | yes | **no — VOID: −25 % path error** |

### 5.1 THE CONTROL PASSES. DR reproduces the answer we already knew.

| `h8bbar` h0 = 0.5 | q/q_exact | source |
|---|---|---|
| merged R3 gate (#722) | 0.9938 | static Newton, plateau |
| note 82 driver | 0.9977 | static Newton, plateau |
| in-tree static leg `qpd_h8bbar_h5.csv` | 0.9977 | static Newton, s/B 0.0744 |
| **this note, tangent-free DR** | **0.9968** | **continuous DR, final HOLD at s/B 0.05** |

`MODE = TARGET`, tail dq/ds = **−0.896 %** of initial ⇒ **PLATEAU**, free advance
⇒ **CAPACITY** on all three clauses. And it is not just the headline that agrees
— **all ten hold checks land on the static Newton curve, and the agreement
tightens monotonically as the plateau forms:**

| hold s/B | 0.005 | 0.010 | 0.015 | 0.020 | 0.025 | 0.030 | 0.035 | 0.040 | 0.045 | 0.050 |
|---|---|---|---|---|---|---|---|---|---|---|
| relaxed q (kPa) | 81.578 | 105.395 | 120.650 | 130.804 | 136.614 | 137.730 | 138.162 | 138.319 | 138.378 | 138.462 |
| **vs static Newton** | −0.45 % | −0.39 % | −0.36 % | −0.35 % | −0.14 % | −0.12 % | −0.08 % | −0.08 % | −0.08 % | **−0.05 %** |
| inertial lag | +0.99 % | +0.37 % | +0.95 % | +0.66 % | +0.66 % | +0.71 % | +0.47 % | +0.28 % | +1.01 % | +0.48 % |
| DR steps to settle | 4000 | 4000 | 4000 | 4000 | 4000 | 4000 | 4000 | 4000 | 4000 | 4000 |

Every hold settled in **4000 steps**, flat across the whole path — the linear
element's relaxation does not get harder as the mechanism forms. **DR is
validated on this problem**, and everything below is therefore interpretable.

### 5.2 THE BINARY ANSWER: the wall at s/B ≈ 0.010 is NOT reproduced by the tangent-free path

The static Newton leg on this exact mesh (`qpd_h20uri_h5.csv`, note 82 §5) dies
at **s/B = 0.01009, q = 95.084 kPa**, on the step floor, still hardening. Under
DR, the same element on the same mesh walks through that settlement without
anything happening there at all:

| hold s/B | 0.00875 | **0.01000** | **0.01125** | **0.01250** |
|---|---|---|---|---|
| relaxed q (kPa) | 89.445 | **93.467** | **97.011** | **97.552** |
| settled? | yes, 4000 steps | yes, 6000 | yes, 6000 | yes, 14000 |
| vs static Newton | −1.50 % | −1.31 % | *(past the wall)* | *(past the wall)* |

> **The quadratic element reaches genuine, gate-passing STATIC EQUILIBRIA at
> 1.11× and 1.24× the settlement where the Newton path died, and nothing
> singular happens at the wall.** So the s/B ≈ 0.010 wall was **in the Newton
> path**, not in the discretization.

**An independent confirmation across a 2× knob change.** The withdrawn
`dtfac = 0.50` leg (§5.4) — a different safety factor *and* a different ramp
rate — reaches the same settlement independently and settles at **93.274 kPa**
against this leg's **93.467 kPa** at s/B = 0.01000: **0.21 %**. Two DR
configurations that disagree about almost every controller setting agree to two
parts in a thousand on the answer, and both agree with the static Newton leg to
~1.4 % right up to the point where Newton stopped being available.

### 5.3 What this note DOES NOT claim, and the observation that must not be over-read

**The quadratic leg did not finish the problem.** It stopped at s/B = 0.0138 on
the **wall clock**, against the linear control's 0.0500 and its plateau. By the
capacity rule that is a **CONTROLLER ALLOWANCE and the allowance is wall-clock
time**; `0.7023` is not an element ceiling and must not be quoted as one. The
full binary question — *does a quadratic element FINISH this problem* — is
therefore answered only in its first half: **the wall is gone; whether the
element plateaus at the exact answer is untested.**

**And there is a real signal in the last three holds that this note deliberately
declines to interpret.** The relaxed sequence is 93.467 → 97.011 → 97.552 kPa:
the increment per hold falls from **+3.54 to +0.54 kPa**, a ~6.5× drop in slope
over one step of the hold ladder, while the cost of settling rose from 6000 to
14000 DR steps and the inertial lag on the ramped sample rose from +2.3 % to
+9.4 %. Two readings fit all of it:

- the element is approaching **its own limit load** at ≈ 0.70 of exact — which
  would make the static campaign's 0.6846 much closer to a real ceiling than
  note 82 was willing to say; or
- the **relaxation is degrading** for the same reason §5.4 documents, and the
  flattening is the leading edge of that, not mechanics.

A third measurement leans, if anything, toward the second reading and is
reported for that reason: the **residual carried by the ramped samples**
(`res_kN`, the true static unbalance during the continuous push) is 0.10 kN at
s/B = 0.010, 0.56 at 0.01125, 1.10 at 0.0125 and 1.71 at 0.01375 — an order of
magnitude of growth over a quarter of a percent of settlement, on a rate that
was quasi-static everywhere before it. The linear control's ramp residual is
**flat at 0.03–0.07 kN across its entire path**, plateau included. So the ramp
rate that was adequate for the quadratic element up to the wall is *not*
adequate past it, which is exactly the confound that makes the flattening
unreadable at three points.

**Three points cannot separate those**, the second reading has a measured
precedent 0.0025 of s/B further along, and note 82's standing rule is that two
runs stopped by the same mechanism agree for reasons unrelated to physics. So it
is recorded as the sharpest open item (§7 item 1) and nothing is concluded from
it. What *is* established is the part that does not depend on the reading: at
s/B = 0.01125 and 0.01250 the state is a settled static equilibrium, and the
Newton path could not get there.

### 5.4 The withdrawn quadratic leg, and why it is reported rather than deleted

The first quadratic leg ran at `dtfac = 0.50` — a safety factor **measured
adequate on the elastic state** (§3: the h20uri zero-push hold is clean at 0.50
and fails at 0.75). It produced settled, static, static-Newton-agreeing holds at
s/B = 0.005 (−1.92 %) and 0.010 (−1.51 %), and then at s/B = 0.015 its hold
**failed to settle in 40 000 steps** while q fell from 125.863 to 91.625 and was
still falling — below its own value at s/B = 0.010. That is not slow
convergence; it is §3's instability, appearing only once the material is deeply
plastic.

**The lesson, and it is the second time this note pays for it:** the DR safety
factor must be validated on the **plastic** state, not the elastic one. A
zero-push hold at the surcharge (control DR-0) certifies the elastic operator
and says nothing about the operator 15 mm of settlement later, when the tangent
has softened and `M*` — rebuilt from it at every KE peak — has softened with it.
Everything after that leg's s/B = 0.0125 is withdrawn; everything before it is
retained, and is what §5.2 uses as its independent cross-check.

## 6. Knob sensitivity — the brief's requirement, four ways

**6.1 `dt` — a predicted EXACT null, and it holds.** With `-mass gershgorin`,
`m ∝ dt²`, so the leap-frog update `du = dt²·a = dt²·R/m = 4R/Σ_j|K_ij|` has
`dt` cancel exactly, and `KE = ½vᵀM*v ∝ du²` is dt-free as well. **The march
must be bit-identical for any `dt`.** Measured (`h8bbar` h0 = 1.0, `dtfac` 0.5,
6 increments, `dt = 1.0` vs `dt = 0.1`): the `s`, `s/B`, `R`, `q`, `dr_iters`
and `settled` columns are **bit-identical**; `ke` and `res` agree to ~7
significant figures (floating-point non-associativity in the `dt` multiply).
**`dt` is not a knob at all — only the ratio `dtfac` is**, which is why §3 is
expressed in `dtfac` throughout.

**6.2 `-noAutoRefresh` — a measured null.** Auto-refresh rebuilds `M*` from the
current (softened) tangent at each KE peak, which was the leading suspect for
§3's element-dependence. Measured (`h20uri` h0 = 0.5, ramped, `dtfac` 0.5,
3 increments): q = 26.420 / 38.710 / 46.043 with refresh on, 26.420 / 38.716 /
46.059 with it off — **identical to four significant figures**. It is not the
mechanism, and it is not a lever on the answer.

**6.3 `dtfac` — NOT a null, and the reason matters.** Measured (`h20uri`
h0 = 0.5, ramped 50 × 100, first two increments): `dtfac` 0.50 gives 26.420 /
38.710; `dtfac` 0.25 gives 26.172 / 37.703 — **−0.94 % and −2.60 %**. The sign is
the one §4.2 predicts and the mechanism is the same one: relaxation progress per
step scales as `dtfac²`, so at a fixed ramp *in steps* a smaller `dtfac` is a
**faster** ramp relative to the dynamics, and buys more spurious softening.
`dtfac` and ramp rate are therefore **not independent knobs**, and a `dtfac`
sweep at fixed `--rampevery` is not a clean stability check. Reported as a
coupled pair, not as two.

**6.4 ramp rate — the binding knob, quantified in §4.2** and carried as a named
~0.6 % systematic rather than treated as converged.

**6.5 The knob check that matters most: does the QUADRATIC answer move?** The two
quadratic legs differ in **three** controller settings at once — safety factor
(0.50 vs 0.25), ramp rate (5.0e-7 vs 2.5e-7 m/step) and hold spacing — and agree
at the one settlement they share to **0.21 %** (93.274 vs 93.467 kPa at
s/B = 0.01000). Both also agree with the static Newton curve to ~1.4 % there.
**On the far side of the DR-controller knobs the answer does not move**, which is
what makes §5.2's claim about the wall a measurement rather than a setting.

**6.6 The knobs that DO move it, and by how much**, collected so no number in
this note is quoted without its exposure:

| knob | effect on q | verdict |
|---|---|---|
| `dt` (at fixed `dtfac`) | **exactly zero — bit-identical** | not a knob |
| `-noAutoRefresh` | < 0.02 % | null |
| ramp: **stepped vs ramped** | **−24 % to −29 %** | dominant; §4 |
| ramp rate 5.0e-7 → 2.0e-7 (sub-jump) | −1.7 % → −0.6 % (linear) | converging |
| ramp rate 2.0e-6 → 5.0e-7 (continuous) | ±2.6 % → ±1 % (linear); ±16 % → ±1 % (quadratic) | converging; **quadratic far more sensitive** |
| `dtfac` 0.50 → 0.25 at fixed ramp-in-steps | −0.9 %, −2.6 % | **coupled to rate, not independent** |
| `dtfac` above the stability threshold | unbounded — the run walks away | fatal; §3 |

## 7. What came back inconclusive

Reported deliberately; three retractions on this problem have cost more than the
caution would have.

1. **THE BIGGEST ONE: does the quadratic element plateau, and where?** §5.3's
   last three holds flatten sharply (93.467 → 97.011 → 97.552 kPa) while the
   relaxation cost more than doubles. Real ceiling near 0.70, or degrading
   relaxation? **Test:** continue the `dtfac = 0.25` leg to s/B ≈ 0.05 with a
   hold every 0.00125 and the safety factor re-validated *on the plastic state*
   (§5.4) — i.e. a zero-*increment* hold repeated at each settlement, which is
   control DR-0 moved down the path. Cost is the honest problem: at the rate
   this note used, ~6 leg-hours for the quadratic alone, and the cost scales as
   `s_total / rate`. **This is the single experiment that would close the binary
   question**, and it is a wall-clock commitment, not a method question.
2. **Whether the linear control's own plateau is rate-converged.** Its holds
   agree with static Newton to −0.05 % at the plateau and −0.45 % early, so the
   bias is small and *shrinking*; but only one rate was run to the plateau.
   **Test:** repeat the control at 2.5e-7 and compare the final hold. Cheap
   relative to item 1 and it would put an error bar on 0.9968.
3. **Why the Gershgorin bound is violated at `dtfac` 0.85 at all.** The triangle
   inequality makes the element-wise row-sum accumulation a *conservative*
   bound on the assembled `Σ_j|K_ij|`, so `dtfac < 1` should be sufficient and
   measurably is not (0.85 fails, 0.75 passes on the linear hex). Auto-refresh is
   excluded (§6.2). Not identified. **Test:** assemble `K` once through
   `FullGeneral`, compute `λ_max(M*⁻¹K)·dt²` directly and compare against 4.
   Cheap, and it would turn a measured safety factor into a derived one.
4. **The residual ramp-rate bias.** Direction and mechanism are
   established and it converges monotonically, but the extrapolated
   zero-rate limit was not measured. **Test:** three rates at one matched
   settlement, Richardson-extrapolated. A few leg-hours.
5. **Whether the DR result carries note 82 §7.1.1's run-to-run scatter.** DR with
   `system Diagonal` involves no threaded factorization, so the FP-determinism
   argument that explained a 5.1 % spread between two bit-identical *static*
   walled legs does not apply — DR should be deterministic. **Untested here.**
   **Test:** repeat one DR leg and diff the CSV; a bit-identical repeat would be
   a genuinely stronger reproducibility footing than the static campaign has.
6. **`h20std` (the locked quadratic) was not run under DR.** Depth was spent on
   the relieved leg. The static campaign has it at 0.4587 (FLOOR).
7. **The identity of note 82 §7.3's abrupt tangent event is untouched here.** DR
   removes the tangent from the path; it does not explain what the tangent was
   doing. That remains note 82 §8 item 1.

<!-- SCOPE-START -->
## 8. Scope

One problem (Prandtl–Reissner), one flow rule (non-associated, ψ = 0), one
material (`DruckerPrager`, φ_ps = 27.47°, σ_y = 0.2 kPa regularizer), one mesh
family (structured graded strip, plane strain), one resolution (h₀ = 0.5), one
displacement-controlled path family, one fictitious-mass mode
(`-mass gershgorin`) and one damping mode (Cundall kinetic). Nothing here
carries to associated flow, to dynamics, or to bending-dominated problems.

**Two scope limits specific to this note.** (i) The quadratic leg reaches
s/B = 0.0138 against the control's 0.0500 — every quadratic statement here is
about the first quarter of the path, and the only *element* claim it supports is
the negative one in §5.2. (ii) The static Newton curve, which anchors every
tracking check in §4.3 and §5, **exists only up to s/B = 0.01009 for the
quadratic element** — because that is where it died. Beyond it the quadratic DR
leg is unanchored by construction, and no cross-check available in this campaign
can validate it there. That is not a gap this note could have closed; it is the
shape of the problem.

## 9. Reproducing

```bash
cd Ladruno_files/testbed/hypo_bearing
py -3.12 dr_elastic_check.py                                   # sec 4.1
py -3.12 quad_dr_probe.py --elem h8bbar --h0 0.5 --sfrac 0.06 --ds 1e-3 \
        --ramp 50 --rampevery 100 --dtfac 0.5 --suffix _RAMP   # THE CONTROL
py -3.12 quad_dr_probe.py --elem h20uri --h0 0.5 --sfrac 0.06 --ds 1e-3 \
        --ramp 50 --rampevery 100 --dtfac 0.5 --suffix _RAMP   # the quadratic
py -3.12 quad_dr_probe.py --elem h8bbar --h0 1.0 --dt 0.1 --suffix _dtnull0.1
py -3.12 quad_dr_summary.py
```

Every leg must be quoted with its termination mode, its `dtfac` and its ramp
rate. A leg quoted as a capacity must have reached `TARGET` with a flat tail and
free advance.
<!-- SCOPE-END -->
