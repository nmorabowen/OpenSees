# Note 81 — What is the TIMs limit-load ceiling a property of? A quadratic-hex measurement

**Kind:** MEASUREMENT. This note reports a number and a verdict. It ships no
feature and changes no element.
**Harness:** `Ladruno_files/testbed/hypo_bearing/h20_prandtl.py` (+ `h20_summary.py`)
**Engine:** `ladrunoBuild() = 1db3394b8c4d4c23e12f2629d5e7ef92a106acfa`
**Problem:** Prandtl–Reissner, `q_u = q0·N_q` exact, φ_txc = 20°, q₀ = 10 kPa,
γ = 0, ν = 0.45, Chen & Han **plane-strain** cone match, non-associated (ρ̄ = 0)
gate + associated control.

---

## 1. The question this was run to answer

The TIMs element-technology campaign measures the one bearing-capacity problem
whose limit analysis has coincident upper and lower bounds — a strip footing on
**weightless** soil under a surcharge, `q_u = q0·N_q`, no shape factor, no
`N_γ`, nothing to calibrate — with different element technologies on the same
problem. Its interim numbers:

| element | basis | order | volumetric treatment | q_max/q_exact |
|---|---|---|---|---|
| `LadrunoBrick -formulation bbar` | Lagrange hex | linear | B-bar | **≈ 0.99**, plateau |
| `TenNodeTetrahedron` | Lagrange tet | quadratic | none | 0.2828 |
| `BezierTet10` std | Bernstein tet | quadratic | none | 0.3075 |
| `BezierTet10 -bbar` | Bernstein tet | quadratic | B-bar | ≈ 0.49–0.55, climbing |

B-bar is worth roughly 0.31 → 0.50 on the quadratic tet — real, and far short
of the linear hex's 0.40 → 0.99. Three explanations survive that data set:

- **(a) simplex geometry** — one volumetric constraint per element is simply
  weaker on tets (constraint counting is adverse on tet meshes), so *any* tet
  tops out, and the fix direction is an F-bar-Patch macro for simplices;
- **(b) quadratic order** — something about quadratic elements under
  near-incompressible perfectly-plastic flow;
- **(c) the Bézier/Bernstein basis** specifically.

The campaign's own data nearly excludes (c) already: at equal order and equal
geometry the Bernstein tet (0.3075) and the Lagrange tet (0.2828) agree to
within 8 %, an order of magnitude smaller than the gap being explained.

**The cell nobody had run is a quadratic HEX.** `LadrunoBrick20` is quadratic,
hexahedral, and Lagrange (serendipity), so it changes exactly one factor at a
time against both the linear hex and the quadratic tets.

## 2. Design

`LadrunoBrick20` has no B-bar, and does not need one. Its two formulations are
`std` (full 27-pt Gauss) and `uri` (uniform reduced 2×2×2) — and H20 + 2×2×2 is
the classical near-incompressible workhorse, the C3D20R of Abaqus. Reduced
integration at quadratic order plays the role B-bar plays at linear order, so
`uri` is the leg technology-matched to the ≈ 0.99 baseline and `std` is the
locked control. ADR 72 §3.4 already computed the constraint ratios, which is
the standard way to predict this:

| element | rule | added DOF / element | constraints / element | ratio r | ideal |
|---|---|---|---|---|---|
| H8 | B-bar | 3 | 1 | **3.0** | 3 |
| H8 | full 2×2×2 | 3 | 8 | 0.38 | 3 |
| H20 | `uri` 2×2×2 | 12 | 8 | **1.5** | 3 |
| H20 | `std` 3×3×3 | 12 | 27 | 0.44 | 3 |

So the design is a 2×2 (order × volumetric treatment) rather than a single
point, and it comes with a prediction: `std` should lock hard, `uri` should
relieve most of it but sit below the linear hex's ratio-optimal 3.0.

Every problem parameter is the reference driver's
(`fork_bundle/r3_prandtl_tet10.py`), unchanged: the cone and its plane-strain
match, q₀, γ = 0, ν = 0.45 (**not** lowered — at a low ν this cone starts at
95 % of yield and the leg measures its own initial condition), the
displacement-controlled push to s/B, the ADR-63 D16 adaptive ladder and its
three guards, the 2 %-drop truncation, and the plateau test (tail dq/ds < 2 %
of initial). The mesh is `dp_strip.py`'s own graded strip grid — 14.5 B of side
clearance, 10 B of depth — so it is the grid the linear-hex baseline was
measured on.

The harness builds that one grid at either order from a half-index lattice (a
site is a node iff at most one of its three half-indices is odd), so
`--order 1` gives the `LadrunoBrick` H8 baseline and `--order 2` the
`LadrunoBrick20` H20 **on the same grid, through the same code path**. The
linear-hex baseline was therefore re-run inside this harness rather than quoted
from the earlier campaign, which removes the driver as a variable.

## 3. Controls

The campaign that produced the numbers in §1 was burned three times by silent
artifacts, so this one opens with a control battery and every leg is refused if
any control fails.

**Control 0 — provenance.** `ops.ladrunoBuild()` is queried and recorded before
anything is built: `1db3394b8c4d4c23e12f2629d5e7ef92a106acfa`.

**Control 2 — the load is consistent with the basis.** A uniform pressure on a
quadratic serendipity face has **negative corner weights**: −A/12 at the four
corners, +A/3 at the four mid-edges. The tributary-area lumping that every mesh
generator in this testbed hands you is correct for an H8 face and *wrong* here.
The surcharge is therefore obtained by 3×3 Gauss integration of the Q8 shape
functions over the real face geometry and checked against the analytic values
(agreement to 3.1e-16 of A), and its sum against the physical top area
(30.000000000 m² vs 30.000000000 m²).

**Control 3 — the base reaction reproduces the applied total.** +0.0000000 %,
i.e. below the 1e-6 gate.

**Control 1 — the 1-D elastic patch test.** Controls 2 and 3 test the *sum*;
this is the one that tests the *distribution*. Roller sides, fixed base and a
consistent uniform surcharge admit the exact 1-D state, and that state lives in
the H20 space, so a consistent load vector must reproduce it at every Gauss
point. Measured: `max |σ_zz/(−q0) − 1| = 7.75e-14` (uri, 8 GP), `1.06e-13`
(std, 27 GP), `4.49e-14` (H8 bbar); shear below 1e-13 throughout.

> **The negative control, which is the reason to trust the above.** Feeding the
> same model a *tributary-lumped* surcharge reproduces the applied total in the
> base reactions to **+0.0000000 %** — so the usual sum identity passes it
> without complaint — while putting **190 %** error into σ_zz. This is the exact
> defect class that made a Bézier deck diverge at first yield, and the sum
> identity cannot see it.

**Control 4 — the spurious-mode census.** H20 at 2×2×2 carries 6 spurious modes
per element and ships with **no hourglass control by design** (ADR 72 §3.2,
following C3D20R). The contract is that they are non-communicable in a solid
assembly — and the usage guidance shipped with it says `uri` wants a mesh **at
least two elements thick**, while this plane-strain strip is **one** element
thick in y. The leg therefore sits exactly on the documented caveat, so the
rank is measured rather than argued: the restrained elastic tangent is
assembled (`FullGeneral` + `printA -ret`) and its eigenvalues taken.

Measured on both meshes, `uri`:

| mesh | free DOF | λ_min / (tr K / n) | modes below 1e-10 |
|---|---|---|---|
| h₀ = 1.0 | 2800 | 2.386e-04 | **0** |
| h₀ = 0.5 | 5460 | 1.061e-04 | **0** |

Full rank on both. The per-element modes do not communicate through the
single-element-thick strip, so the `uri` numbers are measurements and not
mechanisms. (The h₀ = 0.5 patch test runs to 2.48e-13 as well.)

**Control 5 — the subdivision budget is pinned and reported.** A budget-starved
leg reads "still hardening" and is not a ceiling (the reference campaign's note
71: 0.8988 at budget 24 against a 0.9513 plateau once raised). Every leg below
reports the subdivisions it spent against the budget it was given.

## 4. Results

All legs: non-associated (ρ̄ = 0), the gate flow rule, unless marked. Every H20
leg terminated on the **step floor**, none on the subdivision budget — the
`subdiv` column is against a budget of 800.

| leg | element | rule | r | el/B | q_max/q_exact | plateau | tail % | s_end/B | ended by | subdiv | steps | wall s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `h8_nonassoc_bbar_h5` | H8 | B-bar | 3.0 | 4 | **0.9977** | **yes** | −0.00 | **0.2500** | reached the target | 315/500 | 1909 | 2740 |
| `h20_nonassoc_uri_h1.0` | H20 | `uri` | 1.5 | 2 | **0.7715** | NO | 11.6 | 0.0112 | step floor | 49/800 | 260 | 266 |
| `h20_nonassoc_uri_h5` | H20 | `uri` | 1.5 | 4 | **0.6894** | NO | 10.4 | 0.0103 | step floor | 67/800 | 399 | 994 |
| `h20_nonassoc_std_h1.0` | H20 | `std` | 0.44 | 2 | **0.5894** | NO | 27.2 | 0.0045 | step floor | 47/800 | 252 | 310 |
| `..._uri_h1.0_strong` | H20 | `uri` | 1.5 | 2 | 0.7902 (100× rung, §4.2) | NO | 11.2 | 0.0120 | step floor | 46/800 | 305 | 1014 |
| `h20_assoc_uri_h1.0` | H20 | `uri` | 1.5 | 2 | 1.1422 (assoc control) | see §4.5 | 0.02 | 0.0359 | stopped by hand | 0/800 | 2265 | 2265 |

**The baseline ran the full push.** `h8_nonassoc_bbar_h5` reached the s/B = 0.25
target — not a wall, not a cap — with q_max = 138.59 kPa against the exact
138.91, a tail slope of −0.00 kPa/m, and its peak already reached at
s/B = 0.144. That is what a collapse load looks like on this problem, and it
reproduces the campaign's own 1.0020 / ≈0.99 inside this harness. No H20 leg
looks like it.

The associated control was stopped by hand once it had made its point
(§4.5); its row is scored from its CSV by `h20_summary.py`, which applies the
identical truncation and plateau rules.

Wall times are on a **contended box** (three other agents were running); treat
them as an order of magnitude, not a benchmark.

### 4.1 The H20 numbers are not step-size artifacts, and not budget artifacts

Each H20 leg was run on **two independent step ladders** a factor of 10 apart
in both base and floor — `dp_strip.py`'s 2e-4/2e-6/1e-3 and the quadratic
driver's 2e-5/2e-7/2e-4:

| leg | coarse ladder | fine ladder | agreement |
|---|---|---|---|
| H20 `uri`, 2 el/B | 0.7706 at s/B = 0.0112 | 0.7715 at s/B = 0.0112 | **0.1 %** |
| H20 `std`, 2 el/B | 0.5864 at s/B = 0.0044 | 0.5894 at s/B = 0.0045 | **0.5 %** |

Same wall, same settlement, three-digit agreement. This is the control the
brief demanded (note 71's 0.8988-vs-0.9513): the number does not move when the
stepping is loosened or tightened, and the budget was never the binding
constraint (49, 67 and 47 subdivisions spent of 800).

### 4.2 The wall is the discretization, not Newton

The legs fail by **residual stagnation**: 0.045 kN on a 300 kN problem, 1.5e-4
relative, and only 1.5× above the tolerance of the 10× rung that refused it. At
that margin "Newton cannot get there" and "there is nothing there to get to"
look identical, so a fourth rung at **100× tolerance** was added (`--strong`)
and the leg re-run.

It bought **0.7715 → 0.7902, i.e. +2.4 %**, ended at s/B = 0.0120 (against
0.0112) still hardening at 11.2 % of its initial tangent, and again died on the
step floor with 46 of 800 subdivisions. It paid for those 2.4 % with **95 of
305 steps taken on the 100× rung** — nearly a third of the curve at a
tolerance a hundred times looser than the one the headline legs were held to.

A third of the curve on a 100× relaxation moves the answer by under three
percent and still does not produce a plateau. **The wall is a property of the
discretization; a better solver does not walk through it.**

### 4.3 What the volumetric treatment is worth, by order

This is the measurement the 2×2 design was built for:

| order | locked | relieved | gain | ends on a plateau? |
|---|---|---|---|---|
| linear hex | 0.40 (TIMs) | **0.9977** (B-bar, r = 3.0) | ×2.5 | **yes** |
| quadratic hex | **0.5894** (`std`, r = 0.44) | **0.7715** (`uri`, r = 1.5) | ×1.31 | no |
| quadratic tet | 0.2828 / 0.3075 (TIMs) | ≈ 0.49–0.55 (B-bar, TIMs) | ×1.7 | no (climbing) |

Volumetric relief helps at every order. **It only finishes the job at linear
order**, and that is exactly what the constraint ratio predicts: B-bar on an H8
lands on the ideal r = 3.0, while the best relief available on a quadratic hex
(2×2×2 reduced) lands at r = 1.5, half of ideal — ADR 72 §3.4 computed both
before this measurement was taken.

### 4.4 Mesh refinement

| element | 2 el/B | 4 el/B | 8 el/B | 16 el/B |
|---|---|---|---|---|
| H20 `uri` (here) | 0.7715, wall at s/B = 0.0112 | 0.6894, wall at s/B = 0.0103 | — | — |
| H8 B-bar (here / COLLAPSE.md) | — | **0.9977** / 1.0020 | 0.9517 | 0.9260 |

Both families approach from **above** and fall with refinement, so neither
number is a converged limit load and the H20 values are not conservative in the
mesh sense. But the two behave qualitatively differently: the H8 series
*plateaus at every resolution* (a proper convergence sequence against an exact
answer, COLLAPSE.md), whereas the H20 series loses reach as it refines — the
wall arrives *earlier* in settlement, 0.0112 → 0.0103. That is the same
"monotone loss of reach" signature COLLAPSE.md attributes to
Rudnicki–Rice localization under non-associated perfect plasticity.

### 4.5 The falsification control behaves as the campaign said it would

The associated leg runs *past* the exact answer — 1.1422 — which is the
documented behaviour (dilating at ψ = φ must push the surrounding soil aside; a
bounded mesh resists, and the leg hardens through its own oracle). Its
`plateau = yes` flag should **not** be read as a collapse load: the test looks
at the last 10 % of a curve that only reached s/B = 0.036, and the value it
flattens at is 14 % *above* a bound that is exact from both sides. It is
reported here as what it is — a control that fired, confirming the
non-associated leg is the one carrying the measurement.

## 5. Verdict

**The pre-registered rule was: quadratic hex ≈ 1.0 ⇒ (a) simplex geometry;
quadratic hex tops out well below 1.0 ⇒ (b) quadratic order. The quadratic hex
tops out at 0.69–0.78. The data supports (b) QUADRATIC ORDER.**

Three strands, in order of how much weight they carry.

1. **The qualitative split is on ORDER, not on geometry.** The linear hex
   plateaus — cleanly, at every resolution, on the same problem, the same mesh,
   the same script, 0.9977 of an exact answer, running the full s/B = 0.25
   push to completion. *Every* quadratic leg fails to:
   hex and tet, Lagrange and Bernstein, relieved and locked, all of them wall
   while still hardening at 10–27 % of their initial tangent. The plateau, not
   the number, is the robust observable here, and it partitions the elements by
   order.
2. **Volumetric relief does not close the gap at quadratic order.** ×2.5 and a
   plateau at linear order; ×1.31 and a wall at quadratic order. The
   constraint-ratio account predicted this in advance (r = 3.0 vs r = 1.5) and
   it holds.
3. **Geometry is real but second-order.** At matched order and matched relief
   the hex beats the tet by roughly 2× (0.5894 vs 0.2828 locked; 0.77 vs ≈ 0.52
   relieved). That is a genuine simplex penalty — it is simply much smaller than
   the 0.28 → 0.99 gap being explained, and it does not change the plateau
   verdict.

**(c) the Bézier/Bernstein basis is excluded.** It was already nearly excluded
by the 8 % gap between the Bernstein and Lagrange tets; this note reproduces
the entire ceiling inside pure Lagrange elements, with no Bézier anywhere.

### 5.1 The honest limitation, stated plainly

**The H20 numbers are lower bounds, not ceilings.** A leg that walls while
hardening at 10 % of its initial tangent has not shown you where it would have
stopped. What is established is (i) that the wall is reproducible to three
digits across a 10× change in step size, (ii) that it is not the subdivision
budget, and (iii) that a 100×-relaxed solver moves it by 2.4 %. What is *not*
established is that 0.77 is the quadratic hex's asymptote. If a regularized
formulation walked an H20 to 1.0 on this problem, verdict (b) would weaken and
(a) would come back into play.

Two things would settle it, neither cheap: a **regularized** run (rate-
dependent/viscoplastic, gradient, or Cosserat) or an **explicit dynamic** solve,
which removes the singular-tangent wall entirely; and the **H27 + selective**
element, which is the ratio-optimal quadratic brick (r = 3.0, ADR 72 §3.4) and
so is the direct test of whether quadratic order *per se* is the problem or
whether it is merely that no quadratic element in the fork reaches r = 3.

That second one is the sharper experiment, and it is worth naming: **the
present result cannot separate "quadratic order is intrinsically bad here" from
"we have no quadratic element with an optimal constraint ratio."** H27 +
selective would separate them. Until it is run, verdict (b) should be read as
*quadratic order as currently implemented*, which is the engineering-relevant
statement but not the deeper one.

### 5.2 Recommendation on an F-bar-Patch tet lane

The brief asked for this recommendation conditional on verdict (a). The verdict
is (b), so the conditional does not fire in its stated form — but the data
speaks to the question directly, and the answer is **not the lane as
originally scoped**:

- **Do not open F-bar-Patch on QUADRATIC tets.** The evidence says the
  achievable prize is small. Raising the constraint ratio at quadratic order is
  exactly what `uri` does on the H20, and it bought 0.59 → 0.77 without a
  plateau. A quadratic-simplex macro is predicted to land in the same
  neighbourhood, and it would cost a new element to get there.
- **Do open it on LINEAR tets, if it is opened at all.** Linear + optimal
  volumetric relief is the *only* combination this problem has ever been solved
  by — 0.9977 with a plateau — and the fork already owns the 2D form of the
  macro (`LadrunoCSTPair`, dSNPO §15.1.9). An F-bar-Patch Tet4 is the 3D
  analog of a thing that already works here, aimed at the one cell of the 2×2
  that is measured to succeed. That is a far better-supported bet than the
  quadratic-simplex version, and it is the lane that would let unstructured
  meshing reach the accuracy the structured hex mesh already has.
- **Higher value than either, for the fork's own elements:** ADR 72 already
  identifies **H27 + selective/B-bar** as the ratio-optimal quadratic brick.
  It is both the missing scientific control (§5.1) and the only route to a
  quadratic element that could plausibly plateau on this problem. If one
  element gets built off this note, that is the one.

### 5.3 Scope

One problem (Prandtl–Reissner), one flow rule (non-associated, ψ = 0), one
material (`DruckerPrager`, φ_ps = 27.47°, σ_y = 0.2 kPa regularizer), one
mesh family (structured graded strip, plane strain), two resolutions on the
H20 and one on the H8. Nothing here carries to associated flow, to severe
non-associativity (the φ_ps = 53.7° cone, where COLLAPSE.md shows *every*
element loses reach), to dynamics, or to bending-dominated problems — where
quadratic elements are expected to and do win, which is why they exist.

## 6. Reproducing

```bash
cd Ladruno_files/testbed/hypo_bearing
py -3.12 h20_prandtl.py patch modes --h0 0.5 --order 2 --form uri
py -3.12 h20_prandtl.py nonassoc --h0 0.5 --order 2 --form uri --budget 500
py -3.12 h20_prandtl.py nonassoc --h0 0.5 --order 1 --form bbar --budget 500
py -3.12 h20_summary.py
```
