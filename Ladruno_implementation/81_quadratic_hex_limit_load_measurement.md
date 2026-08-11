# Note 81 — Quadratic elements lose the Newton path before their mechanism forms: a quadratic-hex measurement on the Prandtl oracle

**Kind:** MEASUREMENT. Ships no feature, changes no element.
**Status:** **CORRECTED.** The measurements below stand unchanged. The
*explanation* in the first version — that a constraint-ratio (`r`) deficit at
quadratic order was the mechanism — is **withdrawn and recorded as a rejected
hypothesis** (§5.3). Two independent findings from the TIMs campaign forced it,
and re-checking my own legs against them turned up an over-claim of mine as
well (§4.2). What changed is in §0.
**Harness:** `Ladruno_files/testbed/hypo_bearing/h20_prandtl.py` (+ `h20_summary.py`)
**Engine:** `ladrunoBuild() = 1db3394b8c4d4c23e12f2629d5e7ef92a106acfa`
**Problem:** Prandtl–Reissner, `q_u = q0·N_q` exact, φ_txc = 20°, q₀ = 10 kPa,
γ = 0, ν = 0.45, Chen & Han **plane-strain** cone match, non-associated (ρ̄ = 0)
gate + associated control.

> **ATTRIBUTION AND HOLD.** Results credited below to the **TIMs
> element-technology campaign** are theirs, cited here with credit. TIMs
> material is **internal-only**: it is under the same hold as the unfiled UW
> report. Do **not** carry their numbers or their phrasing into anything
> public or upstream-facing — that includes upstream PR bodies, source
> comments destined for upstream files, and any doc ported out of
> `Ladruno_implementation/`. The fork-only measurements in §4 are ours and
> carry no such hold.

---

## 0. What changed in this correction

Three things, in descending order of how much they move the conclusion.

1. **`r` does not predict the behaviour — it is contradicted by a case that
   was already measured.** The TIMs `BezierTet10 -bbar` leg counts at
   **r ≈ 4.4** (≈24·V DOF against ≈5.5·V elements at one volumetric constraint
   each) — *above* the linear hex's curing 3.0 — and it behaves like my H20
   `std` at r = 0.44: no plateau, walls low. Their premise checks out at
   source: `BezierTet10::computeVolumeAveragedDerivatives`
   (`BezierTet10.cpp:1398-1431`) builds `B̄_a,i = (1/V_e)∫∂N_a/∂x_i`, one
   element-wide constant on the full 4-point rule, and `computeBBarMatrix`
   (`:1363-1395`) is textbook mean dilatation — exactly one volumetric
   constraint per element, deviatoric part untouched at full quadratic
   resolution. So the ordering I read as an `r` correlation is three points
   fitted by an ordering with a different cause. §5.3.
2. **A wall moves when a pure solver parameter changes, so a wall is not a
   ceiling.** TIMs ran two legs differing *only* in the adaptive controller's
   subdivision budget — same mesh, element, loads, material, tolerances,
   algorithm ladder. Budget 24 → 0.4930 of exact, end s/B 0.00399, tail 20.8 %
   of initial. Budget 48 → **0.5779**, end s/B **0.00635**, tail **13.3 %**.
   Doubling one integer moved capacity **+17 %**, reach **+59 %**, and softened
   the tail by a third. Neither leg plateaued; both ended on "budget spent".
   That is the falsification test the first version of this note should have
   applied to its own numbers, and my legs all ended on the analogous limit —
   the step floor. §4.1.
3. **My own two-ladder control was narrower than I reported.** I claimed
   three-digit agreement across two step ladders on three legs. Only **one** of
   those pairs is a genuine wall-against-wall comparison; on the other two the
   coarse-ladder run was **stopped by hand** when I relaunched, so its endpoint
   is my `TaskStop`, not a wall, and the "agreement" is not evidence of
   anything. Corrected in §4.2.

None of this touches the measurements, the control battery, or the two
eliminations that were never mechanism-dependent. It changes what the numbers
are numbers *of*.

## 1. The question this was run to answer

The TIMs element-technology campaign measures the one bearing-capacity problem
whose limit analysis has coincident upper and lower bounds — a strip footing on
**weightless** soil under a surcharge, `q_u = q0·N_q`, no shape factor, no
`N_γ`, nothing to calibrate — with different element technologies on the same
problem. Its numbers (TIMs, internal):

| element | basis | order | volumetric treatment | q_max/q_exact |
|---|---|---|---|---|
| `LadrunoBrick -formulation bbar` | Lagrange hex | linear | B-bar | **0.9934**, plateau |
| `TenNodeTetrahedron` | Lagrange tet | quadratic | none | 0.2828 |
| `BezierTet10` std | Bernstein tet | quadratic | none | 0.3075 |
| `BezierTet10 -bbar` | Bernstein tet | quadratic | mean dilatation | ≈ 0.49–0.58, no plateau |

Three explanations survived that data set:

- **(a) simplex geometry** — one volumetric constraint per element is weaker on
  tets, so *any* tet tops out, and the fix direction is an F-bar-Patch macro
  for simplices;
- **(b) quadratic order** — something about quadratic elements under
  near-incompressible perfectly-plastic flow;
- **(c) the Bézier/Bernstein basis** specifically.

**The cell nobody had run is a quadratic HEX.** `LadrunoBrick20` is quadratic,
hexahedral, and Lagrange (serendipity), so it changes exactly one factor at a
time against both the linear hex and the quadratic tets.

## 2. Design

`LadrunoBrick20` has no B-bar. Its two formulations are `std` (full 27-pt
Gauss) and `uri` (uniform reduced 2×2×2) — and H20 + 2×2×2 is the classical
near-incompressible workhorse, the C3D20R of Abaqus. Reduced integration at
quadratic order plays the role B-bar plays at linear order, so `uri` is the leg
technology-matched to the ≈ 0.99 baseline and `std` is the locked control. That
makes the design a 2×2 (order × volumetric treatment) rather than a single
point.

**The hypothesis I went in with, now rejected (§5.3):** ADR 72 §3.4 computes
constraint ratios, and they order the four cells the way the results came out:

| element | rule | added DOF / element | constraints / element | ratio r | ideal |
|---|---|---|---|---|---|
| H8 | B-bar | 3 | 1 | **3.0** | 3 |
| H8 | full 2×2×2 | 3 | 8 | 0.38 | 3 |
| H20 | `uri` 2×2×2 | 12 | 8 | **1.5** | 3 |
| H20 | `std` 3×3×3 | 12 | 27 | 0.44 | 3 |

The prediction — `std` locks hard, `uri` relieves most of it but sits below the
linear hex — is what happened, and I read that as confirmation. It was not:
§5.3 has a measured counter-example at r ≈ 4.4 that behaves like r = 0.44. The
table is retained because the counting is correct and the design was built on
it; the *inference* from it is withdrawn.

Every problem parameter is the reference driver's
(`fork_bundle/r3_prandtl_tet10.py`), unchanged: the cone and its plane-strain
match, q₀, γ = 0, ν = 0.45 (**not** lowered — at a low ν this cone starts at
95 % of yield and the leg measures its own initial condition), the
displacement-controlled push to s/B, the ADR-63 D16 adaptive ladder and its
three guards, the 2 %-drop truncation, and the plateau test (tail dq/ds < 2 %
of initial). The mesh is `dp_strip.py`'s own graded strip grid — 14.5 B of side
clearance, 10 B of depth.

The harness builds that one grid at either order from a half-index lattice (a
site is a node iff at most one of its three half-indices is odd), so
`--order 1` gives the `LadrunoBrick` H8 baseline and `--order 2` the
`LadrunoBrick20` H20 **on the same grid, through the same code path**. The
linear-hex baseline was therefore re-run inside this harness rather than quoted,
which removes the driver as a variable.

## 3. Controls

These are unchanged by the correction and all of them still hold. The campaign
that produced the numbers in §1 was burned three times by silent artifacts, so
every leg is refused if any control fails.

**Control 0 — provenance.** `ops.ladrunoBuild()` queried and recorded before
anything is built: `1db3394b8c4d4c23e12f2629d5e7ef92a106acfa`.

**Control 2 — the load is consistent with the basis.** A uniform pressure on a
quadratic serendipity face has **negative corner weights**: −A/12 at the four
corners, +A/3 at the four mid-edges. The tributary-area lumping that every mesh
generator in this testbed hands you is correct for an H8 face and *wrong* here.
The surcharge is obtained by 3×3 Gauss integration of the Q8 shape functions
over the real face geometry, checked against the analytic values (3.1e-16 of A)
and its sum against the physical top area (30.000000000 m² vs 30.000000000 m²).

**Control 3 — the base reaction reproduces the applied total.** +0.0000000 %.

**Control 1 — the 1-D elastic patch test.** Controls 2 and 3 test the *sum*;
this tests the *distribution*. Roller sides, fixed base and a consistent uniform
surcharge admit the exact 1-D state, which lives in the H20 space, so a
consistent load vector must reproduce it at every Gauss point. Measured:
`max |σ_zz/(−q0) − 1| = 7.75e-14` (uri, 8 GP), `1.06e-13` (std, 27 GP),
`4.49e-14` (H8 bbar); shear below 1e-13 throughout.

> **The negative control, which is why the above is worth trusting — and which
> TIMs have now independently reproduced.** Feeding the same model a
> *tributary-lumped* surcharge reproduces the applied total in the base
> reactions to **+0.0000000 %** — the usual sum identity passes it without
> complaint — while putting **190 %** error into σ_zz on this mesh. TIMs
> measured the same failure at **343 %** on theirs. The sum identity is
> structurally blind to it; the 1-D elastic stress patch test is the
> discriminator that catches it. This is the defect class that made a Bézier
> deck diverge at first yield.

**Control 4 — the spurious-mode census.** H20 at 2×2×2 carries 6 spurious modes
per element and ships with **no hourglass control by design** (ADR 72 §3.2,
following C3D20R), and its usage guidance wants a mesh **at least two elements
thick** while this plane-strain strip is **one**. So the rank is measured, not
argued: the restrained elastic tangent is assembled (`FullGeneral` +
`printA -ret`) and its eigenvalues taken.

| mesh | free DOF | λ_min / (tr K / n) | modes below 1e-10 |
|---|---|---|---|
| h₀ = 1.0 | 2800 | 2.386e-04 | **0** |
| h₀ = 0.5 | 5460 | 1.061e-04 | **0** |

Full rank on both — the `uri` results are not mechanisms. (The h₀ = 0.5 patch
test runs to 2.48e-13 as well.)

**Control 5 — the path-controller limits are pinned and reported.** Every leg
reports *which* guard ended it and what it spent against its allowance. §4.1 is
what that reporting turned out to be for.

## 4. Results

Non-associated (ρ̄ = 0), the gate flow rule, unless marked.

| leg | element | rule | el/B | q_max/q_exact | plateau | tail % | s_end/B | ended by | subdiv | steps | wall s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `h8_nonassoc_bbar_h5` | H8 | B-bar | 4 | **0.9977** | **yes** | −0.00 | **0.2500** | **reached the target** | 315/500 | 1909 | 2740 |
| `h20_nonassoc_uri_h1.0` | H20 | `uri` | 2 | **0.7715** | NO | 11.6 | 0.0112 | step floor | 49/800 | 260 | 266 |
| `h20_nonassoc_uri_h5` | H20 | `uri` | 4 | **0.6894** | NO | 10.4 | 0.0103 | step floor | 67/800 | 399 | 994 |
| `h20_nonassoc_std_h1.0` | H20 | `std` | 2 | **0.5894** | NO | 27.2 | 0.0045 | step floor | 47/800 | 252 | 310 |
| `..._uri_h1.0_strong` | H20 | `uri` | 2 | 0.7902 (100× rung, §4.3) | NO | 11.2 | 0.0120 | step floor | 46/800 | 305 | 1014 |
| `h20_assoc_uri_h1.0` | H20 | `uri` | 2 | 1.1422 (assoc control) | see §4.6 | 0.02 | 0.0359 | stopped by hand | 0/800 | 2265 | 2265 |

**The baseline ran the full push.** `h8_nonassoc_bbar_h5` reached the s/B = 0.25
target — not a wall, not a cap — with q_max = 138.59 kPa against the exact
138.91, tail slope −0.00 kPa/m, peak already reached at s/B = 0.144. That is
what a collapse load looks like on this problem, and it reproduces the
campaign's 0.9934 independently. **No H20 leg looks like it.**

Wall times are on a **contended box** (three other agents were running); treat
them as an order of magnitude, not a benchmark.

### 4.1 Every quadratic leg ended on a path-controller limit — that is the fact that reframes all of them

The `ended by` column is the important one. The H8 leg ran out of *problem*: it
reached the settlement target with a flat tangent. Every H20 leg ran out of
*path*: the adaptive controller bottomed out at its step floor with the
load–settlement curve still rising at 10–27 % of its initial slope.

Those are different kinds of ending, and only the first produces a number that
means "collapse load". The second produces a number that means "this is where
the solver stopped being able to continue", and such a number is only a lower
bound on the element's capacity **if the solver limit is not itself movable**.

It is movable. TIMs demonstrated it directly by changing one integer:

| leg (TIMs, internal) | budget | q/q_exact | end s/B | tail % of initial | ended by |
|---|---|---|---|---|---|
| tet10 + mean-dilatation | 24 | 0.4930 | 0.00399 | 20.8 | budget spent |
| tet10 + mean-dilatation | 48 | **0.5779** | **0.00635** | **13.3** | budget spent |

Same mesh, same element, same loads, same material, same tolerances, same
algorithm ladder. Doubling the subdivision budget bought **+17 % capacity**,
**+59 % reach**, and cut the tail slope by a third. The curve did not converge
to anything — it just got further along before the controller gave up, and the
tail was still falling, which means the sequence had not finished moving.

My step floor is the same class of limit as their budget. So **0.5894, 0.6894
and 0.7715 are not the quadratic hex's ceilings and are not measurements of a
collapse load.** They are the points at which the path was lost.

### 4.2 Step-ladder sensitivity — corrected, and narrower than I first reported

Each H20 configuration was run on two step ladders a factor of 10 apart in base
and floor: `dp_strip.py`'s 2e-4/2e-6/1e-3 and the quadratic driver's
2e-5/2e-7/2e-4. I originally presented three pairs as agreeing to three digits.
Re-checking terminations, only one pair is a legitimate comparison:

| configuration | coarse ladder | fine ladder | is this a wall-vs-wall pair? |
|---|---|---|---|
| H20 `uri`, 2 el/B | 0.7706, walled at s/B = 0.0112 | 0.7715, walled at s/B = 0.0112 | **yes** — 0.1 % |
| H20 `std`, 2 el/B | 0.5864 — **stopped by hand** | 0.5894, walled | no |
| H20 `uri`, 4 el/B | 0.5373 — **stopped by hand** | 0.6894, walled | no |

The coarse-ladder `std` and 4 el/B runs were killed when I relaunched on the
fine ladder; their endpoints are my `TaskStop`, not a convergence wall.
Comparing a truncated run to a walled one and calling the closeness "agreement"
was an error, and in the 4 el/B row the two numbers are not even close (0.5373
against 0.6894).

**What the surviving pair does and does not show.** It shows that scaling this
particular knob — base and floor together, by 10× — leaves the 2 el/B `uri`
wall where it is. It does *not* show the wall is immovable in general, because
TIMs moved an *independent* knob on the same class of leg and got +17 %. One
knob's insensitivity is not a ceiling.

### 4.3 The 100× tolerance rung

The legs fail by residual stagnation: 0.045 kN on a 300 kN problem, 1.5e-4
relative, only 1.5× above the tolerance of the 10× rung that refused it. A
fourth rung at **100× tolerance** (`--strong`) was added and the leg re-run.

It bought **0.7715 → 0.7902, +2.4 %**, ended at s/B = 0.0120 (against 0.0112)
still hardening at 11.2 %, and again died on the step floor. It paid for that
with **95 of 305 steps on the 100× rung**.

I originally read this as "a better solver does not walk through it, therefore
the wall is the discretization". **That inference was too strong.** A relaxed
*tolerance* is not a stronger path algorithm — it is the same algorithm
permitted to accept worse answers, and it moved the number in the same
direction and the same order of magnitude as TIMs' budget change. Taken with
§4.1 it reads as one more knob that shifts the wall a little, not as evidence
that the wall is intrinsic.

### 4.4 What the volumetric treatment is worth, by order

| order | locked | relieved | gain | ends on a plateau? |
|---|---|---|---|---|
| linear hex | 0.40 (TIMs) | **0.9977** (B-bar) | ×2.5 | **yes** |
| quadratic hex | **0.5894** (`std`) | **0.7715** (`uri`) | ×1.31 | no |
| quadratic tet | 0.2828 / 0.3075 (TIMs) | ≈ 0.49–0.58 (mean dilatation, TIMs) | ×1.7 | no |

Volumetric relief moves the number at every order. **It only produces a plateau
at linear order.** That much is observation. The former reading of this table —
that the pattern *is* the constraint ratio at work — is withdrawn; see §5.3.
Note also that the two quadratic rows are made of numbers that §4.1 says are
controller-dependent, so the ×1.31 and ×1.7 are differences between two
non-converged path failures and should not be quoted as element properties.

### 4.5 Mesh refinement

| element | 2 el/B | 4 el/B | 8 el/B | 16 el/B |
|---|---|---|---|---|
| H20 `uri` (here) | 0.7715, wall at s/B = 0.0112 | 0.6894, wall at s/B = 0.0103 | — | — |
| H8 B-bar (here / COLLAPSE.md) | — | **0.9977** / 1.0020 | 0.9517 | 0.9260 |

The H8 series plateaus at every resolution — a proper convergence sequence
against an exact answer. The H20 pair does not, and its wall arrives *earlier*
in settlement as the mesh refines (0.0112 → 0.0103). COLLAPSE.md reads that
"monotone loss of reach" as the Rudnicki–Rice localization signature under
non-associated perfect plasticity; that remains a plausible reading, but with
§4.1 in hand it is equally consistent with the path simply being harder to
follow on a bigger system, and this note cannot separate the two.

### 4.6 The falsification control fired

The associated leg runs *past* the exact answer — 1.1422 — the documented
behaviour (dilating at ψ = φ must push the surrounding soil aside; a bounded
mesh resists, and the leg hardens through its own oracle). Its
`plateau = yes` flag is **not** a collapse load: the test looks at the last
10 % of a curve that only reached s/B = 0.036, and it flattens 14 % *above* a
bound that is exact from both sides. The control fired, confirming the
non-associated leg is the one carrying the measurement.

## 5. What survives

### 5.1 The empirical observable

This is the part that does not depend on any mechanism, and it is the
deliverable:

> **Six quadratic legs, across two campaigns, two harnesses, two element
> families (hex and tet), two bases (Lagrange and Bernstein), relieved and
> locked, ALL wall while still hardening at 10–27 % of their initial tangent.
> One linear relieved element plateaus and reaches its target — twice,
> independently (0.9977 here, 0.9934 TIMs).**

The legs: H20 `uri` at 2 el/B and 4 el/B, H20 `std` (this note);
`TenNodeTetrahedron`, `BezierTet10` std, `BezierTet10` mean-dilatation (TIMs).
The split is perfectly clean and it falls on **order**.

Two eliminations come with it, and both are basis- and geometry-level rather
than mechanism-level, so neither is touched by the correction:

- **(c) the Bézier/Bernstein basis is exonerated.** It was already nearly
  excluded by the 8 % gap between the Bernstein and Lagrange tets; this note
  reproduces the entire behaviour inside a **pure Lagrange hex**, with no
  Bézier anywhere.
- **(a) simplex geometry is not the dominant term.** At matched order the hex
  beats the tet by roughly **2×** (0.5894 vs 0.2828 locked; ≈0.77 vs ≈0.52
  relieved). Real, and much smaller than the ~3.5× gap being explained.

### 5.2 What was NOT measured

**No collapse load was measured for any quadratic element.** Every quadratic
leg in both campaigns terminated on a path failure — my step floor, their
subdivision budget — with the curve still climbing. An independent leg has now
been shown to move **+17 % in capacity and +59 % in reach** under a change to a
single solver integer with the discretization held fixed (§4.1).

So the honest statement of the quadratic result is:

> Quadratic elements **lose the Newton path before their mechanism forms**,
> somewhere between **0.59 and 0.77** of exact on this problem. Those bounds
> describe where the solver stopped, not where the element would have stopped.

Consequently the first version's "the quadratic hex tops out at 0.69–0.78" is
wrong as phrased, and any downstream use of 0.77 as an element ceiling should
be corrected.

### 5.3 The `r` / constraint-ratio hypothesis: REJECTED

I presented the constraint ratio as the mechanism. It is not, and the
disconfirming case was already in the campaign's own data.

**The counter-example.** TIMs' `BezierTet10 -bbar` leg counts at **r ≈ 4.4**
(≈24·V DOF against ≈5.5·V elements, one volumetric constraint each) — *above*
the linear hex's 3.0, the ratio that cures the problem — and it behaves like my
H20 `std` at **r = 0.44**: low, and no plateau. I verified their premise at
source rather than taking the count on trust:
`BezierTet10::computeVolumeAveragedDerivatives` (`BezierTet10.cpp:1398-1431`)
forms `B̄_a,i = (1/V_e)∫∂N_a/∂x_i`, a single element-wide constant on the full
4-point rule, and `computeBBarMatrix` (`:1363-1395`) is textbook mean
dilatation — one volumetric constraint per element, deviatoric part left at
full quadratic resolution. The count is right, and it points the wrong way.

Three points ordered correctly and a fourth that breaks the ordering is not a
mechanism; it is a correlation fitted to an ordering with some other cause.

**Why `r` was the wrong tool here — the reason worth keeping** (TIMs' framing,
and the most useful thing to come out of this correction): `r` counts degrees
of freedom against constraints, and it was built for the **elastic
near-incompressible limit**, where the disease genuinely is *too few DOFs per
constraint* — the element is over-constrained and locks. A **collapse
mechanism** is a different question. What matters there is whether the
**available isochoric velocity fields can represent the slip pattern** — the
*span* of the admissible set in the neighbourhood of the mechanism, not a
*count* of it. A space can have plenty of DOFs per constraint and still not
contain the field the mechanism needs. Those two questions coincide often
enough that `r` is a decent heuristic, and evidently they do not coincide here.

That reframing is what makes §5.5's instrumentation ordering follow.

### 5.4 The counterweight: this is not "the driver is inadequate"

The obvious deflation of §5.2 is that the harness simply cannot follow a
perfectly plastic collapse and the quadratic numbers say nothing about
elements. That is ruled out by the baseline: **the linear hex reaches a genuine
plateau in this same harness** — same script, same mesh, same ladder, same
tolerances, same material — running the full s/B = 0.25 push with a flat
tangent. The path algorithm is demonstrably capable of taking *an* element
through *this* collapse.

So the question is not "can the driver do it" but **"why do quadratic elements
lose the path"**. That is a narrower and better-posed question than the one
this note started with, and it is the right question to point at H27.

### 5.5 H27 + selective — and it must be instrumented

H27 + selective is the ratio-optimal quadratic brick and remains the single
most valuable element to build off this note, but §5.3 changes *why*: not to
test the `r` story (which is rejected) but because it is the quadratic element
most likely to have the richer isochoric span, and because it is the cleanest
way to ask whether *anything* quadratic can complete this collapse.

It is only worth building if it is instrumented to separate **"cannot form the
mechanism"** from **"the path failed"**, since those produce identical CSVs.
In priority order:

1. **The mobilisation field at termination** — first, because §5.3 says the
   question is span, not count. Is a *localised slip surface* forming, or is
   plasticity smeared across the domain? A partly-formed band at the wall says
   the kinematics are available and the path was lost; diffuse plasticity says
   the element cannot represent the mechanism at all. The harness already
   computes a mobilisation field (`r3_prandtl_tet10.py` writes one); it needs
   to be dumped at the wall, not only on success.
2. **Arc-length or displacement control taken well past the wall** — a path
   method that does not depend on the tangent staying invertible. If the curve
   continues and flattens beyond where load control died, the wall was the path.
3. **Tangent conditioning / least eigenvalue as the wall approaches** — does
   the operator degenerate smoothly toward a mechanism (real limit point), or
   abruptly (path/ill-conditioning)?
4. **The decisive one: does the wall MOVE when the path algorithm changes but
   the mesh does not?** This is the generalisation of TIMs' budget experiment
   and my §4.2 pair, and it should be run as a designed sweep rather than
   noticed after the fact. A wall that holds still across genuinely different
   path algorithms is an element property; one that moves is not.

**Every future leg on this problem should report the budget-sensitivity pair
by default** — two runs differing only in the controller allowance. It costs
one extra run and it is the difference between a measurement and an artifact.

### 5.6 Recommendation on an F-bar-Patch tet lane

The brief asked for this conditional on verdict (a). Simplex geometry is *not*
the dominant term (§5.1), so the lane as originally scoped is not supported —
but the reasoning has changed since the first version and the conclusion is now
held more weakly:

- **Do not open F-bar-Patch on QUADRATIC tets.** Previously I argued this from
  the `r` account, which is withdrawn. The surviving argument is the empirical
  one: mean dilatation on a quadratic tet has already been tried by TIMs, at a
  constraint ratio well above the curing value, and it does not plateau —
  it walls like everything else quadratic. A different quadratic-simplex macro
  might behave differently, but nothing measured suggests it, and there is no
  longer a mechanism on the table that predicts it would.
- **Open it on LINEAR tets, if at all.** Linear + volumetric relief is the only
  combination that has ever completed this collapse — twice, independently —
  and the fork already owns the 2D form (`LadrunoCSTPair`, dSNPO §15.1.9). An
  F-bar-Patch Tet4 aims at the one cell of the 2×2 that is measured to succeed,
  and it is the lane that would let unstructured meshing reach what the
  structured hex mesh already reaches.
- **Neither should be opened before the H27 instrumentation in §5.5 has said
  whether quadratic elements are failing on kinematics or on the path.** That
  answer changes which macro is worth building, and it is much cheaper than
  building either.

### 5.7 Scope

One problem (Prandtl–Reissner), one flow rule (non-associated, ψ = 0), one
material (`DruckerPrager`, φ_ps = 27.47°, σ_y = 0.2 kPa regularizer), one mesh
family (structured graded strip, plane strain), two resolutions on the H20 and
one on the H8. Nothing here carries to associated flow, to severe
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

Any leg quoted as a capacity must be run at **two controller allowances**
(§5.5) and both reported.
