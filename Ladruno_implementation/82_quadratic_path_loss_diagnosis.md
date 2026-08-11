# Note 82 — Why quadratic elements lose the Newton path: a diagnosis, a termination-mode table, and the H27 blockage

**Kind:** MEASUREMENT / DIAGNOSIS. Ships no feature, changes no element, and
**deliberately reports no collapse load.**
**Predecessor:** [[81_quadratic_hex_limit_load_measurement]] — read its §0 first;
its retractions are load-bearing here and this note inherits them.
**Harness:** `Ladruno_files/testbed/hypo_bearing/quad_path_diag.py`, built on
`h20_prandtl.py` (mesh, consistent Q8 surcharge, all four of its controls,
imported rather than copied) and classified by the merged R3 gate's three-clause
capacity rule (`tests/test_r3_prandtl_collapse_gate.py`, PR #722).
**Engine:** `ladrunoBuild() = 89353260616ae73ac3bfacf37261b1027defcc70`
**Problem:** Prandtl–Reissner, `q_u = q0·N_q` exact, φ_txc = 20°, q₀ = 10 kPa,
γ = 0, ν = 0.45, Chen & Han plane-strain cone match, non-associated (ρ̄ = 0).
Unchanged from note 81 and from the gate.

> **ATTRIBUTION AND HOLD.** Results credited to the **TIMs element-technology
> campaign** are theirs, cited with credit. TIMs material is internal-only,
> under the same hold as the unfiled UW report: do not carry their numbers or
> phrasing into anything public or upstream-facing.

---

## 0. The scope correction that has to come first: THE FORK OWNS NO 27-NODE HEX

This work was scoped as "H27 + selective integration". **That element does not
exist in this repository and could not be run.** Verified at source, not
assumed:

- `LadrunoBrick20` is **20-node serendipity** (Abaqus C3D20 ordering; see the
  `LadrunoHex20Shape.h` ORDERING MEMO). Its two formulations are `std` (full
  3×3×3) and `uri` (uniform reduced 2×2×2).
- Every "27" in `SRC/element/brick/` is a **Gauss-point count**, not a node
  count — e.g. `Twenty_Node_Brick.h:158-161` declares `shgu[4][20][27]`,
  `wu[27]`, `dvolu[27]`: 20 nodes, 27 integration points.
- A class tag **is reserved** — `SRC/classTags.h:717`,
  `#define ELE_TAG_TwentySevenNodeBrick 45` — but **there is no implementing
  class anywhere in the tree**. `grep -rn "class TwentySevenNodeBrick" SRC/`
  returns nothing. The only other references are four recorder type-maps
  (`PVDRecorder.cpp:1857`, `VTK_Recorder.cpp:1052`,
  `VTKHDF_Recorder.cpp:3154`, `GmshRecorder.cpp:1314`) which map the tag to a
  placeholder `POLY_VERTEX`, i.e. "we do not know how to draw this", not to
  VTK's real `TRIQUADRATIC_HEXAHEDRON`. There is no Tcl or Python command, no
  CMake entry, and no object-broker case.

Building one is **element development, not a measurement**, and it is out of
scope for this note. Nothing below is a substitute for it and nothing below is
presented as if an H27 ran. The recommendation on whether to open that lane —
now informed by what was actually measured — is §7.

---

## 1. The question, restated

Note 81 §5.4 narrowed it correctly and this note takes it verbatim:

> The linear relieved element reaches a genuine plateau **in this same harness**
> — same script, same mesh, same ladder, same tolerances, same material. So the
> question is not "can the driver do it" but **"why do quadratic elements lose
> the path"**.

Two families of answer, and they produce identical CSVs:

- **KINEMATIC (span).** The element's admissible isochoric velocity fields
  cannot represent the collapse mechanism, so there is no nearby equilibrium
  branch to converge onto and Newton is being asked to find something that is
  not there. Note 81 §5.3 reframed the question this way after the constraint
  ratio `r` was rejected: what matters is the **span** of the admissible set
  near the mechanism, not a **count** of DOFs against constraints.
- **NUMERICAL (path).** The mechanism is representable, and the path algorithm
  simply fails to follow it.

The whole design below exists to separate those two.

## 2. Instrumentation, and the controls that make it quotable

Every problem parameter, the mesh, the consistent Q8 surcharge and note 81's
controls 0/1/2/3/4 are **imported from `h20_prandtl.py`, not copied**, so the
expensive load algebra cannot fork. Two new instruments were added, and each one
was calibrated before it was used.

### 2.1 The mobilisation field (ported from the reference driver)

`||dev σ|| / (√(2/3)·σ_y − ρ·I₁)`, exactly 1.0 at yield — the quantity
`r3_prandtl_tet10.py` writes for its mechanism figure, ported and **widened**:
the reference samples material point 1, this samples **every Gauss point** and
reports the element mean and max separately.

**Control M0 — and it fired on the first run, which is why it is worth having.**
At the surcharge step the state is the exact 1-D state, so the mobilisation has
a closed form. Comparing against `h20_prandtl.py`'s printed `m0` gave a
**7.67e-3** discrepancy. That is not a bug in the field: `m0` as the harness
prints it *drops the `√(2/3)·σ_y` apex term*, which is a fine screen for "is the
leg void?" and is 2.9 % wrong as a field reference. Against the exact form

    m0* = q0(1−K0)√6/3 / (√(2/3)·σ_y + ρ·q0(1+2K0)) = 0.260306509

the field agrees to **2.6e-14** on every Gauss point of every element, and
**0 of 200 elements** are above the yield threshold before the push.

### 2.2 The isochoric span ratio χ — an exact volume balance

This is the instrument that carries the mechanism verdict, and it is a
**closed identity on this domain**, not a heuristic. At ψ = 0 plastic flow is
exactly isochoric, so the total volume change is

    ΔV = ∮ u·n dA

and every boundary of this box contributes zero except the top: the base is
fixed, the sides are rollers with `u_x = 0`, and `u_y = 0` on **every** node
(that is what makes it plane strain). So `ΔV = ∫_top u_z dA = Σ_n w_n·Δu_z,n`,
computed with the **consistent** surface weights `w` — the same `∫N_a dΓ` the
surcharge is built from, which is exactly the right operator for integrating an
FE field over that surface. The rigid footing pushes its footprint down by `ds`,
so

    **χ = 1 + ΔV / (A_foot · ds)**,  A_foot = B·THICK,

with **χ = 1** meaning every bit of the footing increment came back up beside it
(a fully isochoric mechanism) and **χ = 0** meaning it was swallowed as
volumetric compaction.

Two things had to be got right and the first cut got both wrong:

- `A_foot` is the **geometric** constant. Summing `w` over the footing *node
  set* gives 1.1667 m² against a true footprint of 1.0000 m², because the nodes
  on the footing edge collect `∫N dΓ` from the element outside it too — a 17 %
  error straight into the headline.
- A plain `max(u_z)` beside the footing is **not** a usable proxy: measured, its
  argmax sits at **x = 30 m, the far boundary**, where the box breathes
  elastically. A surface integral cannot be fooled that way.

**Control CHI — the two-sided calibration.** χ's meaning is a claim about the
measure, and it is checkable with no plasticity at all: a linear elastic body at
ν → ½ is exactly incompressible and must give χ → 1; at ν = 0 it must give χ
well below 1. Measured on the H20 `uri` mesh at h₀ = 1.0:

| ν | 0.0 | 0.25 | 0.45 | 0.49 | 0.4999 |
|---|---|---|---|---|---|
| χ | −3.902 | −3.473 | **−0.883** | +0.534 | **+0.995** |

Monotone, and both endpoints hit. **The ν = 0.45 entry is the operating point
and it is the reason χ is interpretable**: at ν = 0.45 an elastic body is
already fairly incompressible, so "no mechanism" is **not** χ = 0, it is
χ = −0.883 measured on that exact mesh and element. Reporting χ against 0 would
read a purely elastic increment as a half-formed mechanism. Every leg therefore
runs its own elastic probe first and reports the normalised

    **M = (χ − χ_el) / (1 − χ_el)** — 0 = the increment is elastic, 1 = fully isochoric.

**What would have to be true for M to be high while no mechanism is forming?**
M tests **span**, not **localisation**: a completely *diffuse* isochoric plastic
flow also returns M ≈ 1. That is why M and the localisation metrics are reported
separately and neither substitutes for the other — M answers "can this element
deliver an isochoric mechanism at all", the localisation metrics answer "is what
it delivers a slip surface". The converse combination — a large yielded set with
M ≈ 0 — is the unambiguous signature of **locking**.

### 2.3 The localisation metrics, and why they are volume-weighted

All fractions are **volume** weighted. This mesh is graded at 1.35×, so element
volumes span **32.5×** on the h₀ = 1.0 grid; a *count* fraction is dominated by
far-field elements that are an order of magnitude larger and always elastic, and
would read "localised" for a completely diffuse field. Counts are printed beside
the volumes as the control that shows the two disagree.

Connectivity uses **8-connectivity** (face *or* edge neighbours) on the
structured (i,k) element grid. This matters: a Prandtl slip surface runs
*diagonally*, so a one-element-thick diagonal band is **not 4-connected at all**.
Measured on a synthetic one-element-thick band on this mesh: **5 components
under 4-connectivity, 1 under 8**. With 8-connectivity, "the yielded set is not
connected" becomes a strong negative result instead of an artefact of the
stencil.

**The metrics were validated by discrimination, not by inspection.** A synthetic
localised band and a synthetic diffuse blob **of identical plastic volume
fraction** (0.23 % vs 0.25 %) were fed through the same code:

| synthetic field | equivalent thickness | cells/column | bounding-box fill |
|---|---|---|---|
| localised band | **1.30 el** | **1.25** | **32.6 %** |
| diffuse blob | 2.00 el | 2.00 | 66.7 % |

A metric that gave the same answer for both would be useless for the question
being asked; these separate.

## 3. Diagnostic 2, part one: arc-length is not the experiment, and here is why

The brief asked for "arc-length or displacement control taken well beyond where
the current controller dies". **The leg is already displacement-controlled**, and
this is worth stating precisely because it changes what the remaining
experiments can be.

Reading the harness: the push is a `Plain` pattern carrying a `Linear` time
series and `ops.sp(node, 3, 1.0)` on every footing node, driven by
`LoadControl(-ds)` on pseudo-time. The `Transformation` handler applies the
constraint as `u_z = λ · 1.0`, and the reported settlement is `s = uz0 − λ`. So
**the load parameter λ IS the footing settlement**, prescribed, and the
increment is a prescribed displacement increment of size `ds`.

Consequences:

- **λ is monotone by construction. There is no limit point in this
  parameterisation for arc-length to pass.** Arc-length exists to follow a
  branch where the load factor turns over; here the "load factor" is a
  settlement that only increases. `LadrunoArcLength`'s constraint
  (`dUᵀdU + α²dλ² = ℓ²`) would, on this model, be a norm-based step-size
  adaptation on a path that is already single-valued — not a different
  branch-following capability.
- The failure is therefore **not** "the path parameterisation lost the branch".
  It is **Newton failing to converge inside a prescribed-displacement increment**
  — and, at termination, an increment of order 10⁻⁷ m. That is a different
  disease from the one arc-length treats.

So running arc-length here and reporting the result would have produced a number
with no bearing on the question. **It was not run, and that is a deliberate
negative result, not an omission.** What replaces it, as genuinely different
path experiments on an identical mesh:

1. **the controller's ALLOWANCE** (§5) — the direct analogue of the TIMs
   subdivision-budget sweep, and the one note 81 §4.2 thought it had run;
2. **a genuinely different set of nonlinear algorithms** within the increment
   (`--ladder wide`: adds `ModifiedNewton`, `Broyden`, `BFGS` to the ladder).
   This is *not* the tolerance relaxation of note 81 §4.3, which that note
   correctly stopped treating as a stronger path algorithm.

**The path method this note did NOT try, and which is the genuinely decisive
one:** a *dynamic* path — quasi-static transient / dynamic relaxation, which
does not need the tangent to stay invertible at all. The fork already owns the
machinery (`LadrunoFictitiousMass`, `LadrunoDynamicRelaxation`, ADR 21). It is
the single most valuable remaining experiment and it is named in §7.

## 4. What note 81 §4.2's step-ladder pair actually tested — a correction

Note 81 §4.2 ran the H20 `uri` leg on "two step ladders a factor of 10 apart in
base and floor" (2e-4/2e-6 and 2e-5/2e-7), found the wall in the same place
(0.7706 vs 0.7715), and read that as weak evidence for an element property.

**That pair could not have moved the wall, and the reason is arithmetic.** The
controller's actual allowance is not the base step and not the floor — it is
**how many halvings it is permitted before it hits the floor**, i.e.
`log2(ds_base / ds_min)`. Both of note 81's ladders are 100:1, so both allow
**6.64 halvings**. The two runs had *identical* freedom; scaling base and floor
together is a **rescaling, not an allowance change**. The agreement to 0.1 % is
what a rescaling should produce and carries no information about whether the
wall is movable.

This is the fourth instance of the rule note 81 §5.3 credits to TIMs — *what
would have to be true for this check to pass while the thing it checks is
wrong?* — and it was hiding inside the correction that established the rule.

The allowance experiment in §5 therefore moves **the floor alone**, which
changes the halving count, and separately raises the budget.

## 5. The termination-mode table — the deliverable

Build `89353260616ae73ac3bfacf37261b1027defcc70`. Exact `q_u = 138.907 kPa`.
Non-associated (ρ̄ = 0) throughout. Classification is the merged R3 gate's.

| leg | h0 | DOF | q/q_exact | tail % | **MODE** | ds_end (mm) | × floor | subdiv | halvings | plateau | free | **CAPACITY** | wall s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `h8bbar_h5` | 0.50 | 2592 | **0.9977** | 0.004 | **BUDGET** | 5.0e-01 | 125 | 81/80 | 6.64 | yes | yes | **YES** | 1246 |
| `h8bbar_h5_at0103` | 0.50 | 2592 | 0.7698 | 11.858 | TARGET | 6.0e-01 | 250 | 4/80 | 6.64 | NO | yes | NO | 92 |
| `h20uri_h5` | 0.50 | 8814 | 0.6846 | 10.654 | **FLOOR** | 3.1e-04 | **1.6** | 72/200 | 6.64 | NO | NO | NO | 1234 |
| `h20uri_h5_loose` (run A) | 0.50 | 8814 | 0.7306 | 10.346 | **FLOOR** | 5.0e-03 | 78 | 118/800 | **13.29** | NO | NO | NO | 2652 |
| `h20uri_h5_loosefloor` (run B, **same config as A**) | 0.50 | 8814 | 0.6951 | 11.771 | **FLOOR** | 1.3e-03 | 20 | 90/800 | **13.29** | NO | NO | NO | 1731 |
| `h20std_h5` | 0.50 | 8814 | 0.4587 | 31.305 | **FLOOR** | 2.5e-03 | 12.5 | 65/200 | 6.64 | NO | NO | NO | 1056 |
| `h8bbar_h1.0_at0112` | 1.00 | 1386 | 0.8845 | 10.632 | TARGET | 1.0e+00 | 500 | 0/80 | 6.64 | NO | yes | NO | 23 |
| `h8bbar_h1.0_traj` | 1.00 | 1386 | **1.0849** | 0.002 | **BUDGET** | 1.0e+00 | 250 | 81/80 | 6.64 | yes | yes | **YES** | 845 |
| `h20uri_h1.0_traj` | 1.00 | 4659 | 0.7713 | 11.642 | **FLOOR** | 3.1e-03 | 7.8 | 41/200 | 6.64 | NO | NO | NO | 891 |
| `h20uri_h1.0_wide` | 1.00 | 4659 | 0.7718 | 11.623 | **FLOOR** | 3.1e-03 | 15.6 | 38/200 | 6.64 | NO | NO | NO | 688 |
| `h20std_h1.0` | 1.00 | 4659 | 0.5894 | 27.238 | **FLOOR** | 2.5e-03 | 12.5 | 47/200 | 6.64 | NO | NO | NO | 681 |

`_at0103` / `_at0112` are **matched-settlement mechanism controls** (§6), deliberately
stopped early; their `TARGET` is a truncation, not a capacity, and the table says so.

**Only the two LINEAR full legs are capacities.** Every quadratic leg terminates
`FLOOR` — seizure — and none may be quoted as an element ceiling. **Runs A and B
of `h20uri_h5_loose*` are the SAME configuration and differ by 5.1 %** — see
§7.1.1 before reading any quadratic figure to more than one digit.

**This independently reproduces note 81 and the merged gate**, on a different
build and through a different driver:

| leg | here | note 81 / gate #722 | Δ |
|---|---|---|---|
| H8 bbar, h0 = 0.5 | 0.9977 | 0.9977 | 0.00 % |
| H8 bbar, h0 = 1.0 | 1.0849 | 1.0849 | 0.00 % |
| H20 `uri`, h0 = 1.0 | 0.7713 | 0.7715 | −0.03 % |
| H20 `std`, h0 = 1.0 | 0.5894 | 0.5894 | 0.00 % |
| H20 `uri`, h0 = 0.5 | 0.6846 | 0.6894 | −0.70 % |

Note in passing: at h0 = 0.5 the H20 legs carry **8814 DOF against the linear
hex's 2592** — 3.4× the unknowns for a worse answer on this problem.

## 6. Diagnostic 1 — mechanism formation. THE SPAN HYPOTHESIS IS REFUTED

| leg | s/B end | χ_el | χ | **M** | yielded vol % | thickness (el) | comps | reaches |x| |
|---|---|---|---|---|---|---|---|---|
| `h8bbar_h5` (completes) | 0.0744 | −0.889 | 1.000 | **1.000** | 25.15 | 22.96 | 1 | 11.54 m |
| `h8bbar_h5_at0103` (matched) | 0.0103 | −0.889 | 0.592 | **0.784** | 12.50 | 20.12 | 1 | 6.51 m |
| `h20uri_h5` (WALL) | 0.0101 | −0.860 | 0.593 | **0.781** | 9.63 | 20.76 | 1 | 4.84 m |
| `h20uri_h5_loose` (WALL) | 0.0121 | −0.860 | 0.630 | 0.801 | 11.54 | 18.57 | 1 | 6.51 m |
| `h20std_h5` (WALL) | 0.0030 | −0.871 | 0.250 | 0.599 | 4.47 | 13.02 | 1 | 3.56 m |
| `h8bbar_h1.0_traj` (completes) | 0.1457 | −0.937 | 1.000 | **1.000** | 24.91 | 11.46 | 1 | 11.35 m |
| `h8bbar_h1.0_at0112` (matched) | 0.0112 | −0.937 | 0.619 | **0.803** | 15.41 | 13.22 | 1 | 5.95 m |
| `h20uri_h1.0_traj` (WALL) | 0.0112 | −0.883 | 0.600 | **0.788** | 11.89 | 10.21 | 1 | 5.95 m |

### 6.1 The matched-settlement control, which is the whole answer

Comparing a leg that died at s/B = 0.011 with one that ran to s/B = 0.074 tells
you only that one ran further. So the linear control was **re-run and stopped at
exactly the settlement where the quadratic leg dies**:

| at s/B = 0.0112 (h0 = 1.0) | H8 `bbar` | H20 `uri` |
|---|---|---|
| q/q_exact | 0.8845 | 0.7713 |
| **mechanism index M** | **0.803** | **0.788** |
| tail dq/ds, % of initial | **10.63 %** | **11.64 %** |
| yielded volume | 15.4 % | 11.9 % |

| at s/B ≈ 0.0103 (h0 = 0.5) | H8 `bbar` | H20 `uri` |
|---|---|---|
| q/q_exact | 0.7698 | 0.6846 |
| **mechanism index M** | **0.784** | **0.781** |
| tail dq/ds, % of initial | **11.86 %** | **10.65 %** |
| yielded volume | 12.5 % | 9.6 % |

**M agrees to 0.4 % at h0 = 0.5 and 2 % at h0 = 1.0.** At the moment the
quadratic element loses the path, its collapse mechanism is *at the same stage
of formation* as the linear element's — and the linear element continues from
that identical state to a complete, perfectly isochoric mechanism (M = 1.000)
and the exact Prandtl answer.

> **Verdict: quadratic elements are not failing for want of kinematics.** The
> "the available isochoric velocity fields cannot represent the slip pattern"
> reading of note 81 §5.3 is **refuted** at the point of failure. The fields are
> there, the mechanism is forming, and it is forming at the same rate as in the
> element that succeeds.

### 6.2 A correction to note 81 §5.1's headline observable

Note 81's empirical observable was:

> Six quadratic legs ... ALL wall while still hardening at **10–27 %** of their
> initial tangent. ... The split is perfectly clean and it falls on **order**.

**It does not fall on order.** The linear B-bar element, at s/B = 0.0103, is
hardening at **11.86 %** of its initial tangent — squarely inside that band —
and at s/B = 0.0112 at **10.63 %**. It then goes on to 0.9977 with a 0.004 %
tail at s/B = 0.0744, seven times further along.

So "still hardening at 10–27 % of initial" is not a quadratic property. It is
simply **what this problem's load–settlement curve looks like at s/B ≈ 0.01**,
for every element tested. What distinguishes the quadratic legs is only that
they *stop* there. The measurements in note 81 stand; the inference that the
observable "falls on order" does not.

### 6.3 Localised slip surface, or smeared plasticity? Neither — and the question does not discriminate

The brief asked whether a localised slip surface is forming or plasticity is
smeared. Measured, on this mesh family, **no leg forms a slip surface, and the
element that COMPLETES the collapse has the most smeared plastic zone of all**:
25.15 % of the domain by volume, ~23 elements thick, one connected component,
reaching |x| = 11.5 m. Its mobilisation field at the plateau is *saturated* —
every element in the core at mob = 1.000:

```
  [h8bbar_h5]  s/B = 0.0744, q/q_exact = 0.9977, mode BUDGET, M = +1.000
  z= -0.25 |################|      [h20uri_h5]  s/B = 0.0101, WALL, M = +0.781
  z= -0.75 |################|      z= -0.25 |   .+######+.   |
  z= -1.25 |################|      z= -0.75 |...:########:...|
  z= -1.75 |################|      z= -1.25 |..::########::..|
  z= -2.25 |################|      z= -1.75 |.::+#+####+#+::.|
  z= -2.75 |################|      z= -2.25 |.::#++####+##+:.|
  z= -3.25 |################|      z= -2.75 |::+##########+::|
  z= -3.83 |################|      z= -3.25 |:++##########++:|
  z= -4.59 |################|      z= -3.83 |:+##+######+##+:|
  z= -5.59 |################|      z= -4.59 |:+##++#++#++##+:|
  z= -6.90 |################|      z= -6.90 |:+#+#######++#+:|
          |------FFFF------|              |------FFFF------|
   '#' mob>=0.99  '+' >=0.90  ':' >=0.70  '.' >=0.40   F = footing, x = ±6.5 m
```

That is exactly what the merged gate's docstring predicts for displacement-based
elements at a limit load — "too few DOF to form the true velocity
discontinuity, so the cheapest available mechanism is stiffer and stronger" —
and it means the localisation framing is the wrong instrument here. The
quantity that carries information is **M** and the **extent** of the plastic
zone, not its thinness.

**A hypothesis that the field maps raised and the measurement DID NOT SUPPORT.**
The H20 map above is visibly *mottled* — adjacent elements alternating between
`#` and `+` — where the H8 map is smooth, which reads as a checkerboard
(pressure-mode) instability, and H20 `uri` does carry 6 spurious modes per
element with no hourglass control. It was quantified
(`quad_path_osc.py`: RMS of `mob_e − mean(4 neighbours)` over the plastic core)
and **it does not hold up**:

| | H8 `bbar` (matched) | H20 `uri` (wall) |
|---|---|---|
| h0 = 1.00, neighbour residual | **0.00 %** | 3.96 % |
| h0 = 0.50, neighbour residual | **3.99 %** | **2.08 %** |

At the coarse mesh it looks decisive; **at the finer mesh it reverses**, with
the *linear* element reading as the more oscillatory one. Both readings have the
same cause and it is not a checkerboard: the residual cannot separate
high-frequency noise from a high-**curvature** plastic front, and it is
identically zero on a saturated core because that field is *constant*, not
because it is smooth. The visible mottling is a **thresholding artefact** — a
core sitting at 0.985–0.995 straddles the 0.99 render threshold on a ~2 %
variation. Reported as **inconclusive**; §8 names the test that would settle it.

## 7. Diagnostics 2 and 3 — the wall moves, and the operator is healthy until it isn't

### 7.1 The allowance experiment: THE WALL MOVES

Two legs identical in element, mesh, loads, material, tolerances, algorithm
ladder and settlement target. **Only the controller's allowance differs.**

| H20 `uri`, h0 = 0.5 | floor 2e-7 m (6.64 halvings), budget 200 | floor **2e-9 m** (**13.29** halvings), budget **800** |
|---|---|---|
| q/q_exact | 0.6846 | **0.7306** |
| reach s/B | 0.01009 | **0.01205** |
| tail % of initial | 10.65 | 10.35 |
| mode | FLOOR | FLOOR |
| M | 0.781 | 0.801 |

That single pair reads **+6.7 % capacity, +19.4 % reach** for a 100× lower step
floor and nothing else — the TIMs subdivision-budget result (+17 % / +59 %)
reproduced inside the fork's harness on a quadratic **hex**, with the allowance
varied along the axis §4 shows note 81's own control failed to vary.

### 7.1.1 …but a REPEAT of the loose leg does not reproduce it, and that changes the claim

The loose-floor leg was run a **second time, identical in every argument**
(`--elem h20uri --h0 0.5 --budget 800 --floor 2e-9 --sfrac 0.15`). It did not
reproduce the first:

| H20 `uri`, h0 = 0.5 | reference allowance | loose, run A | loose, run B |
|---|---|---|---|
| halvings / budget | 6.64 / 200 | 13.29 / 800 | 13.29 / 800 |
| q/q_exact | 0.6846 | **0.7306** | **0.6951** |
| reach s/B | 0.01009 | **0.01205** | **0.01051** |
| tail % of initial | 10.65 | 10.35 | 11.77 | 
| M | 0.781 | 0.801 | 0.802 |
| mode | FLOOR | FLOOR | FLOOR |

**Two bit-identical configurations differ by 5.1 % in capacity and 14.7 % in
reach.** The mechanism is not mysterious and the merged gate already names it in
its own band justification: a threaded PARDISO factorisation is not
FP-deterministic, and **an adaptive controller amplifies that without limit —
one step converging differently re-sequences every step after it**. A leg that
ends by *seizing* is far more exposed to this than one that ends on a plateau,
because its endpoint is precisely a convergence failure.

**What survives, and what does not.**

- **Does not survive:** the magnitude. "+6.7 % capacity, +19.4 % reach" is one
  draw from a distribution whose spread is of the same order as the effect. The
  allowance effect measured against the two loose runs is **+1.5 % to +6.7 %**
  (mean +4.1 %) in capacity and **+4.2 % to +19.4 %** in reach. With n = 1 at
  the reference allowance and n = 2 at the loose one, **this design does not
  resolve the size of the allowance effect** and no single number from it should
  be quoted.
- **Survives, and is strengthened:** the direction is consistent — **both** loose
  runs went further and higher than the reference run (2/2). And the headline
  conclusion is now over-determined, because it no longer even needs the
  allowance experiment:

> **No quadratic figure on this problem is a ceiling — and none of them is
> stable to three digits.** If two runs of the *same* configuration on the *same
> binary* differ by 5 % in capacity and 15 % in reach, then quoting 0.6846,
> 0.6951, 0.7306, 0.7713 or 0.5894 as an element property is indefensible
> regardless of what the allowance does.

That also re-frames note 81 §4.2's surviving ladder pair ("0.7706 vs 0.7715,
agreement to 0.1 %"): at h0 = 1.0 that was either a fortunate draw or a genuinely
more stable leg, but it was *never* evidence about the element, and §4 already
showed the two ladders had identical allowances anyway.

**A third qualification that does survive both runs.** The wall **recedes but
does not disappear**: every loose leg still ended on `FLOOR`, still hardening at
10–12 %, still with no plateau. Nothing here is a sequence converging to
anything, and it must not be extrapolated.

**What would settle the magnitude:** 4–6 repeats at each of 3 floors, reported
as distributions rather than points. That is ~10 leg-hours and is the honest
price of a number here; note 82 does not pay it and therefore does not quote one.

### 7.2 A different algorithm set does NOT move the wall

Same mesh, same allowance, a genuinely different ladder (adding `ModifiedNewton`
and full `Newton` rungs at the same tolerances — *not* note 81 §4.3's tolerance
relaxation):

| H20 `uri`, h0 = 1.0 | reference ladder | wide ladder |
|---|---|---|
| q/q_exact | 0.7713 | 0.7718 |
| reach s/B | 0.01122 | 0.01124 |
| tail % | 11.64 | 11.62 |

**+0.06 %.** Nothing. Taken with §7.1: the wall is set by **how small a step the
controller is allowed to take**, not by how good the algorithm is at taking it.

### 7.3 Tangent conditioning: the operator is well conditioned until an ABRUPT event

The tangent is unsymmetric (ρ̄ ≠ ρ), so conditioning is a **singular-value**
question; the symmetric part's spectrum is reported separately.

**The LINEAR control, all the way through a completed collapse** (h0 = 1.0,
sampled from s/B = 0.008 to 0.146, ending at M = 0.9999 on a dead-flat plateau):

| s/B | 0.0077 | 0.0202 | 0.0477 | 0.0792 | 0.1167 | 0.1457 (wall) |
|---|---|---|---|---|---|---|
| σ_min/scale | 1.37e-3 | 1.45e-3 | 1.31e-3 | 1.33e-3 | 1.13e-3 | **1.20e-3** |
| cond | 6.9e3 | 6.5e3 | 7.3e3 | 7.1e3 | 8.5e3 | **7.9e3** |
| n_neg(sym) | 8 | 2 | 6 | 5 | 9 | 6 |

**Flat. The linear element reaches a genuine, fully isochoric, perfectly
plateaued collapse mechanism WITHOUT its tangent ever going singular** — σ_min
at the plateau equals σ_min in the elastic range to within 25 %.

That falsifies the assumption written into the harness's own ladder docstring:

> "A perfectly plastic collapse is where Newton is worst: the tangent goes
> singular in the mechanism mode exactly when the answer arrives."

Not on this problem as discretised. The mechanism is a smeared plastic zone, not
a velocity discontinuity, and it does not produce a singular global operator.

**The QUADRATIC leg** (h0 = 1.0, same sampler, same mesh):

| s/B | 0.00014 | 0.0040 | 0.0092 | 0.01088 | 0.01115 | **0.01118** | 0.01122 (wall) |
|---|---|---|---|---|---|---|---|
| σ_min/scale | 2.38e-4 | 2.35e-4 | 2.22e-4 | 2.18e-4 | **2.16e-4** | **2.29e-7** | 2.41e-7 |
| cond | 8.5e4 | 8.7e4 | 9.2e4 | 9.4e4 | **9.5e4** | **5.2e9** | 5.0e9 |
| n_neg(sym) | 0 | 3 | 16 | 12 | 12 | 27 | 32 |

**σ_min is flat to within 10 % over the entire path — including at the sample
immediately before the wall — and then collapses by a factor of 1000, with the
condition number rising 55 000×, between two adjacent samples separated by
Δ(s/B) = 3e-5 (0.06 mm of settlement).** After the jump it does not continue
toward zero; it sits on a new plateau. The locked leg (`h20std_h1.0`) shows the
same signature at its own wall: σ_min/scale 4.6e-7, cond 3.0e9.

**Cross-check.** An earlier, independent run of the same configuration sampled
the tangent by a *different* mechanism (the whole analysis in `FullGeneral`,
sampled in-loop, no zero-increment re-form) and measured σ_min/scale = 2.209e-07,
cond = 5.405e+09 at s/B = 0.01122. Two sampling paths, same answer.

**Two readings this rules out.**

- *"The operator is going singular in the mechanism mode."* No: it is
  well-conditioned at 0.01115, and the failure is a step change, not an
  approach. A limit point degenerates smoothly.
- *"Loss of positive definiteness marks the collapse."* No: λ_min(sym) goes
  negative at **s/B = 0.0012, q = 33 % of exact**, and the **linear** element —
  which completes the collapse — carries 2–11 negative eigenvalues throughout.
  For a non-associated (unsymmetric) tangent this is ordinary, not diagnostic.

> **The quadratic element's operator goes bad ABRUPTLY, at a settlement where
> the curve is still hardening, the mechanism index equals the linear
> element's, and the linear element on the same mesh has a flat, healthy
> tangent. That is the signature of a numerical event, not a physical limit
> point.** What the event *is* was not identified — see §8.

## 8. What came back inconclusive, and what would settle it

Reported deliberately; two prior retractions on this problem cost more than the
caution would have.

1. **The identity of the abrupt tangent event — the single biggest open item.**
   Measured that it happens, that it is abrupt, that it is not a limit point.
   *Not* measured: what changes. The leading candidate is a **return-map branch
   switch at the Drucker–Prager APEX** — `SY = 0.2 kPa` is a tiny apex
   regularizer, and the apex branch has a qualitatively different consistent
   tangent. **Test:** instrument `DruckerPrager` to count apex-branch Gauss
   points per step and re-run the trajectory; a step change in that count at
   s/B = 0.01118 would settle it. Raising `SY` and watching the wall move is the
   confirmatory half.
2. **The checkerboard / pressure-mode hypothesis (§6.3)** — raised by the field
   maps, not supported by the deviatoric measure, and *reverses* with mesh
   refinement. **Test:** dump `I₁` per Gauss point and run the same neighbour
   residual on the **pressure** field, which is the classical instrument. The
   present result neither supports nor excludes it.
3. **A path method that does not need the tangent.** Arc-length is inapplicable
   here for the structural reason in §3, and was not run. The genuinely
   different path is **dynamic relaxation / a quasi-static transient**, for which
   the fork already owns `LadrunoFictitiousMass` and `LadrunoDynamicRelaxation`
   (ADR 21). Untried, and the most valuable remaining experiment: it would
   determine whether a quadratic element can complete this collapse *at all*.
4. **The SIZE of the allowance effect, and where the recession ends.** §7.1.1
   measured a 5.1 % capacity / 14.7 % reach spread between two **identical**
   runs, which is the same order as the effect being measured. Direction is
   established (2/2); magnitude is not. **Test:** 4–6 repeats at each of 3
   floors, reported as distributions. ~10 leg-hours, and it is the honest price
   of a number here.
5. **How far the run-to-run scatter of §7.1.1 extends.** It was measured on one
   configuration. Whether *plateaued* legs (the linear ones, which are the fork's
   gate) carry the same exposure is untested here — the merged gate's ±3 % bands
   assume they carry less, which is plausible (a plateau is not a convergence
   failure) but is now an assumption with a measured counter-example nearby.
   **Test:** repeat the R3 gate's h0 = 0.5 leg 5× and report the spread. Cheap,
   and it directly protects a standing CI gate.
6. **Resolution of the localisation metrics.** At h0 = 1.0 the core holds 18
   elements and at h0 = 0.5 it holds 96; neither can resolve a 1–3-element slip
   band even if one existed. The "no slip surface anywhere" finding is therefore
   a statement about *these meshes*, not about the elements.

## 9. Recommendation on a 27-node / H20-selective element lane: **DO NOT OPEN IT ON THIS EVIDENCE**

Note 81 §5.5 proposed H27 + selective integration as "the quadratic element most
likely to have the richer isochoric span". That proposal rests on a premise this
note measured and **refuted**: at the point of failure the quadratic element's
isochoric span is not the binding constraint (§6.1 — M = 0.788 vs the linear
element's 0.803 at the identical settlement, on the way to the exact answer).

Against that, the three measured facts:

- the mechanism forms in the quadratic element at the same rate as in the linear
  one (§6.1);
- the wall moves on a pure controller parameter (direction consistent, magnitude
  unresolved — §7.1.1) and **not at all** on a different algorithm (§7.2), and
  it is not even reproducible to 5 % between identical runs;
- the operator is well conditioned until an abrupt, unidentified numerical event
  (§7.3).

Every one of those points at the **path and the material integration**, not at
the element's function space. Building a 27-node hex — a new element, its shape
functions, ordering memo, class tag, broker case, recorder type maps, Tcl and
Python commands, mass and body-force paths, and its own control battery — is
weeks of element development aimed at a deficiency that has not been
demonstrated to exist.

**Recommendation:**

1. **Do not open the lane now.** The three experiments in §8 items 1–3 are
   instrumentation, cost days rather than weeks, and each of them could remove
   the motivation for the lane entirely — or, if the apex/branch hypothesis is
   refuted and dynamic relaxation *also* fails to complete the collapse with a
   quadratic element, supply the first real evidence *for* it.
2. **Run §8 item 3 (dynamic relaxation) first.** It is the only experiment that
   can answer "can any quadratic element complete this collapse", which is the
   question H27 was proposed to answer, at a small fraction of the cost.
3. **If a quadratic-hex lane is opened later, justify it somewhere else.**
   Quadratic elements exist to win on bending- and contact-dominated problems,
   and note 81 §5.7 already scopes this measurement to one problem, one flow
   rule, one material and one mesh family. A H27 built for *those* reasons is a
   defensible lane; a H27 built to fix this collapse leg is aimed at a target
   this note could not find.
4. **Note 81 §5.6's F-bar-Patch verdict is unchanged and, if anything,
   strengthened** — for the same reason: it too was a span-side remedy.

## 10. Scope

One problem (Prandtl–Reissner), one flow rule (non-associated, ψ = 0), one
material (`DruckerPrager`, φ_ps = 27.47°, σ_y = 0.2 kPa regularizer), one mesh
family (structured graded strip, plane strain), two resolutions, one
displacement-controlled path family. Nothing here carries to associated flow, to
severe non-associativity, to dynamics, or to bending-dominated problems. The
conditioning statements are for **this** discretisation: a mesh fine enough to
resolve a true velocity discontinuity might well produce the singular tangent
this one does not.

## 11. Reproducing

```bash
cd Ladruno_files/testbed/hypo_bearing
py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --calib          # control CHI
py -3.12 quad_path_diag.py --elem h8bbar --h0 0.5 --ladder linear --budget 80
py -3.12 quad_path_diag.py --elem h20uri --h0 0.5 --budget 200
py -3.12 quad_path_diag.py --elem h20uri --h0 0.5 --budget 800 --floor 2e-9 \
                           --suffix _loose                          # the allowance pair
py -3.12 quad_path_diag.py --elem h8bbar --h0 0.5 --ladder linear --sfrac 0.0103 \
                           --suffix _at0103                         # matched settlement
py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --cond --cond-at 1.5e-4 --suffix _traj
py -3.12 quad_path_summary.py
py -3.12 quad_path_osc.py
```

Every leg must be quoted with its termination mode. Any leg quoted as a capacity
must be run at **two controller allowances** and both reported.
