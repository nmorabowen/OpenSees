# ADR-75b L3-0 — the Lane-3 gate table: MEASURED

**Verdict, in one line: the Lane-1 solver win is what makes Lane 3 worth doing.** On the same
Lane-B model, the same binary, and the same mesh, `system Pardiso -matrixType 2` at 4 threads moves
the element-kernel fraction from **35.8% → 74.9%** and the solve fraction from **55.9% → 16.9%**.
ADR-75 §3's Lane-3 gate (>40% element fraction) **fails under UmfPack and passes with enormous
margin under PARDISO**. Lane B was a solver-bound lane; PARDISO made it an element-bound one.

Date: 2026-07-25 · worktree `patricio-palacios-prs-17c1e2` · `dist/bin/opensees.pyd` built this
session from `ladruno` @ `bdf8adf9f` (PARDISO present, 3 lanes measured).
Sub-ADR: [[75b_ladruno_threaded_assembly_adr]] (stage L3-0). **No threading code exists; this is
measurement only.**

---

## 0. What this measures, and why the shipped profiler alone was not enough

ADR-75 §3 gates every Lane-3 stage on a **">40% element/`update` fraction"**. The shipped profiler
reports per-**phase** totals — but a phase is not a threadable loop. `formTangent` contains **both**
the element kernel (threadable, if the element is re-entrant) **and** the `theSOE->addA` scatter (the
race ADR-75b §4 has to solve). Threading buys the first and not the second, so a phase percentage
**systematically overstates the Amdahl headroom.**

So L3-0's deliverable is the split, per lane, for each of the three loops of ADR-75b §2:

| loop | site | kernel scope | reduction? |
|---|---|---|---|
| **A** | `Domain::update()` `Domain.cpp:2396-2406` | `elem.update` | **none** |
| **B** | `IncrementalIntegrator::formTangent` `:112-125` | `elem.tangent` | `addA` |
| **C** | `IncrementalIntegrator::formElementResidual` `:323-337` | `elem.residual` | `addB` |

recovered as `kernel = Σ(elem_by_type wall under the loop scope)` and
`scatter = loop-scope wall − kernel`, summed over **every** site at which the loop scope appears.

> **Correction to ADR-75 §3 / §6-P4(a), which says the first Lane-3 move is to *build* the
> `elem.update` instrument: it was already built and shipped 2026-07-06** (`Domain.cpp:2394`,
> `:2403`; macro at `ProfilerMacros.h:114`), and ADR-40b's own addenda already used it. The real gap
> was that the >40% gate had never been evaluated **per threadable loop**. That is what this document
> closes.

### Method
- Models: the **unmodified ADR-40b lane models** (`phase0/laneA_fiber_frame.py`,
  `laneB_model.py`, `laneD_model.py`), so numbers are comparable to ADR-40b. The only change is that
  `OPS_PYD` / `OPS_PROF_OUT` / `OPS_SYSTEM` / `OPS_PROF_FLAGS` are now env-overridable; **defaults
  reproduce the original runs verbatim.**
- **Median of 3 interleaved rounds** per configuration (the three Lane-B solver configs are
  interleaved within each round — this box swings ±30% on background load,
  `phase0/PHASE0_RECIPE.md`). Lane-B rows are the interleaved set; fractions reproduced to ±0.2 pp.
- Threads: `MKL_NUM_THREADS` = 1 except the `pard4` row (4). Lane 3 does not exist yet, so no
  assembly threads anywhere.
- Full logs kept, never grepped down (banked P1d trap: a silent `ProfileSPDLinSOE` fallback hides in
  exactly the lines you would have dropped). `OPS_SYSTEM` numeric values are converted to **ints or
  floats**, never left as strings — a string fails `OPS_GetIntInput`/`OPS_GetDoubleInput` and the
  factory falls back silently.

---

## 1. The gate table (medians of 3)

| lane / config | step ms | solve % | **kernel %** | loop % | loop A % | loop B % | loop C % | `addA` % of loop B | Amdahl 4T (kernel) |
|---|---|---|---|---|---|---|---|---|---|
| **A** fiber frame, UmfPack | 7 919 | 5.65 ⚠ **(true 7.86 — see F8)** | **81.94** | 83.32 | **80.79** | 0.45 | 0.79 | 56.05 | 2.59× |
| **B** 3D solid, UmfPack 1T | 19 191 | 55.89 | **35.84** | 43.04 | 2.55 | 29.42 | 3.90 | 19.09 | 1.37× |
| **B** 3D solid, Pardiso `-matrixType 2` 1T | 11 879 | 36.39 | **57.22** | 62.53 | 4.04 | 46.81 | 6.20 | 9.47 | 1.75× |
| **B** 3D solid, Pardiso `-matrixType 2` 4T | 9 200 | 16.87 | **74.85** | 81.73 | 5.24 | **61.56** | 8.04 | 9.45 | 2.28× |
| **D** explicit CDL, Diagonal | 29 468 | 0.03 | **86.84** | 89.66 | **54.35** | ~0.00 | 32.34 | n/a | 2.87× |

`kernel %` = element kernel calls only (the threadable part). `loop %` = kernel + scatter.
`Amdahl 4T` = `1/((1−p)+p/4)` with `p` = kernel fraction: a **ceiling** assuming perfect kernel
scaling and zero fork-join cost, not a prediction.

### Per-element cost (median) — the input to the fork-join regression check

| lane | classTag | loop A `update` | loop B `getTangent` | loop C `getResidual` |
|---|---|---|---|---|
| A fiber frame | 73 `ForceBeamColumn2d` | **20.53 µs/ele** (n=311 610) | **0.12 µs/ele** (n=285 090) | 0.19 µs/ele (n=311 610) |
| B 3D solid | 33002 `LadrunoBrick` | 3.25–3.27 µs/ele | **37.28–38.14 µs/ele** | 3.70–3.75 µs/ele |
| D explicit | 33002 `LadrunoBrick` | 3.20 µs/ele (n=5 000 000) | — | 3.83 µs/ele (n=2 500 000) |

---

## 2. Findings

### F1 — PARDISO flips Lane B's Lane-3 gate. This is the headline.
Same model, same binary, same mesh; only the solver differs:

| | solve % | kernel % | gate (>40%) | Amdahl 4T |
|---|---|---|---|---|
| UmfPack 1T | 55.89 | 35.84 | **FAIL** | 1.37× |
| Pardiso `-matrixType 2` 1T | 36.39 | 57.22 | PASS | 1.75× |
| Pardiso `-matrixType 2` 4T | 16.87 | **74.85** | **PASS (huge margin)** | 2.28× |

ADR-75b §1 argued this would happen; it is now measured. Two consequences:
1. **Lane 3's value is created by Lane 1.** Threading assembly against UmfPack would have been a
   marginal, gate-failing project. Against the shipped desktop solver it is the dominant remaining
   lever on the fork's *primary* lane.
2. **Any Lane-3 gate must be evaluated on the production solver configuration.** Quoting the
   ADR-40b UmfPack-era 30.1% element fraction would have wrongly killed the lane.

### F2 — Kernel time is solver-invariant, which validates the whole measurement
> **Weakened by the adversarial review (§7 finding 2).** The first draft quoted "within 2.3%" over the
> three configs it had run. Adding the `-matrixType 0` config the review introduced widens the spread
> to **5.8%**. The check still passes — but 2.3% was a selective figure, so quote 5.8%.

Absolute loop-B kernel wall across **all four** configurations: **5 646.6 (UmfPack) / 5 536.5
(Pardiso-1T) / 5 663.4 (Pardiso-4T) / 5 857.9 (Pardiso unsym-4T) ms — a 5.8% spread**
(`LadrunoBrick::getTangent` = 38.02 / 37.28 / 38.14 / ~39.5 µs/ele).

The element kernel *cannot* depend on which solver factorizes afterwards, so the residual spread is
this box's measurement noise, not signal — and 5.8% is therefore a **useful noise floor for every
other number here**. The check does its job: the F1 fraction shift (35.8% → 74.9%, a factor of 2.1)
is an order of magnitude outside that floor, so it is the *denominator* moving, exactly as claimed,
and not an instrumentation artifact.

### F3 — `-matrixType 2` cuts the `addA` scatter by 65%: a third, previously unmeasured benefit of symmetric storage
> **This finding was corrected by the adversarial review (§7, finding 1).** The first draft
> attributed the win to *PARDISO's SOE* being a faster assembler than UmfPack's. That was wrong, and
> the review killed it by measuring the one configuration the sweep had omitted: **unsymmetric
> PARDISO.**

Absolute median `addA` scatter in loop B (n=3 each, all at the same mesh):

| config | kernel ms | `addA` scatter ms | scatter / loop |
|---|---|---|---|
| UmfPack 1T | 5 646.6 | 1 322.8 | 18.98% |
| **Pardiso `-matrixType 0`** (unsym) 4T | 5 857.9 | **1 699.2** | 22.48% |
| Pardiso `-matrixType 2` (sym) 4T | 5 663.4 | **591.0** | 9.45% |

Two conclusions, both different from the draft:

1. **The 65% scatter cut (1 699.2 → 591.0 ms) belongs to half-storage, not to PARDISO.**
   `-matrixType 2` scatters only the `col >= row` entries, so `addA` does about half the work — and
   the searches it does perform run over half-length rows. This is a **third** measured benefit of
   `-matrixType 2`, after P1d's ~1.25× time and −41.8% peak memory. P1d never measured it.
2. **PARDISO's *unsymmetric* `addA` is the slowest of the three — 1.28× slower than UmfPack's.**
   That is the O(idSize² × rowlen) linear search of `PARDISOGenLinSOE::addA` (`:439-455`) paying for
   itself on full-length rows. At 1 699.2 ms of an ~11.9 s step it is **~14% of wall on the
   default-configured desktop solver path**, spent locating entries in a structure that has not
   changed since `setSize`.

So ADR-75b §11's open question 5 — replace that search — is no longer a stylistic note: it is a
**~14%-of-wall serial optimization on the default path**, with no determinism cost and no threading
required. On the present evidence it is a better next move than threading loops B/C at all.

**And it is bigger than "a PARDISO fix" — the review checked the other SOEs.** `UmfpackGenLinSOE::addA`
has the *identical* shape on frozen CSC (`for k in Ap[col]..Ap[col+1]: if (Ai[k]==row) { Ax[k] += …;
break; }`), so **both** desktop solvers pay the same O(idSize² × rowlen) search — UmfPack's costs
1 322.8 ms on the same model. Replacing it is a shared win, not a PARDISO patch.
(`DiagonalSOE::addA` is `A[pos] += m(i,i)` into a dense diagonal — no search, nothing to fix, which
is one more reason Lane D is the clean L3-1 target.)

*De-confounding note:* the UmfPack row ran at 1 thread and the unsym-PARDISO row at 4, so the 1.28×
could in principle be a thread-count artifact. It is not: `-matrixType 2` scatter is 579.3 ms @1T vs
591.0 ms @4T — a 2% difference — so thread count does not materially move `addA`, as expected for a
serial routine.

### F4 — Loop A has essentially no non-kernel overhead, confirming the §2.1 determinism claim structurally
Non-kernel share of loop A: **0.63%** (lane A), **4.5–4.8%** (lane B), and there is **no scatter term
at all** — `Domain::update()` has no `addA`/`addB`. That is the measured face of ADR-75b §2.1: loop A
performs no floating-point reduction, so a correct threaded loop A is **bit-identical to serial at
every thread count**, and the existing byte-identical oracles gate it unchanged. The highest-payoff
loop carries zero determinism cost.

### F5 — Payoff and threadability are inverted, and now quantified
Lane A: loop A is **80.79% of step** at **20.53 µs/ele** — a fraction and a per-element cost that
both say "thread this". But the element is `ForceBeamColumn2d`, a **vanilla upstream** file whose
`update()` runs its entire interior Newton on **14 function-scope statics** (ADR-75b §5.2), and
ADR-40b already ruled lane-A optimization out of the fork's change budget.

Meanwhile Lane D's loop A is **54.35% of step** at 3.20 µs/ele on `LadrunoBrick` — **fork-owned**,
already authorized for de-statication (ADR-40a C16 / ADR-68 T3), and with **no SOE interaction
whatsoever** (`system Diagonal`, `linearSolve` = 0.03%).

**So Lane-3 staging cannot be ordered by element fraction alone.** Ordered by
*fork ownership × re-entrancy cost × measured fraction*, the ranking is §3 below.

### F6 — For cheap elements the assembly loop is scatter-bound, not element-bound
On Lane A's tangent loop the **scatter (44.1 ms) exceeds the kernel (35.4 ms)** — 56.05% of the loop
is `addA`. `ForceBeamColumn2d::getTangent` measures **0.12 µs/ele**, confirming the number ADR-75 §3
warned about ("≈0.1 µs, can *regress* from fork-join overhead").

But note the precise reason not to thread it: loops B and C on Lane A are **0.45% and 0.79% of
step** — Amdahl-irrelevant. That, not fork-join, is the decisive argument. Fork-join is the secondary
hazard, and it cannot be quantified yet (see §5).

### F7 — Loop scopes appear at multiple sites; summing one undercounts ~2×
Lane A's `elem.update` appears at **three** sites — `newStep/elem.update` 941.7 ms,
`solveCurrentStep/elem.update` 581.0 ms (the `DisplacementControl` direct call), and
`solveCurrentStep/update/elem.update` 3 952.7 ms. Lane D's appears at two
(`newStep` 8 064.1 ms + `solveCurrentStep/update` 8 534.2 ms — the double constitutive pass ADR-40b
found, confirmed here by `n = 5 000 000 = 2 × 2500 steps × 1000 elements`). `parse_lane3.py` sums all
sites; this is the same trap ADR-40b hit with the hidden second `soe.factor`.

### F8 — Lane A's `solve %` in §1 is UNDERSTATED: it has a hidden factorization, exactly like ADR-40b's lane E
*(Found by the adversarial review, §7 finding 3 — a wrong published number.)*

The §1 `solve %` column reports the `linearSolve` **phase**. On Lane A that is not the solver cost:

| lane A phase | ms | % step |
|---|---|---|
| `linearSolve` (the algorithm's visible solve) | 449.4 | **5.65** |
| `soe.factor` (**all sites**) | 505.7 | 6.39 |
| `soe.trisolve` (all sites) | 120.4 | 1.52 |
| **true solver work = factor + trisolve** | **626.1** | **7.86** |

`soe.factor` (505.7 ms) **exceeds** `linearSolve` (449.4 ms), so factorization is running outside the
visible solve — Lane A uses **`DisplacementControl`**, whose `update()`/`newStep()` solve for the
reference displacement `dUhat` against the same `K`. This is precisely the pathology ADR-40b found on
lane E ("59% of the step is UmfPack numeric factorization, two-thirds of it booked outside
`linearSolve`"), and the first draft walked into it despite that finding being in the cross-validation
table two sections below.

**Nothing in the verdicts changes** — lane A is 81.94% element-kernel either way, and 7.86% is still
negligible against it. But the column was wrong and is now marked.

**Checked for every other lane, and Lane A is the only one affected:** Lane B/UmfPack has
factor+trisolve = 55.32% against `linearSolve` 55.89% (i.e. essentially all factorization is inside
the visible solve — LoadControl + Newton, no hidden `dUhat` solve), and Lane D is ~0 both ways.

**Related instrumentation gap, worth recording:** the `soe.factor` / `soe.trisolve` scopes exist
**only in `UmfpackGenLinSolver.cpp`** (`:195`, `:211`, `:234`, `:249`). **PARDISO has none.** So the
`0.00%` in the factor+trisolve column for every `Pardiso` row above means *"not instrumented"*, not
*"no cost"* — and, more importantly, **ADR-40b's rank-8/10 factor-vs-triangular-solve split does not
exist for the fork's shipped desktop solver.** Adding it is cheap and would make the same
hidden-solve check possible on PARDISO decks.

---

## 3. What this changes in ADR-75b §7 (recommended stage order)

| rank | target | fraction | element | why |
|---|---|---|---|---|
| **1** | **Lane D loop A** (`Domain::update`, explicit) | 54.35% | `LadrunoBrick` (fork-owned) | bit-identical by construction (F4); no SOE interaction at all; de-statication already authorized; ADR-40b already passed rank 7's gate here |
| **2** | **Lane B loop A** (implicit) | 5.24% | `LadrunoBrick` | same code path as rank 1 — near-free once rank 1 lands, though small on its own |
| **3** | **Lane B loop B** (`formTangent`) | **61.56%** | `LadrunoBrick` | the single largest kernel bucket in the fork (38 µs/ele). Needs the §4 scatter work (fast + ordered), so it is gated on the §4.2 gather-memory question |
| **4** | Lane B loop C (`formUnbalance`) | 8.04% | `LadrunoBrick` | rides on rank 3's machinery |
| **defer** | Lane A loop A | 80.79% | `ForceBeamColumn2d` (**vanilla**) | biggest fraction, worst re-entrancy, outside the change budget (F5) — upstream-facing |
| **never** | Lane A/C loops B, C | 0.45% / 0.79% | — | Amdahl-irrelevant and scatter-bound (F6) |

This **confirms** ADR-75b §7's ordering (L3-1 loop A on fork-owned elements first, `ForceBeamColumn2d`
last) and **refines** it: Lane **D** — not Lane B — is the cleanest first target, because its
`system Diagonal` path removes the SOE from the picture entirely.

---

## 4. Cross-validation against ADR-40b

| quantity | ADR-40b | here | reconciliation |
|---|---|---|---|
| Lane A classTag-73 state determination, % of step | ~82% | **81.94%** | reproduces |
| Lane A `ForceBeamColumn2d::update` µs/ele | 57.6 | **20.53** | fraction reproduces, absolute does not: ADR-40b's lane-A step was ~22 s vs **7.9 s** here — a ~3× faster box/build. µs/ele is a wall-derived rate, so it scales with the machine; the *fraction* is the portable quantity |
| Lane D element work, % of step | ~89% | **86.84%** | reproduces |
| Lane D `formTangent` | 12.9 s (16.2%) | **1.3 ms (0.00%)** | **not a regression — the number is stale.** ADR-67 P-NEW-1's constant-mass tangent cache ships `massCache = true` by default (`CentralDifferenceLadruno.cpp:89`, `:154`, override `:259`), which is the "`-factorOnce` with safe invalidation" fix ADR-40b's Finding-3 item 1 asked for. Banked in `LEDGER_quirks` |
| Lane B `linearSolve` (UmfPack), % of step | 66.4% | **55.89%** | same model and mesh; the gap is machine/build and the phases ADR-40b did not yet split. Directionally identical; the *cross-solver* comparison in F1 is within-session and interleaved, so it does not depend on this |

---

## 5. Caveats — stated, and one of them measured rather than argued

**Deep-scope instrumentation tax — measured.** ADR-40b banked ~0.5 µs per deep-scope instance, which
would predict a large tax here (Lane D carries ~15.25 M scope instances: 7.5 M `elem.*` + 7.75 M
`brick.geo`). Rather than argue about it, the same three lanes were re-run with
`OPS_PROF_FLAGS="-perStep"` (coarse — no `-deep`, so no per-element scopes), interleaved with the
deep runs, 6 rounds each:

| lane | deep median s | coarse median s | **tax** | extra s | scope instances | implied µs/instance | kernel % raw | **kernel % corrected** |
|---|---|---|---|---|---|---|---|---|
| A fiber frame | 8.13 | 7.78 | **+4.5%** | 0.35 | 908 310 | 0.385 | 81.94 | **81.13** |
| B 3D solid (pard4) | 9.91 | 9.62 | **+3.0%** | 0.29 | 941 625 | 0.303 | 74.85 | **74.10** |
| D explicit | 28.89 | 26.74 | **+8.0%** | 2.15 | 15 250 000 | 0.141 | 86.84 | **85.78** |

**The correction moves every kernel fraction by ≤1.1 percentage points, so no verdict in this document
changes.** ("corrected" attributes the entire deep-vs-coarse delta to the element buckets, which is the
worst case for these fractions, since that is where the per-element scopes live.)

Two honest qualifications:
- **Only Lane D's tax is statistically established.** Paired sign test over the 6 rounds
  (one-sided, H₀ = deep is no slower than coarse):

  | lane | deep > coarse | p | verdict |
  |---|---|---|---|
  | A fiber frame | 4/6 | 0.344 | **not established** |
  | B 3D solid | 5/6 | 0.109 | **not established** |
  | D explicit | **6/6** | **0.016** | **established** |

  So the honest reading is: **the tax is demonstrated only on Lane D (+8.0%), the lane with by far the
  most scope instances (15.25 M) — which is exactly where it should be largest.** The lane-A and
  lane-B medians (+4.5%, +3.0%) are the right order of magnitude but are **not** separable from noise
  at n=6; individual rounds even show coarse slower than deep (lane A r3: 8.18 coarse vs 7.31 deep).
  The noise floor is large: max/min spread was 1.22×/1.28× (lane A deep/coarse), 1.32×/1.20×
  (lane B), 1.62×/1.34× (lane D) — and F2 independently puts this box's floor near 5.8%.
  **This does not affect any verdict**: the correction is ≤1.1 pp against gate margins of tens of
  points, and Lane D — the one lane where the tax *is* established — still reads 85.8% corrected.
- **Implied cost is 0.14–0.39 µs per scope instance, i.e. at or somewhat below ADR-40b's banked
  ~0.5 µs.** Lane D's 0.141 µs is the figure with the most instances behind it. ADR-40b's caveat
  stands in substance — per-**GP** scopes (8–16 per call) would still be ruinous, and its
  timer-bound-leaf warning is why `getTangent` at 0.12 µs/ele must not be read as a true cost.

**Other caveats.**
- **Single box, single build.** Absolute walls are not portable across sessions (±20–30%); every
  comparison that carries a verdict here is *within one interleaved run*.
- **Lane A's pushover terminates early — deterministically — at step 86 of 150.** All three rounds
  print `FAILED at step 86, stopping` and end at `roof disp = 29.722 in` (0.69% drift vs the 1.2%
  target), i.e. the `KrylovNewton`+substep retry does not clear it. Identical work every round, so
  the fractions are self-consistent and reproducible, but **this is a partial pushover**, and it
  explains part of the gap to ADR-40b's "93 steps". Pre-existing testbed behaviour, not introduced
  here and not fixed here; worth a follow-up for whoever next needs lane A to reach target drift.
- **Lane B is 11 520 DOF.** ADR-75 §1's standing correction applies: the production regime is far
  larger, and the solve fraction grows with N (measured in P1c: 1.61× → 3.40× UmfPack→PARDISO from
  11.5k → 51k DOF). At production size the solve fraction is **higher** than measured here, so the
  element fraction — and hence Lane 3's gate margin — is **lower** than the 74.85% headline. **The
  F1 gate verdict should be re-checked on a production-size deck before L3-3/L3-4 are committed.**
  It is not in doubt for lanes A and D, whose fractions are ~81–87% with no solver in the picture.
- **Amdahl numbers are ceilings, not predictions** — perfect kernel scaling, zero fork-join, no
  memory-bandwidth saturation. P1 already measured PARDISO flattening past 4 threads on this box
  (memory-bandwidth-bound); element kernels will hit the same wall.
- **`getTangent` at 0.12 µs/ele is timer-bound** (ADR-40b: a deep scope instance costs more than
  that), so the true cost is *lower* and Lane A's 0.45% loop-B fraction is an **over**estimate —
  which only strengthens F6.
- **Lane C (shell) and Lane E (IMK) were not re-measured.** ADR-40b covers them (element 55.9% and a
  refuted hypothesis respectively); neither changes the staging in §3.

---

## 6. What L3-0 does NOT answer (feeds ADR-75b §11)

1. **The fork-join / barrier cost on this toolchain, in µs.** It sets the per-element-cost threshold
   in every stage gate. **Not measurable yet: the build has no OpenMP** — `find_package(OpenMP)` is
   absent from the fork's CMake and there is no `/openmp` flag, so the 7 `#pragma omp` lines in
   `SRC/` (all PFEM) are compiler-ignored no-ops. First measurement available at L3-1.
2. **The §4.2 gather memory (Σ idSize² doubles) at production scale.** ~15.5 MB at Lane B,
   ~1.5 GB projected at ~1M DOF. This single number decides whether the byte-identical class-B CI
   gate of ADR-75b §3 survives. Highest-leverage unknown in the lane.
3. **Whether atomic `A[k] +=` contends.** Only measurable once L3-4 exists.
4. **How much of P1's ~34% Amdahl residual is loop A vs B/C vs neither.** Partially answered — on
   Lane B at Pardiso-4T the three loops account for **81.73%** of step, so the residual is
   overwhelmingly *element* work, not "neither". Good news for the lane.
5. **Production-deck re-validation of F1** (see §5).

---

## 7. What the adversarial review changed

Two scoped passes ran against this document, both targeting the claims that would be **wrong without
being obviously wrong**. Between them: **five real defects**, one of them a wrong published number.

### Pass 2 (post-PR, deeper — verification by execution rather than re-reading)

| # | Finding | Status |
|---|---|---|
| 1 | **A published number was wrong.** §1's `solve %` for Lane A reported the `linearSolve` *phase* (5.65%). But `soe.factor` **exceeds** `linearSolve` there (505.7 vs 449.4 ms) — `DisplacementControl` solves for `dUhat` outside the visible solve, so true solver work is **7.86%**. This is precisely ADR-40b's lane-E pathology, cited in this document's own cross-validation table two sections below | **F8 added**, §1 marked; verified Lane A is the *only* affected lane |
| 2 | **F2's "solver-invariant to 2.3%" was a selective figure.** Including the `-matrixType 0` config that pass 1 itself introduced widens the spread to **5.8%** | F2 rewritten; 5.8% is now used as this box's stated noise floor |
| 3 | **The instrumentation tax was over-claimed as a 3-lane result.** A paired sign test over the 6 rounds establishes it **only for Lane D** (6/6, p=0.016); Lanes A and B are 4/6 (p=0.344) and 5/6 (p=0.109) — not separable from noise | §5 now reports the sign test and scopes the claim to Lane D |
| 4 | **The `addA` search fix was under-sold as PARDISO-specific.** `UmfpackGenLinSOE::addA` has the identical frozen-CSC search, so **both** desktop solvers pay it | F3 + ADR §11 q5 corrected to a both-solvers fix |
| 5 | **A latent trap was re-introduced in the harness.** The `OPS_SYSTEM` int-only converter left float option values as strings (`Mumps -BLR 1e-8` → `'1e-8'`), which would fail `OPS_GetDoubleInput` and trigger the *exact* silent `ProfileSPDLinSOE` fallback the harness comment warns about | fixed (ints **and** floats) and re-run end-to-end |
| — | **Arithmetic:** gather projection at 1M DOF is **~1.5 GB**, not 1.4 GB; `static Matrix` count is **1,686**, not ≈1,650 | corrected here, in the ADR, and in `LEDGER_quirks` |

**Cleared in pass 2, by execution rather than argument:**
- **The parser** — `parse_lane3.py` was run against **200 randomized synthetic h5 files** with known
  ground truth (1–3 sites per loop, 1–2 classTags per site): `loop_ms`, `kernel_ms`, `scatter_ms` and
  `step_ms` matched exactly, **0 mismatches**. Multi-site summing and the kernel/scatter split are
  correct by construction, not by inspection.
- **Every published median** re-derived from the raw JSON through an independent code path — all
  five configs' kernel %, loop-A %, and Amdahl 4T reproduce to the last quoted digit.
- **No negative scatter** (`kernel > loop`) in any of the 21 real runs, i.e. the `elem_by_type` bucket
  never over-counts its enclosing scope.
- **§4.1's frozen-sparsity generalization holds beyond PARDISO** — `UmfpackGenLinSOE` (frozen CSC +
  one `Ax[k] +=`) and `DiagonalSOE` (`A[pos] += m(i,i)`, dense, no search) both fit.
- **The 1.28× UmfPack-vs-unsym-PARDISO scatter gap is not a thread-count artifact** — `-matrixType 2`
  scatter is 579.3 ms @1T vs 591.0 ms @4T (2%), so thread count does not move `addA`.

### Pass 1 (pre-PR)

It found two real defects and cleared the one most likely to be an artifact.

| # | Finding | Status |
|---|---|---|
| 1 | **F3 was mis-attributed.** The draft credited PARDISO's SOE with a 2.24× faster `addA`. The sweep had measured only `-matrixType 2`, so "PARDISO" and "half-storage" were confounded. Measuring `-matrixType 0` separates them: unsym PARDISO `addA` = **1 699.2 ms, i.e. 1.28× SLOWER than UmfPack's 1 322.8 ms**. The win is half-storage's, and PARDISO's unsymmetric scatter is the *worst* of the three | **F3 rewritten**; the O(rowlen) search is now quantified at ~14% of wall and promoted in ADR-75b §11 |
| 2 | **ADR-75b §2.1's claim that each `update()` "writes ONLY this element's own state" is too strong.** At least 13 element files write **node trial state**, and `LadrunoRigidBody::update()` (`:490` → `imposeSlaveKinematics()` `:500`, writes at `:518`/`:536`/`:540`) does so **on the loop-A path, to nodes shared with other elements**. That is an *ordering/correctness* race, not merely an FP-order one — so loop A's bit-identicality holds only for allowlisted elements, which is stronger than "the reduction is absent" | ADR-75b §2.1 qualified, §7 exclusion list made concrete, §9 risk added |
| 3 | **Cleared: is the "scatter" figure really `addA`, or is `zeroA` inside the scope?** If `zeroA()` were inside `elem.tangent`, the whole kernel-vs-scatter split — and F3 — would be an artifact of UmfPack storing a full unsymmetric `A`. It is not: `theSOE->zeroA()` is at `IncrementalIntegrator.cpp:100`, the scope opens at `:109`. Scatter = `addA` + loop iteration only | verified sound, no change |
| 4 | **Not claimed rather than asserted:** the draft reasoned that materials are `getCopy()`'d per element so cannot alias across elements. That holds for the standard library but was not exhaustively verified, so it is now an explicit verification item for L3-2 rather than a premise | ADR-75b §11 open question added |

## 8. Reproduce

```bash
cd Ladruno_files/testbed/perf/lane3
./sweep_l3a.sh 3                      # all lanes, 3 interleaved rounds, parses to JSON
L3_LANES="B" ./sweep_l3a.sh 5          # Lane B solver sweep only, 5 rounds
```

Artifacts per configuration: `l3a_<tag>_r<N>.{h5,json,log}`. Single-file inspection:

```bash
PYTHONPATH='C:\Users\nmb\venv\opensees_env\Lib\site-packages' \
  /c/Users/nmb/venv/opensees_env/Scripts/python.exe -S parse_lane3.py l3a_laneB_pard4_r1.h5
```

Two interpreters on purpose: `python3.12` runs OpenSees (no `h5py`); `opensees_env` (py3.12 +
h5py 3.15.1) parses the h5. Do not "simplify" that to one.
