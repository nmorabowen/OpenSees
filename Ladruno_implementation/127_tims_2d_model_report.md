---
title: "Why the Ring Walls — the fork's reply to the TIMs 2D-model requests F18–F23"
project: Ladruno
status: living draft (updated as the pending results land; see the revision log at the end)
date: 2026-09-28
audience: TIMs project team (2D-model act)
answers: _tims_2d_model_requests_2026-09-25.md (F18–F23)
tags:
  - report
  - evidence
  - sanisand
  - integrator
  - tims
---

# Why the Ring Walls

**Reply to the TIMs 2D-model intake of 25 September 2026 (items F18–F23). Living document, started
28 September 2026.**

This is written for the people who run the strip-footing deck. It answers every item you asked for,
says what shipped and how to use it, and keeps apart four things: what is **shipped** (merged on
`ladruno`), what is **measured** (a number with a committed script behind it), what is **pending**
(an open PR or a run still going), and what is **your decision**. Anything pending is marked as
pending. Placeholders for results that have not landed yet look like this:

> **[PENDING — name of the result, where it will come from]**

| | |
|---|---|
| Fork | Ladruno / OpenSees, branch `ladruno` |
| Merged so far | #863 (WP-127), #864 (WP-132), #865 (WP-131 step 1), #866 (WP-133), #869 (WP-128), #870 (WP-136), #871 (WP-129), #872 (WP-134), #885 (WP-144/145 plans), #846 (WP-109, OpenMP on gcc) |
| Still open | #868 (WP-130, CPPM under Newton), #874 (WP-135, PDMY hang), #876 (WP-132 guide follow-up), #878 (WP-138, footing A/B) |
| Plan and findings A–D | `Ladruno_implementation/127_tims_2d_requests_plan.md` |

Every line of §2 of your intake was checked against the source before any work started
(`127_tims_2d_requests_plan.md`, preamble: all citations held on `fb1afe58b`; no commit had touched
`SRC/material/nD/UWmaterials/` since your `79e062367`). One reading was refined later, not refuted:
the "abrupt switch at 0.5 kPa" in the error norm is a continuous 1 kPa floor (§1, F18(a)).

---

## 0. The short version

1. **The wall is an integration failure, not a model failure and not your deck.** `ModifiedEuler`
   (IntScheme 1) has a set of discrete defects that, together, let the back-stress α leave the
   bounding surface in one accepted substep and then commit states the model cannot reach. We traced
   it (WP-128, #869) and confirmed it with an independent integrator written from the paper
   (WP-134, #872): **from admissible ring starts the exact equations never escape (0 escapes);
   `ModifiedEuler` as you run it escapes 25 times on the same starts.** The mechanism is in §2.
2. **Your b8 worst point (element 1950, Gauss point 3) cannot be integrated by anything.** Its α is
   6–7× outside the bounding surface. It is not a hard point; it is an inadmissible state that an
   earlier bad increment produced. The right answer is to refuse it, and the new integrator does.
3. **The error floor you asked for in F18(a) is the wrong fix.** Today's norm already has a 1 kPa
   floor, and at a few kPa the substep count is set by stability, not accuracy. The flag exists
   (`-errFloor`, in SAS-ME), but it will not buy what you hoped.
4. **What shipped:** a new integrator, **SAS-ME (`IntScheme 129`)**, that reproduces the paper
   equations to 2e-7 relative and refuses what it cannot integrate (WP-129, #871); per-point
   post-mortem counters and a material-point replay command (WP-127, #863); a deterministic PARDISO
   mode (WP-132, #864); PDMY03's critical-state constants as flags, plus a verdict on the PDMY
   "dilation brake" (WP-133, #866).
5. **A caution about results you already have.** Every `ModifiedEuler` SANISAND result carries an
   integration error of **6–15 % of the stress increment on 1e-4 strain increments**, measured on
   benign 20–100 kPa states (WP-129). On your footing that has so far moved the load–settlement curve
   by less than about 0.6 % past s/B 0.002, but by up to 5.5 % in the first steps (§4, interim).
6. **Still pending:** CPPM under a global Newton, with a sign fix to its tangent (WP-130, #868);
   threading SANISAND (F19 step 2); the full footing A/B and its verdict (WP-138, #878).
7. **Yours to decide** (§5): the p'-floor rule, the low-p dilatancy sigmoid (`D_factor`), the mesh pair,
   and what "limit load" means for a dense dilatant sand.

---

## 1. Item by item

Status key: **SHIPPED** (merged on `ladruno`), **PENDING** (open PR or run), **THEIR DECISION**,
**NOT DONE** (with the reason).

### F18(a) — an error norm with an absolute floor

**Asked.** `err = ‖dσ₂ − dσ₁‖ / max(2‖σ‖, σ_ref)` as a flag, and substeps and error against a
reference for σ_ref ∈ {0, 0.1, 1, 5} kPa, and where the admitted error is below the Newton tolerance.

**Answer.** The norm you proposed is what `ModifiedEuler` already computes, with σ_ref = 1 kPa. The
"switch at 0.5 kPa" is continuous: `‖dσ₂−dσ₁‖` below ‖σ‖ = 0.5 and `/(2‖σ‖)` above it is exactly
`/max(2‖σ‖, 1)`. A Python port with σ_ref = 1 reproduced today's C++ bit for bit on 640/640 ring
increments. Because ‖σ‖ ≥ √3·p, a floor of 0.1 or 1 kPa only acts below p ≈ 0.29 kPa, and 5 kPa
only below p ≈ 1.4 kPa. To touch your ring at p' ≈ 3.5 kPa, σ_ref would have to exceed about
12 kPa.

More important, **the cost is stability-limited, not accuracy-limited**. At constant p = 2 kPa and
η/M^b = 1.00 the substep count is 60 / 54 / 52 at TolE 1e-6 / 1e-5 / 1e-4. An accuracy-limited Heun
scheme would change 10× over that range; this changes 1.15×. A 20 kPa floor saves 2 of 52 substeps
there and doubles the error. The floor is inert where the error is small and cannot reach the tail
where it is large.

Where the admitted error sits below your Newton tolerance: on every constant-p and active-path state
at p ≥ 2 kPa, today and at σ_ref = 20 (errors 1e-4 to 2e-3 kPa per 1e-5 increment). It fails on the
ring tail (p95 2.9e-2 kPa, max 0.84 kPa), and those tail numbers are understated about 2–3× because
the reference used there was itself α-blind. The conversion from your `NormUnbalance` tolerance to a
stress error at one Gauss point (≈ 2e-3 kPa if the reference load is the footing weight, ≈ 0.12 kPa
if it is the bearing load) is an order-of-magnitude element estimate: **which vector you call "the
reference load" decides it, and only you can say.**

**Status.** SHIPPED as a flag of SAS-ME only; NOT DONE in `ModifiedEuler` (it would change nothing
measurable there).

**How to use it.** `-errFloor $sigRef` on an `IntScheme 129` material; default `P_atm/101` (1 kPa at
`P_atm = 101`), i.e. exactly `ModifiedEuler`'s implicit floor. Leave it at the default.

**Evidence.** `128_sanisand_ring_trace.md` §5 (tables §5.2, §5.3, §5.5); quirks row "The
ModifiedEuler error norm already HAS a 1 kPa floor…" (WP-128, #869); guide
`LadrunoSANISAND_implex_guide.md` §13.1.

### F18(b) — rate-form stages instead of two 6×6 tangents per substep

**Asked.** Compute `dσ = C:dε − Λ·C:m` inside the substep loop, form the tangent once at the end,
match today's results to round-off.

**Answer.** Not done in `ModifiedEuler`, for two reasons. First, it cannot be round-off neutral for
`TanType 2`: the chained "consistent" tangent consumes the per-stage 6×6s
(`127_tims_2d_requests_plan.md` finding D), and that chain has its own defect (it accumulates `T`
where the recurrence needs `dT`; quirks row, WP-129). Second, once the defects in §2 were found,
restructuring `ModifiedEuler` for speed would have made a wrong answer faster. SAS-ME is the
replacement: its stages compute increments at their own state and it forms **one continuum tangent
at the end state** for `TanType 1` and `2`.

**Status.** NOT DONE in `ModifiedEuler`; superseded by SAS-ME (SHIPPED). The allocation-free kernel
that would make the per-substep cost small is on the roadmap (§6).

**Evidence.** Plan finding D; quirks row "`ModifiedEuler`'s `TanType 2` 'consistent' tangent chain
accumulates `T` where the recurrence needs `dT`"; `LEDGER_implementations.md` WP-129 row.

### F18(c) — make `IntScheme 2` (CPPM) usable under a global Newton

**Asked.** Refuse at once instead of 2⁹ recursive halvings, a line search or better start, remove the
static work arrays, rerun F12's bearing deck.

**Answer, and a finding you did not ask for.** **The CPPM's `TanType 2` tangent had the wrong sign in
vanilla `ManzariDafalias`.** `NewtonSol` ends `Cep = -1.0 * CSigma`; the algorithmic tangent is
`+CSigma`. The local return is correct; only the matrix handed to the element is negated, so the
global Newton **diverges from its first iteration** and only the relaxed Krylov rung ever commits a
step. Checked against a finite difference of the return map: the vanilla sign is off by 2.0 relative;
the flipped sign by 1.24e-3 (one local iterate of staleness). That, more than the 2⁹ ladder, is why
F12 found scheme 2 475× shallower. F12 had read the code and called it a genuine algorithmic
tangent; nobody had compared it with a finite difference.

**A sign fix is not a consistent tangent.** Three error sources remain: it is one local iterate
stale (up to 0.27–0.53 relative at the default TolR 1e-7, because the local norm mixes strain and
stress units), the void-ratio dependence is missing from dR/dε (1e-4 to 1e-3), and after a halving
the second half-increment's tangent is handed out. With the fix, the global Newton converges in
about 3 iterations per step but is **superlinear, not quadratic** (median observed order 1.14–1.24).

Measured on F12's bearing deck (same leg, 1200 s budget, driver unchanged), with the recommended
recipe below, run back to back with an `IntScheme 1` control on the same loaded machine:

| arm | s/B at 300 / 600 / 900 / 1200 s | global iterations per committed step | load–settlement vs IntScheme 1 |
|---|---|---|---|
| recommended CPPM recipe | 0.00293 / 0.00421 / 0.00523 / 0.00626 | 3.2 (max 6) | 1.43 / 0.66 / 0.20 / 1.03 % at s/B 0.001 / 0.002 / 0.004 / 0.006 |
| `IntScheme 1`, same load | 0.00138 / 0.00250 / 0.00442 / 0.00698 | 16.8 (max 37) | — |

The recipe is ahead of `IntScheme 1` for the first 900 s (2.1×, 1.7×, 1.2×) and behind at 1200 s,
as refusals cut steps and the driver's 80-subdivision budget ran out (76/80 spent). Vanilla
`IntScheme 2` on the same deck: 3 steps, s/B 0.00002.

**Status.** **PENDING** — PR #868 (WP-130), draft, not merged. The static-array item is done there as
groundwork for F19 (the live shared state was `Matrix::Invert`'s scratch, now a stack-local LU;
`NewtonIter`'s statics are dead code).

**How to use it (once #868 merges).** Recommended recipe for `IntScheme 2` under a global Newton:

```tcl
nDMaterial LadrunoSANISAND $tag <18 params> 2 2 $JacoType $TolF $TolR \
    -cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on
```

`-cppmTangent fixed` is the **default on `LadrunoSANISAND`** (owner decision); `-cppmTangent vanilla`
reproduces the old binary bit for bit. Vanilla `nDMaterial ManzariDafalias` keeps the wrong sign.
`-cppmStart explicit` is **not** in the recipe: on the review's oracle set it doubled the error on 74
of 171 increments. Use a forwarding element (your `LadrunoQuad` is one).

**Evidence.** PR #868 body (tables "F18(c) refusal timing" and "F12 bearing deck, RECOMMENDED
recipe"); `origin/wp/130-sanisand-cppm-under-newton`: guide §9 "IntScheme 2 under a global Newton",
`Ladruno_files/testbed/hypo_bearing/wp130_f18c/tables_recipe.md`.

### F18(d) — per-point fallback from `ModifiedEuler` to CPPM

**Asked.** When `ModifiedEuler` hits `-maxSubsteps`, hand that point's increment to CPPM; refuse only
if CPPM also fails; one-element test.

**Answer.** Built as `-meFallback cppm` (needs `IntScheme 1` and `-maxSubsteps > 0`). One-element
test: a leg that `-maxSubsteps 20` refuses at step 1 runs all 10 steps with the fallback, 20 of 20
capped updates returned by CPPM, stress within 1.3 % of the uncapped integration; where CPPM also
fails the update is refused and nothing is integrated explicitly.

On your footing, though, this lever is small: the cost is **not concentrated in the ring** (§4: the
ring holds about 6–8 % of the substeps), so a per-point fallback would recover under 10 % of the
material time. It also inherits `ModifiedEuler`'s defects for every point it does not rescue. It is
not ported to SAS-ME yet (§6).

**Status.** **PENDING** — PR #868.

**Evidence.** PR #868 body ("F18(d) one-element fallback"); WP-138 census (§4).

### F18(e) — can any integrator take the b8 ring point?

**Asked.** Say plainly if the attached b8 point (p' 0.352 kPa, η 12.87) is one no integrator can
take, and trace how a committed η/M^b ≈ 6 arises.

**Answer. No integrator can take it, and none should.** The dumped state is **inadmissible**: its α
lies 6.26× past the model's bounding surface measured with the Lode angle of n (WP-128), 7.3× with
α's own Lode angle (WP-134), with `b:n = −8.19`. Driven by small probes:

- loading-type probes are "taken" only in the sense that the stress rides the cone around the bad α;
  the state stays inadmissible;
- unloading-type probes either commit **f > 0 as success** (`ModifiedEuler`, at both tolerances), or
  **teleport η from 12.9 to 1.33** in one increment through the force-accept clamp to `Mc` (tight
  `ModifiedEuler`, CPPM's fallback, RK45 even on a 1e-7 increment). A jump of that size in one
  increment is not an integration;
- the exact reference integrates gp 2 and gp 3 honestly under compression (α stays far outside,
  f ≈ 0) and **stops** on gp 3 under shear, where the rate equations are singular (0/0 in the loading
  index) and says so.

How the state arises is §2. In one line: after the low-p clamp sets α = 0, a reversal sets α_in = α = 0;
on the next compression increment h is the 1e10 sentinel, the first Heun stage moves α by about
Δs/p evaluated at the floor pressure while the increment raises p about 20×, and the stress-only
error test passes it at `dT = 1`. **Smallest reproducer:** σ = 0.0101·I, α = α_in = z = 0, one
plane-strain dε_yy = +1e-4. `ModifiedEuler` returns η 10.79 and α at 5.14× the bounding surface in one
substep with rc = 0; the exact equations give η 0.531–0.534 and 0.25×. The ring carries the
signature: all four b8 rows with α outside the bounding surface have α_in ≡ 0 exactly, and no row with
α_in ≠ 0 is outside.

SAS-ME **refuses** both b8 1950 points on entry (`startAlphaOutsideBounding`) and integrates the other
78 rows. That is the behaviour we recommend: refuse with a named code so your step controller cuts the
step, never project α back (a projection would silently rewrite history and hide the upstream defect).

**Status.** Answered (WP-128, WP-134). The refusal is SHIPPED in SAS-ME (WP-129).

**Evidence.** `128_sanisand_ring_trace.md` §0, §2, §4; `134_sanisand_reference_integrator.md` §0.4–0.6,
§6.4–6.6; guide §13.3 (ring row).

### F19 — SANISAND in the threaded state-determination loop

**Asked.** First the inventory of shared mutable state, then per-instance / thread-local / locked
state, then identity and speed-up at 1/2/4/8 threads.

**Answer.** The inventory is done. **The prime suspect for the `IntScheme 1` segfault is not a data
race: it is a print.** Under openseespy, `opserr` goes through `PythonStream` into CPython
(`PySys_FormatStderr`). An OpenMP worker that prints (the `-maxSubsteps` cap warning was the site on
WP-107's crashing run) calls into CPython with no Python thread state. This explains every row of
WP-107's evidence, including "a mutex and `omp critical` still crash", which rules out every
data-race explanation. The fix is a deferred per-thread message buffer flushed in element order after
the loop (which also makes the printed warnings identical at any thread count). The inventory also
found real races for step 2: the process-wide `LadrunoImplexGlobals` counters (on the plain path too,
not only under `-implex`), the warning budgets, `Matrix::Invert`'s shared scratch on the CPPM path, and
class-static return buffers in the plane-strain wrappers.

**Status.** Step 1 (inventory) SHIPPED (#865). Step 2 (code) **PENDING**; it waits on WP-130 (for
`IntScheme 2`) and on the deferred message path. Until then SANISAND stays refused from the threaded
loop.

**Build note for Esmeralda.** PR #846 (WP-109) merged on 2026-09-28 and flipped the CMake default
`LADRUNO_OPENMP` to ON, including gcc builds: until then a bare-cmake Linux build compiled the loop out.
An Esmeralda build from `ladruno` at or after `eeb7847d4` has the threaded loop; SANISAND still runs
serially until step 2.

**Expected gain** (inference, not measured): material update is 85.6 % of your wall time (your §1.2),
so Amdahl bounds 8 threads at about 3.9×; dynamic scheduling and an allocation-free kernel are needed
to approach it.

**Evidence.** `131_sanisand_threaded_inventory.md` §0–§6; PR #846; your intake §1.2.

### F20(a) — a cumulative per-point substep counter

**Answer.** `substepStats`: 17 columns **per integration point** (none process-wide), cumulative since
`revertToStart`, **not** reset by `revertToLastCommit`, carried by `getCopy` and the wire. So a
post-mortem after a failed `analyze` reads the real history instead of the zero you saw. It also counts
what used to be invisible: substeps **force-accepted at `dT_min` after failing the error test**, how
many of those fired the clamp to `Mc`, low-p abandons (the integrator returning at T < 1, silently),
and cap hits. Reading it never changes a number.

```python
s = ops.eleResponse(ele, "material", ip, "substepStats")
substeps, forced, abandoned, caps = s[2], s[5], s[8], s[9]
```

Under SAS-ME the census is `sasStats` (per point, since `revertToStart`, including the last refusal
code). WP-130 (pending) extends `substepStats` to 28 columns for the CPPM.

**Status.** SHIPPED (#863; `sasStats` #871).

**Evidence.** Guide §6.2 (column table); `tests/test_ladruno_sanisand_replay_counters.py`.

### F20(b) — profile scopes inside the integration

**Answer.** Added inside SAS-ME: `sanisand.sasME.predictor`, `.stateDependent`, `.stageArithmetic`,
`.stages`, `.drift`, `.alphaCheck`, `.substeps`, `.update`, `.tangent`. Measured split: about 0.44 ms
per update at a ring state and 0.025 ms at a deep state; at a ring state the stages take roughly half
to 60 % (the state-dependent quantities about 15 % of the total), drift correction about 10 %, the α
check 6–10 %; at a deep state the tangent and drift are about 11 % each.

Not added inside `ModifiedEuler` or the CPPM Newton: `ModifiedEuler` is the integrator we recommend you
leave, and its scopes would not survive the byte-identity constraint cheaply.

**Status.** SHIPPED for SAS-ME (#871); NOT DONE for `ModifiedEuler`/CPPM.

**Evidence.** PR #871 body ("F20(b) profile split"); guide §13.3; `Ladruno_files/testbed/wp129_sasme/out/`.

### F20(c) — a `"tangentEP"` response

**Answer.** `tangentEP` returns the 6×6 continuum elastoplastic tangent at the **committed** state,
whatever `TanType` the deck uses. Checked against a one-sided finite difference at a plastic state
(50 kPa, three strain directions): relative difference below 1e-4.

```python
C = ops.eleResponse(ele, "material", ip, "tangentEP")   # 36 values, row-major
```

**Status.** SHIPPED (#871).

**Evidence.** `tests/test_ladruno_sanisand_sasme.py::test_tangentEP_matches_finite_difference`.

### F21 — replay a dumped material state

**Answer.** `ladrunoSANISANDReplay` puts a private copy of a `LadrunoSANISAND` prototype into a given
(σ, α, α_in, z, e) and drives one strain increment through the same `setTrialStrain` an element uses.
It returns rc, the census of that one update, the returned state (σ, α, α_in, z, e, p, q, f before and
after, the path code) and a per-substep trace (`T, dT, err`, outcome code). Replaying step k+1 from the
committed state of step k reproduces an analysis step to 1e-9.

**Finding A — your CSVs are compression-positive.** The README says "compression negative as
OpenSees stores it", but on all 80 rows `p_kPa = +tr(σ)/3` and every normal stress is ≥ 0: the columns
are the model's internal `mSigma`. So the replay has **no default convention**; you must say which:

```tcl
ladrunoSANISANDReplay $matTag -convention compressionPositive \
    -sigma s11 s22 s33 s12 s23 s31 -alpha ... -alphaIn ... -fabric ... \
    -voidRatio $e -dStrain d11 d22 d33 g12 g23 g31 <-type 3D|PlaneStrain> <-trace 10000>
```

Shear strain is engineering (γ). α, α_in and z are projected to their deviatoric parts with a warning
(b8 row 1859/2 has tr α = 2.3e-3). A Python helper reads your CSVs and runs the standard probes:
`Ladruno_scripts/sanisand_replay.py` (`replay`, `read_ring_csv`, `probes`). For your own dumps from
now on, read `substepStats` in the same dump.

**Status.** SHIPPED (#863).

**Evidence.** Guide §6.3; quirks row "The TIMs ring-point CSVs carry the INTERNAL,
compression-POSITIVE `mSigma`"; `test_replay_reproduces_an_analysis_step`.

### F22 — a deterministic mode

**Asked.** MKL conditional numerical reproducibility for PARDISO, a list of what else is
order-dependent, byte-identical curves twice on 8 threads.

**Answer.**

```tcl
system Pardiso -deterministic              ;# MKL CNR on the AUTO branch + iparm(34)
system Pardiso -cbwr COMPATIBLE            ;# an explicit branch every x86 node can run
```

The first solve prints what MKL actually has in force, e.g.
`PARDISO deterministic mode: MKL CNR branch AUTO, iparm(34)=8 thread(s), CNR ACTIVE`. Measured on a
~22k-DOF push at 8 MKL threads, 5 runs each: mode on, 1 distinct result; mode off, 5 distinct
displacement fields.

Four things to know:

- **Across nodes with different CPUs** (your §1.6 case), AUTO picks a code path per CPU. Pin a branch
  every node can run: **`-cbwr COMPATIBLE`**. The instruction-set branches (`AVX2`, `AVX512`, …) exist
  only on Intel CPUs; on an AMD machine every one of them was refused and only `AUTO` and `COMPATIBLE`
  worked. The thread count must match too.
- **The mode is process-wide and sticky**: it stays on for every later model in the same interpreter.
  MKL refuses to set it once its BLAS/LAPACK dispatch has started (an `eigen` before the `system` line
  triggers this; an earlier PARDISO solve does not). The reliable route is the `MKL_CBWR` environment
  variable set before the process starts.
- **Serial targets only.** `OpenSees.exe` and the sequential `opensees.pyd`; MUMPS and MPI reductions
  are not covered.
- **What else is order-dependent**: the threaded element loop reduces only an integer and is
  bit-identical at 1/2/4/8 threads; SANISAND is refused from it; the `-implex` counters are process-wide
  but serial. See the table in the PARDISO recipe.

**Repeatable is not reliable.** A deterministic mode makes two runs agree; it does not make either run
more correct. Every threaded run is equally correct to machine precision. When a last-bit difference
grows into a 30 % shift in where the wall sits, the model is on a knife edge (a limit point, a yield
state that can flip, a Newton that converges right at its tolerance, an adaptive cut that can go
either way), and a different tolerance, step size or mesh would move it too. Your deck has this
character for a measured reason: with `ModifiedEuler` the stress–strain map is non-smooth (the err = 0
path, §2), so there may be no equilibrium for Newton to converge to. On the fork's own
flip-determinism deck, the first push step has **no reachable equilibrium under any tangent**; a
10⁴× tighter TolR only halves the Newton residual floor (0.20–0.35 kN → 0.12–0.13 kN), and the
"converged" first-step load under `NormDispIncr` moves by about 25 % (9.66 → 7.24) (WP-136). Use
`-deterministic` for regression tests, for reproducing a failure, and for comparing nodes; do not use
it to settle a result.

**Status.** Mode SHIPPED (#864). The "repeatable is not reliable" guide paragraph is **PENDING**
(#876). Measured on Windows (AMD); not yet run on Esmeralda.

**Evidence.** `75c_pardiso_solver_recipe.md` Trap 7, "The deterministic mode"; quirks rows WP-132 (CNR
process-wide and sticky; `mkl_cbwr_set` returning -8; `-cbwr AVX2` refused on AMD; `iparm` zeroed at
every symbolic phase); `tests/test_wp132_deterministic_pardiso.py`; `136_flip_test_drift.md`.

### F23(a) — PDMY03's critical-state constants

**Answer.**

```tcl
nDMaterial PressureDependMultiYield03 $tag ... <-ei $e0> <-cs1 $v> <-cs2 $v> <-cs3 $v>
```

Flags, after every positional argument, in any order; defaults 0.6 / 0.9 / 0.02 / 0.7 (the former
hard-coded values), byte-identical when omitted (Python and Tcl baselines captured before any edit).
Found along the way and fixed: the per-material reallocation every 20 materials overwrote every
existing material's constants with the newest one's. Harmless while they were hard-coded; a silent
cross-material leak once they are user-set.

Also found, **not fixed**: `pAtm` is a static member of PDMY01/02/03, so the last material created sets
the atmospheric pressure for every material of that class. Keep one `$pa` per class per process.

**Status.** SHIPPED (#866).

### F23(b) — the PDMY "dilation brake"

**Answer.** Your reading is right that the brake is keyed to void ratio and that a dense sand never
reaches it with the default constants (from e = 0.6 it must dilate 17.5 % volumetrically at 100 kPa,
9.9 % at 1 652 kPa). It is incomplete in a way that matters: **reaching it would not help**.
`isCriticalState()` is a **crossing detector**. It returns 1 only for the increment whose start and
end lie on opposite sides of the line; past the line both are on the same side again and the full
dilatancy rule resumes (measured: the volumetric rate dips at the crossing step and is back within
0.5 % ten steps later). So retuning `ei`/`cs1..3`, now possible on PDMY03 too, moves *when* one
increment loses its dilatancy; **no choice of constants yields a plateau**. That is consistent with
your candidates with a retuned line failing the saturation gate as well. The route to a plateau is a
model whose dilatancy vanishes at critical state by construction (SANISAND's D ∝ M^d(ψ) − η, PM4Sand).
Your ten-candidate and strip numbers were not re-run.

A related defect found in WP-133 and fixed in the pending WP-135: a wild Newton iterate makes PDMY's
substep count `|Δε|/1e-5` explode to about 1e9 per call, which is why a two-element model "hung" in
`analyze`. With the fix PDMY refuses such a trial in milliseconds. Under `SSPquad` (a host that
discards the refusal) the call is bounded but the step can still be accepted.

**Status.** Note SHIPPED (#866). Hang fix **PENDING** (#874).

**Evidence.** `133_pdmy_notes.md` (b); quirks row "PDMY's 'dilation brake' `isCriticalState()` fires
only on the increment that CROSSES the critical-state line"; PR #874.

---

## 2. What actually walls the deck

### 2.1 The ring, briefly

Just outside the footing edge, the top row of Gauss points sits at p' of a few kPa with η on the
bounding surface. There SANISAND's plastic modulus scales with p and its elastic moduli with √p, so the
rate equations are stiff, and an explicit scheme takes substeps sized by stability, not accuracy. That
part is physics and would cost time under any explicit integrator. It is not what stops the run.

### 2.2 What stops it: a chain of discrete defects in `ModifiedEuler`

Each link is a quirks row; the ranking is the one the independent reference integrator (WP-134)
established, which corrected WP-128's first ranking of F.

| role | mechanism | what it does |
|---|---|---|
| **trigger** | **G** — α_in is re-seated once per increment | inside the substeps (α − α_in):n reaches 0, so h is the 1e10 sentinel, and then goes negative, so h < 0 and the α law becomes a repelling relaxation. 37 of 38 crossing substeps have it. The paper resets α_in at the start of each new loading process, which makes h < 0 impossible (0 of 960 runs of the exact reference). |
| **enabler** | **E** — the substep error looks at stress only | a substep that throws α 5–16× outside the bounding surface passes, because both Heun stages have the same stress increment. Adding α to the error alone keeps α inside. |
| **enabler** | **F** — a loading stage with a negative denominator is taken as elastic, and the step factor has no upper cap | both stages then agree exactly, **the error is exactly 0**, and the next substep swallows the rest of the increment. No tolerance can see it. It accounts for all 25 escapes of your `ModifiedEuler` from admissible ring starts, and for 20–65 % stress errors on benign 20–100 kPa states even at TolE 1e-8. |
| adds error | **U9** — K and G are frozen at the committed state for the whole increment | 0.6 / 6 / 24 % of the stress increment at δ = 1e-5 / 1e-4 / 1e-3. Both stages share the same wrong moduli, so the error test cannot see it and a tighter TolR does not shrink it. |
| adds error | **U10** — the loading test uses n:Δσ, not the yield-function gradient | it ignores the −(n:r)dp term, so an isotropic compression that lowers η can be read as plastic. |
| commits it | **C** — at `dT_min` a substep that failed the error test is accepted anyway, with a clamp to `Mc` | the η 12.9 → 1.33 "teleport". Uncounted until WP-127. |
| commits it | `Stress_Correction` gives up silently | when neither correction direction reduces f, it returns the uncorrected state with f > 0 and rc = 0. Worst measured: f = 11.2 kPa at p = 0.58 kPa. |

The error estimate is the reason all of this stayed hidden. It measured only stress, it read zero on
the err = 0 path, and it compared two stages that shared the same frozen moduli. So every failure mode
above produced an increment the estimator called accurate. Tightening TolR, flooring the norm, or
raising `-maxSubsteps` all act on that estimator, which is why none of them moved your wall.

### 2.3 What SAS-ME does instead

SAS-ME (`IntScheme 129`) is a Sloan–Abbo–Sheng-style explicit modified Euler written against the
oracle, not a patch of `ModifiedEuler`:

- exact elastic path (closed form in √p) for the predictor and the intersection;
- every Heun stage evaluates K, G and every state-dependent quantity at its own state (U9);
- stages classified from the true yield gradient; a loading stage with H ≤ 0 is refused or cut, never
  called elastic (F, U10);
- the error covers σ, α **and** z; TolR is always honoured; the step factor is capped at 1.1 with no
  growth after a rejection (E, F);
- the paper's α_in rule inside the increment (G);
- refusal instead of force-accept, with named codes; its own drift correction fails rather than
  returning f > TolF (C);
- a bound check on α after every substep, and a refusal of inadmissible starts.

Measured against the oracle: benign 20–100 kPa states within 2e-7 relative at TolR 1e-7 (5e-5 at
TolR 1e-4), where `ModifiedEuler` is 6–15 % off on 1e-4 increments; the smallest reproducer at
ρ 0.252, η 0.531 in 251 substeps (oracle 0.252 / 0.531; `ModifiedEuler` 5.14 in 1 substep); the ring,
624 of 640 increments integrated, the 16 from b8 1950/2–3 refused, max f at exit 1e-7, no escape. It
costs more per increment: median 17 substeps on the ring against 4 for `ModifiedEuler`, and 4–6 per
1e-5 increment on smooth monotonic chains against 1–3. That is the honest cost the α-blind test was
hiding.

**How to use it.**

```tcl
nDMaterial LadrunoSANISAND $tag $G0 $nu $e_init $Mc $c $lambda_c $e0 $ksi $P_atm $m $h0 $ch $nb \
    $A0 $nd $z_max $cz $Rho  129 $TanType $JacoType $TolF $TolR \
    <-errFloor 1.0> <-alphaBoundTol 0.1> <-alphaEntryTol 2> <-alphaProject 0> \
    <-sasAlphaIn reseat> <-sasErrorVars full> <-maxSubsteps $n> <-Pmin ...> <-Presidual ...>
```

- `TolR` **is** the substep tolerance. Recommended 1e-4 to 1e-7; default 1e-7. Below about 1e-8, large
  low-p increments cannot meet it above `dT_min` and the update refuses. `-honorTolR` is inert (warned).
- The defaults shown are the shipped defaults. `-sasAlphaIn stale` and `-sasErrorVars stress`
  reproduce `ModifiedEuler`'s defects G and E, for attribution only. `-alphaProject 1` projects α
  instead of refusing; it rewrites history and is off by default.
- `-implex` is refused with 129.
- `TanType 1` and `2` both return the continuum tangent at the end state; `0` returns Ce.
- Refusal codes (in `sasStats` and the warning): 1 startOutsideYield, 2 startAlphaOutsideBounding,
  3 startInadmissible, 4 errorAtDTmin, 5 loadingNonPosH, 6 tensionAtDTmin, 7 driftFailed,
  8 alphaOutsideAtDTmin, 9 maxSubsteps.
- Use a **forwarding** element. `LadrunoQuad` (your element) forwards the refusal, so your step
  controller cuts the step.

The configuration we ran on our copy of your deck (WP-138) was
`129 0 1 1e-7 1e-7 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0`, on
`LadrunoQuad -bbar` at B/8. That is a measured configuration, not yet a recommendation: see §4.

**Behaviour change you should know about.** Since #871, `LadrunoSANISAND::commitState` refuses to
commit a trial whose last update was refused, and this includes a **`ModifiedEuler` `-maxSubsteps` cap
hit**. Under a forwarding element nothing changes (the step already failed). Under a **discarding**
element (`SSPquad`, `stdBrick`, `BbarBrick`, …) such a deck used to commit the strain without the
stress, silently; it now fails the step and the point latches until `revertToStart`. Also, database
and restart files written by an older build will not load (the wire vector grew).

---

## 3. Cautions for results you already have

1. **Every `ModifiedEuler` SANISAND result carries an integration error of this size.** Per increment:
   6–15 % of the stress increment on 1e-4 strain increments at benign 20–100 kPa states, 20–65 % on the
   err = 0 path, and U9 alone at 0.6 / 6 / 24 % for δ = 1e-5 / 1e-4 / 1e-3, none of it visible to the
   error test or shrinking with TolR. This includes your campaign curves. How much it moves a
   load–settlement curve depends on the deck; on your footing it is small past s/B 0.002 and up to
   5.5 % in the first steps (§4, interim). Treat any `ModifiedEuler` ring-point state, and any quantity
   read from ring points, as unreliable.
2. **The explicit lane's failure dumps contain inadmissible states.** The b8 1950/2–3 rows are not a
   hard point of the material; they are the product of the defects in §2. Do not calibrate or test
   anything against them except a refusal.
3. **Your §1.5 tangent comparison had a cause.** Under `IntScheme 1`, `TanType 2`'s chained tangent
   accumulates `T` where it needs `dT`, and the stress–strain map is non-smooth (the err = 0 path).
   `TanType 0` was the only dependable choice under `ModifiedEuler` for that reason. Under SAS-ME,
   `TanType 1` and `2` are the continuum tangent. Under CPPM, the vanilla `TanType 2` had the wrong sign
   (fixed by default in #868, pending).
4. **Committed steps on the relaxed rung.** Your ladder's last rung (`KrylovNewton` at 10× the
   tolerance) is not a small print item on this deck. On our copy, 47 of 79 committed `ModifiedEuler`
   steps and 26 of 61 SAS-ME steps were committed there, and a `TanType 1` run committed almost every
   step there once the cap started refusing off-path iterates. That acceptance alone moved the curve by
   about 1–2 % in our interim runs. Report the rung of every committed step, and the share committed at
   the relaxed tolerance.
5. **The 30 % run-to-run shift is a signal about the deck, not the solver** (F22). Deterministic mode
   will make the two runs agree; it will not tell you which is right.
6. **The `-Presidual` 1.01 / 5.05 kPa comparison of your §1.4 was made with the defective integrator.**
   Its conclusion ("does not move the wall") should be re-measured under SAS-ME before it is used to
   justify a floor (§5).
7. **PDMY under `SSPquad`**: a refused trial is discarded by the host (WP-135 bounds the time, not the
   acceptance). Prefer a forwarding element (`quad`, `LadrunoQuad`, the u-p family) for any deck that
   depends on a material refusal.

---

## 4. The footing A/B (WP-138) — INTERIM

> **Status: interim. PR #878 is a draft and the runs are not finished. Numbers below come from the
> committed run records on `origin/wp/138-footing-sas-me-ab`
> (`Ladruno_files/testbed/footing_sas_me_ab/runs/`) and will be replaced by the WP-138 report
> `138_footing_sas_me_ab.md` when it lands.**

**Setup.** The fork's own copy of your deck, built from the §1 spec of your intake (nothing in the
Workbench was run): plane-strain strip, `LadrunoQuad -bbar`, B/8, `system Pardiso`, `NormUnbalance`
1e-5, Newton → NewtonLineSearch → KrylovNewton (tol × 10), MKL 1 thread with `MKL_CBWR=COMPATIBLE`.
Arm A: `ModifiedEuler` (`IntScheme 1`, TanType 0). Arm B: SAS-ME (`IntScheme 129`, TanType 0, TolR 1e-7,
build `beb6d8333`). A `DruckerPrager` 38° control reached s/B 0.15 in 249 s, as on your side.

**Where the runs are.**

| arm | reached | how it ended | committed steps (rung N / NLS / K×10) |
|---|---|---|---|
| A, `ModifiedEuler` | s/B 0.0174, q 441.6 kPa | the 36 000 s wall-clock cap, **not** the step floor (the machine paused about 2.75 h inside that budget) | 79 (23 / 9 / 47) |
| B, SAS-ME | s/B 0.0142, q 379.8 kPa at the last committed record | running / not concluded | 61 (24 / 11 / 26) |

**No wall yet** on either arm.

**Load–settlement agreement** (same settlement, matched steps):

| s/B | 1.3e-5 | 1e-4 | 5e-4 | 1e-3 | 2e-3 | 0.0038–0.0046 | 0.005–0.0139 | 0.0140–0.0142 |
|---|---|---|---|---|---|---|---|---|
| \|q_A − q_B\| / q_A | 5.5 % | 3.6 % | 1.4 % | 0.72 % | 0.22 % | ≈ 0.59 % | ≤ 0.31 % | 0.37–0.43 % |

SAS-ME carries less load in the first steps (ModifiedEuler's per-increment error is largest there,
relative to the load). Past s/B 0.002 the two curves agree to about 0.6 %.

**Where the cost is.** Over SAS-ME steps 25–50, the ring points hold **7.7 % of all substeps** (3.1 %
by the stricter criterion), and the 100 costliest points 6.8 %. The cost is spread over the whole
plastic zone, not concentrated in the ring. So a per-point fallback (F18(d)) would recover less than
10 % of the material time, and the levers that matter are the ones that act everywhere: fewer global
iterations, and a cheaper substep (§6).

**Replays at `ModifiedEuler`'s last converged step**, the 44 worst points (highest ρ and the low-p
ones), each replayed on its real next increment against the oracle: both integrators land at about
1e-4 of p' from the oracle (medians: `ModifiedEuler` ≈ 7e-5, SAS-ME ≈ 1e-4), except one low-p point that
`ModifiedEuler` took in a single substep, 2.0e-3 of p' off (SAS-ME 22 substeps, 9e-7).

> **[PENDING — WP-138 final: reconcile the replay statistics.** The orchestrator's summary quotes
> "ModifiedEuler 1.2e-3·p' vs SAS-ME 6e-5·p'" on the footing's real increments; the committed
> wall-pair table gives the medians above. Which set of increments that figure covers will be stated
> when the WP-138 report lands.**]**

**TanType 1 under SAS-ME (arm C, stopped).** It cut the global iterations early (step 22: 3 iterations
and 0.25 M substeps against 20 and 1.66 M under TanType 0). But with the 2000-substep cap, off-path
Newton iterates get refused (up to 782 cap hits per step from s/B 0.0016 on), and from then on almost
every step was committed on the KrylovNewton rung at 10× the tolerance; that acceptance cost about 1–2 %
in the curve (C against A: +1.7 % at s/B 0.0021, +1.0 % at 0.0033). The run was killed externally at
s/B 0.0040, and its restart was stopped by an owner re-plan. So TanType 1 is **not yet** the
free 5–7× lever an early-step comparison suggests.

> **[PENDING — WP-138 final verdict: does SAS-ME reach a wall, a peak or a plateau on this deck, and
> where; the q–s curves to the end of both arms; wall-clock per unit s/B.]**
>
> **[PENDING — the default decision: whether the fork recommends `IntScheme 129` as the default for
> the TIMs campaign, with which TanType, TolR, `-maxSubsteps` and step-size policy.]**
>
> **[PENDING — B/16, if WP-138 runs it.]**

---

## 5. Your decisions

These are recommendations. Nobody outside the calibration and the project can make them.

| # | decision | our recommendation | why |
|---|---|---|---|
| D1 | **The p'-floor rule** | (a) `-Pmin` ≤ 0.5 kPa; (b) `-Presidual 0` (at most 1 kPa if used, and then declared as a regularisation, not physics); (c) report the limit load at floor F and at F/2, and accept the floor if the load moves by less than about 2 %; (d) report the number of Gauss points at the floor at the limit state | precedent: PM4Sand and numgeo floor p at 0.5 kPa; HySand's authors needed a 10 kPa surcharge in their own 3D FE. An apparent cohesion p_r·tan φ times N_c ≈ 30–75 means 1 kPa of p_r can move a ~650 kPa load by 3–8 % |
| D2 | **`D_factor`, the UW low-p dilatancy sigmoid** | decide it explicitly; do not inherit it | it is not in Dafalias & Manzari (2004); it acts below p < 0.05·P_atm = 5.05 kPa, i.e. across most of the ring; with `-Presidual 0` it can suppress dilatancy by up to ~900× (the 86 report §6). **There is no deck flag to switch it off today**; if you want the on/off comparison, ask and we add one |
| D3 | **Mesh pair** | run and report B/8 and B/16 | ψ-softening localises; treat any post-peak plateau as mesh-dependent unless regularised. The earlier regularization report (the 90 report) showed the matched-settlement load converges without a regulariser |
| D4 | **The definition of "limit load" for a dense dilatant sand** | say which you mean: the peak, a plateau, or q at a fixed s/B (your ADR 65 D6) before the runs, not after | a dense sand under SANISAND may peak and soften, and the softening branch is mesh-dependent; the fixed-s/B criterion is the one that converged in the 90 report |
| D5 | **Re-running campaign curves** | re-run under SAS-ME the ones that feed a reported number, starting with anything read from ring points | §3.1 |
| D6 | **The reference load for `NormUnbalance`** | name the vector | it decides whether the integrator's per-point error sits under your Newton tolerance (F18(a)) |

---

## 6. Roadmap

**Near term (the pending items above).**

| item | what | where |
|---|---|---|
| CPPM under Newton + tangent sign | merge, then qualify on the strip | #868 (WP-130) |
| Footing A/B verdict | finish both arms; the default decision; B/16 | #878 (WP-138) |
| F19 step 2 | the deferred message buffer, the counters, stack scratch; allowlist `IntScheme 1`, then 2 and 129; identity and speed-up at 1/2/4/8 threads on a deck that prints | after #868 |
| PDMY hang | refuse a wild trial instead of ~1e9 substeps | #874 (WP-135) |
| Guide paragraph | "repeatable is not reliable" | #876 |
| SAS-ME follow-ups from its confirmation review | a 64-sample unload-then-reload locator can miss a very short elastic excursion (up to 5.5e-4 on 5 ring cases); the ψ-driven dead end is relocated to ~30 MPa, not removed; positional TolF/TolR are not range-checked | follow-up WP |

**Performance.** Material cost is (global iterations) × (cost per point per iteration), everywhere in
the plastic zone: every OpenSees iteration re-integrates the whole increment from the committed state.
So the levers, in the order we would pull them:

1. **Fewer global iterations**: a step-size policy that targets a few iterations per step, and the
   continuum tangent where it holds (§4: TanType 1 is not yet qualified on the deck).
2. **Threading** once F19 step 2 lands: bounded near 3.9× at 8 threads by your own 85.6 % share
   (inference).
3. **TolR 1e-4 → 1e-3 as a declared sensitivity**, not a silent default: SAS-ME at 1e-4 is already two
   orders of magnitude inside its accuracy headroom (5e-5 against the oracle), and accuracy-limited
   substeps scale roughly as TolR^−½. It does little in the stability-limited ring.
4. **An allocation-free kernel**: SAS-ME builds dozens of heap vectors per substep; a fixed-size stack
   kernel is typically several times faster for 6-vector arithmetic. Unmeasured; it must be
   benchmarked before a number is quoted, and it is also a precondition for threads to scale.
5. **A stiff-point fallback** SAS-ME → CPPM: a smaller lever on this deck than it first looked (§4).

GPU or SIMD batching and surrogate models are **not** levers for this problem: the per-point work is
wildly unbalanced (tens to 10⁵ substeps per point per step), and a validation study must not trade the
full-order answer for speed.

**Longer term: a model that is consistent by construction.** DM04 SANISAND has no energy function for
α, so non-negative dissipation is not guaranteed, and the α_in memory and the h ∝ 1/((α − α_in):n)
singularity are exactly where the bugs lived. Two classes are planned (tags reserved, no code yet):

- **WP-144 `LadrunoNORSAND`** — NorSand in the Borja & Andrade (2006) form, with DM04's power-law
  critical-state line (so ψ stays bounded as p' → 0) and a Lode-angle dependence for plane strain;
  hyperelastic, implicit, closed-form tangent, proven non-negative dissipation. About 5–6
  engineer-weeks. Start condition: after WP-138 shows what SAS-ME/CPPM reaches on your deck.
- **WP-145 `LadrunoHySAND`** — the 2026 hyperplastic multisurface sand model: the most rigorous
  thermodynamics and the best cyclic behaviour, the least evidence and no public code; a research
  build of about 10–14 weeks, for the later cyclic SSI work.

One expectation to set now: **no dilatant sand model is truly variational.** Non-associated
dilatancy gives a non-symmetric tangent and no incremental minimum principle. "Thermodynamically
admissible, robust implicit return, consistent tangent" is achievable; "variational" is not. And a
p'-floor or a surcharge is established practice even for the rigorous models, so D1 does not go away.

---

## 7. References

**Intake and plan**
- `Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md` + attachments
  `_tims_2d_model_requests_2026-09-25/` (README, `ring_points_b8.csv`, `ring_points_b16.csv`)
- `Ladruno_implementation/127_tims_2d_requests_plan.md` (findings A–D, work packages)

**Evidence documents**
- `128_sanisand_ring_trace.md` — WP-128, #869 (with its 2026-09-27 correction note after WP-134)
- `134_sanisand_reference_integrator.md` — WP-134, #872 (the oracle; U1–U10)
- `131_sanisand_threaded_inventory.md` — WP-131 step 1, #865
- `133_pdmy_notes.md` — WP-133, #866
- `136_flip_test_drift.md` — WP-136, #870
- `_sand_model_survey_2026-09-27.md`, `_sanisand_external_survey_2026-09-27.md`,
  `144_ladruno_norsand_plan.md`, `145_ladruno_hysand_plan.md` — #885

**Guides**
- `LadrunoSANISAND_implex_guide.md` §6.2 (`substepStats`), §6.3 (replay), §13 (choosing an IntScheme;
  SAS-ME); §9 "IntScheme 2 under a global Newton" on `origin/wp/130-sanisand-cppm-under-newton` (#868)
- `75c_pardiso_solver_recipe.md` Trap 7 and "The deterministic mode" (#864; follow-up paragraph #876)

**Ledger rows** (`LEDGER_quirks.md`): the ring CSV convention (finding A); force-accept at `dT_min`
(finding C); F (negative denominator as elastic); the uncapped step factor; G (α_in once per increment);
`ModifiedEuler` TanType-2 chain (T vs dT); RK45 `dAlpha3/4`; IntScheme 4 non-determinism; U9; U10; the
err = 0 path; stress-only error (E); `Stress_Correction`'s silent give-up; RK45 is IntScheme 45 and not a
reference; the 1 kPa floor and stability-limited cost; the flip-determinism pins; the WP-132 CNR rows;
PDMY03 constants and reallocation; the PDMY crossing detector; static `pAtm`; a refusal under a
discarding element was committed. `LEDGER_implementations.md` rows WP-127, WP-129, WP-132, WP-133,
WP-134.

**PRs.** Merged: #863, #864, #865, #866, #869, #870, #871, #872, #884 (Windows-only CI gap), #885, #846.
Open: #868, #874, #876, #878.

**Earlier replies to the TIMs team.** `86_ladruno_sanisand_tims_report.md` (the hidden cohesion),
`90_ladruno_regularization_tims_report.md` (regularization, `-maxSubsteps`).

---

## Revision log

| date | change |
|---|---|
| 2026-09-28 | First issue. F18(a), (b), (e), F20, F21, F22 (mode), F23(a), (b) answered from merged work; F18(c), (d), F19 step 2, the F22 guide paragraph and the WP-138 footing A/B pending; placeholders in §4. |
