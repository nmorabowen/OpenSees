---
title: ADR 94 R2 red lane — tests/process
project: Ladruno
status: draft
owner: nmora
tags: [implementation, material, review, redblue]
---

# ADR-94 R2 — RED (tests/process lane)

Stance: the R1 evidence is real but narrower than the hlist docs read, and the
bookkeeping ADR-87 requires for a warrant is not yet done. Measured on build
`52314165a` in this worktree; no C++ edited.

## 1. Sentinel-rule violations

No CONFIRMED-H test is provably vacuous, but **H13's two tests do not gate the
most likely minimal fix**. Both assertions are behavioural-equivalence checks
(`codes_typo == codes_no_flag`, `f_with_zero_phi <= tol`) — a fix that adds an
`opserr` WARNING on the unmatched token/parameter but still silently falls
through (the same style as every other diagnostic in this file, and the
ADR-84 precedent) changes **nothing** either assertion checks. The test only
flips red if the fix aborts construction or actually applies a fuzzy match.
Recommend adding an assertion on stderr content (child-process pattern,
already built for H4/H12) as a second gate so a warn-only fix is visible too.
No other CONFIRMED test has this gap — H1/H4/H5/H7/H8/H9 all assert a
directly-measured number or a structural absence that a real fix necessarily
changes.

## 2. Docstring naming (literal check)

Every test function *name* embeds its H (`test_H10_...`), but **0/22
docstrings spell the H-number inside the docstring body itself** (verified by
regex over all three files) — traceability from the hlist docs to pytest
node-ids works only through `-k`/node-id matching, not through docstring
grep. Minor/pedantic (the ledger rows already cite exact test names), but the
plan's own phrasing ("each test named after its H") is satisfied only by the
function name, not the docstring, so a docstring-only search (e.g. `pytest
--collect-only -q | grep H7`) works but `grep -A5 '"""' | grep H7` does not.

## 3. Order independence — measured, no leak

Ran the three files together, reverse order, and each alone (`python3.12 -m
pytest`, PYTHONPATH=`dist/bin`): **27/27 pass in every order**, and each file
alone reproduces its documented count (11, 11, 5). No shared-state leak
observed across files. This is consistent with, but does not disprove, H1's
own finding that `ASDPlasticMaterial3D`'s tangent/stress statics are shared
across *instances of one specialization* — the three files mostly use
different YF combos (VM / MC / HB+DP) so cross-file static leakage would not
necessarily show up in these particular assertions even if present.

## 4. Flake check

3 consecutive full-battery runs (all three adr94 files + all four baseline
ASDP files together): **47 passed, 2 skipped, 0 failed, all three runs**,
timings 1.77–1.85s. No `WinError 6` observed. This is a *small* sample
against a documented ~1-in-3 failure rate on the pre-fix code path — 3 runs
has roughly 30% chance of missing a real 1-in-3 flake if the `stdin=DEVNULL`
fix were ever silently reverted. Not a finding against the current fix, but
the claimed "12 consecutive runs, 0 failures" in R1-B is itself the same
class of thin evidence — nobody has run it 50+ times.

## 5. Reading-only CONFIRMED verdicts (not warrant-grade on their own)

- **H2** — UB claim, no test by design (stated explicitly, not hidden).
- **H4 core defect** (`revertToLastCommit` no-op) — pinned by source regex on
  the function body, not by observing a runtime state divergence (the
  `TenNodeTetrahedron` self-heal blocks the runtime channel entirely, per the
  doc's own caveat). Only `revertToStart()`'s silent-swallow half is runtime.
- **H5**, 2 of 6 sites (`Backward_Euler_LineSearch`, `Runge_Kutta_45_Error_Control`
  non-`_old`) — reading only; both integrators fail to converge globally on
  every rig tried, so the claim "unguarded" is untested in the one state
  (a successful commit) where it would matter.
- **H6 classification** ("cutting plane not closest point") — reading only;
  only the *accuracy* and *tangent* claims built on it are measured.
- **H9**, RK45 (non-`_old`) half — reading + "live-code proof", no runtime pin
  (same global-convergence-failure blocker as H5).
- **H10a mechanism** (the `arg=0` branch discontinuity) — reading only; the
  *drift outcome* is measured, the *cause* is asserted from source.
- **H14** — structural regex only; explicitly labeled "practically latent",
  correctly disclosed, but still a reading-only verdict for the underlying
  post-construction scenario.
- **H15** — structural (`et()` call-count in source); the intended StiffSoil
  runtime comparison NaN'd on step 1 for every parameter combination tried,
  so the "genuinely different constitutive operators" claim rests on reading,
  not on two integrators actually disagreeing at runtime.

Roughly a third of the 15 R1 verdicts lean on reading for at least one half
of their claim. That is disclosed in the docs (good practice) but the
plan's acceptance line 159 ("'read the code' is not a verdict") is not met
literally for these — R1 treats "reading + a source-regex test" as
sufficient, which is weaker than a runtime pin.

## 6. Mutation gate — H1, H5, H7, H13

- **H1**: drop `static` from `VoigtMatrix Stiffness;` (`ASDPlasticMaterial3D.h:4119`,
  + its matching out-of-class definition). Read-through: `Stiffness` becomes
  per-instance, `ComputeTangentStiffness`'s writes and `getTangent()`'s read
  are on the same object again, so `test_H1_one_static_tangent_is_shared_by_
  every_element`'s `_rel(blk_pl, blk_el) < 1e-9` assertion fails → flips red.
  Confirmed by reading; not applied.
- **H5**: at `ASDPlasticMaterial3D.h:1423` (Forward_Euler's elastic shortcut),
  insert one line before its `return 0;`: guard on `strict_convergence &&
  yf_val_end > f_absolute_tol` returning `LADRUNO_MATERIAL_REFUSED`. Blocks
  the bad commit for the `Forward_Euler` parametrized case →
  `test_H5_strict_convergence_does_not_gate_other_integrators[Forward_Euler]`'s
  `max(f_MC) > tol` assertion fails → flips red. Confirmed by reading.
- **H7**: at `ASDPlasticMaterial3D.h:2298-2302`, insert one line after
  `Stiffness = Eelastic;` guarding on `strict_convergence` to return
  `LADRUNO_MATERIAL_REFUSED` instead of `return 0;`. Flips
  `test_H7_strict_convergence_does_not_gate_the_inconsistency_branch`'s
  "flag-on bit-identical to flag-off" assertion. Confirmed by reading.
- **H13**: no defensible ONE-LINE fix exists that is guaranteed to flip both
  tests (see §1) — the minimal one-line fix (append `else { opserr <<
  "WARNING..."; }` after the if-chain / tuple base case) does **not** flip
  either test, since both check only `analyze()` codes and committed
  physics, never diagnostic output. A fix that flips them must change
  control flow (throw / return a construction failure), which is not a
  one-line change at either site (`OPS_AllASDPlasticMaterial3Ds.cpp` if-chain
  is independent `if`s, not `else if`; `utuple_storage.h`'s base case has no
  error channel to raise through). **This is the mutation gate's actual
  finding: H13's tests currently only gate the most drastic fix, not the
  most likely one.**

## 7. Ledger obligations already incurred

**Row 337 marker debt is NOT closed**, contrary to the R0 phase description
("confirm... and add the missing markers"). Measured: `HardeningFunction.h`
and `CMakeLists.txt` in `SRC/material/nD/ASDPlasticMaterial3D/` still carry
**zero** `Ladruno` mentions; the four `All*.h` registries
(`AllASDHardeningFunctions.h`, `AllASDInternalVariableTypes.h`,
`AllASDModelParameterTypes.h`, `AllASDPlasticMaterial3Ds.h`) also carry zero.
The `Ladruno` hits found in `ASD_material_definitions.cpp` (1),
`gen_ASD_material_definitions_CPP.py` (3) and `OPS_AllASDPlasticMaterial3Ds.cpp`
(4) are all from the *later* ADR-84 P0/P2a PRs, not markers for the original
2026-07-22 HB/StiffSoil extension the row complains about. The ledger row
text itself still reads "UNMARKED... add them" in the present tense. This
review's own R0 acceptance line was not fulfilled and should not be marked
done until the markers are actually added in a commit.

Draft `LEDGER_quirks.md` entries (for the owner to paste in, not added here):

> ## `capfd` cannot see a native `.pyd`'s `cout`/`cerr` on this Windows build
> **Bites:** a test using pytest's `capfd` fixture to assert on native
> extension stdout/stderr silently sees nothing, even when the exact code
> prints normally under a piped shell or `subprocess`.
> **Why:** the `.pyd`'s own linked CRT writes to fd 1/2 through a stream the
> mid-process `dup2` swap `capfd` performs does not reach.
> **Rule:** to assert on native `cout`/`cerr` content, run the model in a
> child process (`subprocess.run([sys.executable, "-c", script], ...)`) and
> capture its real OS-level stdout/stderr. See `_run_child()` in
> `tests/test_adr94_hlist_mechanical.py`.

> ## `subprocess.run` under pytest can raise `WinError 6` non-deterministically
> **Bites:** `subprocess.run(...)` without an explicit `stdin=` intermittently
> (~1/3 runs observed) raises `OSError: [WinError 6] The handle is invalid`
> from `_winapi.DuplicateHandle`.
> **Why:** pytest's own stdio setup does not always leave the parent's stdin
> as a duplicable handle.
> **Rule:** always pass `stdin=subprocess.DEVNULL` to any `subprocess.run`
> hosted under pytest on Windows.

> ## `printA('-ret')` (dense) is empty for every SOE except `FullGeneral`
> **Bites:** `ops.printA('-ret')` silently returns nothing for `UmfPack` and
> other sparse solvers — looks like "no tangent" rather than "wrong API".
> **Why:** `getA()` is null for non-dense SOEs (`OpenSeesCommands.cpp:2718`).
> **Rule:** use `printA('-sparse', '-ret')` (works with `UmfPack`, returns
> `{rowIndices, colIndices, values}`); note it calls `formTangent()` itself,
> so it reads the tangent state *after* the last `update()`.

> ## `TenNodeTetrahedron::eleResponse` self-heals from nodal trial displacement
> **Bites:** "stresses"/"forces"/"material" eleResponse queries always
> re-derive stress from the CURRENT nodal trial displacement, with or
> without a preceding "forces" call — so a material-level Trial-state
> corruption after a domain-level revert is invisible from this channel.
> **Why:** the element does not cache/read the material's own state; it
> recomputes on every query.
> **Rule:** a source-level structural check (regex on the material's own
> `revertToLastCommit`/`revertToStart` body) is the only observation channel
> for this class of defect on this element; do not trust an eleResponse
> query to prove a revert worked at the material level.

## 8. Coverage gaps remaining after R1

Of 46 registered specializations, R1 exercised at runtime **one member each**
of VonMises (H1/H6/H7), MohrCoulomb (H4/H5/H9/H12/H13/H14/H15), DruckerPrager
(H10b, hydrostatic-only), and HoekBrown (H10a) — 4 of 46. **StiffSoilCap (2
combos) has zero runtime evidence of any kind** (not attempted). **StiffSoilShear
(1 combo) has only a NaN failure** on step 1 at every parameter/`InitialP0`
combination tried — no working runtime evidence, and per H11 (noted in R1-B)
no pre-existing Zone-A coverage exists for any StiffSoil combo either. The
IV/hardening-function combinatorics within each family (Tensor/Scalar
Linear/Null/Armstrong-Frederick) are entirely unexercised — R1 never varied
the hardening-function axis. Of 7 `integration_method` tokens,
`Runge_Kutta_45_Error_Control` (non-`_old`) and `Backward_Euler_LineSearch`
never produced a single successful commit in this lane — both have runtime
evidence only of *failure*, not of the behaviour their H claims about. The
R4 integrator × tangent matrix (35+ cells per YF) has not started; R1's
handful of measured cells (mostly VM × {Continuum, Secant, Elastic, Num1,
Num2} on Backward_Euler) is not a substitute for it.
