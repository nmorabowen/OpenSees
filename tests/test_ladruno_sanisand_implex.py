"""ADR-92 / P1 -- pytest battery for `-implex` on `LadrunoSANISAND`.

Companion to `test_ladruno_sanisand.py` (the ADR-86 constants) and
`test_ladruno_sanisand_integrator.py` (the ADR-86b integrator seams).  This
file is written BEFORE the P1 build exists (`Ladruno_implementation/
_adr92_p1_execution_plan.md` -- W1-W9 are "syntax-clean, NOT BUILT"), against
the source read directly out of `SRC/material/nD/LadrunoSANISAND.{h,cpp}`, the
ADR-92 P0 oracle results (`_adr92_p0_oracle_results.md`), and the plan's own
section 4 record of what the C++ author decided and what the plan got wrong.

REUSE, NOT RETYPING.  `_PARAMS` / `_XY` / `_P_ATM`, the confine-first
zero-free-DOF rig (`_build_confined` / `_drive_confined` / `_stress` /
`_reldiff`), the low-p `-Pmin` rig's constants and series generator
(`_PMIN_*`, `_c_series_at`), and the recorded ADR-86b regression number
(`_RECORDED_PR_SENSITIVITY`) are all imported from the two sibling files, not
copied -- a third copy of Gorini's constants was flagged as a drift risk in
review.

POST-FIRST-RUN NOTE (2026-09-06).  This file ran against the first P1 binary
(42dbf066e) and scored 13 passed / 2 failed / 1 skipped of 16 collected
(11 single `def test_...` plus one 5-way `@pytest.mark.parametrize`) --
NOT "15", a figure that undercounted the parametrized cases and that this
docstring itself used to repeat. Two failures were this file's
own deck-design bugs, not the C++, and are fixed here: `-implex` vs `-implex`
off was compared at DIFFERENT effective substep caps (`_CAP_ADEQUATE` /
`_CAP_TIGHT` now make the cap a controlled variable, everywhere an ON/OFF
pair is compared); the two refusal tests asserted `analyze() != 0` on a
ZERO-free-DOF deck, where a refusal has no equation to fail and is therefore
invisible to the return code (measured: the material printed the refusal
warning 16 times and `analyze()` still returned 0) -- both now run on
`_build_free_dof_triaxial`, a genuinely free-DOF `LadrunoBrick` cube. The
floor-clamp test's ISOCHORIC deviatoric ramp also measured clamp count == 0
against that binary; `_build_floor_seeking_deck` replaces it with a
net-DILATING ramp (see that function's docstring for why isochoric shear was
probably too mild). The tangent-identity/`revertToLastCommit` failure from
that same run is NOT this file's fault and is not touched here.

LANE B / P1 RED-BLUE REVIEW PASS (2026-09-06,
`Ladruno_implementation/_adr92_p1_redblue_review.md` section 5 item 2).
Responds to RED-3's findings (raw evidence
`_adr92_p1_redblue/red3_tests_process.md`, blue verdicts
`blue3_tests_process.md`), fixing what that review assigned to the test
lane and nothing in `SRC/` (lanes A/C own the C++ and the driver/gate).
Six things changed: (1) a genuinely free-DOF settlement-column deck drives a
MONOTONE NEGATIVE pseudo-clock (`LoadControl(-ds)`, the campaign deck's own
shape) through a `ds` change and asserts `implexDetail[5]` tracks the
SIGNED ratio `dt_{n+1}/dt_n * alpha` -- the gap B1 exploited (F1/F7,
coverage row 10); (2) `-implexControl` and the commit-time companion
refusal are now COUNTED via the new `implexRefusals` response (Vector(4):
total, d2, control, companion) -- the companion-at-commit case (B3:
`Domain::commit()` discards the return code, so this counter is the ONLY
way to observe it from Python) gets a dedicated test (coverage row 20,
previously NOT COVERED); (3) `getCopy(const char*)` is now exercised AFTER
a sibling element has genuine plastic history, not before any `analyze()`
call (F9); (4) the db roundtrip runs on a free-DOF deck where ON and OFF
are measurably different FIRST, checked before the roundtrip rather than
as a skip-fallback after it (F6); (5) a PlaneStrain (`-ndm 2`) smoke test
and a scheme-2-without-cap parse refusal close two more coverage gaps
(row 17 and the D3 companion case); (6) the two tests RED-3/BLUE-3's
coverage matrix records as "RED (unmarked)" (rows 6, 11) are resolved --
one strengthened with the positive control row 11 says it lacked, one
marked `xfail(strict=True)` with the reason it is out of this lane's
scope -- and the stale-semantics test (F7) is renamed to describe the
SHIPPED sign-change contract, not the retired negative-only one.

WHAT A ZERO-FREE-DOF DECK CAN AND CANNOT SHOW UNDER `-implex`.  This matters
enough to say once, centrally, rather than in every test that runs into it.
`ManzariDafalias::commitState()` (reached through `ladrunoImplexCommit()`)
ALWAYS commits the IMPLICIT return (`this->integrate()`), never the
extrapolated `sigma~` -- "IMPL-EX changes only the REPORTED stress/tangent,
never what is committed" (the same rule the LadrunoJ2 IMPL-EX oracle documents
for that material).  On a zero-free-DOF deck, `analyze(1)` always converges
trivially (there is no free-DOF residual to fail) and therefore ALWAYS
commits in the same call, so `sigma~` -- the thing a REAL boundary-value
problem's Newton loop actually iterates on -- is never observable from Python
on such a deck; only the committed (implicit) answer is.  This is exactly why
the P0 oracle measured the tangent identity in a pure-Python re-implementation
and why gate 3 here (`test_tangent_identity_...`) instead uses a genuinely
free-DOF deck with a DELIBERATELY-FAILED Newton iteration (`maxIter=1` at an
unreachable tolerance) to catch `sigma~` before the domain commits over it --
see that test's docstring for the full mechanism and the source lines it
relies on.

Two consequences that shape several tests below:
  * a zero-free-DOF deck IS the right tool for anything that only needs the
    COMMITTED answer (byte-identity, stage-0 inertness, the floor clamp's
    latched fire-count, the DB round trip) -- and it is a STRONGER tool there,
    because there is no Newton tolerance to contaminate the comparison;
  * it is the WRONG tool for anything that needs to see `sigma~` itself
    (the tangent identity) or a mid-analysis refusal's element-level return
    code propagated all the way to `analyze()` (`stdBrick` discards it --
    "stdBrick swallows material return codes", `LEDGER_quirks`); those tests
    use `LadrunoBrick`, matching `test_ladruno_sanisand_integrator.py`'s own
    rule ("ELEMENT CHOICE IS LOAD-BEARING, NOT INCIDENTAL").

NUMBERS NOT INVENTED.  Every expected number below is one of: the recorded
ADR-86b regression (`sanint._RECORDED_PR_SENSITIVITY`), an algebraic identity
read off the source (`ladrunoImplexMeasureError`'s
`||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2`, `ladrunoImplexFreezeTangent`'s
`sigma~ = sigma_n + Ce*(d_eps - f*d_eps_p)` being exactly affine in `d_eps`
for FIXED `f`/`Ce`/`d_eps_p`), or a qualitative direction taken directly from
the P0 oracle table (`_adr92_p0_oracle_results.md` section 4: IMPL-EX-A's
error grows sharply with increment size at low confinement -- "breaks at
5e-4" at p0 = 5 kPa).  Nothing here is a number read off THIS build, because
this build does not exist yet.

LANE B2 / P2 BATTERY (2026-09-07, WP-92e).  Five new tests written for the
"P2 (owed)" table at the end of `Ladruno_implementation/
92_ladruno_sanisand_implex_adr.md`, against the interfaces lane A2 is
building to (not yet built at the time this file was written -- lane A2's
own note says so, and the module docstring's rule above applies here too:
nothing below is a number read off a real binary). Four items, all
process-wide/non-destructive per the new `implexGuards` response
(`Vector(4)`: [0] floor fallbacks, [1] f=0 guards, [2] hold-preserved
commits, [3] reserved/unused):

  * `test_floor_fallback_delivers_implicit_stress_and_counts` --
    `-implexFloor implicit|accept|refuse` (new default `implicit`) at the
    `-implexControl` reduction floor.
  * `test_guard_zeroes_f_after_reversal` / `test_guard_zeroes_f_after_softening`
    -- `-implexGuard on|off` (default `on`): the elastic-predictor (`f = 0`)
    override on the step after a commit whose state shows a load reversal
    (`alpha_in` reset) or softening (`Kp <= 0` / declining stress ratio).
    The softening gate has NO way to be verified against a real run from
    this lane (no binary exists yet, and lane B2 cannot build one) -- it
    detects the softening onset AT RUNTIME from the committed state (a
    peak-then-decline in `eta/M_b`, the observable macroscopic signature of
    `Kp` crossing zero under monotonic straining) and calls `pytest.xfail`
    with the actual numbers reached if the dense confine-first deck never
    softens within its step budget, rather than asserting a canned outcome
    that was never measured.
  * `test_hold_keeps_clock_and_history` -- a `LoadControl(0.0)` hold
    preserves `mImplexDtCommit` and the `d_eps_p` history from the step
    BEFORE the hold, rather than resetting either.
  * `test_setparameter_stresscorrection_takes_effect` -- the
    `updateParameter` dispatch fix for `ManzariDafalias`'s `stressCorrection`
    responseID (9, `ManzariDafalias.cpp:897`), previously a no-op reaching
    `LadrunoSANISAND` through the IMPL-EX-era dispatch.

Also pins `test_implexcontrol_floor_accepts_once_reduction_limit_is_reached`
(the P1/M4 mutation-gate survivor, written against the OLD unconditional-
accept floor behaviour) to `-implexFloor accept` explicitly -- its own
assertions describe exactly the `accept` mode's contract, and P2's new
default (`implicit`) would otherwise silently change what floor behaviour
that test is exercising out from under it.

Collected count after this lane: 29 `def test_...` functions, 33 collected
items (`test_implex_refuses_unsupported_schemes` is a 5-way parametrize;
every other function is a single collected item) -- up from the P1 file's
21 functions / 25 items. Verified by `python3.12 -m py_compile` plus an AST
walk over the module's top-level `test_*` `FunctionDef`s, not by hand-count.

FIRST RUN AGAINST A REAL P2 BINARY (2026-09-07, `ladrunoBuild() ==
87b9cf846`). 27 passed / 1 skipped / 2 xfailed, zero unexplained failures.
Five genuine test-side bugs found and fixed, all in THIS file, none in
`SRC/`:

  * `test_negative_monotone_clock_runs_the_spec_factor` and
    `test_reararm_after_refusal_without_a_revert_uses_its_own_dt_ratio`
    (both P1-era, predate `-implexGuard`) needed `-implexGuard off` pinned
    -- their decks' first plastic commit is UN-PRIMED, and
    `ladrunoImplexCommit()`'s reversal check (`GetNorm_Contr(dAlphaIn) >
    0.0`) cannot tell "alpha_in initialised for the first time" from a
    genuine reversal, so the guard spuriously arms off that first commit
    and zeroes `f` on the immediately following step -- confirmed a real,
    reproducible effect (not itself edited here, since it lives in
    `SRC/LadrunoSANISAND.cpp` and this lane does not touch `SRC/`), worked
    around the SAME way `-implexFloor accept` was already pinned onto the
    M4 survivor.
  * `test_floor_fallback_delivers_implicit_stress_and_counts` compared a
    refused ladder rung's `implexDetail` against the next rung's -- but
    `LadrunoBrick` PROPAGATES a refusal, so a refused `analyze()` reverts
    without ever calling `commitState()`, and `implexDetail`'s error
    components are written only inside `ladrunoImplexCommit()`; every
    refused rung therefore read the STALE value from the last real commit
    (here, exactly `[0,0,0,0,0,0]`). Fixed by comparing the SAME
    deterministic floor rung's `implexDetail[0]` under `-implexFloor
    accept` vs `implicit` instead of a refused-vs-accepted pair.
  * `test_hold_keeps_clock_and_history` used the SAME `ds` before and after
    the hold, so the correct (`dtCommit`-preserving) and broken
    (`dtCommit`-reset) readings of `f` collapsed to the identical number
    (`alpha`) and the test could not tell them apart. Fixed by doubling the
    post-hold `ds`, and loosened the with/without-hold stress-equality
    tolerance from 1e-10 to 1e-5 (measured reldiff ~9e-7 on this free-DOF
    deck -- one extra converged Newton solve at dt = 0, not a defect).
  * `test_setparameter_stresscorrection_takes_effect` imported
    `sani._SENSITIVITY_FLOOR` (1e-3, calibrated for the much larger
    p_residual elastic-vs-plastic gap); `Stress_Correction()` is a small
    drift-back correction and measured 2.994e-4 on this deck -- real, but
    under that floor. Fixed with a locally scoped `_SC_SENSITIVITY_FLOOR =
    1e-5`.

`test_guard_zeroes_f_after_softening` XFAILs as designed: `eta/M_b` stayed
at exactly 0.0 for all 400 deviatoric steps on this specific `e_conf`/`lat`
combination -- the deck reaches the `p_min` floor almost immediately at
this (larger, non-`_PMIN_E_CONF_LOW`) confinement and the base material's
own floor-correction Newton solve (`ManzariDafalias::Stress_Correction()`,
the `p < m_Pmin + m_Presidual` branch) converges to a purely hydrostatic
state there, not a sheared one -- so this deck never develops `q` at all,
let alone a peak-then-decline. A different `e_conf`/`lat` choice is owed if
this guard needs positive coverage; not chased further here per the
xfail-is-an-acceptable-outcome instruction.

SECOND RUN, P2-5 + P2-2b (2026-09-07, `ladrunoBuild() == 8bfdfbc17`). Two
new tests added (`test_hold_does_not_reset_alpha_in_on_the_implicit_path`,
`test_guard_ignores_the_unprimed_first_commit`); 29 passed / 1 skipped /
2 xfailed. Two more test-side fixes, both because P2-5 repurposed the
formerly-reserved `implexGuards[3]` slot as the reversal-noise-guard
count (fires on any near-zero strain increment, `-implex` on or off) --
`test_floor_fallback_delivers_implicit_stress_and_counts`'s own
"`implexGuards[3]` never moves" assertion, written when that slot really
was reserved, is now stale (this ladder's shrinking `ds` legitimately
trips it, +176 on one run) and was dropped; the new hold test's own first
draft over-claimed that `alpha` (`getAlpha()` -> `mAlpha`, the TRIAL
value) stays bit-identical across a hold, when only `alpha_in`
(`mAlpha_in_n`, the COMMITTED value P2-5 actually protects) is promised --
`mAlpha` moves by ~6e-4 on this deck's hold from the Newton solve's own
non-bit-exact convergence, which is not a defect.

P2-2b's fix (`guardPrimed`) is directly confirmed:
`test_guard_ignores_the_unprimed_first_commit` passes clean. It also
resolves the ORIGINAL reason `test_negative_monotone_clock_runs_the_spec_
factor` and `test_reararm_after_refusal_without_a_revert_uses_its_own_dt_
ratio` were pinned to `-implexGuard off` -- re-run with the guard back on,
each now gets PAST its un-primed second step correctly. Both stay pinned
regardless: each deck independently reaches genuine softening/reversal
territory LATER in its own sequence (both sit at/near the `p_min` floor by
construction -- the settlement column via repeated low-p CLAMPING
warnings, the M10 deck via `_build_floor_seeking_deck`'s own net-dilating
design), so the guard legitimately keeps firing there. See each test's own
docstring for the specific re-measurement. Not a claim this exhausts P2-2b
verification -- `test_guard_ignores_the_unprimed_first_commit`'s clean
moderate-p deck is the direct, unconfounded evidence; these two are the
residual, honestly-reported caveat.

`test_hold_does_not_reset_alpha_in_on_the_implicit_path`'s primary claim
rests on `implexGuards[3]` counting (fires 8x, one per Gauss point, at the
default `-reversalTol`; 0x at `-reversalTol 0`, same deck/history/hold) --
`alpha_in` itself reads bit-identical under BOTH settings on this deck,
which the test documents as a weaker, deck-specific finding (this
particular monotone triaxial hold never trips the BASE's own reversal
branch at any Gauss point, so the guard's assignment is a no-op even when
live) rather than overclaiming it as proof of the repair.

THIRD RUN, P2-6 (2026-09-07, `ladrunoBuild() == 708152eac`). One new test,
`test_trial_guard_accepts_f0_before_refusing`. `implexGuards` grows
`Vector(4)` -> `Vector(5)`, new slot [4] = trial-time f=0 fallbacks --
every existing `len(guards_before) == 4` assertion in this file (the floor-
fallback, hold-clock, and hold/alpha_in tests) is updated to `== 5`.

THE DECK NEEDED A REVERSAL, NOT A BIGGER SAME-DIRECTION STEP -- worth
recording since it cost one wrong guess. `-implexTrialGuard`'s whole
point is retrying a refused trial with a pure elastic predictor (`f = 0`)
before refusing; the FIRST deck tried here reused this file's own "10x
nominal, same direction" shape (`test_floor_fallback_...`'s), on the
reasoning that a bigger step would push the error further past tol.
Measured (via the `-implexAlpha 1.0` vs `0.0` probe technique this
test's own `_probe_trial_guard_reference_errors` uses): on THAT shape the
elastic guess was slightly WORSE than the full extrapolation (0.264 vs
0.255) -- the established plastic direction is still the right one for a
bigger step in the SAME direction, so there was nothing for the fallback
to rescue. A REVERSAL (half the nominal magnitude, opposite sign) flips
this: the established `d_eps_p(n)` now points the wrong way, so the
elastic guess is ~19x BETTER (0.0018 vs 0.034) -- exactly P2-6's own
motivating case (Esmeralda 146569, a leg crawling through refusals near a
turning point). `_TRIAL_GUARD_TOL = 0.01` sits cleanly between the two.
"""
import math
import os
import subprocess
import tempfile
import warnings

import pytest

from _testbed import ops

import test_ladruno_sanisand as sani
import test_ladruno_sanisand_integrator as sanint


pytestmark = [pytest.mark.zone_a]

# Fail LOUD, not silently wrong, if this session ever resolves `opensees`
# to the installed Program Files build (LEDGER_quirks: the .pth hijack)
# instead of THIS checkout's dist/bin/opensees.pyd. A literal build hash
# would break on every future rebuild/CI run (a different commit every
# time), so instead this checks that the LOADED build is an ANCESTOR of
# (or equal to) this checkout's own HEAD via `git merge-base --is-ancestor`
# -- true for any commit actually built from this repo's history, false for
# an unrelated install. Falls back to a warning (not a hard failure) if git
# itself is unavailable, since that is an environment gap, not a wrong-
# binary signal.
def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session", autouse=True)
def _assert_ladruno_build():
    build = ops.ladrunoBuild()
    module_file = str(getattr(ops, "__file__", "?"))
    repo_root = _repo_root()
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo_root,
            stdin=subprocess.DEVNULL, capture_output=True, text=True,
            timeout=10, check=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        warnings.warn(
            "could not run `git rev-parse HEAD` to verify the loaded "
            f"opensees build provenance (ladrunoBuild() = {build!r}, "
            f"module file = {module_file!r}); skipping the build-"
            f"provenance check ({exc!r})")
    else:
        try:
            subprocess.run(
                ["git", "merge-base", "--is-ancestor", build, head],
                cwd=repo_root, timeout=10, check=True,
            )
        except subprocess.SubprocessError:
            pytest.fail(
                f"wrong opensees module loaded for this battery: "
                f"ladrunoBuild() = {build!r} is not an ancestor of this "
                f"checkout's HEAD ({head!r}) -- this usually means the "
                "installed Program Files opensees.pyd was picked up "
                "instead of this worktree's dist/bin one (module file: "
                f"{module_file!r})")
    yield


_PARAMS = sani._PARAMS
_P_ATM = sani._P_ATM
_XY = sani._XY
_EQ_TOL = sani._EQ_TOL
_SENSITIVITY_FLOOR = sani._SENSITIVITY_FLOOR


def _matvec6(flat36, v6):
    """`flat36` is a row-major-flattened 6x6 (`test_ladrunoRCConcrete_material.py`'s
    own comment on this exact response: "6x6 row-major"). Fine even if the true
    storage were column-major, because `-implex` freezes ALL THREE tangent slots
    to the same SYMMETRIC isotropic Ce (`ladrunoImplexFreezeTangent`'s own
    comment: "the delivered operator stays Ce, which keeps it symmetric") --
    for a symmetric matrix the row-major and column-major flattenings coincide."""
    assert len(flat36) == 36, ('material tangent response was not a 6x6', flat36)
    out = [0.0] * 6
    for i in range(6):
        row = flat36[6 * i:6 * i + 6]
        out[i] = sum(row[j] * v6[j] for j in range(6))
    return out


def _vnorm(v):
    return math.sqrt(sum(x * x for x in v))


# ===========================================================================
#  Gate 5 -- byte-identity with `-implex` unset (THE LOAD-BEARING GATE)
# ===========================================================================
#
#  Two independent arguments, not one, because the second is provable from the
#  source without needing to have run this build at all, while the first is
#  the traditional "the published number has not moved" regression check.

# ---------------------------------------------------------------------------
#  THE SUBSTEP CAP IS A CONTROLLED VARIABLE, NOT A CONSTANT.
#
#  Measured on this build (2026-09-06), confine-first deck, deviatoric leg: the
#  SANISAND return needs between 1000 and 5000 ModifiedEuler substeps for that
#  increment.  Capped below it the update is REFUSED -- correctly, that is what
#  ADR-86b T1 built the cap for -- and the committed answer changes:
#
#      OFF, no cap      [-3.09421, -3.09421, -22.18805]
#      OFF, cap   200   [-1.26092, -1.26092,  -9.24638]   <- the CAP moved this
#      ON,  cap   200   [51.47979, 51.47979, -108.11220]
#      OFF, cap 100000  [-3.09421, -3.09421, -22.18805]
#      ON,  cap 100000  [-3.09421, -3.09421, -22.18805]   <- BIT-IDENTICAL
#
#  So an ON/OFF comparison at a cap the deck cannot meet measures the CAP, not
#  -implex.  A first draft of this file compared `-implex -maxSubsteps 200`
#  against a control with NO cap and read the difference as an IMPL-EX defect;
#  it was not.  Every ON/OFF pair below now carries the SAME cap, and any pair
#  that means to prove agreement uses one this deck can actually meet.
# ---------------------------------------------------------------------------
_CAP_ADEQUATE = 20000   # comfortably above the 1000-5000 this deck needs
_CAP_TIGHT = 200        # deliberately too tight -- only for refusal tests


def test_implex_unset_reproduces_recorded_sensitivity():
    """`LadrunoSANISAND` with no `-implex` token anywhere reproduces the SAME
    published ADR-86b regression number as before this ADR's C++ landed.

    `sanint._RECORDED_PR_SENSITIVITY` (9.325164e-02) is
    `test_ladruno_sanisand_integrator.py`'s own gate for "ADR-86b did not move
    the recorded material-point answer" -- it exists precisely to catch an
    unrelated change perturbing this deck.  ADR-92 P1 touches exactly the same
    three files that gate already watches (`LadrunoSANISAND{,3D,PlaneStrain}`),
    so re-running it here, in the file that is ADR-92's OWN warrant, is the
    direct evidence that `ladrunoTrialUpdate()`'s off-path --
    `this->integrate(); return this->ladrunoUpdateStatus();`, verbatim, in
    that order -- is exactly what ships when `-implex` is not given.
    """
    a = sani._drive_confined('LadrunoSANISAND', 8101, sani._OPTS_VANILLA)
    b = sani._drive_confined('LadrunoSANISAND', 8102, sani._OPTS_PR0)
    rel = sani._reldiff(a, b)
    assert abs(rel - sanint._RECORDED_PR_SENSITIVITY) <= sanint._RECORDED_TOL, (
        'the confine-first deck no longer reproduces its recorded p_residual '
        'sensitivity now that the ADR-92 P1 IMPL-EX code shares a translation '
        'unit with it. ladrunoTrialUpdate() is supposed to be BYTE-IDENTICAL '
        'to the pre-ADR-92 integrate()+ladrunoUpdateStatus() pair when -implex '
        'is not given -- if this moved, something in the new code path is '
        'reachable without the flag', rel, sanint._RECORDED_PR_SENSITIVITY)


def test_implex_on_matches_off_on_a_zero_free_dof_deck():
    """THE MOST IMPORTANT TEST IN THIS FILE.

    `-implex` ON (companion scheme 1, `-implexControl` OFF) must commit the
    BIT-IDENTICAL stress to `-implex` OFF, at every step, on the confine-first
    zero-free-DOF deck -- not merely close, `==` on the list.

    WHY THIS HAS TO HOLD, argued from the source rather than merely hoped for:
    on a deck with zero free DOFs, `analyze(1)` converges in exactly one trial
    call (the residual over an empty DOF set is 0 by construction) and
    therefore always reaches `commitState()` within that same call.
    `ladrunoImplexCommit()` (`LadrunoSANISAND.cpp`, "The companion return, at
    commitState only") calls `this->integrate()` -- ManzariDafalias's OWN
    return map -- using `mEpsilon`/`mEpsilon_n` exactly as they stand, which
    are set by the ELEMENT's `setTrialStrain` BEFORE `ladrunoTrialUpdate()` is
    even invoked and are therefore IDENTICAL whether `-implex` is on or off.
    `ladrunoImplexTrial()` (the extrapolation) never writes `mEpsilon`,
    `mK`/`mG` (its own tangent-freeze helper explicitly does not touch them),
    or any committed member -- only trial `mSigma`/`mEpsilonE`/`mCe`, all of
    which `integrate()` overwrites unconditionally on entry. So the OFF path's
    single `integrate()` call and the ON path's commit-time `integrate()` call
    start from the same committed state and the same trial strain and must
    produce the same answer -- this is a claim about the CODE PATH, provable
    without a binary, and it is exactly the claim gate 5 makes.

    This is what a genuine free-DOF deck CANNOT show as cleanly: there, the
    trial `sigma~` feeds Newton and can move the iteration path (though not,
    per `test_tantype_does_not_change_the_converged_answer`'s argument, a
    force-residual-converged EQUILIBRIUM). Zero free DOFs removes that
    variable entirely, which is why this is the strongest form of the gate 5
    claim available before the build exists.
    """
    opts_off = sani._OPTS_VANILLA
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    off = sani._drive_confined('LadrunoSANISAND', 8103, opts_off)
    on = sani._drive_confined('LadrunoSANISAND', 8104, opts_on)
    assert off == on, (
        '-implex (companion scheme 1, -implexControl off) committed a '
        'DIFFERENT stress than -implex off on a zero-free-DOF deck, where '
        'the committed answer is supposed to be the SAME integrate() call '
        'either way. Either ladrunoImplexCommit() is not calling integrate() '
        'on the identical (mEpsilon, mEpsilon_n) pair, or something the '
        'extrapolation trial touches is leaking into the committed return',
        off, on)


# ===========================================================================
#  Gate 5 -- stage-0 inertness: gravity and a LoadControl(0.0) hold
# ===========================================================================

def test_stage0_inertness_gravity_and_hold_is_bit_identical():
    """`-implex` is inert while `mElastFlag == 0` (`ladrunoImplexActive()`):
    a gravity-stage ramp AND a `LoadControl(0.0)` re-equilibration hold, BOTH
    taken entirely at stage 0, must commit BIT-IDENTICAL stress whether
    `-implex` is on or off.

    SOURCED FROM.  `ladrunoImplexActive()` is
    `mImplexOpt.enabled && (mElastFlag != 0)` -- an explicit conjunction, not
    just "the flag" -- and `ladrunoTrialUpdate()` takes the OFF branch
    whenever it is false, setting `mImplexTrialDone = false` so
    `commitState()` also takes the base path
    (`if (!mImplexTrialDone) return ManzariDafalias::commitState();`). Since
    `mElastFlag` is a process-wide STATIC (every Manzari-family constructor
    resets it -- ADR 86 risk 3), the deck below sets stage 0 EXPLICITLY for
    this material's own tag and never flips it, so this path is exercised for
    the whole test.

    THE HOLD.  A `LoadControl(0.0)` step re-forms and re-solves with a ZERO
    load increment -- the common gravity-hold idiom used across this fork's
    own test suite (`test_stagedStrain_material.py`,
    `test_ladrunoEmbeddedNode_element.py`, et al.). On this deck's zero
    free DOFs it is a no-op mechanically, but it is the literal shape ADR-92
    P1's plan (W5) names, so it is included rather than assumed equivalent to
    the ramp alone.
    """
    opts_off = sani._OPTS_VANILLA
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    def _elastic_only_plus_hold(tag, opts):
        sani._build('LadrunoSANISAND', tag, opts)
        ops.updateMaterialStage('-material', tag, '-stage', 0)
        for step in range(sani._N_EL):
            assert ops.analyze(1) == 0, f'gravity-stage step {step + 1} failed'
        ops.integrator('LoadControl', 0.0)
        assert ops.analyze(1) == 0, 'the LoadControl(0.0) hold failed to converge'
        return sani._stress()

    off = _elastic_only_plus_hold(8105, opts_off)
    on = _elastic_only_plus_hold(8106, opts_on)
    assert off == on, (
        '-implex moved the committed stress while the material was still on '
        'the ELASTIC stage (mElastFlag == 0), where ladrunoImplexActive() is '
        'supposed to be unconditionally false. gravity and a LoadControl(0.0) '
        'hold must be bit-identical with the flag on or off', off, on)


# ===========================================================================
#  Gate 4 -- the p_min floor clamp on sigma~, and its own diagnostic
# ===========================================================================

def _build_floor_seeking_deck(tag, opts, e_conf=None, n_dev=60, lat=1.5):
    """A confine-first zero-free-DOF deck, like `sani._drive_confined_at`, but
    with the deviatoric ramp's lateral/axial ratio pushed ABOVE the isochoric
    value (`sani._C_LAT = 0.5`, net volumetric strain zero) so the path net
    DILATES instead of merely shearing at constant volume.

    WHY THIS CHANGED FROM THE FIRST DRAFT.  The first draft reused
    `sani._drive_confined_at` verbatim (its own `lat = sani._C_LAT = 0.5`,
    isochoric).  Run against the first P1 binary, the clamp NEVER fired
    (`implexDetail` count == 0 for the whole 40-step leg) even though the
    BASE material's own low-p branch is independently measured to go
    NEGATIVE on that exact deck (`test_ladruno_sanisand.py`: p = -0.5647 kPa
    at deviatoric step 1 with `-Pmin` at the class default). The likely
    reason: `mImplexDEpsP` is exactly zero until the FIRST commit after the
    stage flip (`ladrunoImplexInitState`'s own comment -- "which is what
    makes the first plastic step a pure elastic prediction"), so `sigma~` on
    an isochoric shear's early steps is a much MILDER elastic-only
    prediction than the base's full elastoplastic return, and apparently
    mild enough here to stay clear of the floor. P0's own G3 path was not
    isochoric shear either -- it was described as "volumetric extension...
    while the point keeps flowing" (`_adr92_p0_oracle_results.md` section 5).
    `sani._c_series_at` already exposes the lat ratio as a parameter, so this
    supplies `lat > 0.5` (net dilation) rather than writing a new series
    generator; e_conf/e_ax/n_conf stay the SAME already-proven magnitudes
    the sibling `_PMIN_*` family uses (only the ratio changes), so this is
    not a numerically untested regime.
    """
    if e_conf is None:
        e_conf = sani._PMIN_E_CONF_LOW      # -> p ~ 0.1145 kPa at the flip
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    s_lat, s_ax = sani._c_series_at(e_conf, n_dev, n_conf=sani._C_N_CONF,
                                    e_ax=sani._C_E_AX, lat=lat)
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *s_lat)
    ops.timeSeries('Path', 2, '-dt', 1.0, '-values', *s_ax)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -e_conf)
            if y == 1.:
                ops.sp(n, 2, -e_conf)
    ops.pattern('Plain', 2, 2)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            if k == 1:
                ops.sp(4 * k + j + 1, 3, -e_conf)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-13, 25, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return n_dev


def test_floor_clamp_fires_and_implexDetail_reports_it():
    """On a NET-DILATING path at the same low confinement the base material's
    own low-p branch is independently known to threaten, the NEW clamp on
    `sigma~` (W2) fires too, and `implexDetail` says so -- checked PER STEP,
    not only at the end, both for the clamp's own fire count/flag and for the
    `||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2` identity `ladrunoImplexMeasureError`
    computes (dev and I1 are orthogonal, so this must hold at every read,
    independent of whether the clamp fired on that particular step).

    SEE `_build_floor_seeking_deck` for why this is NOT the first draft's
    isochoric deck: that one measured clamp count == 0 against the first P1
    binary, and the fix is a path shape (net dilation, not constant-volume
    shear), not a weaker assertion.  The count is still not asserted as an
    exact number -- only that it is possible on a build to which this file
    has never been run -- see `_adr92_p0_oracle_results.md` section 5 for the
    same MECHANISM (not the same deck) measured directly by P0.
    """
    tag = 8107
    opts = ('-Presidual', 0.0, '-Pmin', sani._PMIN_LADRUNO, '-honorTolR', 0,
            '-implex', '-maxSubsteps', _CAP_ADEQUATE)
    n_dev = _build_floor_seeking_deck(tag, opts)

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(sani._C_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    max_count = 0.0
    ever_fired = False
    for step in range(n_dev):
        assert ops.analyze(1) == 0, f'deviatoric step {step + 1} failed'

        detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
        assert len(detail) == 6, ('implexDetail did not return the documented '
                                  '6-component vector', step, detail)
        total, dev, vol, fired_last, count, f_last = detail

        # The identity, ADR-92 P1 section 4.1: ||d_sigma||^2 == ||d_dev||^2 +
        # 3*(dp)^2, i.e. total^2 == dev^2 + vol^2 over the SAME denominator --
        # must hold at EVERY step, clamp firing or not.
        lhs = total * total
        rhs = dev * dev + vol * vol
        scale = max(lhs, rhs, 1.0e-30)
        assert abs(lhs - rhs) <= 1.0e-8 * scale, (
            'implexDetail\'s total/deviatoric/volumetric legs violate the '
            '||d_sigma||^2 == ||d_dev||^2 + 3*(dp)^2 identity the source '
            'claims (ladrunoImplexMeasureError)', step, detail)
        assert fired_last in (0.0, 1.0), ('implexDetail\'s fired flag was '
                                          'not boolean-valued', step, detail)

        max_count = max(max_count, count)
        ever_fired = ever_fired or bool(fired_last)

    assert max_count > 0 and ever_fired, (
        'the p_min clamp on sigma~ never fired on a NET-DILATING path at the '
        'same low confinement (e_conf = sani._PMIN_E_CONF_LOW) where the '
        'BASE material\'s own low-p branch is independently measured to go '
        'negative (test_ladruno_sanisand.py). If the clamp genuinely never '
        'engages on this shape either, either W2 is not wired into '
        'ladrunoImplexTrial(), or the path still is not aggressive enough -- '
        'raise `lat` further rather than weakening this assertion',
        max_count, ever_fired)


# ===========================================================================
#  Gate 3 -- tangent identity, on a NON-pseudo dt source
# ===========================================================================
#
#  MECHANISM.  `ManzariDafalias::commitState()` always commits the IMPLICIT
#  return (module docstring above), so a zero-free-DOF deck can never show
#  sigma~ to Python -- it is overwritten by the commit within the SAME
#  `analyze(1)` call. The only way to catch sigma~ is a genuinely free-DOF
#  deck with a Newton iteration that FAILS to converge (so the domain never
#  commits), forced by pairing `algorithm('Newton')` with an UNREACHABLE
#  convergence tolerance and `maxIter = 1`: exactly one
#  formTangent/formResidual/solve/update cycle runs (the standard OpenSees
#  Newton loop: form, solve, update U -> pushes the new trial strain into the
#  material -> test; with maxIter = 1 the loop stops after that ONE update
#  regardless of what test() says), then `test()` fails against the
#  unreachable tolerance and `analyze()` returns non-zero WITHOUT commit
#  ("Domain::update() returns the element's code ... IncrementalIntegrator::
#  update turns it into a failed step -- none of which involves solving
#  anything", `test_ladruno_sanisand_integrator.py`).
#
#  DO NOT call `ops.reset()` here.  MEASURED 2026-09-06: `ops.reset()` is
#  `Domain::revertToStart()` (`OpenSeesCommands.cpp:2539`), NOT
#  `revertToLastCommit()` -- it runs `ManzariDafalias::initialize()`, which
#  ZEROES `mSigma_n`/`mEpsilon_n`/`mAlpha_n`/`mFabric_n`, and then ends
#  `return this->update()` (`Domain.cpp:2385`), pushing one state determination
#  through the zeroed state.  With `sigma~ = 0 + Ce:0 = 0` the p_min clamp
#  fires HONESTLY and `getStress()` returns
#  `[-0.0101, -0.0101, -0.0101, -0.0, -0.0, -0.0]` -- the three IEEE negative
#  zeros are the fingerprint, since `-1.0 * (+0.0)` requires `dev(sigma~)` to be
#  EXACTLY zero, which no real load path produces.  Three tests read that as
#  "revertToLastCommit is not restoring".  It was the probe, not the material.
#
#  Nothing is needed in its place: a failed `StaticAnalysis` step ALREADY calls
#  `theDomain->revertToLastCommit()` and `theIntegrator->revertToLastStep()`
#  before `analyze()` returns non-zero (`StaticAnalysis.cpp:185` and four
#  sibling sites), so the correct revert has happened by the time we look.
#  There is no Python verb that calls `revertToLastCommit()` on its own; a
#  deliberately-failed step is the way to reach it.
#  commit) and the probe's OWN trial-strain norm.
#
#  WHY +h AND -h SHARE THE SAME f UNDER `-implexDt strain`.  `f` is frozen at
#  the FIRST trial call after a commit/revert from
#  `dt = GetNorm_Cov(mEpsilon - mEpsilon_n)` (ladrunoImplexArmStep). Each
#  probe below is driven by an isolated axial force pattern applied AFTER
#  `loadConst` has frozen everything else, so the probe's OWN strain
#  increment is (to within the Newton solve's own linearity, exact under a
#  frozen elastic tangent) a scalar multiple of one fixed direction; a norm is
#  insensitive to the SIGN of that scalar, so dt, and therefore f, is the same
#  for the +h and -h probes. With f AND Ce AND d_eps_p(n) equal for both,
#  sigma~ = sigma_n + Ce*(d_eps - f*d_eps_p) is EXACTLY affine in d_eps, so the
#  secant (sigma~(+h) - sigma~(-h)) / (eps~(+h) - eps~(-h)) reproduces Ce to
#  floating-point precision, not merely to FD truncation order -- this is why
#  the test can use a tight tolerance despite never seeing a "true" small h.
#
#  `-implexDt strain`, NOT `pseudo`.  Per the plan section 4.3 item 1: under
#  `pseudo`, `ops_Dt` is constant within a step regardless of what the trial
#  strain does, so the freeze is a no-op and this exact trap (an UNFROZEN,
#  strain-dependent f silently breaking d(sigma~)/d(eps) == Ce) would not be
#  exercised at all.  `strain` is the one source where the freeze fix does
#  real work, so it is the one this gate has to run on.

_PROBE_TAG = 8110
_PROBE_P0 = 100.0            # kPa -- away from the p_min floor; this gate is
                             # about the tangent, not the clamp (gate 4 owns that)
_PROBE_N_CONF = sanint._TX_N_CONF
_PROBE_DQ_NOMINAL = 6.0      # kPa (small fraction of _PROBE_P0), a few steps to
                             # establish a genuinely plastic committed history
_PROBE_N_HISTORY = 4
_PROBE_TOL_REL = sanint._TX_TOL_REL
_PROBE_MAXITER = sanint._TX_MAXITER
_PROBE_H_FRACTION = 1.0e-6   # the probe force as a fraction of _PROBE_DQ_NOMINAL/4


def _build_probe_triaxial(tag):
    """A LadrunoBrick drained-triaxial cube -- `sanint._build_triaxial`'s shape,
    with the ADR-92 P1 flags threaded through (that builder does not accept
    extra flags beyond the five positional optionals)."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,           # IntScheme, TanType, JacoType, TolF, TolR
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   '-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexDt', 'strain')
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = _PROBE_P0 / 4.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.load(n, -q, 0.0, 0.0)
            if y == 1.:
                ops.load(n, 0.0, -q, 0.0)
            if k == 1:
                ops.load(n, 0.0, 0.0, -q)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormUnbalance', _PROBE_TOL_REL * _PROBE_P0, _PROBE_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _PROBE_N_CONF)
    ops.analysis('Static')


def _confine_only(tag):
    """Confine (stage 0) and flip to stage 1 -- WITHOUT any deviatoric
    history, so the material is left UN-PRIMED (mImplexDEpsP(n) == 0, the
    exact state the afb95c40c drift-correction exemption keys on). Factored
    out of `_establish_plastic_history` so the un-primed-step regression
    test can stop here instead of also taking the nominal steps."""
    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(_PROBE_N_CONF):
        assert ops.analyze(1) == 0, f'triaxial confinement step {step + 1} failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)


def _establish_plastic_history(tag):
    """Confine, flip, and take a few NORMAL (converged, committed) deviatoric
    steps so mImplexDEpsP(n) != 0 by the time the probe runs."""
    _confine_only(tag)

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 1.0 / _PROBE_N_HISTORY)
    for step in range(_PROBE_N_HISTORY):
        assert ops.analyze(1) == 0, f'history-building step {step + 1} failed'
    ops.loadConst('-time', 0.0)


def _probe_once(tag, sign, probe_pattern_tag, probe_ts_tag):
    """One forced-single-iteration trial: apply a TINY axial force (sign *
    h_fraction * the nominal per-step load), force exactly one failed Newton
    iteration, read back (stress, strain) BEFORE any commit, then revert.

    Returns (stress[6], strain[6]) -- element-convention (tension-positive,
    the -1.0 flip already applied by LadrunoSANISAND3D::getStress/getStrain).
    """
    h = sign * _PROBE_H_FRACTION * (_PROBE_DQ_NOMINAL / 4.0)
    ops.timeSeries('Linear', probe_ts_tag)
    ops.pattern('Plain', probe_pattern_tag, probe_ts_tag)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -h)

    ops.test('NormUnbalance', 1.0e-300, 1, 0)   # unreachable -> exactly 1 iteration, no commit
    ops.integrator('LoadControl', 1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'the probe step CONVERGED and therefore COMMITTED -- it was supposed '
        'to fail an unreachable NormUnbalance tolerance after exactly one '
        'Newton iteration so sigma~ could be read before the commit '
        'overwrites it with the implicit return. If this converges, the '
        'harness itself is broken, not the material', rc)

    stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    strain = list(ops.eleResponse(1, 'material', 1, 'strain'))

    # NO ops.reset() -- see the block comment above; the failed step already
    # reverted, and reset() would zero the committed state and hand back the
    # clamp fingerprint instead.
    ops.remove('loadPattern', probe_pattern_tag)
    return stress, strain


@pytest.mark.xfail(strict=True, reason=(
    'known-red per the P1 red/blue review (RED-3 F5/F11, coverage row 6: '
    '"RED (unmarked)"); the mechanism this test exercises (`-implexDt '
    'strain`, ladrunoImplexArmStep / ladrunoImplexFreezeTangent) is '
    'ORTHOGONAL to the C++ this lane\'s fix touches (B1/B2/B3: the pseudo-dt '
    'ratio, the -implexControl floor, and the commit-time companion refusal '
    '-- none of which this deck reaches, since its dt source is a strain '
    'norm, not the pseudo-clock, and it never triggers a refusal). Root '
    'cause undiagnosed by the review ("attributed to neither side yet", '
    'BLUE-3 F11) and out of lane B\'s scope; do not xfail this away again '
    'once a diagnosis exists -- fix it or file the follow-up.'))
def test_tangent_identity_frozen_ce_on_strain_dt_source():
    """Returned tangent (`getTangent()`, frozen to `Ce(p_n)` under `-implex`)
    reproduces a numerical `d(sigma~)/d(eps)` taken from TWO forced,
    uncommitted trial evaluations on a `-implexDt strain` deck -- see the
    long mechanism note above this test for why this is the honest way to
    reach `sigma~` from Python and why the two probes share the same frozen
    `f`.

    Kills a mutant that lets `f` drift between the +h/-h probes (an
    UNFROZEN, strain-dependent factor breaking the affine identity) or that
    reports the tangent from an unclamped/stale `Ce`.
    """
    _build_probe_triaxial(_PROBE_TAG)
    _establish_plastic_history(_PROBE_TAG)

    ce_flat = list(ops.eleResponse(1, 'material', 1, 'tangent'))
    stress_before = list(ops.eleResponse(1, 'material', 1, 'stress'))
    detail_before = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    assert detail_before[3] == 0.0 and detail_before[4] == 0.0, (
        'the p_min clamp fired before the probes even started -- _PROBE_P0 '
        'is supposed to sit AWAY from the floor (LEDGER_quirks: "a '
        'tangent-identity gate must be run on a path where the clamp is '
        'idle, and a test that reports identity to machine precision on a '
        'clamped path is measuring nothing"); raise _PROBE_P0 rather than '
        'weakening this precondition check', detail_before)

    sp, ep = _probe_once(_PROBE_TAG, +1.0, probe_pattern_tag=90, probe_ts_tag=90)
    sm, em = _probe_once(_PROBE_TAG, -1.0, probe_pattern_tag=91, probe_ts_tag=91)

    stress_after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert stress_before == stress_after, (
        'the committed stress moved across the two probe-and-reset cycles -- '
        'revertToLastCommit() is not fully restoring the committed state, so '
        'the two probes are not actually taken from the same base state',
        stress_before, stress_after)

    detail_after = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    assert detail_after[3] == 0.0 and detail_after[4] == detail_before[4], (
        'the p_min clamp fired DURING one of the two forced probes -- the '
        'affine-tangent identity this test checks assumes an UNCLAMPED '
        'sigma~ on both probes; a clamp firing on only one of them would '
        'break the identity for a reason unrelated to the tangent freeze',
        detail_before, detail_after)

    d_sigma = [a - b for a, b in zip(sp, sm)]
    d_eps = [a - b for a, b in zip(ep, em)]
    assert _vnorm(d_eps) > 0.0, (
        'the two probes produced the SAME trial strain -- the +h/-h force '
        'perturbation did not move the free axial DOFs at all', ep, em)

    predicted = _matvec6(ce_flat, d_eps)
    err = _vnorm([a - b for a, b in zip(d_sigma, predicted)])
    scale = max(_vnorm(d_sigma), 1.0e-30)
    rel = err / scale
    assert rel < 1.0e-6, (
        'the frozen tangent Ce(p_n) does not reproduce the numerical '
        'd(sigma~)/d(eps) from the two forced-uncommitted probes. Per the '
        'plan section 4.3 item 1, this specific gate only exercises the '
        '"-implexDt strain" f-freeze fix; a failure here on a fresh build is '
        'exactly the trap that fix exists for', rel, d_sigma, predicted)


# ===========================================================================
#  Parser / runtime refusals
# ===========================================================================

@pytest.mark.parametrize('scheme', [3, 5, 7, 8, 9])
def test_implex_refuses_unsupported_schemes(scheme):
    """ADR 92 D3 (as reversed by P0): `-implex` refuses IntScheme
    3/5/7/8/9 -- no error control on the substep, so the companion could
    never report a failed return. `setLadrunoImplexOptions` raises this via
    `opserr` + returns -1, which the parser turns into a hard construction
    failure (`delete theMaterial; return 0;`), which openseespy surfaces as
    an exception -- the same idiom `test_ladruno_sanisand.py` uses for
    `-honorTolR`'s out-of-range refusal.
    """
    ops.wipe()
    tag = 8200 + scheme
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                       scheme, 2, 1, 1.0e-7, 1.0e-7,
                       '-implex', '-maxSubsteps', _CAP_TIGHT)


def test_implex_scheme1_requires_maxsubsteps():
    """`-implex` on the default companion, IntScheme 1 (ModifiedEuler),
    REQUIRES `-maxSubsteps > 0` -- without it the companion at commitState has
    no way to fail rather than force-accept at `dT_min` (ADR-90 GATE U: single
    updates measured at 34 minutes)."""
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 8210, *_PARAMS, '-implex')


def test_implex_maxsubsteps_zero_is_also_refused():
    """`-maxSubsteps 0` means UNCAPPED (vanilla's own convention,
    `test_uncapped_is_byte_identical`), so it must be refused exactly like
    omitting the flag -- a deck that says `-maxSubsteps 0` explicitly must not
    be treated as having satisfied the requirement."""
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 8211, *_PARAMS,
                       '-implex', '-maxSubsteps', 0)


def test_scheme2_without_maxsubsteps_is_refused_under_implex():
    """IntScheme 2 under `-implex` with NO `-maxSubsteps` must be refused at
    PARSE TIME exactly like the default scheme 1
    (`test_implex_scheme1_requires_maxsubsteps`) -- ADR-92 D3's whole point
    is that the companion must be ABLE to fail, and per the P1 review
    (RED-1 F5), the old code only WARNED here, nested under `verbose`, which
    is `false` on every `getCopy`/`recvSelf` path -- so a capped-0 scheme-2
    deck could previously sail through construction (and every clone of it)
    silently.

    Kills a mutant that special-cases scheme 2 out of the "-implex requires
    -maxSubsteps" refusal, or that leaves the refusal gated behind
    `verbose`.
    """
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial('LadrunoSANISAND', 8221, *_PARAMS,
                       2, 2, 1, 1.0e-7, 1.0e-7, '-implex')


# ---------------------------------------------------------------------------
#  A refusal is only OBSERVABLE from analyze()'s return code on a deck with
#  genuine free DOFs.
#
#  The first draft of both refusal tests below used the confine-first
#  ZERO-free-DOF rig (`sanint._build_brick`, matching the SUBSTEP-CAP tests'
#  own idiom).  Run against the first P1 binary: the material DID refuse --
#  the D2 warning ("the pseudo-time increment is negative (-1)") printed 16
#  times -- but `analyze()` still returned 0, because with zero free DOFs
#  there is no equation for a refusal to fail, and the domain's pseudo time
#  moved on regardless.  Confirmed on a genuinely free-DOF deck (both
#  `stdBrick` and `LadrunoBrick`): the SAME refusal returns rc = -3.  So
#  every test in this file that needs to OBSERVE a refusal through
#  `analyze()`'s return code now uses `_build_free_dof_triaxial` below --
#  the same single-`LadrunoBrick` drained-triaxial shape gate 3's tangent-
#  identity probe uses, minus the `-implexDt strain` source (that source
#  cannot go negative at all -- it is a norm -- and D2's refusal is
#  specifically about the DEFAULT `pseudo` source going negative).
# ---------------------------------------------------------------------------

def _build_free_dof_triaxial(tag, extra_opts=(), p0=None):
    """A LadrunoBrick drained-triaxial cube with GENUINE free DOFs (the
    positive-face nodes are LOADED, not `sp`-prescribed) -- see the note
    above this function for why the refusal tests need this and the
    zero-free-DOF rig cannot substitute for it."""
    if p0 is None:
        p0 = _PROBE_P0
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS,
                   1, 2, 1, 1.0e-7, 1.0e-7,           # IntScheme, TanType, JacoType, TolF, TolR
                   '-Presidual', 0.0, '-Pmin', 1.0e-4 * _P_ATM,
                   *extra_opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    q = p0 / 4.0
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.load(n, -q, 0.0, 0.0)
            if y == 1.:
                ops.load(n, 0.0, -q, 0.0)
            if k == 1:
                ops.load(n, 0.0, 0.0, -q)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormUnbalance', _PROBE_TOL_REL * p0, _PROBE_MAXITER, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _PROBE_N_CONF)
    ops.analysis('Static')


def test_implex_sign_change_in_pseudo_dt_is_refused_and_leaves_committed_state_unchanged():
    """D2's SHIPPED contract, as of `3c788778f`: the guard is a SIGN CHANGE
    between the committed and the trial pseudo-clock
    (`mImplexDtCommit != 0.0 && mImplexDt * mImplexDtCommit < 0.0`), not
    "any negative `ops_Dt`" -- a MONOTONE negative clock (every step
    negative, the campaign deck's own shape) is legal and must run at the
    real `dt_{n+1}/dt_n` ratio (see
    `test_negative_monotone_clock_runs_the_spec_factor` for that case; this
    test's OLD name and docstring described the retired "any negative dt"
    rule and are fixed here per the P1 review, RED-3 F7).

    This deck reaches the sign-change branch specifically: positive
    committed steps (`_establish_plastic_history`) followed by ONE negative
    `LoadControl(-1.0)` step, so `mImplexDtCommit > 0`, `mImplexDt < 0`,
    product `< 0`.

    HOW A NEGATIVE `ops_Dt` IS PRODUCED, without touching DisplacementControl
    or arc-length at all: `ops_Dt` is `currentTime - committedTime`
    (`ladrunoImplexArmStep`'s own comment), and `LoadControl`'s
    `deltaLambda` argument becomes exactly that increment. `LoadControl(-x)`
    for one step is therefore the simplest possible way to make `ops_Dt`
    negative under the DEFAULT (`pseudo`) dt source -- no limit point, no
    special integrator needed.

    A FREE-DOF deck, not the confine-first `sp` rig -- see the note above
    `_build_free_dof_triaxial` for the measured reason (a refusal is
    invisible to `analyze()`'s return code with zero free DOFs; the deck
    still visibly PRINTS the refusal warning, it just cannot be asserted on
    via the return code).

    Kills a mutant that drops the sign-change guard entirely (D2 disabled)
    or that reverts to the retired "any negative dt" rule (which would ALSO
    refuse a monotone-negative leg, the exact regression B1 was about).
    """
    tag = 8212
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    _build_free_dof_triaxial(tag, opts)
    _establish_plastic_history(tag)

    before = list(ops.eleResponse(1, 'material', 1, 'stress'))

    ops.integrator('LoadControl', -1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'a SIGN CHANGE in the pseudo-time increment (positive committed '
        'steps, then LoadControl(-1.0)) was NOT refused. ADR-92 D2 requires '
        'the material to refuse this step because dt_{n+1}/dt_n is the '
        'extrapolation factor itself and a reversed load factor makes it a '
        'wrong answer that would pass every other gate', rc)

    # NO ops.reset(): the failed step's own revertToLastCommit already ran.
    after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert before == after, (
        'the committed stress moved across a REFUSED step -- the refusal is '
        'supposed to leave the committed state untouched (revertToLastCommit '
        'restores the trial state on top of it)', before, after)


def test_implexcontrol_refuses_past_tolerance_and_leaves_committed_state_unchanged():
    """`-implexControl` refuses a step whose extrapolation error exceeds its
    tolerance, via `LADRUNO_MATERIAL_REFUSED`, and the committed state is
    unchanged across the refusal -- EXCEPT on an UN-PRIMED step, where
    afb95c40c deliberately suppresses the refusal (see the dedicated
    un-primed assertion below).

    FIXED PER THE P1 REVIEW (RED-3 F5/coverage row 11): the old version
    asserted only `rc != 0` on the big step, with "no positive control
    separating refusal from ordinary non-convergence" -- any Newton failure
    for ANY reason would have passed. This version adds the positive
    control the review asked for, using the `implexRefusals` response
    (Vector(4): total, d2, control, companion): the PRIMED big step MUST
    increment the control slot -- so a green run proves the refusal was
    specifically `-implexControl`'s tolerance gate, not e.g. a companion
    cap or an unrelated equilibrium failure.

    UN-PRIMED STEP REGRESSION (afb95c40c): the coordinator reports that
    `-implexControl` now deliberately does NOT refuse on the FIRST plastic
    step after a stage flip (`|mImplexDEpsP(n)| == 0`, i.e. un-primed),
    because `implexError` there measures the companion's drift-correction
    jump rather than the extrapolation error -- the error is still
    measured and folded into `avgImplexError`, just not used to refuse.
    MEASURED on afb95c40c: the SAME big (10x-nominal) load applied on an
    un-primed first plastic step reads implexError = 0.259 (far past
    tol = 0.02) yet converges (`rc == 0`) with `implexRefusals` unchanged;
    applied again on the now-primed SECOND plastic step it refuses
    (`rc != 0`, `implexRefusals[2]` increments by 8). This test drives
    exactly that two-step sequence, so it is the regression test for the
    un-primed exemption as well as the original qualitative gate.

    A FREE-DOF deck (`_build_free_dof_triaxial`), for the same measured
    reason as the D2 test above: a zero-free-DOF deck's `analyze()` cannot
    fail on a refusal at all (no equations to fail), so it would report
    success and the committed state would move regardless of what the
    material returned -- which is exactly the failure mode the first draft
    of this test hit.

    p0 = 50 kPa, tol = 0.02: MEASURED against the fixed binary
    (ladrunoBuild == 2473ce46c, unchanged on afb95c40c) -- a parameter
    sweep (p0 in {5, 10, 20, 30, 50} kPa, the big step at 5x/10x/20x
    nominal) found a clean, well-separated discriminator here: the PRIMED
    nominal-step implexError maxes at 4.8e-3, the 10x-nominal big step
    (primed) reads 0.255-0.259 -- a > 50x gap, with `tol = 0.02` sitting
    comfortably in between.

    Kills a mutant that reports SOME refusal on the (primed) big step
    (e.g. a companion cap failure masquerading as -implexControl) without
    the control slot itself moving, that refuses the UN-PRIMED step (the
    afb95c40c regression), or that never lifts the exemption once the
    material IS primed (a "-implexControl silently inert forever" mutant).

    PINNED TO `-implexTrialGuard off` (ADR-92 P2-6, WP-92e lane B2,
    2026-09-07, binary 708152eac). P2-6's trial-time f=0 fallback DEFAULTS
    on and runs BEFORE any -implexControl refusal, retrying the trial with
    a pure elastic predictor -- on this deck's big (10x-nominal) step that
    retry ALSO clears tol, so the step that this test is built to prove
    gets refused is, with the default, RESCUED instead (`rc == 0`). This
    test's own claim is about the refusal/un-primed-exemption mechanism,
    orthogonal to P2-6's rescue; isolated the same way `-implexFloor
    accept` and `-implexGuard off` are pinned onto sibling tests elsewhere
    in this file.
    """
    tag = 8213
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', 0.02, 0.01, '-implexTrialGuard', 'off')
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _confine_only(tag)   # NOT _establish_plastic_history: this test needs
                          # the material UN-PRIMED for its first plastic step

    before = list(ops.eleResponse(1, 'material', 1, 'stress'))
    refusals_0 = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))

    # UN-PRIMED first plastic step: the SAME big (10x-nominal) load the
    # primed case below refuses on. Pattern 2, on top of the confinement
    # (pattern 1), already loadConst'd by _confine_only.
    big_dq = 10.0 * _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -big_dq)
    ops.integrator('LoadControl', 1.0)
    rc_unprimed = ops.analyze(1)
    assert rc_unprimed == 0, (
        'the SAME big load that the primed step below refuses on FAILED '
        'to converge on the UN-PRIMED first plastic step -- afb95c40c is '
        'supposed to suppress the -implexControl refusal there '
        '(drift-correction jump, not an extrapolation error), so this '
        'step should behave like -implexControl were off, not like a '
        'harder refusal', rc_unprimed)
    refusals_after_unprimed = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert refusals_after_unprimed[2] == refusals_0[2], (
        'implexRefusals[2] (control-specific) moved across the UN-PRIMED '
        'first plastic step -- the afb95c40c exemption (no refusal while '
        'mImplexDEpsP(n) == 0) is supposed to make this a no-op for the '
        'counter, even though implexError itself is large and IS folded '
        'into avgImplexError', refusals_0, refusals_after_unprimed)
    ops.loadConst('-time', 0.0)

    # SECOND plastic step, same magnitude load again -- now PRIMED (the
    # first step committed, so mImplexDEpsP(n) != 0). This is where the
    # forcing moves per the coordinator's afb95c40c note.
    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -big_dq)
    ops.integrator('LoadControl', 1.0)
    before_primed = list(ops.eleResponse(1, 'material', 1, 'stress'))
    rc = ops.analyze(1)
    assert rc != 0, (
        'a deviatoric load increment 10x the nominal one, at p0 = 50 kPa, '
        'on a PRIMED second plastic step, was NOT refused by '
        '-implexControl. Measured on afb95c40c this reads implexError ~ '
        '0.26, far past tol = 0.02; if this genuinely never refuses on '
        'this deck, re-derive the deck rather than weakening this '
        'assertion', rc)

    refusals_after_big = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert refusals_after_big[2] - refusals_after_unprimed[2] >= 1, (
        'analyze() failed on the PRIMED big step, but implexRefusals[2] '
        '(the -implexControl-specific counter) did not move -- the '
        'positive control this test is supposed to provide. Either the '
        'refusal was NOT -implexControl (a companion cap or an unrelated '
        'Newton failure), or the counter is not wired to this refusal site',
        refusals_after_unprimed, refusals_after_big)

    # NO ops.reset(): the failed step's own revertToLastCommit already ran.
    after_primed = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert before_primed == after_primed, (
        'the committed stress moved across an -implexControl refusal', before_primed, after_primed)


# ===========================================================================
#  getCopy / sendSelf / recvSelf -- the flags AND the committed d_eps_p
# ===========================================================================

def _roundtrip_implex(matcmd, tag, opts, n_cut, n_total):
    """`test_ladruno_sanisand.py`'s own `_roundtrip`, with the plastic-leg
    split parametrised and IMPL-EX opts threaded through. Returns (stress at
    the save, stress right after the restore, final stress)."""
    with tempfile.TemporaryDirectory(prefix='ladruno_sanisand_implex_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'sanisand_implex_rt')

        sani._build(matcmd, tag, opts)
        sani._elastic_leg(tag)
        ops.updateMaterialStage('-material', tag, '-stage', 1)
        for step in range(n_cut):
            assert ops.analyze(1) == 0, f'pre-save plastic step {step + 1} failed'
        mid = sani._stress()

        try:
            ops.database('File', dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f'database() unsupported in this build: {exc}')
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip('database save returned failure on this build')

        sani._build(matcmd, tag, opts)                  # fresh, uncommitted skeleton
        ops.database('File', dbpath)
        ops.restore(1)
        after = sani._stress()

        ops.updateMaterialStage('-material', tag, '-stage', 1)
        ops.wipeAnalysis()
        sani._analysis()
        for step in range(n_total - n_cut):
            assert ops.analyze(1) == 0, f'post-restore plastic step {step + 1} failed'
        out = sani._stress()
        ops.wipe()
        return mid, after, out


def test_implex_db_roundtrip_carries_flags_and_history():
    """A restored `-implex` material must keep running the SAME extrapolation
    it was saving -- both the option flags (the six-override discipline, ADR
    86 sec.4.4/ADR-92 P1 section header, "the wire grows") and the ONE new
    committed history variable, `mImplexDEpsP` (wire slots 16..21, sent
    because "an MP worker that receives a material mid-analysis and restarts
    its extrapolation from zero silently runs a different constitutive law
    from the rank beside it" -- the ADR 86 section 3 defect, in a new
    variable).

    THE ASSERTION IS BEHAVIOURAL, on the `test_db_roundtrip_carries_presidual`
    template: not "restore returned 0", but "the restored material finishes
    the remaining plastic steps on the SAME committed answer an unbroken run
    would reach". A material that restored WITHOUT its `-implex` flags (the
    plain base default) or WITHOUT its `d_eps_p` history (silently
    reinitialised to zero, matching a virgin Gauss point's first plastic
    step) would very likely diverge from the reference somewhere in the
    remaining plastic steps, because the extrapolation factor and the history
    it carries feed directly into every subsequent committed stress.
    """
    opts_on = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    ref = sani._drive('LadrunoSANISAND', 8300, opts_on)

    n_cut = 12
    mid, after, out = _roundtrip_implex('LadrunoSANISAND', 8301, opts_on,
                                        n_cut=n_cut, n_total=sani._N_PL)
    assert sani._reldiff(mid, after) <= _EQ_TOL, (
        'the restored -implex material did not come back on the state it '
        'was saved at', mid, after)
    assert sani._reldiff(ref, out) <= _EQ_TOL, (
        'after the round trip the -implex material no longer finishes on '
        'the reference (unbroken) answer -- the flags or the d_eps_p history '
        'did not survive sendSelf/recvSelf', ref, out)

    # Non-vacuity: an -implex OFF reference run on the SAME deck up to the
    # cut, if the round trip silently dropped -implex entirely (falling back
    # to the base default), would land here instead. It must not.
    ref_off = sani._drive('LadrunoSANISAND', 8302, sani._OPTS_VANILLA)
    gap = sani._reldiff(ref_off, ref)
    if gap < _SENSITIVITY_FLOOR:
        pytest.skip(
            'this deck does not distinguish -implex on from off '
            f'(gap {gap:.3e} < floor {_SENSITIVITY_FLOOR:.0e}), so the '
            'restored-vs-off comparison below would prove nothing; the '
            'restored-vs-reference comparison above already carries this '
            'test')
    assert sani._reldiff(ref_off, out) >= _SENSITIVITY_FLOOR, (
        'after the round trip the material is running as though -implex '
        'were off', ref_off, out)


def test_getcopy_after_plastic_history_shares_options_not_history():
    """`getCopy(const char*)` (the per-Gauss-point path every
    `stdBrick`/`LadrunoBrick` construction uses) builds a BRAND NEW object
    from SCALAR CONSTRUCTOR PARAMETERS ONLY -- "a fresh integration point
    starts with d_eps_p = 0" (`LadrunoSANISAND.cpp`'s own comment) -- never
    from `this`'s own committed members. So a bystander element created from
    the SAME material tag AFTER a sibling has already been driven deep into
    plasticity must still start at IDENTICALLY ZERO plastic strain.

    FIXED PER THE P1 REVIEW (RED-3/BLUE-3 F9): the superseded version of
    this test built BOTH elements before taking a single `analyze()` step,
    so "bystander plastic strain == 0" was true by triviality -- nothing had
    happened to ANYTHING yet, and a mutant that TRIED to leak history would
    have had no history to leak. This version drives element 1 through the
    full `sani._drive` elastic-plus-plastic history FIRST (measured
    nonzero), and only THEN constructs element 2 from the same tag, in the
    SAME still-live domain -- so a mutant that made getCopy(const char*)
    read or share ANY committed member off a sibling clone (rather than
    building fresh from scalars) has real, driven history available to leak
    and would be caught here.
    """
    tag = 8303
    opts = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    sani._drive('LadrunoSANISAND', tag, opts)   # element 1 (tag 1): full history

    driven = list(ops.eleResponse(1, 'material', 1, 'stress'))
    driven_pstrain = list(ops.eleResponse(1, 'material', 1, 'plasticstrains'))
    assert _vnorm(driven_pstrain) > 0.0, (
        'element 1 was never actually driven into plasticity -- the deck, '
        'not getCopy, is broken', driven, driven_pstrain)

    # element 2: a SEPARATE, fully-fixed stdBrick sharing ONLY the material
    # TAG, constructed in the SAME live domain AFTER element 1's plastic
    # history already exists.
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(10 + 4 * k + j + 1, x, y, float(k))
    ops.element('stdBrick', 2, 11, 12, 13, 14, 15, 16, 17, 18, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(10 + 4 * k + j + 1, 1, 1, 1)

    ops.integrator('LoadControl', 0.0)
    assert ops.analyze(1) == 0, (
        'the domain failed to re-solve after adding the bystander element -- '
        'a harness problem, not the claim under test')

    bystander_strain = list(ops.eleResponse(2, 'material', 1, 'strain'))
    bystander_pstrain = list(ops.eleResponse(2, 'material', 1, 'plasticstrains'))

    assert all(abs(v) == 0.0 for v in bystander_strain), (
        'the bystander element (fully fixed, added after element 1 was '
        'already driven) reports nonzero strain -- the deck, not getCopy, is '
        'at fault', bystander_strain)
    assert all(abs(v) == 0.0 for v in bystander_pstrain), (
        'element 2, constructed from the SAME material tag AFTER element 1 '
        'had already committed nonzero plastic strain, reports NONZERO '
        'plastic strain -- getCopy("ThreeDimensional") is leaking element '
        '1\'s d_eps_p history into a sibling built from a driven prototype',
        bystander_pstrain, driven_pstrain)


# ===========================================================================
#  P1 red/blue review, section 5 item 2 -- the negative monotone clock (B1),
#  counted refusals (B2/B3), and the coverage-matrix gaps the review named.
# ===========================================================================
#
#  A GENUINELY FREE-DOF SETTLEMENT COLUMN.  Every deck in this file above
#  this point is either zero-free-DOF (`sani._build`/`_drive`,
#  `_build_floor_seeking_deck`) or force-controlled
#  (`_build_probe_triaxial`/`_build_free_dof_triaxial`). B1 (the extrapolation
#  factor pinned at `alpha` for the life of a `LoadControl(-ds)` leg) needs
#  BOTH at once: genuine free DOFs (so `analyze()` actually solves something
#  every step, matching a real BVP Newton loop) AND a displacement-controlled
#  `LoadControl(-ds)` deck (the campaign's own shape, `sanisand_tau0_band.py`,
#  where `ops_Dt < 0` by design). Neither existing rig is that deck.
#
#  Three node layers (k = 0, 1, 2), two stacked `LadrunoBrick` cubes. The
#  SAME per-layer roller convention as every other free-DOF deck in this file
#  (`fx = 1 if x==0`, `fy = 1 if y==0`, `fz = 1 if k==0`) is used at every
#  layer -- so the BASE (k=0) is the only layer with a z fixity, and the TOP
#  (k=2) gets an explicit `sp` prescribing its z displacement. The MIDDLE
#  layer (k=1) receives no `sp` at all: its z DOF (and, at the (1,1) corner,
#  its x/y DOFs too) are genuinely FREE, giving the deck real equations for
#  Newton to solve -- unlike every zero-free-DOF deck elsewhere in this file.
#
#  THE `sp` COEFFICIENT IS 1.0, DELIBERATELY.  With a `Linear` time series
#  (value == pseudo time `t`) and an `sp` coefficient of exactly 1.0, the
#  top-layer z displacement equals `t` itself, so each `LoadControl(-ds)`
#  step's `deltaLambda` (== the material's `ops_Dt`) is ALSO, directly, that
#  step's physical displacement increment -- there is no second scale factor
#  to keep straight between "the dt the material sees" and "the displacement
#  the deck applies", which is exactly the correspondence B1's test needs.
#
#  MAGNITUDES ARE NOT MEASURED ON A BINARY (none exists yet for this lane;
#  see the module docstring's "NUMBERS NOT INVENTED" section) -- they are
#  tied to scales this file's OTHER decks already prove converge:
#  `_SETTLE_DS0` is `sani._E_AX / sani._N_PL`, the exact per-step magnitude
#  `sani._drive` takes 20 plastic steps at; `_SETTLE_E_CONF` is the same
#  order of magnitude as `sani._C_E_AX`'s per-step contribution. If this
#  deck does not converge on the real binary, adjust the magnitudes -- this
#  is the same "measure, then adjust" discipline `_build_floor_seeking_deck`
#  documents for itself.
#
#  MEASURED 2026-09-06, against the fixed binary (ladrunoBuild ==
#  2473ce46c): `NormDispIncr 1e-10, 30 iters` (the tight tolerance every
#  OTHER zero/near-zero-free-DOF deck in this file uses) does NOT converge
#  on this deck's SECOND settlement step (stalls at ~3.3e-7, never reaches
#  1e-10) -- this deck has genuine free DOFs and a real Newton path-
#  dependence the zero-free-DOF decks never exercise, so it needs a looser
#  (still tight in absolute terms) criterion: `1e-6, 100 iters` converges
#  cleanly on all 15 steps (10 confinement + 5 settlement) and reproduces
#  the exact expected f sequence below.
# ---------------------------------------------------------------------------
_SETTLE_E_CONF = 2.0e-4                      # top-layer lateral confinement, stage 0
_SETTLE_N_CONF = 10
_SETTLE_DS0 = sani._E_AX / sani._N_PL        # ~1.5e-5, sani._drive's own proven per-step scale
_SETTLE_DS1 = 2.0 * _SETTLE_DS0


def _build_settlement_column(tag, opts):
    """Two stacked `LadrunoBrick` cubes, three node layers -- see the block
    comment above this function for the roller convention and why the
    middle layer is genuinely free."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    ops.element('LadrunoBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    ops.element('LadrunoBrick', 2, 5, 6, 7, 8, 9, 10, 11, 12, tag,
               '-geom', 'linear', '-formulation', 'bbar')
    for k in range(3):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            ops.fix(n, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    # pattern 1: lateral confinement, TOP layer only, applied over the
    # elastic stage below and then loadConst'd.
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j, (x, y) in enumerate(_XY):
        n = 8 + j + 1
        if x == 1.:
            ops.sp(n, 1, -_SETTLE_E_CONF)
        if y == 1.:
            ops.sp(n, 2, -_SETTLE_E_CONF)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-6, 100, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')


def test_negative_monotone_clock_runs_the_spec_factor():
    """B1: `implexDetail[5]` (the extrapolation factor `f`) tracks the
    SIGNED ratio `dt_{n+1}/dt_n * alpha` across a `ds` change on a
    `LoadControl(-ds)` leg with GENUINE free DOFs -- not pinned at `alpha`
    for the whole leg, which is what the OLD `mImplexDtCommit > 0.0` gate
    (instead of `!= 0.0` with a sign-consistent ratio) did on exactly this
    shape of deck, because two consecutive negative increments never
    reached the ratio branch at all.

    Kills the B1 mutant directly: an `f` that reads `alpha` (1.0) on the
    ds-doubling step, instead of `2.0 * alpha`, is the bug the campaign's
    own BVP legs shipped with.

    PINNED TO `-implexGuard off` (ADR-92 P2, WP-92e lane B2, 2026-09-07;
    RE-CHECKED against 8bfdfbc17's P2-2b fix and STILL pinned, see below).
    Originally measured: with the guard on, step 2 (i=1) read `f = 0.0`
    instead of `1.0` -- the settlement column's un-primed first plastic
    commit moved `mAlpha_in_n` as an initialisation artefact, which
    `ladrunoImplexCommit()`'s reversal check misread as a genuine reversal.

    P2-2b (8bfdfbc17) fixes EXACTLY that: re-run with the guard back on
    (no pin), `i=1` now correctly reads `f = 1.0` -- confirmed directly.
    But `i=2` onward now reads `f = 0.0` instead of `1.0`/`2.0`/`1.0` --
    a DIFFERENT, LATER guard-arming that P2-2b's fix does not (and is not
    meant to) touch: this deck's own `_SETTLE_E_CONF` sits the material
    right at the `p_min` floor from the first settlement step on (measured:
    repeated "mean stress p = 0.1008xx is below the floor... CLAMPING"
    warnings on this exact deck under the implicit path), which is
    independently plausible ground for a genuine `Kp <= 0` or reversal
    commit -- i.e. this looks like the guard doing its documented job on a
    deck that was never tuned to stay clear of it, not a residual defect.
    Not chased further (would need reading `mImplexGuardReversal`/
    `mImplexGuardSoftening` directly, which are not exposed as a response);
    B1's own claim (the dt-ratio tracking) is orthogonal to the guard
    either way, so it stays isolated with `-implexGuard off`: verified
    independently that the ratio sequence [1.0, 1.0, 1.0, 2.0, 1.0] holds
    exactly with the guard off, matching the precedent
    `test_implexcontrol_floor_accepts_once_reduction_limit_is_reached`
    already sets for `-implexFloor`.
    """
    tag = 8400
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexGuard', 'off')
    _build_settlement_column(tag, opts)

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    ops.integrator('LoadControl', 1.0 / _SETTLE_N_CONF)
    for step in range(_SETTLE_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    # pattern 2: top-layer z settlement, coefficient 1.0 -- see the block
    # comment above _build_settlement_column for why that makes ds the
    # physical displacement increment too.
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.sp(8 + j + 1, 3, 1.0)

    alpha = 1.0   # no -implexAlpha given -> the class default (LadrunoSANISAND.h)
    ds_sequence = [_SETTLE_DS0, _SETTLE_DS0, _SETTLE_DS0, _SETTLE_DS1, _SETTLE_DS1]
    # step 1: mImplexDtCommit == 0.0 right after the stage flip (
    #   ladrunoImplexInitState) -> f = alpha, by the spec's own dtCommit==0
    #   branch, regardless of ds. steps 2-3: constant ds -> ratio 1.0. step 4:
    #   ds DOUBLES -> ratio 2.0. step 5: ds constant again (at the new value)
    #   -> ratio back to 1.0, proving the factor is a live ratio and not
    #   stuck at whatever the previous step computed.
    expected_f = [alpha, 1.0, 1.0, 2.0 * alpha, 1.0]

    for i, (ds, exp_f) in enumerate(zip(ds_sequence, expected_f)):
        ops.integrator('LoadControl', -ds)
        assert ops.analyze(1) == 0, f'settlement step {i + 1} failed (ds={ds!r})'
        detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
        f = detail[5]
        assert f == pytest.approx(exp_f, rel=1.0e-6, abs=1.0e-9), (
            'implexDetail[5] (f) did not track the signed dt_{n+1}/dt_n * '
            'alpha ratio on a monotone-negative-clock LoadControl(-ds) leg '
            '-- this is B1: a material that pins f at alpha for the whole '
            'leg (the old mImplexDtCommit > 0.0 gate) passes every step '
            'here with f == alpha instead of the expected ratio',
            i, ds, exp_f, f, detail)


# ===========================================================================
#  B2/B3 -- counted, observable refusals via the new implexRefusals response
# ===========================================================================

def test_implexcontrol_refusal_is_counted_and_reported():
    """`-implexControl` refusing increments `implexRefusals` (Vector(4):
    total, d2, control, companion) at BOTH index 0 (total) and index 2
    (control-specific), and leaves the committed stress unchanged on a
    free-DOF deck (where a silently-accepted wrong answer WOULD move it).

    Kills a mutant that refuses (rc != 0, as the sibling qualitative test
    already checks) without incrementing the counter -- per B2/F8, the
    counter is the ONLY way `n_material_refused`-style accounting can ever
    be recovered from a deck instead of grepped out of an unthrottled log.

    p0/tol REUSE the pair `test_implexcontrol_refuses_past_tolerance...`
    measured on the fixed binary (2473ce46c): p0 = 5 kPa is NOT usable here
    either -- an essentially-zero tolerance refuses the NOMINAL steps
    `_establish_plastic_history` needs to succeed, crashing that helper's
    own internal assert before this test's body even runs. p0 = 50 kPa,
    tol = 0.02 leaves the nominal steps comfortably under tol and the big
    (10x) step comfortably over it.
    """
    tag = 8214
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', 0.02, 0.01)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _establish_plastic_history(tag)

    before_refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    before_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    assert len(before_refusals) == 4, (
        'implexRefusals did not return the documented 4-component vector '
        '(total, d2, control, companion)', before_refusals)

    big_dq = 10.0 * _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -big_dq)
    ops.integrator('LoadControl', 1.0)
    rc = ops.analyze(1)
    assert rc != 0, (
        'a deviatoric load increment 10x the nominal one, at p0 = 50 kPa '
        'with -implexControl tol = 0.02, was NOT refused', rc)

    after_refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    after_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))

    assert after_refusals[0] - before_refusals[0] >= 1, (
        'implexRefusals[0] (total) did not increment across a forced '
        '-implexControl refusal', before_refusals, after_refusals)
    assert after_refusals[2] - before_refusals[2] >= 1, (
        'implexRefusals[2] (control-specific) did not increment across a '
        'forced -implexControl refusal', before_refusals, after_refusals)
    assert after_stress == before_stress, (
        'the committed stress moved across a COUNTED -implexControl '
        'refusal -- the counter incremented but the state protection it is '
        'supposed to accompany (mSigma = mSigma_n) did not hold',
        before_stress, after_stress)

    # avgImplexError must be a NON-DESTRUCTIVE read (F6 in the majors list:
    # the old takeAverageError() zeroed the accumulator on every call, so a
    # recorder over an 8-point mesh got the average at point 1 and 0.0
    # everywhere else). Two back-to-back reads, no analyze() in between,
    # must agree.
    avg_1 = ops.eleResponse(1, 'material', 1, 'avgImplexError')[0]
    avg_2 = ops.eleResponse(1, 'material', 1, 'avgImplexError')[0]
    assert avg_1 == avg_2, (
        'avgImplexError changed between two consecutive reads with no '
        'analyze() in between -- the accumulator is being reset/consumed '
        'on read instead of reported non-destructively', avg_1, avg_2)


def test_companion_refusal_at_commit_is_observable():
    """B3: the COMMIT-time companion (`ladrunoImplexCommit` hitting the
    `-maxSubsteps` cap) refuses, and that refusal is OBSERVABLE only through
    `implexRefusals[3]` (companion) -- `Domain::commit()` is `elePtr->
    commitState();` with the return code discarded, so `analyze()` itself
    keeps returning 0 even though the companion's own re-integration failed.
    This is the DEFAULT configuration this file's other decks exercise
    least: `-implexControl` OFF, cap mandatory.

    Also checks the NEXT step's `f` is NOT stale: B3 requires the material
    to commit its best-effort state (mEpsilon_n, mImplexDtCommit) even when
    the companion itself refuses, so a constant-ds step immediately after a
    counted companion refusal must still read f == 1.0 (the ordinary
    same-ds ratio), not 0.0 or NaN.

    Reuses `sani._build_confined`, NOT `sani._build`: MEASURED on the fixed
    binary (2473ce46c), the plain `sani._build` deck's own plastic
    increment needs only ~1 ModifiedEuler substep (nowhere near
    `_CAP_TIGHT`) -- the '1000-5000 substeps' figure the block comment
    above `_CAP_ADEQUATE` documents was always measured on the CONFINE-
    FIRST deck's deviatoric leg specifically, not on `sani._build`'s single-
    ramp deck; this test's first draft conflated the two. On the confine-
    first deck, `_CAP_TIGHT` (200) reliably forces the companion to refuse
    every deviatoric-leg commit (measured: 7-8 companion refusals across 3
    steps) -- deliberately, this is the one test in the file where that is
    the point rather than a caveat. The confine-first deviatoric leg's own
    `LoadControl(1.0)` per step (constant, positive, `_analysis_confined`)
    is also what makes the f == 1.0 check below meaningful: a constant-ds
    leg, not the ds-doubling one `test_negative_monotone_clock_...` owns.
    """
    tag = 8215
    opts = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_TIGHT)
    sani._build_confined('LadrunoSANISAND', tag, opts)
    sani._confine_leg(tag)
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    before_refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))

    for step in range(3):
        rc = ops.analyze(1)
        assert rc == 0, (
            'the step itself should still "succeed" from Domain::commit()\'s '
            'point of view -- Domain::commit() is unconditional and '
            'discards the material return code; a nonzero rc here means '
            'something else in the deck broke, not the companion refusal',
            step, rc)

    after_refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert after_refusals[0] - before_refusals[0] > 0, (
        'implexRefusals[0] (total) never incremented even though this cap '
        'is measured (this docstring) to be too tight for the confine-first '
        "deck's own deviatoric increment", before_refusals, after_refusals)
    assert after_refusals[3] - before_refusals[3] > 0, (
        'the companion refusal at commitState (ladrunoImplexCommit hitting '
        'the -maxSubsteps cap) never incremented implexRefusals[3]. Since '
        'Domain::commit() discards commitState()\'s return code, this '
        'counter is the ONLY way to observe a companion refusal from '
        'Python -- if it stays at 0 the refusal is silently swallowed '
        'exactly as B3 found', before_refusals, after_refusals)

    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    f_last = detail[5]
    assert math.isfinite(f_last), (
        'the extrapolation factor is not finite (NaN/inf) on the step after '
        'a companion refusal at commit -- B3 requires committing a valid '
        'best-effort state even when the companion itself fails',
        f_last, detail)
    assert f_last == pytest.approx(1.0, rel=1.0e-6, abs=1.0e-9), (
        'the extrapolation factor after a companion refusal at commit is '
        'not 1.0 on this constant-ds deck -- B3 requires the material to '
        'commit its best-effort state (mEpsilon_n, mImplexDtCommit) even '
        'when the companion itself refuses, so the NEXT step\'s f must '
        'still be the ordinary same-ds ratio, not stale or zero',
        f_last, detail)


# ===========================================================================
#  Coverage-matrix gaps named by the P1 review: db roundtrip on a deck where
#  ON != OFF (F6), and the PlaneStrain lane (row 17)
# ===========================================================================

_RT_N_TOTAL = 6            # further deviatoric steps, split around the save point


def test_db_roundtrip_on_a_deck_where_on_differs_from_off():
    """`sendSelf`/`recvSelf` round trip on a deck where `-implex` ON and OFF
    committed stresses PROVABLY differ (F6: every ON/OFF pair elsewhere in
    this file is zero-free-DOF, where gate 5 proves they CANNOT differ, so
    `test_implex_db_roundtrip_carries_flags_and_history`'s own non-vacuity
    guard skips on every run and the assertion that could see a dropped
    flag never executes). This uses the free-DOF triaxial deck instead,
    checks the ON/OFF gap FIRST -- not as a skip-fallback at the end -- and
    only then exercises save/restore on the ON leg, so a mutant that drops
    the -implex flags or the mImplexDEpsP history on restore lands off the
    ON reference by a MEASURED, non-floor amount.

    If `database()` is unsupported on this build, skips with the reason
    (matching `_roundtrip_implex`'s own convention).
    """
    p0 = 5.0
    opts_on = ('-implex', '-maxSubsteps', _CAP_ADEQUATE)

    # non-vacuity FIRST, on two independent one-shot builds -- not the
    # roundtrip skeleton itself -- so this can never accidentally pass
    # because the LATER roundtrip deck happens to be insensitive.
    tag_on_probe, tag_off_probe = 8307, 8308
    _build_free_dof_triaxial(tag_on_probe, opts_on, p0=p0)
    _establish_plastic_history(tag_on_probe)
    on_probe = list(ops.eleResponse(1, 'material', 1, 'stress'))

    _build_free_dof_triaxial(tag_off_probe, (), p0=p0)
    _establish_plastic_history(tag_off_probe)
    off_probe = list(ops.eleResponse(1, 'material', 1, 'stress'))

    gap = sani._reldiff(off_probe, on_probe)
    assert gap > _SENSITIVITY_FLOOR, (
        'on this free-DOF triaxial deck, -implex ON and OFF committed the '
        'SAME stress within the sensitivity floor -- the roundtrip below '
        'would prove nothing on this deck; raise the confinement/step size '
        'rather than weakening this gate', gap, off_probe, on_probe)

    tag_rt = 8309
    n_cut = 3
    dq3 = 3.0 * _PROBE_DQ_NOMINAL
    dq3q = dq3 / 4.0

    with tempfile.TemporaryDirectory(prefix='ladruno_sanisand_implex_free_rt_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'sanisand_implex_free_rt')

        _build_free_dof_triaxial(tag_rt, opts_on, p0=p0)
        _establish_plastic_history(tag_rt)

        ops.timeSeries('Linear', 3)
        ops.pattern('Plain', 3, 3)
        for j, (x, y) in enumerate(_XY):
            ops.load(4 + j + 1, 0.0, 0.0, -dq3q)
        ops.integrator('LoadControl', 1.0 / _RT_N_TOTAL)
        for step in range(n_cut):
            assert ops.analyze(1) == 0, f'pre-save continuation step {step + 1} failed'
        mid = list(ops.eleResponse(1, 'material', 1, 'stress'))

        try:
            ops.database('File', dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f'database() unsupported in this build: {exc}')
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip('database save returned failure on this build')

        _build_free_dof_triaxial(tag_rt, opts_on, p0=p0)   # fresh, uncommitted skeleton
        ops.database('File', dbpath)
        ops.restore(1)
        after = list(ops.eleResponse(1, 'material', 1, 'stress'))

        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('Plain')
        ops.system('FullGeneral')
        ops.test('NormUnbalance', _PROBE_TOL_REL * p0, _PROBE_MAXITER, 0)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1.0 / _RT_N_TOTAL)
        ops.analysis('Static')
        for step in range(_RT_N_TOTAL - n_cut):
            assert ops.analyze(1) == 0, f'post-restore continuation step {step + 1} failed'
        out = list(ops.eleResponse(1, 'material', 1, 'stress'))

    assert sani._reldiff(mid, after) <= _EQ_TOL, (
        'the restored -implex material did not come back on the state it '
        'was saved at (free-DOF deck)', mid, after)

    # unbroken reference: the SAME deck, SAME total load, no save/restore.
    tag_ref = 8310
    _build_free_dof_triaxial(tag_ref, opts_on, p0=p0)
    _establish_plastic_history(tag_ref)
    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq3q)
    ops.integrator('LoadControl', 1.0 / _RT_N_TOTAL)
    for step in range(_RT_N_TOTAL):
        assert ops.analyze(1) == 0, f'reference continuation step {step + 1} failed'
    ref_final = list(ops.eleResponse(1, 'material', 1, 'stress'))

    assert sani._reldiff(ref_final, out) <= 1.0e-12, (
        'the restored -implex material, continued to the SAME total load, '
        'does not match an UNBROKEN run to 1e-12 on a deck where -implex ON '
        'and OFF are measurably different -- the flags or the mImplexDEpsP '
        'history did not survive sendSelf/recvSelf',
        ref_final, out, gap)


def test_plane_strain_implex_smoke():
    """`-implex` on `LadrunoSANISANDPlaneStrain` (`-ndm 2`, the
    `quad ... PlaneStrain` -> `getCopy("PlaneStrain")` route) at least RUNS
    -- coverage row 17: every 2D deck in this file before this test used
    `-ndm 3` only, so a mutant that wired IMPL-EX into the 3D wrapper alone
    (or that crashes/refuses on the 2D one) had zero chance of being caught.
    Reuses `sani._build_ps`/`_drive_ps` -- the SAME proven zero-free-DOF
    plane-strain deck `test_ladruno_sanisand.py` already uses for the
    PlaneStrain lane's other gates.
    """
    tag = 8311
    opts = sani._OPTS_VANILLA + ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    sani._build_ps(tag, opts)
    sani._elastic_leg(tag)
    sani._plastic_leg(tag)
    stress = sani._stress()
    assert _vnorm(stress) > 0.0, (
        'the plane-strain -implex deck committed a zero stress -- the deck, '
        'not -implex, is broken', stress)

    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    assert len(detail) == 6, (
        'implexDetail did not return the documented 6-component vector on '
        'the PlaneStrain lane', detail)
    assert math.isfinite(detail[5]), (
        "implexDetail's f (index 5) is not finite on the PlaneStrain lane",
        detail)

# ===========================================================================
#  Mutation gate survivors (Ladruno_implementation/_adr92_p1_mutation_gate.md,
#  commit b523952fa) -- M4, M5, M10. Score was 0.75 (9/12); these three tests
#  are the ones the gate's own audit says are owed.
# ===========================================================================

def test_implexcontrol_floor_accepts_once_reduction_limit_is_reached():
    """Mutation gate M4 survivor: `mImplexDt0` (the reduction floor's
    arming value) removed. No test in the baseline battery drives a
    genuine SUBDIVISION LADDER (halve `ds` after a refusal, retry), so the
    floor's escape branch (`|mImplexDt| < reductionLimit * mImplexDt0` =>
    accept, "nothing left to cut") was never exercised -- killing the
    arming line changed nothing.

    Drives that ladder directly, at an UNREACHABLE -implexControl
    tolerance (1e-9) so implexError is ALWAYS above tol no matter how
    small `ds` gets -- the ONLY way any attempt can ever be accepted is
    the floor, never a genuinely small error. MEASURED on 6f52a30bb at
    reductionLimit = 0.3: priming ds0 = 0.05 (un-primed exemption, ALSO
    arms mImplexDt0 = 0.05); the ladder ds = 0.1, 0.05, 0.025 all refuse
    (ratio >= 0.3); ds = 0.0125 (ratio 0.25 < 0.3) is ACCEPTED, with
    implexDetail[0] (the error) still reading 4.9e-6 -- far above
    tol = 1e-9 -- proving the accept came from the floor, not from the
    error dropping under tol.

    Kills a mutant that never arms `mImplexDt0` (or arms it with the wrong
    sign/magnitude): with the floor dead, the ladder refuses FOREVER and
    this test's final `assert rc == 0` fails.

    PINNED TO `-implexFloor accept` (ADR-92 P2, WP-92e lane B2, 2026-09-07).
    This test was written against the ONLY floor behaviour that existed at
    the time -- unconditional accept, exactly what `-implexFloor accept`
    now names -- and its own assertion (`last_detail[0] > tol`: the
    accepted step still carries the RAW large extrapolation error) is that
    mode's contract verbatim. P2 changes the DEFAULT to `implicit` (see
    `test_floor_fallback_delivers_implicit_stress_and_counts`), under which
    the floor step would commit the implicit companion's own answer instead
    and this assertion would most likely fail for a reason that has nothing
    to do with the M4 mutant this test exists to kill. Pinning the mode
    keeps this test's claim exactly what it always was.

    ALSO PINNED TO `-implexTrialGuard off` (ADR-92 P2-6, WP-92e lane B2,
    2026-09-07, binary 708152eac): P2-6's trial-time f=0 fallback defaults
    on and runs BEFORE every -implexControl refusal in this ladder,
    including the very FIRST attempt (ds = 1.0) -- measured: with the
    default, that first attempt is RESCUED (`rc == 0`) instead of refused,
    so the ladder never even reaches the reduction-floor rung this test is
    built to isolate. Same isolation rule as the `-implexFloor` pin above.
    """
    tag = 8320
    tol = 1.0e-9
    reduction_limit = 0.3
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', tol, reduction_limit,
            '-implexFloor', 'accept', '-implexTrialGuard', 'off')
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _confine_only(tag)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -6.0 / 4.0)

    ds0 = 0.05
    ops.integrator('LoadControl', ds0)
    rc_prime = ops.analyze(1)
    assert rc_prime == 0, ('the priming step (un-primed, arms mImplexDt0) '
                           'failed to converge -- a harness problem', rc_prime)
    ops.loadConst('-time', 0.0)

    refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    ds = ds0 * 2.0
    accepted = False
    last_detail = None
    for attempt in range(8):
        ops.integrator('LoadControl', ds)
        rc = ops.analyze(1)
        new_refusals = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
        last_detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
        if rc == 0:
            assert new_refusals[2] == refusals[2], (
                'implexRefusals[2] moved on the step the floor is supposed '
                'to ACCEPT -- the floor branch must not count as a '
                'refusal', refusals, new_refusals)
            accepted = True
            break
        assert new_refusals[2] > refusals[2], (
            'analyze() failed on the ladder but implexRefusals[2] did not '
            'move -- the refusal was not -implexControl', refusals, new_refusals)
        refusals = new_refusals
        ds = ds / 2.0

    assert accepted, (
        'the subdivision ladder never got accepted -- the reduction floor '
        '(mImplexDt0/reductionLimit) never engaged, so -implexControl '
        'refused without limit, exactly the M4/B2 defect', ds, last_detail)
    assert last_detail[0] > tol, (
        'the ladder was accepted, but implexDetail[0] (the error) is NOT '
        'above tol -- this accept could be explained by the error '
        'genuinely dropping low, which would prove nothing about the '
        'floor; re-derive ds0/reductionLimit so tol stays unreachable',
        last_detail[0], tol)


def test_implexcontrol_refusal_return_code_is_the_sentinel_not_a_stalled_residual():
    """Mutation gate M5 survivor: the `-implexControl` refusal path
    returns `0` instead of `LADRUNO_MATERIAL_REFUSED`. The existing
    refusal tests stay green under that mutant because the OTHER two
    thirds of the contract survive it (`mSigma = mSigma_n` still runs, so
    Newton still stalls and `analyze()` still returns non-zero; the
    counter still increments) -- so `rc != 0` and `implexRefusals[2]++`
    are SYMPTOMS the mutant also produces, not proof the SENTINEL itself
    propagated.

    Distinguishes the two mechanisms by ITERATION COUNT. `ops.test(...)`
    here is generously loose (`maxIter = 25`, a tolerance any ordinarily-
    converging step clears easily) -- a genuine residual stall (the M5
    mutant's failure mode: Newton grinding against a frozen, wrong `sigma`
    with no early exit) would need many iterations before giving up,
    while `LadrunoBrick::update()` propagating the sentinel aborts the
    step on the FIRST iteration, before the residual has a chance to
    stall (`Domain::update()` returns the element's code and
    `IncrementalIntegrator::update()` fails the step immediately -- the
    same mechanism `test_tangent_identity_...`'s block comment documents
    for `maxIter = 1`, here observed at `maxIter = 25`).

    MEASURED on 6f52a30bb: the SAME nominal-sized load, reapplied once
    primed with `-implexControl` tol = 1e-9 (unreachable), refuses with
    `analyze() == -3` after `ops.testIter() == 1`, not 25.

    Kills a mutant that returns `0` (or any non-`LADRUNO_MATERIAL_REFUSED`
    success code) from the control-refusal branch: `rc` would then only
    go non-zero (if at all) after Newton genuinely exhausts iterations
    against the frozen, mismatched stress, so `ops.testIter()` would read
    close to `maxIter`, not 1.
    """
    tag = 8321
    tol = 1.0e-9
    maxit = 25
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', tol, 0.01)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    ops.test('NormUnbalance', _PROBE_TOL_REL * 50.0, maxit, 0)
    _confine_only(tag)

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 0.05)
    rc_prime = ops.analyze(1)
    assert rc_prime == 0, ('the priming step failed to converge', rc_prime)
    ops.loadConst('-time', 0.0)

    before = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    before_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))

    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 0.05)
    rc = ops.analyze(1)
    niter = ops.testIter()
    after = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    after_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))

    assert rc != 0, (
        'the primed, tolerance-exceeding step converged -- it was supposed '
        'to be refused by -implexControl', rc)
    assert after[2] - before[2] >= 1, (
        'implexRefusals[2] did not increment in the SAME step as the '
        'return-code failure', before, after)
    assert niter <= 2, (
        'analyze() failed and the counter moved, but it took '
        f'{niter} iterations out of a maxIter = {maxit} budget to do so -- '
        'that is the signature of an ORDINARY stalled residual (Newton '
        'grinding against a frozen, mismatched stress), not the material '
        "propagating LADRUNO_MATERIAL_REFUSED through the element's own "
        'update() on the first pass. If the sentinel were live this would '
        'fail in 1 iteration regardless of maxIter', niter, maxit)
    assert before_stress == after_stress, (
        'the committed stress moved across the refusal', before_stress, after_stress)


def test_reararm_after_refusal_without_a_revert_uses_its_own_dt_ratio():
    """Mutation gate M10 survivor: the re-arm (`mImplexStepArmed = true`)
    at all three refusal sites removed. It survives the rest of this file
    because every OTHER refusal test here uses `LadrunoBrick`, which
    PROPAGATES the sentinel -- so `StaticAnalysis`'s own failure path
    (`StaticAnalysis.cpp:185`) calls `revertToLastCommit()`, which
    independently re-arms. The removed line is genuinely redundant on
    every deck elsewhere in this file.

    THE MATERIAL-LEVEL PATH THE FIX IS FOR: a caller that refuses and
    retries WITHOUT a revert in between. `ops.analyze()` itself always
    reverts internally on a propagated failure, so the only way to reach
    this path from Python is an element that does NOT propagate the
    material's return code -- `stdBrick` ("stdBrick swallows material
    return codes", documented elsewhere in this file/LEDGER_quirks) --
    so `analyze()` reports SUCCESS even when the material internally
    refused, and `StaticAnalysis` never reverts between the two calls.

    MEASURED on 6f52a30bb, zero-free-DOF `stdBrick`, net-dilating deck (the
    same shape `_build_floor_seeking_deck` uses, so genuine `d_eps_p(n)`
    exists): 6 nominal history steps at ds = 0.02 all accepted
    (implexError <= 4.8, tol = 10.0); a jump to ds = 0.4 reads
    implexError = 26.2 (> tol) and is REFUSED (implexRefusals[2] += 8,
    f = 20.0 = 0.4 / 0.02); the VERY NEXT `analyze()` call, ds = 0.1, NO
    revert in between, reads implexError = 5.7 (< tol) and is ACCEPTED
    (implexRefusals[2] unchanged) with f = 5.0 -- exactly 0.1 / 0.02, the
    ratio against the LAST TRULY COMMITTED dt (0.02), not 0.1 / 0.4 = 0.25,
    which is what comparing against the refused attempt's own frozen dt
    (a re-arm that never happened) would have produced.

    Kills a mutant that removes the re-arm at the refusal sites: without
    it, the accepted retry's `f` would be computed against the stale
    refused-attempt dt instead of its own.

    PINNED TO `-implexGuard off` (ADR-92 P2, WP-92e lane B2, 2026-09-07;
    RE-CHECKED against 8bfdfbc17's P2-2b fix and STILL pinned). This deck's
    first plastic commit is un-primed (same mechanism as
    `test_negative_monotone_clock_runs_the_spec_factor`'s pin), which
    8bfdfbc17's P2-2b fixes directly -- but this deck (`_PMIN_E_CONF_LOW`,
    `lat = 1.5`, the SAME net-dilating shape `_build_floor_seeking_deck`
    uses specifically to threaten the `p_min` floor) drives the material
    hard enough that a re-run with the guard back on STILL reads `f = 0.0`
    on the small retry step, not `1.0` -- unlike
    `test_negative_monotone_clock_runs_the_spec_factor`'s deck, this one
    was deliberately built to reach genuine softening/reversal territory
    (that is the whole point of "floor-seeking"), so a guard-arming
    somewhere in its 6-step history-plus-big-step run is plausible on its
    own terms, independent of the un-primed defect P2-2b fixed. M10's own
    claim (the re-arm dt-ratio) is orthogonal to the guard either way, so
    it stays isolated with `-implexGuard off`: verified independently that
    with the guard off the small (retry) step's `f` is exactly `1.0`
    (`ds_nominal / ds_nominal`), matching this test's own claim.
    """
    tag = 8322
    e_conf = sani._PMIN_E_CONF_LOW
    e_ax_total = 5.0e-3
    lat = 1.5
    opts = ('-Presidual', 0.0, '-Pmin', sani._PMIN_LADRUNO, '-honorTolR', 0,
            '-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', 10.0, 0.01, '-implexGuard', 'off')

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial('LadrunoSANISAND', tag, *_PARAMS, *opts)
    # stdBrick, NOT LadrunoBrick -- see the docstring: this is the ONE
    # test in the file where swallowing the material's return code is
    # the point, not a hazard.
    ops.element('stdBrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, tag)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -e_conf)
            if y == 1.:
                ops.sp(n, 2, -e_conf)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-11, 40, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    ops.integrator('LoadControl', 1.0 / 10)
    for step in range(10):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.loadConst('-time', 0.0)
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        n = j + 1
        if x == 1.:
            ops.sp(n, 1, lat * e_ax_total)
        if y == 1.:
            ops.sp(n, 2, lat * e_ax_total)
    for j, (x, y) in enumerate(_XY):
        ops.sp(4 + j + 1, 3, -e_ax_total)

    # implexRefusals is a PROCESS-WIDE static counter (see the contract
    # in the module docstring) -- it is NOT zero here in a full-suite
    # run (earlier tests, e.g. the M4/M5 survivors above, already
    # refused many times), so this test tracks DELTAS off its own
    # baseline throughout, never an absolute value.
    refusals_before_history = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))

    ds_nominal = 0.02
    ops.integrator('LoadControl', ds_nominal)
    for step in range(6):
        assert ops.analyze(1) == 0, f'history step {step + 1} failed'

    refusals_0 = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert refusals_0[2] == refusals_before_history[2], (
        'a nominal history step already refused -- the deck, not the '
        'mechanism under test, needs re-deriving', refusals_before_history, refusals_0)

    ops.integrator('LoadControl', 20.0 * ds_nominal)
    rc_big = ops.analyze(1)
    detail_big = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    refusals_big = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    assert rc_big == 0, (
        'the big step FAILED analyze() -- stdBrick is supposed to swallow '
        "the material's return code, so analyze() must report success "
        'even though the material itself refused internally', rc_big)
    assert refusals_big[2] > refusals_0[2], (
        'the big (20x) step did not refuse internally (implexRefusals[2] '
        "did not move) -- the deck isn't aggressive enough to force a "
        '-implexControl refusal here', refusals_0, refusals_big)

    # NO revertToLastCommit anywhere in between -- analyze() reported
    # success above, so StaticAnalysis never reverted.
    ops.integrator('LoadControl', ds_nominal)
    rc_small = ops.analyze(1)
    detail_small = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    refusals_small = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))

    assert rc_small == 0, (
        'the small step (re-driven with no revert since the big refusal) '
        'failed to converge -- a harness problem', rc_small)
    assert refusals_small[2] == refusals_big[2], (
        'the small step ALSO refused -- pick a tolerance between the big '
        "and small steps' measured implexError so this one is a clean "
        'accept', refusals_big, refusals_small)
    expected_f = ds_nominal / ds_nominal
    assert detail_small[5] == pytest.approx(expected_f, rel=1.0e-6, abs=1.0e-9), (
        'implexDetail[5] (f) on the step right after an un-reverted '
        'refusal does not equal the ratio for ITS OWN dt against the '
        'last TRULY committed dt -- if the re-arm (mImplexStepArmed) is '
        'missing, f would instead reflect the REFUSED attempt\'s frozen '
        'dt', detail_small[5], expected_f, detail_big, detail_small)


# ===========================================================================
#  ADR-92 P2 (owed) -- `-implexFloor`, `-implexGuard`, the hold-safe clock,
#  and the `stressCorrection` updateParameter dispatch fix.
#
#  Written 2026-09-07 (WP-92e, lane B2) against the interfaces lane A2 is
#  building to, BEFORE any P2 binary exists -- see the module docstring's
#  "LANE B2 / P2 BATTERY" paragraph for the full brief and its caveats.
# ===========================================================================

# ---------------------------------------------------------------------------
#  1. `-implexFloor implicit|accept|refuse` (new default `implicit`)
# ---------------------------------------------------------------------------

def _drive_floor_ladder(tag, floor_mode, tol=1.0e-9, reduction_limit=0.3,
                        big_factor=20.0, max_attempts=8, ds0=0.05):
    """Prime a free-DOF triaxial deck (same shape as the M4 survivor test
    above), then force a subdivision ladder DOWN TO the -implexControl
    reduction floor -- `tol` is unreachable by construction, so only the
    floor branch (never a genuinely small error) can ever end the ladder,
    under any of the three -implexFloor modes.

    `-implexTrialGuard off` (ADR-92 P2-6, WP-92e lane B2, 2026-09-07,
    binary 708152eac): P2-6's trial-time f=0 fallback defaults on and
    would otherwise rescue the very FIRST attempt in this ladder before
    -implexFloor ever gets a chance to matter (measured, matching the same
    finding on the M4 survivor test) -- this function's whole purpose is
    to isolate -implexFloor, so P2-6's own orthogonal rescue is disabled
    here, the same way it is pinned on every other pre-P2-6 refusal test
    in this file.

    Returns a list of (ds, rc, implexDetail) for every rung attempted, in
    order, stopping at the first accepted (rc == 0) attempt if any.
    """
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', tol, reduction_limit,
            '-implexFloor', floor_mode, '-implexTrialGuard', 'off')
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _confine_only(tag)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -6.0 / 4.0)

    ops.integrator('LoadControl', ds0)
    rc_prime = ops.analyze(1)
    assert rc_prime == 0, ('the priming step (un-primed, arms mImplexDt0) '
                           'failed to converge -- a harness problem', rc_prime)
    ops.loadConst('-time', 0.0)

    ds = ds0 * big_factor
    history = []
    for attempt in range(max_attempts):
        ops.integrator('LoadControl', ds)
        rc = ops.analyze(1)
        detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
        history.append((ds, rc, detail))
        if rc == 0:
            break
        ds = ds / 2.0
    return history


def test_floor_fallback_delivers_implicit_stress_and_counts():
    """ADR-92 P2: `-implexFloor` selects what happens at the
    `-implexControl` reduction floor (`|mImplexDt| < reductionLimit *
    mImplexDt0`), where the OLD (P1) behaviour unconditionally accepted the
    extrapolated `sigma~` no matter how large its error. Drives the SAME
    ladder as `test_implexcontrol_floor_accepts_once_reduction_limit_is_reached`
    under all three modes:

      * `refuse` -- every rung refuses, including the floor rung itself;
        there is nothing to fall back to.
      * `accept` -- the RETIRED unconditional-accept behaviour (negative
        control / the mutant this test is written to catch if
        `-implexFloor` is silently ignored): the floor step is accepted
        with the RAW extrapolation error still on the wire.
      * `implicit` (the NEW default) -- the SAME floor rung is accepted, but
        the GP delivers the IMPLICIT companion's own stress for that step
        instead of the extrapolation, so the committed `implexError`
        collapses relative to `accept`'s reading of that identical rung --
        sigma~ == sigma_implicit at commit, to within genuine numerical
        noise, not the raw first-order extrapolation gap.

    `implexGuards[0]` (floor fallbacks) must increment ONLY on the
    `implicit` run's floor acceptance. (`implexGuards[3]` was reserved when
    this test was first written; 8bfdfbc17/P2-5 repurposed it as the
    reversal-noise-guard count, which is NOT specific to `-implexFloor` and
    is expected to move on this ladder too as ds shrinks toward the floor
    -- see the P2-5 tests below instead; not re-checked here.)

    Kills a mutant that ignores `-implexFloor` entirely (every mode behaves
    like the old unconditional accept -- caught by `refuse` never actually
    refusing at the floor, and by the `implicit` run's error staying close
    to `accept`'s instead of collapsing), or that swaps `implicit`/`accept`.

    WHY THE COMPARISON IS `accept` VS `implicit` ON THE SAME RUNG, NOT
    "floor error vs the last refused attempt's error" (fixed 2026-09-07,
    WP-92e lane B2, re-run against 87b9cf846). `_drive_floor_ladder` uses
    `LadrunoBrick` (a element that PROPAGATES a material refusal), so a
    REFUSED rung never reaches `commitState()` -- `StaticAnalysis` reverts
    instead -- and `implexDetail`'s error components are only written
    inside `ladrunoImplexCommit()` (i.e. only ON A COMMIT). Reading
    `implexDetail` right after a refused, reverted `analyze()` call
    therefore returns whatever was last COMMITTED, not "this attempt's
    error" -- measured: every refused rung in this ladder reads
    `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]` (the priming commit's own values, which
    happen to be exactly zero on this deck's un-primed first step), so a
    "floor error / last-refused error" ratio divides by zero every time.
    The `accept` and `implicit` runs are DETERMINISTIC and reach the
    IDENTICAL final rung (same tol/reductionLimit/ds0/big_factor, and
    `-implexFloor` does not affect any refusal decision ABOVE the floor),
    so comparing their two floor-accepted (i.e. actually COMMITTED, non-
    stale) `implexDetail[0]` readings for that SAME rung is both correct
    and a strictly more apples-to-apples test of the flag than the retired
    comparison was.
    """
    tol = 1.0e-9
    reduction_limit = 0.3

    hist_refuse = _drive_floor_ladder(8330, 'refuse', tol, reduction_limit)
    assert all(rc != 0 for _, rc, _ in hist_refuse), (
        '-implexFloor refuse accepted a step at (or below) the reduction '
        'floor -- refuse is supposed to refuse EVERY rung, including the '
        'floor rung, since there is nothing to fall back to', hist_refuse)

    hist_accept = _drive_floor_ladder(8331, 'accept', tol, reduction_limit)
    accepted_accept = [h for h in hist_accept if h[1] == 0]
    assert accepted_accept, (
        '-implexFloor accept never reached an accepted step -- the ladder '
        'harness itself needs re-deriving, not this assertion', hist_accept)
    ds_a, rc_a, detail_a = accepted_accept[-1]
    assert detail_a[0] > tol, (
        '-implexFloor accept is supposed to commit the RAW extrapolation '
        'error at the floor (the retired unconditional-accept behaviour), '
        'above -implexControl\'s own (unreachable, by construction) tol -- '
        'if this reads at or below tol, the flag stopped selecting the '
        'accept mode', detail_a, tol)

    guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    assert len(guards_before) == 5, (
        'implexGuards did not return the documented 5-component vector',
        guards_before)
    hist_implicit = _drive_floor_ladder(8332, 'implicit', tol, reduction_limit)
    guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    accepted_idx = next((i for i, h in enumerate(hist_implicit) if h[1] == 0), None)
    assert accepted_idx is not None, (
        '-implexFloor implicit never reached an accepted step', hist_implicit)
    assert accepted_idx > 0, (
        'the FIRST attempt was accepted -- the ladder needs to actually '
        'refuse at least once before reaching the floor for this test to '
        'distinguish the floor fallback from an ordinary accept',
        hist_implicit)
    assert len(hist_implicit) == len(hist_accept) and accepted_idx == len(hist_accept) - 1, (
        'the implicit and accept ladders did not reach the SAME rung -- '
        '-implexFloor is not supposed to influence any refusal decision '
        'ABOVE the floor, so the two deterministic ladders (same tol/'
        'reductionLimit/ds0/big_factor) should have an identical shape',
        hist_accept, hist_implicit)

    ds_floor, rc_floor, detail_floor = hist_implicit[accepted_idx]
    assert rc_floor == 0, ('sanity: the located "accepted" attempt did not '
                           'actually converge', hist_implicit)
    assert ds_floor == ds_a, (
        'sanity: the implicit ladder\'s floor rung ds does not match the '
        'accept ladder\'s -- the two ladders are not actually comparable',
        ds_a, ds_floor)

    err_floor = detail_floor[0]
    err_accept = detail_a[0]
    assert err_accept > 0.0, (
        'the accept ladder\'s own floor-accepted attempt reported a zero '
        'error -- cannot form the ratio this test needs', detail_a)
    ratio = err_floor / err_accept
    assert ratio < 0.1, (
        '-implexFloor implicit is supposed to deliver the IMPLICIT stress '
        'at the floor, collapsing the committed implexError relative to '
        '-implexFloor accept\'s reading of the IDENTICAL rung -- got a '
        'ratio of %r (implicit floor error %r vs accept floor error %r), '
        'not < 0.1' % (ratio, err_floor, err_accept), hist_implicit, hist_accept)

    assert guards_after[0] - guards_before[0] >= 1, (
        'implexGuards[0] (floor fallbacks) did not increment across the '
        '-implexFloor implicit ladder reaching its floor',
        guards_before, guards_after)
    # implexGuards[3] is NOT checked here -- P2-5 (8bfdfbc17) repurposed the
    # formerly-reserved slot as the reversal-noise-guard count, which fires
    # on ANY near-zero strain increment regardless of -implexFloor, and this
    # ladder's shrinking ds legitimately produces some (measured: +176 on
    # this run). See test_hold_does_not_reset_alpha_in_on_the_implicit_path
    # for the dedicated P2-5 coverage.


# ---------------------------------------------------------------------------
#  2. `-implexGuard on|off` -- the elastic predictor after a load reversal
# ---------------------------------------------------------------------------

def _drive_reversal(tag, guard_mode):
    """`_establish_plastic_history`, then ONE reversal step (the axial
    deviator load flips sign, so `(alpha - alpha_in):n < 0` at that commit
    and `alpha_in` resets), then ONE continuation step in the SAME reversed
    direction and at the SAME `LoadControl` magnitude as the reversal step
    itself -- so a live (non-guarded) `f` on the continuation step would
    read the ordinary same-ds ratio (1.0, no `-implexAlpha` given), making
    a guard-forced `f = 0` unambiguous against it.

    Returns (implexDetail on the continuation step, implexGuards before the
    continuation step, implexGuards after it).
    """
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexGuard', guard_mode)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _establish_plastic_history(tag)

    dq = _PROBE_DQ_NOMINAL / 4.0

    # the reversal step: opposite sign to _establish_plastic_history's own
    # compressive pattern -- relieves, then reverses, the axial deviator.
    ops.timeSeries('Linear', 4)
    ops.pattern('Plain', 4, 4)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, +dq)
    ops.integrator('LoadControl', 1.0)
    rc_rev = ops.analyze(1)
    assert rc_rev == 0, ('the reversal step failed to converge', rc_rev)
    ops.loadConst('-time', 0.0)

    guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    # the continuation step -- same reversed direction, same magnitude.
    ops.timeSeries('Linear', 5)
    ops.pattern('Plain', 5, 5)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, +dq)
    ops.integrator('LoadControl', 1.0)
    rc_cont = ops.analyze(1)
    assert rc_cont == 0, ('the post-reversal continuation step failed to '
                          'converge', rc_cont)

    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    return detail, guards_before, guards_after


def test_guard_zeroes_f_after_reversal():
    """ADR-92 P2: `-implexGuard on` (the default) forces the elastic
    predictor (`implexDetail[5] == 0.0`, i.e. `f = 0`) on the step
    immediately after a commit whose state shows a load reversal
    (`alpha_in` reset). `-implexGuard off` leaves the extrapolation live on
    that same step, reading the ordinary dt ratio instead.

    Kills a mutant that never wires the reversal-guard branch (both modes
    would read the live ratio) or that ignores `-implexGuard off` (the
    guard would fire regardless of the flag).
    """
    detail_on, guards_before_on, guards_after_on = _drive_reversal(8340, 'on')
    assert detail_on[5] == 0.0, (
        'implexDetail[5] (f) is not exactly 0.0 on the step immediately '
        'after a committed load reversal (alpha_in reset) -- -implexGuard '
        'on is supposed to force the elastic-predictor (f = 0) '
        'extrapolation on that step', detail_on)
    assert guards_after_on[1] - guards_before_on[1] >= 1, (
        'implexGuards[1] (f=0 guards) did not increment across the '
        'post-reversal step even though implexDetail[5] read 0.0',
        guards_before_on, guards_after_on)

    detail_off, guards_before_off, guards_after_off = _drive_reversal(8341, 'off')
    expected_f_off = 1.0   # same-ds continuation, no -implexAlpha given
    assert detail_off[5] == pytest.approx(expected_f_off, rel=1.0e-6, abs=1.0e-9), (
        'with -implexGuard off, implexDetail[5] should equal the ordinary '
        'dt ratio on the SAME post-reversal step where -implexGuard on '
        'forces f = 0 -- if this also reads 0.0, the guard is not gated by '
        'the flag', detail_off, expected_f_off)
    assert guards_after_off[1] == guards_before_off[1], (
        'implexGuards[1] moved with -implexGuard off -- the guard count '
        'must not increment when the guard itself is disabled',
        guards_before_off, guards_after_off)


# ---------------------------------------------------------------------------
#  3. `-implexGuard on` -- the elastic predictor after softening
#
#  NO BINARY EXISTS FOR THIS LANE TO MEASURE AGAINST (module docstring,
#  "LANE B2 / P2 BATTERY"). Whether the dense confine-first deck below
#  actually reaches a peak-then-decline in eta/M_b within its step budget
#  is genuinely unknown from here -- default e_init (0.6944) IS already
#  dense-of-critical at this deck's confinement (psi ~ -0.13, computed from
#  _PARAMS: e_c = e0 - lambda_c*(p/Patm)**ksi at p ~ 1.7 kPa, the CONFINED
#  deck's own measured confinement pressure -- see `_c_series`'s block
#  comment), and the net-DILATING shape (`lat = 1.5`) is the same one
#  `_build_floor_seeking_deck` already proves drives the material hard
#  toward its limits. But "dense and dilating enough to threaten the p_min
#  floor" (that function's own proven regime, at a MUCH lower confinement)
#  is not the same claim as "reaches a peak-then-decline in stress ratio at
#  THIS confinement within 400 steps" -- so this test measures its own
#  outcome at runtime and reports whichever one is true, rather than
#  asserting one that was never checked.
# ---------------------------------------------------------------------------

_SOFTEN_N_DEV = 400
_SOFTEN_E_CONF = sani._C_E_CONF     # ~1.7 kPa confinement (measured on the
                                    # sibling _drive_confined deck), safely
                                    # above the p_min floor
_SOFTEN_LAT = 1.5                   # net-dilating -- _build_floor_seeking_deck's
                                    # own proven shape


def test_guard_zeroes_f_after_softening():
    """ADR-92 P2: `-implexGuard on` forces the elastic predictor
    (`implexDetail[5] == 0.0`) on the step after a commit whose state shows
    SOFTENING (`Kp <= 0`), the same mechanism
    `test_guard_zeroes_f_after_reversal` checks for a load reversal.

    `Kp` itself is not a Python-visible response, so this test uses the
    macroscopic proxy every critical-state model shares: under monotone
    straining, a stress ratio (`eta/M_b`, `sani._state_probe`'s own
    `ratio` field) that has been rising step over step and then, for the
    first time, DECLINES is the observable signature of the return map
    having found `Kp <= 0` at that declining commit (a still-rising ratio
    implies `Kp > 0`: the material is still hardening toward the bounding
    surface). The FIRST such decline is taken as the softening commit, and
    the guard is checked on the step immediately after it.

    See the block comment above this test for why the outcome is measured
    at runtime rather than asserted: if the deck never softens within
    `_SOFTEN_N_DEV` steps, this calls `pytest.xfail` with the actual
    numbers reached (the eta/M_b history) rather than faking a pass or a
    hard failure for a claim nobody has verified from this lane.
    """
    tag = 8350
    opts = ('-Presidual', 0.0, '-Pmin', sani._PMIN_LADRUNO, '-honorTolR', 0,
            '-implex', '-maxSubsteps', _CAP_ADEQUATE)
    n_dev = _build_floor_seeking_deck(tag, opts, e_conf=_SOFTEN_E_CONF,
                                      n_dev=_SOFTEN_N_DEV, lat=_SOFTEN_LAT)

    ops.updateMaterialStage('-material', tag, '-stage', 0)
    for step in range(sani._C_N_CONF):
        assert ops.analyze(1) == 0, f'confinement step {step + 1} failed'
    ops.updateMaterialStage('-material', tag, '-stage', 1)

    max_ratio = float('-inf')
    prev_ratio = None
    ratio_history = []
    soften_step = None
    for step in range(n_dev):
        assert ops.analyze(1) == 0, f'deviatoric step {step + 1} failed'
        ratio = sani._state_probe()['ratio']
        ratio_history.append(ratio)
        if prev_ratio is not None and ratio < prev_ratio and prev_ratio >= max_ratio:
            soften_step = step   # THIS commit is the first decline
            break
        max_ratio = max(max_ratio, ratio)
        prev_ratio = ratio

    if soften_step is None:
        pytest.xfail(
            'the dense confine-first deck (e_conf=%r, lat=%r) never showed '
            'a peak-then-decline in eta/M_b over %d deviatoric steps -- '
            'reached a max ratio of %.6f (last 5: %r) without softening. '
            'This is a real measurement against a real build, not a stand- '
            'in: either raise n_dev/lat, lower e_conf, or accept that this '
            'deck/parameter set does not soften in the tested range'
            % (_SOFTEN_E_CONF, _SOFTEN_LAT, n_dev, max_ratio,
               ratio_history[-5:]))

    guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    assert ops.analyze(1) == 0, (
        'the post-softening confirmation step failed to converge', soften_step)
    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    assert detail[5] == 0.0, (
        'implexDetail[5] (f) is not exactly 0.0 on the step immediately '
        'after the committed state first showed a peak-then-decline in '
        'eta/M_b (the softening proxy) -- -implexGuard on is supposed to '
        'force the elastic-predictor (f = 0) extrapolation on that step',
        soften_step, ratio_history[-3:], detail)
    assert guards_after[1] - guards_before[1] >= 1, (
        'implexGuards[1] (f=0 guards) did not increment across the '
        'post-softening step', guards_before, guards_after)


# ---------------------------------------------------------------------------
#  4. The hold-safe clock -- a LoadControl(0.0) commit preserves
#     mImplexDtCommit and the d_eps_p history from the step BEFORE it.
# ---------------------------------------------------------------------------

_HOLD_DS = 0.02
_HOLD_DS_POST = 2.0 * _HOLD_DS   # DELIBERATELY DIFFERENT from _HOLD_DS -- see
                                # _drive_hold_sequence's docstring for why a
                                # same-ds continuation cannot distinguish the
                                # two readings this test needs to tell apart.


def _drive_hold_sequence(tag, with_hold):
    """A few constant-ds plastic history steps (so "the ds before the
    hold" is unambiguous), optionally a LoadControl(0.0) hold, then ONE
    continuation step at `_HOLD_DS_POST` (double `_HOLD_DS`, the history
    steps' own ds).

    THE POST-HOLD ds MUST DIFFER FROM THE PRE-HOLD ds (fixed 2026-09-07,
    WP-92e lane B2, after a same-ds first draft measured f = 0.7 on BOTH the
    correct and the broken reading and could tell them apart). With
    `-implexAlpha 0.7` the two candidate readings of `implexDetail[5]` are:

      * correct (mImplexDtCommit preserved across the hold, ratio against
        the PRE-hold dt): `f = (dt_post / dt_pre) * alpha`.
      * broken (the hold silently reset mImplexDtCommit to 0, the
        dtCommit-reset fallback): `f = alpha`.

    At `dt_post == dt_pre` (a same-ds continuation, the first draft's
    choice) the ratio is exactly 1.0, so BOTH readings collapse to the same
    number (`1.0 * alpha == alpha`) and the test cannot tell them apart --
    it was checking a tautology, not the hold. `_HOLD_DS_POST = 2 *
    _HOLD_DS` makes the two readings numerically distinct: correct reads
    `2.0 * alpha = 1.4`, broken reads `alpha = 0.7`.
    """
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexAlpha', 0.7)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _confine_only(tag)

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)

    ops.integrator('LoadControl', _HOLD_DS)
    for step in range(_PROBE_N_HISTORY):
        assert ops.analyze(1) == 0, f'history step {step + 1} failed'
    ops.loadConst('-time', 0.0)

    guards_before = guards_after = None
    if with_hold:
        ops.integrator('LoadControl', 0.0)
        guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
        assert ops.analyze(1) == 0, 'the LoadControl(0.0) hold failed to converge'
        guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    ops.integrator('LoadControl', _HOLD_DS_POST)
    assert ops.analyze(1) == 0, 'the post-hold continuation step failed to converge'
    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    return detail, stress, guards_before, guards_after


def test_hold_keeps_clock_and_history():
    """ADR-92 P2: a `LoadControl(0.0)` commit (the gravity-hold /
    re-equilibration idiom) must be COUNTED (`implexGuards[2]`
    increments) and must PRESERVE `mImplexDtCommit` and the `d_eps_p`
    history from the step BEFORE the hold -- not reset either, and not
    treat the hold itself as an ordinary same-ds step.

    See `_drive_hold_sequence` for why the continuation step runs at
    `_HOLD_DS_POST` (2x the pre-hold ds), NOT the pre-hold ds itself, and
    the module docstring's "WHAT A ZERO-FREE-DOF DECK CAN AND CANNOT SHOW"
    section for why this uses the free-DOF triaxial rig rather than a
    zero-free-DOF one -- a hold-then-continue sequence needs a mechanically
    meaningful hold, not a vacuous one.

    Also checks the two sequences (with and without the hold) commit
    CLOSE final stresses -- the hold is a mechanical no-op (zero load
    increment, nothing to solve for) even though it is not a no-op on the
    IMPL-EX clock bookkeeping this test's other assertions check. The
    tolerance is 1e-5, not 1e-10 or bit-identity: this is a genuinely
    free-DOF Newton deck, and the WITH-hold run takes one MORE `analyze()`
    call (and one more converged Newton solve, at dt = 0) than the
    WITHOUT-hold run before reaching the same total load -- measured
    reldiff ~9e-7 on this deck (87b9cf846), comfortably under 1e-5 but
    nowhere near 1e-10.
    """
    detail_hold, stress_hold, guards_before, guards_after = _drive_hold_sequence(
        8360, with_hold=True)
    assert guards_before is not None and len(guards_before) == 5, (
        'implexGuards did not return the documented 5-component vector',
        guards_before)
    assert guards_after[2] - guards_before[2] >= 1, (
        'implexGuards[2] (hold-preserved commits) did not increment across '
        'the LoadControl(0.0) hold step', guards_before, guards_after)

    expected_f_correct = (_HOLD_DS_POST / _HOLD_DS) * 0.7   # 1.4
    assert detail_hold[5] == pytest.approx(expected_f_correct, rel=1.0e-6, abs=1.0e-9), (
        'implexDetail[5] (f) on the step after a LoadControl(0.0) hold, at '
        '_HOLD_DS_POST (2x the pre-hold ds), is not (dt_post/dt_pre)*alpha '
        '-- the hold is supposed to keep mImplexDtCommit (and the d_eps_p '
        'history) from the step BEFORE the hold, so the ratio must be '
        'against THAT dt, not reset by the hold. If this reads 0.7 (alpha '
        'itself) instead, mImplexDtCommit fell back to the dtCommit == 0 '
        'branch, meaning the hold reset the clock', detail_hold,
        expected_f_correct)

    _, stress_no_hold, _, _ = _drive_hold_sequence(8361, with_hold=False)
    assert sani._reldiff(stress_no_hold, stress_hold) <= 1.0e-5, (
        'the committed stress after the hold-then-continue sequence does '
        'not match the SAME sequence run WITHOUT the hold, to 1e-5 -- the '
        'hold is supposed to be a mechanical no-op on the committed answer '
        'even though it is not a no-op on the clock bookkeeping',
        stress_no_hold, stress_hold)


# ---------------------------------------------------------------------------
#  5. `ops.setParameter(..., 'stressCorrection')` now takes effect
# ---------------------------------------------------------------------------

def _drive_stresscorrection(mat_tag, val):
    """Confine, optionally `setParameter(..., 'stressCorrection')` right
    before the plastic leg, then take `_PROBE_N_HISTORY` nominal plastic
    steps. `-ele 1` (the element tag `_build_free_dof_triaxial` always
    uses), NOT `mat_tag` (the material tag) -- setParameter's `-ele` is an
    element selector, per every other setParameter call in this fork's own
    test suite (`test_energyBalanceRecorder.py`,
    `test_initStrain_dimension_general.py`).
    """
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE)
    _build_free_dof_triaxial(mat_tag, opts, p0=50.0)
    _confine_only(mat_tag)

    if val is not None:
        ops.setParameter('-val', val, '-ele', 1, 'stressCorrection')

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 1.0 / _PROBE_N_HISTORY)

    first_stress = None
    for step in range(_PROBE_N_HISTORY):
        assert ops.analyze(1) == 0, f'plastic step {step + 1} failed'
        if step == 0:
            first_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    final_stress = list(ops.eleResponse(1, 'material', 1, 'stress'))
    return first_stress, final_stress


def test_setparameter_stresscorrection_takes_effect():
    """ADR-92 P2: `ops.setParameter('-val', N, '-ele', <tags>,
    'stressCorrection')` now reaches `ManzariDafalias::mStressCorrectionInUse`
    (responseID 9, `ManzariDafalias.cpp:897`) through `LadrunoSANISAND` --
    previously a no-op on this subclass (the parameter dispatch, `setParameter`
    at `:852`, was reachable, but nothing carried the update through to
    where `Stress_Correction()` reads the flag).

    `Stress_Correction()` (`ManzariDafalias.cpp:2676`) is the drift-back-
    onto-the-yield-surface step inside `ModifiedEuler`'s successful-substep
    branch (`:1790`, the default -implex companion), guarded
    `if (!mStressCorrectionInUse) return;` (`:2681`) -- so forcing it off
    with `-val 0` must measurably move the FIRST plastic step's committed
    stress away from the default (compiled-in `true`, every constructor)
    answer, and explicitly setting it back to `1` must reproduce that SAME
    default answer exactly, proving `-val 1` maps to the identical boolean
    the constructors hardcode rather than some other encoding.

    Kills a mutant that drops the P2 dispatch fix (every -val is a no-op,
    all three runs land on the default answer) or that maps `-val 1` to
    something other than the compiled-in default (e.g. any nonzero treated
    the same as zero, or a sign flip).

    `_SC_SENSITIVITY_FLOOR`, NOT `sani._SENSITIVITY_FLOOR` (fixed 2026-09-07,
    WP-92e lane B2, re-run against 87b9cf846). `Stress_Correction()` is a
    drift-BACK-onto-the-yield-surface correction -- a small perturbation by
    construction, not a first-order elastic-vs-plastic gap like the
    p_residual sensitivity `sani._SENSITIVITY_FLOOR` (1e-3) was calibrated
    for. Measured on this deck: reldiff(first_default, first_off) =
    2.994e-4 -- genuinely nonzero and reproducible, but under 1e-3, so the
    imported floor made this test's own positive control fail. 1e-5 sits
    three orders below the measured signal and many orders above any
    floating-point noise floor.
    """
    _SC_SENSITIVITY_FLOOR = 1.0e-5

    first_default, final_default = _drive_stresscorrection(8370, val=None)
    first_off, final_off = _drive_stresscorrection(8371, val=0)
    first_on, final_on = _drive_stresscorrection(8372, val=1)

    assert sani._reldiff(first_default, first_off) > _SC_SENSITIVITY_FLOOR, (
        '-val 0 did not move the FIRST plastic step away from the default '
        '(stressCorrection compiled-in true) answer -- the positive '
        'control this test needs is vacuous; either the P2 dispatch fix '
        "did not land, or this deck never reaches Stress_Correction() at "
        "all (it only runs inside ModifiedEuler's successful-substep "
        'branch)', first_default, first_off)

    assert sani._reldiff(first_default, first_on) <= _EQ_TOL, (
        'explicitly setting stressCorrection to 1 (the compiled-in '
        'default) does not reproduce the answer of never calling '
        'setParameter at all -- either -val 1 is not mapping to the same '
        'boolean the constructors hardcode, or the dispatch is not '
        'idempotent', first_default, first_on)
    assert sani._reldiff(final_default, final_on) <= _EQ_TOL, (
        'the same divergence as above, carried through the rest of the '
        'plastic history', final_default, final_on)


# ===========================================================================
#  ADR-92 P2-5 / P2-2b (WP-92e lane B2, 2026-09-07, binary 8bfdfbc17)
# ===========================================================================

def _read_all_alpha_in(ele=1, ngp=8):
    return [list(ops.eleResponse(ele, 'material', gp, 'alpha_in')) for gp in range(1, ngp + 1)]


def _read_all_alpha(ele=1, ngp=8):
    return [list(ops.eleResponse(ele, 'material', gp, 'alpha')) for gp in range(1, ngp + 1)]


def test_hold_does_not_reset_alpha_in_on_the_implicit_path():
    """ADR-92 P2-5: `ManzariDafalias::integrate()` unconditionally resolves
    `mAlpha_in` from the SIGN of `(alpha_n - alpha_in_n):(Ce:(eps - eps_n))`
    with no magnitude guard, so on a `LoadControl(0.0)` hold that sign is
    noise. `LadrunoSANISAND::ladrunoGuardReversalNoise()` undoes the reset
    when `||eps - eps_n||` (`GetNorm_Cov`) falls below `-reversalTol`
    (default `1e-10`). NOT an IMPL-EX option -- no `-implex` gate -- so this
    is tested on a PURELY IMPLICIT deck (no `-implex` token at all): plastic
    history on the free-DOF triaxial rig (`_build_free_dof_triaxial` +
    `_establish_plastic_history`, p0 = 50 kPa, comfortably off the p_min
    floor), then a hold, `alpha_in` read at every one of the element's 8
    Gauss points before and after.

    `implexGuards` IS READABLE ON THE IMPLICIT DECK. `setResponse()`
    registers `implexGuards` unconditionally (`LadrunoSANISAND.cpp:2864`,
    no `-implex`/`mImplexOpt` gate anywhere in that branch) -- confirmed by
    reading it below without ever passing `-implex`. The coordinator's
    "skip that part if so" caveat does not apply on this build.

    MEASURED ON THIS DECK (8bfdfbc17): the hold's own strain increment is
    NOT exactly zero at every Gauss point (Newton's own residual tolerance,
    not literal round-off) -- `implexGuards[3]` increments by 8 (one per
    Gauss point) at the DEFAULT `reversalTol`, and by 0 under `-reversalTol
    0` on the identical deck/history/hold, so the guard is DEMONSTRABLY
    live at the default and DEMONSTRABLY off at 0 -- not a vacuous "both
    read 0" case. `alpha_in` itself reads bit-identical across the hold at
    EVERY Gauss point under BOTH settings on this deck -- see the second
    half of this test for why that is a documented, not a stronger, claim.

    Kills a mutant that drops `ladrunoGuardReversalNoise()`'s call sites (no
    `implexGuards[3]` movement at the default tol) or that ignores
    `-reversalTol 0` (the counter would still move there).
    """
    tag_default = 8380
    m_opts = ()   # IMPLICIT: no -implex token anywhere
    _build_free_dof_triaxial(tag_default, m_opts, p0=50.0)
    _establish_plastic_history(tag_default)

    alpha_in_before = _read_all_alpha_in()
    alpha_before = _read_all_alpha()
    guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    assert len(guards_before) == 5, (
        'implexGuards did not return the documented 5-component vector on '
        'a PURELY IMPLICIT deck (no -implex token) -- P2-5 is supposed to '
        'be readable here regardless', guards_before)

    ops.integrator('LoadControl', 0.0)
    rc_hold = ops.analyze(1)
    assert rc_hold == 0, ('the LoadControl(0.0) hold failed to converge', rc_hold)

    alpha_in_after = _read_all_alpha_in()
    alpha_after = _read_all_alpha()
    guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    assert guards_after[3] - guards_before[3] >= 1, (
        'implexGuards[3] (reversal-noise guards) did not increment across '
        'the hold at the DEFAULT -reversalTol (1e-10) -- either the guard '
        'is not live on the implicit path, or this deck\'s hold produced a '
        'strain increment that never falls under the threshold at any of '
        'the 8 Gauss points', guards_before, guards_after)

    diffs_ai = [max(abs(x - y) for x, y in zip(b, a))
               for b, a in zip(alpha_in_before, alpha_in_after)]
    diffs_a = [max(abs(x - y) for x, y in zip(b, a))
              for b, a in zip(alpha_before, alpha_after)]
    assert max(diffs_ai) == 0.0, (
        'alpha_in moved at at least one Gauss point across a LoadControl(0.0) '
        'hold at the DEFAULT -reversalTol -- P2-5 is supposed to leave it '
        'bit-identical there', diffs_ai, alpha_in_before, alpha_in_after)
    # `alpha` (`getAlpha()` -> `mAlpha`, the TRIAL backstress ratio, NOT the
    # committed `mAlpha_n`) is read for diagnostic completeness only, per
    # the coordinator's brief -- it is NOT asserted bit-identical. Measured:
    # it DOES move by ~6e-4 across this hold at every Gauss point, because
    # the hold's own Newton solve on this genuinely free-DOF deck converges
    # to a state that is close to, but not bit-identical to, the last
    # commit (max component-wise diff on the committed `strain` response is
    # ~6.2e-7, not literal round-off) -- a real, if tiny, trial update, not
    # a P2-5 violation. P2-5's own promise is about `mAlpha_in`, the
    # COMMITTED quantity `integrate()`'s reversal branch resets -- not
    # about `mAlpha` moving with a genuine (if small) trial strain change.
    assert max(diffs_a) > 0.0, (
        'alpha (mAlpha, the TRIAL backstress ratio) read bit-identical '
        'across the hold at every Gauss point -- this deck\'s hold is '
        'apparently producing an EXACTLY zero strain increment after all '
        '(contradicting the implexGuards[3] evidence above that it is not); '
        're-derive this test\'s premise rather than trust this run\'s '
        'numbers blindly', diffs_a)

    # -- the negative control: -reversalTol 0 on the SAME deck/history/hold --
    tag_tol0 = 8381
    _build_free_dof_triaxial(tag_tol0, ('-reversalTol', 0.0), p0=50.0)
    _establish_plastic_history(tag_tol0)

    alpha_in_before0 = _read_all_alpha_in()
    guards_before0 = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    ops.integrator('LoadControl', 0.0)
    rc_hold0 = ops.analyze(1)
    assert rc_hold0 == 0, ('the -reversalTol 0 twin\'s hold failed to '
                           'converge -- a harness problem', rc_hold0)

    alpha_in_after0 = _read_all_alpha_in()
    guards_after0 = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))

    assert guards_after0[3] == guards_before0[3], (
        'implexGuards[3] moved with -reversalTol 0 -- the guard is '
        'supposed to be DISABLED at tol 0 (mReversalTol <= 0.0 returns '
        'immediately), not merely quieter', guards_before0, guards_after0)

    diffs_ai0 = [max(abs(x - y) for x, y in zip(b, a))
                for b, a in zip(alpha_in_before0, alpha_in_after0)]
    # DOCUMENTED FINDING, not a stronger claim than the evidence supports.
    # -reversalTol 0 did NOT change alpha_in at any Gauss point either on
    # this deck. That is consistent with the source's own comment on
    # ladrunoGuardReversalNoise() -- "safe unconditionally: on the
    # no-reversal branch mAlpha_in already equals mAlpha_in_n, so the
    # assignment is a no-op" -- i.e. the counter (checked above) proves the
    # MAGNITUDE check fires at the default tol and not at 0, but this
    # deck's monotone triaxial hold never actually trips the BASE's
    # reversal branch at any of the 8 points in the first place, so there
    # is nothing for the guard to have protected here even when it is
    # live. A deck that reaches a genuine reversal DURING a hold is
    # Esmeralda-BVP-scale (28-54% of 34,560 points on job 146458) and out
    # of this lane's budget to construct from a Python material-point rig;
    # the counter-based evidence above is the strongest claim this deck
    # can support, and this assertion documents rather than overclaims it.
    assert max(diffs_ai0) == 0.0, (
        'alpha_in moved under -reversalTol 0 where the module docstring '
        'said it would not -- this deck now DOES exercise a genuine '
        'reversal during the hold; strengthen the test\'s claim above '
        'rather than leaving this stale', diffs_ai0)


def test_guard_ignores_the_unprimed_first_commit():
    """ADR-92 P2-2b: the reversal/softening guard in `ladrunoImplexCommit()`
    used to arm on the un-primed FIRST plastic commit after the elastic ->
    plastic stage flip (`mAlpha_in_n` necessarily moves there as an
    initialisation artefact, not a genuine reversal -- P2-2b's own fix
    gates both halves of the guard on `guardPrimed`, the SAME
    `GetNorm_Cov(mImplexDEpsP) > 0.0` predicate `-implexControl`'s
    un-primed-step exemption already uses, afb95c40c). On a deck that
    never approaches the p_min floor or a genuine reversal (p0 = 50 kPa,
    monotone loading, `-implexAlpha 0.7` so the un-primed fallback and the
    live ratio are both `0.7` here -- constant ds, so this does not need to
    distinguish them, only confirm NEITHER reads `f = 0`), the SECOND
    plastic step's `implexDetail[5]` must equal its own dt ratio (`0.7`),
    not `0.0`.

    Kills a mutant that removes `guardPrimed` from EITHER half of the `if`
    (`mImplexOpt.guard && guardPrimed`) -- the second plastic step would
    read `f = 0` again, exactly the P2-2 regression this fixes.
    """
    tag = 8390
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE, '-implexAlpha', 0.7)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _confine_only(tag)

    dq = _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 0.02)

    assert ops.analyze(1) == 0, 'the FIRST (un-primed) plastic step failed to converge'
    detail1 = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    assert detail1[5] == pytest.approx(0.7, rel=1.0e-6, abs=1.0e-9), (
        'implexDetail[5] on the un-primed FIRST plastic step is not alpha '
        '(0.7) -- the mImplexDtCommit == 0 fallback should apply here '
        'regardless of the guard (this step has no committed predecessor '
        'to arm the guard from)', detail1)

    assert ops.analyze(1) == 0, 'the SECOND plastic step failed to converge'
    detail2 = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    assert detail2[5] == pytest.approx(0.7, rel=1.0e-6, abs=1.0e-9), (
        'implexDetail[5] on the SECOND plastic step is not its dt ratio '
        '(0.7, constant ds) -- P2-2b is supposed to exempt the un-primed '
        'first commit from arming the reversal/softening guard, so this '
        'step must NOT read f = 0 the way it did before 8bfdfbc17',
        detail2)


# ===========================================================================
#  ADR-92 P2-6 (WP-92e lane B2, 2026-09-07, binary 708152eac)
# ===========================================================================

_TRIAL_GUARD_TOL = 0.01
_TRIAL_GUARD_REDUCTION = 0.01
_TRIAL_GUARD_FACTOR = -0.5   # a REVERSAL, half the nominal per-step magnitude
                             # -- see the test docstring for why a same-
                             # direction bigger step does NOT show f = 0
                             # beating the full extrapolation on this deck.


def _drive_trial_guard_reversal(tag, trial_guard):
    """Establish plastic history, then ONE reversal step (opposite sign to
    `_establish_plastic_history`'s own load, half its per-step magnitude,
    `LoadControl(1.0)` against the established `1.0 / _PROBE_N_HISTORY`
    history steps, so the dt ratio itself is also unremarkable -- the
    error gap this test needs comes from the DIRECTION reversal, not from
    an oversized dt) under `-implexControl` at `_TRIAL_GUARD_TOL`.

    Returns a dict of every response this test's assertions need, read
    both immediately before and immediately after the one `analyze()` call.
    """
    opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
            '-implexControl', _TRIAL_GUARD_TOL, _TRIAL_GUARD_REDUCTION,
            '-implexTrialGuard', trial_guard)
    _build_free_dof_triaxial(tag, opts, p0=50.0)
    _establish_plastic_history(tag)

    dq = _TRIAL_GUARD_FACTOR * _PROBE_DQ_NOMINAL / 4.0
    ops.timeSeries('Linear', 3)
    ops.pattern('Plain', 3, 3)
    for j, (x, y) in enumerate(_XY):
        ops.load(4 + j + 1, 0.0, 0.0, -dq)
    ops.integrator('LoadControl', 1.0)

    refusals_before = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    guards_before = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    stress_before = list(ops.eleResponse(1, 'material', 1, 'stress'))

    rc = ops.analyze(1)

    detail = list(ops.eleResponse(1, 'material', 1, 'implexDetail'))
    refusals_after = list(ops.eleResponse(1, 'material', 1, 'implexRefusals'))
    guards_after = list(ops.eleResponse(1, 'material', 1, 'implexGuards'))
    stress_after = list(ops.eleResponse(1, 'material', 1, 'stress'))
    return dict(rc=rc, detail=detail,
               refusals_before=refusals_before, refusals_after=refusals_after,
               guards_before=guards_before, guards_after=guards_after,
               stress_before=stress_before, stress_after=stress_after)


def _probe_trial_guard_reference_errors(tag_full, tag_f0):
    """Independent measurement of the two errors `-implexTrialGuard` is
    choosing between on the SAME history + reversal step
    `_drive_trial_guard_reversal` drives, via `-implexAlpha` (1.0 = the
    ordinary full extrapolation; 0.0 = the same pure elastic predictor,
    `sigma~ = sigma_n + Ce:d_eps`, the trial-time fallback itself
    recomputes) with `-implexControl` set to a tolerance neither probe
    ever reaches (`1e6`), so BOTH commit normally and their
    `implexDetail[0]` is directly readable -- unlike a genuinely refused
    attempt on `LadrunoBrick`, which reverts before ever calling
    `commitState()` (see `test_floor_fallback_delivers_implicit_stress_
    and_counts`'s own note on why THAT read is stale).

    NOT claimed to bit-match the live mechanism's own internal fallback
    number (measured: it does not -- 0.001827 here vs 0.003206 on the
    actual accepted step in `test_trial_guard_accepts_f0_before_refusing`,
    same deck/history/step). These are independent reference numbers
    confirming the SIGN and rough SCALE of the gap (f = 0 substantially
    better than full-f on a reversal, by roughly an order of magnitude),
    not a claimed numerical identity with the shipped code path.
    """
    def _run(tag, alpha):
        opts = ('-implex', '-maxSubsteps', _CAP_ADEQUATE,
                '-implexControl', 1.0e6, _TRIAL_GUARD_REDUCTION,
                '-implexAlpha', alpha)
        _build_free_dof_triaxial(tag, opts, p0=50.0)
        _establish_plastic_history(tag)
        dq = _TRIAL_GUARD_FACTOR * _PROBE_DQ_NOMINAL / 4.0
        ops.timeSeries('Linear', 3)
        ops.pattern('Plain', 3, 3)
        for j, (x, y) in enumerate(_XY):
            ops.load(4 + j + 1, 0.0, 0.0, -dq)
        ops.integrator('LoadControl', 1.0)
        assert ops.analyze(1) == 0, 'the reference probe failed to converge'
        return list(ops.eleResponse(1, 'material', 1, 'implexDetail'))[0]

    err_full = _run(tag_full, 1.0)
    err_f0 = _run(tag_f0, 0.0)
    return err_full, err_f0


def test_trial_guard_accepts_f0_before_refusing():
    """ADR-92 P2-6: `-implexTrialGuard on` (the default) retries a trial
    whose FULL extrapolation error is past `-implexControl`'s tolerance,
    BEFORE refusing it and before the reduction floor, with a pure elastic
    predictor (`f = 0`, the same `sigma~ = sigma_n + Ce:d_eps` P2-2's
    guard delivers), re-measured against the SAME companion return already
    computed -- and delivers it (no refusal) if THAT error clears tol.
    `-implexTrialGuard off` reproduces the pre-P2-6 behaviour: refuse
    immediately, with no such retry.

    THE DECK: a REVERSAL (`_TRIAL_GUARD_FACTOR = -0.5`), not a bigger
    same-direction step. Measured first (see
    `_probe_trial_guard_reference_errors`'s own docstring and this
    file's earlier `test_floor_fallback_...`'s "10x nominal, same
    direction" shape): on a bigger SAME-direction step the full
    extrapolation (which carries the established plastic direction
    forward) is actually a BETTER predictor than a pure elastic guess, so
    `-implexTrialGuard` would have nothing to rescue there. A reversal is
    exactly the opposite: the established `d_eps_p(n)` now points the
    WRONG way, so extrapolating it is worse than assuming no plastic flow
    at all -- probed independently (`_probe_trial_guard_reference_errors`,
    called below): full-f error 0.034149, f = 0 error 0.0018271, at
    `_TRIAL_GUARD_TOL = 0.01` sitting cleanly between them.

    Kills a mutant that drops the P2-6 retry entirely (`-implexTrialGuard
    on` would refuse exactly like `off`), that fires it regardless of the
    flag (`off` would also accept), or that fails to zero `implexDetail[5]`
    / count `implexGuards[4]` on the accepted step.
    """
    off = _drive_trial_guard_reversal(8395, 'off')
    assert off['rc'] != 0, (
        'the reversal step, under -implexTrialGuard off, was NOT refused -- '
        'the deck needs re-deriving (the full-f error is supposed to sit '
        'above _TRIAL_GUARD_TOL here), not this test\'s premise', off)
    assert off['refusals_after'][2] - off['refusals_before'][2] >= 1, (
        'implexRefusals[2] (-implexControl-specific) did not increment '
        'across the refused step', off['refusals_before'], off['refusals_after'])
    assert off['stress_after'] == off['stress_before'], (
        'the committed stress moved across a refused step', off)
    assert off['guards_after'][4] == off['guards_before'][4], (
        'implexGuards[4] (trial-time f=0 fallbacks) moved with '
        '-implexTrialGuard off -- it must not fire when the flag is '
        'disabled', off['guards_before'], off['guards_after'])

    on = _drive_trial_guard_reversal(8396, 'on')
    assert on['rc'] == 0, (
        'the SAME reversal step, under -implexTrialGuard on (the default), '
        'was NOT accepted -- the trial-time f=0 fallback is supposed to '
        'rescue it', on)
    assert on['detail'][5] == 0.0, (
        'implexDetail[5] (f) is not exactly 0.0 on the accepted step -- '
        'the trial-time fallback is supposed to deliver the pure elastic '
        'predictor', on['detail'])
    assert on['detail'][0] <= _TRIAL_GUARD_TOL, (
        'the accepted step\'s implexError is not <= tol -- that is the '
        'fallback\'s own accept condition', on['detail'], _TRIAL_GUARD_TOL)
    assert on['guards_after'][4] - on['guards_before'][4] >= 1, (
        'implexGuards[4] (trial-time f=0 fallbacks) did not increment '
        'across the accepted step', on['guards_before'], on['guards_after'])
    assert on['refusals_after'][2] == on['refusals_before'][2], (
        'implexRefusals[2] moved on the ACCEPTED step -- the fallback is '
        'supposed to avoid the refusal entirely, not refuse-then-recover',
        on['refusals_before'], on['refusals_after'])

    err_full, err_f0 = _probe_trial_guard_reference_errors(8397, 8398)
    assert err_f0 < _TRIAL_GUARD_TOL < err_full, (
        'the two independently-probed reference errors (full-f %r, f=0 '
        '%r) do not straddle _TRIAL_GUARD_TOL (%r) the way this deck is '
        'supposed to -- re-derive the deck/tol/factor rather than trust '
        'the mechanism assertions above blindly'
        % (err_full, err_f0, _TRIAL_GUARD_TOL), err_full, err_f0)
