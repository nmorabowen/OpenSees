"""ADR 90 / WP-A2 — does the tau = 0 (UNREGULARIZED) collapse load of a strip
footing on SOFTENING SANISAND already converge under mesh refinement?

WHY THIS EXISTS (the strategy critic's un-descope gate for ADR 90)
------------------------------------------------------------------
ADR 90 proposes a Duvaut-Lions viscoplastic wrapper whose deliverable is a
mesh-convergent localization WIDTH, and (C5) a collapse load that converges in
`h` at a fixed, declared Deborah number.  Its named consumer is the TIMs / APE
ultimate-surface campaign, which fits loci from radial pushovers on
`LadrunoSANISAND`.  The whole warrant rests on a premise nobody has measured:
that the UNREGULARIZED collapse load is NOT already mesh-converged inside the
campaign's own tolerance.

The fork has that measurement for `DruckerPrager`, and there the answer is
"converges, from above": `tests/test_r3_prandtl_collapse_gate.py` reads
1.0842 / 0.9938 / 0.9513 of the exact Prandtl-Reissner answer at
h0 = 1.0 / 0.5 / 0.25, every leg a genuine plateau.  (Those are that gate's
`_MEASURED` band CENTRES -- the TIMs/APE reference sequence it asserts against.
Its own docstring table logs 1.0849 / 0.9977 / 0.9514 for the fork's own run of
it; the 0.06 / 0.39 / 0.01 % difference is the independent-driver spread its
+/- 3 % band exists to absorb.)  Perfect plasticity is
elliptic; there is nothing to regularize.  **SANISAND is not perfectly plastic.**
It softens from `M^b` back to `M_c` as the fabric evolves, and a NON-ASSOCIATED
softening continuum is where the rate-independent boundary-value problem loses
ellipticity and the answer starts depending on the mesh.

So this driver asks the question directly, at tau = 0, on the campaign's own
element:

    if q_u converges under refinement -> the wrapper is width research, not a
        P1 deliverable, and ADR 90's C5 has no warrant;
    if it does not converge (no plateau/peak, or a three-mesh band wider than
        the campaign's tolerance, or legs that seize) -> the wrapper has one.

**This driver contains no regularization.**  There is no tau, no wrapper, no
`LadrunoDuvautLions`; nothing in `SRC/` is needed and nothing was built for it.
It measures the STATUS QUO, which is ADR 90 §3.1 option (d), "the negative
control of the whole ADR".

WHAT IS REUSED, AND WHAT IS NEW
-------------------------------
Reused **by import**, not by copy, from `tests/test_r3_prandtl_collapse_gate.py`:
the graded plane-strain strip mesh (`_strip_mesh`, `_graded`, `_n_graded`), the
domain (B = 2 m, 14.5 B of clearance, 10 B of depth, grading 1.35), and the
ADR-63 D16 adaptive-stepping guards *including the PINNED `SUBDIV_BUDGET = 80`*
whose calibration history that module documents (at D16's default of 24 the
h0 = 0.25 DruckerPrager leg read 0.8988 "still hardening"; the budget, not the
mesh, was the measurement).  Reused from `tests/test_ladruno_sanisand.py`:
Gorini's calibrated `_PARAMS` (ADR 86 §5, kPa).

New here: the weighted / K0-geostatic staging SANISAND needs (§"THE DECK"), the
softening-aware termination classifier (§"CAPACITY"), and the localization-width
instrumentation (§"WIDTH").

THE DECK
--------
Plane-strain strip footing, `LadrunoBrick -formulation bbar` (ADR 90 S3 freezes
the element; tetrahedra are prohibited on failure legs), one element thick,
u_y = 0 on every node.  Three resolutions h0 = 1.0 / 0.5 / 0.25, the same
sequence the R3 gate runs.  Two densities: Gorini's `e_init = 0.6944` and a
DENSER `e_init = 0.60`.  The second is the case that matters: at p ~ 20 kPa the
state parameter is psi = -0.12 (Gorini) vs -0.22 (dense), so `M^b/M_c` is 1.54
vs 2.14 and `M^d/M_c` is 0.49 vs 0.29 -- the dense leg dilates hard and softens
hard, which is the ill-posed case.

**Weighted soil, not weightless.**  The R3 gate runs `gamma = 0` under a
surcharge because Prandtl-Reissner is a weightless problem with an exact
answer.  SANISAND cannot be run that way: its elastic moduli go as sqrt(p) and
its whole state (`psi`, `M^b`, `M^d`) is pressure-dependent, so p -> 0 is a
pathology, not a simplification.  Self weight through
`eleLoad -type -selfWeight` with `-b 0 0 -gamma`, `gamma = rho*g = 19.62`
kN/m^3 (`rho = 2.0` is `_PARAMS[17]`; dry sand, so effective = total).

**Confine first, then flip, then push** -- the ADR-86 emitter guide's three deck
rules (`86_ladruno_sanisand_apegmsh_emitter_guide.md` §2), all three of which
this deck would otherwise trip:

  1. *Do not shear during the elastic stage.*  Gravity is ramped with the
     material at stage 0; the K0 state it produces has eta = q/p = 0.855
     against `M_c` = 1.3309 (64 %), so `Elastic2Plastic` finds nothing to
     repair and does not inflate the calibrated friction constant.  MEASURED,
     and asserted: zero `Outside Bounding` events, max eta/M_c reported per leg.
  2. *`updateMaterialStage` explicitly, per tag, at every boundary* -- because
     `mElastFlag` is a `static` member of `ManzariDafalias`.
  3. *A proportional strain ramp with p_r = 0 never yields.*  The gravity ramp
     IS proportional, but it runs at stage 0 where no yielding is expected; the
     footing push that follows is strongly non-proportional.  `-Presidual` is
     left at the fork default 0.0 and `-Pmin` at the fork default 1e-3*P_atm
     (ADR-86 defaults, as the WP brief requires); the shallowest Gauss point
     sits at p = 6.2 kPa (h0 = 1.0), three decades above that floor, and zero
     low-p CLAMP events are MEASURED.

Rigid ROUGH footing: the footing nodes carry `u_x = 0` in addition to the
prescribed `u_z`.  (The R3 gate's footing is smooth, because Prandtl-Reissner's
exact answer is the smooth one.  There is no exact answer here, so the deck uses
the rough footing the campaign's problem actually has.  This is a deliberate
divergence from R3 and is the reason no number here may be compared to R3's.)

Push exactly on the R3 pattern: `sp` inside a load pattern under one-argument
`LoadControl`, pseudo-time IS the imposed footing displacement, adaptive
subdivision with the pinned budget.

CAPACITY -- the three-clause rule, adapted for SOFTENING
--------------------------------------------------------
`00_canonical_testbed.md` §1c: a limit load is a capacity only if (1) the tail
is flat, (2) the termination mode is admissible, (3) the controller was still
advancing freely.  Clause 1 is written for a PERFECTLY PLASTIC plateau.  A
softening material does not plateau -- it PEAKS and then descends, and R3's
`DROP_TRUNCATE` rule (cut the curve where q first falls 2 % below its max)
would delete exactly the branch this study exists to see.  So:

  * `DROP_TRUNCATE` is **DISABLED** here, and this is the one R3 guard that is
    not inherited.  It is a bottomed-out-run detector on a plateau problem and a
    signal deleter on a softening one.
  * clause 1 is satisfied by EITHER a plateau (|tail dq/ds| < 2 % of initial)
    OR a **passed peak** (q has fallen `POSTPEAK_DROP` below q_max and stayed
    there for `POSTPEAK_STEPS` steps).  The classifier records which.
  * `PEAK` joins `TARGET` and `BUDGET` as an admissible termination mode: the
    run stopped because the mechanics finished, not because a guard bound.
    `FLOOR`, `WALL` and `TRUNCATED` remain seizure and are never a capacity.

WIDTH
-----
Because this is a localization question, every leg also dumps the end-of-run
deviatoric plastic strain field (per element, CSV) and reports a crude band
width by ADR 90 §7.3's threshold-free metric

    w2 = sqrt(12 * Var),  Var = sum_e p_e[(x_e - xbar)^2 + hx_e^2/12] / sum_e p_e

over the `eps_q^p` profile along a HORIZONTAL line at a FIXED PHYSICAL DEPTH
under the footing edge, restricted to x >= 0 (the mesh is symmetric; a two-sided
profile would measure the footing width, not the band width).  The `hx^2/12`
term is the within-element variance of a piecewise-constant profile: it makes a
one-element band read exactly `hx` and a k-element top hat read exactly `k*hx`,
so the number is comparable across meshes.  The count of yielding elements is
reported alongside, and so is their total VOLUME -- a count is not a
mesh-convergent quantity and must not be read as one.

CONTROLS (all cheap, all mandatory, all per leg)
------------------------------------------------
  * `ops.ladrunoBuild()` and the module path, asserted against a pinned hash;
  * resultant identity `sum R_z(base) == gamma * V` after gravity;
  * the **1-D geostatic patch**, which the resultant is structurally blind to
    (`00_canonical_testbed.md` §1d).  A laterally-rollered, base-fixed box of a
    material with constant Poisson ratio under self weight has the exact
    solution `sigma_zz = gamma*z`, `sigma_xx = sigma_yy = K0 sigma_zz`,
    `K0 = nu/(1-nu)` -- exact even though SANISAND's moduli vary with p, because
    `nu` is constant so K0 is depth-independent.  MEASURED under `bbar`: the
    stress is piecewise constant at the ELEMENT CENTROID value to 1.1e-13, which
    is the b-bar volumetric average making itself visible.  A load-distribution
    or body-force-convention error moves this O(1) while leaving the resultant
    at 1e-16.
  * the per-Gauss-point mobilisation `eta/M_c` field at the stage flip, plus a
    grep of the engine log for `Outside Bounding` and for the low-p `CLAMPING`
    warning -- the M_c-inflation guard;
  * mesh sensitivity of the GRAVITY state: the patch error above is itself that
    measurement (it is a per-GP comparison against a mesh-independent closed
    form), and the summary prints it across the three meshes.

BUDGET
------
`WALL_BUDGET_S` = 90 min per leg, 6 legs.  Coarse meshes run first so a partial
answer exists early.  A leg that ends on the wall clock is mode `WALL`, is
reported **INADMISSIBLE**, and its q_max is not a result.

USAGE
-----
    python3.12 Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py \
        --out <dir> [--legs h1.0_e0.6944,...] [--wall 5400]

Then `sanisand_tau0_summary.py <dir>` for the table.

TWO ENVIRONMENT VARIABLES, and you will meet them the first time you run this
on a box other than the one that produced the report:

    LADRUNO_DIST_BIN          Directory holding the `opensees.pyd` to test.
                              DEFAULT: this repo's own `dist/bin`.  The reported
                              campaign set it to a sibling worktree that held the
                              pinned build; that path was ephemeral and is
                              deliberately not baked in.

    LADRUNO_A2_EXPECT_BUILD   The 40-char git hash the engine must report.
                              DEFAULT: the hash the reported campaign ran on
                              (`548fe911...`).  Point it at YOUR build's hash to
                              re-run this deck elsewhere -- deliberately, with
                              the new hash named, so the numbers stay
                              attributable.  Set it to `any` to skip the check;
                              the run then says so loudly and its numbers must
                              not be quoted alongside the report's.

Both assertions name these variables in their failure text, so a cold
reproduction gets an instruction rather than a bare `AssertionError`.

PROVENANCE: every leg writes `<tag>.json` carrying the engine hash, this
driver's path, the deck constants, the controls and the measurements; the CSVs
carry the same hash in a header comment.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import os
import sys
import time

# --------------------------------------------------------------------------
# Engine binding.  ADR-90 WP-A2 runs against a PINNED build; a collapse load
# with no engine hash attached cannot be re-attributed later (ADR 90 S4).
# --------------------------------------------------------------------------
#
#   LADRUNO_DIST_BIN    where to find `opensees.pyd`.  Defaults to this repo's
#                       own `dist/bin`, which is the conventional location and
#                       the one a cold clone will have.  The campaign reported
#                       in `_adr90_tau0_qu_band.md` set it to a sibling
#                       worktree's `dist/bin` that held the pinned build; that
#                       path was ephemeral and is deliberately NOT baked in here.
#   LADRUNO_A2_EXPECT_BUILD
#                       the 40-char git hash the engine must report.  Defaults
#                       to the hash the reported campaign ran on.  Override it
#                       to re-run this deck on a DIFFERENT build -- deliberately,
#                       with the new hash named, so the numbers stay attributable
#                       (ADR 90 S4: provenance is output).  Set to `any` to skip
#                       the check entirely; the run then prints whatever hash it
#                       found and its numbers must not be compared to the ones in
#                       the report.
#
_HERE_BOOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BIN = os.path.abspath(
    os.path.join(_HERE_BOOT, "..", "..", "..", "dist", "bin"))
EXPECTED_BUILD = os.environ.get(
    "LADRUNO_A2_EXPECT_BUILD", "548fe911427e90a2edfead05cb3672a738d25b6d")

_BIN = os.environ.get("LADRUNO_DIST_BIN", DEFAULT_BIN)
if os.path.isdir(_BIN):
    os.add_dll_directory(_BIN)
    sys.path.insert(0, _BIN)
else:
    # Say it HERE.  Falling through to `import opensees` gives a bare
    # `ModuleNotFoundError: No module named 'opensees'`, which tells a cold
    # reproduction nothing about the variable it needed to set -- and that is
    # the same defect as an unexplained assertion, just moved earlier.
    raise SystemExit(
        f"\nADR-90 WP-A2 driver: no OpenSees build directory at\n"
        f"    {_BIN}\n"
        f"({'from LADRUNO_DIST_BIN' if 'LADRUNO_DIST_BIN' in os.environ else
            'the default, this repo/worktree' + chr(39) + 's own dist/bin'})\n\n"
        f"Set LADRUNO_DIST_BIN to the directory holding the `opensees.pyd` you\n"
        f"mean to test, e.g.\n"
        f"    LADRUNO_DIST_BIN=/path/to/worktree/dist/bin python3.12 {__file__}\n\n"
        f"See this file's USAGE docstring for that variable and for\n"
        f"LADRUNO_A2_EXPECT_BUILD, which pins the engine hash.\n")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import numpy as np                                              # noqa: E402
import opensees as ops                                          # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.append(os.path.join(_REPO, "tests"))

import test_r3_prandtl_collapse_gate as r3                      # noqa: E402
from test_ladruno_sanisand import _PARAMS as _SANISAND_PARAMS   # noqa: E402

# The localization-width metric and its probe depths live in the SUMMARY module
# (which imports no engine), so there is exactly one implementation and an old
# campaign's field CSVs can be re-reduced on any box.
sys.path.insert(0, _HERE)
from sanisand_tau0_summary import (                                 # noqa: E402
    Z_PROBES, widths_from_field)


def assert_engine(strict: bool = True) -> str:
    """Fail loud on the wrong binary.  Printed at the top of every driver.

    The failure messages name the two environment variables, because the whole
    point of a fail-loud pin is defeated if a cold reproduction hits an
    assertion it cannot interpret.
    """
    path = os.path.abspath(ops.__file__)
    build = ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "unknown"
    if strict:
        assert path.lower().startswith(os.path.abspath(_BIN).lower()), (
            f"opensees resolved to {path}, not the expected build dir {_BIN}.\n"
            f"    Set LADRUNO_DIST_BIN to the directory holding the "
            f"opensees.pyd you mean to test.")
        if EXPECTED_BUILD.lower() != "any":
            assert build == EXPECTED_BUILD, (
                f"ladrunoBuild() = {build}, expected {EXPECTED_BUILD}.\n"
                f"    This run has no provenance and its numbers must not be "
                f"quoted alongside the reported campaign's.\n"
                f"    To re-run this deck on a different build, set "
                f"LADRUNO_A2_EXPECT_BUILD to that build's hash (deliberately, "
                f"so the numbers stay attributable), or to 'any' to skip the "
                f"check and forgo comparability.")
        else:
            print(f"    *** LADRUNO_A2_EXPECT_BUILD=any: build hash NOT pinned "
                  f"(running {build}); these numbers are not comparable to the "
                  f"reported campaign's ***")
    return build


# --------------------------------------------------------------------------
# Deck constants.  Everything reused from R3 is referenced through `r3.` so
# that a change there cannot silently diverge from a copy here.
# --------------------------------------------------------------------------
GRAV = 9.81                       # m/s^2
RHO = _SANISAND_PARAMS[17]        # 2.0 Mg/m^3 (kPa-m-s-Mg units)
GAMMA = RHO * GRAV                # 19.62 kN/m^3, dry => effective = total
NU = _SANISAND_PARAMS[1]          # 0.3129
K0 = NU / (1.0 - NU)              # 0.4554
P_ATM = _SANISAND_PARAMS[8]       # 101 kPa
M_C = _SANISAND_PARAMS[3]         # 1.3309

# ADR-86 defaults, explicitly emitted (the guide's rule: never rely on a
# default, and `-Pmin`'s fork default is 10x vanilla's).
OPT_PRESIDUAL = 0.0
OPT_PMIN = 1.0e-3 * P_ATM

# THE POSITIONAL OPTIONALS -- and why this deck must not leave them defaulted.
#
#   <IntScheme TanType JacoType TolF TolR>
#
# `TanType` selects what `getTangent()` returns (`ManzariDafalias3D.cpp:134-142`):
# 0 = `mCe`, the ELASTIC tangent; 1 = `mCep`, the continuum elastoplastic; 2 =
# `mCep_Consistent`.  **The parser default is 0**
# (`LadrunoSANISAND.cpp:117`, matching `OPS_ManzariDafaliasMaterial`), while the
# null / parallel constructors default to 2 (`ManzariDafalias.cpp:365`, `:426`).
# So a deck that emits only the 18 positional parameters -- which is what every
# existing fork SANISAND deck does, `tests/test_ladruno_sanisand.py` included --
# runs `algorithm Newton` against an ELASTIC tangent, i.e. de-facto modified
# Newton, with linear convergence.
#
# On a zero-free-DOF material-point deck that costs nothing.  On this
# boundary-value problem it was the difference between a measurement and a
# non-measurement: MEASURED at TanType 0, the h0 = 0.25 leg spent 350 s to
# advance 11 steps to s/B = 1.6e-4, ~40-65 state-determination passes per step.
# The tangent changes the ITERATION PATH, not the converged answer -- the stress
# update is untouched -- so this is a deck defect being fixed, not a knob being
# turned to make a result appear.  The A/B is reported in
# `_adr90_tau0_qu_band.md`.
#
# `mCep_Consistent` is UNSYMMETRIC for a non-associated flow rule, which is why
# the solver is `Pardiso -matrixType 0`.
INT_SCHEME = 1                    # ModifiedEuler with substepping (the default)
TAN_TYPE = 2                      # consistent elastoplastic -- NOT the default 0
JACO_TYPE = 1
TOL_F, TOL_R = 1.0e-7, 1.0e-7     # the parser defaults, emitted explicitly

N_GRAV = 10                       # gravity ramp steps (stage 0)
GRAV_TOL = 1.0e-9                 # NormDispIncr, metres
GRAV_ITER = 40

# The push.  DS_* and GROW_AFTER and SUBDIV_BUDGET are R3's, imported.
DS_BASE, DS_MIN = r3.DS_BASE, r3.DS_MIN
GROW_AFTER = r3.GROW_AFTER
SUBDIV_BUDGET = r3.SUBDIV_BUDGET          # PINNED at 80 -- see r3's docstring

# DS_MAX is NOT inherited.  R3's 1e-3 m is sized for a weightless deck pushed to
# s/B = 0.15 = 0.30 m.  This deck's target is s/B = 0.25 = 0.50 m on a material
# whose per-Gauss-point cost is two orders above DruckerPrager's, so at R3's cap
# the leg cannot reach its own peak inside any sane wall budget (MEASURED in the
# pilot: 145 s bought s/B = 0.0009).  Pinned here as its own named constant with
# that history, per the slow-tier contract, and reported in every result.
DS_MAX = 5.0e-3

SFRAC = 0.25                      # push target s/B

# THE CONVERGENCE TEST -- and a correction to the WP brief.
#
# The brief says "`test NormDispIncr` per the R3 gate".  The R3 gate does not
# use `NormDispIncr`: it uses **`NormUnbalance` at `1e-5 x the applied load`**
# (`test_r3_prandtl_collapse_gate.py`, `tol = 1.0e-5 * want`).  Both forms were
# MEASURED here and the force form is the one this deck must use, for a reason
# that is specific to this material:
#
#   * `NormDispIncr` with an ABSOLUTE tolerance in metres is tightest where the
#     steps are smallest -- i.e. tightest exactly where they are cheap and
#     loosest where they are expensive -- which is backwards for an adaptive
#     controller.
#   * SANISAND's stress update is a SUBSTEPPED `ModifiedEuler` return with a
#     hardcoded substep tolerance `TolE = 1e-4`, so the discrete stress-strain
#     map is only piecewise smooth and Newton STALLS rather than converging
#     quadratically.  MEASURED on the h0 = 0.5 leg: over 47 failed attempts the
#     residual displacement norm stalls at a median of 6.6e-7 m (min 3.4e-8,
#     max 1.3e-4).  A 1e-8 m target is therefore UNREACHABLE, and the run that
#     nominally used it was in fact carried by the relaxed third ladder rung on
#     18 of its 26 steps -- 65 of every 125 state-determination passes spent
#     failing two rungs that could not succeed.  That is the note-71 failure
#     mode in its purest form: the study was measuring its own convergence
#     test.
#   * A FORCE tolerance scaled to a deck-intrinsic force -- the model's own
#     weight -- is the same number on all three meshes by construction, which
#     is what a mesh-convergence study needs.  A displacement-norm tolerance is
#     not: the norm is taken over the free DOFs, and there are 3.6x more of
#     them at h0 = 0.25 than at h0 = 1.0, so the same nominal tolerance is a
#     different physical requirement on each mesh.
#
# PUSH_TOL is therefore RELATIVE TO THE MODEL WEIGHT `gamma*V`, and the
# `--test` switch keeps the displacement form available for the A/B that
# justified this choice.
# MEASURED, h0 = 0.5, 400 s of wall each, same box, same DS_MAX (the A/B table
# in `_adr90_tau0_qu_band.md`):
#
#   test / tol                reach (s/B)   steps   subdivisions   relaxed
#   NormDispIncr  1e-8 m        0.00106       26         1          18/26
#   NormUnbalance 1e-6 gamma*V  0.00218       31         0          11/31
#   NormUnbalance 1e-5 gamma*V  0.00442       37         0           -
#
# and the ANSWER moves by a MEDIAN of 0.3-0.75 % between all three at matched
# settlement (max excursions of 2.5-8 % are interpolation artefacts of
# comparing two adaptive step sequences on a curve whose slope is still ~25 kPa
# per 0.001 of s/B -- an s/B mismatch of 4e-5 is worth 3 % there.  That artefact
# VANISHES at the peak, where dq/ds -> 0, so q_u is far less tolerance-sensitive
# than these pre-peak numbers look).
#
# 1e-5 is pinned because it is the only setting that lets the step reach
# DS_MAX at all, and reach is the binding constraint on whether this study can
# answer its own question.  One leg is re-run at 1e-6 as the tolerance control.
PUSH_TEST = "NormUnbalance"
PUSH_TOL = 1.0e-5                 # x gamma*V; see the A/B table above
PUSH_TOL_DISP = 1.0e-8            # the metres form, for the A/B only
                                  # The ladder's third rung runs at 10x this and
                                  # every step that needed it is FLAGGED in the
                                  # `relaxed` column of the curve CSV and
                                  # counted in `nrelax`, so a leg that was
                                  # carried by the relaxed rung cannot pass
                                  # unnoticed.
PLATEAU_FRAC = r3.PLATEAU_FRAC    # 2 % of the initial tangent
FREE_ADVANCE_FLOOR_FACTOR = r3.FREE_ADVANCE_FLOOR_FACTOR

# Softening clause-1 (see the module docstring).  A peak counts as "passed" only
# when the curve has come down POSTPEAK_DROP below q_max and STAYED there for
# POSTPEAK_STEPS steps -- one dip is a solver artefact, a sustained descent is a
# post-peak branch.  5 % / 15 steps: the smallest drop that cannot be a single
# relaxed step, over a stretch longer than the 6-step recovery latch.
POSTPEAK_DROP = 0.05
POSTPEAK_STEPS = 15
PEAK_MIN_SFRAC = 0.002            # ignore blips inside the first 0.2 % of s/B

WALL_BUDGET_S = float(os.environ.get("LADRUNO_A2_WALL_BUDGET", 5400.0))

_CAPACITY_MODES = ("TARGET", "PEAK", "BUDGET")
# SEIZURE modes.  `00_canonical_testbed` §1b: "a leg that ends on that stop is
# reported INADMISSIBLE, never as a result".  These are checked explicitly and
# the word is printed, rather than being left to the reader to infer from a
# `CAPACITY = NO` column -- a leg can also fail to be a capacity for an ordinary
# reason (it simply has not peaked yet), and the two must not read the same.
# `KILLED` is here because the summariser assigns it to a leg whose process died
# mid-run: also a machine stopping, not the soil.
_SEIZURE_MODES = ("FLOOR", "WALL", "TRUNCATED", "KILLED")

# Width instrumentation: `Z_PROBES` and `widths_from_field` are
# imported from `sanisand_tau0_summary` (see the import block above).
#
# MATCHED-SETTLEMENT CHECKPOINTS.  Legs on different meshes do NOT terminate at
# the same s/B -- the fine ones are slower and stop earlier -- so a width or a
# load compared only at each leg's own end point compares three different
# problems.  Every leg therefore snapshots the plastic-strain field and the
# reaction at the SAME list of s/B values, and the summary compares those.
# This is the deliverable that survives a leg being stopped by the clock.
CHECKPOINTS = (0.002, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.25)

RESOLUTIONS = (1.0, 0.5, 0.25)
DENSITIES = (("e0.6944", _SANISAND_PARAMS[2]),   # Gorini's calibrated e_init
             ("e0.60", 0.60))                     # denser: the ill-posed case


def leg_tag(h0: float, ename: str) -> str:
    return f"h{h0}_{ename}"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _pick_system():
    """The tangent is UNSYMMETRIC (SANISAND's flow rule is non-associated)."""
    try:
        ops.system("Pardiso", "-matrixType", 0)
        return "Pardiso(-matrixType 0)"
    except Exception:
        ops.system("UmfPack")
        return "UmfPack"


def _dev_and_p(st):
    """(p, q, eta) from a 6-vector in element sign convention (compression -)."""
    p = -(st[0] + st[1] + st[2]) / 3.0
    s = (st[0] + p, st[1] + p, st[2] + p, st[3], st[4], st[5])
    j2 = 0.5 * (s[0] ** 2 + s[1] ** 2 + s[2] ** 2) + s[3] ** 2 + s[4] ** 2 + s[5] ** 2
    q = math.sqrt(3.0 * j2)
    return p, q, (q / p if p > 1.0e-9 else float("nan"))


def _eps_q_p(ep):
    """Deviatoric plastic strain invariant eps_q^p = sqrt(2/3 e^p:e^p).

    `plasticstrains` returns engineering shears, so the off-diagonals are
    halved back to tensor components before the contraction."""
    tr = (ep[0] + ep[1] + ep[2]) / 3.0
    e = (ep[0] - tr, ep[1] - tr, ep[2] - tr, 0.5 * ep[3], 0.5 * ep[4], 0.5 * ep[5])
    ee = e[0] ** 2 + e[1] ** 2 + e[2] ** 2 + 2.0 * (e[3] ** 2 + e[4] ** 2 + e[5] ** 2)
    return math.sqrt(2.0 / 3.0 * ee)


def _plastic_field(n_hex):
    """Element-mean eps_q^p over the 8 Gauss points, for every element."""
    epsq = np.zeros(n_hex)
    for e in range(1, n_hex + 1):
        acc, ngp = 0.0, 0
        for gp in range(1, 9):
            ep = ops.eleResponse(e, "material", gp, "plasticstrains")
            if not ep or len(ep) < 6:
                continue
            acc += _eps_q_p(ep)
            ngp += 1
        epsq[e - 1] = acc / ngp if ngp else 0.0
    return epsq


def _write_leg_json(out_dir, tag, payload):
    """Write `a2_<tag>.json` atomically (tmp + replace).

    Called at every checkpoint with `partial=True` and once more at the end
    with the full record, so a killed process still leaves a leg record the
    summariser can assemble.  Atomic because a half-written JSON read by the
    summariser is worse than no JSON at all."""
    final = os.path.join(out_dir, f"a2_{tag}.json")
    tmp = final + ".tmp"
    with open(tmp, "w") as jf:
        json.dump(payload, jf, indent=2, sort_keys=True, default=float)
    os.replace(tmp, final)


def _dump_field(path, header, xc, zc, hx, hz, vol, epsq, eta_field):
    with open(path, "w", newline="") as ffh:
        ffh.write(f"# {header}\n")
        fw = csv.writer(ffh)
        fw.writerow(["ele", "xc", "zc", "hx", "hz", "vol", "eps_q_p", "eta_grav"])
        for e in range(len(epsq)):
            fw.writerow([e + 1, f"{xc[e]:.9g}", f"{zc[e]:.9g}", f"{hx[e]:.9g}",
                         f"{hz[e]:.9g}", f"{vol[e]:.9g}", f"{epsq[e]:.9g}",
                         f"{eta_field[e]:.9g}"])


# --------------------------------------------------------------------------
# one leg
# --------------------------------------------------------------------------
def run_leg(h0, ename, e_init, out_dir, wall_budget=None, sfrac=SFRAC,
            ds_max=DS_MAX, tol=PUSH_TOL, tan_type=TAN_TYPE,
            test_type=PUSH_TEST, verbose=True):
    wall_budget = WALL_BUDGET_S if wall_budget is None else wall_budget
    tag = leg_tag(h0, ename)
    log_path = os.path.join(out_dir, f"a2_{tag}_engine.log")
    ops.logFile(log_path, "-noEcho")

    t_start = time.time()
    nodes, hexes, trib, sets = r3._strip_mesh(h0)
    n_nodes, n_hex = len(nodes), len(hexes)

    # element geometry, once (rectangular tensor grid)
    p_nod = nodes[hexes]
    xc = p_nod[:, :, 0].mean(axis=1)
    zc = p_nod[:, :, 2].mean(axis=1)
    hx = p_nod[:, :, 0].max(axis=1) - p_nod[:, :, 0].min(axis=1)
    hz = p_nod[:, :, 2].max(axis=1) - p_nod[:, :, 2].min(axis=1)
    vol = hx * hz * r3.THICK

    # ---- model ----------------------------------------------------------
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))

    par = list(_SANISAND_PARAMS)
    par[2] = e_init
    # Positional block first, THEN the flags -- the ADR-86 parser rejects a
    # positional after a flag by design.
    ops.nDMaterial("LadrunoSANISAND", 1, *par,
                   INT_SCHEME, tan_type, JACO_TYPE, TOL_F, TOL_R,
                   "-Presidual", OPT_PRESIDUAL, "-Pmin", OPT_PMIN,
                   "-honorTolR", 0)
    for e, conn in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                    "-geom", "linear", "-b", 0.0, 0.0, -GAMMA,
                    "-formulation", "bbar")

    # ONE fix() per node -- a second fix() ADDS a duplicate SP_Constraint.
    # u_y = 0 everywhere is what makes this plane strain.  The footing is ROUGH,
    # so its nodes also carry u_x = 0 (folded into the same single call).
    bt = set(sets["bottom"].tolist())
    xf = set(sets["xface"].tolist())
    ft = set(sets["footing"].tolist())
    for n in range(n_nodes):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if (n in xf or n in ft) else 0, 1, 0)

    # ---- stage 0: K0 gravity -------------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *range(1, n_hex + 1), "-type", "-selfWeight",
                0.0, 0.0, 1.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    solver = _pick_system()
    ops.test("NormDispIncr", GRAV_TOL, GRAV_ITER, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / N_GRAV)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 0)     # explicit, per tag
    for k in range(N_GRAV):
        assert ops.analyze(1) == 0, f"{tag}: gravity step {k + 1} failed"
    t_grav = time.time() - t_start

    # ---- controls on the gravity state ---------------------------------
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want = GAMMA * float(vol.sum())
    resultant_err = abs(rz / want - 1.0)
    assert resultant_err < 1.0e-6, (
        f"{tag}: initial equilibrium identity violated, {rz} vs {want} kN -- the "
        "body-force convention is wrong and nothing downstream is valid")

    # the field check the identity CANNOT make: the 1-D geostatic patch.
    err_zz = err_xx = 0.0
    eta_max = 0.0
    p_min = float("inf")
    eta_field = np.zeros(n_hex)
    for e in range(1, n_hex + 1):
        z_e = float(zc[e - 1])
        szz_ex = GAMMA * z_e                    # b-bar: piecewise constant at zc
        sxx_ex = K0 * szz_ex
        for gp in range(1, 9):
            st = ops.eleResponse(e, "material", gp, "stress")
            if not st or len(st) < 6:
                continue
            err_zz = max(err_zz, abs(st[2] - szz_ex) / abs(szz_ex))
            err_xx = max(err_xx, abs(st[0] - sxx_ex) / abs(sxx_ex))
            pp, qq, eta = _dev_and_p(st)
            p_min = min(p_min, pp)
            if eta == eta:                       # not nan
                eta_max = max(eta_max, eta)
                eta_field[e - 1] = max(eta_field[e - 1], eta)
    patch_err = max(err_zz, err_xx)

    # ---- stage flip ------------------------------------------------------
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    assert eta_max < M_C, (
        f"{tag}: max eta = {eta_max:.5f} at the stage flip is at or above "
        f"M_c = {M_C} -- Elastic2Plastic will INFLATE the calibrated friction "
        "constant and the leg does not run Gorini's soil")

    # ---- the push --------------------------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    base = [int(n) + 1 for n in sets["bottom"]]
    area = r3.B_FOOT * r3.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = sfrac * r3.B_FOOT

    ops.reactions()
    r0_foot = sum(ops.nodeReaction(t, 3) for t in foot)
    r0_base = sum(ops.nodeReaction(t, 3) for t in base)

    # Pseudo-time IS the imposed footing settlement (R3's idiom): a Linear
    # series makes the pattern factor equal lambda, the sp value is 1.0, and
    # starting the clock at uz0 means the first increment carries no jump.
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    _pick_system()
    ops.analysis("Static")

    # NormUnbalance takes an absolute FORCE; scale the relative tolerance by
    # the model's own weight so the number is identical on all three meshes.
    tol_abs = tol * want if test_type == "NormUnbalance" else tol
    ladder = [("Newton", tol_abs, 25, 0),
              ("NewtonLineSearch", tol_abs, 40, 0),
              ("KrylovNewton", 10.0 * tol_abs, 60, 1)]

    csv_path = os.path.join(out_dir, f"a2_{tag}_curve.csv")
    fh = open(csv_path, "w", newline="")
    w = csv.writer(fh)
    fh.write(f"# ADR90 WP-A2 tau=0 band, build {EXPECTED_BUILD}, leg {tag}\n")
    w.writerow(["s_m", "s_over_B", "q_foot_kPa", "q_base_kPa", "ds_mm",
                "relaxed", "wall_s"])

    rows, ds, good = [], DS_BASE, 0
    nfail = nrelax = nsub = 0
    mode, verdict = "TARGET", "reached the target settlement"
    checkpoints, cp_left = [], [c for c in CHECKPOINTS if c <= sfrac + 1e-12]
    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1.0e-12:
            break
        if time.time() - t0 > wall_budget:
            mode = "WALL"
            verdict = f"WALL-CLOCK budget spent at s/B = {s_now / r3.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test(test_type, tl, it, 0)
            ops.algorithm(algo)
            if ops.analyze(1) == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good = 0
            nsub += 1
            ds *= 0.5
            if nsub > SUBDIV_BUDGET:
                mode = "BUDGET"
                verdict = (f"pinned subdivision budget of {SUBDIV_BUDGET} spent "
                           f"at s/B = {s_now / r3.B_FOOT:.5f}")
                break
            if ds < DS_MIN:
                mode = "FLOOR"
                verdict = (f"step collapsed to the DS_MIN floor at s/B = "
                           f"{s_now / r3.B_FOOT:.5f} (every ladder rung failed "
                           f"at ds = {ds * 1000:.6f} mm)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= GROW_AFTER and ds < ds_max:
            ds, good = min(2.0 * ds, ds_max), 0
        ops.reactions()
        qf = -(sum(ops.nodeReaction(t, 3) for t in foot) - r0_foot) / area
        qb = -(sum(ops.nodeReaction(t, 3) for t in base) - r0_base) / area
        s = uz0 - ops.getTime()
        rows.append((s, s / r3.B_FOOT, qf, qb, ds * 1000.0, relaxed,
                     time.time() - t0))
        w.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
        if verbose and len(rows) % 25 == 0:
            print(f"    [{tag}] {len(rows):5d} s/B={s / r3.B_FOOT:.5f} "
                  f"q={qf:9.3f} ds={ds * 1000:.4f}mm "
                  f"t={time.time() - t0:7.1f}s", flush=True)

        # MATCHED-SETTLEMENT CHECKPOINT.  Taken on the first step at or past the
        # checkpoint (never interpolated: the field is a state, not a number),
        # so the summary can compare three meshes at the SAME s/B even when they
        # terminate at different ones.
        while cp_left and s / r3.B_FOOT >= cp_left[0]:
            cp = cp_left.pop(0)
            fld = _plastic_field(n_hex)
            snap = dict(cp_target=cp, s_over_B=s / r3.B_FOOT, q_foot=qf,
                        q_base=qb, wall_s=time.time() - t0, step=len(rows))
            snap.update(widths_from_field(fld, xc, zc, hx, hz, vol))
            checkpoints.append(snap)
            # Write the leg record NOW, marked partial.  A leg's JSON is
            # otherwise only written when `run_leg` returns, so a process that
            # is killed mid-leg -- which HAPPENED to this campaign's first
            # attempt, four processes at once, none of them past its first
            # leg -- leaves hours of curve and field CSVs that the summariser
            # cannot assemble.  The cost is one small file write per
            # checkpoint; the benefit is that the study survives its own
            # infrastructure.
            _write_leg_json(out_dir, tag, dict(
                tag=tag, h0=h0, e_name=ename, e_init=e_init,
                build=EXPECTED_BUILD, driver=os.path.abspath(__file__),
                date=datetime.datetime.now().isoformat(timespec="seconds"),
                partial=True, solver=solver, nodes=n_nodes, hexes=n_hex,
                dof=3 * n_nodes, gamma=GAMMA, K0=K0, M_c=M_C,
                presidual=OPT_PRESIDUAL, pmin=OPT_PMIN, push_tol=tol,
                push_test=test_type, tan_type=tan_type, ds_max=ds_max,
                subdiv_budget=SUBDIV_BUDGET, sfrac_target=sfrac,
                wall_budget_s=wall_budget, resultant_err=resultant_err,
                patch_err=patch_err, eta_max_grav=eta_max,
                eta_over_Mc=eta_max / M_C, p_min_grav=p_min,
                mode="RUNNING", verdict="leg still running at this snapshot",
                q_u=float(max(r[2] for r in rows)),
                s_end_over_B=s / r3.B_FOOT, steps=len(rows), nsub=nsub,
                nfail=nfail, nrelax=nrelax, wall_s=time.time() - t0,
                csv=csv_path, log=log_path, checkpoints=list(checkpoints)))
            _dump_field(
                os.path.join(out_dir, f"a2_{tag}_field_sB{cp:g}.csv"),
                f"ADR90 WP-A2 plastic field at CHECKPOINT s/B={cp:g} "
                f"(actual {s / r3.B_FOOT:.6f}), build {EXPECTED_BUILD}, "
                f"leg {tag}, q_foot={qf:.4f} kPa",
                xc, zc, hx, hz, vol, fld, eta_field)
            if verbose:
                print(f"    [{tag}] CHECKPOINT s/B={cp:g}: q={qf:.3f} kPa, "
                      f"w2(z={Z_PROBES[0]})="
                      f"{snap[f'w2_z{Z_PROBES[0]}']:.4f} m, "
                      f"yield {snap['n_yield_ele']} ele / "
                      f"{snap['vol_yield']:.3f} m3, "
                      f"t={time.time() - t0:.0f}s", flush=True)

        # softening clause 1: has a peak been PASSED and stayed passed?
        a = np.array([r[2] for r in rows])
        i_pk = int(a.argmax())
        if (rows[i_pk][1] >= PEAK_MIN_SFRAC
                and len(a) - i_pk > POSTPEAK_STEPS
                and float(a[i_pk:].max()) == float(a[i_pk])
                and bool((a[i_pk + 1:] < (1.0 - POSTPEAK_DROP) * a[i_pk]).all())
                and len(a) - i_pk - 1 >= POSTPEAK_STEPS):
            mode = "PEAK"
            verdict = (f"peak passed at s/B = {rows[i_pk][1]:.5f} and held "
                       f"{len(a) - i_pk - 1} steps more than "
                       f"{100 * POSTPEAK_DROP:.0f} % below it")
            break
    fh.close()
    wall = time.time() - t0

    assert len(rows) > 8, (
        f"{tag}: only {len(rows)} converged steps -- the push never started "
        f"({verdict})")

    arr = np.array([r[:5] for r in rows])
    s, qf, qb, dsm = arr[:, 0], arr[:, 2], arr[:, 3], arr[:, 4]
    qmax = float(qf.max())
    i_pk = int(qf.argmax())
    s_pk = float(s[i_pk])

    msk = s >= 0.9 * s[-1]
    t_last = (float(np.polyfit(s[msk], qf[msk], 1)[0]) if msk.sum() > 2
              else float("nan"))
    n0 = max(4, len(s) // 50)
    t_init = float(np.polyfit(s[:n0], qf[:n0], 1)[0])
    plateau = bool(abs(t_last) < PLATEAU_FRAC * abs(t_init))
    peaked = bool(i_pk < len(qf) - POSTPEAK_STEPS
                  and s[i_pk] / r3.B_FOOT >= PEAK_MIN_SFRAC
                  and float(qf[i_pk + 1:].max()) < (1.0 - POSTPEAK_DROP) * qmax)
    floor_mm = DS_MIN * 1000.0
    ds_tail_min_mm = float(dsm[msk].min()) if msk.sum() else float(dsm[-1])
    headroom = ds_tail_min_mm / floor_mm
    free = bool(headroom >= FREE_ADVANCE_FLOOR_FACTOR)
    capacity = bool((plateau or peaked) and free and mode in _CAPACITY_MODES)
    # Clause 2, stated positively AND negatively.  `seized` is not the negation
    # of `capacity`: a leg can fail to be a capacity for the ordinary reason
    # that it has not peaked yet, which is a MEASUREMENT that ran out of road,
    # while a seizure is a RUN that stopped for a reason outside the mechanics.
    # Both are inadmissible as a collapse load; only one of them says the number
    # is where a machine stopped.
    seized = mode in _SEIZURE_MODES

    # ---- end-of-run plastic strain field --------------------------------
    epsq = _plastic_field(n_hex)
    field_path = os.path.join(out_dir, f"a2_{tag}_field.csv")
    _dump_field(field_path,
                f"ADR90 WP-A2 END-OF-RUN plastic field, build "
                f"{EXPECTED_BUILD}, leg {tag}, s/B={s[-1] / r3.B_FOOT:.6f}, "
                f"mode {mode}",
                xc, zc, hx, hz, vol, epsq, eta_field)
    widths = widths_from_field(epsq, xc, zc, hx, hz, vol)
    n_yield = widths.pop("n_yield_ele")
    v_yield = widths.pop("vol_yield")
    epsq_max = widths.pop("epsq_max")

    # ---- engine-log guards ----------------------------------------------
    ops.logFile(os.path.join(out_dir, f"a2_{tag}_engine.after.log"), "-noEcho")
    try:
        txt = open(log_path, errors="ignore").read()
    except OSError:
        txt = ""
    n_outside = txt.count("Outside Bounding")
    n_clamp = txt.count("CLAMPING")

    res = dict(
        tag=tag, h0=h0, e_name=ename, e_init=e_init,
        build=EXPECTED_BUILD, driver=os.path.abspath(__file__),
        date=datetime.datetime.now().isoformat(timespec="seconds"),
        solver=solver, nodes=n_nodes, hexes=n_hex, dof=3 * n_nodes,
        gamma=GAMMA, K0=K0, M_c=M_C, presidual=OPT_PRESIDUAL, pmin=OPT_PMIN,
        ds_max=ds_max, ds_base=DS_BASE, ds_min=DS_MIN, push_tol=tol,
        push_test=test_type, push_tol_abs=tol_abs, force_ref=want,
        int_scheme=INT_SCHEME, tan_type=tan_type, jaco_type=JACO_TYPE,
        tol_f=TOL_F, tol_r=TOL_R,
        subdiv_budget=SUBDIV_BUDGET, sfrac_target=sfrac,
        wall_budget_s=wall_budget,
        # controls
        resultant_err=resultant_err, patch_err=patch_err,
        patch_err_zz=err_zz, patch_err_xx=err_xx,
        eta_max_grav=eta_max, eta_over_Mc=eta_max / M_C, p_min_grav=p_min,
        n_outside_bounding=n_outside, n_clamping=n_clamp,
        # measurements
        q_u=qmax, s_peak=s_pk, s_peak_over_B=s_pk / r3.B_FOOT,
        q_base_at_peak=float(qb[i_pk]),
        base_foot_mismatch=float(np.max(np.abs(qf + qb)) / max(qmax, 1e-12)),
        t_init=t_init, t_last=t_last,
        tail_pct=100.0 * t_last / t_init if t_init else float("nan"),
        plateau=plateau, peaked=peaked, free=free, capacity=capacity,
        seized=seized, admissible_as_capacity=bool(capacity),
        mode=mode, verdict=verdict,
        s_end_over_B=float(s[-1]) / r3.B_FOOT,
        ds_end_mm=float(dsm[-1]), ds_tail_min_mm=ds_tail_min_mm,
        ds_floor_mm=floor_mm, headroom=headroom,
        steps=len(rows), nfail=nfail, nrelax=nrelax, nsub=nsub,
        budget_used_frac=nsub / float(SUBDIV_BUDGET),
        n_yield_ele=n_yield, vol_yield=v_yield, epsq_max=epsq_max,
        wall_s=wall, wall_grav_s=t_grav, csv=csv_path, field=field_path,
        log=log_path, checkpoints=checkpoints,
    )
    res.update(widths)
    res["partial"] = False
    _write_leg_json(out_dir, tag, res)
    return res


def _why(r):
    return (f"[{r['tag']}] mode={r['mode']} ({r['verdict']}); ended s/B="
            f"{r['s_end_over_B']:.5f} of target {r['sfrac_target']:.3f}; "
            f"plateau={'yes' if r['plateau'] else 'NO'} / "
            f"peaked={'yes' if r['peaked'] else 'NO'} "
            f"(tail {r['tail_pct']:.3f} % of initial); "
            f"free-advance={'yes' if r['free'] else 'NO'} "
            f"({r['headroom']:.1f}x floor); subdivisions {r['nsub']}/"
            f"{r['subdiv_budget']}; {r['steps']} steps, {r['nfail']} failed, "
            f"{r['nrelax']} relaxed, {r['wall_s']:.0f} s")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--legs", default="",
                    help="comma list of h<h0>_<ename>; default all 6, coarse first")
    ap.add_argument("--wall", type=float, default=WALL_BUDGET_S)
    ap.add_argument("--sfrac", type=float, default=SFRAC)
    ap.add_argument("--dsmax", type=float, default=DS_MAX)
    ap.add_argument("--tol", type=float, default=None,
                    help="push ladder tolerance: relative to gamma*V for NormUnbalance, absolute metres for NormDispIncr")
    ap.add_argument("--test", default=PUSH_TEST,
                    choices=("NormUnbalance", "NormDispIncr"),
                    help="convergence test for the push ladder")
    ap.add_argument("--tantype", type=int, default=TAN_TYPE,
                    choices=(0, 1, 2),
                    help="ManzariDafalias TanType: 0 elastic (the PARSER DEFAULT, and a trap), 1 continuum ep, 2 consistent ep")
    args = ap.parse_args(argv)
    tol = args.tol if args.tol is not None else (
        PUSH_TOL if args.test == "NormUnbalance" else PUSH_TOL_DISP)
    os.makedirs(args.out, exist_ok=True)

    build = assert_engine()
    print("=== ADR 90 WP-A2 -- tau = 0 q_u band on softening SANISAND ===")
    print(f"    engine build (ladrunoBuild) : {build}")
    print(f"    openseespy module           : {os.path.abspath(ops.__file__)}")
    print(f"    driver                      : {os.path.abspath(__file__)}")
    print(f"    output                      : {os.path.abspath(args.out)}")
    print(f"    date                        : {datetime.datetime.now().isoformat(timespec='seconds')}")
    print(f"    gamma = {GAMMA} kN/m3, K0 = {K0:.4f}, M_c = {M_C}, "
          f"-Presidual {OPT_PRESIDUAL}, -Pmin {OPT_PMIN}")
    print(f"    SUBDIV_BUDGET (pinned, R3)  : {SUBDIV_BUDGET}")
    print(f"    DS_MAX (pinned here)        : {args.dsmax} m")
    print(f"    convergence test (pinned)   : {args.test} @ {tol}"
          + (" x gamma*V" if args.test == "NormUnbalance" else " m"))
    print(f"    TanType (0=elastic default) : {args.tantype}")
    print(f"    wall budget per leg         : {args.wall:.0f} s")

    want = [(h0, en, ei) for h0 in RESOLUTIONS for en, ei in DENSITIES]
    if args.legs:
        keep = {t.strip() for t in args.legs.split(",") if t.strip()}
        want = [t for t in want if leg_tag(t[0], t[1]) in keep]

    for h0, en, ei in want:
        print(f"\n--- leg {leg_tag(h0, en)} (e_init = {ei}) ---", flush=True)
        try:
            r = run_leg(h0, en, ei, args.out, wall_budget=args.wall,
                        sfrac=args.sfrac, ds_max=args.dsmax, tol=tol,
                        tan_type=args.tantype, test_type=args.test)
        except AssertionError as exc:
            print(f"    LEG FAILED A CONTROL: {exc}", flush=True)
            with open(os.path.join(args.out,
                                   f"a2_{leg_tag(h0, en)}.FAILED.txt"), "w") as f:
                f.write(str(exc))
            continue
        if r["seized"]:
            verdict = (f"*** INADMISSIBLE: terminated in SEIZURE mode "
                       f"{r['mode']} -- q_max = {r['q_u']:.3f} kPa is where the "
                       f"RUN stopped, not where the soil failed, and must not "
                       f"be quoted as a capacity ***")
        elif r["capacity"]:
            verdict = f"CAPACITY: q_u = {r['q_u']:.3f} kPa"
        else:
            verdict = (f"not a capacity (no plateau and no passed peak) -- "
                       f"q_max = {r['q_u']:.3f} kPa is a lower bound the leg "
                       f"reached, not a collapse load")
        print(f"    {verdict}", flush=True)
        print(f"    at s/B = {r['s_peak_over_B']:.5f}; {_why(r)}", flush=True)
        print(f"    controls: resultant {r['resultant_err']:.2e}, patch "
              f"{r['patch_err']:.2e}, eta/M_c {r['eta_over_Mc']:.4f}, "
              f"OutsideBounding {r['n_outside_bounding']}, "
              f"CLAMPING {r['n_clamping']}, base-foot {r['base_foot_mismatch']:.2e}",
              flush=True)
        print(f"    width: w2(z={Z_PROBES[0]}) = "
              f"{r[f'w2_z{Z_PROBES[0]}']:.4f} m, "
              f"w2(z={Z_PROBES[1]}) = {r[f'w2_z{Z_PROBES[1]}']:.4f} m, "
              f"yielding {r['n_yield_ele']} ele / {r['vol_yield']:.3f} m3",
              flush=True)
    print("\ndone; run sanisand_tau0_summary.py on the output directory")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
