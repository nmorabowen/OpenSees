"""ADR-97 P1 (wp/97b) — gate 6: every new failure path is LOUD, and every
refusal the ADR promised is actually enforced.

The class of defect ADR-94 found is not "wrong answer" but "wrong answer that
converges": a silent default to a different integrator, a dropped option, an
un-dispatched enum value, a return code no host element looks at.  So every new
code path in `Closest_Point` ends in `LADRUNO_MATERIAL_REFUSED` (never a bare
-1: `LadrunoBrick` compares against the sentinel and `stdBrick` drops
everything, ADR-94 B2), and every combination the ADR declared unsupported is
refused AT PARSE TIME, naming the tokens.

Covered here:

* a starved `n_max_iterations` on a plastic `Closest_Point` step fails the
  analysis instead of committing a non-converged state;
* `stdBrick` swallows that same refusal — pinned as a NEGATIVE control, because
  it is the reason every refusal test in this suite uses `LadrunoBrick`;
* `tangent_type Algorithmic` with any integrator other than `Closest_Point` is
  refused (ADR-97 D2 — a consistent tangent is defined only relative to a
  specific committed map);
* `integration_method Closest_Point` on a family P1 has not converted
  (a MIXED MohrCoulomb_YF x VonMises_PF pairing) is refused,
  naming the ADR phase that will deliver it -- MohrCoulomb x MohrCoulomb
  itself is SHIPPED by ADR-97 P2 and the assertion was inverted there;
* unknown tokens are still rejected (the ADR-94 contract);
* the ADR-94 B4 reproducer — Drucker-Prager driven into hydrostatic TENSION,
  past its apex — never commits a NaN under `Closest_Point`.

Zone-A, ~10 s.
"""
import subprocess  # Ladruno (ADR-97 wp/97f, D5): child-process stderr checks
import sys

import numpy as np
import pytest

from _testbed import ops

import test_adr97_p1_smooth as S  # noqa: E402
import test_asdplastic_mctc as M  # noqa: E402

pytestmark = [pytest.mark.zone_a]

MC_PARAMS = ["YoungsModulus", M.E, "PoissonsRatio", M.NU,
             "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
             "MassDensity", 0.0]


def _mat_mc(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters", *MC_PARAMS, "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


# Ladruno (ADR-97 wp/97c): the family gate is now about MIXED pairings, not about
# MohrCoulomb itself, and wp/97d ships HoekBrown.  `MohrCoulomb_YF` x `VonMises_PF` is one
# specializations the generator really registers and that P2 must still refuse.
def _mat_mc_mixed(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "VonMises_PF", "LinearIsotropic3D_EL",
            "BackStress(TensorLinearHardeningFunction):",
            "Begin_Model_Parameters", *MC_PARAMS,
            "TensorLinearHardeningParameter", 0.0, "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


def _mat_hb(tag, opts=None):
    args = ["ASDPlasticMaterial3D", tag,
            "HoekBrown_YF", "HoekBrown_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            # Ladruno (ADR-97 wp/97d): the parameter is `HB_sigci`, not
            # `HB_sigma_ci`.  With the wrong spelling this deck was rejected for
            # a MISSING PARAMETER, so the refusal assertion below passed without
            # ever exercising the family gate -- a vacuous test in both P1 and
            # P2.  Fixed here, and the assertion inverted (P3 ships HoekBrown).
            "HB_sigci", 50000.0, "HB_mb", 2.396510364418,
            "HB_s", 0.011743628457, "HB_a", 0.502840500848,
            "HB_mb_psi", 2.396510364418, "HB_ds", 0.0, "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables"]
    if opts is not None:
        args += ["Begin_Integration_Options", *opts, "End_Integration_Options"]
    ops.nDMaterial(*args)


@pytest.fixture(scope="module")
def cp_available():
    if not S._constructible(lambda t: S.mat_vm(t)):
        pytest.skip("ASDPlasticMaterial3D Closest_Point not available")


@pytest.fixture(scope="module")
def mc_available():
    if not S._constructible(lambda t: _mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D MohrCoulomb not available")


# ===========================================================================
# 1. parser refusals
# ===========================================================================
def test_closest_point_with_algorithmic_is_the_positive_control(cp_available):
    """The combination the ADR ships must build -- otherwise every refusal
    below would pass for the wrong reason."""
    assert S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Point", tangent="Algorithmic"))


@pytest.mark.parametrize("method", ["Backward_Euler", "Forward_Euler",
                                    "Modified_Euler_Error_Control",
                                    "Runge_Kutta_45_Error_Control"])
def test_algorithmic_is_refused_with_any_other_integrator(cp_available, method):
    """ADR-97 D2.  `tangent_type Algorithmic`'s enum value has existed since
    upstream with NO dispatch case anywhere -- it was dead, which is why nobody
    has been silently getting it.  Offering it on the cutting-plane
    `Backward_Euler` would ship a FOURTH almost-right tangent, which is exactly
    the class of defect ADR-94 M3 found."""
    assert not S._constructible(
        lambda t: S.mat_vm(t, method=method, tangent="Algorithmic")), (
        "tangent_type Algorithmic was accepted with integration_method %s -- "
        "the ADR-97 D2 cross-refusal is gone" % method)


def test_algorithmic_refusal_does_not_break_the_other_tangents(cp_available):
    """The refusal must be surgical: every tangent_type the chosen integrator
    does define still builds."""
    for tg in ("Secant", "Continuum", "Elastic",
               "Numerical_Algorithmic_FirstOrder",
               "Numerical_Algorithmic_SecondOrder"):
        assert S._constructible(
            lambda t, g=tg: S.mat_vm(t, method="Backward_Euler", tangent=g)), tg
        assert S._constructible(
            lambda t, g=tg: S.mat_vm(t, method="Closest_Point", tangent=g)), tg


def test_closest_point_is_refused_for_unconverted_families(mc_available):
    """ADR-97 D3: the closest-point map is added family by family.  A
    specialization whose yield function, plastic flow direction or hardening law
    has not opted in must be refused at parse time -- NOT silently run on the
    inert zero / finite-difference defaults, which is the ADR-94 M3 failure mode.

    Ladruno (ADR-97 wp/97c): this test used to assert that MohrCoulomb ITSELF was
    refused.  P2 ships it (a principal-stress-space multi-surface return; see
    tests/test_adr97_p2_principal.py), so the MohrCoulomb half of the assertion
    is INVERTED here -- deliberately, and recorded in the ADR-97 P2 report.  What
    the family gate must still refuse is a MIXED pairing: `MohrCoulomb_YF` with a
    non-Mohr-Coulomb plastic flow direction is covered by NO oracle (the
    principal-space return assumes both the surface and the potential are
    piecewise linear, and P1's smooth 6D map cannot use MohrCoulomb's Lode-angle
    gradient), so it stays refused.

    Ladruno (ADR-97 wp/97d): the HoekBrown half is now INVERTED too -- P3 ships
    `HoekBrown_YF` x `HoekBrown_PF` (a principal-space return to the CURVED
    surface; see tests/test_adr97_p3_hoekbrown.py).  The mixed HoekBrown
    pairings stay refused and are covered there."""
    assert S._constructible(lambda t: _mat_mc(t)), "control MC deck must build"
    assert S._constructible(
        lambda t: _mat_mc(t, opts=["integration_method", "Closest_Point"])), (
        "integration_method Closest_Point is refused for MohrCoulomb x "
        "MohrCoulomb -- ADR-97 P2 ships it")
    assert not S._constructible(
        lambda t: _mat_mc_mixed(t, opts=["integration_method",
                                         "Closest_Point"])), (
        "integration_method Closest_Point was accepted for the MIXED pairing "
        "MohrCoulomb_YF x VonMises_PF, which no oracle covers -- the ADR-97 D3 "
        "family gate is gone")
    assert S._constructible(lambda t: _mat_hb(t)), (
        "the control HoekBrown deck does not build at all -- fix its parameter "
        "list before reading anything into the assertion below")
    assert S._constructible(
        lambda t: _mat_hb(t, opts=["integration_method", "Closest_Point"])), (
        "integration_method Closest_Point is refused for HoekBrown_YF x "
        "HoekBrown_PF -- ADR-97 P3 ships it")


def test_unknown_tokens_are_still_rejected(cp_available):
    """The ADR-94 wp/94a contract: an unrecognised `integration_method` or
    `tangent_type` is an ERROR, not a silent default."""
    assert not S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Pont", tangent="Secant"))
    assert not S._constructible(
        lambda t: S.mat_vm(t, method="Closest_Point", tangent="Algorythmic"))


# ===========================================================================
# 2. runtime refusals
# ===========================================================================
def _starved(ele):
    """A plastic Closest_Point step with `n_max_iterations` far too small."""
    return S.drive(lambda t: S.mat_vm(t, hiso=7000.0, niter=1),
                   S.VM_PATHS["triaxial"], nstep=10, ele=ele)


def test_starved_newton_is_refused_not_committed(cp_available):
    """`Closest_Point` returns LADRUNO_MATERIAL_REFUSED on exhaustion; the host
    element propagates it and the analysis fails.  Upstream's habit -- and
    `Backward_Euler`'s before ADR-84 P2a -- was to fall out of the loop and
    commit the non-converged state as success."""
    r = _starved("LadrunoBrick")
    print("starved n_max_iterations=1 codes:", r["codes"])
    assert any(c != 0 for c in r["codes"]), (
        "a starved Closest_Point Newton was committed as success")


def test_stdbrick_swallows_the_refusal_negative_control(cp_available):
    """PINNED negative control (ADR-94 B2).  `stdBrick` discards the material's
    return code by design, so the SAME starved deck 'succeeds' on it.  This is
    why every refusal test in the ADR-97 suite runs on `LadrunoBrick`."""
    r = _starved("stdBrick")
    print("starved on stdBrick codes:", r["codes"])
    assert all(c == 0 for c in r["codes"]), (
        "stdBrick now propagates material return codes -- good news, but the "
        "ADR-94 B2 pin and the ADR-97 host-choice rationale both need updating")


def test_strict_convergence_is_accepted_and_inert_on_a_converging_deck(
        cp_available):
    """`Closest_Point` ends every commit path in the existing
    `ladruno_strict_rejects` gate, on the same contract `Backward_Euler` has.
    On a deck that converges, turning the flag on must change nothing."""
    off = S.drive(lambda t: S.mat_vm(t, hiso=7000.0),
                  S.VM_PATHS["triaxial"], nstep=10)
    on = S.drive(lambda t: S.mat_vm(t, hiso=7000.0, strict=1),
                 S.VM_PATHS["triaxial"], nstep=10)
    assert off["codes"] == on["codes"] == [0] * 10
    d = S._rel(on["sigma"][-1], off["sigma"][-1])
    print("strict_convergence on/off rel difference = %.3e" % d)
    assert d == 0.0, "strict_convergence changed a converging Closest_Point deck"


def test_strict_convergence_still_refuses_a_starved_deck(cp_available):
    r = S.drive(lambda t: S.mat_vm(t, hiso=7000.0, niter=1, strict=1),
                S.VM_PATHS["triaxial"], nstep=10)
    assert any(c != 0 for c in r["codes"])


def test_dp_hydrostatic_tension_never_commits_nan_under_closest_point():
    """ADR-94 B4's reproducer, re-run on the new integrator.  Pure hydrostatic
    TENSION drives Drucker-Prager straight past its apex, where the flank return
    would have to make `sqrt(J2)` negative.  `Closest_Point` classifies the apex
    region in the ELASTIC metric and returns to the vertex; nothing NaN may
    reach the committed state, and nothing inadmissible either."""
    if not S._constructible(lambda t: S.mat_dp(t)):
        pytest.skip("DruckerPrager Closest_Point unavailable")
    legs = [np.array([2.4e-3, 2.4e-3, 2.4e-3, 0., 0., 0.])]   # 2x the apex strain
    r = S.drive(lambda t: S.mat_dp(t), legs, nstep=10, want=("stresses",
                                                             "BackStress"))
    assert all(c == 0 for c in r["codes"]), r["codes"]
    assert np.all(np.isfinite(r["sigma"])), r["sigma"]
    p_apex = S.DP_XI_C / S.DP_ETA
    print("DP hydrostatic tension: final sigma = %s (apex p = %.6f)"
          % (np.array2string(r["sigma"][-1], precision=8), p_apex))
    assert abs(float(np.mean(r["sigma"][-1][:3])) - p_apex) < 1e-6 * p_apex
    for i in range(len(r["sigma"])):
        f = S._f_dp(r["sigma"][i], r["BackStress"][i])
        assert f <= 1e-6, (i, f)


# ===========================================================================
# ADR-97 wp/97f (D5): experimental_integrator gate
# ===========================================================================
#: The four integrators D5 gates. `Backward_Euler_LineSearch` and
#: `Runge_Kutta_45_Error_Control_old` are OUT of scope here -- they stay
#: refused outright (ADR-94 M7/M8), with or without the flag.
EXPLICIT_METHODS = ["Forward_Euler", "Forward_Euler_Subincrement",
                     "Modified_Euler_Error_Control",
                     "Runge_Kutta_45_Error_Control"]


def _mat_vm_raw(tag, method="Backward_Euler", tangent="Secant",
                experimental=None, niter=100):
    """Same VonMises deck as ``test_adr97_p1_smooth.mat_vm``, but with a
    direct ``experimental_integrator`` knob -- ``S.mat_vm`` has no such
    parameter."""
    extra = (["experimental_integrator", int(experimental)]
             if experimental is not None else [])
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", S.IV_VM_LIN,
        "Begin_Model_Parameters",
        "YoungsModulus", S.E_VM, "PoissonsRatio", S.NU_VM,
        "ScalarLinearHardeningParameter", 0.0,
        "TensorLinearHardeningParameter", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables", "YieldStress", S.SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter), *extra,
        "End_Integration_Options",
    )


@pytest.mark.parametrize("method", EXPLICIT_METHODS)
def test_explicit_integrator_refused_without_the_flag(cp_available, method):
    """ADR-97 D5.  Without ``experimental_integrator 1`` the four explicit
    integrators are refused at parse time -- ADR-94 M2/M8 found real defects
    in them (empty drift-correction bodies, an error controller that accepts
    unconditionally at ``rk45_dT_min``) and rewriting them is a different ADR;
    this is the gate, not the fix."""
    assert not S._constructible(
        lambda t: _mat_vm_raw(t, method=method, tangent="Secant")), (
        "integration_method %s was accepted with experimental_integrator "
        "unset -- the ADR-97 D5 gate is gone" % method)


@pytest.mark.parametrize("method", EXPLICIT_METHODS)
def test_explicit_integrator_accepted_with_the_flag(cp_available, method):
    """The opt-in half of the same gate: setting the flag must let every one
    of the four explicit integrators build (D5 gates SELECTION, it does not
    remove the integrators)."""
    assert S._constructible(
        lambda t: _mat_vm_raw(t, method=method, tangent="Secant",
                              experimental=1)), (
        "integration_method %s was refused even WITH "
        "experimental_integrator 1 -- the opt-in half of ADR-97 D5 is "
        "broken" % method)


@pytest.mark.parametrize("method", EXPLICIT_METHODS)
def test_explicit_integrator_experimental_flag_zero_still_refuses(
        cp_available, method):
    """Both directions, explicitly: `experimental_integrator 0` (the
    documented default value, spelled out rather than omitted) must refuse
    exactly like leaving the option out entirely."""
    assert not S._constructible(
        lambda t: _mat_vm_raw(t, method=method, tangent="Secant",
                              experimental=0)), (
        "integration_method %s was accepted with experimental_integrator 0 "
        "spelled out explicitly" % method)


def test_backward_euler_and_closest_point_unaffected_by_the_flag(
        cp_available):
    """ADR-97 D1/D5: the gate is scoped to the four explicit integrators
    only.  Backward_Euler and Closest_Point must build identically whether
    experimental_integrator is left unset, 0, or 1 -- byte-identity itself is
    gate 4 (tests/test_adr97_p4_inertness.py); this is the narrower
    constructibility half."""
    for method, tangent in (("Backward_Euler", "Secant"),
                            ("Closest_Point", "Algorithmic")):
        for experimental in (None, 0, 1):
            assert S._constructible(
                lambda t, m=method, g=tangent, e=experimental:
                    _mat_vm_raw(t, method=m, tangent=g, experimental=e)), (
                "%s was refused with experimental_integrator=%r -- the D5 "
                "gate is no longer scoped to the four explicit integrators "
                "only" % (method, experimental))


def test_algorithmic_still_refused_on_an_opted_in_explicit_integrator(
        cp_available):
    """Closes a coverage gap the D5 gate would otherwise open: with
    experimental_integrator=1,
    ``test_algorithmic_is_refused_with_any_other_integrator`` above stops
    isolating the ADR-97 D2 tangent/integrator cross-check for three of its
    four parametrized methods, because the D5 gate now refuses those decks
    FIRST, for an unrelated reason -- exactly the class of "refused, but for
    the wrong reason" defect the ADR-97 P3 report's HB_sigma_ci finding
    warns about. Opting in here removes the D5 refusal so the assertion
    below is unambiguously exercising D2."""
    for method in EXPLICIT_METHODS:
        assert not S._constructible(
            lambda t, m=method: _mat_vm_raw(t, method=m, tangent="Algorithmic",
                                            experimental=1)), (
            "tangent_type Algorithmic was accepted with integration_method "
            "%s (even opted in via experimental_integrator 1) -- the ADR-97 "
            "D2 cross-refusal is gone" % method)


def test_experimental_integrator_typo_is_still_rejected(cp_available):
    """The ADR-94 wp/94a unknown-token contract extends to the new option
    name: a misspelling is an ERROR, not a silently-ignored, still-refused
    deck that happens to look the same."""
    def _typo_deck(t):
        ops.nDMaterial(
            "ASDPlasticMaterial3D", t,
            "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", S.IV_VM_LIN,
            "Begin_Model_Parameters",
            "YoungsModulus", S.E_VM, "PoissonsRatio", S.NU_VM,
            "ScalarLinearHardeningParameter", 0.0,
            "TensorLinearHardeningParameter", 0.0, "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables", "YieldStress", S.SY_VM,
            "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
            "Begin_Integration_Options",
            "integration_method", "Forward_Euler",
            "experimental_integratr", 1,   # deliberate typo
            "End_Integration_Options")
    assert not S._constructible(_typo_deck), (
        "a misspelled 'experimental_integratr' was silently accepted -- the "
        "ADR-94 wp/94a unknown-token contract no longer covers the new D5 "
        "option")


def _child_script(tests_dir, method, experimental):
    """Build a standalone script that constructs a VonMises ASDPlasticMaterial3D
    deck with the given integration_method (and, if not None, an explicit
    ``experimental_integrator`` value), then prints whether construction
    succeeded. Run in a FRESH interpreter -- see ``_run_child``'s docstring
    for why this cannot be done with pytest's own capfd."""
    lines = [
        'import sys; sys.path.insert(0, %r)' % tests_dir,
        'from _testbed import ops',
        'ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)',
        'CUBE = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]',
        'for k, c in enumerate(CUBE):',
        '    ops.node(k + 1, *map(float, c))',
        'rc = None',
        'try:',
        '    ops.nDMaterial(',
        '        "ASDPlasticMaterial3D", 1,',
        '        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL",',
        '        "BackStress(TensorLinearHardeningFunction):",',
        '        "Begin_Model_Parameters",',
        '        "YoungsModulus", 70000.0, "PoissonsRatio", 0.3,',
        '        "ScalarLinearHardeningParameter", 0.0,',
        '        "TensorLinearHardeningParameter", 0.0, "MassDensity", 0.0,',
        '        "End_Model_Parameters",',
        '        "Begin_Internal_Variables", "YieldStress", 30.0,',
        '        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",',
        '        "Begin_Integration_Options",',
        '        "integration_method", %r, "tangent_type", "Secant",' % method,
    ]
    if experimental is not None:
        lines.append('        "experimental_integrator", %d,' % int(experimental))
    lines += [
        '        "End_Integration_Options",',
        '    )',
        '    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)',
        '    rc = "CONSTRUCTED"',
        'except Exception as e:',
        '    rc = "RAISED " + str(e)',
        'print("RESULT", rc)',
    ]
    return chr(10).join(lines) + chr(10)


def _run_child(script, timeout=60):
    """Run ``script`` in a fresh interpreter and return the completed run.

    Same MEASURED TRAP as ``test_adr94_hlist_mechanical._run_child``: pytest's
    ``capfd`` cannot see anything the .pyd writes via ``cout``/``opserr`` on
    this build (a DLL-boundary fd-duplication issue), so the D5 refusal
    message can only be asserted on a CHILD process's real OS-level
    stdout/stderr.
    """
    import os
    env = dict(os.environ)
    dist_dir = os.path.dirname(os.path.abspath(ops.__file__))
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    env["PYTHONPATH"] = dist_dir + os.pathsep + tests_dir + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = dist_dir + os.pathsep + env.get("PATH", "")
    return subprocess.run([sys.executable, "-c", script], capture_output=True,
                          text=True, env=env, timeout=timeout,
                          cwd=tests_dir, stdin=subprocess.DEVNULL)


def test_explicit_integrator_refusal_message_on_real_stderr(cp_available):
    """Child-process check that the refusal is actually LOUD: the real
    process stdout/stderr (not pytest's capfd, which cannot see this .pyd's
    output -- see ``_run_child``) must name ADR-97 D5 and the supported
    implicit pair when an explicit integrator is picked without the flag."""
    import os
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    script = _child_script(tests_dir, "Forward_Euler", None)
    proc = _run_child(script)
    combined = proc.stdout + proc.stderr
    assert "RESULT RAISED" in combined, (
        f"Forward_Euler without experimental_integrator was NOT refused in "
        f"a fresh process (stdout={proc.stdout!r}, stderr={proc.stderr!r})")
    assert "ADR-97 D5" in combined, (
        f"the D5 refusal message did not reach real stdout/stderr "
        f"(stdout={proc.stdout!r}, stderr={proc.stderr!r})")
    assert "Backward_Euler" in combined and "Closest_Point" in combined, (
        f"the D5 refusal message no longer names the supported implicit "
        f"pair (stdout={proc.stdout!r}, stderr={proc.stderr!r})")


def test_explicit_integrator_opt_in_reaches_construction_on_real_process(
        cp_available):
    """No-regression twin, same child-process rig: WITH the flag, the same
    deck actually constructs (real process, not just the in-process
    ``_constructible`` helper)."""
    import os
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    script = _child_script(tests_dir, "Forward_Euler", 1)
    proc = _run_child(script)
    combined = proc.stdout + proc.stderr
    assert "RESULT CONSTRUCTED" in combined, (
        f"Forward_Euler WITH experimental_integrator 1 did not construct in "
        f"a fresh process (stdout={proc.stdout!r}, stderr={proc.stderr!r})")
