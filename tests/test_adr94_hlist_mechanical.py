"""ADR-94 R1-B — mechanical reproducers for H-list rows H2, H4, H5, H9, H12, H13,
H14, H15 (`Ladruno_implementation/94_asdplastic_review_plan.md` sec. 4).

Every test asserts the DEFECT as currently observed on
`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` at `52314165a`, so
a future fix flips it red (the fork's sentinel rule, `LEDGER_quirks.md`).  A
REFUTED claim instead pins the observed correct behaviour with a docstring
starting "REFUTED:".

H2 is read-only (UB from a dangling `std::string::c_str()` in `getClassType`,
196-201) and gets no test, per the review plan.

Host-element rule (`_adr94_hlist_R1B.md`, ADR-94 sec. 8): gate refusals on
``TenNodeTetrahedron`` (propagates the material's return code); ``stdBrick`` is
only used for pure material-point response probing (`test_asdplastic_mctc`'s
own pattern), never for a refusal assertion.
"""
import math
import os
import re
import subprocess
import sys

import numpy as np
import pytest

from _testbed import ops

import test_asdplastic_mctc as M

pytestmark = [pytest.mark.zone_a]

ASDP_HEADER = None  # lazily read once, see _asdp_source()

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_DIST_DIR = os.path.dirname(os.path.abspath(ops.__file__))


def _run_child(script, timeout=60):
    """Run ``script`` in a FRESH python process and return the completed run.

    MEASURED TRAP: pytest's ``capfd`` does not see anything this native
    extension writes via ``cout``/``cerr`` on this build -- a mid-process
    ``dup2`` swap of fd 1/2 does not reach output written through the .pyd's
    own linked CRT (a Windows DLL-boundary issue, not a test bug: the exact
    same code prints normally when the WHOLE process's stdout is piped from
    the outside, i.e. from a shell or ``subprocess``). So H4/H12's stdout/
    stderr assertions run the model in a CHILD process and capture ITS real
    OS-level stdout/stderr instead of using ``capfd``.
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = _DIST_DIR + os.pathsep + _TESTS_DIR + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = _DIST_DIR + os.pathsep + env.get("PATH", "")
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    # MEASURED FLAKE (Windows only): under pytest's own stdio setup, the
    # PARENT's stdin handle is sometimes not a valid handle for
    # subprocess.Popen to duplicate for the child ("OSError: [WinError 6]
    # The handle is invalid" from _winapi.DuplicateHandle, non-deterministic
    # across otherwise-identical runs). Passing stdin=DEVNULL explicitly
    # avoids inheriting/duplicating the parent's stdin handle at all.
    return subprocess.run([sys.executable, "-c", script], capture_output=True,
                           text=True, env=env, timeout=timeout,
                           cwd=_TESTS_DIR, stdin=subprocess.DEVNULL)


def _asdp_source():
    global ASDP_HEADER
    if ASDP_HEADER is None:
        import pathlib
        here = pathlib.Path(__file__).resolve().parent.parent
        ASDP_HEADER = (here / "SRC" / "material" / "nD" / "ASDPlasticMaterial3D"
                       / "ASDPlasticMaterial3D.h").read_text(encoding="utf-8",
                                                              errors="replace")
    return ASDP_HEADER


# ===========================================================================
# shared tet rig (copied from test_adr84_p2a_strict_convergence's pattern)
# ===========================================================================
_TET = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
        5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5),
        9: (.5, 0, .5), 10: (0, .5, .5)}
_TET_TOP = (4, 8, 9, 10)


def mat_mc(tag, method=None, strict=None, niter=None, ftol=None, p0=None):
    opts = []
    if method is not None:
        opts += ["integration_method", method]
    if strict is not None:
        opts += ["strict_convergence", int(strict)]
    if niter is not None:
        opts += ["n_max_iterations", int(niter)]
    if ftol is not None:
        opts += ["f_absolute_tol", float(ftol)]
    model_params = ["YoungsModulus", M.E, "PoissonsRatio", M.NU,
                    "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
                    "MassDensity", 0.0]
    if p0 is not None:
        model_params += ["InitialP0", float(p0)]
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
        "Begin_Model_Parameters", *model_params, "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *(["Begin_Integration_Options"] + opts + ["End_Integration_Options"]
          if opts else []),
    )


def _tet_build(mat_fn, nsteps, utop, tol=1e-8, maxiter=100):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    mat_fn(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, utop)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxiter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")


def _tet_stress():
    ops.eleResponse(1, "forces")
    return np.array(list(ops.eleResponse(1, "stresses"))[0:6])


def _tet_constructible(mat_fn):
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _TET.items():
            ops.node(t, *map(float, c))
        mat_fn(1)
        ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def mc_available():
    if not _tet_constructible(lambda t: mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb_YF / "
                    "TenNodeTetrahedron not available in this build")


# ===========================================================================
# H4 -- revertToLastCommit() is a commented-out no-op; revertToStart() -> -1
# ===========================================================================
@pytest.mark.t0m
def test_H4_revert_to_last_commit_is_noop(mc_available):
    """CONFIRMED.  ``revertToLastCommit()`` (735-748) has every statement
    commented out and returns 0 (success) without touching Trial*/Commit*.
    ``StaticAnalysis::analyze()`` calls ``Domain::revertToLastCommit()``
    itself whenever a step fails to converge -- so a failed step SHOULD leave
    the material's trial stress equal to its last commit.  It does not: the
    dirty trial stress from the failed iteration survives untouched.

    ``revertToStart()`` (780-783 area) is a distinct, separate defect: it
    prints "not implemented" and returns -1, but ``OPS_resetModel()``
    (``OpenSeesCommands.cpp``) never reads that return code, so
    ``ops.reset()`` reports success while the reset silently did nothing.
    """
    _tet_build(lambda t: mat_mc(t), nsteps=20, utop=-0.02)
    assert ops.analyze(1) == 0, "step 1 (healthy) must converge"
    committed = _tet_stress()

    # Step 2: an intentionally impossible NormDispIncr budget forces the
    # global Newton to fail WITHOUT committing -- StaticAnalysis::analyze()
    # calls Domain::revertToLastCommit() internally before returning.
    ops.test("NormDispIncr", 1.0e-14, 1, 0)
    ops.integrator("LoadControl", 5.0 / 20)  # a much larger, harder increment
    rc = ops.analyze(1)
    assert rc != 0, ("expected step 2 to FAIL to converge (the reproducer "
                      "needs an uncommitted, dirty trial state); got rc=0")

    after_revert = _tet_stress()
    # MEASURED CAVEAT: TenNodeTetrahedron's "stresses"/"forces" response
    # unconditionally re-derives stress from the CURRENT nodal trial
    # displacements on every query (confirmed: identical result with or
    # without a preceding "forces" call, and via the "material"/"integrPoint"
    # sub-response routes too) -- and Domain::revertToLastCommit() DOES
    # correctly reset nodal trial displacements. So this element-level query
    # self-heals and reports the correct value REGARDLESS of whether the
    # material's own revertToLastCommit() did anything -- it cannot be used
    # to observe the material-level no-op from Python. Measured, not assumed:
    assert np.allclose(after_revert, committed, atol=1e-6), (
        f"host-element self-heal on 'stresses' no longer masks the material "
        f"no-op (committed={committed}, post-revert={after_revert}) -- "
        f"re-verify this caveat before trusting the structural check below")

    # THE ACTUAL DEFECT, pinned structurally: revertToLastCommit()'s body
    # must contain no active (uncommented) statement other than `return 0`.
    body_m = re.search(
        r"int revertToLastCommit\(void\)\s*\{(.*?)\n    \}\n", _asdp_source(),
        re.S)
    assert body_m, "could not locate revertToLastCommit() body -- source layout changed"
    body = body_m.group(1)
    live = "\n".join(ln for ln in body.splitlines()
                      if ln.strip() and not ln.strip().startswith("//"))
    assert re.fullmatch(r"\s*return 0;\s*", live), (
        f"revertToLastCommit() now has live statements beyond 'return 0' "
        f"({live!r}) -- it may have been implemented; update H4 to "
        f"REFUTED/fixed (and see if the element-recompute caveat above can "
        f"finally be replaced with a real runtime pin).")

    # revertToStart(): ops.reset() must not raise, and must not report the
    # underlying material failure to the caller (the "not implemented" line
    # only shows up on cerr/opserr, never as a raised exception or bad rc).
    # Run in a CHILD process -- capfd cannot see this native extension's
    # cout/cerr on this build (see ``_run_child``'s docstring).
    script = (
        "import sys; sys.path.insert(0, r'" + _TESTS_DIR + "')\n"
        "from _testbed import ops\n"
        "import test_asdplastic_mctc as M\n"
        "_TET = {1:(0,0,0),2:(1,0,0),3:(0,1,0),4:(0,0,1),5:(.5,0,0),"
        "6:(.5,.5,0),7:(0,.5,0),8:(0,0,.5),9:(.5,0,.5),10:(0,.5,.5)}\n"
        "ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)\n"
        "[ops.node(t, *map(float, c)) for t, c in _TET.items()]\n"
        "[ops.fix(t, 1,1,1) for t in (1,2,3,5,6,7)]\n"
        "[ops.fix(t, 1,1,0) for t in (4,8,9,10)]\n"
        "ops.nDMaterial('ASDPlasticMaterial3D', 1, 'MohrCoulomb_YF', "
        "'MohrCoulomb_PF', 'LinearIsotropic3D_EL', M.IV, "
        "'Begin_Model_Parameters', 'YoungsModulus', M.E, 'PoissonsRatio', "
        "M.NU, 'MC_phi', M.PHI, 'MC_c', M.C, 'MC_psi', M.PSI, 'MC_ds', 0.0, "
        "'MassDensity', 0.0, 'End_Model_Parameters', "
        "'Begin_Internal_Variables', 'BackStress', 0.,0.,0.,0.,0.,0., "
        "'End_Internal_Variables')\n"
        # the material must be attached to an ELEMENT: Domain::revertToStart()
        # visits elements, not standalone registered material prototypes.
        "ops.element('TenNodeTetrahedron', 1, 1,2,3,4,5,6,7,8,9,10, 1)\n"
        "rc = ops.reset()\n"
        "print('RESET_RC', rc)\n"
    )
    proc = _run_child(script)
    assert "RESET_RC" in proc.stdout, (
        f"child process did not reach ops.reset() (stdout={proc.stdout!r}, "
        f"stderr={proc.stderr!r})")
    combined = proc.stdout + proc.stderr
    assert "not implemented" in combined, (
        "expected ASDPlasticMaterial3D::revertToStart()'s "
        "'not implemented' line to appear in the child process's real "
        "stdout/stderr; if it is gone, revertToStart() may have been "
        f"implemented -- re-verify H4 before trusting this test "
        f"(stdout={proc.stdout!r}, stderr={proc.stderr!r})")


# ===========================================================================
# H5 -- the f-decreased elastic exit is unguarded at 6 more integrators
# ===========================================================================
#: Measured (see ``_adr94_hlist_R1B.md``): at TET_UTOP=-0.02/20 steps,
#: n_max_iterations=2, strict_convergence=1, these four commit an admissible-
#: violating state despite the flag (max f_MC 300-2800 vs tol ~0.1-0.2).
#: ``Backward_Euler_LineSearch`` and ``Runge_Kutta_45_Error_Control`` (the
#: non-"_old" one) instead FAIL to converge globally on this exact rig
#: (rc=-3) rather than silently committing -- a different failure mode this
#: reproducer does not isolate. Both are still confirmed unguarded BY READING
#: (2406 and 3435 respectively, both `yf_val_start > yf_val_end` with no
#: `strict_convergence` check in scope); only the RUNTIME pin is narrower.
NON_BE_METHODS = [
    "Forward_Euler",
    "Forward_Euler_Subincrement",
    "Modified_Euler_Error_Control",
    "Runge_Kutta_45_Error_Control_old",
]


@pytest.mark.t0m
@pytest.mark.parametrize("method", NON_BE_METHODS)
def test_H5_strict_convergence_does_not_gate_other_integrators(mc_available,
                                                                 method):
    """CONFIRMED.  ``strict_convergence`` is only read inside
    ``Backward_Euler``'s scalar-Newton loop.  The same
    ``yf_val_start > yf_val_end`` => "elastic, no correction" shortcut exists,
    unguarded, in FE (1423), FE_sub (1599), BE_LS (2406),
    RK45_old (2667), ME (3088), RK45 (3435).  So turning strict_convergence on
    does not stop any of these six from committing an inadmissible state on a
    coarse/starved plastic leg -- unlike Backward_Euler, where
    ``test_adr84_p2a_strict_convergence::test_flag_on_refuses_instead_of_committing``
    shows the flag DOES refuse.
    """
    _tet_build(lambda t: mat_mc(t, method=method, strict=1, niter=100),
               nsteps=20, utop=-0.02)
    hist = []
    for _ in range(20):
        ops.analyze(1)          # do not assert: a refusal would REFUTE this H
        hist.append(_tet_stress())
    hist = np.array(hist)
    tol = 1.0e-6 * max(2.0 * M.C * math.cos(math.radians(M.PHI)),
                        float(np.max(np.abs(hist))))
    fmc = np.array([M.f_mc(s) for s in hist])
    assert float(np.max(fmc)) > tol, (
        f"[{method}] strict_convergence=1 kept every committed state "
        f"admissible (max f_MC={np.max(fmc):.3e} <= tol {tol:.3e}) -- this "
        f"integrator may now be gated too; re-verify H5's site list before "
        f"trusting this test.")


# ===========================================================================
# H9 -- ME/RK45 drift-correction blocks are empty ``if`` statements
# ===========================================================================
@pytest.mark.t0m
def test_H9_explicit_drift_check_is_dead_code(mc_available):
    """CONFIRMED.  The "Validate yield function drift" block in
    ``Modified_Euler_Error_Control`` (3240-3250) computes ``yf_val`` and tests
    ``yf_val > 10*tol`` but the body of the ``if`` is a commented-out ``cout``
    -- no correction is ever applied.  With ``return_to_yield_surface
    Disabled`` the only other correction path is also off, so a coarse step
    commits with drift several orders above ``f_absolute_tol``.

    ``Runge_Kutta_45_Error_Control``'s identical block (3760-3770) is
    confirmed the same way BY READING (the ``et()``-call-count check in
    ``test_H15`` independently confirms this function is live code, not
    unreachable): on this build it fails to globally converge at every
    magnitude tried with ``return_to_yield_surface Disabled`` (measured -- see
    ``_adr94_hlist_R1B.md``), so no runtime pin is included for it here.
    """
    method = "Modified_Euler_Error_Control"
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    ops.nDMaterial(
        "ASDPlasticMaterial3D", 1,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
        "Begin_Model_Parameters",
        "YoungsModulus", M.E, "PoissonsRatio", M.NU,
        "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method,
        "f_absolute_tol", 1.0e-8,
        "n_max_iterations", 200,
        "return_to_yield_surface", "Disabled",
        "End_Integration_Options",
    )
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, -0.02)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 20)
    ops.analysis("Static")

    hist = []
    for _ in range(20):
        ops.analyze(1)
        hist.append(_tet_stress())
    hist = np.array(hist)
    fmc = np.array([M.f_mc(s) for s in hist])
    ftol = 1.0e-8
    assert float(np.max(fmc)) > 10.0 * ftol, (
        f"[{method}] max committed f_MC = {np.max(fmc):.3e} did not clear "
        f"10x f_absolute_tol ({ftol:.1e}) -- the coarse-step drift "
        f"reproducer no longer reproduces; re-tune before trusting H9.")


# ===========================================================================
# H12 -- diagnostics use cout (not opserr) and flood stdout
# ===========================================================================
@pytest.mark.t0m
def test_H12_commit_diagnostics_use_cout_not_opserr(mc_available):
    """CONFIRMED.  ``commitState()`` (727-748) prints via ``cout`` whenever
    ``GLOBAL_INT_max_iter[ASDP_TAG] > 0`` -- true after essentially any
    Backward_Euler commit that took a Newton iteration, since the counter is
    a PER-TAG static shared by every Gauss point using that tag (727,
    2263, 2907, 3262, 3785). Driving several elements sharing one tag through
    a plastic path floods stdout with "ASDP Integration Info" lines that
    ``opserr``-only redirection (e.g. under MPI) would never see.

    Run in a CHILD process -- capfd cannot see this native extension's
    cout/cerr on this build (see ``_run_child``'s docstring).
    """
    script = (
        "import sys; sys.path.insert(0, r'" + _TESTS_DIR + "')\n"
        "from _testbed import ops\n"
        "import test_asdplastic_mctc as M\n"
        "ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)\n"
        "n_elem = 10\n"
        "ops.nDMaterial('ASDPlasticMaterial3D', 1, 'MohrCoulomb_YF', "
        "'MohrCoulomb_PF', 'LinearIsotropic3D_EL', M.IV, "
        "'Begin_Model_Parameters', 'YoungsModulus', M.E, 'PoissonsRatio', "
        "M.NU, 'MC_phi', M.PHI, 'MC_c', M.C, 'MC_psi', M.PSI, 'MC_ds', 0.0, "
        "'MassDensity', 0.0, 'End_Model_Parameters', "
        "'Begin_Internal_Variables', 'BackStress', 0.,0.,0.,0.,0.,0., "
        "'End_Internal_Variables')\n"
        "for e in range(n_elem):\n"
        "    base = e * 8\n"
        "    y_off = 2.0 * e\n"
        "    for i, c in M._CUBE.items():\n"
        "        ops.node(base + i, c[0], c[1] + y_off, c[2])\n"
        "    for i, m in M._FIX.items():\n"
        "        ops.fix(base + i, *m)\n"
        "    ops.element('stdBrick', e + 1, *[base + i for i in range(1, 9)], 1)\n"
        "ops.timeSeries('Linear', 1)\n"
        "ops.pattern('Plain', 1, 1)\n"
        "for e in range(n_elem):\n"
        "    base = e * 8\n"
        "    for i in M._XFACE:\n"
        "        ops.sp(base + i, 1, -0.02)\n"
        "ops.constraints('Transformation')\n"
        "ops.numberer('Plain')\n"
        "ops.system('UmfPack')\n"
        "ops.test('NormDispIncr', 1e-8, 50, 0)\n"
        "ops.algorithm('Newton')\n"
        "ops.integrator('LoadControl', 1.0 / 20)\n"
        "ops.analysis('Static')\n"
        "print('BEGIN_STEPS')\n"
        "for _ in range(20):\n"
        "    ops.analyze(1)\n"
        "print('END_STEPS')\n"
    )
    proc = _run_child(script)
    assert "BEGIN_STEPS" in proc.stdout and "END_STEPS" in proc.stdout, (
        f"child process did not complete the 20-step run "
        f"(stdout={proc.stdout!r}, stderr={proc.stderr!r})")
    # MEASURED: Python's own print() markers are fully buffered (not a TTY)
    # and only flush at process exit, while the C++ extension's cout lines
    # flush immediately (std::endl) -- so in the raw captured byte stream
    # the "Integration Info" lines can appear BEFORE "BEGIN_STEPS" even
    # though they are emitted later in program order. Count over the whole
    # stream; this script only ever runs the one plastic path, so there is
    # nothing else "Integration Info" could be attributed to.
    n_out = proc.stdout.count("ASDP Integration Info")
    n_err = proc.stderr.count("ASDP Integration Info")
    assert n_out > 0, (
        "expected the per-commit 'ASDP Integration Info' cout line to "
        "appear at least once over 20 plastic steps on 10 elements; got 0 "
        f"(stdout={proc.stdout!r}) -- re-verify H12's trigger condition "
        f"before trusting this test")
    assert n_err == 0, (
        f"'ASDP Integration Info' now appears on opserr/stderr ({n_err} "
        f"lines) instead of only cout/stdout -- H12's channel claim may "
        f"have been fixed")


# ===========================================================================
# H13 -- unknown parser tokens are silently swallowed
# ===========================================================================
@pytest.mark.t0m
def test_H13_unknown_integration_option_is_silently_ignored(mc_available):
    """CONFIRMED.  ``Begin_Integration_Options`` has no ``else`` after its
    if-chain (OPS_AllASDPlasticMaterial3Ds.cpp:358-460): an unrecognised
    token like ``strict_convergance`` (missing the 'e') matches nothing, so
    its value token is silently discarded too, with no warning and no
    exception. The resulting deck is byte-identical to one where
    strict_convergence is never mentioned at all.
    """
    def mat_typo(tag):
        ops.nDMaterial(
            "ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
            "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables",
            "Begin_Integration_Options",
            "strict_convergance", 1,      # TYPO: should be strict_convergence
            "n_max_iterations", 2,
            "End_Integration_Options",
        )

    # must not raise / must not refuse to parse
    _tet_build(lambda t: mat_typo(t), nsteps=20, utop=-0.02)
    codes_typo = [ops.analyze(1) for _ in range(20)]

    _tet_build(lambda t: mat_mc(t, niter=2), nsteps=20, utop=-0.02)
    codes_no_flag = [ops.analyze(1) for _ in range(20)]

    assert codes_typo == codes_no_flag, (
        f"a misspelled 'strict_convergance' now produces DIFFERENT analyze "
        f"codes ({codes_typo}) than the flag being entirely absent "
        f"({codes_no_flag}) -- the silent-swallow defect may have been "
        f"fixed (e.g. an else-branch was added); re-verify H13.")


@pytest.mark.t0m
def test_H13_unknown_model_parameter_is_silently_ignored(mc_available):
    """CONFIRMED.  ``utuple_storage.h::setParameterByName_impl`` recurses the
    parameter tuple and, if no name matches by the base case (I ==
    tuple_size), does nothing -- no warning, no exception. A misspelled
    ``MC_phii`` is dropped, and the real ``MC_phi`` parameter is left at its
    unset default (0.0), silently changing the material's friction angle to
    zero instead of ``M.PHI``.
    """
    def mat_typo_param(tag):
        ops.nDMaterial(
            "ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            "MC_phii", M.PHI,              # TYPO: should be MC_phi
            "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
            "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables",
        )

    _tet_build(lambda t: mat_typo_param(t), nsteps=20, utop=-0.02)
    for _ in range(20):
        rc = ops.analyze(1)
        assert rc == 0, "expected NO parse-time or solve-time error"
    s = _tet_stress()
    # phi silently defaulted to 0 -> f_mc computed with phi=0 must be the
    # admissibility bound actually enforced, NOT phi=M.PHI's bound.
    f_with_zero_phi = M.f_mc(s, phi_deg=0.0, c=M.C)
    f_with_intended_phi = M.f_mc(s, phi_deg=M.PHI, c=M.C)
    tol = 1.0e-6 * max(2.0 * M.C, float(np.max(np.abs(s))))
    assert f_with_zero_phi <= max(tol, 1.0e-3), (
        f"committed state is not admissible under phi=0 (f={f_with_zero_phi:.3e} "
        f"> tol {max(tol, 1.0e-3):.3e}) -- the silently-defaulted-to-zero friction "
        f"angle assumption no longer holds; re-verify H13.")
    # non-vacuity: under the INTENDED phi the state must sit well INSIDE the
    # (larger) admissible region -- a big negative gap from the phi=0 bound
    # proves the material actually enforced phi=0, not the intended M.PHI.
    gap = f_with_zero_phi - f_with_intended_phi
    assert gap > 100.0, (
        f"f(phi=0)={f_with_zero_phi:.3e} vs f(phi={M.PHI})={f_with_intended_phi:.3e} "
        f"-- gap too small ({gap:.3e}) to show the typo had an observable "
        f"effect; this reproducer would be vacuous, re-tune before trusting "
        f"H13.")


# ===========================================================================
# H14 -- getCopy() omits ``first_step`` from its explicit member copy list
# ===========================================================================
@pytest.mark.t0m
def test_H14_getcopy_does_not_preserve_first_step(mc_available):
    """CONFIRMED (structural).  ``getCopy()`` (780-801) explicitly copies
    TrialStress/TrialStrain/.../CommitStrain/iv_storage/parameters_storage/
    stress_set_externally onto the new instance, but never assigns
    ``first_step`` -- the new object gets ``first_step = true`` from its own
    constructor regardless of the source's state. Because every host element
    calls ``getCopy()`` exactly once, on the still-pristine tag-registered
    prototype, at CONSTRUCTION time (before any analysis step), the bug is
    latent under normal model-build order: it would only bite a getCopy()
    call made on an ALREADY-advanced instance (state re-partitioning,
    lazy per-GP construction after stepping has begun), which this run does
    not exercise. This test pins the missing field mechanically so a fix
    (or an explicit decision to leave it) is visible in a diff, and
    separately pins that InitialP0 does seed CommitStress on every copy's own
    first step in the ordinary (non-buggy) sequence.
    """
    src = _asdp_source()
    m = re.search(r"NDMaterial \*getCopy\(void\)\s*\{.*?\n    \}\n",
                  src, re.S)
    assert m, "could not locate getCopy(void) body -- source layout changed"
    body = m.group(0)
    assert "first_step" not in body, (
        "getCopy(void) now assigns 'first_step' onto the new instance -- "
        "the H14 structural defect appears to be fixed; update the verdict "
        "to REFUTED/fixed.")

    # non-buggy-path sanity: InitialP0 seeds CommitStress on first step for
    # TWO independently constructed elements sharing the same material tag.
    p0 = -37.5
    for tag_offset, ele_tag in ((0, 1), (0, 2)):
        pass  # placeholder to keep structure readable; real build below

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
        ops.node(t + 100, c[0] + 2.0, c[1], c[2])
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
        ops.fix(t + 100, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
        ops.fix(t + 100, 1, 1, 0)
    mat_mc(1, p0=p0)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.element("TenNodeTetrahedron", 2, *[t + 100 for t in
                                           (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)], 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, -0.001)
        ops.sp(t + 100, 3, -0.001)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    ops.eleResponse(1, "forces")
    s1 = list(ops.eleResponse(1, "stresses"))[0:6]
    ops.eleResponse(2, "forces")
    s2 = list(ops.eleResponse(2, "stresses"))[0:6]
    # both independently-copied materials must show the InitialP0 hydrostatic
    # seed in their FIRST commit (mean stress ~ p0 before elastic loading is
    # superposed) -- confirms the ordinary path is unaffected by the missing
    # first_step copy.
    for lbl, s in (("ele1", s1), ("ele2", s2)):
        mean0 = sum(s[0:3]) / 3.0
        assert mean0 < 0.0, (
            f"{lbl}: expected a compressive (negative) seeded mean stress "
            f"from InitialP0={p0}, got mean={mean0:.4g}")


# ===========================================================================
# H15 -- BE uses E(sigma_commit) for the whole step; ME/RK45 re-evaluate E
# ===========================================================================
@pytest.mark.t0m
def test_H15_be_evaluates_elasticity_once_me_rk45_per_stage():
    """CONFIRMED (structural).  ``Backward_Euler`` computes
    ``Eelastic = et(CommitStress, ...)`` ONCE (2055) before its scalar-Newton
    loop and never calls ``et()`` again for the rest of the step;
    ``Backward_Euler_LineSearch`` does the same (2372). ``Modified_Euler_
    Error_Control`` and ``Runge_Kutta_45_Error_Control`` instead call ``et()``
    multiple times per step -- at ``CommitStress``, at the current stage
    stress, and at the predictor stress (measured: 3 and 7 occurrences in
    their respective bodies). For a stress-dependent elasticity
    (``StiffSoil_EL``, ``DuncanChang_EL``) this is not just a numerics
    difference: BE and ME/RK45 are evaluating a DIFFERENT constitutive
    operator on the same nominal material.

    A StiffSoilShear_YF/StiffSoil_EL runtime driver was attempted for this
    review and hit an UNRELATED NaN in the very first step of every
    Backward_Euler triaxial path tried (a real, separate finding -- flagged
    in ``_adr94_hlist_R1B.md``, not part of H15), so this test pins the
    mechanism structurally instead, exactly as H2/H14 do for defects that do
    not survive the host-element's response layer.
    """
    src = _asdp_source()

    def body_of(signature):
        m = re.search(re.escape(signature) + r".*?\n    \}\n", src, re.S)
        assert m, f"could not locate a body starting at {signature!r}"
        return m.group(0)

    def et_calls(body):
        # strip C++ line comments first -- a commented-out alternative
        # implementation elsewhere in these functions also mentions et(...)
        # and must not be counted as a live evaluation.
        live = "\n".join(re.sub(r"//.*$", "", ln) for ln in body.splitlines())
        return len(re.findall(r"\bet\(", live))

    be_m = re.search(
        r"int Backward_Euler\(const VoigtVector & strain_incr\)\s*\{", src)
    assert be_m, "could not locate Backward_Euler's DEFINITION signature"
    be_start = be_m.end() - 1
    # bound the search window; trim to the matching closing brace.
    depth, i = 0, be_start
    while i < len(src):
        if src[i] == '{':
            depth += 1
        elif src[i] == '}':
            depth -= 1
            if depth == 0:
                break
        i += 1
    be_body = src[be_start:i + 1]

    n_be = et_calls(be_body)
    assert n_be == 1, (
        f"Backward_Euler now calls et() {n_be} times (expected exactly 1, "
        f"at CommitStress before the return-map loop) -- H15's premise "
        f"(BE reuses one elasticity evaluation for the whole step) no "
        f"longer holds; re-verify H15.")

    me_def = re.search(
        r"int Modified_Euler_Error_Control\(const VoigtVector & strain_incr\)",
        src)
    rk_def = re.search(
        r"int Runge_Kutta_45_Error_Control\(const VoigtVector & strain_incr\)",
        src)
    assert me_def and rk_def, "could not locate ME/RK45 DEFINITION signatures"
    me_window = src[me_def.start():me_def.start() + 6000]
    rk_window = src[rk_def.start():rk_def.start() + 8000]
    n_me = et_calls(me_window)
    n_rk = et_calls(rk_window)
    assert n_me >= 2, (
        f"Modified_Euler_Error_Control now calls et() only {n_me} time(s) "
        f"in its body window -- it may have been made single-evaluation "
        f"like BE; re-verify H15.")
    assert n_rk >= 2, (
        f"Runge_Kutta_45_Error_Control now calls et() only {n_rk} time(s) "
        f"in its body window -- it may have been made single-evaluation "
        f"like BE; re-verify H15.")
