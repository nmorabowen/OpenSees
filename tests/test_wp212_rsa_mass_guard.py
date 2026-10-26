"""WP-212: responseSpectrumAnalysis mass-shortfall guard and -combine warning.

Feature A. ASCE 7-22 §12.9.1.1 asks the modes of a response spectrum analysis
to carry 100% of the mass (90% per horizontal direction as the exception);
NEC-SE-DS §6.2.2e asks for 90%. Stock RSA never looked: a run with 2 of 10
modes, or with no participating mass at all (the diaphragm-mass bug fixed in
WP-211), returned displacements without a word. Every RSA path now reads the
cumulative effective-mass ratio of the excited direction after
getModalProperties, prints it once per call, warns below 90%, and errors when
the direction has no free mass or no participating mass. `-requireMass r`
turns the warning into an error below r.

Feature B. `-combine` combines nodal DISPLACEMENTS and commits them; a force
read afterwards is derived from that field and is not the CQC of the force.
The call now says so.

Oracles (pure numpy, no reference outputs):
  * the cumulative ratio is sum_k (phi_k' M 1)^2 / (phi_k' M phi_k) / sum m,
    from scipy-free numpy eigh on the same K, M;
  * "nothing changes when the mass is complete" is pinned against values
    recorded with the pre-WP-212 engine (build f2a3a4b6d, see BASELINE_*),
    and the WP-207 -increment path against the stock path;
  * the column shear after -combine CQC is compared with the CQC of the
    per-mode column shears, using the Der Kiureghian rho of
    LadrunoModalCombination.h.
"""
import os
import subprocess

import numpy as np
import pytest

from _testbed import ops
from _testbed.subprocess_run import pinned_child

pytestmark = [pytest.mark.zone_a]

TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]
SPEC = ['-Tn', *TN, '-Sa', *SA]

REPORT = "cumulative effective mass"
LOW_WARN = "below 90%"
COMBINE_WARN = "-combine combines nodal displacements only"
INCR_DISP = 4  # NodeResponseType IncrDisp = trial - committed

# The report line is hidden whenever LADRUNO_OPENSEES_QUIET is set (even empty).
# The engine reads its own copy of the environment, so a monkeypatch cannot
# undo it: in-process cases expect the line only when the var is unset, and
# test_report_line_in_clean_child checks the ratio with the var removed.
QUIET = "LADRUNO_OPENSEES_QUIET" in os.environ
N_REPORT = 0 if QUIET else 1


# --------------------------------------------------------------------------
# 10-DOF shear building. A heavy, stiff podium (DOF 1) keeps ~70% of the mass
# out of the first two modes; the uniform variant puts 93.9% in them.
# --------------------------------------------------------------------------
PODIUM_M = [20.0] + [1.0] * 9
PODIUM_K = [4000.0] + [100.0] * 9
UNIFORM_M = [1.0] * 10
UNIFORM_K = [100.0] * 10


def _shear_building(masses, stiff, nmodes):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(0, 0.0)
    ops.fix(0, 1)
    for i, (mm, kk) in enumerate(zip(masses, stiff), start=1):
        ops.node(i, float(i))
        ops.mass(i, mm)
        ops.uniaxialMaterial('Elastic', i, kk)
        ops.element('Truss', i, i - 1, i, 1.0, i)
    ops.eigen('-fullGenLapack', nmodes)
    ops.modalProperties()


def _numpy_ratio(masses, stiff, nmodes):
    n = len(masses)
    K = np.zeros((n, n))
    for i, kk in enumerate(stiff):
        K[i, i] += kk
        if i > 0:
            K[i - 1, i - 1] += kk
            K[i - 1, i] -= kk
            K[i, i - 1] -= kk
    m = np.asarray(masses)
    # M^-1/2 K M^-1/2 is symmetric; its eigenvectors map back to phi = M^-1/2 v
    s = 1.0 / np.sqrt(m)
    lam, v = np.linalg.eigh(s[:, None] * K * s[None, :])
    phi = s[:, None] * v
    L = phi.T @ m
    meff = L ** 2 / np.einsum('ik,i,ik->k', phi, m, phi)
    return float(np.sum(meff[:nmodes]) / m.sum())


def _reported_ratio(err):
    """The percentage printed on the report line."""
    for line in err.splitlines():
        if REPORT in line:
            return float(line.split('=')[1].split('%')[0]) / 100.0
    raise AssertionError("no mass report line in stderr: %r" % err[-600:])


def test_two_of_ten_modes_warns_with_right_ratio(capfd):
    _shear_building(PODIUM_M, PODIUM_K, 2)
    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1)
    err = capfd.readouterr().err
    ref = _numpy_ratio(PODIUM_M, PODIUM_K, 2)
    assert ref < 0.5
    if not QUIET:
        assert _reported_ratio(err) == pytest.approx(ref, abs=1e-4)
    assert LOW_WARN in err
    assert err.count(REPORT) == N_REPORT


_CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import test_wp212_rsa_mass_guard as t
t._shear_building(t.PODIUM_M, t.PODIUM_K, 2)
t.ops.responseSpectrumAnalysis(1, *t.SPEC, '-mode', 1)
"""


def test_report_line_in_clean_child():
    # the report line itself, in a fresh engine started without the quiet var
    argv, env = pinned_child(['-u', '-c', _CHILD, os.path.dirname(__file__)])
    env.pop("LADRUNO_OPENSEES_QUIET", None)
    p = subprocess.run(argv, env=env, stdin=subprocess.DEVNULL, capture_output=True,
                       encoding='utf-8', errors='replace', timeout=300)
    assert p.returncode == 0, p.stdout + p.stderr
    assert p.stderr.count(REPORT) == 1, p.stderr[-600:]
    assert _reported_ratio(p.stderr) == pytest.approx(
        _numpy_ratio(PODIUM_M, PODIUM_K, 2), abs=1e-4)
    assert LOW_WARN in p.stderr


@pytest.mark.parametrize("path", ["mode", "all", "combine", "increment"])
def test_every_path_reports_once(path, capfd):
    _shear_building(PODIUM_M, PODIUM_K, 2)
    extra = {"mode": ['-mode', 2], "all": [],
             "combine": ['-combine', 'SRSS'], "increment": ['-increment']}[path]
    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, *extra)
    err = capfd.readouterr().err
    assert err.count(REPORT) == N_REPORT, err
    assert err.count(LOW_WARN) == 1, err


def test_require_mass_errors_below_threshold(capfd):
    _shear_building(PODIUM_M, PODIUM_K, 2)
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *SPEC, '-requireMass', 0.9)
    assert "-requireMass" in capfd.readouterr().err


def test_require_mass_passes_above_threshold(capfd):
    # uniform building: 2 modes carry 93.9% -> no warning, 0.9 passes, 0.95 fails
    _shear_building(UNIFORM_M, UNIFORM_K, 2)
    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1, '-requireMass', 0.9)
    err = capfd.readouterr().err
    if not QUIET:
        assert _reported_ratio(err) == pytest.approx(
            _numpy_ratio(UNIFORM_M, UNIFORM_K, 2), abs=1e-4)
    assert LOW_WARN not in err
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1, '-requireMass', 0.95)


@pytest.mark.parametrize("bad", [0.0, -0.5, 1.5])
def test_require_mass_range(bad):
    _shear_building(UNIFORM_M, UNIFORM_K, 2)
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *SPEC, '-requireMass', bad)


# --------------------------------------------------------------------------
# Zero participation. One 2-DOF node on two springs; mass only where stated.
# --------------------------------------------------------------------------
def _two_dof_node(mx, my, nmodes):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.mass(2, mx, my)
    ops.uniaxialMaterial('Elastic', 1, 400.0)   # X: stiff
    ops.uniaxialMaterial('Elastic', 2, 100.0)   # Y: soft
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 2)
    ops.eigen('-fullGenLapack', nmodes)
    ops.modalProperties()


@pytest.mark.parametrize("path", [['-mode', 1], [], ['-combine', 'SRSS'], ['-increment']])
def test_no_free_mass_in_direction_errors(path, capfd):
    # mass in Y only; spectrum in X -> free mass in X is zero
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.mass(2, 0.0, 2.0)
    ops.uniaxialMaterial('Elastic', 2, 100.0)
    ops.element('zeroLength', 1, 1, 2, '-mat', 2, '-dir', 2)
    ops.fix(2, 1, 0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *SPEC, *path)
    err = capfd.readouterr().err
    assert "no free mass in direction 1" in err, err
    assert ops.nodeDisp(2, 1) == 0.0 and ops.nodeDisp(2, 2) == 0.0


def test_zero_participation_errors(capfd):
    # mass in X and Y, but only the (soft) Y mode is extracted: every Gamma_x = 0
    _two_dof_node(2.0, 2.0, 1)
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1)
    err = capfd.readouterr().err
    assert "no participating mass in direction 1" in err, err
    # the same modes are fine for the direction they excite
    ops.responseSpectrumAnalysis(2, *SPEC, '-mode', 1)
    assert ops.nodeDisp(2, 2) != 0.0


# --------------------------------------------------------------------------
# 2-story, 1-bay 2-D frame; mass in X, Y and rotation at every floor node so
# fullGenLapack extracts all 12 modes (complete mass, ratio = 100%).
# --------------------------------------------------------------------------
FRAME_E, FRAME_A, FRAME_I = 2.0e7, 0.09, 6.75e-4
FLOOR_NODES = (3, 4, 5, 6)
COLUMN = 1   # ground-story column 1-3


def _frame(nmodes=12):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    xy = {1: (0.0, 0.0), 2: (5.0, 0.0), 3: (0.0, 3.0), 4: (5.0, 3.0),
          5: (0.0, 6.0), 6: (5.0, 6.0)}
    for n, (x, y) in xy.items():
        ops.node(n, x, y)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    for n, mm in zip(FLOOR_NODES, (3.0, 3.0, 2.0, 2.0)):
        ops.mass(n, mm, mm, 0.01)
    ops.geomTransf('Linear', 1)
    conn = [(1, 3), (2, 4), (3, 5), (4, 6), (3, 4), (5, 6)]   # 4 columns, 2 beams
    for e, (i, j) in enumerate(conn, start=1):
        ops.element('elasticBeamColumn', e, i, j, FRAME_A, FRAME_E, FRAME_I, 1)
    ops.eigen('-fullGenLapack', nmodes)
    ops.modalProperties()


def _frame_state():
    disp = [ops.nodeDisp(n, d) for n in FLOOR_NODES for d in (1, 2, 3)]
    ops.reactions()
    reac = [ops.nodeReaction(n, d) for n in (1, 2) for d in (1, 2, 3)]
    return disp, reac


# Recorded with the pre-WP-212 engine (Ladruno build f2a3a4b6d, which predates
# WP-202/207 too; at the default scale of 1 those change nothing): _frame(),
# dir 1, SPEC.
# BASELINE_MODE1: -mode 1; BASELINE_ALL: no -mode (mode 12 left committed);
# BASELINE_CQC: -combine CQC -damp 0.05.
BASELINE_MODE1 = ([0.0017140390370379052, 9.731344858201036e-06, -0.0005592139292682678, 0.0017140390370379045, -9.731344858201205e-06, -0.0005592139292682702, 0.003583583742366932, 1.3444461067232632e-05, -0.0003474145967942189, 0.0035835837423669322, -1.3444461067232344e-05, -0.0003474145967942203],
                  [-5.25130885881302, -5.838806914920622, 10.393425969926735, -5.251308858812995, 5.838806914920723, 10.393425969926708])
BASELINE_ALL = ([-2.186677022974927e-12, 3.993515543765721e-12, 9.811596716583045e-09, -2.1866770229748127e-12, -3.993515543765716e-12, 9.81159671658301e-09, 1.1171085190869637e-11, 2.2633533712387762e-12, 4.067336198283303e-09, 1.1171085190868159e-11, -2.2633533712387754e-12, 4.0673361982832625e-09],
                [-8.829125038710955e-05, -2.396109326259432e-06, 8.828469035604064e-05, -8.829125038710922e-05, 2.3961093262594296e-06, 8.828469035604032e-05])
BASELINE_CQC = ([0.00171712266058645, 9.74165311305712e-06, 0.0005592051024975392, 0.0017171226605864492, 9.741653113057289e-06, 0.0005592051024975415, 0.0035838066141538126, 1.3473175472287763e-05, 0.00034959641969704316, 0.003583806614153813, 1.3473175472287479e-05, 0.0003495964196970445],
                [-15.335581885996554, -5.844991867834271, 20.486949867755904, -15.335581885996564, -5.8449918678343735, 20.486949867755914])


def _assert_same(state, base):
    disp, reac = state
    bdisp, breac = base
    assert disp == pytest.approx(bdisp, rel=1e-12, abs=1e-15)
    assert reac == pytest.approx(breac, rel=1e-12, abs=1e-12)


def test_complete_mass_reports_100_percent_no_warning(capfd):
    _frame()
    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1, '-requireMass', 1.0 - 1e-9)
    err = capfd.readouterr().err
    if not QUIET:
        assert _reported_ratio(err) == pytest.approx(1.0, abs=1e-4)
    assert LOW_WARN not in err


def test_complete_mass_results_unchanged_per_mode():
    _frame()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1)
    _assert_same(_frame_state(), BASELINE_MODE1)


def test_complete_mass_results_unchanged_all_modes():
    _frame()
    ops.responseSpectrumAnalysis(1, *SPEC)
    _assert_same(_frame_state(), BASELINE_ALL)


def test_complete_mass_results_unchanged_combine():
    _frame()
    ops.responseSpectrumAnalysis(1, *SPEC, '-combine', 'CQC', '-damp', 0.05)
    _assert_same(_frame_state(), BASELINE_CQC)


def test_complete_mass_results_unchanged_increment():
    # about a zero committed state, -increment must reproduce the stock mode
    # (nodeDisp is the committed state; the mode is the trial increment)
    _frame()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1, '-increment')
    disp = [ops.nodeResponse(n, d, INCR_DISP) for n in FLOOR_NODES for d in (1, 2, 3)]
    assert disp == pytest.approx(BASELINE_MODE1[0], rel=1e-12, abs=1e-15)


def test_require_mass_flag_does_not_change_results():
    _frame()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1,
                                 '-requireMass', 0.9)
    _assert_same(_frame_state(), BASELINE_MODE1)


# --------------------------------------------------------------------------
# Feature B: -combine is a displacement combination, not a force combination.
# --------------------------------------------------------------------------
def _rho_cqc(wi, wj, xi):
    # LadrunoModalCombination.h, equal damping
    r = wj / wi
    return (8.0 * xi * xi * (1.0 + r) * r ** 1.5 /
            ((1.0 - r * r) ** 2 + 4.0 * xi * xi * r * (1.0 + r) ** 2))


def _column_shear():
    return ops.eleForce(COLUMN)[0]   # global X force at node i of column 1-3


def test_combine_force_is_not_cqc_of_forces(capfd):
    xi = 0.05
    _frame()
    w = np.sqrt(np.asarray(ops.eigen('-fullGenLapack', 12)))
    ops.modalProperties()
    per_mode = []
    for k in range(1, 13):
        ops.responseSpectrumAnalysis(1, *SPEC, '-mode', k)
        per_mode.append(_column_shear())
    v = np.asarray(per_mode)
    rho = np.array([[_rho_cqc(a, b, xi) for b in w] for a in w])
    v_cqc = float(np.sqrt(v @ rho @ v))

    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, '-combine', 'CQC', '-damp', xi)
    err = capfd.readouterr().err
    v_comb = abs(_column_shear())

    assert COMBINE_WARN in err, err
    assert err.count(COMBINE_WARN) == 1
    assert "-mode k" in err
    assert v_cqc > 0.0
    assert abs(v_comb - v_cqc) / v_cqc > 0.05, (v_comb, v_cqc)


def test_per_mode_paths_do_not_print_combine_warning(capfd):
    _frame()
    capfd.readouterr()
    ops.responseSpectrumAnalysis(1, *SPEC, '-mode', 1)
    ops.responseSpectrumAnalysis(1, *SPEC)
    assert COMBINE_WARN not in capfd.readouterr().err
