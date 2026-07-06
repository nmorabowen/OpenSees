"""EnergyBalanceRecorder (per-region energy sidecar, RECORDER classTag 33000) —
Zone-A T1 check.

Zone-A port of the EnergyBalance check (test 9) from
``Ladruno_scripts/_verify_explicit.py`` (band-remediated 26 -> 33000). The
recorder writes a per-step energy ledger; with ``-time`` the columns are

    time  KE  IE  DW  ULW  RES  ERR%

(confirmed against SRC/recorder/EnergyBalanceRecorder.cpp, column header
{"KE","IE","DW","ULW","RES","ERR%"}, time prepended when ``-time`` is given).

The primary check is the getMass()/getDamp() shared-``theMatrix`` aliasing-fix
regression: with ELEMENT mass (Truss ``-rho``) the recorder must read kinetic
energy from the element MASS matrix, not the (zeroed) damping matrix. Before the
fix the element KE came out ~0; afterwards it tracks 0.5 * M * v^2. This recorder
already runs green on Zone-A CI (CDL-10 exercises the same KE/IE columns).
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


# ==========================================================================
# element-mass KE: the getMass()/getDamp() aliasing-fix regression (test 9)
# ==========================================================================
def test_energybalance_element_mass_ke(tmp_path):
    N = 20
    L = 10.0
    E = 100.0
    A = 1.0
    rho = 1.0
    v0 = 2.0
    efile = str(tmp_path / "bar_energy.txt")
    build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):                 # initial x-velocity on every free node
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    ops.analyze(2, 1e-4)
    ops.wipe()                                # flush/close the recorder

    d = np.loadtxt(efile)
    d = np.atleast_2d(d)
    ke0 = d[0, 1]                             # col 0 = time, col 1 = KE
    M_total = rho * A * L                     # total bar mass (element-based)
    ke_expected = 0.5 * M_total * v0 ** 2     # ~ (fixed end excludes ~rho*A*le/2)
    assert ke0 > 0.5 * ke_expected, (
        "element-mass KE=%.3f vs ~0.5*M*v^2=%.3f — before the aliasing fix this "
        "was ~0 (read from the zeroed damping matrix)" % (ke0, ke_expected)
    )


# ==========================================================================
# column layout: with -time, exactly 7 columns (time + KE IE DW ULW RES ERR%)
# ==========================================================================
def _run_lnvd_bar(tmp_path, fname, v2, alpha):
    """Free-vibration bar with ExplicitBathe [-lnvd alpha]; returns the data."""
    N = 10
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    efile = str(tmp_path / fname)
    le = build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, 1.0, "-commit")
    # NOTE (ADR-69): integrator BEFORE recorder — the -v2 channel columns are
    # fixed at the recorder's first record from the registry's declared set.
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if alpha > 0.0:
        ops.integrator("ExplicitBathe", 0.54, "-lnvd", alpha)
    else:
        ops.integrator("ExplicitBathe", 0.54)
    if v2:
        ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2")
    else:
        ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    ops.analyze(400, 0.2 * le / c)
    ops.wipe()
    d = np.loadtxt(efile)
    return np.atleast_2d(d)


# ==========================================================================
# ADR-69 v2: no-producer layout — time + KE_ele KE_nod IE DW_ele DW_nod ULW
# RES ERR% (no channel columns), element/nodal KE split sane.
# NOTE: must run BEFORE any -lnvd test in this module — channel declarations
# are process-sticky (a later -v2 run would gain a benign zero E_lnvd column).
# ==========================================================================
def test_energybalance_v2_layout_no_producers(tmp_path):
    d = _run_lnvd_bar(tmp_path, "v2_plain.txt", v2=True, alpha=0.0)
    assert d.shape[1] == 9, (
        "expected 9 columns (time + 6 split + RES ERR%%) with no producers; got %d"
        % d.shape[1]
    )
    ke_ele, ke_nod = d[0, 1], d[0, 2]
    assert ke_ele > 0.0 and abs(ke_nod) < 1e-12, (
        "element-mass bar: KE must sit in KE_ele (got ele=%g nod=%g)" % (ke_ele, ke_nod)
    )


# ==========================================================================
# ADR-69 v2: E_lnvd channel closes what v1 leaks into RES. Free vibration +
# FLAC local damping: legacy RES drifts by exactly the (unaccounted) LNVD
# dissipation; v2 publishes it (E_lnvd column) and RES stays put.
# ==========================================================================
def test_energybalance_v2_lnvd_closure(tmp_path):
    alpha = 0.4
    d_legacy = _run_lnvd_bar(tmp_path, "lnvd_legacy.txt", v2=False, alpha=alpha)
    d_v2 = _run_lnvd_bar(tmp_path, "lnvd_v2.txt", v2=True, alpha=alpha)

    # legacy: time KE IE DW ULW RES ERR% -> RES col 5
    res_l = d_legacy[:, 5]
    drift_legacy = res_l[-1] - res_l[0]

    # v2 with LNVD declared: time KE_ele KE_nod IE DW_ele DW_nod ULW E_lnvd RES ERR%
    assert d_v2.shape[1] == 10, (
        "expected 10 columns (time + 6 split + E_lnvd + RES ERR%%); got %d"
        % d_v2.shape[1]
    )
    e_lnvd = d_v2[:, 7]
    res_v2 = d_v2[:, 8]
    drift_v2 = res_v2[-1] - res_v2[0]

    assert e_lnvd[-1] > 0.0, "LNVD dissipation must accumulate positive work"
    assert np.all(np.diff(e_lnvd) > -1e-14), "E_lnvd must be monotone"
    # the legacy RES drift IS the unaccounted LNVD dissipation
    assert drift_legacy > 0.0
    assert 0.7 * drift_legacy < e_lnvd[-1] < 1.3 * drift_legacy, (
        "published E_lnvd=%g should match the legacy RES drift=%g"
        % (e_lnvd[-1], drift_legacy)
    )
    # closure: v2 RES drift shrinks by at least 5x
    assert abs(drift_v2) < 0.2 * abs(drift_legacy), (
        "v2 RES drift %g should be <<%g (legacy)" % (drift_v2, drift_legacy)
    )


# ==========================================================================
# column layout: with -time, exactly 7 columns (time + KE IE DW ULW RES ERR%)
# ==========================================================================
def test_energybalance_time_column_layout(tmp_path):
    N = 5
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    efile = str(tmp_path / "cols.txt")
    le = build_bar(N, L, E, A, rho)
    ops.setNodeVel(N + 1, 1, 1.0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    ops.analyze(5, 0.5 * le / c)
    ops.wipe()

    rows = []
    with open(efile) as fh:
        for line in fh:
            vals = line.split()
            if vals:
                rows.append(vals)
    assert rows, "EnergyBalance recorder produced no output"
    # time + model block (KE IE DW ULW RES ERR%) = 1 + 6 = 7 columns
    assert all(len(r) == 7 for r in rows), (
        "expected 7 columns per row (time + KE IE DW ULW RES ERR%%); got widths %s"
        % sorted({len(r) for r in rows})
    )
