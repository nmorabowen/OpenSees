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
