"""ExplicitBatheLNVD (Noh-Bathe + FLAC local-non-viscous damping, classTag
33002) — Zone-A dynamics battery.

Zone-A port of the LNVD check (test 6) from ``Ladruno_scripts/_verify_explicit.py``
(jaabell's tag, band-remediated 63 -> 33002). LNVD adds FLAC mass-proportional
local damping (alpha, default 0.8) on top of the explicit Noh-Bathe step, so its
defining use is *dynamic relaxation*: damping the transient out drives the
solution to the static equilibrium of the applied load.

The single proven source check (``test_lnvd_static``) is parametrized here over
load magnitude. The bar is linear-elastic, so the relaxation dynamics (and step
count to settle) are load-independent — only the amplitude scales — and each
case self-validates against (a) its own LoadControl/Newton static solve and
(b) the closed-form tip displacement PL/EA.

The integrator signature is ``integrator('ExplicitBatheLNVD', p, alpha)`` with
p=0.54 (sub-step) and alpha=0.8 (classic FLAC default) — confirmed against
SRC/analysis/integrator/ExplicitBatheLNVD.cpp.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

LNVD = "ExplicitBatheLNVD"
P_SUB = 0.54   # Noh-Bathe sub-step parameter
ALPHA = 0.8    # FLAC local-damping coefficient (classic default)


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


@pytest.mark.parametrize("P_load", [2.0, 1.0, 3.0])
def test_lnvd_dynamic_relaxation_to_static(P_load):
    N = 20
    L = 10.0
    E = 100.0
    A = 1.0
    rho = 1.0

    # --- static reference (LoadControl + Newton) ---
    build_bar(N, L, E, A, rho)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P_load, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandSPD")
    ops.test("NormDispIncr", 1e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.analyze(1)
    u_static = np.array([ops.nodeDisp(i + 1, 1) for i in range(N + 1)])

    # --- LNVD dynamic relaxation ---
    build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P_load, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(LNVD, P_SUB, ALPHA)
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    dt = 0.5 * (L / N) / c
    for _ in range(6000):
        ops.analyze(1, dt)
    u_lnvd = np.array([ops.nodeDisp(i + 1, 1) for i in range(N + 1)])

    tip_exact = P_load * L / (E * A)
    err_field = np.abs(u_lnvd - u_static).max() / abs(u_static[-1])
    err_tip = abs(u_lnvd[-1] - tip_exact) / tip_exact
    assert err_field < 0.02, (
        "LNVD field vs static mismatch %.2f%% (expect < 2%%)" % (100 * err_field)
    )
    assert err_tip < 0.02, (
        "LNVD tip %.4f vs PL/EA %.4f (err %.2f%%, expect < 2%%)"
        % (u_lnvd[-1], tip_exact, 100 * err_tip)
    )
