"""LadrunoBrick — NONLINEAR-MATERIAL validation (J2 plasticity).

The elastic gates (test_ladrunoBrick_element.py / _bending.py) prove the element
is correct for linear-elastic small strain. They say nothing about path-dependent
behaviour: the material trial-strain history, the per-step commit cycle, and — for
`ssp` — the stabilization built from the INITIAL tangent while the membrane block
uses the CURRENT tangent. This file closes that gap the same rigorous way the
elastic gates do: **element-vs-reference regression under a real J2 material.**

  * `std`  must track upstream `stdBrick`   step-for-step
  * `bbar` must track upstream `bbarBrick`
  * `ssp`  must track `SSPbrick` (the production stabilized single-point hex we ported)

These need NO analytic oracle — two elements driving the same J2 material through
the same boundary-value problem must produce the same displacement history. The
`ssp`↔`SSPbrick` case specifically verifies that our frozen-initial-tangent
stabilization reproduces the proven element in the inelastic regime, not just the
elastic one. Plus: J2 plastic flow is isochoric, so developed plasticity drives
the same volumetric over-constraint that `bbar`/`ssp` exist to relieve — `std`
must lock relative to them.

Material: J2Plasticity (q(xi) = sigma0 + H*xi with sigInf=sigma0, delta=0 -> linear
isotropic hardening).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# --- a deliberately distorted (positive-Jacobian) unit-ish hex (Brick order) ---
_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.00),
    3: (1.10, 1.00, 0.10),
    4: (0.05, 0.95, 0.00),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_BASE = [1, 2, 3, 4]    # z- face, fully fixed
_TOP = [5, 6, 7, 8]     # z+ face, loaded

# steel-ish, consistent N/mm-style units
_E, _NU = 200000.0, 0.3
_K = _E / (3.0 * (1.0 - 2.0 * _NU))
_G = _E / (2.0 * (1.0 + _NU))
_SIG0, _H = 250.0, 2000.0


def _solve_hex_J2(ele_name, extra_args, ftot, nsteps=25, direction=(0.0, 0.0, 1.0)):
    """Ramp a total nodal force `ftot` (split over the top face, along `direction`)
    over `nsteps` Newton load steps on a single J2 hex, base fully fixed. Returns
    the final 24-vector of nodal displacements. `ftot` is sized to push the element
    well past yield so the comparison exercises plastic flow + the commit cycle."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    # J2Plasticity: tag K G sig0 sigInf delta H  (sigInf=sig0, delta=0 -> linear H)
    ops.nDMaterial("J2Plasticity", 1, _K, _G, _SIG0, _SIG0, 0.0, _H)

    for n in _BASE:
        ops.fix(n, 1, 1, 1)
    ops.element(ele_name, 1, *_CONN, 1, *extra_args)

    ops.timeSeries("Linear", 1)         # factor ramps 0 -> 1 over the steps
    ops.pattern("Plain", 1, 1)
    dx, dy, dz = direction
    fnode = ftot / len(_TOP)
    for n in _TOP:
        ops.load(n, fnode * dx, fnode * dy, fnode * dz)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"{ele_name} {extra_args}: J2 analysis did not converge"

    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    return disps


def _assert_close(a, b, label, rtol=1e-5, atol=1e-7):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x-y):.3e} > {tol:.3e})"
        )


def _mag(u):
    return math.sqrt(sum(v * v for v in u))


# --------------------------------------------------------------------------
# element-vs-reference regression UNDER PLASTICITY (the rigorous, oracle-free gate)
# --------------------------------------------------------------------------
def test_std_matches_stdBrick_J2():
    """-formulation std tracks upstream stdBrick step-for-step under J2."""
    ours = _solve_hex_J2("LadrunoBrick", ["-formulation", "std"], ftot=450.0)
    ref = _solve_hex_J2("stdBrick", [], ftot=450.0)
    _assert_close(ours, ref, "std/J2 disp")


def test_bbar_matches_bbarBrick_J2():
    """-formulation bbar tracks upstream bbarBrick under J2 (the path-dependent
    analogue of the elastic ~1e-9 regression)."""
    ours = _solve_hex_J2("LadrunoBrick", ["-formulation", "bbar"], ftot=450.0)
    ref = _solve_hex_J2("bbarBrick", [], ftot=450.0)
    _assert_close(ours, ref, "bbar/J2 disp")


def test_ssp_matches_sspbrick_J2():
    """THE inelastic gate for ssp: our SSPbrick port must reproduce SSPbrick under
    a path-dependent material, not just elastically. Both freeze Kstab at the
    initial tangent and drive the membrane with the current J2 tangent, so the
    displacement history must agree closely."""
    ours = _solve_hex_J2("LadrunoBrick", ["-formulation", "ssp"], ftot=450.0)
    ref = _solve_hex_J2("SSPbrick", [], ftot=450.0)
    _assert_close(ours, ref, "ssp/J2 disp", rtol=2e-5, atol=1e-7)


# --------------------------------------------------------------------------
# plastic incompressibility: developed J2 flow is isochoric -> std locks
# --------------------------------------------------------------------------
def test_plastic_incompressibility_std_locks():
    """J2 plastic flow is volume-preserving, so once the element yields it sees the
    same incompressibility over-constraint as elastic nu->0.5. Full integration
    (std) locks; bbar and ssp (mean-dilatation / condensed stabilization) relieve
    it. Push well past yield and compare displacement magnitudes."""
    f = 600.0
    std = _solve_hex_J2("LadrunoBrick", ["-formulation", "std"], ftot=f)
    bbar = _solve_hex_J2("LadrunoBrick", ["-formulation", "bbar"], ftot=f)
    ssp = _solve_hex_J2("LadrunoBrick", ["-formulation", "ssp"], ftot=f)
    std_m, bbar_m, ssp_m = _mag(std), _mag(bbar), _mag(ssp)
    assert bbar_m > std_m, f"bbar |u|={bbar_m:.4e} should exceed std |u|={std_m:.4e} (std plastic-locks)"
    assert ssp_m > std_m, f"ssp |u|={ssp_m:.4e} should exceed std |u|={std_m:.4e} (std plastic-locks)"


def test_j2_actually_yields():
    """Guard: confirm the load level genuinely drives the element plastic (so the
    regressions above are exercising the inelastic path, not a glorified elastic
    run). The committed centroid stress should exceed the initial yield sigma0."""
    _solve_hex_J2("LadrunoBrick", ["-formulation", "std"], ftot=450.0)
    ops.eleResponse(1, "forces")                       # ensure strain/stress is set
    sig = list(ops.eleResponse(1, "stresses"))[0:6]    # centroid Gauss point
    svm = math.sqrt(0.5 * ((sig[0] - sig[1]) ** 2 + (sig[1] - sig[2]) ** 2
                           + (sig[2] - sig[0]) ** 2)
                    + 3.0 * (sig[3] ** 2 + sig[4] ** 2 + sig[5] ** 2))
    assert svm > _SIG0, f"von Mises {svm:.1f} should exceed yield {_SIG0} (element must yield)"
