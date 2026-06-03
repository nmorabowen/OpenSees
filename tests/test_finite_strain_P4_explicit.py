"""Finite-strain trifecta — **P4 explicit dynamics** (Taylor bar + energy + dt_cr).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §9
phase **P4** (gate: *Taylor bar; dt_cr caveat documented*). This is the explicit
counterpart of the implicit P1–P3 work: `LadrunoBrick -geom finite` (F-bar) +
`LogStrain(LadrunoJ2)` driven by the fork's **`CentralDifferenceLadruno`**
leap-frog integrator with a **lumped** mass matrix and a `Diagonal` system.

  * **C2 — Taylor-bar impact** (Taylor 1948; Kamoulakos benchmark): a copper
    cylinder at 227 m/s into a frictionless rigid wall mushrooms plastically. The
    final length and the mushroom (footprint) radius match the literature, and
    the deformation localizes at the impact end (the free end is undeformed).
  * **Energy balance**: through the impact the kinetic energy is converted to
    internal (elastic + plastic) energy — `|ΔKE| = |IE|` to engineering tolerance,
    with most of KE₀ absorbed plastically and a small elastic rebound.
  * **dt_cr caveat**: `criticalTimeStep()` is computed from the **reference**
    configuration — it does NOT shrink as elements compress (review GEOM-2), so an
    explicit `-geom finite` run must carry a safety factor (< 1) to stay stable
    through strong compression. The run uses 0.3·dt_cr.

Model — copper, N–mm–tonne–s units: E=117 GPa, ν=0.35, ρ=8.93e-9, σ_y=400 MPa,
linear hardening H=100 MPa; L₀=32.4 mm, r₀=3.2 mm, v₀=227 m/s. Round section as a
structured squircle hex lattice (no gmsh). Half-model: the impact face (z=0) is a
frictionless rigid wall (uz=0, lateral free); the whole bar is given the initial
axial velocity.
"""
import math
import os
import tempfile

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t3]

# Copper, N–mm–tonne–s
_E, _NU, _RHO = 117000.0, 0.35, 8.93e-9
_SY, _H = 400.0, 100.0
_L0, _R0, _V0 = 32.4, 3.2, 227000.0          # mm, mm, mm/s (227 m/s)
_N, _NZ = 4, 12                              # squircle cross n×n, nz axial


def _squircle(u, v, R):
    return (R * u * math.sqrt(max(0.0, 1.0 - 0.5 * v * v)),
            R * v * math.sqrt(max(0.0, 1.0 - 0.5 * u * u)))


def _run_taylor():
    """Build + run the Taylor bar to rebound. Returns a dict of QoIs and the
    energy time-history (from an EnergyBalanceRecorder)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n, nz, c = _N, _NZ, _N // 2
    idx = {}
    tag = 1
    for k in range(nz + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                x, y = _squircle(-1 + 2 * i / n, -1 + 2 * j / n, _R0)
                ops.node(tag, x, y, _L0 * k / nz)
                idx[(i, j, k)] = tag
                tag += 1
    K = _E / 3.0 / (1 - 2 * _NU)
    G = _E / 2.0 / (1 + _NU)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", _SY, 0.0, 1.0, _H, "-kin", 0, "-rho", _RHO)
    ops.nDMaterial("LogStrain", 2, 1)
    e = 1
    for k in range(nz):
        for j in range(n):
            for i in range(n):
                conn = [idx[(i, j, k)], idx[(i + 1, j, k)], idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                        idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)], idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                ops.element("LadrunoBrick", e, *conn, 2, "-formulation", "bbar", "-geom", "finite", "-lumped")
                e += 1
    # frictionless rigid wall at z=0 (uz=0, lateral free); initial axial velocity
    for j in range(n + 1):
        for i in range(n + 1):
            ops.fix(idx[(i, j, 0)], 0, 0, 1)
    for k in range(nz + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.setNodeVel(idx[(i, j, k)], 3, -_V0, "-commit")

    efile = os.path.join(tempfile.gettempdir(), "ladruno_taylor_energy.txt")
    if os.path.exists(efile):
        os.remove(efile)
    ops.recorder("EnergyBalance", "-file", efile, "-time")

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-9)                    # prime → triggers the dt_cr compute
    dtcr_start = ops.criticalTimeStep()
    dt = 0.3 * dtcr_start                     # GEOM-2: reference-config dt_cr ⇒ margin

    topn = idx[(c, c, nz)]
    rebound_vel = None
    for st in range(8000):
        if ops.analyze(1, dt) != 0:
            break
        if st > 300 and st % 50 == 0 and ops.nodeVel(topn, 3) > 0.0:
            rebound_vel = ops.nodeVel(topn, 3)
            break
    dtcr_end = ops.criticalTimeStep()

    def zc(k):
        return np.mean([ops.nodeCoord(idx[(i, j, k)], 3) + ops.nodeDisp(idx[(i, j, k)], 3)
                        for i in range(n + 1) for j in range(n + 1)])

    def rad(k):
        bnd = [(i, j) for j in range(n + 1) for i in range(n + 1) if i in (0, n) or j in (0, n)]
        return np.mean([math.hypot(ops.nodeCoord(idx[(i, j, k)], 1) + ops.nodeDisp(idx[(i, j, k)], 1),
                                   ops.nodeCoord(idx[(i, j, k)], 2) + ops.nodeDisp(idx[(i, j, k)], 2))
                        for (i, j) in bnd])

    res = dict(Lf=float(zc(nz) - zc(0)), r_impact=float(rad(0)), r_free=float(rad(nz)),
               dtcr_start=dtcr_start, dtcr_end=dtcr_end, rebound_vel=rebound_vel,
               t_end=ops.getTime())
    ops.wipe()                                # flush/close the recorder
    d = np.atleast_2d(np.loadtxt(efile))      # cols: time KE IE DW ULW RES ERR%
    res["KE"] = d[:, 1]
    res["IE"] = d[:, 2]
    return res


@pytest.fixture(scope="module")
def taylor():
    r = _run_taylor()
    assert r["rebound_vel"] is not None, "Taylor bar never rebounded — run too short or unstable"
    return r


# =========================================================================== #
#  C2.1 — final length & mushroom radius match the Taylor copper benchmark      #
# =========================================================================== #
def test_C2_taylor_final_length_and_mushroom(taylor):
    Lf_ratio = taylor["Lf"] / _L0
    mush = taylor["r_impact"] / _R0
    # Classic copper Taylor bar @ 227 m/s: L_f ≈ 21.4 mm (L_f/L0 ≈ 0.66), footprint
    # radius ≈ 7 mm (r_f/r0 ≈ 2.2). Engineering band (≤5–7 %, coarse squircle mesh).
    assert 0.62 < Lf_ratio < 0.72, f"Taylor final length L_f/L0 = {Lf_ratio:.3f} (ref ≈ 0.66)"
    assert 2.0 < mush < 2.35, f"Taylor mushroom r_f/r0 = {mush:.3f} (ref ≈ 2.2)"


def test_C2_taylor_deformation_localizes_at_impact(taylor):
    # the mushroom is at the wall; the free end keeps its original radius
    assert taylor["r_free"] == pytest.approx(_R0, rel=0.02), (
        f"free end deformed (r={taylor['r_free']:.3f}, r0={_R0}) — impact should localize")
    assert taylor["r_impact"] > 1.5 * taylor["r_free"], "no mushroom localization at the impact face"


# =========================================================================== #
#  C2.2 / energy — kinetic energy is converted to internal (plastic) energy     #
# =========================================================================== #
def test_P4_energy_balance_through_impact(taylor):
    KE, IE = taylor["KE"], taylor["IE"]
    KE0, KEf = KE[0], KE[-1]
    assert KE0 > 0
    # most of the kinetic energy is absorbed plastically; a small elastic rebound
    assert KEf / KE0 < 0.15, f"too little KE absorbed (KEf/KE0 = {KEf / KE0:.3f}) — not plastic?"
    assert KEf > 0, "no rebound kinetic energy"
    # energy conservation (magnitude): the internal energy gained equals the
    # kinetic energy lost, to engineering tolerance. NOTE the EnergyBalance
    # recorder reports IE with a flipped SIGN for the finite-strain element
    # (LEDGER_quirks) — we compare magnitudes.
    absorbed = KE0 - KEf
    assert abs(abs(IE[-1]) - absorbed) / KE0 < 0.10, (
        f"energy not conserved: |IE|={abs(IE[-1]):.3e} vs KE absorbed={absorbed:.3e}")


# =========================================================================== #
#  dt_cr caveat — criticalTimeStep() is reference-config (GEOM-2)               #
# =========================================================================== #
def test_P4_dtcr_is_reference_config(taylor):
    # the bar compressed ~33 % axially and mushroomed >2×, yet the reported
    # critical time step is UNCHANGED — it is computed from the reference
    # configuration, NOT the deformed one. So an explicit -geom finite run must
    # use a safety factor (< 1) to remain stable as elements compress; this run
    # uses 0.3·dt_cr (documented in the validation report / LEDGER_quirks).
    assert taylor["dtcr_start"] > 0 and math.isfinite(taylor["dtcr_start"])
    assert taylor["dtcr_end"] == pytest.approx(taylor["dtcr_start"], rel=1.0e-6), (
        f"dt_cr changed with deformation ({taylor['dtcr_start']:.4e} → "
        f"{taylor['dtcr_end']:.4e}) — expected reference-config (constant)")
