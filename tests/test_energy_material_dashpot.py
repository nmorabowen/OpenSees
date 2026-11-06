"""WP-184: EnergyBalance books a material dashpot once (IE), not twice (IE + DW).

A rate-dependent material (Viscous, Parallel(..., Viscous)) puts its viscous
force in getStress(), hence in the element resisting force F, so its power is
already in IE = int F.v dt. The element's getDamp() ALSO returns the material
damping tangent; before WP-184 the recorder booked v^T C v of it in DW as well,
so RES = -(dashpot dissipation). DW is now the power of the damping force the
element adds outside F (g = F_incInertia - F - M a); Rayleigh-only elements keep
the original v^T C v.

Classic-Tcl gates (OpenSees.exe), SDOF m = 1 on zeroLength / truss, Newmark
average acceleration, load pulse then free vibration:

  T1  Elastic + Viscous (one zeroLength, two materials), Parallel(Elastic,
      Viscous) -factors (the SSI-spring form), and a separate dashpot element:
      DW = 0, IE = trapz(k u v) + trapz(c v^2), RES ~ 0.
  T2  Rayleigh damping instead of the dashpot (region betaK on the spring):
      DW = trapz(betaK k v^2) (the original v^T C v path), RES ~ 0.
  T3  Mixed: a truss whose material is Parallel(Elastic, Viscous) with Rayleigh
      betaK on the same element: DW books only the Rayleigh part, IE the dashpot.
  T4  Nonlinear dashpot c|v|^a, a = 0.5: DW = 0, RES ~ 0 (before: DW = a x D).
  T5  -v2: DW_ele follows the v1 DW; RES identical.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")

M, K, C, DT, N = 1.0, 100.0, 2.0, 0.001, 4000
PULSE = "timeSeries Path 1 -time {0.0 0.05 0.1 100.0} -values {0.0 10.0 0.0 0.0}"


def run(tmp, name, model, ndm=1):
    d = tmp / name
    d.mkdir()
    if ndm == 1:
        head = ["wipe", "model basic -ndm 1 -ndf 1", "node 1 0.0", f"node 2 0.0 -mass {M!r}", "fix 1 1"]
    else:
        head = ["wipe", "model basic -ndm 2 -ndf 2", "node 1 0.0 0.0", f"node 2 1.0 0.0 -mass {M!r} {M!r}",
                "fix 1 1 1", "fix 2 0 1"]
    deck = head + model + [
        PULSE, "pattern Plain 1 1 { load 2 " + ("1.0" if ndm == 1 else "1.0 0.0") + " }",
        "recorder Node -file u.txt -time -precision 17 -node 2 -dof 1 disp",
        "recorder Node -file v.txt -time -precision 17 -node 2 -dof 1 vel",
        "recorder EnergyBalance -time -precision 17 -file e.txt",
        "recorder EnergyBalance -time -precision 17 -v2 -file e2.txt",
        "constraints Plain", "numberer Plain", "system BandGeneral", "test NormDispIncr 1e-14 30",
        "algorithm Newton", "integrator Newmark 0.5 0.25", "analysis Transient",
        f"set ok [analyze {N} {DT!r}]", 'puts "done $ok"', "wipe",
    ]
    (d / "main.tcl").write_text("\n".join(deck) + "\n")
    r = subprocess.run([TCL_EXE, "main.tcl"], cwd=d, capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    assert "done 0" in out, out[-2000:]
    u = np.r_[0.0, np.loadtxt(d / "u.txt")[:, 1]]
    v = np.r_[0.0, np.loadtxt(d / "v.txt")[:, 1]]
    return u, v, np.loadtxt(d / "e.txt"), np.loadtxt(d / "e2.txt")


def trapz(f):
    return float(np.sum(0.5 * (f[:-1] + f[1:]) * DT))


def check_closure(e, d_visc):
    """RES is the trapezoidal-quadrature residual of the Newmark run (O(dt^2),
    ~1e-7 of the energy scale here; the Rayleigh-only model shows the same).
    Before WP-184 a dashpot made RES = -D, i.e. ~0.5 of the scale."""
    ref = np.max(np.abs(e[:, 1:5]).sum(axis=1))
    assert abs(e[-1, 5]) < 1e-6 * ref, (e[-1], d_visc)


RAYLEIGH_TWIN = [f"uniaxialMaterial Elastic 1 {K!r}", "element zeroLength 1 1 2 -mat 1 -dir 1 -doRayleigh 1",
                 f"region 1 -ele 1 -rayleigh 0.0 {C / K!r} 0.0 0.0"]

DASHPOT_VARIANTS = {
    "two_materials": [f"uniaxialMaterial Elastic 1 {K!r}", f"uniaxialMaterial Viscous 2 {C!r} 1.0",
                      "element zeroLength 1 1 2 -mat 1 2 -dir 1 1 -doRayleigh 0"],
    "parallel_factors": ["uniaxialMaterial Elastic 11 1.0", "uniaxialMaterial Viscous 12 1.0 1.0",
                         f"uniaxialMaterial Parallel 3 11 12 -factors {K!r} {C!r}",
                         "element zeroLength 1 1 2 -mat 3 -dir 1"],
    "separate_dashpot": [f"uniaxialMaterial Elastic 1 {K!r}", f"uniaxialMaterial Viscous 2 {C!r} 1.0",
                         "element zeroLength 1 1 2 -mat 1 -dir 1", "element zeroLength 2 1 2 -mat 2 -dir 1"],
}


@needs_exe
@pytest.mark.parametrize("variant", sorted(DASHPOT_VARIANTS))
def test_T1_linear_dashpot_booked_once(tmp_path, variant):
    u, v, e, e2 = run(tmp_path, variant, DASHPOT_VARIANTS[variant])
    d_visc = trapz(C * v * v)
    ie_expect = trapz(K * u * v) + d_visc
    assert d_visc > 1e-3
    assert abs(e[-1, 3]) < 1e-14 * d_visc, ("DW must be 0", e[-1])
    assert abs(e[-1, 2] - ie_expect) < 1e-9 * d_visc, (e[-1, 2], ie_expect)
    check_closure(e, d_visc)
    # Oracle: the same dynamics with the dashpot as Rayleigh damping (betaK k = c)
    # gives the same RES (quadrature residual) and DW = the dashpot dissipation.
    ur, vr, er, _ = run(tmp_path, variant + "_rayleigh_twin", RAYLEIGH_TWIN)
    assert np.max(np.abs(vr - v)) < 1e-12 * np.max(np.abs(v))
    assert abs(e[-1, 5] - er[-1, 5]) < 1e-12 * d_visc, (e[-1, 5], er[-1, 5])
    assert abs((e[-1, 2] + e[-1, 3]) - (er[-1, 2] + er[-1, 3])) < 1e-12 * d_visc


@needs_exe
def test_T2_rayleigh_spring_unchanged_path(tmp_path):
    beta = 0.02  # c_eq = beta k = 2.0, the same as the dashpot
    model = [f"uniaxialMaterial Elastic 1 {K!r}", "element zeroLength 1 1 2 -mat 1 -dir 1 -doRayleigh 1",
             f"region 1 -ele 1 -rayleigh 0.0 {beta!r} 0.0 0.0"]
    u, v, e, e2 = run(tmp_path, "rayleigh", model)
    dw_expect = trapz(beta * K * v * v)
    assert abs(e[-1, 3] - dw_expect) < 1e-10 * dw_expect, (e[-1, 3], dw_expect)
    assert abs(e[-1, 2] - trapz(K * u * v)) < 1e-10 * dw_expect
    check_closure(e, dw_expect)


@needs_exe
def test_T3_truss_dashpot_plus_rayleigh(tmp_path):
    beta = 0.01
    model = ["uniaxialMaterial Elastic 11 1.0", "uniaxialMaterial Viscous 12 1.0 1.0",
             f"uniaxialMaterial Parallel 3 11 12 -factors {K!r} {C!r}",
             "element truss 1 1 2 1.0 3 -doRayleigh 1",
             f"region 1 -ele 1 -rayleigh 0.0 {beta!r} 0.0 0.0"]
    u, v, e, e2 = run(tmp_path, "truss_mixed", model, ndm=2)
    d_visc = trapz(C * v * v)
    dw_expect = trapz(beta * K * v * v)
    assert abs(e[-1, 3] - dw_expect) < 1e-8 * dw_expect, ("DW = Rayleigh only", e[-1, 3], dw_expect, d_visc)
    assert abs(e[-1, 2] - (trapz(K * u * v) + d_visc)) < 1e-8 * d_visc
    check_closure(e, d_visc)


@needs_exe
def test_T4_nonlinear_dashpot(tmp_path):
    a = 0.5
    model = [f"uniaxialMaterial Elastic 1 {K!r}", f"uniaxialMaterial Viscous 2 {C!r} {a!r}",
             "element zeroLength 1 1 2 -mat 1 2 -dir 1 1"]
    u, v, e, e2 = run(tmp_path, "nonlinear", model)
    d_visc = trapz(C * np.abs(v) ** a * np.abs(v))
    assert abs(e[-1, 3]) < 1e-14 * d_visc, e[-1]
    assert abs(e[-1, 2] - (trapz(K * u * v) + d_visc)) < 1e-8 * d_visc
    check_closure(e, d_visc)


@needs_exe
def test_T5_v2_columns_follow(tmp_path):
    u, v, e, e2 = run(tmp_path, "v2", DASHPOT_VARIANTS["parallel_factors"])
    # v2: time KE_ele KE_nod IE DW_ele DW_nod ULW [channels] RES ERR
    assert abs(e2[-1, 4]) == 0.0 and abs(e2[-1, 5]) == 0.0, e2[-1]
    assert abs(e2[-1, 3] - e[-1, 2]) < 1e-12 * abs(e[-1, 2])
    assert abs(e2[-1, -2] - e[-1, 5]) < 1e-12 * abs(e[-1, 2])
