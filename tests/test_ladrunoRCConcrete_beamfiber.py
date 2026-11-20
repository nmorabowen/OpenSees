"""LadrunoRCConcrete (classTag 33015) -- Ladruno WP-243: the native BeamFiber view.

``getCopy("BeamFiber")`` (what NDFiberSection2d/3d ask for) returns an order-3 view: the fibre
imposes (eps_xx, gamma_xy, gamma_xz) and the lateral block (eps_yy, eps_zz, gamma_yz) is solved so
the concrete balances smeared ties of its own direction, ``-ties rhoY rhoZ fyt Es`` (bilinear,
``-tieHardening b``). No ties is the free-lateral fibre. The view reads the MCFT add-ons from the
full strain (``-addonPlane spatial``) and resolves an unset crack gate to ``damage``; the x-y
interlock is switched off there with a warning. A lateral balance that does not converge REFUSES
the trial (stock BeamFiberMaterial returns 0 there).

The driver's physics and tangent are gated in g++ (tests/_testbed/rc_beamfiber_gpp.cpp: closed
forms, two-leg vs nested bisection, condensed FD, invariance, crack gate, tie law;
test_ladrunoRCConcrete_beamfiber_cpp.py). This file is the wiring: the view as the stock fibre
section drives it, in openseespy and the classic Tcl exe, plus the v8 wire format and the parser.

The fibre is one NDFiber fibre of area A at the centroid of a zeroLengthSection, so the section
force is A times the fibre stress: P = A sig_xx, Vy = A tau_xy, Vz = A tau_xz.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
G = E / (2.0 * (1.0 + NU))
LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]
MCFT = ["-beta", "-tensStiff", "vc"]
FY, ES = 420.0, 200000.0
A = 100.0


def _rc(tag, *flags):
    ops.nDMaterial("LadrunoRCConcrete", tag, E, NU,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC, *flags)


def _build(flags, path):
    """One fibre on a zeroLengthSection; node 2 follows path = [(eps, gy, gz), ...] (one entry per
    step, time = step number). Rotations fixed (a centroidal fibre has no bending stiffness)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 0, 0, 0, 1, 1, 1)
    _rc(1, *flags)
    ops.section("NDFiber", 1, "-GJ", 1.0e12)   # torsion is fixed; classic Tcl demands -GJ in 3D
    ops.fiber(0.0, 0.0, A, 1)
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    t = [0.0] + [float(k + 1) for k in range(len(path))]
    for d in range(3):
        v = [0.0] + [p[d] for p in path]
        ops.timeSeries("Path", d + 1, "-time", *t, "-values", *v)
        ops.pattern("Plain", d + 1, d + 1)
        ops.sp(2, d + 1, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _probe():
    f = list(ops.eleResponse(1, "section", "force"))          # (P, Mz, My, Vy, Vz, T)
    tie = list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 1, "tieStress"))
    lat = list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 1, "lateralStrain"))
    return dict(s=f[0] / A, txy=f[3] / A, txz=f[4] / A, tie=tie, lat=lat)


def _run(flags, path):
    _build(flags, path)
    out = []
    for _ in path:
        if ops.analyze(1) != 0:
            return out, False
        out.append(_probe())
    return out, True


def _ramp(n, eps=0.0, gy=0.0, gz=0.0):
    return [(eps * k / n, gy * k / n, gz * k / n) for k in range(1, n + 1)]


# ----------------------------------------------------------------------------- elastic closed forms
def test_elastic_free_lateral_is_uniaxial_stress():
    """No ties: sig = E eps, tau = G gamma, eps_yy = eps_zz = -nu eps (closed form)."""
    out, ok = _run(MCFT, [(-2.0e-5, 1.0e-5, -0.5e-5)])
    assert ok
    r = out[-1]
    assert r["s"] == pytest.approx(E * -2.0e-5, rel=1e-9)
    assert r["txy"] == pytest.approx(G * 1.0e-5, rel=1e-9)
    assert r["txz"] == pytest.approx(G * -0.5e-5, rel=1e-9)
    assert r["lat"][0] == pytest.approx(NU * 2.0e-5, rel=1e-9)
    assert r["lat"][1] == pytest.approx(NU * 2.0e-5, rel=1e-9)
    assert r["tie"] == [0.0, 0.0]


def test_elastic_ties_stiffen_the_fibre():
    """rhoY = rhoZ = rho: eps_lat = -lam eps / (2(lam+mu) + rho Es), E* = lam + 2 mu - 2 lam^2/(...)."""
    rho = 0.01
    k = 2.0 * (LAM + G) + rho * ES
    elat, estar = -LAM * -2.0e-5 / k, LAM + 2 * G - 2 * LAM ** 2 / k
    out, ok = _run(MCFT + ["-ties", rho, rho, FY, ES], [(-2.0e-5, 0.0, 0.0)])
    assert ok
    r = out[-1]
    assert r["s"] == pytest.approx(estar * -2.0e-5, rel=1e-9)
    assert r["lat"][0] == pytest.approx(elat, rel=1e-9)
    assert r["tie"][0] == pytest.approx(ES * elat, rel=1e-9)
    assert estar > E


# ----------------------------------------------------------------------------- cracked shear
def test_cracked_shear_yields_the_ties():
    """Monotonic gamma_xy with rhoY fyt (0.42) below the cracking stress: the fibre cracks, the tie
    yields and stays at fyt (b = 0), and every step converges."""
    out, ok = _run(MCFT + ["-ties", 0.001, 0.0, FY, ES], _ramp(60, gy=0.012))
    assert ok
    assert out[-1]["tie"][0] == pytest.approx(FY, rel=1e-9)
    assert out[-1]["lat"][0] > FY / ES
    assert out[-1]["tie"][1] == 0.0


def test_xy_and_xz_shear_are_the_same_fibre():
    """Equal ties, gamma_xy vs gamma_xz: identical shear and axial histories (spatial plane)."""
    flags = MCFT + ["-ties", 0.004, 0.004, FY, ES]
    a, oka = _run(flags, _ramp(40, eps=-2.0e-4, gy=5.0e-3))
    b, okb = _run(flags, _ramp(40, eps=-2.0e-4, gz=5.0e-3))
    assert oka and okb
    for ra, rb in zip(a, b):
        assert rb["txz"] == pytest.approx(ra["txy"], rel=1e-7, abs=1e-9)
        assert rb["s"] == pytest.approx(ra["s"], rel=1e-7, abs=1e-9)
        assert rb["tie"][1] == pytest.approx(ra["tie"][0], rel=1e-6, abs=1e-6)


def test_cyclic_shear_with_ties_converges_and_dissipates():
    """Three shear cycles under axial compression: every step converges, the loop has positive
    area, and the tie keeps a plastic strain (the lateral strain does not return to zero).
    MCFT beta only: -tensStiff is a monotonic floor (guide section 7) that re-inflates on
    unloading, so under cyclic shear the lateral balance can lose its root and the view refuses
    (test_cyclic_tension_stiffening_is_refused_not_faked)."""
    amp, n = 4.0e-3, 20
    path = []
    for c in range(3):
        for k in range(1, n + 1):
            path.append((-3.0e-4, amp * k / n, 0.0))
        for k in range(1, 2 * n + 1):
            path.append((-3.0e-4, amp - amp * k / n, 0.0))
        for k in range(1, n + 1):
            path.append((-3.0e-4, -amp + amp * k / n, 0.0))
    out, ok = _run(["-beta", "-ties", 0.002, 0.002, FY, ES], path)
    assert ok, f"failed at step {len(out) + 1} of {len(path)}"
    work = 0.0
    for k in range(1, len(path)):
        work += 0.5 * (out[k]["txy"] + out[k - 1]["txy"]) * (path[k][1] - path[k - 1][1])
    assert work > 0.0
    assert abs(out[-1]["lat"][0]) > 1.0e-4


def test_cyclic_tension_stiffening_is_refused_not_faked():
    """The same first cycle with -tensStiff vc: the monotonic floor carries the loading branch,
    and once unloading leaves the lateral balance without a root the material refuses and
    analyze fails (stock BeamFiberMaterial would return the unbalanced stress)."""
    amp, n = 4.0e-3, 20
    path = [(-3.0e-4, amp * k / n, 0.0) for k in range(1, n + 1)]
    path += [(-3.0e-4, amp - amp * k / n, 0.0) for k in range(1, 2 * n + 1)]
    out, ok = _run(MCFT + ["-ties", 0.002, 0.002, FY, ES], path)
    assert len(out) >= n, f"the monotonic branch did not converge (step {len(out) + 1})"
    assert not ok


# ----------------------------------------------------------------------------- crack gate in the view
def test_unconfined_compression_uses_the_damage_gate():
    """No ties, axial compression to 0.9 fc. The view's default gate (damage) never opens the
    cracked-nu latch. An explicit -crackGate strain opens it from Poisson expansion (~0.5 fc), the
    lateral balance has no root, and the REFUSAL reaches the solver (analyze fails)."""
    path = _ramp(60, eps=-0.0012)
    out, ok = _run(MCFT + ["-crackedNu", 0.0], path)
    assert ok
    assert -out[-1]["s"] > 24.0
    out_s, ok_s = _run(MCFT + ["-crackedNu", 0.0, "-crackGate", "strain"], path)
    assert not ok_s
    assert -out_s[-1]["s"] < 20.0


def test_interlock_is_switched_off_with_a_warning(capfd):
    capfd.readouterr()
    out, ok = _run(MCFT + ["-interlock"], _ramp(10, gy=2.0e-3))
    err = capfd.readouterr().err
    assert ok
    assert "switched off in the BeamFiber view" in err, err


# ----------------------------------------------------------------------------- wire format
def test_serialization_v8_carries_ties_and_tie_history():
    """Schema v8: the ties and the committed tie plastic strains survive a database round trip.
    The path yields the tie, then ends on an elastic unloading step."""
    path = _ramp(30, gy=0.008) + [(0.0, 0.0078, 0.0)]

    def build():
        _build(MCFT + ["-ties", 0.001, 0.001, FY, ES, "-tieHardening", 0.01], path)
        for _ in path:
            assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="rc_wp243",
                       probe_fn=lambda: [_probe()["txy"], *_probe()["tie"], *_probe()["lat"]])


# ----------------------------------------------------------------------------- Model C crack check
CC = ["-crackCheck", "sum", 150.0, 100.0]
SHEAR = [(-3.0e-4, 4.0e-3 * k / 40, 0.0) for k in range(1, 41)]


def _cc_flag():
    return float(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 1, "crackCheck")[0])


def test_crack_check_sum_caps_the_principal_tension():
    """Model C `sum` on a TS shear ramp under axial compression: every step converges, the cap
    binds (recorder crackCheck = 1) and the shear ends lower than the unchecked fibre's (the cap
    trims the tension-stiffening floor, as in the rcbf-0 oracle)."""
    out0, ok0 = _run(MCFT + ["-ties", 0.002, 0.002, FY, ES], SHEAR)
    assert ok0
    assert _cc_flag() == 0.0
    out1, ok1 = _run(MCFT + ["-ties", 0.002, 0.002, FY, ES] + CC, SHEAR)
    assert ok1
    assert _cc_flag() == 1.0
    assert out1[-1]["txy"] < 0.95 * out0[-1]["txy"], (out1[-1]["txy"], out0[-1]["txy"])


def test_serialization_v8_carries_the_crack_check():
    """Schema v8: -crackCheck and its spacings survive a database round trip. The path ends on a
    small shear unload where the cap still binds (the TS floor re-inflates): re-evaluation there
    is a fixed point, while on a loading step the plastic-damage return is not idempotent (xt
    advances), with or without the check. The unchecked fibre's shear is 4.18, so a restore that
    lost the check would not reproduce 3.36."""
    path = SHEAR + [(-3.0e-4, 3.98e-3, 0.0)]

    def build():
        _build(MCFT + ["-ties", 0.002, 0.002, FY, ES] + CC, path)
        for _ in path:
            assert ops.analyze(1) == 0
        assert _cc_flag() == 1.0

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="rc_wp243_cc",
                       probe_fn=lambda: [_probe()["txy"], *_probe()["tie"], *_probe()["lat"], _cc_flag()])


# ----------------------------------------------------------------------------- parser
@pytest.mark.parametrize("flags", [
    ["-ties", -0.01, 0.0, FY, ES],
    ["-ties", 0.01, 1.5, FY, ES],
    ["-ties", 0.01, 0.01, 0.0, ES],
    ["-ties", 0.01, 0.01, FY, 0.0],
    ["-ties", 0.01, 0.01, FY],
    ["-tieHardening", 1.0],
    ["-tieHardening", -0.1],
    ["-crackCheck", "sum", 0.0, 100.0],
    ["-crackCheck", "sum", 150.0],
    ["-crackCheck", "min", 150.0, 100.0],
])
def test_parser_refuses(flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    with pytest.raises(Exception):
        _rc(1, *flags)


def test_parser_accepts():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _rc(1, "-ties", 0.004, 0.002, FY, ES, "-tieHardening", 0.01)
    _rc(2, *MCFT, "-ties", 0.0, 0.0, FY, ES)
    _rc(3, *MCFT, "-ties", 0.002, 0.002, FY, ES, "-crackCheck", "sum", 150.0, 100.0)


# ----------------------------------------------------------------------------- classic Tcl engine
def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_engine_gives_the_same_fibre(tmp_path):
    flags = MCFT + ["-ties", 0.001, 0.0, FY, ES]
    path = _ramp(40, eps=-1.0e-4, gy=0.01)
    py, ok = _run(flags, path)
    assert ok
    bb = " ".join(str(x) for x in ["-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC])
    fl = " ".join(str(x) for x in flags)
    t = " ".join(str(x) for x in [0.0] + [float(k + 1) for k in range(len(path))])
    lines = [
        "model basic -ndm 3 -ndf 6",
        "node 1 0 0 0", "node 2 0 0 0",
        "fix 1 1 1 1 1 1 1", "fix 2 0 0 0 1 1 1",
        f"nDMaterial LadrunoRCConcrete 1 {E} {NU} {bb} {fl}",
        "section NDFiber 1 -GJ 1.0e12 {", f"  fiber 0.0 0.0 {A} 1", "}",
        "element zeroLengthSection 1 1 2 1",
    ]
    for d in range(3):
        v = " ".join(str(x) for x in [0.0] + [p[d] for p in path])
        lines += [f"timeSeries Path {d + 1} -time {{{t}}} -values {{{v}}}",
                  f"pattern Plain {d + 1} {d + 1} {{ sp 2 {d + 1} 1.0 }}"]
    lines += ["system FullGeneral", "numberer Plain", "constraints Transformation",
              "test NormDispIncr 1.0e-12 30 0", "algorithm Newton", "integrator LoadControl 1.0",
              "analysis Static",
              f"for {{set k 0}} {{$k < {len(path)}}} {{incr k}} {{",
              "  if {[analyze 1] != 0} { puts \"FAILED $k\"; exit 1 }",
              "  set f [eleResponse 1 section force]",
              "  puts \"F [lindex $f 0] [lindex $f 3]\"",
              "}"]
    deck = tmp_path / "beamfiber.tcl"
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True, timeout=300,
                          env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    rows = [l.split() for l in (proc.stdout + proc.stderr).splitlines() if l.startswith("F ")]
    assert len(rows) == len(path), proc.stdout + proc.stderr
    for r, p in zip(rows, py):
        assert float(r[1]) / A == pytest.approx(p["s"], rel=1e-9, abs=1e-9)
        assert float(r[2]) / A == pytest.approx(p["txy"], rel=1e-9, abs=1e-9)
