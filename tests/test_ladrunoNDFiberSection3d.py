"""LadrunoNDFiberSection3d (SEC 33001) and LadrunoUniaxialFiber (nD 33030) -- Ladruno WP-245
(program #1131, rcbf-3).

``section LadrunoNDFiber`` is the stock NDFiber 3D section with a chosen shear-strain profile:
``-shearProfile uniform`` is the stock sqrt(alpha) rule; ``-shearProfile parabolic -shearDepth hy hz``
imposes gamma_xy = 1.5 (1 - (2y/hy)^2) gamma_y (and the same in z), as Model C does, with the
resultant ``-shearResultant modelC`` (V = sum tau A, the default, unsymmetric) or ``conjugate``
(V = sum tau phi A, symmetric, k = 6/5 in the continuum). ``nDMaterial LadrunoUniaxialFiber``
puts a uniaxial bar in the same section, carrying no shear.

The shear map is gated in g++ (tests/_testbed/rc_section_gpp.cpp via
test_ladrunoNDFiberSection3d_cpp.py). This file is the wiring through OpenSees: the section on a
zeroLengthSection, the Model C assembly against independently driven strips, the bar adapter, the
parser, the wire format and the classic Tcl exe.

Model C's layout (06 §1.1): b x h rectangle, 40 strips through the depth, fibre at each centroid.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
E, NU = 30000.0, 0.2
G = E / (2.0 * (1.0 + NU))
B = H = 457.2
NSTRIP = 40
DY = H / NSTRIP
YS = [-0.5 * H + DY * (i + 0.5) for i in range(NSTRIP)]
AS = B * DY
PARA = ["-shearProfile", "parabolic", "-shearDepth", H, B]

# RC strip law (as the WP-243 tests): plastic-damage + MCFT, tension stiffening, BeamFiber view
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]
FY, ES = 420.0, 200000.0


def _rc(tag, *flags):
    ops.nDMaterial("LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", 2.0 / 3.0,
                   "-beta", "-tensStiff", "vc", *flags)


# section deformation (P, Mz, My, Vy, Vz, T) <- zeroLengthSection node-2 DOF
DOF_OF = {0: 1, 1: 6, 2: 5, 3: 2, 4: 3, 5: 4}


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def _drive(pairs, path=None):
    """pairs: [(eleTag, nodeI, nodeJ, secTag, rows)]; node J of each zeroLengthSection follows
    rows[k] (6 section deformations per step; rows=None uses `path`). One Path series per
    (element, component); the analysis is defined after every element exists."""
    ts = 1
    for ele, ni, nj, sec, rows in pairs:
        rows = path if rows is None else rows
        t = [0.0] + [float(k + 1) for k in range(len(rows))]
        ops.node(ni, 0.0, 0.0, 0.0)
        ops.node(nj, 0.0, 0.0, 0.0)
        ops.fix(ni, 1, 1, 1, 1, 1, 1)
        ops.element("zeroLengthSection", ele, ni, nj, sec)
        for c in range(6):
            ops.timeSeries("Path", ts, "-time", *t, "-values", *([0.0] + [r[c] for r in rows]))
            ops.pattern("Plain", ts, ts)
            ops.sp(nj, DOF_OF[c], 1.0)
            ts += 1
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _force(ele):
    return list(ops.eleResponse(ele, "section", "force"))


def _defo(ele):
    return list(ops.eleResponse(ele, "section", "deformation"))


def _strips(sec, mat, *opts, bars=()):
    ops.section("LadrunoNDFiber", sec, *opts)
    for y in YS:
        ops.fiber(y, 0.0, AS, mat)
    for (y, a, m) in bars:
        ops.fiber(y, 0.0, a, m)


# ----------------------------------------------------------------------------- elastic closed forms
def _elastic_k(opts, comp):
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    _strips(1, 1, *opts)
    e = [0.0] * 6
    e[comp] = 1.0e-5
    _drive([(1, 1, 2, 1, None)], [e])
    assert ops.analyze(1) == 0
    return _force(1)[comp] / _defo(1)[comp]


def test_elastic_shear_stiffness_closed_forms():
    """40 strips, ElasticIsotropic: modelC K_VV = G A (1 + 1/(2n^2)) (Model C's k = 1);
    conjugate K_VV = 2.25 G A (1 - 2 mu2 + mu4) (-> 6/5); uniform -alpha a gives a G A; EI and EA
    are the midpoint sums."""
    n, A = NSTRIP, B * H
    mu2 = 1.0 / 3.0 - 1.0 / (3.0 * n * n)
    mu4 = 0.2 - 2.0 / (3.0 * n * n) + 7.0 / (15.0 * n ** 4)
    assert _elastic_k(PARA, 3) == pytest.approx(G * A * (1 + 1.0 / (2 * n * n)), rel=1e-12)
    kc = _elastic_k(PARA + ["-shearResultant", "conjugate"], 3)
    assert kc == pytest.approx(2.25 * G * A * (1 - 2 * mu2 + mu4), rel=1e-12)
    assert kc / (G * A) == pytest.approx(1.2, abs=1e-4)
    assert _elastic_k(["-shearProfile", "uniform", "-alpha", 0.83], 3) == pytest.approx(0.83 * G * A, rel=1e-12)
    assert _elastic_k(PARA, 0) == pytest.approx(E * A, rel=1e-12)
    assert _elastic_k(PARA, 1) == pytest.approx(E * B * H ** 3 / 12 * (1 - 1.0 / n ** 2), rel=1e-12)


# ----------------------------------------------------------------------------- Model C assembly
def _rc_section_path():
    """Axial compression, a little curvature, then a shear ramp through cracking."""
    return [(-3.0e-4, 1.0e-6, 0.0, 4.0e-3 * k / 40.0, 0.0, 0.0) for k in range(1, 41)]


def test_model_c_assembly_through_opensees():
    """The parabolic/modelC section of 40 RC strips (ties in the core only, Model C's layout) vs 40
    stock single-fibre NDFiber sections, each driven DIRECTLY at its strip strain
    (eps0 - y kappa, phi(y) gamma): N = sum sig A, M = -sum sig A y, V = sum tau A (06 §1.2, §1.4),
    every step of a cracking shear path, to roundoff."""
    path = _rc_section_path()
    cover = 63.5
    _model()
    _rc(1)                                                   # cover: no ties
    _rc(2, "-ties", 0.002, 0.0, FY, ES)                      # core: ties in y
    core = [abs(y) <= H / 2 - cover for y in YS]
    ops.section("LadrunoNDFiber", 1, *PARA)
    for y, c in zip(YS, core):
        ops.fiber(y, 0.0, AS, 2 if c else 1)
    pairs = [(1, 1, 2, 1, path)]
    for i, (y, c) in enumerate(zip(YS, core)):
        ops.section("NDFiber", 100 + i, "-GJ", 1.0e12)
        ops.fiber(0.0, 0.0, AS, 2 if c else 1)
        phi = 1.5 * (1.0 - (2.0 * y / H) ** 2)
        strip = [[p[0] - y * p[1], 0.0, 0.0, phi * p[3], 0.0, 0.0] for p in path]
        pairs.append((100 + i, 1000 + 2 * i, 1001 + 2 * i, 100 + i, strip))
    _drive(pairs)
    worst, vmax = 0.0, 0.0
    for k in range(len(path)):
        assert ops.analyze(1) == 0, k
        f = _force(1)
        N = sum(_force(100 + i)[0] for i in range(NSTRIP))
        M = -sum(_force(100 + i)[0] * y for i, y in enumerate(YS))
        V = sum(_force(100 + i)[3] for i in range(NSTRIP))
        for a, b in ((f[0], N), (f[1], M), (f[3], V)):
            worst = max(worst, abs(a - b) / (abs(b) + 1.0))
        vmax = max(vmax, f[3])
    assert worst < 1e-10, worst
    assert vmax > 0.5e6                      # the path cracks and the section carries real shear (N)


def test_uniform_profile_reproduces_stock_ndfiber():
    """-shearProfile uniform (alpha 1) gives the stock NDFiber section on the same RC strips."""
    path = _rc_section_path()
    _model()
    _rc(1, "-ties", 0.002, 0.0, FY, ES)
    _strips(1, 1, "-shearProfile", "uniform")
    ops.section("NDFiber", 2, "-GJ", 1.0e12)
    for y in YS:
        ops.fiber(y, 0.0, AS, 1)
    _drive([(1, 1, 2, 1, None), (2, 3, 4, 2, None)], path)
    worst = 0.0
    for k in range(len(path)):
        assert ops.analyze(1) == 0, k
        f1, f2 = _force(1), _force(2)
        worst = max(worst, max(abs(a - b) / (abs(b) + 1.0) for a, b in zip(f1, f2)))
    assert worst < 1e-11, worst


def test_parabolic_profile_is_more_shear_flexible_at_the_faces():
    """Same RC strips and shear path: the parabolic section cracks later at the faces and earlier at
    mid-depth than the uniform one (Model C's reason for the profile, rcbf-0 finding 2), so the two
    shear histories differ."""
    path = _rc_section_path()
    _model()
    _rc(1, "-ties", 0.002, 0.0, FY, ES)
    _strips(1, 1, *PARA)
    _strips(2, 1, "-shearProfile", "uniform")
    _drive([(1, 1, 2, 1, None), (2, 3, 4, 2, None)], path)
    vp, vu = [], []
    for _ in path:
        assert ops.analyze(1) == 0
        vp.append(_force(1)[3])
        vu.append(_force(2)[3])
    assert max(abs(a - b) / abs(b) for a, b in zip(vp, vu)) > 0.02


# ----------------------------------------------------------------------------- bar fibre
def test_bar_fibre_carries_axial_force_and_no_shear():
    """LadrunoUniaxialFiber(Steel01) bars: the axial and bending resultants gain Es A_bar (and its
    lever arm); the shear resultant and shear stiffness are unchanged; the bar's fibre stress is
    (sigma, 0, 0)."""
    abar, yb = 2000.0, 160.0
    res = {}
    for with_bars in (False, True):
        _model()
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
        ops.uniaxialMaterial("Steel01", 5, FY, ES, 0.01)
        ops.nDMaterial("LadrunoUniaxialFiber", 6, 5)
        bars = [(-yb, abar, 6), (yb, abar, 6)] if with_bars else []
        _strips(1, 1, *PARA, bars=bars)
        e = (-1.0e-4, 0.0, 0.0, 2.0e-4, 0.0, 0.0)
        _drive([(1, 1, 2, 1, None)], [e])
        assert ops.analyze(1) == 0
        res[with_bars] = (_force(1), _defo(1))
        if with_bars:
            sb = list(ops.eleResponse(1, "section", "fiber", yb, 0.0, 6, "stress"))
    f0, f1 = res[False][0], res[True][0]
    eps = res[True][1][0]
    assert f1[0] - f0[0] == pytest.approx(2 * abar * ES * eps, rel=1e-10)
    assert f1[3] == pytest.approx(f0[3], rel=1e-13)
    assert sb[0] == pytest.approx(ES * eps, rel=1e-10)
    assert sb[1] == 0.0 and sb[2] == 0.0


def test_bar_fibre_follows_the_uniaxial_law():
    """A bar fibre in axial tension past yield: the section force is Steel01's, fibre by fibre."""
    _model()
    ops.uniaxialMaterial("Steel01", 5, FY, ES, 0.01)
    ops.nDMaterial("LadrunoUniaxialFiber", 6, 5)
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "uniform")
    ops.fiber(0.0, 0.0, 100.0, 6)
    path = [(4.0e-3 * k / 10.0, 0, 0, 0, 0, 0) for k in range(1, 11)]
    _drive([(1, 1, 2, 1, None)], path)
    ey = FY / ES
    for p in path:
        assert ops.analyze(1) == 0
        sig = ES * p[0] if p[0] <= ey else FY + 0.01 * ES * (p[0] - ey)
        assert _force(1)[0] == pytest.approx(100.0 * sig, rel=1e-10)


# ----------------------------------------------------------------------------- element
@pytest.mark.parametrize("res", ["modelC", "conjugate"])
def test_force_beam_column_cantilever_shear_deflection(res):
    """forceBeamColumn cantilever, elastic strips (two fibres per strip at z = +-b/4 so My is not
    singular): tip deflection = P L^3/(3 EI) + P L/(k G A) with the section's own k (modelC
    1 + 1/(2n^2), conjugate 2.25 (1 - 2 mu2 + mu4)) -- the VY row reaches the element. The
    unsymmetric modelC tangent runs on UmfPack."""
    L, P, n, A = 1500.0, 1.0e5, NSTRIP, B * H
    EI = E * B * H ** 3 / 12 * (1 - 1.0 / n ** 2)
    mu2 = 1.0 / 3.0 - 1.0 / (3.0 * n * n)
    mu4 = 0.2 - 2.0 / (3.0 * n * n) + 7.0 / (15.0 * n ** 4)
    k = 1 + 1.0 / (2 * n * n) if res == "modelC" else 2.25 * (1 - 2 * mu2 + mu4)
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.section("LadrunoNDFiber", 1, *PARA, "-shearResultant", res)
    for y in YS:
        for z in (-B / 4, B / 4):
            ops.fiber(y, z, AS / 2, 1)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0, 0, 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, P, 0.0, 0.0, 0.0, 0.0)
    ops.system("UmfPack")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ref = P * L ** 3 / (3 * EI) + P * L / (k * G * A)
    assert ops.nodeDisp(2, 2) == pytest.approx(ref, rel=1e-12)


# ----------------------------------------------------------------------------- parser
def _section_refused(*args):
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    with pytest.raises(Exception):
        ops.section("LadrunoNDFiber", 1, *args)


@pytest.mark.parametrize("args", [
    (),                                                       # -shearProfile is required
    ("-shearProfile", "cubic"),
    ("-shearProfile", "parabolic"),                           # needs -shearDepth
    ("-shearProfile", "parabolic", "-shearDepth", 0.0, 1.0),
    ("-shearProfile", "parabolic", "-shearDepth", H, B, "-alpha", 0.8),
    ("-shearProfile", "parabolic", "-shearDepth", H, B, "-shearResultant", "energy"),
    ("-shearProfile", "uniform", "-shearDepth", H, B),
    ("-shearProfile", "uniform", "-shearResultant", "conjugate"),
    ("-shearProfile", "uniform", "-alpha", -1.0),
    ("-shearProfile", "uniform", "-GJ", 1.0),                 # an ND section takes no -GJ
    ("-shearProfile", "uniform", "-bogus"),
])
def test_parser_refuses(args):
    _section_refused(*args)


def test_parser_accepts():
    _model()
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "uniform", "-alpha", 0.9, "-noCentroid")
    ops.section("LadrunoNDFiber", 2, *PARA, "-shearResultant", "conjugate")
    ops.section("LadrunoNDFiber", 3, *PARA, "-shearResultant", "modelC")


def test_section_refuses_2d_models():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    with pytest.raises(Exception):
        ops.section("LadrunoNDFiber", 1, "-shearProfile", "uniform")


def test_adapter_parser():
    _model()
    ops.uniaxialMaterial("Elastic", 5, ES)
    ops.nDMaterial("LadrunoUniaxialFiber", 6, 5, "-rho", 7.85e-9)
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoUniaxialFiber", 7, 99)         # no such uniaxial material
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoUniaxialFiber", 8, 5, "-bogus", 1.0)


def test_fibre_outside_the_parabolic_band_is_refused():
    """A fibre past the face (|y| > hy/2 from the centroid) would get a negative phi: the first
    state determination refuses."""
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "parabolic", "-shearDepth", 200.0, 200.0)
    ops.fiber(-150.0, 0.0, 100.0, 1)
    ops.fiber(150.0, 0.0, 100.0, 1)
    _drive([(1, 1, 2, 1, None)], [(1e-5, 0, 0, 1e-5, 0, 0)])
    assert ops.analyze(1) != 0


# ----------------------------------------------------------------------------- wire format
def test_serialization_carries_the_profile_and_the_bars():
    """The section's profile, resultant and depths, and the bar adapter's uniaxial state, survive
    a database round trip (two rounds: a mistagged null object would show on the second). Elastic
    concrete under shear with the conjugate resultant (a restore that fell back to modelC or to
    the uniform profile changes V by 20 % or more), bars yielded in tension and then unloaded
    (Steel01 keeps the plastic strain). The RC strips' own wire is WP-243's test: an RC path is
    not used here because the BeamFiber view is not a fixed point at a held strain on some
    cracked states (re-solving the committed strain moves the stress), with or without a round
    trip."""
    path = [(3.0e-3 * k / 10.0, 0.0, 0.0, 1.0e-3 * k / 10.0, 4.0e-4 * k / 10.0, 0.0) for k in range(1, 11)]
    path += [(2.8e-3, 0.0, 0.0, 0.9e-3, 3.6e-4, 0.0)]

    def build():
        _model()
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
        ops.uniaxialMaterial("Steel01", 5, FY, ES, 0.01)
        ops.nDMaterial("LadrunoUniaxialFiber", 6, 5)
        _strips(1, 1, *PARA, "-shearResultant", "conjugate", bars=[(-160.0, 1500.0, 6), (160.0, 1500.0, 6)])
        _drive([(1, 1, 2, 1, None)], path)
        for _ in path:
            assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="rc_wp245", rounds=2,
                       probe_fn=lambda: _force(1) + list(ops.eleResponse(1, "section", "fiber", 160.0, 0.0, 6, "stress")))


# ----------------------------------------------------------------------------- classic Tcl
def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_engine_gives_the_same_section(tmp_path):
    path = _rc_section_path()[:20]
    bars = [(-160.0, 1500.0, 6), (160.0, 1500.0, 6)]

    _model()
    _rc(1, "-ties", 0.002, 0.0, FY, ES)
    ops.uniaxialMaterial("Steel01", 5, FY, ES, 0.01)
    ops.nDMaterial("LadrunoUniaxialFiber", 6, 5)
    _strips(1, 1, *PARA, bars=bars)
    _drive([(1, 1, 2, 1, None)], path)
    py = []
    for _ in path:
        assert ops.analyze(1) == 0
        py.append(_force(1))

    bb = " ".join(str(x) for x in ["-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                                   "-Kc", 2.0 / 3.0, "-beta", "-tensStiff", "vc", "-ties", 0.002, 0.0, FY, ES])
    fib = [f"  fiber {y} 0.0 {AS} 1" for y in YS] + [f"  fiber {y} 0.0 {a} {m}" for (y, a, m) in bars]
    t = " ".join(str(x) for x in [0.0] + [float(k + 1) for k in range(len(path))])
    lines = [
        "model basic -ndm 3 -ndf 6",
        f"nDMaterial LadrunoRCConcrete 1 {E} {NU} {bb}",
        f"uniaxialMaterial Steel01 5 {FY} {ES} 0.01",
        "nDMaterial LadrunoUniaxialFiber 6 5",
        f"section LadrunoNDFiber 1 -shearProfile parabolic -shearDepth {H} {B} {{", *fib, "}",
        "node 1 0 0 0", "node 2 0 0 0", "fix 1 1 1 1 1 1 1",
        "element zeroLengthSection 1 1 2 1",
    ]
    for c in range(6):
        v = " ".join(str(x) for x in [0.0] + [p[c] for p in path])
        lines += [f"timeSeries Path {c + 1} -time {{{t}}} -values {{{v}}}",
                  f"pattern Plain {c + 1} {c + 1} {{ sp 2 {DOF_OF[c]} 1.0 }}"]
    lines += ["system FullGeneral", "numberer Plain", "constraints Transformation",
              "test NormDispIncr 1.0e-12 30 0", "algorithm Newton", "integrator LoadControl 1.0",
              "analysis Static",
              f"for {{set k 0}} {{$k < {len(path)}}} {{incr k}} {{",
              "  if {[analyze 1] != 0} { puts \"FAILED $k\"; exit 1 }",
              "  puts \"F [eleResponse 1 section force]\"",
              "}"]
    deck = tmp_path / "section.tcl"
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True, timeout=300,
                          env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    rows = [l.split()[1:] for l in (proc.stdout + proc.stderr).splitlines() if l.startswith("F ")]
    assert len(rows) == len(path), proc.stdout + proc.stderr
    for r, p in zip(rows, py):
        for a, b in zip(r, p):
            assert float(a) == pytest.approx(b, rel=1e-9, abs=1e-6)
