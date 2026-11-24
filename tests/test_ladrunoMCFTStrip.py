"""LadrunoMCFTStrip (nD 33031) and LadrunoMCFTSteel (uniaxial 33004) -- Ladruno WP-254 (program #1131).

``nDMaterial LadrunoMCFTStrip`` is the Model C MCFT strip as an NDFiber beam fibre: a rotating 2D
membrane with one condensed transverse strain, balanced against smeared ties, with Model C's root
rule; ``uniaxialMaterial LadrunoMCFTSteel`` is Model C's bar law, put in the section through
``nDMaterial LadrunoUniaxialFiber``. The law itself is gated in g++ against the numpy reference
(tests/_testbed/mcft_strip_gpp.cpp via test_ladrunoMCFTStrip_cpp.py, bit-identical strip paths).

This file is the wiring through OpenSees:
- a Model C column section (``section LadrunoNDFiber -shearProfile parabolic -noCentroid``, 40 strips,
  8 bars) on a zeroLengthSection, driven on a cracking cyclic path with commits, against the reference
  Section.forces to round-off;
- the bar law against the reference on a cyclic path;
- the strip's responses, its views (BeamFiber, BeamFiber2d; others refused), the parsers (accept, refuse
  every bad input), the database round trip (two rounds), a wild trial strain;
- the classic Tcl exe: the same section, and `ndTest` (WP-248) views.
"""
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

sys.path.insert(0, str(Path(__file__).resolve().parent / "_testbed"))
import mcft_strip_ref as ref

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
SEC, DATA = ref.column_section()
CORE = int(np.nonzero(DATA["core"])[0][0])
COVER = 0
TOL = 1e-11          # section forces, relative to fc b h (b h^2 for M)

# section deformation (P, Mz, My, Vy, Vz, T) <- zeroLengthSection node-2 DOF (as WP-245's tests)
DOF_OF = {0: 1, 1: 6, 2: 5, 3: 2, 4: 3, 5: 4}


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def _strip_args(i):
    c = SEC.conc
    a = [c.fc, float(DATA["fp"][i]), float(DATA["ep"][i]), float(DATA["e20"][i]), "-crackSpacing", c.smx, c.smy]
    if SEC.rho[i] > 0:
        a += ["-ties", float(SEC.rho[i]), SEC.tie_steel.fy, SEC.tie_steel.fu]
    return a


def _materials():
    ops.nDMaterial("LadrunoMCFTStrip", 1, *_strip_args(CORE))
    ops.nDMaterial("LadrunoMCFTStrip", 2, *_strip_args(COVER))
    ops.uniaxialMaterial("LadrunoMCFTSteel", 10, SEC.bar_steel.fy, SEC.bar_steel.fu)
    ops.nDMaterial("LadrunoUniaxialFiber", 3, 10)


def _section(tag=1):
    ops.section("LadrunoNDFiber", tag, "-shearProfile", "parabolic", "-shearDepth", SEC.h, SEC.b,
                "-shearResultant", "modelC", "-noCentroid")
    for i in range(SEC.nstrip):
        ops.fiber(float(SEC.y[i]), 0.0, float(SEC.Ac[i]), 1 if DATA["core"][i] else 2)
    for yb, A in zip(SEC.bars_y, SEC.bars_A):
        ops.fiber(float(yb), 0.0, float(A), 3)


def _drive(rows):
    t = [0.0] + [float(k + 1) for k in range(len(rows))]
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    for c in range(6):
        ops.timeSeries("Path", c + 1, "-time", *t, "-values", *([0.0] + [r[c] for r in rows]))
        ops.pattern("Plain", c + 1, c + 1)
        ops.sp(2, DOF_OF[c], 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _path():
    """(eps0, kappa, gamma) cycles through cracking, yield and shear softening."""
    rows = []
    amp = [(4e-6, 0.002), (-4e-6, -0.002), (1e-5, 0.005), (-1e-5, -0.005), (2e-5, 0.01), (-2e-5, -0.01)]
    k0, g0 = 0.0, 0.0
    for ka, ga in amp:
        for j in range(1, 21):
            k, g = k0 + (ka - k0) * j / 20, g0 + (ga - g0) * j / 20
            rows.append((-3e-4 + 0.4 * abs(k) * 228.6, k, 0.0, g, 0.0, 0.0))
        k0, g0 = ka, ga
    return rows


# ----------------------------------------------------------------------------- section vs reference
def test_section_matches_the_reference():
    _model()
    _materials()
    _section()
    rows = _path()
    _drive(rows)
    st = SEC.new_state()
    scale = np.array([SEC.conc.fc * SEC.b * SEC.h, SEC.conc.fc * SEC.b * SEC.h ** 2, SEC.conc.fc * SEC.b * SEC.h])
    worst = np.zeros(3)
    vmax = 0.0
    for _ in rows:
        assert ops.analyze(1) == 0
        d = ops.eleResponse(1, "section", "deformation")
        f = ops.eleResponse(1, "section", "force")
        e = (d[0], d[1], d[3])
        s = SEC.forces(e, st)
        SEC.commit(e, st)
        got = np.array([f[0], f[1], f[3]])
        worst = np.maximum(worst, np.abs(got - s) / scale)
        vmax = max(vmax, abs(s[2]))
    assert worst.max() < TOL, f"N, M, V max rel. diff {worst}"
    assert vmax > 4e5          # the path reaches the shear capacity (a few hundred kN)


# ----------------------------------------------------------------------------- bar law
def test_bar_law_matches_the_reference():
    _model()
    ops.uniaxialMaterial("LadrunoMCFTSteel", 10, 414.0, 621.0, "-esh", 0.01, "-esu", 0.1)
    ops.testUniaxialMaterial(10)
    S = ref.Steel(414.0, 621.0, esh=0.01, esu=0.1)
    ep = 0.0
    path = list(np.linspace(0, 0.04, 30)) + list(np.linspace(0.04, -0.01, 30)) + list(np.linspace(-0.01, 0.11, 40))
    worst = 0.0
    for e in path:
        ops.setStrain(float(e))                      # setStrain commits (quirk), as the reference does
        f = S.stress(e, ep)
        ep = float(e - f / S.Es)
        worst = max(worst, abs(ops.getStress() - float(f)))
    assert worst < 1e-10, worst


# ----------------------------------------------------------------------------- responses
def test_strip_responses():
    _model()
    _materials()
    _section()
    rows = _path()[:60]
    _drive(rows)
    for _ in rows:
        assert ops.analyze(1) == 0
    y = float(SEC.y[20])
    st = ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "mcft")
    assert len(st) == 10
    ey, e1, _e2, w, fst, _f1, _f2, beta, _th, found = st
    assert e1 > 10 * SEC.conc.ecr and w > 0.0 and 0.0 < beta < 1.0 and found == 1.0
    assert ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "eps_y") == pytest.approx([ey], abs=0)
    assert ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "crackWidth") == pytest.approx([w], abs=0)
    assert ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "tieStress") == pytest.approx([fst], abs=0)
    assert len(ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "stress")) == 3
    assert len(ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "tangent")) == 9


def test_elastic_section_stiffness_is_the_uncracked_strip():
    """At a tiny deformation from the virgin state the section is elastic: K_PP = Ec Ac + Es As, and the
    gamma_xz rows carry Gxz = Ec/2 (the strip's own -Gxz default)."""
    _model()
    _materials()
    _section()
    _drive([(1e-7, 0.0, 0.0, 1e-7, 1e-7, 0.0)])
    assert ops.analyze(1) == 0
    f = ops.eleResponse(1, "section", "force")
    Ec = 4700 * math.sqrt(SEC.conc.fc)
    EA = Ec * SEC.Ac.sum() + SEC.bar_steel.Es * SEC.bars_A.sum()
    # Popovics is linear only to O((c/ep)^(r-1)) at the tiny compressive e2 the shear makes: rel 1e-3
    assert f[0] / 1e-7 == pytest.approx(EA, rel=1e-3)
    # Vy (modelC resultant, parabolic): V = sum txy A, txy = (Ec/2) phi(y) gamma
    assert f[3] / 1e-7 == pytest.approx(0.5 * Ec * float((SEC.gy * SEC.Ac).sum()), rel=1e-3)
    # Vz: txz = Gxz phi_z(z = 0) gamma_z, phi_z(0) = 1.5, elastic and exact
    assert f[4] / 1e-7 == pytest.approx(0.5 * Ec * 1.5 * SEC.Ac.sum(), rel=1e-12)


# ----------------------------------------------------------------------------- views
def test_views():
    """BeamFiber2d: a 2D NDFiber section of strips, axial tension below cracking is linear (Ec A)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.nDMaterial("LadrunoMCFTStrip", 1, *_strip_args(CORE))
    ops.section("NDFiber", 1)
    ops.fiber(-100.0, 0.0, 1000.0, 1)
    ops.fiber(100.0, 0.0, 1000.0, 1)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, 1e-5); ops.sp(2, 2, 0.0); ops.sp(2, 3, 0.0)       # axial tension below cracking: linear
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-12, 10, 0); ops.algorithm("Newton"); ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    Ec = 4700 * math.sqrt(SEC.conc.fc)
    assert ops.eleResponse(1, "section", "force")[0] == pytest.approx(Ec * 2000.0 * 1e-5, rel=1e-6)
    # The refusal of other views (getCopy("PlaneStress") returns null) is checked through `ndTest` in
    # test_classic_tcl_ndtest_views: the stock quad element calls exit() on a null copy.


# ----------------------------------------------------------------------------- parsers
GOOD = [25.6, 26.4, 0.0023, 0.025, "-crackSpacing", 152.4, 457.2]


def test_parser_accepts():
    _model()
    ops.nDMaterial("LadrunoMCFTStrip", 1, *GOOD)
    ops.nDMaterial("LadrunoMCFTStrip", 2, *GOOD, "-ties", 0.001, 400.0, 600.0, "-tieSteel", 2e5, 0.008, 0.12,
                   "-agg", 19.0, "-betaC", 170.0, "-cvci", 0.18, "-Gxz", 1.0e4, "-tangentStep", 1e-8)
    ops.uniaxialMaterial("LadrunoMCFTSteel", 3, 414.0, 621.0)
    ops.uniaxialMaterial("LadrunoMCFTSteel", 4, 414.0, 621.0, "-Es", 2e5, "-esh", 0.01, "-esu", 0.1)


@pytest.mark.parametrize("args", [
    (25.6, 26.4, 0.0023),                                         # too few
    (25.6, 26.4, 0.0023, 0.025),                                  # -crackSpacing is required
    (-25.6, 26.4, 0.0023, 0.025, "-crackSpacing", 152.4, 457.2),  # fc <= 0
    (25.6, 26.4, 0.0023, 0.002, "-crackSpacing", 152.4, 457.2),   # e20 <= ep
    (25.6, 26.4, 0.0005, 0.025, "-crackSpacing", 152.4, 457.2),   # fp/ep >= Ec (Popovics)
    (*GOOD[:4], "-crackSpacing", 0.0, 457.2),
    (*GOOD, "-ties", 0.001, 400.0, 300.0),                        # fut < fyt
    (*GOOD, "-ties", -0.001, 400.0, 600.0),
    (*GOOD, "-tieSteel", 2e5, 0.008, 0.12),                       # -tieSteel without -ties
    (*GOOD, "-ties", 0.001, 400.0, 600.0, "-tieSteel", 2e5, 0.001, 0.12),   # esh <= fy/Es
    (*GOOD, "-Gxz", -1.0),
    (*GOOD, "-tangentStep", 0.0),
    (*GOOD, "-cvci", 0.0),
    (*GOOD, "-bogus", 1.0),
    (*GOOD, "-ties", 0.001, 400.0),                               # too few values
])
def test_strip_parser_refuses(args):
    _model()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoMCFTStrip", 1, *args)


@pytest.mark.parametrize("args", [(414.0,), (-414.0, 621.0), (414.0, 300.0), (414.0, 621.0, "-esh", 0.001),
                                  (414.0, 621.0, "-esu", 0.005), (414.0, 621.0, "-bogus", 1.0)])
def test_steel_parser_refuses(args):
    _model()
    with pytest.raises(Exception):
        ops.uniaxialMaterial("LadrunoMCFTSteel", 1, *args)


# ----------------------------------------------------------------------------- robustness
def test_wild_trial_strain_is_bounded():
    """One far-off trial (a Newton iterate can be ~1e4 off the path) costs a fixed amount: the root scan has
    81 candidates and 60 halvings whatever the strain, and the stress stays finite."""
    _model()
    _materials()
    _section()
    _drive([(-3e-4, 1e-6, 0.0, 1e-3, 0.0, 0.0), (5e3, -2e1, 0.0, -3e3, 0.0, 0.0)])
    assert ops.analyze(1) == 0
    t0 = time.perf_counter()
    ops.analyze(1)
    assert time.perf_counter() - t0 < 5.0
    assert all(math.isfinite(v) for v in ops.eleResponse(1, "section", "force"))


# ----------------------------------------------------------------------------- wire format
def test_serialization_round_trip():
    """Parameters and the committed memory (eps_y reference, em1, em2, tie and bar plastic strains, the
    first-commit flag) survive a database round trip, two rounds. The strip is a fixed point at a held
    strain (a trial is a function of the committed memory and the strain), so the restored state re-solves
    to the same forces and strip state bit for bit."""
    rows = _path()[:70]
    y = float(SEC.y[20])

    def build():
        _model()
        _materials()
        _section()
        _drive(rows)
        for _ in rows:
            assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="mcft_wp254", rounds=2,
                       probe_fn=lambda: list(ops.eleResponse(1, "section", "force"))
                       + list(ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "mcft"))
                       + list(ops.eleResponse(1, "section", "fiber", float(SEC.bars_y[2]), 0.0, 3, "stress")))


# ----------------------------------------------------------------------------- classic Tcl
def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()


def _tcl(tmp_path, lines):
    deck = tmp_path / "deck.tcl"
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True, timeout=300,
                          env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL, check=False)
    return proc.stdout + proc.stderr


def _tcl_args(args):
    return " ".join(str(a) for a in args)


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_engine_gives_the_same_section(tmp_path):
    rows = _path()[:40]
    _model()
    _materials()
    _section()
    _drive(rows)
    py = []
    for _ in rows:
        assert ops.analyze(1) == 0
        py.append(list(ops.eleResponse(1, "section", "force")))

    fib = [f"  fiber {float(SEC.y[i])!r} 0.0 {float(SEC.Ac[i])!r} {1 if DATA['core'][i] else 2}" for i in range(SEC.nstrip)]
    fib += [f"  fiber {float(yb)!r} 0.0 {float(A)!r} 3" for yb, A in zip(SEC.bars_y, SEC.bars_A)]
    t = " ".join(str(x) for x in [0.0] + [float(k + 1) for k in range(len(rows))])
    lines = [
        "model basic -ndm 3 -ndf 6",
        f"nDMaterial LadrunoMCFTStrip 1 {_tcl_args(_strip_args(CORE))}",
        f"nDMaterial LadrunoMCFTStrip 2 {_tcl_args(_strip_args(COVER))}",
        f"uniaxialMaterial LadrunoMCFTSteel 10 {SEC.bar_steel.fy!r} {SEC.bar_steel.fu!r}",
        "nDMaterial LadrunoUniaxialFiber 3 10",
        (f"section LadrunoNDFiber 1 -shearProfile parabolic -shearDepth {SEC.h!r} {SEC.b!r} -shearResultant modelC "
         "-noCentroid {"), *fib, "}",
        "node 1 0 0 0", "node 2 0 0 0", "fix 1 1 1 1 1 1 1",
        "element zeroLengthSection 1 1 2 1",
    ]
    for c in range(6):
        v = " ".join(repr(float(x)) for x in [0.0] + [r[c] for r in rows])
        lines += [f"timeSeries Path {c + 1} -time {{{t}}} -values {{{v}}}",
                  f"pattern Plain {c + 1} {c + 1} {{ sp 2 {DOF_OF[c]} 1.0 }}"]
    lines += ["system FullGeneral", "numberer Plain", "constraints Transformation",
              "test NormDispIncr 1.0e-12 30 0", "algorithm Newton", "integrator LoadControl 1.0",
              "analysis Static",
              f"for {{set k 0}} {{$k < {len(rows)}}} {{incr k}} {{",
              "  if {[analyze 1] != 0} { puts \"FAILED $k\"; exit 1 }",
              "  puts \"F [eleResponse 1 section force]\"",
              "}"]
    out = _tcl(tmp_path, lines)
    tcl = [[float(x) for x in l.split()[1:]] for l in out.splitlines() if l.startswith("F ")]
    assert len(tcl) == len(py), out[-2000:]
    assert np.abs(np.array(tcl) - np.array(py)).max() <= 1e-9 * SEC.conc.fc * SEC.b * SEC.h ** 2


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_ndtest_views(tmp_path):
    """`ndTest` (WP-248): the BeamFiber view has order 3 and an uncracked tangent diag(Ec, Ec/2, Gxz); the
    BeamFiber2d view order 2; PlaneStress is refused."""
    args = _tcl_args(_strip_args(CORE))
    lines = [
        "model basic -ndm 3 -ndf 6",
        f"nDMaterial LadrunoMCFTStrip 1 {args} -Gxz 5000.0",
        "model testNDMaterial",
        "puts \"O3 [ndTest material 1 BeamFiber]\"",
        "ndTest strain 1e-6 0.0 0.0",
        "puts \"K [ndTest getTangent]\"",
        "puts \"O2 [ndTest material 1 BeamFiber2d]\"",
        "if {[catch {ndTest material 1 PlaneStress} msg]} { puts \"REFUSED\" } else { puts \"ACCEPTED\" }",
    ]
    out = _tcl(tmp_path, lines)
    got = {l.split()[0]: l.split()[1:] for l in out.splitlines() if l.split() and l.split()[0] in ("O3", "O2", "K")}
    assert got["O3"] == ["3"] and got["O2"] == ["2"], out[-2000:]
    K = np.array([float(x) for x in got["K"]]).reshape(3, 3)
    Ec = 4700 * math.sqrt(SEC.conc.fc)
    assert K[0, 0] == pytest.approx(Ec, rel=1e-6) and K[1, 1] == pytest.approx(0.5 * Ec, rel=1e-6)
    assert K[2, 2] == 5000.0
    assert "REFUSED" in out and "ACCEPTED" not in out
