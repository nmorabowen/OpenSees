"""WP-232i: the B-spline basis through the shipped binaries — rows G1, G1-K4, G7 and G8 of
tests/mpm_bspline_cases.txt (gates copied from the pre-registered table, wp-briefs "232i").

* G1     `mpm debug shape` at >= 1e4 random admissible positions per basis and ndm (every face's
         2h band included): S >= 0, partition of unity, sum grad S, linear reproduction; the
         cardinal values and C1 continuity through `mpm debug shape1d -basis`.
* G1-K4  `mpm debug shape -folded` within 2h of each `fix` face (tests/mpm_bspline_interface.md
         B5, B8): the odd / even reproduction rules and S' >= 0.
* G7     `-basis bspline2` round-trips through `mpm info basis`, the `.ladruno` output (C2-GRID,
         232a §2 as amended) and a mode-A database restart; `mpm particles -file` into a B-spline
         grid is refused (K2); an inadmissible particle is refused as today; the new refusals
         (fewer than 2 cells, fix + wall, stale body-force length) are loud.
* G8     images against the full domain (mode F, bspline2), small and large amplitude.
Classic Tcl: G1's cardinal values and partition of unity, and G7's refusals and `info basis`.

Not in the MPM mutation suite (Ladruno_scripts/mutation_gate.py): these are unit and equality
checks (G8 compares the binary with itself), physics-free for ZERO and SCALE by design.
"""
from __future__ import annotations

import copy
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_bspline_harness as B  # noqa: E402

pytestmark = pytest.mark.zone_a

REPO = Path(__file__).resolve().parents[1]
ROWS = B.parse_rows()
G1 = ROWS["G1"]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _grid(ops, nd, basis, x0=(0.3, -0.2, 0.1), h=0.25, n=(6, 5, 5)):
    ops.wipe()
    ops.model("basic", "-ndm", nd, "-ndf", nd)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e3, 0.3)
    ops.mpm("grid", 232, "-origin", *x0[:nd], "-h", h, "-cells", *n[:nd], "-basis", basis)
    return list(x0[:nd]), h, list(n[:nd])


# --------------------------------------------------------------------------- #
#  G1                                                                          #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("nd", [2, 3])
@pytest.mark.parametrize("basis", ["bspline2", "bspline3"])
def test_g1_basis_properties(basis, nd):
    ops = _ops()
    x0, h, n = _grid(ops, nd, basis)
    rng = np.random.default_rng(232 + nd + (basis == "bspline3"))
    pts = B.random_positions(rng, x0, h, n, nd, 10000 if nd == 2 else 6000, 1000)
    assert len(pts) >= 10000
    w = {"minS": 0.0, "pu": 0.0, "grad": 0.0, "lin": 0.0, "lingrad": 0.0}
    for x in pts:
        out = B.shape(ops, x, nd)
        stride = 2 + nd
        sS, sG, sx, xg = 0.0, np.zeros(nd), np.zeros(nd), np.zeros((nd, nd))
        for k in range(0, len(out), stride):
            c = B.node_coord(int(round(out[k])), x0, h, n, nd)
            S, dS = out[k + 1], np.array(out[k + 2:k + stride])
            w["minS"] = min(w["minS"], S)
            sS += S
            sG += dS
            sx += S * c
            xg += np.outer(c, dS)
        w["pu"] = max(w["pu"], abs(sS - 1.0))
        w["grad"] = max(w["grad"], float(np.max(np.abs(sG))))
        w["lin"] = max(w["lin"], float(np.max(np.abs(sx - x))))
        w["lingrad"] = max(w["lingrad"], float(np.max(np.abs(xg - np.eye(nd)))))
    assert w["minS"] >= 0.0, w
    assert w["pu"] <= G1["tol_pu"], w
    assert w["grad"] <= G1["tol_grad"] / h, w
    assert w["lin"] <= G1["tol_lin"] * h, w
    assert w["lingrad"] <= G1["tol_lingrad"], w


@pytest.mark.parametrize("deg", [2, 3])
def test_g1_cardinal_values_and_c1(deg):
    ops = _ops()
    _grid(ops, 2, f"bspline{deg}")
    h = 0.37
    card = {2: [(0.0, 0.75), (0.5, 0.5), (-0.5, 0.5), (1.0, 0.125), (-1.0, 0.125), (1.5, 0.0)],
            3: [(0.0, 2.0 / 3.0), (1.0, 1.0 / 6.0), (-1.0, 1.0 / 6.0), (2.0, 0.0), (-2.0, 0.0)]}[deg]
    for r, S in card:
        got = B._pl(ops.mpm("debug", "shape1d", h, 0.1, r * h, "-basis", f"bspline{deg}"))
        assert abs(got[0] - S) <= G1["tol_cardinal"], (r, got, S)
    knots = [-1.5, -0.5, 0.5, 1.5] if deg == 2 else [-2.0, -1.0, 0.0, 1.0, 2.0]
    eps = 1e-10 * h
    jump = 0.0
    for k in knots:
        a = B._pl(ops.mpm("debug", "shape1d", h, 0.1, k * h - eps, "-basis", f"bspline{deg}"))[1]
        b = B._pl(ops.mpm("debug", "shape1d", h, 0.1, k * h + eps, "-basis", f"bspline{deg}"))[1]
        jump = max(jump, abs(b - a))
    assert jump <= G1["tol_c1"] / h, jump


# --------------------------------------------------------------------------- #
#  G1-K4                                                                       #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("clamped", [False, True])
@pytest.mark.parametrize("nd", [2, 3])
@pytest.mark.parametrize("basis", ["bspline2", "bspline3"])
def test_g1k4_folded_support(basis, nd, clamped):
    ops = _ops()
    x0, h, n = _grid(ops, nd, basis, x0=(-0.4, 0.15, 0.05), h=0.2)
    for d in range(nd):
        comps = "xyz"[:nd] if clamped else "xyz"[d]
        ops.mpm("fix", B.FACES[2 * d], comps)
        ops.mpm("fix", B.FACES[2 * d + 1], comps)
    rng = np.random.default_rng(99 + nd + clamped)
    pts = B.random_positions(rng, x0, h, n, nd, 0, 1500 if nd == 2 else 800)
    lo = np.array(x0)
    hi = lo + np.array(n) * h
    w = {"neg": 0.0, "odd": 0.0, "pu": 0.0, "tang": 0.0}
    n_odd = n_even = 0
    for x in pts:
        # faces folding this particle: its ghost carries weight (quadratic: closer than h/2 to
        # the face; cubic: closer than h)
        reach = 0.5 * h if basis == "bspline2" else 1.0 * h
        side = [(-1 if x[d] - lo[d] < reach - 1e-12 else (1 if hi[d] - x[d] < reach - 1e-12 else 0))
                for d in range(nd)]
        out = B.shape(ops, x, nd, folded=True)
        stride = 1 + nd * (1 + nd)
        nodes = []
        for k in range(0, len(out), stride):
            c = B.node_coord(int(round(out[k])), x0, h, n, nd)
            S = [out[k + 1 + a * (1 + nd)] for a in range(nd)]
            nodes.append((c, S))
        for a in range(nd):
            odd_axes = [d for d in range(nd) if side[d] != 0 and (clamped or d == a)]
            vals = [S[a] for _, S in nodes]
            w["neg"] = min(w["neg"], min(vals))
            if len(odd_axes) == 1:
                d = odd_axes[0]
                xf = lo[d] if side[d] < 0 else hi[d]
                s = sum(S[a] * (c[d] - xf) for c, S in nodes)
                w["odd"] = max(w["odd"], abs(s - (x[d] - xf)))
                n_odd += 1
            elif not odd_axes:
                w["pu"] = max(w["pu"], abs(sum(vals) - 1.0))
                for d in range(nd):
                    if side[d] == 0:
                        w["tang"] = max(w["tang"], abs(sum(S[a] * c[d] for c, S in nodes) - x[d]))
                n_even += 1
    assert n_odd > 0 and (clamped or n_even > 0)
    assert w["neg"] >= 0.0, w
    assert w["odd"] <= G1["tol_lin"] * h, w
    assert w["pu"] <= G1["tol_pu"], w
    assert w["tang"] <= G1["tol_lin"] * h, w


# --------------------------------------------------------------------------- #
#  G7                                                                          #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("basis,code", [("gimp", 0), ("linear", 1), ("bspline2", 2), ("bspline3", 3)])
def test_g7_info_basis(basis, code):
    ops = _ops()
    _grid(ops, 2, basis)
    assert int(ops.mpm("info", "basis")) == code


def test_g7_refuses_a_bspline_grid_with_one_cell():
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    with pytest.raises(Exception):
        ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 30, 1, "-basis", "bspline2")
    # gimp keeps its one-cell strips
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 30, 1, "-basis", "gimp")


def test_g7_ladruno_output_round_trip(tmp_path):
    h5py = pytest.importorskip("h5py")
    sys.path.insert(0, str(REPO / "Ladruno_scripts" / "ladruno_recorder_tests"))
    import mpm_contracts as mc
    ops = _ops()
    _grid(ops, 2, "bspline2", x0=(0.0, 0.0), h=0.25, n=(4, 4))
    ops.mpm("body", 1, "box", 0.25, 0.25, 0.75, 0.75, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("fix", "ymin", "y")
    ops.mpm("velocity", 1, "uniform", 0.1, -0.05)
    ops.mpm("units", "m", "N", "s")
    out = tmp_path / "bs.ladruno"
    ops.mpm("output", str(out), "-every", 2, "-fields", "POSITION", "STRESS")
    ops.mpm("analyze", 6, 1.0e-3)
    ops.mpm("output", "-close")
    assert mc.validate_particle_results(str(out)) == []
    with h5py.File(out, "r") as f:
        b = f["MODEL_STAGE[0]/MODEL/PARTICLES/GRID"].attrs["BASIS"]
        assert (b.decode() if isinstance(b, bytes) else str(b)) == "bspline2"


def test_g7_contract1_set_refused_into_a_bspline_grid(capfd):
    pytest.importorskip("h5py")
    fixture = REPO / "tests" / "data" / "mpm" / "particle_set_v1.h5"
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("ElasticIsotropic", 2, 2.0e4, 0.3)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "bspline2")
    with pytest.raises(Exception):
        ops.mpm("particles", "-file", str(fixture))
    err = capfd.readouterr().err
    assert "differs from the deck's 'bspline2'" in err, err
    assert int(ops.mpm("info", "numParticles")) == 0


def test_g7_inadmissible_particle_refused():
    ops = _ops()
    _grid(ops, 2, "bspline2", x0=(0.0, 0.0), h=0.25, n=(4, 4))
    ops.mpm("particle", 0, 1, 1, 1.0, 0.01, "-x", 1.2, 0.5, "-lp", 0.05, 0.05)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)


def test_g7_fix_and_wall_on_one_face_refused():
    ops = _ops()
    _grid(ops, 2, "bspline2", x0=(0.0, 0.0), h=0.25, n=(4, 4))
    ops.mpm("body", 1, "box", 0.25, 0.25, 0.75, 0.75, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("fix", "ymin", "y")
    ops.mpm("wall", "ymin", 0.3)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)


def test_g7_body_force_lifecycle():
    ops = _ops()
    _grid(ops, 2, "bspline2", x0=(0.0, 0.0), h=0.25, n=(4, 4))
    ops.mpm("body", 1, "box", 0.25, 0.25, 0.75, 0.75, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    N = int(ops.mpm("info", "numParticles"))
    with pytest.raises(Exception):
        ops.mpm("particleBodyForce", *([0.0] * (2 * N - 1)))       # wrong length
    with pytest.raises(Exception):
        ops.mpm("particleVelocity", *([0.0] * 3))
    ops.mpm("particleBodyForce", *([0.0, -9.81] * N))
    ops.mpm("analyze", 10, 1.0e-3)
    v_bf = B._pl(ops.mpm("particleData", "v"))
    ops.mpm("reset")                      # clears the body force
    ops.mpm("analyze", 10, 1.0e-3)
    assert max(abs(v) for v in B._pl(ops.mpm("particleData", "v"))) == 0.0
    ops.mpm("reset")
    ops.mpm("gravity", 0.0, -9.81)        # the same loading through gravity: bit for bit (B11)
    ops.mpm("analyze", 10, 1.0e-3)
    assert B._pl(ops.mpm("particleData", "v")) == v_bf
    with pytest.raises(Exception):        # a velocity is an initial condition
        ops.mpm("particleVelocity", *([0.0] * (2 * N)))


def test_g7_mode_a_restart_round_trip():
    """Mode-A database restart on a B-spline grid with rollers (images): the restored run equals
    the uninterrupted one (rel 1e-10), and the restored host still says bspline2. The 232f
    restart protocol (tests/ladruno_mpm_implicit_harness.py) on its bar, two cells high (a
    B-spline grid needs >= 2 cells per axis), 60 steps, saved at 30."""
    import ladruno_mpm_implicit_harness as H2
    ops = _ops()

    def bspline(c):
        c = copy.deepcopy(c)
        c["basis"] = "bspline2"
        c["cells"] = [c["cells"][0], 2]
        for b in c["bodies"]:
            b[5] = "2.0"                   # body <id> box x0 y0 x1 y1 ...
        c["steps"] = 60
        return c

    ref = bspline(H2.case_by_name("reseed_bar_normal"))
    rr = H2.run_python(ref, ops)
    assert rr.get("stopped") is None
    c = bspline(H2.case_by_name("restart_bar"))
    c["restart"] = 30
    r = H2.run_python(c, ops)
    assert r.get("restart_mode") == "saved" and r.get("stopped") is None
    assert int(ops.mpm("info", "basis")) == 2
    worst = max(H2.field_relerr(r, rr, f) for f in H2.PROTOCOL_FIELDS)
    assert worst <= 1e-10, worst


# --------------------------------------------------------------------------- #
#  G8                                                                          #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("amp", ["small", "large"])
def test_g8_images_equal_the_full_domain(amp):
    row = ROWS["G8"]
    V = float(row["v_small"] if amp == "small" else row["v_large"])
    rx, rs, strain = B.g8(_ops(), V, int(row["steps"]))
    assert rx <= float(row["tol"]) and rs <= float(row["tol"]), (rx, rs)
    if amp == "small":
        assert strain < 0.01, strain
    else:
        assert strain > 0.05, strain


# --------------------------------------------------------------------------- #
#  Classic Tcl                                                                 #
# --------------------------------------------------------------------------- #
def test_classic_tcl_g1_and_g7():
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    deck = """
model basic -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 1000.0 0.3
if {[catch {mpm grid 232 -origin 0 0 -h 1.0 -cells 30 1 -basis bspline2}]} { puts "REFUSED_ONECELL 1" }
mpm grid 232 -origin 0.3 -0.2 -h 0.25 -cells 6 5 -basis bspline2
puts "BASIS [mpm info basis]"
puts "CARD [mpm debug shape1d 0.37 0.1 0.0 -basis bspline2] [mpm debug shape1d 0.37 0.1 0.37 -basis bspline2]"
puts "CARD3 [mpm debug shape1d 0.37 0.1 0.0 -basis bspline3] [mpm debug shape1d 0.37 0.1 0.37 -basis bspline3]"
puts "SHAPE [mpm debug shape 0.31 -0.17 -lp 0.1 0.1]"
mpm body 1 box 0.3 -0.2 1.8 1.05 -ppc 2 -mat 1 -rho 1.0
mpm fix ymin y
mpm wall ymin 0.2
if {[catch {mpm analyze 1 1.0e-4}]} { puts "REFUSED_FIXWALL 1" }
set n [mpm info numParticles]
if {[catch {mpm particleBodyForce 1.0 2.0 3.0}]} { puts "REFUSED_BFLEN 1" }
set b {}
for {set i 0} {$i < $n} {incr i} { lappend b 0.0 -1.0 }
mpm particleBodyForce {*}$b
puts "BF_OK 1"
puts "DECK-DONE"
"""
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "g1g7.tcl"
        p.write_text(deck, encoding="utf-8")
        proc = subprocess.run([exe, str(p)], capture_output=True, text=True, timeout=300, cwd=td,
                              stdin=subprocess.DEVNULL,
                              env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"))
    out = proc.stdout + proc.stderr
    assert "DECK-DONE" in out, out[-3000:]
    vals = {}
    for line in out.splitlines():
        t = line.split()
        if t and t[0] in ("BASIS", "CARD", "CARD3", "SHAPE", "REFUSED_ONECELL", "REFUSED_FIXWALL",
                          "REFUSED_BFLEN", "BF_OK"):
            vals[t[0]] = [float(v) for v in t[1:]]
    assert vals["BASIS"] == [2.0]
    c = vals["CARD"]
    assert abs(c[0] - 0.75) <= 1e-14 and abs(c[2] - 0.125) <= 1e-14
    c3 = vals["CARD3"]
    assert abs(c3[0] - 2.0 / 3.0) <= 1e-14 and abs(c3[2] - 1.0 / 6.0) <= 1e-14
    sh = vals["SHAPE"]
    assert abs(sum(sh[k + 1] for k in range(0, len(sh), 4)) - 1.0) <= 1e-13
    for key in ("REFUSED_ONECELL", "REFUSED_FIXWALL", "REFUSED_BFLEN", "BF_OK"):
        assert key in vals, (key, out[-3000:])
