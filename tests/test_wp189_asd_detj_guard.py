"""WP-189: ASDAbsorbingBoundary3D is independent of the handedness of its node order
and refuses a degenerate brick.

The element integrates its free-field column (addKff / addRff) over the H8 built
from the nodes in the written order. The free-field DOF mapping and the mass use
the coordinate sort, so the only order-dependent term was the signed volume
measure `detJ * gw`: a mirrored (left-handed) write order gave a negative
free-field stiffness and residual (Paper 1 rev3 Tier 3: 161 such skin elements
per cell, unstable model). The volume measure is now |det J|, and setDomain
refuses only a genuinely bad brick: det J changing sign between the 8 gauss
points (folded / inverted) or min |det J| < 1e-6 max |det J| (collapsed).

Checks (classic OpenSees.exe, one-element decks written here):
  * the mirrored write order (n0 n3 n2 n1 n4 n7 n6 n5, the paper generator's
    flip) gives the same element force (eleResponse) and displacements as the
    original order, to round-off, for every boundary type with a vertical face;
  * right-handed decks are byte-identical to LADRUNO_WP189_BASE_EXE (a binary
    built from the base commit) when that variable is set;
  * a collapsed face and a folded (bow-tie) brick are refused with a message
    naming the element tag and nodes; a pure bottom boundary is left alone.

Measured wall time: about 5 s.
"""
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
EXE = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
needs_exe = pytest.mark.skipif(not EXE.exists(), reason="dist OpenSees not built")

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
RIGHT = [1, 2, 3, 4, 5, 6, 7, 8]
LEFT = [5, 6, 7, 8, 1, 2, 3, 4]
MIRROR = [1, 4, 3, 2, 5, 8, 7, 6]
FOLDED = [1, 2, 3, 4, 5, 7, 6, 8]


def deck(order, btype, coords=NODES):
    s = "model BasicBuilder -ndm 3 -ndf 3\n"
    for i, c in enumerate(coords, 1):
        s += "node %d %s\n" % (i, " ".join(str(x) for x in c))
    s += "element ASDAbsorbingBoundary3D 1 %s 1000.0 0.3 2.0 %s\n" % (
        " ".join(map(str, order)), btype)
    s += """for {set i 1} {$i <= 4} {incr i} { fix $i 1 1 1 }
timeSeries Linear 1
pattern Plain 1 1 {
  foreach {n ux uy uz} {5 1e-3 5e-4 -2e-3  6 2e-4 3e-4 -1e-3  7 5e-4 1e-4 -1e-3  8 2e-4 1e-4 -3e-3} {
    sp $n 1 $ux; sp $n 2 $uy; sp $n 3 $uz
  }
}
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1e-10 20
algorithm Linear
integrator LoadControl 1.0
analysis Static
puts "analyze=[analyze 1]"
setParameter -val 1 -ele 1 stage
puts "force [eleResponse 1 force]"
wipeAnalysis
remove loadPattern 1
timeSeries Path 2 -dt 0.01 -values {0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0} -factor 1
pattern Plain 2 2 {
  load 5 1.0 0.5 -2.0
  load 6 0.0 0.3 -1.0
  load 7 0.5 0.0 -1.0
  load 8 0.2 0.1 -3.0
}
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1e-12 30
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient
for {set k 0} {$k < 30} {incr k} {
  set ok [analyze 1 0.01]
  puts "step $k ok=$ok [format %.12e [nodeDisp 5 1]] [format %.12e [nodeDisp 6 2]] [format %.12e [nodeDisp 7 3]] [format %.12e [nodeDisp 8 3]]"
}
foreach n {5 6 7 8} { puts "disp $n [format %.15e [nodeDisp $n 1]] [format %.15e [nodeDisp $n 2]] [format %.15e [nodeDisp $n 3]]" }
"""
    return s


def run(exe, text, tmp_path):
    f = tmp_path / "d.tcl"
    f.write_text(text)
    return subprocess.run([str(exe), str(f)], cwd=tmp_path, capture_output=True, text=True, timeout=120)


@needs_exe
@pytest.mark.parametrize("btype", ["L", "LB", "R", "F", "K", "LF", "LK", "RF", "RK", "LFB"])
def test_mirrored_order_same_response(tmp_path, btype):
    import numpy as np
    out = {}
    for name, order in (("right", RIGHT), ("mirror", MIRROR), ("swapped", LEFT)):
        d = tmp_path / name
        d.mkdir()
        r = run(EXE, deck(order, btype), d)
        text = r.stdout + r.stderr
        assert r.returncode == 0, text
        assert "degenerate" not in text
        disp = np.array([[float(x) for x in l.split()[2:]] for l in text.splitlines() if l.startswith("disp ")])
        K = np.array([float(x) for l in text.splitlines() if l.startswith("force ") for x in l.split()[1:]])
        assert K.size == 24
        hist = np.array([[float(x) for x in l.split()[3:]] for l in text.splitlines() if l.startswith("step ")])
        assert hist.shape == (30, 4)
        # element force is ordered like the written nodes: index it by node id
        byid = np.zeros((9, 3))
        for k, nid in enumerate(order):
            byid[nid] = K[3 * k:3 * k + 3]
        out[name] = (hist, byid)
    assert np.abs(out["right"][1]).max() > 0.0
    for name in ("mirror", "swapped"):
        np.testing.assert_allclose(out[name][1], out["right"][1], rtol=1e-10, atol=1e-6)
        np.testing.assert_allclose(out[name][0], out["right"][0], rtol=1e-8, atol=1e-18)
    assert np.abs(out["right"][0]).max() > 0.0


@needs_exe
def test_collapsed_refused(tmp_path):
    coords = list(NODES)
    coords[4:] = [(0, 0, 0.0), (1, 0, 0.0), (1, 1, 0.0), (0, 1, 0.0)]  # zero height
    r = run(EXE, deck(RIGHT, "L", coords), tmp_path)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    assert "element 1 is a degenerate brick" in out
    assert "1 2 3 4 5 6 7 8" in out


@needs_exe
@pytest.mark.parametrize("btype", ["L", "LB", "R"])
def test_folded_refused(tmp_path, btype):
    r = run(EXE, deck(FOLDED, btype), tmp_path)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    assert "element 1 is a degenerate brick" in out
    assert "1 2 3 4 5 7 6 8" in out
    assert "analyze=" not in out


@needs_exe
def test_pure_bottom_unchanged(tmp_path):
    # a bottom boundary never integrates det J: any order still loads
    for order in (RIGHT, LEFT, MIRROR):
        r = run(EXE, deck(order, "B"), tmp_path)
        assert r.returncode == 0, r.stdout + r.stderr


@needs_exe
@pytest.mark.parametrize("btype", ["L", "LB", "R", "F", "K"])
def test_right_handed_runs(tmp_path, btype):
    r = run(EXE, deck(RIGHT, btype), tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "analyze=" in r.stdout + r.stderr


@needs_exe
@pytest.mark.skipif(not os.environ.get("LADRUNO_WP189_BASE_EXE"), reason="no base-commit binary given")
@pytest.mark.parametrize("btype", ["L", "LB", "B"])
def test_valid_deck_matches_base_binary(tmp_path, btype):
    base = Path(os.environ["LADRUNO_WP189_BASE_EXE"])
    a = run(EXE, deck(RIGHT, btype), tmp_path)
    b = run(base, deck(RIGHT, btype), tmp_path)
    assert a.returncode == b.returncode == 0
    assert a.stdout + a.stderr == b.stdout + b.stderr
