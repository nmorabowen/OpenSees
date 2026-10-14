"""WP-189: ASDAbsorbingBoundary3D refuses a left-handed / degenerate node ordering.

The element integrates its free-field column (addKff / addKffToSoil) with the
node coordinates in the order the user gave, so det J <= 0 at a gauss point
means a negative (or zero) volume weight: a silently negative stiffness. In the
Paper 1 rev3 Tier-3 decks 161 skin elements per cell were left-handed and the
model was unstable. setDomain now prints the element tag, the nodes as given
and the minimum det J, and stops (exit(-1), like the element's other setup
refusals).

Checks (classic OpenSees.exe, one-element decks written here):
  * right-handed element: runs, exit 0;
  * the same element with top and bottom faces swapped (det J = -0.125):
    refused for the vertical boundary types "L" and "LB", message names the
    tag and the nodes, exit != 0 and no `analyze` output;
  * degenerate (collapsed top face): refused;
  * a pure bottom boundary "B" never uses det J, so it is left alone (no change);
  * valid decks are byte-identical to LADRUNO_WP189_BASE_EXE (a binary built
    from the base commit) when that variable is set.

Measured wall time: about 3 s.
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


def deck(order, btype, coords=NODES):
    s = "model BasicBuilder -ndm 3 -ndf 3\n"
    for i, c in enumerate(coords, 1):
        s += "node %d %s\n" % (i, " ".join(str(x) for x in c))
    s += "element ASDAbsorbingBoundary3D 1 %s 1000.0 0.3 2.0 %s\n" % (
        " ".join(map(str, order)), btype)
    s += """for {set i 1} {$i <= 4} {incr i} { fix $i 1 1 1 }
timeSeries Linear 1
pattern Plain 1 1 {
  load 5 1.0 0.5 -2.0
  load 6 0.0 0.3 -1.0
  load 7 0.5 0.0 -1.0
  load 8 0.2 0.1 -3.0
}
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1e-10 20
algorithm Linear
integrator LoadControl 1.0
analysis Static
puts "analyze=[analyze 1]"
foreach n {5 6 7 8} { puts "disp $n [format %.15e [nodeDisp $n 1]] [format %.15e [nodeDisp $n 2]] [format %.15e [nodeDisp $n 3]]" }
"""
    return s


def run(exe, text, tmp_path):
    f = tmp_path / "d.tcl"
    f.write_text(text)
    return subprocess.run([str(exe), str(f)], cwd=tmp_path, capture_output=True, text=True, timeout=120)


@needs_exe
@pytest.mark.parametrize("btype", ["L", "LB", "R", "F"])
def test_right_handed_runs(tmp_path, btype):
    r = run(EXE, deck(RIGHT, btype), tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "analyze=" in r.stdout + r.stderr
    assert "det J" not in (r.stdout + r.stderr)


@needs_exe
@pytest.mark.parametrize("btype", ["L", "LB", "R", "F", "LF"])
def test_left_handed_refused(tmp_path, btype):
    r = run(EXE, deck(LEFT, btype), tmp_path)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    assert "element 1 has det J <= 0" in out
    assert "-0.125" in out
    assert "5 6 7 8 1 2 3 4" in out
    assert "analyze=" not in out


@needs_exe
def test_degenerate_refused(tmp_path):
    coords = list(NODES)
    coords[4:] = [(0, 0, 0.0), (1, 0, 0.0), (1, 1, 0.0), (0, 1, 0.0)]  # zero height
    r = run(EXE, deck(RIGHT, "L", coords), tmp_path)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    assert "det J <= 0" in out


@needs_exe
def test_pure_bottom_unchanged(tmp_path):
    # a bottom boundary never integrates det J: both orderings still load
    for order in (RIGHT, LEFT):
        r = run(EXE, deck(order, "B"), tmp_path)
        assert r.returncode == 0, r.stdout + r.stderr
        assert "det J" not in (r.stdout + r.stderr)


@needs_exe
@pytest.mark.skipif(not os.environ.get("LADRUNO_WP189_BASE_EXE"), reason="no base-commit binary given")
@pytest.mark.parametrize("btype", ["L", "LB", "B"])
def test_valid_deck_matches_base_binary(tmp_path, btype):
    base = Path(os.environ["LADRUNO_WP189_BASE_EXE"])
    a = run(EXE, deck(RIGHT, btype), tmp_path)
    b = run(base, deck(RIGHT, btype), tmp_path)
    assert a.returncode == b.returncode == 0
    assert a.stdout + a.stderr == b.stdout + b.stderr
