"""WP-197: `system UmfPack` re-analysis after a DOF-graph change (fork issue #981).

ADR: Ladruno_implementation/197_ladruno_umfpack_setsize_adr.md.

UmfpackGenLinSOE::setSize used to `reserve` and then `push_back` into Ap/Ai
without clearing them, so a second setSize after the DOF graph changed appended
the new CSC pattern to the old one. addA then searched the OLD column ranges
("an element entry has no slot in the sparsity pattern; its contribution was
DISCARDED"), symbolic factorization could reject the mixed pattern (every later
solve: "setSize has not been called"), and a pattern that shrank left addA
writing past the end of Ax.

The deck is a 2-storey, 1-bay elastic frame (classic OpenSees.exe), Newton +
NormUnbalance, two LoadControl steps, a graph change, two more steps:

  * sp_add      -- `sp` on a free DOF after the first analyze (constraints
                   Transformation): n and nnz shrink (the reported case);
  * elem_add    -- a diagonal truss between two unconnected free nodes: n is the
                   same, nnz grows (the DISCARDED case);
  * elem_remove -- the diagonal removed: n is the same, nnz shrinks (the Ax
                   overrun case);
  * none        -- no change, and `domainChange` forcing a setSize on the same
                   graph: the two must print the same bytes.

Each changed-graph run must print no DISCARDED / setSize / symbolic warning,
converge, and give the same displacements as `system FullGeneral` on the same
deck and as UmfPack after `wipeAnalysis` (a fresh SOE). With
LADRUNO_WP197_BASE_DIST pointing at a dist/ built from the base commit, the
paths that already worked (unchanged graph, with or without `domainChange`, and the
wipeAnalysis workaround) are byte-compared with it.

Measured wall time: about 3 s for the module (each run ~0.1 s).
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DIST = Path(os.environ.get("LADRUNO_DIST", REPO / "dist"))
EXE_NAME = "OpenSees.exe" if os.name == "nt" else "OpenSees"
SER_EXE = DIST / "bin" / EXE_NAME
needs_serial = pytest.mark.skipif(not SER_EXE.exists(), reason="dist OpenSees not built")

BAD = re.compile(r"DISCARDED|setSize has not been called|symbolic analysis returns|"
                 r"numeric analysis returns|solving returns")

DECK = r"""
wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 6.0 0.0
node 3 0.0 3.0
node 4 6.0 3.0
node 5 0.0 6.0
node 6 6.0 6.0
fix 1 1 1 1
fix 2 1 1 1
geomTransf Linear 1
set E 2.5e7
set A 0.09
set I 6.75e-4
element elasticBeamColumn 1 1 3 $A $E $I 1
element elasticBeamColumn 2 2 4 $A $E $I 1
element elasticBeamColumn 3 3 5 $A $E $I 1
element elasticBeamColumn 4 4 6 $A $E $I 1
element elasticBeamColumn 5 3 4 $A $E $I 1
element elasticBeamColumn 6 5 6 $A $E $I 1
uniaxialMaterial Elastic 1 2.0e8
if {"@CASE@" == "elem_remove"} { element truss 7 3 6 0.002 1 }
pattern Plain 1 Linear {
    load 3 40.0 -10.0 0.0
    load 5 80.0 -10.0 5.0
}
proc mkAnalysis {} {
    constraints Transformation
    numberer RCM
    system @SYSTEM@
    test NormUnbalance 1.0e-9 25
    algorithm Newton
    integrator LoadControl 0.5
    analysis Static
}
proc report {tag} {
    set s $tag
    foreach n {3 4 5 6} {
        foreach d {1 2 3} { append s " " [format %.17g [nodeDisp $n $d]] }
    }
    puts $s
}
mkAnalysis
puts "OK1 [analyze 2]"
report D1
switch -- "@CASE@" {
    sp_add      { pattern Plain 2 Linear { sp 4 1 0.002 } }
    elem_add    { element truss 7 3 6 0.002 1 }
    elem_remove { remove element 7 }
    none        { if {@DC@} { domainChange } }
}
if {@WIPE@} { wipeAnalysis; mkAnalysis }
puts "OK2 [analyze 2]"
report D2
puts "DONE"
"""


def _run(tmp, case, system="UmfPack", wipe=0, dc=0, dist=DIST):
    tmp = Path(tmp)
    name = f"{case}_{system.split()[0]}_w{wipe}_dc{dc}"
    deck = tmp / (name + ".tcl")
    deck.write_text(DECK.replace("@CASE@", case).replace("@SYSTEM@", system)
                    .replace("@WIPE@", str(wipe)).replace("@DC@", str(dc)))
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    tcl = Path(dist) / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    p = subprocess.run([str(Path(dist) / "bin" / EXE_NAME), str(deck)],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=120, env=env, check=False,
                       cwd=str(tmp))
    log = p.stdout + p.stderr
    # classic OpenSees.exe sends Tcl `puts` to stderr; parse both streams
    out = {"rc": p.returncode, "log": log}
    for line in log.splitlines():
        tok = line.split()
        if tok and tok[0] in ("OK1", "OK2"):
            out[tok[0]] = int(tok[1])
        elif tok and tok[0] in ("D1", "D2"):
            out[tok[0]] = [float(x) for x in tok[1:]]
            out[tok[0] + "_text"] = line
    return out


def _close(a, b, rtol=1e-9):
    scale = max(abs(x) for x in b)
    return max(abs(x - y) for x, y in zip(a, b)) <= rtol * scale


def _assert_clean(r, what):
    assert r["rc"] == 0 and re.search(r"^DONE$", r["log"], re.MULTILINE), f"{what}: run died (rc={r['rc']})\n{r['log'][-3000:]}"
    m = BAD.search(r["log"])
    assert m is None, f"{what}: solver warning {m.group(0)!r}\n{r['log'][-3000:]}"
    assert r.get("OK1") == 0 and r.get("OK2") == 0, f"{what}: analyze failed\n{r['log'][-3000:]}"


@needs_serial
@pytest.mark.parametrize("case", ["sp_add", "elem_add", "elem_remove"])
def test_graph_change_reanalysis_matches_fresh_and_dense(tmp_path, case):
    umf = _run(tmp_path, case)
    _assert_clean(umf, f"UmfPack {case}")
    dense = _run(tmp_path, case, system="FullGeneral")
    _assert_clean(dense, f"FullGeneral {case}")
    fresh = _run(tmp_path, case, wipe=1)
    _assert_clean(fresh, f"UmfPack+wipeAnalysis {case}")
    # the change really changes the answer (otherwise the test proves nothing)
    assert not _close(dense["D2"], [2.0 * x for x in dense["D1"]], 1e-6)
    assert _close(umf["D1"], dense["D1"])
    assert _close(umf["D2"], dense["D2"]), (umf["D2"], dense["D2"])
    assert _close(umf["D2"], fresh["D2"]), (umf["D2"], fresh["D2"])


@needs_serial
def test_sp_add_imposes_the_displacement(tmp_path):
    r = _run(tmp_path, "sp_add")
    _assert_clean(r, "UmfPack sp_add")
    # node 4 dof 1 is the 4th entry; pattern 2 is Linear in time, t = 2.0 at the end
    assert abs(r["D2"][3] - 0.002 * 2.0) < 1e-14


@needs_serial
def test_unchanged_graph_resetsize_is_byte_identical(tmp_path):
    # `domainChange` forces a second setSize on the same graph. umfpack_*_symbolic
    # reads Ax; before WP-197 Ax still held the last assembled tangent there (resize
    # kept it), while the first setSize saw zeros, so the forced re-setSize differed
    # from the plain run in the last bits. With Ax zeroed on every setSize the two
    # runs print the same bytes.
    plain = _run(tmp_path, "none", dc=0)
    forced = _run(tmp_path, "none", dc=1)
    _assert_clean(plain, "none dc=0")
    _assert_clean(forced, "none dc=1")
    assert plain["D1_text"] == forced["D1_text"]
    assert plain["D2_text"] == forced["D2_text"]


BASE = os.environ.get("LADRUNO_WP197_BASE_DIST")


@needs_serial
@pytest.mark.skipif(not BASE, reason="LADRUNO_WP197_BASE_DIST not set")
@pytest.mark.parametrize("case,wipe,dc", [("none", 0, 0), ("sp_add", 1, 0),
                                          ("elem_add", 1, 0), ("elem_remove", 1, 0)])
def test_paths_that_worked_are_byte_identical_to_base(tmp_path, case, wipe, dc):
    """Unchanged-graph re-analysis and the wipeAnalysis workaround print the same
    bytes as a dist/ built from the base. (A FORCED re-setSize on an unchanged graph
    is not in the list: it changes in the last bits, to equal the plain run; see
    test_unchanged_graph_resetsize_is_byte_identical.)"""
    (tmp_path / "new").mkdir()
    (tmp_path / "old").mkdir()
    new = _run(tmp_path / "new", case, wipe=wipe, dc=dc)
    old = _run(tmp_path / "old", case, wipe=wipe, dc=dc, dist=Path(BASE))
    _assert_clean(new, "new")
    _assert_clean(old, "base")
    assert new["D1_text"] == old["D1_text"] and new["D2_text"] == old["D2_text"]
