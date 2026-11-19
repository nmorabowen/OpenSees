"""WP-239: LadrunoMembrane holds its supporting frames by TAG, not by pointer.

Before WP-239 the membrane resolved its frames once (when they all existed) and
cached a raw `Element*`. `remove element` deletes the beam, so a beam removed
and re-created with the same tag -- the hinge step of a nonlinear model, which
moves each beam onto hinge-twin nodes -- left the membrane writing its slab
share through a freed pointer (undefined behaviour; it "worked" only when the
allocator reused the freed block for the new beam). And a model restored from a
database re-resolved the frames against the twins, which are not membrane
nodes: the frames lost their link lines and printed one "lies on the panel but
is not one of its nodes" warning per end.

Now:
  * addLoad looks every frame up by tag through the domain;
  * a frame end on a node coincident with a membrane external node (a hinge
    twin) ties to that external node, the ORIGINAL one.

The model (tests/_testbed/membrane_twin_model.py) is a 2 x 1 bay floor on
columns with two membranes on seven beams; the replay re-creates every beam on
twins after the membranes exist.
"""
import json
import math
import os
import subprocess
import tempfile

import pytest
from _testbed import membrane_twin_model as M
from _testbed import ops
from _testbed.subprocess_run import pinned_child

pytestmark = [pytest.mark.zone_a]

SHEAR_IDX = (1, 2, 7, 8)   # Vy, Vz at end i and end j of localForce


def _solve(hinge=None, decoy=False, lateral=True):
    M.build(hinge=hinge, decoy=decoy)
    M.apply_loads(lateral=lateral)
    M.analysis()
    assert ops.analyze(1) == 0
    return M.snapshot()


def _close(a, b, rtol, atol=1e-9):
    return abs(a - b) <= rtol * max(abs(a), abs(b)) + atol


def _assert_same(s0, s1, rtol, what):
    for b in M.BEAMS:
        f0, f1 = s0["forces"][b], s1["forces"][b]
        for i in SHEAR_IDX:
            assert _close(f0[i], f1[i], rtol), (
                f"{what}: beam {b} localForce[{i}] {f0[i]} vs {f1[i]}")
    for n in M.SLAB:
        for d in range(6):
            a, c = s0["disp"][n][d], s1["disp"][n][d]
            assert _close(a, c, rtol, 1e-12), f"{what}: node {n} dof {d} {a} vs {c}"


def _sum_rz(snap):
    return sum(r[2] for r in snap["reac"].values())


@pytest.mark.parametrize("decoy", [False, True])
@pytest.mark.parametrize("hinge", ["equalDOF", "zeroLength"])
def test_replayed_frames_keep_their_slab_share(hinge, decoy):
    """Rigid twins: the replayed model must equal the model with no replay,
    beam by beam, and the reactions must carry the whole applied load. The
    decoy beams (allocated between each remove and re-add) carry nothing."""
    ref = _solve(hinge=None)
    got = _solve(hinge=hinge, decoy=decoy)
    _assert_same(ref, got, 1e-6, f"{hinge}/decoy={decoy}")
    total = M.total_applied_vertical()
    assert math.isclose(_sum_rz(got), total, rel_tol=1e-9), (_sum_rz(got), total)
    for n in M.DECOY_NODES:
        if n in got["reac"]:
            assert max(abs(v) for v in got["reac"][n]) < 1e-9, got["reac"][n]


def test_membrane_declared_after_the_replay_ties_to_the_original_nodes(capfd):
    """Resolving against frames that already sit on twins (what a restore
    does) must give the same slab as resolving before the replay, with no
    'not one of its nodes' warning: twin ends map to the coincident
    membrane node."""
    ref = _solve(hinge="zeroLengthPinned")
    capfd.readouterr()

    M.build(hinge="zeroLengthPinned", membranes_after_hinge=True)
    M.apply_loads()
    M.analysis()
    assert ops.analyze(1) == 0
    late = M.snapshot()
    out = capfd.readouterr()
    text = out.out + out.err
    assert "not one of its nodes" not in text, text
    _assert_same(ref, late, 1e-9, "membrane after replay")


def test_removed_frame_is_reported_not_dereferenced(capfd):
    """A frame removed and never re-created: the membrane warns once and
    skips its share (it used to write through the freed pointer)."""
    M.build(hinge=None)
    ops.remove("element", 6)            # the interior beam shared by both bays
    ops.uniaxialMaterial("Elastic", 3, 1.0)   # churn the allocator a little
    M.apply_loads(lateral=False)
    M.analysis()
    assert ops.analyze(1) == 0
    out = capfd.readouterr()
    text = out.out + out.err
    assert "supporting frame 6 is no longer in the domain" in text, text
    ops.reactions()
    rz = sum(ops.nodeReaction(10 + n, 3) for n in M.SLAB)
    assert rz < M.total_applied_vertical() - 1.0   # the share of beam 6 is gone


_CHILD = r"""
import json, sys
from _testbed import ops
from _testbed import membrane_twin_model as M
mode, hinge, db, out = sys.argv[1:5]
hinge = None if hinge == "none" else hinge
res = {}
if mode == "save":
    M.build(hinge=hinge)
    M.apply_loads()
    M.analysis()
    assert ops.analyze(1) == 0
    ops.database("File", db)
    ops.save(1)
    res["raw_disp"] = {n: list(ops.nodeDisp(n)) for n in M.SLAB}
else:
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.database("File", db)
    ops.restore(1)
    res["raw_disp"] = {n: list(ops.nodeDisp(n)) for n in M.SLAB}
    M.analysis()
# re-evaluate at the same load factor: the forces include the element loads
ops.integrator("LoadControl", 0.0)
assert ops.analyze(1) == 0
res.update(M.snapshot())
ops.wipe()
with open(out, "w") as f:
    json.dump(res, f)
print("CHILD-DONE", flush=True)
"""


def _child(mode, hinge, db, out):
    argv, env = pinned_child(["-u", "-c", _CHILD, mode, hinge, db, out])
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    env["PYTHONPATH"] = tests_dir + os.pathsep + env.get("PYTHONPATH", "")
    p = subprocess.run(argv, check=False, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                       errors="replace", timeout=300, env=env, cwd=os.path.dirname(db))
    assert p.returncode == 0 and "CHILD-DONE" in p.stdout, p.stdout[-4000:]
    with open(out) as f:
        return json.load(f), p.stdout


@pytest.mark.parametrize("hinge", ["zeroLengthPinned", "equalDOF"])
def test_save_then_restore_in_a_fresh_process(hinge):
    """database File + save in one process, restore in another: the same
    displacements and beam forces, no 'not one of its nodes' warning, and the
    five Plain patterns with 21 -beamUniform eleLoads come back."""
    with tempfile.TemporaryDirectory(prefix="wp239_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "model")
        live, _ = _child("save", hinge, db, os.path.join(td, "live.json"))
        back, log = _child("restore", hinge, db, os.path.join(td, "back.json"))
    assert "not one of its nodes" not in log, log
    assert "EleLoad" not in log, log
    # the committed state comes back bit for bit
    assert live["raw_disp"] == back["raw_disp"]
    # one more step at the same load factor, in each process: the same model
    # gives the same answer (round-off only; the 1e12 / 1e-6 springs make K
    # ill-conditioned, so the step itself moves u by ~1e-8 relative)
    for n, d in live["disp"].items():
        for k in range(6):
            assert _close(d[k], back["disp"][n][k], 1e-10, 1e-14), (n, k)
    for b, f in live["forces"].items():
        for i, v in enumerate(f):
            assert _close(v, back["forces"][b][i], 1e-8, 1e-8), (b, i, v, back["forces"][b][i])
    sum_live = sum(r[2] for r in live["reac"].values())
    sum_back = sum(r[2] for r in back["reac"].values())
    assert math.isclose(sum_back, M.total_applied_vertical(), rel_tol=1e-9)
    assert math.isclose(sum_live, sum_back, rel_tol=1e-12)
