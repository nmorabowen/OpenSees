"""WP-258 -- soil materials on the WP-107 / WP-192 threaded FE element loops.

DruckerPrager, LadrunoSANISAND (ManzariDafalias) and LadrunoHS (HSSmall) answer
`ladrunoThreadSafeUpdate()` / `ladrunoThreadSafeForm()` true for the options
the WP-258 audit proved (Ladruno_implementation/258_thread_safe_soil_materials.md
section 4). This file is the determinism gate for the FE half: the threaded
result must be BIT-IDENTICAL to one thread, at 2, 4 and 8 threads, and the
elements must really run on the threads.

THE DECK. A self-weight footing problem, after the WP-112 strip
(test_ladruno_sanisand_flip_determinism.py): gravity plus a surcharge ramped at
stage 0, `loadConst`, `updateMaterialStage 1` (which reaches ONE element: the
rest pick the stage up lazily, inside the threaded loops -- the HS tag registry
and the SANISAND class-wide flag are both exercised there), then a rigid
footing pushed by `DisplacementControl` with `FixedNumIter` (a determinism gate
needs identical arithmetic, not convergence). Two meshes:

  * 2D: `LadrunoQuad -formulation bbar -type PlaneStrain` (the update loop
    threads; LadrunoQuad does not opt into the WP-192 form loops);
  * 3D: `stdBrick` (update, tangent, residual and commit loops all thread).

Elements are created in a SHUFFLED tag order (a fixed permutation), so the
dynamic schedule hands each thread a different spatial mix than the 2D/3D
natural order does. The digest is every node's displacement after every push
step, by `float.hex`, plus the footing load factor; `system BandGeneral`, so
the solve is serial and deterministic.

Each configuration runs in a child process per thread count (thread state and
the once-per-process latches stay isolated).
"""
import json
import os
import sys

import pytest

_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops                                  # noqa: E402
from _testbed.subprocess_run import run_python_script     # noqa: E402

_TESTS = os.path.dirname(os.path.abspath(__file__))

_OPENMP_PROBE = (
    "import sys; sys.path.insert(0, %r)\n"
    "import opensees as ops\n"
    "print('LADRUNO_PROBE_HW', ops.ladrunoThreads(1 << 20))\n"
)
_OPENMP_ABSENT_MARKER = "built WITHOUT LADRUNO_OPENMP"


def _probe():
    try:
        code, out = run_python_script(_OPENMP_PROBE % _DIST_BIN)
    except Exception:
        return True, 0
    if code != 0 or "LADRUNO_PROBE_HW" not in out:
        return True, 0
    if _OPENMP_ABSENT_MARKER in out:
        return False, 0
    for line in out.splitlines():
        if line.startswith("LADRUNO_PROBE_HW"):
            try:
                return True, int(line.split()[1])
            except (IndexError, ValueError):
                return True, 0
    return True, 0


_OPENMP_IN, _HW = _probe()

pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(not _OPENMP_IN, reason="binary built with LADRUNO_OPENMP=OFF"),
]

SANISAND = [264.32, 0.3129, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3,
            0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
HS = [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0, "-nu", 0.2, "-pref", 100.0,
      "-Rf", 0.9, "-OCR", 1.5]

# name -> (nDMaterial args after the tag, stage-1 id or None). The SANISAND tail is
# `IntScheme TanType JacoType TolF TolR`.
MATERIALS = {
    # weak enough to yield under the footing (secant ratio ~0.7 over the push)
    "dp": (["DruckerPrager", 20000.0, 9000.0, 1.0, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            2.0], None),
    "sanisand_is1": (["LadrunoSANISAND", *SANISAND, 1, 2, 1, 1e-7, 1e-7,
                      "-Presidual", 0.0, "-Pmin", 1.0e-4 * 101.0], 1),
    "sanisand_is2": (["LadrunoSANISAND", *SANISAND, 2, 2, 1, 1e-7, 1e-7,
                      "-Presidual", 0.0, "-Pmin", 1.0e-4 * 101.0], 1),
    "sanisand_is129": (["LadrunoSANISAND", *SANISAND, 129, 2, 1, 1e-7, 1e-7,
                        "-Presidual", 0.0, "-Pmin", 1.0e-4 * 101.0], 1),
    "sanisand_implex": (["LadrunoSANISAND", *SANISAND, 2, 2, 1, 1e-7, 1e-7,
                         "-Presidual", 0.0, "-Pmin", 1.0e-4 * 101.0, "-implex",
                         "-maxSubsteps", 1000], 1),
    # vanilla ManzariDafalias (no -Pmin; the deck keeps p > 0), IntScheme 1. Vanilla IntScheme 2
    # has no substep cap (`-maxSubsteps` is a SANISAND option) and one serial run of this deck
    # exceeds the child timeout -- the uncapped-substepping quirk, not a threading one; the base
    # class's IntScheme-2 path is gated through LadrunoSANISAND (sanisand_is2 / _implex).
    "manzari_is1": (["ManzariDafalias", *SANISAND, 1, 2, 1, 1e-7, 1e-7], 1),
    "hs": (["LadrunoHS", *HS, "-stage", 0], 1),
    "hs_small": (["LadrunoHS", *HS, "-stage", 0, "-smallStrain", 1.0e5, 1.0e-4], 1),
    "hs_implex": (["LadrunoHS", *HS, "-stage", 0, "-implex"], 1),
    "hs_sigmat": (["LadrunoHS", *HS, "-stage", 0, "-sigmat", 5.0], 1),
}

_CHILD = r'''
import json, random, sys
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[2])
from _testbed import ops
cfg = json.loads(sys.argv[3])
mat, stage1, dim, nthreads = cfg["mat"], cfg["stage1"], cfg["dim"], cfg["threads"]
n_push, iters = cfg["n_push"], cfg["iters"]
GAMMA, Q_SUR, NG = 18.0, 20.0, 8
W, H, B, DS = 6.0, 3.0, 1.0, 2.0e-4

ops.wipe()
if dim == 2:
    NX, NZ = 12, 6
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    dx, dz = W / NX, H / NZ
    nid = lambda i, k: 1 + i + (NX + 1) * k
    for k in range(NZ + 1):
        for i in range(NX + 1):
            ops.node(nid(i, k), -0.5 * W + i * dx, -H + k * dz)
    ops.nDMaterial(mat[0], 1, *mat[1:])
    cells = [(i, k) for k in range(NZ) for i in range(NX)]
    order = list(range(len(cells)))
    random.Random(258).shuffle(order)
    for tag0, c in enumerate(order):
        i, k = cells[c]
        ops.element("LadrunoQuad", tag0 + 1, nid(i, k), nid(i + 1, k), nid(i + 1, k + 1),
                    nid(i, k + 1), 1, "-thick", 1.0, "-type", "PlaneStrain",
                    "-formulation", "bbar")
    for k in range(1, NZ + 1):
        ops.fix(nid(0, k), 1, 0)
        ops.fix(nid(NX, k), 1, 0)
    for i in range(NX + 1):
        ops.fix(nid(i, 0), 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    fz = GAMMA * dx * dz / 4.0
    for (i, k) in cells:
        for (a, c) in ((i, k), (i + 1, k), (i + 1, k + 1), (i, k + 1)):
            ops.load(nid(a, c), 0.0, -fz)
    for i in range(NX + 1):
        ops.load(nid(i, NZ), 0.0, -Q_SUR * (dx if 0 < i < NX else 0.5 * dx))
    foot = [nid(i, NZ) for i in range(NX + 1) if abs(-0.5 * W + i * dx) <= 0.5 * B + 1e-9]
    master, vdof = nid(NX // 2, NZ), 2
    nodes = [nid(i, k) for k in range(NZ + 1) for i in range(NX + 1)]
else:
    NX, NY, NZ = 8, 3, 4
    D = 1.0
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx, dy, dz = W / NX, D / NY, H / NZ
    nid = lambda i, j, k: 1 + i + (NX + 1) * (j + (NY + 1) * k)
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                ops.node(nid(i, j, k), -0.5 * W + i * dx, j * dy, -H + k * dz)
    ops.nDMaterial(mat[0], 1, *mat[1:])
    cells = [(i, j, k) for k in range(NZ) for j in range(NY) for i in range(NX)]
    order = list(range(len(cells)))
    random.Random(258).shuffle(order)
    for tag0, c in enumerate(order):
        i, j, k = cells[c]
        ops.element("stdBrick", tag0 + 1,
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                    nid(i, j + 1, k + 1), 1)
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                n = nid(i, j, k)
                if k == 0:
                    ops.fix(n, 1, 1, 1)
                else:
                    ops.fix(n, 1 if i in (0, NX) else 0, 1 if j in (0, NY) else 0, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    fz = GAMMA * dx * dy * dz / 8.0
    for (i, j, k) in cells:
        for a in (i, i + 1):
            for b in (j, j + 1):
                for c in (k, k + 1):
                    ops.load(nid(a, b, c), 0.0, 0.0, -fz)
    for j in range(NY + 1):
        wy = dy if 0 < j < NY else 0.5 * dy
        for i in range(NX + 1):
            wx = dx if 0 < i < NX else 0.5 * dx
            ops.load(nid(i, j, NZ), 0.0, 0.0, -Q_SUR * wx * wy)
    foot = [nid(i, j, NZ) for j in range(NY + 1) for i in range(NX + 1)
            if abs(-0.5 * W + i * dx) <= 0.5 * B + 1e-9]
    master, vdof = nid(NX // 2, 0, NZ), 3
    nodes = [nid(i, j, k) for k in range(NZ + 1) for j in range(NY + 1)
             for i in range(NX + 1)]

ops.ladrunoThreads(nthreads)
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1.0e-8, 60, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0 / NG)
ops.analysis("Static")
if stage1 is not None:      # DruckerPrager has no updateMaterialStage (returns -1)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
rcs = []
for s in range(NG):
    rcs.append(ops.analyze(1))
ops.loadConst("-time", 0.0)
if stage1 is not None:
    ops.updateMaterialStage("-material", 1, "-stage", stage1)
ops.timeSeries("Linear", 2)
ops.pattern("Plain", 2, 2)
for n in foot:
    if n != master:
        ops.equalDOF(master, n, vdof)
ops.load(master, *([0.0] * (vdof - 1)), -1.0)
ops.test("FixedNumIter", iters, 0)
ops.integrator("DisplacementControl", master, vdof, -DS)
steps = []
for s in range(n_push):
    r = ops.analyze(1)
    q = ops.getLoadFactor(2)
    steps.append([r, q.hex(), [x.hex() for n in nodes for x in ops.nodeDisp(n)]])
ledger = None
if "-implex" in mat:   # the process-wide IMPL-EX ledger, read after the threaded commits
    rd = lambda what: [float(x) for x in ops.eleResponse(1, "material", 1, what)]
    ledger = {"max": [x.hex() for x in rd("implexError")], "avg": rd("avgImplexError"),
              "refusals": rd("implexRefusals")}
ops.ladrunoThreads(1)
print("RESULT " + json.dumps({"build": ops.ladrunoBuild(), "gravity_rc": rcs, "steps": steps,
                              "q": [float.fromhex(s[1]) for s in steps], "ledger": ledger}))
'''

_THREADS = [n for n in (1, 2, 4, 8) if _HW == 0 or n <= max(_HW, 1)]


def _run(mat, dim, nthreads, n_push=6, iters=12):
    args, stage1 = MATERIALS[mat]
    cfg = {"mat": args, "stage1": stage1, "dim": dim, "threads": nthreads,
           "n_push": n_push, "iters": iters}
    code, out = run_python_script(_CHILD, argv=[_DIST_BIN, _TESTS, json.dumps(cfg)])
    assert code == 0, out[-3000:]
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    assert res, out[-3000:]
    return json.loads(res[-1][len("RESULT "):]), out


def _threaded_text(dim, n, nele):
    """The loops that must have announced themselves on n threads."""
    want = ["element update loop THREADED on %d threads" % n, "(%d elements)" % nele]
    if dim == 3:
        want.append("%d of %d FE_Elements computed in parallel" % (nele, nele))
    return want


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("mat", sorted(MATERIALS))
def test_threaded_equal_serial(mat, dim):
    if len(_THREADS) < 2:
        pytest.skip("single hardware thread")
    nele = 72 if dim == 2 else 96
    base, _ = _run(mat, dim, 1)
    assert all(r == 0 for r in base["gravity_rc"]), (mat, dim, base["gravity_rc"])
    assert len(base["steps"]) > 0
    # the gate means something only if the push ran and was plastic (critic of record): every
    # push step returned 0, and the footing's secant stiffness fell by 5 % or more over the push
    assert all(st[0] == 0 for st in base["steps"]), (mat, dim, [st[0] for st in base["steps"]])
    q = base["q"]                       # the load factor starts the push at 0
    ratio = (q[-1] - q[-2]) / q[0]      # last increment's stiffness over the first's
    assert 0.0 < ratio < 0.95, (mat, dim, "not plastic enough to gate anything", q)
    for n in _THREADS[1:]:
        got, text = _run(mat, dim, n)
        assert got["build"] == base["build"]
        assert got["gravity_rc"] == base["gravity_rc"], (mat, dim, n)
        for s, (a, b) in enumerate(zip(base["steps"], got["steps"])):
            assert a[0] == b[0] and a[1] == b[1], (mat, dim, n, s, a[:2], b[:2])
            assert a[2] == b[2], (mat, dim, n, "nodal displacements differ at push step %d" % s)
        for w in _threaded_text(dim, n, nele):
            assert w in text, (mat, dim, n, w, text[-4000:])
        assert "PARTLY threaded" not in text, text[-4000:]
        if base["ledger"] is not None:
            # max and every counter exact; the average is an FP sum folded per thread, so under
            # a threaded commit loop it may differ in its last bits (258 doc section 2.2)
            assert got["ledger"]["max"] == base["ledger"]["max"], (mat, dim, n)
            assert got["ledger"]["refusals"] == base["ledger"]["refusals"], (mat, dim, n)
            for a, b in zip(base["ledger"]["avg"], got["ledger"]["avg"]):
                assert abs(a - b) <= 1.0e-12 * max(abs(a), 1.0e-300), (mat, dim, n, a, b)


# A purely cohesive DP (rho = 0) on the 2D deck reaches the apex, where the
# consistent tangent vanishes and DruckerPrager prints its throttled "NormCep"
# floor message from inside the threaded update loop -- i.e. from OpenMP
# WORKERS. Under openseespy `opserr` is a PythonStream, and before WP-258 its
# err_out() called PySys_FormatStderr from a thread with no Python thread
# state: a segfault (the WP-107 "ManzariDafalias IntScheme 1 crashes on a
# worker" was the same defect, behind ModifiedEuler's substep-cap warning).
# The deck's numbers are not meaningful (FixedNumIter on a softened apex); the
# gate is that the child survives, prints, and still matches 1 thread bit for bit.
APEX = ["DruckerPrager", 20000.0, 9000.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]


def test_worker_prints_under_openseespy():
    if len(_THREADS) < 2:
        pytest.skip("single hardware thread")
    MATERIALS["_apex"] = (APEX, None)
    try:
        base, text1 = _run("_apex", 2, 1)
        assert "NormCep" in text1, text1[-3000:]
        for n in _THREADS[1:]:
            got, text = _run("_apex", 2, n)
            assert "NormCep" in text, (n, text[-3000:])
            assert "element update loop THREADED on %d threads" % n in text, text[-3000:]
            assert got["steps"] == base["steps"], n
    finally:
        del MATERIALS["_apex"]
