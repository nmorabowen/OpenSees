"""WP-259 -- the Ladruno concretes on the WP-107 / WP-192 threaded element loops.

LadrunoConcrete3D and LadrunoRCConcrete answer `ladrunoThreadSafeUpdate()` /
`ladrunoThreadSafeForm()` true on the views the WP-259 audit proved
(Ladruno_implementation/259_thread_safe_concretes.md section 4). This file is the
determinism gate: the threaded result must be BIT-IDENTICAL to one thread, at 2, 4
and 8 threads, the elements must really run on the threads, and the deck must
really crack.

FE DECKS. Plain (unreinforced) concrete, N / mm / MPa, pushed by
`DisplacementControl` with `FixedNumIter` (a determinism gate needs identical
arithmetic, not convergence) until the root cracks (or, `3dc`, the prism crushes):

  * `3d`:    a 400 x 100 x 100 cantilever of `stdBrick` (update, tangent,
             residual and commit loops all thread), the tip face tied in z;
  * `3dc`:   the same prism crushed axially (tip face tied in x), so the
             compressive hardening / damage branch runs;
  * `3dcr`:  LadrunoConcrete3D pushed past the crush: return maps fail, steps are
             refused, and the process-wide warning budget saturates;
  * `pe`/`ps`: the same cantilever in 2D, `LadrunoQuad -type PlaneStrain` /
             `PlaneStress` (the update loop threads; the PlaneStress view runs
             the sigma_22 / sigma_33 = 0 condensation Newton);
  * `shell`: a 1000 x 1000 x 150 wall of `ASDShellQ4` + `LayeredShell` (the
             PlateFiber view, condensed; all four loops thread), pushed
             laterally at the top.

Elements are created in a SHUFFLED tag order. The digest is every node's
displacement after every push step and the load factor (`float.hex`), every
step's return code, and every Gauss point's own diagnostics at the end (damage;
LadrunoConcrete3D's `substeps` / `returnFailures`; LadrunoRCConcrete's
`implexError` under `-implex`). Neither material keeps a process-wide
floating-point accumulator, so everything is compared exactly; the one
process-wide counter (LadrunoConcrete3D's return-map warning budget, now
atomic) is compared by the number of lines it printed.

MPM mode F is not gated: it refuses both concretes at `mpm body` (ADR 232
section 6 admits ElasticIsotropic, LadrunoJ2, DruckerPrager, LadrunoSANISAND and
LadrunoHS only), so the opt-in cannot reach it (259 doc section 2.3).

Each configuration runs in a child process per thread count.
"""
import json
import os
import sys

import pytest

_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

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

_THREADS = [n for n in (1, 2, 4, 8) if _HW == 0 or n <= max(_HW, 1)]

# ---------------------------------------------------------------- materials (N, mm, MPa)
C3D = ["LadrunoConcrete3D", 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0]
E, NU, FT = 30000.0, 0.2, 3.0
RC = ["LadrunoRCConcrete", E, NU,
      "-Ce", 0.0, 0.0007, 0.0020, 0.0100, "-Cs", 0.0, 24.0, 30.0, 8.0,
      "-Cd", 0.0, 0.0, 0.25, 1.0 - 8.0 / (E * 0.0100),
      "-Te", 0.0, FT / E, 2.0e-3, 1.5e-2, "-Ts", 0.0, FT, 1.8, 0.15,
      "-Td", 0.0, 0.0, 1.0 - 1.8 / (E * 2.0e-3), 1.0 - 0.15 / (E * 1.5e-2),
      "-Kc", 2.0 / 3.0, "-beta"]

# name -> (nDMaterial args after the tag, hosts it runs on)
MATERIALS = {
    "c3d": (C3D + ["-lch", 50.0], ("3d", "3dc", "3dcr", "pe", "ps", "shell")),
    "c3d_implex": (C3D + ["-lch", 50.0, "-implex"], ("3d", "shell")),
    "c3d_eta": (C3D + ["-lch", 50.0, "-eta", 0.05], ("3d",)),
    "c3d_autoreg": (C3D + ["-autoRegularization"], ("3d", "pe")),
    # -verbose: every final failure prints (no budget) and every commit prints a summary
    "c3d_verbose": (C3D + ["-lch", 50.0, "-verbose"], ("3dcr",)),
    "rc": (RC, ("3d", "3dc", "ps", "shell")),
    "rc_implex": (RC + ["-implex"], ("3d", "shell")),
    "rc_autoreg": (RC + ["-autoRegularization", 50.0], ("3d", "shell")),
    "rc_interlock": (RC + ["-interlock", "-agg", 16.0, "-cyclic"], ("shell",)),
    # -implexControl with a tolerance every commit exceeds: the commit loop prints from workers
    "rc_implexctl": (RC + ["-implex", "-implexControl", 1.0e-14, 0.1], ("shell",)),
    # the WP-196 EnergyBalance element sweep (a fifth threaded loop, on ladrunoThreadSafeForm)
    "c3d_energy": (C3D + ["-lch", 50.0], ("3d+eb",)),
    "rc_energy": (RC, ("shell+eb",)),
}

_FE_CHILD = r'''
import json, random, sys
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[2])
from _testbed import ops
cfg = json.loads(sys.argv[3])
mat, host, nthreads = cfg["mat"], cfg["host"], cfg["threads"]
energy = host.endswith("+eb")
host = host.split("+")[0]
n_push, iters, ds = cfg["n_push"], cfg["iters"], cfg["ds"]
rc_mat = mat[0] == "LadrunoRCConcrete"
diag = ["damage"] + (["implexError"] if (rc_mat and "-implex" in mat) else []) \
       + ([] if rc_mat else ["substeps", "returnFailures"])
NL = 4                                                    # shell layers

ops.wipe()
if host in ("3d", "3dc", "3dcr"):
    L, Bw, Hd, NX, NY, NZ = 400.0, 100.0, 100.0, 12, 2, 4
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx, dy, dz = L / NX, Bw / NY, Hd / NZ
    nid = lambda i, j, k: 1 + i + (NX + 1) * (j + (NY + 1) * k)
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                ops.node(nid(i, j, k), i * dx, j * dy, k * dz)
    ops.nDMaterial(mat[0], 1, *mat[1:])
    cells = [(i, j, k) for k in range(NZ) for j in range(NY) for i in range(NX)]
    order = list(range(len(cells)))
    random.Random(259).shuffle(order)
    for tag0, c in enumerate(order):
        i, j, k = cells[c]
        ops.element("stdBrick", tag0 + 1,
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                    nid(i, j + 1, k + 1), 1)
    for k in range(NZ + 1):
        for j in range(NY + 1):
            ops.fix(nid(0, j, k), 1, 1, 1)
    tip = [nid(NX, j, k) for k in range(NZ + 1) for j in range(NY + 1)]
    master, vdof, ndf = (nid(NX, 0, NZ), 3, 3) if host == "3d" else (nid(NX, 0, 0), 1, 3)
    gps = [(e, ("material", gp)) for e in range(1, len(cells) + 1) for gp in range(1, 9)]
    nodes = [nid(i, j, k) for k in range(NZ + 1) for j in range(NY + 1) for i in range(NX + 1)]
elif host in ("pe", "ps"):
    L, Hd, NX, NZ = 400.0, 100.0, 12, 6
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    dx, dz = L / NX, Hd / NZ
    nid = lambda i, k: 1 + i + (NX + 1) * k
    for k in range(NZ + 1):
        for i in range(NX + 1):
            ops.node(nid(i, k), i * dx, k * dz)
    ops.nDMaterial(mat[0], 1, *mat[1:])
    cells = [(i, k) for k in range(NZ) for i in range(NX)]
    order = list(range(len(cells)))
    random.Random(259).shuffle(order)
    ptype = "PlaneStrain" if host == "pe" else "PlaneStress"
    for tag0, c in enumerate(order):
        i, k = cells[c]
        ops.element("LadrunoQuad", tag0 + 1, nid(i, k), nid(i + 1, k), nid(i + 1, k + 1),
                    nid(i, k + 1), 1, "-thick", 100.0, "-type", ptype)
    for k in range(NZ + 1):
        ops.fix(nid(0, k), 1, 1)
    tip = [nid(NX, k) for k in range(NZ + 1)]
    master, vdof, ndf = nid(NX, NZ), 2, 2
    gps = [(e, ("material", gp)) for e in range(1, len(cells) + 1) for gp in range(1, 5)]
    nodes = [nid(i, k) for k in range(NZ + 1) for i in range(NX + 1)]
else:                                                     # shell wall, PlateFiber view
    Lw, Hw, T, NX, NY = 1000.0, 1000.0, 150.0, 6, 6
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    dx, dy = Lw / NX, Hw / NY
    nid = lambda i, k: 1 + i + (NX + 1) * k
    for k in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, k), i * dx, k * dy, 0.0)
    ops.nDMaterial(mat[0], 1, *mat[1:])
    layers = []
    for _ in range(NL):
        layers += [1, T / NL]
    ops.section("LayeredShell", 1, NL, *layers)
    cells = [(i, k) for k in range(NY) for i in range(NX)]
    order = list(range(len(cells)))
    random.Random(259).shuffle(order)
    for tag0, c in enumerate(order):
        i, k = cells[c]
        ops.element("ASDShellQ4", tag0 + 1, nid(i, k), nid(i + 1, k), nid(i + 1, k + 1),
                    nid(i, k + 1), 1)
    for k in range(NY + 1):
        for i in range(NX + 1):
            if k == 0:
                ops.fix(nid(i, k), 1, 1, 1, 1, 1, 1)
            else:
                ops.fix(nid(i, k), 0, 0, 1, 1, 1, 1)   # planar membrane
    tip = [nid(i, NY) for i in range(NX + 1)]
    master, vdof, ndf = nid(0, NY), 1, 6
    gps = [(e, ("material", gp, "fiber", ly)) for e in range(1, len(cells) + 1)
           for gp in range(1, 5) for ly in range(1, NL + 1)]
    nodes = [nid(i, k) for k in range(NY + 1) for i in range(NX + 1)]

ops.ladrunoThreads(nthreads)
if energy:
    ops.region(1, "-eleRange", 1, len(cells))
    ops.recorder("EnergyBalance", "-v2", "-region", 1, "-time", "-precision", 17,
                 "-file", cfg["eb_file"])
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for n in tip:
    if n != master:
        ops.equalDOF(master, n, vdof)
ld = [0.0] * ndf
ld[vdof - 1] = -1.0 if host != "shell" else 1.0
ops.load(master, *ld)
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("FixedNumIter", iters, 0)
ops.algorithm("Newton")
ops.integrator("DisplacementControl", master, vdof, ds if host == "shell" else -ds)
ops.analysis("Static")
steps = []
for s in range(n_push):
    r = ops.analyze(1)
    q = ops.getLoadFactor(1)
    steps.append([r, q.hex(), [x.hex() for n in nodes for x in ops.nodeDisp(n)]])
pts = {}
for what in diag:
    pts[what] = [[float(v).hex() for v in ops.eleResponse(e, *path, what)] for e, path in gps]
dmg = [[float.fromhex(v) for v in p] for p in pts["damage"]]
ops.ladrunoThreads(1)
eb = None
if energy:
    ops.remove("recorders")
    eb = open(cfg["eb_file"]).read()
print("RESULT " + json.dumps({"build": ops.ladrunoBuild(), "steps": steps,
                              "q": [float.fromhex(s[1]) for s in steps], "pts": pts,
                              "dmax": [max(p[c] for p in dmg if len(p) > c) for c in range(2)],
                              "eb": eb}))
'''

_FE_CASES = [(m, h) for m in sorted(MATERIALS) for h in MATERIALS[m][1]]


def _base(host):
    return host.split("+")[0]
_NELE = {"3d": 96, "3dc": 96, "3dcr": 96, "pe": 72, "ps": 72, "shell": 36}
_DS = {"3d": 0.03, "3dc": 0.085, "3dcr": 0.12, "pe": 0.03, "ps": 0.03, "shell": 0.05}
_NPUSH = {"3d": 8, "3dc": 10, "3dcr": 9, "pe": 8, "ps": 8, "shell": 8}
_WARN = "return map did not converge"
_COUNTED = ("sigma_22 condensation", "sigma_33 condensation", "commit: ", "IMPL-EX error")
_SUPPRESSED = "further return-map warnings suppressed"


def _run_fe(mat, host, nthreads, iters=8):
    import tempfile
    fd, eb_file = tempfile.mkstemp(suffix=".eb.txt")
    os.close(fd)
    cfg = {"mat": MATERIALS[mat][0], "host": host, "threads": nthreads,
           "n_push": _NPUSH[_base(host)], "iters": iters, "ds": _DS[_base(host)],
           "eb_file": eb_file}
    code, out = run_python_script(_FE_CHILD, argv=[_DIST_BIN, _TESTS, json.dumps(cfg)],
                                  timeout=900)
    try:
        os.remove(eb_file)
    except OSError:
        pass
    assert code == 0, out[-3000:]
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    assert res, out[-3000:]
    return json.loads(res[-1][len("RESULT "):]), out


def _threaded_text(host, n, nele):
    """The loops that must have announced themselves on n threads."""
    want = ["element update loop THREADED on %d threads" % n, "(%d elements)" % nele]
    if _base(host) in ("3d", "3dc", "3dcr", "shell"):   # stdBrick / ASDShellQ4: the WP-192 loops
        for loop in ("tangent", "residual"):
            want.append("element %s loop threaded on %d threads: %d of %d FE_Elements computed "
                        "in parallel" % (loop, n, nele, nele))
    if host.endswith("+eb"):
        want.append("EnergyBalanceRecorder: element sweep threaded on %d threads: %d of %d"
                    % (n, nele, nele))
    return want


@pytest.mark.parametrize("mat,host", _FE_CASES)
def test_fe_threaded_equal_serial(mat, host):
    if len(_THREADS) < 2:
        pytest.skip("single hardware thread")
    base, text1 = _run_fe(mat, host, 1)
    rcs = [st[0] for st in base["steps"]]
    nwarn = text1.count(_WARN)
    if host == "3dcr":
        # pushed past the crush: return maps fail, steps are refused, and the process-wide
        # warning budget (the one shared counter, atomic since WP-259) saturates. The gate is
        # that the refusals, every point's failure count and the printed lines match 1 thread.
        assert any(r < 0 for r in rcs), (mat, host, rcs)
        if "-verbose" in MATERIALS[mat][0]:   # no budget: every failure prints, every commit sums
            assert nwarn > 20 and _SUPPRESSED not in text1 and "commit: " in text1, (mat, nwarn)
        else:
            assert nwarn == 20 and text1.count(_SUPPRESSED) == 1, (mat, host, nwarn)
        fails = [float.fromhex(p[0]) for p in base["pts"]["returnFailures"]]
        assert sum(fails) > 20, (mat, host, sum(fails))
    else:
        # the gate means something only if the push ran and the concrete cracked / crushed.
        # `3dc` sits at a sharp crush (LadrunoConcrete3D on nmb: no compressive damage after step
        # 9, omega_c 0.17 after step 10, refusals from step 11), so only its first 8 steps must
        # converge: libm differences (Zone-A is GCC / glibc) may move the crush by a step.
        conv = rcs[:8] if host == "3dc" else rcs
        assert all(r == 0 for r in conv), (mat, host, rcs)
        if "-implexControl" in MATERIALS[mat][0]:
            assert "IMPL-EX error" in text1, (mat, host, "the commit-loop print never ran")
        if host == "3dc":
            assert max(base["dmax"]) > 0.1, (mat, host, "no point cracked or crushed", base["dmax"])
        else:
            assert base["dmax"][0] > 0.3, (mat, host, "no point cracked", base["dmax"])
            q = base["q"]
            ratio = (q[-1] - q[-2]) / q[0]   # last increment's stiffness over the first's
            assert ratio < 0.95, (mat, host, "no softening over the push", q)
    for n in _THREADS[1:]:
        got, text = _run_fe(mat, host, n)
        assert got["build"] == base["build"]
        for s, (a, b) in enumerate(zip(base["steps"], got["steps"])):
            assert a[0] == b[0] and a[1] == b[1], (mat, host, n, s, a[:2], b[:2])
            assert a[2] == b[2], (mat, host, n, "nodal displacements differ at push step %d" % s)
        assert got["pts"] == base["pts"], (mat, host, n, "Gauss-point diagnostics differ")
        assert text.count(_WARN) == nwarn, (mat, host, n, "warning budget", nwarn)
        assert text.count(_SUPPRESSED) == text1.count(_SUPPRESSED), (mat, host, n)
        for w in _COUNTED:              # every other worker print: same number of lines
            assert text.count(w) == text1.count(w), (mat, host, n, w)
        assert got["eb"] == base["eb"], (mat, host, n, "EnergyBalance output differs")
        for w in _threaded_text(host, n, _NELE[_base(host)]):
            assert w in text, (mat, host, n, w, text[-4000:])
        assert "PARTLY threaded" not in text, text[-4000:]
