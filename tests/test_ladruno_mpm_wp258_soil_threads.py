"""WP-258 (program #1075): soil particles thread their material calls in MPM mode F.

WP-253 threads mode F's transfer loops and, when every particle material answers
`ladrunoThreadSafeUpdate()` AND `ladrunoThreadSafeForm()`, the material calls too
(`setTrialStrain`, `getStress`, `commitState`, `getInitialTangent` on workers). WP-258 makes
DruckerPrager, LadrunoSANISAND (IntScheme 1 / 2 / 129, `-implex`) and LadrunoHS (HSSmall,
`-implex`) answer both (Ladruno_implementation/258_thread_safe_soil_materials.md section 4),
lifting ADR 232 wall W14. This file is the determinism gate for the MPM half: every particle
array after K steps is EXACTLY equal (float64 ==) at 2, 4 and 8 threads to the 1-thread run,
and `mpm debug threads` reports the material loop threaded.

The soil decks are staged, as in 232k: a laterally confined column settles under gravity at
stage 0 (elastic), `mpm stage <tag> 1` flips every particle (SANISAND's class-wide flag and
`Elastic2Plastic` reset, HS's per-tag registry -- the HS particles then pick the stage up
lazily, on the workers), and a rigid `gridvel` plate then pushes the top down, so the plastic
branch runs at every particle. SANISAND runs on `bspline2` only (232k's K0 admission). DP runs
WP-253's collapsing column (`gimp`, a frictional floor). One SANISAND case adds its particles
in a shuffled order. Steps are taken at 0.25 dt_crit of each step's own start, so the dt
sequence is itself bit-identical; nothing is frozen.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS))
from _testbed.subprocess_run import run_python_script  # noqa: E402

_PROBE = (
    "import sys; sys.path.insert(0, %r)\n"
    "from _testbed import ops\n"
    "print('LADRUNO_PROBE_HW', ops.ladrunoThreads(1 << 20))\n"
) % str(TESTS)


def _openmp():
    try:
        code, out = run_python_script(_PROBE)
    except Exception:
        return True, ""
    if code != 0 or "LADRUNO_PROBE_HW" not in out:
        return True, ""
    if "built WITHOUT LADRUNO_OPENMP" in out:
        return False, "binary built with LADRUNO_OPENMP=OFF: nothing to thread"
    return True, ""


_OPENMP_IN, _SKIP_WHY = _openmp()
pytestmark = [pytest.mark.zone_a, pytest.mark.skipif(not _OPENMP_IN, reason=_SKIP_WHY)]

FIELDS = ("x", "v", "V", "F", "R", "epsHat", "sigHat", "sig")

SANISAND = [264.32, 0.3129, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3, 0.968,
            3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
HS = [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0, "-nu", 0.2, "-pref", 100.0, "-Rf", 0.9,
      "-OCR", 1.5, "-stage", 0]
DP = ["DruckerPrager", 1, 1.6666666666666666e7, 7.6923076923076922e6, 50.0, 0.21, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

SAN_TAIL = ["-Presidual", 0.0, "-Pmin", 1.0e-4 * 101.0]
SOILS = {
    "sanisand_is1": ["LadrunoSANISAND", 3, *SANISAND, 1, 2, 1, 1e-7, 1e-7, *SAN_TAIL],
    "sanisand_is2": ["LadrunoSANISAND", 3, *SANISAND, 2, 2, 1, 1e-7, 1e-7, *SAN_TAIL],
    "sanisand_is129": ["LadrunoSANISAND", 3, *SANISAND, 129, 2, 1, 1e-7, 1e-7, *SAN_TAIL],
    "sanisand_implex": ["LadrunoSANISAND", 3, *SANISAND, 2, 2, 1, 1e-7, 1e-7, *SAN_TAIL,
                        "-implex", "-maxSubsteps", 1000],
    "hs": ["LadrunoHS", 3, *HS],
    "hs_small": ["LadrunoHS", 3, *HS, "-smallStrain", 1.0e5, 1.0e-4],
    "hs_implex": ["LadrunoHS", 3, *HS, "-implex"],
}

# The child builds one deck, runs it at one thread count, and prints the particle arrays.
_CHILD = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
from _testbed import ops
cfg = json.loads(sys.argv[2])


def l(x):
    return [float(v) for v in x] if isinstance(x, (list, tuple)) else [float(x)]


def steps(k):
    for _ in range(k):
        ops.mpm("analyze", 1, 0.25 * float(ops.mpm("criticalTimeStep")))


def particles(dim, hs, seed):
    # the shuffled column, generated here (288 commands would overflow a Windows command line)
    import random
    if dim == 3:
        pts = [((i + 0.5) * hs, (j + 0.5) * hs, (k + 0.5) * hs)
               for k in range(12) for j in range(4) for i in range(6)]
    else:
        pts = [((i + 0.5) * hs, (k + 0.5) * hs) for k in range(16) for i in range(10)]
    random.Random(seed).shuffle(pts)
    for p, x in enumerate(pts):
        ops.mpm("particle", p, 1, 3, 2.0, hs ** dim, "-x", *x, "-lp", *([hs / 2] * dim))


for name, args in cfg["cmds"]:
    if name == "_particles":
        particles(*args)
    else:
        getattr(ops, name)(*args)
ops.ladrunoThreads(cfg["threads"])
n = int(ops.mpm("info", "numParticles"))
x0 = l(ops.mpm("particleData", "x"))
if cfg["stage"] is not None:
    steps(cfg["k0"])
    ops.mpm("stage", 3, 1)
    for name, args in cfg["after_stage"]:
        getattr(ops, name)(*args)
steps(cfg["k1"])
data = {f: l(ops.mpm("particleData", f)) for f in cfg["fields"]}
data["moved"] = max(abs(a - b) for a, b in zip(data["x"], x0))
thr = [int(v) for v in l(ops.mpm("debug", "threads"))]
ops.ladrunoThreads(1)
print("RESULT " + json.dumps({"build": ops.ladrunoBuild(), "n": n, "thr": thr, "data": data}))
'''


def _soil_column(dim, mat, shuffled=False):
    """A laterally confined column on a fixed base (bspline2, images on the fix faces); kPa / t."""
    h = 0.25
    cmds = [("wipe", []), ("model", ["basic", "-ndm", dim, "-ndf", dim]),
            ("nDMaterial", list(SOILS[mat]))]
    if dim == 3:
        cmds.append(("mpm", ["grid", 232, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", 3, 2, 8,
                             "-basis", "bspline2"]))
        box = [0.0, 0.0, 0.0, 0.75, 0.5, 1.5]
        faces = [("xmin", "x"), ("xmax", "x"), ("ymin", "y"), ("ymax", "y"), ("zmin", "xyz")]
        grav = [0.0, 0.0, -9.81]
        # inside the fix faces (a moving box must not reach their image ghosts, WP-232i B10)
        plate = ["box", 0.001, 0.001, 1.4, 0.749, 0.499, 3.0, "v", "free", "free", -0.05]
    else:
        cmds.append(("mpm", ["grid", 232, "-origin", 0.0, 0.0, "-h", h, "-cells", 5, 9,
                             "-basis", "bspline2"]))
        box = [0.0, 0.0, 1.25, 2.0]
        faces = [("xmin", "x"), ("xmax", "x"), ("ymin", "xy")]
        grav = [0.0, -9.81]
        plate = ["box", 0.001, 1.9, 1.249, 3.0, "v", "free", -0.05]
    if shuffled:
        cmds.append(("_particles", [dim, h / 2.0, 258]))   # built in the child, shuffled
    else:
        cmds.append(("mpm", ["body", 1, "box", *box, "-ppc", 2, "-mat", 3, "-rho", 2.0]))
    for face, dofs in faces:
        cmds.append(("mpm", ["fix", face, dofs]))
    cmds.append(("mpm", ["gravity", *grav]))
    after = [("mpm", ["gridvel", "plate", *plate])]
    return cmds, after


def _dp_column(dim):
    """WP-253's DruckerPrager column (gimp, frictional floor), in 2D and 3D."""
    if dim == 3:
        return [("wipe", []), ("model", ["basic", "-ndm", 3, "-ndf", 3]), ("nDMaterial", DP),
                ("mpm", ["grid", 232, "-origin", 0.0, 0.0, 0.0, "-h", 0.25, "-cells", 6, 4, 10,
                         "-basis", "gimp"]),
                ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, 0.75, 0.5, 1.5, "-ppc", 2, "-mat", 1,
                         "-rho", 2000.0]),
                ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "ymin", "y"]),
                ("mpm", ["wall", "zmin", 0.3]), ("mpm", ["gravity", 0.0, 0.0, -10.0])]
    return [("wipe", []), ("model", ["basic", "-ndm", 2, "-ndf", 2]), ("nDMaterial", DP),
            ("mpm", ["grid", 232, "-origin", 0.0, 0.0, "-h", 0.25, "-cells", 10, 12,
                     "-basis", "gimp"]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2000.0]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["wall", "ymin", 0.3]),
            ("mpm", ["gravity", 0.0, -10.0])]


def _case(name):
    if name.startswith("dp_"):
        dim = int(name[3])
        return dict(cmds=_dp_column(dim), after_stage=[], stage=None, k0=0, k1=24)
    mat, tail = name.rsplit("_", 1)          # "<mat>_<dim>d" or "<mat>_<dim>ds" (shuffled)
    dim, shuffled = int(tail[0]), tail.endswith("s")
    cmds, after = _soil_column(dim, mat, shuffled)
    return dict(cmds=cmds, after_stage=after, stage=1, k0=8, k1=16)


CASES = (["dp_2d", "dp_3d"]
         + ["%s_%dd" % (m, d) for m in sorted(SOILS) for d in (2, 3)]
         + ["sanisand_is2_3ds"])


def _hw():
    code, out = run_python_script(_PROBE)
    for line in out.splitlines():
        if line.startswith("LADRUNO_PROBE_HW"):
            return int(line.split()[1])
    return 1


def _run(case, nthreads):
    import json
    cfg = dict(case, threads=nthreads, fields=list(FIELDS))
    code, out = run_python_script(_CHILD, argv=[str(TESTS), json.dumps(cfg)], timeout=900)
    assert code == 0, out[-3000:]
    res = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    assert res, out[-3000:]
    return json.loads(res[-1][len("RESULT "):])


@pytest.mark.parametrize("name", CASES)
def test_soil_particles_threaded_equal_serial(name):
    hw = _hw()
    ns = sorted({n for n in (2, 4, min(8, hw)) if 2 <= n <= hw})
    if not ns:
        pytest.skip("this machine has one hardware thread")
    case = _case(name)
    ref = _run(case, 1)
    assert ref["thr"] == [1, 1], ref["thr"]
    assert ref["n"] > 50, ref["n"]
    assert ref["data"]["moved"] > 0.0, "the model must move for the comparison to mean anything"
    for nt in ns:
        got = _run(case, nt)
        assert got["build"] == ref["build"]
        # the material loop threads too: [transfer-loop threads, material-loop threads]
        assert got["thr"] == [nt, nt], (name, nt, got["thr"])
        for f in FIELDS:
            a, b = ref["data"][f], got["data"][f]
            assert len(a) == len(b), (name, nt, f)
            bad = [i for i, (u, v) in enumerate(zip(a, b)) if not (u == v or (math.isnan(u) and math.isnan(v)))]
            assert not bad, f"{name}: {f} differs at {nt} threads in {len(bad)} entries (first {bad[0]})"
