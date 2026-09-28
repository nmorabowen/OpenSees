"""WP-147 -- the class-wide tangent/residual pools leave slot [MAX_NUM_DOF] uninitialized.

`FE_Element`, `TransformationFE`, `DOF_Group` and `TransformationDOF_Group` each keep a
class-wide array of `MAX_NUM_DOF+1` Matrix/Vector pointers indexed by the DOF count, and
take the POOLED branch for `numDOF <= MAX_NUM_DOF`. Upstream zeroed only `i < MAX_NUM_DOF`,
so the last slot -- the one an object with EXACTLY `MAX_NUM_DOF` DOFs uses -- held whatever
the heap left there:

  * nonzero garbage  -> taken as an existing Vector/Matrix -> the element's tangent and
    residual live at a wild address (access violation, or silent corruption);
  * zero             -> works, and the destructor (same off-by-one) leaks that slot.

Which one you get is heap luck, which is how it survived. For FE_Element and
TransformationFE `MAX_NUM_DOF` is 64, reachable by an ordinary coupling element; for
TransformationDOF_Group (16) and DOF_Group (256) the fix is by inspection only.

How the test makes the luck deterministic: right before the analysis objects are built
(the first FE_Element of a model allocates the pool), the child FILLS the allocator's
free list for the pool array's size (65 pointers = 520 bytes) with a nonzero byte pattern
and frees it again, so the pool array is carved from dirty memory. Unfixed, slot [64] then
reads 0xA5A5... and the element dereferences it; fixed, it is zeroed like every other slot.
Same heap on both platforms: MSVC's UCRT (static or DLL) allocates from GetProcessHeap(),
glibc's operator new from malloc. The poison is best-effort by construction (it cannot
prove the pool array landed on a poisoned block), so the mutation record in the PR body
is what shows it bites; the green assertions below are exact regardless.

Two decks, one 64-DOF element each -- an RBE3 (LadrunoDistributingCoupling) in 2D:
reference node ndf 3 + one independent with ndf 3 + 29 independents with ndf 2 = 64.

  * `plain`          -- `constraints Plain`, independents fixed: the RBE3's FE_Element
                        takes FE_Element pool slot [64].
  * `transformation` -- `constraints Transformation`, independents on ground springs, one
                        of them tied by `equalDOF` so the RBE3 touches an MP-constrained
                        node (-> a TransformationFE) whose transformed DOF count is STILL
                        64 (equalDOF keeps a 2-DOF node at 2) -> TransformationFE slot [64]
                        as well as the FE_Element slot of its base.

Each deck runs several build -> poison -> analyze -> check -> wipe cycles in a child
interpreter (a wild pointer kills the process, which no in-process assert can observe).
The physics check is the RBE3 force distribution, exact by equilibrium for a load at the
centroid of a symmetric ring: a force splits equally, a moment becomes equal tangential
forces M / (N a).
"""
import os

import pytest

from _testbed import ops
from _testbed.subprocess_run import run_python_script

pytestmark = [pytest.mark.zone_a]

# The child must load the SAME binary the parent resolved (see
# test_adr85_contact2d_t0_refusals.py for why only this directory is pinned).
ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))

N_INDEP = 30          # 1 x ndf-3 + 29 x ndf-2 independents; + ndf-3 reference = 64 DOF
CYCLES = 4

CHILD = r'''
import ctypes, math, os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    _add = getattr(os, "add_dll_directory", None)   # Windows-only; AttributeError elsewhere
    if _add is not None:
        try:
            _add(_D)
        except OSError:
            pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# a site-packages boot .pth can preload `opensees` from ANOTHER build before this
# script runs; then the pin above is silently ignored. Refuse rather than test it.
_got = os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__)))
if os.path.isdir(_D) and _got != os.path.normcase(_D):
    print("WP147_WRONG_BINARY", _got, flush=True)
    raise SystemExit(7)

DECK = sys.argv[1]
N = %(N_INDEP)d
CYCLES = %(CYCLES)d
A = 1.0                       # ring radius
FX, FY, MZ = 120.0, -45.0, 30.0
KS = 1.0e5                    # ground springs (transformation deck)


def _heap_api():
    if sys.platform == "win32":
        k32 = ctypes.windll.kernel32
        k32.GetProcessHeap.restype = ctypes.c_void_p
        k32.HeapAlloc.restype = ctypes.c_void_p
        k32.HeapAlloc.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_size_t]
        k32.HeapFree.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_void_p]
        h = k32.GetProcessHeap()
        return (lambda n: k32.HeapAlloc(h, 0, n)), (lambda p: k32.HeapFree(h, 0, p))
    libc = ctypes.CDLL(None)
    libc.malloc.restype = ctypes.c_void_p
    libc.malloc.argtypes = [ctypes.c_size_t]
    libc.free.argtypes = [ctypes.c_void_p]
    return libc.malloc, libc.free


_ALLOC, _FREE = _heap_api()


def poison(size=65 * 8, count=4096, byte=0xA5):
    """Dirty the free list for `size`-byte blocks: allocate, fill, free."""
    ptrs = []
    for _ in range(count):
        p = _ALLOC(size)
        if p:
            ctypes.memset(p, byte, size)
            ptrs.append(p)
    for p in ptrs:
        _FREE(p)


def ring():
    return [(A * math.cos(2 * math.pi * k / N), A * math.sin(2 * math.pi * k / N))
            for k in range(N)]


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0, "-ndf", 3)                       # reference, at the centroid
    xy = ring()
    indep = []
    for k, (x, y) in enumerate(xy):
        tag = 2 + k
        if k == 0:
            ops.node(tag, x, y, "-ndf", 3)                 # the one ndf-3 independent
        else:
            ops.node(tag, x, y)
        indep.append(tag)
    ops.element("LadrunoDistributingCoupling", 1, 1, N, *indep, "-k", 1.0e8)

    ground = {}                                            # indep -> node whose reaction = -f_i
    if DECK == "plain":
        ops.fix(2, 1, 1, 1)
        for tag in indep[1:]:
            ops.fix(tag, 1, 1)
        ground = {t: t for t in indep}
    elif DECK == "transformation":
        ops.uniaxialMaterial("Elastic", 1, KS)
        for k, tag in enumerate(indep):
            x, y = xy[k]
            spring_node = tag
            if k == 1:
                # indep 3 is tied by equalDOF to a retained twin that carries its spring:
                # node 3 becomes MP-constrained (-> the RBE3 gets a TransformationFE) while
                # its transformed DOF count stays 2, so the RBE3's stays 64.
                spring_node = 2000 + tag
                ops.node(spring_node, x, y)
                ops.equalDOF(spring_node, tag, 1, 2)
            g = 1000 + tag
            if k == 0:
                ops.node(g, x, y, "-ndf", 3)
                ops.fix(g, 1, 1, 1)
                ops.element("zeroLength", 100 + tag, g, spring_node, "-mat", 1, 1, 1,
                            "-dir", 1, 2, 3)
            else:
                ops.node(g, x, y)
                ops.fix(g, 1, 1)
                ops.element("zeroLength", 100 + tag, g, spring_node, "-mat", 1, 1,
                            "-dir", 1, 2)
            ground[tag] = g
    else:
        raise SystemExit("unknown deck " + DECK)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, FX, FY, MZ)
    return xy, indep, ground


def check(xy, indep, ground):
    ndof = sum(len(ops.nodeDisp(n)) for n in ops.eleNodes(1))
    print("WP147_NDOF", ndof, flush=True)
    ops.reactions()
    worst = 0.0
    sx = sy = mz = 0.0
    for k, tag in enumerate(indep):
        x, y = xy[k]
        # RBE3 at the centroid of a symmetric ring, equal weights: f_i = F/N + M (-y, x)/(N a^2)
        fx = FX / N - MZ * y / (N * A * A)
        fy = FY / N + MZ * x / (N * A * A)
        r = ops.nodeReaction(ground[tag])
        sx += r[0]
        sy += r[1]
        mz += x * r[1] - y * r[0]
        worst = max(worst, abs(r[0] + fx), abs(r[1] + fy))
    scale = abs(FX) / N
    print("WP147_DIST %%.3e" %% (worst / scale), flush=True)
    print("WP147_SUMS %%.9e %%.9e %%.9e" %% (sx, sy, mz), flush=True)


for cycle in range(CYCLES):
    xy, indep, ground = build()
    poison()                          # the pool array is allocated below, from this memory
    ops.constraints("Plain" if DECK == "plain" else "Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    rc = ops.analyze(1)
    print("WP147_RC", cycle, rc, flush=True)
    if rc != 0:
        raise SystemExit(3)
    check(xy, indep, ground)
    ops.wipe()                        # last FE_Element dies -> pool freed (and slot [64] with it)
    print("WP147_CYCLE_OK", cycle, flush=True)

print("WP147_DONE", flush=True)
'''


def _run(deck):
    script = CHILD % {"ENGINE_DIR": ENGINE_DIR, "N_INDEP": N_INDEP, "CYCLES": CYCLES}
    return run_python_script(script, argv=(deck,), timeout=240)


@pytest.mark.parametrize("deck", ["plain", "transformation"])
def test_64dof_element_survives_a_dirty_pool_array(deck):
    rc, out = _run(deck)
    tail = out[-3000:]
    assert "WP147_WRONG_BINARY" not in out, tail
    # the process must SURVIVE every cycle -- a wild slot-[64] pointer kills it
    assert rc == 0, "child died (rc=%s) -- slot [MAX_NUM_DOF] garbage?\n%s" % (rc, tail)
    assert "WP147_DONE" in out, tail
    assert out.count("WP147_CYCLE_OK") == CYCLES, tail

    lines = out.splitlines()
    # the deck must actually exercise slot [64]: guard against deck drift making this
    # test observe some other pool slot (a green test on the wrong path proves nothing)
    ndofs = [int(l.split()[1]) for l in lines if l.startswith("WP147_NDOF")]
    assert ndofs == [64] * CYCLES, ndofs

    # the RBE3 distribution is exact by equilibrium; penalty only moves the reference
    for l in lines:
        if l.startswith("WP147_DIST"):
            assert float(l.split()[1]) < 1.0e-6, l
        if l.startswith("WP147_SUMS"):
            sx, sy, mz = (float(v) for v in l.split()[1:])
            assert sx == pytest.approx(-120.0, rel=1e-8)
            assert sy == pytest.approx(45.0, rel=1e-8)
            assert mz == pytest.approx(-30.0, rel=1e-8)
