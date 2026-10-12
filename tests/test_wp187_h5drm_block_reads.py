# ============================================================================
# test_wp187_h5drm_block_reads.py  (WP-187)
# ----------------------------------------------------------------------------
# The H5DRM pattern now reads its motion datasets through a chunk cache sized
# per rank and serves the per-step samples from an in-memory block of time
# steps (H5DRMLoadPattern::ensure_block / block_value). The values are the same
# file bytes and the interpolation is unchanged, so every response must be
# BIT-IDENTICAL to the old per-node path, which the binary keeps behind
# `-prefetch 0 -chunkCacheMB 0`.
#
# Fixture: a 6 x 5 x 4 station grid (all stations matched to model nodes, the
# outer shell external, the rest internal), gzip-chunked motion datasets with
# small chunks (9 rows x 7 steps) so a step crosses chunk rows and time blocks,
# and a SHUFFLED data_location (station p's rows are not at 3p), so a rank's
# rows form many short runs. Transient run through tend and past it (hold
# branch). Variants compared with np.array_equal (no tolerance):
#   legacy  : -prefetch 0 -chunkCacheMB 0   (the pre-WP-187 read path)
#   default : no flags                       (chunk-length block + sized cache)
#   p1, p5  : -prefetch 1 / 5                (block boundaries every step / mid-chunk)
#   big     : -prefetch 1000 -chunkCacheMB 3 (one block for the whole record)
#   verbose : default + -verbose             (log only)
# The flags follow the 12 T / x0 values: equality with the no-flag run also
# proves the positional parse still takes T and x0.
#
# Runs through openseespy (worktree dist/bin/opensees.pyd) and, when present,
# the classic Tcl OpenSees.exe (TclPatternCommand's own flag filter).
# ============================================================================
import os
import subprocess
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_DIST = os.path.normpath(os.path.join(_HERE, "..", "dist", "bin"))

h5py = pytest.importorskip("h5py")
if not os.path.exists(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip("worktree dist/bin/opensees.pyd not built", allow_module_level=True)
os.add_dll_directory(_DIST)
sys.path.insert(0, _DIST)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees as ops  # noqa: E402

NX, NY, NZ = 6, 5, 4
H = 0.01            # km
CRD = 1000.0        # km -> m
DT = 0.01
NT = 40             # samples; tend = (NT - 1) * DT
TEND = (NT - 1) * DT
NSTEPS = 52         # runs past tend
CHUNK = (9, 7)
T = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
X0 = (0.0, 0.0, 0.0)

VARIANTS = {
    "legacy": ["-prefetch", 0, "-chunkCacheMB", 0],
    "default": [],
    "p1": ["-prefetch", 1],
    "p5": ["-prefetch", 5],
    "big": ["-prefetch", 1000, "-chunkCacheMB", 3],
    "verbose": ["-verbose"],
}


def _grid():
    pts = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                pts.append((i, j, k))
    return pts


def _internal(i, j, k):
    return 0 < i < NX - 1 and 0 < j < NY - 1 and 0 < k < NZ - 1


def write_h5drm(path):
    st = _grid()
    n = len(st)
    rng = np.random.default_rng(187)
    xyz = np.array([(i * H, j * H, k * H) for i, j, k in st], dtype=float)
    perm = rng.permutation(n)                     # station p -> slot perm[p]
    t = np.arange(NT) * DT
    dis = np.zeros((3 * n, NT))
    acc = np.zeros((3 * n, NT))
    for p in range(n):
        for c in range(3):
            w = 2 * np.pi * (1.5 + rng.random() * 3.0)
            ph = rng.random() * 2 * np.pi
            amp = 1e-3 * (0.5 + rng.random())
            row = 3 * perm[p] + c
            dis[row] = amp * np.sin(w * t + ph) + 1e-4 * rng.standard_normal(NT)
            acc[row] = -amp * w * w * np.sin(w * t + ph) + 1e-2 * rng.standard_normal(NT)
    with h5py.File(path, "w") as f:
        d = f.create_group("DRM_Data")
        d["xyz"] = xyz
        d["internal"] = np.array([_internal(*s) for s in st], dtype=bool)
        d["data_location"] = (3 * perm).astype(np.int32)
        d.create_dataset("displacement", data=dis, chunks=CHUNK, compression="gzip", compression_opts=4)
        d.create_dataset("acceleration", data=acc, chunks=CHUNK, compression="gzip", compression_opts=4)
        d.create_dataset("velocity", data=np.gradient(dis, DT, axis=1), chunks=CHUNK,
                         compression="gzip", compression_opts=4)
        m = f.create_group("DRM_Metadata")
        m["dt"] = DT
        m["tstart"] = 0.0
        m["tend"] = TEND
        m["drmbox_x0"] = np.array([0.0, 0.0, 0.0])
        for kname, v in (("drmbox_xmin", 0.0), ("drmbox_xmax", (NX - 1) * H),
                         ("drmbox_ymin", 0.0), ("drmbox_ymax", (NY - 1) * H),
                         ("drmbox_zmin", 0.0), ("drmbox_zmax", (NZ - 1) * H)):
            m[kname] = float(v)
        m["h"] = np.array([H, H, H])


def _tag(i, j, k):
    return 1 + (k * NY + j) * NX + i        # node tags NOT in station order


def _bricks():
    out = []
    e = 1
    for i in range(NX - 1):
        for j in range(NY - 1):
            for k in range(NZ - 1):
                ns = [_tag(i, j, k), _tag(i + 1, j, k), _tag(i + 1, j + 1, k), _tag(i, j + 1, k),
                      _tag(i, j, k + 1), _tag(i + 1, j, k + 1), _tag(i + 1, j + 1, k + 1), _tag(i, j + 1, k + 1)]
                out.append((e, ns))
                e += 1
    return out


def run_py(path, flags):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                ops.node(_tag(i, j, k), i * H * CRD, j * H * CRD, k * H * CRD)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e8, 0.3, 2.0)
    for e, ns in _bricks():
        ops.element("stdBrick", e, *ns, 1)
    for i in range(NX):
        for j in range(NY):
            ops.fix(_tag(i, j, 0), 1, 1, 1)
    ops.pattern("H5DRM", 1, path, 1.0, CRD, 1.0e-3, 1, *T, *X0, *flags)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    tags = sorted(ops.getNodeTags())
    hist = []
    for _ in range(NSTEPS):
        assert ops.analyze(1, DT) == 0
        hist.append([ops.nodeDisp(nt) for nt in tags])
    ops.wipe()
    return np.array(hist)


@pytest.fixture(scope="module")
def drmfile(tmp_path_factory):
    p = str(tmp_path_factory.mktemp("wp187") / "grid.h5drm")
    write_h5drm(p)
    return p


@pytest.fixture(scope="module")
def py_runs(drmfile):
    return {name: run_py(drmfile, fl) for name, fl in VARIANTS.items()}


def test_drm_moves_the_model(py_runs):
    h = py_runs["legacy"]
    assert np.all(np.isfinite(h))
    assert np.max(np.abs(h)) > 1e-6, "DRM applied no motion: the fixture does not test anything"


@pytest.mark.parametrize("name", [k for k in VARIANTS if k != "legacy"])
def test_block_reads_bit_identical_py(py_runs, name):
    a, b = py_runs["legacy"], py_runs[name]
    assert a.shape == b.shape
    assert np.array_equal(a, b), (
        f"{name}: max |diff| {np.max(np.abs(a - b)):.3e} -- block path is not bit-identical")


def test_unknown_flag_is_refused(drmfile):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    with pytest.raises(Exception):
        ops.pattern("H5DRM", 1, drmfile, 1.0, CRD, 1.0e-3, 1, *T, *X0, "-noSuchFlag")
    ops.wipe()


# ---------------------------------------------------------------------------
# classic Tcl OpenSees.exe (TclPatternCommand flag filter)
# ---------------------------------------------------------------------------
_EXE = os.path.join(_DIST, "OpenSees.exe")


def _tcl_deck(path, flags, out):
    L = ["model basic -ndm 3 -ndf 3"]
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                L.append(f"node {_tag(i, j, k)} {i * H * CRD!r} {j * H * CRD!r} {k * H * CRD!r}")
    L.append("nDMaterial ElasticIsotropic 1 2.0e8 0.3 2.0")
    for e, ns in _bricks():
        L.append(f"element stdBrick {e} {' '.join(map(str, ns))} 1")
    for i in range(NX):
        for j in range(NY):
            L.append(f"fix {_tag(i, j, 0)} 1 1 1")
    f = " ".join(str(x) for x in flags)
    tx = " ".join(repr(x) for x in (*T, *X0))
    L.append(f"pattern H5DRM 1 {{{path}}} 1.0 {CRD!r} 1.0e-3 1 {tx} {f}")
    L += ["constraints Transformation", "numberer RCM", "system UmfPack",
          "test NormDispIncr 1e-10 10", "algorithm Linear", "integrator Newmark 0.5 0.25",
          "analysis Transient",
          f"set fo [open {{{out}}} w]",
          "set tags [lsort -integer [getNodeTags]]",
          f"for {{set s 0}} {{$s < {NSTEPS}}} {{incr s}} {{",
          f"  if {{[analyze 1 {DT!r}] != 0}} {{ puts $fo FAIL; break }}",
          "  foreach t $tags { foreach v [nodeDisp $t] { puts -nonewline $fo [format \"%.17g \" $v] } }",
          "  puts $fo \"\"",
          "}",
          "close $fo"]
    return "\n".join(L) + "\n"


@pytest.mark.skipif(not os.path.exists(_EXE), reason="worktree dist/bin/OpenSees.exe not built")
def test_block_reads_bit_identical_classic_tcl(drmfile, tmp_path):
    outs = {}
    for name in ("legacy", "default", "p5", "big"):
        out = str(tmp_path / f"{name}.txt").replace("\\", "/")
        deck = tmp_path / f"{name}.tcl"
        deck.write_text(_tcl_deck(drmfile.replace("\\", "/"), VARIANTS[name], out))
        r = subprocess.run([_EXE, str(deck)], cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr[-2000:]
        txt = open(out).read()
        assert "FAIL" not in txt and txt.count("\n") == NSTEPS, txt[:200]
        outs[name] = txt
    vals = np.array(outs["legacy"].split(), dtype=float)
    assert np.max(np.abs(vals)) > 1e-6
    for name in ("default", "p5", "big"):
        assert outs[name] == outs["legacy"], f"classic Tcl {name}: output differs from the legacy read path"
    # the per-step log lines are gone without -verbose, back with it
    deck = tmp_path / "verbose.tcl"
    out = str(tmp_path / "verbose.txt").replace("\\", "/")
    deck.write_text(_tcl_deck(drmfile.replace("\\", "/"), ["-verbose"], out))
    rv = subprocess.run([_EXE, str(deck)], cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    assert rv.returncode == 0
    assert open(out).read() == outs["legacy"]
    deck = tmp_path / "quiet.tcl"
    out = str(tmp_path / "quiet.txt").replace("\\", "/")
    deck.write_text(_tcl_deck(drmfile.replace("\\", "/"), [], out))
    rq = subprocess.run([_EXE, str(deck)], cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    logv, logq = rv.stdout + rv.stderr, rq.stdout + rq.stderr
    nv = logv.count("H5DRM(0) : t = ")
    nq = logq.count("H5DRM(0) : t = ")
    assert nv >= NSTEPS and nq == 0, (nv, nq)
    # the block path really served the default run (else the equality above is vacuous) ...
    assert "WP-187 block 1:" in logv and "WP-187 reads:" in logq, logv[-1500:]
    # ... and -prefetch 0 really is the per-node path
    deck = tmp_path / "legacy_verbose.tcl"
    out = str(tmp_path / "legacy_verbose.txt").replace("\\", "/")
    deck.write_text(_tcl_deck(drmfile.replace("\\", "/"), VARIANTS["legacy"] + ["-verbose"], out))
    rl = subprocess.run([_EXE, str(deck)], cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    assert rl.returncode == 0
    assert open(out).read() == outs["legacy"]
    assert "WP-187 block" not in (rl.stdout + rl.stderr)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
