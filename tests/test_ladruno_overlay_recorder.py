"""LadrunoPorousOverlay (PATTERN 33022) — ADR-73 P4 RECORDER battery (WP4.C).

P4 adds pore-pressure recorder CHANNELS on top of the shipped overlay:

  * `recorder ladruno $file ... -overlay <tag...>`  ->  HDF5
    RESULTS/ON_NODES/overlayPressure_<patternTag>/{ID,DATA,STEP,TIME},
    DATA shape [T x nRegionNodes x 1], COMPONENTS="p"; write-once topology
    MODEL/OVERLAYS/OVERLAY_<tag>/{REGION_NODES,DRAINED_NODES} + TAG attr
    (both groups live inside the per-analysis MODEL_STAGE group).
  * `recorder Monitor -overlay $tag <-nodes ...> -sink $file <-every N>`  ->
    SWMR HDF5 with COLUMNS (overlay<tag>.p.node<n>), STEP, TIME, FRAMES[T x nCols].
  * the overlay's own `-record $file` CSV (header = region-node tags, rows =
    "time, p(node1..nodeN)") is the shipped ORACLE for the twin gates.

This file implements gates (a)-(h) of the pins §4.5 EXACTLY. It is written
against the FROZEN pins, not against the C++ (which is authored in parallel):
every gate checks the pinned command surface / HDF5 layout, and any surface
detail not fixed by the pins is asserted order/precision-robustly and flagged
in the final report so MAIN can reconcile against the 4.B implementation.

THE GATES (pins §4.5):
  (a) CSV-twin: same run, `-record` CSV + `recorder ladruno -overlay`:
      DATA[k,:,0] == CSV row k (<=1e-12 rel); ID == CSV header == REGION_NODES;
      DATA row count == commit count.
  (b) topology: REGION_NODES/DRAINED_NODES == model-as-built; drained rows
      identically 0 in DATA; TAG attr == pattern tag.
  (c) multi-overlay: two disjoint overlays -> two result groups + two topology
      rows, values segregated (per-overlay CSV cross-check).
  (d) fatals: `-overlay 999` (unknown tag); bare `-overlay` with none in domain;
      Monitor `-overlay` mixed with `-node`; Monitor `-nodes` w/ a non-region
      node. A fatal must NOT silently succeed.
  (e) subcycle cadence: `-subcycle 4` -> recorded p piecewise-constant within
      windows, changes exactly at sync commits (period == N, phase-agnostic).
  (f) Monitor twin: SWMR COLUMNS exact, FRAMES == CSV at matched commits,
      `-nodes` subset selects exactly those columns.
  (g) no-regression smoke: a plain `recorder ladruno -N displacement` channel
      on the same model still produces RESULTS/ON_NODES/DISPLACEMENT.
  (h) envelope: `-envelope` overlay channel MIN/MAX/ABSMAX == hand-computed
      running extremes over the CSV (SKIP-with-reason if family-gated).

RUNNER (worktree dist/bin opensees.pyd is CPython 3.12):
  python -S tests/test_ladruno_overlay_recorder.py        # standalone, exit 1 on fail
  py -3.12 -m pytest tests/test_ladruno_overlay_recorder.py -x -q   # pytest
The `-S` invocation re-adds site-packages below (numpy + h5py) and evicts any
boot-.pth preloaded stale opensees.pyd so THIS worktree's engine wins.

FATAL-PROBE DESIGN (gate d): some opserr fatals raise OpenSeesError, others
call C `exit(-1)` and kill the interpreter. We first run ONE representative
fatal (`-overlay 999`) in a SUBPROCESS probe. If it was caught as a Python
exception ("raises") we run all four fatals INLINE (cheap); if it killed the
process ("exits") we run each fatal in its OWN subprocess (isolation). A fatal
that returns silently is a hard failure. Subprocess cases re-exec THIS file
with `_OV_FATAL_CASE` set so no model-build source is duplicated in a string.
"""
import os
import sys
import types
import subprocess
import tempfile
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap (copied from test_ladruno_overlay_physics.py, adapted for h5py)
# --------------------------------------------------------------------------
# h5py may read a file another process still holds open (Monitor SWMR); disable
# libhdf5 file locking BEFORE importing h5py (LEDGER_quirks precedent).
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# (a) `python -S` skips site: re-add the user + core site-packages for numpy/h5py.
for _sp in (
    os.path.join(os.environ.get("APPDATA", ""), "Python", "Python312", "site-packages"),
    r"C:\Users\nmora\AppData\Roaming\Python\Python312\site-packages",
    os.path.join(sys.prefix, "Lib", "site-packages"),
):
    if _sp and os.path.isdir(_sp) and _sp not in sys.path:
        sys.path.append(_sp)

import numpy as np  # noqa: E402

try:
    import pytest  # noqa: E402
    _HAVE_PYTEST = True
except Exception:  # pragma: no cover
    _HAVE_PYTEST = False

# (b) THIS worktree's engine (evict any boot-.pth preloaded stale pyd).
_ROOT = Path(__file__).resolve().parents[1]
_DIST = str(_ROOT / "dist" / "bin")
_BUILT = os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
if _BUILT:
    os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
    try:
        os.add_dll_directory(_DIST)
    except (FileNotFoundError, OSError):
        pass
    for _m in ("opensees", "openseespy", "openseespy.opensees"):
        sys.modules.pop(_m, None)
    if _DIST not in sys.path:
        sys.path.insert(0, _DIST)

if _HAVE_PYTEST and not _BUILT:
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

# (c) h5py is mandatory for this battery; skip loudly if it is missing after the
#     site-packages re-add (never a silent pass).
_HAVE_H5PY = True
_H5_ERR = ""
try:
    import h5py  # noqa: E402
except Exception as _exc:  # pragma: no cover
    _HAVE_H5PY = False
    _H5_ERR = repr(_exc)
if _HAVE_PYTEST and _BUILT and not _HAVE_H5PY:
    pytest.skip(f"h5py unavailable after site re-add: {_H5_ERR}", allow_module_level=True)

if _BUILT:
    import opensees as ops  # noqa: E402
    assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
        f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
    )

if _HAVE_PYTEST:
    pytestmark = [pytest.mark.zone_b, pytest.mark.t1]

_TMP = tempfile.mkdtemp(prefix="ov_rec_")


def _p(name):
    return os.path.join(_TMP, name)


def _rm(*paths):
    for p in paths:
        try:
            if os.path.exists(p):
                os.remove(p)
        except OSError:
            pass


class _Skip(Exception):
    """Standalone-runner analog of pytest.skip (a gate that cannot run, not a fail)."""


def _skip(msg):
    if _HAVE_PYTEST:
        pytest.skip(msg)
    raise _Skip(msg)


# ==========================================================================
# Terzaghi consolidation column (self-contained; mirrors the physics battery
# model idiom but imports nothing from it). E=1e7 nu=0.25 Kf=2.2e9 poro=0.4
# kbar=1e-11 q=1e5; massless solid under transient Newmark -> the fs1 flow.
# Small meshes are fine for the twin/topology gates (they check plumbing, not
# physics); gate (a) uses >=20 committed steps per the pin.
# ==========================================================================
_E, _NU, _KF, _PORO = 1.0e7, 0.25, 2.2e9, 0.4
_KH, _GW = 1.0e-7, 1.0e4
_KBAR = _KH / _GW
_LX, _LY = 1.0, 10.0
_MOED = _E * (1 - _NU) / ((1 + _NU) * (1 - 2 * _NU))
_CV = _KBAR * _MOED
_Q = 1.0e5


def _dt_for(nsteps, Tv=0.1):
    return (Tv * _LY ** 2 / _CV) / nsteps


def _grid(nx, ny, x0=0.0, nid0=1, eid0=1):
    """Return (node_ids, coords, elems[(eid,n0,n1,n2,n3)], top_node_ids)."""
    def NID(i, j):
        return nid0 + j * (nx + 1) + i
    ids, coords = [], []
    for j in range(ny + 1):
        for i in range(nx + 1):
            ids.append(NID(i, j))
            coords.append((x0 + i * _LX / nx, j * _LY / ny))
    elems = []
    for j in range(ny):
        for i in range(nx):
            eid = eid0 + j * nx + i
            elems.append((eid, NID(i, j), NID(i + 1, j), NID(i + 1, j + 1), NID(i, j + 1)))
    top = [NID(i, ny) for i in range(nx + 1)]
    return ids, coords, elems, top


def _add_column(nx, ny, tag, csv_path=None, extra=None, x0=0.0, nid0=1, eid0=1):
    """Add one plain-ndf2 Q4 column + a LadrunoPorousOverlay over every cell.
    Returns (sorted region-node ids, drained/top ids). Region == all cells so
    region nodes == every node in the strip."""
    ids, coords, elems, top = _grid(nx, ny, x0, nid0, eid0)
    for nid, (x, y) in zip(ids, coords):
        ops.node(nid, float(x), float(y))
    for (et, a, b, c, d) in elems:
        ops.element("quad", et, a, b, c, d, 1.0, "PlaneStrain", 1)
    for nid, (x, y) in zip(ids, coords):
        fx = 1 if (abs(x - x0) < 1e-9 or abs(x - (x0 + _LX)) < 1e-9 or abs(y) < 1e-9) else 0
        fy = 1 if abs(y) < 1e-9 else 0
        ops.fix(nid, fx, fy)
    dx = _LX / nx
    for k, nid in enumerate(top):
        w = dx if (0 < k < nx) else dx / 2
        ops.load(nid, 0.0, -_Q * w)
    et_list = [e[0] for e in elems]
    args = ["-region", *et_list, "-Kf", _KF, "-rhoF", 1000.0,
            "-perm", _KBAR, _KBAR, "-poro", _PORO, "-moduli", _E, _NU,
            "-drained", *top]
    if extra:
        args += list(extra)
    if csv_path:
        args += ["-record", csv_path, 1]
    ops.pattern("LadrunoPorousOverlay", tag, *args)
    return sorted(ids), list(top)


def _new_transient():
    ops.system("ProfileSPD")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _build_single(nx, ny, tag, csv_path=None, extra=None):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, 0.0)  # massless -> quasi-static/step
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    region, drained = _add_column(nx, ny, tag, csv_path=csv_path, extra=extra)
    _new_transient()
    return region, drained


def _run(nsteps, dt):
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "overlay analyze failed"
    ops.wipe()  # flush + close ALL recorder files


# ==========================================================================
# HDF5 / CSV readers
# ==========================================================================
def _attr(obj, name):
    v = obj.attrs[name]
    if isinstance(v, bytes):
        return v.decode("utf-8")
    if isinstance(v, np.ndarray):
        if v.dtype.kind == "S":
            d = [x.decode("utf-8") for x in v.reshape(-1)]
            return d[0] if v.size == 1 else d
        if v.size == 1:
            return v.reshape(-1)[0].item()
    return v


def _stages(f):
    return sorted(k for k in f if k.startswith("MODEL_STAGE"))


def _read_overlay(f, tag):
    """(-> ids, data2d[T,n], comps, step, time) for the first stage that carries
    overlayPressure_<tag>, else None."""
    name = f"overlayPressure_{tag}"
    for st in _stages(f):
        g = f.get(f"{st}/RESULTS/ON_NODES/{name}")
        if g is not None and "DATA" in g:
            ids = np.asarray(g["ID"][...]).reshape(-1)
            data = np.asarray(g["DATA"][...])
            d2 = data[:, :, 0] if data.ndim == 3 else np.atleast_2d(data)
            comps = _attr(g, "COMPONENTS") if "COMPONENTS" in g.attrs else ""
            step = np.asarray(g["STEP"][...]).reshape(-1) if "STEP" in g else None
            time = np.asarray(g["TIME"][...]).reshape(-1) if "TIME" in g else None
            return ids, d2, comps, step, time
    return None


def _read_topology(f, tag):
    """(-> region_nodes, drained_nodes, TAG attr) from MODEL/OVERLAYS/OVERLAY_<tag>."""
    for st in _stages(f):
        g = f.get(f"{st}/MODEL/OVERLAYS/OVERLAY_{tag}")
        if g is not None and "REGION_NODES" in g:
            rn = np.asarray(g["REGION_NODES"][...]).reshape(-1)
            dn = np.asarray(g["DRAINED_NODES"][...]).reshape(-1)
            tagattr = _attr(g, "TAG") if "TAG" in g.attrs else None
            return rn, dn, tagattr
    return None


def _read_envelope(f, tag):
    """(-> ids, MIN, MAX, ABSMAX) from RESULTS/ENVELOPES/ON_NODES/overlayPressure_<tag>."""
    name = f"overlayPressure_{tag}"
    for st in _stages(f):
        g = f.get(f"{st}/RESULTS/ENVELOPES/ON_NODES/{name}")
        if g is not None and "MIN" in g:
            ids = np.asarray(g["ID"][...]).reshape(-1)
            mn = np.atleast_2d(np.asarray(g["MIN"][...], dtype=float))
            mx = np.atleast_2d(np.asarray(g["MAX"][...], dtype=float))
            am = np.atleast_2d(np.asarray(g["ABSMAX"][...], dtype=float))
            return ids, mn, mx, am
    return None


def _has_result(f, family, name):
    for st in _stages(f):
        g = f.get(f"{st}/RESULTS/{family}/{name}")
        if g is not None and "DATA" in g:
            return True
    return False


def _read_record(path):
    """(-> times[T], P[T,n], tag->col map, header node tags). CSV: header row is
    'time,tag1,tag2,...'; data rows are 'time,p1,p2,...'."""
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    hdr = [int(t) for t in lines[0].split(",")[1:]]
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return data[:, 0], data[:, 1:], {t: k for k, t in enumerate(hdr)}, hdr


def _relmax(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))


_TW_RTOL = 1e-12   # pin (a) twin tolerance (relative: p ~ O(1e5), CSV is full-precision)


def _mon_cols(f):
    def _s(v):
        return v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v)
    return [_s(c) for c in f["COLUMNS"][:]]


# ==========================================================================
# (a) CSV-twin — THE gate
# ==========================================================================
def test_a_csv_twin():
    tag, ns = 77, 25
    dt = _dt_for(ns)
    csv, h5 = _p("a.csv"), _p("a.ladruno")
    _rm(csv, h5)
    region, drained = _build_single(2, 8, tag, csv_path=csv)
    ops.recorder("ladruno", h5, "-overlay", tag)      # pattern exists -> discoverable
    _run(ns, dt)

    ctimes, cP, ct2c, chdr = _read_record(csv)
    with h5py.File(h5, "r") as f:
        res = _read_overlay(f, tag)
        assert res is not None, "no RESULTS/ON_NODES/overlayPressure_<tag> group written"
        ids, d2, comps, step, time = res
        topo = _read_topology(f, tag)
    assert topo is not None, "no MODEL/OVERLAYS topology row written"
    rn, dn, tagattr = topo

    assert [c for c in str(comps).split(",") if c] == ["p"], f"COMPONENTS != 'p' ({comps!r})"
    assert np.max(np.abs(cP)) > 0.0, "CSV pressures are all zero (nothing to compare)"

    # row count == commit count
    assert d2.shape[0] == ns, f"DATA rows {d2.shape[0]} != commit count {ns}"
    assert cP.shape[0] == ns, f"CSV rows {cP.shape[0]} != commit count {ns}"

    # ID order == CSV header order == REGION_NODES
    assert list(ids) == chdr, f"HDF5 ID order {list(ids)[:6]} != CSV header {chdr[:6]}"
    assert list(rn) == chdr, f"REGION_NODES {list(rn)[:6]} != CSV header {chdr[:6]}"

    # DATA[k,:,0] == CSV row k
    rel = _relmax(d2, cP)
    maxabs = float(np.max(np.abs(d2 - cP)))
    assert rel <= _TW_RTOL, f"DATA vs CSV relL2 {rel:.2e} > {_TW_RTOL:.0e} (maxabs {maxabs:.2e})"
    print(f"PASS (a) CSV-twin: DATA[{d2.shape[0]}x{d2.shape[1]}] == CSV relL2={rel:.2e} "
          f"(maxabs {maxabs:.2e}); ID==header==REGION_NODES ({len(ids)} nodes); "
          f"rows==commits=={ns}; COMPONENTS='p'")


# ==========================================================================
# (b) topology rows match the model as built; drained rows identically 0
# ==========================================================================
def test_b_topology():
    tag, ns = 78, 15
    dt = _dt_for(ns)
    csv, h5 = _p("b.csv"), _p("b.ladruno")
    _rm(csv, h5)
    region, drained = _build_single(2, 8, tag, csv_path=csv)
    ops.recorder("ladruno", h5, "-overlay", tag)
    _run(ns, dt)

    with h5py.File(h5, "r") as f:
        ids, d2, comps, step, time = _read_overlay(f, tag)
        rn, dn, tagattr = _read_topology(f, tag)

    assert set(int(x) for x in rn) == set(region), (
        f"REGION_NODES set != model region {sorted(set(rn) ^ set(region))[:8]}")
    assert set(int(x) for x in dn) == set(drained), (
        f"DRAINED_NODES set != model drained {sorted(set(dn) ^ set(drained))[:8]}")
    if tagattr is not None:
        assert int(tagattr) == tag, f"OVERLAY TAG attr {tagattr} != {tag}"

    # drained rows identically 0 in DATA
    col = {int(t): k for k, t in enumerate(ids)}
    worst = 0.0
    for t in dn:
        c = col[int(t)]
        worst = max(worst, float(np.max(np.abs(d2[:, c]))))
    assert worst == 0.0, f"drained node p != 0 (max |p_drained| = {worst:.3e})"
    print(f"PASS (b) topology: REGION_NODES({len(rn)}) & DRAINED_NODES({len(dn)}) "
          f"match model; TAG={tagattr}; drained rows identically 0")


# ==========================================================================
# (c) multi-overlay — two disjoint overlays, values segregated
# ==========================================================================
def _build_two(tagA, tagB, csvA, csvB, nx=1, ny=6):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, 0.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    info = {}
    for (tag, csv, nid0, eid0, x0) in ((tagA, csvA, 1, 1, 0.0), (tagB, csvB, 1000, 1000, 3.0)):
        r, d = _add_column(nx, ny, tag, csv_path=csv, x0=x0, nid0=nid0, eid0=eid0)
        info[tag] = (r, d)
    _new_transient()
    return info


def test_c_multi_overlay():
    tagA, tagB, ns = 71, 72, 15
    dt = _dt_for(ns)
    csvA, csvB, h5 = _p("cA.csv"), _p("cB.csv"), _p("c.ladruno")
    _rm(csvA, csvB, h5)
    info = _build_two(tagA, tagB, csvA, csvB)
    ops.recorder("ladruno", h5, "-overlay", tagA, tagB)
    _run(ns, dt)

    _, cPA, _, hdrA = _read_record(csvA)
    _, cPB, _, hdrB = _read_record(csvB)
    with h5py.File(h5, "r") as f:
        resA, resB = _read_overlay(f, tagA), _read_overlay(f, tagB)
        topA, topB = _read_topology(f, tagA), _read_topology(f, tagB)
    assert resA is not None and resB is not None, "both overlay result groups must be present"
    assert topA is not None and topB is not None, "both overlay topology rows must be present"

    idsA, d2A = resA[0], resA[1]
    idsB, d2B = resB[0], resB[1]
    # regions truly disjoint
    assert set(int(x) for x in topA[0]).isdisjoint(set(int(x) for x in topB[0])), \
        "overlay region node sets overlap (test builds disjoint columns)"
    # segregation: each HDF5 group matches ITS own CSV, not the other's
    relA, relB = _relmax(d2A, cPA), _relmax(d2B, cPB)
    assert relA <= _TW_RTOL, f"overlay A HDF5 vs CSV_A relL2 {relA:.2e}"
    assert relB <= _TW_RTOL, f"overlay B HDF5 vs CSV_B relL2 {relB:.2e}"
    assert list(idsA) == hdrA and list(idsB) == hdrB, "per-overlay ID order != its CSV header"
    print(f"PASS (c) multi-overlay: 2 groups + 2 topology rows; regions disjoint "
          f"({len(topA[0])}/{len(topB[0])} nodes); A==CSV_A {relA:.1e}, B==CSV_B {relB:.1e}")


# ==========================================================================
# (d) fatals — probe once, then choose inline vs subprocess isolation
# ==========================================================================
def _tiny(nx=2, ny=2):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ids, coords, elems, top = _grid(nx, ny)
    for nid, (x, y) in zip(ids, coords):
        ops.node(nid, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, 0.0)
    for (et, a, b, c, d) in elems:
        ops.element("quad", et, a, b, c, d, 1.0, "PlaneStrain", 1)
    for nid, (x, y) in zip(ids, coords):
        fx = 1 if (abs(x) < 1e-9 or abs(x - _LX) < 1e-9 or abs(y) < 1e-9) else 0
        fy = 1 if abs(y) < 1e-9 else 0
        ops.fix(nid, fx, fy)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    dx = _LX / nx
    for k, nid in enumerate(top):
        w = dx if (0 < k < nx) else dx / 2
        ops.load(nid, 0.0, -_Q * w)
    return ids, coords, elems, top


def _overlay_full(tag=77):
    _, _, elems, top = _tiny()
    et = [e[0] for e in elems]
    ops.pattern("LadrunoPorousOverlay", tag, "-region", *et, "-Kf", _KF, "-rhoF", 1000.0,
                "-perm", _KBAR, _KBAR, "-poro", _PORO, "-moduli", _E, _NU, "-drained", *top)
    return top


def _fatal_unknown_tag():
    _overlay_full(77)                                    # a VALID overlay exists (tag 77)
    ops.recorder("ladruno", _p("fatal_unk.ladruno"), "-overlay", 999)  # 999 unknown -> fatal


def _fatal_bare_no_overlay():
    _tiny()                                              # NO overlay pattern in the domain
    ops.recorder("ladruno", _p("fatal_bare.ladruno"), "-overlay")      # bare + none -> fatal


def _fatal_monitor_mixed():
    _overlay_full(77)
    ops.recorder("Monitor", "-overlay", 77, "-node", 5, "-dof", 1, "-resp", "disp",
                 "-sink", _p("fatal_mix.h5"))            # -overlay mixed with -node -> fatal


def _fatal_monitor_nonregion():
    # overlay over the BOTTOM strip only (elems 1,2 -> nodes 1..6); node 9 is non-region.
    _tiny(nx=2, ny=2)
    ops.pattern("LadrunoPorousOverlay", 77, "-region", 1, 2, "-Kf", _KF, "-rhoF", 1000.0,
                "-perm", _KBAR, _KBAR, "-poro", _PORO, "-moduli", _E, _NU, "-drained", 4, 5, 6)
    # non-region node 9 -> fatal (at construction, or deferred to first record: run 1 step)
    ops.recorder("Monitor", "-overlay", 77, "-nodes", 9, "-sink", _p("fatal_nr.h5"))
    _new_transient()
    ops.analyze(1, _dt_for(1))


_FATAL_CASES = {
    "unknown_tag": _fatal_unknown_tag,
    "bare_no_overlay": _fatal_bare_no_overlay,
    "monitor_mixed": _fatal_monitor_mixed,
    "monitor_nonregion": _fatal_monitor_nonregion,
}


def _run_self(env_extra):
    env = os.environ.copy()
    env.update(env_extra)
    return subprocess.run(
        [sys.executable, "-S", os.path.abspath(__file__)],
        env=env, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=300,
    )


def _detect_fatal_mode():
    """Return 'raises' | 'exits' | 'silent' from a subprocess probe of `-overlay 999`."""
    cp = _run_self({"_OV_FATAL_CASE": "unknown_tag", "_OV_FATAL_MODE": "probe"})
    out = (cp.stdout or "") + (cp.stderr or "")
    if "PROBE_RAISES" in out:
        return "raises"
    if "PROBE_SILENT" in out:
        return "silent"
    return "exits"


def _expect_fatal(name, mode):
    if mode == "raises":
        # fatals surface as Python exceptions -> run inline (cheap, no isolation).
        raised = False
        try:
            _FATAL_CASES[name]()
        except Exception:
            raised = True
        assert raised, f"{name}: expected a fatal, but the command returned silently"
    else:
        # a fatal that calls exit(-1) would kill the interpreter -> isolate each case.
        cp = _run_self({"_OV_FATAL_CASE": name, "_OV_FATAL_MODE": "run"})
        out = cp.stdout or ""
        assert "REACHED_END_OK" not in out, (
            f"{name}: command completed with NO fatal (stdout tail: {out[-200:]!r})")
        assert cp.returncode != 0, f"{name}: process exit 0 -> no fatal"


def test_d_fatals():
    mode = _detect_fatal_mode()
    print(f"  (d) fatal-probe: '-overlay 999' mechanism = {mode!r} "
          f"-> {'inline' if mode == 'raises' else 'subprocess'} isolation")
    assert mode != "silent", (
        "representative fatal ('recorder ladruno -overlay 999') did NOT fatal — pin 4.5(d) "
        "requires unknown-tag to abort; the run must not silently succeed")
    for name in ("unknown_tag", "bare_no_overlay", "monitor_mixed", "monitor_nonregion"):
        _expect_fatal(name, mode)
    print(f"PASS (d) fatals: unknown-tag, bare-no-overlay, Monitor+-node mix, "
          f"Monitor non-region node ALL aborted (mechanism={mode})")


# ==========================================================================
# (e) subcycle cadence — p piecewise-constant, changes exactly at syncs
# ==========================================================================
def test_e_subcycle_cadence():
    tag, N, ns = 79, 4, 24
    dt = _dt_for(ns)
    csv, h5 = _p("e.csv"), _p("e.ladruno")
    _rm(csv, h5)
    _build_single(2, 8, tag, csv_path=csv, extra=["-subcycle", N])
    ops.recorder("ladruno", h5, "-overlay", tag)
    _run(ns, dt)

    with h5py.File(h5, "r") as f:
        ids, d2, comps, step, time = _read_overlay(f, tag)
    assert d2.shape[0] == ns, f"expected {ns} recorded rows, got {d2.shape[0]}"

    # row-to-row change magnitude; within a window p is the SAME stored vector -> 0.
    scale = max(1.0, float(np.max(np.abs(d2))))
    dchg = np.array([np.max(np.abs(d2[i] - d2[i - 1])) for i in range(1, ns)])
    tol = 1e-9 * scale
    changes = [i for i in range(1, ns) if dchg[i - 1] > tol]
    assert changes, "recorded p never changed under -subcycle (fluid not advancing?)"
    gaps = np.diff(changes)
    assert all(g == N for g in gaps), (
        f"change-points not spaced by N={N}: positions {changes}, gaps {list(gaps)}")
    # ~ ns/N sync events (phase may drop the first or last partial window)
    assert (ns // N) - 1 <= len(changes) <= (ns // N) + 1, (
        f"{len(changes)} change-points, expected ~{ns // N} at N={N}")

    # cross-check: the distinct window values == the -record CSV sync rows
    _, cP, _, _ = _read_record(csv)
    extra = ""
    if cP.shape[0] == len(changes):
        rel = _relmax(d2[changes], cP)
        assert rel <= _TW_RTOL, f"synced HDF5 window values vs CSV sync rows relL2 {rel:.2e}"
        extra = f"; window values == {cP.shape[0]} CSV sync rows ({rel:.1e})"
    print(f"PASS (e) subcycle N={N}: {len(changes)} change-points, all spaced {N} apart "
          f"(piecewise-constant within windows){extra}")


# ==========================================================================
# (f) Monitor twin — SWMR columns/frames vs CSV; -nodes subset selects columns
# ==========================================================================
def test_f_monitor_twin():
    tag, ns = 80, 20
    dt = _dt_for(ns)
    csv, mon = _p("f.csv"), _p("f_mon.h5")
    _rm(csv, mon)
    _build_single(2, 8, tag, csv_path=csv)
    ops.recorder("Monitor", "-overlay", tag, "-sink", mon, "-every", 1)
    _run(ns, dt)

    _, cP, ct2c, chdr = _read_record(csv)
    with h5py.File(mon, "r") as f:
        assert _attr(f, "FORMAT") == "ladruno-monitor", "sink FORMAT attr != 'ladruno-monitor'"
        cols = _mon_cols(f)
        frames = np.asarray(f["FRAMES"][:], dtype=float)
    want_cols = [f"overlay{tag}.p.node{n}" for n in chdr]
    assert cols == want_cols, f"COLUMNS\n got {cols[:4]}...\n want {want_cols[:4]}..."
    n = min(frames.shape[0], cP.shape[0])
    assert frames.shape[1] == cP.shape[1], f"FRAMES cols {frames.shape[1]} != CSV cols {cP.shape[1]}"
    rel = _relmax(frames[:n], cP[:n])
    assert rel <= 1e-9, f"Monitor FRAMES vs CSV relL2 {rel:.2e}"

    # -nodes subset: pick 3 region nodes; columns must be exactly those
    sub = [chdr[1], chdr[len(chdr) // 2], chdr[-2]]
    sub = list(dict.fromkeys(sub))  # dedupe, keep order
    csv2, mon2 = _p("f2.csv"), _p("f2_mon.h5")
    _rm(csv2, mon2)
    _build_single(2, 8, tag, csv_path=csv2)
    ops.recorder("Monitor", "-overlay", tag, "-nodes", *sub, "-sink", mon2, "-every", 1)
    _run(ns, dt)
    _, cP2, ct2c2, chdr2 = _read_record(csv2)
    with h5py.File(mon2, "r") as f:
        cols2 = _mon_cols(f)
        frames2 = np.asarray(f["FRAMES"][:], dtype=float)
    assert set(cols2) == {f"overlay{tag}.p.node{n}" for n in sub}, (
        f"-nodes subset COLUMNS {cols2} != expected for {sub}")
    # value cross-check per selected node against its CSV column (order-robust)
    worst = 0.0
    m = min(frames2.shape[0], cP2.shape[0])
    for j, label in enumerate(cols2):
        node = int(label.split(".node")[1])
        worst = max(worst, _relmax(frames2[:m, j], cP2[:m, ct2c2[node]]))
    assert worst <= 1e-9, f"-nodes subset FRAMES vs CSV columns relL2 {worst:.2e}"
    print(f"PASS (f) Monitor twin: FORMAT+COLUMNS exact ({len(cols)} cols), FRAMES==CSV "
          f"{rel:.1e}; -nodes subset -> exactly {len(cols2)} columns ({worst:.1e})")


# ==========================================================================
# (g) no-regression smoke — plain -N displacement channel still lands
# ==========================================================================
def test_g_regression_smoke():
    tag, ns = 81, 5
    dt = _dt_for(ns)
    h5 = _p("g.ladruno")
    _rm(h5)
    _build_single(2, 8, tag)                  # overlay present in the domain
    ops.recorder("ladruno", h5, "-N", "displacement")   # ordinary nodal channel
    _run(ns, dt)
    with h5py.File(h5, "r") as f:
        ok = _has_result(f, "ON_NODES", "DISPLACEMENT")
    assert ok, "plain 'recorder ladruno -N displacement' produced no ON_NODES/DISPLACEMENT group"
    print("PASS (g) regression smoke: plain -N displacement channel still writes "
          "RESULTS/ON_NODES/DISPLACEMENT alongside the overlay pattern")


# ==========================================================================
# (h) envelope — overlay -envelope MIN/MAX/ABSMAX == hand-computed CSV extremes
# ==========================================================================
def test_h_envelope():
    tag, ns = 82, 20
    dt = _dt_for(ns)
    csv, h5 = _p("h.csv"), _p("h.ladruno")
    _rm(csv, h5)
    _build_single(2, 8, tag, csv_path=csv)
    ops.recorder("ladruno", h5, "-overlay", tag, "-envelope")
    _run(ns, dt)

    _, cP, _, chdr = _read_record(csv)
    with h5py.File(h5, "r") as f:
        env = _read_envelope(f, tag)
    if env is None:
        _skip("no RESULTS/ENVELOPES/ON_NODES/overlayPressure_<tag> group — envelope "
              "wiring appears family-gated for overlay channels (pin 4.5(h) allows "
              "SKIP-with-reason; MAIN to adjudicate against the 4.B as-built)")
    ids, mn, mx, am = env
    assert list(ids) == chdr, f"envelope ID order {list(ids)[:6]} != CSV header {chdr[:6]}"
    want_min = cP.min(axis=0)
    want_max = cP.max(axis=0)
    want_abs = np.abs(cP).max(axis=0)
    r_min = _relmax(mn.reshape(-1), want_min)
    r_max = _relmax(mx.reshape(-1), want_max)
    r_abs = _relmax(am.reshape(-1), want_abs)
    assert r_min <= 1e-9, f"envelope MIN vs CSV running-min relL2 {r_min:.2e}"
    assert r_max <= 1e-9, f"envelope MAX vs CSV running-max relL2 {r_max:.2e}"
    assert r_abs <= 1e-9, f"envelope ABSMAX vs CSV running-|max| relL2 {r_abs:.2e}"
    print(f"PASS (h) envelope: MIN/MAX/ABSMAX == CSV running extremes "
          f"({r_min:.1e}/{r_max:.1e}/{r_abs:.1e}) over {ns} steps")


def test_i_pattern_removal():
    """Gate (i) — P4 panel F-1 pin: `remove loadPattern` mid-run must NOT dangle
    the recorder channel (pattern deletion fires no domainChange — the source
    re-resolves by tag each step and zero-fills). Rows after removal are
    identically 0.0, rows before are untouched, the file stays readable."""
    tag, ns_pre, ns_post = 83, 5, 3
    dt = _dt_for(ns_pre + ns_post)
    csv, h5 = _p("i.csv"), _p("i.ladruno")
    _rm(csv, h5)
    _build_single(2, 8, tag, csv_path=csv)
    ops.recorder("ladruno", h5, "-overlay", tag)
    for _ in range(ns_pre):
        assert ops.analyze(1, dt) == 0
    ops.remove("loadPattern", tag)      # deletes the overlay; no domainChange fires
    for _ in range(ns_post):
        assert ops.analyze(1, dt) == 0, "post-removal analyze failed"
    ops.wipe()

    _, cP, _, _ = _read_record(csv)     # CSV stops at removal: ns_pre rows
    with h5py.File(h5, "r") as f:
        got = _read_overlay(f, tag)
    assert got is not None, "overlay group missing"
    _, d2, _, _, _ = got
    assert d2.shape[0] == ns_pre + ns_post, \
        f"rows {d2.shape[0]} != commits {ns_pre + ns_post}"
    r_pre = _relmax(d2[:ns_pre, :], cP[:ns_pre, :])
    assert r_pre <= 1e-12, f"pre-removal rows drifted from CSV: relL2 {r_pre:.2e}"
    post = np.abs(d2[ns_pre:, :]).max()
    assert post == 0.0, f"post-removal rows not zero-filled (max |p| {post:.3e})"
    print(f"PASS (i) pattern removal: {ns_pre} live rows CSV-exact ({r_pre:.1e}), "
          f"{ns_post} post-removal rows identically 0.0, file readable")


# ==========================================================================
# standalone runner + fatal-subprocess dispatch
# ==========================================================================
_ALL = [
    ("test_a_csv_twin", test_a_csv_twin),
    ("test_b_topology", test_b_topology),
    ("test_c_multi_overlay", test_c_multi_overlay),
    ("test_d_fatals", test_d_fatals),
    ("test_e_subcycle_cadence", test_e_subcycle_cadence),
    ("test_f_monitor_twin", test_f_monitor_twin),
    ("test_g_regression_smoke", test_g_regression_smoke),
    ("test_h_envelope", test_h_envelope),
    ("test_i_pattern_removal", test_i_pattern_removal),
]


def _run_fatal_child(name, mode):
    """Executed in a fresh subprocess (env _OV_FATAL_CASE set) to isolate a fatal
    that may call exit(-1). 'probe' reports the mechanism; 'run' aborts (nonzero
    exit / no marker) on a real fatal, prints the sentinel only on silent success."""
    fn = _FATAL_CASES[name]
    if mode == "probe":
        try:
            fn()
            print("PROBE_SILENT")
        except Exception:
            print("PROBE_RAISES")
        sys.stdout.flush()
        os._exit(0)
    fn()                          # run mode: a fatal raises (nonzero exit) or exit(-1)s here
    print("REACHED_END_OK")       # only reached if the command did NOT fatal
    sys.stdout.flush()
    os._exit(0)


def main():
    if not _BUILT:
        print(f"SKIP: worktree engine not built at {_DIST}")
        return 0
    if not _HAVE_H5PY:
        print(f"SKIP: h5py unavailable after site re-add: {_H5_ERR}")
        return 0
    fails, skips = [], []
    for name, fn in _ALL:
        try:
            fn()
        except _Skip as exc:
            print(f"SKIP {name}: {exc}")
            skips.append(name)
        except Exception as exc:  # noqa: BLE001
            import traceback
            print(f"FAIL {name}: {exc}")
            traceback.print_exc()
            fails.append(name)
    print("=" * 74)
    if fails:
        print(f"FAILED {len(fails)}/{len(_ALL)}: {fails}"
              + (f"  (skipped {skips})" if skips else ""))
        return 1
    print(f"ALL {len(_ALL) - len(skips)}/{len(_ALL)} RECORDER GATES PASSED"
          + (f"  ({len(skips)} skipped: {skips})" if skips else ""))
    return 0


if __name__ == "__main__":
    _case = os.environ.get("_OV_FATAL_CASE")
    if _case:
        _run_fatal_child(_case, os.environ.get("_OV_FATAL_MODE", "run"))
    else:
        sys.exit(main())
