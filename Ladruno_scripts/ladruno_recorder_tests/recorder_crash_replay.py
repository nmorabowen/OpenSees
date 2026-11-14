"""Crash-consistency replay model of the LadrunoRecorder streaming write path.

What it shows: whether a ``.ladruno`` written by ``recorder ladruno`` (time-series mode, the
default ``StreamingSink``) survives a kill at EVERY point the process can die, and what a
killed run keeps. It is the recorder's twin of
``Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py`` (the MPM particle writer) and reuses
its logging file object. It does NOT run the C++ recorder: it mirrors its HDF5 call pattern in
h5py, through the h5py ``fileobj`` driver:

- ``LadrunoRecorder::initialize``: ``H5Fcreate`` with link creation order tracked + indexed
  and the DEFAULT file access list (no libver bounds, so v1 B-tree chunk indexes; no SWMR,
  ``LADRUNO_USE_SWMR`` is defined nowhere); ``INFO`` + ``INFO/PROVENANCE``;
- ``writeModel`` (abridged): ``MODEL_STAGE[1]/MODEL/NODES/{ID,COORDINATES}`` contiguous;
- ``StreamingSink::begin`` per result: ``RESULTS/ON_NODES/<name>`` with its attributes, a
  contiguous ``ID``, and ``DATA`` [0 x nIds x nComp] chunked by ``planTimeSeries`` (~1 MiB
  chunks; shuffle + deflate at ``-compress``, default 1; a dataset chunk cache sized to a row
  of chunks, w0 = 1), ``TIME`` f64 / ``STEP`` i32 chunked 1024, unfiltered, default cache;
- ``StreamingSink::accept`` per step: extend + write one slab of ``DATA``, then ``TIME``,
  ``STEP``; handles stay open for the stage;
- ``LadrunoRecorder::record``: ``H5Fflush`` when the ``-flush`` cadence is due. The first
  recorded step always flushes (``last_flush`` starts at the epoch). The cadence is wall
  clock (default 10 s); here it is ``--flush-every K`` recorded steps (0 = every step, the
  ``-flush 0`` behaviour).

Every byte HDF5 hands the driver is logged in order; the file is then rebuilt as it stands
after each write (the state a TerminateProcess leaves) and read the way a reader does: open
plainly, read every result's DATA / TIME / STEP in full, and compare them EXACTLY with what
was written. Outcomes per kill point:

- ``OK``: opens, every result has DATA, TIME, STEP of one length n, all equal to steps 1..n,
  and n is at least the steps of the last COMPLETED flush (nothing flushed was lost);
- ``RAGGED``: values are right, but DATA / TIME / STEP of a result disagree in length (the
  format validator rejects it; a reader must cut to the shortest);
- ``LOST``: readable and right, but fewer steps than the last completed flush committed;
- ``MISMATCH``: readable, but a value differs from what was written (silent corruption);
- ``READ_FAIL`` / ``OPEN_FAIL``: h5py raises.

``--chunk-scale s`` multiplies the plan's 1 MiB chunk target and 32 MiB row budget (and the
1024-element TIME / STEP chunk) by s, so a chunk index splits after s times fewer steps; the
mechanism is unchanged (the MPM finding: bigger chunks only make splits rarer). Use 1 for the
real layout. ``--check-from`` / ``--check-to`` restrict the replay to a window of steps (the
image at the window start is snapshotted, so a long run is cheap to check near one split).

Run from anywhere (no OpenSees binary needed):

    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 400 --flush-every 0
    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 400 --flush-every 20
    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 400 --flush-every 0 --swmr
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "Ladruno_scripts" / "mpm_crash_tests"))
import mpm_contracts as mc
from hdf5_crash_replay import LogFile

STAGE = "MODEL_STAGE[1]"


# ---------------------------------------------------------------- layout (mirrors C++)
def plan_time_series(n_ids, n_comp, elem_bytes, scale=1.0):
    """Ladruno_Hdf5.h planTimeSeries, with the two byte budgets scaled by ``scale``."""
    target = max(1, int((1 << 20) * scale))
    row_budget = max(1, int((32 << 20) * scale))
    row_bytes = n_comp * elem_bytes
    slab = n_ids * row_bytes
    if slab <= target:
        chunk_ids = n_ids
        ct = max(1, min(1024, target // slab if slab else 1))
    else:
        ct = max(1, min(16, row_budget // slab))
        chunk_ids = max(1, min(n_ids, target // (ct * row_bytes)))
    tiles = (n_ids + chunk_ids - 1) // chunk_ids
    chunk_bytes = ct * chunk_ids * row_bytes
    cache = tiles * chunk_bytes + chunk_bytes
    cache = max(1 << 20, min(256 << 20, cache))
    slots = (tiles * 10 + 521) | 1
    return ct, chunk_ids, cache, slots


def parse_results(spec):
    """'DISPLACEMENT:3,VELOCITY:3' -> [(name, ncomp)]"""
    out = []
    for tok in spec.split(","):
        name, nc = tok.split(":")
        out.append((name.strip(), int(nc)))
    return out


def values(ri, k, n_ids, nc):
    """Slab of result ri at step k (1-based): full-mantissa floats, deterministic."""
    i = np.arange(n_ids, dtype=np.float64)[:, None]
    c = np.arange(nc, dtype=np.float64)[None, :]
    return np.sin(0.0137 * k + 0.113 * i + 1.31 * c + 0.7 * ri) * (1e-3 * (ri + 1))


def step_time(k):
    return k * 1e-4


class NoLog:
    """Stand-in for LogFile when writing a real file (the kill stress)."""
    buf = b""
    log = ()

    def mark(self, label):
        pass


def write_run(fobj, a, results, real_path=None, flush_due=None):
    """Mirror of the recorder's HDF5 calls. ``real_path``: write that file with the
    default (sec2) driver instead of the logging fileobj. ``flush_due(k)``: overrides the
    step-count cadence (the kill stress passes a wall-clock one, as the recorder has)."""
    n_ids = a.nodes
    fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    if real_path is None:
        fapl.set_fileobj_driver(h5py.h5fd.fileobj_driver, fobj)
    if a.swmr:
        fapl.set_libver_bounds(h5py.h5f.LIBVER_LATEST, h5py.h5f.LIBVER_LATEST)
    fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
    name = os.fsencode(real_path) if real_path else b"replay.ladruno"
    fid = h5py.h5f.create(name, h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl)
    f = h5py.File(fid)

    # initialize(): INFO (+ PROVENANCE)
    info = f.create_group("INFO", track_order=True)
    for k_, v in (("GENERATOR", "Ladruno"), ("FORMAT_VERSION", 1), ("FORMAT_MINOR_VERSION", 1),
                  ("SOLVER_NAME", "OpenSees"), ("SOLVER_VERSION", np.array([3, 5, 1], "<i4")),
                  ("SPATIAL_DIM", 3), ("STORED_PRECISION", a.precision), ("PARTITION_ID", 0),
                  ("PARTITIONED", 0), ("NUM_PARTITIONS", 1),
                  ("RUN_ID", "00000000-0000-4000-8000-000000000000"), ("RUN_ID_SCOPE", "process")):
        info.attrs[k_] = v
    info.create_group("PROVENANCE", track_order=True)

    # writeModel() (abridged: the node set; elements add more contiguous datasets)
    st = f.create_group(STAGE, track_order=True)
    st.attrs["STEP"] = 0
    st.attrs["TIME"] = 0.0
    st.attrs["KIND"] = "transient"
    model = st.create_group("MODEL", track_order=True)
    nodes = model.create_group("NODES", track_order=True)
    ids = np.arange(1, n_ids + 1, dtype="<i4")
    nodes.create_dataset("ID", data=ids.reshape(-1, 1), track_times=True)
    nodes.create_dataset("COORDINATES", data=np.random.default_rng(1).random((n_ids, 3)),
                         track_times=True)

    # StreamingSink::begin() per result
    on_nodes = st.create_group("RESULTS", track_order=True).create_group("ON_NODES", track_order=True)
    disk = "<f4" if a.precision == "f32" else "<f8"
    eb = 4 if a.precision == "f32" else 8
    axis_chunk = max(1, int(1024 * a.chunk_scale)) if a.scale_axis else 1024
    sinks = []
    for ri, (name, nc) in enumerate(results):
        g = on_nodes.create_group(name, track_order=True)
        for at in ("DISPLAY_NAME", "COMPONENTS", "DIMENSION", "DESCRIPTION"):
            g.attrs[at] = name
        g.attrs["TYPE"] = 0
        g.attrs["DATA_TYPE"] = 1
        g.attrs["NUM_COMPONENTS"] = nc
        g.attrs["PARTITION_REDUCTION"] = "NONE"
        g.create_dataset("ID", data=ids.reshape(-1, 1), track_times=True)
        ct, cids, cache, slots = plan_time_series(n_ids, nc, eb, a.chunk_scale)
        kw = {"compression": "gzip", "compression_opts": a.compress, "shuffle": True} if a.compress > 0 else {}
        d = g.create_dataset("DATA", shape=(0, n_ids, nc), maxshape=(None, n_ids, nc),
                             chunks=(ct, cids, nc), dtype=disk, track_times=True,
                             rdcc_nbytes=cache, rdcc_nslots=slots, rdcc_w0=1.0, **kw)
        dt = g.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(axis_chunk,),
                              dtype="<f8", track_times=True)
        ds = g.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(axis_chunk,),
                              dtype="<i4", track_times=True)
        sinks.append((ri, name, nc, d, dt, ds))
        if a.verbose:
            print(f"  {name}: DATA chunk ({ct}, {cids}, {nc}) cache {cache} B; TIME/STEP chunk {axis_chunk}")
    if a.swmr:
        f.flush()
        f.swmr_mode = True        # the fix under study: arm after the stage layout
    fobj.mark("stage-open")

    snapshot = None
    for k in range(1, a.steps + 1):
        if k == a.check_from:
            snapshot = (bytes(fobj.buf), len(fobj.log))
        fobj.mark(f"s{k}:write")
        t = k - 1
        for ri, name, nc, d, dt, ds in sinks:          # StreamingSink::accept
            d.resize(t + 1, axis=0)
            d[t] = values(ri, k, n_ids, nc)
            dt.resize(t + 1, axis=0)
            dt[t] = step_time(k)
            ds.resize(t + 1, axis=0)
            ds[t] = k
        if flush_due is not None:
            due = flush_due(k)
        else:
            due = (k == 1) or a.flush_every <= 0 or (k - 1) % a.flush_every == 0
        if due:
            fobj.mark(f"s{k}:flush")
            f.flush()
            fobj.mark(f"s{k}:flushed")
        if a.check_to and k >= a.check_to:
            break
    fobj.mark("end")
    # no close: the run is killed somewhere inside the logged window
    return snapshot


# ---------------------------------------------------------------- reader
def reference(ri, n_ids, nc, n_max, precision):
    """What result ri holds after n_max steps, as read back (f32 rounded if stored f32)."""
    arr = np.stack([values(ri, k, n_ids, nc) for k in range(1, n_max + 1)]) if n_max else \
        np.zeros((0, n_ids, nc))
    if precision == "f32":
        arr = arr.astype("<f4").astype(np.float64)
    return arr


def expected(results, n_ids, n_max, precision):
    return {name: reference(ri, n_ids, nc, n_max, precision) for ri, (name, nc) in enumerate(results)}


def check(path, results, exp, swmr, precision="f64"):
    """Returns (status, n_min, detail)."""
    try:
        f = h5py.File(path, "r")
    except Exception as e:  # noqa: BLE001
        if not swmr:
            return "OPEN_FAIL", -1, str(e).splitlines()[0][:160]
        try:
            if "already open for write" in str(e):
                mc.set_superblock(path, flags=0x05)   # sec2's SWMR flags; fileobj writes 0x04
            f = h5py.File(path, "r", swmr=True)
        except Exception as e2:  # noqa: BLE001
            return "OPEN_FAIL", -1, (str(e).splitlines()[0][:100] + " | swmr: " + str(e2).splitlines()[0][:100])
    with f:
        ns = []
        ragged = False
        try:
            base = f[f"{STAGE}/RESULTS/ON_NODES"]
            for ri, (name, nc) in enumerate(results):
                g = base[name]
                try:
                    data = g["DATA"][...]
                except Exception as e:  # noqa: BLE001
                    # how much is still readable: whole chunk rows from the start
                    dd = g["DATA"]
                    ct = dd.chunks[0]
                    ok = 0
                    for j in range(0, dd.shape[0], ct):
                        try:
                            dd[j:j + ct]
                        except Exception:  # noqa: BLE001
                            break
                        ok = min(j + ct, dd.shape[0])
                    return "READ_FAIL", -1, (f"{name}/DATA: {str(e).splitlines()[0][:90]}; "
                                             f"readable rows {ok} of {dd.shape[0]} (chunk {ct} steps)")
                tim = g["TIME"][...]
                stp = g["STEP"][...]
                n = min(len(data), len(tim), len(stp))
                if not (len(data) == len(tim) == len(stp)):
                    ragged = True
                if exp is not None:
                    ref = exp[name]
                else:   # the kill stress: build the reference for the steps found
                    ref = reference(ri, data.shape[1], nc, len(data), precision)
                if n > len(ref):
                    return "MISMATCH", n, f"{name}: {n} steps > written"
                if not np.array_equal(np.asarray(data[:n], dtype=np.float64), ref[:n]):
                    bad = int(np.argmax([not np.array_equal(np.asarray(data[j], np.float64), ref[j])
                                         for j in range(n)]))
                    return "MISMATCH", n, f"{name}/DATA step {bad + 1} of {n}"
                if not np.array_equal(stp[:n], np.arange(1, n + 1)) or \
                        not np.array_equal(tim[:n], step_time(np.arange(1, n + 1))):
                    return "MISMATCH", n, f"{name}/TIME|STEP"
                ns.append(n)
        except Exception as e:  # noqa: BLE001
            return "READ_FAIL", -1, str(e).splitlines()[0][:160]
    n_min = min(ns) if ns else 0
    return ("RAGGED" if ragged else "OK"), n_min, (f"lengths {ns}" if ragged else "")


# ---------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--nodes", type=int, default=200)
    ap.add_argument("--results", default="DISPLACEMENT:3,VELOCITY:3")
    ap.add_argument("--flush-every", type=int, default=0,
                    help="flush every K recorded steps (0 = every step, -flush 0); first step always")
    ap.add_argument("--compress", type=int, default=1, help="-compress level (default 1; 0 = no filter)")
    ap.add_argument("--precision", choices=("f64", "f32"), default="f64")
    ap.add_argument("--chunk-scale", type=float, default=1.0,
                    help="scale of the 1 MiB chunk target / 32 MiB row budget (1 = real layout)")
    ap.add_argument("--scale-axis", action="store_true",
                    help="also scale the 1024-element TIME/STEP chunk by --chunk-scale")
    ap.add_argument("--check-from", type=int, default=1, help="first step whose kill points are checked")
    ap.add_argument("--check-to", type=int, default=0, help="stop the run after this step (0 = --steps)")
    ap.add_argument("--swmr", action="store_true", help="latest format + SWMR after the layout (the fix)")
    ap.add_argument("--show", type=int, default=15)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    results = parse_results(a.results)

    t0 = time.time()
    fobj = LogFile()
    print(f"h5py {h5py.__version__} / HDF5 {h5py.version.hdf5_version}; nodes={a.nodes} "
          f"results={a.results} steps={a.steps} flush_every={a.flush_every} compress={a.compress} "
          f"precision={a.precision} chunk_scale={a.chunk_scale} scale_axis={a.scale_axis} swmr={a.swmr}")
    snapshot = write_run(fobj, a, results)
    log = fobj.log
    n_steps = a.check_to or a.steps
    if snapshot is None:
        snapshot = (b"", 0)
    image = bytearray(snapshot[0])
    start = snapshot[1]
    nw = sum(1 for e in log[start:] if e[0] == "w")
    print(f"  wrote in {time.time() - t0:.1f} s; ops in window={len(log) - start} writes={nw} "
          f"image at window start={len(image)} B, final={len(fobj.buf)} B")

    exp = expected(results, a.nodes, n_steps, a.precision)
    tmp = tempfile.mkdtemp(prefix="rec_replay_")
    path = os.path.join(tmp, "k.ladruno")
    phase = "start"
    # steps committed by the last COMPLETED flush (before the window: the cadence's last one)
    flushed = max([k for k in range(1, a.check_from)
                   if k == 1 or a.flush_every <= 0 or (k - 1) % a.flush_every == 0], default=0)
    completed = a.check_from - 1   # steps whose accept() finished
    stats = collections.Counter()
    by_phase = collections.Counter()
    fail_steps = collections.defaultdict(set)
    fails = []
    lag = collections.Counter()
    for i in range(start, len(log)):
        e = log[i]
        if e[0] == "m":
            phase = e[1]
            if phase.startswith("s"):
                k = int(phase[1:].split(":")[0])
                kind = phase.split(":")[1]
                if kind == "flushed":
                    flushed = k
                if kind in ("flush",) or (kind == "write" and k > 1):
                    completed = k if kind == "flush" else k - 1
            continue
        if e[0] == "w":
            _, off, data = e
            end = off + len(data)
            if end > len(image):
                image.extend(b"\0" * (end - len(image)))
            image[off:end] = data
        else:
            size = e[1]
            if size < len(image):
                del image[size:]
            else:
                image.extend(b"\0" * (size - len(image)))
        if not phase.startswith("s"):
            continue
        with open(path, "wb") as fh:
            fh.write(image)
        status, n, detail = check(path, results, exp, a.swmr)
        if flushed == 0 and status == "OPEN_FAIL":
            status = "NO_FILE"   # before the first flush ends: nothing was committed yet
        if status in ("OK", "RAGGED") and n < flushed:
            status, detail = "LOST", f"read {n} < flushed {flushed}"
        stats[status] += 1
        kind = phase.split(":")[1]
        by_phase[(kind, status)] += 1
        if status in ("OK", "NO_FILE"):
            lag[completed - max(n, 0)] += 1
        else:
            fail_steps[status].add(int(phase[1:].split(":")[0]))
            fails.append((i, phase, status, n, flushed, completed, detail))
    n_fl = sum(1 for e in log[start:] if e[0] == "m" and e[1].endswith(":flush"))
    # the first flush (step 1) writes the file's first superblock: a kill before it ends
    # leaves no openable file but loses nothing committed, so it is reported apart
    bad_fl = len({int(ph[1:].split(":")[0]) for _, ph, *_ in fails
                  if ph.endswith(":flush") and ph != "s1:flush"})
    in_fl = sum(c for (kind, st_), c in by_phase.items() if kind == "flush" and st_ != "NO_FILE")
    bad_in_fl = sum(c for (kind, st_), c in by_phase.items()
                    if kind == "flush" and st_ not in ("OK", "NO_FILE"))
    print(f"  replayed in {time.time() - t0:.1f} s")
    print(f"flushes in window: {n_fl}; flushes with a failing kill point: {bad_fl}; "
          f"kill points inside flushes: {in_fl}, failing: {bad_in_fl}")
    print("kill points checked:", sum(stats.values()), dict(stats))
    print("by phase (the mark before the kill point: write = between flushes, inside accept(); "
          "flush = inside H5Fflush; flushed = after the flush, before the next step):")
    for (kind, status), c in sorted(by_phase.items()):
        print(f"  {kind:8s} {status:10s} {c}")
    for status, ks in sorted(fail_steps.items()):
        ks = sorted(ks)
        print(f"  steps with a {status} kill point ({len(ks)}): {ks[:40]}{' ...' if len(ks) > 40 else ''}")
    print("OK points: completed steps minus steps read:", dict(sorted(lag.items())[:12]),
          "..." if len(lag) > 12 else "")
    for i, phase, status, n, fl, done, detail in fails[: a.show]:
        print(f"  FAIL op#{i} after-mark={phase} {status} read={n} flushed={fl} completed={done} {detail}")


if __name__ == "__main__":
    main()
