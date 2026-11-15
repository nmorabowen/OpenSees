"""WP-233 (#1082): file size and write cost of the candidate streaming layouts under SWMR.

SWMR never reuses file space that is freed while it is armed. The WP-164 layout deflates a
partially filled chunk at every flush; its compressed size grows, so HDF5 moves it and the old
copy is freed and lost. This bench writes the same run, on a real file (sec2 driver), with the
recorder's call pattern (``recorder_crash_replay.write_run`` semantics: extend + write one slab
per step, ``H5Fflush`` on a cadence of K recorded steps), for each candidate layout:

- ``base``      WP-164 layout, deflate ``-compress 1`` + shuffle, no SWMR (today's file);
- ``swmr``      the same layout, latest format, SWMR armed after the layout (the bloat);
- ``raw``       SWMR, ``DATA`` unfiltered (same chunk plan);
- ``rawNN``     as ``raw`` with the plan's 1 MiB chunk target scaled to NN KiB (an unfiltered
                chunk is allocated in full, so the last, partial chunk of every DATA costs up
                to one whole chunk: a short run with many results pays that per result);
- ``step``      SWMR, deflate, one step per chunk along T (ct = 1, ids tiled to <= 1 MiB):
                every chunk is complete when written, so it is written once;
- ``edge``      SWMR, deflate + ``H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS`` (HDF5 >= 1.10): the
                partial last chunk is kept unfiltered and deflated once, when it is complete.

Reports the file size, the write wall time (steps + flushes) and the time to read one id's
whole history (the post-processing access the WP-164 stacking optimises).

    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_swmr_layout_bench.py --nodes 200 --steps 300 --flush-every 1
"""
from __future__ import annotations

import argparse
import ctypes
import os
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from recorder_crash_replay import plan_time_series, values  # noqa: E402

_H5 = ctypes.CDLL(os.path.join(os.path.dirname(h5py.__file__), "hdf5.dll")) \
    if os.name == "nt" else ctypes.CDLL(h5py.h5.__file__)
_H5.H5Pset_chunk_opts.argtypes = (ctypes.c_int64, ctypes.c_uint)
_H5.H5Pset_chunk_opts.restype = ctypes.c_int
DONT_FILTER_PARTIAL = 0x0002


def make_data(g, n_ids, nc, layout, compress):
    eb = 8
    scale = 1.0
    if layout.startswith("raw") and layout != "raw":
        scale = int(layout[3:]) / 1024.0           # rawNN: chunk target NN KiB
    ct, cids, cache, slots = plan_time_series(n_ids, nc, eb, scale)
    if layout == "step":
        ct = 1
        cids = max(1, min(n_ids, (1 << 20) // (nc * eb)))
    dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
    dcpl.set_chunk((ct, cids, nc))
    if layout in ("base", "swmr", "step", "edge") and compress > 0:
        dcpl.set_shuffle()
        dcpl.set_deflate(compress)
    if layout == "edge":
        if _H5.H5Pset_chunk_opts(dcpl.id, DONT_FILTER_PARTIAL) < 0:
            raise RuntimeError("H5Pset_chunk_opts failed")
    dapl = h5py.h5p.create(h5py.h5p.DATASET_ACCESS)
    dapl.set_chunk_cache(slots, cache, 1.0)
    space = h5py.h5s.create_simple((0, n_ids, nc), (h5py.h5s.UNLIMITED, n_ids, nc))
    did = h5py.h5d.create(g.id, b"DATA", h5py.h5t.IEEE_F64LE, space, dcpl=dcpl, dapl=dapl)
    return h5py.Dataset(did), (ct, cids, nc)


def run(path, a, layout):
    results = [(f"R{i}", a.comp) for i in range(a.results)]
    swmr = layout != "base"
    fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    if swmr:
        fapl.set_libver_bounds(h5py.h5f.LIBVER_LATEST, h5py.h5f.LIBVER_LATEST)
    fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
    f = h5py.File(h5py.h5f.create(os.fsencode(path), h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl))
    on = f.create_group("MODEL_STAGE[1]/RESULTS/ON_NODES")
    sinks = []
    chunk = None
    for ri, (name, nc) in enumerate(results):
        g = on.create_group(name)
        d, chunk = make_data(g, a.nodes, nc, layout, a.compress)
        dt = g.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<f8")
        ds = g.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<i4")
        sinks.append((ri, nc, d, dt, ds))
    f.flush()
    if swmr:
        f.swmr_mode = True
    slabs = {(ri, k): values(ri, k, a.nodes, nc) for ri, nc, *_ in sinks for k in range(1, a.steps + 1)} \
        if a.nodes * a.steps <= 4_000_000 else None
    t0 = time.perf_counter()
    for k in range(1, a.steps + 1):
        t = k - 1
        for ri, nc, d, dt, ds in sinks:
            d.resize(t + 1, axis=0)
            d[t] = slabs[(ri, k)] if slabs is not None else values(ri, k, a.nodes, nc)
            dt.resize(t + 1, axis=0)
            dt[t] = k * 1e-4
            ds.resize(t + 1, axis=0)
            ds[t] = k
        if k == 1 or a.flush_every <= 1 or (k - 1) % a.flush_every == 0:
            f.flush()
    f.flush()
    t_write = time.perf_counter() - t0
    f.close()
    size = os.path.getsize(path)
    with h5py.File(path, "r") as f2:
        d = f2["MODEL_STAGE[1]/RESULTS/ON_NODES/R0/DATA"]
        t1 = time.perf_counter()
        h = d[:, a.nodes // 2, :]
        t_read = time.perf_counter() - t1
        assert h.shape == (a.steps, a.comp)
        assert np.array_equal(h[-1], values(0, a.steps, a.nodes, a.comp)[a.nodes // 2])
    return size, t_write, t_read, chunk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", type=int, default=200)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--flush-every", type=int, default=1, help="flush every K steps (1 = -flush 0)")
    ap.add_argument("--compress", type=int, default=1)
    ap.add_argument("--results", type=int, default=2, help="number of results (one DATA each)")
    ap.add_argument("--comp", type=int, default=3, help="components per result")
    ap.add_argument("--layouts", default="base,swmr,raw,step,edge")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="swmr_layout_")
    raw_bytes = a.steps * a.nodes * a.results * a.comp * 8
    print(f"h5py {h5py.__version__} / HDF5 {h5py.version.hdf5_version}; nodes={a.nodes} "
          f"steps={a.steps} flush_every={a.flush_every} compress={a.compress}; "
          f"raw DATA payload {raw_bytes / 1e6:.2f} MB")
    print(f"{'layout':7s} {'chunk':>16s} {'size MB':>9s} {'x raw':>6s} {'write s':>8s} {'read 1 id ms':>12s}")
    for lay in a.layouts.split(","):
        size, tw, tr, chunk = run(os.path.join(tmp, f"{lay}.ladruno"), a, lay)
        print(f"{lay:7s} {str(chunk):>16s} {size / 1e6:9.2f} {size / raw_bytes:6.2f} {tw:8.2f} {tr * 1e3:12.1f}")


if __name__ == "__main__":
    main()
