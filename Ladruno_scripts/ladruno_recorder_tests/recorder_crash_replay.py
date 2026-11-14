"""Crash-consistency replay model of the LadrunoRecorder write path (#1082, WP-233).

What it shows: whether a ``.ladruno`` written by ``recorder ladruno`` survives a kill at EVERY
point the process can die, and what a killed run keeps. It is the recorder's twin of
``Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py`` (the MPM particle writer) and reuses
its logging file object. It does NOT run the C++ recorder: it mirrors its HDF5 call pattern in
h5py, through the h5py ``fileobj`` driver. Every byte HDF5 hands the driver is logged in order;
the file is then rebuilt as it stands after each write (the state a TerminateProcess leaves)
and read the way a reader does.

Two writers:

``--writer swmr`` (default; the recorder since WP-233):

- ``initialize``: ``H5Fcreate`` with the latest format (``H5Pset_libver_bounds(LATEST,
  LATEST)``), link creation order tracked + indexed; ``INFO`` + ``INFO/PROVENANCE``;
- every LAYOUT STEP (a stage's first recorded step; with ``--kind-rewrite`` the stage's second
  step, which rewrites the provisional ``KIND``; with ``--modes-every M`` every M-th step, a
  modes snapshot ``<MODES>/DATA/STEP_<k>/MODE_<j>``): if SWMR is armed, flush, close (which
  clears the SWMR flags) and reopen read-write; create the objects (``writeModel`` abridged
  to the node set, every result's group + ``ID`` + empty ``DATA``/``TIME``/``STEP``, the
  envelope group with the first sample); hold the step's rows; flush; ``H5Fstart_swmr_write``;
  append the held rows; flush. No streaming row is written outside SWMR;
- steady steps: extend + write one slab of ``DATA``, then ``TIME``, ``STEP`` per result (the
  envelope accumulators rewritten in place with ``--envelope``); ``H5Fflush`` when the cadence
  is due (``--flush-every K`` recorded steps; 0 = every step, ``-flush 0``);
- ``DATA`` layout = ``Ladruno_Hdf5.h::planTimeSeries``: ``-compress N > 0`` (the default 1)
  one step per chunk, ids tiled to 1 MiB, shuffle + deflate N; ``-compress 0`` unfiltered,
  ~64 KiB chunks stacking up to 1024 steps. ``TIME``/``STEP`` chunk 1024, unfiltered.

``--writer wp164`` (the recorder before WP-233, kept for the evidence on #1082): default file
access list (v1 B-tree chunk indexes), no SWMR, one stage, the WP-164 ~1 MiB chunk plan with
deflated partial chunks.

The reader (``--writer swmr``) applies schema §7b.2: the reader rule (plain open, then
SWMR-read, then a repaired COPY, ``ladruno_format.open_ladruno`` / ``repair_copy``) and the
crash rule (every stage cut to ``ladruno_format.committed_steps``). The fileobj driver has no
lock method, so the consistency flags it writes are not those of the real (sec2) file; the
checker sets them to what sec2 writes in each phase, as the MPM model does: 0x01 from a create
or read-write reopen until SWMR is armed, 0x05 while armed, 0x00 after a clean close.

Outcomes per kill point:

- ``OK``: every stage reads; every result's rows ``0 .. T'-1`` equal what was written; every
  stage that ended keeps all its steps; the open stage keeps at least the steps of its last
  COMPLETED flush (and of the close that ended a layout window). ``OK+REPAIRED``: the same,
  through reader rule step 3;
- ``TORN_LAYOUT``: step 3 was needed and the LAST stage does not open (a stage being laid out
  holds no committed step, ``lost=0``; a mid-stage layout step can lose the open stage's
  steps, reported as ``lost``);
- ``LOST``: readable and right, but a stage holds fewer steps than were committed;
- ``MISMATCH``: readable, but a value differs from what was written (silent corruption);
- ``READ_FAIL`` / ``OPEN_FAIL``: the rule does not get a readable file.

Run from anywhere (no OpenSees binary needed):

    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 300 --flush-every 0
    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 300 --flush-every 0 \\
        --stage-every 50 --kind-rewrite --modes-every 40
    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py --steps 300 --writer wp164
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ladruno_format as lf  # noqa: E402
from hdf5_crash_replay import LogFile  # noqa: E402

STAGE = "MODEL_STAGE[1]"
NMODES = 2

# H5Pset_link_phase_change / H5Pset_est_link_info (h5py does not wrap them)
import ctypes  # noqa: E402
_H5 = ctypes.CDLL(os.path.join(os.path.dirname(h5py.__file__), "hdf5.dll")) \
    if os.name == "nt" else ctypes.CDLL(h5py.h5.__file__)
_H5.H5Pset_link_phase_change.argtypes = (ctypes.c_int64, ctypes.c_uint, ctypes.c_uint)
_H5.H5Pset_est_link_info.argtypes = (ctypes.c_int64, ctypes.c_uint, ctypes.c_uint)


def presize_links(plist, n_links, name_len):
    """WP-233 (Ladruno_Hdf5.h::plist::presizeLinks): keep the group's links COMPACT for
    good (max_compact 65535) and size its object header for ``n_links`` up front, so adding
    a link is ONE in-place write of that header. Otherwise HDF5 goes dense after 8 links
    (fractal heap + v2 B-trees rewritten in place) or adds a continuation chunk, and a
    kill inside the layout flush can leave the group unreadable (measured: the root)."""
    if _H5.H5Pset_link_phase_change(plist.id, 65535, 65534) < 0 or \
            _H5.H5Pset_est_link_info(plist.id, n_links, name_len) < 0:
        raise RuntimeError("cannot presize the group's links")


# ---------------------------------------------------------------- layout (mirrors C++)
def plan_time_series(n_ids, n_comp, elem_bytes, scale=1.0):
    """The WP-164 plan (the pre-WP-233 writer), byte budgets scaled by ``scale``."""
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


def plan_swmr(n_ids, n_comp, elem_bytes, compressed, scale=1.0):
    """``Ladruno_Hdf5.h::planTimeSeries`` since WP-233, byte budgets scaled by ``scale``."""
    row_bytes = n_comp * elem_bytes
    slab = n_ids * row_bytes
    if compressed:
        ct = 1
        chunk_ids = max(1, min(n_ids, int((1 << 20) * scale) // row_bytes))
    else:
        target = max(1, int((64 << 10) * scale))
        row_budget = max(1, int((2 << 20) * scale))
        if slab <= target:
            chunk_ids = n_ids
            ct = max(1, min(1024, target // slab if slab else 1))
        else:
            ct = max(1, min(16, row_budget // slab))
            chunk_ids = max(1, min(n_ids, target // (ct * row_bytes)))
    tiles = (n_ids + chunk_ids - 1) // chunk_ids
    chunk_bytes = ct * chunk_ids * row_bytes
    cache = max(1 << 20, min(256 << 20, tiles * chunk_bytes + chunk_bytes))
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


def mode_values(j, k, n_ids):
    i = np.arange(n_ids, dtype=np.float64)[:, None]
    c = np.arange(3, dtype=np.float64)[None, :]
    return np.cos(0.05 * k + 0.3 * j + 0.07 * i + c)


def step_time(k):
    return k * 1e-4


class NoLog:
    """Stand-in for LogFile when writing a real file (the kill stress)."""
    buf = b""
    log = ()

    def mark(self, label):
        pass


def _mark(fobj, label):
    fobj.mark(label)


# ---------------------------------------------------------------- the WP-233 writer
class SwmrWriter:
    """Mirror of LadrunoRecorder's HDF5 calls since WP-233 (see the module docstring)."""

    def __init__(self, fobj, a, results, real_path=None):
        self.fobj, self.a, self.results, self.real_path = fobj, a, results, real_path
        self.name = os.fsencode(real_path) if real_path else b"replay.ladruno"
        self.f = None
        self.armed = False
        self.stage = 0
        self.stage_name = None
        self.stage_step = 0            # recorded steps in the open stage
        self.kind_pending = False
        self.sinks = []                # (ri, nc, path)
        self.t = 0
        self.env = None                # envelope accumulators {ri: (min, max, absmax, arg)}

    def _fapl(self):
        fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
        if self.real_path is None:
            fapl.set_fileobj_driver(h5py.h5fd.fileobj_driver, self.fobj)
        fapl.set_libver_bounds(h5py.h5f.LIBVER_LATEST, h5py.h5f.LIBVER_LATEST)
        return fapl

    def create(self):
        fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
        fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
        if self.a.presize:
            presize_links(fcpl, self.a.presize_root, 20)    # INFO + MODEL_STAGE[...]
        _mark(self.fobj, "phase:layout")
        self.f = h5py.File(h5py.h5f.create(self.name, h5py.h5f.ACC_TRUNC, fcpl=fcpl,
                                           fapl=self._fapl()))
        info = self.f.create_group("INFO", track_order=True)
        for k_, v in (("GENERATOR", "Ladruno"), ("FORMAT_VERSION", 1),
                      ("FORMAT_MINOR_VERSION", 1), ("SOLVER_NAME", "OpenSees"),
                      ("SOLVER_VERSION", np.array([3, 5, 1], "<i4")), ("SPATIAL_DIM", 3),
                      ("STORED_PRECISION", self.a.precision), ("PARTITION_ID", 0),
                      ("PARTITIONED", 0), ("NUM_PARTITIONS", 1),
                      ("RUN_ID", "00000000-0000-4000-8000-000000000000"),
                      ("RUN_ID_SCOPE", "process")):
            info.attrs[k_] = v
        info.create_group("PROVENANCE", track_order=True)

    # -- lifecycle
    def open_layout(self, k):
        if not self.armed:
            return
        _mark(self.fobj, f"s{k}:lflush")
        self.f.flush()
        _mark(self.fobj, f"s{k}:lflushed")
        _mark(self.fobj, "phase:closing")
        self.f.close()                         # clears the SWMR flags (a clean close)
        _mark(self.fobj, "phase:closed")
        _mark(self.fobj, "phase:layout")       # sec2 writes WRITE_ACCESS at the RW open
        self.f = h5py.File(h5py.h5f.open(self.name, h5py.h5f.ACC_RDWR, self._fapl()))
        self.armed = False

    def close_layout(self, k, pending):
        self.f.flush()
        self.f.swmr_mode = True                # H5Fstart_swmr_write
        _mark(self.fobj, "phase:armed")
        self.armed = True
        self._append(pending)
        _mark(self.fobj, f"s{k}:flush")
        self.f.flush()
        _mark(self.fobj, f"s{k}:flushed")

    # -- objects
    def new_stage(self, k):
        self.stage += 1
        self.stage_name = f"MODEL_STAGE[{self.stage}]"
        _mark(self.fobj, f"stage:{self.stage}:{k}")
        n_ids = self.a.nodes
        st = self.f.create_group(self.stage_name, track_order=True)
        st.attrs["STEP"] = k
        st.attrs["TIME"] = step_time(k)
        self.kind_pending = bool(self.a.kind_rewrite)
        if self.a.presize:
            # fixed-length 9 (the longest KIND), NULLTERM: rewritten in place later
            st.attrs.create("KIND", np.array(b"static" if self.kind_pending else b"transient",
                                             dtype="S9"))
        else:
            st.attrs["KIND"] = "static" if self.kind_pending else "transient"
        nodes = st.create_group("MODEL", track_order=True).create_group("NODES", track_order=True)
        ids = np.arange(1, n_ids + 1, dtype="<i4")
        nodes.create_dataset("ID", data=ids.reshape(-1, 1))
        nodes.create_dataset("COORDINATES", data=np.random.default_rng(1).random((n_ids, 3)))
        res = st.create_group("RESULTS", track_order=True)
        on = res.create_group("ON_NODES", track_order=True)
        disk = "<f4" if self.a.precision == "f32" else "<f8"
        eb = 4 if self.a.precision == "f32" else 8
        self.sinks = []
        self.t = 0
        self.stage_step = 0
        self.env = None
        if self.a.envelope:
            return
        if self.a.modes_every and self.a.presize:
            self._modes_group(on)              # created with the stage, not mid-stage
        for ri, (name, nc) in enumerate(self.results):
            g = on.create_group(name, track_order=True)
            for at in ("DISPLAY_NAME", "COMPONENTS", "DIMENSION", "DESCRIPTION"):
                g.attrs[at] = name
            g.attrs["TYPE"] = 0
            g.attrs["DATA_TYPE"] = 1
            g.attrs["PARTITION_REDUCTION"] = "NONE"
            g.create_dataset("ID", data=ids.reshape(-1, 1))
            ct, cids, cache, slots = plan_swmr(n_ids, nc, eb, self.a.compress > 0, self.a.chunk_scale)
            kw = ({"compression": "gzip", "compression_opts": self.a.compress, "shuffle": True}
                  if self.a.compress > 0 else {})
            g.create_dataset("DATA", shape=(0, n_ids, nc), maxshape=(None, n_ids, nc),
                             chunks=(ct, cids, nc), dtype=disk, rdcc_nbytes=cache,
                             rdcc_nslots=slots, rdcc_w0=1.0, **kw)
            axis = max(1, int(1024 * self.a.chunk_scale)) if self.a.scale_axis else 1024
            g.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(axis,), dtype="<f8")
            g.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(axis,), dtype="<i4")
            self.sinks.append((ri, nc, f"{self.stage_name}/RESULTS/ON_NODES/{name}"))
            if self.a.verbose and self.stage == 1:
                print(f"  {name}: DATA chunk ({ct}, {cids}, {nc}); TIME/STEP chunk {axis}")

    def envelope_accept(self, k):
        n_ids = self.a.nodes
        if self.env is None:
            self.env = {}
            for ri, (name, nc) in enumerate(self.results):
                v = values(ri, k, n_ids, nc)
                self.env[ri] = [v.copy(), v.copy(), np.abs(v), np.full(v.shape, k, "<i4")]
            return True       # first sample: the group must be created (a layout step)
        for ri, (name, nc) in enumerate(self.results):
            v = values(ri, k, n_ids, nc)
            mn, mx, ab, arg = self.env[ri]
            np.minimum(mn, v, out=mn)
            np.maximum(mx, v, out=mx)
            upd = np.abs(v) > ab
            ab[upd] = np.abs(v)[upd]
            arg[upd] = k
        return False

    def envelope_create(self):
        base = self.f.require_group(f"{self.stage_name}/RESULTS/ENVELOPES/ON_NODES")
        for ri, (name, nc) in enumerate(self.results):
            g = base.create_group(name, track_order=True)
            g.attrs["COMPONENTS"] = name
            g.create_dataset("ID", data=np.arange(1, self.a.nodes + 1, dtype="<i4").reshape(-1, 1))
            for key, arr in zip(("MIN", "MAX", "ABSMAX", "ARG_STEP"), self.env[ri]):
                g.create_dataset(key, data=arr)

    def envelope_rewrite(self):
        for ri, (name, nc) in enumerate(self.results):
            g = self.f[f"{self.stage_name}/RESULTS/ENVELOPES/ON_NODES/{name}"]
            for key, arr in zip(("MIN", "MAX", "ABSMAX", "ARG_STEP"), self.env[ri]):
                g[key][...] = arr

    def kind_rewrite(self):
        st = self.f[self.stage_name]
        if self.a.presize:
            st.attrs.modify("KIND", np.array(b"transient", dtype="S9"))   # one in-place write
        else:
            del st.attrs["KIND"]
            st.attrs["KIND"] = "transient"
        self.kind_pending = False

    def _modes_group(self, on):
        g = on.create_group("MODES", track_order=True)
        g.attrs["COMPONENTS"] = "U1,U2,U3"
        g.create_dataset("ID", data=np.arange(1, self.a.nodes + 1, dtype="<i4").reshape(-1, 1))
        gcpl = h5py.h5p.create(h5py.h5p.GROUP_CREATE)
        gcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
        if self.a.presize:
            presize_links(gcpl, self.a.presize_modes, 16)
        h5py.h5g.create(g.id, b"DATA", gcpl=gcpl)

    def modes_snapshot(self, k):
        on = self.f[f"{self.stage_name}/RESULTS/ON_NODES"]
        if "MODES" not in on:
            self._modes_group(on)
        sg = on["MODES/DATA"].create_group(f"STEP_{k}", track_order=True)
        sg.attrs["STEP"] = k
        sg.attrs["TIME"] = step_time(k)
        for j in range(NMODES):
            d = sg.create_dataset(f"MODE_{j}", data=mode_values(j, k, self.a.nodes))
            d.attrs["MODE"] = j
            d.attrs["PERIOD"] = 1.0 / (j + 1)

    def _append(self, k_list):
        n_ids = self.a.nodes
        for k in k_list:
            t = self.t
            for ri, nc, path in self.sinks:
                g = self.f[path]
                d, dt, ds = g["DATA"], g["TIME"], g["STEP"]
                d.resize(t + 1, axis=0)
                d[t] = values(ri, k, n_ids, nc)
                dt.resize(t + 1, axis=0)
                dt[t] = step_time(k)
                ds.resize(t + 1, axis=0)
                ds[t] = k
            self.t += 1

    # -- one recorded step (LadrunoRecorder::record)
    def record(self, k, due):
        a = self.a
        new_stage = (self.stage == 0) or (a.stage_every > 0 and (k - 1) % a.stage_every == 0)
        kind = (not new_stage) and self.kind_pending and self.stage_step == 1
        modes = a.modes_every > 0 and k % a.modes_every == 0
        env_first = False
        if a.envelope and not new_stage:
            env_first = self.envelope_accept(k)
        layout = new_stage or kind or modes or env_first
        _mark(self.fobj, f"s{k}:write")
        if layout:
            if new_stage and self.stage and a.envelope:
                pass                               # the close-out rewrite happened in place
            self.open_layout(k)
            if new_stage:
                self.new_stage(k)
                if a.envelope:
                    self.envelope_accept(k)
                    self.envelope_create()
            if kind:
                self.kind_rewrite()
            if modes:
                self.modes_snapshot(k)
            if env_first and not new_stage:
                self.envelope_create()
            self.stage_step += 1
            self.close_layout(k, [k] if not a.envelope else [])
            return
        self._append([k])
        if a.envelope:
            self.envelope_rewrite()                # small envelopes: every recorded step
        self.stage_step += 1
        if due:
            _mark(self.fobj, f"s{k}:flush")
            self.f.flush()
            _mark(self.fobj, f"s{k}:flushed")


def write_run_swmr(fobj, a, results, real_path=None, flush_due=None):
    w = SwmrWriter(fobj, a, results, real_path)
    w.create()
    snapshot = None
    for k in range(1, a.steps + 1):
        if k == a.check_from:
            snapshot = (bytes(fobj.buf), len(fobj.log))
        due = flush_due(k) if flush_due is not None else (
            k == 1 or a.flush_every <= 0 or (k - 1) % a.flush_every == 0)
        w.record(k, due)
        if a.check_to and k >= a.check_to:
            break
    _mark(fobj, "end")
    return snapshot


# ---------------------------------------------------------------- the WP-164 writer
def write_run_wp164(fobj, a, results, real_path=None, flush_due=None):
    """The recorder before WP-233: default fapl (v1 B-trees), no SWMR, one stage."""
    n_ids = a.nodes
    fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    if real_path is None:
        fapl.set_fileobj_driver(h5py.h5fd.fileobj_driver, fobj)
    fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
    name = os.fsencode(real_path) if real_path else b"replay.ladruno"
    f = h5py.File(h5py.h5f.create(name, h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl))
    info = f.create_group("INFO", track_order=True)
    info.attrs["GENERATOR"] = "Ladruno"
    info.create_group("PROVENANCE", track_order=True)
    st = f.create_group(STAGE, track_order=True)
    st.attrs["KIND"] = "transient"
    ids = np.arange(1, n_ids + 1, dtype="<i4")
    nodes = st.create_group("MODEL", track_order=True).create_group("NODES", track_order=True)
    nodes.create_dataset("ID", data=ids.reshape(-1, 1))
    nodes.create_dataset("COORDINATES", data=np.random.default_rng(1).random((n_ids, 3)))
    on_nodes = st.create_group("RESULTS", track_order=True).create_group("ON_NODES", track_order=True)
    disk = "<f4" if a.precision == "f32" else "<f8"
    eb = 4 if a.precision == "f32" else 8
    axis_chunk = max(1, int(1024 * a.chunk_scale)) if a.scale_axis else 1024
    sinks = []
    for ri, (rname, nc) in enumerate(results):
        g = on_nodes.create_group(rname, track_order=True)
        g.attrs["COMPONENTS"] = rname
        g.create_dataset("ID", data=ids.reshape(-1, 1))
        ct, cids, cache, slots = plan_time_series(n_ids, nc, eb, a.chunk_scale)
        kw = {"compression": "gzip", "compression_opts": a.compress, "shuffle": True} if a.compress > 0 else {}
        d = g.create_dataset("DATA", shape=(0, n_ids, nc), maxshape=(None, n_ids, nc),
                             chunks=(ct, cids, nc), dtype=disk, rdcc_nbytes=cache,
                             rdcc_nslots=slots, rdcc_w0=1.0, **kw)
        dt = g.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(axis_chunk,), dtype="<f8")
        ds = g.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(axis_chunk,), dtype="<i4")
        sinks.append((ri, nc, d, dt, ds))
    _mark(fobj, "stage:1:1")
    snapshot = None
    for k in range(1, a.steps + 1):
        if k == a.check_from:
            snapshot = (bytes(fobj.buf), len(fobj.log))
        _mark(fobj, f"s{k}:write")
        t = k - 1
        for ri, nc, d, dt, ds in sinks:
            d.resize(t + 1, axis=0)
            d[t] = values(ri, k, n_ids, nc)
            dt.resize(t + 1, axis=0)
            dt[t] = step_time(k)
            ds.resize(t + 1, axis=0)
            ds[t] = k
        due = flush_due(k) if flush_due is not None else (
            k == 1 or a.flush_every <= 0 or (k - 1) % a.flush_every == 0)
        if due:
            _mark(fobj, f"s{k}:flush")
            f.flush()
            _mark(fobj, f"s{k}:flushed")
        if a.check_to and k >= a.check_to:
            break
    _mark(fobj, "end")
    return snapshot


def write_run(fobj, a, results, real_path=None, flush_due=None):
    if a.writer == "wp164":
        return write_run_wp164(fobj, a, results, real_path, flush_due)
    return write_run_swmr(fobj, a, results, real_path, flush_due)


# ---------------------------------------------------------------- reader
def _ref(ri, k0, n, n_ids, nc, precision):
    arr = np.stack([values(ri, k, n_ids, nc) for k in range(k0, k0 + n)]) if n else \
        np.zeros((0, n_ids, nc))
    if precision == "f32":
        arr = arr.astype("<f4").astype(np.float64)
    return arr


def _check_stage(st, a, results, k0, n_written, n_must, inflight=None):
    """One stage: (status, n_read, detail). k0 = the stage's first step. ``inflight``: the
    step being recorded, whose modes snapshot may be torn (skipped, as a reader does)."""
    n = lf.committed_steps(st)
    if a.envelope:
        return _check_envelope(st, a, results, k0, n_written)
    if n is None:
        return "READ_FAIL", -1, f"{st.name}: no streaming result"
    if n > n_written:
        return "MISMATCH", n, f"{st.name}: {n} steps > {n_written} written"
    for ri, (name, nc) in enumerate(results):
        g = st[f"RESULTS/ON_NODES/{name}"]
        data = np.asarray(g["DATA"][:n], dtype=np.float64)
        if not np.array_equal(data, _ref(ri, k0, n, a.nodes, nc, a.precision)):
            return "MISMATCH", n, f"{st.name}/{name}/DATA"
        ks = np.arange(k0, k0 + n)
        if not np.array_equal(g["STEP"][:n], ks) or not np.array_equal(g["TIME"][:n], step_time(ks)):
            return "MISMATCH", n, f"{st.name}/{name}/TIME|STEP"
    if a.modes_every:
        on = st["RESULTS/ON_NODES"]
        if "MODES" in on:
            for sname in on["MODES/DATA"]:
                k = int(sname[5:])
                try:
                    snap = [on[f"MODES/DATA/{sname}/MODE_{j}"][...] for j in range(NMODES)]
                except Exception as e:  # noqa: BLE001
                    if k == inflight:
                        continue           # the snapshot being written: a reader skips it
                    return "READ_FAIL", n, f"{st.name}/MODES/{sname}: " + str(e).splitlines()[0][:80]
                for j in range(NMODES):
                    if not np.array_equal(snap[j], mode_values(j, k, a.nodes)):
                        if k == inflight:
                            break
                        return "MISMATCH", n, f"{st.name}/MODES/{sname}/MODE_{j}"
    kind = st.attrs.get("KIND")
    if kind not in ("static", "transient", b"static", b"transient"):
        return "MISMATCH", n, f"{st.name}/KIND={kind!r}"
    if n < n_must:
        return "LOST", n, f"{st.name}: read {n} < committed {n_must}"
    return "OK", n, ""


def _check_envelope(st, a, results, k0, n_written):
    """Every accumulator entry is the true extreme over steps k0..k for some k <= written
    (an in-place rewrite torn between two steps mixes neighbours, both true)."""
    if "RESULTS/ENVELOPES/ON_NODES" not in st:
        return "OK", 0, ""
    for ri, (name, nc) in enumerate(results):
        g = st[f"RESULTS/ENVELOPES/ON_NODES/{name}"]
        mn = g["MIN"][...]
        seq = np.stack([values(ri, k, a.nodes, nc) for k in range(k0, k0 + max(1, n_written))])
        prefix_min = np.minimum.accumulate(seq, axis=0)
        if not np.all(np.any(prefix_min == mn[None], axis=0)):
            return "MISMATCH", -1, f"{st.name}/ENVELOPES/{name}/MIN"
    return "OK", n_written, ""


def check_swmr(path, a, results, phase, stages, committed, inflight=None, flags=None):
    """stages: [(name, k0, n_written)]; committed: {name: steps a completed flush holds}.
    ``flags``: the superblock flags to set (default: what sec2 writes in ``phase``)."""
    if flags is None:
        flags = {"layout": 0x01, "armed": 0x05, "closing": 0x05, "closed": 0x00}[phase]
    try:
        lf.set_superblock(path, flags=flags)
    except ValueError:
        pass                                    # no superblock yet
    repaired = False
    try:
        f = lf.open_ladruno(path)
    except lf.StageLayoutKilledError:
        lf.repair_copy(path, path + ".rep")
        repaired = True
        try:
            f = h5py.File(path + ".rep", "r")
        except Exception as e:  # noqa: BLE001
            return "OPEN_FAIL", -1, "repaired copy: " + str(e).splitlines()[0][:120], 0
    except Exception as e:  # noqa: BLE001
        return "OPEN_FAIL", -1, str(e).splitlines()[0][:160], 0
    suffix = "+REPAIRED" if repaired else ""
    with f:
        n_last = 0
        for i, (name, k0, n_written) in enumerate(stages):
            last = i == len(stages) - 1
            try:
                st = f[name]
                status, n, detail = _check_stage(st, a, results, k0, n_written,
                                                 committed.get(name, 0), inflight)
            except Exception as e:  # noqa: BLE001
                status, n, detail = "READ_FAIL", -1, f"{name}: " + str(e).splitlines()[0][:120]
            if status in ("READ_FAIL",) and last and repaired:
                return "TORN_LAYOUT", -1, f"{name}: {detail}", committed.get(name, 0)
            if status == "READ_FAIL" and last and phase == "layout":
                return "TORN_LAYOUT", -1, f"{name}: {detail}", committed.get(name, 0)
            if status != "OK":
                return status + suffix, n, detail, 0
            n_last = n
    return "OK" + suffix, n_last, "", 0


def check(path, results, exp, swmr, precision="f64"):
    """The pre-WP-233 checker (one stage, plain open then SWMR-read); used by the
    h5py-mirror kill stress for ``--writer wp164``. Returns (status, n_min, detail)."""
    try:
        f = lf.open_ladruno(path)
    except Exception as e:  # noqa: BLE001
        return "OPEN_FAIL", -1, str(e).splitlines()[0][:160]
    with f:
        ns, ragged = [], False
        try:
            base = f[f"{STAGE}/RESULTS/ON_NODES"]
            for ri, (name, nc) in enumerate(results):
                g = base[name]
                try:
                    data = g["DATA"][...]
                except Exception as e:  # noqa: BLE001
                    return "READ_FAIL", -1, f"{name}/DATA: {str(e).splitlines()[0][:90]}"
                tim, stp = g["TIME"][...], g["STEP"][...]
                n = min(len(data), len(tim), len(stp))
                ragged |= not (len(data) == len(tim) == len(stp))
                ref = exp[name] if exp is not None else _ref(ri, 1, len(data), data.shape[1], nc, precision)
                if n > len(ref):
                    return "MISMATCH", n, f"{name}: {n} steps > written"
                if not np.array_equal(np.asarray(data[:n], dtype=np.float64), ref[:n]):
                    return "MISMATCH", n, f"{name}/DATA"
                if not np.array_equal(stp[:n], np.arange(1, n + 1)):
                    return "MISMATCH", n, f"{name}/STEP"
                ns.append(n)
        except Exception as e:  # noqa: BLE001
            return "READ_FAIL", -1, str(e).splitlines()[0][:160]
    return ("RAGGED" if ragged else "OK"), (min(ns) if ns else 0), ""


# ---------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--writer", choices=("swmr", "wp164"), default="swmr")
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--nodes", type=int, default=200)
    ap.add_argument("--results", default="DISPLACEMENT:3,VELOCITY:3")
    ap.add_argument("--flush-every", type=int, default=0,
                    help="flush every K recorded steps (0 = every step, -flush 0); first step always")
    ap.add_argument("--compress", type=int, default=1, help="-compress level (default 1; 0 = no filter)")
    ap.add_argument("--precision", choices=("f64", "f32"), default="f64")
    ap.add_argument("--chunk-scale", type=float, default=1.0,
                    help="scale of the chunk byte budgets (1 = the real layout)")
    ap.add_argument("--scale-axis", action="store_true",
                    help="also scale the 1024-element TIME/STEP chunk by --chunk-scale")
    ap.add_argument("--stage-every", type=int, default=0, help="a new MODEL_STAGE every K steps")
    ap.add_argument("--kind-rewrite", action="store_true",
                    help="each stage opens with a provisional KIND, rewritten at its 2nd step")
    ap.add_argument("--modes-every", type=int, default=0, help="a modes snapshot every M steps")
    ap.add_argument("--envelope", action="store_true", help="-envelope mode (no time series)")
    ap.add_argument("--no-presize", dest="presize", action="store_false",
                    help="the first WP-233 cut: no presized link headers, KIND delete + create, "
                         "MODES created at the first snapshot")
    ap.add_argument("--presize-root", type=int, default=1024, help="root links presized for")
    ap.add_argument("--presize-modes", type=int, default=1024, help="modes snapshots presized for")
    ap.add_argument("--check-from", type=int, default=1, help="first step whose kill points are checked")
    ap.add_argument("--check-to", type=int, default=0, help="stop the run after this step (0 = --steps)")
    ap.add_argument("--show", type=int, default=15)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    if a.writer == "wp164":
        for opt in ("stage_every", "kind_rewrite", "modes_every", "envelope"):
            if getattr(a, opt):
                ap.error(f"--{opt.replace('_', '-')} needs --writer swmr")
    results = parse_results(a.results)

    t0 = time.time()
    fobj = LogFile()
    print(f"h5py {h5py.__version__} / HDF5 {h5py.version.hdf5_version}; writer={a.writer} "
          f"nodes={a.nodes} results={a.results} steps={a.steps} flush_every={a.flush_every} "
          f"compress={a.compress} precision={a.precision} chunk_scale={a.chunk_scale} "
          f"stage_every={a.stage_every} kind_rewrite={a.kind_rewrite} "
          f"modes_every={a.modes_every} envelope={a.envelope}")
    snapshot = write_run(fobj, a, results)
    log = fobj.log
    if snapshot is None:
        snapshot = (b"", 0)
    image = bytearray(snapshot[0])
    start = snapshot[1]
    nw = sum(1 for e in log[start:] if e[0] == "w")
    print(f"  wrote in {time.time() - t0:.1f} s; ops in window={len(log) - start} writes={nw}; "
          f"final image {len(fobj.buf)} B")

    # replay the marks before the window to know the state at its start
    phase = "layout"
    stages = []              # [name, k0, n_written]
    committed = {}           # name -> steps covered by the last completed flush / close
    cur_step = 0
    region = "start"
    for e in log[:start]:
        if e[0] == "m":
            phase, region, cur_step = _apply_mark(e[1], phase, region, cur_step, stages, committed)

    tmp = tempfile.mkdtemp(prefix="rec_replay_")
    path = os.path.join(tmp, "k.ladruno")
    stats = collections.Counter()
    by_region = collections.Counter()
    fail_steps = collections.defaultdict(set)
    torn = collections.Counter()
    torn_detail = []
    fails = []
    for i in range(start, len(log)):
        e = log[i]
        if e[0] == "m":
            phase, region, cur_step = _apply_mark(e[1], phase, region, cur_step, stages, committed)
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
        if cur_step < a.check_from:
            continue
        with open(path, "wb") as fh:
            fh.write(image)
        if a.writer == "wp164":
            status, n, detail = check(path, results, None, False, a.precision)
            lost = 0
            if status in ("OK", "RAGGED") and n < committed.get(STAGE, 0):
                status, detail = "LOST", f"read {n} < flushed {committed.get(STAGE, 0)}"
        else:
            status, n, detail, lost = check_swmr(
                path, a, results, phase, [tuple(s) for s in stages], committed, cur_step)
            if phase == "layout":
                # sec2 writes WRITE_ACCESS into the superblock only when it next flushes
                # it, so a layout step's flags may still read 0x00 (a plain open): check
                # that image too, and keep the worse outcome
                with open(path, "wb") as fh:
                    fh.write(image)
                s2 = check_swmr(path, a, results, phase, [tuple(s) for s in stages],
                                committed, cur_step, flags=0x00)
                if not s2[0].startswith("OK") and (status.startswith("OK") or s2[0] != "TORN_LAYOUT"):
                    status, n, detail, lost = s2
        if status == "OPEN_FAIL" and not committed:
            status = "NO_FILE"       # before the first flush ends: nothing committed yet
        stats[status] += 1
        by_region[(region, phase, status)] += 1
        if status == "TORN_LAYOUT":
            torn[lost] += 1
            if lost:
                torn_detail.append((i, cur_step, region, phase, status, n, detail))
        if not status.startswith("OK") and status not in ("NO_FILE", "TORN_LAYOUT"):
            fail_steps[status].add(cur_step)
            fails.append((i, cur_step, region, phase, status, n, detail))
    print(f"  replayed in {time.time() - t0:.1f} s")
    print("kill points checked:", sum(stats.values()), dict(stats))
    if torn:
        print("TORN_LAYOUT (last stage does not open after the step-3 repair): committed steps "
              "lost -> kill points:", dict(sorted(torn.items())))
    print("by region (the mark before the kill point) / sec2 phase / status:")
    for (reg, ph, status), c in sorted(by_region.items()):
        print(f"  {reg:9s} {ph:8s} {status:22s} {c}")
    for status, ks in sorted(fail_steps.items()):
        ks = sorted(ks)
        print(f"  steps with a {status} kill point ({len(ks)}): {ks[:40]}{' ...' if len(ks) > 40 else ''}")
    for i, k, reg, ph, status, n, detail in torn_detail[: a.show]:
        print(f"  TORN op#{i} step={k} region={reg} phase={ph} read={n} {detail}")
    for i, k, reg, ph, status, n, detail in fails[: a.show]:
        print(f"  FAIL op#{i} step={k} region={reg} phase={ph} {status} read={n} {detail}")
    bad = sum(c for s, c in stats.items() if not s.startswith("OK") and s not in ("NO_FILE", "TORN_LAYOUT"))
    print(f"RESULT: {bad} bad kill points of {sum(stats.values())}")
    return 1 if bad else 0


def _apply_mark(label, phase, region, cur_step, stages, committed):
    """Advance the replay's view of the writer at a mark."""
    if label.startswith("phase:"):
        return label[6:], region, cur_step
    if label.startswith("stage:"):
        _, sid, k = label.split(":")
        stages.append([f"MODEL_STAGE[{sid}]", int(k), 1])   # its first step is in flight
        return phase, region, cur_step
    if label == "end":
        return phase, "end", cur_step
    if label.startswith("s"):
        k = int(label[1:].split(":")[0])
        kind = label.split(":")[1]
        if kind == "write":
            # step k is being recorded: the open stage may show its rows (an upper bound)
            if stages:
                stages[-1][2] = k - stages[-1][1] + 1
            return phase, "write", k
        if kind in ("flush", "lflush"):
            return phase, kind, k
        if kind == "flushed":
            # a completed flush after step k's rows: the open stage holds k0 .. k
            if stages:
                committed[stages[-1][0]] = k - stages[-1][1] + 1
            return phase, kind, k
        if kind == "lflushed":
            # the flush that opens a layout step (before step k's rows): k0 .. k-1
            if stages:
                committed[stages[-1][0]] = k - stages[-1][1]
            return phase, kind, k
    return phase, region, cur_step


if __name__ == "__main__":
    sys.exit(main())
