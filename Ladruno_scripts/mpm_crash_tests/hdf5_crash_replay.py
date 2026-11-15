"""Crash-consistency replay model of the MPM particle writer's HDF5 call pattern.

What it shows: whether a ``.ladruno`` particle file survives a kill at EVERY point the
process can die, not only at the few points a real kill happens to land. It mirrors the
HDF5 call pattern of ``SRC/domain/pattern/ladrunoMPM/LadrunoMPMIO.cpp`` (createFile /
openStage / writeOutput at the default ``-flush 0``: the rows, the write barrier flush,
the OFFSETS entry, the flush after the output) in h5py, through the h5py ``fileobj``
driver. It does NOT run the C++ writer itself: it models its sequence of HDF5 calls and
layout (12 particles, POSITION / VELOCITY / STRESS / VOLUME, gzip one-chunk-per-output
ID / DATA, uncompressed STEP / TIME / OFFSETS). Every byte HDF5 hands the driver is
logged in order; the file is then rebuilt as it stands after each write (the state a
TerminateProcess leaves) and read the way a crash reader does: open by the reader rule,
232a section 2.5 ``committed_outputs()`` (``tests/mpm_contracts.py``), then the committed
rows, compared EXACTLY with what was written. Kill points are checked from output 0's
first write on: before it the file and stage are still being laid out, and nothing is
committed (a reader skips such a stage: no OFFSETS, or crash-rule step 0).

Measured (2026-11-14, h5py 3.16 / HDF5 2.0.0):

- without SWMR (the pre-cbe1b7149 writer): 30 failing kill points of ~3,300 in 300
  outputs, all inside the BARRIER flush at outputs 64 and 121, where a dataset's v1 B-tree
  chunk index splits (an ID node is 2096 B, a DATA node 2616 B). H5Fflush writes the
  in-place parent (user ring) before the superblock's new end-of-allocation (superblock
  ring); a kill between them leaves a parent pointing past the stored EOA ("addr
  overflow"). Identical with metadata-cache eviction on or off (``--evict-off``); chunks
  of R outputs (``--rpc R``) only move the same window R times later;
- with SWMR (``--swmr``, the writer since cbe1b7149): 0 failing kill points of 17,439 in
  600 outputs; every completed output recovers bit-exact.

The ``mpm reset`` window (``--resets N``, WP-232h / #1085; measured 2026-11-15, 240 outputs, a
reset every 10, 24 stages, every kill point from output 0 on): each reset closes the file,
reopens it read-write, lays out a stage WITHOUT SWMR and re-arms. Before the fix
(``--no-presize``: the root's header unsized) 634 of 9,376 kill points lost the ROOT, and with
it every earlier stage ("addr overflow" on ``/``); after it (``plist::presizeLinks`` on the
fcpl) 0 of 9,249, every earlier stage exact; only the stage being laid out may be torn.

HDF5 2.0 + fileobj quirk: under SWMR the fileobj driver writes only the 0x04 consistency
flag into the superblock, where the sec2 driver (the real writer) writes 0x05, and an
SWMR-read open requires 0x05. So when a rebuilt image is refused as "already open for
write", the model sets the flags to 0x05 and re-checksums the version-3 superblock with
HDF5's lookup3 hash before the SWMR-read open (``mc.set_superblock``).

Run from anywhere (no OpenSees binary needed; ~minutes per few hundred outputs):

    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 300
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 300 --evict-off
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 600 --swmr
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 240 --resets 10
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 240 --resets 10 --no-presize
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "Ladruno_scripts" / "ladruno_recorder_tests"))
import mpm_contracts as mc

M = 12
FIELDS = [("POSITION", 2), ("VELOCITY", 2), ("STRESS", 3), ("VOLUME", 1)]
IDS = np.array([3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610], dtype=np.int64)


class LogFile:
    """In-memory file that logs every write / truncate in order."""

    def __init__(self):
        self.buf = bytearray()
        self.pos = 0
        self.log = []          # ('w', off, bytes) | ('t', size) | ('m', label)

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else len(self.buf) + off)
        return self.pos

    def tell(self):
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = len(self.buf) - self.pos
        b = bytes(self.buf[self.pos:self.pos + n])
        self.pos += len(b)
        return b

    def readinto(self, mv):
        b = self.read(len(mv))
        mv[:len(b)] = b
        return len(b)

    def write(self, data):
        data = bytes(data)
        end = self.pos + len(data)
        if end > len(self.buf):
            self.buf.extend(b"\0" * (end - len(self.buf)))
        self.buf[self.pos:end] = data
        self.log.append(("w", self.pos, data))
        self.pos = end
        return len(data)

    def truncate(self, size=None):
        size = self.pos if size is None else size
        if size < len(self.buf):
            del self.buf[size:]
        else:
            self.buf.extend(b"\0" * (size - len(self.buf)))
        self.log.append(("t", size))
        return size

    def flush(self):
        pass

    def mark(self, label):
        self.log.append(("m", label))

    def readable(self):
        return True

    def writable(self):
        return True

    def seekable(self):
        return True


def rows(k, nc):
    i = np.arange(M)[:, None]
    c = np.arange(nc)[None, :]
    return (k * 1000.0 + i * 10.0 + c + 0.5).astype(np.float64)


def append(ds, data):
    n0 = ds.shape[0]
    ds.resize(n0 + data.shape[0], axis=0)
    ds[n0:] = data


def write_run(fobj, n_out, evict_off, cache_init, rpc=0, swmr=False):
    fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    fapl.set_fileobj_driver(h5py.h5fd.fileobj_driver, fobj)
    if swmr:
        fapl.set_libver_bounds(h5py.h5f.LIBVER_LATEST, h5py.h5f.LIBVER_LATEST)
    cfg = fapl.get_mdc_config()
    if evict_off:
        cfg.evictions_enabled = False
        cfg.incr_mode = 0
        cfg.flash_incr_mode = 0
        cfg.decr_mode = 0
        cfg.set_initial_size = True
        cfg.initial_size = cache_init
        cfg.max_size = max(cfg.max_size, cache_init)
        cfg.min_size = min(cfg.min_size, cache_init)
        fapl.set_mdc_config(cfg)
    fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
    fid = h5py.h5f.create(b"replay.ladruno", h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl)
    f = h5py.File(fid)
    f.attrs["FORMAT_VERSION"] = 1
    f.attrs["FORMAT_MINOR_VERSION"] = 2
    f.attrs["KIND"] = "transient"
    f.flush()
    fobj.mark("file-created")

    cmp = {"compression": "gzip", "compression_opts": 1, "shuffle": True, "track_times": True}
    rows_cmp = cmp if rpc == 0 else {"track_times": True}   # multi-output chunks: no filter
    crows = M if rpc == 0 else rpc * M
    st = f.create_group("MODEL_STAGE[0]", track_order=True)
    st.attrs["PARTICLE_EPOCH"] = 0
    st.attrs["STEP"] = 0
    model = st.create_group("MODEL", track_order=True)
    mp = model.create_group("PARTICLES", track_order=True)
    for name, dt, val in (("ID", "<i8", IDS), ("BODY_ID", "<i4", np.ones(M, "<i4")),
                          ("MATERIAL_TAG", "<i4", np.ones(M, "<i4")),
                          ("MASS", "<f8", np.full(M, 2.5))):
        d = mp.create_dataset(name, shape=(0,), maxshape=(None,), chunks=(M,), dtype=dt, **cmp)
        append(d, np.asarray(val, dtype=dt))
    res = st.create_group("RESULTS", track_order=True)
    op = res.create_group("ON_PARTICLES", track_order=True)
    op.attrs["PARTITION_REDUCTION"] = "NONE"
    raw = {"track_times": True}
    ds_step = op.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<i4", **raw)
    ds_time = op.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<f8", **raw)
    ds_off = op.create_dataset("OFFSETS", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<i8", **raw)
    append(ds_off, np.array([0], "<i8"))
    ds_id = op.create_dataset("ID", shape=(0,), maxshape=(None,), chunks=(crows,), dtype="<i8", **rows_cmp)
    ds_data = {}
    for name, nc in FIELDS:
        g = op.create_group(name, track_order=True)
        for a in ("DISPLAY_NAME", "COMPONENTS", "DIMENSION", "DESCRIPTION"):
            g.attrs[a] = name
        g.attrs["TYPE"] = 0
        g.attrs["DATA_TYPE"] = 1
        g.attrs["PARTITION_REDUCTION"] = "NONE"
        ds_data[name] = g.create_dataset("DATA", shape=(0, nc), maxshape=(None, nc),
                                         chunks=(crows, nc), dtype="<f8", **rows_cmp)
    f.flush()
    if swmr:
        f.swmr_mode = True
    fobj.mark("stage-open")

    S = 0
    for k in range(n_out):
        fobj.mark(f"out{k}:rows")
        for name, nc in FIELDS:
            append(ds_data[name], rows(k, nc))
        append(ds_id, IDS)
        append(ds_step, np.array([k + 1], "<i4"))
        append(ds_time, np.array([(k + 1) * 1e-6], "<f8"))
        fobj.mark(f"out{k}:barrier")
        f.flush()                                   # the write barrier
        fobj.mark(f"out{k}:offsets")
        S += M
        append(ds_off, np.array([S], "<i8"))
        fobj.mark(f"out{k}:commit")
        f.flush()                                   # flush after every output (-flush 0)
        fobj.mark(f"out{k}:done")
    f.close()
    fobj.mark("closed")


def check(path):
    """Returns (status, n_committed, detail)."""
    try:
        f = h5py.File(path, "r")
    except Exception as e:  # noqa: BLE001
        try:
            if "already open for write" in str(e):
                mc.set_superblock(path, flags=0x05)   # sec2's SWMR flags; fileobj writes 0x04
            f = h5py.File(path, "r", swmr=True)
        except Exception as e2:  # noqa: BLE001
            return "OPEN_FAIL", -1, (str(e).splitlines()[0][:120] + " | swmr: " + str(e2).splitlines()[0][:120])
    with f:
        try:
            if "MODEL_STAGE[0]/RESULTS/ON_PARTICLES/OFFSETS" not in f:
                return "NO_STAGE", 0, ""
            stg = f["MODEL_STAGE[0]"]
            n = mc.committed_outputs(stg)
        except Exception as e:  # noqa: BLE001
            return "RULE_FAIL", -1, repr(e)[:200]
        op = stg["RESULTS/ON_PARTICLES"]
        try:
            offs = op["OFFSETS"][: n + 1]
            steps = op["STEP"][:n]
            times = op["TIME"][:n]
            end = int(offs[n]) if n else 0
            ids = op["ID"][:end]
            data = {name: op[name]["DATA"][:end] for name, _ in FIELDS}
        except Exception as e:  # noqa: BLE001
            return "READ_FAIL", n, str(e).splitlines()[0][:200]
        for k in range(n):
            a, b = int(offs[k]), int(offs[k + 1])
            if b - a != M or not np.array_equal(ids[a:b], IDS):
                return "MISMATCH", n, f"ids of output {k}"
            if steps[k] != k + 1 or times[k] != (k + 1) * 1e-6:
                return "MISMATCH", n, f"step/time of output {k}"
            for name, nc in FIELDS:
                if not np.array_equal(data[name][a:b], rows(k, nc)):
                    return "MISMATCH", n, f"{name} of output {k}"
        return "OK", n, ""


# ---------------------------------------------------------------- the mpm reset window
# WP-232h (#1085). The writer with ``mpm reset``: a new MODEL_STAGE every ``resets``
# outputs. This mirrors LadrunoMPMIO.cpp's createFile / openStage / writeOutput since
# cbe1b7149 (SWMR), including the INFO group the root holds, so the root's links are
# INFO + one MODEL_STAGE per stage:
#
# - createFile: H5Fcreate (latest format; link creation order tracked + indexed; WITH
#   ``plist::presizeLinks(fcpl, 1024, 20)`` unless ``--no-presize``), INFO + its
#   attributes + INFO/PROVENANCE;
# - openStage: at a reset, flush, close (clears the SWMR flags), reopen read-write (and,
#   since WP-232h, flush at once); then the stage's layout (attributes, MODEL/NODES,
#   ELEMENTS, PARTICLES + GRID + BODIES + the registry, RESULTS/ON_* and ON_PARTICLES
#   with its fields, OFFSETS[0]); flush; ``H5Fstart_swmr_write``;
# - the stage's output 0 at once, then every output: rows, the barrier flush, the
#   OFFSETS entry, the commit flush (``-flush 0``).
#
# The fileobj driver has no lock method, so the consistency flags it writes are not
# sec2's; the checker sets what sec2 writes in each phase (as recorder_crash_replay.py
# does): 0x01 from a create or read-write reopen until SWMR is armed, 0x05 while armed or
# closing, 0x00 once a reset's close completed. In a layout the image is ALSO checked
# with flags 0x00, because sec2 writes WRITE_ACCESS into the superblock lazily (only
# when it next flushes it), so a plain open can get in; the worse outcome counts.

ROOT_LINKS = 1024


def _stage_rows(k):
    """Output k's position etc. (global output index: every stage's rows differ)."""
    return {name: rows(k, nc) for name, nc in FIELDS}


def write_run_resets(fobj, n_out, resets, presize):
    def fapl_():
        fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
        fapl.set_fileobj_driver(h5py.h5fd.fileobj_driver, fobj)
        fapl.set_libver_bounds(h5py.h5f.LIBVER_LATEST, h5py.h5f.LIBVER_LATEST)
        return fapl

    fcpl = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    fcpl.set_link_creation_order(h5py.h5p.CRT_ORDER_TRACKED | h5py.h5p.CRT_ORDER_INDEXED)
    if presize:
        from recorder_crash_replay import presize_links   # ctypes H5Pset_* (WP-233)
        presize_links(fcpl, ROOT_LINKS, 20)
    fobj.mark("phase:layout")
    f = h5py.File(h5py.h5f.create(b"replay.ladruno", h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl_()))
    info = f.create_group("INFO", track_order=True)
    for k_, v in (("GENERATOR", "Ladruno"), ("FORMAT_VERSION", 1), ("FORMAT_MINOR_VERSION", 2),
                  ("SOLVER_NAME", "OpenSees"), ("SOLVER_VERSION", np.array([3, 5, 1], "<i4")),
                  ("SPATIAL_DIM", 2), ("STORED_PRECISION", "f64"), ("PARTITION_ID", 0),
                  ("PARTITIONED", 0), ("NUM_PARTITIONS", 1),
                  ("RUN_ID", "proc-0123456789abcdef-12345678"), ("RUN_ID_SCOPE", "process")):
        info.attrs[k_] = v
    info.create_group("PROVENANCE", track_order=True)

    cmp = {"compression": "gzip", "compression_opts": 1, "shuffle": True, "track_times": True}
    raw = {"track_times": True}
    stage = -1
    S = 0
    ds = None
    for k in range(n_out):
        if k == 0 or k % resets == 0:
            if stage >= 0:                              # an mpm reset: openStage's reopen
                f.flush()
                fobj.mark("phase:closing")
                f.close()                               # clears the SWMR flags
                fobj.mark("phase:closed")
                fobj.mark("phase:layout")               # sec2 writes WRITE_ACCESS at the RW open
                f = h5py.File(h5py.h5f.open(b"replay.ladruno", h5py.h5f.ACC_RDWR, fapl_()))
                if presize:
                    f.flush()                           # WP-232h: WRITE_ACCESS on disk at once
            stage += 1
            fobj.mark(f"stage:{stage}:{k}")
            st = f.create_group(f"MODEL_STAGE[{stage}]", track_order=True)
            st.attrs["STEP"] = k
            st.attrs["TIME"] = k * 1e-6
            st.attrs["KIND"] = "transient"
            model = st.create_group("MODEL", track_order=True)
            res = st.create_group("RESULTS", track_order=True)
            for g in ("ON_NODES", "ON_ELEMENTS", "ON_DOMAIN", "ON_REGIONS"):
                res.create_group(g, track_order=True)
            st.attrs["EMPTY_PARTITION"] = 1
            nodes = model.create_group("NODES", track_order=True)
            nodes.create_dataset("ID", shape=(0, 1), dtype="<i4")
            nodes.create_dataset("COORDINATES", shape=(0, 2), dtype="<f8")
            model.create_group("ELEMENTS", track_order=True)
            st.attrs["PARTICLE_EPOCH"] = stage
            mp = model.create_group("PARTICLES", track_order=True)
            for a_, v in (("PATTERN_TAG", 1), ("MODE", "explicit"), ("PLANE_STRAIN_THICKNESS", 1.0),
                          ("UNIT_LENGTH", "m"), ("UNIT_FORCE", "kN"), ("UNIT_TIME", "s")):
                mp.attrs[a_] = v
            grid = mp.create_group("GRID", track_order=True)
            grid.attrs["ORIGIN"] = np.zeros(2)
            grid.attrs["CELL_SIZE"] = np.ones(2)
            grid.attrs["NUM_CELLS"] = np.array([4, 4], "<i4")
            grid.attrs["BASIS"] = "gimp"
            bodies = mp.create_group("BODIES", track_order=True)
            bodies.create_dataset("ID", data=np.array([1], "<i4"))
            bodies.create_dataset("NAME", data=np.array([b"soil"]))
            for name, dt, val in (("ID", "<i8", IDS), ("BODY_ID", "<i4", np.ones(M, "<i4")),
                                  ("MATERIAL_TAG", "<i4", np.ones(M, "<i4")),
                                  ("MASS", "<f8", np.full(M, 2.5))):
                d = mp.create_dataset(name, shape=(0,), maxshape=(None,), chunks=(M,), dtype=dt, **cmp)
                append(d, np.asarray(val, dtype=dt))
            op = res.create_group("ON_PARTICLES", track_order=True)
            op.attrs["PARTITION_REDUCTION"] = "NONE"
            ds_step = op.create_dataset("STEP", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<i4", **raw)
            ds_time = op.create_dataset("TIME", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<f8", **raw)
            ds_off = op.create_dataset("OFFSETS", shape=(0,), maxshape=(None,), chunks=(1024,), dtype="<i8", **raw)
            append(ds_off, np.array([0], "<i8"))
            ds_id = op.create_dataset("ID", shape=(0,), maxshape=(None,), chunks=(M,), dtype="<i8", **cmp)
            ds_data = {}
            for name, nc in FIELDS:
                g = op.create_group(name, track_order=True)
                for a_ in ("DISPLAY_NAME", "COMPONENTS", "DIMENSION", "DESCRIPTION"):
                    g.attrs[a_] = name
                g.attrs["TYPE"] = 0
                g.attrs["DATA_TYPE"] = 1
                g.attrs["PARTITION_REDUCTION"] = "NONE"
                ds_data[name] = g.create_dataset("DATA", shape=(0, nc), maxshape=(None, nc),
                                                 chunks=(M, nc), dtype="<f8", **cmp)
            f.flush()
            f.swmr_mode = True                          # H5Fstart_swmr_write
            fobj.mark("phase:armed")
            S = 0
            ds = (ds_step, ds_time, ds_off, ds_id, ds_data)
        ds_step, ds_time, ds_off, ds_id, ds_data = ds
        fobj.mark(f"out{k}:rows")
        for name, nc in FIELDS:
            append(ds_data[name], rows(k, nc))
        append(ds_id, IDS)
        append(ds_step, np.array([k + 1], "<i4"))
        append(ds_time, np.array([(k + 1) * 1e-6], "<f8"))
        fobj.mark(f"out{k}:barrier")
        f.flush()
        fobj.mark(f"out{k}:offsets")
        S += M
        append(ds_off, np.array([S], "<i8"))
        fobj.mark(f"out{k}:commit")
        f.flush()
        fobj.mark(f"out{k}:done")
    fobj.mark("phase:closing")
    f.close()
    fobj.mark("phase:closed")


def _check_one_stage(st, k0, n_must, n_max):
    """(status, n, detail) of one stage against what was written."""
    n = mc.committed_outputs(st)
    op = st["RESULTS/ON_PARTICLES"]
    if n < n_must:
        return "LOST", n, f"{n} < committed {n_must}"
    if n > n_max:
        return "MISMATCH", n, f"{n} > written {n_max}"
    offs = op["OFFSETS"][: n + 1]
    end = int(offs[n]) if n else 0
    ids = op["ID"][:end]
    data = {name: op[name]["DATA"][:end] for name, _ in FIELDS}
    steps, times = op["STEP"][:n], op["TIME"][:n]
    for j in range(n):
        k = k0 + j
        a, b = int(offs[j]), int(offs[j + 1])
        if b - a != M or not np.array_equal(ids[a:b], IDS):
            return "MISMATCH", n, f"ids of output {k}"
        if steps[j] != k + 1 or times[j] != (k + 1) * 1e-6:
            return "MISMATCH", n, f"step/time of output {k}"
        for name, nc in FIELDS:
            if not np.array_equal(data[name][a:b], rows(k, nc)):
                return "MISMATCH", n, f"{name} of output {k}"
    return "OK", n, ""


def check_resets(path, phase, stages, flags=None):
    """stages: [[name, k0, n_written, n_committed]] in order. Returns (status, detail).

    OK / OK+REPAIRED: every stage reads, every earlier stage keeps all its outputs, the
    last keeps at least its committed ones, every value exact. TORN_LAYOUT: only the
    stage being laid out (the last; it holds no committed output) does not read.
    EARLIER_LOST: the root or an EARLIER stage does not read (the #1085 failure).
    OPEN_FAIL / READ_FAIL / LOST / MISMATCH: other bad outcomes. NO_FILE: nothing
    committed yet and nothing opens."""
    if flags is None:
        flags = {"layout": 0x01, "armed": 0x05, "closing": 0x05, "closed": 0x00}[phase]
    try:
        mc.set_superblock(path, flags=flags)
    except (ValueError, OSError):
        pass                                        # no (whole) superblock yet
    repaired = False
    try:
        f = mc.open_particle_file(path)
    except mc.StageLayoutKilledError:
        mc.lf.repair_copy(path, path + ".rep")
        repaired = True
        try:
            f = h5py.File(path + ".rep", "r")
        except Exception as e:  # noqa: BLE001
            return "OPEN_FAIL", "repaired copy: " + str(e).splitlines()[0][:140]
    except Exception as e:  # noqa: BLE001
        if not any(s[3] for s in stages):
            return "NO_FILE", ""
        return "OPEN_FAIL", str(e).splitlines()[0][:160]
    suffix = "+REPAIRED" if repaired else ""
    with f:
        try:
            names = set(f.keys())                   # reads the root's links
        except Exception as e:  # noqa: BLE001
            return "EARLIER_LOST", "root: " + str(e).splitlines()[0][:140]
        for i, (name, k0, n_written, n_committed) in enumerate(stages):
            last = i == len(stages) - 1
            layout_of_last = last and (phase == "layout" or repaired) and n_committed == 0
            if name not in names:
                if layout_of_last:
                    return "TORN_LAYOUT", f"{name} not linked yet"
                return ("EARLIER_LOST" if not last else "LOST"), f"{name} missing"
            try:
                status, n, detail = _check_one_stage(f[name], k0, n_committed if last else n_written,
                                                     n_written)
            except Exception as e:  # noqa: BLE001
                status, detail = "READ_FAIL", str(e).splitlines()[0][:120]
            if status == "READ_FAIL" and layout_of_last:
                return "TORN_LAYOUT", f"{name}: {detail}"
            if status == "READ_FAIL" and not last:
                return "EARLIER_LOST", f"{name}: {detail}"
            if status != "OK":
                return status + suffix, f"{name}: {detail}"
    return "OK" + suffix, ""


def replay_resets(a):
    fobj = LogFile()
    write_run_resets(fobj, a.outputs, a.resets, a.presize)
    log = fobj.log
    nw = sum(1 for e in log if e[0] == "w")
    print(f"h5py {h5py.__version__} / HDF5 {h5py.version.hdf5_version}; outputs={a.outputs} "
          f"resets every {a.resets} outputs; presize={a.presize}; ops={len(log)} writes={nw} "
          f"final size={len(fobj.buf)}")
    tmp = tempfile.mkdtemp(prefix="replay_reset_")
    path = os.path.join(tmp, "k.ladruno")
    image = bytearray()
    phase, region = "layout", "create"
    stages = []
    started = False
    stats = collections.Counter()
    by_region = collections.Counter()
    fails = []
    for i, e in enumerate(log):
        if e[0] == "m":
            lab = e[1]
            if lab.startswith("phase:"):
                phase = lab[6:]
                if phase in ("closing", "closed", "layout") and stages:
                    region = "reset"
            elif lab.startswith("stage:"):
                _, s, k = lab.split(":")
                stages.append([f"MODEL_STAGE[{s}]", int(k), 0, 0])
                region = "layout"
            elif lab.startswith("out"):
                k = int(lab[3:].split(":")[0])
                kind = lab.split(":")[1]
                region = kind
                if kind == "rows":
                    stages[-1][2] = k - stages[-1][1] + 1
                    if k >= a.from_output:
                        started = True
                elif kind == "done":
                    stages[-1][3] = k - stages[-1][1] + 1
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
        if not started:
            continue
        snap = [list(s) for s in stages]
        with open(path, "wb") as fh:
            fh.write(image)
        status, detail = check_resets(path, phase, snap)
        if phase == "layout":
            # sec2 writes WRITE_ACCESS lazily: the image may still carry 0x00 (a plain open)
            with open(path, "wb") as fh:
                fh.write(image)
            s2, d2 = check_resets(path, phase, snap, flags=0x00)
            if not s2.startswith("OK") and (status.startswith("OK") or s2 != "TORN_LAYOUT"):
                status, detail = s2 + "(0x00)", d2
        stats[status] += 1
        by_region[(region, phase, status)] += 1
        # TORN_LAYOUT(0x00): the last stage torn after a PLAIN open (lazy WRITE_ACCESS);
        # as acceptable as after the step-3 repair: it holds no committed output
        good = status.startswith("OK") or status.split("(")[0] in ("NO_FILE", "TORN_LAYOUT")
        if not good:
            fails.append((i, len(stages) - 1, region, phase, status, detail))
    print("kill points checked:", sum(stats.values()), dict(stats))
    print("by region (the mark before the kill point) / sec2 phase / status:")
    for (reg, ph, status), c in sorted(by_region.items()):
        print(f"  {reg:9s} {ph:8s} {status:22s} {c}")
    for i, s, reg, ph, status, detail in fails[: a.show]:
        print(f"  FAIL op#{i} stage={s} region={reg} phase={ph} {status} {detail}")
    bad = len(fails)
    print(f"RESULT: {bad} bad kill points of {sum(stats.values())}")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", type=int, default=300)
    ap.add_argument("--evict-off", action="store_true",
                    help="disable metadata-cache eviction (build 8's change; no effect measured)")
    ap.add_argument("--cache-init", type=int, default=8 << 20, help="initial cache with --evict-off")
    ap.add_argument("--from-output", type=int, default=0,
                    help="check only kill points at or after this output's first write "
                         "(default 0: every point after the stage layout)")
    ap.add_argument("--show", type=int, default=12)
    ap.add_argument("--rpc", type=int, default=0, help="outputs per ID/DATA chunk, uncompressed (0 = current layout)")
    ap.add_argument("--swmr", action="store_true", help="latest format + SWMR (the current writer)")
    ap.add_argument("--resets", type=int, default=0,
                    help="WP-232h: an `mpm reset` (a new MODEL_STAGE, laid out without SWMR) "
                         "every N outputs; the writer since cbe1b7149 (SWMR), INFO group included")
    ap.add_argument("--no-presize", dest="presize", action="store_false",
                    help="with --resets: the writer before WP-232h (root header not presized, "
                         "no flush after the reset's reopen)")
    a = ap.parse_args()
    if a.resets:
        return replay_resets(a)

    fobj = LogFile()
    write_run(fobj, a.outputs, a.evict_off, a.cache_init, a.rpc, a.swmr)
    log = fobj.log
    nw = sum(1 for e in log if e[0] == "w")
    print(f"h5py {h5py.__version__} / HDF5 {h5py.version.hdf5_version}; outputs={a.outputs} "
          f"evict_off={a.evict_off} rpc={a.rpc} swmr={a.swmr}; ops={len(log)} writes={nw} "
          f"final size={len(fobj.buf)}")

    tmp = tempfile.mkdtemp(prefix="replay_")
    path = os.path.join(tmp, "k.ladruno")
    image = bytearray()
    phase, done = "start", 0
    started = False   # kill points while the file and stage are created hold nothing committed
    stats = collections.Counter()
    by_phase = collections.Counter()
    fails = []
    lag = collections.Counter()
    for i, e in enumerate(log):
        if e[0] == "m":
            phase = e[1]
            if phase.endswith(":done"):
                done = int(phase[3:].split(":")[0]) + 1
            if phase == f"out{a.from_output}:rows":
                started = True
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
        if not started:
            continue
        with open(path, "wb") as fh:
            fh.write(image)
        status, n, detail = check(path)
        stats[status] += 1
        kind = phase.split(":")[-1] if phase.startswith("out") else phase
        by_phase[(kind, status)] += 1
        if status == "OK":
            lag[done - n] += 1
        else:
            fails.append((i, phase, status, n, done, detail))
    print("kill points checked:", sum(stats.values()), dict(stats))
    print("by phase (the mark before the kill point: rows = between flushes before the barrier, "
          "barrier = inside the barrier flush, offsets = between flushes, commit = inside the "
          "commit flush):")
    for (kind, status), c in sorted(by_phase.items()):
        print(f"  {kind:14s} {status:10s} {c}")
    print("OK points: completed outputs minus committed:", dict(sorted(lag.items())))
    for i, phase, status, n, done, detail in fails[: a.show]:
        print(f"  FAIL op#{i} after-mark={phase} {status} committed={n} completed={done} {detail}")


if __name__ == "__main__":
    sys.exit(main())
