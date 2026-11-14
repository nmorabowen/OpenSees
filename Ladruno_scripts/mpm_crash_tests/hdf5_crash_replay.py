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

HDF5 2.0 + fileobj quirk: under SWMR the fileobj driver writes only the 0x04 consistency
flag into the superblock, where the sec2 driver (the real writer) writes 0x05, and an
SWMR-read open requires 0x05. So when a rebuilt image is refused as "already open for
write", the model sets the flags to 0x05 and re-checksums the version-3 superblock with
HDF5's lookup3 hash before the SWMR-read open (``mc.set_superblock``).

Run from anywhere (no OpenSees binary needed; ~minutes per few hundred outputs):

    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 300
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 300 --evict-off
    python3.12 Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py --outputs 600 --swmr
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
    a = ap.parse_args()

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
    main()
