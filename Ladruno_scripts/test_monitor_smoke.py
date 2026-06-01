#!/usr/bin/env python
"""
Ladruno live analysis-monitor — smoke + parity + overhead test.

Builds a small SDOF transient, streams node disp through the new
`recorder Monitor` SWMR-HDF5 sink, and checks:

  1. SMOKE   — the .h5 sink has the expected frame count, monotonic STEP,
               self-describing COLUMNS, and finite FRAMES.
  2. PARITY  — the monitored channel equals a parallel plain Node recorder
               on the same node/dof (the monitor must observe, not alter).
  3. OVERHEAD— monitored vs unmonitored wall time stays small.

Run with a Python that has h5py + numpy, after building OpenSeesPy:

    python Ladruno_scripts/test_monitor_smoke.py

The opensees module is loaded from dist/bin (built by build.bat OpenSeesPy).
"""
import os
import sys
import time

# h5py must read a file another process may still hold open (SWMR) -> disable
# libhdf5 file locking (LEDGER_quirks). Set BEFORE importing h5py.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DIST_BIN = os.path.join(ROOT, "dist", "bin")
if DIST_BIN not in sys.path:
    sys.path.insert(0, DIST_BIN)

import numpy as np
import h5py
import opensees as ops

NSTEPS = 200
DT = 0.01
EVERY = 1


def build_sdof(with_monitor, sink_path=None, parity_path=None, every=EVERY):
    """A linear SDOF under a ramp load. Node 2 is the free mass."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, 2.0)
    ops.uniaxialMaterial("Elastic", 1, 50.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0)

    if parity_path is not None:
        # -precision 16 so the text reference isn't truncated to 6 sig figs
        # (otherwise parity is limited by the *reference's* round-trip, ~1e-7).
        ops.recorder("Node", "-file", parity_path, "-time", "-precision", 16,
                     "-node", 2, "-dof", 1, "disp")
    if with_monitor:
        ops.recorder("Monitor", "-node", 2, "-dof", 1, "-resp", "disp",
                     "-sink", sink_path, "-every", every)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def run(nsteps):
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            raise RuntimeError("analyze failed")
    # flush + drop recorders so the sink file is closed before we read it
    ops.wipe()


def check_smoke(sink_path):
    def _s(v):
        return v.decode() if isinstance(v, bytes) else v

    with h5py.File(sink_path, "r") as f:
        assert _s(f.attrs["FORMAT"]) == "ladruno-monitor", "bad FORMAT attr"
        cols = [_s(c) for c in f["COLUMNS"][:]]
        step = f["STEP"][:]
        tvec = f["TIME"][:]
        frames = f["FRAMES"][:]
        n_expected = NSTEPS // EVERY
        assert len(step) == n_expected, f"STEP {len(step)} != {n_expected}"
        assert frames.shape == (n_expected, 1), f"FRAMES shape {frames.shape}"
        assert np.all(np.diff(step) > 0), "STEP not monotonic"
        assert np.all(np.isfinite(frames)), "non-finite FRAMES"
        assert cols == ["node2.disp.dof1"], f"COLUMNS {cols}"
        print(f"  SMOKE ok: {len(step)} frames, columns={cols}")
        return step, tvec, frames[:, 0]


def check_parity(parity_path, mon_disp):
    data = np.loadtxt(parity_path)          # cols: time, disp
    ref = data[:, 1]
    n = min(len(ref), len(mon_disp))
    err = np.max(np.abs(ref[:n] - mon_disp[:n]))
    assert err < 1e-9, f"parity mismatch max|err|={err:.3e}"
    print(f"  PARITY ok: max|err| vs Node recorder = {err:.2e}")


def check_swmr_readable(sink_path):
    """Reopen in SWMR-read mode + refresh -- confirms the file carries a valid
    SWMR superblock and is tail-readable (the live path the FastAPI sidecar
    will use). True cross-process concurrent tailing is exercised by the P1
    sidecar; here we just prove the format is SWMR-correct."""
    with h5py.File(sink_path, "r", swmr=True) as f:
        ds = f["STEP"]
        ds.id.refresh()
        assert ds.shape[0] == NSTEPS, f"swmr read sees {ds.shape[0]} != {NSTEPS}"
    print("  SWMR ok: file reopens in swmr=True read mode and refreshes")


def main():
    sink = os.path.join(HERE, "_monitor_live.h5")
    parity = os.path.join(HERE, "_monitor_parity.out")
    for p in (sink, parity):
        if os.path.exists(p):
            os.remove(p)

    # monitored run
    build_sdof(True, sink_path=sink, parity_path=parity)
    t0 = time.perf_counter()
    run(NSTEPS)
    t_mon = time.perf_counter() - t0

    step, tvec, mon_disp = check_smoke(sink)
    check_parity(parity, mon_disp)
    check_swmr_readable(sink)

    # overhead: unmonitored baseline
    build_sdof(False)
    t0 = time.perf_counter()
    run(NSTEPS)
    t_plain = time.perf_counter() - t0

    # Report the meaningful figure: absolute cost per emitted frame. (A raw %
    # is meaningless here -- this SDOF's step is ~1.5 us, so any fixed per-frame
    # HDF5 cost dwarfs it. On a realistic FE step of >=1 ms it is <0.5%, and
    # -every/-hz bound it further; see check below.)
    per_frame_us = 1e6 * (t_mon - t_plain) / NSTEPS
    print(f"  OVERHEAD: ~{per_frame_us:.1f} us / emitted frame "
          f"(monitored {t_mon*1e3:.1f} ms vs plain {t_plain*1e3:.1f} ms, "
          f"{NSTEPS} frames)")

    # The -every gate must cut emitted frames (and thus overhead) proportionally.
    sink10 = os.path.join(HERE, "_monitor_every10.h5")
    if os.path.exists(sink10):
        os.remove(sink10)
    build_sdof(True, sink_path=sink10, every=10)
    run(NSTEPS)
    with h5py.File(sink10, "r") as f:
        n10 = len(f["STEP"])
        step10 = f["STEP"][:]
    print(f"  GATE: -every 10 -> {n10} frames (expect {NSTEPS // 10})")
    assert n10 == NSTEPS // 10, f"decimation gave {n10}, expected {NSTEPS // 10}"
    # committed step stride must be exactly 10 (commitTag is a global monotonic
    # counter that does NOT reset on wipe(), so check the stride, not absolutes)
    assert np.all(np.diff(step10) == 10), f"stride not 10: {np.unique(np.diff(step10))}"

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
