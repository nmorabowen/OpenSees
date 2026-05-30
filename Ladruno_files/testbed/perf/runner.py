#!/usr/bin/env python3
"""Axis-3 performance micro-bench (Zone B). See the canonical doc §2 Axis 3.

Policy (locked):
  * metric  : state-determination work, reported as seconds/iteration (the
              portable Python proxy for ns/call; a C++ micro-driver is the
              eventual exact instrument).
  * env     : MKL/OMP threads pinned to 1, warmup discarded, MEDIAN-of-7.
  * gate    : baseline JSON keyed by (machine, bench). WARN at +10%, FAIL at
              +25% of the baseline median. Re-baseline ONLY via an explicit
              `perf: re-baseline <reason>` commit editing the JSON.

Usage:
  python runner.py --baseline        # record/overwrite this machine's baseline
  python runner.py --check           # compare vs baseline, exit 1 on FAIL
  python runner.py                    # just print medians
"""
import argparse
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

# pin threads BEFORE importing the solver
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

BASELINE_DIR = Path(__file__).resolve().parent / "baselines"
WARN_FRAC, FAIL_FRAC = 0.10, 0.25
REPEATS = 7
MACHINE = platform.node() or "unknown"


def _ops():
    try:
        import opensees as ops
    except ModuleNotFoundError:
        import openseespy.opensees as ops
    return ops


# --------------------------------------------------------------- benches
# Each bench returns seconds for one timed unit of work. Add a new entry per
# feature whose state-determination cost you want tracked. Keep the model
# fixed across versions so the number is comparable.

def bench_forceBeamColumn_chain(n_ele=200, inner=50):
    """Assemble+factor a cantilever chain `inner` times; seconds per pass."""
    ops = _ops()

    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.fix(1, 1, 1, 1)
        for i in range(2, n_ele + 2):
            ops.node(i, float(i - 1), 0.0)
        ops.section("Elastic", 1, 29000.0, 20.0, 800.0)
        ops.beamIntegration("Lobatto", 1, 1, 3)
        ops.geomTransf("Linear", 1)
        for i in range(1, n_ele + 1):
            ops.element("forceBeamColumn", i, i, i + 1, 1, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(n_ele + 1, 0.0, -10.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Newton")
        ops.analysis("Static", "-noWarnings")

    build()
    t0 = time.perf_counter()
    for _ in range(inner):
        ops.analyze(1)
    return (time.perf_counter() - t0) / inner


BENCHES = {"forceBeamColumn_chain": bench_forceBeamColumn_chain}


# --------------------------------------------------------------- harness
def measure(name, fn):
    fn()  # warmup, discarded
    samples = sorted(fn() for _ in range(REPEATS))
    return statistics.median(samples)


def baseline_path():
    return BASELINE_DIR / f"{MACHINE}.json"


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", action="store_true")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv[1:])

    medians = {name: measure(name, fn) for name, fn in BENCHES.items()}
    for name, m in medians.items():
        print(f"{name}: {m*1e3:.4f} ms/pass (median of {REPEATS}, machine {MACHINE})")

    if args.baseline:
        BASELINE_DIR.mkdir(exist_ok=True)
        baseline_path().write_text(json.dumps(medians, indent=2))
        print(f"baseline written: {baseline_path()}")
        return 0

    if args.check:
        if not baseline_path().exists():
            print(f"perf: no baseline for {MACHINE} — run --baseline first")
            return 0  # not a failure: a new machine has no reference yet
        base = json.loads(baseline_path().read_text())
        worst = 0
        for name, m in medians.items():
            b = base.get(name)
            if b is None:
                print(f"WARN  {name}: no baseline entry")
                continue
            frac = (m - b) / b
            tag = "OK" if frac < WARN_FRAC else ("WARN" if frac < FAIL_FRAC else "FAIL")
            print(f"{tag}  {name}: {frac*100:+.1f}% vs baseline {b*1e3:.4f} ms")
            if frac >= FAIL_FRAC:
                worst = 1
        return worst
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
