"""MPCO modes recorder after `eigen; wipe` (apeGmsh #1983) — runner + check.

Runs mpco_eigen_wipe_model.tcl under the classic Tcl exe and asserts:

1. the process exits 0 AND prints the deck's trailing marker. OpenSees.exe
   exits 0 even when a Tcl error ends the script, so the marker is the proof;
   the unfixed build exit(-1)s inside `analyze` in phase B and prints neither;
2. phase A (no wipe): eig_nowipe.mpco holds one modal step with MODE_0, and
   its LAMBDA is the closed form 3EI/(m L^3) = 6.0e7 (rel 1e-8), with a purely
   lateral eigenvector (Uy = 0 at the tip);
3. phase B (eigen -> wipe -> new model): eig_wipe.mpco holds no modal step.
   The new domain has no spectrum, so the recorder records no modes, as the
   Ladruno recorder already does (WP-163 R1).

Usage: python mpco_eigen_wipe_check.py <OpenSees.exe> <out_dir>

TCL_LIBRARY is set to <exe dir>/../lib/tcl8.6 when that directory exists and
the caller did not set it.
"""

from __future__ import annotations

import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DECK = os.path.join(HERE, "mpco_eigen_wipe_model.tcl")
MARKER = "MPCO_EIGEN_WIPE_MODEL OK"
MODES = "RESULTS/ON_NODES/MODES_OF_VIBRATION(U)"
LAMBDA_EXACT = 3.0 * 200.0e9 * 1.0e-4 / (1.0 * 1.0**3)   # 3EI/(m L^3)


def _modal_steps(path):
    """{stage/STEP_n: {MODE_k: (LAMBDA, data)}} for every modal step in the file."""
    out = {}
    with h5py.File(path, "r") as f:
        for stage in sorted(k for k in f if k.startswith("MODEL_STAGE")):
            data = f.get(f"{stage}/{MODES}/DATA")
            if data is None:
                continue
            for step in data:
                out[f"{stage}/{step}"] = {
                    # MPCO writes scalar attributes as 1-element arrays
                    m: (float(np.ravel(data[step][m].attrs["LAMBDA"])[0]), data[step][m][()])
                    for m in data[step]
                }
    return out


def main():
    exe, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    env = dict(os.environ)
    tcl_lib = os.path.normpath(os.path.join(os.path.dirname(exe), "..", "lib", "tcl8.6"))
    if "TCL_LIBRARY" not in env and os.path.isdir(tcl_lib):
        env["TCL_LIBRARY"] = tcl_lib

    proc = subprocess.run([exe, DECK, out_dir], capture_output=True, text=True, env=env)
    log = proc.stdout + proc.stderr
    with open(os.path.join(out_dir, "mpco_eigen_wipe_model.log"), "w") as fh:
        fh.write(log)

    ok = True
    if proc.returncode != 0:
        print(f"  FAIL exit code {proc.returncode} (expected 0)")
        ok = False
    if MARKER not in log:
        print(f"  FAIL deck did not reach '{MARKER}'")
        ok = False
    if not ok:
        print(log[-2000:])
        print("MPCO_EIGEN_WIPE CHECK FAILED")
        return 1
    print("  OK   deck ran to its marker, exit 0")

    steps = _modal_steps(os.path.join(out_dir, "eig_nowipe.mpco"))
    if len(steps) != 1:
        print(f"  FAIL phase A: expected 1 modal step, got {sorted(steps)}")
        ok = False
    else:
        (step, modes), = steps.items()
        if sorted(modes) != ["MODE_0"]:
            print(f"  FAIL phase A: expected [MODE_0] in {step}, got {sorted(modes)}")
            ok = False
        else:
            lam, vec = modes["MODE_0"]
            rel = abs(lam - LAMBDA_EXACT) / LAMBDA_EXACT
            # rows follow the ID dataset (nodes 1, 2); node 2 is the free tip
            tip = vec[-1]
            if rel > 1e-8 or abs(tip[0]) == 0.0 or abs(tip[1]) > 1e-12 * abs(tip[0]):
                print(f"  FAIL phase A: LAMBDA {lam!r} (exact {LAMBDA_EXACT!r}, rel {rel:.2e}), tip {tip}")
                ok = False
            else:
                print(f"  OK   phase A: {step}/MODE_0 LAMBDA {lam!r} = 3EI/(mL^3), lateral tip mode")

    steps = _modal_steps(os.path.join(out_dir, "eig_wipe.mpco"))
    if steps:
        print(f"  FAIL phase B: expected no modal step after wipe, got {sorted(steps)}")
        ok = False
    else:
        print("  OK   phase B: no modal step recorded after eigen -> wipe")

    if not ok:
        print("MPCO_EIGEN_WIPE CHECK FAILED")
        return 1
    print("MPCO_EIGEN_WIPE CHECK OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
