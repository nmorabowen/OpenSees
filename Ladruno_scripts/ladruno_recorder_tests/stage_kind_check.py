"""Per-stage KIND + stage-boundary gate (apeGmsh #2034) -- driver and checker.

Runs every case of stage_kind_model.tcl through the classic Tcl exe (and, when
--buildpy is given, the openseespy twin stage_kind_model.py), asserts each deck
reached its trailing marker (OpenSees.exe exits 0 on a Tcl error), and checks the
written .ladruno with h5py: the MODEL_STAGE list (KIND and recorded rows per
stage, in stamp order), strictly increasing stage ids, the shared validator, and
INFO/FORMAT_MINOR_VERSION == 1. With --buildpy, every case that runs in both
interpreters must also give IDENTICAL stage lists in Tcl and openseespy (the
classic `eigen` keeps the transient analysis it builds, openseespy deletes it:
the two drifted apart once, apeGmsh#2034 review).

    python stage_kind_check.py <dist_bin> <out_dir> [--buildpy <python.exe>]

<dist_bin> holds OpenSees.exe (and opensees.pyd for --buildpy); TCL_LIBRARY is
set to <dist_bin>\\..\\lib\\tcl8.6. Exit 0 iff every case passes. Run it against
a build without the #2034 fix and it must fail (dc_static_transient,
flip_static_transient, dc_static_static, transient_only, record_before_analysis).
"""

from __future__ import annotations

import os
import subprocess
import sys

import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ladruno_format as lf  # noqa: E402

# case -> expected [(KIND, recorded rows), ...] in MODEL_STAGE order
EXPECT = {
    "dc_static_transient": [("static", 4), ("transient", 5)],
    "flip_static_transient": [("static", 4), ("transient", 5)],
    "dc_static_static": [("static", 4), ("static", 3)],
    "repeat_static": [("static", 7)],
    "internal_restamp": [("static", 6)],
    "dc_first": [("static", 4)],
    "transient_only": [("transient", 5)],
    "record_before_analysis": [("transient", 6)],
    "override_kind": [("static", 4), ("static", 5)],
    "override_kind_dc": [("transient", 4), ("transient", 3)],
    "eigen_snapshot": [("static", 5)],
    "eigen_snapshot_transient": [("transient", 6)],
    "flip_then_topology": [("static", 4), ("transient", 5), ("transient", 2)],
}
PY_CASES = ("dc_static_transient", "flip_static_transient",
            "dc_static_static", "internal_restamp", "override_kind",
            "eigen_snapshot", "eigen_snapshot_transient")
MARKER = "STAGE_KIND_MODEL OK"


def _s(v):
    import numpy as np
    v = np.atleast_1d(v)[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def stages_of(path):
    with h5py.File(path, "r") as f:
        names = [k for k in f if k.startswith("MODEL_STAGE[")]
        ids = sorted(int(k[len("MODEL_STAGE["):-1]) for k in names)
        out = []
        for i in ids:
            g = f[f"MODEL_STAGE[{i}]"]
            d = g.get("RESULTS/ON_NODES/DISPLACEMENT/DATA")
            out.append((_s(g.attrs["KIND"]), 0 if d is None else int(d.shape[0])))
        minor = f["INFO"].attrs.get("FORMAT_MINOR_VERSION")
        minor = None if minor is None else int(minor[0])
    return ids, out, minor


def run_case(label, cmd, env, out, case, log):
    if os.path.exists(out):
        os.remove(out)
    with open(log, "w", encoding="utf-8", errors="replace") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env,
                       cwd=HERE, timeout=300)
    text = open(log, encoding="utf-8", errors="replace").read()
    probs = []
    if f"{MARKER} {case}" not in text:
        return [f"deck did not reach '{MARKER} {case}' (see {log})"]
    if not os.path.exists(out):
        return [f"no output file {out}"]
    ids, got, minor = stages_of(out)
    if got != EXPECT[case]:
        probs.append(f"stages {got} != expected {EXPECT[case]} (ids {ids})")
    if ids != sorted(set(ids)):
        probs.append(f"stage ids not strictly increasing: {ids}")
    if minor != 1:
        probs.append(f"INFO/FORMAT_MINOR_VERSION = {minor}, expected 1")
    v = lf.validate(out)
    if v:
        probs.append(f"validator: {v[:3]}")
    return probs


def main():
    args = sys.argv[1:]
    buildpy = None
    if "--buildpy" in args:
        k = args.index("--buildpy")
        buildpy = args[k + 1]
        del args[k:k + 2]
    dist_bin, out_dir = args[0], args[1]
    os.makedirs(out_dir, exist_ok=True)
    env = dict(os.environ)
    env["TCL_LIBRARY"] = os.path.join(dist_bin, "..", "lib", "tcl8.6")
    exe = os.path.join(dist_bin, "OpenSees.exe")
    fails = 0
    results = {}
    runs = [("tcl", c) for c in EXPECT]
    if buildpy:
        runs += [("py", c) for c in PY_CASES]
    for lang, case in runs:
        # forward slashes: Tcl reads a backslash in a path as an escape
        out = os.path.join(out_dir, f"stage_kind_{lang}_{case}.ladruno").replace("\\", "/")
        log = os.path.join(out_dir, f"stage_kind_{lang}_{case}.log")
        if lang == "tcl":
            cmd = [exe, os.path.join(HERE, "stage_kind_model.tcl"), out, case]
        else:
            cmd = [buildpy, os.path.join(HERE, "stage_kind_model.py"), dist_bin, out, case]
        probs = run_case(lang, cmd, env, out, case, log)
        if not probs:
            results[(lang, case)] = stages_of(out)[1]
        if probs:
            fails += 1
            print(f" FAIL [{lang}] {case}: " + "; ".join(probs))
        else:
            print(f"  ok  [{lang}] {case}: {EXPECT[case]}")
    if buildpy:
        for case in PY_CASES:
            t, p = results.get(("tcl", case)), results.get(("py", case))
            if t is None or p is None or t != p:
                fails += 1
                print(f" FAIL [tcl==py] {case}: tcl {t} vs py {p}")
            else:
                print(f"  ok  [tcl==py] {case}")
    if fails:
        print(f"STAGE_KIND CHECK FAILED ({fails} case(s))")
        return 1
    print("STAGE_KIND CHECK OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
