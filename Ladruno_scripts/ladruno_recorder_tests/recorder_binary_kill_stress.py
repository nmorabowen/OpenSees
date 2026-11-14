"""Kill stress of ``recorder ladruno`` on the REAL binary (WP-233, #1082).

What it shows: whether a run of ``dist/bin/opensees.pyd`` killed at a random moment always
leaves a ``.ladruno`` that a reader can read by schema §7b.2, holding exactly the values the
run computed. Each trial starts a child deck, waits a random time, kills it with
TerminateProcess (``Popen.kill``), and checks the file against a REFERENCE: the same deck run
to completion and closed (the analysis is deterministic, so a killed run's rows are a prefix
of the reference's, bit for bit).

The deck: a 2-D cantilever of 60 ``elasticBeamColumn`` (61 nodes, lumped masses, a sine tip
load, Newmark), ``recorder ladruno -N displacement velocity -E force``, one step per
``analyze(1)``. ``--mode`` drives the recorder's lifecycle windows (the points where it leaves
SWMR to create objects), so a random kill lands in them often:

- ``steady``: one stage, only appends (the steady state);
- ``stages``: ``domainChange`` every ``--every`` steps: a new MODEL_STAGE each time (a close,
  a read-write reopen, a layout, SWMR armed again);
- ``kind``: every ``--every`` steps ``domainChange``, ``wipeAnalysis``, an explicit
  ``record`` (the stage opens with no analysis: a provisional KIND), then the analysis is
  defined again: the next step rewrites KIND (a mid-stage layout step);
- ``modes``: ``-N modesOfVibration`` too, and ``eigen 2`` every ``--every`` steps: a modes
  snapshot (``DATA/STEP_<k>/MODE_<j>``) each time (a mid-stage layout step);
- ``envelope``: ``-envelope`` and ``domainChange`` every ``--every`` steps.

The check (``check``):

1. the reader rule: ``ladruno_format.open_ladruno`` (a plain open, else SWMR-read); when it
   raises ``StageLayoutKilledError`` (flags 0x01: killed inside a layout step) the file is
   recovered on a COPY with ``ladruno_format.repair_copy``;
2. every stage of the killed file is a stage of the reference, in order;
3. per stage, the crash rule's ``committed_steps`` T': every streaming result's DATA, TIME and
   STEP rows ``0 .. T'-1`` equal the reference's exactly; a stage that ENDED (not the last)
   holds all of the reference's steps; modes snapshots equal the reference's; KIND is one of
   static / transient / eigen; envelopes equal the reference's envelope over SOME prefix of
   the stage's steps, entry by entry;
4. only the LAST stage may fail to open (it was being laid out and holds no committed step;
   the open may be any step of the rule, since sec2 writes WRITE_ACCESS to the superblock
   lazily): ``TORN_LAST``, not a failure.

Statuses: ``OK`` (plain / swmr / repaired), ``TORN_LAST``, ``CHILD_EXITED`` (the child
finished before the kill; harmless), and the failures ``UNREADABLE``, ``MISMATCH``,
``LOST``, ``EXTRA_STAGE``. Failing trials keep their file under ``<out>/t<nnnn>/``.

Run from the repo root after a full ``build.bat`` (the pyd must be current):

    python3.12 Ladruno_scripts/ladruno_recorder_tests/recorder_binary_kill_stress.py \\
        --mode stages --every 20 --flush 0 --trials 300 --workers 4 --out <dir>
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import random
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DISTBIN = ROOT / "dist" / "bin"
sys.path.insert(0, str(HERE))
import ladruno_format as lf  # noqa: E402

CHILD = r"""
import os, sys
pyd, out, mode, every, flush, compress, nsteps = sys.argv[1:8]
every, nsteps = int(every), int(nsteps)
os.add_dll_directory(pyd); sys.path.insert(0, pyd)
import opensees as ops
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)
N = 60
for i in range(N + 1):
    ops.node(i + 1, 0.0, 0.1 * i)
    if i:
        ops.mass(i + 1, 1.0, 1.0, 0.01)
ops.fix(1, 1, 1, 1)
ops.geomTransf("Linear", 1)
for i in range(N):
    ops.element("elasticBeamColumn", i + 1, i + 1, i + 2, 0.01, 2.0e8, 1.0e-5, 1)
req = ["-N", "displacement", "velocity", "-E", "force"]
if mode == "modes":
    req += ["-N", "modesOfVibration"]
opt = ["-flush", float(flush), "-compress", int(compress)]
if mode == "envelope":
    opt += ["-envelope"]
ops.recorder("ladruno", out, *req, *opt)
ops.timeSeries("Sine", 1, 0.0, 1.0e9, 0.37)
ops.pattern("Plain", 1, 1)
ops.load(N + 1, 10.0, 0.0, 0.0)
def chain():
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
chain()
print("CHILD START", flush=True)
for k in range(1, nsteps + 1):
    if ops.analyze(1, 0.001) != 0:
        raise SystemExit("analyze failed")
    if every and k % every == 0:
        if mode in ("stages", "envelope"):
            ops.domainChange()
        elif mode == "kind":
            ops.domainChange(); ops.wipeAnalysis(); ops.record(); chain()
        elif mode == "modes":
            ops.eigen(2)
ops.wipe()
print("CHILD DONE", flush=True)
"""


def child_env():
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    env.pop("HDF5_USE_FILE_LOCKING", None)
    return env


def stage_names(f):
    return sorted((k for k in f if k.startswith("MODEL_STAGE[")), key=lambda s: int(s[12:-1]))


def _streams(st):
    return {g.name.split("/RESULTS/", 1)[1]: g for g in lf._streaming_results(st)}


def _envelopes(st):
    out = {}
    if "RESULTS/ENVELOPES" in st:
        def visit(name, obj):
            if isinstance(obj, h5py.Group) and "MIN" in obj:
                out[name] = obj
        st["RESULTS/ENVELOPES"].visititems(visit)
    return out


def check_stage(st, ref_st, last, a):
    """(status, T') for one stage, against the reference stage."""
    kind = lf._attr(st, "KIND") if "KIND" in st.attrs else None
    if kind not in ("static", "transient", "eigen"):
        return "MISMATCH", f"KIND={kind!r}"
    n = lf.committed_steps(st)
    if a.mode == "envelope":
        ref_streams = _streams(ref_st)
        for name, g in _envelopes(st).items():
            # the reference stream of the same result: ENVELOPES/<fam>/<path> -> <fam>/<path>
            r = ref_streams.get(name)
            if r is None:
                return "MISMATCH", f"envelope {name}: no reference stream"
            data = np.asarray(r["DATA"][...], dtype=np.float64)
            pmin = np.minimum.accumulate(data, axis=0)
            pmax = np.maximum.accumulate(data, axis=0)
            for key, pre in (("MIN", pmin), ("MAX", pmax)):
                v = g[key][...]
                if not np.all(np.any(pre == v[None], axis=0)):
                    return "MISMATCH", f"envelope {name}/{key} is no prefix extreme"
        return "OK", n
    if n is None:
        return "MISMATCH", "no streaming result"
    ref = _streams(ref_st)
    got = _streams(st)
    if set(got) != set(ref):
        return "MISMATCH", f"results {sorted(set(got) ^ set(ref))[:4]}"
    n_ref = lf.committed_steps(ref_st)
    if n > n_ref:
        return "MISMATCH", f"{n} steps > reference {n_ref}"
    if not last and n != n_ref:
        return "LOST", f"ended stage holds {n} of {n_ref} steps"
    for name, g in got.items():
        r = ref[name]
        for ds in ("DATA", "TIME", "STEP"):
            if not np.array_equal(g[ds][:n], r[ds][:n]):
                return "MISMATCH", f"{name}/{ds}"
    for mname in [k for k in st.get("RESULTS/ON_NODES", {}) if "MODES" in k.upper()]:
        mg, rg = st[f"RESULTS/ON_NODES/{mname}/DATA"], ref_st[f"RESULTS/ON_NODES/{mname}/DATA"]
        for sname in mg:
            try:
                snap = {m: mg[sname][m][...] for m in mg[sname]}
            except Exception:  # noqa: BLE001
                if last:
                    continue   # the snapshot being written when killed
                raise
            for m, v in snap.items():
                if not np.array_equal(v, rg[sname][m][...]):
                    return "MISMATCH", f"{mname}/{sname}/{m}"
    return "OK", n


def check(path, ref, a, work):
    """Returns (status, detail, how, n_stages, T_last)."""
    flags = lf.superblock_flags(path)
    how = "plain" if flags in (None, 0) else "swmr"
    try:
        f = lf.open_ladruno(path)
    except lf.StageLayoutKilledError:
        rep = os.path.join(work, "repaired.ladruno")
        lf.repair_copy(path, rep)
        f = h5py.File(rep, "r")
        how = "repaired"
    with f, h5py.File(ref, "r") as fr:
        names, rnames = stage_names(f), stage_names(fr)
        if names != rnames[: len(names)]:
            return "EXTRA_STAGE", f"{names} vs {rnames[:len(names) + 1]}", how, len(names), -1
        t_last = -1
        for i, name in enumerate(names):
            last = i == len(names) - 1
            try:
                status, detail = check_stage(f[name], fr[name], last, a)
            except Exception as e:  # noqa: BLE001
                if last:   # schema 7b.2: the stage being laid out (any open path)
                    return "TORN_LAST", str(e).splitlines()[0][:160], how, len(names), -1
                return "UNREADABLE", f"{name}: " + str(e).splitlines()[0][:160], how, len(names), -1
            if status != "OK":
                return status, f"{name}: {detail}", how, len(names), -1
            t_last = detail
        return "OK", "", how, len(names), t_last


def run_child(out, a, nsteps, wait=None, mode=None):
    p = subprocess.Popen([sys.executable, "-c", CHILD, str(DISTBIN), out, mode or a.mode, str(a.every),
                          str(a.flush), str(a.compress), str(nsteps)],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL,
                         text=True, env=child_env())
    if wait is None:
        o, e = p.communicate()
        if p.returncode != 0 or "CHILD DONE" not in o:
            raise SystemExit(f"reference run failed:\n{o[-800:]}\n{e[-1500:]}")
        return "DONE"
    line = p.stdout.readline()
    if "CHILD START" not in line:
        p.kill()
        return "CHILD_FAILED: " + line + p.stderr.read()[-400:]
    time.sleep(wait)
    exited = p.poll() is not None
    p.kill()
    p.wait()
    return "CHILD_EXITED" if exited else "KILLED"


def trial(i, wait, a, ref):
    d = os.path.join(a.out, f"t{i:04d}")
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "killed.ladruno")
    st = run_child(out, a, a.steps, wait)
    if st != "KILLED":
        shutil.rmtree(d, ignore_errors=True)
        return {"i": i, "wait": round(wait, 3), "status": st.split(":")[0], "detail": st[:300]}
    try:
        status, detail, how, ns, tl = check(out, ref, a, d)
    except Exception as e:  # noqa: BLE001
        status, detail, how, ns, tl = "UNREADABLE", "".join(
            traceback.format_exception_only(type(e), e)).strip()[-400:], "-", -1, -1
    flags = lf.superblock_flags(out)
    r = {"i": i, "wait": round(wait, 3), "status": status, "how": how, "flags": flags,
         "stages": ns, "T_last": tl, "size": os.path.getsize(out), "detail": detail}
    if status in ("OK", "TORN_LAST"):
        shutil.rmtree(d, ignore_errors=True)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=("steady", "stages", "kind", "modes", "envelope"), default="steady")
    ap.add_argument("--every", type=int, default=20, help="steps between lifecycle events")
    ap.add_argument("--flush", type=float, default=0.0, help="-flush seconds (0 = every step)")
    ap.add_argument("--compress", type=int, default=1)
    ap.add_argument("--steps", type=int, default=20000, help="steps of the reference / the child")
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--min-wait", type=float, default=0.3)
    ap.add_argument("--max-wait", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=233)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if not (DISTBIN / "opensees.pyd").exists():
        sys.exit(f"no binary in {DISTBIN}: run Ladruno_scripts\\build.bat first")
    os.makedirs(a.out, exist_ok=True)
    ref = os.path.join(a.out, "reference.ladruno")
    t0 = time.time()
    # the reference is always a time series (an envelope run is checked against the
    # prefix extremes of the same deck's streaming rows)
    run_child(ref, a, a.steps, mode="stages" if a.mode == "envelope" else None)
    with h5py.File(ref, "r") as fr:
        rn = stage_names(fr)
    print(f"reference: {len(rn)} stages, {time.time() - t0:.1f} s for {a.steps} steps "
          f"({os.path.getsize(ref) / 1e6:.1f} MB)", flush=True)
    rng = random.Random(a.seed)
    waits = [rng.uniform(a.min_wait, a.max_wait) for _ in range(a.trials)]
    results = []
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(trial, i, w, a, ref) for i, w in enumerate(waits)]
        for k, fu in enumerate(cf.as_completed(futs), 1):
            r = fu.result()
            results.append(r)
            if r["status"] not in ("OK", "TORN_LAST", "CHILD_EXITED"):
                print("FAIL", json.dumps(r), flush=True)
            if k % 50 == 0:
                bad = sum(1 for x in results if x["status"] not in ("OK", "TORN_LAST", "CHILD_EXITED"))
                print(f"progress {k}/{a.trials} failures={bad} elapsed={time.time() - t0:.0f}s", flush=True)
    results.sort(key=lambda r: r["i"])
    with open(os.path.join(a.out, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    c = Counter(r["status"] for r in results)
    how = Counter(r.get("how") for r in results if r["status"] in ("OK", "TORN_LAST"))
    bad = sum(v for s, v in c.items() if s not in ("OK", "TORN_LAST", "CHILD_EXITED"))
    print(f"SUMMARY mode={a.mode} every={a.every} flush={a.flush} compress={a.compress}: "
          f"{dict(c)}; opened {dict(how)}; failures {bad} of {len(results)}; "
          f"{time.time() - t0:.0f} s", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
