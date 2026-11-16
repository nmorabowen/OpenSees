"""WP-232g (M0b): run the measurement scripts against one build, in sequence.

    python3.12 Ladruno_scripts/mpm_m0b/run_all.py --smoke                 # N = 10^3, everything
    python3.12 Ladruno_scripts/mpm_m0b/run_all.py --plan full --dry-run   # the run list + estimate
    python3.12 Ladruno_scripts/mpm_m0b/run_all.py --plan full             # the orchestrator's plan
    python3.12 Ladruno_scripts/mpm_m0b/run_all.py --only m1,m3            # each script's own default

Each plan item is one script invocation (its own results directory under <out>/), each run in it
a fresh child. `--plan full` is the orchestrator's run decision of 2026-11-15 (it moves no
threshold of the pre-registered table). The orchestrator launches it; nothing else heavy may run
on the machine meanwhile.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import glob
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C       # noqa: E402

SCRIPTS = {"m1": "m1_rebuild.py", "m2": "m2_explicit.py", "m3": "m3_memory.py", "m4": "m4_cantilever.py"}

ALL_N = "1000,10000,100000"
PAIRS = "qs,transient"
MATS = "elastic,dp"

# (name, script, arguments) -- the full plan, item 4 of the round-2 message
PLAN_FULL = [
    ("m1_pardiso_renormalize", "m1", ["--solvers", "Pardiso", "--treatments", "renormalize",
                                      "--N", ALL_N, "--modes", PAIRS, "--mats", MATS, "--repeats", "3"]),
    ("m1_umfpack_renormalize", "m1", ["--solvers", "UmfPack", "--treatments", "renormalize",
                                      "--N", ALL_N, "--modes", PAIRS, "--mats", MATS, "--repeats", "3",
                                      "--repeats-by-N", "100000=1", "--timeout-s", "2700"]),
    ("m1_pardiso_ghost", "m1", ["--solvers", "Pardiso", "--treatments", "ghost", "--N", ALL_N,
                                "--modes", PAIRS, "--mats", MATS, "--repeats", "3",
                                "--repeats-by-N", "100000=1"]),
    ("m1_pardiso_exclude", "m1", ["--solvers", "Pardiso", "--treatments", "exclude",
                                  "--N", "1000,10000", "--modes", PAIRS, "--mats", MATS, "--repeats", "1"]),
    ("m2", "m2", ["--repeats", "3"]),
    ("m3", "m3", ["--N", ALL_N, "--modes", "A,F", "--mats", MATS, "--repeats", "2"]),
    ("m4", "m4", ["--cases", "r1,r2,r3", "--cutoffs", "5e-2,1e-6,1e-10", "--repeats", "1"]),
]

# Wall-time scaling for --dry-run, from the smoke's per-run child wall at N = 10^3:
#   t(N) = t_fixed + (t_smoke - t_fixed) * (N / N_smoke)^p * (steps / smoke steps)
# p: implicit with PARDISO 1.5 (on the pre-232g plumbing build, the quasi-static DP block went
# from ~1.5 s at 10^3 to 41 s at 10^4: x27 for x10.6 in N); UmfPack 2.0 (3D fill of a direct
# solver without nested dissection); explicit 1.0; memory runs as their mode; the M2 driver and
# the cantilever have no N (the cantilever scales by its own particle count, 640 / 2560 / 10240,
# with p = 1.5). Only the stepping scales (the Python wall of the smoke's analyze calls); the rest
# of the smoke child's wall is fixed. Where the smoke has no stepping wall (M3), t_fixed = 1 s.
EXPONENT = {("m1", "Pardiso"): 1.5, ("m1", "UmfPack"): 2.0, ("m2", None): 1.0,
            ("m3", "A"): 1.5, ("m3", "F"): 1.0, ("m4", None): 1.5}
T_FIXED = 1.0
CANTILEVER_N = {"cantilever_qs_r1": 640, "cantilever_qs_r2": 2560, "cantilever_qs_r3": 10240}


def script_cmd(w, args, out, extra):
    cmd = [sys.executable, str(HERE / SCRIPTS[w]), "--pyd", args.pyd, "--out", str(out)] + extra
    if args.smoke:
        cmd.append("--smoke")
    if args.repeats and "--repeats" not in extra:
        cmd += ["--repeats", str(args.repeats)]
    return cmd


def list_of(cmd):
    out = subprocess.run(cmd + ["--list"], capture_output=True, text=True, stdin=subprocess.DEVNULL)
    for line in out.stdout.splitlines():
        if line.startswith("M0B-LIST "):
            return json.loads(line[len("M0B-LIST "):])
    raise RuntimeError("no run list from " + " ".join(cmd) + "\n" + out.stdout + out.stderr)


# --------------------------------------------------------------------------- #
#  the estimate                                                               #
# --------------------------------------------------------------------------- #
def smoke_reference(path):
    """{(script, key): child wall at the smoke} from a run_all smoke directory."""
    ref = {}
    for w in SCRIPTS:
        f = Path(path) / w / f"{w}_results.json"
        if not f.exists():
            continue
        doc = json.loads(f.read_text(encoding="utf-8"))
        steps = doc.get("pinned", {})
        for r in doc.get("runs", []):
            c = r.get("config") or {}
            ref[(w, _key(w, c))] = {"wall": r.get("child_wall_s"), "compute": _compute_s(w, r),
                                    "N": c.get("N"),
                                    "warmup": steps.get("warmup_steps"), "timed": steps.get("timed_steps"),
                                    "case": c.get("case")}
    return ref


def _compute_s(w, r):
    """The part of a smoke run that scales with N: the stepping (Python wall of the analyze
    calls); the rest of the child's wall (interpreter, import, build) is taken as fixed."""
    if w in ("m1", "m4"):
        v = sum(s.get("py_wall_s") or 0.0 for s in r.get("steps", []))
    elif w == "m2" and (r.get("config") or {}).get("kind") == "column":
        v = (r.get("warmup_py_wall_s") or 0.0) + ((r.get("timed") or {}).get("py_wall_s") or 0.0)
    else:
        v = None
    return v


def _key(w, c):
    drop = {"N", "pyd", "threads", "warmup", "timed", "increment", "increment_source", "repeats"}
    if w == "m4":
        drop |= {"case"}
    return json.dumps({k: v for k, v in sorted(c.items()) if k not in drop})


def estimate(w, row, ref):
    c = row["config"]
    r = ref.get((w, _key(w, c)))
    if r is None or r["wall"] is None:
        # fall back to any smoke run of the script at the same mode / kind
        cands = [v for (ww, k), v in ref.items() if ww == w and v["wall"]]
        if not cands:
            return None, "no smoke reference"
        r = max(cands, key=lambda v: v["wall"])
        note = "nearest smoke run (no exact match)"
    else:
        note = "smoke match"
    if w == "m1":
        p = EXPONENT[("m1", c["solver"])]
        n_ratio = c["N"] / (r["N"] or 1000)
        s_ratio = (row.get("warmup", 2) + max(row.get("timed", 6), 4 if c["N"] >= 100000 else 0)) / \
            ((r["warmup"] or 2) + (r["timed"] or 6))
    elif w == "m2":
        if c["kind"] == "driver":
            return r["wall"], note + " (driver: no N)"
        p, n_ratio = EXPONENT[("m2", None)], c["N"] / (r["N"] or 1000)
        s_ratio = (row.get("warmup", 5) + row.get("timed", 20)) / ((r["warmup"] or 2) + (r["timed"] or 5))
    elif w == "m3":
        p, n_ratio, s_ratio = EXPONENT[("m3", c["mode"])], c["N"] / (r["N"] or 1000), 1.0
    else:
        p, s_ratio = EXPONENT[("m4", None)], 1.0
        n_ratio = CANTILEVER_N[c["case"]] / CANTILEVER_N.get(r.get("case") or "cantilever_qs_r1", 640)
    comp = r.get("compute")
    if comp is None:
        comp = max(r["wall"] - T_FIXED, 0.05)
        fixed = T_FIXED
    else:
        fixed = max(r["wall"] - comp, 0.0)
    t = fixed + comp * (n_ratio ** p) * s_ratio
    to = row.get("timeout_s")
    if to and t > to:
        return to, note + f" (capped at the {to:.0f} s timeout)"
    return t, note


def latest_smoke():
    c = sorted(glob.glob(str(C.RESULTS / "run_all_smoke_*" / "index.json")), key=os.path.getmtime)
    return str(Path(c[-1]).parent) if c else None


# --------------------------------------------------------------------------- #
#  main                                                                       #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pyd", default=str(C.DEFAULT_PYD))
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--plan", choices=["full"], default=None,
                    help="run the orchestrator's full plan (see PLAN_FULL)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan's run list and a wall-time estimate; start nothing")
    ap.add_argument("--smoke-dir", default=None,
                    help="the smoke results the estimate scales from (default: the latest run_all_smoke_*)")
    ap.add_argument("--only", default=None, help="comma-separated subset of m1,m2,m3,m4 (or plan item names)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--repeats", type=int, default=None)
    args = ap.parse_args()
    args.pyd = str(Path(args.pyd).resolve())          # the children run in their results directories
    tag = "plan_full" if args.plan else ("run_all_smoke" if args.smoke else "run_all")
    root = Path(args.out).resolve() if args.out else \
        C.RESULTS / f"{tag}_{_dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"

    if args.plan:
        items = PLAN_FULL
    else:
        items = [(w, w, []) for w in SCRIPTS]
    if args.only:
        sel = {x.strip() for x in args.only.split(",")}
        items = [it for it in items if it[0] in sel or it[1] in sel]

    if args.dry_run:
        sd = args.smoke_dir or latest_smoke()
        ref = smoke_reference(sd) if sd else {}
        print(f"[dry-run] smoke reference: {sd}")
        print("[dry-run] scaling: t(N) = t_fixed + t_step_smoke (N / 10^3)^p (steps / smoke steps); "
              "p = 1.5 PARDISO / mode A, 2.0 UmfPack, 1.0 explicit, cantilever 1.5 by particle count; "
              "UmfPack capped at its timeout")
        grand, n_runs, plan_rows = 0.0, 0, []
        for name, w, extra in items:
            rows = list_of(script_cmd(w, args, root / name, extra))
            sub = 0.0
            for row in rows:
                t, note = estimate(w, row, ref)
                row["estimate_s"], row["estimate_note"] = t, note
                sub += t or 0.0
            grand += sub
            n_runs += len(rows)
            plan_rows.append({"item": name, "script": w, "args": extra, "runs": len(rows),
                              "estimate_h": sub / 3600.0, "rows": rows})
            print(f"[dry-run] {name:26s} {len(rows):4d} runs  ~{sub / 3600.0:7.2f} h")
            for row in rows:
                print(f"            {row['tag']:58s} ~{(row['estimate_s'] or 0) / 60.0:8.1f} min")
        print(f"[dry-run] TOTAL {n_runs} runs, ~{grand / 3600.0:.1f} h")
        return 0

    root.mkdir(parents=True, exist_ok=True)
    index = {"pyd": args.pyd, "smoke": args.smoke, "plan": args.plan, "items": {}}
    rc_all = 0
    for name, w, extra in items:
        cmd = script_cmd(w, args, root / name, extra)
        print(f"[run_all] {' '.join(cmd)}", flush=True)
        t0 = time.time()
        rc = subprocess.call(cmd, stdin=subprocess.DEVNULL)
        index["items"][name] = {"script": w, "args": extra, "rc": rc, "wall_s": time.time() - t0,
                                "dir": str(root / name)}
        rc_all = rc_all or rc
        (root / "index.json").write_text(json.dumps(index, indent=1), encoding="utf-8")
    print(f"[run_all] done: {root}")
    return rc_all


if __name__ == "__main__":
    sys.exit(main())
