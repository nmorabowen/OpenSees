"""WP-232g (M0b): what every measurement script shares.

* the binary: `--pyd` (default: this worktree's dist/bin/opensees.pyd), imported in a CHILD
  process only, and checked to be the file that was asked for;
* every run is one fresh child process (`<script> --child <config.json>`), started with the
  thread environment pinned (OMP_NUM_THREADS / MKL_NUM_THREADS) before opensees is imported;
* the child writes a result JSON and, when profiled, one profiler HDF5 report per step
  (`profiler report <h5> -run s####`); the PARENT reads the HDF5 with h5py. The child never
  imports h5py: opensees.pyd and h5py each bring an HDF5 library, and two in one process is a
  version clash waiting to happen;
* memory through ctypes GetProcessMemoryInfo (PrivateUsage), since psutil is not installed;
* the profiler tree helpers and the M1 decomposition (t_rb, t_it, R, k, f) exactly as the
  232g brief defines them (.claude/skills/ladruno-mpm/references/wp-briefs.md, "232g").
"""
from __future__ import annotations

import argparse
import ctypes
import datetime as _dt
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_PYD = REPO / "dist" / "bin" / "opensees.pyd"
RESULTS = HERE / "results"
HOST = 232                     # the MPM host pattern tag (as tests/ladruno_mpm_implicit_harness.py)

# The pinned machine (232g "Pinned conditions"); the harness also RECORDS what it sees.
PINNED_MACHINE = "i7-13700KF (8 P + 8 E cores), 64 GB"


# --------------------------------------------------------------------------- #
#  Parent side: arguments, children, results                                  #
# --------------------------------------------------------------------------- #
P_CORES_MASK = 0xFFFF          # i7-13700KF: logical 0-15 are the 8 P-cores (HT); 16-23 E-cores


def _csv(t):
    return lambda s: [t(x.strip()) for x in s.split(",") if x.strip()]


def base_parser(desc):
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--pyd", default=str(DEFAULT_PYD),
                    help="the opensees.pyd to measure (default: this worktree's dist/bin); a relative "
                         "path is taken from the current directory")
    ap.add_argument("--smoke", action="store_true",
                    help="N = 10^3, one repeat, short runs (the round-1 smoke)")
    ap.add_argument("--repeats", type=int, default=None, help="repeats per configuration (default 3)")
    ap.add_argument("--repeats-by-N", dest="repeats_by_N", default=None,
                    help="per-N repeats overriding --repeats, e.g. 100000=1")
    ap.add_argument("--out", default=None, help="results directory (default results/<script>_<stamp>)")
    ap.add_argument("--only", default=None,
                    help="comma-separated key=value filters on the configuration (e.g. N=1000,solver=Pardiso)")
    # subset flags (comma lists); a flag a script's configurations do not carry is ignored
    ap.add_argument("--N", type=_csv(int), default=None, help="particle counts, e.g. 1000,10000")
    ap.add_argument("--modes", type=_csv(str), default=None, help="M1: qs,transient; M3: A,F")
    ap.add_argument("--mats", type=_csv(str), default=None, help="elastic,dp[,j2,sanisand]")
    ap.add_argument("--solvers", type=_csv(str), default=None, help="Pardiso,UmfPack")
    ap.add_argument("--treatments", type=_csv(str), default=None, help="renormalize,exclude,ghost")
    ap.add_argument("--timeout-s", "--timeout", dest="timeout", type=float, default=6 * 3600.0,
                    help="per-run timeout [s]; a timeout is recorded (ok=false, reason=timeout)")
    ap.add_argument("--affinity", default=hex(P_CORES_MASK),
                    help="process affinity mask of every child (default 0xFFFF, the P-cores); 'none' "
                         "leaves it")
    ap.add_argument("--list", action="store_true",
                    help="print the selected runs as JSON and exit (run_all --dry-run uses it)")
    ap.add_argument("--child", default=None, help=argparse.SUPPRESS)
    return ap


def parse_args(ap):
    args = ap.parse_args()
    # resolve in the PARENT: the child runs in the results directory (review fix, round 2)
    args.pyd = str(Path(args.pyd).resolve())
    if args.out:
        args.out = str(Path(args.out).resolve())
    args.affinity_mask = None if str(args.affinity).lower() == "none" else int(str(args.affinity), 0)
    return args


_FLAG_KEYS = (("N", "N"), ("modes", "mode"), ("mats", "mat"), ("solvers", "solver"),
              ("treatments", "treatment"))


def select(cfgs, args):
    """Apply the subset flags and --only to a configuration list."""
    out = []
    for c in cfgs:
        ok = True
        for flag, key in _FLAG_KEYS:
            vals = getattr(args, flag, None)
            if key == "N" and c.get("N") == 0:      # N-free runs (the M2 single-point driver)
                continue
            if vals is not None and key in c and c[key] not in vals:
                ok = False
        if ok and keep(c, args.only):
            out.append(c)
    return out


def repeats_for(c, args, default):
    """--repeats, overridden per N by --repeats-by-N (e.g. 100000=1)."""
    n = args.repeats or default
    for kv in (args.repeats_by_N or "").split(","):
        if "=" in kv:
            k, _, v = kv.partition("=")
            if "N" in c and int(k) == c["N"]:
                n = int(v)
    return n


def list_runs(cfgs, args, script, tag_of, default_reps, extra=None):
    """--list: every run this invocation would start, as one JSON line on stdout."""
    rows = [{"script": script, "tag": tag_of(c, r), "config": c, "timeout_s": args.timeout,
             **(extra or {})}
            for c in cfgs for r in range(repeats_for(c, args, default_reps))]
    print("M0B-LIST " + json.dumps(rows))
    return 0


def out_dir(args, name):
    if args.out:
        d = Path(args.out)
    else:
        d = RESULTS / f"{name}_{_dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def keep(cfg, only):
    """--only filter: every key=value must match str(cfg[key])."""
    if not only:
        return True
    for kv in only.split(","):
        k, _, v = kv.partition("=")
        if str(cfg.get(k.strip())) != v.strip():
            return False
    return True


def child_env(threads):
    """The child's environment with the thread pins. threads=None leaves both unset (the
    'default' row), recording that it did."""
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env.pop(k, None)
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
        env["MKL_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["LADRUNO_OPENSEES_QUIET"] = "1"          # no splash banner in the child logs
    for k in ("LADRUNO_MPM_EMPTY_NODES", "LADRUNO_MPM_MASSCUT", "LADRUNO_MPM_ALGORITHM"):
        env.pop(k, None)                         # the 232f harness's data-run switches: never inherited
    return env


def run_child(script, cfg, workdir, tag, timeout, affinity=P_CORES_MASK):
    """One configuration in a fresh process. Returns (result dict, log path). A timeout or a
    crash is a result too (ok=false with its reason and wall time), never a dropped run."""
    workdir = Path(workdir)
    cfg = dict(cfg)
    cfg["_result"] = str(workdir / f"{tag}.json")
    cfg["_h5"] = str(workdir / f"{tag}.h5")
    cfg["_affinity"] = affinity
    cfgp = workdir / f"{tag}.cfg.json"
    cfgp.write_text(json.dumps(cfg, indent=1), encoding="utf-8")
    for p in (cfg["_result"], cfg["_h5"]):
        if os.path.exists(p):
            os.remove(p)
    log = workdir / f"{tag}.log"
    t0 = time.perf_counter()
    timed_out = False
    with open(log, "w", encoding="utf-8", errors="replace") as fh:
        try:
            proc = subprocess.run([sys.executable, "-X", "utf8", str(script), "--child", str(cfgp)],
                                  stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                  env=child_env(cfg.get("threads")), cwd=str(workdir),
                                  timeout=timeout)
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            rc, timed_out = "timeout", True
    wall = time.perf_counter() - t0
    res = None
    if os.path.exists(cfg["_result"]):
        try:
            res = json.loads(Path(cfg["_result"]).read_text(encoding="utf-8"))
        except ValueError:
            res = None
    if res is None:
        res = {"ok": False, "config": {k: v for k, v in cfg.items() if not k.startswith("_")},
               "error": f"child exited {rc} without a result (see {log.name})"}
    if timed_out:
        res["ok"] = False
        res["reason"] = "timeout"
        res["timeout_s"] = timeout
    elif not res.get("ok") and "reason" not in res:
        res["reason"] = "error"
    res["child_rc"] = rc
    res["child_wall_s"] = wall
    res["log"] = log.name
    return res, log


def safe_read_profile(path):
    """read_profile, but a file a killed child left unreadable is recorded, not fatal."""
    try:
        return read_profile(path), None
    except Exception as e:                    # noqa: BLE001 -- an unreadable file is data
        return {}, repr(e)


def parent_conditions(args):
    return {
        "pinned_machine": PINNED_MACHINE,
        "machine_seen": {"processor": platform.processor(), "cpu_count": os.cpu_count(),
                         "total_ram_bytes": total_ram_bytes(), "platform": platform.platform()},
        "python": sys.version.split()[0],
        "pyd": str(Path(args.pyd).resolve()),
        "smoke": bool(args.smoke),
        "affinity_requested": None if args.affinity_mask is None else hex(args.affinity_mask),
        "timeout_s": args.timeout,
        "started": _dt.datetime.now().isoformat(timespec="seconds"),
        "note": "no build or other heavy job may run during a timed run (232g pinned conditions)",
    }


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default), encoding="utf-8")


def _json_default(o):
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, set):
        return sorted(o)
    raise TypeError(type(o))


def median_spread(vals):
    """{median, min, max, spread = (max-min)/|median|, n} over the finite values."""
    v = [x for x in vals if x is not None and isinstance(x, (int, float)) and math.isfinite(x)]
    if not v:
        return {"median": None, "min": None, "max": None, "spread": None, "n": 0}
    m = statistics.median(v)
    return {"median": m, "min": min(v), "max": max(v),
            "spread": (max(v) - min(v)) / abs(m) if m else None, "n": len(v)}


def loglog_slope(xs, ys):
    """Least-squares slope of log y against log x (None with fewer than two usable points)."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys)
           if x and y and x > 0 and y > 0 and math.isfinite(y)]
    if len(pts) < 2:
        return None
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx == 0:
        return None
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx


# --------------------------------------------------------------------------- #
#  Memory (both sides; ctypes, no psutil)                                     #
# --------------------------------------------------------------------------- #
class _PMC_EX(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t), ("PrivateUsage", ctypes.c_size_t)]


def process_memory():
    """{private, working_set, peak_working_set, peak_pagefile} in bytes (Windows), else {}."""
    if os.name != "nt":
        return {}
    c = _PMC_EX()
    c.cb = ctypes.sizeof(c)
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    try:
        fn = k32.K32GetProcessMemoryInfo
    except AttributeError:
        fn = ctypes.WinDLL("psapi").GetProcessMemoryInfo
    fn.argtypes = [ctypes.c_void_p, ctypes.POINTER(_PMC_EX), ctypes.c_ulong]
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    if not fn(k32.GetCurrentProcess(), ctypes.byref(c), c.cb):
        return {}
    return {"private": int(c.PrivateUsage), "working_set": int(c.WorkingSetSize),
            "peak_working_set": int(c.PeakWorkingSetSize), "peak_pagefile": int(c.PeakPagefileUsage)}


def total_ram_bytes():
    if os.name != "nt":
        return None

    class MEMSTAT(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    m = MEMSTAT()
    m.dwLength = ctypes.sizeof(m)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
        return None
    return int(m.ullTotalPhys)


# --------------------------------------------------------------------------- #
#  Child side                                                                 #
# --------------------------------------------------------------------------- #
def load_cfg(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


_AFFINITY = {"requested": None, "process_mask": None, "system_mask": None, "set_ok": None}


def set_affinity(mask):
    """SetProcessAffinityMask on this (child) process, BEFORE opensees / MKL start threads;
    records the request and the mask the OS reports back."""
    if os.name != "nt":
        return dict(_AFFINITY)
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    h = k32.GetCurrentProcess()
    if mask is not None:
        k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        _AFFINITY["requested"] = hex(mask)
        _AFFINITY["set_ok"] = bool(k32.SetProcessAffinityMask(h, ctypes.c_size_t(mask)))
    pm, sm = ctypes.c_size_t(0), ctypes.c_size_t(0)
    k32.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t),
                                           ctypes.POINTER(ctypes.c_size_t)]
    if k32.GetProcessAffinityMask(h, ctypes.byref(pm), ctypes.byref(sm)):
        _AFFINITY["process_mask"], _AFFINITY["system_mask"] = hex(pm.value), hex(sm.value)
    return dict(_AFFINITY)


def import_ops(pyd, affinity=None):
    """Pin the affinity, then import opensees from exactly `pyd` (refuse anything else: an
    installed .pth or a stale tests/ copy must never be measured by mistake)."""
    set_affinity(affinity)
    pyd = Path(pyd).resolve()
    if not pyd.exists():
        raise SystemExit(f"no binary at {pyd}")
    sys.path.insert(0, str(pyd.parent))
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(pyd.parent))
    import opensees as ops      # noqa: E402  (the path above decides which one)
    got = Path(getattr(ops, "__file__", "")).resolve()
    if os.path.normcase(str(got)) != os.path.normcase(str(pyd)):
        raise SystemExit(f"imported opensees from {got}, asked for {pyd}")
    return ops


def child_conditions(ops):
    cond = {"ladrunoBuild": str(ops.ladrunoBuild()),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "default (unset)"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS", "default (unset)"),
            "python": sys.version.split()[0], "pid": os.getpid(), "affinity": dict(_AFFINITY)}
    try:                                            # the WP-107 element-loop threads (informative)
        cond["ladrunoThreads"] = ops.ladrunoThreads()
    except Exception as e:
        cond["ladrunoThreads"] = f"unavailable: {e}"
    return cond


class StepProfiler:
    """Coarse (or -deep / -memory) profiler run around exactly one call, reported as its own
    HDF5 run (`s0000`, `s0001`, ...). The report and the reset happen outside the call."""

    def __init__(self, ops, h5, flags=()):
        self.ops, self.h5, self.flags, self.k = ops, str(h5), list(flags), 0

    def __call__(self, fn):
        ops = self.ops
        ops.profiler("reset")
        ops.profiler("start", *self.flags)
        t0 = time.perf_counter()
        try:
            out = fn()
        finally:
            wall = time.perf_counter() - t0
            ops.profiler("stop")
            run = f"s{self.k:04d}"
            ops.profiler("report", self.h5, "-run", run)
            ops.profiler("reset")
            self.k += 1
        return out, wall, run


def write_child_result(cfg, res):
    write_json(cfg["_result"], res)


def particle_field(ops, f):
    v = ops.mpm("particleData", f)
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def info(ops, what):
    return int(ops.mpm("info", what))


# --------------------------------------------------------------------------- #
#  Profiler HDF5 (parent side)                                                #
# --------------------------------------------------------------------------- #
def read_profile(h5path):
    """{run id: tree}; tree = {"calls", "wall_ns", "children": {name: tree}, "elem": [...]}.
    Keeps the per-step rollup verbatim (minus min/max/cpu): the 'raw per-step totals'."""
    import h5py
    out = {}
    if not os.path.exists(h5path):
        return out

    def node(g):
        n = {"calls": int(g.attrs.get("calls", 0)), "wall_ns": int(g.attrs.get("wall_ns", 0)),
             "children": {}}
        for k in g:
            obj = g[k]
            if isinstance(obj, h5py.Group):
                n["children"][k] = node(obj)
            elif k == "elem_by_type":
                n["elem"] = [{"classTag": int(r["classTag"]), "count": int(r["count"]),
                              "wall_ns": int(r["wall_ns"])} for r in obj[()]]
        return n

    with h5py.File(h5path, "r") as f:
        runs = f.get("runs")
        if runs is None:
            return out
        for rid in runs:
            g = runs[rid]
            meta = {k: (v.decode() if isinstance(v, bytes) else
                        (v.item() if hasattr(v, "item") else v)) for k, v in g.attrs.items()}
            tree = node(g["rollup"]["root"]) if "rollup" in g and "root" in g["rollup"] else None
            mem = None
            if "memory" in g:
                m = g["memory"]
                mem = {k: int(v) for k, v in m.attrs.items()}
                if "components_live" in m:
                    mem["components_live"] = {int(r["classTag"]): int(r["count"])
                                              for r in m["components_live"][()]}
            series = None
            if "series" in g:          # `profiler start -perStep`: one row per committed step
                sg = g["series"]
                series = {k: [float(x) for x in sg[k][()].ravel()] for k in sg}
                if "wall_ms" in sg:
                    series["phases"] = [p.decode() if isinstance(p, bytes) else str(p)
                                        for p in sg["wall_ms"].attrs.get("phases", [])]
            out[rid] = {"meta": meta, "tree": tree, "memory": mem, "series": series}
    return out


def sub(tree, *path):
    n = tree
    for p in path:
        if n is None:
            return None
        n = n["children"].get(p)
    return n


def wall(n):
    return 0 if n is None else n["wall_ns"]


def calls(n):
    return 0 if n is None else n["calls"]


def total(n, name, exclude_under=()):
    """Σ wall_ns of every node called `name` in the subtree of n (n itself included), not
    descending into nodes named in exclude_under."""
    if n is None:
        return 0
    return _walk(n, name, set(exclude_under), "wall_ns")


def total_calls(n, name, exclude_under=()):
    if n is None:
        return 0
    return _walk(n, name, set(exclude_under), "calls")


def _walk(n, name, excl, key):
    s = 0
    for k, c in n["children"].items():
        if k == name:
            s += c[key]
        elif k in excl:
            continue
        else:
            s += _walk(c, name, excl, key)
    return s


def top_level_names(tree):
    return sorted(tree["children"]) if tree else []


# --------------------------------------------------------------------------- #
#  The M1 decomposition (232g "Definitions")                                  #
# --------------------------------------------------------------------------- #
PREP_PHASES = ("mpm.a.support", "mpm.a.sp", "mpm.a.churn", "mpm.a.p2g", "mpm.a.repoint",
               "mpm.a.domainChange")
DC_PHASES = ("dc.clearAll", "dc.handle", "dc.numberDOF", "dc.setSize", "dc.integrator",
             "dc.algorithm")
NEWTON_PARTS = ("formTangent", "linearSolve", "formUnbalance", "update", "convTest")
# Top-level scopes a serial mode-A / mode-F step may open. Anything else at the rollup root
# would be a scope opened on another thread (each thread's tree has its own root and the
# report merges them by name): flagged, never silently summed.
EXPECTED_TOP = {"mpm.a.step", "step"}


def m1_step(tree):
    """Decompose one mode-A step (one `mpm analyze 1 dt` profiled run) into the brief's
    quantities. Every number is ns except the dimensionless ones."""
    st = sub(tree, "mpm.a.step")
    if st is None:
        return {"ok": False, "why": "no mpm.a.step scope in this run"}
    an = sub(st, "step")                                   # the analysis's own step
    prep = sub(st, "mpm.a.prepare")
    dc = sub(an, "domainChanged")
    solve = sub(an, "solveCurrentStep")
    sym_all = total(st, "soe.symbolic")
    sym_dc = total(dc, "soe.symbolic")
    sym_solve = total(solve, "soe.symbolic")
    sym_other = sym_all - sym_dc - sym_solve
    t_prep, t_dc = wall(prep), wall(dc)
    t_rb = t_prep + t_dc + (sym_all - sym_dc)               # symbolic counted once
    lin = total(solve, "linearSolve")
    parts = {p: total(solve, p) for p in NEWTON_PARTS}
    parts["linearSolve_numeric"] = lin - sym_solve
    parts["soe.factor"] = total(solve, "soe.factor")
    parts["soe.trisolve"] = total(solve, "soe.trisolve")
    k = total_calls(solve, "linearSolve")
    it_sum = (parts["formTangent"] + parts["linearSolve_numeric"] + parts["formUnbalance"]
              + parts["update"] + parts["convTest"])
    t_it = it_sum / k if k else None
    step_wall = wall(st)
    # the rebuild's sub-phases, no remainder hidden
    phases = {p: wall(sub(prep, p)) for p in PREP_PHASES}
    phases["mpm.a.prepare.other"] = t_prep - sum(phases[p] for p in PREP_PHASES)
    dcp = {p: wall(sub(dc, p)) for p in DC_PHASES}
    setsize = sub(dc, "dc.setSize")
    graph = wall(sub(setsize, "dc.s.graph"))
    sym_in_setsize = total(setsize, "soe.symbolic")
    phases["dc.clearAll"] = dcp["dc.clearAll"]
    phases["dc.handle"] = dcp["dc.handle"]
    phases["dc.numberDOF"] = dcp["dc.numberDOF"]
    phases["dc.s.graph"] = graph
    phases["dc.setSize.soe"] = dcp["dc.setSize"] - graph - sym_in_setsize
    phases["dc.integrator"] = dcp["dc.integrator"]
    phases["dc.algorithm"] = dcp["dc.algorithm"]
    # symbolic run under domainChanged but outside dc.setSize is counted in "soe.symbolic"
    # below, so it is taken out of the remainder here (review, PR #1089)
    phases["dc.other"] = t_dc - sum(dcp.values()) - (sym_dc - sym_in_setsize)
    phases["soe.symbolic"] = sym_all              # wherever it ran (setSize or first solve)
    largest = max(phases, key=lambda p: phases[p]) if phases else None
    commit = sub(an, "commit")
    return {
        "ok": True,
        "t_rb_ns": t_rb, "t_it_ns": t_it, "k": k,
        "R": (t_rb / t_it) if t_it else None,
        "f": (t_rb / step_wall) if step_wall else None,
        "step_wall_ns": step_wall,
        "t_prepare_ns": t_prep, "t_domainChanged_ns": t_dc,
        "symbolic_ns": {"total": sym_all, "in_domainChanged": sym_dc, "in_solve": sym_solve,
                        "elsewhere": sym_other},
        "symbolic_share_of_t_rb": (sym_all / t_rb) if t_rb else None,
        "rebuild_phases_ns": phases, "largest_rebuild_phase": largest,
        "newton_ns": parts, "formTangent_share_of_t_it":
            (parts["formTangent"] / it_sum) if it_sum else None,
        "newStep_ns": wall(sub(an, "newStep")), "commit_ns": wall(commit),
        "g2p_ns": total(commit, "mpm.a.g2p"),
        "unattributed_step_ns": step_wall - t_prep - wall(an),
        # > 1: the step was cut and retried (its attempts are all in these totals); a cut step
        # also finishes its `mpm analyze 1 dt` in several accepted sub-steps (substeps > 1)
        "attempts": calls(prep), "analysis_steps": calls(an), "substeps": calls(st),
        "foreign_top_level": sorted(set(top_level_names(tree)) - EXPECTED_TOP),
    }


def m1_run_summary(steps):
    """One run's numbers over its timed steps: means for t_rb and k, t_it as total iteration
    time over total iterations, R = t_rb / t_it, f = Σ t_rb / Σ step wall."""
    ok = [s for s in steps if s.get("ok")]
    if not ok:
        return None
    n = len(ok)
    t_rb = sum(s["t_rb_ns"] for s in ok) / n
    ksum = sum(s["k"] for s in ok)
    it_tot = sum((s["t_it_ns"] or 0.0) * s["k"] for s in ok)
    t_it = it_tot / ksum if ksum else None
    wall_tot = sum(s["step_wall_ns"] for s in ok)
    phases = {}
    for s in ok:
        for p, v in s["rebuild_phases_ns"].items():
            phases[p] = phases.get(p, 0.0) + v / n
    # formTangent / t_it as a ratio of sums over the timed steps (the document's form), not a
    # mean of per-step ratios (review, PR #1089)
    ft_num = sum(s["newton_ns"]["formTangent"] for s in ok)
    ft_den = sum((s["t_it_ns"] or 0.0) * s["k"] for s in ok)
    # R per rebuild cycle and per accepted sub-step: a cut harness step runs several cycles
    # (one per prepare call) and several accepted sub-steps (review, PR #1089)
    cycles = sum(s.get("attempts") or 1 for s in ok)
    subs = sum(s.get("substeps") or 1 for s in ok)
    t_rb_cycle = sum(s["t_rb_ns"] for s in ok) / cycles
    # the C-cheaper-first row: rebuild on steps whose connectivity changed against steps whose
    # did not (`mpm info connectivityChanges`); a rebuild of an unchanged graph is avoidable
    ch = [s for s in ok if s.get("connectivity_changed")]
    un = [s for s in ok if not s.get("connectivity_changed")]

    def mean(rows, f):
        return sum(f(s) for s in rows) / len(rows) if rows else None
    split = {"n_changed": len(ch), "n_unchanged": len(un),
             "t_rb_changed_ns": mean(ch, lambda s: s["t_rb_ns"]),
             "t_rb_unchanged_ns": mean(un, lambda s: s["t_rb_ns"]),
             "symbolic_changed_ns": mean(ch, lambda s: s["symbolic_ns"]["total"]),
             "symbolic_unchanged_ns": mean(un, lambda s: s["symbolic_ns"]["total"]),
             "rebuild_phases_unchanged_ns": {p: mean(un, lambda s, p=p: s["rebuild_phases_ns"][p])
                                             for p in (un[0]["rebuild_phases_ns"] if un else {})}}
    return {**split, "n_timed_steps": n, "t_rb_ns": t_rb, "t_it_ns": t_it, "k": ksum / n,
            "R": (t_rb / t_it) if t_it else None,
            "f": (sum(s["t_rb_ns"] for s in ok) / wall_tot) if wall_tot else None,
            "step_wall_ns": wall_tot / n,
            "symbolic_share_of_t_rb": (sum(s["symbolic_ns"]["total"] for s in ok)
                                       / sum(s["t_rb_ns"] for s in ok)) if t_rb else None,
            "rebuild_phases_ns": phases,
            "largest_rebuild_phase": max(phases, key=lambda p: phases[p]) if phases else None,
            "formTangent_share_of_t_it": (ft_num / ft_den) if ft_den else None,
            "rebuild_cycles": cycles, "accepted_substeps": subs,
            "t_rb_per_cycle_ns": t_rb_cycle,
            "R_per_cycle": (t_rb_cycle / t_it) if t_it else None,
            "R_per_accepted_substep": ((sum(s["t_rb_ns"] for s in ok) / subs) / t_it)
                                      if t_it else None}
