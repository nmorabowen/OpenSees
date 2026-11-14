"""WP-232f harness: run a tests/mpm_implicit_cases.txt case (mode A) through the shipped binaries.

Written from tests/mpm_implicit_interface.md (the frozen mode-A interface) and the cases
file's own grammar only; it never reads tests/mpm_implicit_reference.py.

One translator turns a case into interpreter commands (`model`, `nDMaterial`, `timeSeries`,
`mpm ...`, and the user's analysis: `constraints`, `numberer`, `system`, `test`, `algorithm`,
`integrator`, `analysis`). The Python engine (opensees.pyd) executes them through `ops`; the
classic Tcl engine (OpenSees.exe) gets a deck whose `puts` lines are parsed back. Both give the
same result dict, and one evaluator turns it into the cases file's observables.

An expected value written `TBD` (232f-V still filling it) is skipped LOUDLY: `evaluate`
reports it as a skip with its label, never as a pass.
"""
from __future__ import annotations

import copy
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "tests" / "mpm_implicit_cases.txt"
HOST = 232            # the MPM host pattern tag
TS_TAG = 9001         # the time series of `series linear`

# keys a protocol case inherits from its `base`
_INHERIT = ("mode", "ndm", "basis", "origin", "h", "cells", "mats", "bodies", "pforces", "fixes",
            "gridvels", "gravity", "series", "dt", "steps", "schedule", "newton", "integrator",
            "velocities", "particles_file", "masscut")


# --------------------------------------------------------------------------- #
#  The cases file                                                              #
# --------------------------------------------------------------------------- #
def _new(name):
    return {"name": name, "kind": "run", "mode": None, "status": "gate", "ndm": 2, "basis": "gimp",
            "origin": None, "h": None, "cells": None, "mats": [], "bodies": [], "pforces": [],
            "fixes": [], "gridvels": [], "gravity": None, "series": None, "dt": None, "steps": None,
            "schedule": [], "newton": None, "integrator": None, "velocities": [],
            "particles_file": None, "masscut": None, "expects": [], "compares": [], "base": None,
            "injects": [], "output": None, "restart": None, "engines": None,
            "integrator_command": None, "analysis_command": None, "own": set()}


def _parse_expect(t, raw):
    e = {"obs": t[1], "args": [], "op": None, "vals": [], "tbd": False, "tol": None,
         "report": t[-1] == "report", "line": raw.strip()}
    toks = t[2:-1] if e["report"] else t[2:]
    i = 0
    while i < len(toks) and toks[i] not in ("eq", "le", "ge"):
        e["args"].append(toks[i])
        i += 1
    if i < len(toks):
        e["op"] = toks[i]
        i += 1
    while i < len(toks) and toks[i] != "tol":
        if toks[i] == "TBD":
            e["tbd"] = True
        else:
            try:
                e["vals"].append(float(toks[i]))
            except ValueError:
                e["args"].append(toks[i])        # e.g. `relerr_ge <case> <field> <r>`
        i += 1
    if i < len(toks):
        e["tol"] = (float(toks[i + 1]), float(toks[i + 2]))
    return e


def parse_cases(path=CASES):
    raw_cases, cur = {}, None
    order = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        t = raw.split()
        if not t or t[0].startswith("#"):
            continue
        k = t[0]
        if k == "case":
            cur = _new(t[1])
            continue
        if cur is None:
            continue
        if k == "end":
            raw_cases[cur["name"]] = cur
            order.append(cur["name"])
            cur = None
            continue
        cur["own"].add(k)
        if k in ("desc", "ref", "gap", "note"):
            pass
        elif k in ("kind", "mode", "status", "basis"):
            cur[k] = t[1]
        elif k == "ndm":
            cur["ndm"] = int(t[1])
        elif k == "origin":
            cur["origin"] = [float(v) for v in t[1:]]
        elif k == "h":
            cur["h"] = float(t[1])
        elif k == "cells":
            cur["cells"] = [int(v) for v in t[1:]]
        elif k == "material":
            cur["mats"].append((int(t[1]), t[2], [float(v) for v in t[3:]]))
        elif k == "body":
            cur["bodies"].append(t[1:])
        elif k == "pforce":
            cur["pforces"].append(t[1:])
        elif k == "fix":
            cur["fixes"].append((t[1], t[2]))
        elif k == "gridvel":
            cur["gridvels"].append(t[1:])
        elif k == "gravity":
            cur["gravity"] = [float(v) for v in t[1:]]
        elif k == "series":
            cur["series"] = float(t[2]) if t[1] == "linear" else None
        elif k == "dt":
            cur["dt"] = float(t[1])
        elif k == "steps":
            cur["steps"] = int(t[1])
        elif k == "schedule":
            cur["schedule"].append((float(t[1]), int(t[2])))
        elif k == "newton":
            cur["newton"] = (t[1], float(t[2]), int(t[3]))
        elif k == "integrator":
            cur["integrator"] = (float(t[2]), float(t[3]))
        elif k == "velocity":
            cur["velocities"].append(t[1:])
        elif k == "particles":
            cur["particles_file"] = t[1]
        elif k == "masscut":
            cur["masscut"] = float(t[1])
        elif k == "expect":
            cur["expects"].append(_parse_expect(t, raw))
        elif k == "compare":
            cur["compares"].append(t[1:])
        elif k == "base":
            cur["base"] = t[1]
        elif k == "inject":
            cur["injects"].append(t[1:])
        elif k == "output":
            cur["output"] = t[1:]
        elif k == "restart":
            cur["restart"] = int(t[2])
        elif k == "engines":
            cur["engines"] = t[1:]
        elif k == "integrator_command":
            cur["integrator_command"] = t[1:]
        elif k == "analysis_command":
            cur["analysis_command"] = t[1:]
        else:
            raise ValueError(f"mpm_implicit_cases.txt: unknown key '{k}' in case {cur['name']}")
    out = []
    for name in order:
        out.append(_resolve(raw_cases, name))
    return out


def _resolve(raw_cases, name):
    c = raw_cases[name]
    if c["base"] is None:
        return c
    b = _resolve(raw_cases, c["base"])
    r = copy.deepcopy(c)
    for k in _INHERIT:
        own_key = {"mats": "material", "bodies": "body", "pforces": "pforce", "fixes": "fix",
                   "gridvels": "gridvel", "velocities": "velocity", "particles_file": "particles",
                   "integrator": "integrator"}.get(k, k)
        if own_key not in c["own"]:
            r[k] = copy.deepcopy(b[k])
    # a protocol that sets `steps` (or `schedule`) replaces the base's stepping
    if "steps" in c["own"] and "schedule" not in c["own"]:
        r["schedule"] = []
    return r


def case_by_name(name):
    for c in parse_cases():
        if c["name"] == name:
            return c
    raise KeyError(name)


def schedule_of(c):
    """[(dt, n)]: the case's schedule lines, or one line of `steps` steps of `dt`."""
    if c["schedule"]:
        return list(c["schedule"])
    return [(c["dt"], c["steps"])]


def num_steps(c):
    return sum(n for _, n in schedule_of(c))


# --------------------------------------------------------------------------- #
#  Contract-1 seeding (the plate with a hole)                                  #
# --------------------------------------------------------------------------- #
def seed_contract1(c, path):
    """Seed the case's bodies by the cases file's rules and write ONE Contract-1 file."""
    import numpy as np
    import mpm_contract_fixtures as fx
    import mpm_contracts as mc
    nd, h = c["ndm"], c["h"]
    xs, vols, lps, bodies, mats, rhos = [], [], [], [], [], []
    for b in c["bodies"]:
        bid, shape = int(b[0]), b[1]
        opts = dict(zip(b[2:][::2], b[2:][1::2]))
        if shape == "holeplate":
            a, hw, nth, nr = float(b[2]), float(b[3]), int(b[4]), int(b[5])
            opts = dict(zip(b[6:][::2], b[6:][1::2]))
            gx, gw = np.polynomial.legendre.leggauss(8)
            dth = (math.pi / 2) / nth
            ds = 1.0 / nr

            def rho_o(th):
                return hw / np.maximum(np.cos(th), np.sin(th))
            for j in range(nr):
                s = (j + 0.5) * ds
                for i in range(nth):
                    th = (i + 0.5) * dth
                    ro = rho_o(th)
                    r = a + s * (ro - a)
                    px, py = r * math.cos(th), r * math.sin(th)
                    xs.append((px, py))
                    tq = th + 0.5 * dth * gx
                    rq = rho_o(tq)
                    rp = a + (s + 0.5 * ds) * (rq - a)
                    rm = a + (s - 0.5 * ds) * (rq - a)
                    V = float(np.sum(gw * 0.5 * (rp ** 2 - rm ** 2)) * 0.5 * dth)
                    vols.append(V)
                    # clipped so no GIMP domain leaves [0, W/2]^2 (C1-POS; cases amendment)
                    l0 = math.sqrt(V) / 2
                    lps.append((min(l0, px, hw - px), min(l0, py, hw - py)))
                    bodies.append(bid)
                    mats.append(int(opts["mat"]))
                    rhos.append(float(opts["rho"]))
        elif shape == "box":
            lo = [float(v) for v in b[2:2 + nd]]
            hi = [float(v) for v in b[2 + nd:2 + 2 * nd]]
            opts = dict(zip(b[2 + 2 * nd:][::2], b[2 + 2 * nd:][1::2]))
            ppc = int(opts["ppc"])
            hs = h / ppc
            nx = [ppc * c["cells"][d] for d in range(nd)]
            for ay in range(nx[1]):
                y = c["origin"][1] + (ay + 0.5) * hs
                if not lo[1] <= y <= hi[1]:
                    continue
                for ax in range(nx[0]):
                    x = c["origin"][0] + (ax + 0.5) * hs
                    if not lo[0] <= x <= hi[0]:
                        continue
                    xs.append((x, y))
                    vols.append(hs ** nd)
                    lps.append((h / (2 * ppc),) * nd)
                    bodies.append(bid)
                    mats.append(int(opts["mat"]))
                    rhos.append(float(opts["rho"]))
        else:
            raise ValueError(shape)
    n = len(xs)
    model = fx.contract1_model()
    model["attrs"].update({"NUM_PARTICLES": np.int64(n), "SET_ID": "6f1c2a3b-4d5e-4f60-8a7b-9c0d1e2f3a4b",
                           "UNIT_LENGTH": "m", "UNIT_FORCE": "N", "UNIT_TIME": "s",
                           "UNIT_SYSTEM": "N_m_s", "MASS_SOURCE": "density",
                           "VOIGT_ORDER": mc.VOIGT[nd], "NDM": np.int32(nd)})
    model["grid"] = {"ORIGIN": np.array(c["origin"], dtype=float),
                     "CELL_SIZE": np.array([h] * nd, dtype=float),
                     "NUM_CELLS": np.array(c["cells"], dtype=np.int32), "BASIS": c["basis"]}
    model["bodies"] = [(int(b[0]), "body" + b[0]) for b in c["bodies"]]
    model.pop("seeding", None)
    model["particles"] = {
        "ID": np.arange(n, dtype="<i8"),
        "POSITION": np.array(xs, dtype="<f8"),
        "VOLUME": np.array(vols, dtype="<f8"),
        "DENSITY": np.array(rhos, dtype="<f8"),
        "GIMP_HALF_LENGTH": np.array(lps, dtype="<f8"),
        "STRESS": np.zeros((n, 3 if nd == 2 else 6)),
        "STRESS33": np.zeros(n),
        "BODY_ID": np.array(bodies, dtype="<i4"),
        "MATERIAL_TAG": np.array(mats, dtype="<i4"),
    }
    fx.write_contract1(path, model)
    return n


# --------------------------------------------------------------------------- #
#  Translation into interpreter commands                                       #
# --------------------------------------------------------------------------- #
def _has_prescribed(c):
    for g in c["gridvels"]:
        vi = g.index("v")
        for tok in g[vi + 1:vi + 1 + c["ndm"]]:
            if tok != "free" and float(tok) != 0.0:
                return True
    return False


def model_commands(c, particle_path=None):
    """[(command, [args])] up to (not including) the analysis."""
    nd = c["ndm"]
    cmds = [("wipe", []), ("model", ["basic", "-ndm", nd, "-ndf", nd])]
    for tag, kind, v in c["mats"]:
        if kind == "elastic":
            cmds.append(("nDMaterial", ["ElasticIsotropic", tag, *v]))
        elif kind == "druckerprager":
            cmds.append(("nDMaterial", ["DruckerPrager", tag, *v]))
        else:
            raise ValueError(kind)
    cmds.append(("mpm", ["grid", HOST, "-origin", *c["origin"], "-h", c["h"], "-cells",
                         *c["cells"], "-basis", c["basis"]]))
    cmds.append(("mpm", ["mode", "implicit" if c["mode"] == "transient" else "quasiStatic"]))
    if c["particles_file"]:
        cmds.append(("mpm", ["particles", "-file", ("PATH", particle_path)]))
    else:
        for b in c["bodies"]:
            bid, shape = int(b[0]), b[1]
            ng = 2 * nd if shape == "box" else nd + 1
            geom = [float(v) for v in b[2:2 + ng]]
            o = dict(zip(b[2 + ng:][::2], b[2 + ng:][1::2]))
            cmds.append(("mpm", ["body", bid, shape, *geom, "-ppc", int(o["ppc"]),
                                 "-mat", int(o["mat"]), "-rho", float(o["rho"])]))
    for v in c["velocities"]:
        if v[1] == "uniform":
            cmds.append(("mpm", ["velocity", int(v[0]), "uniform", *[float(w) for w in v[2:2 + nd]]]))
        else:
            cmds.append(("mpm", ["velocity", int(v[0]), "sine", v[2], v[3], float(v[4]),
                                 float(v[5]), float(v[6])]))
    for p in c["pforces"]:
        cmds.append(("mpm", ["pforce", int(p[0]), "box", *[float(w) for w in p[2:2 + 3 * nd]]]))
    for face, comps in c["fixes"]:
        cmds.append(("mpm", ["fix", face, comps]))
    for g in c["gridvels"]:
        name = g[0]
        lo = [float(w) for w in g[2:2 + nd]]
        hi = [float(w) for w in g[2 + nd:2 + 2 * nd]]
        vi = g.index("v")
        vals = [w if w == "free" else float(w) for w in g[vi + 1:vi + 1 + nd]]
        rest = g[vi + 1 + nd:]
        args = ["gridvel", name, "box", *lo, *hi, "v", *vals]
        i = 0
        while i < len(rest):
            if rest[i] == "start":
                args += ["start", float(rest[i + 1])]
            elif rest[i] == "track":
                args += ["track", rest[i + 1]]
            i += 2
        cmds.append(("mpm", args))
    if c["gravity"] is not None:
        cmds.append(("mpm", ["gravity", *c["gravity"]]))
    if c["series"] is not None:
        Tr = c["series"]
        cmds.append(("timeSeries", ["Path", TS_TAG, "-time", (0.0, Tr, 1.0e9),
                                    "-values", (0.0, 1.0, 1.0)]))
        cmds.append(("mpm", ["timeSeries", TS_TAG]))
    if c["masscut"] is not None:
        cmds.append(("mpm", ["massCutoff", c["masscut"]]))
    # data runs only (never set in CI): the empty-node treatment and cutoff of every case
    en = os.environ.get("LADRUNO_MPM_EMPTY_NODES")
    if en:
        mode, _, g = en.partition(":")
        cmds.append(("mpm", ["emptyNodes", mode] + ([float(g)] if g else [])))
    mc = os.environ.get("LADRUNO_MPM_MASSCUT")
    if mc and c["masscut"] is None:
        cmds.append(("mpm", ["massCutoff", float(mc)]))
    return cmds


def analysis_commands(c, override_integrator=True):
    test, tol, it = c["newton"] or ("NormDispIncr", 1e-10, 30)
    cmds = [("constraints", ["Transformation" if _has_prescribed(c) else "Plain"]),
            ("numberer", ["RCM"]), ("system", ["UmfPack"]), ("test", [test, tol, it]),
            ("algorithm", ["Newton"])]
    if c["integrator_command"] and override_integrator:
        integ = c["integrator_command"]
        cmds.append(("integrator", [_num(x) for x in integ]))
    elif c["mode"] == "transient":
        g, b = c["integrator"] or (0.5, 0.25)
        cmds.append(("integrator", ["Newmark", g, b]))
    else:
        cmds.append(("integrator", ["LoadControl", schedule_of(c)[0][0]]))
    if c["analysis_command"]:
        cmds.append(("analysis", list(c["analysis_command"])))
    else:
        cmds.append(("analysis", ["Transient" if c["mode"] == "transient" else "Static"]))
    return cmds


def _num(x):
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


def _tcl_arg(a):
    if isinstance(a, tuple):
        if a and a[0] == "PATH":
            return "{" + str(a[1]).replace("\\", "/") + "}"
        return "{" + " ".join(_tcl_arg(x) for x in a) + "}"
    if isinstance(a, float):
        return repr(a)
    return str(a)


def _py_args(args):
    out = []
    for a in args:
        if isinstance(a, tuple):
            if a and a[0] == "PATH":
                out.append(str(a[1]))
            else:
                out.extend(a)
        else:
            out.append(a)
    return out


# --------------------------------------------------------------------------- #
#  What a run records                                                          #
# --------------------------------------------------------------------------- #
FIELDS_END = ("x", "v", "a", "F", "V", "sig", "sigHat", "epsHat")
FIELDS_INIT = ("x", "m", "id", "body", "mat", "V")


def tip_args(c):
    return sorted({float(e["args"][0]) for e in c["expects"]
                   if e["obs"].startswith(("tip_disp_y", "first_peak_tip", "first_period_tip"))}
                  | {float(cp[1]) for cp in c["compares"] if cp[0].startswith("tip_disp_y")})


def gridvel_names(c):
    return [g[0] for g in c["gridvels"]]


def _tip(nd, X, x, m, xf):
    sel = [p for p in range(len(m)) if X[nd * p] >= xf - 1e-12]
    mt = sum(m[p] for p in sel)
    return sum(m[p] * (x[nd * p + 1] - X[nd * p + 1]) for p in sel) / mt if mt > 0 else 0.0


def _inject_commands(c):
    out = []
    for inj in c["injects"]:
        if inj[0] == "fail_attempt":
            o = dict(zip(inj[1:][::2], inj[1:][1::2]))
            out.append(("mpm", ["debug", "failAttempt", int(o["step"]), int(o.get("attempts", 1))]))
        elif inj[0] == "force_domainchange_every_step":
            out.append(("mpm", ["debug", "reseed", "always"]))
        elif inj[0] == "skip_domainchange_without_crossing":
            out.append(("mpm", ["debug", "reseed", "skip"]))
        else:
            raise ValueError("unknown inject " + inj[0])
    return out


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def run_python(c, ops, workdir=None, mass_cutoff=None, empty_nodes=None):
    """Run a case through opensees.pyd. Returns the result dict."""
    td = tempfile.TemporaryDirectory(ignore_cleanup_errors=True) if workdir is None else None
    wd = Path(workdir) if workdir else Path(td.name)
    try:
        ppath = None
        if c["particles_file"]:
            ppath = wd / c["particles_file"]
            seed_contract1(c, ppath)
        for name, args in model_commands(c, ppath):
            getattr(ops, name)(*_py_args(args))
        if mass_cutoff is not None:
            ops.mpm("massCutoff", mass_cutoff)
        if empty_nodes is not None:
            ops.mpm("emptyNodes", *empty_nodes)
        for name, args in _inject_commands(c):
            getattr(ops, name)(*_py_args(args))
        r = _record_init(c, lambda f: _as_list(ops.mpm("particleData", f)))
        refusal_pre = c["kind"] == "refusal" and c["integrator_command"] and \
            c["integrator_command"][0] == "DisplacementControl"
        if refusal_pre:
            # one admitted step first, so a grid node exists to name
            for name, args in analysis_commands(c, override_integrator=False):
                getattr(ops, name)(*_py_args(args))
            ops.mpm("analyze", 1, schedule_of(c)[0][0])
            tag = int(_as_list(ops.mpm("info", "gridNodeTags"))[0])
            integ = [tag if x == "<grid" else x for x in c["integrator_command"] if x != "node>"]
            ops.integrator(*[_num(x) if not isinstance(x, int) else x for x in integ])
        else:
            for name, args in analysis_commands(c):
                getattr(ops, name)(*_py_args(args))
        out = None
        if c["output"]:
            o = c["output"]
            every = int(o[o.index("every") + 1])
            flds = [f.upper() for f in o[o.index("fields") + 1:]]
            out = wd / (c["name"] + ".ladruno")
            ops.mpm("output", str(out), "-every", every, "-fields", *flds)
        r["x_before"] = _as_list(ops.mpm("particleData", "x"))
        steps_cb = _make_recorder(c, r, lambda *a: ops.mpm(*a), sync=bool(out))
        steps_cb(initial=True)
        done = 0
        restart_at = c["restart"]
        for dt, n in schedule_of(c):
            per_step = not c["schedule"]
            calls = [(1, dt)] * n if per_step else [(n, dt)]
            for k, (m, h) in enumerate(calls):
                try:
                    ops.mpm("analyze", m, h)
                except Exception:
                    r["stopped"] = done
                    break
                done += m
                steps_cb()
                if restart_at is not None and done == restart_at and r.get("restart_mode") is None:
                    r["restart_mode"] = "saved"
                    ops.database("File", str(wd / "mpmdb"))
                    ops.save(restart_at)
                    # restore into a fresh model in THIS process, rebuild the analysis
                    ops.wipe()
                    ops.model("basic", "-ndm", c["ndm"], "-ndf", c["ndm"])
                    ops.database("File", str(wd / "mpmdb"))
                    ops.restore(restart_at)
                    for name, args in analysis_commands(c):
                        getattr(ops, name)(*_py_args(args))
            if r.get("stopped") is not None:
                break
        if out is not None:
            ops.mpm("output", "-close")
            r["output_path"] = str(out)
            r["sync_ok"] = _check_sync(r, out)
        for f in FIELDS_END:
            r["end_" + f] = _as_list(ops.mpm("particleData", f))
        r["x_after"] = r["end_x"]
        for w in ("g2p", "acceptedSteps", "retries", "domainChanges", "connectivityChanges"):
            r["info_" + w] = int(ops.mpm("info", w))
        r["dtHistory"] = _as_list(ops.mpm("info", "dtHistory"))
        return r
    finally:
        if td is not None:
            td.cleanup()


def _record_init(c, get):
    r = {"engine": "python", "stopped": None, "hist_t": [], "hist_tip": {}, "hist_cf": {},
         "summary": [], "sync": []}
    for f in FIELDS_INIT:
        r["init_" + f] = get(f)
    return r


def _make_recorder(c, r, mpm, sync=False):
    nd = c["ndm"]
    tips = tip_args(c)
    names = gridvel_names(c)
    X, m = r["init_x"], r["init_m"]

    def cb(initial=False):
        t = float(mpm("info", "time"))
        x = _as_list(mpm("particleData", "x"))
        r["hist_t"].append(t)
        for xf in tips:
            r["hist_tip"].setdefault(xf, []).append(_tip(nd, X, x, m, xf))
        for nm in names:
            r["hist_cf"].setdefault(nm, []).append(_as_list(mpm("condForce", nm)))
        r["summary"].append(_as_list(mpm("info", "summary")))
        if sync:
            r["sync"].append((t, x, _as_list(mpm("particleData", "sig"))))
    return cb


def _check_sync(r, path):
    """Every output record's positions and stresses equal the in-memory state of the same
    accepted step, exactly (matched by TIME)."""
    import numpy as np
    import h5py
    ids = [int(round(v)) for v in r["init_id"]]
    order = np.argsort(ids)
    with h5py.File(path, "r") as f:
        stages = [s for s in f if s.startswith("MODEL_STAGE[")]
        op = f[stages[0] + "/RESULTS/ON_PARTICLES"]
        off, times = op["OFFSETS"][...], op["TIME"][...]
        pos, sig = op["POSITION/DATA"], op["STRESS/DATA"]
        mem = {t: (np.array(x), np.array(s)) for t, x, s in r["sync"]}
        if len(times) < 2:
            return False
        for k, t in enumerate(times):
            if float(t) not in mem:
                return False
            x, s = mem[float(t)]
            a, b = off[k], off[k + 1]
            if not np.array_equal(pos[a:b], x.reshape(-1, 2)[order]):
                return False
            if not np.array_equal(sig[a:b], s.reshape(-1, 3)[order]):
                return False
    return True


# --------------------------------------------------------------------------- #
#  Classic Tcl                                                                 #
# --------------------------------------------------------------------------- #
def tcl_deck(c, particle_path=None, phase="all", db=None):
    """phase: all | save (run to the restart step, save) | restore (restore, continue)."""
    nd = c["ndm"]
    lines = ["proc dump {tag f} { puts \"$tag $f [mpm particleData $f]\" }"]
    tips = tip_args(c)
    names = gridvel_names(c)
    if phase != "restore":
        for name, args in model_commands(c, particle_path):
            lines.append(name + " " + " ".join(_tcl_arg(a) for a in args))
        for name, args in _inject_commands(c):
            lines.append(name + " " + " ".join(_tcl_arg(a) for a in args))
        for f in FIELDS_INIT:
            lines.append("dump INIT %s" % f)
    else:
        lines.append("model basic -ndm %d -ndf %d" % (nd, nd))
        lines.append("database File {%s}" % str(db).replace("\\", "/"))
        lines.append("restore %d" % c["restart"])
    dc = c["kind"] == "refusal" and c["integrator_command"] and         c["integrator_command"][0] == "DisplacementControl"
    lines.append("set ok [catch {")
    for name, args in analysis_commands(c, override_integrator=not dc):
        lines.append("  " + name + " " + " ".join(_tcl_arg(a) for a in args))
    lines.append("}]")
    if dc:      # one admitted step first, so a grid node exists to name
        lines.append("mpm analyze 1 %r" % schedule_of(c)[0][0])
        lines.append("set gt [lindex [mpm info gridNodeTags] 0]")
        rest = [x for x in c["integrator_command"][1:] if x not in ("<grid", "node>")]
        lines.append("integrator DisplacementControl [expr {int($gt)}] %s" % " ".join(rest))
        lines.append("puts \"XB [mpm particleData x]\"")
    rec = ["puts \"T [mpm info time]\"", "puts \"SUM [mpm info summary]\""]
    rec.append("puts \"X [mpm particleData x]\"" if tips else "")
    for nm in names:
        rec.append("puts \"CF %s [mpm condForce %s]\"" % (nm, nm))
    rec = [x for x in rec if x]
    if phase != "restore":
        lines += rec
    if phase == "all":
        sched = schedule_of(c)
        per_step = not c["schedule"]
        for dt, n in sched:
            if per_step:
                lines.append("for {set s 0} {$s < %d} {incr s} {" % n)
                lines.append("  if {[catch {mpm analyze 1 %r}]} { puts \"STOP\"; break }" % dt)
                lines += ["  " + x for x in rec]
                lines.append("}")
            else:
                lines.append("if {[catch {mpm analyze %d %r}]} { puts \"STOP\" }" % (n, dt))
                lines += rec
    elif phase in ("save", "inproc"):
        dt = schedule_of(c)[0][0]
        lines.append("for {set s 0} {$s < %d} {incr s} {" % c["restart"])
        lines.append("  if {[catch {mpm analyze 1 %r}]} { puts \"STOP\"; break }" % dt)
        lines.append("}")
        lines.append("database File {%s}" % str(db).replace("\\", "/"))
        lines.append("save %d" % c["restart"])
        if phase == "inproc":          # wipe + restore in the same OpenSees.exe
            lines.append("wipe")
            lines.append("model basic -ndm %d -ndf %d" % (nd, nd))
            lines.append("database File {%s}" % str(db).replace("\\", "/"))
            lines.append("restore %d" % c["restart"])
            for name, args in analysis_commands(c):
                lines.append(name + " " + " ".join(_tcl_arg(a) for a in args))
            lines.append("for {set s %d} {$s < %d} {incr s} {" % (c["restart"], num_steps(c)))
            lines.append("  if {[catch {mpm analyze 1 %r}]} { puts \"STOP\"; break }" % dt)
            lines.append("}")
    else:
        dt = schedule_of(c)[0][0]
        lines.append("for {set s %d} {$s < %d} {incr s} {" % (c["restart"], num_steps(c)))
        lines.append("  if {[catch {mpm analyze 1 %r}]} { puts \"STOP\"; break }" % dt)
        lines.append("}")
    if phase in ("all", "restore", "inproc"):
        for f in FIELDS_END:
            lines.append("dump END %s" % f)
        for w in ("g2p", "acceptedSteps", "retries"):
            lines.append("puts \"INFO %s [mpm info %s]\"" % (w, w))
        lines.append("puts \"DTH [mpm info dtHistory]\"")
    lines.append("puts \"DECK-DONE\"")
    return "\n".join(lines) + "\n"


def _run_deck(text, exe, cwd, name, timeout):
    deck = Path(cwd) / (name + ".tcl")
    deck.write_text(text, encoding="utf-8")
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    proc = subprocess.run([exe, str(deck)], capture_output=True, text=True, timeout=timeout,
                          env=env, cwd=cwd, stdin=subprocess.DEVNULL)
    out = proc.stdout + "\n" + proc.stderr
    if "DECK-DONE" not in out:
        raise RuntimeError(f"Tcl deck {name} did not finish (rc={proc.returncode}):\n" + out[-4000:])
    return out


def run_tcl(c, exe, timeout=3600, inproc=False):
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        ppath = None
        if c["particles_file"]:
            ppath = Path(td) / c["particles_file"]
            seed_contract1(c, ppath)
        if c["restart"] is not None and inproc:
            db = Path(td) / "mpmdb"
            return _parse_tcl(c, _run_deck(tcl_deck(c, ppath, "inproc", db), exe, td,
                                           c["name"] + "_inproc", timeout))
        if c["restart"] is not None:
            db = Path(td) / "mpmdb"
            out1 = _run_deck(tcl_deck(c, ppath, "save", db), exe, td, c["name"] + "_save", timeout)
            out = _run_deck(tcl_deck(c, ppath, "restore", db), exe, td, c["name"] + "_restore", timeout)
            init = _parse_tcl(c, out1)
            r = _parse_tcl(c, out)
            for f in FIELDS_INIT:
                r["init_" + f] = init.get("init_" + f, [])
            return r
        out = _run_deck(tcl_deck(c, ppath), exe, td, c["name"], timeout)
        return _parse_tcl(c, out)


def _parse_tcl(c, out):
    nd = c["ndm"]
    r = {"engine": "tcl", "stopped": None, "hist_t": [], "hist_tip": {}, "hist_cf": {},
         "summary": [], "log": out}
    xs = []
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "T":
            r["hist_t"].append(float(t[1]))
        elif t[0] == "SUM":
            r["summary"].append([float(v) for v in t[1:]])
        elif t[0] == "X":
            xs.append([float(v) for v in t[1:]])
        elif t[0] == "CF":
            r["hist_cf"].setdefault(t[1], []).append([float(v) for v in t[2:]])
        elif t[0] == "STOP":
            r["stopped"] = len(r["hist_t"]) - 1
        elif t[0] in ("INIT", "END") and len(t) >= 2:
            r[t[0].lower() + "_" + t[1]] = [float(v) for v in t[2:]]
        elif t[0] == "INFO":
            r["info_" + t[1]] = int(float(t[2]))
        elif t[0] == "XB":
            r["x_before_dc"] = [float(v) for v in t[1:]]
        elif t[0] == "DTH":
            r["dtHistory"] = [float(v) for v in t[1:]]
    if "init_x" in r:
        for xf in tip_args(c):
            r["hist_tip"][xf] = [_tip(nd, r["init_x"], x, r["init_m"], xf) for x in xs]
        r["x_before"] = r.get("x_before_dc", r["init_x"])
        r["x_after"] = r.get("end_x", [])
    return r


# --------------------------------------------------------------------------- #
#  Observables and expectations                                                #
# --------------------------------------------------------------------------- #
def _minima(t, d):
    """Local minima of the history, each refined by the parabola through its neighbours."""
    out = []
    for k in range(1, len(d) - 1):
        if d[k] < d[k - 1] and d[k] <= d[k + 1]:
            t0, t1, t2 = t[k - 1], t[k], t[k + 1]
            y0, y1, y2 = d[k - 1], d[k], d[k + 1]
            den = (t0 - t1) * (t0 - t2) * (t1 - t2)
            A = (t2 * (y1 - y0) + t1 * (y0 - y2) + t0 * (y2 - y1)) / den
            B = (t2 * t2 * (y0 - y1) + t1 * t1 * (y2 - y0) + t0 * t0 * (y1 - y2)) / den
            C = (t1 * t2 * (t1 - t2) * y0 + t2 * t0 * (t2 - t0) * y1 + t0 * t1 * (t0 - t1) * y2) / den
            if A > 0:
                tv = -B / (2 * A)
                out.append((tv, C - B * B / (4 * A)))
            else:
                out.append((t1, y1))
    return out


def field_relerr(ra, rb, field):
    a, b = ra["end_" + field], rb["end_" + field]
    if len(a) != len(b) or not b:
        return float("inf")
    den = max(abs(v) for v in b) or 1.0
    return max(abs(x - y) for x, y in zip(a, b)) / den


PROTOCOL_FIELDS = ("x", "v", "F", "sig")


def observable(c, r, e, other=None):
    """The value list of one expectation (`other`: the run of a case it names)."""
    nd = c["ndm"]
    nv = 6 if nd == 3 else 3
    ob, args = e["obs"], e["args"]
    if ob == "tip_disp_y":
        return [r["hist_tip"][float(args[0])][-1]]
    if ob in ("first_peak_tip_disp_y", "first_period_tip"):
        mins = _minima(r["hist_t"], r["hist_tip"][float(args[0])])
        if ob == "first_peak_tip_disp_y":
            return [mins[0][1]] if mins else [float("nan")]
        return [mins[1][0] - mins[0][0]] if len(mins) >= 2 else [float("nan")]
    if ob in ("cond_pressure_max", "cond_pressure_plateau"):
        name, comp, width = args[0], "xyz".index(args[1]), float(args[2])
        q = [-f[comp] / width for f in r["hist_cf"][name]]
        if ob == "cond_pressure_max":
            return [max(q)]
        tf = float(args[3])
        sel = [qq for tt, qq in zip(r["hist_t"], q) if tt >= tf - 1e-9]
        return [sum(sel) / len(sel)] if sel else [float("nan")]
    if ob == "peak_sthth_over_snet":
        xc, yc, lig = (float(a) for a in args[:3])
        snet = r["hist_cf"]["load"][-1][1] / lig
        x, s = r["end_x"], r["end_sig"]
        best = -float("inf")
        for p in range(len(x) // nd):
            th = math.atan2(x[nd * p + 1] - yc, x[nd * p] - xc)
            sx, sy, sxy = s[nv * p], s[nv * p + 1], s[nv * p + 2]
            stt = sx * math.sin(th) ** 2 + sy * math.cos(th) ** 2 - 2 * sxy * math.sin(th) * math.cos(th)
            best = max(best, stt / snet)
        return [best]
    if ob == "seed_count":
        return [float(len(r["init_V"]))]
    if ob == "seed_volume_sum":
        return [sum(r["init_V"])]
    if ob == "seed_moment":
        V, X = r["init_V"], r["init_x"]
        return [sum(V[p] * X[nd * p + d] for p in range(len(V))) for d in range(nd)]
    if ob == "selfweight_error":
        rho0, g, l0, yb = (float(a) for a in args)
        X, sig = r["init_x"], r["end_sig"]
        res = [sig[nv * p + 1] - (-rho0 * g * (l0 - (X[nd * p + 1] - yb))) for p in range(len(X) // nd)]
        return [max(abs(v) for v in res) / (rho0 * g * l0)]
    if ob == "period_comx":
        t = r["hist_t"]
        cv = [s[1] / s[0] for s in r["summary"]]
        cr = []
        for k in range(len(cv) - 1):
            if cv[k] < 0.0 <= cv[k + 1]:
                cr.append(t[k] + (t[k + 1] - t[k]) * (-cv[k]) / (cv[k + 1] - cv[k]))
        return [(cr[-1] - cr[0]) / (len(cr) - 1) if len(cr) >= 2 else float("nan")]
    if ob == "g2p_count":
        return [float(r["info_g2p"])]
    if ob == "accepted_steps":
        return [float(r["info_acceptedSteps"])]
    if ob == "relerr_ge":
        return [field_relerr(r, other, args[1])]
    if ob == "recorder_sync":
        return [1.0 if r.get("sync_ok") else 0.0]
    if ob == "refused":
        unchanged = r.get("x_before") == r.get("x_after")
        return [1.0 if (r["stopped"] is not None and unchanged) else 0.0]
    raise KeyError("unknown observable " + ob)


def evaluate(c, r, runs=None):
    """[(label, verdict, detail)] with verdict 'pass' | 'fail' | 'tbd' | 'report'.

    `runs(name)` returns the run of another case (relerr_ge, compare lines)."""
    res = []
    for e in c["expects"]:
        label = f"{c['name']} {e['obs']} {' '.join(e['args'])}".strip()
        other = None
        if e["obs"] == "relerr_ge":
            other = runs(e["args"][0])
        if e["obs"] == "g2p_count":
            got, want = observable(c, r, e), [float(r["info_acceptedSteps"])]
            ok = got == want
            res.append((label, "pass" if ok else "fail", f"g2p {got} accepted {want}"))
            continue
        if e["tbd"]:
            res.append((label, "tbd", "expected value TBD (232f-V): skipped, not passed"))
            continue
        got = observable(c, r, e, other)
        if e["obs"] in ("refused", "recorder_sync"):
            ok = got[0] == 1.0
        elif e["obs"] == "relerr_ge":
            ok = got[0] >= float(e["args"][2])
        elif e["op"] == "eq":
            at, rt = e["tol"] or (0.0, 0.0)
            ok = len(got) == len(e["vals"]) and all(abs(g - x) <= at + rt * abs(x)
                                                    for g, x in zip(got, e["vals"]))
        elif e["op"] == "le":
            ok = got[0] <= e["vals"][0]
        else:
            ok = got[0] >= e["vals"][0]
        if c["kind"] == "run" and r.get("stopped") is not None:
            ok = False
        verdict = "pass" if ok else "fail"
        if e["report"] or c["status"] == "report":
            verdict = "report"
        res.append((label, verdict, f"got {got} {e['op'] or ''} {e['vals']}"))
    for cp in c["compares"]:
        if cp[0] == "particle_fields":
            other = runs(cp[2])
            tol = float(cp[cp.index("tol") + 2]) if "tol" in cp else 0.0
            worst = max(field_relerr(r, other, f) for f in PROTOCOL_FIELDS)
            res.append((f"{c['name']} particle_fields vs {cp[2]}", "pass" if worst <= tol else "fail",
                        f"worst field relerr {worst:.3g} (tol {tol})"))
        elif "relerrle" in cp:
            i = cp.index("relerrle")
            oc = case_by_name(cp[i + 1])
            e = {"obs": cp[0], "args": cp[1:i]}

            def own(cc, rr):
                ex = [x for x in cc["expects"] if x["obs"] == e["obs"] and x["args"] == e["args"]][0]
                if ex["tbd"]:
                    return None
                v = ex["vals"][0]
                return abs(observable(cc, rr, e)[0] - v) / abs(v)
            a = own(c, r)
            b = own(oc, runs(cp[i + 1]))
            if a is None or b is None:
                res.append((f"{c['name']} compare {' '.join(cp)}", "tbd", "reference value TBD"))
                continue
            res.append((f"{c['name']} compare {' '.join(cp)}", "pass" if a <= b else "fail",
                        f"relerr {a:.4g} vs {cp[i + 1]} {b:.4g}"))
        elif "errle" in cp:
            i = cp.index("errle")
            if cp[-1] == "TBD":
                res.append((f"{c['name']} compare {' '.join(cp)}", "tbd", "reference value TBD"))
                continue
            val = float(cp[-1])
            e = {"obs": cp[0], "args": cp[1:i]}
            a = observable(c, r, e)[0]
            oc = case_by_name(cp[i + 1])
            b = observable(oc, runs(cp[i + 1]), e)[0]
            ok = abs(a - val) <= abs(b - val)
            res.append((f"{c['name']} compare {' '.join(cp)}", "pass" if ok else "fail",
                        f"|{a}-{val}| vs |{b}-{val}|"))
    return res


# --------------------------------------------------------------------------- #
#  A fresh second process (the restart test's other form)                      #
# --------------------------------------------------------------------------- #
_RESTORE_SCRIPT = r'''
import json, sys
sys.path.insert(0, {tests!r})
import ladruno_mpm_implicit_harness as H
from _testbed import ops
c = H.case_by_name({name!r})
ops.wipe()
ops.model("basic", "-ndm", c["ndm"], "-ndf", c["ndm"])
ops.database("File", {db!r})
ops.restore(c["restart"])
for n, a in H.analysis_commands(c):
    getattr(ops, n)(*H._py_args(a))
dt = H.schedule_of(c)[0][0]
for k in range(c["restart"], H.num_steps(c)):
    ops.mpm("analyze", 1, dt)
out = {{f: H._as_list(ops.mpm("particleData", f)) for f in H.FIELDS_END}}
print("JSON " + json.dumps(out))
'''


def run_python_two_process(c, ops):
    """Save at the restart step in this process, restore and continue in a fresh one."""
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        db = str(Path(td) / "mpmdb")
        for name, args in model_commands(c):
            getattr(ops, name)(*_py_args(args))
        for name, args in analysis_commands(c):
            getattr(ops, name)(*_py_args(args))
        dt = schedule_of(c)[0][0]
        for _ in range(c["restart"]):
            ops.mpm("analyze", 1, dt)
        ops.database("File", db)
        ops.save(c["restart"])
        script = _RESTORE_SCRIPT.format(tests=str(REPO / "tests"), name=c["name"], db=db)
        proc = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                              timeout=3600, env=dict(os.environ), cwd=str(REPO / "tests"))
        line = [x for x in proc.stdout.splitlines() if x.startswith("JSON ")]
        if not line:
            raise RuntimeError("the restore process failed:\n" + (proc.stdout + proc.stderr)[-4000:])
        data = json.loads(line[-1][5:])
        return {"engine": "python-2proc", **{"end_" + k: v for k, v in data.items()}}
