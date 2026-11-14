"""WP-232c harness: run a tests/mpm_cases.txt case through the shipped binaries.

One translator turns a case into a list of interpreter commands (`model`,
`nDMaterial`, `timeSeries`, `mpm ...`). The Python engine (opensees.pyd)
executes them through `ops`; the classic Tcl engine (OpenSees.exe) gets them as a
deck whose `puts` lines this module parses back. Both produce the same result
dict, and one evaluator turns it into the interface's §10 observables, so the
two engines are held to the same expectations and to each other (ADR 232 §8,
"engine parity").

Written from tests/mpm_interface.md and tests/mpm_cases.txt only; it never
reads tests/mpm_reference.py.
"""
from __future__ import annotations

import math
import os
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "tests" / "mpm_cases.txt"
TS_TAG = 9001          # the timeSeries tag a case with `series linear` uses
HOST_TAG = 232         # the MPM host pattern tag


# --------------------------------------------------------------------------- #
#  The case file                                                               #
# --------------------------------------------------------------------------- #
def parse_cases(path=CASES):
    cases, cur = [], None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        t = raw.split()
        if not t or t[0].startswith("#"):
            continue
        k = t[0]
        if k == "case":
            cur = {"name": t[1], "ndm": 2, "basis": "gimp", "mats": {}, "bodies": [],
                   "velocities": [], "stresses": [], "particles": [], "fixes": [], "walls": [],
                   "gravity": None, "series": None, "damping": 0.0, "flip": 1.0,
                   "masscut": 1e-10, "safety": 0.5, "dt": 0.0, "steps": 0, "tol": (0.0, 0.0),
                   "expects": [], "compares": [], "fields": [], "kind": None}
            continue
        if cur is None:
            continue
        if k == "end":
            cases.append(cur)
            cur = None
        elif k == "kind":
            cur["kind"] = t[1]
        elif k in ("desc", "note"):
            pass
        elif k == "ndm":
            cur["ndm"] = int(t[1])
        elif k == "basis":
            cur["basis"] = t[1]
        elif k == "origin":
            cur["origin"] = [float(v) for v in t[1:]]
        elif k == "h":
            cur["h"] = float(t[1])
        elif k == "cells":
            cur["cells"] = [int(v) for v in t[1:]]
        elif k == "material":
            cur["mats"][int(t[1])] = (float(t[3]), float(t[4]))
        elif k == "body":
            cur["bodies"].append(t)
        elif k == "velocity":
            cur["velocities"].append(t)
        elif k == "stress":
            cur["stresses"].append(t)
        elif k == "particle":
            cur["particles"].append(t)
        elif k == "point":
            nd = cur["ndm"]
            cur["point"] = [float(v) for v in t[1:1 + nd]]
            cur["plp"] = [float(v) for v in t[2 + nd:2 + 2 * nd]]
        elif k == "fix":
            cur["fixes"].append((t[1], t[2]))
        elif k == "wall":
            cur["walls"].append((t[1], float(t[2])))
        elif k == "gravity":
            cur["gravity"] = [float(v) for v in t[1:]]
        elif k == "series":
            cur["series"] = float(t[2]) if t[1] == "linear" else None
        elif k in ("damping", "flip", "masscut", "safety", "dt"):
            cur[k] = float(t[1])
        elif k == "steps":
            cur["steps"] = int(t[1])
        elif k == "tol":
            cur["tol"] = (float(t[1]), float(t[2]))
        elif k == "compare":
            cur["compares"].append(t)
        elif k == "field":
            cur["fields"].append(t)
        elif k == "expect":
            e = {"obs": t[1], "args": [], "op": None, "vals": [], "pending": False,
                 "fields": False, "tol": None, "line": raw}
            i = 2
            while i < len(t) and t[i] not in ("eq", "le", "ge"):
                e["args"].append(t[i])
                i += 1
            if i < len(t):
                e["op"] = t[i]
                i += 1
            while i < len(t) and t[i] != "tol":
                if t[i] == "pending-oracle":
                    e["pending"] = True
                elif t[i] == "fields":
                    e["fields"] = True
                else:
                    e["vals"].append(float(t[i]))
                i += 1
            if i < len(t):
                e["tol"] = (float(t[i + 1]), float(t[i + 2]))
            cur["expects"].append(e)
    return cases


def case_by_name(name):
    for c in parse_cases():
        if c["name"] == name:
            return c
    raise KeyError(name)


# --------------------------------------------------------------------------- #
#  Translation into interpreter commands                                       #
# --------------------------------------------------------------------------- #
def setup_commands(c):
    """[(command, [args...])]; a tuple arg is a Tcl list / flattened in Python."""
    nd = c["ndm"]
    cmds = [("wipe", []), ("model", ["basic", "-ndm", nd, "-ndf", nd])]
    for tag, (E, nu) in sorted(c["mats"].items()):
        cmds.append(("nDMaterial", ["ElasticIsotropic", tag, E, nu]))
    if c["kind"] == "shape1d":
        cmds.append(("mpm", ["grid", HOST_TAG, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 1, 1]))
        return cmds
    grid = ["grid", HOST_TAG, "-origin", *c["origin"], "-h", c["h"], "-cells", *c["cells"],
            "-basis", c["basis"]]
    cmds.append(("mpm", grid))
    for b in c["bodies"]:
        bid, shape = int(b[1]), b[2]
        ng = 2 * nd if shape == "box" else nd + 1
        geom = [float(v) for v in b[3:3 + ng]]
        rest = b[3 + ng:]
        ppc, mat, rho = int(rest[1]), int(rest[3]), float(rest[5])
        cmds.append(("mpm", ["body", bid, shape, *geom, "-ppc", ppc, "-mat", mat, "-rho", rho]))
    nv = 6 if nd == 3 else 3
    for p in c["particles"]:
        pid, body, mat, rho, V0 = int(p[1]), int(p[2]), int(p[3]), float(p[4]), float(p[5])
        i = 7
        x = [float(v) for v in p[i:i + nd]]; i += nd + 1
        lp = [float(v) for v in p[i:i + nd]]; i += nd + 1
        v = [float(w) for w in p[i:i + nd]]; i += nd + 1
        s = [float(w) for w in p[i:i + nv]]
        cmds.append(("mpm", ["particle", pid, body, mat, rho, V0, "-x", *x, "-lp", *lp,
                             "-v", *v, "-s", *s]))
    for vl in c["velocities"]:
        if vl[2] == "uniform":
            cmds.append(("mpm", ["velocity", int(vl[1]), "uniform", *[float(w) for w in vl[3:3 + nd]]]))
        else:
            cmds.append(("mpm", ["velocity", int(vl[1]), "sine", vl[3], vl[4],
                                 float(vl[5]), float(vl[6]), float(vl[7])]))
    for st in c["stresses"]:
        cmds.append(("mpm", ["stress", int(st[1]), "geostatic", *[float(w) for w in st[3:7]]]))
    for face, comps in c["fixes"]:
        cmds.append(("mpm", ["fix", face, comps]))
    for face, mu in c["walls"]:
        cmds.append(("mpm", ["wall", face, mu]))
    if c["gravity"] is not None:
        cmds.append(("mpm", ["gravity", *c["gravity"]]))
    if c["series"] is not None:
        Tr = c["series"]
        cmds.append(("timeSeries", ["Path", TS_TAG, "-time", (0.0, Tr, 1.0e9),
                                    "-values", (0.0, 1.0, 1.0)]))
        cmds.append(("mpm", ["timeSeries", TS_TAG]))
    cmds.append(("mpm", ["damping", "local", c["damping"]]))
    cmds.append(("mpm", ["flip", c["flip"]]))
    cmds.append(("mpm", ["massCutoff", c["masscut"]]))
    cmds.append(("mpm", ["safety", c["safety"]]))
    return cmds


def _tcl_arg(a):
    if isinstance(a, tuple):
        return "{" + " ".join(_tcl_arg(x) for x in a) + "}"
    if isinstance(a, float):
        return repr(a)
    return str(a)


def _py_args(args):
    out = []
    for a in args:
        if isinstance(a, tuple):
            out.extend(a)
        else:
            out.append(a)
    return out


FIELDS_END = ("x", "v", "F", "V", "sig", "sigHat", "epsHat")
FIELDS_INIT = ("x", "m", "id", "body", "mat", "V")


# --------------------------------------------------------------------------- #
#  Running                                                                     #
# --------------------------------------------------------------------------- #
def run_python(c, ops):
    """Run a case through opensees.pyd. Returns the result dict."""
    for name, args in setup_commands(c):
        getattr(ops, name)(*_py_args(args))
    r = {"engine": "python"}
    kind = c["kind"]
    if kind == "shape1d":
        r["s1"] = [list(ops.mpm("debug", "shape1d", *[float(a) for a in e["args"][:3]]))
                   for e in c["expects"]]
        return r
    if kind == "shape":
        r["shape"] = list(ops.mpm("debug", "shape", *c["point"], "-lp", *c["plp"]))
        return r
    if kind == "critdt":
        r["dt_crit"] = float(ops.mpm("criticalTimeStep"))
        return r
    if kind == "masscut":
        r["active"] = [int(round(v)) for v in _as_list(ops.mpm("debug", "activeNodes"))]
        return r
    for f in FIELDS_INIT:
        r["init_" + f] = _as_list(ops.mpm("particleData", f))
    r["summary"] = [_as_list(ops.mpm("info", "summary"))]
    r["stopped"] = None
    for s in range(c["steps"]):
        try:
            ops.mpm("analyze", 1, c["dt"])
        except Exception:
            r["stopped"] = s
            break
        r["summary"].append(_as_list(ops.mpm("info", "summary")))
    for f in FIELDS_END:
        r["end_" + f] = _as_list(ops.mpm("particleData", f))
    r["time"] = float(ops.mpm("info", "time"))
    return r


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [float(x) for x in v]
    return [float(v)]


def tcl_deck(c):
    lines = ["proc dump {tag f} { puts \"$tag $f [mpm particleData $f]\" }"]
    for name, args in setup_commands(c):
        lines.append(name + " " + " ".join(_tcl_arg(a) for a in args))
    kind = c["kind"]
    if kind == "shape1d":
        for e in c["expects"]:
            lines.append("puts \"S1 [mpm debug shape1d %s]\"" % " ".join(e["args"][:3]))
    elif kind == "shape":
        lines.append("puts \"SHAPE [mpm debug shape %s -lp %s]\"" % (
            " ".join(repr(v) for v in c["point"]), " ".join(repr(v) for v in c["plp"])))
    elif kind == "critdt":
        lines.append("puts \"DTC [mpm criticalTimeStep]\"")
    elif kind == "masscut":
        lines.append("puts \"ACTIVE [mpm debug activeNodes]\"")
    else:
        for f in FIELDS_INIT:
            lines.append("dump INIT %s" % f)
        lines.append("puts \"SUM [mpm info summary]\"")
        lines.append("for {set s 0} {$s < %d} {incr s} {" % c["steps"])
        lines.append("  if {[catch {mpm analyze 1 %r}]} { puts \"STOP $s\"; break }" % c["dt"])
        lines.append("  puts \"SUM [mpm info summary]\"")
        lines.append("}")
        for f in FIELDS_END:
            lines.append("dump END %s" % f)
        lines.append("puts \"TIME [mpm info time]\"")
    lines.append("puts \"DECK-DONE\"")
    return "\n".join(lines) + "\n"


def run_tcl(c, exe, timeout=900):
    """Run a case through classic Tcl (OpenSees.exe). Returns the result dict."""
    with tempfile.TemporaryDirectory() as td:
        deck = Path(td) / (c["name"] + ".tcl")
        deck.write_text(tcl_deck(c), encoding="utf-8")
        env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
        proc = subprocess.run([exe, str(deck)], capture_output=True, text=True, timeout=timeout,
                              env=env, cwd=td, stdin=subprocess.DEVNULL)
    out = proc.stdout + "\n" + proc.stderr
    if "DECK-DONE" not in out:
        raise RuntimeError(f"Tcl deck for {c['name']} did not finish (rc={proc.returncode}):\n"
                           + out[-4000:])
    r = {"engine": "tcl", "summary": [], "stopped": None, "s1": [], "log": out}
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "S1":
            r["s1"].append([float(v) for v in t[1:]])
        elif t[0] == "SHAPE":
            r["shape"] = [float(v) for v in t[1:]]
        elif t[0] == "DTC":
            r["dt_crit"] = float(t[1])
        elif t[0] == "ACTIVE":
            r["active"] = [int(round(float(v))) for v in t[1:]]
        elif t[0] == "SUM":
            r["summary"].append([float(v) for v in t[1:]])
        elif t[0] == "STOP":
            r["stopped"] = int(t[1])
        elif t[0] in ("INIT", "END") and len(t) >= 2:
            r[t[0].lower() + "_" + t[1]] = [float(v) for v in t[2:]]
        elif t[0] == "TIME":
            r["time"] = float(t[1])
    return r


# --------------------------------------------------------------------------- #
#  Observables (interface §10) and expectations                                #
# --------------------------------------------------------------------------- #
def _lame(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return lam, mu


def strain_energy(c, r, prefix="end_"):
    nd = c["ndm"]
    nv = 6 if nd == 3 else 3
    eps = r.get(prefix + "epsHat")
    if eps is None:
        return 0.0
    V = r[prefix + "V"]
    mats = r.get("init_mat")
    se = 0.0
    n = len(V)
    # every case's particles share one material unless they declare more
    tags = sorted(c["mats"])
    for p in range(n):
        E, nu = c["mats"][tags[0] if len(tags) == 1 else int(mats[p])]
        lam, mu = _lame(E, nu)
        e = eps[nv * p:nv * p + nv]
        if nd == 3:
            tr = e[0] + e[1] + e[2]
            s = [lam * tr + 2 * mu * e[k] for k in range(3)] + [mu * e[k] for k in range(3, 6)]
        else:
            tr = e[0] + e[1]
            s = [lam * tr + 2 * mu * e[0], lam * tr + 2 * mu * e[1], mu * e[2]]
        se += 0.5 * V[p] * sum(a * b for a, b in zip(e, s))
    return se


def observables(c, r):
    """{(obs, args...): [values]} for every expectation of the case."""
    nd = c["ndm"]
    nv = 6 if nd == 3 else 3
    out = {}
    kind = c["kind"]
    for e in c["expects"]:
        key = (e["obs"], *e["args"])
        ob = e["obs"]
        if kind == "shape1d":
            idx = [x for x in c["expects"]].index(e)
            out[("s1", idx)] = r["s1"][idx]
            continue
        if kind == "shape":
            sh, stride = r["shape"], 2 + nd
            got = [0.0] * (1 + nd)
            for k in range(0, len(sh), stride):
                if int(round(sh[k])) == int(e["args"][0]):
                    got = sh[k + 1:k + stride]
            out[key] = got
            continue
        if ob == "dt_crit":
            out[key] = [r["dt_crit"]]
        elif ob == "active_nodes":
            out[key] = [float(v) for v in r["active"]]
        elif ob in ("refused", "accepted"):
            out[key] = [1.0 if (r["stopped"] == 0 if ob == "refused" else r["stopped"] is None) else 0.0]
        elif ob == "particle_fields":
            out[key] = r
        elif ob == "translation_error":
            v0 = [float(a) for a in e["args"]]
            X, x, t = r["init_x"], r["end_x"], r["time"]
            n = len(X) // nd
            out[key] = [max(math.sqrt(sum((x[nd * p + d] - (X[nd * p + d] + v0[d] * t)) ** 2
                                          for d in range(nd))) for p in range(n))]
        elif ob == "F_error":
            F = r["end_F"]
            out[key] = [max(abs(F[k] - (1.0 if (k % (nd * nd)) // nd == (k % nd) else 0.0))
                            for k in range(len(F)))]
        elif ob == "max_abs_stress":
            out[key] = [max(abs(s) for s in r["end_sig"])]
        elif ob == "period_comx":
            dt = c["dt"]
            cv = [s[1] / s[0] for s in r["summary"]]
            cr = []
            for k in range(len(cv) - 1):
                if cv[k] < 0.0 <= cv[k + 1]:
                    cr.append(k * dt + dt * (-cv[k]) / (cv[k + 1] - cv[k]))
            out[key] = [(cr[-1] - cr[0]) / (len(cr) - 1) if len(cr) >= 2 else float("nan")]
        elif ob == "energy_ratio_end":
            ke0 = r["summary"][0][-1]
            keN = r["summary"][-1][-1]
            out[key] = [(keN + strain_energy(c, r)) / ke0]
        elif ob in ("selfweight_error", "cell_crossing_meter"):
            rho0, g, l0, yb = (float(a) for a in e["args"])
            X, sig = r["init_x"], r["end_sig"]
            n = len(X) // nd
            res = [sig[nv * p + 1] - (-rho0 * g * (l0 - (X[nd * p + 1] - yb))) for p in range(n)]
            if ob == "selfweight_error":
                out[key] = [max(abs(v) for v in res) / (rho0 * g * l0)]
            else:
                out[key] = [math.sqrt(sum(v * v for v in res) / n) / (rho0 * g * l0)]
        elif ob == "max_speed":
            out[key] = [max(s[1 + nd] for s in r["summary"])]
        elif ob == "momentum":
            out[key] = [max(abs(s[1 + d]) for s in r["summary"]) for d in range(nd)]
        elif ob in ("com_disp_x", "com_v_x"):
            b = int(e["args"][0])
            body, m = r["init_body"], r["init_m"]
            n = len(m)
            sel = [p for p in range(n) if int(round(body[p])) == b]
            mt = sum(m[p] for p in sel)
            if ob == "com_v_x":
                out[key] = [sum(m[p] * r["end_v"][nd * p] for p in sel) / mt]
            else:
                out[key] = [sum(m[p] * (r["end_x"][nd * p] - r["init_x"][nd * p]) for p in sel) / mt]
        else:
            raise KeyError("unknown observable " + ob)
    return out


def check_fields(c, r, tol):
    """Particle fields after the last step against the case's `field` lines."""
    nd = c["ndm"]
    nv = 6 if nd == 3 else 3
    nn = nd * nd
    ids = [int(round(v)) for v in r["init_id"]]
    bad = []
    for fl in c["fields"]:
        pid = int(fl[1])
        p = ids.index(pid)
        got = (r["end_x"][nd * p:nd * p + nd] + r["end_v"][nd * p:nd * p + nd] +
               r["end_F"][nn * p:nn * p + nn] + [r["end_V"][p]] +
               r["end_sig"][nv * p:nv * p + nv] + r["end_sigHat"][nv * p:nv * p + nv] +
               r["end_epsHat"][nv * p:nv * p + nv])
        exp = []
        for tok in fl[2:]:
            try:
                exp.append(float(tok))
            except ValueError:
                pass
        if len(exp) != len(got):
            bad.append(f"particle {pid}: {len(got)} values, {len(exp)} expected")
            continue
        for i, (g, x) in enumerate(zip(got, exp)):
            if abs(g - x) > tol[0] + tol[1] * abs(x):
                bad.append(f"particle {pid} value {i}: got {g!r}, expected {x!r}")
    return bad


def evaluate(c, r):
    """[(label, ok, detail)] for every expectation of the case."""
    obs = observables(c, r)
    res = []
    for idx, e in enumerate(c["expects"]):
        label = f"{c['name']} {e['obs']} {' '.join(e['args'])}".strip()
        if e["pending"]:
            continue
        if c["kind"] == "shape1d":
            got = obs[("s1", idx)]
        else:
            got = obs[(e["obs"], *e["args"])]
        if e["obs"] == "particle_fields":
            tol = e["tol"] or c["tol"]
            bad = check_fields(c, r, tol)
            res.append((label, not bad, "; ".join(bad[:5])))
            continue
        if e["obs"] in ("refused", "accepted"):
            res.append((label, got[0] == 1.0, f"stopped at {r['stopped']}"))
            continue
        at, rt = e["tol"] or c["tol"]
        if e["op"] == "eq":
            ok = len(got) == len(e["vals"]) and all(
                abs(g - x) <= at + rt * abs(x) for g, x in zip(got, e["vals"]))
        elif e["op"] == "le":
            ok = got[0] <= e["vals"][0]
        else:
            ok = got[0] >= e["vals"][0]
        if r.get("stopped") is not None and c["kind"] == "run" and e["obs"] != "refused":
            ok = False
        res.append((label, ok, f"got {got} {e['op']} {e['vals']}"))
    return res


def scalar(c, r, obs, args):
    return observables(c, r)[(obs, *args)][0]
