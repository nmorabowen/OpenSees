"""WP-232j harness: run a tests/mpm_antilocking_cases.txt case through the shipped binaries.

Written from tests/mpm_antilocking_interface.md (frozen as J1-J13 with the "232j" freeze
conditions) and the cases files' own grammar; it never reads tests/mpm_reference.py or
tests/mpm_implicit_reference.py.

The cases file extends tests/mpm_implicit_cases.txt's grammar (its translator,
tests/ladruno_mpm_implicit_harness.py, builds the model and the analysis) with:
  - `antilocking off|fbar` (default fbar in this file, freeze Q2; inherited through `base`);
  - `ladder <alg> <tol> <it> ... relaxed <k> halve [<h_min>]` (K6; h_min default 5e-4, freeze Q1);
  - `material <tag> ladrunoj2 <K> <G> <s0>` (L3: `nDMaterial LadrunoJ2 tag K G -iso voce s0 0 0 0`);
  - kinds `affine` (L1), `tangent` (L2), `mms` (L5);
  - observables cond_pressure_window / _argmax, reach_t, the op `out`, mms_p / mms_nonincr,
    tangent_rel_undeformed / tangent_ratio_strained.
Every footing run reports freeze F1's anti-locking counts: averaged and fallback particle-steps,
and the fallback share in the footing's plastic zone (x <= 3B, y >= -2B; B = 2 m).
"""
from __future__ import annotations

import copy
import math
import os
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_implicit_harness as H  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
CASES = REPO / "tests" / "mpm_antilocking_cases.txt"
EXTRA = ("antilocking", "ladder", "affine", "row", "tangent_eps", "stress", "tangent_strain",
         "energy_mode", "dt_factor", "sample", "v0", "ledger_case")
H_MIN_DEFAULT = 5e-4          # freeze Q1: the FE's dsc < 1e-6 m at v = 0.002 per unit pseudo-time
PLASTIC_ZONE = (6.0, -4.0)    # freeze F1: x <= 3B, y >= -2B with B = 2 m


# --------------------------------------------------------------------------- #
#  The cases file                                                              #
# --------------------------------------------------------------------------- #
def _parse_expect(t, raw):
    if "out" in t:
        i = t.index("out")
        return {"obs": t[1], "args": t[2:i], "op": "out", "vals": [float(t[i + 1]), float(t[i + 2])],
                "tol": None, "report": t[-1] == "report", "tbd": False, "line": raw.strip()}
    return H._parse_expect(t, raw)


def parse_cases(path=CASES):
    """The cases, resolved: the implicit grammar plus this file's keys (see the docstring)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    kept, extras, cur, order = [], {}, None, []
    for raw in lines:
        t = raw.split()
        if t and t[0] == "case":
            cur = t[1]
            order.append(cur)
            extras[cur] = {"antilocking": None, "ladder": None, "affine": [], "row": None,
                           "tangent_eps": None, "expects": [], "base": None, "j2": [],
                           "stress": [], "tangent_strain": None, "energy_mode": None,
                           "dt_factor": None, "sample": None, "v0": None, "ledger_case": None}
            kept.append(raw)
            continue
        if cur is not None and t:
            k = t[0]
            if k in EXTRA:
                if k in ("affine", "stress"):
                    extras[cur][k].append(t[1:])
                else:
                    extras[cur][k] = t[1:]
                continue
            if k == "expect":
                extras[cur]["expects"].append(_parse_expect(t, raw))
                continue
            if k == "material" and t[2] == "ladrunoj2":
                extras[cur]["j2"].append((int(t[1]), [float(v) for v in t[3:]]))
                continue
            if k == "base":
                extras[cur]["base"] = t[1]
            if k == "end":
                cur = None
        kept.append(raw)
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "cases.txt"
        p.write_text("\n".join(kept) + "\n", encoding="utf-8")
        base_cases = {c["name"]: c for c in H.parse_cases(p)}

    def ext(name, key):
        e = extras[name]
        if key == "j2" and e["j2"]:
            return e["j2"]
        if key != "j2" and e[key] is not None:
            return e[key]
        return ext(e["base"], key) if e["base"] else ([] if key == "j2" else None)

    out = []
    for name in order:
        c = base_cases[name]
        e = extras[name]
        c["expects"] = e["expects"]
        al = ext(name, "antilocking")
        c["antilocking"] = (al[0] if al else "fbar")
        c["ladder"] = _ladder(ext(name, "ladder"))
        c["affine"] = e["affine"]
        c["row"] = e["row"][0] if e["row"] else None
        te = ext(name, "tangent_eps")
        c["tangent_eps"] = [float(v) for v in te] if te else None
        c["j2"] = ext(name, "j2")
        c["stress"] = e["stress"]
        for k in ("tangent_strain", "dt_factor", "v0"):
            c[k] = float(e[k][0]) if e[k] else None
        c["sample"] = int(e["sample"][0]) if e["sample"] else None
        c["energy_mode"] = e["energy_mode"][0] if e["energy_mode"] else None
        c["ledger_case"] = e["ledger_case"][0] if e["ledger_case"] else None
        out.append(c)
    return out


def _ladder(t):
    if not t:
        return None
    tries, i, relaxed, hmin = [], 0, None, None
    while i < len(t):
        if t[i] == "relaxed":
            relaxed = int(t[i + 1])
            i += 2
        elif t[i] == "halve":
            if i + 1 < len(t):
                hmin = float(t[i + 1])
                i += 2
            else:
                i += 1
        else:
            tries.append((t[i], float(t[i + 1]), int(t[i + 2])))
            i += 3
    return {"tries": tries, "relaxed": relaxed or len(tries) + 1,
            "h_min": hmin if hmin is not None else H_MIN_DEFAULT}


def case_by_name(name):
    for c in parse_cases():
        if c["name"] == name:
            return c
    raise KeyError(name)


# --------------------------------------------------------------------------- #
#  Model commands                                                              #
# --------------------------------------------------------------------------- #
def model_commands(c):
    """The implicit translator's commands, LadrunoJ2 materials inserted, then the setting."""
    cc = copy.deepcopy(c)
    cmds = H.model_commands(cc)
    if c.get("mode") == "explicit":           # the implicit translator knows modes A only
        cmds = [("mpm", ["mode", "explicit"]) if (n == "mpm" and a and a[0] == "mode") else (n, a)
                for n, a in cmds]
    j2 = [("nDMaterial", ["LadrunoJ2", tag, v[0], v[1], "-iso", "voce", v[2], 0.0, 0.0, 0.0])
          for tag, v in c.get("j2") or []]
    if j2:
        k = [i for i, (n, _) in enumerate(cmds) if n == "model"][0]
        cmds[k + 1:k + 1] = j2
    for st in c.get("stress") or []:
        cmds.append(("mpm", ["stress", int(st[0]), st[1], *[float(v) for v in st[2:]]]))
    cmds.append(("mpm", ["antiLocking", c.get("antilocking", "fbar")]))
    return cmds


def _apply(ops, cmds):
    for name, args in cmds:
        getattr(ops, name)(*H._py_args(args))


def _as_list(v):
    return H._as_list(v)


def ladruno_build(ops):
    try:
        return str(ops.ladrunoBuild())
    except Exception as exc:  # pragma: no cover
        return f"unknown ({exc})"


# --------------------------------------------------------------------------- #
#  L1: affine fields                                                           #
# --------------------------------------------------------------------------- #
def run_affine(c, ops):
    """Every `affine` line of the case through `mpm debug affineCheck`. Returns {name: out}."""
    _apply(ops, model_commands(c))
    nd = c["ndm"]
    N = int(ops.mpm("info", "numParticles"))
    res = {"N": N}
    for line in c["affine"]:
        name = line[0]
        ia, ib = line.index("A"), line.index("b")
        A = [float(v) for v in line[ia + 1:ia + 1 + nd * nd]]
        b = [float(v) for v in line[ib + 1:ib + 1 + nd]]
        gates = dict(zip(line[ib + 1 + nd:][::2], [float(v) for v in line[ib + 1 + nd:][1::2]]))
        out = _as_list(ops.mpm("debug", "affineCheck", *A, *b))
        res[name] = {"out": out, "gates": gates}
    return res


def evaluate_affine(c, res):
    rows = []
    N = res["N"]
    s0 = max([abs(float(v)) for st in c.get("stress") or [] for v in st[2:]] or [0.0])
    for name, r in res.items():
        if name == "N":
            continue
        rel, deps, navg, nfb, nbad, spur = r["out"]
        ok = navg == N and nfb == 0 and nbad == 0
        det = (f"rel {rel:.3e}, deps {deps:.3e}, spurious stress {spur:.3e} (|sigma_0| {s0:.3g}), "
               f"averaged {navg:.0f}/{N}, fallback {nfb:.0f}, refused {nbad:.0f}")
        if "rel_le" in r["gates"]:
            ok = ok and rel <= r["gates"]["rel_le"]
        if "deps_le" in r["gates"]:
            ok = ok and deps <= r["gates"]["deps_le"]
        if "spur_rel_le" in r["gates"]:
            ok = ok and s0 > 0.0 and spur <= r["gates"]["spur_rel_le"] * s0
        rows.append((f"{c['name']} {name}", "pass" if ok else "fail", det))
    return rows


# --------------------------------------------------------------------------- #
#  L2: the tangent                                                             #
# --------------------------------------------------------------------------- #
def run_tangent(c, ops):
    _apply(ops, model_commands(c))
    _apply(ops, H.analysis_commands(c))
    eps = c["tangent_eps"] or [1e-8, 1e-9, 1e-10]

    def tc(e, al=None):
        args = ["debug", "tangentCheck", 1e-6, e] + (["-antiLocking", al] if al else [])
        return _as_list(ops.mpm(*args))

    r = {"undeformed_on": [tc(e) for e in eps], "undeformed_off": [tc(e, "off") for e in eps]}
    for _ in range(H.num_steps(c)):
        ops.mpm("analyze", 1, c["dt"])
    r["strain_max"] = max(abs(v) for v in _as_list(ops.mpm("particleData", "epsHat")))
    r["N"] = int(ops.mpm("info", "numParticles"))
    r["strained_on"] = [tc(e) for e in eps]
    r["strained_off"] = [tc(e, "off") for e in eps]
    r["eps"] = eps
    return r


def _converged(rows):
    """The first eps whose max rel agrees with the next to 10 %, else the first."""
    for a, b in zip(rows, rows[1:]):
        if abs(a[0] - b[0]) <= 0.1 * max(a[0], 1e-300) + 1e-12:
            return a[0]
    return rows[0][0]


# --------------------------------------------------------------------------- #
#  The footing: K6's ladder per accepted step, F1's counts                     #
# --------------------------------------------------------------------------- #
def run_ladder(c, ops, log=None, dump=None):
    """A footing case through `mpm analyze 1 h` with K6's ladder (interface §5): per step, each
    try sets `algorithm` and `test NormUnbalance`; `mpm cut 2 0` so a failed try reverts and the
    harness moves on; all tries failing halves h (regrown x2 per accepted sub-step, to dt); h
    below h_min stops the run. Records at the step ends only."""
    t0w = time.time()
    _apply(ops, model_commands(c))
    r = H._record_init(c, lambda f: _as_list(ops.mpm("particleData", f)))
    _apply(ops, H.analysis_commands(c))
    ops.mpm("cut", 2.0, 0)
    if dump:                                  # the L4 diagnosis (`mpm debug failDump`)
        ops.mpm("debug", "failDump", str(dump[0]), int(dump[1]))
    lad = c["ladder"] or {"tries": [("Newton", c["newton"][1], c["newton"][2])], "relaxed": 99,
                          "h_min": H_MIN_DEFAULT}
    rec = H._make_recorder(c, r, lambda *a: ops.mpm(*a))
    rec(initial=True)
    dt, n = c["dt"], H.num_steps(c)
    nd = c["ndm"]
    t = float(ops.mpm("info", "time"))
    t_end = t + n * dt
    h = dt
    failed = relaxed = halvings = 0
    al_avg = al_fb = zone_avg = zone_fb = 0
    r["stopped"] = None
    k = 1
    while k <= n:
        target = t_end - (n - k) * dt
        x_prev = _as_list(ops.mpm("particleData", "x"))
        hh = min(h, target - t)
        ok_try = None
        for i, (alg, tol, it) in enumerate(lad["tries"]):
            ops.algorithm(alg)
            ops.test("NormUnbalance", tol, it)
            try:
                with H._stderr_capture():
                    ops.mpm("analyze", 1, hh)
                ok_try = i
                break
            except Exception:
                failed += 1
        if ok_try is None:
            h = hh / 2.0
            halvings += 1
            if h < lad["h_min"]:
                r["stopped"] = k - 1
                if log:
                    log(f"{c['name']}: stopped at t = {t:.6g} (h {h:.3g} < {lad['h_min']:.3g})")
                break
            continue
        if ok_try + 1 >= lad["relaxed"]:
            relaxed += 1
        flags = _as_list(ops.mpm("particleData", "antiLocking"))
        for p, f in enumerate(flags):
            if f == 0:
                continue
            inzone = x_prev[nd * p] <= PLASTIC_ZONE[0] + 1e-12 and \
                x_prev[nd * p + 1] >= PLASTIC_ZONE[1] - 1e-12
            if f == 1:
                al_avg += 1
                zone_avg += inzone
            else:
                al_fb += 1
                zone_fb += inzone
        t = float(ops.mpm("info", "time"))
        h = min(2.0 * hh, dt) if hh < dt else dt
        if abs(t - target) <= 1e-9 * max(1.0, abs(target)):
            t = target
            rec()
            if log and k % 10 == 0:
                q = -r["hist_cf"]["footing"][-1][1] if "footing" in r["hist_cf"] else float("nan")
                log(f"{c['name']}: t = {t:g}, q = {q:.2f}, failed {failed}, relaxed {relaxed}, "
                    f"halvings {halvings}")
            k += 1
    if dump:
        ops.mpm("debug", "failDump", "off")
    counts = _as_list(ops.mpm("info", "antiLockingCounts"))
    r.update({"failed_attempts": failed, "relaxed_steps": relaxed, "halvings": halvings,
              "al_averaged": al_avg, "al_fallback": al_fb, "al_counts_engine": counts,
              "zone_averaged": zone_avg, "zone_fallback": zone_fb,
              "zone_fallback_share": zone_fb / (zone_avg + zone_fb) if zone_avg + zone_fb else 0.0,
              "wall_s": time.time() - t0w, "antilocking": c["antilocking"],
              "ladrunoBuild": ladruno_build(ops)})
    return r


def footing_series(r, name="footing", comp=1, width=1.0):
    ts = r["hist_t"][1:]
    qs = [-v[comp] / width for v in r["hist_cf"][name][1:]]
    return ts, qs


def observable(c, r, e):
    ts, qs = footing_series(r)
    o = e["obs"]
    if o == "cond_pressure_max":
        return max(qs) if qs else float("nan")
    if o == "cond_pressure_argmax":
        return ts[qs.index(max(qs))] if qs else float("nan")
    if o == "cond_pressure_plateau":
        t_from = float(e["args"][3])
        v = [q for t, q in zip(ts, qs) if t >= t_from - 1e-9]
        return sum(v) / len(v) if v else float("nan")
    if o == "cond_pressure_window":
        a, b = float(e["args"][3]), float(e["args"][4])
        v = [q for t, q in zip(ts, qs) if a - 1e-9 <= t <= b + 1e-9]
        return sum(v) / len(v) if v else float("nan")
    if o == "reach_t":
        return r["hist_t"][-1] if r["hist_t"] else float("nan")
    raise ValueError("unknown observable " + o)


def judge(e, got):
    if e["op"] is None:
        return "report", f"{got:.6g}"
    if not (isinstance(got, (int, float)) and math.isfinite(got)):
        return "fail", f"{got} (not a finite number)"
    if e["op"] == "out":
        lo, hi = e["vals"]
        ok = not (lo <= got <= hi)
        return ("pass" if ok else "fail"), f"{got:.6g} outside [{lo}, {hi}]"
    exp = e["vals"][0]
    if e["op"] == "eq":
        a, rr = e["tol"] if e["tol"] else (0.0, 0.0)
        ok = abs(got - exp) <= a + rr * abs(exp)
    elif e["op"] == "le":
        ok = got <= exp
    else:
        ok = got >= exp
    return ("pass" if ok else "fail"), f"{got:.6g} {e['op']} {exp}"


def evaluate(c, r):
    rows = []
    for e in c["expects"]:
        got = observable(c, r, e)
        v, det = judge(e, got)
        if e["report"] or c.get("status") == "report":
            v = "report"
        rows.append((e["line"], v, det))
    return rows


def run_record(c, r):
    """The JSON record of one footing run (tests/data/mpm/antilocking/)."""
    ts, qs = footing_series(r)
    rec = {k: r.get(k) for k in ("ladrunoBuild", "antilocking", "failed_attempts", "relaxed_steps",
                                 "halvings", "al_averaged", "al_fallback", "al_counts_engine",
                                 "zone_averaged", "zone_fallback", "zone_fallback_share",
                                 "wall_s", "stopped")}
    rec.update({"case": c["name"], "t": ts, "q": qs,
                "s": [0.002 * max(t - 10.0, 0.0) for t in ts]})
    for e in c["expects"]:
        rec.setdefault("observables", {})[e["line"]] = observable(c, r, e)
    rec["verdicts"] = [list(x) for x in evaluate(c, r)]
    return rec


# --------------------------------------------------------------------------- #
#  L9: energy                                                                  #
# --------------------------------------------------------------------------- #
def _energy_run(c, ops, al):
    cc = dict(c)
    cc["antilocking"] = al
    cmds = [x for x in model_commands(cc) if not (x[0] == "mpm" and x[1] and x[1][0] == "mode")]
    _apply(ops, cmds)
    if c["energy_mode"] == "A":
        ops.mpm("mode", "implicit")
    if c["v0"]:
        ops.mpm("velocity", 1, "uniform", c["v0"], *([0.0] * (c["ndm"] - 1)))
    dt = c["dt_factor"] * float(ops.mpm("criticalTimeStep"))
    if c["energy_mode"] == "A":
        _apply(ops, [("constraints", ["Plain"]), ("numberer", ["RCM"]), ("system", ["UmfPack"]),
                     ("test", ["NormDispIncr", 1e-10, 30]), ("algorithm", ["Newton"]),
                     ("integrator", ["Newmark", 0.5, 0.25]), ("analysis", ["Transient"])])
    ke = [_as_list(ops.mpm("info", "summary"))[-1]]
    done, n, blk = 0, c["steps"], c["sample"] or 20
    try:
        with H._stderr_capture():
            while done < n:
                m = min(blk, n - done)
                ops.mpm("analyze", m, dt)
                done += m
                ke.append(_as_list(ops.mpm("info", "summary"))[-1])
    except Exception:
        pass
    return {"ke": ke, "done": done, "dt": dt, "ke_max": max(ke)}


def run_energy(c, ops):
    """L9: the case off and on; ke_ratio = KE_max(on) / KE_max(off); completed = both ran."""
    t0 = time.time()
    off = _energy_run(c, ops, "off")
    on = _energy_run(c, ops, "fbar")
    return {"off": off, "on": on, "ke_ratio": on["ke_max"] / max(off["ke_max"], 1e-300),
            "completed": 1.0 if (off["done"] == c["steps"] and on["done"] == c["steps"]) else 0.0,
            "wall_s": time.time() - t0}


def run_ledger(c, ops):
    """L9: an S1 energy-ledger case of tests/mpm_cases.txt with anti-locking on (slice-1 harness)."""
    import ladruno_mpm_harness as H1
    base = [x for x in H1.parse_cases() if x["name"] == c["ledger_case"]][0]
    orig = H1.setup_commands
    H1.setup_commands = lambda cc: orig(cc) + [("mpm", ["antiLocking", "fbar"])]
    try:
        r = H1.run_python(base, ops)
    finally:
        H1.setup_commands = orig
    obs = H1.observables(base, r)
    return {"obs": obs, "stopped": r.get("stopped"),
            "counts": _as_list(ops.mpm("info", "antiLockingCounts"))}


def evaluate_values(c, values):
    """[(line, verdict, detail)] for expects over a dict obs -> value."""
    rows = []
    for e in c["expects"]:
        v, det = judge(e, values.get(e["obs"], float("nan")))
        rows.append((e["line"], "report" if e["report"] else v, det))
    return rows


def flushprint(*a):
    print(*a, flush=True)


def finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)
